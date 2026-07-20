//! winding 기반 단일 arrangement 엔진 (M5) — 설계 `docs/winding-engine.md`.
//!
//! 현 case-routed 엔진(detector + 생존표)을 대체할 문헌식 winding 엔진을 **격리·미배선**으로
//! 증분 구축한다(cold-start 금지). 이 모듈의 코드는 커토버(M-E)까지 프로덕션에서 안 쓰이며,
//! 매 마일스톤 격리 테스트로 검증한다. 모든 판정은 exact plane-triple substrate(좌표 미읽음).
//!
//! **M-A:** per-face propagation-winding **분류** 프리미티브 — 원본 정점 anchor
//! (`point_in_solid_idx`) + seam 교차에서 alternation 전파(`run_classes`)로, 횡단 절단된 면의
//! 각 boundary run이 keep되는지 판정. 순수 횡단·hole-free 면만.
//!
//! **M-B(현재):** **unified per-face 분류** — 글로벌 detector 없이 로컬 집합-멤버십
//! (`other_classes.contains(&pi)`)으로 면당 공면-vs-횡단을 가른다. 공면이면 on±
//! (`coincident_overlap_faces`+`same_normal`+`coplanar_survival`), 횡단이면 M-A. 분류만
//! (sub-face 방출·조립은 M-C+). 이게 case-routed detector 폭발을 죽이는 winding 엔진의 핵심.
#![cfg_attr(not(test), allow(dead_code))]

use super::*;

/// per-run 분류 결과: `(seam 교차 triple들, ∂f run별 원본정점, run별 kept 여부)`.
type FaceRunClasses = (Vec<[usize; 3]>, Vec<Vec<usize>>, Vec<bool>);

/// 면 `f`의 winding 분류 결과 (M-B, 분류만 — 방출은 M-C).
#[derive(Debug)]
enum FaceClass {
    /// 통째 keep(`true`)/drop(`false`): 미절단 횡단(균일 정점 클래스) 또는 disjoint-coplanar.
    Whole(bool),
    /// 횡단 절단: M-A의 `(crossings, runs, kept)` — run별 kept.
    Transversal(FaceRunClasses),
    /// 공면 on±: 겹치는 상대 면 `cf`, 법선 동측 여부, Requicha 생존 선택자.
    Coplanar {
        cf: Handle<Face>,
        same_normal: bool,
        survive: PSurvive,
    },
}

/// **unified per-face 분류(M-B)** — 면 `f`(A의)를 상대 solid `other`에 대해 로컬로 분류한다.
/// 글로벌 detector 없이 **`other_classes.contains(&pi)`** 로만 공면-vs-횡단을 가른다(pi=F의 평면
/// 클래스, other_classes=상대 solid의 평면 클래스 집합):
/// - **공면**(pi가 상대 클래스): `coincident_overlap_faces`로 겹치는 상대 면을 찾아 on±
///   (`same_normal`+`coplanar_survival`)로 분류. 겹침 없으면 disjoint-coplanar → 정점 하나로 통째.
///   상대 면이 여럿(다중 cf)이면 정직 거절(후속 마일스톤).
/// - **횡단**(그 외): M-A `classify_face_runs_winding`. 미절단(crossings 없음)이면 `Whole`.
///
/// `keep`=이 op·side가 남기는 쪽(Fuse a-side=Outside), `kind`=op(공면 생존표용). 방출·조립은 M-C+.
#[allow(clippy::too_many_arguments)]
fn classify_face_winding(
    model: &Model,
    f: Handle<Face>,
    other: Handle<Solid>,
    other_classes: &HashSet<usize>,
    kind: BoolKind,
    keep: Side,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    canon: &[usize],
    inc_f: &arrange::EdgePlanes,
    inc_o: &arrange::EdgePlanes,
) -> Result<FaceClass, BoolError> {
    let f_idx = surf_ix[&f];
    let pi = canon[f_idx];
    if other_classes.contains(&pi) {
        // F's plane coincides with a plane class of `other` → coplanar (on±) or disjoint-coplanar.
        let cfs =
            coincident_overlap_faces(model, f, other, pi, planes, surf_ix, canon, inc_f, inc_o)?;
        match cfs.as_slice() {
            [] => {
                // Disjoint coplanar: the shared plane meets `other` but the footprints miss, so no
                // F-vertex lies on `other`'s boundary — F is wholly in/out by any one of them.
                let v = he_start(model, model.faces.get(f).outer.half_edges[0]);
                let side = arrange::point_in_solid_idx(
                    model, v, inc_f, other, inc_o, planes, surf_ix, canon,
                )?;
                Ok(FaceClass::Whole(side == keep))
            }
            [cf] => {
                let cf = *cf;
                let same_normal = planes[f_idx].n_out.dot(planes[surf_ix[&cf]].n_out) > 0.0;
                let (survive, _flip) = coplanar_survival(kind, same_normal);
                Ok(FaceClass::Coplanar {
                    cf,
                    same_normal,
                    survive,
                })
            }
            // Several `other` faces seat on this plane (bar-top hosts two legs) — later milestone.
            _ => Err(reject(tag::COPLANAR_OVERLAP_MULTI)),
        }
    } else {
        // Transversal: no coplanar contribution on F's plane. M-A per-run classification.
        let (crossings, runs, kept) = classify_face_runs_winding(
            model, f, f_idx, other, keep, planes, surf_ix, canon, inc_f, inc_o,
        )?;
        if crossings.is_empty() {
            Ok(FaceClass::Whole(kept[0]))
        } else {
            Ok(FaceClass::Transversal((crossings, runs, kept)))
        }
    }
}

/// 솔리드 `other`의 모든 면의 평면 클래스 집합(`canon[surf_ix[cf]]`) — 로컬 공면 dispatch용.
fn solid_plane_classes(
    model: &Model,
    other: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
    canon: &[usize],
) -> HashSet<usize> {
    solid_shell_handles(model, other)
        .into_iter()
        .flat_map(|sh| model.shells.get(sh).faces.clone())
        .filter_map(|cf| surf_ix.get(&cf).map(|&ix| canon[ix]))
        .collect()
}

/// 솔리드 `qs`의 모든 정점을 `other`에 대해 `point_in_solid_idx`로 분류한 결과(정점 순서대로).
/// `Err`는 그 정점이 exact 인덱스-substrate로 판정 불가함을 뜻한다 — winding 엔진에서 **B 경계
/// 위(on)** 정점(공유평면 위 코너 등)이 여기 걸린다. propagation의 anchor는 `Ok`인 정점만 쓰고,
/// `Err`(on-boundary)는 on± 분류(M-B)가 담당한다.
fn anchor_vertices(
    model: &Model,
    qs: Handle<Solid>,
    other: Handle<Solid>,
) -> Vec<(Handle<Vertex>, Result<Side, BoolError>)> {
    let (planes, surf_ix, inc_q, inc_o, canon) = match plane_index_setup(model, qs, other) {
        Ok(t) => t,
        Err(e) => {
            return solid_vertex_handles(model, qs)
                .into_iter()
                .map(|v| (v, Err(e)))
                .collect();
        }
    };
    solid_vertex_handles(model, qs)
        .into_iter()
        .map(|vh| {
            let side = arrange::point_in_solid_idx(
                model, vh, &inc_q, other, &inc_o, &planes, &surf_ix, &canon,
            );
            (vh, side)
        })
        .collect()
}

/// 횡단 절단된 **hole-free** 면 `f`(평면 인덱스 `plane_idx`, A의 면)의 각 boundary run이 keep되는지
/// 판정한다 — propagation-winding: 원본 정점 anchor(`point_in_solid_idx`) + seam 교차에서
/// alternation 전파(`run_classes`). 반환 `(crossings, runs, kept)`: run별 kept 여부.
///
/// `keep` = 이 op·side가 남기는 쪽(Fuse a-side=Outside 등). 조립·on±·공면은 이 함수 밖(M-B/M-C).
/// inner ring(hole) 있는 면은 `SEAM_BRANCH`류로 정직 거절(M-A 스코프 밖).
#[allow(clippy::too_many_arguments)]
fn classify_face_runs_winding(
    model: &Model,
    f: Handle<Face>,
    plane_idx: usize,
    other: Handle<Solid>,
    keep: Side,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    canon: &[usize],
    inc_f: &arrange::EdgePlanes,
    inc_o: &arrange::EdgePlanes,
) -> Result<FaceRunClasses, BoolError> {
    let face = model.faces.get(f);
    if !face.inner.is_empty() {
        return Err(reject(tag::SEAM_BRANCH)); // M-A scope: hole-free only
    }
    let hes = &face.outer.half_edges;
    let verts: Vec<Handle<Vertex>> = hes.iter().map(|&he| he_start(model, he)).collect();
    // anchor: each outer vertex classified vs `other`; kept iff its side is the keep side.
    let kept_vert: Vec<bool> = verts
        .iter()
        .map(|&v| {
            arrange::point_in_solid_idx(model, v, inc_f, other, inc_o, planes, surf_ix, canon)
                .map(|s| s == keep)
        })
        .collect::<Result<_, _>>()?;

    let paths = arrange::seam_paths_on(model, f, other, planes, surf_ix, inc_f, inc_o)?;
    // ∂f crossings grouped by the outer edge each rides (on_edge = Some).
    let edge_pos: HashMap<Handle<Edge>, usize> =
        hes.iter().enumerate().map(|(i, he)| (he.edge, i)).collect();
    let mut by_edge: Vec<Vec<[usize; 3]>> = vec![Vec::new(); hes.len()];
    for nd in paths.iter().flat_map(|p| p.nodes()) {
        if let Some(e) = nd.on_edge {
            let ix = *edge_pos
                .get(&e)
                .ok_or_else(|| reject(tag::SEAM_COUNT_MISMATCH))?;
            by_edge[ix].push(nd.triple);
        }
    }
    // The seam misses this face (no ∂f crossing): it is wholly kept or dropped by its uniform vertex
    // class. Every vertex must agree — a disagreement means an uncrossed seam (out of M-A scope).
    if by_edge.iter().all(|v| v.is_empty()) {
        if kept_vert.iter().any(|&k| k != kept_vert[0]) {
            return Err(reject(tag::HOLE_CLASS_SPLIT));
        }
        return Ok((
            Vec::new(),
            vec![(0..verts.len()).collect()],
            vec![kept_vert[0]],
        ));
    }
    // `bnd` (face vertex triples) is read by boundary_runs only to order two crossings on one edge.
    let bnd = if by_edge.iter().any(|v| v.len() > 1) {
        arrange::face_vertex_triples(model, f, plane_idx, inc_f)?
    } else {
        Vec::new()
    };
    let br = arrange::boundary_runs(planes, plane_idx, &bnd, &by_edge)?;
    let kept = arrange::run_classes(&br.runs, &kept_vert, &[br.runs.len()])?;
    Ok((br.crossings, br.runs, kept))
}

#[cfg(test)]
mod tests {
    use super::*;

    // M-A anchor 계약(스파이크 발견 잠금). off-boundary 정점은 anchor(propagation seed),
    // B 경계 위(공유평면) 정점은 실패(=on± 대상, M-B). M-B가 딛는 계약.
    #[test]
    fn anchor_on_coplanar_case_leaves_on_boundary_unanchored() {
        // same_ground: a=[0,1]³, b=[0.5,1.5]²×[0,1] — share z=0/z=1. a's (1,1,0)/(1,1,1) corners
        // lie ON b's shared planes AND inside b's footprint = "on" → must NOT anchor.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.0]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        m.rebuild_adjacency();
        let anchors = anchor_vertices(&m, a, b);
        let ok = anchors.iter().filter(|(_, s)| s.is_ok()).count();
        let err: Vec<[i64; 3]> = anchors
            .iter()
            .filter(|(_, s)| s.is_err())
            .map(|(v, _)| {
                let p = m.vertices.get(*v).point.as_array();
                [p[0] as i64, p[1] as i64, p[2] as i64]
            })
            .collect();
        assert_eq!(ok, 6, "6 off-boundary vertices anchor (propagation seeds)");
        assert_eq!(err.len(), 2, "2 on-boundary corners do not (on± target)");
        assert!(
            err.contains(&[1, 1, 0]) && err.contains(&[1, 1, 1]),
            "{err:?}"
        );
        // Every anchored vertex is Outside b (a's off-boundary corners are all outside b's box).
        assert!(
            anchors
                .iter()
                .filter_map(|(_, s)| s.as_ref().ok())
                .all(|s| *s == Side::Outside)
        );
    }

    // M-A hand-known: propagation-winding classifies a's faces vs b for a pure transversal corner
    // overlap (a=[0,1]³, b=[0.5,1.5]³, no shared plane). Fuse a-side keeps Outside b. The 3 faces
    // facing away from b's corner (z=0, x=0, y=0) are wholly kept; the 3 facing it (z=1, x=1, y=1)
    // are each cut into an outside-b run (kept) and an inside-b run (dropped). Exact substrate only.
    #[test]
    fn classify_transversal_corner_overlap() {
        let run = || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let b = m.add_cuboid(
                Point3::from_array([0.5, 0.5, 0.5]),
                Point3::from_array([1.5, 1.5, 1.5]),
            );
            m.rebuild_adjacency();
            let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
            let shell = m.shells.get(m.solids.get(a).outer);
            let mut whole = 0;
            let mut cut = 0;
            let mut sig: Vec<(usize, Vec<bool>)> = Vec::new();
            for &fh in &shell.faces {
                let pidx = surf_ix[&fh];
                let (cx, _runs, kept) = classify_face_runs_winding(
                    &m,
                    fh,
                    pidx,
                    b,
                    Side::Outside,
                    &planes,
                    &surf_ix,
                    &canon,
                    &inc_a,
                    &inc_b,
                )
                .unwrap();
                if cx.is_empty() {
                    whole += 1;
                    assert_eq!(kept, vec![true], "uncut a-face (outside b) is wholly kept");
                } else {
                    cut += 1;
                    assert_eq!(cx.len(), 2, "each cut face has two crossings");
                    assert_eq!(
                        kept.iter().filter(|&&k| k).count(),
                        1,
                        "cut face: one run kept (outside b), one dropped (inside)"
                    );
                }
                sig.push((cx.len(), kept));
            }
            assert_eq!((whole, cut), (3, 3), "3 uncut + 3 cut faces");
            sig
        };
        assert_eq!(run(), run(), "deterministic");
    }

    // Run every a-face of a boolean through the unified per-face classifier, returning
    // `(n_out, Ok(FaceClass) | Err(reject-tag))` per face — the M-B dispatch, no global detector.
    fn classify_solid_faces(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        kind: BoolKind,
        keep: Side,
    ) -> Vec<([i64; 3], Result<FaceClass, &'static str>)> {
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(m, a, b).unwrap();
        let other_classes = solid_plane_classes(m, b, &surf_ix, &canon);
        let shell = m.shells.get(m.solids.get(a).outer);
        shell
            .faces
            .iter()
            .map(|&fh| {
                let n = planes[surf_ix[&fh]].n_out.as_array();
                let key = [n[0] as i64, n[1] as i64, n[2] as i64];
                LAST_REJECT.with(|c| c.take());
                let fc = classify_face_winding(
                    m,
                    fh,
                    b,
                    &other_classes,
                    kind,
                    keep,
                    &planes,
                    &surf_ix,
                    &canon,
                    &inc_a,
                    &inc_b,
                )
                .map_err(|_| LAST_REJECT.with(|c| c.get()).unwrap());
                (key, fc)
            })
            .collect()
    }

    // M-B core (measured 2026-07-20): the unified per-face dispatch routes each a-face by the LOCAL
    // check `pi ∈ other_classes` — NO global detector. same_ground Fuse (a=[0,1]³, b=[0.5,1.5]²×[0,1],
    // shared z=0/z=1 caps): the two caps land on the coplanar on± path (same_normal → Fuse=Whole);
    // the four walls land on the transversal path (dispatch correct) but hit M-A's coplanar-adjacent
    // FRONTIER — the +axis walls' (1,1,·) corners are on b's boundary (anchor `no_clear_ray`), the
    // −axis walls' z-edges lie on the shared caps (seam `vertex_on_face_plane`). Closing that
    // frontier (partial-anchor propagation + coplanar-aware seam) is a later milestone; here we lock
    // that the dispatch + on± are correct and the walls reject through the transversal path.
    #[test]
    fn same_ground_dispatch_caps_coplanar_walls_frontier() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.0]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        m.rebuild_adjacency();
        let faces = classify_solid_faces(&m, a, b, BoolKind::Fuse, Side::Outside);
        let mut cap_cfs: Vec<Handle<Face>> = Vec::new();
        let mut wall_tags: Vec<&str> = Vec::new();
        for (n, fc) in &faces {
            if n[2] != 0 {
                // cap: coplanar on± — shared caps point the same way, Fuse/same = Whole. Each maps to
                // its own overlapping b-cap.
                match fc {
                    Ok(FaceClass::Coplanar {
                        cf,
                        same_normal: true,
                        survive: PSurvive::Whole,
                    }) => cap_cfs.push(*cf),
                    _ => panic!("cap {n:?} -> {fc:?}"),
                }
            } else {
                // wall: dispatched to the transversal path, rejects at the M-A frontier.
                wall_tags.push(*fc.as_ref().unwrap_err());
            }
        }
        assert_eq!(cap_cfs.len(), 2, "two coplanar caps");
        assert_ne!(cap_cfs[0], cap_cfs[1], "each cap overlaps a distinct b-cap");
        wall_tags.sort_unstable();
        assert_eq!(
            wall_tags,
            [
                "no_clear_ray",
                "no_clear_ray",
                "vertex_on_face_plane",
                "vertex_on_face_plane"
            ],
            "four walls reject through the transversal path (coplanar-adjacent frontier)"
        );
        // determinism
        let again = classify_solid_faces(&m, a, b, BoolKind::Fuse, Side::Outside);
        let tag_of = |v: &[([i64; 3], Result<FaceClass, &'static str>)]| {
            v.iter()
                .map(|(n, fc)| (*n, fc.as_ref().err().copied()))
                .collect::<Vec<_>>()
        };
        assert_eq!(tag_of(&faces), tag_of(&again), "deterministic");
    }

    // Unified entry on a PURE transversal (no shared plane): every a-face routes to the transversal
    // path (pi ∉ other_classes) and classifies — 3 uncut faces `Whole(true)`, 3 cut `Transversal`.
    // Proves the dispatch feeds M-A correctly when there is no coplanar adjacency (the classifiable
    // regime), complementing same_ground where the walls hit the frontier.
    #[test]
    fn unified_dispatch_pure_transversal_classifies() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        m.rebuild_adjacency();
        let faces = classify_solid_faces(&m, a, b, BoolKind::Fuse, Side::Outside);
        let mut whole = 0;
        let mut cut = 0;
        for (n, fc) in &faces {
            match fc {
                Ok(FaceClass::Whole(true)) => whole += 1,
                Ok(FaceClass::Transversal((cx, _, kept))) => {
                    assert_eq!(cx.len(), 2, "{n:?}: cut face has two crossings");
                    assert_eq!(
                        kept.iter().filter(|&&k| k).count(),
                        1,
                        "{n:?}: one run kept"
                    );
                    cut += 1;
                }
                other => panic!("{n:?}: unexpected {other:?}"),
            }
        }
        assert_eq!((whole, cut), (3, 3), "3 uncut Whole + 3 cut Transversal");
    }

    // Disjoint coplanar: a=[0,1]³, b=[2,3]²×[0,1] share the z=0/z=1 planes but their footprints
    // miss. a's caps dispatch to the coplanar path, find no overlapping b-face (empty cfs), and are
    // classified whole by a single vertex (all outside b) → `Whole(true)` for Fuse. Locks the
    // disjoint-coplanar branch (`cfs == []`).
    #[test]
    fn disjoint_coplanar_caps_whole() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([2.0, 2.0, 0.0]),
            Point3::from_array([3.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let faces = classify_solid_faces(&m, a, b, BoolKind::Fuse, Side::Outside);
        let caps: Vec<_> = faces.iter().filter(|(n, _)| n[2] != 0).collect();
        assert_eq!(caps.len(), 2);
        for (n, fc) in caps {
            assert!(
                matches!(fc, Ok(FaceClass::Whole(true))),
                "disjoint-coplanar cap {n:?} -> {fc:?}"
            );
        }
    }

    // Pure transversal (no shared plane): every vertex anchors — propagation has full seeds.
    #[test]
    fn anchor_on_transversal_case_classifies_all() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        m.rebuild_adjacency();
        let anchors = anchor_vertices(&m, a, b);
        assert!(
            anchors.iter().all(|(_, s)| s.is_ok()),
            "all anchor (pure transversal)"
        );
        // Only a's (1,1,1) corner is inside b=[0.5,1.5]³; the rest Outside.
        for (vh, s) in &anchors {
            let p = m.vertices.get(*vh).point.as_array();
            let want = if p.iter().all(|&c| c > 0.5 && c < 1.5) {
                Side::Inside
            } else {
                Side::Outside
            };
            assert_eq!(*s.as_ref().unwrap(), want, "{p:?}");
        }
    }
}
