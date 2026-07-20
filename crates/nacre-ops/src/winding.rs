//! winding 기반 단일 arrangement 엔진 (M5) — 설계 `docs/winding-engine.md`.
//!
//! 현 case-routed 엔진(detector + 생존표)을 대체할 문헌식 winding 엔진을 **격리·미배선**으로
//! 증분 구축한다(cold-start 금지). 이 모듈의 코드는 커토버(M-E)까지 프로덕션에서 안 쓰이며,
//! 매 마일스톤 격리 테스트로 검증한다. 모든 판정은 exact plane-triple substrate(좌표 미읽음).
//!
//! **M-A(현재):** per-face propagation-winding **분류** 프리미티브 — 원본 정점 anchor
//! (`point_in_solid_idx`) + seam 교차에서 alternation 전파(`run_classes`)로, 횡단 절단된 면의
//! 각 boundary run이 keep되는지 판정. 순수 횡단·hole-free 면만(공면 on±·조립은 M-B/M-C).
#![cfg_attr(not(test), allow(dead_code))]

use super::*;

/// per-run 분류 결과: `(seam 교차 triple들, ∂f run별 원본정점, run별 kept 여부)`.
type FaceRunClasses = (Vec<[usize; 3]>, Vec<Vec<usize>>, Vec<bool>);

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
