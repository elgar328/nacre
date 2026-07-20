#![cfg_attr(not(test), allow(dead_code))]
//! Segment trace producer for the plane-class arrangement (the B engine). Isolated and unwired:
//! nothing outside this module and its tests calls it, following the `winding.rs` precedent.
//!
//! Unlike [`crate::section_of_solid`], which rejects the whole solid the moment any vertex lies
//! on the cut plane and then assembles **closed loops** (needing per-face even parity, degree-2,
//! and an alternating walk), this produces **line segments** — the pieces where each face of the
//! solid meets a plane class — and leaves assembly to the arrangement. Geometry lying in the
//! plane is *input*, not a degeneracy.
//!
//! Two face cases are handled. **Seated** — a face whose own class is the cut class — lies wholly
//! in the plane, so its whole boundary is trace. **Transversal / mixed** — a face crossing the
//! plane, possibly with an edge lying in it — is 2D polygon-vs-line clipping on `L = W ∩ fp`, with
//! a three-valued (−/0/+) scan that collapses on-line edges instead of rejecting them. Faces with
//! **holes** are declined this brick (a forced-covered run inside a hole would claim its void as
//! material, and the corpus has no case to verify it). Declining is a per-face record, never a
//! whole-solid abort — a non-empty `declined` means "this trace is incomplete, do not conclude
//! from it", which is the defect `section_of_solid` had.

use super::*;

/// Which operand a segment came from — the boolean's per-cell label needs both solids' material
/// above and below, so provenance cannot be merged away.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SolidSide {
    A,
    B,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SegKind {
    /// A face crossing the plane transversally. Crossing this segment in the arrangement flips
    /// **both** the above-label and the below-label together. `mat` is which side of the directed
    /// line `W ∩ fp` the solid's material lies on (`+1`/`-1`, relative to `d = n_W × n_fp`).
    Transversal { mat: i8 },
    /// A boundary edge of a face lying in the plane. Crossing it flips **one** label — the side
    /// the seated face's body occupies (`body_above`).
    Seated { body_above: bool },
}

/// One segment of a solid's trace on a plane class, named entirely in plane-class triples.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Seg {
    /// Canon plane class the segment rides besides the cut plane (its line is `W ∩ wall`).
    pub wall: usize,
    /// Endpoint triples, canonized.
    pub end: [[usize; 3]; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
}

/// A solid's trace on one plane class. `declined` non-empty ⇒ incomplete: a consumer must not
/// read "no segments" as "the plane misses the solid".
#[derive(Default, Debug)]
pub(crate) struct Trace {
    pub segs: Vec<Seg>,
    /// Single-point tangential contacts — real arrangement vertices, but not segments (a
    /// zero-length chord would abort at `Line::through_points`).
    pub touches: Vec<[usize; 3]>,
    /// `(face plane index, reason)` for every face this brick could not trace.
    pub declined: Vec<(usize, &'static str)>,
}

/// Canonize a raw triple and sort it.
fn canon3(t: [usize; 3], canon: &[usize]) -> [usize; 3] {
    let mut c = [canon[t[0]], canon[t[1]], canon[t[2]]];
    c.sort_unstable();
    c
}

/// One feature node on the line `L = W ∩ fp`: a point named by its third plane `r` (raw index).
struct Node {
    /// Third plane naming this point on `L` (the point is `{W, fp, r}`).
    r: usize,
    /// Whether crossing this node toggles inside/outside. Strict crossings and single-vertex runs
    /// set it in Phase A; a two-vertex run's *hi* node has it set after the sort.
    flip: bool,
    /// A two-vertex on-line edge's run id; the two nodes share it, and the gap between them is
    /// forced-covered (it is boundary of `F`, hence of `F ∩ W`).
    run: Option<usize>,
    /// For run nodes: whether the flanks (the off-line neighbours) sit on opposite sides.
    flanks_differ: bool,
    /// A single on-line vertex with equal flanks — a tangential touch, not a crossing.
    single_touch: bool,
}

/// Trace one **transversal or mixed** face `f` (plane `fp ≠ W`) onto plane class `wc`, as the
/// segments where `f` meets `W`. This is 2D polygon-vs-line clipping on `L = W ∩ fp`, handling
/// on-line edges (a run of `side == 0` vertices) instead of rejecting them the way
/// `section_of_solid` does. Holes are declined this brick (a forced-covered run inside a hole
/// would wrongly claim the hole's void as material, and the corpus has no case to verify it).
///
/// Every emitted segment rides `wall = fp` (its line is `W ∩ fp`); the side walls appear only as
/// the third plane naming each endpoint. `mat = orient_sign(fp)` is a per-face constant.
#[allow(clippy::too_many_arguments)]
fn trace_transversal_face(
    model: &Model,
    fh: Handle<Face>,
    fp: usize,
    which: SolidSide,
    wc: usize,
    planes: &[PlaneInfo],
    inc: &arrange::EdgePlanes,
    canon: &[usize],
    out: &mut Trace,
) {
    if arrange::hole_rings(model, fh, fp, inc).is_ok_and(|h| !h.is_empty()) {
        out.declined.push((fp, "has-holes"));
        return;
    }
    let ring = match arrange::face_vertex_triples(model, fh, fp, inc) {
        Ok(r) => r,
        Err(_) => {
            out.declined.push((fp, "outer-ring"));
            return;
        }
    };
    let n = ring.len();
    let side: Vec<i8> = (0..n)
        .map(|i| arrange::side_of(planes, ring[i], wc))
        .collect();
    let Some(start) = side.iter().position(|&s| s != 0) else {
        // A face off `W` cannot have every vertex on `W`; if it does, do not guess.
        out.declined.push((fp, "all-on-plane"));
        return;
    };

    // `L`'s third plane naming a point of the on-line edge: the vertex triple `{fp, W-class, r}`.
    let third_on_l = |t: [usize; 3]| -> Option<usize> {
        let (mut r, mut has_fp, mut has_w) = (None, false, false);
        for &x in &t {
            if x == fp {
                has_fp = true;
            } else if canon[x] == wc {
                has_w = true;
            } else if r.replace(x).is_some() {
                return None; // two off-planes: not a clean point on L
            }
        }
        (has_fp && has_w).then_some(r).flatten()
    };

    // Phase A — one pass around the ring collecting feature nodes.
    let mut nodes: Vec<Node> = Vec::new();
    let mut run_counter = 0usize;
    let mut declined: Option<&'static str> = None;
    let mut j = 0;
    while j < n {
        let i = (start + j) % n;
        if side[i] != 0 {
            let ni = (i + 1) % n;
            if side[ni] != 0 && side[ni] != side[i] {
                // Strict crossing on edge i; its wall is the plane the edge rides besides fp.
                match arrange::ring_edge(fp, &ring, i) {
                    Ok((wall, _, _)) => nodes.push(Node {
                        r: wall,
                        flip: true,
                        run: None,
                        flanks_differ: false,
                        single_touch: false,
                    }),
                    Err(_) => declined = Some("crossing-name"),
                }
            }
            j += 1;
        } else {
            // A maximal run of side==0 vertices (at most two, since three would be a straight
            // angle already rejected by `loop_triples`).
            let run_start = i;
            let mut m = 0;
            while j < n && side[(start + j) % n] == 0 {
                m += 1;
                j += 1;
            }
            let before = side[(run_start + n - 1) % n];
            let after = side[(start + j) % n];
            let flanks_differ = before != after;
            let name = |k: usize| third_on_l(ring[(run_start + k) % n]);
            if m == 1 {
                match name(0) {
                    Some(r) => nodes.push(Node {
                        r,
                        flip: flanks_differ,
                        run: None,
                        flanks_differ,
                        single_touch: !flanks_differ,
                    }),
                    None => declined = Some("run-name"),
                }
            } else if m == 2 {
                let id = run_counter;
                run_counter += 1;
                match (name(0), name(1)) {
                    (Some(ra), Some(rb)) => {
                        for r in [ra, rb] {
                            nodes.push(Node {
                                r,
                                flip: false,
                                run: Some(id),
                                flanks_differ,
                                single_touch: false,
                            });
                        }
                    }
                    _ => declined = Some("run-name"),
                }
            } else {
                declined = Some("long-run");
            }
        }
        if declined.is_some() {
            break;
        }
    }
    if let Some(reason) = declined {
        out.declined.push((fp, reason));
        return;
    }

    // Phase B — order the nodes along L and fix run structure.
    nodes.sort_by(
        |a, b| match arrange::order_along(planes, wc, fp, a.r, b.r) {
            -1 => std::cmp::Ordering::Less,
            1 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        },
    );
    for w in nodes.windows(2) {
        if arrange::order_along(planes, wc, fp, w[0].r, w[1].r) == 0 {
            out.declined.push((fp, "coincident-features"));
            return;
        }
    }
    // Each run's two nodes must be adjacent after the sort; the hi node carries the flip.
    for k in 0..nodes.len() {
        if let Some(id) = nodes[k].run {
            let prev = k > 0 && nodes[k - 1].run == Some(id);
            let next = k + 1 < nodes.len() && nodes[k + 1].run == Some(id);
            if !prev && !next {
                out.declined.push((fp, "run-split"));
                return;
            }
            if prev {
                nodes[k].flip = nodes[k].flanks_differ;
            }
        }
    }

    // Phase C — sweep left to right; covered gaps merge into segments.
    let mut parity = 0i8;
    let mut seg_start: Option<usize> = None;
    let emit = |a: usize, b: usize, out: &mut Trace| {
        out.segs.push(Seg {
            wall: canon[fp],
            end: [canon3([wc, fp, a], canon), canon3([wc, fp, b], canon)],
            solid: which,
            kind: SegKind::Transversal {
                mat: arrange::orient_sign(planes, fp),
            },
        });
    };
    for k in 0..nodes.len() {
        if nodes[k].single_touch && parity == 0 {
            out.touches.push(canon3([wc, fp, nodes[k].r], canon));
        }
        if nodes[k].flip {
            parity ^= 1;
        }
        let forced =
            k + 1 < nodes.len() && nodes[k].run.is_some() && nodes[k].run == nodes[k + 1].run;
        let covered = parity == 1 || forced;
        match (seg_start, covered) {
            (None, true) => seg_start = Some(nodes[k].r),
            (Some(a), false) => {
                emit(a, nodes[k].r, out);
                seg_start = None;
            }
            _ => {}
        }
    }
    if seg_start.is_some() || parity != 0 {
        // A line enters and leaves a bounded region equally; an unbalanced sweep is degenerate.
        out.declined.push((fp, "odd-parity"));
    }
}

/// Trace one solid on plane class `wc` (a canon root, i.e. an index into `planes`). This brick:
/// seated faces → their boundary as segments; every other face → `declined`.
#[allow(clippy::too_many_arguments)]
fn trace_one(
    model: &Model,
    solid: Handle<Solid>,
    which: SolidSide,
    wc: usize,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &arrange::EdgePlanes,
    canon: &[usize],
    out: &mut Trace,
) {
    let w_normal = planes[wc].plane.normal();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let fp = surf_ix[&fh];
            if canon[fp] != wc {
                trace_transversal_face(model, fh, fp, which, wc, planes, inc, canon, out);
                continue;
            }
            // Seated: the face lies in W, so its whole boundary is trace. The body lies on one
            // side of W — `n_out` points away from the body, so the body is above W exactly when
            // `n_out · n_W < 0`. (f64 sign is exact for axis-aligned normals; rotated inputs are
            // outside this corpus and this brick does not claim them.)
            let body_above = planes[fp].n_out.dot(w_normal) < 0.0;
            let kind = SegKind::Seated { body_above };
            let mut emit_ring = |tris: &[[usize; 3]]| {
                let n = tris.len();
                for i in 0..n {
                    // Edge i runs vertex i → vertex i+1; the wall it rides is the plane the two
                    // endpoint triples share besides `fp`.
                    let (t0, t1) = (tris[i], tris[(i + 1) % n]);
                    let shared: Vec<usize> = t0
                        .iter()
                        .copied()
                        .filter(|&x| x != fp && t1.contains(&x))
                        .collect();
                    let [wall] = shared[..] else {
                        out.declined.push((fp, "seated-edge-naming"));
                        continue;
                    };
                    out.segs.push(Seg {
                        wall: canon[wall],
                        end: [canon3(t0, canon), canon3(t1, canon)],
                        solid: which,
                        kind,
                    });
                }
            };
            match arrange::face_vertex_triples(model, fh, fp, inc) {
                Ok(ts) => emit_ring(&ts),
                Err(_) => {
                    out.declined.push((fp, "outer-ring"));
                    continue;
                }
            }
            if let Ok(rings) = arrange::hole_rings(model, fh, fp, inc) {
                for r in rings {
                    emit_ring(&r);
                }
            }
        }
    }
}

/// Both operands' traces on plane class `wc`, merged into one `Trace` (segments keep their
/// `solid` tag).
#[allow(clippy::too_many_arguments)]
fn trace_on_class(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    wc: usize,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &arrange::EdgePlanes,
    inc_b: &arrange::EdgePlanes,
    canon: &[usize],
) -> Trace {
    let mut out = Trace::default();
    trace_one(
        model,
        a,
        SolidSide::A,
        wc,
        planes,
        surf_ix,
        inc_a,
        canon,
        &mut out,
    );
    trace_one(
        model,
        b,
        SolidSide::B,
        wc,
        planes,
        surf_ix,
        inc_b,
        canon,
        &mut out,
    );
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Point of a canon triple, for asserting geometry by hand.
    fn pt(t: [usize; 3], planes: &[PlaneInfo]) -> [f64; 3] {
        three_planes(
            &planes[t[0]].plane,
            &planes[t[1]].plane,
            &planes[t[2]].plane,
        )
        .unwrap()
        .as_array()
    }

    /// A seated face's whole boundary is emitted, with `body_above` matching the geometry.
    ///
    /// Two stacked cubes sharing the plane z=1: `a=[0,1]³`, `b=[0,1]²×[1,2]`. On the shared
    /// class, `a`'s top cap is seated with its body **below** (body_above=false) and `b`'s bottom
    /// cap is seated with its body **above** (body_above=true). Each emits its 4 boundary edges.
    #[test]
    fn a_seated_cap_emits_its_boundary_with_the_right_body_side() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();

        // The shared class z=1: the class both solids seat a cap on.
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                let seats =
                    |s: Handle<Solid>| {
                        solid_shell_handles(&m, s).into_iter().any(|sh| {
                            m.shells.get(sh).faces.iter().any(|fh| {
                                canon[surf_ix[fh]] == c && face_on_z1(*fh, &surf_ix, &planes)
                            })
                        })
                    };
                seats(a) && seats(b)
            })
            .expect("a shared cap class");

        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        // Two seated caps × 4 edges each. (The side walls now also trace as transversal chords;
        // the two far caps are parallel to z=1 and contribute nothing — misses, not declines.)
        let seated: Vec<&Seg> = tr
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Seated { .. }))
            .collect();
        assert_eq!(seated.len(), 8, "two seated caps, four edges each: {tr:?}");
        // Four side walls per cube each cross z=1 in one chord.
        let trans = tr
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .count();
        assert_eq!(trans, 8, "eight side-wall chords: {tr:?}");
        assert!(
            tr.declined.is_empty(),
            "axis-aligned cubes decline nothing: {tr:?}"
        );
        assert!(tr.touches.is_empty());
        for s in &seated {
            // On z=1, riding a side wall (not W itself), body on the geometrically correct side.
            for e in &s.end {
                assert!((pt(*e, &planes)[2] - 1.0).abs() < 1e-12, "endpoint on z=1");
            }
            assert_ne!(s.wall, wc, "a seated edge rides a side wall, not W");
            match (s.solid, s.kind) {
                (SolidSide::A, SegKind::Seated { body_above }) => assert!(!body_above, "a below"),
                (SolidSide::B, SegKind::Seated { body_above }) => assert!(body_above, "b above"),
                _ => unreachable!(),
            }
        }
    }

    /// Round a point for set membership by hand.
    fn rp(p: [f64; 3]) -> [i64; 3] {
        [
            (p[0] * 1e6).round() as i64,
            (p[1] * 1e6).round() as i64,
            (p[2] * 1e6).round() as i64,
        ]
    }

    /// The transversal chord across a non-convex reflex plane — the cap chord — is produced whole,
    /// collapsing the on-line edge, which the seated brick provably cannot do (it declines this
    /// face). l_prism (profile (0,0),(2,0),(2,1),(1,1),(1,2),(0,2), extruded z∈[0,1]) cut at y=1:
    /// the z=0 cap's trace is the single segment x∈[0,2], **not** [0,1]+[1,2].
    #[test]
    fn the_cap_chord_collapses_the_on_line_edge() {
        let profile = Profile2d {
            points: vec![
                Point2::from_array([0.0, 0.0]),
                Point2::from_array([2.0, 0.0]),
                Point2::from_array([2.0, 1.0]),
                Point2::from_array([1.0, 1.0]),
                Point2::from_array([1.0, 2.0]),
                Point2::from_array([0.0, 2.0]),
            ],
        };
        let mut m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        // A far cube so plane_index_setup has two solids; it never touches the y=1 class.
        m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([11.0, 11.0, 11.0]),
        );
        m.rebuild_adjacency();
        let a2 = m.live_solids[0];
        let b2 = m.live_solids[1];
        let (planes, surf_ix, inc_a, _inc_b, canon) = plane_index_setup(&m, a2, b2).unwrap();

        // The y=1 class: a plane through all-y=1 points that a's reflex face sits on.
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                planes.iter().enumerate().any(|(i, p)| {
                    canon[i] == c && p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
                })
            })
            .expect("a y=1 class");

        let mut out = Trace::default();
        trace_one(
            &m,
            a2,
            SolidSide::A,
            wc,
            &planes,
            &surf_ix,
            &inc_a,
            &canon,
            &mut out,
        );

        // Collect transversal segments as rounded endpoint-pairs.
        let chords: Vec<([i64; 3], [i64; 3])> = out
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .map(|s| {
                let (mut u, mut v) = (rp(pt(s.end[0], &planes)), rp(pt(s.end[1], &planes)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u, v)
            })
            .collect();

        // The cap chords: the z=0 (and z=1) cap ∩ (y=1) is the single segment x∈[0,2]. A cap is a
        // horizontal segment (both endpoints share z), which distinguishes it from the reflex wall
        // x=1's own vertical chord (1,1,0)-(1,1,1), a legitimate different face's trace.
        let horiz_z0: Vec<_> = chords
            .iter()
            .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
            .collect();
        // ★ Falsifiable core: a failed collapse yields TWO horizontal chords ([0,1]+[1,2], an x=1
        // endpoint appears); a correct collapse yields exactly ONE spanning [0,2].
        assert_eq!(
            horiz_z0.len(),
            1,
            "the on-line edge collapsed to one chord: {chords:?}"
        );
        assert_eq!(
            *horiz_z0[0],
            (rp([0.0, 1.0, 0.0]), rp([2.0, 1.0, 0.0])),
            "cap chord spans x∈[0,2]: {chords:?}"
        );
        // The seated brick provably cannot produce this: it declines the z=0 cap as `not-seated`
        // (its class is z=0, not y=1). So a transversal chord riding fp=z=0 is new capability.
        // mat is a per-face constant, never 0.
        for s in &out.segs {
            if let SegKind::Transversal { mat } = s.kind {
                assert!(mat == 1 || mat == -1, "mat is a clean side: {mat}");
            }
        }
    }

    /// The flanks-equal (tangential on-line edge) branch: u_prism's notch bottom on y=1. The run
    /// x∈[1,2] is flanked by material on **both** sides (the two prongs), so it does not toggle
    /// parity yet must still appear inside the trace. The cap ∩ y=1 is the single span x∈[0,3].
    #[test]
    fn a_tangential_on_line_edge_still_spans() {
        let profile = Profile2d {
            points: vec![
                Point2::from_array([0.0, 0.0]),
                Point2::from_array([3.0, 0.0]),
                Point2::from_array([3.0, 2.3]),
                Point2::from_array([2.0, 2.3]),
                Point2::from_array([2.0, 1.0]),
                Point2::from_array([1.0, 1.0]),
                Point2::from_array([1.0, 2.0]),
                Point2::from_array([0.0, 2.0]),
            ],
        };
        let mut m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([11.0, 11.0, 11.0]),
        );
        m.rebuild_adjacency();
        let a2 = m.live_solids[0];
        let b2 = m.live_solids[1];
        let (planes, surf_ix, inc_a, _inc_b, canon) = plane_index_setup(&m, a2, b2).unwrap();
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                planes.iter().enumerate().any(|(i, p)| {
                    canon[i] == c && p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
                })
            })
            .expect("a y=1 class");
        let mut out = Trace::default();
        trace_one(
            &m,
            a2,
            SolidSide::A,
            wc,
            &planes,
            &surf_ix,
            &inc_a,
            &canon,
            &mut out,
        );

        let horiz_z0: Vec<_> = out
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .map(|s| {
                let (mut u, mut v) = (rp(pt(s.end[0], &planes)), rp(pt(s.end[1], &planes)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u, v)
            })
            .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
            .collect();
        assert_eq!(
            horiz_z0.len(),
            1,
            "tangential run bridged into one span: {horiz_z0:?}"
        );
        assert_eq!(
            horiz_z0[0],
            (rp([0.0, 1.0, 0.0]), rp([3.0, 1.0, 0.0])),
            "spans x∈[0,3]"
        );
    }

    fn face_on_z1(
        fh: Handle<Face>,
        surf_ix: &HashMap<Handle<Face>, usize>,
        planes: &[PlaneInfo],
    ) -> bool {
        let p = &planes[surf_ix[&fh]];
        // A cap in the z=1 plane: all three defining points at z=1.
        p.tri.iter().all(|q| (q.as_array()[2] - 1.0).abs() < 1e-12)
    }

    /// Axis-aligned cubes never decline: every face is seated, a clean transversal chord, or a
    /// parallel miss. This is falsifiable — a `declined` entry would mean the producer hit a
    /// degeneracy it should not on this input.
    #[test]
    fn axis_aligned_cubes_decline_nothing() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        for wc in {
            let mut c: Vec<usize> = canon.clone();
            c.sort_unstable();
            c.dedup();
            c
        } {
            let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
            assert!(tr.declined.is_empty(), "class {wc} declined: {tr:?}");
        }
    }
}
