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
//! **This first brick handles `seated` faces exactly** — a face whose own plane class is the cut
//! class, so the face lies wholly in the plane and its whole boundary is trace. Transversal and
//! mixed faces (a face that crosses the plane, possibly with an edge lying in it) are **declined
//! with a reason**; the three-valued run collapse they need is the next brick. Declining is a
//! per-face record, never a whole-solid abort — a non-empty `declined` means "this trace is
//! incomplete, do not conclude from it", which is the defect `section_of_solid` had.

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
    /// **both** the above-label and the below-label together. (Emitted by the next brick.)
    #[cfg_attr(test, allow(dead_code))]
    Transversal,
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
                // Not seated on this class. Transversal / mixed tracing is the next brick.
                out.declined.push((fp, "not-seated"));
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
        // Two seated caps × 4 edges each; the other 10 faces (8 side walls + 2 far caps) decline.
        assert_eq!(tr.segs.len(), 8, "two caps, four edges each: {tr:?}");
        assert_eq!(
            tr.declined.len(),
            10,
            "the non-seated faces are recorded: {tr:?}"
        );
        assert!(tr.touches.is_empty());
        for s in &tr.segs {
            // On z=1, riding a side wall (not W itself), body on the geometrically correct side.
            for e in &s.end {
                assert!((pt(*e, &planes)[2] - 1.0).abs() < 1e-12, "endpoint on z=1");
            }
            assert_ne!(s.wall, wc, "a seated edge rides a side wall, not W");
            match (s.solid, s.kind) {
                (SolidSide::A, SegKind::Seated { body_above }) => assert!(!body_above, "a below"),
                (SolidSide::B, SegKind::Seated { body_above }) => assert!(body_above, "b above"),
                _ => panic!("seated segments only"),
            }
        }
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

    /// No face is silently dropped: on every class, (seated faces) + (declined faces) equals the
    /// 12 faces of the two cubes. This is falsifiable — a face that fell through both arms would
    /// make the sum 11.
    #[test]
    fn every_face_is_traced_or_declined_never_dropped() {
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
            // Axis-aligned cubes never hit `seated-edge-naming`, so seated faces = segs/4.
            assert_eq!(tr.segs.len() % 4, 0, "class {wc}: whole caps only");
            let seated = tr.segs.len() / 4;
            assert_eq!(
                seated + tr.declined.len(),
                12,
                "class {wc}: no face dropped: {tr:?}"
            );
        }
    }
}
