//! Toleranced boolean predicates (overhaul stage 3): the geom sign predicates routed
//! through TIP so they stay exact under rotation.
//!
//! A rotated face's plane coefficients and `tri` coordinates are rounded irrationals, so
//! the axis-aligned geom predicates (`nacre_geom::intersect`) are exact only w.r.t. the
//! *rounded* geometry. When the boolean's operands are rotated, these wrappers rebuild
//! each plane from the three exact `Pt3` its face carries ([`PlaneInfo::tri_pt3`]) and
//! decide the sign with the `nacre-scalar::frame3` indirect judges instead. The routing
//! is a per-boolean `rotated` flag ([`solid_is_rotated`](crate::solid_is_rotated)) — the
//! axis-aligned hot path is unchanged and never builds a `Pt3`.
//!
//! Stage 3a-i wires only [`t_orient3d`] (the `order_along` comparator); the live
//! arrangement still calls the geom predicates until stage 3b.

use crate::PlaneInfo;
use nacre_geom::intersect::three_plane_orient3d;
use nacre_scalar::Orient;
use nacre_scalar::frame3::indirect_orient3d_judge;

/// The sign of `orient3d(V, tri_j)` where `V = ∩(planes p, q, r)` is an implicit point —
/// the toleranced twin of [`three_plane_orient3d`], matching `order_along`'s shape so it
/// is a drop-in (`+1`/`-1`/`0`).
///
/// - `!rotated`: the exact axis-aligned path — `three_plane_orient3d` on the stored plane
///   coefficients and `tri` coordinates (unchanged, fast).
/// - `rotated`: each of `p, q, r` and the explicit triangle `j` is taken as the three
///   exact `Pt3` in [`PlaneInfo::tri_pt3`], and [`indirect_orient3d_judge`] decides the
///   sign from the definitions — never materializing `V` or reading the rounded `tri`.
///
/// Winding-invariant, so the three points defining each plane may be in any order.
// `#[allow(dead_code)]`: the wrapper lands in 3a-i (unit-tested); `order_along` calls it
// (retiring this allow) in stage 3b.
#[allow(dead_code)]
pub(crate) fn t_orient3d(
    planes: &[PlaneInfo],
    p: usize,
    q: usize,
    r: usize,
    j: usize,
    rotated: bool,
) -> i8 {
    if !rotated {
        return three_plane_orient3d(
            &planes[p].plane,
            &planes[q].plane,
            &planes[r].plane,
            planes[j].tri[0],
            planes[j].tri[1],
            planes[j].tri[2],
        );
    }
    let def = |k: usize| {
        let t = planes[k]
            .tri_pt3
            .as_ref()
            .expect("a rotated plane must carry tri_pt3 (collect_planes builds it)");
        (&t[0], &t[1], &t[2])
    };
    let jd = planes[j]
        .tri_pt3
        .as_ref()
        .expect("a rotated plane must carry tri_pt3 (collect_planes builds it)");
    match indirect_orient3d_judge(def(p), def(q), def(r), &jd[0], &jd[1], &jd[2]) {
        Orient::Positive => 1,
        Orient::Negative => -1,
        Orient::Zero => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Operation, apply, collect_planes};
    use nacre_math::Point3;
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation as SRot};
    use nacre_store::Handle;
    use nacre_topo::{Model, Solid};

    fn cuboid() -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        (m, s)
    }

    /// Rotate solid `s` by 30° about Z through a rational pivot; return the new solid.
    fn rotated(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
        let iso = Isometry::rotation(SRot {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        let out = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: iso,
            },
        )
        .unwrap();
        m.rebuild_adjacency();
        match out {
            crate::OpOutput::Transform { solid } => solid,
            _ => panic!("expected Transform"),
        }
    }

    /// The signed volume of the normal triple `(n_p, n_q, n_r)` — nonzero iff the three
    /// planes meet in a single point (so `three_plane_orient3d` is well-defined).
    fn normals_independent(planes: &[PlaneInfo], p: usize, q: usize, r: usize) -> bool {
        planes[p]
            .n_out
            .dot(planes[q].n_out.cross(planes[r].n_out))
            .abs()
            > 0.5
    }

    /// Core: `orient3d` is rigid-rotation invariant, so the frame3 path over a rotated
    /// cuboid's exact `Pt3` definitions must agree, on every definite config, with the
    /// geom path over the same cuboid unrotated. Validates the ops-side assembly + routing
    /// (predicate soundness itself is `frame3`'s H-c).
    #[test]
    fn t_orient3d_rotation_invariant() {
        let (mut m, s) = cuboid();
        let pu = collect_planes(&m, s).unwrap();
        let r = rotated(&mut m, s);
        let pr = collect_planes(&m, r).unwrap();
        assert_eq!(pu.len(), pr.len(), "rotation preserves the face list");
        let n = pu.len();
        let mut checked = 0usize;
        for p in 0..n {
            for q in (p + 1)..n {
                for rr in (q + 1)..n {
                    if !normals_independent(&pu, p, q, rr) {
                        continue; // a parallel triple has no meeting point
                    }
                    for j in 0..n {
                        if j == p || j == q || j == rr {
                            continue;
                        }
                        let su = t_orient3d(&pu, p, q, rr, j, false);
                        if su == 0 {
                            continue; // indefinite — skip
                        }
                        let sr = t_orient3d(&pr, p, q, rr, j, true);
                        assert_eq!(su, sr, "rotation-invariant at ({p},{q},{rr},{j})");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 0, "no definite config in the corpus");
    }

    /// `collect_planes` on a rotated solid fills each plane's exact `Pt3` triple, whose
    /// coordinates match the stored `tri` and whose definition is rotated (tol > 0).
    #[test]
    fn plane_def_from_face() {
        let (mut m, s) = cuboid();
        let r = rotated(&mut m, s);
        let planes = collect_planes(&m, r).unwrap();
        for pi in &planes {
            let def = pi.tri_pt3.as_ref().expect("rotated plane carries tri_pt3");
            for (d, t) in def.iter().zip(pi.tri.iter()) {
                assert_eq!(d.coord, t.as_array(), "tri_pt3 matches tri");
            }
            assert!(
                def.iter().any(|p| p.tol.iter().any(|&t| t > 0.0)),
                "a rotated plane's def carries rotation tol"
            );
        }
        // The axis-aligned original carries no def.
        let pu = collect_planes(&m, s).unwrap();
        assert!(pu.iter().all(|pi| pi.tri_pt3.is_none()));
    }

    /// `rotated = false` forwards to the geom predicate bit-for-bit (axis-aligned hot path
    /// unchanged).
    #[test]
    fn t_orient3d_unrotated_forwards_geom() {
        let (m, s) = cuboid();
        let planes = collect_planes(&m, s).unwrap();
        let n = planes.len();
        for p in 0..n {
            for q in (p + 1)..n {
                for rr in (q + 1)..n {
                    if !normals_independent(&planes, p, q, rr) {
                        continue;
                    }
                    for j in 0..n {
                        if j == p || j == q || j == rr {
                            continue;
                        }
                        let want = three_plane_orient3d(
                            &planes[p].plane,
                            &planes[q].plane,
                            &planes[rr].plane,
                            planes[j].tri[0],
                            planes[j].tri[1],
                            planes[j].tri[2],
                        );
                        assert_eq!(t_orient3d(&planes, p, q, rr, j, false), want);
                    }
                }
            }
        }
    }
}
