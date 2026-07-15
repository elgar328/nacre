//! Toleranced boolean predicates (overhaul stage 3): the geom sign predicates routed
//! through TIP so they stay exact under rotation.
//!
//! A rotated face's plane coefficients and `tri` coordinates are rounded irrationals, so
//! the axis-aligned geom predicates (`nacre_geom::intersect`) are exact only w.r.t. the
//! *rounded* geometry. When the boolean's operands are rotated, these wrappers rebuild
//! each plane from the three exact `Pt3` its face carries ([`PlaneInfo::tri_pt3`], or —
//! for the axis-aligned operand of a *mixed*-rotation boolean, whose `tri_pt3` is `None` —
//! exactly from its `tri` coordinates, see [`plane_def`]) and decide the sign with the
//! `nacre-scalar::frame3` judges instead. The routing is a per-boolean `rotated` flag
//! ([`solid_is_rotated`](crate::solid_is_rotated)) — the axis-aligned hot path is unchanged
//! and never builds a `Pt3`.
//!
//! [`t_orient3d`] (order_along, 3a-i), [`t_cmp_coord`] (loop_winding) and [`t_plane_side`]
//! (straddle, 3a-ii) are here; the live arrangement still calls the geom predicates until
//! stage 3b wires these in.

use crate::PlaneInfo;
use nacre_geom::intersect::{plane_side, three_plane_cmp_coord, three_plane_orient3d};
use nacre_math::Point3;
use nacre_scalar::frame3::{
    Pt3, indirect_cmp_coord_judge, indirect_orient3d_judge, orient3d_judge,
};
use nacre_scalar::{Orient, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Vertex};

/// The three exact `Pt3` defining plane `k`. A rotated plane carries them cached in
/// [`PlaneInfo::tri_pt3`] (clone — a shallow copy of the rotation chain, no forest walk);
/// an axis-aligned plane (`None`, e.g. the unrotated operand of a *mixed*-rotation
/// boolean) is built exactly from its `tri` coordinates, which are already exact f64.
fn plane_def(planes: &[PlaneInfo], k: usize) -> [Pt3; 3] {
    match &planes[k].tri_pt3 {
        Some(t) => t.clone(),
        None => planes[k].tri.map(pt3_from_exact),
    }
}

/// An exact axis-aligned point as a tol-0 `Pt3` (its f64 coordinates are exact rationals).
fn pt3_from_exact(p: Point3) -> Pt3 {
    let a = p.as_array();
    let rat = |x: f64| Rat::try_from_f64(x).expect("axis-aligned coordinate is an exact rational");
    Pt3::at([rat(a[0]), rat(a[1]), rat(a[2])])
}

fn to_i8(o: Orient) -> i8 {
    match o {
        Orient::Positive => 1,
        Orient::Negative => -1,
        Orient::Zero => 0,
    }
}

/// Borrow an owned plane def as the `&Pt3` tuple the judges take.
fn borrow3(d: &[Pt3; 3]) -> (&Pt3, &Pt3, &Pt3) {
    (&d[0], &d[1], &d[2])
}

/// Borrow three owned plane defs as the tuples `indirect_cmp_coord_judge` takes.
fn borrow_triple(d: &[[Pt3; 3]; 3]) -> [(&Pt3, &Pt3, &Pt3); 3] {
    [borrow3(&d[0]), borrow3(&d[1]), borrow3(&d[2])]
}

/// The sign of `orient3d(V, tri_j)` where `V = ∩(planes p, q, r)` is an implicit point —
/// the toleranced twin of [`three_plane_orient3d`], matching `order_along`'s shape so it
/// is a drop-in (`+1`/`-1`/`0`).
///
/// - `!rotated`: the exact axis-aligned path — `three_plane_orient3d` on the stored plane
///   coefficients and `tri` coordinates (unchanged, fast).
/// - `rotated`: each of `p, q, r` and the explicit triangle `j` is taken as its three
///   exact `Pt3` ([`plane_def`]), and [`indirect_orient3d_judge`] decides the sign from
///   the definitions — never materializing `V` or reading the rounded `tri`.
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
    let (dp, dq, dr, dj) = (
        plane_def(planes, p),
        plane_def(planes, q),
        plane_def(planes, r),
        plane_def(planes, j),
    );
    to_i8(indirect_orient3d_judge(
        borrow3(&dp),
        borrow3(&dq),
        borrow3(&dr),
        &dj[0],
        &dj[1],
        &dj[2],
    ))
}

/// The sign of `a[axis] − b[axis]` between the two implicit points `a = ∩(planes a…)` and
/// `b = ∩(planes b…)` — the toleranced twin of [`three_plane_cmp_coord`] (loop_winding),
/// matching its shape (`+1` = `a[axis] > b[axis]`). `!rotated` → the geom predicate on the
/// stored coefficients; `rotated` → each triple's three planes as exact `Pt3` →
/// [`indirect_cmp_coord_judge`].
#[allow(dead_code)]
pub(crate) fn t_cmp_coord(
    planes: &[PlaneInfo],
    a: [usize; 3],
    b: [usize; 3],
    axis: usize,
    rotated: bool,
) -> i8 {
    if !rotated {
        let tri = |t: [usize; 3]| {
            [
                &planes[t[0]].plane,
                &planes[t[1]].plane,
                &planes[t[2]].plane,
            ]
        };
        return three_plane_cmp_coord(tri(a), tri(b), axis);
    }
    let da = a.map(|k| plane_def(planes, k));
    let db = b.map(|k| plane_def(planes, k));
    to_i8(indirect_cmp_coord_judge(
        borrow_triple(&da),
        borrow_triple(&db),
        axis,
    ))
}

/// Which side of plane `plane_idx` the explicit vertex `p` lies on — the toleranced twin
/// of [`plane_side`] (an all-explicit `orient3d(p, tri0, tri1, tri2)`). `!rotated` → the
/// geom predicate on `p`'s coordinate; `rotated` → `p` and the plane's three points as
/// exact `Pt3` → [`orient3d_judge`] (the direct toleranced orient3d).
#[allow(dead_code)]
pub(crate) fn t_plane_side(
    model: &Model,
    planes: &[PlaneInfo],
    plane_idx: usize,
    p: Handle<Vertex>,
    rotated: bool,
) -> i8 {
    if !rotated {
        return plane_side(planes[plane_idx].tri, model.vertices.get(p).point);
    }
    let pp = nacre_tip::vertex_pt3(model, p).expect("an original vertex is a direct point");
    let d = plane_def(planes, plane_idx);
    to_i8(orient3d_judge(&pp, &d[0], &d[1], &d[2]))
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

    fn boundary_verts(m: &Model, s: Handle<Solid>) -> Vec<Handle<nacre_topo::Vertex>> {
        let mut seen = std::collections::HashSet::new();
        let mut vs = Vec::new();
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh) {
                            vs.push(vh);
                        }
                    }
                }
            }
        }
        vs
    }

    /// The independent-normal plane triples of a cuboid (each meets at one corner).
    fn corner_triples(planes: &[PlaneInfo]) -> Vec<[usize; 3]> {
        let n = planes.len();
        let mut out = Vec::new();
        for p in 0..n {
            for q in (p + 1)..n {
                for r in (q + 1)..n {
                    if normals_independent(planes, p, q, r) {
                        out.push([p, q, r]);
                    }
                }
            }
        }
        out
    }

    /// `plane_side` is an orient3d sign → rigid-rotation invariant: the frame3 path over a
    /// rotated cuboid agrees with the geom path over the same cuboid unrotated.
    #[test]
    fn t_plane_side_rotation_invariant() {
        let (mut m, s) = cuboid();
        let pu = collect_planes(&m, s).unwrap();
        let vu = boundary_verts(&m, s);
        let r = rotated(&mut m, s);
        let pr = collect_planes(&m, r).unwrap();
        let vr = boundary_verts(&m, r);
        assert_eq!((pu.len(), vu.len()), (pr.len(), vr.len()));
        let mut checked = 0usize;
        for pi in 0..pu.len() {
            for vk in 0..vu.len() {
                let su = t_plane_side(&m, &pu, pi, vu[vk], false);
                if su == 0 {
                    continue; // vertex on the plane — skip
                }
                let sr = t_plane_side(&m, &pr, pi, vr[vk], true);
                assert_eq!(
                    su, sr,
                    "plane_side rotation-invariant at plane {pi}, vert {vk}"
                );
                checked += 1;
            }
        }
        assert!(checked > 0, "no definite side test in the corpus");
    }

    /// `rotated = false` forwards to `plane_side` bit-for-bit.
    #[test]
    fn t_plane_side_unrotated_forwards_geom() {
        let (m, s) = cuboid();
        let planes = collect_planes(&m, s).unwrap();
        let vs = boundary_verts(&m, s);
        for pi in 0..planes.len() {
            for &v in &vs {
                let want = plane_side(planes[pi].tri, m.vertices.get(v).point);
                assert_eq!(t_plane_side(&m, &planes, pi, v, false), want);
            }
        }
    }

    /// `t_cmp_coord` orders two implicit points (corner triples) by an axis exactly as
    /// their `three_planes` coordinates do (cmp is not rotation-invariant, so the f64
    /// coordinate is the oracle) — never a wrong sign, resolves some, antisymmetric.
    #[test]
    fn t_cmp_coord_matches_coord() {
        use nacre_geom::intersect::three_planes;
        let (mut m, s) = cuboid();
        let r = rotated(&mut m, s);
        let planes = collect_planes(&m, r).unwrap();
        let triples = corner_triples(&planes);
        let meet = |t: [usize; 3]| {
            three_planes(
                &planes[t[0]].plane,
                &planes[t[1]].plane,
                &planes[t[2]].plane,
            )
        };
        let (mut checked, mut resolved) = (0usize, 0usize);
        for i in 0..triples.len() {
            for jx in (i + 1)..triples.len() {
                let (a, b) = (triples[i], triples[jx]);
                let (Some(va), Some(vb)) = (meet(a), meet(b)) else {
                    continue;
                };
                for axis in 0..3 {
                    let diff = va.as_array()[axis] - vb.as_array()[axis];
                    if diff.abs() < 1e-6 {
                        continue; // near-tie: f64 oracle unreliable
                    }
                    checked += 1;
                    let want = if diff > 0.0 { 1 } else { -1 };
                    let got = t_cmp_coord(&planes, a, b, axis, true);
                    assert!(got == want || got == 0, "cmp {got} vs coord {want}");
                    if got == want {
                        resolved += 1;
                    }
                    assert_eq!(
                        t_cmp_coord(&planes, b, a, axis, true),
                        -got,
                        "cmp is antisymmetric"
                    );
                }
            }
        }
        assert!(checked > 0, "no definite pair");
        assert!(resolved > 0, "bridge must resolve some orderings");
    }

    /// `rotated = false` forwards to `three_plane_cmp_coord` bit-for-bit.
    #[test]
    fn t_cmp_coord_unrotated_forwards_geom() {
        let (m, s) = cuboid();
        let planes = collect_planes(&m, s).unwrap();
        let triples = corner_triples(&planes);
        for i in 0..triples.len() {
            for jx in (i + 1)..triples.len() {
                let (a, b) = (triples[i], triples[jx]);
                for axis in 0..3 {
                    let tri = |t: [usize; 3]| {
                        [
                            &planes[t[0]].plane,
                            &planes[t[1]].plane,
                            &planes[t[2]].plane,
                        ]
                    };
                    let want = three_plane_cmp_coord(tri(a), tri(b), axis);
                    assert_eq!(t_cmp_coord(&planes, a, b, axis, false), want);
                }
            }
        }
    }

    /// A mixed-rotation boolean (one rotated operand, one axis-aligned) leaves the
    /// axis-aligned operand's `tri_pt3 = None`; the wrappers build it exactly from the tri
    /// coordinates on demand, so a predicate spanning both never panics.
    #[test]
    fn mixed_rotation_handled() {
        let mut m = Model::new();
        let a0 = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([12.0, 13.0, 14.0]),
        );
        let a = rotated(&mut m, a0);
        let mut planes = collect_planes(&m, a).unwrap();
        let na = planes.len();
        let pb = collect_planes(&m, b).unwrap();
        assert!(
            planes.iter().all(|p| p.tri_pt3.is_some()),
            "rotated operand cached"
        );
        assert!(
            pb.iter().all(|p| p.tri_pt3.is_none()),
            "axis-aligned operand not cached"
        );
        planes.extend(pb);
        // plane_def on an axis-aligned (None) plane builds an exact, tol-0 def.
        let bdef = plane_def(&planes, na);
        assert!(
            bdef.iter().all(|d| d.tol == [0.0; 3]),
            "axis-aligned def is tol 0"
        );
        // A predicate spanning the rotated operand (Some) and axis-aligned operand (None)
        // must not panic — plane_def resolves the None plane on demand.
        let mut exercised = false;
        for p in 0..na {
            for q in (p + 1)..na {
                for r in (q + 1)..na {
                    if normals_independent(&planes, p, q, r) {
                        let _ = t_orient3d(&planes, p, q, r, na, true); // j = first B plane
                        exercised = true;
                    }
                }
            }
        }
        assert!(exercised, "no mixed predicate exercised");
    }
}
