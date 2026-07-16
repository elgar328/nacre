//! Toleranced boolean predicates (overhaul stage 3): the geom sign predicates routed
//! through TIP so they stay exact under rotation.
//!
//! A rotated face's plane coefficients and `tri` coordinates are rounded irrationals, so
//! the axis-aligned geom predicates (`nacre_geom::intersect`) are exact only w.r.t. the
//! *rounded* geometry. When a predicate's planes are rotated, these wrappers rebuild each
//! plane from the three exact `Pt3` its face carries ([`PlaneInfo::tri_pt3`], or — for the
//! axis-aligned operand of a *mixed*-rotation boolean, whose `tri_pt3` is `None` — exactly
//! from its `tri` coordinates, see [`plane_def`]) and decide the sign with the
//! `nacre-scalar::frame3` judges instead.
//!
//! **Routing is per-predicate, derived — no `rotated` flag is threaded.** Each wrapper asks
//! [`any_rotated`] of just the planes it touches: all-axis-aligned → the exact geom hot path
//! (never builds a `Pt3`); any rotated → the frame3 judge. So a mixed-rotation boolean keeps
//! the axis-aligned operand's own predicates on the fast path, finer than a per-boolean flag.
//!
//! All four are wired into the live arrangement: [`t_orient3d`] (order_along, side_of),
//! [`t_cmp_coord`] (loop_winding) and [`t_plane_pair_dir_sign`] (dir_sign, turn_at,
//! point_on_ring, every_ray) in stage 3b-i; [`t_plane_side`] (edge_crosses_face straddle,
//! with a vertex handle + `model` threaded down) in 3b-ii.

use crate::PlaneInfo;
use crate::arrange::orient_sign;
use nacre_geom::intersect::{
    plane_pair_dir_sign, plane_side, three_plane_cmp_coord, three_plane_orient3d,
};
use nacre_math::Point3;
use nacre_scalar::frame3::{
    Pt3, dir_sign_judge, indirect_cmp_coord_judge, indirect_orient3d_judge, orient3d_judge,
};
use nacre_scalar::{Orient, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Origin, Vertex};

/// Whether any of the named planes is rotated — the per-predicate routing signal.
///
/// A plane carries [`PlaneInfo::tri_pt3`] `Some` exactly when it came from a rotated solid
/// (`collect_planes` fills it via [`solid_is_rotated`](crate::solid_is_rotated), whose truth
/// is owned by the classifier: `transform` marks a vertex `Origin::Rotated` only when
/// `Isometry::is_exact` is false, i.e. the realization is irrational). A predicate must
/// escalate to frame3 if **any** — not all — of its planes is irrational: a single rounded
/// coordinate can flip an f64 `orient3d`/`cmp`, whereas all-rational planes are exact on the
/// geom path. Widening what counts as exact (e.g. more angle families) is a classifier-layer
/// change; this consumer only reads the flag.
fn any_rotated(planes: &[PlaneInfo], idx: &[usize]) -> bool {
    idx.iter().any(|&k| planes[k].tri_pt3.is_some())
}

/// The three exact `Pt3` defining plane `k`. A rotated plane carries them cached in
/// [`PlaneInfo::tri_pt3`] (clone — a shallow copy of the rotation chain, no forest walk);
/// an axis-aligned plane (`None`, e.g. the unrotated operand of a *mixed*-rotation
/// boolean) is built exactly from its `tri` coordinates, which are already exact f64.
pub(crate) fn plane_def(planes: &[PlaneInfo], k: usize) -> [Pt3; 3] {
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
///
/// Routes on [`any_rotated`] of `p, q, r, j`: all-axis-aligned → geom, any rotated → frame3.
pub(crate) fn t_orient3d(planes: &[PlaneInfo], p: usize, q: usize, r: usize, j: usize) -> i8 {
    if !any_rotated(planes, &[p, q, r, j]) {
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
/// [`indirect_cmp_coord_judge`]. Routes on [`any_rotated`] of the six planes in `a` and `b`.
pub(crate) fn t_cmp_coord(planes: &[PlaneInfo], a: [usize; 3], b: [usize; 3], axis: usize) -> i8 {
    if !any_rotated(planes, &[a[0], a[1], a[2], b[0], b[1], b[2]]) {
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
///
/// Routes on the plane *or* the vertex being rotated: unlike the plane-only predicates, an
/// axis-aligned plane may still be tested against a rotated (irrational) vertex, so `p`'s
/// `Origin::Rotated` also forces frame3.
///
/// `p` is an edge endpoint — an original (`Constructed`/`Rotated`) vertex. The rotated path's
/// `vertex_pt3` errors on a `Discovered` (three-plane) definition; a rotated boolean whose
/// operand is a prior boolean result carries such corners, so 3c must route those endpoints
/// through the indirect `orient3d` (nacre-tip dispatch) before the rotation guard retires.
pub(crate) fn t_plane_side(
    model: &Model,
    planes: &[PlaneInfo],
    plane_idx: usize,
    p: Handle<Vertex>,
) -> i8 {
    let rotated = any_rotated(planes, &[plane_idx])
        || matches!(model.vertices.get(p).origin, Origin::Rotated { .. });
    if !rotated {
        return plane_side(planes[plane_idx].tri, model.vertices.get(p).point);
    }
    let pp = nacre_tip::vertex_pt3(model, p).expect("an original vertex is a direct point");
    let d = plane_def(planes, plane_idx);
    to_i8(orient3d_judge(&pp, &d[0], &d[1], &d[2]))
}

/// `sign(det[n_p; n_a; n_b])` over the three planes' stored normals — the toleranced twin
/// of [`plane_pair_dir_sign`] (how the line `p ∩ a` runs relative to plane `b`), matching
/// its shape (index-based drop-in, `+1`/`-1`/`0`).
///
/// `!rotated` → the geom predicate. `rotated` → the frame3 `D` (det of the *outward*
/// `tri` normals, [`dir_sign_judge`]) bridged to the *stored*-normal convention by the
/// per-plane [`orient_sign`]: `det(stored) = orient_sign(p)·orient_sign(a)·orient_sign(b)·
/// det(outward)`. `orient_sign` is an f64 dot of two parallel unit vectors (`|·| ≈ 1`),
/// robust under rotation. Routes on [`any_rotated`] of `p, a, b`.
pub(crate) fn t_plane_pair_dir_sign(planes: &[PlaneInfo], p: usize, a: usize, b: usize) -> i8 {
    if !any_rotated(planes, &[p, a, b]) {
        return plane_pair_dir_sign(&planes[p].plane, &planes[a].plane, &planes[b].plane);
    }
    let (dp, da, db) = (
        plane_def(planes, p),
        plane_def(planes, a),
        plane_def(planes, b),
    );
    orient_sign(planes, p)
        * orient_sign(planes, a)
        * orient_sign(planes, b)
        * to_i8(dir_sign_judge(borrow3(&dp), borrow3(&da), borrow3(&db)))
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
                        let su = t_orient3d(&pu, p, q, rr, j);
                        if su == 0 {
                            continue; // indefinite — skip
                        }
                        let sr = t_orient3d(&pr, p, q, rr, j);
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
                        assert_eq!(t_orient3d(&planes, p, q, rr, j), want);
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
                let su = t_plane_side(&m, &pu, pi, vu[vk]);
                if su == 0 {
                    continue; // vertex on the plane — skip
                }
                let sr = t_plane_side(&m, &pr, pi, vr[vk]);
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
                assert_eq!(t_plane_side(&m, &planes, pi, v), want);
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
                    let got = t_cmp_coord(&planes, a, b, axis);
                    assert!(got == want || got == 0, "cmp {got} vs coord {want}");
                    if got == want {
                        resolved += 1;
                    }
                    assert_eq!(
                        t_cmp_coord(&planes, b, a, axis),
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
                    assert_eq!(t_cmp_coord(&planes, a, b, axis), want);
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
                        let _ = t_orient3d(&planes, p, q, r, na); // j = first B plane
                        exercised = true;
                    }
                }
            }
        }
        assert!(exercised, "no mixed predicate exercised");
    }

    /// `plane_pair_dir_sign` is a determinant of normals → rigid-rotation invariant: the
    /// frame3 path (`orient_sign` · `D`) over a rotated cuboid agrees with the geom path
    /// over the same cuboid unrotated, on every definite ordered triple.
    #[test]
    fn t_dir_sign_rotation_invariant() {
        let (mut m, s) = cuboid();
        let pu = collect_planes(&m, s).unwrap();
        let r = rotated(&mut m, s);
        let pr = collect_planes(&m, r).unwrap();
        let n = pu.len();
        let mut checked = 0usize;
        for p in 0..n {
            for a in 0..n {
                for b in 0..n {
                    if p == a || p == b || a == b {
                        continue;
                    }
                    let su = t_plane_pair_dir_sign(&pu, p, a, b);
                    if su == 0 {
                        continue; // coplanar normals — skip
                    }
                    let sr = t_plane_pair_dir_sign(&pr, p, a, b);
                    assert_eq!(su, sr, "dir_sign rotation-invariant at ({p},{a},{b})");
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no definite triple in the corpus");
    }

    /// `rotated = false` forwards to `plane_pair_dir_sign` bit-for-bit.
    #[test]
    fn t_dir_sign_unrotated_forwards_geom() {
        let (m, s) = cuboid();
        let planes = collect_planes(&m, s).unwrap();
        let n = planes.len();
        for p in 0..n {
            for a in 0..n {
                for b in 0..n {
                    if p == a || p == b || a == b {
                        continue;
                    }
                    let want =
                        plane_pair_dir_sign(&planes[p].plane, &planes[a].plane, &planes[b].plane);
                    assert_eq!(t_plane_pair_dir_sign(&planes, p, a, b), want);
                }
            }
        }
    }
}
