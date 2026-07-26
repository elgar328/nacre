//! The b-rep side of the toleranced predicates: `nacre-ops`'s arrangement tables
//! ([`PlaneGeom`], [`FaceInfo`]) implement the [`Witness`]/[`PlaneWitness`] ports so the
//! rotation-general sign predicates in [`nacre_cip::predicate`] can run over them, and the
//! `t_*` predicates are re-exported here so existing `crate::tolerant::t_*` call sites are
//! unchanged. The predicate logic itself (exact-vs-kernel routing, the frame3 judges) lives in
//! `nacre-cip`.

use crate::planes::{FaceInfo, PlaneGeom};
use nacre_cip::Pt3;
use nacre_cip::predicate::{PlaneWitness, Witness};
use nacre_math::Point3;

// Re-export the toleranced predicates from cip so existing call sites keep working.
pub(crate) use nacre_cip::predicate::{
    any_rotated, plane_def, t_cmp_coord, t_orient3d, t_plane_pair_dir_sign, t_planes_coplanar,
};

impl Witness for PlaneGeom {
    fn tri(&self) -> [Point3; 3] {
        self.tri
    }
    fn tri_pt3(&self) -> &[Pt3; 3] {
        &self.tri_pt3
    }
    fn is_rotated(&self) -> bool {
        self.rotated
    }
}

impl Witness for FaceInfo {
    fn tri(&self) -> [Point3; 3] {
        self.tri
    }
    fn tri_pt3(&self) -> &[Pt3; 3] {
        &self.tri_pt3
    }
    fn is_rotated(&self) -> bool {
        self.rotated
    }
}

impl PlaneWitness for PlaneGeom {
    fn coeffs(&self) -> [f64; 4] {
        self.plane.coefficients()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planes::collect_planes;
    use crate::{Operation, apply};
    // The exact geom predicates, used here as independent oracles for the `t_*` wrappers.
    use nacre_geom::intersect::{plane_pair_dir_sign, three_plane_cmp_coord, three_plane_orient3d};

    /// The plane table of one solid, built through the real path so these tests exercise the same
    /// `PlaneGeom` the engine does. A single convex operand has no coplanar pair, so the numbering
    /// is the identity — indices below name a face and its plane interchangeably.
    fn plane_table(m: &Model, s: Handle<Solid>) -> Vec<PlaneGeom> {
        let faces = collect_planes(m, s).unwrap();
        let canon = crate::planes::plane_classes(&faces);
        crate::planes::dense_planes(&faces, &canon).0
    }
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
    fn normals_independent(planes: &[PlaneGeom], p: usize, q: usize, r: usize) -> bool {
        let n = |k: usize| planes[k].plane.normal();
        n(p).dot(n(q).cross(n(r))).abs() > 0.5
    }

    /// Core: `orient3d` is rigid-rotation invariant, so the frame3 path over a rotated
    /// cuboid's exact `Pt3` definitions must agree, on every definite config, with the
    /// geom path over the same cuboid unrotated. Validates the ops-side assembly + routing
    /// (predicate soundness itself is `frame3`'s H-c).
    #[test]
    fn t_orient3d_rotation_invariant() {
        let (mut m, s) = cuboid();
        let pu = plane_table(&m, s);
        let r = rotated(&mut m, s);
        let pr = plane_table(&m, r);
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
        let planes = plane_table(&m, r);
        for pi in &planes {
            assert!(pi.rotated, "a rotated solid's planes are flagged rotated");
            let def = &pi.tri_pt3;
            for (d, t) in def.iter().zip(pi.tri.iter()) {
                assert_eq!(d.coord, t.as_array(), "tri_pt3 matches tri");
            }
            assert!(
                def.iter().any(|p| p.tol.iter().any(|&t| t > 0.0)),
                "a rotated plane's def carries rotation tol"
            );
        }
        // The axis-aligned original also carries a def now — that is the change this test used
        // to forbid. What distinguishes it is the flag, not the presence of a definition, and its
        // def is exact: tol 0 and coordinates equal to `tri`.
        let pu = plane_table(&m, s);
        assert!(pu.iter().all(|pi| !pi.rotated));
        for pi in &pu {
            for (d, t) in pi.tri_pt3.iter().zip(pi.tri.iter()) {
                assert_eq!(d.coord, t.as_array(), "axis-aligned def matches tri");
                assert_eq!(d.tol, [0.0; 3], "axis-aligned def is exact");
            }
        }
    }

    /// `rotated = false` forwards to the geom predicate bit-for-bit (axis-aligned hot path
    /// unchanged).
    #[test]
    fn t_orient3d_unrotated_forwards_geom() {
        let (m, s) = cuboid();
        let planes = plane_table(&m, s);
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

    /// The independent-normal plane triples of a cuboid (each meets at one corner).
    fn corner_triples(planes: &[PlaneGeom]) -> Vec<[usize; 3]> {
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

    /// `t_planes_coplanar` decides plane identity from the faces' own coordinates, so it must
    /// (a) refuse to conclude anything from a degenerate `tri` — three equal points lie on *every*
    /// plane, and merging on that evidence would fuse genuinely different planes — and
    /// (b) answer the same for a rotated solid as for the unrotated one.
    #[test]
    fn t_planes_coplanar_guards_degeneracy_and_survives_rotation() {
        let (mut m, s) = cuboid();
        let pu = plane_table(&m, s);
        // (a) A hand-built degenerate pair: same-normal parallel planes, but `tri` is a point.
        let degenerate: Vec<PlaneGeom> = (0..2)
            .map(|k| PlaneGeom {
                surf: pu[0].surf,
                plane: pu[0].plane,
                tri: [Point3::from_array([k as f64, 0.0, 0.0]); 3],
                tri_pt3: std::array::from_fn(|_| Pt3::exact([k as f64, 0.0, 0.0]).expect("exact")),
                rotated: false,
                frame_sign: pu[0].frame_sign,
            })
            .collect();
        assert!(
            !t_planes_coplanar(&degenerate, 0, 1),
            "a degenerate tri is no evidence"
        );
        // (b) Rotation invariance over every pair of the cuboid's faces.
        let r = rotated(&mut m, s);
        let pr = plane_table(&m, r);
        assert_eq!(pu.len(), pr.len());
        let mut same = 0usize;
        for i in 0..pu.len() {
            for j in (i + 1)..pu.len() {
                assert_eq!(
                    t_planes_coplanar(&pu, i, j),
                    t_planes_coplanar(&pr, i, j),
                    "rotation-invariant at ({i},{j})"
                );
                if t_planes_coplanar(&pu, i, j) {
                    same += 1;
                }
            }
        }
        // A cuboid has six distinct planes, so no pair may merge — the answer is not vacuously true.
        assert_eq!(same, 0, "a cuboid has no two coplanar faces");
    }

    /// `t_cmp_coord` orders two implicit points (corner triples) by an axis exactly as
    /// their `three_planes` coordinates do (cmp is not rotation-invariant, so the f64
    /// coordinate is the oracle) — never a wrong sign, resolves some, antisymmetric.
    #[test]
    fn t_cmp_coord_matches_coord() {
        use nacre_geom::intersect::three_planes;
        let (mut m, s) = cuboid();
        let r = rotated(&mut m, s);
        let planes = plane_table(&m, r);
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
        let planes = plane_table(&m, s);
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
        let mut planes = plane_table(&m, a);
        let na = planes.len();
        let pb = plane_table(&m, b);
        assert!(planes.iter().all(|p| p.rotated), "rotated operand flagged");
        assert!(
            pb.iter().all(|p| !p.rotated),
            "axis-aligned operand not flagged"
        );
        planes.extend(pb);
        // The axis-aligned operand's stored def is exact (tol 0) — it is cached like the
        // rotated one, and only the flag tells them apart.
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
        let pu = plane_table(&m, s);
        let r = rotated(&mut m, s);
        let pr = plane_table(&m, r);
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
        let planes = plane_table(&m, s);
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
