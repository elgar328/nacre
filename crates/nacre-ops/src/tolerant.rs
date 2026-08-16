//! The b-rep side of the toleranced predicates: `nacre-ops`'s arrangement tables
//! ([`WorkingPlane`], [`FaceInfo`]) implement the [`Witness`]/[`PlaneWitness`] ports so the
//! rotation-general sign predicates in [`nacre_cip::predicate`] can run over them.
//!
//! The tables are **pure description** — geometry and provenance, nothing about how this
//! operation judges. That belongs to [`Judge`], which the engine builds once per boolean
//! (`plane_index_setup` supplies the standard and the collector) and passes down; the predicates
//! are its methods. The predicate logic itself (exact-vs-kernel routing, the frame3 judges) lives
//! in `nacre-cip`.

use crate::planes::{FaceInfo, WorkingPlane};
use nacre_cip::WitnessPoint;
use nacre_cip::predicate::{PlaneWitness, Witness};
use nacre_math::Point3;

// The judging context and the helpers the engine reaches for by name. The predicates themselves
// are methods on `Judge`, so there is nothing else to re-export.
#[cfg(test)]
use nacre_cip::predicate::plane_def;
pub(crate) use nacre_cip::predicate::{ImplicitPoint, Judge};

impl Witness for WorkingPlane {
    fn base_coeffs_rat(&self) -> Option<[nacre_scalar::Rat; 4]> {
        self.base_rat
    }
    fn tri(&self) -> [Point3; 3] {
        self.tri
    }
    fn chain_id(&self) -> u64 {
        self.base.chain_id
    }
    fn base_tri(&self) -> Option<[Point3; 3]> {
        self.base.tri
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.tri_pt3
    }
    fn is_rotated(&self) -> bool {
        self.rotated
    }
}

impl Witness for FaceInfo {
    fn base_coeffs_rat(&self) -> Option<[nacre_scalar::Rat; 4]> {
        self.base_rat
    }
    fn tri(&self) -> [Point3; 3] {
        self.tri
    }
    // A face table exists before plane classes do, and the only predicate that runs on it is
    // `Judge::planes_coplanar` during class discovery. Opting out here keeps that path unchanged.
    fn chain_id(&self) -> u64 {
        0
    }
    fn base_tri(&self) -> Option<[Point3; 3]> {
        None
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.tri_pt3
    }
    fn is_rotated(&self) -> bool {
        self.rotated
    }
}

impl PlaneWitness for WorkingPlane {
    fn coeffs(&self) -> [f64; 4] {
        self.plane.coefficients()
    }
    fn exact_coeffs(&self) -> Option<[f64; 4]> {
        self.exact_coeffs
    }
    fn exact_normal(&self) -> Option<[f64; 3]> {
        self.exact_normal
    }
    fn frame_sign(&self) -> i8 {
        self.frame_sign
    }
    fn base_coeffs(&self) -> Option<[f64; 4]> {
        self.base.coeffs
    }
    fn name_ints(&self) -> Option<&nacre_cip::predicate::NameInts> {
        self.name_ints.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
    /// when the plane is not one the model already holds (a seed, or a face's).
    fn datum_frame(m: &mut Model, plane: crate::SketchPlane) -> crate::SketchFrame {
        match crate::apply(
            m,
            &crate::Operation::DatumPlane {
                def: crate::DatumDef::Stated(plane),
            },
        ) {
            Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
            other => panic!("stating a plane: {other:?}"),
        }
    }
    use crate::planes::collect_planes;
    use crate::{Operation, apply};
    // The exact geom predicates, used here as independent oracles for the `t_*` wrappers.
    use nacre_geom::intersect::{plane_pair_dir_sign, three_plane_cmp_coord, three_plane_orient3d};

    /// The plane table of one solid, built through the real path so these tests exercise the same
    /// `WorkingPlane` the engine does. A single convex operand has no coplanar pair, so the numbering
    /// is the identity — indices below name a face and its plane interchangeably.
    fn plane_table(m: &Model, s: Handle<Solid>) -> Vec<WorkingPlane> {
        let faces = collect_planes(m, s).unwrap();
        let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces));
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
    fn normals_independent(planes: &[WorkingPlane], p: usize, q: usize, r: usize) -> bool {
        let n = |k: usize| planes[k].plane.normal();
        n(p).dot(n(q).cross(n(r))).abs() > 0.5
    }

    /// Core: `orient3d` is rigid-rotation invariant, so the frame3 path over a rotated
    /// cuboid's exact `WitnessPoint` definitions must agree, on every definite config, with the
    /// geom path over the same cuboid unrotated. Validates the ops-side assembly + routing
    /// (predicate soundness itself is `frame3`'s H-c).
    #[test]
    fn orient3d_is_rotation_invariant() {
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
                        let su = crate::planes::test_judge(&pu).orient3d(p, q, rr, j);
                        if su == 0 {
                            continue; // indefinite — skip
                        }
                        let sr = crate::planes::test_judge(&pr).orient3d(p, q, rr, j);
                        assert_eq!(su, sr, "rotation-invariant at ({p},{q},{rr},{j})");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 0, "no definite config in the corpus");
    }

    /// `collect_planes` on a rotated solid fills each plane's exact `WitnessPoint` triple, whose
    /// coordinates **span the face's own plane** and whose definition is rotated (tol > 0).
    ///
    /// ★ The witness is the surface's recorded triple, which is the *first pusher's* — since S9
    /// an origin cuboid's axis faces intern onto the seeded world planes, so the witness is the
    /// seed's `[0, u, v]` (rotated along), not this face's corners. "Matches `tri`" was a
    /// first-push coincidence, never the contract; on-plane is.
    #[test]
    fn plane_def_from_face() {
        let on_plane = |pi: &crate::planes::WorkingPlane| {
            for d in &pi.tri_pt3 {
                let dist = pi.plane.distance(Point3::from_array(d.coord));
                assert!(dist < 1e-9, "a witness point sits {dist:e} off its plane");
            }
        };
        let (mut m, s) = cuboid();
        let r = rotated(&mut m, s);
        let planes = plane_table(&m, r);
        for pi in &planes {
            assert!(pi.rotated, "a rotated solid's planes are flagged rotated");
            on_plane(pi);
            assert!(
                pi.tri_pt3.iter().any(|p| p.tol.iter().any(|&t| t > 0.0)),
                "a rotated plane's def carries rotation tol"
            );
        }
        // The axis-aligned original also carries a def now — that is the change this test used
        // to forbid. What distinguishes it is the flag, not the presence of a definition, and
        // its def is exact: tol 0, on the plane.
        let pu = plane_table(&m, s);
        assert!(pu.iter().all(|pi| !pi.rotated));
        for pi in &pu {
            on_plane(pi);
            for d in &pi.tri_pt3 {
                assert_eq!(d.tol, [0.0; 3], "axis-aligned def is exact");
            }
        }
    }

    /// `rotated = false` forwards to the geom predicate bit-for-bit (axis-aligned hot path
    /// unchanged).
    #[test]
    fn orient3d_unrotated_forwards_geom() {
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
                        assert_eq!(
                            crate::planes::test_judge(&planes).orient3d(p, q, rr, j),
                            want
                        );
                    }
                }
            }
        }
    }

    /// The independent-normal plane triples of a cuboid (each meets at one corner).
    fn corner_triples(planes: &[WorkingPlane]) -> Vec<[usize; 3]> {
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

    /// `Judge::planes_coplanar` decides plane identity from the faces' own coordinates, so it must
    /// (a) refuse to conclude anything from a degenerate `tri` — three equal points lie on *every*
    /// plane, and merging on that evidence would fuse genuinely different planes — and
    /// (b) answer the same for a rotated solid as for the unrotated one.
    #[test]
    fn planes_coplanar_guards_degeneracy_and_survives_rotation() {
        let (mut m, s) = cuboid();
        let pu = plane_table(&m, s);
        // (a) A hand-built degenerate pair: same-normal parallel planes, but `tri` is a point.
        let degenerate: Vec<WorkingPlane> = (0..2)
            .map(|k| WorkingPlane {
                // A hand-built table has no recorded coefficients; the composed-rotation route
                // declines and the fixture takes the same escalating path it always did.
                base_rat: None,
                name_ints: None,
                base: crate::planes::BaseFrame::none(),
                surf: pu[0].surf,
                plane: pu[0].plane,
                tri: [Point3::from_array([k as f64, 0.0, 0.0]); 3],
                tri_pt3: std::array::from_fn(|_| {
                    WitnessPoint::exact([k as f64, 0.0, 0.0]).expect("exact")
                }),
                rotated: false,
                frame_sign: pu[0].frame_sign,
                // Derived by the same rule the arrangement uses -- a fixture that routed
                // differently would be testing a different engine.
                exact_coeffs: WorkingPlane::reconcile(
                    &pu[0].plane,
                    [Point3::from_array([k as f64, 0.0, 0.0]); 3],
                    false,
                )
                .0,
                exact_normal: WorkingPlane::reconcile(
                    &pu[0].plane,
                    [Point3::from_array([k as f64, 0.0, 0.0]); 3],
                    false,
                )
                .1,
            })
            .collect();
        assert!(
            !crate::planes::test_judge(&degenerate).planes_coplanar(0, 1),
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
                    crate::planes::test_judge(&pu).planes_coplanar(i, j),
                    crate::planes::test_judge(&pr).planes_coplanar(i, j),
                    "rotation-invariant at ({i},{j})"
                );
                if crate::planes::test_judge(&pu).planes_coplanar(i, j) {
                    same += 1;
                }
            }
        }
        // A cuboid has six distinct planes, so no pair may merge — the answer is not vacuously true.
        assert_eq!(same, 0, "a cuboid has no two coplanar faces");
    }

    /// `Judge::cmp_coord` orders two implicit points (corner triples) by an axis exactly as
    /// their `three_planes` coordinates do (cmp is not rotation-invariant, so the f64
    /// coordinate is the oracle) — never a wrong sign, resolves some, antisymmetric.
    #[test]
    fn cmp_coord_matches_coord() {
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
                    let got = crate::planes::test_judge(&planes).cmp_coord(a, b, axis);
                    assert!(got == want || got == 0, "cmp {got} vs coord {want}");
                    if got == want {
                        resolved += 1;
                    }
                    assert_eq!(
                        crate::planes::test_judge(&planes).cmp_coord(b, a, axis),
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
    fn cmp_coord_unrotated_forwards_geom() {
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
                    assert_eq!(
                        crate::planes::test_judge(&planes).cmp_coord(a, b, axis),
                        want
                    );
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
                        let _ = crate::planes::test_judge(&planes).orient3d(p, q, r, na); // j = first B plane
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
                    let su = crate::planes::test_judge(&pu).plane_pair_dir_sign(p, a, b);
                    if su == 0 {
                        continue; // coplanar normals — skip
                    }
                    let sr = crate::planes::test_judge(&pr).plane_pair_dir_sign(p, a, b);
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
                    assert_eq!(
                        crate::planes::test_judge(&planes).plane_pair_dir_sign(p, a, b),
                        want
                    );
                }
            }
        }
    }

    /// The same geometry told two ways: as a constructed solid, and as the reflection of its
    /// own mirror image. A reflection in `x = 1` sends an integer to an integer, so the second
    /// spelling reproduces the first **bit for bit** — same triangles, same planes,
    /// same everything a predicate reads directly. What differs is only the *provenance*: one
    /// table carries no motion (`chain_id = 0`, so every shortcut declines and the predicates run
    /// on the coordinates), the other carries a reflection (`chain_id ≠ 0`, so the shortcuts fire
    /// and answer in the pre-motion frame).
    ///
    /// **That is the whole point.** A reflection is improper, and the shortcuts' licence is that
    /// a motion preserves the determinant they take — which an improper one does not. The base
    /// frame is canonicalised for that (`BaseFrame::of`), and this test is what says the
    /// canonicalisation is right: it puts the shortcut path and the no-shortcut path side by
    /// side on identical geometry, with no bypass switch in between. A sign that survives one
    /// too many or one too few reflections shows up here as a disagreement.
    fn two_spellings() -> (Vec<WorkingPlane>, Vec<WorkingPlane>) {
        use crate::planes::{FaceInfo, dense_planes, plane_classes};
        use nacre_topo::Motion;

        let (mut m, s) = cuboid();
        // A **non-zero** mirror plane on purpose: with `offset = 0` the base frame would come out
        // equal to the moved one and every comparison below would pass for free. At `x = 1` the
        // canonicalised base is the moved triangle displaced by `(−2, 0, 0)` — a different frame,
        // still exact.
        let leaf = m.push_motion(
            Motion::Mirror {
                axis: Axis::X,
                offset: Rat::from_int(1),
            },
            None,
        );
        let plain = collect_planes(&m, s).unwrap();
        // The same faces, described as the image of their own reflection: base = `−coord`, one
        // `Mirror` node, and the realized coordinate lands back on `coord`.
        let mirrored: Vec<FaceInfo> = plain
            .iter()
            .map(|f| FaceInfo {
                // A hand-built table has no surface to have recorded coefficients on, so the base
                // frame derives as it always did. Filling this by hand is how a fixture and the
                // engine come to route differently (see `WorkingPlane::reconcile`).
                base_rat: None,
                name: None,
                surf: f.surf,
                face: f.face,
                plane: f.plane,
                tri: f.tri,
                n_out: f.n_out,
                orient_sign: f.orient_sign,
                tri_pt3: std::array::from_fn(|i| {
                    let c = f.tri_pt3[i].coord;
                    let r = |v: f64| Rat::try_from_f64(v).expect("integer cuboid coordinate");
                    WitnessPoint::at([r(2.0 - c[0]), r(c[1]), r(c[2])])
                        .mirror(Axis::X, Rat::from_int(1))
                }),
                motion: Some(leaf),
                rotated: true,
            })
            .collect();
        for (a, b) in plain.iter().zip(&mirrored) {
            for i in 0..3 {
                assert_eq!(a.tri_pt3[i].coord, b.tri_pt3[i].coord, "same coordinates");
                // The *operation* is exact on these integers, but `tol` is an a-priori bound
                // (`ε·(|2c| + |x|)`), not the realized error — so it is a couple of ulp, not
                // zero. Sound, and it costs nothing here: the base frame is chosen on the
                // rational `base`, not on `tol`.
                for k in 0..3 {
                    let t = b.tri_pt3[i].tol[k];
                    let c = b.tri_pt3[i].coord[k].abs();
                    assert!(
                        t <= 4.0 * f64::EPSILON * c.max(1.0),
                        "reflection tol stays within a few ulp"
                    );
                }
            }
        }
        let ca = plane_classes(&crate::planes::test_judge(&plain));
        let cb = plane_classes(&crate::planes::test_judge(&mirrored));
        assert_eq!(
            ca, cb,
            "the two spellings agree on which faces are coplanar"
        );
        let pa = dense_planes(&plain, &ca).0;
        let pb = dense_planes(&mirrored, &cb).0;
        assert!(
            pb.iter()
                .all(|g| g.base.chain_id != 0 && g.base.tri.is_some() && g.base.coeffs.is_some()),
            "the reflected spelling really does offer a base frame — else this proves nothing"
        );
        assert!(
            pa.iter().all(|g| g.base.chain_id == 0),
            "the plain spelling really does decline"
        );
        assert!(
            pb.iter()
                .any(|g| g.base.tri.unwrap() != g.tri || g.base.coeffs.unwrap() != g.coeffs()),
            "the base frame differs from the moved one — else the comparison is free"
        );
        (pa, pb)
    }

    /// The three tests above compare a table that *offers* a base frame against one that does
    /// not — and agreement proves nothing unless the base frame is what answered. So corrupt it
    /// and require the answers to **change**: negating a base plane reverses both the side
    /// `indirect_plane_side` reports and the determinant `plane_pair_dir_sign` takes. If either
    /// sweep ever stops disagreeing, the shortcut has stopped firing and the tests above have
    /// gone vacuous.
    ///
    /// ★ **And the base *triangle* must now be inert.** `orient3d` used to ask the fourth plane
    /// with its pre-motion triangle while the other three spoke in coefficients, so a plane whose
    /// `d` was a rounded product was two planes inside one judgement. It asks all four in
    /// coefficients now, and reversing the triangles has to change nothing — that is the property,
    /// not an accident, so it is asserted rather than left untested.
    ///
    /// A caution the first draft of this test walked into: not every perturbation is detectable.
    /// Reflecting the base triangles *and* re-deriving their planes changes nothing at all — the
    /// two sign flips cancel, which is the same pseudovector identity this whole cell is about.
    #[test]
    fn the_base_frame_is_actually_consulted() {
        let (pa, pb) = two_spellings();
        let n = pa.len();

        let mut reversed = two_spellings().1;
        for g in &mut reversed {
            let [a, b, c] = g.base.tri.unwrap();
            g.base.tri = Some([b, a, c]);
        }
        let mut negated_planes = two_spellings().1;
        for g in &mut negated_planes {
            g.base.coeffs = g.base.coeffs.map(|c| c.map(|v| -v));
        }
        let mut orient_changed = 0usize;
        let mut tri_changed = 0usize;
        for p in 0..n {
            for q in (p + 1)..n {
                for r in (q + 1)..n {
                    if !normals_independent(&pa, p, q, r) {
                        continue;
                    }
                    for j in 0..n {
                        if j == p || j == q || j == r {
                            continue;
                        }
                        let want = crate::planes::test_judge(&pb).orient3d(p, q, r, j);
                        if want != crate::planes::test_judge(&negated_planes).orient3d(p, q, r, j) {
                            orient_changed += 1;
                        }
                        if want != crate::planes::test_judge(&reversed).orient3d(p, q, r, j) {
                            tri_changed += 1;
                        }
                    }
                }
            }
        }
        assert!(
            orient_changed > 0,
            "`orient3d` never read the base coefficients"
        );
        assert_eq!(
            tri_changed, 0,
            "`orient3d` still reads the base triangle — the fourth plane must be asked in \
             coefficients, or one judgement holds two descriptions of it again"
        );

        let mut negated = two_spellings().1;
        for g in &mut negated {
            g.base.coeffs = g.base.coeffs.map(|c| c.map(|v| -v));
        }
        let mut dir_changed = 0usize;
        for p in 0..n {
            for a in 0..n {
                for b in 0..n {
                    if p == a || p == b || a == b {
                        continue;
                    }
                    if crate::planes::test_judge(&pb).plane_pair_dir_sign(p, a, b)
                        != crate::planes::test_judge(&negated).plane_pair_dir_sign(p, a, b)
                    {
                        dir_changed += 1;
                    }
                }
            }
        }
        assert!(
            dir_changed > 0,
            "`plane_pair_dir_sign` never read the base plane"
        );
    }

    /// `orient3d` — the shortcut ([`nacre_cip::frame3::shared_base`] via `Judge::orient3d`) is
    /// answering the same question as the coordinates. See [`two_spellings`].
    #[test]
    fn a_reflected_spelling_orients_as_the_plain_one() {
        let (pa, pb) = two_spellings();
        let n = pa.len();
        let mut checked = 0usize;
        for p in 0..n {
            for q in (p + 1)..n {
                for r in (q + 1)..n {
                    if !normals_independent(&pa, p, q, r) {
                        continue;
                    }
                    for j in 0..n {
                        if j == p || j == q || j == r {
                            continue;
                        }
                        let sa = crate::planes::test_judge(&pa).orient3d(p, q, r, j);
                        if sa == 0 {
                            continue;
                        }
                        let sb = crate::planes::test_judge(&pb).orient3d(p, q, r, j);
                        assert_eq!(sa, sb, "orient3d at ({p},{q},{r},{j})");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 0, "no definite config in the corpus");
    }

    /// `plane_pair_dir_sign` — the one shortcut that reads **normals**, where a pseudovector sign
    /// would hide. See [`two_spellings`].
    #[test]
    fn a_reflected_spelling_takes_the_same_direction_signs() {
        let (pa, pb) = two_spellings();
        let n = pa.len();
        let mut checked = 0usize;
        for p in 0..n {
            for a in 0..n {
                for b in 0..n {
                    if p == a || p == b || a == b {
                        continue;
                    }
                    let sa = crate::planes::test_judge(&pa).plane_pair_dir_sign(p, a, b);
                    if sa == 0 {
                        continue;
                    }
                    let sb = crate::planes::test_judge(&pb).plane_pair_dir_sign(p, a, b);
                    assert_eq!(sa, sb, "dir_sign at ({p},{a},{b})");
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no definite triple in the corpus");
    }

    /// `planes_coplanar` and `cmp_coord`. The latter's shortcut is switched **off** for a
    /// reflection (a reflection negates the axis it compares along, which is not a handedness
    /// question), so what this pins is that the fallback still answers correctly.
    /// See [`two_spellings`].
    #[test]
    fn a_reflected_spelling_compares_coordinates_the_same_way() {
        let (pa, pb) = two_spellings();
        let n = pa.len();
        for i in 0..n {
            for j in 0..n {
                assert_eq!(
                    crate::planes::test_judge(&pa).planes_coplanar(i, j),
                    crate::planes::test_judge(&pb).planes_coplanar(i, j),
                    "coplanar at ({i},{j})"
                );
            }
        }
        let triples: Vec<[usize; 3]> = (0..n)
            .flat_map(|p| ((p + 1)..n).flat_map(move |q| ((q + 1)..n).map(move |r| [p, q, r])))
            .filter(|t| normals_independent(&pa, t[0], t[1], t[2]))
            .collect();
        let mut checked = 0usize;
        for &t in &triples {
            for &u in &triples {
                if t == u {
                    continue;
                }
                for axis in 0..3 {
                    let sa = crate::planes::test_judge(&pa).cmp_coord(t, u, axis);
                    if sa == 0 {
                        continue;
                    }
                    let sb = crate::planes::test_judge(&pb).cmp_coord(t, u, axis);
                    assert_eq!(sa, sb, "cmp_coord {t:?} vs {u:?} on axis {axis}");
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no definite comparison in the corpus");
    }

    // ---- stage 0b: does the two-caps defect reproduce, and is it wrong or conservative? ----

    /// Two bosses on one tilted face, raised to the same height by different arithmetic —
    /// `7.7` against `1.1` then `6.6`. Their caps are the same plane.
    ///
    /// Returns the model, the solid, and the two cap faces.
    /// How each column reaches `7.7` — the only thing that varies between a run and its control.
    ///
    /// ★ **A shape-matched control is not available.** Raising *both* columns as `1.1 + 6.6` would
    /// hold the shape fixed and change only the coincidence question, but it cannot be built: the
    /// two columns' shoulders are then one plane and the second column's upper pad is refused with
    /// `Rejected { SeamAlias }` — the same family of defect, one level down. `Single` differs by
    /// one shoulder at height `3.1`, which is nowhere near the caps at `9.7`, so the third control
    /// (the plate over one cap only) carries the isolation instead.
    #[derive(Clone, Copy, Debug, PartialEq)]
    enum Recipe {
        /// A: one `7.7`. B: `1.1` then `6.6`. The caps land one ulp apart.
        Split,
        /// Both columns: one `7.7`. Caps agree, but B has no shoulder.
        Single,
    }

    fn two_caps_on_a_tilted_face(
        recipe: Recipe,
    ) -> (Model, Handle<Solid>, [Handle<nacre_topo::Face>; 2]) {
        use crate::{OpOutput, Profile2d};
        use nacre_math::Vector3;
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let mut m = Model::new();
        let mut s = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        m.rebuild_adjacency();
        // Two turns, so the face's normal is off every world axis and its frame is not exact.
        for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid: s,
                    isometry: Isometry::rotation(SRot {
                        axis,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                    }),
                },
            )
            .expect("turn") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            s = solid;
        }
        let (sy, cy) = (53f64).to_radians().sin_cos();
        let (sz, cz) = (17f64).to_radians().sin_cos();
        let up = Vector3::from_array([cz * sy, sz * sy, cy]);
        // The **original** tilted face — lowest along `up` among the faces that point that way,
        // so the second boss is placed beside the first rather than on top of it.
        let facing = |m: &Model, s: Handle<Solid>| -> Handle<nacre_topo::Face> {
            *m.shells
                .get(m.solids.get(s).outer)
                .faces
                .iter()
                .filter(|&&f| {
                    crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99)
                })
                .min_by(|&&a, &&b| {
                    let h = |f: Handle<nacre_topo::Face>| {
                        (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                    };
                    h(a).partial_cmp(&h(b)).unwrap()
                })
                .expect("a face along up")
        };
        // A footprint in the face's own frame, centred on `(cu + lo … cu + hi, cv ± half)`.
        let footprint = |m: &Model, f: Handle<nacre_topo::Face>, lo: f64, hi: f64, half: f64| {
            let sp = crate::ops::face_plane(m, f).expect("planar");
            let d = nacre_props::face_props(m, f).unwrap().centroid - sp.origin;
            let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
            Profile2d::polygon(vec![
                p(cu + lo, cv - half),
                p(cu + hi, cv - half),
                p(cu + hi, cv + half),
                p(cu + lo, cv + half),
            ])
            .unwrap()
        };
        let pad = |m: &mut Model,
                   f: Handle<nacre_topo::Face>,
                   profile: Profile2d,
                   dist: f64|
         -> (Handle<Solid>, Handle<nacre_topo::Face>) {
            let OpOutput::PadOnFace { solid, top_face } = apply(
                m,
                &Operation::PadOnFace {
                    face: f,
                    profile,
                    dist,
                },
            )
            .unwrap_or_else(|e| panic!("pad dist={dist} on {f:?}: {e:?}")) else {
                unreachable!()
            };
            m.rebuild_adjacency();
            (solid, top_face)
        };
        // One column, raised either in one pad or in two.
        let column = |m: &mut Model,
                      s: &mut Handle<Solid>,
                      lo: f64,
                      hi: f64,
                      two_tier: bool|
         -> Handle<nacre_topo::Face> {
            let f = facing(m, *s);
            let prof = footprint(m, f, lo, hi, 0.5);
            if two_tier {
                let (s1, low) = pad(m, f, prof, 1.1);
                *s = s1;
                let prof = footprint(m, low, -0.25, 0.25, 0.25);
                let (s2, cap) = pad(m, low, prof, 6.6);
                *s = s2;
                cap
            } else {
                let (s1, cap) = pad(m, f, prof, 7.7);
                *s = s1;
                cap
            }
        };
        let (a_two_tier, b_two_tier) = match recipe {
            Recipe::Split => (false, true),
            Recipe::Single => (false, false),
        };
        let cap_a = column(&mut m, &mut s, -1.0, -0.4, a_two_tier);
        let cap_b = column(&mut m, &mut s, 0.4, 1.0, b_two_tier);
        (m, s, [cap_a, cap_b])
    }

    /// The world direction the tilted face points, for both fixtures.
    fn tilt_up() -> nacre_math::Vector3 {
        let (sy, cy) = (53f64).to_radians().sin_cos();
        let (sz, cz) = (17f64).to_radians().sin_cos();
        nacre_math::Vector3::from_array([cz * sy, sz * sy, cy])
    }

    /// The two highest faces along the tilt — the caps — as indices into `faces`.
    fn cap_pair(m: &Model, faces: &[FaceInfo]) -> (usize, usize) {
        let up = tilt_up();
        let mut along: Vec<(usize, f64)> = faces
            .iter()
            .enumerate()
            .filter_map(|(i, fi)| {
                let p = nacre_props::face_props(m, fi.face.expect("a face")).ok()?;
                (p.normal?.dot(up) > 0.99).then(|| (i, (p.centroid - Point3::origin()).dot(up)))
            })
            .collect();
        along.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        assert!(along.len() >= 2, "two caps point along the tilt");
        (along[0].0, along[1].0)
    }

    /// ★★★★★ **Two bosses that reach one height by different arithmetic land on one plane.**
    ///
    /// The heights are `2.0 + 7.7` against `2.0 + 1.1 + 6.6`. In the rationals those are the same
    /// number — `77/10` either way — but on a tilted face the cap is reached through an irrational
    /// frame, and the judge compares the plane's **witness**.
    ///
    /// ★★★ **This used to answer `Sign(Negative)` — two planes, with full confidence.** The
    /// witness was the exact ring points *realized*, and a third of those do not survive f64
    /// (`11/10` is not one), so the judge was describing the plane through three rounded points
    /// and answering correctly about the wrong plane. Now `Model::surface_points` carries the
    /// exact triple and the two caps meet.
    ///
    /// The `Single` control (both columns raised in one pad) merged before and merges now.
    #[test]
    fn two_caps_reached_by_different_arithmetic_are_one_plane() {
        for (recipe, want_same) in [(Recipe::Split, true), (Recipe::Single, true)] {
            let (m, s, _) = two_caps_on_a_tilted_face(recipe);
            let faces = collect_planes(&m, s).unwrap();
            let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces));
            let (ia, ib) = cap_pair(&m, &faces);
            let ha = nacre_props::face_props(&m, faces[ia].face.unwrap())
                .unwrap()
                .centroid;
            let hb = nacre_props::face_props(&m, faces[ib].face.unwrap())
                .unwrap()
                .centroid;
            let (ha, hb) = (
                (ha - Point3::origin()).dot(tilt_up()),
                (hb - Point3::origin()).dot(tilt_up()),
            );
            assert_eq!(
                canon[ia] == canon[ib],
                want_same,
                "{recipe:?}: caps at {ha:.17} and {hb:.17}"
            );
            assert_eq!(
                crate::planes::test_judge(&faces).planes_coplanar(ia, ib),
                want_same,
                "{recipe:?}: the judge and the class assignment must agree"
            );
            assert!(nacre_validate::validate(&m).is_empty());
        }
    }

    /// ★★★★★ **A plate laid across both caps builds a well-formed solid.**
    ///
    /// This is what the split cost, and it is why the split was not merely conservative. With the
    /// caps in different classes the result came back carrying a face whose stored surface normal
    /// was **76° off its own outward normal** — the invariant `collect_planes` asserts, so the
    /// solid could not be used as an operand again, and in release `orient_sign` would have been
    /// read off that dot. **Nothing reported it**: the volume was right and `validate` was clean.
    ///
    /// The cause was isolated with three controls, and the malformed face appeared only in the one
    /// cell where a plate spanned two caps in different classes. All four are kept: they are what
    /// says the fix closed the cause rather than the symptom.
    #[test]
    fn a_plate_across_both_caps_builds_a_well_formed_solid() {
        for (recipe, span, want_bad) in [
            (Recipe::Split, true, 0),
            (Recipe::Split, false, 0),
            (Recipe::Single, true, 0),
            (Recipe::Single, false, 0),
        ] {
            let (bad, volume_ok, valid) = plate_across_caps(recipe, span);
            assert_eq!(bad, want_bad, "{recipe:?} span={span}");
            // ★ And the two things a caller would look at say nothing is wrong.
            assert!(volume_ok, "{recipe:?} span={span}: volume");
            assert!(valid, "{recipe:?} span={span}: validate");
        }
    }

    /// Raise a plate from the taller cap — across both columns when `span`, over one otherwise —
    /// and report `(faces whose stored normal is not parallel to their own outward normal,
    /// volume is right, validate is clean)`.
    fn plate_across_caps(recipe: Recipe, span: bool) -> (usize, bool, bool) {
        use crate::{OpOutput, Profile2d};
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let (mut m, s, _) = two_caps_on_a_tilted_face(recipe);
        let up = tilt_up();
        let before = nacre_props::mass_props(&m, s).unwrap().volume;
        let cap = *m
            .shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .filter(|&&f| crate::ops::face_plane(&m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99))
            .max_by(|&&a, &&b| {
                let h = |f: Handle<nacre_topo::Face>| {
                    (nacre_props::face_props(&m, f).unwrap().centroid - Point3::origin()).dot(up)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a cap");
        let sp = crate::ops::face_plane(&m, cap).expect("planar");
        let d = nacre_props::face_props(&m, cap).unwrap().centroid - sp.origin;
        let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
        // The columns sit at u ≈ ∓0.7 of the base face's centre; `-1.6` reaches the far one.
        let lo = if span { -1.6 } else { -0.2 };
        let profile = Profile2d::polygon(vec![
            p(cu + lo, cv - 0.2),
            p(cu + 0.4, cv - 0.2),
            p(cu + 0.4, cv + 0.2),
            p(cu + lo, cv + 0.2),
        ])
        .unwrap();
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: cap,
                profile,
                dist: 0.5,
            },
        )
        .unwrap_or_else(|e| panic!("{recipe:?} span={span}: the plate must build, got {e:?}")) else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let after = nacre_props::mass_props(&m, solid).unwrap().volume;
        let volume_ok = (after - before - (0.4 - lo) * 0.4 * 0.5).abs() < 1e-9;
        let valid = nacre_validate::validate(&m).is_empty();
        // `collect_planes` would assert on this; count it instead, so the test can say how many.
        let bad = m
            .shells
            .get(m.solids.get(solid).outer)
            .faces
            .iter()
            .filter(|&&f| {
                let face = m.faces.get(f);
                let Some((tri, _)) = crate::planes::outer_tri(&m, face) else {
                    return false;
                };
                let Some(n_out) = (tri[1] - tri[0]).cross(tri[2] - tri[0]).normalize() else {
                    return false;
                };
                let nacre_geom::Surface::Plane(pl) = m.surface(face.surface) else {
                    return false;
                };
                pl.normal().dot(n_out).abs() <= 0.5
            })
            .count();
        (bad, volume_ok, valid)
    }

    // ---- stage 0g: the prediction the whole plan rests on ----

    /// ★★★★★ **Are the two caps *really* one plane?**
    ///
    /// The plan's premise is that describing each cap by **exact points** — rather than by the
    /// rounded f64 witness — makes the judge see one plane. That is an algebraic claim about the
    /// composed frames, and nothing has measured it.
    ///
    /// This measures it without changing any production type. Each cap already carries what the
    /// plan would store: exact **frame** coefficients (`Model::surface_name`) and the motion that
    /// carries them out (`SurfaceDef::Moved`). Three exact points on that frame plane, replayed
    /// through that motion, are the definition the plan proposes — for the *plane* question any
    /// non-collinear triple on the plane is equivalent, so the synthetic triple answers it.
    ///
    /// The judge is then run with a **deliberately unreachable** coincidence limit, so it never
    /// short-circuits and instead reports the bound it achieved. Reading that across precisions:
    ///
    /// * falling like `2^-prec` ⇒ the true value is **0** ⇒ the caps are one plane ⇒ premise holds;
    /// * flattening at some nonzero value ⇒ they are genuinely different planes and the plan's
    ///   diagnosis is wrong.
    ///
    /// ★★★★★ **Measured: it falls one bit per bit**, and the mantissa stays at `0.9338…` — the
    /// bound is `C · 2^-prec` with a constant `C`, exactly the model `escalate` derives its jump
    /// from. So the caps *are* one plane, and what splits them today is only the rounded witness.
    ///
    /// ```text
    /// prec=128  within=2^-114
    /// prec=192  within=2^-178   (-64)
    /// prec=256  within=2^-242   (-64)
    /// prec=384  within=2^-370   (-128)
    /// prec=512  within=2^-498   (-128)
    /// ```
    #[test]
    fn two_caps_described_exactly_are_one_plane() {
        use nacre_cip::predicate::{Judge, Notes};
        use nacre_cip::{Standard, WitnessPoint};
        use nacre_scalar::Mag;

        let (m, s, _) = two_caps_on_a_tilted_face(Recipe::Split);
        let faces = collect_planes(&m, s).unwrap();
        let (ia, ib) = cap_pair(&m, &faces);

        // The exact definition each cap already carries: frame coefficients + the motion.
        // `nudge` offsets the cap's height inside its own frame — the negative control.
        let define = |i: usize, nudge: Option<Rat>| -> [WitnessPoint; 3] {
            let surf = faces[i].surf;
            let c = *m
                .surface_name
                .get(&surf)
                .expect("a frame cap records its name")
                .narrow()
                .expect("a frame cap's name is narrow");
            let motion = match m.surface_truth(surf) {
                nacre_topo::SurfaceTruth::Plane {
                    motion: Some(motion),
                    ..
                } => *motion,
                other => panic!("cap {i} is not moved: {other:?}"),
            };
            // A cap's frame coefficients are `[0, 0, ±1, ∓h]`, so `z = -d/c` and any two in-frame
            // directions complete the triple.
            assert!(
                c[0] == Rat::from_int(0) && c[1] == Rat::from_int(0),
                "cap {i} is not a frame cap: {c:?}"
            );
            // Canonicalisation divides the content, so `c[2]` is any nonzero integer — the height
            // is `−d / c` exactly, as a rational.
            let (dn, dd) = (c[3].numer(), c[3].denom());
            let (cn, cd) = (c[2].numer(), c[2].denom());
            let z = Rat::new(
                dn.checked_mul(cd)
                    .expect("no overflow")
                    .checked_neg()
                    .unwrap(),
                dd.checked_mul(cn).expect("no overflow"),
            )
            .expect("a cap's height");
            let z = match nudge {
                Some(n) => z.checked_add(n).expect("nudge"),
                None => z,
            };
            let zero = Rat::from_int(0);
            let one = Rat::from_int(1);
            let chain = crate::rotated_vertex::motion_chain(&m, motion).expect("chain");
            [[zero, zero, z], [one, zero, z], [zero, one, z]].map(|b| {
                crate::rotated_vertex::replay(WitnessPoint::at(b), &chain).expect("replay")
            })
        };
        let (da, db) = (define(ia, None), define(ib, None));

        // A two-plane table over those definitions. `base_rat: None` keeps the composed-rotation
        // shortcut out of it, so what runs is the interval route the plan's stage C1 exercises.
        let mk = |d: [WitnessPoint; 3]| {
            let tri = d.clone().map(|p| Point3::from_array(p.coord));
            WorkingPlane {
                base_rat: None,
                name_ints: None,
                base: crate::planes::BaseFrame::none(),
                surf: faces[ia].surf,
                plane: nacre_geom::Plane::through_points(tri[0], tri[1], tri[2])
                    .expect("non-collinear"),
                tri,
                tri_pt3: d,
                rotated: true,
                frame_sign: 1,
                exact_coeffs: None,
                exact_normal: None,
            }
        };
        let planes = vec![mk(da), mk(db)];

        let mut prev: Option<(usize, i64)> = None;
        for prec in [128usize, 192, 256, 384, 512] {
            let notes = Notes::new();
            // ★ A limit no realization can reach, so the judge never answers `Coincident` and has
            // to report the bound it actually achieved.
            let standard = Standard {
                prec,
                coincidence: Mag::pow2(-4000),
                scale: Mag::of(16.0),
                cap: prec,
            };
            let same = Judge::new(&planes, standard, &notes).planes_coplanar(0, 1);
            let ev = notes.sorted();
            let within = ev.first().and_then(|e| match e.outcome {
                nacre_cip::Decision::Coincident { within } => Some(within),
                nacre_cip::Decision::Exhausted { within, .. } => within,
                _ => None,
            });
            assert!(same, "prec={prec}: no definite separation may be found");
            // `Mag` is `m · 2^e`; the exponent is the reading that matters — a true zero makes it
            // fall one per bit of precision, a real separation makes it flatten.
            let w = within
                .and_then(|b| b.exp2())
                .unwrap_or_else(|| panic!("prec={prec}: no bound to read, outcome {ev:?}"));
            if let Some((pp, pw)) = prev {
                let gained = pw - w;
                let spent = (prec - pp) as i64;
                assert!(
                    gained >= spent - 4,
                    "prec {pp} -> {prec} spent {spent} bits and gained only {gained}: \
                     the residual is flattening, so the two caps are not one plane"
                );
            }
            prev = Some((prec, w));
        }

        // ★★★★ **The negative control — without it the loop above proves nothing.**
        //
        // Offset one cap by `1e-12` *inside its own frame* and the same sweep must stop gaining a
        // bit per bit: a real separation is a floor the precision cannot go under. If this half
        // passed too, the assertion above would be measuring the ladder rather than the geometry.
        let off = Rat::new(1, 1_000_000_000_000).expect("1e-12");
        let planes = vec![mk(define(ia, None)), mk(define(ib, Some(off)))];
        let mut prev: Option<(usize, i64)> = None;
        let mut flattened = false;
        for prec in [128usize, 192, 256, 384, 512] {
            let notes = Notes::new();
            let standard = Standard {
                prec,
                coincidence: Mag::pow2(-4000),
                scale: Mag::of(16.0),
                cap: prec,
            };
            let same = Judge::new(&planes, standard, &notes).planes_coplanar(0, 1);
            if !same {
                flattened = true; // a definite separation — even better than a floor
                break;
            }
            let ev = notes.sorted();
            let w = match ev.first().and_then(|e| match e.outcome {
                nacre_cip::Decision::Coincident { within } => Some(within),
                nacre_cip::Decision::Exhausted { within, .. } => within,
                _ => None,
            }) {
                Some(b) => b.exp2().expect("a bound"),
                None => break,
            };
            if let Some((pp, pw)) = prev
                && pw - w < (prec - pp) as i64 - 4
            {
                flattened = true;
                break;
            }
            prev = Some((prec, w));
        }
        assert!(
            flattened,
            "a cap offset by 1e-12 still looked like the same plane at every precision — \
             the sweep is not measuring the geometry"
        );
    }

    /// ★ **How many planes now carry an exact point triple, and do the triples describe the
    /// plane they are filed under?**
    ///
    /// The points are written but nothing reads them yet, so this is what says the producers were
    /// wired — and the residual check is what says they were wired *correctly*: a triple filed
    /// under the wrong surface would still be exact and still be three points.
    #[test]
    fn every_plane_that_can_records_its_three_exact_points() {
        for (name, recipe) in [("split", Recipe::Split), ("single", Recipe::Single)] {
            let (m, s, _) = two_caps_on_a_tilted_face(recipe);
            let faces = collect_planes(&m, s).unwrap();
            let (mut with, mut without) = (0usize, 0usize);
            for fi in &faces {
                let nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Known(pts),
                    ..
                } = m.surface_truth(fi.surf)
                else {
                    without += 1;
                    continue;
                };
                with += 1;
                // ★ The points must satisfy the coefficients they are filed beside — both are
                // stated in the same frame, so this is an exact rational identity with no
                // tolerance anywhere.
                if let Some(c) = m.surface_name.get(&fi.surf).and_then(|n| n.narrow()) {
                    for p in pts {
                        let mut acc = c[3];
                        for k in 0..3 {
                            acc = acc
                                .checked_add(c[k].checked_mul(p[k]).expect("no overflow"))
                                .expect("no overflow");
                        }
                        assert_eq!(
                            acc,
                            Rat::from_int(0),
                            "{name}: a recorded point is not on its own plane {:?}",
                            fi.surf
                        );
                    }
                }
            }
            println!("{name}: {with} planes with points, {without} without");
            assert!(
                with > 0 && without == 0,
                "{name}: {without} of {} faces record no points",
                with + without
            );
        }
    }

    /// ★★★★★ **A boolean creates no surface — so a plane never comes from a discovered point.**
    ///
    /// This is the well-foundedness the whole "a plane's truth is three construction points"
    /// design rests on. If a boolean minted a plane, that plane's points would be seam vertices,
    /// which are themselves defined *by planes* — and the two-level structure would become an
    /// unbounded recursion with no base case.
    ///
    /// It holds because the arrangement reuses its operands' surface handles: the classes it
    /// computes are over the planes it was handed, and a result face carries the handle of
    /// whichever operand face it came from. Only `ops.rs` (construction) and `transform.rs`
    /// (moving one) ever push.
    ///
    /// ★ Asserted rather than argued: the grep that says so today is not a thing the compiler
    /// keeps true, and the failure it would hide is silent.
    #[test]
    fn a_boolean_mints_no_surface() {
        use crate::BoolKind;
        for (name, recipe) in [("split", Recipe::Split), ("single", Recipe::Single)] {
            let (mut m, s, _) = two_caps_on_a_tilted_face(recipe);
            // A second body to cut with, built before the count so its own construction does not
            // land inside the window.
            // The turns pivot on the origin, so the base cuboid's corner there does not move —
            // a box straddling it bites the solid whatever the tilt.
            let tool = m.add_cuboid(
                Point3::from_array([-0.5, -0.5, -0.5]),
                Point3::from_array([0.6, 0.6, 0.6]),
            );
            m.rebuild_adjacency();
            let vol_before = nacre_props::mass_props(&m, s).unwrap().volume;
            let before = m.surface_count();
            let out = crate::boolean(&mut m, BoolKind::Cut, s, tool);
            let after = m.surface_count();
            let bodies = out.unwrap_or_else(|e| panic!("{name}: the cut must solve, got {e:?}"));
            assert_eq!(
                after,
                before,
                "{name}: the boolean minted {} surface(s) — a plane defined by discovered points \
                 would break the design's base case",
                after - before
            );
            // ★ And the cut actually cut: a boolean that missed would keep the arena still for
            // the wrong reason.
            m.rebuild_adjacency();
            let vol_after: f64 = bodies
                .iter()
                .map(|&b| nacre_props::mass_props(&m, b).unwrap().volume)
                .sum();
            assert!(
                vol_after < vol_before - 1e-9,
                "{name}: the tool missed — {vol_before} -> {vol_after}"
            );
        }
    }

    /// ★★★★ **A named plane's triple names the caller's plane, full-width normals included —
    /// and the point order carries the caller's polarity.**
    ///
    /// History: this used to pin *"`c · p` is computable in `i128`"* and chose an axis-solved
    /// triple to keep every point coefficient-sized (a frame-walk triple's third point is
    /// `n × u`, a product, and 1,551 of 83,813 names went out unverified over it). That
    /// computability was a design constraint of the narrow-only pipeline; since S2 (wide names)
    /// and S4 (wide frames) an overflowing narrow derivation takes the arbitrary-precision road
    /// instead of failing, so `PlaneDef` states the *structural* triple `[o, o+u, o+v]` and the
    /// old property is retired. What must still hold, and is pinned here on the very normals
    /// whose widths used to break things:
    ///
    /// - the definition **exists** (no width-based decline is left in `normal_def`),
    /// - the triple **names a plane** (`plane_name_exact` — total, `None ⇔ collinear`),
    /// - the derived name is the same one the retired coefficient route computed
    ///   (`plane_from_point_normal` on the lifted normal — same plane, one name),
    /// - the point order faces the **caller's** normal (`u × v = |u|²·n` — checked through the
    ///   realization, far from degenerate, so the f64 sign is exact here).
    #[test]
    fn a_named_planes_points_name_the_callers_plane() {
        use crate::SketchPlane;
        use nacre_math::Vector3;
        for n in [
            [0.3, -0.2, 1.0],
            [0.3141592653589793, -0.2718281828459045, 1.0],
            [1.0, 0.0, 0.0],
            [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
            [0.1, 0.2, 0.30000000000000004],
        ] {
            let o = Point3::from_array([0.25, -0.5, 1.5]);
            let p = SketchPlane::from_origin_normal(o, Vector3::from_array(n)).expect("a plane");
            let d = p.def.as_ref().expect("an exact definition");
            let pts = d.points();
            let name = nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2])
                .unwrap_or_else(|| panic!("{n:?}: the triple is collinear"));

            // Same plane as the coefficient route would have stated (both from the same lifts).
            let lift = |v: [f64; 3]| v.map(|x| Rat::from_decimal(x).unwrap());
            let coeffs = nacre_scalar::plane_from_point_normal(lift(n), lift(o.as_array()))
                .unwrap_or_else(|| panic!("{n:?}: the reference route overflowed"));
            assert_eq!(
                name,
                nacre_scalar::PlaneName::Narrow(coeffs),
                "{n:?}: the points name a different plane than the caller's normal"
            );

            // Polarity: (p1−p0) × (p2−p0) must face the caller's n, not its negation.
            let f = |q: [Rat; 3]| q.map(|r| r.to_f64());
            let (a, b, c) = (f(pts[0]), f(pts[1]), f(pts[2]));
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let dot: f64 = (0..3).map(|k| cross[k] * n[k]).sum();
            assert!(
                dot > 0.0,
                "{n:?}: the point order faces away from the caller"
            );
        }
    }

    /// ★★★★★ **Did the operation take the exact road, or only arrive at the right answer?**
    ///
    /// Everything else in this suite checks the *answer*: a volume, a face count, a coordinate.
    /// None of them can see a build that silently gave up on exact arithmetic and reached the
    /// same number through f64 — and that is precisely the failure that has to be caught here,
    /// because it does not misbehave until two such builds have to agree with each other.
    ///
    /// ★ **A prism on an arbitrarily tilted plane is the case that exercises it.** An axis-aligned
    /// sketch lifts to exact rationals on its own (`SketchPlane::exact`), so it takes the exact
    /// road even when everything else is broken. A tilted one can only get there through its own
    /// frame — the caller's `PlaneDef`, its coefficients, `frame_chain`, `Motion::Frame` — and
    /// every link in that chain is load-bearing. Break any one and the extrude quietly reverts to
    /// the f64 path, still producing a prism of the right volume.
    ///
    /// ★★★★★ Written because that is what happened: teaching `push_surface_with_coeffs` to drop a
    /// name it could not verify dropped it for *every* tilted plane (the check overflows on wide
    /// coefficients), which unhooked the frame and closed this road — and the whole suite, plus a
    /// bit-identical coordinate census, said nothing at all.
    #[test]
    fn a_prism_on_a_tilted_plane_takes_the_exact_road() {
        use crate::{OpOutput, Profile2d, SketchPlane};
        use nacre_math::{Point2, Vector3};

        // ★★ **A full-width direction, not a tidy one.** `[0.3, -0.2, 1.0]` is also tilted, but its
        // decimals are three digits long, so every rational derived from it stays small and the
        // road stays open however wide the arithmetic gets. The normals a real caller hands over —
        // a slider, a computed draft angle, a proptest — carry seventeen digits, and those are what
        // push the coefficients out to a hundred-odd bits. Both belong here.
        for (name, n) in [
            ("short decimals", [0.3, -0.2, 1.0]),
            ("full-width", [0.3141592653589793, -0.2718281828459045, 1.0]),
        ] {
            let plane = SketchPlane::from_origin_normal(
                Point3::from_array([0.25, -0.5, 1.5]),
                Vector3::from_array(n),
            )
            .expect("a tilted plane");
            let sq = Profile2d::polygon(vec![
                Point2::from_array([-1.0, -1.0]),
                Point2::from_array([1.0, -1.0]),
                Point2::from_array([1.0, 1.0]),
                Point2::from_array([-1.0, 1.0]),
            ])
            .unwrap();
            let mut m = Model::new();
            let __f104 = datum_frame(&mut m, plane);
            let out = apply(
                &mut m,
                &Operation::Extrude {
                    frame: __f104,
                    profile: sq,
                    dist: 0.75,
                },
            )
            .unwrap_or_else(|e| panic!("{name}: the extrude must build, got {e:?}"));
            let OpOutput::Extrude { solid: s, .. } = out else {
                panic!("{name}: extrude returned no solid")
            };
            m.rebuild_adjacency();

            let faces = collect_planes(&m, s).unwrap();
            let (mut exact, mut f64_only) = (0usize, 0usize);
            for fi in &faces {
                // ★ **The points are the whole test.** A plane with an exact triple was reached by
                // exact arithmetic; one without was computed in f64 and rounded. The `def` is a
                // separate axis and both of its exact spellings are right here: the base cap is
                // `Constructed` because it *is* the world plane the caller named exactly, and
                // everything else is `Moved` against the sketch frame. Only `Inexact` — a plane
                // the forest cannot replay — means the road was lost.
                // S6b: the truth is total — a planar surface always carries points, so the
                // old "has points and is not Inexact" test collapses to "is a plane truth",
                // which the type now guarantees. The sweep stays as the retrospective record
                // of what this lock used to have to check.
                let ok = matches!(
                    m.surface_truth(fi.surf),
                    nacre_topo::SurfaceTruth::Plane { .. }
                );
                if ok {
                    exact += 1;
                } else {
                    f64_only += 1;
                    println!(
                        "  {name}: {:?} named={} truth={:?}",
                        fi.surf,
                        m.surface_name.contains_key(&fi.surf),
                        m.surface_truth(fi.surf)
                    );
                }
            }
            println!("{name}: {exact} faces exact, {f64_only} on the f64 path");
            assert_eq!(
                f64_only,
                0,
                "{name}: {f64_only} of {} faces of a tilted prism fell to the f64 path — the \
                 exact road is closed and no answer-checking test can see it",
                exact + f64_only
            );
        }
    }

    /// ★★★★★ **One plane, one handle — the gain, stated as something that can fail.**
    ///
    /// Everything else in this suite says *"nothing broke"*. This says what the wide derivation
    /// **bought**: a U-shaped footprint has two walls on one plane, and until the kernel could name
    /// that plane they were two `Surface` handles for one thing. Handle identity is what
    /// `plane_classes` merges on in O(1) and what a coplanar contact is recognised by, so a plane it
    /// cannot name is a plane the engine has to rediscover geometrically every time.
    ///
    /// ★ **The coordinates have to be full-width.** On short decimals (`1.0`, `0.5`) the narrow
    /// `Rat` route already names every wall and the two already share a handle — a version of this
    /// test written that way passes before the change and proves nothing. These are seventeen
    /// digits, which is what a slider, a computed dimension or a proptest actually produces.
    #[test]
    fn two_walls_of_one_plane_become_one_surface() {
        use crate::{OpOutput, Profile2d, SketchFrame};
        use nacre_math::Point2;
        let t = 1.9384756293847; // the plane the notch's two arms share
        let p = |x: f64, y: f64| Point2::from_array([x, y]);
        let mut m = Model::new();
        let __w0 = SketchFrame::world(&m, Axis::Z);
        let out = apply(
            &mut m,
            &Operation::Extrude {
                frame: __w0,
                profile: Profile2d::polygon(vec![
                    p(0.0, 0.0),
                    p(2.8374652839472, 0.0),
                    p(2.8374652839472, t),
                    p(1.8473625849372, t),
                    p(1.8473625849372, 0.9384756293847),
                    p(0.9937465283947, 0.9384756293847),
                    p(0.9937465283947, t),
                    p(0.0, t),
                ])
                .unwrap(),
                dist: 0.8473625849372,
            },
        )
        .expect("the extrude must build");
        let OpOutput::Extrude { faces, .. } = out else {
            panic!("extrude returned no faces")
        };
        let mut surfaces: Vec<_> = faces.iter().map(|&f| m.faces.get(f).surface).collect();
        let named = faces
            .iter()
            .filter(|&&f| m.surface_name.contains_key(&m.faces.get(f).surface))
            .count();
        let total = surfaces.len();
        surfaces.sort_unstable();
        surfaces.dedup();
        assert_eq!(
            named,
            total,
            "{} of {total} faces sit on an unnamed plane",
            total - named
        );
        assert_eq!(
            surfaces.len(),
            total - 1,
            "the notch's two arms are one plane and must share one surface handle \
             ({total} faces, {} distinct surfaces)",
            surfaces.len()
        );
    }

    /// ★★ Open item 15 — a **wide** name's integer rescue, locked two ways: against the climb
    /// on the real population (shared-motion gate, mirrored included), and by a starved judge
    /// that can answer *only* through the name (world gate, path proved from both sides).
    mod wide_name_rescue {
        use super::*;
        use crate::planes::WorkingPlane;
        use crate::{DatumDef, OpOutput, Profile2d};
        use nacre_cip::Standard;
        use nacre_cip::predicate::{Judge, Notes, name_stored_ints};
        use nacre_math::Point2;
        use nacre_scalar::{Mag, PlaneName};
        use num_bigint::BigInt;

        fn p2(x: f64, y: f64) -> Point2 {
            Point2::from_array([x, y])
        }

        /// The census `wf` family with a pocket — the kernel's one natural producer of a wide
        /// canonical name (three of its discovered vertices name a plane past `i128`;
        /// `tests/wide_datum_cost.rs` measured it, and is duplicated here because an
        /// integration test cannot be imported). A prism on a tilted decimal orthonormal
        /// frame, with a pocket on the frame-general wall.
        fn wf_pocket_model() -> Model {
            let mut m = Model::new();
            let plane = crate::SketchPlane::from_axes(
                Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
                nacre_math::Vector3::from_array([0.6, 0.8, 0.0]),
                nacre_math::Vector3::from_array([-0.48, 0.36, 0.8]),
            );
            let frame = datum_frame(&mut m, plane);
            let OpOutput::Extrude { solid, .. } = apply(
                &mut m,
                &Operation::Extrude {
                    frame,
                    profile: Profile2d::polygon(vec![
                        p2(0.1111111111111111, 0.1234567890123456),
                        p2(4.123456789012345, 0.2345678901234567),
                        p2(3.9876543210987654, 3.1234567890123459),
                        p2(0.2222222222222222, 2.765432109876543),
                    ])
                    .unwrap(),
                    dist: 2.5,
                },
            )
            .expect("the wf base prism") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            let wall = *m
                .shells
                .get(m.solids.get(solid).outer)
                .faces
                .iter()
                .find(|&&f| {
                    let s = m.faces.get(f).surface;
                    m.surface_name
                        .get(&s)
                        .and_then(|n| n.narrow())
                        .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
                })
                .expect("the wf population");
            let sp = crate::face_plane(&m, wall).expect("planar");
            // An interior point of the convex wall: its outer-loop vertex average (the
            // integration test uses the area centroid; any interior point serves).
            let pts: Vec<Point3> = m
                .faces
                .get(wall)
                .outer
                .half_edges
                .iter()
                .map(|&he| m.vertex_point(m.he_start(he)))
                .collect();
            let n = pts.len() as f64;
            let c = Point3::from_array(std::array::from_fn(|k| {
                pts.iter().map(|p| p.as_array()[k]).sum::<f64>() / n
            }));
            let d = c - sp.origin();
            let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
            apply(
                &mut m,
                &Operation::PocketOnFace {
                    face: wall,
                    profile: Profile2d::polygon(vec![
                        p2(cu - 0.3, cv - 0.3),
                        p2(cu + 0.3, cv - 0.3),
                        p2(cu + 0.3, cv + 0.3),
                        p2(cu - 0.3, cv + 0.3),
                    ])
                    .unwrap(),
                    dist: 0.4,
                },
            )
            .expect("the wf pocket");
            m.rebuild_adjacency();
            m
        }

        /// Every live vertex, deduplicated in first-seen order.
        fn live_verts(m: &Model) -> Vec<nacre_store::Handle<nacre_topo::Vertex>> {
            let mut seen = std::collections::HashSet::new();
            let mut out = Vec::new();
            for &s in &m.live_solids {
                let sol = m.solids.get(s);
                for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
                    for &fh in &m.shells.get(sh).faces {
                        let f = m.faces.get(fh);
                        for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                            for &he in &lp.half_edges {
                                let vh = m.he_start(he);
                                if seen.insert(vh) {
                                    out.push(vh);
                                }
                            }
                        }
                    }
                }
            }
            out
        }

        /// A tool raised on a **wide** vertex-named datum of `m` — the first triple whose
        /// canonical name needs the wide vessel and whose datum carries an extrude
        /// (`half`-sized square profile, `dist` deep). The carriers are the pocket's
        /// *discovered* vertices, which carry no motion chain, so the datum names the
        /// **world** — the probe below (`the wide plane names the world`) is what taught this
        /// module that fact, and it is why the populations under test are the world gate and a
        /// later whole-solid motion, not the boolean that first makes the wide plane.
        fn wide_datum_tool(m: &mut Model, half: f64, dist: f64) -> Option<Handle<Solid>> {
            let verts = live_verts(m);
            for i in 0..verts.len() {
                for j in (i + 1)..verts.len() {
                    for k in (j + 1)..verts.len() {
                        let mut t = [verts[i], verts[j], verts[k]];
                        t.sort_by_key(|v| v.index());
                        match m.plane_name_through(t) {
                            Some(n) if n.narrow().is_none() => {}
                            _ => continue,
                        }
                        let Ok(OpOutput::DatumPlane { frame, .. }) = apply(
                            m,
                            &Operation::DatumPlane {
                                def: DatumDef::ThroughVertices([verts[i], verts[j], verts[k]]),
                            },
                        ) else {
                            continue;
                        };
                        if let Ok(OpOutput::Extrude { solid, .. }) = apply(
                            m,
                            &Operation::Extrude {
                                frame,
                                profile: Profile2d::polygon(vec![
                                    p2(-half, -half),
                                    p2(half, -half),
                                    p2(half, half),
                                    p2(-half, half),
                                ])
                                .unwrap(),
                                dist,
                            },
                        ) {
                            m.rebuild_adjacency();
                            return Some(solid);
                        }
                    }
                }
            }
            None
        }

        /// A solid whose faces are **all world-named, one of them wide**: a cuboid sliced by
        /// the wide plane alone — the tool is made so large that its walls and far cap miss
        /// the cuboid entirely, so the result's faces are the trimmed cuboid's (integer world
        /// names) plus the wide cap (the discovered-vertex world name).
        fn cuboid_sliced_by_the_wide_plane(m: &mut Model) -> Handle<Solid> {
            let tool =
                wide_datum_tool(m, 50.0, 50.0).expect("a wide-named datum carries a slab tool");
            let cub = m.add_cuboid(
                Point3::from_array([-1.0, -1.0, -1.0]),
                Point3::from_array([5.0, 5.0, 4.0]),
            );
            let out = crate::boolean(m, crate::BoolKind::Cut, cub, tool)
                .expect("the wide slab slices the cuboid");
            m.rebuild_adjacency();
            *out.first().expect("the slice leaves material")
        }

        /// Whether the name-integer gate opens for these planes — the fixture-validity probe
        /// (mirrors `Judge::name_rescue`'s condition; a fixture where no question passes it
        /// would run the differential against nothing).
        fn gate_opens(t: &[WorkingPlane], idx: &[usize]) -> bool {
            idx.iter().all(|&k| t[k].name_ints.is_some())
                && idx
                    .iter()
                    .any(|&k| t[k].name_ints.as_ref().is_some_and(|n| n.wide))
                && (!idx.iter().any(|&k| t[k].rotated)
                    || (t[idx[0]].base.chain_id != 0
                        && idx
                            .iter()
                            .all(|&k| t[k].base.chain_id == t[idx[0]].base.chain_id)))
        }

        /// ★★★ **The probe that corrected the plan's population claim.** The wide datum's
        /// carriers are the pocket's *discovered* vertices — solved world rationals with no
        /// motion chain — so the wide name is a **world** description (`chain = []`,
        /// `rotated = false`), not a moved one. The boolean that first makes the wide plane
        /// therefore asks mixed-frame questions (world-wide cap × frame-moved prism walls ×
        /// `FrameWide`-raised tool walls) that *neither* gate can speak to; the name's real
        /// populations are the world gate and a later solid moved **whole**. Kept as its own
        /// test because the wrong claim survived investigation and review twice.
        #[test]
        fn the_wide_plane_names_the_world() {
            let mut m = wf_pocket_model();
            let tool = wide_datum_tool(&mut m, 0.4, 1.2).expect("a wide datum carries a tool");
            let table = plane_table(&m, tool);
            let wide = table
                .iter()
                .find(|p| p.name_ints.as_ref().is_some_and(|n| n.wide))
                .expect("the tool's base cap carries the wide name");
            assert!(
                wide.tri_pt3[0].chain.is_empty() && !wide.rotated,
                "the wide name speaks for the world, not for a frame"
            );
            // And its table-mates are moved (the tool's walls live on the derived frame), so
            // no question pairing them with the wide cap passes either gate — the fact that
            // sends the first boolean's judgements to the climb, by design.
            assert!(
                table.iter().any(|p| p.rotated),
                "the first boolean's table is mixed-frame"
            );
        }

        /// The rescue against the climb, over the populations the gates cover:
        ///
        /// - **world** — the wide cap beside an unmoved cuboid, one combined table (the shape
        ///   a boolean of the two would judge): every cap-and-cuboid question is all-world,
        ///   all-named, one wide.
        /// - **rotated** — a solid whose faces are all world-named with the wide slice among
        ///   them, then turned 30° as one: every question shares one chain (even parity).
        /// - **mirrored** — the same through a reflection: odd parity, where the pseudovector
        ///   correction (negate `y`, `z`, `d` — keep `x`) is the part a sign slip would hide
        ///   in.
        ///
        /// Both routes are exact, so any disagreement is a bug in one of them. That the gate
        /// opens at all is asserted per arm (`gate_opens` — a differential over a gate that
        /// never fires would measure nothing), and that only the name route can answer is
        /// proved by `a_starved_judge_answers_only_through_the_name`.
        #[test]
        fn a_wide_name_judges_the_same_signs_the_climb_does() {
            for arm in ["world", "rotated", "mirrored"] {
                let mut m = wf_pocket_model();
                let named: Vec<WorkingPlane> = if arm == "world" {
                    let tool =
                        wide_datum_tool(&mut m, 0.4, 1.2).expect("a wide datum carries a tool");
                    let cub = m.add_cuboid(
                        Point3::from_array([-1.0, -1.0, -1.0]),
                        Point3::from_array([5.0, 5.0, 4.0]),
                    );
                    let faces: Vec<crate::planes::FaceInfo> = collect_planes(&m, cub)
                        .unwrap()
                        .into_iter()
                        .chain(collect_planes(&m, tool).unwrap())
                        .collect();
                    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces));
                    crate::planes::dense_planes(&faces, &canon).0
                } else {
                    let sliced = cuboid_sliced_by_the_wide_plane(&mut m);
                    let moved = if arm == "rotated" {
                        rotated(&mut m, sliced)
                    } else {
                        let Ok(OpOutput::Mirror { solid }) = apply(
                            &mut m,
                            &Operation::Mirror {
                                solid: sliced,
                                axis: Axis::X,
                                offset: Rat::from_int(0),
                            },
                        ) else {
                            panic!("mirror applies");
                        };
                        m.rebuild_adjacency();
                        solid
                    };
                    plane_table(&m, moved)
                };
                let n = named.len();
                let mut stripped = named.clone();
                for p in &mut stripped {
                    p.name_ints = None;
                }
                let jn = crate::planes::test_judge(&named);
                let js = crate::planes::test_judge(&stripped);
                let (mut open3, mut open4) = (0usize, 0usize);
                let mut diffs: Vec<String> = Vec::new();
                for p in 0..n {
                    for q in (p + 1)..n {
                        for r in (q + 1)..n {
                            // The predicates' precondition: the triple meets in a point.
                            if !normals_independent(&named, p, q, r) {
                                continue;
                            }
                            for j in 0..n {
                                if j == p || j == q || j == r {
                                    continue;
                                }
                                open4 += usize::from(gate_opens(&named, &[p, q, r, j]));
                                let (a, b) = (jn.orient3d(p, q, r, j), js.orient3d(p, q, r, j));
                                if a != b {
                                    diffs.push(format!("orient3d({p},{q},{r};{j}): {a} vs {b}"));
                                }
                            }
                        }
                    }
                }
                for p in 0..n {
                    for a in (p + 1)..n {
                        for b in (a + 1)..n {
                            open3 += usize::from(gate_opens(&named, &[p, a, b]));
                            let (x, y) = (
                                jn.plane_pair_dir_sign(p, a, b),
                                js.plane_pair_dir_sign(p, a, b),
                            );
                            if x != y {
                                diffs.push(format!("dir_sign({p},{a},{b}): {x} vs {y}"));
                            }
                        }
                    }
                }
                assert!(
                    open3 > 0 && open4 > 0,
                    "{arm}: the gate never opened — the differential measured nothing"
                );
                if arm == "mirrored" {
                    let wide = named
                        .iter()
                        .find(|p| p.name_ints.as_ref().is_some_and(|w| w.wide))
                        .expect("a wide plane");
                    assert!(
                        nacre_cip::chain_parity(&wide.tri_pt3[0].chain) < 0,
                        "the mirrored arm must carry an odd chain"
                    );
                }
                assert!(
                    diffs.is_empty(),
                    "{arm}: name route vs climb disagree: {diffs:?}"
                );
            }
        }

        /// ★★★ **The rescue path itself, proved from both sides** — over questions separated
        /// by `2⁻¹⁰⁰`, which the f64 filter cannot see and an 8-bit-capped climb cannot reach:
        /// a starved judge with names answers the hand-derived signs (nothing but the
        /// name-integer route can), and the same starved judge without names answers `0` for
        /// every nonzero question (without the name, this judge *cannot* know). The full-budget
        /// climb on the stripped table is the independent oracle, and the expected signs are
        /// stated by hand, so no route is compared only to itself.
        ///
        /// Also under test: the world gate (`!any_rotated`), the σ fold (plane 3's canonical
        /// name opposes its stored orientation, `frame_sign = −1`), and σ's negative control —
        /// negating every name (and then only some) changes no answer.
        #[test]
        fn a_starved_judge_answers_only_through_the_name() {
            let (m, s) = cuboid();
            let donor = plane_table(&m, s)[0].surf;
            let e_den: i128 = 1 << 100;
            let one_plus = Rat::new(e_den + 1, e_den).unwrap(); // 1 + 2⁻¹⁰⁰
            let one_minus = Rat::new(e_den - 1, e_den).unwrap();
            let r = Rat::from_int;

            // Integer spellings of the five planes' canonical names, scaled by a ~201-bit odd
            // factor so every name is genuinely wide (the predicates are row-linear, so the
            // scale is invisible to every answer — C1's lock).
            let s_factor: BigInt = (BigInt::from(1) << 201) + 3;
            let wide_name = |c: [i128; 4]| -> PlaneName {
                PlaneName::Wide(c.map(|x| BigInt::from(x) * &s_factor))
            };

            // (witness triangle wound to the tri normal; stored normal; frame_sign; name)
            type Fixture = ([[Rat; 3]; 3], [f64; 3], i8, PlaneName);
            let fixtures: [Fixture; 5] = [
                // 0: x = 1
                (
                    [[r(1), r(0), r(0)], [r(1), r(1), r(0)], [r(1), r(0), r(1)]],
                    [1.0, 0.0, 0.0],
                    1,
                    wide_name([1, 0, 0, -1]),
                ),
                // 1: x = 1 + 2⁻¹⁰⁰
                (
                    [
                        [one_plus, r(0), r(0)],
                        [one_plus, r(1), r(0)],
                        [one_plus, r(0), r(1)],
                    ],
                    [1.0, 0.0, 0.0],
                    1,
                    wide_name([e_den, 0, 0, -(e_den + 1)]),
                ),
                // 2: y = 1
                (
                    [[r(0), r(1), r(0)], [r(0), r(1), r(1)], [r(1), r(1), r(0)]],
                    [0.0, 1.0, 0.0],
                    1,
                    wide_name([0, 1, 0, -1]),
                ),
                // 3: z = 1, stored the other way round — the σ fold and the frame_sign bridge.
                (
                    [[r(0), r(0), r(1)], [r(1), r(0), r(1)], [r(0), r(1), r(1)]],
                    [0.0, 0.0, -1.0],
                    -1,
                    wide_name([0, 0, 1, -1]),
                ),
                // 4: x + 2⁻¹⁰⁰·y = 1 — nearly parallel to plane 0.
                (
                    [
                        [r(1), r(0), r(0)],
                        [one_minus, r(1), r(0)],
                        [r(1), r(0), r(1)],
                    ],
                    [1.0, 0.0, 0.0],
                    1,
                    wide_name([e_den, 1, 0, -e_den]),
                ),
            ];

            let build = |names_sign: [i8; 5], keep_names: bool| -> Vec<WorkingPlane> {
                fixtures
                    .iter()
                    .zip(names_sign)
                    .map(|((bases, n_stored, fs, name), flip)| {
                        let name = if flip < 0 {
                            let PlaneName::Wide(c) = name else {
                                unreachable!()
                            };
                            PlaneName::Wide(core::array::from_fn(|i| -&c[i]))
                        } else {
                            name.clone()
                        };
                        let tri_pt3: [WitnessPoint; 3] = bases.map(WitnessPoint::at);
                        let tri = std::array::from_fn(|i| Point3::from_array(tri_pt3[i].coord));
                        WorkingPlane {
                            base_rat: None,
                            // ★ Deliberately not `reconcile`: a deep-rational witness's f64
                            // `tri` is a rounded cache, and reconciling against it would label
                            // the *rounded* plane exact — the two-descriptions trap. `None` is
                            // what forces every question here onto the name or the climb.
                            exact_coeffs: None,
                            exact_normal: None,
                            name_ints: keep_names
                                .then(|| name_stored_ints(Some(&name), &tri_pt3, *fs))
                                .flatten(),
                            base: crate::planes::BaseFrame::none(),
                            surf: donor,
                            plane: nacre_geom::Plane::from_point_normal(
                                tri[0],
                                nacre_math::Vector3::from_array(*n_stored),
                            )
                            .unwrap(),
                            tri,
                            tri_pt3,
                            rotated: false,
                            frame_sign: *fs,
                        }
                    })
                    .collect()
            };
            let named = build([1; 5], true);
            assert!(
                named
                    .iter()
                    .all(|p| p.name_ints.as_ref().is_some_and(|n| n.wide)),
                "every fixture name folds and is wide"
            );
            let stripped = build([1; 5], false);
            let negated = build([-1; 5], true);
            let mixed = build([1, -1, 1, 1, -1], true);

            let starved = Standard {
                prec: 8,
                coincidence: Mag::pow2(-4000),
                scale: Mag::of(16.0),
                cap: 8,
            };
            let notes = Notes::new();
            let oracle = crate::planes::test_judge(&stripped);

            // (question, expected) — hand-derived, so no route is its own oracle.
            // q1: ∩(1,2,3) = (1+2⁻¹⁰⁰, 1, 1) vs ∩(0,2,3) = (1, 1, 1) on axis 0 → +1.
            // q2: (1,1,1) against plane 1 (outward +x, a hair above) → −1.
            // q3: (1+2⁻¹⁰⁰,1,1) against plane 0 → +1.
            // q4: det of stored normals [(1,0,0), (1,2⁻¹⁰⁰,0), (0,0,−1)] = −2⁻¹⁰⁰ → −1.
            // q5: det of stored normals [(1,0,0), (1,2⁻¹⁰⁰,0), (0,1,0)] = 0 — the exact-zero
            //     case, where every route must say 0.
            type Q = (
                &'static str,
                Box<dyn Fn(&Judge<'_, WorkingPlane>) -> i8>,
                i8,
            );
            let questions: Vec<Q> = vec![
                ("cmp", Box::new(|j| j.cmp_coord([1, 2, 3], [0, 2, 3], 0)), 1),
                ("orient a", Box::new(|j| j.orient3d(0, 2, 3, 1)), -1),
                ("orient b", Box::new(|j| j.orient3d(1, 2, 3, 0)), 1),
                ("dir a", Box::new(|j| j.plane_pair_dir_sign(0, 4, 3)), -1),
                ("dir zero", Box::new(|j| j.plane_pair_dir_sign(0, 4, 2)), 0),
            ];
            for (label, ask, expected) in &questions {
                // The name route, too starved to climb: only the integers can answer.
                for (t, arm) in [(&named, "named"), (&negated, "negated"), (&mixed, "mixed")] {
                    let jd = Judge::new(t, starved, &notes);
                    assert_eq!(ask(&jd), *expected, "{label}: starved {arm} table");
                }
                // The independent oracle: the full-budget climb, no names anywhere.
                assert_eq!(ask(&oracle), *expected, "{label}: full climb oracle");
                // The other side of the proof: starved and nameless, the judge cannot know.
                let js = Judge::new(&stripped, starved, &notes);
                assert_eq!(
                    ask(&js),
                    0,
                    "{label}: a starved judge without names must come back empty-handed"
                );
            }
        }
    }
}
