//! The b-rep side of the toleranced predicates: `nacre-ops`'s arrangement tables
//! ([`PlaneGeom`], [`FaceInfo`]) implement the [`Witness`]/[`PlaneWitness`] ports so the
//! rotation-general sign predicates in [`nacre_cip::predicate`] can run over them.
//!
//! The tables are **pure description** — geometry and provenance, nothing about how this
//! operation judges. That belongs to [`Judge`], which the engine builds once per boolean
//! (`plane_index_setup` supplies the standard and the collector) and passes down; the predicates
//! are its methods. The predicate logic itself (exact-vs-kernel routing, the frame3 judges) lives
//! in `nacre-cip`.

use crate::planes::{FaceInfo, PlaneGeom};
use nacre_cip::Pt3;
use nacre_cip::predicate::{PlaneWitness, Witness};
use nacre_math::Point3;

// The judging context and the helpers the engine reaches for by name. The predicates themselves
// are methods on `Judge`, so there is nothing else to re-export.
pub(crate) use nacre_cip::predicate::{ImplicitPoint, Judge, any_rotated, plane_def};

impl Witness for PlaneGeom {
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
    fn tri_pt3(&self) -> &[Pt3; 3] {
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
    fn normals_independent(planes: &[PlaneGeom], p: usize, q: usize, r: usize) -> bool {
        let n = |k: usize| planes[k].plane.normal();
        n(p).dot(n(q).cross(n(r))).abs() > 0.5
    }

    /// Core: `orient3d` is rigid-rotation invariant, so the frame3 path over a rotated
    /// cuboid's exact `Pt3` definitions must agree, on every definite config, with the
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

    /// `Judge::planes_coplanar` decides plane identity from the faces' own coordinates, so it must
    /// (a) refuse to conclude anything from a degenerate `tri` — three equal points lie on *every*
    /// plane, and merging on that evidence would fuse genuinely different planes — and
    /// (b) answer the same for a rotated solid as for the unrotated one.
    #[test]
    fn planes_coplanar_guards_degeneracy_and_survives_rotation() {
        let (mut m, s) = cuboid();
        let pu = plane_table(&m, s);
        // (a) A hand-built degenerate pair: same-normal parallel planes, but `tri` is a point.
        let degenerate: Vec<PlaneGeom> = (0..2)
            .map(|k| PlaneGeom {
                // A hand-built table has no recorded coefficients; the composed-rotation route
                // declines and the fixture takes the same escalating path it always did.
                base_rat: None,
                base: crate::planes::BaseFrame::none(),
                surf: pu[0].surf,
                plane: pu[0].plane,
                tri: [Point3::from_array([k as f64, 0.0, 0.0]); 3],
                tri_pt3: std::array::from_fn(|_| Pt3::exact([k as f64, 0.0, 0.0]).expect("exact")),
                rotated: false,
                frame_sign: pu[0].frame_sign,
                // Derived by the same rule the arrangement uses -- a fixture that routed
                // differently would be testing a different engine.
                exact_coeffs: PlaneGeom::reconcile(
                    &pu[0].plane,
                    [Point3::from_array([k as f64, 0.0, 0.0]); 3],
                    false,
                )
                .0,
                exact_normal: PlaneGeom::reconcile(
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
    fn two_spellings() -> (Vec<PlaneGeom>, Vec<PlaneGeom>) {
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
                // engine come to route differently (see `PlaneGeom::reconcile`).
                base_rat: None,
                surf: f.surf,
                face: f.face,
                plane: f.plane,
                tri: f.tri,
                n_out: f.n_out,
                orient_sign: f.orient_sign,
                tri_pt3: std::array::from_fn(|i| {
                    let c = f.tri_pt3[i].coord;
                    let r = |v: f64| Rat::try_from_f64(v).expect("integer cuboid coordinate");
                    Pt3::at([r(2.0 - c[0]), r(c[1]), r(c[2])]).mirror(Axis::X, Rat::from_int(1))
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
    /// `Unsupported { SeamAlias }` — the same family of defect, one level down. `Single` differs by
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

    /// ★★★★★ **The defect, pinned.** Two bosses on a tilted face reach the same height by
    /// different arithmetic, and the kernel calls their caps **different planes**.
    ///
    /// The heights are `2.0 + 7.7` against `2.0 + 1.1 + 6.6`. In the rationals those are the same
    /// number — `77/10` either way, which is why the axis-aligned path merges them — but on a
    /// tilted face the cap is reached through an irrational frame, and what the judge compares is
    /// the **witness triangle**, whose f64 coordinates land one ulp apart.
    ///
    /// ★ The control run (`Single`, both columns raised in one pad) merges, so the split is the
    /// arithmetic and not the tilt. What this test becomes once the witness is exact: `Split`
    /// merges too. See `docs/dev-log.md`.
    #[test]
    fn two_caps_reached_by_different_arithmetic_are_two_planes_today() {
        for (recipe, want_same) in [(Recipe::Split, false), (Recipe::Single, true)] {
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
            // ★ On its own the split costs nothing: the solid is right either way.
            assert!(nacre_validate::validate(&m).is_empty());
        }
    }

    /// ★★★★★ **The split is not conservative — it builds a malformed solid, silently.**
    ///
    /// A plate laid across **both** caps asks the boolean a question it can only answer correctly
    /// if the caps are one plane. With the caps split, the result comes back with a face whose
    /// stored surface normal is **76° off its own outward normal** — `collect_planes` asserts
    /// exactly that invariant, so the solid cannot be used as an operand again, and in release
    /// `orient_sign` would be read off that dot.
    ///
    /// ★ **Nothing reports it**: the volume is right and `validate` is clean.
    ///
    /// The cause is isolated by three controls — the bad face appears only when a plate spans two
    /// caps that are in different classes:
    ///
    /// | recipe | plate | malformed faces |
    /// |---|---|---|
    /// | `Split` | across both caps | **1** |
    /// | `Split` | over one cap | 0 |
    /// | `Single` | across both caps | 0 |
    /// | `Single` | over one cap | 0 |
    #[test]
    fn a_plate_across_two_split_caps_builds_a_malformed_solid_today() {
        for (recipe, span, want_bad) in [
            (Recipe::Split, true, 1),
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
        ]);
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
                let nacre_geom::Surface::Plane(pl) = m.surfaces.get(face.surface) else {
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
    /// plan would store: exact **frame** coefficients (`Model::surface_coeffs`) and the motion that
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
        use nacre_cip::{Pt3, Standard};
        use nacre_scalar::Bound;

        let (m, s, _) = two_caps_on_a_tilted_face(Recipe::Split);
        let faces = collect_planes(&m, s).unwrap();
        let (ia, ib) = cap_pair(&m, &faces);

        // The exact definition each cap already carries: frame coefficients + the motion.
        // `nudge` offsets the cap's height inside its own frame — the negative control.
        let define = |i: usize, nudge: Option<Rat>| -> [Pt3; 3] {
            let surf = faces[i].surf;
            let c = *m
                .surface_coeffs
                .get(&surf)
                .expect("a frame cap records its coefficients");
            let motion = match m.surface_defs.get(&surf) {
                Some(nacre_topo::SurfaceDef::Moved { motion, .. }) => *motion,
                other => panic!("cap {i} is not Moved: {other:?}"),
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
            [[zero, zero, z], [one, zero, z], [zero, one, z]]
                .map(|b| crate::rotated_vertex::replay(Pt3::at(b), &chain).expect("replay"))
        };
        let (da, db) = (define(ia, None), define(ib, None));

        // A two-plane table over those definitions. `base_rat: None` keeps the composed-rotation
        // shortcut out of it, so what runs is the interval route the plan's stage C1 exercises.
        let mk = |d: [Pt3; 3]| {
            let tri = d.clone().map(|p| Point3::from_array(p.coord));
            PlaneGeom {
                base_rat: None,
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
                coincidence: Bound::pow2(-4000),
                scale: Bound::of(16.0),
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
            // `Bound` is `m · 2^e`; the exponent is the reading that matters — a true zero makes it
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
                coincidence: Bound::pow2(-4000),
                scale: Bound::of(16.0),
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
}
