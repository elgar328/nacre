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
    /// and require the answers to **change**: reverse each base triangle (which reverses the
    /// normal `orient3d`'s convention is tied to) and negate each base plane (which reverses the
    /// determinant `plane_pair_dir_sign` takes). If either sweep ever stops disagreeing, the
    /// shortcut has stopped firing and the tests above have gone vacuous.
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
        let mut orient_changed = 0usize;
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
                        if crate::planes::test_judge(&pb).orient3d(p, q, r, j)
                            != crate::planes::test_judge(&reversed).orient3d(p, q, r, j)
                        {
                            orient_changed += 1;
                        }
                    }
                }
            }
        }
        assert!(
            orient_changed > 0,
            "`orient3d` never read the base triangle"
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
}
