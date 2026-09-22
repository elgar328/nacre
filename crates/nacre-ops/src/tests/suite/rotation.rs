//! Booleans on rotated operands: witness invariance, reuse of a rotated result, rotation stress.

use super::*;

// A boolean commutes with a rigid motion, so rotating both operands by the same
// irrational-angle isometry must give the rigid image of the unrotated result — identical
// volume, solid count, and cavity count, and still valid. These are the first live proof
// that the CIP-wired machinery (arrangement, seam, in/out, outer/cavity — 3a–3c-vi) is
// sound end-to-end on rotated (rounded-irrational) geometry.

/// A cutter's convex corner landing exactly on the target's concave corner — three planes
/// (x=1,y=1,z=1) meeting at one point (1,1,1), six faces there — makes a **non-manifold pinch**
/// (two face-fans touching at the point). No valid 2-manifold solid exists, so `boolean` rejects
/// with the clear `NonManifoldVertex` reason (not the incidental `EulerParity`), leaving the
/// live model untouched. **Rotation-independent**: both the axis-aligned and the rotated framings
/// (exact and CIP-kernel paths) hit the same pinch and reject. R = a cube minus a far-corner
/// octant; C = a cube whose +corner is that removed octant's inner corner. (Coplanar contact away
/// from a corner works — [`a_rotated_boolean_result_can_be_cut_again`].)
#[test]
fn a_corner_coincident_cut_is_rejected_not_silently_wrong() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    // `rotate`: None = axis-aligned (exact path); Some = the result and cutter tilted (CIP path).
    let run = |rotate: bool| {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        // C's +corner is (1,1,1) = R's concave corner → three shared planes meet there.
        let c = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([1.0; 3]),
        );
        m.rebuild_adjacency();
        let (r, c) = if rotate {
            let r = transform(&mut m, r, &iso).unwrap();
            m.rebuild_adjacency();
            let c = transform(&mut m, c, &iso).unwrap();
            m.rebuild_adjacency();
            (r, c)
        } else {
            (r, c)
        };
        let live = m.live_solids().to_vec();
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, r, c),
            RejectReason::NonManifoldVertex,
        );
        assert_eq!(
            m.live_solids().to_vec(),
            live,
            "reject must not mutate the live set"
        );
    };
    run(false); // axis-aligned
    run(true); // rotated
}

/// Two cubes touching only at the corner (1,1,1): nothing is joined, so the Fuse comes back as
/// the two bodies it was handed. ★ This used to be `NON_MANIFOLD_VERTEX` — true of the single
/// welded body the reconstruction built then, and beside the point once the pieces are minted
/// per solid. The pinch reject is still there for a body that touches *itself*
/// (`tests/probes/contact_separates.rs`).
#[test]
fn two_cubes_touching_at_a_corner_fuse_to_two_bodies() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("a point contact separates");
    assert_eq!(out.len(), 2);
    m.rebuild_adjacency();
    assert_eq!(nacre_validate::validate(&m), Vec::new());
}

/// The kernel's validator catches the pinch too. Bypass `boolean`'s reject via
/// `arrangement::boolean` to obtain a malformed
/// solid, then confirm `validate` reports a `NonManifoldVertex`.
///
/// ★ The shape has to be one that **cannot** part: two bodies touching at a corner now come
/// back as two bodies and there is no malformed solid to obtain from them. So A and B meet only
/// at `(2,2,1)` while two bridges run around the contact and join them elsewhere — the material
/// loops, the pinch is real, and the vertex check (which lives in `boolean`, not in the
/// arrangement) is the one this bypasses.
#[test]
fn validate_reports_the_non_manifold_pinch() {
    let mut m = Model::new();
    let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
        let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        s
    };
    let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let g1 = cub(&mut m, [1.0, 0.3, 0.2], [4.5, 1.3, 0.8]);
    let g2 = cub(&mut m, [3.5, 0.3, 0.2], [4.5, 3.8, 1.6]);
    let b = cub(&mut m, [2.0, 2.0, 1.0], [4.0, 4.0, 2.0]);
    let t1 = boolean(&mut m, BoolKind::Fuse, a, g1).expect("a and g1 overlap");
    m.rebuild_adjacency();
    let t2 = boolean(&mut m, BoolKind::Fuse, t1[0], g2).expect("g1 and g2 overlap");
    m.rebuild_adjacency();
    crate::arrangement::boolean(&mut m, BoolKind::Fuse, t2[0], b).unwrap();

    m.rebuild_adjacency();
    let issues = nacre_validate::validate(&m);
    assert!(
        issues
            .iter()
            .any(|v| matches!(v, nacre_validate::Violation::NonManifoldVertex { .. })),
        "validate must flag the pinch: {issues:?}"
    );
}

/// Predicates over a rotated result's witness planes are rotation-invariant against the same
/// result unrotated — the provenance witness (with outward winding) defines the exact plane,
/// so `Judge::orient3d` agrees on every definite triple. A regression guard on the witness itself.
#[test]
fn rotated_result_witness_predicates_are_invariant() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let table = |m: &Model, s: Handle<Solid>| {
        let f = collect_planes(m, s).unwrap();
        let c = plane_classes(&crate::planes::test_judge(&f));
        dense_planes(&f, &c).0
    };
    let pu = table(&m, r);
    let r2 = transform(&mut m, r, &iso).unwrap();
    m.rebuild_adjacency();
    let pr = table(&m, r2);
    assert_eq!(pu.len(), pr.len(), "rotation preserves the plane count");
    let n = pu.len();
    let indep = |p: &[WorkingPlane], a: usize, b: usize, c: usize| {
        let nrm = |k: usize| p[k].plane.normal();
        nrm(a).dot(nrm(b).cross(nrm(c))).abs() > 0.3
    };
    let mut disagree = 0;
    for p in 0..n {
        for q in (p + 1)..n {
            for rr in (q + 1)..n {
                if !indep(&pu, p, q, rr) {
                    continue;
                }
                for j in 0..n {
                    if j == p || j == q || j == rr {
                        continue;
                    }
                    let su = crate::planes::test_judge(&pu).orient3d(p, q, rr, j);
                    let sr = crate::planes::test_judge(&pr).orient3d(p, q, rr, j);
                    if su != 0 && sr != 0 && su != sr {
                        eprintln!("DISAGREE orient3d ({p},{q},{rr},{j}): u={su} r={sr}");
                        disagree += 1;
                    }
                }
            }
        }
    }
    assert_eq!(
        disagree, 0,
        "{disagree} predicate disagreements (witness wrong)"
    );
}

/// Rotating a boolean *result* and feeding it back into a boolean: `collect_planes`
/// witnesses each rotated seam face's plane through provenance — its
/// plane is `R(π)` for an operand plane `π`, recovered from the operand face still on `π` and
/// rotated by the face's own chain. A boolean commutes with a rigid motion, so the rotated
/// chain's result matches the unrotated chain's (volume, solid count) and stays valid.
#[test]
fn a_rotated_boolean_result_can_be_cut_again() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    // Chain: R = Cut(A, B) removes a far-corner octant; then Cut(R, C) removes a near one.
    let build = |m: &mut Model| -> (Handle<Solid>, Handle<Solid>) {
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let r = boolean_one(m, BoolKind::Cut, a, b).unwrap();
        // A clean slab that severs R at x = 0.5 — no plane of C coincides with any of R's
        // (avoids the separate rotated-coplanar-contact gap; isolates the witness).
        let c = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([0.5, 4.0, 4.0]),
        );
        (r, c)
    };
    // Unrotated reference (reuse already works when nothing is rotated).
    let mut m0 = Model::new();
    let (r0, c0) = build(&mut m0);
    m0.rebuild_adjacency();
    let ref_out = boolean(&mut m0, BoolKind::Cut, r0, c0).unwrap();
    m0.rebuild_adjacency();
    let ref_vol: f64 = ref_out
        .iter()
        .map(|&s| nacre_props::mass_props(&m0, s).unwrap().volume)
        .sum();
    // Rotated: turn the *result* R (and C) by the same isometry, then reuse R.
    let mut m = Model::new();
    let (r, c) = build(&mut m);
    m.rebuild_adjacency();
    let r = transform(&mut m, r, &iso).unwrap();
    m.rebuild_adjacency();
    let c = transform(&mut m, c, &iso).unwrap();
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Cut, r, c).unwrap();
    m.rebuild_adjacency();
    let issues = nacre_validate::validate(&m);
    assert!(
        issues.is_empty(),
        "rotated-result reuse must be valid: {issues:?}"
    );
    let vol: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert_eq!(
        out.len(),
        ref_out.len(),
        "solid count invariant under rotation"
    );
    assert!(
        (vol - ref_vol).abs() < 1e-6,
        "rotated reuse volume {vol} vs unrotated {ref_vol}"
    );
}

/// Result-reuse rotation stress: build R with a first boolean, then feed R into a second
/// boolean with a fresh cutter C — once unrotated, once with R and C rotated by the same
/// isometry. A boolean commutes with a rigid motion, so the rotated reuse must equal the
/// unrotated one (volume, solid count, cavity count) or be an honest reject — never silently
/// wrong. Every vertex of a boolean result is arrangement-named, so its rotation exercises the
/// provenance witness on every face. `#[ignore]`: rotated booleans escalate to astro-float (~1–3
/// s).
#[test]
#[ignore = "slow: rotated result-reuse booleans (run with --ignored)"]
fn rotated_result_reuse_stress() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    #[derive(PartialEq, Debug)]
    enum Out {
        Rej,
        Ok(f64, usize, usize),
    }
    let rot = |axis: Axis, deg: i128, piv: [i128; 3]| {
        Isometry::rotation(Rotation {
            axis,
            pivot: [
                Rat::from_int(piv[0]),
                Rat::from_int(piv[1]),
                Rat::from_int(piv[2]),
            ],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    };
    // (first kind, second kind, |m| -> (a, b, c)). R = kind1(a, b); out = kind2(R, c).
    type Build = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>, Handle<Solid>)>;
    let cuboid = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
        m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
    };
    let fixtures: Vec<(&str, BoolKind, BoolKind, Build)> = vec![
        (
            "corner_then_slab",
            BoolKind::Cut,
            BoolKind::Cut,
            Box::new(move |m: &mut Model| {
                let a = cuboid(m, [0.0; 3], [2.0; 3]);
                let b = cuboid(m, [1.0; 3], [3.0; 3]);
                let c = cuboid(m, [-1.0, -1.0, -1.0], [0.5, 4.0, 4.0]);
                (a, b, c)
            }),
        ),
        (
            "fuse_then_bite",
            BoolKind::Fuse,
            BoolKind::Cut,
            Box::new(move |m: &mut Model| {
                let a = cuboid(m, [0.0; 3], [2.0, 1.0, 1.0]);
                let b = cuboid(m, [0.0, 0.0, 0.0], [1.0, 2.0, 1.0]);
                let c = cuboid(m, [1.3, 1.3, -1.0], [3.0, 3.0, 2.0]);
                (a, b, c)
            }),
        ),
        (
            "cut_then_fuse",
            BoolKind::Cut,
            BoolKind::Fuse,
            Box::new(move |m: &mut Model| {
                let a = cuboid(m, [0.0; 3], [2.0; 3]);
                let b = cuboid(m, [1.3, 1.3, 1.3], [3.0, 3.0, 3.0]);
                let c = cuboid(m, [-0.7, 0.4, 0.4], [0.3, 1.4, 1.4]);
                (a, b, c)
            }),
        ),
    ];
    let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
        ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
        ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
        (
            "Z50>Y37",
            vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::Y, 37, [0, 0, 1])],
        ),
    ];
    let run = |k1: BoolKind, k2: BoolKind, build: &Build, isos: &[Isometry]| -> Out {
        let mut m = Model::new();
        let (a, b, c) = build(&mut m);
        m.rebuild_adjacency();
        let Ok(r) = boolean(&mut m, k1, a, b) else {
            return Out::Rej;
        };
        assert_eq!(r.len(), 1, "first boolean is a single solid");
        let mut r = r[0];
        let mut c = c;
        m.rebuild_adjacency();
        for iso in isos {
            r = transform(&mut m, r, iso).unwrap();
            m.rebuild_adjacency();
            c = transform(&mut m, c, iso).unwrap();
            m.rebuild_adjacency();
        }
        match boolean(&mut m, k2, r, c) {
            Ok(solids) => {
                m.rebuild_adjacency();
                assert!(
                    nacre_validate::validate(&m).is_empty(),
                    "INVALID rotated reuse result"
                );
                let vol: f64 = solids
                    .iter()
                    .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                    .sum();
                let cav: usize = solids.iter().map(|&s| m.solid(s).cavities.len()).sum();
                Out::Ok(vol, solids.len(), cav)
            }
            Err(_) => Out::Rej,
        }
    };
    let vclose =
        |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());
    let (mut success, mut reject, mut silent) = (0, 0, 0);
    for (fname, k1, k2, build) in &fixtures {
        let base = run(*k1, *k2, build, &[]);
        for (rname, isos) in &isos_list {
            let r = run(*k1, *k2, build, isos);
            match (&base, &r) {
                (_, Out::Rej) => reject += 1,
                (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2))
                    if vclose(*v1, *v2) && s1 == s2 && c1 == c2 =>
                {
                    success += 1
                }
                _ => {
                    silent += 1;
                    eprintln!("SILENT-WRONG {fname} {rname} base={base:?} rot={r:?}");
                }
            }
        }
    }
    eprintln!("REUSE STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent}");
    assert_eq!(
        silent, 0,
        "a rotated result-reuse boolean was silently wrong"
    );
    assert!(
        success >= 1,
        "at least one rotated reuse must actually succeed"
    );
}

/// Adversarial rotation stress (overhaul 3d-iii): many fixtures × kinds × rotations
/// (single-axis, and Euler chains reaching arbitrary orientation) confirm the DNA
/// invariant — a rotated boolean is *never silently wrong*: its result either equals the
/// unrotated one (a boolean commutes with a rigid motion, so volume/solid-count/cavity-count
/// are invariant) or is an honest reject. `#[ignore]`: each rotated boolean escalates its
/// CIP predicates to astro-float and costs ~0.5–2.5 s, so this runs on demand, not per commit
/// (the invariance regression guard is the fast `rotated_*` tests above).
#[test]
#[ignore = "slow: rotated booleans ~2s each (run with --ignored)"]
fn rotation_invariance_stress() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    #[derive(PartialEq, Debug)]
    enum Out {
        Rej,
        Ok(f64, usize, usize),
    }
    let rot = |axis: Axis, deg: i128, piv: [i128; 3]| {
        Isometry::rotation(Rotation {
            axis,
            pivot: [
                Rat::from_int(piv[0]),
                Rat::from_int(piv[1]),
                Rat::from_int(piv[2]),
            ],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    };
    let run = |build: &dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>),
               kind: BoolKind,
               isos: &[Isometry]|
     -> Out {
        let (mut m, mut a, mut b) = build();
        for iso in isos {
            a = transform(&mut m, a, iso).unwrap();
            m.rebuild_adjacency();
            b = transform(&mut m, b, iso).unwrap();
            m.rebuild_adjacency();
        }
        match boolean(&mut m, kind, a, b) {
            Ok(solids) => {
                m.rebuild_adjacency();
                assert!(
                    nacre_validate::validate(&m).is_empty(),
                    "INVALID rotated result"
                );
                let vol: f64 = solids
                    .iter()
                    .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                    .sum();
                let cav: usize = solids.iter().map(|&s| m.solid(s).cavities.len()).sum();
                Out::Ok(vol, solids.len(), cav)
            }
            Err(_) => Out::Rej,
        }
    };
    let vclose =
        |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());
    let matches = |base: &Out, r: &Out| match (base, r) {
        (Out::Rej, Out::Rej) => true,
        (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2)) => vclose(*v1, *v2) && s1 == s2 && c1 == c2,
        _ => false,
    };
    let cube = |lo: [f64; 3], hi: [f64; 3]| {
        move || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
            (m, a, b)
        }
    };
    type Build = Box<dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>)>;
    let fixtures: Vec<(&str, Build)> = vec![
        ("corner", Box::new(l_and_corner_box)),
        (
            "rod",
            Box::new(|| {
                let (m, l, r) = l_and_rod();
                (m, r, l) // sever = Cut(rod, l): a=rod, b=l
            }),
        ),
        ("inner", Box::new(l_and_inner_box)),
        ("cube_corner", Box::new(cube([0.5; 3], [1.5; 3]))),
    ];
    let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
        ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
        ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
        (
            "Z50>X37",
            vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::X, 37, [0, 0, 1])],
        ),
        // A three-axis Euler chain reaches an arbitrary orientation (axes are X/Y/Z only).
        (
            "Z30>X30>Y73",
            vec![
                rot(Axis::Z, 30, [1, 1, 0]),
                rot(Axis::X, 30, [0, 0, 0]),
                rot(Axis::Y, 73, [0, 2, 0]),
            ],
        ),
    ];
    let (mut success, mut reject, mut silent, mut skipped) = (0, 0, 0, 0);
    for (fname, build) in &fixtures {
        for kind in [BoolKind::Cut, BoolKind::Fuse, BoolKind::Common] {
            let base = run(build.as_ref(), kind, &[]);
            for (rname, isos) in &isos_list {
                if matches!(base, Out::Rej) {
                    skipped += 1;
                    continue;
                }
                let r = run(build.as_ref(), kind, isos);
                if matches!(r, Out::Rej) {
                    reject += 1; // an honest reject is acceptable, not a counterexample
                } else if matches(&base, &r) {
                    success += 1;
                } else {
                    silent += 1;
                    eprintln!("SILENT-WRONG {fname} {kind:?} {rname} base={base:?} rot={r:?}");
                }
            }
        }
    }
    eprintln!(
        "ROTATION STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent} skipped(base_rej)={skipped}"
    );
    assert_eq!(
        silent, 0,
        "a rotated boolean was silently wrong (valid but != unrotated)"
    );
}
