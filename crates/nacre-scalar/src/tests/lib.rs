/// **A cylinder's unit frame is its own statement when the frame is already orthonormal** —
/// and it declines rather than projecting when `ref_dir` is not perpendicular, because every
/// consumer's soundness rests on `û·dir = 0` being a fact.
#[test]
fn a_cylinder_frame_is_unit_and_perpendicular_or_it_is_nothing() {
    use crate::{Rat, cyl_unit_frame};
    let q = |n: i128| Rat::from_int(n);
    let v = |a: i128, b: i128, c: i128| [q(a), q(b), q(c)];
    // Orthonormal in, and the two axes come back as the statement itself.
    let (u1, u2) = cyl_unit_frame(&v(0, 0, 1), &v(1, 0, 0)).expect("a world frame");
    assert_eq!(u1, v(1, 0, 0));
    assert_eq!(u2, v(0, 1, 0));
    // A non-unit but perpendicular seed still normalizes exactly when the norms are rational.
    let (u1, u2) = cyl_unit_frame(&v(0, 0, 5), &v(3, 4, 0)).expect("a 3-4-5 seed");
    assert_eq!(u1, [Rat::new(3, 5).unwrap(), Rat::new(4, 5).unwrap(), q(0)]);
    assert_eq!(
        u2,
        [Rat::new(-4, 5).unwrap(), Rat::new(3, 5).unwrap(), q(0)]
    );
    // ★ Not perpendicular: declined, **not** projected. The projecting recipe would answer
    // here, and its `û₁` would carry an axis component — a "rim" point off the cap plane.
    assert!(cyl_unit_frame(&v(0, 0, 1), &v(3, 0, 4)).is_none());
    // Parallel, and the zero vector with it (`parallel_rat` says a zero is parallel to all).
    assert!(cyl_unit_frame(&v(0, 0, 1), &v(0, 0, 2)).is_none());
    assert!(cyl_unit_frame(&v(0, 0, 1), &v(0, 0, 0)).is_none());
    assert!(cyl_unit_frame(&v(0, 0, 0), &v(1, 0, 0)).is_none());
    // An irrational norm is the third cause, and it is a decline like the others.
    assert!(cyl_unit_frame(&v(0, 0, 1), &v(1, 1, 0)).is_none());
}

/// `a + b·π` signed exactly through π's bracket: `22/7` and `355/113` straddle π and both
/// decide, a rational alone is its own sign, and the integer door agrees with the `Rat` one.
#[test]
fn the_sign_of_a_plus_b_pi_is_read_off_pis_bracket() {
    let r = Rat::from_int;
    assert_eq!(sign_a_plus_b_pi(r(-3), r(1)), Some(Orient::Positive)); // π > 3
    assert_eq!(sign_a_plus_b_pi(r(-4), r(1)), Some(Orient::Negative)); // π < 4
    assert_eq!(sign_a_plus_b_pi(r(22), r(-7)), Some(Orient::Positive)); // 22/7 > π
    assert_eq!(sign_a_plus_b_pi(r(-355), r(113)), Some(Orient::Negative)); // 355/113 > π
    assert_eq!(sign_a_plus_b_pi(r(0), r(-2)), Some(Orient::Negative)); // −2π
    assert_eq!(sign_a_plus_b_pi(r(5), r(0)), Some(Orient::Positive));
    assert_eq!(sign_a_plus_b_pi(r(0), r(0)), Some(Orient::Zero));
    assert_eq!(
        sign_a_plus_b_pi(Rat::new(-1, 2).unwrap(), Rat::new(1, 7).unwrap()), // π/7 ≈ 0.449 < 1/2
        Some(Orient::Negative)
    );
    assert_eq!(
        sign_a_plus_b_pi(Rat::new(-2, 5).unwrap(), Rat::new(1, 7).unwrap()), // π/7 > 2/5
        Some(Orient::Positive)
    );
    // A whole circle's winding: no straight steps, one full turn — 2·area = 2πr².
    let circle = QuarterArc {
        center: [r(0), r(0)],
        r2: r(25),
        start: [r(5), r(0)],
        end: [r(5), r(0)],
        ccw: true,
    };
    assert_eq!(
        winding_sign_quarter_arcs(&[], &[circle]),
        Some(Orient::Positive)
    );
    let cw = QuarterArc {
        ccw: false,
        ..circle
    };
    assert_eq!(
        winding_sign_quarter_arcs(&[], &[cw]),
        Some(Orient::Negative)
    );
    // A square, and a slot (two lines, two half circles) drawn counter-clockwise.
    let sq = [
        [[r(0), r(0)], [r(4), r(0)]],
        [[r(4), r(0)], [r(4), r(4)]],
        [[r(4), r(4)], [r(0), r(4)]],
        [[r(0), r(4)], [r(0), r(0)]],
    ];
    assert_eq!(winding_sign_quarter_arcs(&sq, &[]), Some(Orient::Positive));
    let slot_lines = [
        [[r(0), r(-5)], [r(30), r(-5)]],
        [[r(30), r(5)], [r(0), r(5)]],
    ];
    let slot_arcs = [
        QuarterArc {
            center: [r(30), r(0)],
            r2: r(25),
            start: [r(30), r(-5)],
            end: [r(30), r(5)],
            ccw: true,
        },
        QuarterArc {
            center: [r(0), r(0)],
            r2: r(25),
            start: [r(0), r(5)],
            end: [r(0), r(-5)],
            ccw: true,
        },
    ];
    assert_eq!(
        winding_sign_quarter_arcs(&slot_lines, &slot_arcs),
        Some(Orient::Positive)
    );
    // A 60° arc is not a quarter-turn multiple: undecided, not guessed.
    let sixty = QuarterArc {
        center: [r(0), r(0)],
        r2: r(4),
        start: [r(2), r(0)],
        end: [r(1), r(1)], // not even on the circle — the shape of the refusal is the same
        ccw: true,
    };
    assert_eq!(winding_sign_quarter_arcs(&[], &[sixty]), None);
}

use super::*;

/// The depth these tests realize at. **Test-local on purpose**: production has no default
/// precision — it is computed from the model — so a constant here must not be reachable from
/// outside the fixtures that chose it.
const GT_PREC: usize = 160;

fn ints(v: [i128; 4]) -> [Rat; 4] {
    v.map(Rat::from_int)
}

/// [`Isometry::fixes_plane`] — the set-invariance predicate the invariant-plane
/// restatement gates on. Hand-known planes and motions; the overflow case locks the
/// conservative direction (a miss records a node — slower, never wrong).
mod fixes_plane {
    use super::*;

    fn rot(axis: Axis, deg: i128, point: [i128; 3]) -> Isometry {
        Isometry::rotation(Rotation {
            axis,
            pivot: point.map(Rat::from_int),
            angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
        })
    }

    #[test]
    fn an_axis_parallel_normal_rides_any_pivot() {
        let z_plane = ints([0, 0, 1, -3]);
        assert!(rot(Axis::Z, 30, [0, 0, 0]).fixes_plane(&z_plane));
        assert!(rot(Axis::Z, 30, [7, -2, 5]).fixes_plane(&z_plane));
        let x_wall = ints([1, 0, 0, -2]);
        assert!(!rot(Axis::Z, 30, [0, 0, 0]).fixes_plane(&x_wall));
        assert!(rot(Axis::X, 30, [1, 1, 1]).fixes_plane(&x_wall));
        assert!(!rot(Axis::X, 30, [1, 1, 1]).fixes_plane(&z_plane));
    }

    #[test]
    fn a_translation_moves_the_plane_iff_it_leaves_it() {
        let z_plane = ints([0, 0, 1, -3]);
        let mut in_plane = rot(Axis::Z, 30, [0, 0, 0]);
        in_plane.translate = [Rat::from_int(1), Rat::from_int(2), Rat::from_int(0)];
        assert!(in_plane.fixes_plane(&z_plane));
        let mut off_plane = rot(Axis::Z, 30, [0, 0, 0]);
        off_plane.translate = [Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)];
        assert!(!off_plane.fixes_plane(&z_plane));
    }

    #[test]
    fn a_pure_translation_needs_no_rotation_clause() {
        // The population the whole-solid probe loses to magnitude rounding: a huge
        // in-plane offset still fixes every plane it slides along.
        let big = Rat::from_int(1i128 << 100);
        let iso = Isometry::translation([big, big, Rat::from_int(0)]);
        assert!(iso.fixes_plane(&ints([0, 0, 1, -3])));
        assert!(!iso.fixes_plane(&ints([1, 0, 0, -2])));
    }

    #[test]
    fn an_overflowing_dot_is_a_conservative_miss() {
        // a·tx and −b·ty would cancel to a true zero, but each product overflows
        // i128 — the checked arithmetic answers false rather than guessing.
        let wide = 1i128 << 100;
        let plane = [
            Rat::from_int(wide),
            Rat::from_int(wide),
            Rat::from_int(0),
            Rat::from_int(-1),
        ];
        let iso =
            Isometry::translation([Rat::from_int(wide), Rat::from_int(-wide), Rat::from_int(0)]);
        assert!(!iso.fixes_plane(&plane));
    }
}

/// [`plane_name_from_meets`] — the name of a plane through meets of any width,
/// locked on both output widths with hand-known planes so the derivation is
/// never its own oracle.
mod name_from_meets {
    use super::*;
    use num_bigint::BigInt;
    use proptest::prelude::*;

    fn r(n: i128, d: i128) -> Rat {
        Rat::new(n, d).unwrap()
    }

    /// The carriers of one fixture vertex: `x = m·2⁻⁸⁰`, `y = w·5⁻⁴⁰`, and a third plane
    /// `base + λ·(x − m·2⁻⁸⁰)` — which meets the first two exactly where `base` does, so
    /// the meet lies on `base` by construction while no carrier *is* `base`.
    fn meet_on(base: [Rat; 4], m: i128, w: i128, lambda: i128) -> MeetPoint {
        let a = PlaneName::Narrow([r(1 << 80, 1), r(0, 1), r(0, 1), r(-m, 1)]);
        let b = PlaneName::Narrow([r(0, 1), r(5i128.pow(40), 1), r(0, 1), r(-w, 1)]);
        let u = r(m, 1 << 80);
        let c = PlaneName::Narrow([
            base[0].checked_add(r(lambda, 1)).unwrap(),
            base[1],
            base[2],
            base[3]
                .checked_sub(r(lambda, 1).checked_mul(u).unwrap())
                .unwrap(),
        ]);
        three_planes_big([&a, &b, &c]).expect("the fixture carriers meet")
    }

    /// ★ Narrow-name output: three **wide** meets on the hand-known plane
    /// `T: 7x + 11y − 13z + 1 = 0` (already canonical — gcd 1, first coefficient
    /// positive) name exactly `T`.
    #[test]
    fn wide_meets_on_a_narrow_plane_name_it() {
        let t = [r(7, 1), r(11, 1), r(-13, 1), r(1, 1)];
        let meets = [
            meet_on(t, 1, 1, 1),
            meet_on(t, 3, 1, 2),
            meet_on(t, 1, 7, 3),
        ];
        for m in &meets {
            assert!(
                matches!(m, MeetPoint::Wide(_)),
                "fixture validity: every meet must be wide, or this test measures nothing"
            );
        }
        assert_eq!(
            plane_name_from_meets([&meets[0], &meets[1], &meets[2]]),
            Some(PlaneName::Narrow([r(7, 1), r(11, 1), r(-13, 1), r(1, 1)])),
            "the plane the meets lie on, by hand"
        );
    }

    /// ★ Wide-name output — the other side of the derivation's fork. `T_w: 2⁷⁰·x + y/5³⁰
    /// − z + 1 = 0` has all-`Rat` coefficients but its canonical integers are `×5³⁰`:
    /// `(2⁷⁰·5³⁰, 1, −5³⁰, 5³⁰)`, whose largest is ~2¹⁴⁰ — a `Wide` name, still known by
    /// hand.
    #[test]
    fn wide_meets_on_a_wide_plane_name_it() {
        let tw = [r(1 << 70, 1), r(1, 5i128.pow(30)), r(-1, 1), r(1, 1)];
        let meets = [
            meet_on(tw, 1, 1, 1),
            meet_on(tw, 3, 1, 2),
            meet_on(tw, 1, 7, 3),
        ];
        for m in &meets {
            assert!(matches!(m, MeetPoint::Wide(_)), "fixture validity");
        }
        let five30 = BigInt::from(5i128.pow(30));
        let expected = [
            BigInt::from(1i128 << 70) * &five30,
            BigInt::from(1),
            -&five30,
            five30.clone(),
        ];
        assert_eq!(
            plane_name_from_meets([&meets[0], &meets[1], &meets[2]]),
            Some(PlaneName::Wide(expected)),
            "the canonical integers of T_w, by hand"
        );
    }

    /// Collinear meets name nothing — width is never the cause of a `None` here.
    #[test]
    fn collinear_meets_name_no_plane() {
        let p = |x: i128| MeetPoint::Narrow([r(x, 1), r(x, 1), r(0, 1)]);
        let (a, b, c) = (p(0), p(1), p(2));
        assert_eq!(plane_name_from_meets([&a, &b, &c]), None);
    }

    proptest! {
        #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
            proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
        ))]
        /// All-narrow meets: the new derivation answers exactly what [`plane_name_exact`]
        /// answers — the independent oracle for the lift being a pure generalization.
        #[test]
        fn narrow_meets_agree_with_the_point_derivation(
            pts in proptest::array::uniform3(proptest::array::uniform3((-40i128..=40, 1u32..=40))),
        ) {
            let p: [[Rat; 3]; 3] = pts.map(|q| q.map(|(n, e)| r(n, 1i128 << e)));
            let meets = p.map(MeetPoint::Narrow);
            prop_assert_eq!(
                plane_name_from_meets([&meets[0], &meets[1], &meets[2]]),
                plane_name_exact(p[0], p[1], p[2])
            );
        }
    }
}

/// The integer sign predicates, differentially locked.
///
/// Three locks, matching the three ways they can be wrong:
/// - **against the f64 twins** on narrow input — same convention, two vessels;
/// - **scale-invariance into wide width** — the sign predicates are row-linear, so a row
///   times a 200-bit positive integer is a genuine `> i128` input whose true answer the
///   narrow twin still knows. This is the only exact ground truth for wide inputs that does
///   not test `f(x)` against `f(x)`;
/// - **negation as the direction negative-control** — exactly the direction-sensitive
///   slots flip (`plane_side`: `j` only; `cmp_coord`: none; `dir_sign`: every row), which
///   is the per-slot σ analysis the judging layer's rescue relies on.
mod int_predicates {
    use super::*;
    use nacre_predicates::{ThreePlane, det3_sign, indirect_cmp_coord, indirect_plane_side};
    use num_bigint::BigInt;
    use proptest::prelude::*;

    /// Small integer coefficients — exact in `f64`, so the twins speak the same input.
    fn coeffs() -> impl Strategy<Value = [i32; 4]> {
        proptest::array::uniform4(-9i32..=9)
    }

    fn big(c: [i32; 4]) -> [BigInt; 4] {
        c.map(BigInt::from)
    }
    fn as_f64(c: [i32; 4]) -> [f64; 4] {
        c.map(f64::from)
    }
    fn normals(p: [[i32; 4]; 3]) -> [[f64; 3]; 3] {
        p.map(|c| [f64::from(c[0]), f64::from(c[1]), f64::from(c[2])])
    }
    fn rows(p: &[[BigInt; 4]; 3]) -> [&[BigInt; 4]; 3] {
        [&p[0], &p[1], &p[2]]
    }

    proptest! {
        #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
            proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
        ))]
        #[test]
        fn the_integer_predicates_agree_with_the_f64_twins(
            p in proptest::array::uniform3(coeffs()),
            j in coeffs(),
        ) {
            let bp = p.map(big);
            prop_assert_eq!(
                int_plane_side(rows(&bp), &big(j)),
                indirect_plane_side(&ThreePlane(p.map(as_f64)), as_f64(j))
            );
            prop_assert_eq!(int_dir_sign(rows(&bp)), det3_sign(normals(p)));
        }

        #[test]
        fn cmp_coord_agrees_with_the_f64_twin(
            a in proptest::array::uniform3(coeffs()),
            b in proptest::array::uniform3(coeffs()),
            axis in 0usize..3,
        ) {
            let (ba, bb) = (a.map(big), b.map(big));
            // Both twins' precondition: each triple meets in a point.
            prop_assume!(int_dir_sign(rows(&ba)) != 0 && int_dir_sign(rows(&bb)) != 0);
            prop_assert_eq!(
                int_cmp_coord(rows(&ba), rows(&bb), axis),
                indirect_cmp_coord(&ThreePlane(a.map(as_f64)), &ThreePlane(b.map(as_f64)), axis)
            );
        }

        /// Rows scaled by ~200-bit positive integers — genuinely wider than `i128`, with
        /// the unscaled f64 twin as exact ground truth.
        #[test]
        fn a_wide_positive_scale_cannot_move_any_answer(
            p in proptest::array::uniform3(coeffs()),
            j in coeffs(),
            b in proptest::array::uniform3(coeffs()),
            axis in 0usize..3,
        ) {
            // A distinct scale per row, so no cross-row cancellation can hide a bug.
            let scale = |c: [i32; 4], k: u32| -> [BigInt; 4] {
                let s = (BigInt::from(1) << 199) + BigInt::from(12345 + k);
                c.map(|x| BigInt::from(x) * &s)
            };
            let sp: [[BigInt; 4]; 3] = core::array::from_fn(|k| scale(p[k], k as u32));
            let sj = scale(j, 7);
            prop_assert_eq!(
                int_plane_side(rows(&sp), &sj),
                indirect_plane_side(&ThreePlane(p.map(as_f64)), as_f64(j))
            );
            prop_assert_eq!(int_dir_sign(rows(&sp)), det3_sign(normals(p)));
            let sb: [[BigInt; 4]; 3] = core::array::from_fn(|k| scale(b[k], 11 + k as u32));
            if int_dir_sign(rows(&sp)) != 0 && int_dir_sign(rows(&sb)) != 0 {
                prop_assert_eq!(
                    int_cmp_coord(rows(&sp), rows(&sb), axis),
                    indirect_cmp_coord(
                        &ThreePlane(p.map(as_f64)),
                        &ThreePlane(b.map(as_f64)),
                        axis
                    )
                );
            }
        }

        /// The direction negative-control: negating one row flips exactly the
        /// direction-sensitive answers and nothing else.
        #[test]
        fn a_negated_row_flips_exactly_the_direction_sensitive_answers(
            p in proptest::array::uniform3(coeffs()),
            j in coeffs(),
            k in 0usize..3,
            axis in 0usize..3,
        ) {
            let bp = p.map(big);
            let bj = big(j);
            let neg = |c: &[BigInt; 4]| -> [BigInt; 4] { core::array::from_fn(|i| -&c[i]) };
            let mut np = bp.clone();
            np[k] = neg(&bp[k]);
            let nj = neg(&bj);

            // `plane_side`: a negated `p` row negates `D` and the dot together — invariant;
            // a negated `j` negates the dot alone — flipped.
            let side = int_plane_side(rows(&bp), &bj);
            prop_assert_eq!(int_plane_side(rows(&np), &bj), side);
            prop_assert_eq!(int_plane_side(rows(&bp), &nj), -side);

            // `dir_sign`: one negated row negates the determinant.
            prop_assert_eq!(int_dir_sign(rows(&np)), -int_dir_sign(rows(&bp)));

            // `cmp_coord`: orientation-invariant in every row.
            if int_dir_sign(rows(&bp)) != 0 && int_dir_sign(rows(&bj_triple(&bp, &bj))) != 0 {
                let other = bj_triple(&bp, &bj);
                prop_assert_eq!(
                    int_cmp_coord(rows(&np), rows(&other), axis),
                    int_cmp_coord(rows(&bp), rows(&other), axis)
                );
            }
        }
    }

    /// A second triple for `cmp_coord`, made from the first by swapping in `j` — cheap, and
    /// dependent enough on the inputs to exercise real coordinates.
    fn bj_triple(p: &[[BigInt; 4]; 3], j: &[BigInt; 4]) -> [[BigInt; 4]; 3] {
        [p[1].clone(), p[2].clone(), j.clone()]
    }
}

/// The projection lands **exactly on** the plane it came from, tilted ones included — checked
/// as `n·p + d == 0` in rationals, not within a tolerance.
/// A plane pushed along its own normal lands where the distance says — checked by taking the
/// projected origin of the moved plane and measuring it against the original's.
#[test]
fn a_plane_pushed_along_its_normal_moves_by_that_distance() {
    for (raw, name) in [
        ([0, 0, 1, -2], "z = 2"),
        ([0, 0, 10, -21], "z = 2.1"),
        ([3, 0, -4, -5], "3-4-5 tilt"),
        ([5, 12, 0, -13], "5-12-13"),
    ] {
        let c = ints(raw);
        let t = Rat::from_decimal(7.7).unwrap();
        let moved = plane_offset(c, t).unwrap_or_else(|| panic!("{name} has a rational normal"));
        // The two projected origins differ by exactly `t` along the shared unit normal.
        let (p0, p1) = (
            plane_origin_projection(c).unwrap(),
            plane_origin_projection(moved).unwrap(),
        );
        let d2 = (0..3).fold(Rat::from_int(0), |a, i| {
            let d = p1[i].checked_sub(p0[i]).unwrap();
            a.checked_add(d.checked_mul(d).unwrap()).unwrap()
        });
        assert_eq!(
            d2,
            t.checked_mul(t).unwrap(),
            "{name} moved the wrong distance"
        );
    }
}

/// ★★★ **Two steps compose into one, exactly.** This is the property the whole thing exists for:
/// a cap raised `7.7` and a cap raised `1.1` then `6.6` record the *same* plane, because
/// `11/10 + 66/10 = 77/10` holds in rationals where it does not in `f64`.
#[test]
fn pushing_twice_lands_where_pushing_once_does() {
    let d = |x: f64| Rat::from_decimal(x).unwrap();
    assert_ne!(1.1f64 + 6.6, 7.7, "the fixture must discriminate");
    for raw in [[0, 0, 1, -2], [3, 0, -4, -5], [0, 0, 10, -21]] {
        let c = ints(raw);
        let once = plane_offset(c, d(7.7)).unwrap();
        let twice = plane_offset(plane_offset(c, d(1.1)).unwrap(), d(6.6)).unwrap();
        assert_eq!(once, twice, "{raw:?}");
    }
}

/// Pushing by zero is the plane itself, and pushing back undoes it — the sign convention is the
/// one thing here a reader has to trust, so it is pinned in both directions.
#[test]
fn pushing_by_zero_and_pushing_back_are_identities() {
    // Every normal here has a rational length; `[7,-13,5]` does not, and `plane_offset`
    // declines it outright — that is the next test's business.
    for raw in [[0, 0, 1, -2], [3, 0, -4, -5], [5, 12, 0, 91]] {
        let c = canonical_plane_coeffs(ints(raw)).unwrap();
        assert_eq!(plane_offset(c, Rat::from_int(0)), Some(c), "{raw:?}");
        let t = Rat::new(-7, 2).unwrap();
        let there = plane_offset(c, t).unwrap();
        assert_eq!(
            plane_offset(there, Rat::from_int(0).checked_sub(t).unwrap()),
            Some(c),
            "{raw:?} did not come back"
        );
    }
}

/// **An irrational normal has no exact offset** — `|n|` is what must be rational, and `[1,1,1]`
/// is the smallest thing that is not. The caller keeps its f64 path; nothing is approximated.
#[test]
fn a_plane_whose_normal_has_no_rational_length_declines() {
    assert_eq!(plane_offset(ints([1, 1, 1, -7]), Rat::from_int(1)), None);
    assert_eq!(plane_offset(ints([1, 2, 3, 0]), Rat::from_int(1)), None);
    // A 3-4-5 direction does have one, so this is a statement about lengths, not about tilt.
    assert!(plane_offset(ints([3, 4, 0, -5]), Rat::from_int(1)).is_some());
    // Not a plane.
    assert_eq!(plane_offset(ints([0, 0, 0, 1]), Rat::from_int(1)), None);
}

#[test]
fn the_projected_origin_lies_exactly_on_its_plane() {
    for (raw, name) in [
        ([0, 0, 1, -2], "z = 2"),
        ([0, 0, 10, -21], "z = 2.1"),
        ([3, 0, -4, -5], "3-4-5 tilt"),
        ([1, 1, 1, -7], "diagonal"),
        ([5, 12, 0, -13], "5-12-13"),
        ([7, -13, 5, 91], "ugly"),
        ([0, 0, 1, 0], "through the origin"),
    ] {
        let c = ints(raw);
        let p = plane_origin_projection(c).unwrap_or_else(|| panic!("{name} has a projection"));
        let on = (0..3).fold(c[3], |acc, i| {
            acc.checked_add(c[i].checked_mul(p[i]).unwrap()).unwrap()
        });
        assert_eq!(
            on,
            Rat::from_int(0),
            "{name}: p = {p:?} is off its own plane"
        );
    }
}

/// ★★★ **One plane, three spellings, one point.** The design rests on this: the canonical form
/// carries no direction (`push_surface_with_coeffs` returns a `flipped` flag for exactly that
/// reason), and a plane is scale-invariant, so an origin derived from the coefficients would be
/// worthless if it moved when the coefficients were negated or scaled.
#[test]
fn the_projection_does_not_depend_on_how_the_plane_is_spelled() {
    let want = plane_origin_projection(ints([0, 0, 1, -3])).expect("a plane");
    for raw in [[0, 0, -1, 3], [0, 0, 2, -6], [0, 0, -5, 15]] {
        assert_eq!(plane_origin_projection(ints(raw)), Some(want), "{raw:?}");
    }
    // And with denominators: 11/10·x − 7/2 = 0 is 11x − 35 = 0, both giving p = (35/11, 0, 0).
    let fracs = [
        Rat::new(11, 10).unwrap(),
        Rat::from_int(0),
        Rat::from_int(0),
        Rat::new(-7, 2).unwrap(),
    ];
    assert_eq!(
        plane_origin_projection(fracs),
        plane_origin_projection(ints([11, 0, 0, -35]))
    );
}

/// The projection is the **nearest** point of the plane to the origin, which is what makes it a
/// sensible frame origin rather than merely a reproducible one: `p` is parallel to `n`, so no
/// other point of the plane is closer.
#[test]
fn the_projected_origin_is_the_nearest_point_of_its_plane() {
    let c = ints([1, 1, 1, -7]);
    let p = plane_origin_projection(c).expect("a plane");
    let d2 = |q: [Rat; 3]| {
        (0..3).fold(Rat::from_int(0), |a, i| {
            a.checked_add(q[i].checked_mul(q[i]).unwrap()).unwrap()
        })
    };
    // Step along an in-plane direction (n × ê is perpendicular to n) and the distance grows.
    for step in [Rat::from_int(1), Rat::new(-3, 7).unwrap()] {
        let dir = [Rat::from_int(1), Rat::from_int(-1), Rat::from_int(0)]; // ⊥ to (1,1,1)
        let q = [0, 1, 2].map(|i| p[i].checked_add(dir[i].checked_mul(step).unwrap()).unwrap());
        assert!(
            d2(q) > d2(p),
            "stepping by {step:?} did not move away from the origin"
        );
    }
}

/// Coefficients that are not a plane have no projection, and neither does an `i128` overflow —
/// both are the kernel's ordinary demotion, not a failure.
#[test]
fn a_non_plane_and_an_overflow_both_decline() {
    assert_eq!(plane_origin_projection(ints([0, 0, 0, 5])), None);
    assert_eq!(plane_origin_projection(ints([0, 0, 0, 0])), None);
    let huge = i128::MAX / 3;
    assert_eq!(plane_origin_projection(ints([huge, huge, huge, -1])), None);
}

#[test]
fn one_plane_written_at_any_scale_canonicalizes_to_one_vector() {
    let want = ints([1, 2, 0, -3]);
    for scale in [1, 2, 7, -1, -13] {
        let scaled = ints([scale, 2 * scale, 0, -3 * scale]);
        assert_eq!(canonical_plane_coeffs(scaled), Some(want), "scale {scale}");
    }
    // And with denominators: 11/10·x + 3/5·y − 7/2 = 0 is 11x + 6y − 35 = 0.
    let fracs = [
        Rat::new(11, 10).unwrap(),
        Rat::new(3, 5).unwrap(),
        Rat::from_int(0),
        Rat::new(-7, 2).unwrap(),
    ];
    assert_eq!(
        canonical_plane_coeffs(fracs),
        Some(ints([11, 6, 0, -35])),
        "denominators are cleared and the content divided out"
    );
}

#[test]
fn the_sign_convention_picks_the_first_nonzero_component() {
    assert_eq!(
        canonical_plane_coeffs(ints([0, -2, 4, 6])),
        Some(ints([0, 1, -2, -3])),
        "a plane and its negation are one plane, so one of the two spellings wins"
    );
}

/// ★★★ **The rule this function cannot enforce, made executable.**
///
/// Two faces of the plane `x = 3` reach `nacre_geom::Plane::coefficients()` as these two f64
/// vectors, because the un-normalized normal scales with the face's size and `d` is a rounded
/// product. Lifting them is lossless — and still gives two different planes, because the
/// *values* differ: `2.2 × 3` is not `6.6000000000000005`. Canonicalization cannot undo a
/// rounding that already happened, so the coefficients have to be rational from construction.
#[test]
fn lifting_rounded_f64_coefficients_does_not_merge_them() {
    let lift = |v: [f64; 4]| v.map(|x| Rat::try_from_f64(x).expect("finite"));
    let a = canonical_plane_coeffs(lift([2.2, 0.0, 0.0, -6.6000000000000005])).unwrap();
    let b = canonical_plane_coeffs(lift([13.2, 0.0, 0.0, -39.599999999999994])).unwrap();
    assert_ne!(
        a, b,
        "lifting a rounded coefficient carries the rounding in"
    );

    // Built from what the user *wrote*, and multiplied **in the rationals**, the same two
    // faces agree exactly.
    //
    // ★★ The `d` in `−3·k` has to be a `Rat` product. Writing `d(-3.0 * k)` instead puts the
    // multiplication back in f64, `−3.0 × 2.2` comes out `−6.6000000000000005`, and its
    // shortest decimal is that — not `−6.6`. A value the caller computed is not a value the
    // caller wrote, and `from_decimal` cannot tell them apart.
    let d = |x: f64| Rat::from_decimal(x).expect("decimal");
    let by_hand = |k: f64| {
        let k = d(k);
        [
            k,
            Rat::from_int(0),
            Rat::from_int(0),
            k.checked_mul(Rat::from_int(-3)).expect("small"),
        ]
    };
    assert_eq!(
        canonical_plane_coeffs(by_hand(2.2)),
        canonical_plane_coeffs(by_hand(13.2)),
        "rational construction is scale-independent"
    );
}

fn on_plane(c: [Rat; 4], p: [Rat; 3]) -> bool {
    let mut acc = c[3];
    for i in 0..3 {
        acc = acc.checked_add(c[i].checked_mul(p[i]).unwrap()).unwrap();
    }
    acc == Rat::from_int(0)
}

/// ★ **The moved plane must contain the moved points** — checked against a *separate*
/// point transform, so this cannot pass by restating the plane formula.
#[test]
fn a_moved_plane_still_contains_its_moved_points() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // A slanted plane, 2x − 3y + z − 6 = 0, and three points on it.
    let plane = ints([2, -3, 1, -6]);
    let pts = [
        [r(3, 1), r(0, 1), r(0, 1)],
        [r(0, 1), r(0, 1), r(6, 1)],
        [r(1, 2), r(1, 1), r(8, 1)],
    ];
    for p in pts {
        assert!(on_plane(plane, p), "fixture point is on the plane");
    }

    let cases = [
        Isometry::translation([r(1, 2), r(-3, 1), r(7, 5)]),
        Isometry::rotation(Rotation {
            axis: Axis::X,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
        }),
        Isometry::rigid(
            Rotation {
                axis: Axis::Z,
                pivot: [r(1, 1), r(2, 1), r(0, 1)],
                angle: Angle::from_deg(Rat::from_int(270)).unwrap(),
            },
            [r(0, 1), r(5, 1), r(-1, 2)],
        ),
    ];
    for iso in cases {
        let moved = iso.plane_coeffs(plane).expect("exact motion");
        for p in pts {
            assert!(
                on_plane(moved, iso.point_rat(p).expect("exact motion")),
                "moved plane {moved:?} must contain the moved point"
            );
        }
    }
}

#[test]
fn a_mirrored_plane_still_contains_its_mirrored_points() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let plane = ints([2, -3, 1, -6]);
    let (axis, offset) = (Axis::Z, r(5, 4));
    let moved = mirror_plane_coeffs(plane, axis, offset).unwrap();
    for p in [[r(3, 1), r(0, 1), r(0, 1)], [r(0, 1), r(0, 1), r(6, 1)]] {
        let q = mirror_point_rat(p, axis, offset).unwrap();
        assert!(on_plane(moved, q), "mirrored plane must contain {q:?}");
    }
}

/// A non-90° rotation has no rational `cos`/`sin`, so there is nothing exact to return —
/// ★ **for the point as well as the plane**. A caller that could move one description and
/// silently drop the other would leave a plane whose stored points are no longer on it.
#[test]
fn an_inexact_rotation_declines() {
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    assert!(!iso.is_exact());
    assert_eq!(iso.plane_coeffs(ints([0, 0, 1, -3])), None);
    assert_eq!(iso.point_rat([Rat::from_int(1); 3]), None);
}

#[test]
fn overflow_demotes_instead_of_panicking() {
    // Four coprime denominators near 10^10: their lcm is their product, ~10^40, past i128.
    let r = |den| Rat::new(1, den).unwrap();
    let coeffs = [
        r(10_000_000_019),
        r(10_000_000_033),
        r(10_000_000_061),
        r(10_000_000_069),
    ];
    assert_eq!(canonical_plane_coeffs(coeffs), None);
}
use proptest::prelude::*;

/// `try_from_f64` is the exact rational of the f64: it round-trips (`to_f64` gives
/// back the same bits) for finite values in range, handles 0/integers/dyadic
/// fractions exactly, and returns `None` for non-finite or overflowing inputs.
#[test]
fn try_from_f64_is_exact_and_round_trips() {
    assert_eq!(Rat::try_from_f64(0.0), Some(Rat::from_int(0)));
    assert_eq!(Rat::try_from_f64(6.0), Some(Rat::from_int(6)));
    assert_eq!(Rat::try_from_f64(-2.0), Some(Rat::from_int(-2)));
    assert_eq!(Rat::try_from_f64(0.5), Rat::new(1, 2));
    assert_eq!(Rat::try_from_f64(-0.75), Rat::new(-3, 4));
    assert_eq!(Rat::try_from_f64(f64::NAN), None);
    assert_eq!(Rat::try_from_f64(f64::INFINITY), None);
    assert_eq!(Rat::try_from_f64(1e300), None); // exponent overflows i128
    // round-trip over a spread of normal-range values (incl. non-dyadic f64s).
    for &x in &[0.1, 1.0 / 3.0, 2.0, 1000.0, -6.1, 4_503.7, 1e-6, 1e6] {
        assert_eq!(Rat::try_from_f64(x).unwrap().to_f64(), x, "round-trip {x}");
    }
}

/// **The property `from_decimal` rests on**: the shortest decimal of an f64,
/// read as an exact rational and realized again, gives back the same bits.
///
/// It holds *because* `to_f64` is correctly rounded, and only because of that. The
/// decimal is by construction a value whose nearest f64 is `x`; nearest rounding
/// therefore has no choice. Rounding numerator and denominator separately first —
/// what this function used to do — fails **17.4%** of the values below (measured),
/// since a 17-digit decimal has a numerator past 2⁵³. That is the whole reason this
/// cell touches `to_f64` at all.
#[test]
fn a_shortest_decimal_realizes_back_to_its_own_f64() {
    let mut state = 0x243f_6a88_85a3_08d3u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let (mut tried, mut declined) = (0, 0);
    for _ in 0..200_000 {
        // A random sign and mantissa over a CAD-plausible exponent range.
        let bits = next();
        let exp = 1023 - 60 + bits % 121;
        let x = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (exp << 52));
        match Rat::from_decimal(x) {
            Some(r) => {
                tried += 1;
                assert_eq!(r.to_f64(), x, "{x:?} → {r:?} → {:?}", r.to_f64());
            }
            None => declined += 1,
        }
    }
    // Nothing declines over this exponent range; the guard is here so a corpus that
    // drifted out of i128 could not turn this test vacuous.
    assert!(tried > 100_000, "{tried} tried, {declined} declined");
}

/// The same property, but over **every** finite f64 rather than the CAD-plausible
/// band: subnormals, `MIN_POSITIVE`, `MAX`, and the whole exponent range in between.
/// Where the power of ten leaves i128 the answer is `None` and the caller keeps its
/// f64 — what must never happen is a `Some` that realizes back to a *different*
/// number, because that would move a coordinate this cell has no business moving.
#[test]
fn from_decimal_never_lies_anywhere_in_the_finite_range() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut corpus = vec![
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::MIN,
        f64::from_bits(1), // the smallest subnormal
        1.1,
        7.7,
        1e-30,
        1e300,
    ];
    for _ in 0..200_000 {
        let bits = next();
        // Every exponent a finite f64 can have, subnormals included.
        let exp = (bits % 2047) << 52;
        corpus.push(f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | exp));
    }
    let mut declined = 0;
    for x in corpus {
        match Rat::from_decimal(x) {
            Some(r) => assert_eq!(r.to_f64(), x, "{x:e} → {r:?}"),
            None => declined += 1,
        }
    }
    assert!(
        declined > 0,
        "the i128 limit should bite somewhere out here"
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(Rat::from_decimal(bad), None, "{bad}");
    }
}

/// Where the i128 limit actually bites, pinned so the doc comment stays true and a
/// caller can tell whether its dimensions are anywhere near it. They are not: a CAD
/// model in millimetres sits around `1e0`, sixty orders of magnitude inside.
///
/// The two edges differ because a 17-digit mantissa spends 16 of its own powers of
/// ten on the fraction, so it runs out on the small side first. Note that exponent
/// *notation* is no obstacle in itself — `1e-30` is accepted; only the magnitude is.
#[test]
fn the_decimal_window_is_wide_and_its_edges_are_where_they_should_be() {
    let at = |s: &str| Rat::from_decimal(s.parse::<f64>().unwrap()).is_some();
    for e in -22..=38 {
        assert!(at(&format!("1.2345678901234567e{e}")), "17 digits at 1e{e}");
    }
    assert!(!at("1.2345678901234567e-23"), "17 digits at 1e-23");
    assert!(!at("1.2345678901234567e39"), "17 digits at 1e39");
    assert!(
        at("1e-38") && at("1e-30") && !at("1e-39"),
        "a short decimal reaches further down"
    );
}

/// The ends of the `i128` range, where the two operands' bounds differ.
///
/// `i128::MIN.unsigned_abs()` is `2¹²⁷` **exactly**, one past what the denominator
/// may be — an asymmetry easy to assert away and, when it was, a debug-only panic
/// on a value the algorithm handles correctly. The answers here are exact powers of
/// two and their neighbours, so they can be written down rather than approximated.
#[test]
fn the_ends_of_the_range_realize_exactly() {
    assert_eq!(Rat::from_int(i128::MIN).to_f64(), -(2f64.powi(127)));
    assert_eq!(Rat::from_int(i128::MAX).to_f64(), 2f64.powi(127)); // rounds up to 2¹²⁷
    assert_eq!(Rat::new(i128::MIN, 2).unwrap().to_f64(), -(2f64.powi(126)));
    assert_eq!(Rat::new(i128::MIN, i128::MAX).unwrap().to_f64(), -1.0);
    assert_eq!(Rat::new(1, i128::MAX).unwrap().to_f64(), 2f64.powi(-127));
    // Every one is finite and signed the way its numerator is.
    for (n, d) in [
        (i128::MIN, 3),
        (i128::MAX, 7),
        (-1, i128::MAX),
        (i128::MIN + 1, 1),
    ] {
        let q = Rat::new(n, d).unwrap().to_f64();
        assert!(q.is_finite() && (q < 0.0) == (n < 0), "{n}/{d} -> {q:e}");
    }
}

/// **The identity this whole cell exists for.** A dimension split into two and
/// stacked must land where the undivided one does. Lifting the f64 *values* cannot
/// give that — the drift is already inside them — and reading the decimals can.
#[test]
fn a_split_dimension_stacks_back_to_the_whole_one() {
    let split = Rat::from_decimal(1.1)
        .unwrap()
        .checked_add(Rat::from_decimal(6.6).unwrap())
        .unwrap();
    assert_eq!(split, Rat::new(77, 10).unwrap());
    assert_eq!(split, Rat::from_decimal(7.7).unwrap());
    assert_eq!(split.to_f64(), 7.7);
    // The f64 arithmetic this replaces, and the exact-but-binary lift that does not
    // help either — both land an ulp away.
    assert_ne!(1.1 + 6.6, 7.7);
    assert_ne!(
        Rat::try_from_f64(1.1)
            .unwrap()
            .checked_add(Rat::try_from_f64(6.6).unwrap())
            .unwrap(),
        Rat::try_from_f64(7.7).unwrap()
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Nearest-ness, checked against the definition rather than against a second
    /// implementation of it: `q` is the nearest f64 to `n/d` exactly when no
    /// neighbour of `q` is closer, and the comparison `|n/d − a/b| ≤ |n/d − c/e|`
    /// is decidable in exact rationals. Held to the range where those stay in i128.
    #[test]
    fn to_f64_lands_on_the_nearest_f64(n in -(1i128 << 60)..(1i128 << 60), d in 1i128..(1i128 << 60)) {
        prop_assume!(n != 0);
        let r = Rat::new(n, d).unwrap();
        let q = r.to_f64();
        let zero = Rat::from_int(0);
        let err = |y: f64| {
            let e = Rat::try_from_f64(y).and_then(|yr| r.checked_sub(yr))?;
            if e < zero { zero.checked_sub(e) } else { Some(e) }
        };
        // ★ The *check* has a range, and it is narrower than the generator: comparing exactly
        // needs `r`'s denominator times `q`'s power of two, which for a small enough value
        // leaves `i128` (measured: `n = 1, d = 500930446045` — 2^39 · 2^91). The two
        // neighbours below have always skipped on the same limit; `here` said `expect` and
        // so turned a limit of the instrument into a failure of the thing measured.
        prop_assume!(err(q).is_some());
        let here = err(q).expect("just assumed representable");
        for nb in [f64::from_bits(q.to_bits() + 1), f64::from_bits(q.to_bits() - 1)] {
            if let Some(there) = err(nb) {
                prop_assert!(here <= there, "{r:?}: {q:?} is not nearest ({nb:?} is closer)");
            }
        }
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// ★★★★★ **The total predicate must answer what the checked one answered.**
    ///
    /// `parallel_rat` replaced a checked-`Rat` cross that declined on overflow, and the
    /// change is only sound if the two agree wherever the old one could speak at all. The
    /// old spelling is kept here as the oracle — an independent derivation, not a call back
    /// into the code under test — and compared on every input it can answer.
    ///
    /// ★ **The width is drawn, not fixed**, because both halves of the range have to be
    /// visited: at `bits ≈ 30` the checked cross answers everything (so the two must agree),
    /// past `bits ≈ 63` its products leave `i128` and it falls silent (so the total one is
    /// alone). Measured with a counter before this spelling settled — a fixed narrow range
    /// left the oracle answering 100% of the time (differential, no coverage of the new
    /// behaviour) and a fixed wide one left it silent 100% of the time (coverage, no
    /// differential). Either way the test would have been half a test.
    #[test]
    fn the_total_cross_answers_what_the_checked_one_did(
        m in 0u32..40,
        xs in prop::array::uniform6(-(1i128 << 20)..(1i128 << 20)),
        ds in prop::array::uniform6(1i128..(1i128 << 20)),
    ) {
        // ★ The denominator's width is what decides whether the checked route survives, and
        // it is **drawn** so both halves of the range are visited. The scale is a power of
        // **three**, not two: scaling numerator and denominator by the same power of two
        // reduces straight back out (measured — the first spelling did exactly that and left
        // every case narrow), so the fraction has to be widened where the gcd cannot undo it.
        let scale = 3i128.saturating_pow(m);
        let r = |i: usize| {
            Rat::new(xs[i], (2 * ds[i] + 1).saturating_mul(scale)).unwrap()
        };
        let (a, b) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)]);
        // The retired spelling, verbatim: checked `Rat`, `None` on overflow.
        let checked = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<bool> {
            let term = |i: usize, j: usize| -> Option<Rat> {
                x[i].checked_mul(y[j])?.checked_sub(x[j].checked_mul(y[i])?)
            };
            let zero = Rat::from_int(0);
            Some(term(1, 2)? == zero && term(2, 0)? == zero && term(0, 1)? == zero)
        };
        if let Some(want) = checked(&a, &b) {
            prop_assert_eq!(parallel_rat(&a, &b), want, "{:?} x {:?}", a, b);
        }
        // Parallelism is symmetric and scale-invariant — properties the checked route could
        // not always demonstrate, and the total one must. (The scaling is asserted only when
        // every component survived it: a `Rat` whose numerator is already near the ceiling
        // has no ×3, and asserting through an `unwrap_or` fallback would compare a *different
        // vector* — the shape of a test that measures nothing.)
        prop_assert_eq!(parallel_rat(&a, &b), parallel_rat(&b, &a));
        let three = Rat::from_int(3);
        let scaled: Option<Vec<Rat>> = a.iter().map(|c| c.checked_mul(three)).collect();
        if let Some(s) = scaled {
            let s = [s[0], s[1], s[2]];
            prop_assert_eq!(parallel_rat(&s, &b), parallel_rat(&a, &b), "scale-invariance");
        }
    }

    /// **The clearance predicate must not drop the point's scale** — the negative control
    /// for the one decision in this section that is not a matter of taste.
    ///
    /// `(n·p + d)² − r²|n|²` is *not* homogeneous in `p`: clearing `p`'s denominators
    /// multiplies `n·p` while leaving `d` where it was, which states a different
    /// proposition. The naive spelling is written out here and compared with the truth in
    /// exact rationals; on mixed denominators the two part company (measured: 1.3% of
    /// random cases, so this test bites without needing a hand-picked adversary).
    #[test]
    fn clearance_keeps_the_points_scale(
        ns in prop::array::uniform3(-10i128..10),
        ps in prop::array::uniform3(-1000i128..1000),
        pd in prop::array::uniform3(1i128..(1i128 << 24)),
        rn in 1i128..100,
        rd in 1i128..100_000,
    ) {
        // ★ The generator is aimed at the **straddle**: an integer normal through the
        // origin, a point a hair off the plane (denominators up to 2²⁴), and a radius of
        // comparable size — so the answer is genuinely `Negative` about as often as
        // `Positive`. A generator whose points sit far outside every radius would compare
        // two implementations that always say `Positive`, which is how the first spelling
        // of this test passed while the scale-dropping defect was installed (measured: the
        // probe went green here and was caught only by a fixture two crates away).
        let coeffs: [Rat; 4] = [
            Rat::from_int(ns[0]),
            Rat::from_int(ns[1]),
            Rat::from_int(ns[2]),
            Rat::from_int(0),
        ];
        prop_assume!(coeffs[..3].iter().any(|c| *c != Rat::from_int(0)));
        let p: [Rat; 3] = core::array::from_fn(|i| Rat::new(ps[i], pd[i]).unwrap());
        let r = Rat::new(rn, rd).unwrap();

        // The truth, in exact rationals — no lifting, no scales, just the definition.
        let truth = (|| -> Option<Orient> {
            let dot = (0..3).try_fold(Rat::from_int(0), |a, i| {
                a.checked_add(coeffs[i].checked_mul(p[i])?)
            })?.checked_add(coeffs[3])?;
            let nn = (0..3).try_fold(Rat::from_int(0), |a, i| {
                a.checked_add(coeffs[i].checked_mul(coeffs[i])?)
            })?;
            let val = dot.checked_mul(dot)?
                .checked_sub(r.checked_mul(r)?.checked_mul(nn)?)?;
            Some(match val.cmp(&Rat::from_int(0)) {
                core::cmp::Ordering::Less => Orient::Negative,
                core::cmp::Ordering::Equal => Orient::Zero,
                core::cmp::Ordering::Greater => Orient::Positive,
            })
        })();
        if let Some(want) = truth {
            prop_assert_eq!(point_plane_clearance_rat(&coeffs, &p, &BigRat::square_of(r)), want,
                "coeffs {:?} p {:?} r {:?}", coeffs, p, r);
        }
    }

    /// The dot sign is the sign of the rational dot product, wherever that can be formed at
    /// all — and is total where it cannot.
    #[test]
    fn dot_sign_is_the_sign_of_the_dot(
        xs in prop::array::uniform6(-(1i128 << 30)..(1i128 << 30)),
        ds in prop::array::uniform6(1i128..(1i128 << 30)),
    ) {
        let r = |i: usize| Rat::new(xs[i], ds[i]).unwrap();
        let (a, b) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)]);
        let checked = (0..3).try_fold(Rat::from_int(0), |acc, i| {
            acc.checked_add(a[i].checked_mul(b[i])?)
        });
        if let Some(v) = checked {
            let want = match v.cmp(&Rat::from_int(0)) {
                core::cmp::Ordering::Less => Orient::Negative,
                core::cmp::Ordering::Equal => Orient::Zero,
                core::cmp::Ordering::Greater => Orient::Positive,
            };
            prop_assert_eq!(dot_sign_rat(&a, &b), want);
        }
        // Perpendicularity is symmetric, whatever the widths.
        prop_assert_eq!(dot_sign_rat(&a, &b), dot_sign_rat(&b, &a));
    }

    /// A vector is parallel to itself, to its multiples, and to zero — and **not** to a
    /// vector off its line. Without this the test above passes for a predicate that always
    /// answers `false` on the inputs the oracle cannot check.
    #[test]
    fn parallel_knows_a_line_from_a_plane(
        xs in prop::array::uniform3(-(1i128 << 40)..(1i128 << 40)),
        ds in prop::array::uniform3(1i128..(1i128 << 40)),
        k in 1i128..(1i128 << 20),
    ) {
        let v = [
            Rat::new(xs[0], ds[0]).unwrap(),
            Rat::new(xs[1], ds[1]).unwrap(),
            Rat::new(xs[2], ds[2]).unwrap(),
        ];
        let zero = [Rat::from_int(0); 3];
        prop_assert!(parallel_rat(&v, &v));
        prop_assert!(parallel_rat(&v, &zero), "zero is parallel to everything");
        let scaled: [Rat; 3] = core::array::from_fn(|i| {
            Rat::new(xs[i], ds[i]).unwrap().checked_mul(Rat::new(k, 1).unwrap()).unwrap_or(v[i])
        });
        prop_assert!(parallel_rat(&v, &scaled), "a multiple stays on the line");
        // A vector with one coordinate moved off the line is not parallel — unless `v` was
        // itself degenerate in that coordinate's plane, which the cross decides exactly.
        let mut off = v;
        off[0] = v[0].checked_add(Rat::from_int(1)).unwrap_or(v[0]);
        off[1] = v[1].checked_sub(Rat::from_int(1)).unwrap_or(v[1]);
        let cross_zero = {
            let t = |i: usize, j: usize| {
                v[i].checked_mul(off[j])
                    .and_then(|a| a.checked_sub(v[j].checked_mul(off[i])?))
            };
            match (t(1, 2), t(2, 0), t(0, 1)) {
                (Some(a), Some(b), Some(c)) => {
                    let z = Rat::from_int(0);
                    Some(a == z && b == z && c == z)
                }
                _ => None,
            }
        };
        if let Some(want) = cross_zero {
            prop_assert_eq!(parallel_rat(&v, &off), want);
        }
    }

    /// ★★★★★ **The two derivations must be the same function.**
    ///
    /// `plane_name_exact` runs the `Rat` route first and only falls back, so wherever the
    /// narrow one answers, the wide one is never consulted — and an error in it would sit
    /// there unseen until the day it *is* consulted, on inputs no test covers. Calling both on
    /// the same inputs is the only way to say they agree.
    ///
    /// ★ `plane_name_big` is `pub(crate)` for exactly this reason: a fallback hidden behind
    /// its filter cannot be tested against it.
    #[test]
    fn the_wide_derivation_answers_what_the_narrow_one_does(
        xs in prop::array::uniform9(-(1i64 << 20)..(1i64 << 20)),
        ds in prop::array::uniform9(1i64..(1i64 << 20)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let (a, b, c) = (
            [r(0), r(1), r(2)],
            [r(3), r(4), r(5)],
            [r(6), r(7), r(8)],
        );
        let narrow = plane_through_points(a, b, c);
        let wide = plane_name_big(a, b, c);
        if let Some(n) = narrow {
            // ★ The invariant rides along: an answer the narrow route reached fits `i128`
            // by construction, so the wide route must store it `Narrow` — same value, same
            // representation, structural equality.
            prop_assert_eq!(wide.clone(), Some(PlaneName::Narrow(n)),
                "narrow answered but wide disagrees");
        }
        // ★ And the wide one, whenever it answers narrowly at all, answers about a plane
        // these points are actually on — checked in the rationals, no tolerance.
        //
        // ★★★★★ **The residual can overflow even when the name and the points both fit**, and
        // this test found that by asserting it could not. `c · p` multiplies a canonical
        // coefficient by a point coordinate, so it needs the *sum* of their widths — which is
        // exactly the population `Model::push_surface_with_coeffs` cannot verify either. An
        // unevaluable check is not a failed one, here as there: skip it, never fail on it.
        if let Some(w) = wide.as_ref().and_then(|n| n.narrow()) {
            for p in [a, b, c] {
                let residual = (|| {
                    let mut acc = w[3];
                    for k in 0..3 {
                        acc = acc.checked_add(w[k].checked_mul(p[k])?)?;
                    }
                    Some(acc)
                })();
                if let Some(acc) = residual {
                    prop_assert_eq!(acc, Rat::from_int(0), "a point is off the derived plane");
                }
            }
        }
    }
}

/// ★★★★ **The width the narrow route cannot reach, and the wide one can.**
///
/// A triple whose reduced denominators are coprime — one a power of two, one a power of five,
/// which is what decimal arithmetic produces once it reduces — needs their product to state
/// the plane, and `(b − a) × (c − a)` needs it squared. `plane_through_points` gives up there;
/// the answer is small, and this is the case that says so.
#[test]
fn a_plane_the_narrow_route_gives_up_on_is_still_named() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let a = [r(1, 1 << 53), r(0, 1), r(0, 1)];
    let b = [r(0, 1), r(1, 5i128.pow(23)), r(0, 1)];
    let c = [r(0, 1), r(0, 1), r(1, (1 << 40) * 5i128.pow(11))];
    assert_eq!(
        plane_through_points(a, b, c),
        None,
        "the narrow route was expected to overflow on coprime denominators"
    );
    let name = plane_name_exact(a, b, c).expect("the wide route names it");
    // The canonical answer here is small (the doc above says so) — the invariant demands it
    // come back `Narrow`.
    let wide = name.narrow().expect("a small answer must be stored Narrow");
    for p in [a, b, c] {
        let mut acc = wide[3];
        for k in 0..3 {
            acc = acc
                .checked_add(wide[k].checked_mul(p[k]).expect("no overflow"))
                .expect("no overflow");
        }
        assert_eq!(acc, Rat::from_int(0), "a point is off the derived plane");
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// ★★★★★ **The two solves must be the same function** — the point-side twin of
    /// `the_wide_derivation_answers_what_the_narrow_one_does`, and for the same reason: a
    /// wide route that is only ever consulted where the narrow one declined would carry an
    /// error unseen until the day it is consulted. Calling both on the same inputs is the
    /// only way to say they agree.
    #[test]
    fn the_wide_solve_answers_what_the_narrow_one_does(
        xs in prop::array::uniform12(-(1i64 << 20)..(1i64 << 20)),
        ds in prop::array::uniform12(1i64..(1i64 << 20)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let rows: [[Rat; 4]; 3] =
            core::array::from_fn(|k| core::array::from_fn(|j| r(4 * k + j)));
        let names = rows.map(PlaneName::Narrow);
        let wide = three_planes_big([&names[0], &names[1], &names[2]]);

        if let Some(n) = three_planes_rat(rows) {
            // ★ The invariant rides along: a point the narrow route reached fits `Rat` by
            // construction, so the wide route must store it `Narrow` — same value, same
            // representation, structural equality.
            prop_assert_eq!(wide.clone(), Some(MeetPoint::Narrow(n)),
                "narrow answered but wide disagrees");
        }

        // ★★ And wherever the wide one answers **at all**, the point it names is on all three
        // planes — checked in the rationals, no tolerance. This half covers the `Wide` arm,
        // which the differential above cannot reach: the narrow route is silent there by
        // definition, so agreement says nothing and only the residual does.
        if let Some(w) = wide.as_ref() {
            use num_bigint::BigInt;
            use num_rational::Ratio;
            use num_traits::Zero;
            let big = |n: i128, d: i128| Ratio::new(BigInt::from(n), BigInt::from(d));
            let coord = |i: usize| -> Ratio<BigInt> {
                match w {
                    MeetPoint::Narrow(p) => big(p[i].numer(), p[i].denom()),
                    MeetPoint::Wide(p) => Ratio::new(p[i].0.clone(), p[i].1.clone()),
                }
            };
            let point = [coord(0), coord(1), coord(2)];
            for row in &rows {
                let mut acc = big(row[3].numer(), row[3].denom());
                for (j, x) in point.iter().enumerate() {
                    acc += big(row[j].numer(), row[j].denom()) * x;
                }
                prop_assert!(acc.is_zero(), "the solved point is off one of the planes");
            }
        }
    }
}

/// ★★★★★ **An overflowing intermediate no longer costs the answer** —
/// locked on the fixture that documents the overflow.
///
/// These three planes meet at `(1, 1, 1)`, which fits `Rat` with room to spare. But their
/// coefficients carry coprime denominators — a power of two and a power of five, what decimal
/// arithmetic produces once it reduces — and the determinant is their product: `2⁹³·5³⁴`, a
/// ~172-bit denominator. The narrow Cramer multiplies before it can reduce, so it overflows;
/// this fixture used to assert the resulting decline, and now asserts the fallback answers
/// through the integer core instead — the hand-known point, so this is an independent oracle
/// and not the two routes agreeing with each other.
///
/// ★ The private narrow route still declines here (asserted), so the fixture keeps proving
/// the fallback is *reached*, not merely present.
#[test]
fn an_overflowing_intermediate_no_longer_costs_the_answer() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let rows = [
        [r(1, 1 << 53), r(0, 1), r(0, 1), r(-1, 1 << 53)],
        [r(0, 1), r(1, 5i128.pow(23)), r(0, 1), r(-1, 5i128.pow(23))],
        [
            r(0, 1),
            r(0, 1),
            r(1, (1 << 40) * 5i128.pow(11)),
            r(-1, (1 << 40) * 5i128.pow(11)),
        ],
    ];
    assert_eq!(
        three_planes_rat_narrow(rows),
        None,
        "the narrow route was expected to overflow on coprime denominators — \
             without that this fixture no longer exercises the fallback"
    );
    assert_eq!(
        three_planes_rat(rows),
        Some([Rat::from_int(1); 3]),
        "the point fits `Rat`, and the solve now says so"
    );
    let names = rows.map(PlaneName::Narrow);
    let found = three_planes_big([&names[0], &names[1], &names[2]])
        .expect("the wide solve finds the point");
    assert_eq!(
        found,
        MeetPoint::Narrow([Rat::from_int(1); 3]),
        "the point fits `Rat` — the invariant demands it come back Narrow"
    );
    assert_eq!(found.width_bits(), 1, "1/1 is one bit wide");
}

/// The two remaining meanings of [`three_planes_rat`]'s `None`, each still honest:
/// no unique point (parallel planes), and a point that truly does not fit `Rat` — which
/// [`three_planes_big`] tells apart by answering `Wide`.
#[test]
fn the_solve_still_declines_what_it_should() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // Parallel pair: x = 0 and x = 1 — no unique point, both routes say so.
    let parallel = [
        [r(1, 1), r(0, 1), r(0, 1), r(0, 1)],
        [r(1, 1), r(0, 1), r(0, 1), r(-1, 1)],
        [r(0, 1), r(0, 1), r(1, 1), r(0, 1)],
    ];
    assert_eq!(
        three_planes_rat(parallel),
        None,
        "parallel planes meet nowhere"
    );
    // Narrow rows whose meeting point is wider than `Rat`: x + y = 2⁻¹⁰⁰, x − y = 5⁻⁵⁰
    // put a denominator of 2¹⁰¹·5⁵⁰ (~217 bits) on x. `three_planes_rat` declines;
    // the wide twin answers, and answers `Wide` — the causes stay told apart.
    let wide_point = [
        [r(1, 1), r(1, 1), r(0, 1), r(-1, 1 << 100)],
        [r(1, 1), r(-1, 1), r(0, 1), r(-1, 5i128.pow(50))],
        [r(0, 1), r(0, 1), r(1, 1), r(0, 1)],
    ];
    assert_eq!(
        three_planes_rat(wide_point),
        None,
        "a point no `Rat` can hold is a decline, not an answer"
    );
    let names = wide_point.map(PlaneName::Narrow);
    assert!(
        matches!(
            three_planes_big([&names[0], &names[1], &names[2]]),
            Some(MeetPoint::Wide(_))
        ),
        "the twin names the cause: the point exists and is wide"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Totality: wherever the integer core answers `Narrow`, [`three_planes_rat`] answers
    /// the same. ★ A **wiring** lock, stated as such: after the refactor both routes
    /// converge on `three_planes_int`, so this pins the `.or_else` plumbing (lift,
    /// `narrow()` return) rather than serving as an independent oracle — that role belongs
    /// to `an_overflowing_intermediate_no_longer_costs_the_answer`'s hand-known point and
    /// to the agreement direction the ops-side `point_width` invariant keeps.
    #[test]
    fn the_solve_answers_wherever_the_answer_is_narrow(
        rows in proptest::array::uniform3(proptest::array::uniform4((-9i128..=9, 1u32..=60))),
    ) {
        let rat_rows: [[Rat; 4]; 3] =
            rows.map(|row| row.map(|(n, e)| Rat::new(n, 1i128 << e).unwrap()));
        let names = rat_rows.map(PlaneName::Narrow);
        let big = three_planes_big([&names[0], &names[1], &names[2]]);
        let expect = match &big {
            Some(MeetPoint::Narrow(p)) => Some(*p),
            _ => None, // no unique point, or truly wide — the honest declines
        };
        prop_assert_eq!(three_planes_rat(rat_rows), expect);
    }
}

/// ★★★★★ **The negative control for [`MeetPoint::width_bits`]** — without it, a corpus that
/// reports "nothing over 127 bits" is indistinguishable from a dead probe.
///
/// A `Wide` carrier states a plane at `x = 2²⁰⁰`, which no `Rat` can hold; the meeting point
/// inherits that width and the fork has to take its other branch. So the instrument can say
/// the other thing, and `over127 = 0` in a measurement means the population, not the meter.
#[test]
fn the_width_meter_reports_a_point_no_rat_can_hold() {
    use num_bigint::BigInt;
    let far = BigInt::from(1) << 200;
    let wide = PlaneName::Wide([BigInt::from(1), BigInt::from(0), BigInt::from(0), -&far]);
    let r = |v: [i128; 4]| PlaneName::Narrow(v.map(Rat::from_int));
    let (py, pz) = (r([0, 1, 0, 0]), r([0, 0, 1, 0]));
    let found = three_planes_big([&wide, &py, &pz]).expect("three planes meet at (2²⁰⁰, 0, 0)");
    assert!(
        matches!(found, MeetPoint::Wide(_)),
        "a 201-bit coordinate cannot be Narrow"
    );
    assert_eq!(
        found.width_bits(),
        201,
        "the meter reads the coordinate's width"
    );
    assert_eq!(found.narrow(), None);
}

/// ★★★★★ **A canonical answer wider than `i128` is a name now, not a `None`** — and two
/// statements of that plane are one value. The
/// interning consequence is locked on the model side
/// (`a_wide_plane_interns_but_opens_no_shortcut` in nacre-topo).
#[test]
fn a_plane_too_wide_for_i128_is_named_wide() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // Cross-product terms multiply two coordinates' numerators, so two ~2^90 coprime
    // numerators push the canonical coefficients past i128 with no content to divide out.
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let a = [r(big1, 3), r(big2, 7), r(0, 1)];
    let b = [r(-big2, 5), r(big1, 11), r(0, 1)];
    let c = [r(1, 13), r(1, 17), r(1, 19)];
    assert_eq!(
        plane_through_points(a, b, c),
        None,
        "expected the narrow route to overflow"
    );
    let name = plane_name_exact(a, b, c).expect("collinear it is not — it must be named");
    // ★ Wide carries identity only: the arithmetic shortcuts' door stays shut
    // (a wide name must NOT open a frame).
    assert!(
        name.narrow().is_none(),
        "an answer past i128 must be stored Wide"
    );
    // ★ Same plane, two spellings (a permuted triple) — one value, structurally.
    let permuted = plane_name_exact(b, c, a).expect("the same plane, permuted");
    assert_eq!(
        name, permuted,
        "two statements of one wide plane must be one value"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// **The two orientation arms must be the same function** — same argument as
    /// `the_wide_derivation_answers_what_the_narrow_one_does`: `orient2d_rat` runs the `Rat`
    /// route first, so wherever it answers the `BigInt` arm is never consulted, and an error
    /// there would wait for the first overflowing input. `orient2d_big` is `pub(crate)` for
    /// exactly this call.
    #[test]
    fn the_big_orientation_answers_what_the_narrow_one_does(
        xs in prop::array::uniform6(-(1i64 << 20)..(1i64 << 20)),
        ds in prop::array::uniform6(1i64..(1i64 << 20)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let (a, b, c) = ([r(0), r(1)], [r(2), r(3)], [r(4), r(5)]);
        // Small operands: the narrow route always answers, so this compares the arms.
        prop_assert_eq!(orient2d_rat(a, b, c), orient2d_big(a, b, c));
    }
}

/// ★★ **The population the narrow route cannot reach — and the reason the sign is total.**
///
/// Two decimal-window denominators of `10^22` put the determinant's product denominator at
/// `10^44`, past `i128` — precisely what `from_decimal`'d profile coordinates produce. The
/// answers are hand-computable: `(0,0) → (1/d,0) → (0,1/d)` turns counter-clockwise, its
/// mirror clockwise, and a doubled point on the same ray is collinear.
#[test]
fn an_orientation_past_i128_still_gets_its_sign() {
    let d = 10i128.pow(22);
    let r = |n: i128, den: i128| Rat::new(n, den).unwrap();
    let (o, x, y) = (
        [r(0, 1), r(0, 1)],
        [r(1, d), r(0, 1)],
        [r(0, 1), r(1, d + 1)],
    );
    // The narrow route must actually be dead here, or this pins nothing (self-qualification).
    assert!(
        r(1, d).checked_mul(r(1, d + 1)).is_none(),
        "the product denominator was expected to overflow i128"
    );
    assert_eq!(orient2d_rat(o, x, y), 1, "counter-clockwise");
    assert_eq!(orient2d_rat(o, y, x), -1, "clockwise");
    let far = [r(2, d), r(2, d)];
    let near = [r(1, d + 1), r(1, d + 1)];
    assert_eq!(orient2d_rat(o, near, far), 0, "one ray, three points");
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// **The residual's two arms must be the same function** — `plane_residual_sign` runs the
    /// `Rat` substitution first, so wherever it answers the `BigInt` arm is never consulted
    /// (the `orient2d` arrangement, one dimension up).
    #[test]
    fn the_big_residual_answers_what_the_narrow_one_does(
        // Small enough that the canonical name always fits i128 (the lifted integers stay
        // near 2^40 and the derivation's peak near 2^122), so `narrow()` below cannot shrug.
        xs in prop::array::uniform9(-(1i64 << 10)..(1i64 << 10)),
        ds in prop::array::uniform9(1i64..(1i64 << 10)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let (a, b, c) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)], [r(6), r(7), r(8)]);
        if let Some(name) = plane_name_exact(a, b, c) {
            let coeffs = name.narrow().expect("small operands stay narrow");
            let ci = coeffs.map(|x| num_bigint::BigInt::from(x.numer()));
            // The naming points themselves: both arms must call them on-plane...
            for p in [a, b, c] {
                prop_assert_eq!(plane_residual_sign(&name, p), 0);
                prop_assert_eq!(residual_sign_big(&ci, p), 0);
            }
            // ...and an off-plane probe (a naming point pushed along the normal) must get
            // the same nonzero sign from both.
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            if let Some(q) = (|| -> Option<[Rat; 3]> {
                Some([
                    a[0].checked_add(n[0])?,
                    a[1].checked_add(n[1])?,
                    a[2].checked_add(n[2])?,
                ])
            })() {
                let s = plane_residual_sign(&name, q);
                prop_assert_eq!(s, 1, "n·n > 0: the push is to the positive side");
                prop_assert_eq!(residual_sign_big(&ci, q), s);
            }
        }
    }
}

/// ★★ **A `Wide` name still judges its point** — the population `plane_residual_sign` exists
/// for. Same fixture as `a_plane_too_wide_for_i128_is_named_wide`: canonical coefficients
/// past `i128`, no `Rat` route to fall back on.
#[test]
fn a_residual_against_a_wide_name_still_gets_its_sign() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let a = [r(big1, 3), r(big2, 7), r(0, 1)];
    let b = [r(-big2, 5), r(big1, 11), r(0, 1)];
    let c = [r(1, 13), r(1, 17), r(1, 19)];
    let name = plane_name_exact(a, b, c).expect("a genuine plane");
    assert!(name.narrow().is_none(), "the fixture must actually be Wide");
    for p in [a, b, c] {
        assert_eq!(plane_residual_sign(&name, p), 0, "a naming point is on it");
    }
    // A point off the plane along ±z (the naming triangle is not vertical: a and b span
    // z = 0 and c leaves it, so the normal has a z-component): opposite pushes must get
    // opposite, nonzero signs.
    let one = Rat::from_int(1);
    let up = [c[0], c[1], c[2].checked_add(one).unwrap()];
    let down = [c[0], c[1], c[2].checked_sub(one).unwrap()];
    let (su, sd) = (
        plane_residual_sign(&name, up),
        plane_residual_sign(&name, down),
    );
    assert_ne!(su, 0);
    assert_eq!(su, -sd, "opposite sides, opposite signs");
}

/// ★ `three_planes_rat` — the exact corner of three rational planes. Hand-checkable
/// fixture, a degenerate (line-sharing) triple, and an i128-overflow decline.
#[test]
fn three_rational_planes_meet_where_they_should() {
    let r = Rat::from_int;
    // x = 2, y = 3, z = 5  (as a·x + d = 0 rows: [1,0,0,-2] etc.)
    let pt = three_planes_rat([
        [r(1), r(0), r(0), r(-2)],
        [r(0), r(1), r(0), r(-3)],
        [r(0), r(0), r(1), r(-5)],
    ])
    .expect("axis planes meet in one point");
    assert_eq!(pt, [r(2), r(3), r(5)]);
    // A tilted but exact triple: x+y=1, x−y=0, z=7 ⇒ (1/2, 1/2, 7).
    let pt = three_planes_rat([
        [r(1), r(1), r(0), r(-1)],
        [r(1), r(-1), r(0), r(0)],
        [r(0), r(0), r(1), r(-7)],
    ])
    .expect("a fair triple");
    assert_eq!(pt, [Rat::new(1, 2).unwrap(), Rat::new(1, 2).unwrap(), r(7)]);
    // Three planes through one line (z = 0, y = 0, y + z = 0): det = 0 — no unique point.
    assert_eq!(
        three_planes_rat([
            [r(0), r(0), r(1), r(0)],
            [r(0), r(1), r(0), r(0)],
            [r(0), r(1), r(1), r(0)],
        ]),
        None
    );
    // Coefficients near the i128 edge overflow the narrow Cramer — and the answer still
    // fits `Rat` (denominator 2¹²⁶ + 1), so the fallback answers it: declining here would be a defect.
    let big = Rat::from_int(1 << 126);
    let edge = [
        [big, r(1), r(0), r(-1)],
        [r(1), big, r(0), r(-1)],
        [r(0), r(0), big, big],
    ];
    assert_eq!(
        three_planes_rat_narrow(edge),
        None,
        "the narrow route overflows here"
    );
    let inv = Rat::new(1, (1 << 126) + 1).unwrap();
    assert_eq!(three_planes_rat(edge), Some([inv, inv, r(-1)]));
}

/// The core rational-representation property in miniature: exact rational accumulation does not
/// drift, where the f64 control does. `(1/10)` summed ten times is exactly
/// `1`, but `0.1_f64` summed ten times is not `1.0`.
#[test]
fn rational_accumulation_is_exact_where_f64_drifts() {
    let tenth = Rat::new(1, 10).unwrap();
    let mut acc = Rat::from_int(0);
    for _ in 0..10 {
        acc = acc.checked_add(tenth).unwrap();
    }
    assert_eq!(acc, Rat::from_int(1));

    let mut f = 0.0_f64;
    for _ in 0..10 {
        f += 0.1;
    }
    assert_ne!(f, 1.0); // 0.9999999999999999 — the drift Rat avoids
}

/// The "thin film" example: `1.1 × 7` must be exactly `7.7`. In rationals
/// `11/10 × 7 = 77/10`; in f64 `1.1 * 7.0` is not `7.7`.
#[test]
fn one_point_one_times_seven_is_exact() {
    let a = Rat::new(11, 10).unwrap();
    let seven = Rat::from_int(7);
    assert_eq!(a.checked_mul(seven).unwrap(), Rat::new(77, 10).unwrap());

    assert_ne!(1.1_f64 * 7.0, 7.7); // f64 cannot represent 7.7 exactly
}

/// Measurement: with fixed-width i128, chained coprime-denominator
/// accumulation *does* overflow (the finite-precision cliff handled by
/// downgrading). Summing `1/p` over successive primes forces the denominator
/// toward the primorial, which exceeds i128. This test pins two facts:
///   (a) a modest sum (first 8 primes) stays exact — normal use is fine;
///   (b) accumulation eventually overflows within the prime list — proving
///       the downgrade trigger fires, and reporting *where* (onset index).
#[test]
fn coprime_accumulation_overflows_and_reports_onset() {
    const PRIMES: [i128; 30] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113,
    ];

    // (a) first 8 primes stay exact.
    let mut acc = Rat::from_int(0);
    for &p in &PRIMES[..8] {
        acc = acc
            .checked_add(Rat::new(1, p).unwrap())
            .expect("first 8 primes must not overflow");
    }

    // (b) full accumulation eventually overflows; record the onset.
    let mut acc = Rat::from_int(0);
    let mut onset = None;
    let mut last_bits = 0;
    for (i, &p) in PRIMES.iter().enumerate() {
        match acc.checked_add(Rat::new(1, p).unwrap()) {
            Some(next) => {
                acc = next;
                last_bits = acc.bit_width();
            }
            None => {
                onset = Some((i, last_bits));
                break;
            }
        }
    }
    let (idx, bits) = onset.expect("i128 rational accumulation must overflow within 30 primes");
    // Measurement record (run with `-- --nocapture`).
    eprintln!(
        "[rational] coprime 1/p accumulation overflows at prime index {idx} (p={}); \
             denominator bit-width just before onset = {bits}",
        PRIMES[idx]
    );
    // Onset is comfortably past normal use and near the i128 ceiling (~127 bits).
    assert!(idx >= 8, "overflow onset {idx} should be past modest use");
    assert!(
        bits > 100,
        "denominator should be near the i128 ceiling at onset, got {bits} bits"
    );
}

/// Closure: a rational angle accumulates exactly, so a full turn lands back
/// on exactly `0` — while the f64 control drifts. `360/7` degrees added seven
/// times is exactly `360 → 0`; `360.0/7.0` summed seven times in f64 is not
/// `360.0`.
#[test]
fn rational_angle_closes_exactly_where_f64_drifts() {
    let step = Rat::new(360, 7).unwrap();
    let mut a = Angle::from_deg(Rat::from_int(0)).unwrap();
    for _ in 0..7 {
        a = a.checked_add(step).unwrap();
    }
    assert_eq!(a, Angle::from_deg(Rat::from_int(0)).unwrap());

    let mut f = 0.0_f64;
    for _ in 0..7 {
        f += 360.0 / 7.0;
    }
    assert_ne!(f, 360.0); // f64 drifts off the full turn
}

/// Many small rational steps also close: `1/3` degree added 1080 times is
/// exactly one full turn → `0`.
#[test]
fn many_small_steps_close_to_zero() {
    let step = Rat::new(1, 3).unwrap();
    let mut a = Angle::from_deg(Rat::from_int(0)).unwrap();
    for _ in 0..1080 {
        a = a.checked_add(step).unwrap();
    }
    assert_eq!(a, Angle::from_deg(Rat::from_int(0)).unwrap());
}

/// The angle value stays exact, and its cos/sin realization is f64 — but a *correctly rounded*
/// one, so the irrational-realization boundary now costs at most half an ulp rather than
/// whatever the platform's libm happened to do.
///
/// **`cos 45°` is the case to look at**: `√2/2` cannot be an f64, so the realization is
/// genuinely lossy, and it must land on the nearest f64 to the truth. Checked by asking a far
/// deeper realization whether anything is closer.
#[test]
fn realization_is_f64_while_angle_stays_exact() {
    let a = Angle::from_deg(Rat::from_int(45)).unwrap();
    assert_eq!(a.deg(), Rat::from_int(45)); // angle exact
    let (c, _) = a.cos_sin_f64();
    assert!(a.try_exact_cos_sin().is_none()); // √2/2 is not rational

    // Nothing is nearer: both neighbours are further from the deep truth than `c` is.
    let (deep, _) = a.cos_sin_at(512);
    let dist = |f: f64| BigFloat::from_f64(f, 512).sub(&deep, 512, HP_RM).abs();
    let here = dist(c);
    for nb in [
        f64::from_bits(c.to_bits() - 1),
        f64::from_bits(c.to_bits() + 1),
    ] {
        assert!(
            here.cmp(&dist(nb)).is_some_and(|s| s < 0),
            "a neighbour of {c:e} is nearer the truth"
        );
    }
}

/// A 90°-family angle yields exact rational cos/sin, so rotating a rational
/// point stays exact (tol 0): `(x, y)` rotated 90° is `(-y, x)`. A 45° angle
/// has no exact rational cos/sin (√2/2), so it returns `None` and would fall
/// to the f64/dd realization. Note the f64 path is *not* exact here:
/// `cos 90°` realizes to ~6e-17, not `0`.
#[test]
fn exact_cos_sin_only_for_quadrantal_angles() {
    let a90 = Angle::from_deg(Rat::from_int(90)).unwrap();
    let (c, s) = a90.try_exact_cos_sin().unwrap();
    assert_eq!((c, s), (Rat::from_int(0), Rat::from_int(1)));

    // Exact rotation of (3, 5) by 90° → (-5, 3), all rational (tol 0):
    // x' = x·cos − y·sin,  y' = x·sin + y·cos.
    let (x, y) = (Rat::from_int(3), Rat::from_int(5));
    let xr = x
        .checked_mul(c)
        .unwrap()
        .checked_sub(y.checked_mul(s).unwrap())
        .unwrap();
    let yr = x
        .checked_mul(s)
        .unwrap()
        .checked_add(y.checked_mul(c).unwrap())
        .unwrap();
    assert_eq!((xr, yr), (Rat::from_int(-5), Rat::from_int(3)));

    // 45° has no exact rational realization → None (falls to the rounded high-precision path).
    let a45 = Angle::from_deg(Rat::from_int(45)).unwrap();
    assert!(a45.try_exact_cos_sin().is_none());
    // The high-precision realization of 90° is *not* zero either — it is ~2⁻¹²⁸ — which is
    // exactly why the branch above exists rather than being an optimization.
    assert!(!a90.cos_sin_at(128).0.is_zero());
}

/// `cos_sin_f64` snaps the 90°-family to exact `0.0`/`±1.0`, and the exactly-representable
/// values Niven allows off that family come out exact too.
///
/// ★★ **`cos 60° == 0.5` is the visible proof that libm left.** `1/2` is one of the three
/// rational values a rational-degree cosine can take, and it *is* an f64 — but reaching it
/// through `(60.0 * PI / 180.0).cos()` does not land on it. Rounding the high-precision value
/// does.
#[test]
fn cos_sin_f64_is_exact_where_the_true_value_is_representable() {
    let deg = |d| Angle::from_deg(Rat::from_int(d)).unwrap();
    assert_eq!(deg(0).cos_sin_f64(), (1.0, 0.0));
    assert_eq!(deg(90).cos_sin_f64(), (0.0, 1.0));
    assert_eq!(deg(180).cos_sin_f64(), (-1.0, 0.0));
    assert_eq!(deg(270).cos_sin_f64(), (0.0, -1.0));
    for (d, want) in [(60, 0.5), (120, -0.5), (240, -0.5), (300, 0.5)] {
        assert_eq!(deg(d).cos_sin_f64().0, want, "cos {d}°");
    }
    for (d, want) in [(30, 0.5), (150, 0.5), (210, -0.5), (330, -0.5)] {
        assert_eq!(deg(d).cos_sin_f64().1, want, "sin {d}°");
    }
}

/// **The 90°-family realizes with no error at all, and that zero is load-bearing.**
///
/// `WitnessPoint::rotate_about` reads [`Angle::realization_error_of`] to decide whether a rotation
/// contributes any tolerance. For `cos`/`sin` in `{0, ±1}` the f64 values *are* the true ones,
/// and every product and difference downstream is exact too — which is why a quadrantal origin
/// rotation stays at tol 0 and an axis-aligned model never leaves the exact predicate path.
/// A nonzero here would not fail loudly; it would quietly move those models.
///
/// Everything else must report *something* — but only where the true value is *not*
/// representable. Since the realization became correctly rounded, `cos 60°` really is `0.5`
/// exactly, so the corpus below has to avoid the four angles where that happens or it would be
/// asserting a nonzero error that does not exist.
#[test]
fn only_the_quadrantal_family_realizes_exactly() {
    let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
    let err = |a: Angle| {
        let (c, s) = a.cos_sin_f64();
        a.realization_error_of(c, s)
    };
    for d in [0, 90, 180, 270] {
        assert_eq!(err(deg(d, 1)), (0.0, 0.0), "{d} deg");
    }
    // ★ And the zero is of the *pair*, not of the angle: hand a 90°-family angle some other
    // realization and it must be measured like anything else. Otherwise a caller that got its
    // cos/sin from somewhere other than `cos_sin_f64` would be handed an exactness claim that
    // does not hold of what it is holding.
    let a90 = deg(90, 1);
    assert!(a90.realization_error_of(1e-17, 1.0).0 > 0.0);
    for (n, d) in [(1, 1), (37, 1), (45, 1), (337, 1), (1, 3), (359999, 1000)] {
        let (dc, ds) = err(deg(n, d));
        assert!(
            dc > 0.0 && ds > 0.0,
            "{n}/{d} deg claimed an exact realization"
        );
        // ε-scale: a bound this large would mean the measurement, not the platform, is wrong.
        assert!(
            dc < 64.0 * f64::EPSILON && ds < 64.0 * f64::EPSILON,
            "{n}/{d} deg: {dc:e}"
        );
    }
}

/// **The f64 read-out is a re-encoding, not a computation — checked by round trip.**
///
/// `set_precision(53, ToEven)` is where a value loses bits; everything after it is supposed to
/// be pure bookkeeping over the mantissa words. So any `f64` put in must come back out
/// unchanged. **The 32-bit-`Word` path is the one this is really for** — the mantissa is
/// assembled across two words there, wasm is a 32-bit target, and no amount of reading the
/// crate source substitutes for running it on the target that ships.
#[test]
fn the_f64_readout_round_trips() {
    let mut cases = vec![
        1.0,
        0.5,
        -0.5,
        0.9999999999999999,
        1e-300,
        -3.7e17,
        f64::MIN_POSITIVE,
    ];
    let mut st = 0x2545_F491_4F6C_DD1Du64;
    for _ in 0..2000 {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        // Any finite normal double; the exponent is squeezed into the normal range.
        let bits = (st & !(0x7ffu64 << 52)) | ((1 + (st >> 53) % 2045) << 52);
        let x = f64::from_bits(bits);
        if x.is_normal() {
            cases.push(x);
        }
    }
    for x in cases {
        let back = to_f64_exact(&BigFloat::from_f64(x, 128));
        assert_eq!(back, Some(x), "{x:e} did not survive the round trip");
    }
    assert_eq!(to_f64_exact(&BigFloat::from_f64(0.0, 128)), Some(0.0));
}

/// **The rounding check has to be able to say no.**
///
/// [`round_to_f64`] answers `None` when the interval straddles a rounding boundary, and that
/// branch is the whole reason the function is not just "round the midpoint". A check that can
/// only say yes is indistinguishable from no check — and there is a specific way to build one
/// here, by forming `mid ± rad` at the realization's own precision so the radius rounds away.
/// So: a radius wide enough to be undecidable must be refused, and a tight one accepted.
#[test]
fn the_rounding_check_refuses_an_undecidable_interval() {
    let a = Angle::from_deg(Rat::new(37, 1).unwrap()).unwrap();
    let (c, _) = a.cos_sin_bounded(128);
    let (c, rc) = (c.value, c.error);
    assert!(
        round_to_f64(&c, rc, 128).is_some(),
        "a 2^-128 radius is decidable"
    );
    // An ulp-wide radius cannot be: it reaches both neighbours.
    assert!(round_to_f64(&c, Mag::pow2(-52), 128).is_none());
    // And zero is the case no precision resolves — cos 90 is exactly 0, so its interval
    // straddles zero forever. This is why `cos_sin_f64` resolves the family first.
    let a90 = Angle::from_deg(Rat::from_int(90)).unwrap();
    for prec in [128usize, 256, 512] {
        let (c90, _) = a90.cos_sin_bounded(prec);
        let (c90, r90) = (c90.value, c90.error);
        assert!(
            round_to_f64(&c90, r90, prec).is_none(),
            "cos 90 became decidable at {prec}, which would make the quadrantal branch optional"
        );
    }
}

/// The escalation and fallback rungs, counted — over a corpus that reaches for them.
///
/// Near an axis `|cos|` is tiny while its error bound is absolute, so the relative radius grows
/// and 128 bits stops being enough; that is the only place the second rung is reachable. **The
/// count is reported rather than asserted nonzero**: with `Rat` bounded by `i128` an angle
/// cannot get closer to 90° than ~6e-39 degrees, so it is entirely possible that nothing in a
/// finite corpus needs it — but a silent zero and an unreachable branch look identical, and
/// this at least says which corpus produced the zero.
#[test]
fn the_escalation_rung_is_reachable() {
    let before = ROUND_ESCALATED.with_borrow(|c| *c);
    let mut asked = 0usize;
    for k in 1..400i128 {
        // Just off 90 degrees, by ever smaller amounts.
        for d in [
            10i128.pow(9),
            10i128.pow(18),
            10i128.pow(30),
            i128::MAX / 91,
        ] {
            if let Some(a) = Rat::new(90 * d + k, d).and_then(Angle::from_deg) {
                let (c, s) = a.cos_sin_f64();
                assert!(
                    c.is_finite() && s.is_finite(),
                    "{a:?} realized to a non-number"
                );
                // Whatever rung answered, the answer must still bound its own error.
                let (dc, ds) = a.realization_error_of(c, s);
                assert!(dc >= 0.0 && ds >= 0.0);
                asked += 1;
            }
        }
    }
    let after = ROUND_ESCALATED.with_borrow(|c| *c);
    eprintln!(
        "[round_to_f64] {} escalated to 256, {} fell back, over {asked} near-axis angles",
        after.0 - before.0,
        after.1 - before.1
    );
    assert!(asked > 500, "corpus shrank to {asked}");
    assert_eq!(
        after.1, before.1,
        "the undecidable fallback fired, which is derived not to"
    );
}

/// The memo answers the second ask, keys on the angle's *value* rather than its spelling — same
/// failure mode as [`TRIG`]'s, showing up as work done twice rather than a wrong answer — **and
/// keys on the realized pair**, which is the part that is about correctness rather than cost.
#[test]
fn the_realization_error_memo_keys_on_the_angle_and_the_pair() {
    let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
    let entries = || F64_ERR.with_borrow(|m| m.len());
    // An angle no other test in this thread asks for.
    let a = deg(1234567, 9973);
    let (c, s) = a.cos_sin_f64();
    let before = entries();
    let first = a.realization_error_of(c, s);
    assert_eq!(entries(), before + 1, "the first ask must insert");
    assert_eq!(a.realization_error_of(c, s), first);
    assert_eq!(
        entries(),
        before + 1,
        "the second ask must be answered, not recomputed"
    );
    // The same angle, unreduced, with the same pair, is the same entry.
    assert_eq!(deg(2469134, 19946).realization_error_of(c, s), first);
    assert_eq!(
        entries(),
        before + 1,
        "the key is the angle, not its spelling"
    );
    // ★ A neighbouring realization of that same angle is a *different question* and must get
    // its own answer — an entry keyed on the angle alone would hand back `first`, an error
    // measured against a pair this caller is not holding.
    let nudged = f64::from_bits(c.to_bits() + 1);
    assert_ne!(a.realization_error_of(nudged, s), first);
    assert_eq!(entries(), before + 2, "the pair is part of the key");
}

/// `apply_point`/`apply_dir` realize a 90°-family rotation bit-exactly: no ~6e-17
/// spurious cross-term. (3,5,z) about Z by 90° → exactly (-5,3,z); a rational-pivot
/// rotation is exact too; a non-quadrantal angle is unchanged from the f64 path.
#[test]
fn apply_point_exact_for_quadrantal() {
    let iso = |d| {
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(d)).unwrap(),
        })
    };
    assert_eq!(iso(90).apply_point([3.0, 5.0, 7.0]), [-5.0, 3.0, 7.0]);
    assert_eq!(iso(180).apply_point([3.0, 5.0, 7.0]), [-3.0, -5.0, 7.0]);
    assert_eq!(iso(270).apply_point([3.0, 5.0, 7.0]), [5.0, -3.0, 7.0]);
    assert_eq!(iso(90).apply_dir([0.0, 1.0, 0.0]), [-1.0, 0.0, 0.0]);

    // Non-origin pivot (2,2): 90° maps (3,5)→pivot+R(1,3)=(2-3, 2+1)=(-1,3).
    let piv = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(2), Rat::from_int(2), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
    });
    assert_eq!(piv.apply_point([3.0, 5.0, 0.0]), [-1.0, 3.0, 0.0]);

    // Non-quadrantal: the same arithmetic on the same realized pair, so this pins the *route*
    // (`px + u·c − v·s`, in that order) rather than the values.
    let a = Angle::from_deg(Rat::from_int(37)).unwrap();
    let (c, s) = a.cos_sin_f64();
    assert_eq!(
        iso(37).apply_point([3.0, 5.0, 0.0]),
        [3.0 * c - 5.0 * s, 3.0 * s + 5.0 * c, 0.0]
    );
}

/// H1.5: the arbitrary-precision cos/sin realization must be far more accurate
/// than any f64/double-double — the accuracy gate astro-float passes and the
/// double-double `twofloat` failed (its trig degraded to ~1e-16 near zero-
/// crossings). Error at rational angles must beat `2^-100` (~1e-30).
#[test]
fn high_precision_trig_meets_gate() {
    const GATE_EXP: i32 = -100; // 2^-100 ≈ 7.9e-31
    // (deg, exact_cos, exact_sin | None where irrational, e.g. sin 60°)
    let cases = [
        (0i128, 1.0, Some(0.0)),
        (60, 0.5, None),
        (90, 0.0, Some(1.0)),
        (120, -0.5, None),
        (180, -1.0, Some(0.0)),
        (270, 0.0, Some(-1.0)),
    ];
    let mut worst = i32::MIN;
    for (deg, exact_cos, exact_sin) in cases {
        let a = Angle::from_deg(Rat::from_int(deg)).unwrap();
        let ec = hp_err_exp(&a.cos_sin_at(GT_PREC).0, exact_cos);
        assert!(ec < GATE_EXP, "cos {deg}° error 2^{ec} exceeds gate");
        worst = worst.max(ec);
        if let Some(s) = exact_sin {
            let es = hp_err_exp(&a.cos_sin_at(GT_PREC).1, s);
            assert!(es < GATE_EXP, "sin {deg}° error 2^{es} exceeds gate");
            worst = worst.max(es);
        }
    }
    eprintln!("[H1.5] astro-float {GT_PREC}-bit worst cos/sin error at rational angles: 2^{worst}");
}

/// A large angle normalizes in one step, not one step per turn.
///
/// Normalization used to subtract 360° in a loop, so `2⁶⁰` degrees needed ~3·10¹⁵ iterations —
/// the kernel did not reject that input, it stopped responding to it. Found by execution: a
/// trig corpus reached for a big numerator and the test never returned.
#[test]
fn a_huge_angle_normalizes_without_counting_turns() {
    let huge = Rat::new(1i128 << 60, 7).unwrap();
    let a = Angle::from_deg(huge).expect("a large rational angle is representable");
    assert!(a.deg() >= Rat::from_int(0) && a.deg() < Rat::from_int(360));
    // Same residue class as the input, so the reduction is `− 360k`, not a different angle.
    let back = a.deg().checked_sub(huge).unwrap();
    let turns = back.to_f64() / -360.0;
    assert_eq!(
        turns.fract(),
        0.0,
        "the reduction was not a whole number of turns"
    );
    // Negatives land in range too, and exactly on 0 at a full turn.
    assert_eq!(
        Angle::from_deg(Rat::from_int(-720)).unwrap().deg(),
        Rat::from_int(0)
    );
    assert_eq!(
        Angle::from_deg(Rat::new(-1, 2).unwrap()).unwrap().deg(),
        Rat::new(719, 2).unwrap()
    );
}

/// **The memo memoizes, and its key is the *value* of the angle rather than its spelling.**
///
/// Two failure modes, and neither shows up as a wrong answer — [`Angle::cos_sin_bounded`] is a
/// pure function either way, so a broken memo is only slow. That is exactly why it needs a
/// test: a rotated boolean spent 22% of its time re-running Ziv's loop for angles it had
/// already realized, and nothing but a measurement would say so.
///
/// - **It does not memoize** (a rewrite drops the lookup): asking twice would insert twice.
/// - ★ **The key fragments**: `90/1` and `180/2` are the same angle. If they hashed apart the
///   memo would still be *correct* and still show a high hit rate, while paying twice for every
///   angle a caller happened to spell in two ways. `Angle::from_deg` reduces through
///   `Rat::new`, so they are one entry — this is the guard on that.
///
/// Deltas, not absolute counts: the memo is a `thread_local`, and the harness gives each test
/// its own thread, but nothing here should depend on which tests ran first.
#[test]
fn the_trig_memo_keys_on_the_angle_not_its_spelling() {
    let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
    // A precision no other test asks for, so this thread's map cannot be pre-warmed for it.
    let prec = 704;
    let before = trig_entries();
    let first = deg(37, 1).cos_sin_bounded(prec);
    assert_eq!(trig_entries(), before + 1, "the first ask must insert");

    let again = deg(37, 1).cos_sin_bounded(prec);
    assert_eq!(
        trig_entries(),
        before + 1,
        "the second ask must be answered, not recomputed"
    );
    assert_eq!(
        (first.0.value.clone(), first.0.error),
        (again.0.value.clone(), again.0.error),
        "and answered with the same value"
    );

    // ★ The same angle, spelled as an unreduced ratio, is the same entry.
    let spelled = deg(74, 2).cos_sin_bounded(prec);
    assert_eq!(
        trig_entries(),
        before + 1,
        "74/2 is 37/1: a second entry means the key is the spelling, not the angle"
    );
    assert_eq!(
        (first.0.value.clone(), first.0.error),
        (spelled.0.value.clone(), spelled.0.error)
    );

    // …and precision *is* part of the key: a different depth is a different answer, so reusing
    // an entry across depths would hand back coordinates realized at the wrong one.
    let deeper = deg(37, 1).cos_sin_bounded(prec + 64);
    assert_eq!(trig_entries(), before + 2, "precision must key the memo");
    assert_ne!(
        deeper.0.error, again.0.error,
        "a deeper realization has a smaller bound"
    );
}

/// **The seed of every error radius, checked against a realization far deeper than itself.**
///
/// [`Angle::cos_sin_bounded`] derives its bound from two things the crate does not promise in
/// writing: that `Consts::pi` is correctly rounded, and that `cos`/`sin` are too (they run
/// Ziv's loop, which is how one builds a correctly-rounded transcendental — but an
/// implementation detail, not a documented contract). If either weakens, every interval above
/// it is unsound, so the claim is measured: at each rung the ladder uses, the value must sit
/// within its own bound of the same value realized with 512 extra bits.
///
/// The reference is not independent code — it is the same routine at higher precision — so
/// this cannot catch an error that grows with precision in the same shape. What it does catch
/// is the failure that matters here: a bound that is simply too small.
#[test]
fn the_trig_bound_holds_against_a_far_deeper_realization() {
    // Angles spanning the quadrants, plus rationals with awkward denominators and one whose
    // numerator is large enough to exercise the `i128 → f64` term.
    let angles = [
        (0i128, 1i128),
        (30, 1),
        (45, 1),
        (60, 1),
        (90, 1),
        (135, 1),
        (180, 1),
        (271, 1),
        (359, 1),
        (1, 7),
        (22, 7),
        (1000, 3),
        (1, 1_000_000),
        // A denominator past 2⁵³ that is still exact in `f64` (a power of two), so the bound
        // must *not* charge it the conversion term.
        (1, 1i128 << 60),
        // …and one that genuinely does not round-trip, where the angle itself is only known
        // to a relative `2⁻⁵³` and no working precision can recover it.
        (1, (1i128 << 60) + 1),
    ];
    // Slack is tracked per rung so the two stories stay separable: a word-aligned precision
    // is delivered as asked, while `200` is silently rounded up to 256 and the extra bits show
    // up as slack that is astro-float's, not this bound's.
    let mut worst = std::collections::BTreeMap::<usize, (i64, String)>::new();
    for prec in [128usize, 200, 256, 512, 1024] {
        for (num, den) in angles {
            let a = Angle::from_deg(Rat::new(num, den).unwrap()).unwrap();
            let (c, s) = a.cos_sin_bounded(prec);
            let (bc, bs) = (c.error, s.error);
            let (c, s) = (c.value, s.value);
            let deep = prec + 512;
            let (rc, rs) = a.cos_sin_at(deep);
            for (got, reference, bound, what) in [(&c, &rc, bc, "cos"), (&s, &rs, bs, "sin")] {
                let diff = got.sub(reference, deep, HP_RM);
                let Some(de) = (if diff.is_zero() {
                    None
                } else {
                    diff.exponent()
                }) else {
                    continue; // exactly equal — nothing to bound
                };
                // `|diff| < 2^de`; the bound must be at least that.
                let observed = Mag::pow2(de as i64);
                assert!(
                    !bound.lt(observed),
                    "{what} {num}/{den}° at {prec} bits: error 2^{de} exceeds its bound 2^{:?}",
                    bound.exp2()
                );
                // Track how much slack the bound carries, so a bound that is merely
                // enormous does not pass as a bound that is right.
                // Slack is only meaningful where the angle converts exactly. Where it does
                // not, the `2⁻⁵³` term dominates by design and the gap to the observed error
                // is the honest cost of an unrepresentable angle, not looseness.
                let angle_exact = (num as f64) as i128 == num && (den as f64) as i128 == den;
                // And only where the two realizations actually disagree above the
                // *reference's* own resolution. `cos 60° = 1/2` is exact at both precisions,
                // so their difference measures the reference, not this bound.
                let above_reference_noise = (de as i64) > -(deep as i64) + 8;
                if let (Some(be), true) = (bound.exp2(), angle_exact && above_reference_noise) {
                    let slack = be - de as i64;
                    let e = worst.entry(prec).or_insert((i64::MIN, String::new()));
                    if slack > e.0 {
                        *e = (
                            slack,
                            format!("{what} {num}/{den}°, bound 2^{be} vs error 2^{de}"),
                        );
                    }
                }
            }
        }
    }
    for (prec, (slack, at)) in &worst {
        eprintln!("[cip] {prec}-bit trig bound slack: 2^{slack}  ({at})");
    }
    for (&prec, (slack, at)) in &worst {
        if prec % 64 == 0 {
            assert!(
                *slack <= 16,
                "at {prec} bits the bound is 2^{slack} above the worst observed error ({at}) — \
                     that is a fudge factor wearing a derivation's clothes, not a tight bound"
            );
        }
    }
    // …and the odd rung out proves why stage 4's rungs are multiples of 64: asking for 200
    // bits buys 256, so ~56 bits of the result are paid for and then claimed away.
    let (slack_200, _) = &worst[&200];
    assert!(
        (40..=80).contains(slack_200),
        "expected ~56 bits of unclaimed precision at 200 bits (the word-size round-up), got \
             2^{slack_200}"
    );
}

/// Squared lengths a frame normal actually produces, plus awkward ones.
///
/// `(0,0,1)` and `(1,1,0)` give `1` and `2`; a profile edge `(1,2)` gives a wall normal
/// `(2,−1,0)` and so `5`; the fractions are what a normal reduced by its own content leaves;
/// and the last two exercise the `i128 → f64` term at and past 2⁵³.
const INV_SQRT_CASES: [(i128, i128); 12] = [
    (1, 1),
    (2, 1),
    (3, 1),
    (5, 1),
    (4, 1),
    (9, 25),
    (1, 2),
    (13, 7),
    (1_000_003, 3),
    (1, 1_000_000),
    (1, 1i128 << 60),
    (1, (1i128 << 60) + 1),
];

#[test]
fn the_inverse_sqrt_bound_holds_against_a_far_deeper_realization() {
    let mut worst = std::collections::BTreeMap::<usize, (i64, String)>::new();
    for prec in [128usize, 256, 512] {
        for (num, den) in INV_SQRT_CASES {
            let v = Rat::new(num, den).unwrap();
            let HpBounded {
                value: z,
                error: bound,
            } = inv_sqrt_bounded(v, prec).unwrap();
            let deep = prec + 512;
            let rz = inv_sqrt_bounded(v, deep).unwrap().value;
            let diff = z.sub(&rz, deep, HP_RM);
            let Some(de) = (if diff.is_zero() {
                None
            } else {
                diff.exponent()
            }) else {
                continue; // exactly equal — nothing to bound
            };
            // `|diff| < 2^de`; the bound must be at least that.
            assert!(
                !bound.lt(Mag::pow2(de as i64)),
                "1/sqrt({num}/{den}) at {prec} bits: error 2^{de} exceeds its bound 2^{:?}",
                bound.exp2()
            );
            // ★★ Slack is only meaningful where the observation actually *measures* this
            // bound. A `prec`-bit realization of a value of magnitude `2^zexp` must carry an
            // error near `2^(zexp − prec)`; when the observed gap is far below that, the two
            // realizations agreed better than either one's own accuracy warrants — measured,
            // `1/√(1/(2⁶⁰+1))` at 128 bits agrees to 154 bits where 98 is all that is earned,
            // and `1/√(1/10⁶)` is exactly `1000` at every precision. Those samples understate
            // the true error, so counting them would read a *lucky observation* as a loose
            // bound. Same shape as the trig test's reference-noise filter.
            let earned = z.exponent().unwrap_or(0) as i64 - prec as i64 - 8;
            if let (Some(be), true) = (bound.exp2(), (de as i64) >= earned) {
                let slack = be - de as i64;
                let e = worst.entry(prec).or_insert((i64::MIN, String::new()));
                if slack > e.0 {
                    *e = (slack, format!("{num}/{den}, bound 2^{be} vs error 2^{de}"));
                }
            }
        }
    }
    for (prec, (slack, at)) in &worst {
        eprintln!("[cip] {prec}-bit inverse-sqrt bound slack: 2^{slack}  ({at})");
    }
    for (prec, (slack, at)) in &worst {
        assert!(
            *slack <= 8,
            "at {prec} bits the bound is 2^{slack} above the worst observed error ({at}) — \
                 that is a fudge factor wearing a derivation's clothes, not a tight bound"
        );
    }
}

/// ★ The claim the axis-aligned corpus rests on: a frame whose normal is a coordinate
/// direction realizes its `1/|n|` **exactly**, so it never reaches arbitrary precision and
/// carries no realization error at all.
#[test]
fn an_axis_aligned_frame_needs_no_arbitrary_precision() {
    // `n·n` for (0,0,1) and (2,0,0) — an axis-aligned face and a wall — then the Pythagorean
    // ones, where the answer is exactly rational but (`5/3`) need not be an f64.
    for (v, want) in [
        ((1, 1), 1.0),
        ((4, 1), 0.5),
        ((25, 1), 0.2),
        ((9, 25), 5.0 / 3.0),
    ] {
        let r = Rat::new(v.0, v.1).unwrap();
        assert_eq!(
            inv_sqrt_exact(r).map(Rat::to_f64),
            Some(want),
            "1/sqrt{v:?}"
        );
        assert_eq!(inv_sqrt_f64(r), Some(want));
    }
    // …and a tilted one is honestly irrational, so it falls to the ladder.
    for v in [(2, 1), (3, 1), (5, 1), (1, 2)] {
        assert_eq!(inv_sqrt_exact(Rat::new(v.0, v.1).unwrap()), None, "{v:?}");
    }
    // There is no frame with a non-positive squared length; that is a broken premise.
    for v in [(0, 1), (-1, 1)] {
        let r = Rat::new(v.0, v.1).unwrap();
        assert_eq!(inv_sqrt_exact(r), None);
        assert_eq!(inv_sqrt_f64(r), None);
        assert!(inv_sqrt_bounded(r, 128).is_none());
    }
}

/// **Correctly rounded, which is what keeps debug and release identical** — the failure this
/// crate already paid for once in the trig path. Checked against a realization deep enough
/// that its own rounding cannot reach the 53rd bit.
#[test]
fn inv_sqrt_f64_lands_on_the_nearest_f64() {
    // Deltas, not absolute counts: the memo is a `thread_local` and the harness may have run
    // other tests on this thread first.
    let before = INV_SQRT_ESCALATED.with_borrow(|c| *c);
    for (num, den) in INV_SQRT_CASES {
        let v = Rat::new(num, den).unwrap();
        let deep = inv_sqrt_bounded(v, 1024).unwrap().value;
        let want = to_f64_exact(&deep).unwrap();
        assert_eq!(
            inv_sqrt_f64(v),
            Some(want),
            "1/sqrt({num}/{den}) is not the nearest f64"
        );
    }
    let after = INV_SQRT_ESCALATED.with_borrow(|c| *c);
    let (deeper, undecided) = (after.0 - before.0, after.1 - before.1);
    eprintln!("[cip] inverse-sqrt escalations: {deeper} deeper, {undecided} undecided");
    // ★ 128 bits sufficed for every case here, and the second rung is a proven cap rather
    // than a rung anything reaches: `1/√v` is never zero for a positive `v`, so its interval
    // cannot straddle zero the way `cos 90°`'s does and the rounding always decides.
    assert_eq!(
        (deeper, undecided),
        (0, 0),
        "a case escalated — the ladder's first rung no longer covers the population"
    );
}

/// **The reported realization error really does cover the error that is there** — checked
/// against a realization far deeper than the one the measurement itself uses (128 bits), so
/// the ruler and the thing being measured are not the same instrument.
#[test]
fn the_inverse_sqrt_realization_error_covers_the_error_that_is_there() {
    let mut worst = (i32::MIN, String::new());
    for (num, den) in INV_SQRT_CASES {
        let v = Rat::new(num, den).unwrap();
        let f = inv_sqrt_f64(v).unwrap();
        let reported = inv_sqrt_error_of(v, f).unwrap();
        let HpBounded {
            value: deep,
            error: deep_rad,
        } = inv_sqrt_bounded(v, 1024).unwrap();
        let true_err = hp_err_exp(&deep, f);
        // A `BigFloat` is `m · 2^e` with `m ∈ [0.5, 1)`, so the exponent gives
        // `2^(e−1) ≤ |f − deep| < 2^e` — the *lower* end is what a bound has to clear. Using
        // `2^e` would demand the reported value exceed an over-estimate.
        //
        // ★ **And the reference is not exact either**, so its own radius comes off:
        // `|f − true| ≥ |f − deep| − rad`. Without that, `1/√(1/10⁶)` — which is exactly
        // `1000.0`, correctly reported as a zero error — fails against the 1024-bit
        // realization's `2⁻¹⁰¹³` residue, and the test would be measuring its own ruler.
        let ref_rad = deep_rad.exp2().map_or(0.0, |e| 2f64.powi(e as i32));
        let floor = if true_err == i32::MIN {
            0.0
        } else {
            (2f64.powi(true_err - 1) - ref_rad).max(0.0)
        };
        assert!(
            reported >= floor,
            "1/sqrt({num}/{den}): reported {reported:e} is below the true error 2^{true_err}"
        );
        if true_err > worst.0 {
            worst = (true_err, format!("{num}/{den}"));
        }
    }
    // ★★★★ **Exact is not the same as dyadic, and only dyadic earns the zero.** An
    // axis-aligned frame lands on `1`/`½` and must report a literal `0.0` — that is what keeps
    // it on the exact predicate path. A Pythagorean one lands on `1/5` or `5/3`, which are
    // exactly rational and still not f64, so they must report the rounding that really
    // happened. Reporting `0` there (the shape the trig path can safely use) is the bug this
    // pair of assertions exists to keep out.
    // `1/9` is here rather than below because `1/√(1/9)` is `3` — an integer is dyadic too.
    for (num, den) in [(1i128, 1i128), (4, 1), (1, 4), (1, 64), (1, 9)] {
        let v = Rat::new(num, den).unwrap();
        let f = inv_sqrt_f64(v).unwrap();
        assert_eq!(inv_sqrt_error_of(v, f), Some(0.0), "1/sqrt({num}/{den})");
    }
    // `(3,4,0)` gives `n·n = 25` and `1/|n| = 1/5`; `(3,4,0)` scaled gives `9/25` and `5/3`.
    for (num, den) in [(25i128, 1i128), (9, 25), (49, 1)] {
        let v = Rat::new(num, den).unwrap();
        let f = inv_sqrt_f64(v).unwrap();
        assert!(
            inv_sqrt_exact(v).is_some() && inv_sqrt_error_of(v, f).unwrap() > 0.0,
            "1/sqrt({num}/{den}) is exactly rational but not an f64 — it must be charged"
        );
    }
    eprintln!(
        "[cip] worst inverse-sqrt f64 realization error: 2^{} ({})",
        worst.0, worst.1
    );
}

/// ★★ **The naive f64 route is not good enough, measured** — `1.0 / v.to_f64().sqrt()` is
/// three roundings and lands on the wrong `f64` often enough to see. Without a case that
/// actually differs, the exact realization above would be a cost with no effect, and this
/// test would be indistinguishable from one that is not running.
#[test]
fn the_naive_f64_route_gets_the_last_bit_wrong() {
    let (mut n, mut differ) = (0usize, 0usize);
    for num in 1i128..400 {
        for den in 1i128..7 {
            let v = Rat::new(num, den).unwrap();
            let naive = 1.0 / v.to_f64().sqrt();
            n += 1;
            if inv_sqrt_f64(v) != Some(naive) {
                differ += 1;
            }
        }
    }
    eprintln!("[cip] naive 1/sqrt disagrees with the correctly rounded one: {differ}/{n}");
    assert!(
        differ > 0,
        "the naive route agreed everywhere on {n} cases — either the correctly rounded path \
             is not running, or this population cannot tell them apart"
    );
}

/// Error of a high-precision realization `a` against an exact f64 `b`, as a
/// power-of-two exponent (`i32::MIN` when exactly equal).
fn hp_err_exp(a: &BigFloat, b: f64) -> i32 {
    let err = a.sub(&BigFloat::from_f64(b, GT_PREC), GT_PREC, HP_RM);
    if err.is_zero() {
        i32::MIN
    } else {
        err.exponent().unwrap_or(i32::MIN)
    }
}

// Small ranges keep checked arithmetic inside i128, so these exercise the
// algebraic laws, not the overflow path (which its own test above pins).
prop_compose! {
    fn small_rat()(num in -1000i128..=1000, den in 1i128..=1000) -> Rat {
        Rat::new(num, den).unwrap()
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Rational `+` and `×` are commutative (exact, no drift).
    #[test]
    fn add_and_mul_commute(a in small_rat(), b in small_rat()) {
        prop_assert_eq!(a.checked_add(b), b.checked_add(a));
        prop_assert_eq!(a.checked_mul(b), b.checked_mul(a));
    }

    /// `a + b + c` associates regardless of grouping.
    #[test]
    fn add_associates(a in small_rat(), b in small_rat(), c in small_rat()) {
        let left = a.checked_add(b).and_then(|ab| ab.checked_add(c));
        let right = b.checked_add(c).and_then(|bc| a.checked_add(bc));
        prop_assert_eq!(left, right);
    }

    /// `from_deg` always normalizes into `[0, 360)`.
    #[test]
    fn from_deg_normalizes_into_range(deg in -3600i128..=3600) {
        let a = Angle::from_deg(Rat::from_int(deg)).unwrap();
        prop_assert!(a.deg() >= Rat::from_int(0));
        prop_assert!(a.deg() < Rat::from_int(360));
    }
}
// ---- the cylinder's strip on a wall plane ----

/// The running fixture: the wall plane `y = 12`, and a cylinder on the vertical line
/// `x = 8, y = 10` with `r = 3`. The axis stands `d = 2` from the plane, so the strip has
/// half-width `h = √(9 − 4) = √5 ≈ 2.236` about `x = 8`, and `e = n × m = x̂` — which makes
/// `U` simply `p.x − 8`.
mod strip {
    // `crate::*`, not `super::*`: a glob import does not re-export what the parent glob-imported.
    use crate::*;

    fn r3(v: [i128; 3]) -> [Rat; 3] {
        v.map(Rat::from_int)
    }
    fn wall(d: i128) -> [Rat; 4] {
        [0, 1, 0, -d].map(Rat::from_int) // y = d
    }
    fn at(x: i128, y: i128, z: i128) -> MeetPoint {
        MeetPoint::Narrow(r3([x, y, z]))
    }
    fn side(coeffs: &[Rat; 4], p: &MeetPoint, ox: i128) -> StripSide {
        cylinder_strip_side(
            coeffs,
            p,
            &r3([ox, 10, -1]),
            &r3([0, 0, 1]),
            &BigRat::from(Rat::from_int(9)), // r = 3, stated as r²
        )
    }

    #[test]
    fn a_point_each_side_of_the_strip_is_clear_and_says_which_side() {
        assert_eq!(side(&wall(12), &at(11, 12, 0), 8), StripSide::Plus);
        assert_eq!(side(&wall(12), &at(5, 12, 0), 8), StripSide::Minus);
    }

    #[test]
    fn a_point_inside_the_strip_is_inside() {
        // |U| = 2, and 4 < 5 — the point stands 2√2 from the axis, inside r = 3.
        assert_eq!(side(&wall(12), &at(10, 12, 0), 8), StripSide::Inside);
        // On the axis' own in-plane line: U = 0, answered before the magnitude test.
        assert_eq!(side(&wall(12), &at(8, 12, 0), 8), StripSide::Inside);
    }

    /// **Exact tangency, with no irrational to dodge.** Put the plane *through* the axis
    /// (`d = 0`), where the half-width is `r` itself — then a point at `U = r` sits exactly
    /// on the lateral surface, and "boundary included" is a statement this test can make.
    #[test]
    fn a_point_exactly_on_the_surface_is_inside() {
        assert_eq!(side(&wall(10), &at(11, 10, 0), 8), StripSide::Inside);
        assert_eq!(side(&wall(10), &at(12, 10, 0), 8), StripSide::Plus);
    }

    /// ★ **The negative control.** The same point, the same plane — move the axis under it
    /// and the answer must change. Without this, every assertion above could be passing for
    /// a reason that has nothing to do with the cylinder.
    #[test]
    fn moving_the_axis_moves_the_answer() {
        let p = at(11, 12, 0);
        assert_eq!(side(&wall(12), &p, 8), StripSide::Plus); // U = 3, 9 > 5
        assert_eq!(side(&wall(12), &p, 10), StripSide::Inside); // U = 1, 1 < 5
    }

    /// A plane that already clears the axis by more than `r` cuts **no** strip, so every
    /// point off the axis' line is clear — including one the near plane called `Inside`.
    #[test]
    fn a_plane_clear_of_the_axis_cuts_no_strip() {
        assert_eq!(side(&wall(12), &at(10, 12, 0), 8), StripSide::Inside);
        assert_eq!(side(&wall(14), &at(10, 14, 0), 8), StripSide::Plus);
    }

    /// **The consumer-side net for the on-plane precondition.** `cylinder_strip_side` takes
    /// the plane's distance from the axis as the point's, so a point that does not satisfy the
    /// coefficients is judged against geometry that is not there. The gate asks this first.
    #[test]
    fn a_point_is_on_the_plane_only_when_it_satisfies_the_coefficients() {
        let w = wall(12);
        assert!(point_on_plane_exact(&w, &at(11, 12, 0)));
        assert!(!point_on_plane_exact(&w, &at(11, 13, 0)));
        // Scaling the coefficients states the same plane, so it cannot move the answer.
        let scaled = [0, 5, 0, -60].map(Rat::from_int);
        assert!(point_on_plane_exact(&scaled, &at(11, 12, 0)));
        assert!(!point_on_plane_exact(&scaled, &at(11, 13, 0)));
    }

    /// The same question at a width `Rat` cannot hold.
    #[test]
    fn a_wide_point_is_placed_against_the_plane_too() {
        use num_bigint::BigInt;
        let huge: BigInt = BigInt::from(10u8).pow(40);
        let on = MeetPoint::Wide([
            (BigInt::from(11) * &huge, huge.clone()),
            (BigInt::from(12) * &huge, huge.clone()),
            (BigInt::from(0), huge.clone()),
        ]);
        let off = MeetPoint::Wide([
            (BigInt::from(11) * &huge, huge.clone()),
            (BigInt::from(12) * &huge + BigInt::from(1), huge.clone()),
            (BigInt::from(0), huge.clone()),
        ]);
        assert!(point_on_plane_exact(&wall(12), &on));
        assert!(!point_on_plane_exact(&wall(12), &off));
    }

    /// ★★ **The one case the `U = 0` early answer exists for, and the only one that measures
    /// it.** On a plane that clears the axis the strip is *empty*, so `U² > (negative)` holds
    /// even at `U = 0` — and a "which side" read off a zero would **invent** one. Removing the
    /// early answer makes this say `Minus` for a point that has no side at all; every other
    /// test in this module stays green through that change (measured), which is why it is
    /// written out rather than assumed to be covered.
    #[test]
    fn a_point_with_no_side_is_never_given_one() {
        assert_eq!(side(&wall(14), &at(8, 14, 0), 8), StripSide::Inside);
    }

    /// The coefficients' and the direction's magnitudes cancel — the derivation says their
    /// denominators enter only as positive squares — so doubling them cannot move an answer.
    #[test]
    fn the_answer_does_not_depend_on_how_the_plane_and_axis_are_scaled() {
        let p = at(11, 12, 0);
        let plain = side(&wall(12), &p, 8);
        let scaled = cylinder_strip_side(
            &[0, 2, 0, -24].map(Rat::from_int),
            &p,
            &r3([8, 10, -1]),
            &r3([0, 0, 7]),
            &BigRat::from(Rat::from_int(9)),
        );
        assert_eq!(plain, scaled);
    }

    /// ★ **A point too wide for `Rat` gets the same answer.** This is why the predicate takes
    /// a `MeetPoint`: a meet that overflows `Rat` is an ordinary point of the geometry, and a
    /// gate that declined it would put a width limit inside a refusal that claims to be about
    /// shape.
    #[test]
    fn a_point_wider_than_rat_travels_the_same_road() {
        use num_bigint::BigInt;
        let huge: BigInt = BigInt::from(10u8).pow(40);
        let wide = MeetPoint::Wide([
            (BigInt::from(11) * &huge, huge.clone()),
            (BigInt::from(12) * &huge, huge.clone()),
            (BigInt::from(0), huge.clone()),
        ]);
        assert!(wide.narrow().is_none(), "the fixture must really be wide");
        assert_eq!(side(&wall(12), &wide, 8), StripSide::Plus);
        assert_eq!(side(&wall(12), &wide, 10), StripSide::Inside);
    }
}

// ---- the same rectangle's other axis: along the cylinder ----

/// The running fixture: the cylinder of `mod strip`, whose axis starts at `z = −1` and runs
/// `+z`. With the raw direction `(0,0,1)` the parameter is simply `t = z + 1`.
mod axis {
    use crate::*;

    fn r3(v: [i128; 3]) -> [Rat; 3] {
        v.map(Rat::from_int)
    }
    fn at(x: i128, y: i128, z: i128) -> MeetPoint {
        MeetPoint::Narrow(r3([x, y, z]))
    }
    fn side(p: &MeetPoint, t: Rat) -> Orient {
        point_axis_side(p, &r3([8, 10, -1]), &r3([0, 0, 1]), t)
    }

    #[test]
    fn a_point_each_side_of_the_plane_says_which() {
        // z = 6 is t = 7, z = 0 is t = 1.
        assert_eq!(side(&at(8, 10, 6), Rat::from_int(5)), Orient::Positive);
        assert_eq!(side(&at(8, 10, 0), Rat::from_int(5)), Orient::Negative);
    }

    /// ★★ **The `Zero` this predicate exists to distinguish.** A face resting exactly on a
    /// band's cap plane is what the footprint reading calls *clear* — the uniform-slab theorem
    /// speaks of the **open** slab. A predicate that folded this into one of the sides would
    /// decide that question here, out of sight of the rule that owns it.
    #[test]
    fn a_point_exactly_on_the_plane_is_zero() {
        assert_eq!(side(&at(8, 10, 4), Rat::from_int(5)), Orient::Zero);
        // Off the axis, at the same height — the parameter does not care how far out it is.
        assert_eq!(side(&at(100, -7, 4), Rat::from_int(5)), Orient::Zero);
    }

    /// A `t` with a denominator is ordinary: `t = 11/2` is the plane `z = 4.5`.
    #[test]
    fn a_fractional_parameter_is_ordinary() {
        let t = Rat::new(11, 2).unwrap();
        assert_eq!(side(&at(8, 10, 5), t), Orient::Positive);
        assert_eq!(side(&at(8, 10, 4), t), Orient::Negative);
    }

    /// ★ **The direction's scale is part of the parameter's meaning, and the predicate honours
    /// it.** `axis(t) = o + t·m` with the *raw* `m`, so stretching `m` sevenfold divides every
    /// parameter by seven — and the same plane keeps the same answer.
    #[test]
    fn the_parameter_lives_in_the_raw_directions_scale() {
        let p = at(8, 10, 8); // t = 9 with m = (0,0,1); t = 9/7 with m = (0,0,7)
        let o = r3([8, 10, -1]);
        let plain = point_axis_side(&p, &o, &r3([0, 0, 1]), Rat::from_int(7));
        let stretched = point_axis_side(&p, &o, &r3([0, 0, 7]), Rat::from_int(1));
        assert_eq!(plain, stretched, "both name the plane z = 6");
        assert_eq!(plain, Orient::Positive);
        // ★ And reading the stretched axis on the *unstretched* scale gets it wrong — which is
        // what proves the scale is doing work here rather than cancelling out.
        assert_eq!(
            point_axis_side(&p, &o, &r3([0, 0, 7]), Rat::from_int(7)),
            Orient::Negative,
            "t = 7 on a sevenfold direction is z = 48, well above the point"
        );
    }

    /// ★ **The negative control.** The same point and the same `t` — move the origin and the
    /// answer must change, or every assertion above could be passing for some other reason.
    #[test]
    fn moving_the_origin_moves_the_answer() {
        let p = at(8, 10, 4);
        let m = r3([0, 0, 1]);
        let t = Rat::from_int(5);
        assert_eq!(point_axis_side(&p, &r3([8, 10, -1]), &m, t), Orient::Zero);
        assert_eq!(
            point_axis_side(&p, &r3([8, 10, -3]), &m, t),
            Orient::Positive
        );
    }

    /// ★★ **A negative parameter is ordinary, and the sign rides in the numerator.** The
    /// derivation clears `dp·d_o·d_m²·td` as a *positive* factor; `Rat` guarantees that by
    /// normalising a negative denominator into the numerator (`Rat::new(11, -2)` is `-11/2`),
    /// and `axis_param_of_plane` builds parameters by inverting a ratio whose numerator may
    /// well be negative. A wall standing below its cylinder — the overhanging boss's plate, at
    /// `t ∈ [−2, 0]` — is exactly this case.
    #[test]
    fn a_negative_parameter_is_ordinary() {
        // The axis starts at z = −1, so t = −1 is z = −2 and t = −3 is z = −4.
        assert_eq!(side(&at(8, 10, -3), Rat::from_int(-1)), Orient::Negative);
        assert_eq!(side(&at(8, 10, 0), Rat::from_int(-1)), Orient::Positive);
        assert_eq!(side(&at(8, 10, -2), Rat::from_int(-1)), Orient::Zero);
        // Written with the sign on the denominator instead — the same plane, the same answers.
        let t = Rat::new(3, -2).unwrap(); // −3/2, i.e. z = −5/2
        assert_eq!(side(&at(8, 10, -2), t), Orient::Positive);
        assert_eq!(side(&at(8, 10, -3), t), Orient::Negative);
    }

    /// A point too wide for `Rat` is an ordinary point of the geometry — the reason this takes
    /// a [`MeetPoint`] rather than a `[Rat; 3]`.
    #[test]
    fn a_point_wider_than_rat_travels_the_same_road() {
        use num_bigint::BigInt;
        let huge: BigInt = BigInt::from(10u8).pow(40);
        let wide = |z: i128| {
            MeetPoint::Wide([
                (BigInt::from(8) * &huge, huge.clone()),
                (BigInt::from(10) * &huge, huge.clone()),
                (BigInt::from(z) * &huge, huge.clone()),
            ])
        };
        assert!(
            wide(6).narrow().is_none(),
            "the fixture must really be wide"
        );
        assert_eq!(side(&wide(6), Rat::from_int(5)), Orient::Positive);
        assert_eq!(side(&wide(4), Rat::from_int(5)), Orient::Zero);
        assert_eq!(side(&wide(0), Rat::from_int(5)), Orient::Negative);
    }
}
