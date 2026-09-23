//! Planes in exact arithmetic: frames, fixed planes, names from meets, integer predicates, offsets,
//! projection, canonical form.

use super::*;

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
/// carries no direction (`Model::push_plane` returns a `flipped` flag for exactly that
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

/// A triple's sense, exactly — against an independent `Rat` derivation, and across widths.
mod triple_sense {
    use crate::{MeetPoint, Orient, Rat, normal_sense, same_sense, triple_normal};
    use proptest::prelude::*;

    fn r(n: i128, d: i128) -> Rat {
        Rat::new(n, d).expect("a rational")
    }

    /// `(b − a) × (c − a)` in `Rat`, for inputs small enough never to overflow.
    fn rat_cross(p: &[[Rat; 3]; 3]) -> [Rat; 3] {
        let e = |q: [Rat; 3]| -> [Rat; 3] {
            core::array::from_fn(|i| q[i].checked_sub(p[0][i]).expect("small"))
        };
        let (u, v) = (e(p[1]), e(p[2]));
        let t = |i: usize, j: usize| {
            u[i].checked_mul(v[j])
                .expect("small")
                .checked_sub(u[j].checked_mul(v[i]).expect("small"))
                .expect("small")
        };
        [t(1, 2), t(2, 0), t(0, 1)]
    }

    fn rat_dot_sign(a: [Rat; 3], b: [Rat; 3]) -> Orient {
        let d = (0..3).fold(Rat::from_int(0), |acc, i| {
            acc.checked_add(a[i].checked_mul(b[i]).expect("small"))
                .expect("small")
        });
        match d {
            x if x > Rat::from_int(0) => Orient::Positive,
            x if x < Rat::from_int(0) => Orient::Negative,
            _ => Orient::Zero,
        }
    }

    /// Collinear points span no direction.
    #[test]
    fn collinear_points_have_no_normal() {
        let p = |x: i128| MeetPoint::Narrow([r(x, 1), r(2 * x, 1), r(0, 1)]);
        let (a, b, c) = (p(0), p(1), p(3));
        assert_eq!(triple_normal([&a, &b, &c]), None);
    }

    /// Swapping two points reverses the sense; a cyclic turn keeps it.
    #[test]
    fn a_transposition_reverses_and_a_cycle_keeps() {
        let p = |x: i128, y: i128, z: i128| MeetPoint::Narrow([r(x, 1), r(y, 3), r(z, 7)]);
        let (a, b, c) = (p(0, 0, 0), p(1, 0, 2), p(0, 1, 5));
        let n = triple_normal([&a, &b, &c]).expect("a plane");
        assert!(!same_sense(
            &n,
            &triple_normal([&a, &c, &b]).expect("a plane")
        ));
        assert!(same_sense(
            &n,
            &triple_normal([&b, &c, &a]).expect("a plane")
        ));
    }

    proptest! {
        #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
            proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
        ))]
        /// The sign of `n · dir` is the `Rat` derivation's sign, zero included.
        #[test]
        fn normal_sense_is_the_rational_dot_sign(
            pts in proptest::array::uniform3(proptest::array::uniform3((-20i128..=20, 1i128..=9))),
            dir in proptest::array::uniform3((-20i128..=20, 1i128..=9)),
        ) {
            let p: [[Rat; 3]; 3] = pts.map(|q| q.map(|(n, d)| r(n, d)));
            let dir = dir.map(|(n, d)| r(n, d));
            let meets = p.map(MeetPoint::Narrow);
            let cross = rat_cross(&p);
            match triple_normal([&meets[0], &meets[1], &meets[2]]) {
                None => prop_assert!(cross.iter().all(|c| *c == Rat::from_int(0))),
                Some(n) => prop_assert_eq!(normal_sense(&n, dir), rat_dot_sign(cross, dir)),
            }
        }

        /// Two triples on one plane agree exactly when their `Rat` crosses do.
        #[test]
        fn same_sense_is_the_rational_agreement(
            pts in proptest::array::uniform3(proptest::array::uniform3((-20i128..=20, 1i128..=9))),
            w in proptest::array::uniform3((-4i128..=4, -4i128..=4)),
        ) {
            let p: [[Rat; 3]; 3] = pts.map(|q| q.map(|(n, d)| r(n, d)));
            // A second triple on the same plane: affine combinations of the first.
            let q: [[Rat; 3]; 3] = w.map(|(s, t)| core::array::from_fn(|i| {
                let u = p[1][i].checked_sub(p[0][i]).expect("small");
                let v = p[2][i].checked_sub(p[0][i]).expect("small");
                p[0][i]
                    .checked_add(u.checked_mul(Rat::from_int(s)).expect("small"))
                    .expect("small")
                    .checked_add(v.checked_mul(Rat::from_int(t)).expect("small"))
                    .expect("small")
            }));
            let (mp, mq) = (p.map(MeetPoint::Narrow), q.map(MeetPoint::Narrow));
            if let (Some(a), Some(b)) = (
                triple_normal([&mp[0], &mp[1], &mp[2]]),
                triple_normal([&mq[0], &mq[1], &mq[2]]),
            ) {
                prop_assert_eq!(
                    same_sense(&a, &b),
                    rat_dot_sign(rat_cross(&p), rat_cross(&q)) == Orient::Positive
                );
            }
        }

        /// Width does not move a sense: shrinking every point by `2⁻²⁰⁰` (a denominator no `i128`
        /// holds, so the meets are `Wide`) leaves the direction where it was.
        #[test]
        fn a_wide_triple_faces_where_its_narrow_twin_does(
            pts in proptest::array::uniform3(proptest::array::uniform3(-20i128..=20)),
            dir in proptest::array::uniform3(-20i128..=20),
        ) {
            use num_bigint::BigInt;
            let narrow = pts.map(|q| MeetPoint::Narrow(q.map(Rat::from_int)));
            let big = BigInt::from(1) << 200u32;
            let wide = pts.map(|q| MeetPoint::Wide(q.map(|n| {
                // Lowest terms with a positive denominator: strip the shared powers of two.
                let (mut num, mut den) = (BigInt::from(n), big.clone());
                while num != BigInt::from(0) && &num % 2 == BigInt::from(0) && &den % 2 == BigInt::from(0) {
                    num /= 2;
                    den /= 2;
                }
                if num == BigInt::from(0) {
                    den = BigInt::from(1);
                }
                (num, den)
            })));
            let dir = dir.map(Rat::from_int);
            match (
                triple_normal([&narrow[0], &narrow[1], &narrow[2]]),
                triple_normal([&wide[0], &wide[1], &wide[2]]),
            ) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    prop_assert!(same_sense(&a, &b));
                    prop_assert_eq!(normal_sense(&a, dir), normal_sense(&b, dir));
                }
                other => prop_assert!(false, "collinearity disagrees across widths: {:?}", other),
            }
        }
    }
}
