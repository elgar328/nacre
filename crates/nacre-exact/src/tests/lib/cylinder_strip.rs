//! The cylinder's strip on a wall plane, and the same rectangle's other axis.

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
