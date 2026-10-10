use super::*;
use num_bigint::BigInt;

/// A radical realized at precision must agree with the digits of √2, and with itself.
#[test]
fn a_radical_realizes_to_its_own_digits() {
    let two = Rat::from_int(2);
    // √2 = 1.414213562373095048801688724209 6980785696 7187537694...
    //                                        ^places 31-40  ^place 41 is 7, so 40 places
    // round the last digit up: ...5696 -> ...5697. The digit string is public knowledge and
    // the rounding is done here by hand, so this is an oracle and not a restatement.
    let want = "1.4142135623730950488016887242096980785697";
    for p in [256usize, 512, 1024] {
        let HpBounded { value: v, error: e } = sqrt_bounded(two, p).expect("a positive radicand");
        let d = round_to_digits(&v, e, 40).expect("40 places at this precision");
        assert_eq!(d, want, "prec {p}");
    }
    // A rational-valued QuadVal takes the exact road.
    let q = quad::QuadVal::from_rat(Rat::new(1, 8).unwrap());
    let HpBounded { value: v, error: e } = realize_quad(&q, 256).expect("rational");
    assert_eq!(
        round_to_digits(&v, e, 20).as_deref(),
        Some("0.12500000000000000000")
    );
    // A perfect square resolves rather than being approached: 3 + 2·√4 = 7.
    let q = quad::QuadVal::new(Rat::from_int(3), Rat::from_int(2), Rat::from_int(4)).unwrap();
    let HpBounded { value: v, error: e } = realize_quad(&q, 256).expect("resolves");
    assert_eq!(round_to_digits(&v, e, 5).as_deref(), Some("7.00000"));
}

/// A seam point is the centre plus a radius along the normalized perpendicular — the radius
/// stated as its square, exact when that square has a rational root and realized otherwise.
#[test]
fn a_seam_point_lands_on_the_rim() {
    let r = |n: i128| Rat::from_int(n);
    // centre at origin, e1 = +x (already unit), r² = 25 -> (5, 0, 0) exactly.
    let big = |n: i128| BigRat::from(r(n));
    let out = realize_seam_point([r(0); 3], [r(1), r(0), r(0)], &big(25), 256).expect("ok");
    let got: Vec<_> = out
        .iter()
        .map(|b| round_to_digits(&b.value, b.error, 10).expect("decided"))
        .collect();
    assert_eq!(got, ["5.0000000000", "0.0000000000", "0.0000000000"]);
    // e1 = (1,1,0): the seam is at radius/√2 on each of x and y.
    let out = realize_seam_point([r(0); 3], [r(1), r(1), r(0)], &big(1), 512).expect("ok");
    let x = round_to_digits(&out[0].value, out[0].error, 20).expect("decided");
    assert_eq!(x, "0.70710678118654752440");
    // r² = 2 with e1 = +x: no rational radius exists, and the seam is at √2, realized.
    let out = realize_seam_point([r(0); 3], [r(1), r(0), r(0)], &big(2), 512).expect("ok");
    let x = round_to_digits(&out[0].value, out[0].error, 20).expect("decided");
    assert_eq!(x, "1.41421356237309504880");
    // A wide square: `r = 5.000000000000001e-8` has a square no `Rat` holds, and the seam is
    // its exact `r` again — the road a stated radius always took.
    let wide = Rat::from_decimal(5.000000000000001e-8).expect("in the window");
    assert!(
        wide.checked_mul(wide).is_none(),
        "the fixture's square must leave i128"
    );
    let out = realize_seam_point([r(0); 3], [r(1), r(0), r(0)], &BigRat::square_of(wide), 256)
        .expect("ok");
    assert_eq!(
        round_to_digits(&out[0].value, out[0].error, 24).expect("decided"),
        "0.000000050000000000000010"
    );
}

/// `sqrt_f64`: a rational's square comes back as that rational's own `to_f64`, bit for bit —
/// the road every stated radius takes — and a non-square radicand as the correctly rounded
/// root. The oracle for the second half is the hardware `sqrt`, which is correctly rounded
/// **when its input is exact**, so the radicands are dyadic (`n / 2^k`) on purpose: a
/// radicand that itself rounds would compare two roundings against one.
#[test]
fn sqrt_f64_is_exact_on_squares_and_correctly_rounded_otherwise() {
    for k in 1..=50i128 {
        for d in [1i128, 3, 7, 10, 1000] {
            let r = Rat::new(k, d).unwrap();
            assert_eq!(
                sqrt_f64(&BigRat::square_of(r)).map(f64::to_bits),
                Some(r.to_f64().to_bits()),
                "{r:?}"
            );
        }
    }
    for (n, d) in [
        (2i128, 1i128),
        (3, 1),
        (1, 2),
        (5, 4),
        (7, 8),
        (10_000_000_019, 1),
        (3, 1 << 40),
    ] {
        let v = Rat::new(n, d).unwrap();
        assert!(rat_sqrt_exact(v).is_none(), "{v:?} must not be a square");
        let want = (n as f64 / d as f64).sqrt(); // exact quotient, then one correct rounding
        assert_eq!(
            sqrt_f64(&BigRat::from(v)).map(f64::to_bits),
            Some(want.to_bits()),
            "{v:?}"
        );
    }
    assert_eq!(sqrt_f64(&BigRat::zero()), Some(0.0));
    assert_eq!(sqrt_f64(&BigRat::from(Rat::from_int(-1))), None);
    // A wide square's root is the stated radius, bit for bit — the cache's number for a
    // sub-micron cylinder stated to f64's last digit.
    let wide = Rat::from_decimal(5.000000000000001e-8).expect("in the window");
    assert_eq!(
        sqrt_f64(&BigRat::square_of(wide)).map(f64::to_bits),
        Some(wide.to_f64().to_bits())
    );
}

/// A perfect square has its root found at **any** width `Rat` holds — the roots past `2⁵⁵` an
/// `f64` estimate lands too far from to search back to: four radii a user writes to f64's last
/// digit, and `2^k + 3` up to the top of `i128`. `inv_sqrt_exact` reads the same answer.
#[test]
fn a_perfect_square_has_its_root_at_any_width() {
    let decimals = [
        2.7848112834336347,
        0.36681130224225267,
        26.288485480848863,
        27.644747522971063,
    ]
    .map(|x| Rat::from_decimal(x).expect("in the window"));
    let ints = (56..=63).map(|k| Rat::from_int((1i128 << k) + 3));
    for r in decimals.into_iter().chain(ints) {
        let r2 = r.checked_mul(r).expect("the square fits i128");
        assert_eq!(rat_sqrt_exact(r2), Some(r), "{r:?}");
        assert_eq!(
            inv_sqrt_exact(r2),
            Some(Rat::new(r.denom(), r.numer()).unwrap()),
            "{r:?}"
        );
    }
}

/// A sweep for the carry and tie cases the hand-picked table cannot reach.
#[test]
fn the_big_road_agrees_with_the_rat_oracle_under_sweep() {
    // deterministic xorshift
    let mut x: u64 = 0x0020_2609_11c0_ffee;
    let mut next = || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut bad = 0usize;
    let mut carries = 0usize;
    for _ in 0..20000 {
        let num = (next() as i128) - (i64::MAX as i128);
        let den = ((next() % 1_000_000) + 1) as i128;
        let Some(r) = Rat::new(num, den) else {
            continue;
        };
        let got = nearest_f64_big(&BigInt::from(num), &BigInt::from(den)).expect("range");
        if got != r.to_f64() {
            bad += 1;
            if bad < 4 {
                println!("MISMATCH {num}/{den}: big={got:?} rat={:?}", r.to_f64());
            }
        }
    }
    // boundary family: values that round up across a power of two
    for p in [52i32, 53, 54, 60] {
        let two_p = BigInt::from(1i128) << p;
        for off in [-1i128, 0, 1] {
            let num = &two_p * 2 - BigInt::from(1) + BigInt::from(off);
            let den = BigInt::from(2);
            let got = nearest_f64_big(&num, &den).expect("range");
            let want = Rat::new(i128::try_from(&num).expect("fits"), 2)
                .expect("in range")
                .to_f64();
            if got != want {
                bad += 1;
                println!("BOUNDARY p={p} off={off}: big={got:?} rat={want:?}");
            } else {
                carries += 1;
            }
        }
    }
    println!("boundary agreements={carries}");
    assert_eq!(bad, 0, "{bad} disagreements with the Rat oracle");
}

/// A value below `2⁻¹⁰²²` is still a value — the scale must not underflow before it does.
///
/// `2⁻¹⁰⁰⁰` is an ordinary normal `f64`; an earlier spelling answered `0.0` for it and for
/// every subnormal, because `2f64.powi(-k)` flushed before the product did.
#[test]
fn a_tiny_rational_still_names_its_f64() {
    let one = BigInt::from(1);
    let at = |e: usize| nearest_f64_big(&one, &(BigInt::from(1) << e));
    assert_eq!(at(1000), Some(2f64.powi(-500) * 2f64.powi(-500))); // normal
    assert_eq!(at(1040), Some(2f64.powi(-520) * 2f64.powi(-520))); // subnormal
    assert_eq!(at(1074), Some(f64::from_bits(1))); // the smallest subnormal
    assert_eq!(at(1100), Some(0.0)); // genuinely below the range
    // And the ordinary range is unmoved.
    assert_eq!(at(0), Some(1.0));
    assert_eq!(at(10), Some(1.0 / 1024.0));
}

/// The BigInt road must agree with the `Rat` road wherever both can speak.
///
/// An independent oracle: `Rat::to_f64` is correctly rounded by a different route
/// (`nearest_f64` on `u128`), so agreement is evidence, not a restatement.
#[test]
fn the_big_road_agrees_with_the_rat_road() {
    let mut n = 0usize;
    for num in [-97i128, -7, -1, 1, 3, 5, 7, 11, 97, 1234567, -98765432] {
        for den in [1i128, 2, 3, 7, 10, 1024, 999983, 1_000_000_007] {
            let r = Rat::new(num, den).expect("in range");
            let got =
                nearest_f64_big(&BigInt::from(num), &BigInt::from(den)).expect("in f64 range");
            assert_eq!(got, r.to_f64(), "{num}/{den}");
            n += 1;
        }
    }
    // ★★★ **The tie, which a random sweep can never reach.** At `2⁵² + ½` both neighbours
    // are representable, so this is where a tie rule is visible — and where half-away-from-zero
    // (what this shipped first) disagreed with the road the point cache is built on.
    for e in [50u32, 51, 52, 53] {
        let two_e = 1i128 << e;
        let (num, den) = (two_e * 2 + 1, 2i128);
        let big = nearest_f64_big(&BigInt::from(num), &BigInt::from(den)).expect("in range");
        assert_eq!(
            big,
            Rat::new(num, den).expect("in range").to_f64(),
            "2^{e} + 1/2"
        );
    }
    // Wide inputs no `Rat` can hold — the reason this twin exists at all.
    let big = BigInt::from(1i128) << 200;
    let v = nearest_f64_big(&(&big * 3), &big).expect("in range");
    assert_eq!(v, 3.0);
    assert!(n > 50, "oracle swept {n} pairs");
}

/// Does `big_to_ratio` really spell the value? Exact f64s have known ratios.
#[test]
fn a_binary_float_lifts_to_the_ratio_it_is() {
    for (v, n, d) in [
        (0.5f64, 1i64, 2i64),
        (1.0, 1, 1),
        (3.0, 3, 1),
        (-0.25, -1, 4),
        (0.75, 3, 4),
        (1024.0, 1024, 1),
        (-7.0, -7, 1),
    ] {
        let b = BigFloat::from_f64(v, 128);
        let (gn, gd) = big_to_ratio(&b).expect("finite");
        // compare as a reduced fraction
        let lhs = &gn * BigInt::from(d);
        let rhs = &gd * BigInt::from(n);
        assert_eq!(lhs, rhs, "value {v}: got {gn}/{gd}, want {n}/{d}");
    }
}

/// Exact decimals of a rational, including the tie rule and negatives.
#[test]
fn a_rational_prints_its_own_digits() {
    let c = |n: i64, d: i64, p: usize| decimals_of_ratio(&BigInt::from(n), &BigInt::from(d), p);
    assert_eq!(c(1, 3, 10), "0.3333333333");
    assert_eq!(c(2, 3, 10), "0.6666666667");
    assert_eq!(c(-2, 3, 5), "-0.66667");
    assert_eq!(c(1, 2, 0), "0"); // 0.5 ties to even, as IEEE and `Rat::to_f64` do
    assert_eq!(c(3, 2, 0), "2"); // 1.5 ties to even the other way
    assert_eq!(c(7, 1, 3), "7.000");
    assert_eq!(c(0, 1, 4), "0.0000");
    assert_eq!(c(1, 8, 20), "0.12500000000000000000");
}

/// The escalation signal: too few bits must say "undecided", more bits must decide.
#[test]
fn too_few_bits_report_undecided() {
    // A wide interval cannot pin many digits; a zero radius pins all of them.
    let mid = BigFloat::from_f64(0.125, 256);
    assert_eq!(
        round_to_digits(&mid, Mag::ZERO, 10).as_deref(),
        Some("0.1250000000")
    );
    let fuzzy = Mag::pow2(-10); // ±~0.001 — cannot decide the 5th place
    assert_eq!(round_to_digits(&mid, fuzzy, 5), None);
    assert_eq!(round_to_digits(&mid, fuzzy, 1).as_deref(), Some("0.1"));
}

/// **Both ends of the exact door, and the band in the middle where its flag was lying.**
///
/// This door had **no direct test at all** — its two callers only ever hand it coordinates
/// from the middle of the range, so nothing exercised either end. The ends here were *found*
/// by scanning, not computed from the `k` guard, and computing them would have been wrong:
///
/// - **Large end is `2^1024`, not `2^1252`.** The `|k| > 1200` guard never gets to decide it:
///   the scaled value goes infinite first and `is_finite` refuses. Reading the guard and
///   solving for the value describes a branch nothing reaches.
/// - **Small end is `2^-1149`.** There the guard really is what refuses.
/// - ⚠ **Between `2^-1075` and that end, the flag lied.** The value is right — below the
///   smallest subnormal the nearest `f64` genuinely is `0.0`, and
///   `a_tiny_rational_still_names_its_f64` locks exactly that — but it came back marked
///   *exact*, and [`Realized::to_f64`] turns that mark into a proven zero bound. The value
///   stays; the claim is gone.
#[test]
fn the_exact_door_names_both_ends_and_stops_claiming_a_flushed_zero() {
    let one = BigInt::from(1);
    let small = |e: usize| nearest_f64_big_exact(&one, &(BigInt::from(1) << e));
    let big = |e: usize| nearest_f64_big_exact(&(BigInt::from(1) << e), &one);

    // Large end: the last power of two that names an `f64`, and the first that does not.
    assert_eq!(big(1023), Some((8.98846567431158e307, true)));
    assert_eq!(
        big(1024),
        None,
        "2^1024 overflows — finiteness refuses first"
    );
    // Small end: exact all the way down to the smallest subnormal.
    assert_eq!(small(1074), Some((f64::from_bits(1), true)));
    // Below it the nearest `f64` is zero — the value is kept, the exactness claim is not.
    for e in [1075usize, 1100, 1148] {
        assert_eq!(
            small(e),
            Some((0.0, false)),
            "1/2^{e} rounds to zero, and zero does not name it exactly"
        );
    }
    assert_eq!(
        small(1149),
        None,
        "past the guard there is no answer at all"
    );

    // A real zero is still exactly zero, and the flag still tells exact from rounded.
    assert_eq!(
        nearest_f64_big_exact(&BigInt::from(0), &BigInt::from(5)),
        Some((0.0, true))
    );
    assert_eq!(
        nearest_f64_big_exact(&one, &BigInt::from(4)),
        Some((0.25, true))
    );
    let (third, exact) = nearest_f64_big_exact(&one, &BigInt::from(3)).expect("1/3 is in range");
    assert!(!exact, "1/3 is not an f64");
    assert_eq!(third, 1.0 / 3.0);
}
