use super::*;
use num_bigint::BigInt;
use num_rational::Ratio;

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
        let HpBounded { value: v, error: e } = sqrt_hp(two, p).expect("a positive radicand");
        let d = round_to_digits(&v, e, 40).expect("40 places at this precision");
        assert_eq!(d, want, "prec {p}");
    }
    // A rational-valued QuadVal takes the exact road.
    let q = quad::QuadVal::from_rat(Rat::new(1, 8).unwrap());
    let HpBounded { value: v, error: e } = quad_hp(&q, 256).expect("rational");
    assert_eq!(
        round_to_digits(&v, e, 20).as_deref(),
        Some("0.12500000000000000000")
    );
    // A perfect square resolves rather than being approached: 3 + 2·√4 = 7.
    let q = quad::QuadVal::new(Rat::from_int(3), Rat::from_int(2), Rat::from_int(4)).unwrap();
    let HpBounded { value: v, error: e } = quad_hp(&q, 256).expect("resolves");
    assert_eq!(round_to_digits(&v, e, 5).as_deref(), Some("7.00000"));
}

/// A seam point is the centre plus a radius along the normalized perpendicular — the radius
/// stated as its square, exact when that square has a rational root and realized otherwise.
#[test]
fn a_seam_point_lands_on_the_rim() {
    let r = |n: i128| Rat::from_int(n);
    // centre at origin, e1 = +x (already unit), r² = 25 -> (5, 0, 0) exactly.
    let big = |n: i128| BigRat::from(r(n));
    let out = seam_point_hp([r(0); 3], [r(1), r(0), r(0)], &big(25), 256).expect("ok");
    let got: Vec<_> = out
        .iter()
        .map(|b| round_to_digits(&b.value, b.error, 10).expect("decided"))
        .collect();
    assert_eq!(got, ["5.0000000000", "0.0000000000", "0.0000000000"]);
    // e1 = (1,1,0): the seam is at radius/√2 on each of x and y.
    let out = seam_point_hp([r(0); 3], [r(1), r(1), r(0)], &big(1), 512).expect("ok");
    let x = round_to_digits(&out[0].value, out[0].error, 20).expect("decided");
    assert_eq!(x, "0.70710678118654752440");
    // r² = 2 with e1 = +x: no rational radius exists, and the seam is at √2, realized.
    let out = seam_point_hp([r(0); 3], [r(1), r(0), r(0)], &big(2), 512).expect("ok");
    let x = round_to_digits(&out[0].value, out[0].error, 20).expect("decided");
    assert_eq!(x, "1.41421356237309504880");
    // A wide square: `r = 5.000000000000001e-8` has a square no `Rat` holds, and the seam is
    // its exact `r` again — the road a stated radius always took.
    let wide = Rat::from_decimal(5.000000000000001e-8).expect("in the window");
    assert!(
        wide.checked_mul(wide).is_none(),
        "the fixture's square must leave i128"
    );
    let out =
        seam_point_hp([r(0); 3], [r(1), r(0), r(0)], &BigRat::square_of(wide), 256).expect("ok");
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
        assert!(sqrt_rat(v).is_none(), "{v:?} must not be a square");
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

    // The edges, collected so that one wrong row does not hide the next.
    let big = |n: BigInt, d: BigInt| BigRat(Ratio::new(n, d));
    let two = |e: usize| BigInt::from(1) << e;
    let nan = f64::NAN.to_bits();
    let wrong: Vec<_> = [
        (
            // An exact root too wide for `Rat`, sitting exactly on a rounding boundary: the
            // exact branch rounds it once, ties to even.
            "((2⁵³+1)·2²⁰⁰)²",
            big(
                (BigInt::from((1i64 << 53) + 1) * two(200)).pow(2),
                BigInt::from(1),
            ),
            2f64.powi(253).to_bits(),
        ),
        (
            // Irrational, `2⁻⁴⁰¹` above the boundary `1 + 2⁻⁵³`: 256 bits cannot see the gap.
            "(1+2⁻⁵³)² + 2⁻⁴⁰⁰",
            big(
                BigInt::from((1i64 << 53) + 1).pow(2) * two(294) + 1,
                two(400),
            ),
            (1.0 + 2f64.powi(-52)).to_bits(),
        ),
        // An exact root whose nearest `f64` is not normal: `NaN`, as the ladder decides one.
        ("2⁻²²⁰⁰", big(BigInt::from(1), two(2200)), nan),
    ]
    .into_iter()
    .filter_map(|(what, v, want)| {
        let got = sqrt_f64(&v).map(f64::to_bits);
        (got != Some(want)).then(|| format!("{what}: got {got:?}, want {want}"))
    })
    .collect();
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// **A radicand's width bounds its climb.** `(1+2⁻⁵³)² + 2⁻⁴⁰⁰⁰⁰` has its root about `2⁻⁴⁰⁰⁰¹`
/// above the rounding boundary `1 + 2⁻⁵³`, so it needs some forty thousand bits; its own width — a
/// 40,001-bit denominator — is what earns them. A test of its own, so a climb that stops short
/// cannot hide another row.
#[test]
fn a_wide_radicand_climbs_as_deep_as_its_width_needs() {
    let two = |e: usize| BigInt::from(1) << e;
    let v = BigRat(Ratio::new(
        BigInt::from((1i64 << 53) + 1).pow(2) * two(39894) + 1,
        two(40000),
    ));
    assert_eq!(
        sqrt_f64(&v).map(f64::to_bits),
        Some((1.0 + 2f64.powi(-52)).to_bits())
    );
}

/// **The separation bound at its tightest.** `p·2¹⁰⁸ − q·(2⁵⁴−1)² = 1` with `p` a square and `q`
/// not, both 108 bits, so `√(p/q)` sits a relative `2⁻²¹⁷` above the boundary `1 − 2⁻⁵⁴` — two
/// bits inside the bound `2^-(108+111)`. 128 bits cannot see that gap and 256 can, and 256 is
/// exactly the depth 108-bit terms earn; `1/√` of the reciprocal is the same number.
#[test]
fn the_tightest_radicand_needs_the_whole_separation_margin() {
    let p: i128 = 324518553658426654725561982648324;
    let q: i128 = 324518553658426690754359001612287;
    let before = LADDER_CLIMBED.with_borrow(|c| *c);
    let got = [
        sqrt_f64(&BigRat::from(Rat::new(p, q).unwrap())),
        inv_sqrt_f64(Rat::new(q, p).unwrap()),
    ];
    let after = LADDER_CLIMBED.with_borrow(|c| *c);
    assert_eq!(got, [Some(1.0), Some(1.0)]);
    assert_eq!(
        (after.0 - before.0, after.1 - before.1),
        (2, 0),
        "both decided at 256"
    );
}

/// **The bottom of the range is where the nearest `f64` stops being normal**, on every road.
/// `x₀ = 2⁻¹⁰²² − 2⁻¹⁰⁷⁵` is the subnormal-grid midpoint below `MIN_POSITIVE` and ties up to it, so
/// `x₀` and the band above it read `MIN_POSITIVE` and a value below it reads `NaN` — through an
/// exact root (`v = x²`), through `round_to_f64` itself (`None` there), and through the ladder.
#[test]
fn the_bottom_edge_is_where_the_nearest_f64_stops_being_normal() {
    let mp = f64::MIN_POSITIVE;
    let n0: BigInt = ((BigInt::from(1) << 53) - 1) << 25; // x₀ · 2¹¹⁰⁰
    // `x₀` is a boundary itself, so only a zero radius decides it.
    let rows = [
        ("x₀", n0.clone(), Mag::ZERO, Some(mp)),
        ("x₀ + 2⁻¹¹⁰⁰", &n0 + 1, Mag::pow2(-1140), Some(mp)),
        ("x₀ − 2⁻¹¹⁰⁰", &n0 - 1, Mag::pow2(-1140), None),
        ("−x₀", -n0.clone(), Mag::ZERO, Some(-mp)),
    ];
    let at = |n: &BigInt, p: usize| bigint_to_bigfloat(n, p).mul(&pow2(-1100, p), p, HP_RM);
    let mut wrong = Vec::new();
    for (what, n, rad, want) in &rows {
        let read = want.unwrap_or(f64::NAN).to_bits();
        if n.sign() == num_bigint::Sign::Plus {
            let v = BigRat(Ratio::new(n * n, BigInt::from(1) << 2200));
            let got = sqrt_f64(&v).map(f64::to_bits);
            if got != Some(read) {
                wrong.push(format!("sqrt_f64 {what}: got {got:?}, want {read}"));
            }
        }
        let got = round_to_f64(&at(n, 256), *rad);
        if got.map(f64::to_bits) != want.map(f64::to_bits) {
            wrong.push(format!("round_to_f64 {what}: got {got:?}, want {want:?}"));
        }
        let got = climb(|p| Some([HpBounded::new(at(n, p), *rad)]));
        if got != (Some(vec![read]), (0, 0)) {
            wrong.push(format!("ladder {what}: got {got:?}, want {read}"));
        }
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// **A value is rounded once, at its own width.** `1 + 2⁻⁵³ + 2⁻²⁵⁰` with a zero radius lies just
/// above the midpoint between `1` and `1 + 2⁻⁵²`, so its nearest `f64` is `1 + 2⁻⁵²`. Its ends formed
/// at fewer bits than it has (192, a 128-bit caller's guard) are the midpoint itself, which ties to
/// the even side, `1.0`. The value is built at 320 bits; both signs.
#[test]
fn a_value_wider_than_any_caller_is_rounded_at_its_own_width() {
    let p = 320;
    let up = 1.0 + f64::EPSILON;
    let mut wrong = Vec::new();
    for sign in [1.0, -1.0] {
        let signed = |e: i64| pow2(e, p).mul(&BigFloat::from_f64(sign, p), p, HP_RM);
        let x = signed(0)
            .add(&signed(-53), p, HP_RM)
            .add(&signed(-250), p, HP_RM);
        let got = round_to_f64(&x, Mag::ZERO);
        if got != Some(sign * up) {
            wrong.push(format!("sign {sign}: got {got:?}, want {}", sign * up));
        }
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// **A subnormal is rounded once, on its own grid.** `((2j+1)/2 ± 2⁻⁶⁰)·2⁻¹⁰⁷⁴` lies just off a
/// midpoint of the subnormal grid: rounding to 53 bits first would land on the midpoint, and a
/// second rounding would then pick the even side — wrong for half of these rows. The expected value
/// is built from its bits (a subnormal's bits are its grid index), not by rounding again;
/// `j = 2⁵² − 1` rounds up into `MIN_POSITIVE`. A tiny negative keeps its sign, and a value below
/// the grid is `0.0`, never exact.
#[test]
fn a_subnormal_rounds_once_on_its_own_grid() {
    let den = BigInt::from(1) << 1134;
    let mut wrong = Vec::new();
    let bits = |r: Option<(f64, bool)>| r.map(|(v, x)| (v.to_bits(), x));
    for j in [2u64, 3, (1 << 51) + 6, (1 << 51) + 7, (1 << 52) - 1] {
        let mid: BigInt = BigInt::from(2 * j + 1) << 59;
        let even = if j % 2 == 0 { j } else { j + 1 };
        for (what, num, m) in [
            ("+ε", &mid + 1, j + 1),
            ("−ε", &mid - 1, j),
            ("tie", mid.clone(), even),
        ] {
            let got = bits(nearest_f64_big(&num, &den));
            if got != Some((m, false)) {
                wrong.push(format!("j = {j} {what}: got {got:?}, want bits {m}"));
            }
        }
    }
    let one = BigInt::from(1);
    for (what, num, e, want) in [
        ("−2⁻¹¹⁰⁰", -one.clone(), 1100usize, (-0.0f64).to_bits()),
        ("2⁻¹¹⁵⁰", one.clone(), 1150, 0),
        ("2⁻²⁰⁰⁰", one.clone(), 2000, 0),
    ] {
        let got = bits(nearest_f64_big(&num, &(BigInt::from(1) << e)));
        if got != Some((want, false)) {
            wrong.push(format!("{what}: got {got:?}, want bits {want}"));
        }
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// A perfect square has its root found at **any** width `Rat` holds — the roots past `2⁵⁵` an
/// `f64` estimate lands too far from to search back to: four radii a user writes to f64's last
/// digit, and `2^k + 3` up to the top of `i128`. `inv_sqrt_rat` reads the same answer.
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
        assert_eq!(sqrt_rat(r2), Some(r), "{r:?}");
        assert_eq!(
            inv_sqrt_rat(r2),
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
        let got = nearest_f64_big(&BigInt::from(num), &BigInt::from(den))
            .map(|(v, _)| v)
            .expect("range");
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
            let got = nearest_f64_big(&num, &den).map(|(v, _)| v).expect("range");
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
    let at = |e: usize| nearest_f64_big(&one, &(BigInt::from(1) << e)).map(|(v, _)| v);
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
            let got = nearest_f64_big(&BigInt::from(num), &BigInt::from(den))
                .map(|(v, _)| v)
                .expect("in f64 range");
            assert_eq!(got, r.to_f64(), "{num}/{den}");
            n += 1;
        }
    }
    // **The tie, which a random sweep can never reach.** At `2⁵² + ½` both neighbours
    // are representable, so this is where a tie rule is visible — and where half-away-from-zero
    // (what this shipped first) disagreed with the road the point cache is built on.
    for e in [50u32, 51, 52, 53] {
        let two_e = 1i128 << e;
        let (num, den) = (two_e * 2 + 1, 2i128);
        let big = nearest_f64_big(&BigInt::from(num), &BigInt::from(den))
            .map(|(v, _)| v)
            .expect("in range");
        assert_eq!(
            big,
            Rat::new(num, den).expect("in range").to_f64(),
            "2^{e} + 1/2"
        );
    }
    // Wide inputs no `Rat` can hold — the reason this twin exists at all.
    let big = BigInt::from(1i128) << 200;
    let v = nearest_f64_big(&(&big * 3), &big)
        .map(|(v, _)| v)
        .expect("in range");
    assert_eq!(v, 3.0);
    assert!(n > 50, "oracle swept {n} pairs");
}

/// **A wide integer enters a `BigFloat` whole**, from its words: it reads back exactly
/// (`big_to_ratio` is the inverse reading) at `max(bits, floor, 128)` bits in whole words — across
/// word boundaries, for both signs, and for powers of two whose low words are all zero.
#[test]
fn a_wide_integer_enters_a_bigfloat_whole() {
    let wb = astro_float::WORD_BIT_SIZE;
    let mut xs = vec![BigInt::from(0), BigInt::from(1), BigInt::from(3)];
    for b in [31usize, 32, 33, 63, 64, 65, 127, 128, 129, 200, 1000, 4097] {
        let p = BigInt::from(1) << b;
        xs.extend([&p - 1, p.clone(), &p + 12345]);
    }
    let mut wrong = Vec::new();
    for x in xs.iter().flat_map(|x| [x.clone(), -x]) {
        for floor in [0usize, 128, 300] {
            let f = bigint_to_bigfloat(&x, floor);
            let exact = big_to_ratio(&f).is_some_and(|(n, d)| n == &x * d);
            let want = (x.bits() as usize).max(floor).max(128).div_ceil(wb) * wb;
            if !exact || bits_of(&f) != want {
                wrong.push(format!(
                    "{x} at floor {floor}: exact {exact}, {} bits",
                    bits_of(&f)
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
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

/// **The exact door ends only above.** The large end is `2^1024`, where the scaled value stops
/// being finite. Below there is no end: exact down to the smallest subnormal, and under it `0.0` —
/// the nearest `f64` at any depth, never marked exact, because `Realization::to_f64` turns that
/// mark into a proven zero bound. `a_tiny_rational_still_names_its_f64` locks the values; this
/// locks the ends and the flag.
#[test]
fn the_exact_door_ends_only_above_and_never_claims_a_rounded_zero() {
    let one = BigInt::from(1);
    let small = |e: usize| nearest_f64_big(&one, &(BigInt::from(1) << e));
    let big = |e: usize| nearest_f64_big(&(BigInt::from(1) << e), &one);
    let mut rows = vec![
        // The last power of two that names an `f64`, and the first that does not.
        (
            "2^1023".to_string(),
            big(1023),
            Some((8.98846567431158e307, true)),
        ),
        ("2^1024".to_string(), big(1024), None),
        (
            "2^-1074".to_string(),
            small(1074),
            Some((f64::from_bits(1), true)),
        ),
        // A real zero is exactly zero, and the flag tells exact from rounded.
        (
            "0/5".to_string(),
            nearest_f64_big(&BigInt::from(0), &BigInt::from(5)),
            Some((0.0, true)),
        ),
        (
            "1/4".to_string(),
            nearest_f64_big(&one, &BigInt::from(4)),
            Some((0.25, true)),
        ),
        (
            "1/3".to_string(),
            nearest_f64_big(&one, &BigInt::from(3)),
            Some((1.0 / 3.0, false)),
        ),
    ];
    for e in [1075usize, 1100, 1148, 1149, 1150, 4000] {
        rows.push((format!("2^-{e}"), small(e), Some((0.0, false))));
    }
    let bits = |r: Option<(f64, bool)>| r.map(|(v, x)| (v.to_bits(), x));
    let wrong: Vec<_> = rows
        .iter()
        .filter(|(_, got, want)| bits(*got) != bits(*want))
        .map(|(what, got, want)| format!("{what}: got {got:?}, want {want:?}"))
        .collect();
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// What one climb gave, as bits — `NaN` equals itself and `+0.0` differs from `-0.0` — and how much
/// it moved [`LADDER_CLIMBED`]: `(decided at 256, decided past 256)`.
type Climbed = (Option<Vec<u64>>, (usize, usize));

/// The answer [`ladder_f64`] gives for `realize`, with the counts it moved.
fn climb<const N: usize>(realize: impl FnMut(usize) -> Option<[HpBounded; N]>) -> Climbed {
    let before = LADDER_CLIMBED.with_borrow(|c| *c);
    let got = ladder_f64(TRIG_LADDER_CAP, realize).map(|v| v.iter().map(|x| x.to_bits()).collect());
    let after = LADDER_CLIMBED.with_borrow(|c| *c);
    (got, (after.0 - before.0, after.1 - before.1))
}

/// `±2^e` at `prec` bits, exactly, for any `e` — `2f64.powi` is `0.0` or `inf` past `f64`'s range.
fn pow2(e: i64, prec: usize) -> BigFloat {
    let mut v = BigFloat::from_f64(1.0, prec);
    let mut rem = e;
    while rem != 0 {
        let step = rem.clamp(-1000, 1000);
        v = v.mul(
            &BigFloat::from_f64(2f64.powi(step as i32), prec),
            prec,
            HP_RM,
        );
        rem -= step;
    }
    v
}

/// `b + 2^e` for the rounding boundary `b = 1 + 2⁻⁵³`, exact at `prec` bits.
fn above_boundary(e: i64, prec: usize) -> BigFloat {
    pow2(0, prec)
        .add(&pow2(-53, prec), prec, HP_RM)
        .add(&pow2(e, prec), prec, HP_RM)
}

/// **The ladder climbs until each value is decided** — driven through realizations built to sit near
/// a rounding boundary, since nothing in the corpus climbs past 256 bits.
///
/// The boundary is `b = 1 + 2⁻⁵³`, exactly between `1.0` and `1 + 2⁻⁵²`. The 512-bit row is honest
/// but unlucky: each rung's midpoint lies below the boundary by more than the truth lies above it,
/// within its radius, so a ladder that stopped at 256 and read the midpoint would say `1.0`; its
/// midpoints are built at `p + 64` bits so that offset survives. A value outside `f64`'s normal
/// range is decided as `NaN` at the first rung its interval clears the range on one side.
#[test]
fn the_ladder_climbs_until_each_value_is_decided() {
    let hp = |v: BigFloat, radius: i64| HpBounded::new(v, Mag::pow2(radius));
    let bits = |v: &[f64]| Some(v.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    let p2 = |e: i32| 2f64.powi(e);
    let nan = f64::NAN;
    // 128 bits' radius straddles `b`, 256 bits' does not.
    let late = |p: usize| {
        hp(
            above_boundary(-70, p + 64),
            if p == 128 { -60 } else { -100 },
        )
    };
    let cases = [
        (
            "decided at 128",
            climb(|p| Some([hp(pow2(0, p).add(&pow2(-1, p), p, HP_RM), -100)])),
            bits(&[1.5]),
            (0, 0),
        ),
        (
            "decided at 256",
            climb(|p| Some([late(p)])),
            bits(&[1.0 + p2(-52)]),
            (1, 0),
        ),
        (
            "decided at 512",
            climb(|p| {
                let w = p + 64;
                let mid = above_boundary(-300, w).sub(&pow2(-(p as i64) - 4, w), w, HP_RM);
                Some([hp(mid, -(p as i64))])
            }),
            bits(&[1.0 + p2(-52)]),
            (0, 1),
        ),
        (
            "one of three decided at 256",
            climb(|p| Some([hp(pow2(0, p), -100), late(p), hp(pow2(-2, p), -100)])),
            bits(&[1.0, 1.0 + p2(-52), 0.25]),
            (1, 0),
        ),
        (
            "outside the normal range, both sides and both signs",
            climb(|p| {
                let rel = |e: i64, neg: bool| {
                    let v = pow2(e, p);
                    hp(if neg { v.neg() } else { v }, e - p as i64)
                };
                Some([
                    rel(-1100, false),
                    rel(-1100, true),
                    rel(1030, false),
                    rel(1030, true),
                ])
            }),
            bits(&[nan, nan, nan, nan]),
            (0, 0),
        ),
        (
            // Across zero at 128 bits with both ends subnormal, of two signs — refused like a value
            // below the range, but not on one side — and wholly below the range at 256. (A radius
            // past `2⁻¹⁰²²` would make the ends normal and test nothing about the sign.)
            "across zero, then below the range",
            climb(|p| {
                let r = if p == 128 { -1050 } else { -1100 - p as i64 };
                Some([hp(pow2(-1100, p), r)])
            }),
            bits(&[nan]),
            (1, 0),
        ),
        (
            // Across the bottom edge of the normal range at 128 bits (one end subnormal, one
            // normal) — undecided, not outside — and inside it at 256.
            "across the normal range's edge, then inside it",
            climb(|p| {
                let v = pow2(-1022, p).add(&pow2(-1032, p), p, HP_RM);
                let r = if p == 128 { -1030 } else { -1022 - p as i64 };
                Some([hp(v, r)])
            }),
            bits(&[f64::MIN_POSITIVE * (1.0 + p2(-10))]),
            (1, 0),
        ),
        (
            // Past the exponent guard of the radius `round_to_f64` assembles (`2^±10⁶`): the range
            // is read off the exponents, so these decide too.
            "past the radius' exponent guard",
            climb(|p| {
                let rel = |e: i64| hp(pow2(e, p), e - p as i64);
                Some([rel(-1_000_100), rel(1_000_100)])
            }),
            bits(&[nan, nan]),
            (0, 0),
        ),
        ("no realization", climb::<1>(|_| None), None, (0, 0)),
    ];
    let mut wrong: Vec<_> = cases
        .iter()
        .filter(|(_, got, want, counted)| got != &(want.clone(), *counted))
        .map(|(what, got, want, counted)| {
            format!("{what}: got {got:?}, want {:?}", (want, counted))
        })
        .collect();
    // A zero component is an exact zero — radius `0` — and the first rung decides it as `+0.0`.
    let before = LADDER_CLIMBED.with_borrow(|c| *c);
    let n = [BigInt::from(1), BigInt::from(1), BigInt::from(0)];
    let got = unit_vector_f64(&n).map(|v| v[2].to_bits());
    if got != Some(0f64.to_bits()) || LADDER_CLIMBED.with_borrow(|c| *c) != before {
        wrong.push(format!("zero component of (1, 1, 0): got {got:?}"));
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// A value **on** a rounding boundary never decides, and the ladder panics at its cap rather than
/// climbing forever or handing back a midpoint.
#[test]
#[should_panic(expected = "reached its cap (32768 bits) undecided")]
fn a_value_on_a_rounding_boundary_trips_the_cap() {
    let b = |p: usize| pow2(0, p).add(&pow2(-53, p), p, HP_RM);
    let _ = ladder_f64(TRIG_LADDER_CAP, |p| {
        Some([HpBounded::new(b(p), Mag::pow2(-60))])
    });
}
