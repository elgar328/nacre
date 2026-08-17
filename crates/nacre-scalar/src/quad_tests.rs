//! Gates for the quadratic sign tower (M6-1).
//!
//! Two instruments, deliberately different: the case fixtures carry **hand-derived** expected
//! signs (the spec's ladder, exercised arm by arm), and the differential tests compare against
//! a 512-bit astro-float realization — an oracle whose own resolving power is verified by a
//! negative control before it is trusted ([instrument-decides-the-answer]).

use super::*;
use crate::HP_RM;
use astro_float::BigFloat;

fn r(n: i128, d: i128) -> Rat {
    Rat::new(n, d).expect("test rational")
}

fn ri(n: i128) -> Rat {
    Rat::from_int(n)
}

fn qv(a: Rat, b: Rat, c: Rat) -> QuadVal {
    QuadVal::new(a, b, c).expect("test value exists")
}

/// The 512-bit oracle: realize `a + b·√c` and read the sign, calling anything below 1e−120
/// zero. Our fixtures either separate far above that or are constructed to be exactly zero —
/// `the_oracle_actually_resolves_the_adversarial_gap` proves the instrument moves at the
/// smallest separation this file uses.
fn oracle_sign(a: Rat, b: Rat, c: Rat) -> Orient {
    const P: usize = 512;
    let bf = |x: Rat| {
        BigFloat::from_i128(x.numer(), P).div(&BigFloat::from_i128(x.denom(), P), P, HP_RM)
    };
    let val = bf(a).add(&bf(b).mul(&bf(c).sqrt(P, HP_RM), P, HP_RM), P, HP_RM);
    bigfloat_sign(&val)
}

fn oracle_biquad(a: Rat, b: Rat, c: Rat, d: Rat, u: Rat, v: Rat) -> Orient {
    const P: usize = 512;
    let bf = |x: Rat| {
        BigFloat::from_i128(x.numer(), P).div(&BigFloat::from_i128(x.denom(), P), P, HP_RM)
    };
    let (su, sv) = (bf(u).sqrt(P, HP_RM), bf(v).sqrt(P, HP_RM));
    let suv = su.mul(&sv, P, HP_RM);
    let val = bf(a)
        .add(&bf(b).mul(&su, P, HP_RM), P, HP_RM)
        .add(&bf(c).mul(&sv, P, HP_RM), P, HP_RM)
        .add(&bf(d).mul(&suv, P, HP_RM), P, HP_RM);
    bigfloat_sign(&val)
}

fn bigfloat_sign(val: &BigFloat) -> Orient {
    const P: usize = 512;
    // ★ Never compare against a `from_f64(0.0)` zero: astro-float's comparison against zero
    // is a known upstream issue (hit in this project before; measured here again — 0.414
    // answered "< 0" against a from_f64 zero while `sign()` and a nonzero-eps comparison
    // stayed consistent). The sign is read from `sign()` directly.
    let eps = BigFloat::from_f64(1e-120, P);
    if val.abs().cmp(&eps).is_some_and(|s| s < 0) {
        return Orient::Zero;
    }
    match val.sign() {
        Some(astro_float::Sign::Neg) => Orient::Negative,
        _ => Orient::Positive,
    }
}

/// The spec's case ladder, arm by arm, with hand-derived answers — including the two arms the
/// plan reviews had to correct on paper (c = 0 first; the opposite-sign arm is a *product*).
#[test]
fn the_sign_cases_are_exhaustive() {
    let cases: [(Rat, Rat, Rat, Orient); 16] = [
        // c = 0 first — b√0 = 0 whatever b is (the spec trap: sign(b) here would be wrong).
        (ri(0), ri(5), ri(0), Orient::Zero),
        (ri(-3), ri(7), ri(0), Orient::Negative),
        (ri(2), ri(-9), ri(0), Orient::Positive),
        // b = 0 / a = 0.
        (ri(3), ri(0), ri(2), Orient::Positive),
        (ri(-3), ri(0), ri(2), Orient::Negative),
        (ri(0), ri(0), ri(2), Orient::Zero),
        (ri(0), ri(2), ri(2), Orient::Positive),
        (ri(0), ri(-2), ri(2), Orient::Negative),
        // Same sign.
        (ri(1), ri(1), ri(2), Orient::Positive),
        (ri(-1), ri(-1), ri(2), Orient::Negative),
        // Opposite signs — all four (sign(a) × relation) quadrants.
        (ri(2), ri(-1), ri(2), Orient::Positive), // 2 − √2 > 0, a² > b²c
        (ri(1), ri(-1), ri(2), Orient::Negative), // 1 − √2 < 0, a² < b²c
        (ri(-2), ri(1), ri(2), Orient::Negative), // −2 + √2 < 0
        (ri(-1), ri(1), ri(2), Orient::Positive), // −1 + √2 > 0
        // Exact zero needs a rational √c: perfect-square radicand.
        (ri(2), ri(-1), ri(4), Orient::Zero), // 2 − √4
        (ri(-2), ri(1), ri(4), Orient::Zero),
    ];
    for (a, b, c, want) in cases {
        assert_eq!(qv(a, b, c).sign(), want, "sign({a:?} + {b:?}√{c:?})");
        assert_eq!(
            oracle_sign(a, b, c),
            want,
            "oracle disagrees on {a:?} + {b:?}√{c:?}"
        );
    }
    // A rational radicand exercises the √(p/q) = √(pq)/q integerization: −1 + 2·√(1/4) = 0.
    assert_eq!(qv(ri(-1), ri(2), r(1, 4)).sign(), Orient::Zero);
}

/// One value, two spellings — within one radical (√4 beside a rational) and across two
/// (√8 = 2√2, which drives the biquad recursion through its `P² − v·Q² = 0` arm).
#[test]
fn two_spellings_of_one_value_cancel() {
    let x = qv(ri(1), ri(1), ri(4)); // 1 + √4 = 3, kept unsimplified by construction
    let y = QuadVal::from_rat(ri(3));
    assert_eq!(x.checked_sub(&y).expect("compatible").sign(), Orient::Zero);
    // √8 − 2√2 = 0: opposite-sign P/Q, and P² − vQ² = 8 − 2·4 = 0 exactly.
    assert_eq!(
        biquad_sign(ri(0), ri(1), ri(-2), ri(0), ri(8), ri(2)),
        Some(Orient::Zero)
    );
}

#[test]
fn biquad_hand_cases() {
    let cases: [(Rat, Rat, Rat, Rat, Orient); 5] = [
        // 1 + √2 − √3 ≈ 0.68.
        (ri(1), ri(1), ri(-1), ri(0), Orient::Positive),
        // −4 + √2 + √3 ≈ −0.85.
        (ri(-4), ri(1), ri(1), ri(0), Orient::Negative),
        // (√2−1)(√3−1) = 1 − √2 − √3 + √6 ≈ 0.30.
        (ri(1), ri(-1), ri(-1), ri(1), Orient::Positive),
        // Its negation.
        (ri(-1), ri(1), ri(1), ri(-1), Orient::Negative),
        // Pure √6 term.
        (ri(0), ri(0), ri(0), ri(-3), Orient::Negative),
    ];
    for (a, b, c, d, want) in cases {
        assert_eq!(
            biquad_sign(a, b, c, d, ri(2), ri(3)),
            Some(want),
            "sign({a:?} + {b:?}√2 + {c:?}√3 + {d:?}√6)"
        );
        assert_eq!(oracle_biquad(a, b, c, d, ri(2), ri(3)), want);
    }
    // Degenerate radicals collapse to the lower storeys.
    assert_eq!(
        biquad_sign(ri(-1), ri(1), ri(9), ri(9), ri(2), ri(0)),
        Some(Orient::Positive), // v = 0: just −1 + √2
    );
    assert_eq!(biquad_sign(ri(1), ri(1), ri(1), ri(1), ri(-1), ri(2)), None);
}

/// **The adversarial ladder is self-verifying**: Pell convergents `p/q` of √2 satisfy
/// `p² − 2q² = ±1`, so `√2 − p/q` has the sign `−sign(p² − 2q²)` **by construction** — no
/// magic literals — while the separation shrinks like `1/(2√2·q²)`, reaching ~1e−73 before
/// `p, q` leave `i128`. Every rung must agree with both the derived answer and the oracle.
#[test]
fn pell_convergents_walk_the_precision_ladder() {
    let (mut p, mut q) = (1i128, 1i128);
    let mut rungs = 0;
    loop {
        let big = num_bigint::BigInt::from;
        let pell = big(p) * big(p) - big(2) * big(q) * big(q);
        let want = match i32::try_from(pell).expect("±1") {
            -1 => Orient::Positive, // p² − 2q² = −1 → p/q < √2
            1 => Orient::Negative,  // p/q > √2
            _ => unreachable!("Pell residue is ±1"),
        };
        // sign(√2 − p/q) = sign(−p/q + 1·√2).
        let val = qv(r(-p, q), ri(1), ri(2));
        assert_eq!(val.sign(), want, "rung p/q = {p}/{q}");
        assert_eq!(
            oracle_sign(r(-p, q), ri(1), ri(2)),
            want,
            "oracle at {p}/{q}"
        );
        rungs += 1;
        match pell_next(p, q) {
            Some((np, nq)) => (p, q) = (np, nq),
            None => break,
        }
    }
    assert!(
        rungs > 80,
        "the ladder should run deep into i128, got {rungs}"
    );
}

/// The next √2 convergent, `None` once it leaves `i128`.
fn pell_next(p: i128, q: i128) -> Option<(i128, i128)> {
    Some((p.checked_add(q.checked_mul(2)?)?, p.checked_add(q)?))
}

/// The oracle itself resolves the smallest gap the ladder reaches — the negative control that
/// licenses trusting it everywhere else.
#[test]
fn the_oracle_actually_resolves_the_adversarial_gap() {
    // Walk to the deepest in-range convergent, then check the oracle answers ±, not Zero.
    let (mut p, mut q) = (1i128, 1i128);
    while let Some((np, nq)) = pell_next(p, q) {
        (p, q) = (np, nq);
    }
    let got = oracle_sign(r(-p, q), ri(1), ri(2));
    assert_ne!(
        got,
        Orient::Zero,
        "oracle blind at separation ~1/(2√2·q²), q = {q}"
    );
}

#[test]
fn arithmetic_agrees_with_the_oracle() {
    // (a₁+b₁√c)(a₂+b₂√c) and sums, spot-checked against the oracle over a small grid.
    for a1 in -2..=2i128 {
        for b1 in -2..=2i128 {
            for a2 in -2..=2i128 {
                for b2 in -2..=2i128 {
                    let (x, y) = (qv(ri(a1), ri(b1), ri(3)), qv(ri(a2), ri(b2), ri(3)));
                    let prod = x.checked_mul(&y).expect("small values");
                    let want = oracle_sign(prod.a(), prod.b(), prod.c());
                    assert_eq!(prod.sign(), want, "({a1}+{b1}√3)({a2}+{b2}√3)");
                    let sum = x.checked_add(&y).expect("small values");
                    assert_eq!(sum.sign(), oracle_sign(sum.a(), sum.b(), sum.c()));
                }
            }
        }
    }
}

#[test]
fn the_constructor_refuses_a_negative_radicand_and_normalizes_a_dead_radical() {
    assert!(QuadVal::new(ri(1), ri(1), ri(-2)).is_none());
    let dead = QuadVal::new(ri(5), ri(0), ri(7)).expect("b = 0 is fine");
    assert_eq!(dead.c(), ri(0), "b = 0 stores one spelling");
    let dead2 = QuadVal::new(ri(5), ri(7), ri(0)).expect("c = 0 is fine");
    assert_eq!(dead2.b(), ri(0));
}

#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "same-radical arithmetic")]
fn mixing_radicals_is_a_caller_bug() {
    let x = qv(ri(1), ri(1), ri(2));
    let y = qv(ri(1), ri(1), ri(3));
    let _ = x.checked_add(&y);
}

#[test]
fn to_f64_realizes_the_value() {
    let x = qv(ri(1), ri(2), ri(2)); // 1 + 2√2
    assert!((x.to_f64() - (1.0 + 2.0 * 2f64.sqrt())).abs() < 1e-15);
}

mod props {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Random rationals against the oracle — both storeys.
        #[test]
        fn sign1_matches_the_oracle(
            an in -1000i128..1000, ad in 1i128..100,
            bn in -1000i128..1000, bd in 1i128..100,
            c in 0i128..1000,
        ) {
            let (a, b, c) = (r(an, ad), r(bn, bd), ri(c));
            prop_assert_eq!(qv(a, b, c).sign(), oracle_sign(a, b, c));
        }

        #[test]
        fn biquad_matches_the_oracle(
            a in -50i128..50, b in -50i128..50, c in -50i128..50, d in -50i128..50,
            u in 0i128..60, v in 0i128..60,
        ) {
            let got = biquad_sign(ri(a), ri(b), ri(c), ri(d), ri(u), ri(v)).expect("radicands ≥ 0");
            prop_assert_eq!(got, oracle_biquad(ri(a), ri(b), ri(c), ri(d), ri(u), ri(v)));
        }
    }
}
