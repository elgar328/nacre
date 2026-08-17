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
fn oracle_val(a: Rat, b: Rat, c: Rat) -> BigFloat {
    const P: usize = 512;
    let bf = |x: Rat| {
        BigFloat::from_i128(x.numer(), P).div(&BigFloat::from_i128(x.denom(), P), P, HP_RM)
    };
    bf(a).add(&bf(b).mul(&bf(c).sqrt(P, HP_RM), P, HP_RM), P, HP_RM)
}

fn oracle_sign(a: Rat, b: Rat, c: Rat) -> Orient {
    bigfloat_sign(&oracle_val(a, b, c))
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
    // ★ The differential runs against the **inputs' realizations** — the first spelling fed
    // the result's own coefficients back to the oracle, which verifies `sign()` but lets a
    // wrong `checked_mul` coefficient pass unseen (the shared-derivation trap): the oracle
    // must derive the product from x and y, never from the thing under test.
    const P: usize = 512;
    for a1 in -2..=2i128 {
        for b1 in -2..=2i128 {
            for a2 in -2..=2i128 {
                for b2 in -2..=2i128 {
                    let (x, y) = (qv(ri(a1), ri(b1), ri(3)), qv(ri(a2), ri(b2), ri(3)));
                    let (xv, yv) = (
                        oracle_val(ri(a1), ri(b1), ri(3)),
                        oracle_val(ri(a2), ri(b2), ri(3)),
                    );
                    let prod = x.checked_mul(&y).expect("small values");
                    assert_eq!(
                        prod.sign(),
                        bigfloat_sign(&xv.mul(&yv, P, HP_RM)),
                        "({a1}+{b1}√3)({a2}+{b2}√3)"
                    );
                    let sum = x.checked_add(&y).expect("small values");
                    assert_eq!(sum.sign(), bigfloat_sign(&xv.add(&yv, P, HP_RM)));
                    let negated = x.checked_neg().expect("small values");
                    assert_eq!(
                        negated.sign(),
                        bigfloat_sign(&xv.mul(&BigFloat::from_i128(-1, P), P, HP_RM))
                    );
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

// ──────────────────── plane · plane · cylinder (commit 2) ────────────────────

/// The unit z-cylinder at the origin, radius 1 — the degenerate ladder's stage.
fn zcyl() -> ([Rat; 3], [Rat; 3], Rat) {
    ([ri(0), ri(0), ri(0)], [ri(0), ri(0), ri(1)], ri(1))
}

fn plane(a: i128, b: i128, c: i128, d: i128) -> [Rat; 4] {
    [ri(a), ri(b), ri(c), ri(d)]
}

/// Every rung of the degenerate ladder answers with its own name — no silent fallback.
#[test]
fn the_degenerate_ladder_names_every_outcome() {
    let (o, m, r1) = zcyl();
    let meet = |p1: &[Rat; 4], p2: &[Rat; 4]| {
        plane_plane_cylinder(p1, p2, &o, &m, r1).expect("no overflow in fixtures")
    };
    // One plane, two spellings (x = 1 and 2x = 2).
    assert!(matches!(
        meet(&plane(1, 0, 0, -1), &plane(2, 0, 0, -2)),
        CylinderMeet::CoincidentPlanes
    ));
    assert!(matches!(
        meet(&plane(1, 0, 0, -1), &plane(1, 0, 0, -2)),
        CylinderMeet::ParallelPlanes
    ));
    // x = 1 ∧ y = 0: the vertical line ON the surface — a ruling, not points.
    assert!(matches!(
        meet(&plane(1, 0, 0, -1), &plane(0, 1, 0, 0)),
        CylinderMeet::OnRuling(_)
    ));
    // x = 2 ∧ y = 0: vertical, off the surface.
    assert!(matches!(
        meet(&plane(1, 0, 0, -2), &plane(0, 1, 0, 0)),
        CylinderMeet::AxisParallelMiss(_)
    ));
    // x = 1 ∧ z = 0: horizontal, touching at (1, 0, 0) — a rational double root.
    match meet(&plane(1, 0, 0, -1), &plane(0, 0, 1, 0)) {
        CylinderMeet::Tangent { line, s } => {
            let p = line.point_f64(&QuadVal::from_rat(s));
            assert_eq!(p, [1.0, 0.0, 0.0], "the tangency point");
        }
        other => panic!("expected Tangent, got {other:?}"),
    }
    // x = 0 ∧ z = 0: through the axis — two points (0, ±1, 0), s ascending.
    match meet(&plane(1, 0, 0, 0), &plane(0, 0, 1, 0)) {
        CylinderMeet::Pair { line, s } => {
            let (lo, hi) = (line.point_f64(&s[0]), line.point_f64(&s[1]));
            assert_eq!(lo[1].abs(), 1.0);
            assert_eq!(hi[1], -lo[1]);
            // Both roots lie on both cutting planes — the zero cross-check for plane_side,
            // meaningful precisely because the general root is irrational.
            for sv in &s {
                assert_eq!(
                    plane_side(&plane(1, 0, 0, 0), &line, sv),
                    Some(Orient::Zero)
                );
                assert_eq!(
                    plane_side(&plane(0, 0, 1, 0), &line, sv),
                    Some(Orient::Zero)
                );
            }
            // And the membership planes y = ±1 pin which root is which.
            assert_eq!(
                plane_side(&plane(0, 1, 0, -1), &line, &s[0]).unwrap(),
                if lo[1] > 0.0 {
                    Orient::Zero
                } else {
                    Orient::Negative
                }
            );
        }
        other => panic!("expected Pair, got {other:?}"),
    }
    // x = 2 ∧ z = 0: crosses nothing.
    assert!(matches!(
        meet(&plane(1, 0, 0, -2), &plane(0, 0, 1, 0)),
        CylinderMeet::Miss(_)
    ));
}

/// A tilted axis: the roots are genuinely irrational; their realizations must sit on the
/// cylinder and on both cutting planes.
#[test]
fn a_tilted_meet_realizes_onto_its_carriers() {
    let o = [ri(1), r(-1, 2), ri(2)];
    let m = [ri(1), ri(2), ri(2)];
    let radius = r(3, 2);
    let p1 = plane(1, 0, 0, -1); // x = 1 (contains the base point)
    let p2 = plane(0, 0, 1, -2); // z = 2
    match plane_plane_cylinder(&p1, &p2, &o, &m, radius).expect("fits") {
        CylinderMeet::Pair { line, s } => {
            for sv in &s {
                let p = line.point_f64(sv);
                // Plane residuals.
                assert!((p[0] - 1.0).abs() < 1e-12 && (p[2] - 2.0).abs() < 1e-12);
                // Cylinder residual: |w|²|m|² − (w·m)² = r²|m|².
                let w = [p[0] - 1.0, p[1] + 0.5, p[2] - 2.0];
                let mf = [1.0, 2.0, 2.0];
                let dot = |x: &[f64; 3], y: &[f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
                let lhs = dot(&w, &w) * dot(&mf, &mf) - dot(&w, &mf).powi(2);
                let rhs = 1.5f64.powi(2) * dot(&mf, &mf);
                assert!((lhs - rhs).abs() < 1e-9, "off the cylinder: {lhs} vs {rhs}");
                // plane_side against its own carriers is exactly zero.
                assert_eq!(plane_side(&p1, &line, sv), Some(Orient::Zero));
                assert_eq!(plane_side(&p2, &line, sv), Some(Orient::Zero));
            }
        }
        other => panic!("expected Pair, got {other:?}"),
    }
}

/// Cut the unit z-cylinder (seam +x) by x = c planes and sort every root with the exact
/// comparator: the order must equal the f64 angle order wherever f64 can see it.
#[test]
fn the_circular_order_matches_the_realized_angles() {
    let (o, m, radius) = zcyl();
    let ref_dir = [ri(1), ri(0), ri(0)]; // seam at +x, θ right-handed about +z
    let mut points: Vec<(MeetLine, QuadVal)> = Vec::new();
    for (cn, cd) in [(-9, 10), (-1, 2), (0, 1), (1, 2), (9, 10)] {
        let cut = [ri(1), ri(0), ri(0), r(-cn, cd)]; // x = cn/cd
        match plane_plane_cylinder(&cut, &plane(0, 0, 1, 0), &o, &m, radius).expect("fits") {
            CylinderMeet::Pair { line, s } => {
                for sv in s {
                    points.push((line.clone(), sv));
                }
            }
            other => panic!("expected Pair, got {other:?}"),
        }
    }
    let theta = |p: [f64; 3]| -> f64 {
        let t = p[1].atan2(p[0]);
        if t <= 0.0 {
            t + std::f64::consts::TAU
        } else {
            t
        }
    };
    let order = |a: &(MeetLine, QuadVal), b: &(MeetLine, QuadVal)| -> std::cmp::Ordering {
        match circular_order_about_seam(&o, &m, &ref_dir, (&a.0, &a.1), (&b.0, &b.1))
            .expect("no overflow")
        {
            SeamOrder::Ordered(ord) => ord,
            SeamOrder::SeamIncident { .. } => panic!("no fixture point sits on the seam"),
        }
    };
    let mut by_exact = points.clone();
    by_exact.sort_by(order);
    let mut by_angle = points.clone();
    by_angle.sort_by(|a, b| {
        theta(a.0.point_f64(&a.1))
            .partial_cmp(&theta(b.0.point_f64(&b.1)))
            .unwrap()
    });
    let realized: Vec<[f64; 3]> = by_exact.iter().map(|(l, s)| l.point_f64(s)).collect();
    let realized_angle: Vec<[f64; 3]> = by_angle.iter().map(|(l, s)| l.point_f64(s)).collect();
    assert_eq!(
        realized, realized_angle,
        "exact order ≠ realized angle order"
    );
}

/// A point on the seam generator is answered by name, never ranked — and θ = π is a real
/// class of its own (the y = 0 cut hits both at once).
#[test]
fn the_seam_point_is_surfaced_not_ranked() {
    let (o, m, radius) = zcyl();
    let ref_dir = [ri(1), ri(0), ri(0)];
    // y = 0 ∧ z = 0: the diameter through the seam — points (±1, 0, 0).
    let meet =
        plane_plane_cylinder(&plane(0, 1, 0, 0), &plane(0, 0, 1, 0), &o, &m, radius).expect("fits");
    let CylinderMeet::Pair { line, s } = meet else {
        panic!("expected Pair, got {meet:?}");
    };
    // Which root is +x?
    let (seam_s, pi_s) = if line.point_f64(&s[0])[0] > 0.0 {
        (&s[0], &s[1])
    } else {
        (&s[1], &s[0])
    };
    let got = circular_order_about_seam(&o, &m, &ref_dir, (&line, seam_s), (&line, pi_s))
        .expect("no overflow");
    assert_eq!(
        got,
        SeamOrder::SeamIncident {
            first: true,
            second: false
        },
        "θ = 0 must come back by name; θ = π is an ordinary class"
    );
}

/// Two cuts 1e−18 apart: their upper roots realize to the **same** f64 angle (the instrument
/// f64 is blind), and the exact comparator still gives the strict order — the reason this
/// tower exists.
#[test]
fn the_exact_order_outresolves_f64() {
    let (o, m, radius) = zcyl();
    let ref_dir = [ri(1), ri(0), ri(0)];
    let z0 = plane(0, 0, 1, 0);
    let c1 = r(1, 2);
    let c2 = r(500_000_000_000_000_001, 1_000_000_000_000_000_000); // 1/2 + 1e−18
    let upper = |c: Rat| -> (MeetLine, QuadVal) {
        let cut = [ri(1), ri(0), ri(0), ri(0).checked_sub(c).unwrap()];
        match plane_plane_cylinder(&cut, &z0, &o, &m, radius).expect("fits") {
            CylinderMeet::Pair { line, s } => {
                let sv = if line.point_f64(&s[0])[1] > 0.0 {
                    s[0]
                } else {
                    s[1]
                };
                (line, sv)
            }
            other => panic!("expected Pair, got {other:?}"),
        }
    };
    let (l1, s1) = upper(c1);
    let (l2, s2) = upper(c2);
    // f64 is blind: both x-coordinates realize to the same float.
    assert_eq!(
        l1.point_f64(&s1)[0].to_bits(),
        l2.point_f64(&s2)[0].to_bits(),
        "the fixture must sit below f64 resolution"
    );
    // Exact is not: larger x ⇒ smaller θ on the upper half, strictly.
    let got =
        circular_order_about_seam(&o, &m, &ref_dir, (&l1, &s1), (&l2, &s2)).expect("no overflow");
    assert_eq!(got, SeamOrder::Ordered(std::cmp::Ordering::Greater));
}
