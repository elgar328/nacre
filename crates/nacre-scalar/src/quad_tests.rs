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
                assert_eq!(plane_side(&plane(1, 0, 0, 0), &line, sv), Orient::Zero);
                assert_eq!(plane_side(&plane(0, 0, 1, 0), &line, sv), Orient::Zero);
            }
            // And the membership planes y = ±1 pin which root is which.
            assert_eq!(
                plane_side(&plane(0, 1, 0, -1), &line, &s[0]),
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
                assert_eq!(plane_side(&p1, &line, sv), Orient::Zero);
                assert_eq!(plane_side(&p2, &line, sv), Orient::Zero);
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

/// The radial-side sign agrees with the f64 distance oracle across the surface — including the
/// exact-zero shell, which only rational points on radius-r circles can witness.
#[test]
fn the_radial_side_knows_its_shell() {
    let (o, m, radius) = zcyl();
    // (3,4)/5-scaled points: exactly on the unit circle, inside, outside — rational witnesses.
    let cases: [([Rat; 3], Orient); 5] = [
        ([r(3, 5), r(4, 5), ri(7)], Orient::Zero),
        ([r(3, 5), r(4, 5), ri(-2)], Orient::Zero), // z is radially irrelevant
        ([r(1, 2), ri(0), ri(0)], Orient::Negative),
        ([ri(2), ri(0), ri(3)], Orient::Positive),
        ([ri(0), ri(0), ri(5)], Orient::Negative), // the axis itself
    ];
    for (p, want) in cases {
        assert_eq!(
            cylinder_radial_side(&p, &o, &m, radius),
            want,
            "point {p:?}"
        );
    }
    // A tilted rational axis: cross-check against the f64 realization.
    let o2 = [ri(1), r(-1, 2), ri(2)];
    let m2 = [ri(1), ri(2), ri(2)];
    let r2 = r(3, 2);
    for p in [
        [ri(1), ri(1), ri(1)],
        [ri(3), ri(3), ri(3)],
        [ri(1), r(-1, 2), ri(2)],
        [r(5, 2), r(5, 2), ri(5)],
    ] {
        let got = cylinder_radial_side(&p, &o2, &m2, r2);
        let pf: Vec<f64> = p.iter().map(|x| x.to_f64()).collect();
        let w = [pf[0] - 1.0, pf[1] + 0.5, pf[2] - 2.0];
        let mf = [1.0, 2.0, 2.0];
        let dot = |x: &[f64; 3], y: &[f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
        let val = dot(&w, &w) * dot(&mf, &mf) - dot(&w, &mf).powi(2) - 2.25 * dot(&mf, &mf);
        let want = if val.abs() < 1e-9 {
            Orient::Zero
        } else if val < 0.0 {
            Orient::Negative
        } else {
            Orient::Positive
        };
        assert_eq!(got, want, "point {pf:?}");
    }
}

/// **The radial side answers where the checked spelling had to give up.** A point and an axis
/// whose coordinates carry ~10²³ denominators: `|w|²|m|²` alone leaves `i128`, so the retired
/// `Option` version returned `None` and its callers rejected for the arithmetic. The total one
/// answers, and the answer is checked against the geometry by hand — the point sits at radius
/// `2·r` from the axis, so it is outside; halving that puts it inside.
#[test]
fn the_radial_side_answers_where_the_checked_one_could_not() {
    let wide = |x: f64| Rat::from_decimal(x).expect("in the decimal window");
    // A sub-micron cylinder on ẑ: radius 5e-8, axis through (1e-7, 1e-7).
    let o = [
        wide(1.0000000000000002e-7),
        wide(1.0000000000000002e-7),
        ri(0),
    ];
    let m = [ri(0), ri(0), ri(1)];
    let radius = wide(5.000000000000001e-8);
    // The retired spelling, verbatim — it cannot even form `|w|²`.
    let checked = |p: &[Rat; 3]| -> Option<Orient> {
        let w = [
            p[0].checked_sub(o[0])?,
            p[1].checked_sub(o[1])?,
            p[2].checked_sub(o[2])?,
        ];
        let dot = |x: &[Rat; 3], y: &[Rat; 3]| {
            x[0].checked_mul(y[0])?
                .checked_add(x[1].checked_mul(y[1])?)?
                .checked_add(x[2].checked_mul(y[2])?)
        };
        let mm = dot(&m, &m)?;
        let wm = dot(&w, &m)?;
        let val = dot(&w, &w)?
            .checked_mul(mm)?
            .checked_sub(wm.checked_mul(wm)?)?
            .checked_sub(radius.checked_mul(radius)?.checked_mul(mm)?)?;
        Some(match val.cmp(&Rat::from_int(0)) {
            core::cmp::Ordering::Less => Orient::Negative,
            core::cmp::Ordering::Equal => Orient::Zero,
            core::cmp::Ordering::Greater => Orient::Positive,
        })
    };
    let outside = [
        o[0].checked_add(radius)
            .and_then(|x| x.checked_add(radius))
            .expect("small sum"),
        o[1],
        ri(3),
    ];
    let inside = [
        o[0].checked_add(Rat::new(radius.numer(), radius.denom() * 2).expect("half r"))
            .expect("small sum"),
        o[1],
        ri(-4),
    ];
    assert_eq!(
        checked(&outside),
        None,
        "the fixture must actually overflow the old road"
    );
    assert_eq!(checked(&inside), None, "…on both witnesses");
    assert_eq!(
        cylinder_radial_side(&outside, &o, &m, radius),
        Orient::Positive
    );
    assert_eq!(
        cylinder_radial_side(&inside, &o, &m, radius),
        Orient::Negative
    );
    // And exactly on the shell — the case only exact arithmetic can witness.
    let on = [o[0].checked_add(radius).expect("small sum"), o[1], ri(9)];
    assert_eq!(cylinder_radial_side(&on, &o, &m, radius), Orient::Zero);
}

/// Two **parallel** axes: clear when their distance exceeds `r₁ + r₂`, tangent at equality, and
/// overlapping below it — **with the sum formed inside the integer arithmetic**. The radii here
/// are wide enough that `r₁.checked_add(r₂)` still works but `(r₁+r₂)²·|m|²` does not, which is
/// what the retired gate computed.
#[test]
fn parallel_axes_clear_compares_against_the_radius_sum() {
    let wide = |x: f64| Rat::from_decimal(x).expect("in the decimal window");
    let m = [ri(0), ri(0), ri(1)];
    let (ra, rb) = (wide(5.000000000000001e-8), wide(3.000000000000001e-8));
    let o_a = [ri(0), ri(0), ri(0)];
    let far = [wide(1.0000000000000002e-7), ri(0), ri(5)];
    // Short on purpose: `near` only has to sit inside the radius sum — the width that matters
    // is in `ra`, `rb` and `far`, which is where the old arithmetic gave out.
    let near = [wide(7e-8), ri(0), ri(5)];
    assert_eq!(
        crate::cylinders_clear(&o_a, &m, ra, &far, &m, rb),
        Orient::Positive,
        "1e-7 apart, radii summing to 8e-8: clear"
    );
    assert_eq!(
        crate::cylinders_clear(&o_a, &m, ra, &near, &m, rb),
        Orient::Negative,
        "7e-8 apart, radii summing to 8e-8: overlapping"
    );
    let touching = [ra.checked_add(rb).expect("small sum"), ri(0), ri(5)];
    assert_eq!(
        crate::cylinders_clear(&o_a, &m, ra, &touching, &m, rb),
        Orient::Zero,
        "exactly tangent"
    );
}

/// ★★ **Axes that are not parallel are the case the caller used to refuse outright.** The
/// distance between two skew lines is their common perpendicular, and the same comparison
/// against `r₁ + r₂` decides the pair — a drill crossing a bore six apart with radii summing to
/// four does not touch it, whatever the angle between them.
///
/// Tangency is here because it is the boundary a `≥` would swallow: at exactly `r₁ + r₂` the two
/// surfaces meet along a line, which is contact, not clearance.
#[test]
fn skew_axes_are_compared_by_their_common_perpendicular() {
    let z = [ri(0), ri(0), ri(1)];
    let y = [ri(0), ri(1), ri(0)];
    let o_a = [ri(0), ri(0), ri(0)]; // the vertical bore
    let (ra, rb) = (ri(3), ri(1));
    // The user's own model: a horizontal tunnel six away from a vertical bore.
    assert_eq!(
        crate::cylinders_clear(&o_a, &z, ra, &[ri(-6), ri(0), ri(0)], &y, rb),
        Orient::Positive,
        "six apart, radii summing to four: clear"
    );
    assert_eq!(
        crate::cylinders_clear(&o_a, &z, ra, &[ri(-4), ri(0), ri(0)], &y, rb),
        Orient::Zero,
        "exactly four apart: tangent, which is contact"
    );
    assert_eq!(
        crate::cylinders_clear(&o_a, &z, ra, &[ri(-2), ri(0), ri(0)], &y, rb),
        Orient::Negative,
        "two apart: the surfaces cut through each other"
    );
    // Axes that actually meet are distance zero, whatever the radii.
    assert_eq!(
        crate::cylinders_clear(&o_a, &z, ra, &[ri(0), ri(0), ri(0)], &y, rb),
        Orient::Negative,
        "crossing axes cannot be clear"
    );
    // The answer does not depend on how long the direction vectors are, nor on where along its
    // own axis each cylinder is measured from.
    let long_y = [ri(0), ri(17), ri(0)];
    assert_eq!(
        crate::cylinders_clear(&o_a, &z, ra, &[ri(-6), ri(40), ri(0)], &long_y, rb),
        Orient::Positive,
        "scale and axial offset change nothing"
    );
}

// ---- a segment against a cylinder's axis (M6-2b preparation) ----
//
// The running fixture: the axis is the vertical line through the origin with `r = 3`, and every
// segment below lies in the plane `z = 0`, which is perpendicular to it.

fn seg_meets(p0: [i128; 3], p1: [i128; 3], radius: i128) -> bool {
    crate::segment_meets_cylinder(
        &p0.map(ri),
        &p1.map(ri),
        &[ri(0), ri(0), ri(0)],
        &[ri(0), ri(0), ri(1)],
        ri(radius),
    )
}

#[test]
fn an_endpoint_inside_the_cylinder_meets_it() {
    assert!(seg_meets([0, 0, 0], [10, 0, 0], 3), "starts on the axis");
    assert!(seg_meets([2, 0, 0], [10, 0, 0], 3), "starts inside");
    assert!(
        !seg_meets([4, 0, 0], [10, 0, 0], 3),
        "starts outside, runs away"
    );
}

/// A **chord**: both ends outside, the middle through the disk.
#[test]
fn a_chord_meets_though_both_ends_are_outside() {
    assert!(seg_meets([-10, 1, 0], [10, 1, 0], 3));
}

/// ★★ **The case the whole check turns on.** The line passes within `r`, but the segment stops
/// before it gets there — the perpendicular foot lies outside `[0,1]`. Ten edges of today's
/// corpus are exactly this shape (a boss standing far away whose wall plane, extended, crosses a
/// bore), and answering "meets" here would close the family opened for them.
#[test]
fn a_near_line_whose_foot_is_off_the_segment_does_not_meet() {
    // The line y = 1 passes 1 from the axis; this piece of it lives at x ∈ [20, 30].
    assert!(!seg_meets([20, 1, 0], [30, 1, 0], 3));
    // …and the piece that does reach the foot meets it.
    assert!(seg_meets([-1, 1, 0], [30, 1, 0], 3));
}

/// Exact tangency: the closest approach is `r` itself, and the closed disk includes it.
#[test]
fn a_segment_grazing_at_exactly_r_meets_it() {
    assert!(seg_meets([-10, 3, 0], [10, 3, 0], 3));
    assert!(!seg_meets([-10, 4, 0], [10, 4, 0], 3));
}

/// The axis' direction magnitude and the points' denominators are not part of the question.
#[test]
fn the_answer_does_not_depend_on_scale() {
    let long_axis = [ri(0), ri(0), ri(7)];
    let plain = crate::segment_meets_cylinder(
        &[r(-10, 1), r(1, 1), ri(0)],
        &[r(10, 1), r(1, 1), ri(0)],
        &[ri(0), ri(0), ri(0)],
        &[ri(0), ri(0), ri(1)],
        ri(3),
    );
    let scaled = crate::segment_meets_cylinder(
        &[r(-20, 2), r(2, 2), ri(0)],
        &[r(20, 2), r(2, 2), ri(0)],
        &[ri(0), ri(0), ri(0)],
        &long_axis,
        ri(3),
    );
    assert_eq!(plain, scaled);
}

/// ★ **The negative control.** The same segment, the same radius — move the axis and the answer
/// must change. Without it every assertion above could be passing for a reason that has nothing
/// to do with the cylinder.
#[test]
fn moving_the_axis_moves_the_answer() {
    let seg = ([ri(-10), ri(1), ri(0)], [ri(10), ri(1), ri(0)]);
    let z = [ri(0), ri(0), ri(1)];
    assert!(crate::segment_meets_cylinder(
        &seg.0,
        &seg.1,
        &[ri(0), ri(0), ri(0)],
        &z,
        ri(3)
    ));
    assert!(!crate::segment_meets_cylinder(
        &seg.0,
        &seg.1,
        &[ri(0), ri(40), ri(0)],
        &z,
        ri(3)
    ));
}

// ---- One ruler for two kinds of point (M6-2b preparation) ----
//
// ★★ An arrangement that carries arcs holds vertices of two shapes: a three-plane node, whose
// coordinate is rational, and a plane·plane·cylinder node, whose coordinate is `a + b√c`. The
// winding is read at the ring's lexicographically least node, so the two must be comparable —
// and the comparison must never decline, because a "cannot order" would leave a valid solid
// unbuilt for a reason about arithmetic rather than about shape.

/// The unit z-cylinder cut by `x = 0` gives `(0, ±1, z)` — **rational** roots, so the same two
/// points can be named the other way too, by three planes. ★ That is the only place the new road
/// and the old one answer about the *same value*, which makes it the one independent oracle
/// available: everything else about `a + b√c` has no three-plane spelling at all.
#[test]
fn the_two_roads_agree_where_a_crossing_is_rational() {
    let (o, m, radius) = zcyl();
    let cut = plane(1, 0, 0, 0); // x = 0
    let lid = plane(0, 0, 1, -3); // z = 3
    let CylinderMeet::Pair { line, s } =
        plane_plane_cylinder(&cut, &lid, &o, &m, radius).expect("fits")
    else {
        panic!("x = 0 crosses the unit cylinder twice");
    };
    // ★ Which root is which is **derived, not read back**: `plane_plane_cylinder` orders along
    // `d = n_cut × n_lid = (1,0,0) × (0,0,1) = (0,-1,0)`, so ascending `s` runs toward −y and
    // `s[0]` is the `+y` end.
    for (sv, want_y) in s.iter().zip([1.0, -1.0]) {
        let p = branch_point_f64(&line, sv);
        assert!(
            (p[0]).abs() < 1e-12 && (p[1] - want_y).abs() < 1e-12 && (p[2] - 3.0).abs() < 1e-12,
            "{p:?}"
        );
    }
    // The same point, named by three planes: x = 0, z = 3, y = ±1.
    for (sv, y) in s.iter().zip([1i128, -1]) {
        let meet = crate::three_planes_big([
            &crate::PlaneName::Narrow(cut),
            &crate::PlaneName::Narrow(lid),
            &crate::PlaneName::Narrow(plane(0, 1, 0, -y)),
        ])
        .expect("three independent planes");
        for axis in 0..3 {
            assert_eq!(
                cmp_coord_meet_branch(&meet, &line, sv, axis),
                Orient::Zero,
                "the two roads name one point (axis {axis})"
            );
        }
    }
    // ★★ **Coincidence alone cannot see a reversed subtraction** — `Zero` is symmetric. So the
    // same cross-road fixture also asserts a *direction*: the `y = +1` three-plane point against
    // the `y = −1` root, where the answer has a sign to get wrong.
    let top = crate::three_planes_big([
        &crate::PlaneName::Narrow(cut),
        &crate::PlaneName::Narrow(lid),
        &crate::PlaneName::Narrow(plane(0, 1, 0, -1)),
    ])
    .expect("independent");
    assert_eq!(
        cmp_coord_meet_branch(&top, &line, &s[1], 1),
        Orient::Positive,
        "+1 lies above the −1 root"
    );
}

/// **Hand-derived, from the geometry and not from the engine.** `x = 1/2` cuts the unit circle at
/// `y = ±√3/2`, so on the `y` axis the second root is above the first and both straddle any
/// rational between them. The three-plane point `(1/2, 0, 0)` sits between them.
#[test]
fn a_rational_point_is_placed_between_two_irrational_ones() {
    let (o, m, radius) = zcyl();
    let cut = plane(2, 0, 0, -1); // x = 1/2
    let CylinderMeet::Pair { line, s } =
        plane_plane_cylinder(&cut, &plane(0, 0, 1, 0), &o, &m, radius).expect("fits")
    else {
        panic!("x = 1/2 crosses twice");
    };
    let mid = crate::three_planes_big([
        &crate::PlaneName::Narrow(cut),
        &crate::PlaneName::Narrow(plane(0, 0, 1, 0)),
        &crate::PlaneName::Narrow(plane(0, 1, 0, 0)), // y = 0
    ])
    .expect("independent");
    // ★ Which root is which is **derived from the two normals**, not read back: the roots ascend
    // along `d = n_cut × n_lid = (2,0,0) × (0,0,1) = (0,−2,0)`, so `s[0]` is the `+√3/2` end.
    // Hence `y`: `s[1] = −√3/2  <  mid = 0  <  s[0] = +√3/2`.
    assert_eq!(
        cmp_coord_meet_branch(&mid, &line, &s[0], 1),
        Orient::Negative
    );
    assert_eq!(
        cmp_coord_meet_branch(&mid, &line, &s[1], 1),
        Orient::Positive
    );
    // The two roots against each other, on the same radical — first storey.
    assert_eq!(
        cmp_coord_branch((&line, &s[0]), (&line, &s[1]), 1),
        Orient::Positive
    );
    assert_eq!(
        cmp_coord_branch((&line, &s[1]), (&line, &s[0]), 1),
        Orient::Negative
    );
    // x is the same for both — the cut plane pins it.
    assert_eq!(
        cmp_coord_branch((&line, &s[0]), (&line, &s[1]), 0),
        Orient::Zero
    );
    assert_eq!(cmp_coord_meet_branch(&mid, &line, &s[0], 0), Orient::Zero);
}

/// **Two different radicands — the second storey.** `x = 1/2` gives `y = ±√3/2`; `x = 3/5` gives
/// `y = ±4/5`… so pick `x = 1/3`, whose `y = ±√8/3` shares no radicand with `√3/2`. The exact
/// order must match the realized one, and the realization is only the *witness* here: the
/// assertion is the hand-derived inequality `√8/3 > √3/2` (`8/9 > 3/4`).
#[test]
fn two_points_from_different_cuts_order_on_the_second_storey() {
    let (o, m, radius) = zcyl();
    let lid = plane(0, 0, 1, 0);
    let roots = |a: i128, b: i128| {
        let cut = plane(b, 0, 0, -a); // x = a/b
        let CylinderMeet::Pair { line, s } =
            plane_plane_cylinder(&cut, &lid, &o, &m, radius).expect("fits")
        else {
            panic!("crosses twice");
        };
        (line, s)
    };
    let (l_half, s_half) = roots(1, 2); // y = ±√3/2 ≈ ±0.8660
    let (l_third, s_third) = roots(1, 3); // y = ±√8/3 ≈ ±0.9428
    assert_ne!(
        s_half[1].c(),
        s_third[1].c(),
        "the fixture must really reach the second storey"
    );
    // ★ Both lines order toward −y (`d = (b,0,0) × (0,0,1) = (0,−b,0)`), so `s[0]` is the upper
    // root. √8/3 > √3/2 because 8/9 > 3/4 — derived from the radii, not read back.
    assert_eq!(
        cmp_coord_branch((&l_third, &s_third[0]), (&l_half, &s_half[0]), 1),
        Orient::Positive
    );
    assert_eq!(
        cmp_coord_branch((&l_half, &s_half[1]), (&l_third, &s_third[1]), 1),
        Orient::Positive,
        "mirrored below the axis, −√3/2 is the larger of the two"
    );
}

/// ★ **A point too wide for `Rat` still answers.** Three planes whose meet needs a numerator past
/// `i128` give a [`MeetPoint::Wide`]; the road that narrows to `Rat` would have declined here, and
/// declining is what this whole design exists to avoid.
#[test]
fn a_wide_three_plane_point_is_still_placed() {
    use num_bigint::BigInt;
    let (o, m, radius) = zcyl();
    let cut = plane(1, 0, 0, 0);
    let lid = plane(0, 0, 1, 0);
    let CylinderMeet::Pair { line, s } =
        plane_plane_cylinder(&cut, &lid, &o, &m, radius).expect("fits")
    else {
        panic!("crosses twice");
    };
    // x = 10^40 · t / 10^40 — a wide name for the plane x = 7, so the meet is wide by carrier.
    let huge: BigInt = BigInt::from(10u8).pow(40);
    let wide_x = crate::PlaneName::Wide([
        huge.clone(),
        BigInt::from(0),
        BigInt::from(0),
        -(&huge * BigInt::from(7)),
    ]);
    let meet = crate::three_planes_big([
        &wide_x,
        &crate::PlaneName::Narrow(lid),
        &crate::PlaneName::Narrow(plane(0, 1, 0, 0)),
    ])
    .expect("independent");
    // The three-plane point is `(7, 0, 0)`; the branch points are `(0, +1, 0)` and `(0, −1, 0)`
    // (ascending `s` runs toward −y, as derived above).
    assert_eq!(
        cmp_coord_meet_branch(&meet, &line, &s[0], 0),
        Orient::Positive
    );
    assert_eq!(
        cmp_coord_meet_branch(&meet, &line, &s[0], 1),
        Orient::Negative
    );
    assert_eq!(
        cmp_coord_meet_branch(&meet, &line, &s[1], 1),
        Orient::Positive
    );
}

/// The degenerate rung: an axis the line does not move along (`dir_k = 0`) makes the coordinate
/// rational, so `b` vanishes and the comparison drops to a plain integer sign.
#[test]
fn an_axis_the_line_does_not_move_along_is_rational() {
    let (o, m, radius) = zcyl();
    let cut = plane(2, 0, 0, -1); // x = 1/2 — the line runs in y, not x or z
    let lid = plane(0, 0, 1, -5); // z = 5
    let CylinderMeet::Pair { line, s } =
        plane_plane_cylinder(&cut, &lid, &o, &m, radius).expect("fits")
    else {
        panic!("crosses twice");
    };
    for sv in &s {
        let p = branch_point_f64(&line, sv);
        assert!(
            (p[0] - 0.5).abs() < 1e-12 && (p[2] - 5.0).abs() < 1e-12,
            "{p:?}"
        );
    }
    // Both roots share x and z exactly, whatever the radical does to y.
    assert_eq!(
        cmp_coord_branch((&line, &s[0]), (&line, &s[1]), 0),
        Orient::Zero
    );
    assert_eq!(
        cmp_coord_branch((&line, &s[0]), (&line, &s[1]), 2),
        Orient::Zero
    );
}
