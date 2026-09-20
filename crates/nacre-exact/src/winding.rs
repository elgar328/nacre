use super::*;
/// The exact square root of a non-negative rational, or `None` when it is irrational.
///
/// A fraction in lowest terms is a perfect square exactly when its numerator and denominator both
/// are — they share no factor to trade — so this is two integer square roots and two checks.
/// `√v` when it is rational — numerator and denominator both perfect squares — else `None`
/// (also for `v < 0`). The sketch's arcs ask this of `|start − centre|²`: a radius the kernel can
/// state is a rational one, and this is the only place that question is answered.
/// **The sign of `a + b·π`** for integers `a`, `b`, exactly. π is irrational, so the sum is zero
/// only when both are; otherwise `f(x) = a + b·x` is monotone and its sign at π is the sign it has
/// at both ends of a rational bracket around π — here π's decimal expansion to 37 places. If the
/// two ends disagree, π sits where `f` changes sign and the answer is `None`: undecidable at this
/// width, never a guess.
///
/// The consumer is the prism builder's winding question: a ring whose arcs are quarter turns has
/// twice its signed area equal to such a sum (`nacre-ops`, `Ring2d::winding_sign`).
pub fn sign_a_plus_b_pi_int(a: &num_bigint::BigInt, b: &num_bigint::BigInt) -> Option<Orient> {
    use num_bigint::BigInt;
    use num_traits::{Signed, Zero};
    let sign = |x: &BigInt| {
        if x.is_zero() {
            Orient::Zero
        } else if x.is_positive() {
            Orient::Positive
        } else {
            Orient::Negative
        }
    };
    if b.is_zero() {
        return Some(sign(a));
    }
    let den: BigInt = "10000000000000000000000000000000000000"
        .parse()
        .expect("10^37");
    let lo: BigInt = "31415926535897932384626433832795028841"
        .parse()
        .expect("π low");
    let hi: BigInt = "31415926535897932384626433832795028842"
        .parse()
        .expect("π high");
    // f(x) at x = num/den has the sign of a·den + b·num.
    let at = |num: &BigInt| sign(&(a * &den + b * num));
    let (s_lo, s_hi) = (at(&lo), at(&hi));
    (s_lo == s_hi && s_lo != Orient::Zero).then_some(s_lo)
}

/// A circular arc for [`winding_sign_quarter_arcs`]: `start → end` about `center`,
/// counter-clockwise when `ccw`; `start == end` is the whole circle. `r2` is the squared radius —
/// the form the sketch truth holds, and the one the area term wants.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QuarterArc {
    pub center: [Rat; 2],
    pub r2: Rat,
    pub start: [Rat; 2],
    pub end: [Rat; 2],
    pub ccw: bool,
}

/// **Which way a closed ring of straight steps and quarter-turn arcs runs**, exactly: `Positive`
/// is counter-clockwise (`+x → +y`). Twice the signed area is `∮ x dy − y dx`; a straight step
/// `p → q` contributes the shoelace term `p × q`, an arc `S → E` about `C` through the signed
/// angle `φ` contributes `C × (E − S) + r²·φ` (parametrize `p = C + r·u(t)`: `x dy − y dx =
/// r·(C × u′) + r²` dt). With every `φ` a multiple of `π/2` the total is `a + b·π`, whose sign
/// [`sign_a_plus_b_pi_int`] decides.
///
/// ★ **Integers, not `Rat`.** Every coordinate goes over one common denominator and the sums run
/// in `BigInt`: two 17-digit decimals multiplied already leave `i128` (the reason [`orient2d_rat`]
/// carries a `BigInt` arm), and a winding read that failed open on overflow would build a solid
/// inside out. The area is doubled so `r²·k/2` stays integral. `None` for an arc that is not a
/// quarter-turn multiple (no producer makes one today — the kit's fillet is axis-aligned, so its
/// `arc_rat` arcs are quarter-turns; the type allows more) or for a sum that lands inside the π
/// bracket.
pub fn winding_sign_quarter_arcs(lines: &[[[Rat; 2]; 2]], arcs: &[QuarterArc]) -> Option<Orient> {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;
    let mut den = BigInt::from(1);
    let mut fold = |r: &Rat| den = den.lcm(&BigInt::from(r.denom()));
    for [p, q] in lines {
        for c in p.iter().chain(q.iter()) {
            fold(c);
        }
    }
    for a in arcs {
        for c in a.center.iter().chain(a.start.iter()).chain(a.end.iter()) {
            fold(c);
        }
        fold(&a.r2);
    }
    let int = |r: &Rat| -> BigInt { BigInt::from(r.numer()) * (&den / BigInt::from(r.denom())) };
    let pt = |p: &[Rat; 2]| -> [BigInt; 2] { [int(&p[0]), int(&p[1])] };
    let cross2 = |a: &[BigInt; 2], b: &[BigInt; 2]| -> BigInt { &a[0] * &b[1] - &a[1] * &b[0] };
    let dot2 = |a: &[BigInt; 2], b: &[BigInt; 2]| -> BigInt { &a[0] * &b[0] + &a[1] * &b[1] };
    let sub2 = |a: &[BigInt; 2], b: &[BigInt; 2]| -> [BigInt; 2] { [&a[0] - &b[0], &a[1] - &b[1]] };
    // `a_int` collects 2·(shoelace + C × (E − S)), `b_int` collects r²·k — units of den².
    let (mut a_int, mut b_int) = (BigInt::from(0), BigInt::from(0));
    for [p, q] in lines {
        a_int += 2 * cross2(&pt(p), &pt(q));
    }
    for arc in arcs {
        let (c, s0, e0) = (pt(&arc.center), pt(&arc.start), pt(&arc.end));
        a_int += 2 * cross2(&c, &sub2(&e0, &s0));
        // The counter-clockwise quarter turns from S to E.
        let quarters_ccw: i128 = if s0 == e0 {
            4
        } else {
            let (v1, v2) = (sub2(&s0, &c), sub2(&e0, &c));
            match (dot2(&v1, &v2).sign(), cross2(&v1, &v2).sign()) {
                (Sign::NoSign, Sign::Plus) => 1,
                (Sign::Minus, Sign::NoSign) => 2,
                (Sign::NoSign, Sign::Minus) => 3,
                _ => return None,
            }
        };
        // φ = k·π/2 signed by the direction walked — a whole turn is ±4 either way, a partial arc
        // walked clockwise is the complement, negated; 2·r²·φ = (r²·k)·π.
        let k = match (arc.ccw, quarters_ccw) {
            (true, q) => q,
            (false, 4) => -4,
            (false, q) => q - 4,
        };
        // `r²` is stated, so `int` puts it in units of `den`; one more `den` brings it to the
        // `den²` the area terms carry.
        b_int += int(&arc.r2) * &den * BigInt::from(k);
    }
    sign_a_plus_b_pi_int(&a_int, &b_int)
}

/// [`sign_a_plus_b_pi_int`] for rationals: `a + b·π` and `(a·db + b·da·π)` share a sign since
/// `da·db > 0`.
pub fn sign_a_plus_b_pi(a: Rat, b: Rat) -> Option<Orient> {
    use num_bigint::BigInt;
    let (an, ad) = (BigInt::from(a.numer()), BigInt::from(a.denom()));
    let (bn, bd) = (BigInt::from(b.numer()), BigInt::from(b.denom()));
    sign_a_plus_b_pi_int(&(an * &bd), &(bn * &ad))
}
