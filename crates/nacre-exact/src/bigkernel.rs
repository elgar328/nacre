use super::*;
/// `-1 / 0 / +1` of a `BigInt`.
pub(super) fn big_sign(x: &num_bigint::BigInt) -> i8 {
    use num_bigint::Sign;
    match x.sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}

/// `det3` over borrowed `BigInt` entries — the same cofactor expansion
/// `nacre_predicates::det3` evaluates through `Expansion`, exact by being integer arithmetic
/// rather than by tracking roundoff.
pub(super) fn det3_big(m: [[&num_bigint::BigInt; 3]; 3]) -> num_bigint::BigInt {
    let minor = |p: &num_bigint::BigInt,
                 q: &num_bigint::BigInt,
                 r: &num_bigint::BigInt,
                 t: &num_bigint::BigInt| p * q - r * t;
    m[0][0] * minor(m[1][1], m[2][2], m[1][2], m[2][1])
        - m[0][1] * minor(m[1][0], m[2][2], m[1][2], m[2][0])
        + m[0][2] * minor(m[1][0], m[2][1], m[1][1], m[2][0])
}

/// The implicit point `∩(p₀, p₁, p₂)` by Cramer, over integer plane coefficients: numerators
/// `(Dx, Dy, Dz)` over the common denominator `D`. The integer twin of `nacre_predicates`'
/// `cramer`, for coefficients too wide for any `f64`-chunk representation — an `Expansion`'s
/// pieces are `f64`s, so its exponent range ends near `2¹⁰²³` while a wide canonical name is
/// measured out to `~2²²⁹¹`.
pub(super) fn cramer_big(
    p: [&[num_bigint::BigInt; 4]; 3],
) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
    let rhs = [-&p[0][3], -&p[1][3], -&p[2][3]];
    let row = |r: usize, k: usize| {
        let mut m = [&p[r][0], &p[r][1], &p[r][2]];
        m[k] = &rhs[r];
        m
    };
    let col_replaced = |k: usize| det3_big([row(0, k), row(1, k), row(2, k)]);
    (
        [col_replaced(0), col_replaced(1), col_replaced(2)],
        det3_big(p.map(|r| [&r[0], &r[1], &r[2]])),
    )
}

// ---------------------------------------------------------------------------
// Total rational-vector predicates.
//
// A question whose answer is a **sign** has no width: the answer is one of three values, and
// only the road to it can overflow. These clear the denominators once and run in `BigInt`, so
// they cannot decline — which is what lets their callers say `None`/reject for the geometry
// alone. The same move `three_planes_rat` makes for its value ("the caller resolves the
// conflation by asking the integer core") and `plane_name_exact` for its name.
//
// ★ **The lift keeps each vector's denominator.** Clearing denominators multiplies a vector by
// a positive factor, which is harmless for an expression that is *homogeneous* in that vector
// (a cross, a dot) and **wrong** for one that is not — an expression mixing a term in `p` with
// a constant would state a different proposition after scaling `p` alone. Every predicate here
// says which case it is in, so the ones that need the scale can multiply it back in.
// ---------------------------------------------------------------------------

/// A rational 3-vector as **(integer components, positive denominator)**: `v = out.0 / out.1`.
pub(super) fn lift3(v: &[Rat; 3]) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let den: [BigInt; 3] = core::array::from_fn(|i| BigInt::from(v[i].denom()));
    let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
    let num = core::array::from_fn(|i| BigInt::from(v[i].numer()) * (&d / &den[i]));
    (num, d)
}
