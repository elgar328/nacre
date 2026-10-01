use super::*;
impl BigRat {
    pub fn from_rat(r: Rat) -> BigRat {
        BigRat(Ratio::new(
            num_bigint::BigInt::from(r.numer()),
            num_bigint::BigInt::from(r.denom()),
        ))
    }

    /// `r²`, exactly — the one product this type exists for. Never overflows; `r` in lowest terms
    /// makes `r²` so too.
    pub fn square_of(r: Rat) -> BigRat {
        let (n, d) = (
            num_bigint::BigInt::from(r.numer()),
            num_bigint::BigInt::from(r.denom()),
        );
        BigRat(Ratio::new_raw(&n * &n, &d * &d))
    }

    pub fn zero() -> BigRat {
        BigRat(Ratio::from_integer(num_bigint::BigInt::from(0)))
    }

    pub fn numer(&self) -> &num_bigint::BigInt {
        self.0.numer()
    }

    pub fn denom(&self) -> &num_bigint::BigInt {
        self.0.denom()
    }

    /// The `Rat` this is, when both parts fit `i128` — `None` is width, not a value.
    pub fn narrow(&self) -> Option<Rat> {
        use num_traits::ToPrimitive;
        Rat::new(self.0.numer().to_i128()?, self.0.denom().to_i128()?)
    }

    pub fn is_positive(&self) -> bool {
        self.0.numer().sign() == num_bigint::Sign::Plus
    }

    pub fn is_negative(&self) -> bool {
        self.0.numer().sign() == num_bigint::Sign::Minus
    }

    pub fn is_zero(&self) -> bool {
        self.0.numer().sign() == num_bigint::Sign::NoSign
    }

    /// `self · r`, exactly.
    pub fn mul_rat(&self, r: Rat) -> BigRat {
        BigRat(&self.0 * BigRat::from_rat(r).0)
    }
}

impl From<Rat> for BigRat {
    fn from(r: Rat) -> BigRat {
        BigRat::from_rat(r)
    }
}

/// [`rat_sqrt_exact`] for a wide radicand: the root, when it is rational **and** fits `Rat` —
/// which every stated radius does, its square being what widened.
pub fn rat_sqrt_exact_big(v: &BigRat) -> Option<Rat> {
    if v.is_negative() {
        return None;
    }
    let (n, d) = (v.numer(), v.denom());
    let (rn, rd) = (n.sqrt(), d.sqrt());
    if &rn * &rn != *n || &rd * &rd != *d {
        return None;
    }
    BigRat(Ratio::new_raw(rn, rd)).narrow()
}

pub fn rat_sqrt_exact(v: Rat) -> Option<Rat> {
    let (num, den) = (*v.0.numer(), *v.0.denom());
    if num < 0 {
        return None;
    }
    let root = |x: i128| -> Option<i128> {
        let r = (x as f64).sqrt() as i128;
        // `as f64` rounds past 2⁵³, so search a small neighbourhood rather than trusting it.
        (r.saturating_sub(2).max(0)..=r.saturating_add(2)).find(|&c| c.checked_mul(c) == Some(x))
    };
    Rat::new(root(num)?, root(den)?)
}

/// The exact `1/√v` of a positive rational, or `None` when it is irrational.
///
/// `1/√v = √(1/v)`, so this is [`rat_sqrt_exact`] of the reciprocal — and the test is symmetric
/// under inversion (a fraction in lowest terms is a perfect square exactly when both of its parts
/// are), so nothing is lost by asking it that way round.
///
/// ★ **This branch is why an axis-aligned model pays nothing.** A frame normal of `(0, 0, 1)` has
/// `n·n = 1` and lands here with `1.0`, exactly, without touching arbitrary precision — the same
/// role [`Angle::try_exact_cos_sin`] plays for a quadrantal rotation.
pub fn inv_sqrt_exact(v: Rat) -> Option<Rat> {
    if v <= Rat::from_int(0) {
        return None;
    }
    rat_sqrt_exact(Rat::new(*v.0.denom(), *v.0.numer())?)
}

/// `1/√v` at `prec` bits **with an upper bound on how far it may be from the true value** — the
/// seed for every error radius that a frame's realization grows from.
///
/// The bound is *derived*, term by term, in the shape [`Angle::cos_sin_bounded`] uses:
///
/// - `numer`/`denom` enter as `i128`, exactly. Routing them through `f64` would cost a relative
///   `2⁻⁵³` past 2⁵³ that no working precision recovers.
/// - `x = fl(n/d)` is one rounded operation: relative `u = 2⁻ᵖʳᵉᶜ`.
/// - `√` **halves** a relative error (`d√x/√x = ½ · dx/x`), so the argument arrives as `u/2`, and
///   the square root's own realization adds a half-ulp: `u`.
/// - The reciprocal passes relative error through unchanged (`d(1/y)/(1/y) = −dy/y`) and adds its
///   own `u`.
///
/// That totals `2.5u`; the bound below uses `4u`, and
/// `the_inverse_sqrt_bound_holds_against_a_far_deeper_realization` checks it against a far deeper
/// realization rather than trusting the arithmetic or astro-float's rounding contract.
///
/// **Memoized by `(v, prec)`** in [`INV_SQRT`], for the same reason [`TRIG`] exists: a frame's
/// realization asks for the same `n·n` once per coordinate.
///
/// `None` when `v ≤ 0` — there is no frame normal with a non-positive squared length, so that is a
/// broken premise rather than an unusual input.
pub fn inv_sqrt_bounded(v: Rat, prec: usize) -> Option<HpBounded> {
    if v <= Rat::from_int(0) {
        return None;
    }
    if let Some(hit) = INV_SQRT.with_borrow(|t| t.get(&(v, prec)).cloned()) {
        return Some(hit);
    }
    let out = realize_inv_sqrt(v, prec);
    INV_SQRT.with_borrow_mut(|t| t.insert((v, prec), out.clone()));
    Some(out)
}

/// An arbitrary-precision integer as a `BigFloat`, **exactly** — the working precision is the
/// integer's own bit length, so no digit is ever rounded away.
///
/// The wide-frame realization feeds plane data wider than `i128` through this; keeping the
/// conversion exact is what keeps [`inv_sqrt_bigint_bounded`]'s error budget identical to
/// [`inv_sqrt_bounded`]'s: every rounding still happens *after* the value has entered whole,
/// exactly as the narrow twin's `i128`s do.
///
/// Horner over the base-2⁶⁴ digits: the scale multiply is an exponent shift (exact), and each
/// add lands in a mantissa wide enough for the whole running value (exact).
pub fn bigint_to_bigfloat(x: &num_bigint::BigInt, prec_floor: usize) -> BigFloat {
    let bits = x.magnitude().bits() as usize;
    let p = bits.max(prec_floor).max(128);
    let scale = BigFloat::from_i128(1i128 << 64, p);
    let mut acc = BigFloat::from_i128(0, p);
    for d in x.magnitude().iter_u64_digits().rev() {
        acc = acc.mul(&scale, p, HP_RM);
        acc = acc.add(&BigFloat::from_i128(d as i128, p), p, HP_RM);
    }
    if x.sign() == num_bigint::Sign::Minus {
        acc = acc.neg();
    }
    acc
}

/// [`inv_sqrt_bounded`] for a squared length wider than `i128` — the wide-frame twin.
///
/// Same ladder, same derived bound: the integer enters **exactly** ([`bigint_to_bigfloat`] at
/// its own bit length), the square root rounds once (halving the incoming relative error, which
/// is zero here) and the reciprocal rounds once — under the narrow twin's `4u`, which is kept
/// for symmetry rather than tightened.
///
/// **Uncached** — the wide population is a fraction of a percent of pushes and the memo key
/// would be a `BigInt`; measured before optimizing, per the cache philosophy.
pub fn inv_sqrt_bigint_bounded(v: &num_bigint::BigInt, prec: usize) -> Option<HpBounded> {
    if v.sign() != num_bigint::Sign::Plus {
        return None;
    }
    let ip = prec.max(128);
    let x = bigint_to_bigfloat(v, ip);
    let one = BigFloat::from_i128(1, ip);
    let z = one.div(&x.sqrt(prec, HP_RM), prec, HP_RM);
    let u = Mag::pow2(-(prec as i64));
    let rel = u.times(Mag::of(4.0));
    let error = Mag::above(&z).times(rel);
    Some(HpBounded::new(z, error))
}

/// [`inv_sqrt_bounded`] without the memo — the evaluation itself, kept separate so no `INV_SQRT`
/// borrow is held across the arbitrary-precision work.
fn realize_inv_sqrt(v: Rat, prec: usize) -> HpBounded {
    // As `i128`, not through `f64`: the loss would happen before astro-float saw the value.
    let ip = prec.max(128);
    let n = BigFloat::from_i128(*v.0.numer(), ip);
    let d = BigFloat::from_i128(*v.0.denom(), ip);
    let one = BigFloat::from_i128(1, ip);
    let x = n.div(&d, prec, HP_RM);
    let z = one.div(&x.sqrt(prec, HP_RM), prec, HP_RM);
    let u = Mag::pow2(-(prec as i64));
    // `½ + 1 + 1 = 2.5`, rounded up. Relative, so it is scaled by the result's magnitude below.
    let rel = u.times(Mag::of(4.0));
    let error = Mag::above(&z).times(rel); // `|x| < 2^exponent`
    HpBounded::new(z, error)
}

/// **`1/√v` as the `f64` nearest the true value**, or `None` when `v ≤ 0`.
///
/// ★★ **The exact branch runs first** ([`inv_sqrt_exact`]) and covers every axis-aligned frame,
/// plus the Pythagorean ones a CAD user actually draws — a `(3, 4, 0)` normal has `n·n = 25` and
/// `1/|n| = 1/5`. Only a genuinely irrational length reaches the ladder.
///
/// ★ **Exact does not mean free of rounding**: `1/5` is exactly rational and still not an `f64`,
/// so that branch returns the *nearest* f64 to a known-exact value. [`inv_sqrt_error_of`] reports
/// what that rounding cost, and returns a literal zero only where there was none.
///
/// **The ladder is 128 then 256 bits, and it terminates.** Unlike `cos 90°`, `1/√v` is never zero
/// for a positive `v`, so its interval never straddles zero and [`round_to_f64`] cannot answer
/// `None` forever. The rungs match [`Angle::realize_rounded_f64`]'s so a judgement at
/// `nacre_judge`'s trial precision shares this realization instead of paying for a second one.
///
/// ★ **Correct rounding is what keeps debug and release the same.** A faithfully-rounded value
/// would let the two builds disagree by an ulp, which is the failure this crate already paid for
/// once in the trig path.
///
/// If even 256 bits leave the rounding undecided the midpoint is returned rather than a panic, and
/// a test build counts the event (`INV_SQRT_ESCALATED`).
pub fn inv_sqrt_f64(v: Rat) -> Option<f64> {
    if v <= Rat::from_int(0) {
        return None;
    }
    if let Some(r) = inv_sqrt_exact(v) {
        return Some(r.to_f64());
    }
    if let Some(hit) = F64_INV_SQRT.with_borrow(|m| m.get(&v).copied()) {
        return Some(hit);
    }
    // Outside the borrow: the realization below takes `INV_SQRT`'s in turn.
    let out = realize_inv_sqrt_rounded(v);
    F64_INV_SQRT.with_borrow_mut(|m| m.insert(v, out));
    Some(out)
}

/// [`inv_sqrt_f64`]'s general branch without the memo or the exact test.
fn realize_inv_sqrt_rounded(v: Rat) -> f64 {
    for (i, prec) in [128usize, 256].into_iter().enumerate() {
        let HpBounded {
            value: z,
            error: rad,
        } = realize_inv_sqrt_memoized(v, prec);
        if let Some(f) = round_to_f64(&z, rad, prec) {
            #[cfg(test)]
            if i > 0 {
                INV_SQRT_ESCALATED.with_borrow_mut(|(e, _)| *e += 1);
            }
            #[cfg(not(test))]
            let _ = i;
            return f;
        }
    }
    #[cfg(test)]
    INV_SQRT_ESCALATED.with_borrow_mut(|(_, f)| *f += 1);
    let z = realize_inv_sqrt_memoized(v, 256).value;
    to_f64_exact(&z).unwrap_or(f64::NAN)
}

/// [`inv_sqrt_bounded`] for a `v` already known positive.
fn realize_inv_sqrt_memoized(v: Rat, prec: usize) -> HpBounded {
    inv_sqrt_bounded(v, prec).expect("v > 0 checked by the caller")
}

/// **`√v` as the `f64` nearest the true value**, or `None` when `v < 0` — the f64 realization of
/// a radius stated as its square, at any width.
///
/// The twin of [`inv_sqrt_f64`]: the exact branch runs first ([`rat_sqrt_exact`] — every radius a
/// user writes as a decimal lands here, and comes back as the `to_f64` of the rational it always
/// was, bit for bit), and only a genuinely irrational root reaches the 128 → 256 bit ladder, under
/// the same correct-rounding contract. `√v` is nonzero for `v > 0`, so the interval never
/// straddles zero and [`round_to_f64`] cannot decline forever; the documented fallback to the
/// midpoint is kept for the same reason the inverse has it.
pub fn sqrt_f64(v: &BigRat) -> Option<f64> {
    if v.is_negative() {
        return None;
    }
    if let Some(r) = rat_sqrt_exact_big(v) {
        return Some(r.to_f64());
    }
    for prec in [128usize, 256] {
        let HpBounded {
            value: z,
            error: rad,
        } = sqrt_bounded_big(v, prec)?;
        if let Some(f) = round_to_f64(&z, rad, prec) {
            return Some(f);
        }
    }
    let z = sqrt_bounded_big(v, 256)?.value;
    Some(to_f64_exact(&z).unwrap_or(f64::NAN))
}

/// **The unit direction of an integer vector, each component the `f64` nearest `nₖ/|n|`**, or
/// `None` for the zero vector — the realization of a plane's normal from its name.
///
/// The twin of [`sqrt_f64`], under the same correct-rounding contract, so the answer is unique and
/// asking at more bits cannot move it. The exact branch runs first: when `n·n` is a perfect square
/// every component is a rational, rounded once ([`nearest_f64_big_exact`]) — an axis-aligned or
/// Pythagorean normal never touches arbitrary precision. Otherwise `1/|n|` is realized at 128 and
/// then 256 bits — through the memoized [`inv_sqrt_bounded`] where `n·n` fits `Rat`, and
/// [`inv_sqrt_bigint_bounded`] where it does not — and each product rounded; the midpoint stands
/// if even 256 bits leave a component undecided, as in [`sqrt_f64`].
///
/// ★ **A zero component is `+0.0`.** The sign of the vector belongs to the integers handed in, so
/// a caller that turns the direction negates them before asking — negating the `f64` afterwards
/// would write `-0.0`, a different bit pattern for the same value.
pub fn unit_vector_f64(n: &[num_bigint::BigInt; 3]) -> Option<[f64; 3]> {
    use num_traits::ToPrimitive;
    let nn = &n[0] * &n[0] + &n[1] * &n[1] + &n[2] * &n[2];
    if nn.sign() == num_bigint::Sign::NoSign {
        return None;
    }
    let root = nn.sqrt();
    if &root * &root == nn {
        let c = |k: usize| nearest_f64_big_exact(&n[k], &root).map(|(v, _)| v);
        return Some([c(0)?, c(1)?, c(2)?]);
    }
    let inv = |prec: usize| match nn.to_i128() {
        Some(v) => inv_sqrt_bounded(Rat::from_int(v), prec),
        None => inv_sqrt_bigint_bounded(&nn, prec),
    };
    for prec in [128usize, 256] {
        let s = inv(prec)?;
        let c = |k: usize| {
            let v = HpBounded::of_bigint(&n[k], prec).mul(&s, prec);
            round_to_f64(&v.value, v.error, prec)
        };
        if let (Some(x), Some(y), Some(z)) = (c(0), c(1), c(2)) {
            return Some([x, y, z]);
        }
    }
    let s = inv(256)?;
    let c = |k: usize| {
        let v = HpBounded::of_bigint(&n[k], 256).mul(&s, 256).value;
        to_f64_exact(&v).unwrap_or(f64::NAN)
    };
    Some([c(0), c(1), c(2)])
}

/// **How far the `1/√v` the caller was handed sits from the true one** — measured against an
/// arbitrary-precision realization, the twin of [`Angle::realization_error_of`].
///
/// ★★★ **The caller passes the value in rather than letting this re-realize it**, for the same
/// reason the trig one does: the consumer is a frame's `tol`, which must bound the error in the
/// `coord` it wrote from *its* `inv_sqrt_f64` result. An error measured against a second,
/// independent realization would bound a number nobody stored.
///
/// ★★★★ **Exactly zero only when the exact value is also *dyadic*, which is not the same thing
/// as being exact.** [`Angle::realization_error_of`] can return a flat `0` for its exact family
/// because `{0, ±1}` are f64 values; `1/√v` cannot. A **Pythagorean** normal like `(3, 4, 0)` has
/// `n·n = 25` and an exactly rational `1/|n| = 1/5` — which is *not* an f64, so the realization
/// still rounds. Copying the trig zero here reported `0` for a real `2⁻⁵³` error, and
/// `the_inverse_sqrt_realization_error_covers_the_error_that_is_there` caught it. So the exact
/// branch **measures the rational's own rounding** instead, in exact arithmetic, and reaches `0`
/// where it genuinely belongs: an axis-aligned frame, whose `1/|n|` is `1`, `½`, `¼`, …
///
/// Callers rely on that zero: a frame whose realization carries no error also performs no rounding
/// downstream, which is what keeps an axis-aligned sketch at tol 0.
///
/// ★ **This measures rather than assumes even though [`inv_sqrt_f64`] is correctly rounded.**
/// `round_to_f64` returning `Some` *is* a certificate, so `½ ulp` would be defensible — but the
/// ladder has a documented fallback for an undecided rounding at 256 bits, and a term derived from
/// a guarantee that has an escape hatch is exactly the kind that goes quietly wrong. Measuring
/// covers the fallback for free.
///
/// The reading is at octave granularity, so it can sit up to 2× above the true error —
/// conservative in the sound direction, and still a measurement. `None` when `v ≤ 0`.
pub fn inv_sqrt_error_of(v: Rat, f: f64) -> Option<f64> {
    if v <= Rat::from_int(0) {
        return None;
    }
    // ★ The claim is made of *this* value, not of `v`: what a caller depends on is "the number I
    // am holding is the realization of the true `1/√v`", and only comparing says that.
    //
    // The gap is taken in **exact rational arithmetic** — `try_from_f64` is exact, so `f − r` is
    // exact and its being zero is a fact rather than a measurement below some resolution. That is
    // what lets an axis-aligned frame reach a literal `0.0`; going through `inv_sqrt_bounded`
    // would hand it the realization's own `2⁻¹²⁴` radius and no exact route would ever be taken.
    if let Some(r) = inv_sqrt_exact(v) {
        if r.to_f64() == f {
            let d = Rat::try_from_f64(f).and_then(|fr| fr.checked_sub(r))?;
            if d == Rat::from_int(0) {
                return Some(0.0);
            }
            // `to_f64` rounds to nearest, so nudge up to keep the bound above the truth.
            return Some(d.to_f64().abs() * (1.0 + 2.0 * f64::EPSILON));
        }
    }
    const P: usize = 128; // the hp radius is then ~2⁻¹²⁸ against an ε-scale quantity
    let HpBounded {
        value: h,
        error: rad,
    } = inv_sqrt_bounded(v, P)?;
    let diff = BigFloat::from_f64(f, P).sub(&h, P, HP_RM);
    let mag = if diff.is_zero() {
        0.0
    } else {
        2f64.powi(diff.exponent().unwrap_or(0))
    };
    // Rounded up, for the reason `Angle::measure_realization_error` spells out: when `|diff|` is
    // itself a power of two the octave bound has no slack, and `mag + rad` would round back down
    // to `mag` in f64 — short by the radius.
    Some((mag + rad.upper_f64()) * (1.0 + 2.0 * f64::EPSILON))
}

/// Greatest common divisor of two magnitudes, Euclid. `gcd(0, 0) == 0`.
pub(super) fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}
