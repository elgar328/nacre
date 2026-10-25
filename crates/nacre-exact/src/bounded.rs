//! **A value and the error it is caught in** — one machine at two precisions.
//!
//! [`Bounded`] is the `f64` filter: a value with an error radius, propagated through every
//! operation, whose sign is only reported when the interval clears zero. [`HpBounded`] is the same
//! thing in astro-float, differing in exactly one place — the per-operation rounding is `2⁻ᵖʳᵉᶜ`
//! instead of `2·ε`. Named for what the type promises — the truth lies within `value ± error` —
//! not for what the value is (an approximation), because the promise is the whole point.
//!
//! **They live in one file on purpose.** Both are one kind of machine — each propagates a radius —
//! and side by side, an operation that propagates its radius in one and not the other is visible.
//! (The alternative, a determinant held against a hand-picked floor built from a
//! separately-maintained magnitude estimate, collapses when the estimate is read off a value that
//! has already cancelled, and reports a rounding residue as a confident sign.)
//!
//! **And they live in this crate, not in the judge that runs them,** for the same reason: this
//! crate realizes curved coordinates with the same high-precision arithmetic, and two spellings of
//! one machine in two crates is how one of them drifts. One spelling here, and the judge imports
//! it.
//!
//! **The radius is a [`Mag`], not an `f64`.** At a deep rung `2⁻ᵖʳᵉᶜ` underflows an `f64` to
//! zero, and a zero radius claims exactness — the same failure in new clothes. See [`Mag`]'s own
//! doc.

use astro_float::BigFloat;

use crate::{HP_RM, Mag, Rat};

/// A rational as an arbitrary-precision float.
///
/// **Numerator and denominator go in as `i128`, not through `f64`.** Routing them through `f64`
/// costs a relative `2⁻⁵³` the moment either exceeds 2⁵³ — and *no amount of working precision
/// recovers it*, because the loss happens before astro-float ever sees the value. Measured: a
/// denominator of `2⁵⁴+1` took a coordinate's error bound from `2⁻²⁵¹` to `2⁻⁵⁰`. Chained exact
/// rational arithmetic multiplies denominators, so that is not an exotic input; it is what the
/// crate's own bit-growth note is about.
pub fn rat_to_bigfloat(r: Rat, prec: usize) -> BigFloat {
    // `from_i128` needs at least 128 bits to hold the integer before the division rounds it.
    let p = prec.max(128);
    BigFloat::from_i128(r.numer(), p).div(&BigFloat::from_i128(r.denom(), p), prec, HP_RM)
}

/// A value with a symmetric error radius (`value ± error`, `error ≥ 0`). Arithmetic keeps
/// `error` a sound upper bound (worst case), plus a per-op f64-rounding inflation.
#[derive(Clone, Copy, Debug)]
pub struct Bounded {
    pub value: f64,
    pub error: f64,
}

/// The radius's own rounding, charged upward once per operation: `1 + 4ε`, exactly representable
/// (`1 + 2⁻⁵⁰`). See [`Bounded::sub`].
const ERROR_UP: f64 = 1.0 + 4.0 * f64::EPSILON;

// `add`/`sub`/`mul` are spelled like the std traits on purpose and are not them: the same names
// at both precisions is what lets the two machines be read side by side, and `HpBounded`'s take
// a `prec` the trait signatures have no room for. Naming these apart would split the one machine.
#[allow(clippy::should_implement_trait)]
impl Bounded {
    pub fn new(value: f64, error: f64) -> Self {
        Bounded { value, error }
    }
    /// The rounding that brought a value to this nearest `f64`: at most a half-ulp, `|value|·ε/2`,
    /// floored at the smallest normal, as [`Mag::upper_f64`] and [`crate::bigfloat_upper_f64`]
    /// floor: a subnormal's rounding error is a fixed half step, which `|value|·ε/2` would read as
    /// zero. The interval arithmetic below (`add`·`sub`·`mul`) keeps no such floor — it charges
    /// `2ε·|value|` and rests on the values staying in the normal range. `HpBounded::round_off` is
    /// the same charge at `prec` bits. Never zero, `0.0` included — the `f64` alone does not say
    /// its rounding was exact; a caller that knows gives zero itself.
    pub fn round_off(value: f64) -> f64 {
        (value.abs() * (f64::EPSILON * 0.5)).max(f64::MIN_POSITIVE)
    }
    /// A rational as its nearest `f64` **with the error that rounding carries** — [`HpBounded::of_rat`]
    /// at `f64`. `0` exactly when the rational is an `f64` ([`Rat::as_f64`]); otherwise
    /// [`Self::round_off`], because [`Rat::to_f64`] rounds to nearest and so lands within half an
    /// ulp.
    ///
    /// From that contract, not measured: a measurement at a fixed precision reads a rational
    /// nearer its `f64` than that precision resolves as exactly representable.
    pub fn of_rat(r: Rat) -> Self {
        match r.as_f64() {
            Some(f) => Bounded::new(f, 0.0),
            None => {
                let f = r.to_f64();
                Bounded::new(f, Bounded::round_off(f))
            }
        }
    }
    /// **Each result's radius is a proved bound, operation by operation.** Two charges:
    ///
    /// - `2ε·|value|` for the value's own rounding. Round-to-nearest misses the exact result by
    ///   at most `ε/2` of it, which is under `ε·|value|` of the rounded one; the charge is twice that.
    /// - `× ERROR_UP` (`1 + 4ε`, exact in `f64`) on the summed radius, for the radius's own
    ///   arithmetic. That sum is itself rounded to nearest, at most three times in `add`/`sub` (the
    ///   input sum, the value charge, this product) and five in `mul` (a radius product, two sums,
    ///   the value charge, this product), each step at most `ε/2` low — so the computed radius is at
    ///   least `(1 − ε/2)⁵ · (1 + 4ε) > 1` times the exact one. Without it the sum could come out
    ///   below the radius it stands for wherever the radius outweighs the value, which is most of a
    ///   determinant's intermediate steps.
    ///
    /// Both rest on the values staying in the normal range (design 「수의 이름」). Measured slack on
    /// the soundness corpora: the worst actual error is `0.101` of the radius for the `orient3d`
    /// determinant (5,000 rotated four-point configurations, `orient3d_iv_soundness`) and `0.347`
    /// for plane coefficients (10,000, `plane_coefficient_tol_soundness`).
    pub fn sub(self, o: Bounded) -> Bounded {
        let value = self.value - o.value;
        Bounded::new(
            value,
            (self.error + o.error + 2.0 * f64::EPSILON * value.abs()) * ERROR_UP,
        )
    }
    pub fn add(self, o: Bounded) -> Bounded {
        let value = self.value + o.value;
        Bounded::new(
            value,
            (self.error + o.error + 2.0 * f64::EPSILON * value.abs()) * ERROR_UP,
        )
    }
    pub fn mul(self, o: Bounded) -> Bounded {
        let value = self.value * o.value;
        // Worst-case product radius `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, plus the
        // f64 rounding of the product itself.
        let error = self.value.abs() * o.error + o.value.abs() * self.error + self.error * o.error;
        Bounded::new(value, (error + 2.0 * f64::EPSILON * value.abs()) * ERROR_UP)
    }
    /// `Some(true)` if definitely positive, `Some(false)` if definitely negative,
    /// `None` if the interval straddles 0 (escalate).
    pub fn sign(self) -> Option<bool> {
        if self.value > self.error {
            Some(true)
        } else if self.value < -self.error {
            Some(false)
        } else {
            None
        }
    }
}

/// [`Bounded`] in astro-float: a high-precision value with a **computed** error radius.
///
/// Every constructor must supply a radius that actually bounds its value's distance from the
/// truth — for a coordinate that means the rotation chain's realization error
/// ([`crate::Angle::cos_sin_hp`] is where it starts), and for a rational read at enough
/// bits it means [`Mag::ZERO`]. A radius invented for convenience makes every sign above it
/// unearned.
///
/// A consumer holds one without naming astro-float (`nacre-ops` does), and the radius is not
/// decoration: [`crate::round_to_f64`] and [`crate::round_to_digits`] both need it to say whether
/// an answer is determined — a value without it can only be rounded by guessing.
#[derive(Clone, Debug)]
pub struct HpBounded {
    pub value: BigFloat,
    pub error: Mag,
}

impl HpBounded {
    pub fn new(value: BigFloat, error: Mag) -> Self {
        HpBounded { value, error }
    }

    /// A value known exactly at this precision — a rational whose realization did not round.
    pub fn exact(value: BigFloat) -> Self {
        HpBounded {
            value,
            error: Mag::ZERO,
        }
    }

    /// A rational realized at `prec` bits **with the error that realization carries**.
    ///
    /// The integers themselves are exact — [`rat_to_bigfloat`] feeds them in as `i128` — so the
    /// only error is the division, and even that vanishes when it terminates: a power-of-two
    /// denominator with a numerator inside `prec` bits is exact, and then the radius is genuinely
    /// zero.
    pub fn of_rat(r: Rat, prec: usize) -> Self {
        let value = rat_to_bigfloat(r, prec);
        let (n, d) = (r.numer(), r.denom()); // `d > 0` after reduction
        let n_bits = (128 - n.unsigned_abs().leading_zeros()) as usize;
        if d & (d - 1) == 0 && n_bits <= prec {
            return HpBounded::exact(value); // a terminating division: a power-of-two denominator
        }
        // The integers enter exactly (see `rat_to_bigfloat`), so the only error left is the
        // division's own rounding.
        let error = Self::round_off(&value, prec);
        HpBounded::new(value, error)
    }

    /// An `f64` as an **exact** interval: every finite `f64` is a binary fraction, held whole at
    /// `prec` bits (at least 64), so the radius is zero. `None` for an infinite or NaN value.
    ///
    /// It reads the `f64` exactly and nothing more: the `f64` is not the truth it may stand for
    /// (a cache, a figure), so this is a way to measure the very number a cache or a file holds,
    /// never a door from an `f64` to a truth.
    pub fn of_f64(v: f64, prec: usize) -> Option<Self> {
        v.is_finite()
            .then(|| HpBounded::exact(BigFloat::from_f64(v, prec.max(64))))
    }

    /// An arbitrary-precision integer as an **exact** interval — the wide-frame entry point.
    ///
    /// [`crate::bigint_to_bigfloat`] converts at the integer's own bit length, so the value enters
    /// whole and the radius is genuinely zero; downstream operations charge their own rounding,
    /// exactly as [`Self::of_rat`]'s exact branch does.
    pub fn of_bigint(x: &num_bigint::BigInt, prec: usize) -> Self {
        HpBounded::exact(crate::bigint_to_bigfloat(x, prec))
    }

    /// The rounding a `prec`-bit operation adds to its own result: at most a half-ulp,
    /// `|result| · 2⁻ᵖʳᵉᶜ` — [`Mag::ZERO`] for an exact zero, which rounding cannot move.
    pub(crate) fn round_off(value: &BigFloat, prec: usize) -> Mag {
        Mag::above(value).times(Mag::pow2(-(prec as i64)))
    }

    pub fn sub(&self, o: &HpBounded, prec: usize) -> HpBounded {
        let value = self.value.sub(&o.value, prec, HP_RM);
        let error = self.error.plus(o.error).plus(Self::round_off(&value, prec));
        HpBounded::new(value, error)
    }

    pub fn add(&self, o: &HpBounded, prec: usize) -> HpBounded {
        let value = self.value.add(&o.value, prec, HP_RM);
        let error = self.error.plus(o.error).plus(Self::round_off(&value, prec));
        HpBounded::new(value, error)
    }

    pub fn mul(&self, o: &HpBounded, prec: usize) -> HpBounded {
        let value = self.value.mul(&o.value, prec, HP_RM);
        // `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, then the rounding of the product itself.
        let error = Mag::above(&self.value)
            .times(o.error)
            .plus(Mag::above(&o.value).times(self.error))
            .plus(self.error.times(o.error))
            .plus(Self::round_off(&value, prec));
        HpBounded::new(value, error)
    }

    /// Division by an **exact** nonzero divisor — what a wide frame's origin (`num / den`)
    /// realizes through. The dividend's radius scales by `1/|b| ≤ 2^(1−e_b)` (from
    /// `|b| ≥ 2^(e_b−1)`), and the division's own rounding is charged on top. The divisor
    /// being exact is a premise (its producer is [`Self::of_bigint`]), so it is asserted rather
    /// than handled.
    pub fn div_exact(&self, b: &HpBounded, prec: usize) -> Option<HpBounded> {
        debug_assert!(
            b.error.exp2().is_none(),
            "div_exact's divisor must carry a zero radius"
        );
        if b.value.is_zero() {
            return None;
        }
        let e = b.value.exponent()? as i64;
        let value = self.value.div(&b.value, prec, HP_RM);
        let error = self
            .error
            .times(Mag::pow2(1 - e))
            .plus(Self::round_off(&value, prec));
        Some(HpBounded::new(value, error))
    }

    /// **An `f64` at or above every `|x|` the interval holds** — the far end (the radius read as
    /// the power of two at or above it) rounded to its nearest `f64`, one step up: an exact value
    /// comes back one ulp above itself. A far end below the normal range reads as
    /// `f64::MIN_POSITIVE`; `None` past the largest `f64`.
    pub fn abs_upper_f64(&self) -> Option<f64> {
        let (lo, hi) = crate::realize::interval_ends(&self.value, self.error)?;
        let far = if lo.abs().cmp(&hi.abs())? >= 0 {
            lo.abs()
        } else {
            hi.abs()
        };
        match crate::realize::bigfloat_to_f64(&far) {
            Some(0.0) => Some(0.0),
            // A positive finite `f64`'s next one up is its bit pattern plus one.
            Some(f) => Some(f64::from_bits(f.to_bits() + 1)),
            None => far
                .exponent()
                .is_some_and(|e| e < 0)
                .then_some(f64::MIN_POSITIVE),
        }
    }

    /// **An `f64` at or below every `|x|` the interval holds, and above zero** — the near end (the
    /// radius read as for [`Self::abs_upper_f64`]) rounded to its nearest `f64`, one step down.
    /// `None` where the interval may reach zero, or its near end lies below the normal range.
    pub fn abs_lower_f64(&self) -> Option<f64> {
        let (lo, hi) = crate::realize::interval_ends(&self.value, self.error)?;
        if lo.is_zero() || hi.is_zero() || lo.is_positive() != hi.is_positive() {
            return None;
        }
        let near = if lo.abs().cmp(&hi.abs())? <= 0 {
            lo.abs()
        } else {
            hi.abs()
        };
        match crate::realize::bigfloat_to_f64(&near) {
            Some(f) => {
                // …and its next one down, the pattern minus one (`f > 0` here).
                let down = f64::from_bits(f.to_bits() - 1);
                (down >= f64::MIN_POSITIVE).then_some(down)
            }
            None => near.exponent().is_some_and(|e| e > 0).then_some(f64::MAX),
        }
    }

    /// `Some(true)` if definitely positive, `Some(false)` if definitely negative, `None` if the
    /// interval straddles 0.
    ///
    /// `None` is **"not decided at this precision"**, not "proved zero" — no finite precision
    /// proves a transcendental equality. The caller either climbs the ladder or says so.
    pub fn sign(&self) -> Option<bool> {
        let low = Mag::below(&self.value)?;
        self.error.lt(low).then(|| self.value.is_positive())
    }

    /// A positive lower bound on this interval's true magnitude, or `None` when it may reach
    /// zero — [`Mag::below`]'s exponent floor on the computed value, minus the radius, in the
    /// rounded-toward-zero direction [`Mag::minus`] exists for.
    fn mag_lo(&self) -> Option<Mag> {
        Mag::below(&self.value)?.minus(self.error)
    }

    /// Division by an **interval** divisor — what [`HpBounded::div_exact`] refuses, made sound by
    /// bounding the divisor away from zero first.
    ///
    /// This is a *realization* operation (a frame's origin is the foot of a perpendicular,
    /// `(−d·n)/(n·n)`, and for a judged plane `n·n` carries a radius); the kernel's no-division
    /// discipline is about **sign questions**, where a quotient would manufacture an irrational
    /// the predicate then has to trust. A realization hands back the radius it incurred, which is
    /// exactly what this does.
    ///
    /// The bound: with `q̂` the computed quotient and `L ≤ |b|` (true divisor),
    /// `a/b − â/b̂ = (a−â)/b + (â/b̂)·(b̂−b)/b`, so the true-input error is at most
    /// `(r_a + |â/b̂|·r_b) / L`, and `|â/b̂| ≤ |q̂| + round_off`. `None` when the divisor's
    /// interval may reach zero — the caller climbs or rejects by name, never guesses.
    pub fn div(&self, b: &HpBounded, prec: usize) -> Option<HpBounded> {
        let lo = b.mag_lo()?;
        let value = self.value.div(&b.value, prec, HP_RM);
        let ro = Self::round_off(&value, prec);
        let q = Mag::above(&value).plus(ro);
        let error = self.error.plus(q.times(b.error)).over(lo)?.plus(ro);
        Some(HpBounded::new(value, error))
    }

    /// `1/√x` of an **interval** — the judged-frame twin of [`crate::inv_sqrt_hp`] (exact
    /// `Rat` input) and [`crate::inv_sqrt_hp_big`] (exact `BigInt` input): same shape, an
    /// input that carries a radius.
    ///
    /// The input error's amplification comes from the derivative: `|d(1/√x)| = ½·x^(−3/2)`,
    /// largest at the interval's low end, so `error ≤ r / (2·L^(3/2))` with `L ≤ x` a positive
    /// lower bound. `L ≥ 2^(e−1)` from its exponent, so `2·L^(3/2) ≥ 2^(1+(e−1)+⌊(e−1)/2⌋)` —
    /// a pure power of two, which [`Mag::over`] divides by exactly. The two arithmetic roundings
    /// (sqrt, then the reciprocal) are charged on top.
    ///
    /// `None` when `x` may reach zero or is not positive — a degenerate normal has no direction,
    /// and the caller says so by name.
    pub fn inv_sqrt(&self, prec: usize) -> Option<HpBounded> {
        if !self.value.is_positive() {
            return None;
        }
        let lo = self.mag_lo()?;
        let e = lo.exp2()?;
        let s = self.value.sqrt(prec, HP_RM);
        let one = BigFloat::from_f64(1.0, prec);
        let value = one.div(&s, prec, HP_RM);
        let denom = Mag::pow2(1 + (e - 1) + (e - 1).div_euclid(2));
        let ro = Self::round_off(&value, prec);
        let error = self.error.over(denom)?.plus(ro).plus(ro);
        Some(HpBounded::new(value, error))
    }
}

#[cfg(test)]
#[path = "tests/bounded.rs"]
mod tests;
