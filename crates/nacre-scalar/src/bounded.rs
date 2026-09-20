//! **A value and the error it is caught in** — one machine at two precisions.
//!
//! [`Bounded`] is the `f64` filter: a value with an error radius, propagated through every
//! operation, whose sign is only reported when the interval clears zero. [`HpBounded`] is the same
//! thing in astro-float, differing in exactly one place — the per-operation rounding is `2⁻ᵖʳᵉᶜ`
//! instead of `2·ε`. Named for what the type promises — the truth lies within `value ± error` —
//! not for what the value is (an approximation), because the promise is the whole point.
//!
//! **They live in one file on purpose.** The two used to be different kinds of machine: the filter
//! propagated a radius, while the escalation compared a determinant against a hand-picked floor
//! (`FLOOR_K · mag · 2⁻ᵖʳᵉᶜ`) built from a separately-maintained magnitude estimate. One of those
//! estimates was read off a value that had already cancelled, the floor collapsed with it, and a
//! rounding residue was reported as a confident sign. Side by side, an operation that propagates
//! its radius in one and not the other is visible.
//!
//! **And they live in this crate, not in the judge that runs them,** for the same reason. The
//! judge (`nacre-cip`) used to hold both as crate-private types while this crate realized curved
//! coordinates with a second, looser spelling of the high-precision arithmetic — a tuple alias
//! and five free functions whose magnitude reader had no zero guard and whose rational entry had
//! no exact branch. Two spellings of one machine in two crates is how one of them drifts; one
//! spelling here, and the judge imports it.
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
pub fn rat_to_big(r: Rat, prec: usize) -> BigFloat {
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

// `add`/`sub`/`mul` are spelled like the std traits on purpose and are not them: the same names
// at both precisions is what lets the two machines be read side by side, and `HpBounded`'s take
// a `prec` the trait signatures have no room for. Naming these apart would split the one machine.
#[allow(clippy::should_implement_trait)]
impl Bounded {
    pub fn new(value: f64, error: f64) -> Self {
        Bounded { value, error }
    }
    /// `2·ε` per operation, four times the `ε/2` that round-to-nearest can cost — derived, not
    /// picked, and the same charge in `add` and `mul`.
    pub fn sub(self, o: Bounded) -> Bounded {
        let value = self.value - o.value;
        Bounded::new(
            value,
            self.error + o.error + 2.0 * f64::EPSILON * value.abs(),
        )
    }
    pub fn add(self, o: Bounded) -> Bounded {
        let value = self.value + o.value;
        Bounded::new(
            value,
            self.error + o.error + 2.0 * f64::EPSILON * value.abs(),
        )
    }
    pub fn mul(self, o: Bounded) -> Bounded {
        let value = self.value * o.value;
        // Worst-case product radius `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, plus the
        // f64 rounding of the product itself.
        let error = self.value.abs() * o.error + o.value.abs() * self.error + self.error * o.error;
        Bounded::new(value, error + 2.0 * f64::EPSILON * value.abs())
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
/// ([`crate::Angle::cos_sin_bounded`] is where it starts), and for a rational read at enough
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
    /// The integers themselves are exact — [`rat_to_big`] feeds them in as `i128` — so the only
    /// error is the division, and even that vanishes when it terminates: a power-of-two
    /// denominator with a numerator inside `prec` bits is exact, and then the radius is genuinely
    /// zero.
    pub fn of_rat(r: Rat, prec: usize) -> Self {
        let value = rat_to_big(r, prec);
        let (n, d) = (r.numer(), r.denom()); // `d > 0` after reduction
        let n_bits = (128 - n.unsigned_abs().leading_zeros()) as usize;
        if d & (d - 1) == 0 && n_bits <= prec {
            return HpBounded::exact(value); // a terminating division: a power-of-two denominator
        }
        // The integers enter exactly (see `rat_to_big`), so the only error left is the division's
        // own rounding.
        let error = Self::round_off(&value, prec);
        HpBounded::new(value, error)
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

    /// `1/√x` of an **interval** — the judged-frame twin of [`crate::inv_sqrt_bounded`] (exact
    /// `Rat` input) and [`crate::inv_sqrt_bigint_bounded`] (exact `BigInt` input): same shape, an
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
