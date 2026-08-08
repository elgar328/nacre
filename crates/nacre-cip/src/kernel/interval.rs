//! The two intervals the judgment ladder runs on — **one machine at two precisions**.
//!
//! [`Approx`] is the `f64` filter: a value with an error radius, propagated through every operation,
//! whose sign is only reported when the interval clears zero. [`HpApprox`] is the same thing in
//! astro-float, differing in exactly one place — the per-operation rounding is `2⁻ᵖʳᵉᶜ` instead of
//! `2·ε`.
//!
//! **They live in one file on purpose.** The two used to be different kinds of machine: the filter
//! propagated a radius, while the escalation compared a determinant against a hand-picked floor
//! (`FLOOR_K · mag · 2⁻ᵖʳᵉᶜ`) built from a separately-maintained magnitude estimate. One of those
//! estimates was read off a value that had already cancelled, the floor collapsed with it, and a
//! rounding residue was reported as a confident sign. Side by side, an operation that propagates
//! its radius in one and not the other is visible.
//!
//! **The radius is a [`Mag`], not an `f64`.** At a deep rung `2⁻ᵖʳᵉᶜ` underflows an `f64` to
//! zero, and a zero radius claims exactness — the same failure in new clothes. See
//! [`nacre_scalar::bound`].

use astro_float::BigFloat;
use nacre_scalar::{Mag, Rat};

use super::HP_RM;

// ---- realizing an exact value, and what that realization costs ----
//
// These four moved here when the 2D frame that first held them was deleted (nothing consumed it).
// They belong with the intervals: three of them are how an exact rational *enters* this machine,
// and the fourth is the f64 rounding the filter charges per operation.

/// A rational as an arbitrary-precision float.
///
/// **Numerator and denominator go in as `i128`, not through `f64`.** Routing them through `f64`
/// costs a relative `2⁻⁵³` the moment either exceeds 2⁵³ — and *no amount of working precision
/// recovers it*, because the loss happens before astro-float ever sees the value. Measured: a
/// denominator of `2⁵⁴+1` took a coordinate's error bound from `2⁻²⁵¹` to `2⁻⁵⁰`. Chained exact
/// rational arithmetic multiplies denominators, so that is not an exotic input; it is what the
/// crate's own bit-growth note is about.
pub(crate) fn rat_to_big(r: Rat, prec: usize) -> BigFloat {
    // `from_i128` needs at least 128 bits to hold the integer before the division rounds it.
    let p = prec.max(128);
    BigFloat::from_i128(r.numer(), p).div(&BigFloat::from_i128(r.denom(), p), prec, HP_RM)
}

/// A rational realized at `prec` bits **with the error that realization carries**.
///
/// Two things can go wrong and both are bounded here rather than assumed away:
///
/// The integers themselves are exact — [`rat_to_big`] feeds them in as `i128` — so the only
/// error is the division, and even that vanishes when it terminates: a power-of-two denominator
/// with a numerator inside `prec` bits is exact, and then the radius is genuinely zero.
pub(crate) fn rat_to_hp(r: Rat, prec: usize) -> HpApprox {
    let value = rat_to_big(r, prec);
    let (n, d) = (r.numer(), r.denom()); // `d > 0` after reduction
    let n_bits = (128 - n.unsigned_abs().leading_zeros()) as usize;
    if d & (d - 1) == 0 && n_bits <= prec {
        return HpApprox::exact(value); // a terminating division: a power-of-two denominator
    }
    // The integers enter exactly (see `rat_to_big`), so the only error left is the division's
    // own rounding.
    let error = ub(&value).times(Mag::pow2(-(prec as i64)));
    HpApprox::new(value, error)
}

/// An arbitrary-precision integer as an **exact** interval — the wide-frame (S4) entry point.
///
/// [`nacre_scalar::bigint_to_bigfloat`] converts at the integer's own bit length, so the value
/// enters whole and the radius is genuinely zero; downstream operations charge their own
/// rounding, exactly as [`rat_to_hp`]'s exact branch does.
pub(crate) fn bigint_to_hp(x: &num_bigint::BigInt, prec: usize) -> HpApprox {
    HpApprox::exact(nacre_scalar::bigint_to_bigfloat(x, prec))
}

/// Magnitude of a `BigFloat` as an f64 power of two (0 when exactly zero).
pub(crate) fn bf_mag(bf: &BigFloat) -> f64 {
    if bf.is_zero() {
        0.0
    } else {
        2f64.powi(bf.exponent().unwrap_or(0))
    }
}

/// A value with a symmetric error radius (`value ± error`, `error ≥ 0`). Arithmetic keeps
/// `error` a sound upper bound (worst case), plus a per-op f64-rounding inflation.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Approx {
    pub value: f64,
    pub error: f64,
}

impl Approx {
    pub fn new(value: f64, error: f64) -> Self {
        Approx { value, error }
    }
    /// `2·ε` per operation, four times the `ε/2` that round-to-nearest can cost — derived, not
    /// picked, and the same charge in `add` and `mul`.
    pub fn sub(self, o: Approx) -> Approx {
        let value = self.value - o.value;
        Approx::new(
            value,
            self.error + o.error + 2.0 * f64::EPSILON * value.abs(),
        )
    }
    pub fn add(self, o: Approx) -> Approx {
        let value = self.value + o.value;
        Approx::new(
            value,
            self.error + o.error + 2.0 * f64::EPSILON * value.abs(),
        )
    }
    pub fn mul(self, o: Approx) -> Approx {
        let value = self.value * o.value;
        // Worst-case product radius `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, plus the
        // f64 rounding of the product itself.
        let error = self.value.abs() * o.error + o.value.abs() * self.error + self.error * o.error;
        Approx::new(value, error + 2.0 * f64::EPSILON * value.abs())
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

/// An upper bound on `|x|`, from its exponent (`|x| < 2^exponent`).
pub(crate) fn ub(x: &BigFloat) -> Mag {
    match x.exponent() {
        Some(e) if !x.is_zero() => Mag::pow2(e as i64),
        _ => Mag::ZERO,
    }
}

/// A lower bound on `|x|` (`|x| ≥ 2^(exponent−1)`), or `None` when `x` is zero.
///
/// Reading the mantissa would tighten this by up to one bit. It is left at the exponent because
/// the error is in the safe direction — a judgement declines slightly sooner than it must, and
/// the ladder is what recovers those, not a tighter comparison here.
pub(crate) fn lb(x: &BigFloat) -> Option<Mag> {
    match x.exponent() {
        Some(e) if !x.is_zero() => Some(Mag::pow2(e as i64 - 1)),
        _ => None,
    }
}

/// [`Approx`] in astro-float: a high-precision value with a **computed** error radius.
///
/// Every constructor must supply a radius that actually bounds its value's distance from the
/// truth — for a coordinate that means the rotation chain's realization error
/// ([`nacre_scalar::Angle::cos_sin_bounded`] is where it starts), and for a rational read at
/// enough bits it means [`Mag::ZERO`]. A radius invented for convenience makes every sign above
/// it unearned.
#[derive(Clone, Debug)]
pub(crate) struct HpApprox {
    pub value: BigFloat,
    pub error: Mag,
}

impl HpApprox {
    pub fn new(value: BigFloat, error: Mag) -> Self {
        HpApprox { value, error }
    }

    /// A value known exactly at this precision — a rational whose realization did not round.
    pub fn exact(value: BigFloat) -> Self {
        HpApprox {
            value,
            error: Mag::ZERO,
        }
    }

    /// The rounding a `prec`-bit operation adds to its own result: at most a half-ulp,
    /// `|result| · 2⁻ᵖʳᵉᶜ`.
    fn round_off(value: &BigFloat, prec: usize) -> Mag {
        ub(value).times(Mag::pow2(-(prec as i64)))
    }

    pub fn sub(&self, o: &HpApprox, prec: usize) -> HpApprox {
        let value = self.value.sub(&o.value, prec, HP_RM);
        let error = self.error.plus(o.error).plus(Self::round_off(&value, prec));
        HpApprox::new(value, error)
    }

    pub fn add(&self, o: &HpApprox, prec: usize) -> HpApprox {
        let value = self.value.add(&o.value, prec, HP_RM);
        let error = self.error.plus(o.error).plus(Self::round_off(&value, prec));
        HpApprox::new(value, error)
    }

    pub fn mul(&self, o: &HpApprox, prec: usize) -> HpApprox {
        let value = self.value.mul(&o.value, prec, HP_RM);
        // `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, then the rounding of the product itself.
        let error = ub(&self.value)
            .times(o.error)
            .plus(ub(&o.value).times(self.error))
            .plus(self.error.times(o.error))
            .plus(Self::round_off(&value, prec));
        HpApprox::new(value, error)
    }

    /// Division by an **exact** nonzero divisor — what a wide frame's origin (`num / den`)
    /// realizes through (S4). The dividend's radius scales by `1/|b| ≤ 2^(1−e_b)` (from
    /// `|b| ≥ 2^(e_b−1)`), and the division's own rounding is charged on top. The divisor
    /// being exact is a premise (its producer is [`bigint_to_hp`]), so it is asserted rather
    /// than handled.
    pub fn div_exact(&self, b: &HpApprox, prec: usize) -> Option<HpApprox> {
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
        Some(HpApprox::new(value, error))
    }

    /// `Some(true)` if definitely positive, `Some(false)` if definitely negative, `None` if the
    /// interval straddles 0.
    ///
    /// `None` is **"not decided at this precision"**, not "proved zero" — no finite precision
    /// proves a transcendental equality. The caller either climbs the ladder or says so.
    pub fn sign(&self) -> Option<bool> {
        let low = lb(&self.value)?;
        self.error.lt(low).then(|| self.value.is_positive())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_scalar::{Angle, Mag};

    fn big(x: f64, prec: usize) -> BigFloat {
        BigFloat::from_f64(x, prec)
    }

    /// **Does the error `Angle::realization_error_of` reports actually bound the real one?**
    ///
    /// `Pt3::rotate_about` charges that number for the one input here without a rounding contract:
    /// neither Rust nor any libm promises an accuracy for `f64::cos`. It used to charge a constant
    /// measured once and written into a doc — sound only on machines like the one it was taken on,
    /// which a kernel that ships to browsers cannot assume. Now it measures, so **a worse platform
    /// reports a bigger number and the tolerance grows to match** rather than silently under-stating.
    ///
    /// That moves the question. There is no budget left to overrun; what has to hold is that the
    /// measurement is an *upper* bound. So this compares production's answer — taken at 128 bits —
    /// against a realization 384 bits deeper.
    ///
    /// **Not `#[ignore]`d, and that is the point.** The slow astro-float suites are skipped because
    /// they are slow; this one exists so that *every build checks its own platform*, and a few
    /// hundred angles is milliseconds.
    ///
    /// ★★★ **It asks about the pair, because that is what the answer is about.** The values go in
    /// as arguments — the same ones `rotate_about` writes into `coord` — rather than being re-derived
    /// inside. Written the other way, this test asked whether the error of *some* realization bounded
    /// the error of *another*, and it failed at 27°: a debug and a release build realize `sin 27°`
    /// one ulp apart, and so do two call sites within one release build. That divergence is a real
    /// and separate defect in the kernel's determinism story; what this test now checks is the
    /// property `tol` actually needs, which is about the pair in hand.
    ///
    /// ★ **The 90°-family is included, not skipped, and that is the load-bearing part.** There
    /// `realization_error_of` claims *exactly zero*, and `rotate_about` reads that zero to contribute
    /// no tolerance at all — which is what keeps a quadrantal origin rotation at tol 0 and keeps
    /// axis-aligned models on the exact predicate path. A claimed zero that was not really zero
    /// would not fail loudly; it would quietly move those models. Here it has to be real.
    ///
    /// ★★ **It also still reports the worst error, though nothing depends on it.** Charging a
    /// measured value means a degraded platform no longer fails anything — it just escalates more
    /// and runs slower, invisibly. This line is what keeps that visible.
    ///
    /// What it does not do: it checks one input. `frame3::tol_bounds_error_over_random_chains`
    /// checks the conclusion — that `tol` bounds the error, second-order terms and all.
    #[test]
    fn the_measured_trig_error_bounds_the_real_one() {
        const GT: usize = 512; // the reference it is checked against
        let mut st = 0x9E37_79B9_7F4A_7C15u64;
        let mut rng = |lo: i128, hi: i128| {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            lo + ((st >> 33) as i128) % (hi - lo).max(1)
        };
        let mut angles: Vec<(i128, i128)> = (0..360).map(|d| (d, 1)).collect();
        for _ in 0..200 {
            angles.push((rng(0, 360_000), rng(1, 9973)));
            // ★ Near an axis, where the realized cos or sin is close to zero — the region that
            // refuted the tangential tol model, and where a relative error reads worst.
            angles.push((90 * rng(0, 4) * 1000 + rng(-3, 4), 1000));
            // ★ Past `2⁵³`, where `Rat::to_f64` leaves its fast path. The radian argument is what
            // dominates the error, so an angle whose conversion takes the other route has to be here.
            angles.push((rng(1, 1 << 62) * 360, rng(1 << 53, 1 << 62)));
        }
        // ★ **Reported to the nearest ε, which needs a comparison rather than `bf_mag`.** `bf_mag`
        // is `2^exponent` — an upper *octave*, up to 2× above the true magnitude — so a figure read
        // from it blurs exactly the range that matters. The smallest integer `k` with `x < k·ε` is
        // one comparison per candidate and says the thing plainly.
        //
        // ★★★ **The range starts at 1, and it must**: `0·ε` is `0.0`, and astro-float's `cmp`
        // answers `Some(-1)` in *both* directions against a zero — so `x.cmp(&zero)` reports "less
        // than" for a positive `x`, every inexact realization reads `k = 0`, and `worst_k` is then
        // only ever raised by the *exact* ones, which read the floor. That prints a headline of
        // "worst < 1ε" for a platform whose real figure is several times that. A `k` of 0 would be
        // meaningless anyway — nothing is below `0·ε` — which is what makes the widening look free.
        let eps_multiple = |x: &BigFloat| -> usize {
            (1..=64)
                .find(|k| {
                    // astro-float's `cmp` yields a sign as `Option<i128>`, not an `Ordering`.
                    x.abs()
                        .cmp(&big(*k as f64 * f64::EPSILON, GT))
                        .is_some_and(|sign| sign < 0)
                })
                .unwrap_or(usize::MAX)
        };
        // ★★★ **The reading is pinned before it is trusted.** Three separate figures in this
        // kernel's error work turned out to be artefacts of how they were printed rather than of
        // what was measured, and each time the tell was a number that would not move. A reading
        // that cannot be shown to respond to a known input is not evidence.
        assert_eq!(eps_multiple(&big(0.0, GT)), 1, "zero reads as the floor");
        assert_eq!(eps_multiple(&big(3.5 * f64::EPSILON, GT)), 4);
        assert_eq!(eps_multiple(&big(-9.0 * f64::EPSILON, GT)), 10);
        assert_eq!(eps_multiple(&big(64.0 * f64::EPSILON, GT)), usize::MAX);
        let (mut worst_k, mut worst_at) = (0usize, (0i128, 1i128));
        let (mut checked, mut exact_seen) = (0usize, 0usize);
        for (n, d) in angles {
            let Some(a) = Rat::new(n, d).and_then(Angle::from_deg) else {
                continue; // outside Rat's range — `from_deg` declines, and so does the kernel
            };
            let (c, s) = a.cos_sin_f64();
            let (dc, ds) = a.realization_error_of(c, s);
            if dc == 0.0 && ds == 0.0 {
                exact_seen += 1;
            }
            let (hc, hs, rc, rs) = a.cos_sin_bounded(GT);
            for (which, f, h, error, reported) in [("cos", c, &hc, rc, dc), ("sin", s, &hs, rs, ds)]
            {
                // |f64 − true| ≤ |f64 − deep midpoint| + that realization's own radius.
                let truth = big(f, GT).sub(h, GT, HP_RM).abs().add(
                    &big(error.exp2().map_or(0.0, |e| 2f64.powi(e as i32)), GT),
                    GT,
                    HP_RM,
                );
                assert!(
                    truth.cmp(&big(reported, GT)).is_some_and(|sign| sign <= 0),
                    "realization_error_of under-states at {n}/{d} deg [{which}]: \
                     f = {f:e}, reported {reported:e}, true error 2^{:?}",
                    truth.exponent(),
                );
                // ★ The worst is tracked over the *inexact* realizations only. An exact one has
                // error zero, which reads as the floor `k = 1` — and since the corpus opens at 0°,
                // including them let the family claim the headline and named an angle with no
                // error at all as this platform's worst case.
                let k = eps_multiple(&truth);
                if reported != 0.0 && k > worst_k {
                    (worst_k, worst_at) = (k, (n, d));
                }
                checked += 1;
            }
        }
        assert!(checked > 700, "corpus shrank to {checked} realizations");
        assert!(
            exact_seen > 0,
            "no quadrantal angle reached the zero-error claim — the strongest case went untested"
        );
        eprintln!(
            "[trig realization] worst error < {worst_k}ε (at {}/{} deg), {exact_seen} exact, \
             over {checked} realizations",
            worst_at.0, worst_at.1,
        );
    }

    /// The point of the type: a value that is only rounding residue must not report a sign, no
    /// matter how confidently its `value` is nonzero.
    #[test]
    fn a_residue_left_by_cancellation_reports_no_sign() {
        let prec = 200;
        // Two large values differing by far less than the radius they carry. The separation has
        // to stay inside `prec` bits of the operands, or the subtraction is exactly zero and the
        // test proves nothing.
        let error = Mag::of(1.0e-10);
        let a = HpApprox::new(big(1.0e30, prec), error);
        let b = HpApprox::new(
            big(1.0e30, prec).sub(&big(1.0e-20, prec), prec, HP_RM),
            error,
        );
        let d = a.sub(&b, prec);
        assert!(
            !d.value.is_zero(),
            "the midpoint is nonzero — that is the trap"
        );
        assert_eq!(
            d.sign(),
            None,
            "a difference under its own radius claimed a sign"
        );
    }

    /// …and a value genuinely clear of its radius still decides, so the type is not merely
    /// refusing to answer.
    #[test]
    fn a_value_clear_of_its_radius_still_decides() {
        let prec = 200;
        let a = HpApprox::new(big(3.0, prec), Mag::pow2(-100));
        let b = HpApprox::new(big(2.0, prec), Mag::pow2(-100));
        assert_eq!(a.sub(&b, prec).sign(), Some(true));
        assert_eq!(b.sub(&a, prec).sign(), Some(false));
        assert_eq!(a.mul(&b, prec).sign(), Some(true));
    }

    /// The radius survives a rung deep enough to flush an `f64` radius to zero — the reason
    /// [`Mag`] exists. With an `f64` radius this product would report a confident sign.
    #[test]
    fn a_deep_rung_does_not_lose_the_radius() {
        let prec = 2048;
        let tiny = HpApprox::new(big(1.0, prec), Mag::pow2(-(prec as i64)));
        let p = tiny.mul(&tiny, prec);
        assert!(!p.error.is_zero(), "the radius vanished at {prec} bits");
        // A difference of exactly that size is therefore undecided, not positive.
        let q = HpApprox::new(
            big(1.0, prec).add(&BigFloat::from_f64(1.0, prec), prec, HP_RM),
            Mag::pow2(-(prec as i64) + 4),
        );
        let r = HpApprox::new(big(2.0, prec), Mag::pow2(-(prec as i64) + 4));
        assert_eq!(q.sub(&r, prec).sign(), None);
    }

    /// An exact input contributes no radius of its own, so the only radius an exact computation
    /// carries is the half-ulp charged per operation.
    ///
    /// That charge is unconditional: astro-float does not report whether an operation was exact,
    /// so a subtraction of two dyadics — exact in fact — is still charged. Conservative in the
    /// safe direction, and negligible in size, but it does mean **this path never proves zero**.
    /// Judgements that need a proved zero take the exact route above this kernel instead.
    #[test]
    fn exact_inputs_keep_a_radius_only_from_rounding() {
        let prec = 200;
        let a = HpApprox::exact(big(0.5, prec));
        let b = HpApprox::exact(big(0.25, prec));
        assert!(
            a.error.is_zero() && b.error.is_zero(),
            "an exact input carried a radius"
        );
        let d = a.sub(&b, prec);
        assert!(
            d.error.lt(Mag::pow2(-190)),
            "an exact subtraction picked up more than a half-ulp: 2^{:?}",
            d.error.exp2()
        );
        assert_eq!(d.sign(), Some(true));
        // Exactly equal inputs are the case that cannot be settled here: the difference is zero
        // and zero has no sign to read.
        assert_eq!(a.sub(&a, prec).sign(), None);
    }
}
