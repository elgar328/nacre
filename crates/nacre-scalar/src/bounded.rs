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

    /// An arbitrary-precision integer as an **exact** interval — the wide-frame (S4) entry point.
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
    /// realizes through (S4). The dividend's radius scales by `1/|b| ≤ 2^(1−e_b)` (from
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
mod tests {
    use super::*;
    use crate::Angle;

    /// An upper `f64` reading of a [`Mag`], for comparisons in tests only.
    fn mag_f64(m: Mag) -> f64 {
        m.exp2().map_or(0.0, |e| 2f64.powi(e as i32))
    }

    /// |x − y| at the two values' own precision, as f64 — a test-side residual reader.
    fn resid(x: &BigFloat, y: &BigFloat, prec: usize) -> f64 {
        let d = x.sub(y, prec, HP_RM).abs();
        if d.is_zero() {
            return 0.0;
        }
        d.exponent().map_or(f64::INFINITY, |e| 2f64.powi(e))
    }

    /// ★ The differential the two new interval operations are locked by: on **exact** inputs they
    /// must agree with the exact routes that already exist, and on **interval** inputs the radius
    /// must contain the truth — probed at the interval's endpoints, which is where a monotone
    /// function's image is extreme (so the corners are the whole question, not a sample).
    #[test]
    fn the_interval_division_agrees_with_the_exact_route_and_contains_the_corners() {
        let prec = 192;
        let gt = 512;
        for (a, ra, b, rb) in [
            (3.75, 0.0, 1.5, 0.0),
            (-7.0, 1e-20, 0.3, 1e-21),
            (1e12, 1e-6, -2.5e-3, 1e-12),
            (0.1, 1e-30, 12345.678, 1e-18),
        ] {
            let av = HpBounded::new(big(a, prec), Mag::of(ra));
            let bv = HpBounded::new(big(b, prec), Mag::of(rb));
            let q = av.div(&bv, prec).expect("divisor clears zero");
            if ra == 0.0 && rb == 0.0 {
                // The exact route answers the same value, bit for bit — same BigFloat division.
                let qe = av.div_exact(&bv, prec).expect("exact divisor");
                assert!(
                    q.value.sub(&qe.value, prec, HP_RM).is_zero(),
                    "interval and exact division disagree on exact inputs"
                );
            }
            // Soundness at the corners: the true quotient for any (a', b') in the box lies
            // within value ± error. A quotient is monotone in each argument on a sign-constant
            // box, so the four corners bound the image.
            for (da, db) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                let at = big(a + da * ra, gt);
                let bt = big(b + db * rb, gt);
                let truth = at.div(&bt, gt, HP_RM);
                let err = resid(&q.value, &truth, gt);
                assert!(
                    err <= mag_f64(q.error) || err == 0.0,
                    "corner ({da},{db}) of {a}±{ra} / {b}±{rb}: off by {err:e}, radius {:e}",
                    mag_f64(q.error)
                );
            }
        }
    }

    #[test]
    fn the_interval_inv_sqrt_agrees_with_the_exact_route_and_contains_the_endpoints() {
        let prec = 192;
        let gt = 512;
        for (x, r) in [(2.0, 0.0), (0.09, 1e-22), (1e20, 1.0), (5.0e-7, 1e-27)] {
            let xv = HpBounded::new(big(x, prec), Mag::of(r));
            let s = xv.inv_sqrt(prec).expect("bounded away from zero");
            if r == 0.0 {
                // Against the exact-rational route on a value both can state.
                let rat = Rat::from_decimal(x).expect("in the window");
                let exact = crate::inv_sqrt_bounded(rat, prec).expect("positive");
                let err = resid(&s.value, &exact.value, prec);
                assert!(
                    err <= mag_f64(s.error) + mag_f64(exact.error),
                    "exact input {x}: the two routes disagree by {err:e}"
                );
            }
            // 1/√x is monotone, so the endpoints are the extreme truths.
            for d in [1.0, -1.0] {
                let xt = big(x + d * r, gt);
                let truth = big(1.0, gt).div(&xt.sqrt(gt, HP_RM), gt, HP_RM);
                let err = resid(&s.value, &truth, gt);
                assert!(
                    err <= mag_f64(s.error) || err == 0.0,
                    "endpoint {d} of {x}±{r}: off by {err:e}, radius {:e}",
                    mag_f64(s.error)
                );
            }
        }
    }

    /// The refusals: a divisor or radicand whose interval may reach zero is `None`, never a
    /// guessed bound — and a negative radicand is refused outright.
    #[test]
    fn a_quantity_that_may_reach_zero_is_refused_not_bounded() {
        let prec = 128;
        let wide = HpBounded::new(big(1e-10, prec), Mag::of(1.0)); // straddles zero
        let a = HpBounded::exact(big(1.0, prec));
        assert!(a.div(&wide, prec).is_none(), "divisor may reach zero");
        assert!(wide.inv_sqrt(prec).is_none(), "radicand may reach zero");
        let neg = HpBounded::exact(big(-4.0, prec));
        assert!(neg.inv_sqrt(prec).is_none(), "a negative has no real root");
    }

    fn big(x: f64, prec: usize) -> BigFloat {
        BigFloat::from_f64(x, prec)
    }

    /// **Does the error `Angle::realization_error_of` reports actually bound the real one?**
    ///
    /// `WitnessPoint::rotate_about` charges that number for the one input here without a rounding contract:
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
    /// What it does not do: it checks one input. `nacre-cip`'s `tol_bounds_error_over_random_chains`
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
        // ★ **Reported to the nearest ε, which needs a comparison rather than an octave reader.**
        // `2^exponent` is an upper *octave*, up to 2× above the true magnitude — so a figure read
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
            let (hc, hs) = a.cos_sin_bounded(GT);
            for (which, f, h, reported) in [("cos", c, &hc, dc), ("sin", s, &hs, ds)] {
                // |f64 − true| ≤ |f64 − deep midpoint| + that realization's own radius.
                let truth = big(f, GT).sub(&h.value, GT, HP_RM).abs().add(
                    &big(h.error.exp2().map_or(0.0, |e| 2f64.powi(e as i32)), GT),
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
        let a = HpBounded::new(big(1.0e30, prec), error);
        let b = HpBounded::new(
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
        let a = HpBounded::new(big(3.0, prec), Mag::pow2(-100));
        let b = HpBounded::new(big(2.0, prec), Mag::pow2(-100));
        assert_eq!(a.sub(&b, prec).sign(), Some(true));
        assert_eq!(b.sub(&a, prec).sign(), Some(false));
        assert_eq!(a.mul(&b, prec).sign(), Some(true));
    }

    /// The radius survives a rung deep enough to flush an `f64` radius to zero — the reason
    /// [`Mag`] exists. With an `f64` radius this product would report a confident sign.
    #[test]
    fn a_deep_rung_does_not_lose_the_radius() {
        let prec = 2048;
        let tiny = HpBounded::new(big(1.0, prec), Mag::pow2(-(prec as i64)));
        let p = tiny.mul(&tiny, prec);
        assert!(!p.error.is_zero(), "the radius vanished at {prec} bits");
        // A difference of exactly that size is therefore undecided, not positive.
        let q = HpBounded::new(
            big(1.0, prec).add(&BigFloat::from_f64(1.0, prec), prec, HP_RM),
            Mag::pow2(-(prec as i64) + 4),
        );
        let r = HpBounded::new(big(2.0, prec), Mag::pow2(-(prec as i64) + 4));
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
        let a = HpBounded::exact(big(0.5, prec));
        let b = HpBounded::exact(big(0.25, prec));
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

    /// The one place the zero guard shows: rounding an exact zero costs nothing, so a product or
    /// sum that lands exactly on zero from exact inputs reports the radius it earned — none.
    /// [`Mag::above`] answers `ZERO` there where a guard-less exponent reader would not.
    #[test]
    fn an_exact_zero_is_charged_no_rounding() {
        let prec = 128;
        let zero = HpBounded::exact(big(0.0, prec));
        let one = HpBounded::exact(big(1.0, prec));
        assert!(
            zero.mul(&one, prec).error.is_zero(),
            "0·1 charged a rounding"
        );
        assert!(
            one.sub(&one, prec).error.is_zero(),
            "1−1 charged a rounding"
        );
        assert_eq!(HpBounded::round_off(&big(0.0, prec), prec), Mag::ZERO);
    }
}
