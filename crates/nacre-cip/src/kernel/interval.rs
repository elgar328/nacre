//! The two intervals the judgment ladder runs on — **one machine at two precisions**.
//!
//! [`Iv`] is the `f64` filter: a value with an error radius, propagated through every operation,
//! whose sign is only reported when the interval clears zero. [`HpIv`] is the same thing in
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
//! **The radius is a [`Bound`], not an `f64`.** At a deep rung `2⁻ᵖʳᵉᶜ` underflows an `f64` to
//! zero, and a zero radius claims exactness — the same failure in new clothes. See
//! [`nacre_scalar::bound`].

use astro_float::BigFloat;
use nacre_scalar::{Bound, Rat};

use super::HP_RM;

// ---- realizing an exact value, and what that realization costs ----
//
// These four moved here when the 2D frame that first held them was deleted (nothing consumed it).
// They belong with the intervals: three of them are how an exact rational *enters* this machine,
// and the fourth is the f64 rounding the filter charges per operation.

/// f64 trig-realization error bound used by the directional tol formula — a conservative multiple
/// of ulp covering cos/sin rounding, the deg→rad conversion, and the combining arithmetic.
///
/// **The one constant here that cannot be derived, and it is measured instead.** Everything else
/// in this module charges round-to-nearest, which is contracted: `≤ ε/2` relative per operation.
/// `f64::cos` has no such contract — neither Rust nor the platform libm promises an accuracy —
/// so the only sound basis is measurement plus margin. Erring loose costs an escalation now and
/// then; erring tight is a wrong sign on a platform whose trig is a little worse than this one's.
///
/// ★★ **And whether it holds *here* is checked, not assumed** —
/// [`this_platforms_trig_stays_inside_the_charged_budget`] runs on every build and compares this
/// platform's realization against arbitrary-precision truth, angle by angle. That test also prints
/// the margin (`[trig budget]`), which is where the figure belongs: one pasted into this doc went
/// stale the moment the corpus it described grew, and stayed wrong for two commits.
///
/// **The budget it checks is `DA_F64 − ε`, not `DA_F64`** — for an origin pivot `Pt3::rotate_about`
/// has no `piv` term, so this one also carries the arithmetic of `fl(u·c) − fl(v·s)`. The
/// derivation is written out there.
///
/// ★ **And the dominant term is not `f64::cos` itself.** `Angle`'s f64 route is
/// `(deg.to_f64() * PI / 180.0).cos()`, so the *argument* already carries ~3ε relative from three
/// roundings and a rounded `PI`; with `|θ| ≤ 2π` that is ~9ε absolute, and `d cos = −sin·dθ`
/// carries it straight through. **The deg→rad conversion dominates, and libm is the small part** —
/// so "measure libm per angle and tighten this" is aimed at the wrong term *if the goal is a
/// tighter bound*. If the goal is **soundness**, measuring the realization per angle is exactly
/// right: it captures the conversion and libm together, and the kernel would then adapt to a bad
/// platform instead of assuming a good one. Removing the constant outright means realizing cos/sin
/// at arbitrary precision and rounding *that* to f64 (≤ ε/2, contracted), which moves every stored
/// rotated coordinate and needs its own cell.
pub(crate) const DA_F64: f64 = 16.0 * f64::EPSILON;

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
pub(crate) fn rat_to_hp(r: Rat, prec: usize) -> HpIv {
    let mid = rat_to_big(r, prec);
    let (n, d) = (r.numer(), r.denom()); // `d > 0` after reduction
    let n_bits = (128 - n.unsigned_abs().leading_zeros()) as usize;
    if d & (d - 1) == 0 && n_bits <= prec {
        return HpIv::exact(mid); // a terminating division: a power-of-two denominator
    }
    // The integers enter exactly (see `rat_to_big`), so the only error left is the division's
    // own rounding.
    let rad = ub(&mid).times(Bound::pow2(-(prec as i64)));
    HpIv::new(mid, rad)
}

/// Magnitude of a `BigFloat` as an f64 power of two (0 when exactly zero).
pub(crate) fn bf_mag(bf: &BigFloat) -> f64 {
    if bf.is_zero() {
        0.0
    } else {
        2f64.powi(bf.exponent().unwrap_or(0))
    }
}

/// A value with a symmetric error radius (`mid ± rad`, `rad ≥ 0`). Arithmetic keeps
/// `rad` a sound upper bound (worst case), plus a per-op f64-rounding inflation.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Iv {
    pub mid: f64,
    pub rad: f64,
}

impl Iv {
    pub fn new(mid: f64, rad: f64) -> Self {
        Iv { mid, rad }
    }
    /// `2·ε` per operation, four times the `ε/2` that round-to-nearest can cost — derived, not
    /// picked, and the same charge in `add` and `mul`.
    pub fn sub(self, o: Iv) -> Iv {
        let mid = self.mid - o.mid;
        Iv::new(mid, self.rad + o.rad + 2.0 * f64::EPSILON * mid.abs())
    }
    pub fn add(self, o: Iv) -> Iv {
        let mid = self.mid + o.mid;
        Iv::new(mid, self.rad + o.rad + 2.0 * f64::EPSILON * mid.abs())
    }
    pub fn mul(self, o: Iv) -> Iv {
        let mid = self.mid * o.mid;
        // Worst-case product radius `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, plus the
        // f64 rounding of the product itself.
        let rad = self.mid.abs() * o.rad + o.mid.abs() * self.rad + self.rad * o.rad;
        Iv::new(mid, rad + 2.0 * f64::EPSILON * mid.abs())
    }
    /// `Some(true)` if definitely positive, `Some(false)` if definitely negative,
    /// `None` if the interval straddles 0 (escalate).
    pub fn sign(self) -> Option<bool> {
        if self.mid > self.rad {
            Some(true)
        } else if self.mid < -self.rad {
            Some(false)
        } else {
            None
        }
    }
}

/// An upper bound on `|x|`, from its exponent (`|x| < 2^exponent`).
pub(crate) fn ub(x: &BigFloat) -> Bound {
    match x.exponent() {
        Some(e) if !x.is_zero() => Bound::pow2(e as i64),
        _ => Bound::ZERO,
    }
}

/// A lower bound on `|x|` (`|x| ≥ 2^(exponent−1)`), or `None` when `x` is zero.
///
/// Reading the mantissa would tighten this by up to one bit. It is left at the exponent because
/// the error is in the safe direction — a judgement declines slightly sooner than it must, and
/// the ladder is what recovers those, not a tighter comparison here.
pub(crate) fn lb(x: &BigFloat) -> Option<Bound> {
    match x.exponent() {
        Some(e) if !x.is_zero() => Some(Bound::pow2(e as i64 - 1)),
        _ => None,
    }
}

/// [`Iv`] in astro-float: a high-precision value with a **computed** error radius.
///
/// Every constructor must supply a radius that actually bounds its value's distance from the
/// truth — for a coordinate that means the rotation chain's realization error
/// ([`nacre_scalar::Angle::cos_sin_bounded`] is where it starts), and for a rational read at
/// enough bits it means [`Bound::ZERO`]. A radius invented for convenience makes every sign above
/// it unearned.
#[derive(Clone, Debug)]
pub(crate) struct HpIv {
    pub mid: BigFloat,
    pub rad: Bound,
}

impl HpIv {
    pub fn new(mid: BigFloat, rad: Bound) -> Self {
        HpIv { mid, rad }
    }

    /// A value known exactly at this precision — a rational whose realization did not round.
    pub fn exact(mid: BigFloat) -> Self {
        HpIv {
            mid,
            rad: Bound::ZERO,
        }
    }

    /// The rounding a `prec`-bit operation adds to its own result: at most a half-ulp,
    /// `|result| · 2⁻ᵖʳᵉᶜ`.
    fn round_off(mid: &BigFloat, prec: usize) -> Bound {
        ub(mid).times(Bound::pow2(-(prec as i64)))
    }

    pub fn sub(&self, o: &HpIv, prec: usize) -> HpIv {
        let mid = self.mid.sub(&o.mid, prec, HP_RM);
        let rad = self.rad.plus(o.rad).plus(Self::round_off(&mid, prec));
        HpIv::new(mid, rad)
    }

    pub fn add(&self, o: &HpIv, prec: usize) -> HpIv {
        let mid = self.mid.add(&o.mid, prec, HP_RM);
        let rad = self.rad.plus(o.rad).plus(Self::round_off(&mid, prec));
        HpIv::new(mid, rad)
    }

    pub fn mul(&self, o: &HpIv, prec: usize) -> HpIv {
        let mid = self.mid.mul(&o.mid, prec, HP_RM);
        // `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, then the rounding of the product itself.
        let rad = ub(&self.mid)
            .times(o.rad)
            .plus(ub(&o.mid).times(self.rad))
            .plus(self.rad.times(o.rad))
            .plus(Self::round_off(&mid, prec));
        HpIv::new(mid, rad)
    }

    /// `Some(true)` if definitely positive, `Some(false)` if definitely negative, `None` if the
    /// interval straddles 0.
    ///
    /// `None` is **"not decided at this precision"**, not "proved zero" — no finite precision
    /// proves a transcendental equality. The caller either climbs the ladder or says so.
    pub fn sign(&self) -> Option<bool> {
        let low = lb(&self.mid)?;
        self.rad.lt(low).then(|| self.mid.is_positive())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_scalar::{Angle, Bound};

    fn big(x: f64, prec: usize) -> BigFloat {
        BigFloat::from_f64(x, prec)
    }

    /// **Does *this* platform's `f64` trig fit inside what [`DA_F64`] charges for it?**
    ///
    /// That constant is the one place the kernel's soundness rests on a measurement rather than a
    /// theorem: every other term here charges round-to-nearest, which IEEE 754 contracts, while
    /// neither Rust nor any libm promises an accuracy for `f64::cos`. The measurement used to live
    /// only in prose beside the constant — so it was an assumption on every machine but the one it
    /// was taken on, and this kernel ships to browsers, whose libm is not this one.
    ///
    /// **Not `#[ignore]`d, and that is the point.** The slow astro-float suites are skipped by
    /// default because they are slow; this one exists precisely so that *every build checks its own
    /// platform*, and a few hundred angles is milliseconds.
    ///
    /// ## The budget is `DA_F64 − ε`, not `DA_F64`
    ///
    /// `Pt3::rotate_about` charges `rot = (|u| + |v|)·DA_F64`, and for an **origin pivot** that term
    /// is alone — `piv` is exactly zero there, so `rot` also has to cover the arithmetic of
    /// `fl(u·c) − fl(v·s)`:
    ///
    /// | | |
    /// |---|---|
    /// | trig realization | `\|u\|·dc + \|v\|·ds` |
    /// | `fl(u·c)`, `fl(v·s)` | `≤ (ε/2)(\|u\| + \|v\|)` together |
    /// | their difference | `≤ (ε/2)(\|u\| + \|v\|)` |
    ///
    /// Requiring the sum to stay under `(|u| + |v|)·DA_F64` and sending `|v|` (or `|u|`) to zero
    /// gives `dc + ε ≤ DA_F64`. **Asserting `dc ≤ DA_F64` instead would pass a platform at `15.9ε`,
    /// which is unsound** — a guard that only ever passes is the worst kind.
    ///
    /// The 90°-family is excluded and that is not a gap: there `cos`/`sin` are exactly `0`/`±1`, so
    /// `dc = ds = 0` *and* every product and their difference are exact — which is why
    /// `rotate_about` sets `rot = 0` for them. `cos_sin_f64_is_exact_for_quadrantal` pins that.
    ///
    /// ## What this does *not* do
    ///
    /// It checks a **premise**. `frame3::tol_bounds_error_over_random_chains` checks the
    /// **conclusion** — that `tol` really bounds the error, second-order terms and all — and both
    /// are needed: that one is only sensitive to a libm ~4.5× worse (its rotation-only chains sit at
    /// 0.2231 of their bound) while soundness is lost at ~1.9×, so there is a window where the
    /// kernel is wrong and the end-to-end test is still green. This closes that window. It also
    /// only protects *builds* — a released wasm binary on a bad platform is not checked by anything
    /// here, which is what charging a *measured* per-angle error instead of a constant would fix.
    #[test]
    fn this_platforms_trig_stays_inside_the_charged_budget() {
        const P: usize = 128; // rad ≈ 2⁻¹²⁸ against an 8ε ≈ 2⁻⁴⁹ quantity: not a participant
        let budget = Bound::of(DA_F64)
            .minus(Bound::of(f64::EPSILON))
            .expect("DA_F64 exceeds one epsilon");
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
            // dominates `dc`, so an angle whose conversion takes the other route has to be here.
            angles.push((rng(1, 1 << 62) * 360, rng(1 << 53, 1 << 62)));
        }
        // ★ **Reported to the nearest ε, which needs a comparison rather than `bf_mag`.**
        // `bf_mag` is `2^exponent` — an upper *octave*, up to 2× above the true magnitude — and
        // `Bound` keeps that shape, so a ratio printed from either would blur exactly the range
        // that matters here. The smallest integer `k` with `|diff| < k·ε` is one comparison per
        // candidate and says the thing plainly.
        let eps_multiple = |diff: &BigFloat| -> usize {
            (1..=64)
                .find(|k| {
                    // astro-float's `cmp` yields a sign as `Option<i128>`, not an `Ordering`.
                    diff.abs()
                        .cmp(&big(*k as f64 * f64::EPSILON, P))
                        .is_some_and(|sign| sign < 0)
                })
                .unwrap_or(usize::MAX)
        };
        let (mut worst_k, mut worst_at) = (0usize, (0i128, 1i128));
        let mut checked = 0usize;
        for (n, d) in angles {
            let Some(a) = Rat::new(n, d).and_then(Angle::from_deg) else {
                continue; // outside Rat's range — `from_deg` declines, and so does the kernel
            };
            let (c, s, exact) = a.cos_sin_f64();
            if exact {
                continue; // 90°-family: zero realization error, and the products are exact too
            }
            let (hc, hs, rc, rs) = a.cos_sin_bounded(P);
            for (f, h, rad) in [(c, &hc, rc), (s, &hs, rs)] {
                // |f64 − true| ≤ |f64 − hp midpoint| + hp's own radius. The assertion runs on
                // `Bound`, which never under-states; the reported `k` is the readable version.
                let diff = big(f, P).sub(h, P, HP_RM);
                let err = Bound::of(bf_mag(&diff).abs()).plus(rad);
                assert!(
                    err.lt(budget),
                    "this platform's trig overruns the budget DA_F64 charges: \
                     {n}/{d} deg, err {err:?} ≥ {budget:?}"
                );
                let k = eps_multiple(&diff);
                if k > worst_k {
                    (worst_k, worst_at) = (k, (n, d));
                }
                checked += 1;
            }
        }
        assert!(checked > 700, "corpus shrank to {checked} realizations");
        // ★ The margin lives here, not in `DA_F64`'s doc — a figure written beside the constant
        // went stale two commits after it was measured and stayed wrong.
        let budget_k = (DA_F64 - f64::EPSILON) / f64::EPSILON;
        eprintln!(
            "[trig budget] worst realization error < {worst_k}ε (at {}/{} deg) against a \
             {budget_k:.0}ε budget — {:.2}x margin, over {checked} realizations",
            worst_at.0,
            worst_at.1,
            budget_k / worst_k as f64,
        );
    }

    /// The point of the type: a value that is only rounding residue must not report a sign, no
    /// matter how confidently its `mid` is nonzero.
    #[test]
    fn a_residue_left_by_cancellation_reports_no_sign() {
        let prec = 200;
        // Two large values differing by far less than the radius they carry. The separation has
        // to stay inside `prec` bits of the operands, or the subtraction is exactly zero and the
        // test proves nothing.
        let rad = Bound::of(1.0e-10);
        let a = HpIv::new(big(1.0e30, prec), rad);
        let b = HpIv::new(big(1.0e30, prec).sub(&big(1.0e-20, prec), prec, HP_RM), rad);
        let d = a.sub(&b, prec);
        assert!(
            !d.mid.is_zero(),
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
        let a = HpIv::new(big(3.0, prec), Bound::pow2(-100));
        let b = HpIv::new(big(2.0, prec), Bound::pow2(-100));
        assert_eq!(a.sub(&b, prec).sign(), Some(true));
        assert_eq!(b.sub(&a, prec).sign(), Some(false));
        assert_eq!(a.mul(&b, prec).sign(), Some(true));
    }

    /// The radius survives a rung deep enough to flush an `f64` radius to zero — the reason
    /// [`Bound`] exists. With an `f64` radius this product would report a confident sign.
    #[test]
    fn a_deep_rung_does_not_lose_the_radius() {
        let prec = 2048;
        let tiny = HpIv::new(big(1.0, prec), Bound::pow2(-(prec as i64)));
        let p = tiny.mul(&tiny, prec);
        assert!(!p.rad.is_zero(), "the radius vanished at {prec} bits");
        // A difference of exactly that size is therefore undecided, not positive.
        let q = HpIv::new(
            big(1.0, prec).add(&BigFloat::from_f64(1.0, prec), prec, HP_RM),
            Bound::pow2(-(prec as i64) + 4),
        );
        let r = HpIv::new(big(2.0, prec), Bound::pow2(-(prec as i64) + 4));
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
        let a = HpIv::exact(big(0.5, prec));
        let b = HpIv::exact(big(0.25, prec));
        assert!(
            a.rad.is_zero() && b.rad.is_zero(),
            "an exact input carried a radius"
        );
        let d = a.sub(&b, prec);
        assert!(
            d.rad.lt(Bound::pow2(-190)),
            "an exact subtraction picked up more than a half-ulp: 2^{:?}",
            d.rad.exp2()
        );
        assert_eq!(d.sign(), Some(true));
        // Exactly equal inputs are the case that cannot be settled here: the difference is zero
        // and zero has no sign to read.
        assert_eq!(a.sub(&a, prec).sign(), None);
    }
}
