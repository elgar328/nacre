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
use nacre_scalar::Bound;

use super::HP_RM;

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
    use nacre_scalar::Bound;

    fn big(x: f64, prec: usize) -> BigFloat {
        BigFloat::from_f64(x, prec)
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
