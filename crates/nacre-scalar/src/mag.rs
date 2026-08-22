//! [`Mag`] — a non-negative upper bound that does not run out of exponent range.
//!
//! **Why not `f64`.** Every error radius in the judgment path is proportional to `2⁻ᵖʳᵉᶜ`, and the
//! precision ladder climbs. At 2048 bits that factor is `≈1e-617`, which underflows an `f64` to
//! `0` — and a radius of zero says *"this value is exact"*, so the interval declares a sign it has
//! not earned. That is the same failure as a collapsed declare-0 floor, in new clothes. The
//! mantissa needed for a radius is tiny (a handful of bits buys a bound within a percent), but the
//! **exponent** must survive, so this splits them: `m · 2^e` with `e: i64`.
//!
//! **Every operation rounds away from zero.** A bound that rounds to nearest can land just below
//! the quantity it bounds, and one such step is enough to turn "cannot decide" into a wrong
//! answer. Each result is therefore inflated by a relative `2·ε`, which covers the `f64` rounding
//! of the bound arithmetic itself.

/// A non-negative upper bound, `m · 2^e` with `m ∈ [0.5, 1)` when nonzero.
///
/// Comparison and arithmetic are **conservative**: the stored value is never smaller than the
/// quantity it bounds. See the module docs for why an `f64` cannot play this role.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mag {
    m: f64,
    e: i64,
}

/// `x = m · 2^e` with `m ∈ [0.5, 1)`, for finite `x > 0`.
fn frexp(x: f64) -> (f64, i64) {
    debug_assert!(x > 0.0 && x.is_finite());
    let bits = x.to_bits();
    let raw = ((bits >> 52) & 0x7ff) as i64;
    if raw == 0 {
        // Subnormal: scale into the normal range and correct the exponent.
        let (m, e) = frexp(x * 2f64.powi(64));
        (m, e - 64)
    } else {
        // Replace the biased exponent with the one that puts the mantissa in [0.5, 1).
        let m = f64::from_bits((bits & !(0x7ffu64 << 52)) | (1022u64 << 52));
        (m, raw - 1022)
    }
}

/// One step of relative inflation — enough to absorb the `f64` rounding of a single operation on
/// the mantissa, so a rounded-to-nearest result still bounds the true one.
const UP: f64 = 1.0 + 2.0 * f64::EPSILON;

/// Beyond this exponent gap the smaller addend cannot change the larger mantissa, so the sum is
/// the larger one inflated. Keeps `add` free of `2^-de` underflow.
const GAP: i64 = 60;

/// Exactly `2^k`, built from the exponent field.
///
/// `2f64.powi(k)` is what this replaces, and with a **runtime** `k` that is not a constant fold —
/// it is `__powidf2`, a square-and-multiply loop, which a profile of a rotated boolean found at
/// **2.1% of the whole run** (`Mag`'s arithmetic totalled 7.6%, for a type that is two f64s of
/// bookkeeping). A power of two is one exponent field, so this is a shift and a load.
///
/// Callers stay inside `|k| <= GAP`, far from the subnormal edge, so the result is always normal
/// and always exact — the value is identical to `powi`'s, only cheaper.
#[inline]
fn pow2_f64(k: i64) -> f64 {
    debug_assert!(k.abs() <= GAP, "outside the range that keeps this normal");
    f64::from_bits(((1023 + k) as u64) << 52)
}

impl Mag {
    /// The bound `0` — used only where a quantity is *exactly* representable (a rational read at
    /// enough bits), never as a default.
    pub const ZERO: Mag = Mag { m: 0.0, e: 0 };

    /// Exactly `2^e`.
    pub fn pow2(e: i64) -> Mag {
        Mag { m: 0.5, e: e + 1 }
    }

    /// The smallest representable bound at or above `|x| · 2^e`.
    pub fn scaled(x: f64, e: i64) -> Mag {
        let x = x.abs();
        if x == 0.0 || !x.is_finite() {
            return Mag { m: 0.0, e: 0 };
        }
        let (m, me) = frexp(x);
        Mag { m, e: me + e }.inflate()
    }

    /// The smallest representable bound at or above `|x|`.
    pub fn of(x: f64) -> Mag {
        Mag::scaled(x, 0)
    }

    /// Is this bound zero? A zero bound asserts exactness, so callers that could be wrong about
    /// that must not produce one.
    pub fn is_zero(self) -> bool {
        self.m == 0.0
    }

    /// Round the mantissa away from zero by one relative step.
    fn inflate(self) -> Mag {
        if self.m == 0.0 {
            return self;
        }
        let m = self.m * UP;
        // The inflation can carry the mantissa to 1.0; renormalize rather than leave it out of range.
        if m >= 1.0 {
            Mag {
                m: m * 0.5,
                e: self.e + 1,
            }
        } else {
            Mag { m, e: self.e }
        }
    }

    /// An upper bound on `self + other`. Named apart from `std::ops::Add` on purpose: this
    /// rounds **away from zero**, which an operator would not lead a reader to expect.
    pub fn plus(self, other: Mag) -> Mag {
        if self.m == 0.0 {
            return other;
        }
        if other.m == 0.0 {
            return self;
        }
        let (big, small) = if self.e >= other.e {
            (self, other)
        } else {
            (other, self)
        };
        if big.e - small.e > GAP {
            return big.inflate();
        }
        let m = big.m + small.m * pow2_f64(small.e - big.e);
        let (nm, ne) = frexp(m);
        Mag {
            m: nm,
            e: ne + big.e,
        }
        .inflate()
    }

    /// An upper bound on `self · other`, rounded away from zero. See [`plus`](Self::plus).
    pub fn times(self, other: Mag) -> Mag {
        if self.m == 0.0 || other.m == 0.0 {
            return Mag::ZERO;
        }
        let (m, e) = frexp(self.m * other.m);
        Mag {
            m,
            e: e + self.e + other.e,
        }
        .inflate()
    }

    /// An upper bound on `self / other`, rounded away from zero. `None` when `other` is zero —
    /// there is no bound to give, and returning a large one would be a guess.
    pub fn over(self, other: Mag) -> Option<Mag> {
        if other.m == 0.0 {
            return None;
        }
        if self.m == 0.0 {
            return Some(Mag::ZERO);
        }
        let (m, e) = frexp(self.m / other.m);
        Some(
            Mag {
                m,
                e: e + self.e - other.e,
            }
            .inflate(),
        )
    }

    /// A **lower** bound on `self − other`, or `None` when the difference may not be positive.
    ///
    /// This is the one operation here that rounds *toward* zero, and it has to: its result is
    /// used where a quantity must be shown to be at least so large (a norm that has to stay away
    /// from zero before anything is divided by it). Rounding it up the way everything else
    /// rounds would turn a lower bound into a claim, which is the failure this whole type exists
    /// to prevent.
    pub fn minus(self, other: Mag) -> Option<Mag> {
        if other.m == 0.0 {
            return (self.m != 0.0).then_some(self);
        }
        if !other.lt(self) {
            return None; // the difference may be zero or negative — no positive lower bound
        }
        let de = self.e - other.e;
        if de > GAP {
            return Some(self.deflate()); // `other` cannot reach `self`'s mantissa
        }
        let m = self.m - other.m * pow2_f64(-de);
        if m <= 0.0 {
            return None;
        }
        let (nm, ne) = frexp(m);
        Some(
            Mag {
                m: nm,
                e: ne + self.e,
            }
            .deflate(),
        )
    }

    /// Round the mantissa toward zero by one relative step — the mirror of [`Mag::inflate`], for the
    /// lower-bound direction.
    fn deflate(self) -> Mag {
        if self.m == 0.0 {
            return self;
        }
        let m = self.m / UP;
        if m < 0.5 {
            Mag {
                m: m * 2.0,
                e: self.e - 1,
            }
        } else {
            Mag { m, e: self.e }
        }
    }

    /// Is this bound strictly below `other`? Used to ask whether a value's magnitude clears its
    /// own error radius.
    pub fn lt(self, other: Mag) -> bool {
        if self.m == 0.0 {
            return other.m != 0.0;
        }
        if other.m == 0.0 {
            return false;
        }
        (self.e, self.m) < (other.e, other.m)
    }

    /// The base-2 exponent, for reporting. `None` when the bound is zero.
    pub fn exp2(self) -> Option<i64> {
        (self.m != 0.0).then_some(self.e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The reason this type exists: a radius the ladder actually produces must not become zero.
    /// An `f64` cannot hold `2^-2048`; this must, and must still compare as positive.
    #[test]
    fn a_deep_ladder_radius_does_not_underflow_to_zero() {
        let r = Mag::pow2(-2048);
        assert!(!r.is_zero(), "a 2048-bit rung's radius vanished");
        assert!(
            2f64.powi(-2048) == 0.0,
            "the f64 this replaces does underflow"
        );
        // …and it still orders correctly against an even smaller one.
        assert!(Mag::pow2(-4096).lt(r));
        assert!(!r.lt(Mag::pow2(-4096)));
    }

    /// Products and sums stay bounds — never below the exact answer — across the whole exponent
    /// range the ladder uses, where the `f64` arithmetic they are made of has long since flushed
    /// to zero.
    #[test]
    fn arithmetic_stays_above_the_exact_value() {
        // (3/4)·2^-1000 × (5/8)·2^-1000 = (15/32)·2^-2000
        let a = Mag::scaled(0.75, -1000);
        let b = Mag::scaled(0.625, -1000);
        let p = a.times(b);
        // The exact product, one relative step below, must be strictly under the bound.
        let exact_lo = Mag::scaled(0.46875 * (1.0 - 4.0 * f64::EPSILON), -2000);
        assert!(exact_lo.lt(p), "product bound fell below the exact value");
        // Adding a far smaller term never shrinks the sum.
        let s = p.plus(Mag::pow2(-9000));
        assert!(!s.lt(p));
    }

    /// Addition across a small exponent gap is a real sum, not a shortcut to the larger term.
    #[test]
    fn addition_across_a_small_gap_sums_both_terms() {
        let s = Mag::pow2(-100).plus(Mag::pow2(-101)); // = 0.75 · 2^-99
        assert!(Mag::scaled(0.7, -99).lt(s), "the smaller term was dropped");
        assert!(s.lt(Mag::scaled(0.8, -99)), "the sum overshot");
    }

    /// `pow2_f64` is a *replacement*, so what matters is that it replaces exactly — every
    /// exponent a caller can reach must give the bit-identical `powi` value.
    #[test]
    fn the_power_of_two_shortcut_is_bit_identical() {
        for k in -GAP..=GAP {
            assert_eq!(
                pow2_f64(k).to_bits(),
                2f64.powi(k as i32).to_bits(),
                "2^{k}"
            );
        }
    }

    /// `Mag::of` round-trips ordinary magnitudes without collapsing them.
    #[test]
    fn ordinary_magnitudes_survive_the_split_representation() {
        // Including the smallest subnormal, where `f64` has no room left below.
        for x in [1.0, 0.5, 1e-300, 5e-324, 1e300] {
            let b = Mag::of(x);
            assert!(!b.is_zero(), "{x} became a zero bound");
            // Halve via the exponent, so the comparison stays meaningful where `x · 0.9` would
            // round back onto `x`.
            assert!(Mag::scaled(x, -1).lt(b), "{x} lost its ordering");
        }
        assert!(Mag::of(0.0).is_zero());
    }
}
