//! Exact rational scalars — the overhaul's §4 value engine.
//!
//! The truth layer for user-input dimensions and angles: a value the user typed
//! is preserved *exactly*, so `1.1` stays `11/10` and `1.1 × 7` is exactly `7.7`
//! (the "thin film" problem's root fix). This is the exact-*value* counterpart to
//! `nacre-predicates`, which decides exact *signs* of geometric determinants —
//! complementary, not redundant.
//!
//! - [`Rat`] — a rational scalar (tol 0). Fixed-width `Ratio<i128>` with
//!   **checked** arithmetic: overflow is a *signal* (the §4 downgrade trigger),
//!   not a panic or silent wrap. The caller downgrades that value's cache to
//!   f64/double-double and records the tol on its `Origin`; the definition is
//!   never lost.
//! - [`Angle`] — rational degrees, normalized mod-360, with exact accumulation so
//!   a full turn lands back on exactly `0` (no f64 drift — a sketch closes). The
//!   `cos`/`sin` realization crosses into f64 (the irrational boundary), but the
//!   90°-family (`0/90/180/270°`) realizes to exact rationals `{0, ±1}` (Niven),
//!   so those rotations of a rational point stay tol 0.
//!
//! Scope (overhaul port #1): this is the exact **value** engine only. The unified
//! `Scalar { value, tol }` wrapper (§4), the arbitrary-precision trig judgment
//! path (astro-float, H1.5), the 2D-frame `orient2d` pipeline, and the kernel
//! wiring that makes points carry these are later cells. Ported from the verified
//! 2D experiment (`experiments/exact2d`, overhaul checkpoint 1).

use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedMul, CheckedSub};
use std::f64::consts::PI;

/// A rational scalar (exact, tol 0). Arithmetic returns `None` on i128 overflow
/// so the caller sees the §4 downgrade trigger explicitly; on overflow the kernel
/// switches that value's cache to f64/double-double and tags its `Origin` with
/// the resulting tol, while the definition (the op-log of input rationals) is
/// preserved. Overflow is far from normal use — adversarial coprime-denominator
/// accumulation reaches it near the i128 ceiling (~122 bits).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Rat(Ratio<i128>);

impl Rat {
    /// Construct `num/den` in lowest terms. `None` if `den == 0`.
    pub fn new(num: i128, den: i128) -> Option<Self> {
        (den != 0).then(|| Rat(Ratio::new(num, den)))
    }

    /// The integer `n` as `n/1`.
    pub fn from_int(n: i128) -> Self {
        Rat(Ratio::from_integer(n))
    }

    /// Exact addition; `None` on overflow (§4 downgrade trigger).
    pub fn checked_add(self, rhs: Self) -> Option<Self> {
        self.0.checked_add(&rhs.0).map(Rat)
    }

    /// Exact subtraction; `None` on overflow.
    pub fn checked_sub(self, rhs: Self) -> Option<Self> {
        self.0.checked_sub(&rhs.0).map(Rat)
    }

    /// Exact multiplication; `None` on overflow.
    pub fn checked_mul(self, rhs: Self) -> Option<Self> {
        self.0.checked_mul(&rhs.0).map(Rat)
    }

    /// Best-f64 image (for the §4 downgrade path and for measurement/export).
    pub fn to_f64(self) -> f64 {
        *self.0.numer() as f64 / *self.0.denom() as f64
    }

    /// Reduced numerator (denominator is always positive after reduction).
    pub fn numer(self) -> i128 {
        *self.0.numer()
    }

    /// Reduced denominator (`> 0`).
    pub fn denom(self) -> i128 {
        *self.0.denom()
    }

    /// Bit-width of the larger of |numer|, |denom| — the "size" §4 watches for bit
    /// growth under chained rational arithmetic (the downgrade threshold).
    pub fn bit_width(self) -> u32 {
        let n = self.numer().unsigned_abs();
        let d = self.denom().unsigned_abs();
        (128 - n.max(d).leading_zeros()).max(1)
    }
}

/// A *direction* angle in degrees, kept normalized to `[0, 360)` exactly
/// (rational). Accumulation is exact: turning by a rational angle repeatedly and
/// completing a full turn lands back on exactly `0` — no f64 drift. `cos`/`sin`
/// realization crosses into f64 (deg→rad via π): the irrational-realization
/// boundary (§4), where the kernel uses double-double / arbitrary precision for
/// *judgment* (a later cell). The angle stays exact; only its realized coordinate
/// carries tol.
///
/// This is the direction type. A multi-turn *amount* (helix pitch × turns, revolve
/// sweep > 360°) must preserve the turn count, so it belongs to a separate
/// *unnormalized* `Sweep` type — flat 2D sketches never need it, so it is left as
/// a documented companion, not built here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Angle(Rat); // invariant: 0 <= inner < 360

impl Angle {
    /// Normalize `deg` into `[0, 360)` (exact). `None` on overflow during
    /// reduction (the §4 downgrade trigger). Assumes `deg` is within a few turns
    /// of the range — a multi-turn amount is the future `Sweep` type's job, so
    /// this reduces by subtracting whole turns rather than a `360·denom` divide
    /// (which would overflow at large denominators).
    pub fn from_deg(deg: Rat) -> Option<Self> {
        let full = Rat::from_int(360);
        let zero = Rat::from_int(0);
        let mut v = deg;
        while v >= full {
            v = v.checked_sub(full)?;
        }
        while v < zero {
            v = v.checked_add(full)?;
        }
        Some(Angle(v))
    }

    /// Turn by `delta_deg` (exact) and renormalize into `[0, 360)`.
    pub fn checked_add(self, delta_deg: Rat) -> Option<Self> {
        Self::from_deg(self.0.checked_add(delta_deg)?)
    }

    /// The normalized degree value (exact, in `[0, 360)`).
    pub fn deg(self) -> Rat {
        self.0
    }

    /// cos, realized in f64 (deg→rad via π — the irrational boundary).
    pub fn cos(self) -> f64 {
        (self.0.to_f64() * PI / 180.0).cos()
    }

    /// sin, realized in f64.
    pub fn sin(self) -> f64 {
        (self.0.to_f64() * PI / 180.0).sin()
    }

    /// Exact `(cos, sin)` as rationals — `Some` only for the quadrantal angles
    /// (0/90/180/270°), the sole angles where *both* are rational (`{0, ±1}`, by
    /// Niven); `None` otherwise, so the caller falls to the f64/dd realization.
    ///
    /// The values `0`/`±1` are exact in f64 too, so this returns `Rat` not for
    /// representability but to keep a rotated *point* exact: `x·cos − y·sin` must
    /// stay in rational arithmetic, and multiplying a rational coordinate by an
    /// f64 (even an exact `0.0`/`1.0`) would drop the point into f64 and lose the
    /// very exactness this path exists for. So a 90°-family rotation of a rational
    /// point stays tol 0.
    pub fn try_exact_cos_sin(self) -> Option<(Rat, Rat)> {
        let zero = Rat::from_int(0);
        let one = Rat::from_int(1);
        let neg_one = Rat::from_int(-1);
        if self.0 == Rat::from_int(0) {
            Some((one, zero))
        } else if self.0 == Rat::from_int(90) {
            Some((zero, one))
        } else if self.0 == Rat::from_int(180) {
            Some((neg_one, zero))
        } else if self.0 == Rat::from_int(270) {
            Some((zero, neg_one))
        } else {
            None
        }
    }
}
