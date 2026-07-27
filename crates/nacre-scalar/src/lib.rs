//! Exact rational scalars — the overhaul's exact rational value engine (`Rat`/`Angle`) and
//! the axis/rotation/isometry value types the kernel builds on. (The toleranced-sign frame
//! judgment this value engine enables now lives in `nacre-cip`.)
//!
//! The truth layer for user-input dimensions and angles: a value the user typed
//! is preserved *exactly*, so `1.1` stays `11/10` and `1.1 × 7` is exactly `7.7`
//! (the "thin film" problem's root fix). This is the exact-*value* counterpart to
//! `nacre-predicates`, which decides exact *signs* of geometric determinants —
//! complementary, not redundant.
//!
//! - [`Rat`] — a rational scalar (tol 0). Fixed-width `Ratio<i128>` with
//!   **checked** arithmetic: overflow is a *signal* (the downgrade trigger),
//!   not a panic or silent wrap. The caller downgrades that value's cache to
//!   f64/double-double and records the tol on its `Origin`; the definition is
//!   never lost.
//! - [`Angle`] — rational degrees, normalized mod-360, with exact accumulation so
//!   a full turn lands back on exactly `0` (no f64 drift — a sketch closes). The
//!   `cos`/`sin` realization crosses into f64 (the irrational boundary); the
//!   90°-family (`0/90/180/270°`) realizes to exact rationals `{0, ±1}` (Niven),
//!   so those rotations of a rational point stay tol 0; and `cos_hp`/`sin_hp`
//!   realize in arbitrary precision (astro-float) for the judgment path.
//!
//! Scope: the exact value engine (the toleranced-sign frame judgment it enabled now lives
//! in `nacre-cip`). Still deferred to later cells: the unified `Scalar { value, tol }`
//! wrapper, the declare-0 → user-confirmation policy, and the kernel wiring that makes
//! geometry carry these. Ported from the verified 2D experiment (`experiments/exact2d`).

pub mod bound;
pub use bound::Bound;

use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedMul, CheckedSub};
use std::f64::consts::PI;

use astro_float::{BigFloat, Consts, RoundingMode};
use std::cell::RefCell;

/// Precision (bits) and rounding for the high-precision realization layer
/// (astro-float). ~160 bits ≈ 48 decimal digits — far below CAD tolerance and
/// dial-able if amplification ever needs more (the ceiling twofloat lacked, H1.5).
pub(crate) const HP_PREC: usize = 160;
pub(crate) const HP_RM: RoundingMode = RoundingMode::ToEven;

thread_local! {
    /// Transcendental-constant cache (π, …) for the high-precision layer.
    static HP_CONSTS: RefCell<Consts> = RefCell::new(Consts::new().expect("astro-float consts"));
}

/// The result of an orientation judgment (§CIP) — the shared sign vocabulary used by both
/// the 2D and 3D toleranced-sign judges (now in `nacre-cip`) and their downstream consumers.
/// A cross-cutting "judgment result" carried here as a fundamental value (a candidate to
/// split into its own vocabulary type later).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orient {
    Positive,
    Negative,
    /// Declared 0 — collinear within the escalation precision cap. This is the
    /// `declare-0` case where the kernel (later) stops and asks the user.
    Zero,
}

/// A rational scalar (exact, tol 0). Arithmetic returns `None` on i128 overflow
/// so the caller sees the downgrade trigger explicitly; on overflow the kernel
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

    /// Exact addition; `None` on overflow (downgrade trigger).
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

    /// Best-f64 image (for the downgrade path and for measurement/export).
    pub fn to_f64(self) -> f64 {
        *self.0.numer() as f64 / *self.0.denom() as f64
    }

    /// The **exact** rational value of an f64 (`mantissa · 2^exp`). `None` for a
    /// non-finite input or when the exact numerator/denominator overflows i128
    /// (subnormals, extreme exponents — the downgrade trigger). Round-trips:
    /// `to_f64(try_from_f64(x).unwrap()) == x` for finite f64 in the normal CAD range.
    /// The bridge from an f64 coordinate cache to an exact `[Rat]` definition.
    pub fn try_from_f64(x: f64) -> Option<Self> {
        if x == 0.0 {
            return Some(Rat::from_int(0));
        }
        if !x.is_finite() {
            return None;
        }
        let bits = x.to_bits();
        let neg = bits >> 63 == 1;
        let exp_field = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & 0x000f_ffff_ffff_ffff;
        // value = mantissa · 2^exp (implicit leading 1 for a normal; bias 1023, and the
        // 52-bit fraction shifts the exponent by another 52).
        let (mantissa, exp) = if exp_field == 0 {
            (frac, -1074) // subnormal
        } else {
            (frac | 0x0010_0000_0000_0000, exp_field - 1075)
        };
        let m = mantissa as i128;
        let (numer, denom) = if exp >= 0 {
            if exp > 126 {
                return None;
            }
            (m.checked_mul(1i128 << exp)?, 1i128)
        } else {
            let k = (-exp) as u32;
            if k > 126 {
                return None;
            }
            (m, 1i128 << k)
        };
        Rat::new(if neg { -numer } else { numer }, denom)
    }

    /// Reduced numerator (denominator is always positive after reduction).
    pub fn numer(self) -> i128 {
        *self.0.numer()
    }

    /// Reduced denominator (`> 0`).
    pub fn denom(self) -> i128 {
        *self.0.denom()
    }

    /// Bit-width of the larger of |numer|, |denom| — the "size" the kernel watches for bit
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
/// boundary, where the kernel uses double-double / arbitrary precision for
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
    /// reduction (the downgrade trigger). Assumes `deg` is within a few turns
    /// of the range — a multi-turn amount is the future `Sweep` type's job, so
    /// this reduces by subtracting whole turns rather than a `360·denom` divide
    /// (which would overflow at large denominators).
    pub fn from_deg(deg: Rat) -> Option<Self> {
        // Exact reduction mod 360, in one division. Subtracting a turn at a time is the same
        // arithmetic but costs one iteration per turn, so an input like `2⁶⁰` degrees — a value a
        // script can produce without meaning anything unusual by it — does not return.
        let (n, d) = (deg.numer(), deg.denom()); // `d > 0` after reduction
        let full = 360i128.checked_mul(d)?; // one turn, in units of `1/d`
        let rem = n.rem_euclid(full); // `0 ≤ rem < full`, so the result is in `[0, 360)`
        Rat::new(rem, d).map(Angle)
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

    /// `(cos, sin)` realized in arbitrary precision at `prec` bits — the judgment
    /// path. astro-float replaces twofloat here (H1.5: twofloat's trig was
    /// f64-level near zero-crossings). Higher `prec` gives ground truth; the
    /// default [`HP_PREC`] gives the working judgment realization. (numer/denom
    /// pass through f64, exact for the small values used here; a general large-
    /// rational path would build from a string.)
    pub fn cos_sin_at(self, prec: usize) -> (BigFloat, BigFloat) {
        let (c, s, _, _) = self.cos_sin_bounded(prec);
        (c, s)
    }

    /// `(cos, sin)` at `prec` bits **with an upper bound on how far each may be from the true
    /// value** — the seed every error radius in the judgment path grows from.
    ///
    /// The bound is *derived*, term by term, not chosen:
    ///
    /// - `numer`/`denom` enter as `i128`, exactly. Routing them through `f64` would cost a
    ///   relative `2⁻⁵³` past 2⁵³ that no working precision recovers, so it is not done.
    /// - `n/d`, `·π`, `/180` are three round-to-nearest operations at `prec` bits, each a relative
    ///   `2⁻ᵖʳᵉᶜ`, and `π` itself carries one more.
    /// - The argument's absolute error `δθ = |θ|·ρ` passes through the **derivative**:
    ///   `d cos = −sin·dθ` and `d sin = cos·dθ`. Slope 1 would also be sound, but near a zero
    ///   crossing the true slope is what keeps a tiny result from being swamped by its own bound.
    /// - astro-float's `cos`/`sin` run **Ziv's loop** — `cos_series` at a working precision, then
    ///   `try_set_precision(p, rm, p_wrk)`, retrying with more bits when the value sits too close
    ///   to a rounding boundary to decide. That is the standard construction for a *correctly
    ///   rounded* transcendental, so the realization adds at most a half-ulp — `|result| · 2⁻ᵖʳᵉᶜ`,
    ///   *relative* to the value, which is why the two functions get separate bounds. The crate
    ///   does not document this, so `the_trig_bound_holds_against_a_far_deeper_realization` checks
    ///   it rather than trusting it.
    ///
    /// Returns `(cos, sin, |Δcos|, |Δsin|)`.
    pub fn cos_sin_bounded(self, prec: usize) -> (BigFloat, BigFloat, Bound, Bound) {
        HP_CONSTS.with_borrow_mut(|cc| {
            let pi = cc.pi(prec, HP_RM);
            let d180 = BigFloat::from_f64(180.0, prec);
            // As `i128`, not through `f64`: past 2⁵³ the conversion would cost a relative
            // `2⁻⁵³` that no working precision recovers, because the loss happens before
            // astro-float sees the value.
            let ip = prec.max(128);
            let n = BigFloat::from_i128(self.0.numer(), ip);
            let d = BigFloat::from_i128(self.0.denom(), ip);
            let rad = n
                .div(&d, prec, HP_RM)
                .mul(&pi, prec, HP_RM)
                .div(&d180, prec, HP_RM);
            let u = Bound::pow2(-(prec as i64));
            // Relative error of the argument: the two `i128 → f64` conversions, then four
            // rounded high-precision operations (the division, the product, the division, and π).
            // Four rounded operations build the argument (the division, the product, the
            // division, and π itself). The integers contribute nothing — they go in exactly.
            let rel = u.times(Bound::of(4.0));
            // `|θ|` in radians, over-estimated from its exponent (`|x| < 2^exponent`).
            let theta = match rad.exponent() {
                Some(e) if !rad.is_zero() => Bound::pow2(e as i64),
                _ => Bound::ZERO,
            };
            let d_theta = theta.times(rel);
            let (c, s) = (rad.cos(prec, HP_RM, cc), rad.sin(prec, HP_RM, cc));
            // `|x| < 2^exponent` — the slope of the *other* function, and the scale of the
            // half-ulp of this one.
            let ub = |x: &BigFloat| match x.exponent() {
                Some(e) if !x.is_zero() => Bound::pow2(e as i64),
                _ => Bound::ZERO,
            };
            let (uc, us) = (ub(&c), ub(&s));
            let err_cos = us.times(d_theta).plus(uc.times(u));
            let err_sin = uc.times(d_theta).plus(us.times(u));
            (c, s, err_cos, err_sin)
        })
    }

    /// cos realized at the default judgment precision ([`HP_PREC`]).
    pub fn cos_hp(self) -> BigFloat {
        self.cos_sin_at(HP_PREC).0
    }

    /// sin realized at the default judgment precision ([`HP_PREC`]).
    pub fn sin_hp(self) -> BigFloat {
        self.cos_sin_at(HP_PREC).1
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

    /// `(cos, sin)` realized in f64 — **exact** (`0.0`/`±1.0`) for the 90°-family,
    /// plain [`cos`](Self::cos)/[`sin`](Self::sin) otherwise. The single source of
    /// truth for realizing a rotation angle into f64: every path that turns a point or
    /// direction by an angle must go through here, so the quadrantal case never
    /// re-introduces the `cos(90°)≈6e-17` spurious cross-term (an axis-aligned rotation
    /// then lands its coordinates exactly on the grid — a 90°-family rotation is tol 0).
    pub fn cos_sin_f64(self) -> (f64, f64) {
        match self.try_exact_cos_sin() {
            Some((cr, sr)) => (cr.to_f64(), sr.to_f64()),
            None => (self.cos(), self.sin()),
        }
    }
}

/// A coordinate axis — the fixed axis of an axis-aligned rotation (overhaul stage
/// 1b restricts to `X`/`Y`/`Z`, the form `exact3d` validated; arbitrary rational
/// axes via Rodrigues are a later extension).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    /// The two in-plane coordinate indices (the third is the fixed rotation axis).
    /// The order gives a right-handed (CCW-about-the-axis) rotation.
    pub fn plane(self) -> (usize, usize) {
        match self {
            Axis::X => (1, 2), // rotate y,z
            Axis::Y => (2, 0), // rotate z,x
            Axis::Z => (0, 1), // rotate x,y
        }
    }
}

/// An axis-aligned rigid rotation: turn about `axis` (the line through the rational
/// `point`) by the rational `angle`. Exact for the 90°-family (`try_exact_cos_sin`);
/// otherwise the realized coordinate is irrational (cos/sin) and carries tol.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotation {
    pub axis: Axis,
    pub point: [Rat; 3],
    pub angle: Angle,
}

/// A rigid-body isometry (§ Transform): a rotation (optional) then a translation.
/// The exact rational data is the **definition**; the `apply_*`/`offset_f64`
/// realizers give the f64 cache. Math-type independent — operates on plain
/// `[f64; 3]`, mirroring [`frame2::Pt2`]'s `(f64, f64)` (so `nacre-scalar` never
/// depends on `nacre-math`); the caller (`nacre-ops`) applies it to `Point3`/`Plane`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Isometry {
    /// Applied first: an axis-aligned rotation, or `None` (pure translation).
    pub rotate: Option<Rotation>,
    /// Applied second: an exact rational translation.
    pub translate: [Rat; 3],
}

impl Isometry {
    /// A pure translation by the rational vector `translate`.
    pub fn translation(translate: [Rat; 3]) -> Self {
        Isometry {
            rotate: None,
            translate,
        }
    }

    /// A pure rotation (no translation).
    pub fn rotation(rotate: Rotation) -> Self {
        Isometry {
            rotate: Some(rotate),
            translate: [Rat::from_int(0); 3],
        }
    }

    /// A rotation followed by a translation.
    pub fn rigid(rotate: Rotation, translate: [Rat; 3]) -> Self {
        Isometry {
            rotate: Some(rotate),
            translate,
        }
    }

    /// The translation realized in f64.
    pub fn offset_f64(&self) -> [f64; 3] {
        [
            self.translate[0].to_f64(),
            self.translate[1].to_f64(),
            self.translate[2].to_f64(),
        ]
    }

    /// Whether the isometry realizes exactly: no rotation, or a 90°-family rotation
    /// (`try_exact_cos_sin` gives rational cos/sin, so an f64-representable point
    /// stays exact — tol 0). A non-90° rotation realizes to irrational f64 (tol > 0).
    pub fn is_exact(&self) -> bool {
        match self.rotate {
            None => true,
            Some(r) => r.angle.try_exact_cos_sin().is_some(),
        }
    }

    /// Apply the full isometry (rotate about the axis point, then translate) to a
    /// point realized in f64.
    pub fn apply_point(&self, p: [f64; 3]) -> [f64; 3] {
        let mut q = p;
        if let Some(r) = self.rotate {
            let (i, j) = r.axis.plane();
            let (px, py) = (r.point[i].to_f64(), r.point[j].to_f64());
            let (c, s) = r.angle.cos_sin_f64();
            let (dx, dy) = (p[i] - px, p[j] - py);
            q[i] = px + dx * c - dy * s;
            q[j] = py + dx * s + dy * c;
        }
        let off = self.offset_f64();
        [q[0] + off[0], q[1] + off[1], q[2] + off[2]]
    }

    /// Apply only the rotation (no axis point, no translation) to a direction.
    pub fn apply_dir(&self, d: [f64; 3]) -> [f64; 3] {
        let mut q = d;
        if let Some(r) = self.rotate {
            let (i, j) = r.axis.plane();
            let (c, s) = r.angle.cos_sin_f64();
            let (dx, dy) = (d[i], d[j]);
            q[i] = dx * c - dy * s;
            q[j] = dx * s + dy * c;
        }
        q
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// `try_from_f64` is the exact rational of the f64: it round-trips (`to_f64` gives
    /// back the same bits) for finite values in range, handles 0/integers/dyadic
    /// fractions exactly, and returns `None` for non-finite or overflowing inputs.
    #[test]
    fn try_from_f64_is_exact_and_round_trips() {
        assert_eq!(Rat::try_from_f64(0.0), Some(Rat::from_int(0)));
        assert_eq!(Rat::try_from_f64(6.0), Some(Rat::from_int(6)));
        assert_eq!(Rat::try_from_f64(-2.0), Some(Rat::from_int(-2)));
        assert_eq!(Rat::try_from_f64(0.5), Rat::new(1, 2));
        assert_eq!(Rat::try_from_f64(-0.75), Rat::new(-3, 4));
        assert_eq!(Rat::try_from_f64(f64::NAN), None);
        assert_eq!(Rat::try_from_f64(f64::INFINITY), None);
        assert_eq!(Rat::try_from_f64(1e300), None); // exponent overflows i128
        // round-trip over a spread of normal-range values (incl. non-dyadic f64s).
        for &x in &[0.1, 1.0 / 3.0, 2.0, 1000.0, -6.1, 4_503.7, 1e-6, 1e6] {
            assert_eq!(Rat::try_from_f64(x).unwrap().to_f64(), x, "round-trip {x}");
        }
    }

    /// The core rational-representation property in miniature: exact rational accumulation does not
    /// drift, where the f64 control does. `(1/10)` summed ten times is exactly
    /// `1`, but `0.1_f64` summed ten times is not `1.0`.
    #[test]
    fn rational_accumulation_is_exact_where_f64_drifts() {
        let tenth = Rat::new(1, 10).unwrap();
        let mut acc = Rat::from_int(0);
        for _ in 0..10 {
            acc = acc.checked_add(tenth).unwrap();
        }
        assert_eq!(acc, Rat::from_int(1));

        let mut f = 0.0_f64;
        for _ in 0..10 {
            f += 0.1;
        }
        assert_ne!(f, 1.0); // 0.9999999999999999 — the drift Rat avoids
    }

    /// The "thin film" example: `1.1 × 7` must be exactly `7.7`. In rationals
    /// `11/10 × 7 = 77/10`; in f64 `1.1 * 7.0` is not `7.7`.
    #[test]
    fn one_point_one_times_seven_is_exact() {
        let a = Rat::new(11, 10).unwrap();
        let seven = Rat::from_int(7);
        assert_eq!(a.checked_mul(seven).unwrap(), Rat::new(77, 10).unwrap());

        assert_ne!(1.1_f64 * 7.0, 7.7); // f64 cannot represent 7.7 exactly
    }

    /// Measurement: with fixed-width i128, chained coprime-denominator
    /// accumulation *does* overflow (the finite-precision cliff handled by
    /// downgrading). Summing `1/p` over successive primes forces the denominator
    /// toward the primorial, which exceeds i128. This test pins two facts:
    ///   (a) a modest sum (first 8 primes) stays exact — normal use is fine;
    ///   (b) accumulation eventually overflows within the prime list — proving
    ///       the downgrade trigger fires, and reporting *where* (onset index).
    #[test]
    fn coprime_accumulation_overflows_and_reports_onset() {
        const PRIMES: [i128; 30] = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97, 101, 103, 107, 109, 113,
        ];

        // (a) first 8 primes stay exact.
        let mut acc = Rat::from_int(0);
        for &p in &PRIMES[..8] {
            acc = acc
                .checked_add(Rat::new(1, p).unwrap())
                .expect("first 8 primes must not overflow");
        }

        // (b) full accumulation eventually overflows; record the onset.
        let mut acc = Rat::from_int(0);
        let mut onset = None;
        let mut last_bits = 0;
        for (i, &p) in PRIMES.iter().enumerate() {
            match acc.checked_add(Rat::new(1, p).unwrap()) {
                Some(next) => {
                    acc = next;
                    last_bits = acc.bit_width();
                }
                None => {
                    onset = Some((i, last_bits));
                    break;
                }
            }
        }
        let (idx, bits) = onset.expect("i128 rational accumulation must overflow within 30 primes");
        // Measurement record (run with `-- --nocapture`).
        eprintln!(
            "[rational] coprime 1/p accumulation overflows at prime index {idx} (p={}); \
             denominator bit-width just before onset = {bits}",
            PRIMES[idx]
        );
        // Onset is comfortably past normal use and near the i128 ceiling (~127 bits).
        assert!(idx >= 8, "overflow onset {idx} should be past modest use");
        assert!(
            bits > 100,
            "denominator should be near the i128 ceiling at onset, got {bits} bits"
        );
    }

    /// Closure: a rational angle accumulates exactly, so a full turn lands back
    /// on exactly `0` — while the f64 control drifts. `360/7` degrees added seven
    /// times is exactly `360 → 0`; `360.0/7.0` summed seven times in f64 is not
    /// `360.0`.
    #[test]
    fn rational_angle_closes_exactly_where_f64_drifts() {
        let step = Rat::new(360, 7).unwrap();
        let mut a = Angle::from_deg(Rat::from_int(0)).unwrap();
        for _ in 0..7 {
            a = a.checked_add(step).unwrap();
        }
        assert_eq!(a, Angle::from_deg(Rat::from_int(0)).unwrap());

        let mut f = 0.0_f64;
        for _ in 0..7 {
            f += 360.0 / 7.0;
        }
        assert_ne!(f, 360.0); // f64 drifts off the full turn
    }

    /// Many small rational steps also close: `1/3` degree added 1080 times is
    /// exactly one full turn → `0`.
    #[test]
    fn many_small_steps_close_to_zero() {
        let step = Rat::new(1, 3).unwrap();
        let mut a = Angle::from_deg(Rat::from_int(0)).unwrap();
        for _ in 0..1080 {
            a = a.checked_add(step).unwrap();
        }
        assert_eq!(a, Angle::from_deg(Rat::from_int(0)).unwrap());
    }

    /// The angle value stays exact, but its cos/sin *realization* is f64: `cos 90°`
    /// is not exactly `0` (it is ~6e-17), showing the irrational-realization
    /// boundary. `cos 0°`/`sin 0°` happen to be exact in f64.
    #[test]
    fn realization_is_f64_while_angle_stays_exact() {
        let a = Angle::from_deg(Rat::from_int(90)).unwrap();
        assert_eq!(a.deg(), Rat::from_int(90)); // angle exact
        assert!(a.cos().abs() < 1e-15 && a.cos() != 0.0); // realized near 0, not exact

        let z = Angle::from_deg(Rat::from_int(0)).unwrap();
        assert_eq!(z.cos(), 1.0);
        assert_eq!(z.sin(), 0.0);
    }

    /// A 90°-family angle yields exact rational cos/sin, so rotating a rational
    /// point stays exact (tol 0): `(x, y)` rotated 90° is `(-y, x)`. A 45° angle
    /// has no exact rational cos/sin (√2/2), so it returns `None` and would fall
    /// to the f64/dd realization. Note the f64 path is *not* exact here:
    /// `cos 90°` realizes to ~6e-17, not `0`.
    #[test]
    fn exact_cos_sin_only_for_quadrantal_angles() {
        let a90 = Angle::from_deg(Rat::from_int(90)).unwrap();
        let (c, s) = a90.try_exact_cos_sin().unwrap();
        assert_eq!((c, s), (Rat::from_int(0), Rat::from_int(1)));

        // Exact rotation of (3, 5) by 90° → (-5, 3), all rational (tol 0):
        // x' = x·cos − y·sin,  y' = x·sin + y·cos.
        let (x, y) = (Rat::from_int(3), Rat::from_int(5));
        let xr = x
            .checked_mul(c)
            .unwrap()
            .checked_sub(y.checked_mul(s).unwrap())
            .unwrap();
        let yr = x
            .checked_mul(s)
            .unwrap()
            .checked_add(y.checked_mul(c).unwrap())
            .unwrap();
        assert_eq!((xr, yr), (Rat::from_int(-5), Rat::from_int(3)));

        // 45° has no exact rational realization → None (falls to f64/dd).
        let a45 = Angle::from_deg(Rat::from_int(45)).unwrap();
        assert!(a45.try_exact_cos_sin().is_none());
        assert_ne!(a90.cos(), 0.0); // the general f64 path is not exact at 90°
    }

    /// `cos_sin_f64` snaps the 90°-family to exact `0.0`/`±1.0` (unlike `cos()`/`sin()`
    /// which realize ~6e-17 at 90°), and falls through to `cos()`/`sin()` otherwise.
    #[test]
    fn cos_sin_f64_is_exact_for_quadrantal() {
        let deg = |d| Angle::from_deg(Rat::from_int(d)).unwrap();
        assert_eq!(deg(0).cos_sin_f64(), (1.0, 0.0));
        assert_eq!(deg(90).cos_sin_f64(), (0.0, 1.0));
        assert_eq!(deg(180).cos_sin_f64(), (-1.0, 0.0));
        assert_eq!(deg(270).cos_sin_f64(), (0.0, -1.0));
        // non-quadrantal: identical to the plain f64 realization.
        let a45 = deg(45);
        assert_eq!(a45.cos_sin_f64(), (a45.cos(), a45.sin()));
    }

    /// `apply_point`/`apply_dir` realize a 90°-family rotation bit-exactly: no ~6e-17
    /// spurious cross-term. (3,5,z) about Z by 90° → exactly (-5,3,z); a rational-pivot
    /// rotation is exact too; a non-quadrantal angle is unchanged from the f64 path.
    #[test]
    fn apply_point_exact_for_quadrantal() {
        let iso = |d| {
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(d)).unwrap(),
            })
        };
        assert_eq!(iso(90).apply_point([3.0, 5.0, 7.0]), [-5.0, 3.0, 7.0]);
        assert_eq!(iso(180).apply_point([3.0, 5.0, 7.0]), [-3.0, -5.0, 7.0]);
        assert_eq!(iso(270).apply_point([3.0, 5.0, 7.0]), [5.0, -3.0, 7.0]);
        assert_eq!(iso(90).apply_dir([0.0, 1.0, 0.0]), [-1.0, 0.0, 0.0]);

        // Non-origin pivot (2,2): 90° maps (3,5)→pivot+R(1,3)=(2-3, 2+1)=(-1,3).
        let piv = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(2), Rat::from_int(2), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
        });
        assert_eq!(piv.apply_point([3.0, 5.0, 0.0]), [-1.0, 3.0, 0.0]);

        // Non-quadrantal: unchanged from the plain cos/sin realization.
        let a = Angle::from_deg(Rat::from_int(37)).unwrap();
        let (c, s) = (a.cos(), a.sin());
        assert_eq!(
            iso(37).apply_point([3.0, 5.0, 0.0]),
            [3.0 * c - 5.0 * s, 3.0 * s + 5.0 * c, 0.0]
        );
    }

    /// H1.5: the arbitrary-precision cos/sin realization must be far more accurate
    /// than any f64/double-double — the accuracy gate astro-float passes and the
    /// double-double `twofloat` failed (its trig degraded to ~1e-16 near zero-
    /// crossings). Error at rational angles must beat `2^-100` (~1e-30).
    #[test]
    fn high_precision_trig_meets_gate() {
        const GATE_EXP: i32 = -100; // 2^-100 ≈ 7.9e-31
        // (deg, exact_cos, exact_sin | None where irrational, e.g. sin 60°)
        let cases = [
            (0i128, 1.0, Some(0.0)),
            (60, 0.5, None),
            (90, 0.0, Some(1.0)),
            (120, -0.5, None),
            (180, -1.0, Some(0.0)),
            (270, 0.0, Some(-1.0)),
        ];
        let mut worst = i32::MIN;
        for (deg, exact_cos, exact_sin) in cases {
            let a = Angle::from_deg(Rat::from_int(deg)).unwrap();
            let ec = hp_err_exp(&a.cos_hp(), exact_cos);
            assert!(ec < GATE_EXP, "cos {deg}° error 2^{ec} exceeds gate");
            worst = worst.max(ec);
            if let Some(s) = exact_sin {
                let es = hp_err_exp(&a.sin_hp(), s);
                assert!(es < GATE_EXP, "sin {deg}° error 2^{es} exceeds gate");
                worst = worst.max(es);
            }
        }
        eprintln!(
            "[H1.5] astro-float {HP_PREC}-bit worst cos/sin error at rational angles: 2^{worst}"
        );
    }

    /// A large angle normalizes in one step, not one step per turn.
    ///
    /// Normalization used to subtract 360° in a loop, so `2⁶⁰` degrees needed ~3·10¹⁵ iterations —
    /// the kernel did not reject that input, it stopped responding to it. Found by execution: a
    /// trig corpus reached for a big numerator and the test never returned.
    #[test]
    fn a_huge_angle_normalizes_without_counting_turns() {
        let huge = Rat::new(1i128 << 60, 7).unwrap();
        let a = Angle::from_deg(huge).expect("a large rational angle is representable");
        assert!(a.deg() >= Rat::from_int(0) && a.deg() < Rat::from_int(360));
        // Same residue class as the input, so the reduction is `− 360k`, not a different angle.
        let back = a.deg().checked_sub(huge).unwrap();
        let turns = back.to_f64() / -360.0;
        assert_eq!(
            turns.fract(),
            0.0,
            "the reduction was not a whole number of turns"
        );
        // Negatives land in range too, and exactly on 0 at a full turn.
        assert_eq!(
            Angle::from_deg(Rat::from_int(-720)).unwrap().deg(),
            Rat::from_int(0)
        );
        assert_eq!(
            Angle::from_deg(Rat::new(-1, 2).unwrap()).unwrap().deg(),
            Rat::new(719, 2).unwrap()
        );
    }

    /// **The seed of every error radius, checked against a realization far deeper than itself.**
    ///
    /// [`Angle::cos_sin_bounded`] derives its bound from two things the crate does not promise in
    /// writing: that `Consts::pi` is correctly rounded, and that `cos`/`sin` are too (they run
    /// Ziv's loop, which is how one builds a correctly-rounded transcendental — but an
    /// implementation detail, not a documented contract). If either weakens, every interval above
    /// it is unsound, so the claim is measured: at each rung the ladder uses, the value must sit
    /// within its own bound of the same value realized with 512 extra bits.
    ///
    /// The reference is not independent code — it is the same routine at higher precision — so
    /// this cannot catch an error that grows with precision in the same shape. What it does catch
    /// is the failure that matters here: a bound that is simply too small.
    #[test]
    fn the_trig_bound_holds_against_a_far_deeper_realization() {
        // Angles spanning the quadrants, plus rationals with awkward denominators and one whose
        // numerator is large enough to exercise the `i128 → f64` term.
        let angles = [
            (0i128, 1i128),
            (30, 1),
            (45, 1),
            (60, 1),
            (90, 1),
            (135, 1),
            (180, 1),
            (271, 1),
            (359, 1),
            (1, 7),
            (22, 7),
            (1000, 3),
            (1, 1_000_000),
            // A denominator past 2⁵³ that is still exact in `f64` (a power of two), so the bound
            // must *not* charge it the conversion term.
            (1, 1i128 << 60),
            // …and one that genuinely does not round-trip, where the angle itself is only known
            // to a relative `2⁻⁵³` and no working precision can recover it.
            (1, (1i128 << 60) + 1),
        ];
        // Slack is tracked per rung so the two stories stay separable: a word-aligned precision
        // is delivered as asked, while `200` is silently rounded up to 256 and the extra bits show
        // up as slack that is astro-float's, not this bound's.
        let mut worst = std::collections::BTreeMap::<usize, (i64, String)>::new();
        for prec in [128usize, 200, 256, 512, 1024] {
            for (num, den) in angles {
                let a = Angle::from_deg(Rat::new(num, den).unwrap()).unwrap();
                let (c, s, bc, bs) = a.cos_sin_bounded(prec);
                let deep = prec + 512;
                let (rc, rs) = a.cos_sin_at(deep);
                for (got, reference, bound, what) in [(&c, &rc, bc, "cos"), (&s, &rs, bs, "sin")] {
                    let diff = got.sub(reference, deep, HP_RM);
                    let Some(de) = (if diff.is_zero() {
                        None
                    } else {
                        diff.exponent()
                    }) else {
                        continue; // exactly equal — nothing to bound
                    };
                    // `|diff| < 2^de`; the bound must be at least that.
                    let observed = Bound::pow2(de as i64);
                    assert!(
                        !bound.lt(observed),
                        "{what} {num}/{den}° at {prec} bits: error 2^{de} exceeds its bound 2^{:?}",
                        bound.exp2()
                    );
                    // Track how much slack the bound carries, so a bound that is merely
                    // enormous does not pass as a bound that is right.
                    // Slack is only meaningful where the angle converts exactly. Where it does
                    // not, the `2⁻⁵³` term dominates by design and the gap to the observed error
                    // is the honest cost of an unrepresentable angle, not looseness.
                    let angle_exact = (num as f64) as i128 == num && (den as f64) as i128 == den;
                    // And only where the two realizations actually disagree above the
                    // *reference's* own resolution. `cos 60° = 1/2` is exact at both precisions,
                    // so their difference measures the reference, not this bound.
                    let above_reference_noise = (de as i64) > -(deep as i64) + 8;
                    if let (Some(be), true) = (bound.exp2(), angle_exact && above_reference_noise) {
                        let slack = be - de as i64;
                        let e = worst.entry(prec).or_insert((i64::MIN, String::new()));
                        if slack > e.0 {
                            *e = (
                                slack,
                                format!("{what} {num}/{den}°, bound 2^{be} vs error 2^{de}"),
                            );
                        }
                    }
                }
            }
        }
        for (prec, (slack, at)) in &worst {
            eprintln!("[cip] {prec}-bit trig bound slack: 2^{slack}  ({at})");
        }
        for (&prec, (slack, at)) in &worst {
            if prec % 64 == 0 {
                assert!(
                    *slack <= 16,
                    "at {prec} bits the bound is 2^{slack} above the worst observed error ({at}) — \
                     that is a fudge factor wearing a derivation's clothes, not a tight bound"
                );
            }
        }
        // …and the odd rung out proves why stage 4's rungs are multiples of 64: asking for 200
        // bits buys 256, so ~56 bits of the result are paid for and then claimed away.
        let (slack_200, _) = &worst[&200];
        assert!(
            (40..=80).contains(slack_200),
            "expected ~56 bits of unclaimed precision at 200 bits (the word-size round-up), got \
             2^{slack_200}"
        );
    }

    /// Error of a high-precision realization `a` against an exact f64 `b`, as a
    /// power-of-two exponent (`i32::MIN` when exactly equal).
    fn hp_err_exp(a: &BigFloat, b: f64) -> i32 {
        let err = a.sub(&BigFloat::from_f64(b, HP_PREC), HP_PREC, HP_RM);
        if err.is_zero() {
            i32::MIN
        } else {
            err.exponent().unwrap_or(i32::MIN)
        }
    }

    // Small ranges keep checked arithmetic inside i128, so these exercise the
    // algebraic laws, not the overflow path (which its own test above pins).
    prop_compose! {
        fn small_rat()(num in -1000i128..=1000, den in 1i128..=1000) -> Rat {
            Rat::new(num, den).unwrap()
        }
    }

    proptest! {
        /// Rational `+` and `×` are commutative (exact, no drift).
        #[test]
        fn add_and_mul_commute(a in small_rat(), b in small_rat()) {
            prop_assert_eq!(a.checked_add(b), b.checked_add(a));
            prop_assert_eq!(a.checked_mul(b), b.checked_mul(a));
        }

        /// `a + b + c` associates regardless of grouping.
        #[test]
        fn add_associates(a in small_rat(), b in small_rat(), c in small_rat()) {
            let left = a.checked_add(b).and_then(|ab| ab.checked_add(c));
            let right = b.checked_add(c).and_then(|bc| a.checked_add(bc));
            prop_assert_eq!(left, right);
        }

        /// `from_deg` always normalizes into `[0, 360)`.
        #[test]
        fn from_deg_normalizes_into_range(deg in -3600i128..=3600) {
            let a = Angle::from_deg(Rat::from_int(deg)).unwrap();
            prop_assert!(a.deg() >= Rat::from_int(0));
            prop_assert!(a.deg() < Rat::from_int(360));
        }
    }
}
