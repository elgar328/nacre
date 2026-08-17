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
//!   f64 and records the tol on its `Origin`; the definition is never lost, so a
//!   judgement can still realize it at whatever precision it needs.
//! - [`Angle`] — rational degrees, normalized mod-360, with exact accumulation so
//!   a full turn lands back on exactly `0` (no f64 drift — a sketch closes). The
//!   `cos`/`sin` realization crosses into f64 (the irrational boundary); the
//!   90°-family (`0/90/180/270°`) realizes to exact rationals `{0, ±1}` (Niven),
//!   so those rotations of a rational point stay tol 0; and `cos_sin_at`
//!   realize in arbitrary precision (astro-float) for the judgment path.
//!
//! Scope: the exact value engine (the toleranced-sign frame judgment it enabled now lives
//! in `nacre-cip`). [`Rat::from_decimal`] carries this into *construction*: a prism's
//! placement and sweep are done in the rationals the dimensions were **written** as, so
//! `1.1` then `6.6` reaches the same plane as `7.7` (`nacre_ops::exact`). That holds where
//! the sketch frame is exactly orthonormal — a rotated frame's axes are irrational, and
//! there this crate has nothing to offer; `nacre-cip` is what keeps *judgments* sound
//! there. Still deferred: the unified `Scalar { value, tol }` wrapper.
//! **The declare-0 → user-confirmation policy that stood here is retired, not pending**:
//! measurement refuted both halves, and an unprovable sign now leaves as a proved
//! coincidence carrying its evidence, or as a reject named for its cause (`nacre-cip`).
//! Ported from an isolated 2D experiment that verified it first.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
pub mod mag;
pub use mag::Mag;

use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedDiv, CheckedMul, CheckedSub};

use astro_float::{BigFloat, Consts, RoundingMode};
use std::cell::RefCell;
use std::collections::HashMap;

/// Rounding for the high-precision realization layer (astro-float).
///
/// **There is no default precision here, on purpose.** A fixed one used to sit beside this
/// (`HP_PREC = 160`) with `cos_hp`/`sin_hp` reading it, and nothing outside this crate's own
/// tests ever called them — a hand-picked depth waiting to be wired into a judgement whose
/// precision belongs to the *model* (`nacre_cip::judge_precision`). Callers pass `prec`.
pub(crate) const HP_RM: RoundingMode = RoundingMode::ToEven;

thread_local! {
    /// Transcendental-constant cache (π, …) for the high-precision layer.
    ///
    /// **Constants, not results.** `astro_float::Consts` holds `pi/e/ln2/ln10/tenpowers` and
    /// nothing else, and `cos`/`sin` take it only to reach π (degrees → radians, and argument
    /// reduction). So it saves computing π once per call — a real saving, and orthogonal to
    /// [`TRIG`] below, which is what saves the *evaluation*.
    static HP_CONSTS: RefCell<Consts> = RefCell::new(Consts::new().expect("astro-float consts"));

    /// **Realized `cos`/`sin` per `(angle, precision)`** — see [`Angle::cos_sin_bounded`].
    ///
    /// A rotation's realization asks for the *same* angle once per point: a 60-fin fold evaluated
    /// **29,037 of them for 112 distinct `(angle, prec)` pairs**, and a solid turned 4,200 times
    /// re-evaluated its one angle 4,200 times per point. astro-float's `cos`/`sin` run Ziv's loop
    /// (a series, then a retry at more bits when the value sits too close to a rounding boundary),
    /// which is what that repetition was paying for — measured, **22% of a rotated boolean**.
    ///
    /// **The key is a value, not a handle**, which is what makes a process-wide memo sound here:
    /// `Angle` is an exact rational number of degrees, so two models asking for the same angle are
    /// asking the same question. (A cache keyed by a `Handle` could not be shared this way — a
    /// handle is an index into *one* model's store.) And the result is a pure function of the key,
    /// so the memo cannot move an answer; only how often Ziv's loop runs. `prec` is *in* the key,
    /// so a model whose judging precision grows simply lands on a different entry rather than
    /// reading one realized at the wrong depth.
    static TRIG: RefCell<HashMap<(Angle, usize), TrigAt>> = RefCell::new(HashMap::new());

    /// **How far this platform's `f64` cos/sin land from the truth** — see
    /// [`Angle::realization_error_of`].
    ///
    /// ★ **No precision in the key, and that is the point.** [`TRIG`] is keyed by
    /// `(angle, precision)` because it holds a value *realized at* a depth. This holds how wrong
    /// the **f64** realization is, and f64 has one precision. The `P` used to measure it says how
    /// finely the answer is read, not what the answer is.
    ///
    /// ★★★ **The realized pair is in the key, not just the angle.** The quantity is how far *this
    /// `(cos, sin)`* sits from the truth, so the pair that was measured is part of the question.
    /// Keying on the angle alone would hand a caller the error of a realization it did not use —
    /// which is not hypothetical: `(deg.to_f64() * PI / 180.0).sin()` is measured to differ by
    /// 1 ulp between a debug and a release build, and between two call sites within one release
    /// build (LLVM folds the literal-angle site at compile time and leaves the other to libm).
    ///
    /// ★★ **So a second entry under one angle is a signal, not waste**: it says two realizations of
    /// that angle are live in this process. The tests read the count for exactly that reason.
    static F64_ERR: RefCell<HashMap<RealizedAt, (f64, f64)>> = RefCell::new(HashMap::new());

    /// **The f64 realization of an angle, keyed by the angle alone** — see [`Angle::cos_sin_f64`].
    ///
    /// ★ Unlike [`F64_ERR`], the pair is *not* in the key, and it must not be: the value here **is**
    /// the pair, and it is the correctly rounded one, which is unique. This memo therefore cannot
    /// change an answer — it only stops the arbitrary-precision realization from running per vertex.
    ///
    /// **Measured, which is why it is here from the start**: one realization costs ~32µs against
    /// libm's 22ns, and `Isometry::apply_point` runs per vertex.
    static F64_TRIG: RefCell<HashMap<Angle, (f64, f64)>> = RefCell::new(HashMap::new());

    /// How often the f64 realization needed a second, deeper pass, and how often even that left
    /// the rounding undecided — see [`Angle::cos_sin_f64`]. **Read by tests**: a checker that
    /// never fires is indistinguishable from one that is not running.
    static ROUND_ESCALATED: RefCell<(usize, usize)> = const { RefCell::new((0, 0)) };

    /// **Realized `1/√v` per `(rational, precision)`** — see [`inv_sqrt_bounded`].
    ///
    /// The same argument that makes [`TRIG`] sound applies unchanged: the key is a *value*, not a
    /// handle, so two models asking for the same squared length are asking the same question; the
    /// result is a pure function of the key, so the memo cannot move an answer; and `prec` is in
    /// the key, so a model judged more deeply lands on a different entry rather than reading one
    /// realized too shallowly.
    static INV_SQRT: RefCell<HashMap<(Rat, usize), (BigFloat, Mag)>> =
        RefCell::new(HashMap::new());

    /// **The f64 realization of `1/√v`, keyed by the rational alone** — see [`inv_sqrt_f64`].
    ///
    /// As with [`F64_TRIG`], the value is the correctly rounded one, which is unique, so this memo
    /// is a cost question only.
    static F64_INV_SQRT: RefCell<HashMap<Rat, f64>> = RefCell::new(HashMap::new());

    /// How often `1/√v`'s f64 realization needed the deeper rung, and how often even that left the
    /// rounding undecided — the counterpart of [`ROUND_ESCALATED`], read by tests for the same
    /// reason.
    static INV_SQRT_ESCALATED: RefCell<(usize, usize)> = const { RefCell::new((0, 0)) };
}

/// One `(angle, precision)` realization: `(cos, sin, |Δcos|, |Δsin|)`.
type TrigAt = (BigFloat, BigFloat, Mag, Mag);

/// **The `f64` nearest the true value that `mid ± rad` encloses — or `None` when `mid ± rad` is
/// not narrow enough to say.**
///
/// An arbitrary-precision realization is an *interval*, and rounding its midpoint to 53 bits is
/// rounding an approximation: near a rounding boundary the answer would be the approximation's,
/// not the truth's. So both ends are rounded and compared. Round-to-nearest is **monotonic**, so
/// `lo ≤ v ≤ hi` gives `round(lo) ≤ round(v) ≤ round(hi)` — when the ends agree, the value between
/// them cannot round anywhere else, whatever it is. `None` says "realize deeper and ask again",
/// which is the same shape as `precision_for` escalating a judgement.
///
/// ★★★ **The ends are formed at `prec + 64`, and that line is load-bearing.** `rad` is about
/// `2⁻ᵖʳᵉᶜ` of `mid`, so forming `mid ± rad` *at* `prec` rounds the radius away, both ends collapse
/// onto `mid`, and the comparison then passes for every input — a check that runs, reports success,
/// and verifies nothing. The guard digits keep the perturbation alive.
///
/// **Not for a value whose true magnitude may be zero.** `cos 90°` is exactly `0`, so its interval
/// straddles zero, the two ends have opposite signs, and no precision ever makes them agree — this
/// would return `None` forever. Callers must resolve the exactly-representable cases first;
/// [`Angle::cos_sin_f64`] does that with `try_exact_cos_sin`.
///
/// `pub` because **STEP export** will want the same rounding, and two implementations of "round
/// this interval to f64" is exactly the drift this crate keeps deleting. STEP carries `f64`
/// coordinates and nothing else — not the motion chain, not the three-plane definition — so export
/// is the one place where a point's exact truth must be realized as well as `f64` allows, and it
/// can be *verified* as it goes: this returning `Some` **is** the proof that the written coordinate
/// is the nearest `f64` to the exact one.
///
/// ★★ **Not for the coordinates the kernel works with.** Rounding `px + u·c − v·s` as a whole
/// would make every `WitnessPoint::coord` literally `round(compute_hp)` and take `tol` to a half-ulp, and it
/// was **declined**: that is arbitrary precision *per vertex*, where [`Angle::cos_sin_f64`]'s is
/// per *angle* and memoised, so it would pay at construction for a precision that `WitnessPoint`'s lazy
/// `compute_hp` already buys **only where a judgement actually needs it**. (`nacre-cip` depends on
/// this crate, so that type cannot be named here as a link.) Cheaper ways to shrink that
/// arithmetic (an FMA, a compensated
/// evaluation) stay in `f64` and keep the laziness, so they are the candidates if the term ever
/// needs to move.
pub fn round_to_f64(mid: &BigFloat, rad: Mag, prec: usize) -> Option<f64> {
    if mid.is_nan() || mid.is_inf() {
        return None;
    }
    let p = prec + 64;
    let r = rad_upper_big(rad, p)?;
    let (lo, hi) = (mid.sub(&r, p, HP_RM), mid.add(&r, p, HP_RM));
    let (rlo, rhi) = (to_f64_exact(&lo)?, to_f64_exact(&hi)?);
    (rlo == rhi && !rlo.is_nan()).then_some(rlo)
}

/// `2^exp2()` as a `BigFloat` — an **upper bound** on the radius, exactly representable.
///
/// [`Mag`] is `m · 2^e` with `m ∈ [0.5, 1)`, so `2^e` is above it; the mantissa is left out
/// because widening the interval can only cost an escalation, never buy a wrong acceptance, and
/// `Mag` does not expose its mantissa. A power of two is exact in `BigFloat` at any precision.
///
/// `None` when `2^e` is outside `f64`'s range. That cannot happen for the radii this crate
/// produces (`prec ≤ 256` against magnitudes above `2⁻¹³³`), and returning `None` rather than
/// silently flushing to zero is what keeps a broken premise from reading as a *tighter* interval.
fn rad_upper_big(rad: Mag, p: usize) -> Option<BigFloat> {
    let Some(e) = rad.exp2() else {
        return Some(BigFloat::from_f64(0.0, p)); // an exact realization: a zero radius is honest
    };
    if !(-1000..=1000).contains(&e) {
        return None;
    }
    Some(BigFloat::from_f64(2f64.powi(e as i32), p))
}

/// `x` rounded to 53 significant bits and read out as the `f64` with those bits.
///
/// Two steps that must not be confused: `set_precision(53, ToEven)` performs the *rounding* (this
/// is the only place a value loses bits), and the assembly below is a pure re-encoding of what
/// that produced. Splitting them is why the caller can round two interval ends and compare.
///
/// Restricted to the normal range on purpose — every caller here holds a `cos`/`sin` of a
/// non-quadrantal rational-degree angle, whose magnitude is between `2⁻¹³³` and `1`, so a subnormal
/// or an overflow means a premise broke rather than an input being unusual. `None` says so.
fn to_f64_exact(x: &BigFloat) -> Option<f64> {
    if x.is_zero() {
        return Some(0.0);
    }
    let mut v = x.clone();
    v.set_precision(53, HP_RM).ok()?;
    let (words, _bits, sign, e, _inexact) = v.as_raw_parts()?;
    // Most significant word last (`Mantissa::to_u64` reads `m[len - 1]`), and the mantissa is
    // normalized so that top bit is set. Assembled across words because `Word` is `u32` on
    // 32-bit targets — wasm is one, and it is where this kernel actually ships.
    const WB: u32 = astro_float::WORD_BIT_SIZE as u32;
    let mut top: u64 = 0;
    let mut filled = 0u32;
    for w in words.iter().rev() {
        if filled >= 64 {
            break;
        }
        // `checked_shl` rather than `<<`: on a 64-bit target `WB` *is* 64, and the shift would be
        // undefined. It cannot actually run there — `filled` reaches 64 after one word and the
        // loop stops — but the expression still has to be well-formed for the compiler.
        // The widening is a no-op where `Word` is already `u64` and required where it is `u32`;
        // clippy sees only the target it is run on, and dropping it would fail to compile the
        // other one.
        #[allow(clippy::useless_conversion)]
        {
            top = top.checked_shl(WB).unwrap_or(0) | u64::from(*w);
        }
        filled += WB;
    }
    top <<= 64 - filled.min(64);
    // `e` is astro-float's exponent for a mantissa in `[0.5, 1)`; f64's biased exponent for the
    // same value is `e - 1 + 1023`. Anything outside the normal range is a broken premise.
    let biased = i64::from(e) + 1022;
    if !(1..=2046).contains(&biased) {
        return None;
    }
    let sign_bit = u64::from(sign == astro_float::Sign::Neg) << 63;
    // Drop the implicit leading 1, then take the 52 stored bits.
    Some(f64::from_bits(
        sign_bit | ((biased as u64) << 52) | ((top << 1) >> 12),
    ))
}

/// An angle **together with one f64 realization of it** — `(angle, cos.to_bits(), sin.to_bits())`.
/// The bits, not the floats, because the key has to be `Hash` and `Eq`.
type RealizedAt = (Angle, u64, u64);

/// How many entries [`TRIG`] holds — for the tests that pin the memo actually memoizes.
///
/// **The count, not a hit tally.** A hit rate cannot tell "the memo works" from "the memo is
/// fragmenting": a caller that spelled one angle two ways would show high hits while paying twice
/// for every angle. Every miss inserts, so the entry count is the direct reading.
#[cfg(test)]
fn trig_entries() -> usize {
    TRIG.with_borrow(|t| t.len())
}

/// The result of an orientation judgment (§CIP) — the shared sign vocabulary used by both
/// the 2D and 3D toleranced-sign judges (now in `nacre-cip`) and their downstream consumers.
/// A cross-cutting "judgment result" carried here as a fundamental value (a candidate to
/// split into its own vocabulary type later).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orient {
    Positive,
    Negative,
    /// Zero — the two things are on the same side of nothing.
    ///
    /// **Whether that was *proved* is not carried here.** This type is the sign the geometry
    /// consumes; a proved zero, a coincidence established within the judging standard, and a
    /// judgement that ran out of precision all reach it looking alike. `nacre_cip::Decision`
    /// is what keeps them apart — and what turns the last of them into a named reject rather
    /// than a silent merge.
    Zero,
}

/// A rational scalar (exact, tol 0). Arithmetic returns `None` on i128 overflow
/// so the caller sees the downgrade trigger explicitly; on overflow the kernel
/// switches that value's cache to f64 and tags its `Origin` with the resulting
/// tol, while the definition (the op-log of input rationals) is preserved. Overflow is far from normal use — adversarial coprime-denominator
/// accumulation reaches it near the i128 ceiling (~122 bits).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    /// The **nearest** f64 to this rational, ties to even — the same answer IEEE-754
    /// would give if it could divide exactly (for the downgrade path, and for export).
    ///
    /// `numer as f64 / denom as f64` is *not* that answer: past 2⁵³ each conversion
    /// rounds before the division sees it, and two roundings do not compose into one.
    /// The error is a matter of significant digits, not of magnitude — a 17-digit
    /// value lands off-by-an-ulp whether it is 1e-3 or 1e17. That is fine for a cache
    /// but not for the round-trip `from_decimal(x).to_f64() == x`, which is exactly the
    /// property that lets a decimal dimension be recovered without moving the model.
    pub fn to_f64(self) -> f64 {
        let (n, d) = (*self.0.numer(), *self.0.denom());
        // `Ratio` keeps the sign in the numerator and the denominator positive.
        let neg = n < 0;
        let (n, d) = (n.unsigned_abs(), d as u128);
        // Both terms exact as f64, so the hardware division is already the nearest
        // value — and this covers nearly every rational the kernel actually holds.
        if n < (1 << 53) && d < (1 << 53) {
            let q = n as f64 / d as f64;
            return if neg { -q } else { q };
        }
        let q = nearest_f64(n, d);
        if neg { -q } else { q }
    }

    /// The rational **the decimal spelled**, rather than the one the f64 holds.
    ///
    /// These differ, and the difference is the whole point. [`try_from_f64`] lifts the
    /// binary value — `1.1` becomes `2476979795053773/2251799813685248`, not `11/10` —
    /// so exact arithmetic on lifted f64 reproduces the binary drift exactly instead of
    /// removing it: `try_from_f64(1.1) + try_from_f64(6.6) != try_from_f64(7.7)`, while
    /// `11/10 + 66/10 = 77/10` and realizes back to exactly `7.7`. A dimension a person
    /// typed, and one the kernel derived from it, agree again.
    ///
    /// The decimal is the **shortest that reads back as `x`** — a deterministic normal
    /// form, not a guess at intent: distinct f64 map to distinct rationals and equal
    /// f64 always to the same one. Decimal is the choice only because people type
    /// decimal. `None` for non-finite input, or when the power of ten overflows i128,
    /// in which case the caller keeps its f64 path.
    ///
    /// Measured, that limit sits at `1e38` above and — for a value carrying all 17
    /// significant digits — `1e-22` below; a short decimal reaches down to `1e-38`.
    /// CAD dimensions live nowhere near either edge.
    ///
    /// [`try_from_f64`]: Rat::try_from_f64
    pub fn from_decimal(x: f64) -> Option<Self> {
        if !x.is_finite() {
            return None;
        }
        // `{:e}`, not `{}`: both are the shortest round-tripping spelling, but `{}`
        // writes `1e300` out in full, so the mantissa is not bounded and the i128
        // verdict comes only after 301 characters. Here the mantissa is at most 17
        // digits and the exponent arrives separately.
        let s = format!("{x:e}");
        let (mantissa, exp) = s.split_once('e')?;
        let exp: i32 = exp.parse().ok()?;
        let (int, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        let digits: i128 = format!("{int}{frac}").parse().ok()?;
        // value = digits · 10^(exp − frac.len())
        let scale = exp - frac.len() as i32;
        if scale >= 0 {
            Rat::new(digits.checked_mul(10i128.checked_pow(scale as u32)?)?, 1)
        } else {
            Rat::new(digits, 10i128.checked_pow(-scale as u32)?)
        }
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

/// The canonical rational plane with normal `n` through `point` — `n·x − n·point = 0`.
///
/// See [`canonical_plane_coeffs`] for why the result is canonicalized and why the inputs must be
/// rational by construction rather than lifted from f64 coefficients. `None` on `i128` overflow.
pub fn plane_from_point_normal(n: [Rat; 3], point: [Rat; 3]) -> Option<[Rat; 4]> {
    let mut d = Rat::from_int(0);
    for i in 0..3 {
        d = d.checked_sub(n[i].checked_mul(point[i])?)?;
    }
    canonical_plane_coeffs([n[0], n[1], n[2], d])
}

/// The canonical rational plane through three points, normal `(b − a) × (c − a)` — the exact twin
/// of `nacre_geom::Plane::through_points`, so the two describe the same plane with the same
/// orientation. `None` if the points are collinear (no plane) or on `i128` overflow.
pub fn plane_through_points(a: [Rat; 3], b: [Rat; 3], c: [Rat; 3]) -> Option<[Rat; 4]> {
    let d = |p: [Rat; 3], q: [Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            p[0].checked_sub(q[0])?,
            p[1].checked_sub(q[1])?,
            p[2].checked_sub(q[2])?,
        ])
    };
    let (u, v) = (d(b, a)?, d(c, a)?);
    let term = |i: usize, j: usize| u[i].checked_mul(v[j])?.checked_sub(u[j].checked_mul(v[i])?);
    let n = [term(1, 2)?, term(2, 0)?, term(0, 1)?];
    if n.iter().all(|c| *c == Rat::from_int(0)) {
        return None; // collinear
    }
    plane_from_point_normal(n, a)
}

/// The canonical name of a rational plane — primitive integer coefficients, sign-fixed
/// (see [`canonical_plane_coeffs`] for the canonical form). One vessel, two widths.
///
/// ★ **Normalization invariant: a value that fits `i128` is ALWAYS stored `Narrow`** — the only
/// constructor ([`plane_name_exact`]) enforces it, so two statements of one plane are
/// structurally equal (`==`/`Hash`) across representations. That is what makes this the
/// interning key: identity never depends on which route derived the name.
///
/// ★★ `Wide` carries **identity only**. Arithmetic shortcuts (frames, Shewchuk transports,
/// `base_rat`) read [`PlaneName::narrow`] and decline on `None`, exactly as they declined on a
/// missing name before — a wide name does not open a frame (measured refutation in
/// `docs/truth-and-cache.md`'s refuted-answers table).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum PlaneName {
    Narrow([Rat; 4]),
    Wide([num_bigint::BigInt; 4]),
}

impl PlaneName {
    /// The `i128` form when there is one — what the exact shortcuts (Shewchuk integer
    /// predicates, frame derivation, `Isometry` transport) consume. `None` (`Wide`) means those
    /// shortcuts decline and judgment takes the general route: slower, never wrong.
    #[inline]
    pub fn narrow(&self) -> Option<&[Rat; 4]> {
        match self {
            PlaneName::Narrow(c) => Some(c),
            PlaneName::Wide(_) => None,
        }
    }

    /// The canonical coefficients as integers, **whichever width holds them** — the one vessel
    /// the integer sign predicates ([`int_plane_side`], [`int_cmp_coord`], [`int_dir_sign`])
    /// consume. A canonical narrow name is a primitive integer vector (the invariant
    /// [`plane_name_exact`] normalizes to), so `numer()` is the value; a wide name already is
    /// the integers.
    pub fn coeff_ints(&self) -> [num_bigint::BigInt; 4] {
        match self {
            PlaneName::Narrow(c) => {
                debug_assert!(c.iter().all(|r| r.denom() == 1));
                core::array::from_fn(|k| num_bigint::BigInt::from(c[k].numer()))
            }
            PlaneName::Wide(c) => c.clone(),
        }
    }
}

/// **The plane three points name, computed so that the arithmetic on the way cannot lose it.**
///
/// The same answer [`plane_through_points`] gives, in a vessel that always holds it: `Narrow`
/// when the canonical answer fits `i128`, `Wide` (arbitrary-precision integers) when it does
/// not. `None` means exactly one thing — the points are collinear and name no plane. It is
/// never a shrug about an intermediate, and no longer one about the answer's width either.
///
/// ★★★★★ **The distinction is not academic — it was most of the failures.** `plane_through_points`
/// works in `Rat`, so `(b − a) × (c − a)` multiplies the points' denominators together and
/// `−n · a` multiplies once more. Measured on the failing cases, that peak needs **271 bits** while
/// the canonical answer, after the content is divided out, comes back to about **50**. Recomputing
/// those at unbounded precision reproduced the stored name **1,545 times out of 1,545** — the two
/// vectors differ by a scalar factor and `canonical_plane_coeffs` removes exactly that freedom. So
/// what overflowed was the road, not the destination.
///
/// ★★★ **Why a second function rather than widening the first.** `plane_through_points` is what a
/// *caller* uses to state a plane they wrote down ([`crate::plane_frame`]'s callers, a sketch, a
/// box's corners); this is what the kernel uses to **derive** the name of a plane it already holds
/// three exact points for. Widening the shared one would change both at once, and they are
/// different propositions with different evidence. (Making the caller's path exact too is worth
/// doing and is its own measurement.)
///
/// ★ **Cost is paid only on the fallback.** The `Rat` route runs first and is the answer whenever
/// it fits; `BigInt` is reached on the rest. Nothing here is on a boolean's inner loop — a plane is
/// named once per `Model::push_surface_with_coeffs`.
pub fn plane_name_exact(a: [Rat; 3], b: [Rat; 3], c: [Rat; 3]) -> Option<PlaneName> {
    plane_through_points(a, b, c)
        .map(PlaneName::Narrow)
        .or_else(|| plane_name_big(a, b, c))
}

/// [`plane_name_exact`]'s unbounded arm, always taken — the differential test needs to call it on
/// inputs the `Rat` route handles, which it cannot do through the filter.
///
/// **The same plane [`plane_through_points`] computes**, reached by clearing each point's
/// denominators first so the arithmetic is integer throughout. The two agreeing wherever the narrow
/// one answers is the correctness argument, and `the_wide_derivation_answers_what_the_narrow_one_does`
/// is what holds it. `None` is collinearity, or a canonical component that no longer fits `Rat`.
pub(crate) fn plane_name_big(a: [Rat; 3], b: [Rat; 3], c: [Rat; 3]) -> Option<PlaneName> {
    use num_bigint::BigInt;
    use num_integer::Integer;

    // ★★★★ **Integers, not rationals.** The obvious spelling is `Ratio<BigInt>`, mirroring
    // `plane_through_points` term for term — and that is how this started. But `Ratio` reduces by a
    // gcd on *every* multiply and subtract, and there are a dozen of them, so the reduction work
    // dominates. Clearing each point's denominators once up front leaves plain integer arithmetic
    // and exactly one gcd at the end, where the content has to come out anyway.
    //
    // ★ **A plane is scale-free, which is what makes this legal.** Each point is scaled by its own
    // `D_k`, so `u'` and `v'` are the true edges times `D_a·D_b` and `D_a·D_c`; their cross product
    // is the true normal times a positive factor, and the canonical form divides all of it out.
    // The `d` term is `−N·a = −(N·P_a)/D_a`, so scaling the whole 4-vector by `D_a` clears it.
    let lift = |p: [Rat; 3]| -> ([BigInt; 3], BigInt) {
        let den = p.map(|r| BigInt::from(r.denom()));
        let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
        let num = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&d / &den[i]));
        (num, d)
    };
    plane_name_from_lifted([lift(a), lift(b), lift(c)])
}

/// [`plane_name_big`]'s body after the lift — three points as `(integer coordinates, positive
/// denominator scale)` pairs. Separate so [`plane_name_from_meets`] can reach it with points
/// that never were `Rat` triples ([`MeetPoint::Wide`]).
fn plane_name_from_lifted(
    pts: [([num_bigint::BigInt; 3], num_bigint::BigInt); 3],
) -> Option<PlaneName> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{ToPrimitive, Zero};

    let [(pa, da), (pb, db), (pc, dc)] = pts;
    let edge = |q: &[BigInt; 3], dq: &BigInt| -> [BigInt; 3] {
        core::array::from_fn(|i| &q[i] * &da - &pa[i] * dq)
    };
    let (u, v) = (edge(&pb, &db), edge(&pc, &dc));
    let term = |i: usize, j: usize| &u[i] * &v[j] - &u[j] * &v[i];
    let n = [term(1, 2), term(2, 0), term(0, 1)];
    if n.iter().all(Zero::is_zero) {
        return None; // collinear
    }
    let dot: BigInt = (0..3).map(|i| &n[i] * &pa[i]).sum();
    let mut num = [&n[0] * &da, &n[1] * &da, &n[2] * &da, -dot];

    // Divide out the content and fix the sign of the first nonzero — `canonical_plane_coeffs`'
    // steps ② and ③; ① is already done, since these are integers.
    let g = num.iter().fold(BigInt::zero(), |g, x| g.gcd(x));
    if g.is_zero() {
        return None; // the zero vector is not a plane
    }
    for x in &mut num {
        *x /= &g;
    }
    if num
        .iter()
        .find(|x| !x.is_zero())
        .is_some_and(|x| *x < BigInt::zero())
    {
        for x in &mut num {
            *x = -&*x;
        }
    }

    // ★ The normalization invariant lives here: narrow whenever the canonical answer fits
    // `i128`, `Wide` only when it does not — so equal planes are structurally equal whichever
    // route derived them. (This used to be the one honest failure; now it is the fork.)
    let mut out = [Rat::from_int(0); 4];
    for (o, x) in out.iter_mut().zip(&num) {
        match x.to_i128() {
            Some(v) => *o = Rat::from_int(v),
            None => return Some(PlaneName::Wide(num)),
        }
    }
    Some(PlaneName::Narrow(out))
}

/// **The canonical name of the plane through three meeting points, at whatever width the points
/// needed** — what lets a datum through [`MeetPoint::Wide`] vertices keep a name (truth-and-cache
/// open item 17). The same lift-and-join [`plane_name_big`] performs, with each point's
/// denominators cleared from whichever vessel holds them; a plane is scale-free per point, so the
/// per-point scale cannot move the canonical answer.
///
/// `None` means the points are collinear and name no plane — width is never a cause here, on
/// either side: wide points are lifted the same way, and a canonical answer too wide for `Rat`
/// comes back [`PlaneName::Wide`].
pub fn plane_name_from_meets(points: [&MeetPoint; 3]) -> Option<PlaneName> {
    use num_bigint::BigInt;
    use num_integer::Integer;

    let lift = |p: &MeetPoint| -> ([BigInt; 3], BigInt) {
        match p {
            MeetPoint::Narrow(p) => {
                let den = p.map(|r| BigInt::from(r.denom()));
                let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
                let num = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&d / &den[i]));
                (num, d)
            }
            MeetPoint::Wide(p) => {
                let d = p.iter().fold(BigInt::from(1), |l, (_, den)| l.lcm(den));
                let num = core::array::from_fn(|i| &p[i].0 * (&d / &p[i].1));
                (num, d)
            }
        }
    };
    plane_name_from_lifted([lift(points[0]), lift(points[1]), lift(points[2])])
}

/// The sign of the 2-D orientation determinant `(b − a) × (c − a)` — positive when `abc` turns
/// counter-clockwise, zero when collinear. The exact-`Rat` twin of `nacre_predicates::orient2d`
/// (same convention), for coordinates that *are* rational rather than f64.
///
/// ★ **Total.** Unlike every other `Rat` derivation here, this cannot decline: a sign consumed by
/// a simplicity check that failed open on overflow would read "no intersection" where there is
/// one — silent-wrong, the exact shape [`plane_name_exact`] was built to kill. So the `Rat` route
/// runs first and an overflow falls through to integer arithmetic that cannot overflow.
///
/// ★ **The `BigInt` arm is load-bearing, not a corner case.** [`Rat::from_decimal`] denominators
/// reach `10^38`, and the determinant multiplies two coordinate *differences* — two 17-digit
/// dimensions already push the product denominator past `i128`. Same denominator-clearing shape
/// as [`plane_name_big`]: each point is scaled by its own positive denominator, which multiplies
/// the determinant by `Da²·Db·Dc > 0` and therefore cannot move the sign.
pub fn orient2d_rat(a: [Rat; 2], b: [Rat; 2], c: [Rat; 2]) -> i8 {
    let narrow = || -> Option<i8> {
        let u = [b[0].checked_sub(a[0])?, b[1].checked_sub(a[1])?];
        let v = [c[0].checked_sub(a[0])?, c[1].checked_sub(a[1])?];
        let det = u[0]
            .checked_mul(v[1])?
            .checked_sub(u[1].checked_mul(v[0])?)?;
        Some(rat_sign(det))
    };
    narrow().unwrap_or_else(|| orient2d_big(a, b, c))
}

/// `-1 / 0 / +1` of a rational.
fn rat_sign(x: Rat) -> i8 {
    match x.cmp(&Rat::from_int(0)) {
        core::cmp::Ordering::Less => -1,
        core::cmp::Ordering::Equal => 0,
        core::cmp::Ordering::Greater => 1,
    }
}

/// [`orient2d_rat`]'s unbounded arm, always taken — the differential test needs to call it on
/// inputs the `Rat` route handles, which it cannot do through the filter (the 2-D analogue of
/// [`plane_name_big`]'s arrangement).
pub(crate) fn orient2d_big(a: [Rat; 2], b: [Rat; 2], c: [Rat; 2]) -> i8 {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;

    // Clear each point's denominators once (positive by `Ratio`'s reduction invariant), then the
    // arithmetic is plain integers and the determinant is the true one times `Da²·Db·Dc`.
    let lift = |p: [Rat; 2]| -> ([BigInt; 2], BigInt) {
        let den = p.map(|r| BigInt::from(r.denom()));
        let d = den[0].lcm(&den[1]);
        let num = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&d / &den[i]));
        (num, d)
    };
    let ((pa, da), (pb, db), (pc, dc)) = (lift(a), lift(b), lift(c));
    let edge = |q: &[BigInt; 2], dq: &BigInt| -> [BigInt; 2] {
        core::array::from_fn(|i| &q[i] * &da - &pa[i] * dq)
    };
    let (u, v) = (edge(&pb, &db), edge(&pc, &dc));
    match (&u[0] * &v[1] - &u[1] * &v[0]).sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}

/// **The sign of a plane's residual at a rational point** — `sign(a·x + b·y + c·z + d)` for the
/// plane a [`PlaneName`] names. `0` means the point lies exactly on the plane.
///
/// ★ **Total**, like [`orient2d_rat`] and for the same reason: this is what a *validating
/// constructor* consumes (is a caller's stated sketch origin on the plane they picked?), and a
/// check that failed open on overflow would accept an off-plane origin — silent-wrong. The `Rat`
/// route runs first; an overflow, or a `Wide` name, falls through to integer arithmetic that
/// cannot overflow.
///
/// The sign convention is the name's own: which side is positive depends on how the canonical
/// coefficients came out, so callers should compare against `0`, not against each other across
/// planes.
pub fn plane_residual_sign(name: &PlaneName, p: [Rat; 3]) -> i8 {
    match name {
        PlaneName::Narrow(c) => {
            let narrow = || -> Option<i8> {
                let mut acc = c[3];
                for k in 0..3 {
                    acc = acc.checked_add(c[k].checked_mul(p[k])?)?;
                }
                Some(rat_sign(acc))
            };
            narrow().unwrap_or_else(|| {
                // A canonical narrow name is a primitive integer vector (the invariant
                // `plane_name_exact` normalizes to), so `numer()` is the value.
                debug_assert!(c.iter().all(|r| r.denom() == 1));
                residual_sign_big(&c.map(|x| num_bigint::BigInt::from(x.numer())), p)
            })
        }
        PlaneName::Wide(c) => residual_sign_big(c, p),
    }
}

/// [`plane_residual_sign`]'s unbounded arm — integer coefficients (which both name forms reduce
/// to), a rational point. Clearing the point's denominators multiplies the residual by a
/// positive factor, which cannot move the sign.
fn residual_sign_big(ci: &[num_bigint::BigInt; 4], p: [Rat; 3]) -> i8 {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;
    let den = p.map(|r| BigInt::from(r.denom()));
    let l = den.iter().fold(BigInt::from(1), |acc, x| acc.lcm(x));
    let num: [BigInt; 3] = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&l / &den[i]));
    let res: BigInt = (0..3).map(|i| &ci[i] * &num[i]).sum::<BigInt>() + &ci[3] * &l;
    match res.sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}

/// `-1 / 0 / +1` of a `BigInt`.
fn big_sign(x: &num_bigint::BigInt) -> i8 {
    use num_bigint::Sign;
    match x.sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}

/// `det3` over borrowed `BigInt` entries — the same cofactor expansion
/// `nacre_predicates::det3` evaluates through `Expansion`, exact by being integer arithmetic
/// rather than by tracking roundoff.
fn det3_big(m: [[&num_bigint::BigInt; 3]; 3]) -> num_bigint::BigInt {
    let minor = |p: &num_bigint::BigInt,
                 q: &num_bigint::BigInt,
                 r: &num_bigint::BigInt,
                 t: &num_bigint::BigInt| p * q - r * t;
    m[0][0] * minor(m[1][1], m[2][2], m[1][2], m[2][1])
        - m[0][1] * minor(m[1][0], m[2][2], m[1][2], m[2][0])
        + m[0][2] * minor(m[1][0], m[2][1], m[1][1], m[2][0])
}

/// The implicit point `∩(p₀, p₁, p₂)` by Cramer, over integer plane coefficients: numerators
/// `(Dx, Dy, Dz)` over the common denominator `D`. The integer twin of `nacre_predicates`'
/// `cramer`, for coefficients too wide for any `f64`-chunk representation — an `Expansion`'s
/// pieces are `f64`s, so its exponent range ends near `2¹⁰²³` while a wide canonical name is
/// measured out to `~2²²⁹¹`.
fn cramer_big(p: [&[num_bigint::BigInt; 4]; 3]) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
    let rhs = [-&p[0][3], -&p[1][3], -&p[2][3]];
    let row = |r: usize, k: usize| {
        let mut m = [&p[r][0], &p[r][1], &p[r][2]];
        m[k] = &rhs[r];
        m
    };
    let col_replaced = |k: usize| det3_big([row(0, k), row(1, k), row(2, k)]);
    (
        [col_replaced(0), col_replaced(1), col_replaced(2)],
        det3_big(p.map(|r| [&r[0], &r[1], &r[2]])),
    )
}

/// **Which side of the plane `j` the implicit point `∩(p₀, p₁, p₂)` lies on**, over integer
/// coefficients — `sign(j·[Dvec : D]) · sign(D)`, the integer twin of
/// `nacre_predicates::indirect_plane_side` (same convention: `+1` on the side `j`'s normal
/// points to, and the sign is the plane's own, so a face whose stored normal opposes its
/// outward direction applies that relation itself).
///
/// Row-linear in each input, so any per-row **positive** scale — clearing a denominator,
/// un-reducing a canonical vector — cannot move the answer; negating a `p` row cannot either
/// (it negates `D` and the dot together), while negating `j` negates the answer. `D = 0`
/// (the three planes meet in no point) returns `0` with no branch, as in the twin.
pub fn int_plane_side(p: [&[num_bigint::BigInt; 4]; 3], j: &[num_bigint::BigInt; 4]) -> i8 {
    let (dvec, d) = cramer_big(p);
    let side: num_bigint::BigInt =
        (0..3).map(|i| &j[i] * &dvec[i]).sum::<num_bigint::BigInt>() + &j[3] * &d;
    big_sign(&side) * big_sign(&d)
}

/// The exact sign of `a[axis] − b[axis]` for two implicit points over integer coefficients —
/// `sign(Na·Db − Nb·Da) · sign(Da) · sign(Db)`, the integer twin of
/// `nacre_predicates::indirect_cmp_coord`. Orientation-invariant: negating any row negates a
/// numerator and its denominator together.
///
/// **Precondition:** both triples meet in a point (`Da, Db ≠ 0`), as in the twin.
pub fn int_cmp_coord(
    a: [&[num_bigint::BigInt; 4]; 3],
    b: [&[num_bigint::BigInt; 4]; 3],
    axis: usize,
) -> i8 {
    debug_assert!(axis < 3, "int_cmp_coord: axis must be 0, 1 or 2");
    let (na, da) = cramer_big(a);
    let (nb, db) = cramer_big(b);
    debug_assert!(
        big_sign(&da) != 0 && big_sign(&db) != 0,
        "int_cmp_coord: degenerate three-plane input (D = 0)"
    );
    big_sign(&(&na[axis] * &db - &nb[axis] * &da)) * big_sign(&da) * big_sign(&db)
}

/// `sign(det[n₀; n₁; n₂])` over the three planes' integer normals — how the line `p ∩ a` runs
/// relative to plane `b`, the integer twin of `nacre_predicates::det3_sign` on plane normals
/// (`d` never enters). Direction-sensitive in every row: negating one normal negates the
/// answer.
pub fn int_dir_sign(n: [&[num_bigint::BigInt; 4]; 3]) -> i8 {
    big_sign(&det3_big(n.map(|r| [&r[0], &r[1], &r[2]])))
}

/// **The world origin projected onto a rational plane** — `p = (−d / n·n) · n` for
/// `a·x + b·y + c·z + d = 0`.
///
/// This is the point of the plane closest to `(0, 0, 0)`, and it is what a face's sketch frame
/// takes for its origin. The alternative — the face's area centroid — is computed in `f64` from the
/// face's own vertices, and lifting *that* back into a rational makes a rounded cache into the
/// truth, which is the one thing construction here must never do.
///
/// ★★★ **The answer does not depend on how the plane is spelled.** A plane has a family of
/// coefficient vectors and this formula is invariant across all of them:
///
/// ```text
/// sign:  ( −(−d) / n·n )·(−n)   =  ( −d / n·n )·n
/// scale: ( −λd / λ²(n·n) )·λn   =  ( −d / n·n )·n
/// ```
///
/// So it is a function of the *plane*, not of the vector describing it — which matters because
/// [`canonical_plane_coeffs`] deliberately carries no direction (`[0,0,1,−3]` and `[0,0,−1,3]`
/// canonicalize together), and because two faces of one plane may hold it either way round.
/// Canonicalizing first is therefore not needed for correctness, only for overflow headroom.
///
/// `None` when `n·n = 0` — the coefficients are not a plane — or on `i128` overflow, which is the
/// kernel's ordinary demotion signal: the caller keeps its f64 path.
pub fn plane_origin_projection(coeffs: [Rat; 4]) -> Option<[Rat; 3]> {
    let zero = Rat::from_int(0);
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let mut nn = zero;
    for c in n {
        nn = nn.checked_add(c.checked_mul(c)?)?;
    }
    if nn == zero {
        return None;
    }
    let neg_d = zero.checked_sub(coeffs[3])?;
    // `Rat` exposes no division — the sign and scale invariance above is what makes one here
    // legitimate, and `Ratio`'s checked division reports the overflow the rest of the crate does.
    let mut out = [zero; 3];
    for i in 0..3 {
        out[i] = Rat(neg_d.checked_mul(n[i])?.0.checked_div(&nn.0)?);
    }
    Some(out)
}

/// **The point where three rational planes meet**, exactly — the `Rat` Cramer route first, and
/// where an intermediate overflows `i128`, the same answer through [`three_planes_big`]'s
/// integer core, denominators cleared once per row (a plane row is scale-free, and scaling a
/// row of a linear system leaves its solution). Input rows are `[a, b, c, d]` for
/// `a·x + b·y + c·z + d = 0`.
///
/// `None` means exactly two things, **neither of them the arithmetic**: the determinant is zero
/// (no unique point — parallel or line-sharing planes), or the point itself does not fit `Rat`
/// ([`MeetPoint::Wide`]). It used to also mean "an intermediate overflowed" — the rational
/// cofactor expansion builds `a.num·b.den ± b.num·a.den` before it can reduce — and that
/// conflation silently cost a decimal-framed tool its whole class reuse: its constructed
/// corners solve to points that fit `Rat` (measured 8/8), but the road there overflowed
/// (truth-and-cache open item 0). The fallback runs only on the decline path, so the narrow
/// route's cost and answers are untouched.
///
/// This is what replaces a dissolved sketch-frame base vertex: the frame-shared triple of a
/// prism corner, solved in the frame the planes are stated in, is the corner's exact base —
/// measured bit-identical to the stored base-and-replay road (8/8, `docs/dev-log.md`).
pub fn three_planes_rat(p: [[Rat; 4]; 3]) -> Option<[Rat; 3]> {
    three_planes_rat_narrow(p).or_else(|| {
        use num_bigint::BigInt;
        use num_integer::Integer;
        let lift = |row: &[Rat; 4]| -> [BigInt; 4] {
            let den: [BigInt; 4] = core::array::from_fn(|i| BigInt::from(row[i].denom()));
            let l = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
            core::array::from_fn(|i| BigInt::from(row[i].numer()) * (&l / &den[i]))
        };
        three_planes_int([lift(&p[0]), lift(&p[1]), lift(&p[2])])?
            .narrow()
            .copied()
    })
}

/// [`three_planes_rat`]'s narrow route — Cramer over `Rat` with checked arithmetic (S7), whose
/// `None` still conflates "no unique point" with "an intermediate overflowed". That is fine
/// *here*: the caller above resolves the conflation by asking the integer core.
fn three_planes_rat_narrow(p: [[Rat; 4]; 3]) -> Option<[Rat; 3]> {
    let zero = Rat::from_int(0);
    // 3×3 determinant by cofactor expansion, all checked.
    let det3 = |m: [[Rat; 3]; 3]| -> Option<Rat> {
        let minor = |r0: usize, r1: usize, c0: usize, c1: usize| -> Option<Rat> {
            m[r0][c0]
                .checked_mul(m[r1][c1])?
                .checked_sub(m[r0][c1].checked_mul(m[r1][c0])?)
        };
        m[0][0]
            .checked_mul(minor(1, 2, 1, 2)?)?
            .checked_sub(m[0][1].checked_mul(minor(1, 2, 0, 2)?)?)?
            .checked_add(m[0][2].checked_mul(minor(1, 2, 0, 1)?)?)
    };
    let rhs = [
        zero.checked_sub(p[0][3])?,
        zero.checked_sub(p[1][3])?,
        zero.checked_sub(p[2][3])?,
    ];
    let d = det3([
        [p[0][0], p[0][1], p[0][2]],
        [p[1][0], p[1][1], p[1][2]],
        [p[2][0], p[2][1], p[2][2]],
    ])?;
    if d == zero {
        return None;
    }
    let with_col = |c: usize| -> [[Rat; 3]; 3] {
        let mut m = [
            [p[0][0], p[0][1], p[0][2]],
            [p[1][0], p[1][1], p[1][2]],
            [p[2][0], p[2][1], p[2][2]],
        ];
        for r in 0..3 {
            m[r][c] = rhs[r];
        }
        m
    };
    let mut out = [zero; 3];
    for (i, o) in out.iter_mut().enumerate() {
        *o = Rat(det3(with_col(i))?.0.checked_div(&d.0)?);
    }
    Some(out)
}

/// Where three planes meet, at whatever width the answer needs — [`three_planes_big`]'s result.
///
/// The fork mirrors [`PlaneName`]'s: `Narrow` **whenever the answer fits**, so two routes to one
/// point are structurally equal.
///
/// ★ The reduction is **per coordinate**, because that is what `Rat` — a `Ratio<i128>` each — has
/// to hold. Cramer hands back `detᵢ/det` sharing one denominator, and that shared form is
/// systematically wider than the coordinates it names; measuring it would report the width of the
/// *arithmetic* rather than of the point.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MeetPoint {
    Narrow([Rat; 3]),
    /// `(numerator, denominator)` per coordinate, each in lowest terms with a positive
    /// denominator — at least one of which does not fit `i128`.
    Wide([(num_bigint::BigInt, num_bigint::BigInt); 3]),
}

impl MeetPoint {
    /// The `Rat` form when there is one — `None` (`Wide`) means the point cannot be a `Rat`
    /// triple at all, which is a different fact from [`three_planes_rat`] declining.
    #[inline]
    pub fn narrow(&self) -> Option<&[Rat; 3]> {
        match self {
            MeetPoint::Narrow(p) => Some(p),
            MeetPoint::Wide(_) => None,
        }
    }

    /// **The width `Rat` would have to hold** — the widest of the six magnitudes (numerator and
    /// denominator of each coordinate) *after* the per-coordinate reduction. `≤ 127` is exactly
    /// the `Narrow` condition, so this is the quantity, not a proxy for it.
    ///
    /// ★ Reduced, deliberately: the unreduced `detᵢ/det` Cramer produces is systematically wider
    /// and would report how the answer was computed rather than how wide the answer is.
    pub fn width_bits(&self) -> u64 {
        match self {
            MeetPoint::Narrow(p) => p
                .iter()
                .flat_map(|r| [r.numer(), r.denom()])
                .map(|v| (128 - v.unsigned_abs().leading_zeros()) as u64)
                .max()
                .unwrap_or(0),
            MeetPoint::Wide(p) => p
                .iter()
                .flat_map(|(n, d)| [n.bits(), d.bits()])
                .max()
                .unwrap_or(0),
        }
    }
}

/// **The same meeting point [`three_planes_rat`] computes, without the `i128` ceiling** — the
/// wide twin that [`plane_name_exact`] has had on the plane-*name* side all along and the vertex
/// solve did not.
///
/// [`three_planes_rat`] is `checked_*` throughout, so its `None` conflates two different facts:
/// the point does not fit `Rat`, and an **intermediate** of the rational cofactor expansion
/// overflowed while the point itself would have fit. Only a route without the ceiling can tell
/// those apart — and telling them apart is what decides whether a limit belongs to the *type* or
/// to the *arithmetic*.
///
/// **Integers, not rationals** — [`plane_name_big`]'s argument, applied per row instead of per
/// point: clearing each row's denominators once up front leaves plain integer arithmetic, with
/// the reductions at the end where the content has to come out anyway. ★ **A plane row is
/// scale-free**, so scaling row `k` by its own denominators' lcm leaves the same plane, and
/// scaling a row of a linear system leaves the same solution.
///
/// Takes [`PlaneName`]s rather than `[Rat; 4]` rows so a `Wide` carrier — one whose canonical
/// name no longer fits `Rat` — is solvable too. That population is invisible to
/// [`three_planes_rat`], which reads [`PlaneName::narrow`] and so never sees it.
///
/// `None` for a zero determinant only: parallel or line-sharing planes, no unique point.
pub fn three_planes_big(p: [&PlaneName; 3]) -> Option<MeetPoint> {
    use num_bigint::BigInt;
    use num_integer::Integer;

    let row = |name: &PlaneName| -> [BigInt; 4] {
        match name {
            PlaneName::Wide(c) => c.clone(),
            PlaneName::Narrow(c) => {
                let den: [BigInt; 4] = core::array::from_fn(|i| BigInt::from(c[i].denom()));
                let l = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
                core::array::from_fn(|i| BigInt::from(c[i].numer()) * (&l / &den[i]))
            }
        }
    };
    three_planes_int([row(p[0]), row(p[1]), row(p[2])])
}

/// [`three_planes_big`]'s integer core — Cramer over integer plane rows, one reduction per
/// coordinate, the [`MeetPoint`] fork at the end. Separate so [`three_planes_rat`]'s fallback
/// can reach it with rows that are *not* canonical names (wrapping those in
/// [`PlaneName::Narrow`] would make the type say something false — that variant means "a
/// canonical name", and [`PlaneName::coeff_ints`] asserts its denominators are 1).
fn three_planes_int(m: [[num_bigint::BigInt; 4]; 3]) -> Option<MeetPoint> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{Signed, ToPrimitive, Zero};

    let det3 = |a: &[[BigInt; 3]; 3]| -> BigInt {
        let minor = |r0: usize, r1: usize, c0: usize, c1: usize| -> BigInt {
            &a[r0][c0] * &a[r1][c1] - &a[r0][c1] * &a[r1][c0]
        };
        &a[0][0] * minor(1, 2, 1, 2) - &a[0][1] * minor(1, 2, 0, 2) + &a[0][2] * minor(1, 2, 0, 1)
    };

    let base: [[BigInt; 3]; 3] =
        core::array::from_fn(|r| core::array::from_fn(|c| m[r][c].clone()));
    let rhs: [BigInt; 3] = core::array::from_fn(|r| -&m[r][3]);

    let d = det3(&base);
    if d.is_zero() {
        return None; // no unique point — parallel or line-sharing planes
    }

    // Cramer, then **one reduction per coordinate**: `d` is nonzero, so every gcd here is at
    // least 1 and the division is total.
    let out: [(BigInt, BigInt); 3] = core::array::from_fn(|i| {
        let mut a = base.clone();
        for (r, x) in rhs.iter().enumerate() {
            a[r][i] = x.clone();
        }
        let (mut n, mut q) = (det3(&a), d.clone());
        let g = n.gcd(&q);
        n /= &g;
        q /= &g;
        if q.is_negative() {
            n = -n;
            q = -q;
        }
        (n, q)
    });

    // ★ The same normalization invariant `plane_name_big` states: narrow whenever the canonical
    // answer fits `i128`, `Wide` only when it does not — so equal points are structurally equal
    // whichever route derived them.
    let narrow = (|| {
        let mut r = [Rat::from_int(0); 3];
        for (o, (n, q)) in r.iter_mut().zip(&out) {
            *o = Rat::new(n.to_i128()?, q.to_i128()?)?;
        }
        Some(r)
    })();
    Some(match narrow {
        Some(r) => MeetPoint::Narrow(r),
        None => MeetPoint::Wide(out),
    })
}

/// **A plane pushed `t` along its own unit normal**, exactly — `d′ = d − t·|n|`.
///
/// What a prism's far cap *is*: the face it was raised from, moved out by the sweep. Deriving it
/// this way rather than from the realized geometry is what lets two prisms raised to one height
/// record **one plane** — `7.7` in a single step and `1.1` then `6.6` give the same `d′`, because
/// `11/10 + 66/10` is `77/10` in rationals — and what lets the far cap of a face on a turned solid
/// be exact at all: the coefficients live in that face's own pre-motion frame, where the rotation's
/// irrational numbers never appear.
///
/// ★★ **`None` unless `|n|` is rational**, which is the one thing that can stop this: `n·n` must be
/// a perfect square. It is `1` for a box face and `25` for a `3-4-5` normal; it is `3` for `[1,1,1]`,
/// and there the caller keeps whatever path it had. Also `None` on `i128` overflow and for `n = 0`,
/// which is not a plane.
///
/// The result is canonicalized, so it compares by `==` with any other exact description of the same
/// plane — including the one the sketch frame produces when it is exact, which is what makes the
/// two derivations checkable against each other.
pub fn plane_offset(coeffs: [Rat; 4], t: Rat) -> Option<[Rat; 4]> {
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let mut nn = Rat::from_int(0);
    for c in n {
        nn = nn.checked_add(c.checked_mul(c)?)?;
    }
    if nn == Rat::from_int(0) {
        return None; // not a plane
    }
    let len = rat_sqrt_exact(nn)?;
    canonical_plane_coeffs([
        n[0],
        n[1],
        n[2],
        coeffs[3].checked_sub(t.checked_mul(len)?)?,
    ])
}

/// **A plane's own frame, exactly** — the origin it is drawn from and the raw direction its `u`
/// axis runs along, both rational, plus the normal reduced to its primitive direction.
///
/// ```text
/// origin  = the world origin projected onto the plane
/// n       = [a, b, c] divided by its content (sign kept)
/// u_raw   = ẑ × n = (−b, a, 0),  or  ŷ × n = (c, 0, 0) when the normal is vertical
/// ```
///
/// ★★★ **`u_raw` is not projected and need not be a unit vector.** The general recipe for laying
/// a reference direction into a plane is `(n·n)·ref − (ref·n)·n`, and it is unnecessary here: a
/// cross product is perpendicular to both its arguments, so `ẑ × n` is *already* in the plane.
/// Skipping it is the difference between coefficients that grow cubically and ones that do not —
/// measured over the faces a face-based operation actually targets, projecting left 21% of them
/// inside `i128` and this leaves **100%**.
///
/// ★★★ **The sign is kept, unlike [`canonical_plane_coeffs`].** That function answers *"are these
/// the same plane"*, where direction is noise. A frame's `n` **is** a direction: negating it
/// negates `û` and `ŵ` together, which is a half-turn about `v` — a different frame, not the same
/// one written differently.
///
/// ★★ **The branch is exact, not toleranced**, and it matches `nacre-ops`' f64 `frame_axes` term
/// for term. DXF's arbitrary-axis convention switches on `|n_x| < 1/64` because a float-only
/// kernel cannot ask the real question; `a == 0 && b == 0` is the real question.
///
/// `None` when the plane is degenerate, when the origin projection is not rational, or when the
/// squared lengths the realization needs do not fit `i128` — all three are honest declines that
/// leave a caller on the f64 path it was already on, never a reject.
pub fn plane_frame(coeffs: [Rat; 4]) -> Option<PlaneFrame> {
    let (origin, ref_dir) = plane_frame_default(coeffs)?;
    plane_frame_named(coeffs, origin, ref_dir)
}

/// **Where a plane's frame sits when nobody names it** — `(origin, ref_dir)`.
///
/// The origin is the world origin projected onto the plane and `ref_dir` is `ẑ × n`
/// (`ŷ × n` when the normal is vertical), which is the arbitrary-axis convention a *face* takes.
/// A caller who names their own sketch origin and `+u` passes those to [`plane_frame_named`]
/// instead — a named plane can insist on axes no derivation would produce, and the script layer's
/// `ZX` (whose `+u` is `+ẑ`, not `ẑ × n = −x̂`) is exactly such a case.
pub fn plane_frame_default(coeffs: [Rat; 4]) -> Option<([Rat; 3], [Rat; 3])> {
    let origin = plane_origin_projection(coeffs)?;
    let n = reduce_direction([coeffs[0], coeffs[1], coeffs[2]])?;
    let zero = Rat::from_int(0);
    let ref_dir = if n[0] == zero && n[1] == zero {
        [n[2], zero, zero]
    } else {
        [zero.checked_sub(n[1])?, n[0], zero]
    };
    Some((origin, ref_dir))
}

/// **A plane's frame with the origin and `+u` direction its author chose.**
///
/// `ref_dir` must lie in the plane and not be zero; it need **not** be a unit vector, and it is
/// reduced to its primitive form here so that two spellings of one direction (`[10,10,0]` and
/// `[20,20,0]`) name **one** frame. That reduction is what makes a frame node interning-friendly
/// — the document's blocker 6.
///
/// `None` on a degenerate plane or direction, when the origin projection is not rational, or when
/// the squared lengths the realization divides by do not fit `i128`.
pub fn plane_frame_named(
    coeffs: [Rat; 4],
    origin: [Rat; 3],
    ref_dir: [Rat; 3],
) -> Option<PlaneFrame> {
    let zero = Rat::from_int(0);
    let n = reduce_direction([coeffs[0], coeffs[1], coeffs[2]])?;
    let u_raw = reduce_direction(ref_dir)?;
    let dot =
        |a: &[Rat; 3]| (0..3).try_fold(zero, |acc, k| acc.checked_add(a[k].checked_mul(a[k])?));
    let (uu, nn) = (dot(&u_raw)?, dot(&n)?);
    // ★★★ **`v_raw` is exact, and taking it exactly is what makes the axes come out exactly.**
    // Realizing `v̂` as `ŵ × û` in f64 costs two roundings that do not cancel: a plain wall whose
    // `v` is exactly `ẑ` came out `0.999999999999999_7`, which is a frame that is not quite
    // orthonormal and an exact path quietly lost. `n ⊥ u_raw` by construction, so
    // `|v_raw|² = |n|²·|u_raw|²` — one inverse square root of an exact rational, and the same
    // wall lands on `1.0`.
    //
    // ★ Its components are bounded by `|n|·|u_raw|`, so the one `checked_mul` below covers them:
    // if the squared length fits, so does every component.
    // ★★★ `v̂` exactly when it fits, and an honest fallback when it does not.
    //
    // `|v_raw|² = |n|²·|u_raw|²` is a product of two squared lengths, so it needs **twice** the
    // width they do — measured, that halves the per-component budget from ~62 bits to ~31 and
    // declines 774 of 780 planes built from a wide normal. When it does not fit, the realization
    // falls back to `v̂ = ŵ × û`, which costs two roundings that do not cancel instead of one.
    // ★ That is the accuracy this crate had before `v_raw` existed — a graceful step down, never
    // a wrong frame. (Filling `vv` with something else would not be a fallback but a corruption:
    // `v̂` would come out the wrong *length* and the basis would not be orthonormal.)
    let cx = |i: usize, j: usize| {
        n[i].checked_mul(u_raw[j])?
            .checked_sub(n[j].checked_mul(u_raw[i])?)
    };
    let v = (|| {
        let vv = nn.checked_mul(uu)?;
        Some(([cx(1, 2)?, cx(2, 0)?, cx(0, 1)?], vv))
    })();
    Some(PlaneFrame {
        origin,
        u_raw,
        n,
        uu,
        nn,
        v,
    })
}

/// A plane's own frame, exactly — what [`plane_frame`] derives.
///
/// The `*_raw` vectors are rational and **not** unit length; `uu`/`nn` are their squared lengths,
/// carried because they are what the realization divides by and because computing them here is
/// what proves the frame fits `i128` before anything tries to use it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PlaneFrame {
    pub origin: [Rat; 3],
    pub u_raw: [Rat; 3],
    pub n: [Rat; 3],
    pub uu: Rat,
    pub nn: Rat,
    /// `(v_raw, |v_raw|²)` when both fit `i128` — `v̂` is then one inverse square root of an exact
    /// rational, and a wall whose `v` is exactly `ẑ` lands on `1.0`. `None` when the product
    /// `|n|²·|u_raw|²` overflows, and the realization takes `v̂ = ŵ × û` instead.
    pub v: Option<([Rat; 3], Rat)>,
}

/// A direction vector divided by its content — the **primitive** integer vector along it, with
/// its sign kept. `None` for the zero vector or on overflow.
///
/// Steps ① and ② of [`canonical_plane_coeffs`] and deliberately not step ③: see [`plane_frame`]
/// for why a direction may not have its sign normalized.
fn reduce_direction(v: [Rat; 3]) -> Option<[Rat; 3]> {
    let mut lcm: i128 = 1;
    for c in v {
        let d = c.denom();
        let g = gcd_u128(lcm.unsigned_abs(), d.unsigned_abs()) as i128;
        lcm = lcm.checked_div(g)?.checked_mul(d)?;
    }
    let mut num = [0i128; 3];
    for (i, c) in v.iter().enumerate() {
        num[i] = c.numer().checked_mul(lcm.checked_div(c.denom())?)?;
    }
    let g = num.iter().fold(0u128, |g, n| gcd_u128(g, n.unsigned_abs()));
    if g == 0 {
        return None; // the zero vector is not a direction
    }
    for n in &mut num {
        *n /= g as i128;
    }
    Some(num.map(Rat::from_int))
}

/// The exact square root of a non-negative rational, or `None` when it is irrational.
///
/// A fraction in lowest terms is a perfect square exactly when its numerator and denominator both
/// are — they share no factor to trade — so this is two integer square roots and two checks.
fn rat_sqrt_exact(v: Rat) -> Option<Rat> {
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
pub fn inv_sqrt_bounded(v: Rat, prec: usize) -> Option<(BigFloat, Mag)> {
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
/// The wide-frame realization (S4) feeds plane data wider than `i128` through this; keeping the
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

/// [`inv_sqrt_bounded`] for a squared length wider than `i128` — the wide-frame (S4) twin.
///
/// Same ladder, same derived bound: the integer enters **exactly** ([`bigint_to_bigfloat`] at
/// its own bit length), the square root rounds once (halving the incoming relative error, which
/// is zero here) and the reciprocal rounds once — under the narrow twin's `4u`, which is kept
/// for symmetry rather than tightened.
///
/// **Uncached** — the wide population is a fraction of a percent of pushes and the memo key
/// would be a `BigInt`; measured before optimizing, per the cache philosophy.
pub fn inv_sqrt_bigint_bounded(v: &num_bigint::BigInt, prec: usize) -> Option<(BigFloat, Mag)> {
    if v.sign() != num_bigint::Sign::Plus {
        return None;
    }
    let ip = prec.max(128);
    let x = bigint_to_bigfloat(v, ip);
    let one = BigFloat::from_i128(1, ip);
    let z = one.div(&x.sqrt(prec, HP_RM), prec, HP_RM);
    let u = Mag::pow2(-(prec as i64));
    let rel = u.times(Mag::of(4.0));
    let mag = match z.exponent() {
        Some(e) if !z.is_zero() => Mag::pow2(e as i64),
        _ => Mag::ZERO,
    };
    Some((z, mag.times(rel)))
}

/// [`inv_sqrt_bounded`] without the memo — the evaluation itself, kept separate so no `INV_SQRT`
/// borrow is held across the arbitrary-precision work.
fn realize_inv_sqrt(v: Rat, prec: usize) -> (BigFloat, Mag) {
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
    let mag = match z.exponent() {
        Some(e) if !z.is_zero() => Mag::pow2(e as i64), // `|x| < 2^exponent`
        _ => Mag::ZERO,
    };
    (z, mag.times(rel))
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
/// `nacre_cip`'s trial precision shares this realization instead of paying for a second one.
///
/// ★ **Correct rounding is what keeps debug and release the same.** A faithfully-rounded value
/// would let the two builds disagree by an ulp, which is the failure this crate already paid for
/// once in the trig path.
///
/// If even 256 bits leave the rounding undecided the midpoint is returned rather than a panic, and
/// the event is counted — see [`INV_SQRT_ESCALATED`].
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
        let (z, rad) = realize_inv_sqrt_memoized(v, prec);
        if let Some(f) = round_to_f64(&z, rad, prec) {
            if i > 0 {
                INV_SQRT_ESCALATED.with_borrow_mut(|(e, _)| *e += 1);
            }
            return f;
        }
    }
    INV_SQRT_ESCALATED.with_borrow_mut(|(_, f)| *f += 1);
    let (z, _) = realize_inv_sqrt_memoized(v, 256);
    to_f64_exact(&z).unwrap_or(f64::NAN)
}

/// [`inv_sqrt_bounded`] for a `v` already known positive.
fn realize_inv_sqrt_memoized(v: Rat, prec: usize) -> (BigFloat, Mag) {
    inv_sqrt_bounded(v, prec).expect("v > 0 checked by the caller")
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
    let (h, rad) = inv_sqrt_bounded(v, P)?;
    let diff = BigFloat::from_f64(f, P).sub(&h, P, HP_RM);
    let mag = if diff.is_zero() {
        0.0
    } else {
        2f64.powi(diff.exponent().unwrap_or(0))
    };
    // Rounded up, for the reason `Angle::measure_realization_error` spells out: when `|diff|` is
    // itself a power of two the octave bound has no slack, and `mag + rad` would round back down
    // to `mag` in f64 — short by the radius.
    Some((mag + rad.exp2().map_or(0.0, |e| 2f64.powi(e as i32))) * (1.0 + 2.0 * f64::EPSILON))
}

/// Greatest common divisor of two magnitudes, Euclid. `gcd(0, 0) == 0`.
fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// **The canonical representative of a rational plane** `a·x + b·y + c·z + d = 0`.
///
/// A plane is scale-invariant — `2x + 4y − 6 = 0` and `x + 2y − 3 = 0` are the same plane — so a
/// family of coefficient vectors describes it and equality has to pick one of them. Three steps:
/// clear the denominators, divide out the content, and fix the sign of the first nonzero component.
/// What comes back is a **primitive integer vector**, and two vectors describing the same plane
/// canonicalize to **bit-identical** arrays. That is what turns *"are these the same plane?"* from a
/// predicate into `==`.
///
/// ★★ **Scale-independence is the whole point, and it is what f64 coefficients cannot give.**
/// `nacre_geom::Plane` stores an un-normalized `raw` normal whose length follows the *face's size*,
/// so two faces of the plane `x = 3` come out as `[2.2, 0, 0, −6.6000000000000005]` and
/// `[13.2, 0, 0, −39.599999999999994]` — the same plane, not exactly proportional, because `d` was a
/// rounded product. Measured: 18 pairs in the census are merged only because a *second* test looks at
/// the faces' coordinates instead (`docs/dev-log.md`).
///
/// ★★★ **The input must be rational by construction, not lifted from those f64 coefficients.**
/// Lifting is lossless but it preserves the drift: `2.2` and `6.6000000000000005` are different
/// dyadics whose exact ratio is not 3, so canonicalizing them still gives two different vectors.
/// Rationals here come from the dimensions the user *wrote* ([`Rat::from_decimal`]) carried through
/// exact arithmetic — the same rule the construction path already follows.
///
/// `None` on `i128` overflow (the denominators' lcm, or a numerator scaled by it), which is the
/// kernel's ordinary demotion signal: the caller keeps the plain rational and the geometric tests
/// answer instead. Nothing is wrong, one shortcut is unavailable.
///
/// All-zero coefficients are not a plane; they canonicalize to themselves.
pub fn canonical_plane_coeffs(coeffs: [Rat; 4]) -> Option<[Rat; 4]> {
    // ① Clear the denominators: multiply through by their lcm.
    let mut lcm: i128 = 1;
    for c in coeffs {
        let d = c.denom();
        let g = gcd_u128(lcm.unsigned_abs(), d.unsigned_abs()) as i128;
        lcm = lcm.checked_div(g)?.checked_mul(d)?;
    }
    let mut num = [0i128; 4];
    for (i, c) in coeffs.iter().enumerate() {
        // Exact: `lcm` is a multiple of every denominator, so the division has no remainder.
        num[i] = c.numer().checked_mul(lcm.checked_div(c.denom())?)?;
    }

    // ② Divide out the content.
    let g = num.iter().fold(0u128, |g, n| gcd_u128(g, n.unsigned_abs()));
    if g == 0 {
        return Some(coeffs); // the zero vector is not a plane
    }
    let g = g as i128;
    for n in &mut num {
        *n /= g;
    }

    // ③ Fix the sign: the first nonzero component is positive.
    if num.iter().find(|n| **n != 0).is_some_and(|n| *n < 0) {
        for n in &mut num {
            *n = n.checked_neg()?;
        }
    }
    Some(num.map(Rat::from_int))
}

/// The nearest f64 to `n / d`, ties to even. Both arguments are strictly positive and
/// came from an `i128`, which is what closes every shift below — so this is a helper
/// for [`Rat::to_f64`] and nothing else.
///
/// The two operands are bounded differently, and the difference is load-bearing. The
/// **denominator** must be `< 2¹²⁷`, because `rem < d` is doubled in the loop; it is,
/// since `Ratio` keeps it positive and `i128::MAX < 2¹²⁷`. The **numerator** may be
/// `2¹²⁷` exactly — `i128::MIN.unsigned_abs()` is — and that is fine: it appears only
/// in shifts whose width is the *other* operand's, so nothing overflows.
///
/// Textbook restoring long division: emit the quotient's leading 54 bits, keep the
/// remainder to tell a tie from a near-tie, then round once. Doing it in integers
/// rather than in a wide float sidesteps double rounding entirely — there is only
/// ever the one rounding, at the end.
fn nearest_f64(n: u128, d: u128) -> f64 {
    debug_assert!(n > 0 && d > 0 && d < (1 << 127));
    let bits = |x: u128| 128 - x.leading_zeros() as i32;

    // The quotient's binary exponent: `2^e ≤ n/d < 2^(e+1)`. The bit-width difference
    // pins it to two candidates, and one comparison picks between them. Both shifts
    // below stay under 2¹²⁸: the shifted operand's width is the *other* one's.
    let t = bits(n) - bits(d);
    let e = if t >= 0 {
        if n >= (d << t) { t } else { t - 1 }
    } else if (n << -t) >= d {
        t
    } else {
        t - 1
    };

    // `m = ⌊(n/d) · 2^(53−e)⌋`, which lies in `[2⁵³, 2⁵⁴)`: 54 bits, one more than an
    // f64 keeps, so the extra bit is the round bit and `rem` is the sticky bit.
    let s = 53 - e;
    let (mut m, mut rem) = (n / d, n % d);
    let sticky;
    if s >= 0 {
        for _ in 0..s {
            m <<= 1;
            // `rem < d < 2¹²⁷`, so this cannot overflow.
            rem <<= 1;
            if rem >= d {
                rem -= d;
                m += 1;
            }
        }
        sticky = rem != 0;
    } else {
        let drop = (-s) as u32;
        sticky = rem != 0 || (m & ((1 << drop) - 1)) != 0;
        m >>= drop;
    }
    debug_assert!((1 << 53..1 << 54).contains(&m));

    // Round to 53 bits, ties to even.
    let (round, mut mant) = (m & 1, m >> 1);
    let mut e = e;
    if round == 1 && (sticky || mant & 1 == 1) {
        mant += 1;
        if mant == 1 << 53 {
            mant >>= 1;
            e += 1;
        }
    }

    // Exact: a 53-bit integer is an exact f64, and scaling by a power of two is exact
    // as long as the result stays normal — which it does, since `|e| ≤ 127` here.
    mant as f64 * (2.0f64).powi(e - 52)
}

/// A *direction* angle in degrees, kept normalized to `[0, 360)` exactly
/// (rational). Accumulation is exact: turning by a rational angle repeatedly and
/// completing a full turn lands back on exactly `0` — no f64 drift. `cos`/`sin`
/// realization crosses into f64 (deg→rad via π): the irrational-realization
/// boundary, past which a *judgement* realizes at arbitrary precision instead
/// ([`Angle::cos_sin_at`]). The angle stays exact; only its realized coordinate
/// carries tol.
///
/// This is the direction type. A multi-turn *amount* (helix pitch × turns, revolve
/// sweep > 360°) must preserve the turn count, so it belongs to a separate
/// *unnormalized* `Sweep` type — flat 2D sketches never need it, so it is left as
/// a documented companion, not built here.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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

    /// `(cos, sin)` realized in arbitrary precision at `prec` bits — the judgment
    /// path. astro-float replaces twofloat here (H1.5: twofloat's trig was
    /// f64-level near zero-crossings). The depth is the caller's: a judgement's precision is a
    /// property of the model it judges, not of this crate. (numer/denom pass through f64, exact
    /// for the small values used here; a general large-rational path would build from a string.)
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
    /// **Memoized by `(self, prec)`** in [`TRIG`] — the value below is a pure function of those
    /// two, and a rotation's realization asks for the same angle once per point. See [`TRIG`] for
    /// why a process-wide memo is sound here and a handle-keyed one would not be.
    ///
    /// Returns `(cos, sin, |Δcos|, |Δsin|)`.
    pub fn cos_sin_bounded(self, prec: usize) -> (BigFloat, BigFloat, Mag, Mag) {
        if let Some(hit) = TRIG.with_borrow(|t| t.get(&(self, prec)).cloned()) {
            return hit;
        }
        let out = self.realize_cos_sin(prec);
        TRIG.with_borrow_mut(|t| t.insert((self, prec), out.clone()));
        out
    }

    /// [`cos_sin_bounded`](Self::cos_sin_bounded) without the memo — the evaluation itself.
    ///
    /// **Separate so no `TRIG` borrow is held across it.** `HP_CONSTS` is borrowed for the whole
    /// realization and the trig calls are the slow part; nesting the memo's borrow around that is
    /// how a re-entrant call would panic rather than merely be slow.
    fn realize_cos_sin(self, prec: usize) -> (BigFloat, BigFloat, Mag, Mag) {
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
            let u = Mag::pow2(-(prec as i64));
            // Relative error of the argument: the two `i128 → f64` conversions, then four
            // rounded high-precision operations (the division, the product, the division, and π).
            // Four rounded operations build the argument (the division, the product, the
            // division, and π itself). The integers contribute nothing — they go in exactly.
            let rel = u.times(Mag::of(4.0));
            // `|θ|` in radians, over-estimated from its exponent (`|x| < 2^exponent`).
            let theta = match rad.exponent() {
                Some(e) if !rad.is_zero() => Mag::pow2(e as i64),
                _ => Mag::ZERO,
            };
            let d_theta = theta.times(rel);
            let (c, s) = (rad.cos(prec, HP_RM, cc), rad.sin(prec, HP_RM, cc));
            // `|x| < 2^exponent` — the slope of the *other* function, and the scale of the
            // half-ulp of this one.
            let ub = |x: &BigFloat| match x.exponent() {
                Some(e) if !x.is_zero() => Mag::pow2(e as i64),
                _ => Mag::ZERO,
            };
            let (uc, us) = (ub(&c), ub(&s));
            let err_cos = us.times(d_theta).plus(uc.times(u));
            let err_sin = uc.times(d_theta).plus(us.times(u));
            (c, s, err_cos, err_sin)
        })
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

    /// `(cos, sin)` realized in f64 — **exact** (`0.0`/`±1.0`) for the 90°-family, and everywhere
    /// else the arbitrary-precision value **correctly rounded**. The single source of truth for
    /// realizing a rotation angle into f64: every path that turns a point or direction by an angle
    /// goes through here.
    ///
    /// ★★★ **No libm.** This used to be `(deg.to_f64() * PI / 180.0).cos()`, and the error of that
    /// is not something anyone contracts: neither Rust nor any platform promises an accuracy for
    /// `f64::cos`, measured here at `< 5ε`. Worse, it was not a function of the angle — `sin 27°`
    /// came out one ulp apart between a debug and a release build, and between two call sites
    /// *within one release build*, because LLVM evaluates a visible constant angle at compile time
    /// and its answer differs from the runtime library's. Rounding the high-precision realization
    /// instead makes the result **unique**: the same bits on every platform, profile and call site.
    ///
    /// ★★ **The 90°-family branch is not an optimization, it is what makes this terminate.**
    /// `cos 90°` is exactly `0`, so its interval straddles zero and the two ends never round to the
    /// same f64 no matter how deep the realization goes — [`round_to_f64`] would answer `None`
    /// forever. Niven's theorem says the only rational values are `{0, ±1/2, ±1}`, and only the
    /// zeros have this problem; they are exactly the quadrantal ones caught here. (`±1/2` at 60°
    /// and friends is exactly representable and comes out of the general path just fine.)
    ///
    /// **The realization is memoized per angle** ([`F64_TRIG`]), because one costs ~32µs against
    /// libm's 22ns and `Isometry::apply_point` runs per vertex. The memo cannot change an answer —
    /// a correctly rounded value is unique — so it is a cost question only.
    ///
    /// A caller accounting for the realization's *error* hands what it got back to
    /// [`realization_error_of`](Self::realization_error_of), which returns zero for exactly this
    /// family — so nothing here has to report which branch ran.
    pub fn cos_sin_f64(self) -> (f64, f64) {
        if let Some((cr, sr)) = self.try_exact_cos_sin() {
            return (cr.to_f64(), sr.to_f64());
        }
        if let Some(hit) = F64_TRIG.with_borrow(|m| m.get(&self).copied()) {
            return hit;
        }
        // Outside the borrow: the realization below takes `TRIG`'s and `HP_CONSTS`' in turn.
        let out = self.realize_rounded_f64();
        F64_TRIG.with_borrow_mut(|m| m.insert(self, out));
        out
    }

    /// [`cos_sin_f64`](Self::cos_sin_f64)'s general branch without the memo.
    ///
    /// **The ladder is `TRIAL_PREC` then twice that, and both rungs are derived rather than tried.**
    /// The realization's error is dominated by the degrees→radians conversion, not by the cosine:
    /// `≈ 30 · 2⁻ᵖʳᵉᶜ`. Against a result of magnitude `2^e` that has to clear a half-ulp of `2^(e-54)`,
    /// so `prec > 54 - e + 5`.
    ///
    /// - **128** covers every `|cos| > 2⁻⁶⁹`, which is every angle a model has ever held. It is also
    ///   `nacre_cip`'s trial precision, so a model that goes on to be judged **shares this exact
    ///   realization** rather than paying for a second one at a different depth. That sharing is
    ///   why the ladder does not start lower: 64 bits would satisfy the inequality and measured no
    ///   cheaper (31.9µs against 32.3µs — the cost is setup, not bit count), but it would be a
    ///   different `TRIG` key and so pure duplication for anything judged.
    /// - **256** is the proven cap. An `Angle` holds `Ratio<i128>`, so a normalized angle cannot
    ///   come closer to 90° than `1/denominator ≥ 5.9e-39` degrees; `|cos|` is therefore never
    ///   below `~2⁻¹³³`, which needs `prec > 192`.
    ///
    /// ★ **If even 256 leaves it undecided the answer is still returned, not a panic.** The value is
    /// then *faithfully* rounded (within an ulp) instead of correctly rounded, which stays sound
    /// because [`realization_error_of`](Self::realization_error_of) measures the error that is
    /// actually there and the tolerance grows to match — and it stays deterministic, because a
    /// 256-bit midpoint is. It is counted so that "can't happen" does not quietly become "happens".
    fn realize_rounded_f64(self) -> (f64, f64) {
        for (i, prec) in [128usize, 256].into_iter().enumerate() {
            let (c, s, rc, rs) = self.cos_sin_bounded(prec);
            if let (Some(cf), Some(sf)) = (round_to_f64(&c, rc, prec), round_to_f64(&s, rs, prec)) {
                if i > 0 {
                    ROUND_ESCALATED.with_borrow_mut(|(e, _)| *e += 1);
                }
                return (cf, sf);
            }
        }
        ROUND_ESCALATED.with_borrow_mut(|(_, f)| *f += 1);
        let (c, s, _, _) = self.cos_sin_bounded(256);
        (
            to_f64_exact(&c).unwrap_or(f64::NAN),
            to_f64_exact(&s).unwrap_or(f64::NAN),
        )
    }

    /// **How far the `(cos, sin)` the caller was handed sits from the true ones** — measured
    /// against an arbitrary-precision realization, not assumed from a constant.
    ///
    /// `|f64 − true| ≤ |f64 − hp midpoint| + hp's own radius`, which is the ruler this kernel
    /// already uses for a rational's realization in `WitnessPoint`'s `translate`, `mirror` and pivot terms.
    /// Memoized in [`F64_ERR`], keyed by the angle *and the pair* — see there for why the pair.
    ///
    /// ★★★ **The caller passes the values in rather than letting this re-realize them**, and that
    /// is the whole soundness argument. The consumer is `WitnessPoint::rotate_about`, whose `tol` must bound
    /// the error in the `coord` it just wrote from *its* `cos_sin_f64()` result. An error measured
    /// against a second, independent realization would bound a number nobody stored — and those two
    /// realizations are measured to differ (see [`F64_ERR`]). Taking `c` and `s` as arguments makes
    /// "the error describes the value that was used" hold by construction instead of by hope.
    ///
    /// ★★ **This is what lets the error accounting stop guessing.** `f64::cos` has no accuracy
    /// contract — neither Rust nor any libm promises one — so the bound above it used to be a
    /// measured-once constant with margin, sound only on platforms like the one it was taken on.
    /// A kernel that ships to browsers cannot know that. Measuring instead means a worse libm
    /// simply reports a larger error and the tolerance grows to match: **the kernel adapts rather
    /// than assumes.**
    ///
    /// ★ **Exactly zero for the 90°-family**, because `cos_sin_f64` returns `0.0`/`±1.0` there and
    /// those *are* the true values. Callers rely on that zero: a rotation whose realization carries
    /// no error also performs no rounding downstream (`u·(±1)` and `u·0` are exact), which is what
    /// keeps a quadrantal origin rotation at tol 0.
    ///
    /// ★ **An `f64`, not a [`Mag`].** `Mag` exists because a deep ladder's `2⁻ᵖʳᵉᶜ` underflows
    /// `f64` to zero and a zero radius claims exactness; this quantity is always ε-scale, so that
    /// hazard is absent — and the consumer is `WitnessPoint::tol`, which is `f64`.
    ///
    /// The reading is at octave granularity (`bf_mag` is `2^exponent`), so it can sit up to 2×
    /// above the true error. Conservative in the sound direction, and still a measurement.
    pub fn realization_error_of(self, c: f64, s: f64) -> (f64, f64) {
        let key = (self, c.to_bits(), s.to_bits());
        if let Some(hit) = F64_ERR.with_borrow(|m| m.get(&key).copied()) {
            return hit;
        }
        // Outside the borrow: the computation below takes `TRIG`'s and `HP_CONSTS`' in turn.
        let out = self.measure_realization_error(c, s);
        F64_ERR.with_borrow_mut(|m| m.insert(key, out));
        out
    }

    /// [`realization_error_of`](Self::realization_error_of) without the memo — the measurement.
    fn measure_realization_error(self, c: f64, s: f64) -> (f64, f64) {
        const P: usize = 128; // the hp radius is then ~2⁻¹²⁸ against an ε-scale quantity
        // ★ The zero is claimed of *these* values, not of the angle: the exactness that callers
        // depend on is "the pair I am holding is the true cos/sin", and only comparing the pair
        // says that. An angle in the family whose caller somehow realized it the general way is
        // then measured like any other rather than being handed a zero it has not earned.
        if self
            .try_exact_cos_sin()
            .is_some_and(|(cr, sr)| cr.to_f64() == c && sr.to_f64() == s)
        {
            return (0.0, 0.0);
        }
        let (hc, hs, rc, rs) = self.cos_sin_bounded(P);
        let gap = |f: f64, h: &BigFloat, rad: Mag| {
            let diff = BigFloat::from_f64(f, P).sub(h, P, HP_RM);
            let mag = if diff.is_zero() {
                0.0
            } else {
                2f64.powi(diff.exponent().unwrap_or(0))
            };
            // ★ **The sum is rounded *up*, and that is not pedantry.** `mag` is an octave bound, so
            // it usually sits well above the truth — but when `|diff|` is itself a power of two the
            // slack is exactly zero, and then `mag + rad` rounds back down to `mag` in f64 and the
            // "bound" is short by the radius. `sin 30°` is that case: it misses 0.5 by exactly
            // `2⁻⁵⁴`, and the ground-truth test caught the missing ulp the day this was written.
            (mag + rad.exp2().map_or(0.0, |e| 2f64.powi(e as i32))) * (1.0 + 2.0 * f64::EPSILON)
        };
        (gap(c, &hc, rc), gap(s, &hs, rs))
    }
}

/// A coordinate axis — the fixed axis of an axis-aligned rotation (overhaul stage
/// 1b restricts to `X`/`Y`/`Z`, the form `exact3d` validated; arbitrary rational
/// axes via Rodrigues are a later extension).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    /// The coordinate index this axis names (`X → 0`, `Y → 1`, `Z → 2`) — the one a reflection in
    /// a plane perpendicular to it negates.
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }

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

/// Whether **any** rotation about an axis-parallel line fixes the plane
/// `a·x + b·y + c·z + d = 0` as a set: the normal rides the axis (its other two
/// components exactly zero), so the normal coordinate of every point is untouched
/// whatever the pivot or angle. One of the three set-invariance atoms — the others
/// are [`translation_fixes_plane`] and [`mirror_fixes_plane`]; `Isometry::fixes_plane`
/// composes the first two, and a recorded motion chain is checked node by node.
pub fn axis_rotation_fixes_plane(axis: Axis, coeffs: &[Rat; 4]) -> bool {
    let i = axis.index();
    let zero = Rat::from_int(0);
    coeffs[i] != zero && (0..3).all(|k| k == i || coeffs[k] == zero)
}

/// Whether the translation `offset` fixes the plane as a set: its component along the
/// normal is exactly zero (`n · offset == 0`, checked — an overflow answers `false`,
/// a conservative miss).
pub fn translation_fixes_plane(offset: &[Rat; 3], coeffs: &[Rat; 4]) -> bool {
    let zero = Rat::from_int(0);
    let dot = coeffs[..3]
        .iter()
        .zip(offset)
        .try_fold(zero, |acc, (&a, &t)| a.checked_mul(t)?.checked_add(acc));
    dot == Some(zero)
}

/// Whether the reflection in the coordinate plane `axis = offset` fixes the plane as a
/// set: either the normal is perpendicular to the mirror axis (the plane contains the
/// mirrored direction, so the set maps to itself), or the plane **is** the mirror plane
/// (`n ∥ axis` and `n_axis · offset + d == 0`, checked — overflow answers `false`).
/// Set-fixed only: the second case flips the normal's sense, which a canonical name
/// does not carry.
pub fn mirror_fixes_plane(axis: Axis, offset: Rat, coeffs: &[Rat; 4]) -> bool {
    let i = axis.index();
    let zero = Rat::from_int(0);
    if coeffs[i] == zero {
        return true;
    }
    (0..3).all(|k| k == i || coeffs[k] == zero)
        && coeffs[i]
            .checked_mul(offset)
            .and_then(|p| p.checked_add(coeffs[3]))
            == Some(zero)
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
/// `[f64; 3]`, so `nacre-scalar` never depends on `nacre-math`; the caller
/// (`nacre-ops`) applies it to `Point3`/`Plane`.
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

    /// Whether this isometry maps the plane `a·x + b·y + c·z + d = 0` onto itself
    /// **as a set** — exactly, on the rational definition.
    ///
    /// A rotation about an axis parallel to the plane's normal permutes the plane
    /// within itself whatever the pivot (the normal coordinate is untouched), and a
    /// translation moves the plane iff its component along the normal is nonzero.
    /// So the condition is: the rotation, if any, has its axis parallel to
    /// `(a, b, c)`, and `(a, b, c) · translate == 0`. `d` plays no part — every
    /// plane sharing the normal answers alike. Checked arithmetic; an overflow
    /// answers `false`, a conservative miss (the plane is then carried on the
    /// recorded path, slower but never wrong).
    ///
    /// A reflection is not an `Isometry`; [`mirror_fixes_plane`] answers for the mirror
    /// motion, and a recorded chain is checked node by node with the same three atoms.
    pub fn fixes_plane(&self, coeffs: &[Rat; 4]) -> bool {
        self.rotate
            .is_none_or(|r| axis_rotation_fixes_plane(r.axis, coeffs))
            && translation_fixes_plane(&self.translate, coeffs)
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

    /// **This isometry applied to a rational plane**, exactly — the definition-level twin of
    /// [`apply_point`](Isometry::apply_point).
    ///
    /// A plane is not a bag of points, so it does not go through `apply_point`: under
    /// `x ↦ R(x − p) + p + t` the plane `n·x + d = 0` becomes
    ///
    /// ```text
    /// n' = R·n            d' = d + n·p − n'·(p + t)
    /// ```
    ///
    /// (For a pure translation that collapses to the familiar `d' = d − n·t`.)
    ///
    /// `None` unless the rotation is one the rationals can state — the 90°-family, where
    /// [`Angle::try_exact_cos_sin`] gives `cos`/`sin` in `{0, ±1}` — or on `i128` overflow.
    /// That is the same condition [`Isometry::is_exact`] reports, so a caller that already
    /// checked it will not be surprised here.
    ///
    /// The result is canonicalized ([`canonical_plane_coeffs`]), so a plane reached by two
    /// different routes lands on the **same array**.
    pub fn plane_coeffs(&self, c: [Rat; 4]) -> Option<[Rat; 4]> {
        let n = [c[0], c[1], c[2]];
        let (n2, pivot) = match self.rotate {
            None => (n, [Rat::from_int(0); 3]),
            Some(r) => {
                let (cos, sin) = r.angle.try_exact_cos_sin()?;
                let (i, j) = r.axis.plane();
                let mut m = n;
                m[i] = n[i].checked_mul(cos)?.checked_sub(n[j].checked_mul(sin)?)?;
                m[j] = n[i].checked_mul(sin)?.checked_add(n[j].checked_mul(cos)?)?;
                (m, r.point)
            }
        };
        // d' = d + n·p − n'·(p + t)
        let mut d = c[3];
        for k in 0..3 {
            d = d.checked_add(n[k].checked_mul(pivot[k])?)?;
            let q = pivot[k].checked_add(self.translate[k])?;
            d = d.checked_sub(n2[k].checked_mul(q)?)?;
        }
        canonical_plane_coeffs([n2[0], n2[1], n2[2], d])
    }

    /// **This isometry applied to a rational point**, exactly — the rational twin
    /// [`apply_point`](Isometry::apply_point) never had, and the map
    /// [`plane_coeffs`](Isometry::plane_coeffs) is described against.
    ///
    /// `x ↦ R(x − p) + p + t`, with `R` read from [`Angle::try_exact_cos_sin`]. `None` under the
    /// same conditions as `plane_coeffs` — a rotation the rationals cannot state, or `i128`
    /// overflow — so a caller that can move one description exactly can move the other.
    ///
    /// ★ **Why a plane needs it even though a plane is not a bag of points.** A plane's *truth* is
    /// three points on it (`Model::surface_points`): its coefficients are a product of two point
    /// differences and overflow `i128` far sooner than the points do. Moving such a plane means
    /// moving its points, and carrying them through f64 would put a rounded coordinate back into a
    /// definition — the thing storing points was meant to stop.
    pub fn point_rat(&self, x: [Rat; 3]) -> Option<[Rat; 3]> {
        let turned = match self.rotate {
            None => x,
            Some(r) => {
                let (cos, sin) = r.angle.try_exact_cos_sin()?;
                let (i, j) = r.axis.plane();
                // Pivot-relative, exactly as the realization does it: `u = x − p`, turn, shift back.
                let u = x[i].checked_sub(r.point[i])?;
                let v = x[j].checked_sub(r.point[j])?;
                let mut m = x;
                m[i] = r.point[i]
                    .checked_add(u.checked_mul(cos)?.checked_sub(v.checked_mul(sin)?)?)?;
                m[j] = r.point[j]
                    .checked_add(u.checked_mul(sin)?.checked_add(v.checked_mul(cos)?)?)?;
                m
            }
        };
        Some([
            turned[0].checked_add(self.translate[0])?,
            turned[1].checked_add(self.translate[1])?,
            turned[2].checked_add(self.translate[2])?,
        ])
    }
}

/// **A reflection in `axis = offset` applied to a rational point**, exactly — `x_a ↦ 2·offset − x_a`.
///
/// The point twin of [`mirror_plane_coeffs`], and the same one-line map: a reflection is its own
/// inverse and touches one coordinate. `None` only on `i128` overflow.
pub fn mirror_point_rat(x: [Rat; 3], axis: Axis, offset: Rat) -> Option<[Rat; 3]> {
    let a = axis.index();
    let mut out = x;
    out[a] = offset.checked_mul(Rat::from_int(2))?.checked_sub(x[a])?;
    Some(out)
}

/// **A reflection in `axis = offset` applied to a rational plane**, exactly.
///
/// `x_a ↦ 2·offset − x_a` negates that component of the normal and shifts the offset:
/// `n'_a = −n_a`, `d' = d + 2·offset·n_a`. The reflection is its own inverse, which is why the
/// map and its transpose-inverse coincide and no case analysis is needed.
///
/// ★ **The determinant is `−1`.** A caller that relies on orientation being preserved has to
/// account for that itself; this function states where the plane goes, nothing more.
pub fn mirror_plane_coeffs(c: [Rat; 4], axis: Axis, offset: Rat) -> Option<[Rat; 4]> {
    let a = axis.index();
    let mut out = c;
    out[a] = Rat::from_int(0).checked_sub(c[a])?;
    out[3] = c[3].checked_add(offset.checked_mul(Rat::from_int(2))?.checked_mul(c[a])?)?;
    canonical_plane_coeffs(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The depth these tests realize at. **Test-local on purpose**: production has no default
    /// precision — it is computed from the model — so a constant here must not be reachable from
    /// outside the fixtures that chose it.
    const GT_PREC: usize = 160;

    fn ints(v: [i128; 4]) -> [Rat; 4] {
        v.map(Rat::from_int)
    }

    /// [`Isometry::fixes_plane`] — the set-invariance predicate the invariant-plane
    /// restatement gates on. Hand-known planes and motions; the overflow case locks the
    /// conservative direction (a miss records a node — slower, never wrong).
    mod fixes_plane {
        use super::*;

        fn rot(axis: Axis, deg: i128, point: [i128; 3]) -> Isometry {
            Isometry::rotation(Rotation {
                axis,
                point: point.map(Rat::from_int),
                angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
            })
        }

        #[test]
        fn an_axis_parallel_normal_rides_any_pivot() {
            let z_plane = ints([0, 0, 1, -3]);
            assert!(rot(Axis::Z, 30, [0, 0, 0]).fixes_plane(&z_plane));
            assert!(rot(Axis::Z, 30, [7, -2, 5]).fixes_plane(&z_plane));
            let x_wall = ints([1, 0, 0, -2]);
            assert!(!rot(Axis::Z, 30, [0, 0, 0]).fixes_plane(&x_wall));
            assert!(rot(Axis::X, 30, [1, 1, 1]).fixes_plane(&x_wall));
            assert!(!rot(Axis::X, 30, [1, 1, 1]).fixes_plane(&z_plane));
        }

        #[test]
        fn a_translation_moves_the_plane_iff_it_leaves_it() {
            let z_plane = ints([0, 0, 1, -3]);
            let mut in_plane = rot(Axis::Z, 30, [0, 0, 0]);
            in_plane.translate = [Rat::from_int(1), Rat::from_int(2), Rat::from_int(0)];
            assert!(in_plane.fixes_plane(&z_plane));
            let mut off_plane = rot(Axis::Z, 30, [0, 0, 0]);
            off_plane.translate = [Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)];
            assert!(!off_plane.fixes_plane(&z_plane));
        }

        #[test]
        fn a_pure_translation_needs_no_rotation_clause() {
            // The population the whole-solid probe loses to magnitude rounding: a huge
            // in-plane offset still fixes every plane it slides along.
            let big = Rat::from_int(1i128 << 100);
            let iso = Isometry::translation([big, big, Rat::from_int(0)]);
            assert!(iso.fixes_plane(&ints([0, 0, 1, -3])));
            assert!(!iso.fixes_plane(&ints([1, 0, 0, -2])));
        }

        #[test]
        fn an_overflowing_dot_is_a_conservative_miss() {
            // a·tx and −b·ty would cancel to a true zero, but each product overflows
            // i128 — the checked arithmetic answers false rather than guessing.
            let wide = 1i128 << 100;
            let plane = [
                Rat::from_int(wide),
                Rat::from_int(wide),
                Rat::from_int(0),
                Rat::from_int(-1),
            ];
            let iso = Isometry::translation([
                Rat::from_int(wide),
                Rat::from_int(-wide),
                Rat::from_int(0),
            ]);
            assert!(!iso.fixes_plane(&plane));
        }
    }

    /// [`plane_name_from_meets`] — the name of a plane through meets of any width (open
    /// item 17), locked on both output widths with hand-known planes so the derivation is
    /// never its own oracle.
    mod name_from_meets {
        use super::*;
        use num_bigint::BigInt;
        use proptest::prelude::*;

        fn r(n: i128, d: i128) -> Rat {
            Rat::new(n, d).unwrap()
        }

        /// The carriers of one fixture vertex: `x = m·2⁻⁸⁰`, `y = w·5⁻⁴⁰`, and a third plane
        /// `base + λ·(x − m·2⁻⁸⁰)` — which meets the first two exactly where `base` does, so
        /// the meet lies on `base` by construction while no carrier *is* `base`.
        fn meet_on(base: [Rat; 4], m: i128, w: i128, lambda: i128) -> MeetPoint {
            let a = PlaneName::Narrow([r(1 << 80, 1), r(0, 1), r(0, 1), r(-m, 1)]);
            let b = PlaneName::Narrow([r(0, 1), r(5i128.pow(40), 1), r(0, 1), r(-w, 1)]);
            let u = r(m, 1 << 80);
            let c = PlaneName::Narrow([
                base[0].checked_add(r(lambda, 1)).unwrap(),
                base[1],
                base[2],
                base[3]
                    .checked_sub(r(lambda, 1).checked_mul(u).unwrap())
                    .unwrap(),
            ]);
            three_planes_big([&a, &b, &c]).expect("the fixture carriers meet")
        }

        /// ★ Narrow-name output: three **wide** meets on the hand-known plane
        /// `T: 7x + 11y − 13z + 1 = 0` (already canonical — gcd 1, first coefficient
        /// positive) name exactly `T`.
        #[test]
        fn wide_meets_on_a_narrow_plane_name_it() {
            let t = [r(7, 1), r(11, 1), r(-13, 1), r(1, 1)];
            let meets = [
                meet_on(t, 1, 1, 1),
                meet_on(t, 3, 1, 2),
                meet_on(t, 1, 7, 3),
            ];
            for m in &meets {
                assert!(
                    matches!(m, MeetPoint::Wide(_)),
                    "fixture validity: every meet must be wide, or this test measures nothing"
                );
            }
            assert_eq!(
                plane_name_from_meets([&meets[0], &meets[1], &meets[2]]),
                Some(PlaneName::Narrow([r(7, 1), r(11, 1), r(-13, 1), r(1, 1)])),
                "the plane the meets lie on, by hand"
            );
        }

        /// ★ Wide-name output — the other side of the derivation's fork. `T_w: 2⁷⁰·x + y/5³⁰
        /// − z + 1 = 0` has all-`Rat` coefficients but its canonical integers are `×5³⁰`:
        /// `(2⁷⁰·5³⁰, 1, −5³⁰, 5³⁰)`, whose largest is ~2¹⁴⁰ — a `Wide` name, still known by
        /// hand.
        #[test]
        fn wide_meets_on_a_wide_plane_name_it() {
            let tw = [r(1 << 70, 1), r(1, 5i128.pow(30)), r(-1, 1), r(1, 1)];
            let meets = [
                meet_on(tw, 1, 1, 1),
                meet_on(tw, 3, 1, 2),
                meet_on(tw, 1, 7, 3),
            ];
            for m in &meets {
                assert!(matches!(m, MeetPoint::Wide(_)), "fixture validity");
            }
            let five30 = BigInt::from(5i128.pow(30));
            let expected = [
                BigInt::from(1i128 << 70) * &five30,
                BigInt::from(1),
                -&five30,
                five30.clone(),
            ];
            assert_eq!(
                plane_name_from_meets([&meets[0], &meets[1], &meets[2]]),
                Some(PlaneName::Wide(expected)),
                "the canonical integers of T_w, by hand"
            );
        }

        /// Collinear meets name nothing — width is never the cause of a `None` here.
        #[test]
        fn collinear_meets_name_no_plane() {
            let p = |x: i128| MeetPoint::Narrow([r(x, 1), r(x, 1), r(0, 1)]);
            let (a, b, c) = (p(0), p(1), p(2));
            assert_eq!(plane_name_from_meets([&a, &b, &c]), None);
        }

        proptest! {
            /// All-narrow meets: the new derivation answers exactly what [`plane_name_exact`]
            /// answers — the independent oracle for the lift being a pure generalization.
            #[test]
            fn narrow_meets_agree_with_the_point_derivation(
                pts in proptest::array::uniform3(proptest::array::uniform3((-40i128..=40, 1u32..=40))),
            ) {
                let p: [[Rat; 3]; 3] = pts.map(|q| q.map(|(n, e)| r(n, 1i128 << e)));
                let meets = p.map(MeetPoint::Narrow);
                prop_assert_eq!(
                    plane_name_from_meets([&meets[0], &meets[1], &meets[2]]),
                    plane_name_exact(p[0], p[1], p[2])
                );
            }
        }
    }

    /// The integer sign predicates, differentially locked.
    ///
    /// Three locks, matching the three ways they can be wrong:
    /// - **against the f64 twins** on narrow input — same convention, two vessels;
    /// - **scale-invariance into wide width** — the sign predicates are row-linear, so a row
    ///   times a 200-bit positive integer is a genuine `> i128` input whose true answer the
    ///   narrow twin still knows. This is the only exact ground truth for wide inputs that does
    ///   not test `f(x)` against `f(x)`;
    /// - **negation as the direction negative-control** — exactly the direction-sensitive
    ///   slots flip (`plane_side`: `j` only; `cmp_coord`: none; `dir_sign`: every row), which
    ///   is the per-slot σ analysis the judging layer's rescue relies on.
    mod int_predicates {
        use super::*;
        use nacre_predicates::{ThreePlane, det3_sign, indirect_cmp_coord, indirect_plane_side};
        use num_bigint::BigInt;
        use proptest::prelude::*;

        /// Small integer coefficients — exact in `f64`, so the twins speak the same input.
        fn coeffs() -> impl Strategy<Value = [i32; 4]> {
            proptest::array::uniform4(-9i32..=9)
        }

        fn big(c: [i32; 4]) -> [BigInt; 4] {
            c.map(BigInt::from)
        }
        fn as_f64(c: [i32; 4]) -> [f64; 4] {
            c.map(f64::from)
        }
        fn normals(p: [[i32; 4]; 3]) -> [[f64; 3]; 3] {
            p.map(|c| [f64::from(c[0]), f64::from(c[1]), f64::from(c[2])])
        }
        fn rows(p: &[[BigInt; 4]; 3]) -> [&[BigInt; 4]; 3] {
            [&p[0], &p[1], &p[2]]
        }

        proptest! {
            #[test]
            fn the_integer_predicates_agree_with_the_f64_twins(
                p in proptest::array::uniform3(coeffs()),
                j in coeffs(),
            ) {
                let bp = p.map(big);
                prop_assert_eq!(
                    int_plane_side(rows(&bp), &big(j)),
                    indirect_plane_side(&ThreePlane(p.map(as_f64)), as_f64(j))
                );
                prop_assert_eq!(int_dir_sign(rows(&bp)), det3_sign(normals(p)));
            }

            #[test]
            fn cmp_coord_agrees_with_the_f64_twin(
                a in proptest::array::uniform3(coeffs()),
                b in proptest::array::uniform3(coeffs()),
                axis in 0usize..3,
            ) {
                let (ba, bb) = (a.map(big), b.map(big));
                // Both twins' precondition: each triple meets in a point.
                prop_assume!(int_dir_sign(rows(&ba)) != 0 && int_dir_sign(rows(&bb)) != 0);
                prop_assert_eq!(
                    int_cmp_coord(rows(&ba), rows(&bb), axis),
                    indirect_cmp_coord(&ThreePlane(a.map(as_f64)), &ThreePlane(b.map(as_f64)), axis)
                );
            }

            /// Rows scaled by ~200-bit positive integers — genuinely wider than `i128`, with
            /// the unscaled f64 twin as exact ground truth.
            #[test]
            fn a_wide_positive_scale_cannot_move_any_answer(
                p in proptest::array::uniform3(coeffs()),
                j in coeffs(),
                b in proptest::array::uniform3(coeffs()),
                axis in 0usize..3,
            ) {
                // A distinct scale per row, so no cross-row cancellation can hide a bug.
                let scale = |c: [i32; 4], k: u32| -> [BigInt; 4] {
                    let s = (BigInt::from(1) << 199) + BigInt::from(12345 + k);
                    c.map(|x| BigInt::from(x) * &s)
                };
                let sp: [[BigInt; 4]; 3] = core::array::from_fn(|k| scale(p[k], k as u32));
                let sj = scale(j, 7);
                prop_assert_eq!(
                    int_plane_side(rows(&sp), &sj),
                    indirect_plane_side(&ThreePlane(p.map(as_f64)), as_f64(j))
                );
                prop_assert_eq!(int_dir_sign(rows(&sp)), det3_sign(normals(p)));
                let sb: [[BigInt; 4]; 3] = core::array::from_fn(|k| scale(b[k], 11 + k as u32));
                if int_dir_sign(rows(&sp)) != 0 && int_dir_sign(rows(&sb)) != 0 {
                    prop_assert_eq!(
                        int_cmp_coord(rows(&sp), rows(&sb), axis),
                        indirect_cmp_coord(
                            &ThreePlane(p.map(as_f64)),
                            &ThreePlane(b.map(as_f64)),
                            axis
                        )
                    );
                }
            }

            /// The direction negative-control: negating one row flips exactly the
            /// direction-sensitive answers and nothing else.
            #[test]
            fn a_negated_row_flips_exactly_the_direction_sensitive_answers(
                p in proptest::array::uniform3(coeffs()),
                j in coeffs(),
                k in 0usize..3,
                axis in 0usize..3,
            ) {
                let bp = p.map(big);
                let bj = big(j);
                let neg = |c: &[BigInt; 4]| -> [BigInt; 4] { core::array::from_fn(|i| -&c[i]) };
                let mut np = bp.clone();
                np[k] = neg(&bp[k]);
                let nj = neg(&bj);

                // `plane_side`: a negated `p` row negates `D` and the dot together — invariant;
                // a negated `j` negates the dot alone — flipped.
                let side = int_plane_side(rows(&bp), &bj);
                prop_assert_eq!(int_plane_side(rows(&np), &bj), side);
                prop_assert_eq!(int_plane_side(rows(&bp), &nj), -side);

                // `dir_sign`: one negated row negates the determinant.
                prop_assert_eq!(int_dir_sign(rows(&np)), -int_dir_sign(rows(&bp)));

                // `cmp_coord`: orientation-invariant in every row.
                if int_dir_sign(rows(&bp)) != 0 && int_dir_sign(rows(&bj_triple(&bp, &bj))) != 0 {
                    let other = bj_triple(&bp, &bj);
                    prop_assert_eq!(
                        int_cmp_coord(rows(&np), rows(&other), axis),
                        int_cmp_coord(rows(&bp), rows(&other), axis)
                    );
                }
            }
        }

        /// A second triple for `cmp_coord`, made from the first by swapping in `j` — cheap, and
        /// dependent enough on the inputs to exercise real coordinates.
        fn bj_triple(p: &[[BigInt; 4]; 3], j: &[BigInt; 4]) -> [[BigInt; 4]; 3] {
            [p[1].clone(), p[2].clone(), j.clone()]
        }
    }

    /// The projection lands **exactly on** the plane it came from, tilted ones included — checked
    /// as `n·p + d == 0` in rationals, not within a tolerance.
    /// A plane pushed along its own normal lands where the distance says — checked by taking the
    /// projected origin of the moved plane and measuring it against the original's.
    #[test]
    fn a_plane_pushed_along_its_normal_moves_by_that_distance() {
        for (raw, name) in [
            ([0, 0, 1, -2], "z = 2"),
            ([0, 0, 10, -21], "z = 2.1"),
            ([3, 0, -4, -5], "3-4-5 tilt"),
            ([5, 12, 0, -13], "5-12-13"),
        ] {
            let c = ints(raw);
            let t = Rat::from_decimal(7.7).unwrap();
            let moved =
                plane_offset(c, t).unwrap_or_else(|| panic!("{name} has a rational normal"));
            // The two projected origins differ by exactly `t` along the shared unit normal.
            let (p0, p1) = (
                plane_origin_projection(c).unwrap(),
                plane_origin_projection(moved).unwrap(),
            );
            let d2 = (0..3).fold(Rat::from_int(0), |a, i| {
                let d = p1[i].checked_sub(p0[i]).unwrap();
                a.checked_add(d.checked_mul(d).unwrap()).unwrap()
            });
            assert_eq!(
                d2,
                t.checked_mul(t).unwrap(),
                "{name} moved the wrong distance"
            );
        }
    }

    /// ★★★ **Two steps compose into one, exactly.** This is the property the whole thing exists for:
    /// a cap raised `7.7` and a cap raised `1.1` then `6.6` record the *same* plane, because
    /// `11/10 + 66/10 = 77/10` holds in rationals where it does not in `f64`.
    #[test]
    fn pushing_twice_lands_where_pushing_once_does() {
        let d = |x: f64| Rat::from_decimal(x).unwrap();
        assert_ne!(1.1f64 + 6.6, 7.7, "the fixture must discriminate");
        for raw in [[0, 0, 1, -2], [3, 0, -4, -5], [0, 0, 10, -21]] {
            let c = ints(raw);
            let once = plane_offset(c, d(7.7)).unwrap();
            let twice = plane_offset(plane_offset(c, d(1.1)).unwrap(), d(6.6)).unwrap();
            assert_eq!(once, twice, "{raw:?}");
        }
    }

    /// Pushing by zero is the plane itself, and pushing back undoes it — the sign convention is the
    /// one thing here a reader has to trust, so it is pinned in both directions.
    #[test]
    fn pushing_by_zero_and_pushing_back_are_identities() {
        // Every normal here has a rational length; `[7,-13,5]` does not, and `plane_offset`
        // declines it outright — that is the next test's business.
        for raw in [[0, 0, 1, -2], [3, 0, -4, -5], [5, 12, 0, 91]] {
            let c = canonical_plane_coeffs(ints(raw)).unwrap();
            assert_eq!(plane_offset(c, Rat::from_int(0)), Some(c), "{raw:?}");
            let t = Rat::new(-7, 2).unwrap();
            let there = plane_offset(c, t).unwrap();
            assert_eq!(
                plane_offset(there, Rat::from_int(0).checked_sub(t).unwrap()),
                Some(c),
                "{raw:?} did not come back"
            );
        }
    }

    /// **An irrational normal has no exact offset** — `|n|` is what must be rational, and `[1,1,1]`
    /// is the smallest thing that is not. The caller keeps its f64 path; nothing is approximated.
    #[test]
    fn a_plane_whose_normal_has_no_rational_length_declines() {
        assert_eq!(plane_offset(ints([1, 1, 1, -7]), Rat::from_int(1)), None);
        assert_eq!(plane_offset(ints([1, 2, 3, 0]), Rat::from_int(1)), None);
        // A 3-4-5 direction does have one, so this is a statement about lengths, not about tilt.
        assert!(plane_offset(ints([3, 4, 0, -5]), Rat::from_int(1)).is_some());
        // Not a plane.
        assert_eq!(plane_offset(ints([0, 0, 0, 1]), Rat::from_int(1)), None);
    }

    #[test]
    fn the_projected_origin_lies_exactly_on_its_plane() {
        for (raw, name) in [
            ([0, 0, 1, -2], "z = 2"),
            ([0, 0, 10, -21], "z = 2.1"),
            ([3, 0, -4, -5], "3-4-5 tilt"),
            ([1, 1, 1, -7], "diagonal"),
            ([5, 12, 0, -13], "5-12-13"),
            ([7, -13, 5, 91], "ugly"),
            ([0, 0, 1, 0], "through the origin"),
        ] {
            let c = ints(raw);
            let p = plane_origin_projection(c).unwrap_or_else(|| panic!("{name} has a projection"));
            let on = (0..3).fold(c[3], |acc, i| {
                acc.checked_add(c[i].checked_mul(p[i]).unwrap()).unwrap()
            });
            assert_eq!(
                on,
                Rat::from_int(0),
                "{name}: p = {p:?} is off its own plane"
            );
        }
    }

    /// ★★★ **One plane, three spellings, one point.** The design rests on this: the canonical form
    /// carries no direction (`push_surface_with_coeffs` returns a `flipped` flag for exactly that
    /// reason), and a plane is scale-invariant, so an origin derived from the coefficients would be
    /// worthless if it moved when the coefficients were negated or scaled.
    #[test]
    fn the_projection_does_not_depend_on_how_the_plane_is_spelled() {
        let want = plane_origin_projection(ints([0, 0, 1, -3])).expect("a plane");
        for raw in [[0, 0, -1, 3], [0, 0, 2, -6], [0, 0, -5, 15]] {
            assert_eq!(plane_origin_projection(ints(raw)), Some(want), "{raw:?}");
        }
        // And with denominators: 11/10·x − 7/2 = 0 is 11x − 35 = 0, both giving p = (35/11, 0, 0).
        let fracs = [
            Rat::new(11, 10).unwrap(),
            Rat::from_int(0),
            Rat::from_int(0),
            Rat::new(-7, 2).unwrap(),
        ];
        assert_eq!(
            plane_origin_projection(fracs),
            plane_origin_projection(ints([11, 0, 0, -35]))
        );
    }

    /// The projection is the **nearest** point of the plane to the origin, which is what makes it a
    /// sensible frame origin rather than merely a reproducible one: `p` is parallel to `n`, so no
    /// other point of the plane is closer.
    #[test]
    fn the_projected_origin_is_the_nearest_point_of_its_plane() {
        let c = ints([1, 1, 1, -7]);
        let p = plane_origin_projection(c).expect("a plane");
        let d2 = |q: [Rat; 3]| {
            (0..3).fold(Rat::from_int(0), |a, i| {
                a.checked_add(q[i].checked_mul(q[i]).unwrap()).unwrap()
            })
        };
        // Step along an in-plane direction (n × ê is perpendicular to n) and the distance grows.
        for step in [Rat::from_int(1), Rat::new(-3, 7).unwrap()] {
            let dir = [Rat::from_int(1), Rat::from_int(-1), Rat::from_int(0)]; // ⊥ to (1,1,1)
            let q = [0, 1, 2].map(|i| p[i].checked_add(dir[i].checked_mul(step).unwrap()).unwrap());
            assert!(
                d2(q) > d2(p),
                "stepping by {step:?} did not move away from the origin"
            );
        }
    }

    /// Coefficients that are not a plane have no projection, and neither does an `i128` overflow —
    /// both are the kernel's ordinary demotion, not a failure.
    #[test]
    fn a_non_plane_and_an_overflow_both_decline() {
        assert_eq!(plane_origin_projection(ints([0, 0, 0, 5])), None);
        assert_eq!(plane_origin_projection(ints([0, 0, 0, 0])), None);
        let huge = i128::MAX / 3;
        assert_eq!(plane_origin_projection(ints([huge, huge, huge, -1])), None);
    }

    #[test]
    fn one_plane_written_at_any_scale_canonicalizes_to_one_vector() {
        let want = ints([1, 2, 0, -3]);
        for scale in [1, 2, 7, -1, -13] {
            let scaled = ints([scale, 2 * scale, 0, -3 * scale]);
            assert_eq!(canonical_plane_coeffs(scaled), Some(want), "scale {scale}");
        }
        // And with denominators: 11/10·x + 3/5·y − 7/2 = 0 is 11x + 6y − 35 = 0.
        let fracs = [
            Rat::new(11, 10).unwrap(),
            Rat::new(3, 5).unwrap(),
            Rat::from_int(0),
            Rat::new(-7, 2).unwrap(),
        ];
        assert_eq!(
            canonical_plane_coeffs(fracs),
            Some(ints([11, 6, 0, -35])),
            "denominators are cleared and the content divided out"
        );
    }

    #[test]
    fn the_sign_convention_picks_the_first_nonzero_component() {
        assert_eq!(
            canonical_plane_coeffs(ints([0, -2, 4, 6])),
            Some(ints([0, 1, -2, -3])),
            "a plane and its negation are one plane, so one of the two spellings wins"
        );
    }

    /// ★★★ **The rule this function cannot enforce, made executable.**
    ///
    /// Two faces of the plane `x = 3` reach `nacre_geom::Plane::coefficients()` as these two f64
    /// vectors, because the un-normalized normal scales with the face's size and `d` is a rounded
    /// product. Lifting them is lossless — and still gives two different planes, because the
    /// *values* differ: `2.2 × 3` is not `6.6000000000000005`. Canonicalization cannot undo a
    /// rounding that already happened, so the coefficients have to be rational from construction.
    #[test]
    fn lifting_rounded_f64_coefficients_does_not_merge_them() {
        let lift = |v: [f64; 4]| v.map(|x| Rat::try_from_f64(x).expect("finite"));
        let a = canonical_plane_coeffs(lift([2.2, 0.0, 0.0, -6.6000000000000005])).unwrap();
        let b = canonical_plane_coeffs(lift([13.2, 0.0, 0.0, -39.599999999999994])).unwrap();
        assert_ne!(
            a, b,
            "lifting a rounded coefficient carries the rounding in"
        );

        // Built from what the user *wrote*, and multiplied **in the rationals**, the same two
        // faces agree exactly.
        //
        // ★★ The `d` in `−3·k` has to be a `Rat` product. Writing `d(-3.0 * k)` instead puts the
        // multiplication back in f64, `−3.0 × 2.2` comes out `−6.6000000000000005`, and its
        // shortest decimal is that — not `−6.6`. A value the caller computed is not a value the
        // caller wrote, and `from_decimal` cannot tell them apart.
        let d = |x: f64| Rat::from_decimal(x).expect("decimal");
        let by_hand = |k: f64| {
            let k = d(k);
            [
                k,
                Rat::from_int(0),
                Rat::from_int(0),
                k.checked_mul(Rat::from_int(-3)).expect("small"),
            ]
        };
        assert_eq!(
            canonical_plane_coeffs(by_hand(2.2)),
            canonical_plane_coeffs(by_hand(13.2)),
            "rational construction is scale-independent"
        );
    }

    fn on_plane(c: [Rat; 4], p: [Rat; 3]) -> bool {
        let mut acc = c[3];
        for i in 0..3 {
            acc = acc.checked_add(c[i].checked_mul(p[i]).unwrap()).unwrap();
        }
        acc == Rat::from_int(0)
    }

    /// ★ **The moved plane must contain the moved points** — checked against a *separate*
    /// point transform, so this cannot pass by restating the plane formula.
    #[test]
    fn a_moved_plane_still_contains_its_moved_points() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        // A slanted plane, 2x − 3y + z − 6 = 0, and three points on it.
        let plane = ints([2, -3, 1, -6]);
        let pts = [
            [r(3, 1), r(0, 1), r(0, 1)],
            [r(0, 1), r(0, 1), r(6, 1)],
            [r(1, 2), r(1, 1), r(8, 1)],
        ];
        for p in pts {
            assert!(on_plane(plane, p), "fixture point is on the plane");
        }

        let cases = [
            Isometry::translation([r(1, 2), r(-3, 1), r(7, 5)]),
            Isometry::rotation(Rotation {
                axis: Axis::X,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
            }),
            Isometry::rigid(
                Rotation {
                    axis: Axis::Z,
                    point: [r(1, 1), r(2, 1), r(0, 1)],
                    angle: Angle::from_deg(Rat::from_int(270)).unwrap(),
                },
                [r(0, 1), r(5, 1), r(-1, 2)],
            ),
        ];
        for iso in cases {
            let moved = iso.plane_coeffs(plane).expect("exact motion");
            for p in pts {
                assert!(
                    on_plane(moved, iso.point_rat(p).expect("exact motion")),
                    "moved plane {moved:?} must contain the moved point"
                );
            }
        }
    }

    #[test]
    fn a_mirrored_plane_still_contains_its_mirrored_points() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        let plane = ints([2, -3, 1, -6]);
        let (axis, offset) = (Axis::Z, r(5, 4));
        let moved = mirror_plane_coeffs(plane, axis, offset).unwrap();
        for p in [[r(3, 1), r(0, 1), r(0, 1)], [r(0, 1), r(0, 1), r(6, 1)]] {
            let q = mirror_point_rat(p, axis, offset).unwrap();
            assert!(on_plane(moved, q), "mirrored plane must contain {q:?}");
        }
    }

    /// A non-90° rotation has no rational `cos`/`sin`, so there is nothing exact to return —
    /// ★ **for the point as well as the plane**. A caller that could move one description and
    /// silently drop the other would leave a plane whose stored points are no longer on it.
    #[test]
    fn an_inexact_rotation_declines() {
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        assert!(!iso.is_exact());
        assert_eq!(iso.plane_coeffs(ints([0, 0, 1, -3])), None);
        assert_eq!(iso.point_rat([Rat::from_int(1); 3]), None);
    }

    #[test]
    fn overflow_demotes_instead_of_panicking() {
        // Four coprime denominators near 10^10: their lcm is their product, ~10^40, past i128.
        let r = |den| Rat::new(1, den).unwrap();
        let coeffs = [
            r(10_000_000_019),
            r(10_000_000_033),
            r(10_000_000_061),
            r(10_000_000_069),
        ];
        assert_eq!(canonical_plane_coeffs(coeffs), None);
    }
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

    /// **The property `from_decimal` rests on**: the shortest decimal of an f64,
    /// read as an exact rational and realized again, gives back the same bits.
    ///
    /// It holds *because* `to_f64` is correctly rounded, and only because of that. The
    /// decimal is by construction a value whose nearest f64 is `x`; nearest rounding
    /// therefore has no choice. Rounding numerator and denominator separately first —
    /// what this function used to do — fails **17.4%** of the values below (measured),
    /// since a 17-digit decimal has a numerator past 2⁵³. That is the whole reason this
    /// cell touches `to_f64` at all.
    #[test]
    fn a_shortest_decimal_realizes_back_to_its_own_f64() {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let (mut tried, mut declined) = (0, 0);
        for _ in 0..200_000 {
            // A random sign and mantissa over a CAD-plausible exponent range.
            let bits = next();
            let exp = 1023 - 60 + bits % 121;
            let x = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (exp << 52));
            match Rat::from_decimal(x) {
                Some(r) => {
                    tried += 1;
                    assert_eq!(r.to_f64(), x, "{x:?} → {r:?} → {:?}", r.to_f64());
                }
                None => declined += 1,
            }
        }
        // Nothing declines over this exponent range; the guard is here so a corpus that
        // drifted out of i128 could not turn this test vacuous.
        assert!(tried > 100_000, "{tried} tried, {declined} declined");
    }

    /// The same property, but over **every** finite f64 rather than the CAD-plausible
    /// band: subnormals, `MIN_POSITIVE`, `MAX`, and the whole exponent range in between.
    /// Where the power of ten leaves i128 the answer is `None` and the caller keeps its
    /// f64 — what must never happen is a `Some` that realizes back to a *different*
    /// number, because that would move a coordinate this cell has no business moving.
    #[test]
    fn from_decimal_never_lies_anywhere_in_the_finite_range() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut corpus = vec![
            0.0,
            -0.0,
            f64::MIN_POSITIVE,
            f64::MAX,
            f64::MIN,
            f64::from_bits(1), // the smallest subnormal
            1.1,
            7.7,
            1e-30,
            1e300,
        ];
        for _ in 0..200_000 {
            let bits = next();
            // Every exponent a finite f64 can have, subnormals included.
            let exp = (bits % 2047) << 52;
            corpus.push(f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | exp));
        }
        let mut declined = 0;
        for x in corpus {
            match Rat::from_decimal(x) {
                Some(r) => assert_eq!(r.to_f64(), x, "{x:e} → {r:?}"),
                None => declined += 1,
            }
        }
        assert!(
            declined > 0,
            "the i128 limit should bite somewhere out here"
        );
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(Rat::from_decimal(bad), None, "{bad}");
        }
    }

    /// Where the i128 limit actually bites, pinned so the doc comment stays true and a
    /// caller can tell whether its dimensions are anywhere near it. They are not: a CAD
    /// model in millimetres sits around `1e0`, sixty orders of magnitude inside.
    ///
    /// The two edges differ because a 17-digit mantissa spends 16 of its own powers of
    /// ten on the fraction, so it runs out on the small side first. Note that exponent
    /// *notation* is no obstacle in itself — `1e-30` is accepted; only the magnitude is.
    #[test]
    fn the_decimal_window_is_wide_and_its_edges_are_where_they_should_be() {
        let at = |s: &str| Rat::from_decimal(s.parse::<f64>().unwrap()).is_some();
        for e in -22..=38 {
            assert!(at(&format!("1.2345678901234567e{e}")), "17 digits at 1e{e}");
        }
        assert!(!at("1.2345678901234567e-23"), "17 digits at 1e-23");
        assert!(!at("1.2345678901234567e39"), "17 digits at 1e39");
        assert!(
            at("1e-38") && at("1e-30") && !at("1e-39"),
            "a short decimal reaches further down"
        );
    }

    /// The ends of the `i128` range, where the two operands' bounds differ.
    ///
    /// `i128::MIN.unsigned_abs()` is `2¹²⁷` **exactly**, one past what the denominator
    /// may be — an asymmetry easy to assert away and, when it was, a debug-only panic
    /// on a value the algorithm handles correctly. The answers here are exact powers of
    /// two and their neighbours, so they can be written down rather than approximated.
    #[test]
    fn the_ends_of_the_range_realize_exactly() {
        assert_eq!(Rat::from_int(i128::MIN).to_f64(), -(2f64.powi(127)));
        assert_eq!(Rat::from_int(i128::MAX).to_f64(), 2f64.powi(127)); // rounds up to 2¹²⁷
        assert_eq!(Rat::new(i128::MIN, 2).unwrap().to_f64(), -(2f64.powi(126)));
        assert_eq!(Rat::new(i128::MIN, i128::MAX).unwrap().to_f64(), -1.0);
        assert_eq!(Rat::new(1, i128::MAX).unwrap().to_f64(), 2f64.powi(-127));
        // Every one is finite and signed the way its numerator is.
        for (n, d) in [
            (i128::MIN, 3),
            (i128::MAX, 7),
            (-1, i128::MAX),
            (i128::MIN + 1, 1),
        ] {
            let q = Rat::new(n, d).unwrap().to_f64();
            assert!(q.is_finite() && (q < 0.0) == (n < 0), "{n}/{d} -> {q:e}");
        }
    }

    /// **The identity this whole cell exists for.** A dimension split into two and
    /// stacked must land where the undivided one does. Lifting the f64 *values* cannot
    /// give that — the drift is already inside them — and reading the decimals can.
    #[test]
    fn a_split_dimension_stacks_back_to_the_whole_one() {
        let split = Rat::from_decimal(1.1)
            .unwrap()
            .checked_add(Rat::from_decimal(6.6).unwrap())
            .unwrap();
        assert_eq!(split, Rat::new(77, 10).unwrap());
        assert_eq!(split, Rat::from_decimal(7.7).unwrap());
        assert_eq!(split.to_f64(), 7.7);
        // The f64 arithmetic this replaces, and the exact-but-binary lift that does not
        // help either — both land an ulp away.
        assert_ne!(1.1 + 6.6, 7.7);
        assert_ne!(
            Rat::try_from_f64(1.1)
                .unwrap()
                .checked_add(Rat::try_from_f64(6.6).unwrap())
                .unwrap(),
            Rat::try_from_f64(7.7).unwrap()
        );
    }

    proptest! {
        /// Nearest-ness, checked against the definition rather than against a second
        /// implementation of it: `q` is the nearest f64 to `n/d` exactly when no
        /// neighbour of `q` is closer, and the comparison `|n/d − a/b| ≤ |n/d − c/e|`
        /// is decidable in exact rationals. Held to the range where those stay in i128.
        #[test]
        fn to_f64_lands_on_the_nearest_f64(n in -(1i128 << 60)..(1i128 << 60), d in 1i128..(1i128 << 60)) {
            prop_assume!(n != 0);
            let r = Rat::new(n, d).unwrap();
            let q = r.to_f64();
            let zero = Rat::from_int(0);
            let err = |y: f64| {
                let e = Rat::try_from_f64(y).and_then(|yr| r.checked_sub(yr))?;
                if e < zero { zero.checked_sub(e) } else { Some(e) }
            };
            let here = err(q).expect("in range");
            for nb in [f64::from_bits(q.to_bits() + 1), f64::from_bits(q.to_bits() - 1)] {
                if let Some(there) = err(nb) {
                    prop_assert!(here <= there, "{r:?}: {q:?} is not nearest ({nb:?} is closer)");
                }
            }
        }
    }

    proptest! {
        /// ★★★★★ **The two derivations must be the same function.**
        ///
        /// `plane_name_exact` runs the `Rat` route first and only falls back, so wherever the
        /// narrow one answers, the wide one is never consulted — and an error in it would sit
        /// there unseen until the day it *is* consulted, on inputs no test covers. Calling both on
        /// the same inputs is the only way to say they agree.
        ///
        /// ★ `plane_name_big` is `pub(crate)` for exactly this reason: a fallback hidden behind
        /// its filter cannot be tested against it.
        #[test]
        fn the_wide_derivation_answers_what_the_narrow_one_does(
            xs in prop::array::uniform9(-(1i64 << 20)..(1i64 << 20)),
            ds in prop::array::uniform9(1i64..(1i64 << 20)),
        ) {
            let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
            let (a, b, c) = (
                [r(0), r(1), r(2)],
                [r(3), r(4), r(5)],
                [r(6), r(7), r(8)],
            );
            let narrow = plane_through_points(a, b, c);
            let wide = plane_name_big(a, b, c);
            if let Some(n) = narrow {
                // ★ The invariant rides along: an answer the narrow route reached fits `i128`
                // by construction, so the wide route must store it `Narrow` — same value, same
                // representation, structural equality.
                prop_assert_eq!(wide.clone(), Some(PlaneName::Narrow(n)),
                    "narrow answered but wide disagrees");
            }
            // ★ And the wide one, whenever it answers narrowly at all, answers about a plane
            // these points are actually on — checked in the rationals, no tolerance.
            //
            // ★★★★★ **The residual can overflow even when the name and the points both fit**, and
            // this test found that by asserting it could not. `c · p` multiplies a canonical
            // coefficient by a point coordinate, so it needs the *sum* of their widths — which is
            // exactly the population `Model::push_surface_with_coeffs` cannot verify either. An
            // unevaluable check is not a failed one, here as there: skip it, never fail on it.
            if let Some(w) = wide.as_ref().and_then(|n| n.narrow()) {
                for p in [a, b, c] {
                    let residual = (|| {
                        let mut acc = w[3];
                        for k in 0..3 {
                            acc = acc.checked_add(w[k].checked_mul(p[k])?)?;
                        }
                        Some(acc)
                    })();
                    if let Some(acc) = residual {
                        prop_assert_eq!(acc, Rat::from_int(0), "a point is off the derived plane");
                    }
                }
            }
        }
    }

    /// ★★★★ **The width the narrow route cannot reach, and the wide one can.**
    ///
    /// A triple whose reduced denominators are coprime — one a power of two, one a power of five,
    /// which is what decimal arithmetic produces once it reduces — needs their product to state
    /// the plane, and `(b − a) × (c − a)` needs it squared. `plane_through_points` gives up there;
    /// the answer is small, and this is the case that says so.
    #[test]
    fn a_plane_the_narrow_route_gives_up_on_is_still_named() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        let a = [r(1, 1 << 53), r(0, 1), r(0, 1)];
        let b = [r(0, 1), r(1, 5i128.pow(23)), r(0, 1)];
        let c = [r(0, 1), r(0, 1), r(1, (1 << 40) * 5i128.pow(11))];
        assert_eq!(
            plane_through_points(a, b, c),
            None,
            "the narrow route was expected to overflow on coprime denominators"
        );
        let name = plane_name_exact(a, b, c).expect("the wide route names it");
        // The canonical answer here is small (the doc above says so) — the invariant demands it
        // come back `Narrow`.
        let wide = name.narrow().expect("a small answer must be stored Narrow");
        for p in [a, b, c] {
            let mut acc = wide[3];
            for k in 0..3 {
                acc = acc
                    .checked_add(wide[k].checked_mul(p[k]).expect("no overflow"))
                    .expect("no overflow");
            }
            assert_eq!(acc, Rat::from_int(0), "a point is off the derived plane");
        }
    }

    proptest! {
        /// ★★★★★ **The two solves must be the same function** — the point-side twin of
        /// `the_wide_derivation_answers_what_the_narrow_one_does`, and for the same reason: a
        /// wide route that is only ever consulted where the narrow one declined would carry an
        /// error unseen until the day it is consulted. Calling both on the same inputs is the
        /// only way to say they agree.
        #[test]
        fn the_wide_solve_answers_what_the_narrow_one_does(
            xs in prop::array::uniform12(-(1i64 << 20)..(1i64 << 20)),
            ds in prop::array::uniform12(1i64..(1i64 << 20)),
        ) {
            let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
            let rows: [[Rat; 4]; 3] =
                core::array::from_fn(|k| core::array::from_fn(|j| r(4 * k + j)));
            let names = rows.map(PlaneName::Narrow);
            let wide = three_planes_big([&names[0], &names[1], &names[2]]);

            if let Some(n) = three_planes_rat(rows) {
                // ★ The invariant rides along: a point the narrow route reached fits `Rat` by
                // construction, so the wide route must store it `Narrow` — same value, same
                // representation, structural equality.
                prop_assert_eq!(wide.clone(), Some(MeetPoint::Narrow(n)),
                    "narrow answered but wide disagrees");
            }

            // ★★ And wherever the wide one answers **at all**, the point it names is on all three
            // planes — checked in the rationals, no tolerance. This half covers the `Wide` arm,
            // which the differential above cannot reach: the narrow route is silent there by
            // definition, so agreement says nothing and only the residual does.
            if let Some(w) = wide.as_ref() {
                use num_bigint::BigInt;
                use num_rational::Ratio;
                use num_traits::Zero;
                let big = |n: i128, d: i128| Ratio::new(BigInt::from(n), BigInt::from(d));
                let coord = |i: usize| -> Ratio<BigInt> {
                    match w {
                        MeetPoint::Narrow(p) => big(p[i].numer(), p[i].denom()),
                        MeetPoint::Wide(p) => Ratio::new(p[i].0.clone(), p[i].1.clone()),
                    }
                };
                let point = [coord(0), coord(1), coord(2)];
                for row in &rows {
                    let mut acc = big(row[3].numer(), row[3].denom());
                    for (j, x) in point.iter().enumerate() {
                        acc += big(row[j].numer(), row[j].denom()) * x;
                    }
                    prop_assert!(acc.is_zero(), "the solved point is off one of the planes");
                }
            }
        }
    }

    /// ★★★★★ **An overflowing intermediate no longer costs the answer** — the repair of open
    /// item 0, locked on the fixture that used to document the defect.
    ///
    /// These three planes meet at `(1, 1, 1)`, which fits `Rat` with room to spare. But their
    /// coefficients carry coprime denominators — a power of two and a power of five, what decimal
    /// arithmetic produces once it reduces — and the determinant is their product: `2⁹³·5³⁴`, a
    /// ~172-bit denominator. The narrow Cramer multiplies before it can reduce, so it overflows;
    /// this fixture used to assert the resulting decline, and now asserts the fallback answers
    /// through the integer core instead — the hand-known point, so this is an independent oracle
    /// and not the two routes agreeing with each other.
    ///
    /// ★ The private narrow route still declines here (asserted), so the fixture keeps proving
    /// the fallback is *reached*, not merely present.
    #[test]
    fn an_overflowing_intermediate_no_longer_costs_the_answer() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        let rows = [
            [r(1, 1 << 53), r(0, 1), r(0, 1), r(-1, 1 << 53)],
            [r(0, 1), r(1, 5i128.pow(23)), r(0, 1), r(-1, 5i128.pow(23))],
            [
                r(0, 1),
                r(0, 1),
                r(1, (1 << 40) * 5i128.pow(11)),
                r(-1, (1 << 40) * 5i128.pow(11)),
            ],
        ];
        assert_eq!(
            three_planes_rat_narrow(rows),
            None,
            "the narrow route was expected to overflow on coprime denominators — \
             without that this fixture no longer exercises the fallback"
        );
        assert_eq!(
            three_planes_rat(rows),
            Some([Rat::from_int(1); 3]),
            "the point fits `Rat`, and the solve now says so"
        );
        let names = rows.map(PlaneName::Narrow);
        let found = three_planes_big([&names[0], &names[1], &names[2]])
            .expect("the wide solve finds the point");
        assert_eq!(
            found,
            MeetPoint::Narrow([Rat::from_int(1); 3]),
            "the point fits `Rat` — the invariant demands it come back Narrow"
        );
        assert_eq!(found.width_bits(), 1, "1/1 is one bit wide");
    }

    /// The two remaining meanings of [`three_planes_rat`]'s `None`, each still honest:
    /// no unique point (parallel planes), and a point that truly does not fit `Rat` — which
    /// [`three_planes_big`] tells apart by answering `Wide`.
    #[test]
    fn the_solve_still_declines_what_it_should() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        // Parallel pair: x = 0 and x = 1 — no unique point, both routes say so.
        let parallel = [
            [r(1, 1), r(0, 1), r(0, 1), r(0, 1)],
            [r(1, 1), r(0, 1), r(0, 1), r(-1, 1)],
            [r(0, 1), r(0, 1), r(1, 1), r(0, 1)],
        ];
        assert_eq!(
            three_planes_rat(parallel),
            None,
            "parallel planes meet nowhere"
        );
        // Narrow rows whose meeting point is wider than `Rat`: x + y = 2⁻¹⁰⁰, x − y = 5⁻⁵⁰
        // put a denominator of 2¹⁰¹·5⁵⁰ (~217 bits) on x. `three_planes_rat` declines;
        // the wide twin answers, and answers `Wide` — the causes stay told apart.
        let wide_point = [
            [r(1, 1), r(1, 1), r(0, 1), r(-1, 1 << 100)],
            [r(1, 1), r(-1, 1), r(0, 1), r(-1, 5i128.pow(50))],
            [r(0, 1), r(0, 1), r(1, 1), r(0, 1)],
        ];
        assert_eq!(
            three_planes_rat(wide_point),
            None,
            "a point no `Rat` can hold is a decline, not an answer"
        );
        let names = wide_point.map(PlaneName::Narrow);
        assert!(
            matches!(
                three_planes_big([&names[0], &names[1], &names[2]]),
                Some(MeetPoint::Wide(_))
            ),
            "the twin names the cause: the point exists and is wide"
        );
    }

    proptest! {
        /// Totality: wherever the integer core answers `Narrow`, [`three_planes_rat`] answers
        /// the same. ★ A **wiring** lock, stated as such: after the refactor both routes
        /// converge on `three_planes_int`, so this pins the `.or_else` plumbing (lift,
        /// `narrow()` return) rather than serving as an independent oracle — that role belongs
        /// to `an_overflowing_intermediate_no_longer_costs_the_answer`'s hand-known point and
        /// to the agreement direction the ops-side `point_width` invariant keeps.
        #[test]
        fn the_solve_answers_wherever_the_answer_is_narrow(
            rows in proptest::array::uniform3(proptest::array::uniform4((-9i128..=9, 1u32..=60))),
        ) {
            let rat_rows: [[Rat; 4]; 3] =
                rows.map(|row| row.map(|(n, e)| Rat::new(n, 1i128 << e).unwrap()));
            let names = rat_rows.map(PlaneName::Narrow);
            let big = three_planes_big([&names[0], &names[1], &names[2]]);
            let expect = match &big {
                Some(MeetPoint::Narrow(p)) => Some(*p),
                _ => None, // no unique point, or truly wide — the honest declines
            };
            prop_assert_eq!(three_planes_rat(rat_rows), expect);
        }
    }

    /// ★★★★★ **The negative control for [`MeetPoint::width_bits`]** — without it, a corpus that
    /// reports "nothing over 127 bits" is indistinguishable from a dead probe.
    ///
    /// A `Wide` carrier states a plane at `x = 2²⁰⁰`, which no `Rat` can hold; the meeting point
    /// inherits that width and the fork has to take its other branch. So the instrument can say
    /// the other thing, and `over127 = 0` in a measurement means the population, not the meter.
    #[test]
    fn the_width_meter_reports_a_point_no_rat_can_hold() {
        use num_bigint::BigInt;
        let far = BigInt::from(1) << 200;
        let wide = PlaneName::Wide([BigInt::from(1), BigInt::from(0), BigInt::from(0), -&far]);
        let r = |v: [i128; 4]| PlaneName::Narrow(v.map(Rat::from_int));
        let (py, pz) = (r([0, 1, 0, 0]), r([0, 0, 1, 0]));
        let found = three_planes_big([&wide, &py, &pz]).expect("three planes meet at (2²⁰⁰, 0, 0)");
        assert!(
            matches!(found, MeetPoint::Wide(_)),
            "a 201-bit coordinate cannot be Narrow"
        );
        assert_eq!(
            found.width_bits(),
            201,
            "the meter reads the coordinate's width"
        );
        assert_eq!(found.narrow(), None);
    }

    /// ★★★★★ **A canonical answer wider than `i128` is a name now, not a `None`** — and two
    /// statements of that plane are one value. This is the S2 vessel change in miniature; the
    /// interning consequence is locked on the model side
    /// (`a_wide_plane_interns_but_opens_no_shortcut` in nacre-topo).
    #[test]
    fn a_plane_too_wide_for_i128_is_named_wide() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        // Cross-product terms multiply two coordinates' numerators, so two ~2^90 coprime
        // numerators push the canonical coefficients past i128 with no content to divide out.
        let big1 = (1i128 << 90) + 1;
        let big2 = (1i128 << 90) + 3;
        let a = [r(big1, 3), r(big2, 7), r(0, 1)];
        let b = [r(-big2, 5), r(big1, 11), r(0, 1)];
        let c = [r(1, 13), r(1, 17), r(1, 19)];
        assert_eq!(
            plane_through_points(a, b, c),
            None,
            "expected the narrow route to overflow"
        );
        let name = plane_name_exact(a, b, c).expect("collinear it is not — it must be named");
        // ★ Wide carries identity only: the arithmetic shortcuts' door stays shut (S2/S4
        // boundary — a wide name must NOT open a frame).
        assert!(
            name.narrow().is_none(),
            "an answer past i128 must be stored Wide"
        );
        // ★ Same plane, two spellings (a permuted triple) — one value, structurally.
        let permuted = plane_name_exact(b, c, a).expect("the same plane, permuted");
        assert_eq!(
            name, permuted,
            "two statements of one wide plane must be one value"
        );
    }

    proptest! {
        /// **The two orientation arms must be the same function** — same argument as
        /// `the_wide_derivation_answers_what_the_narrow_one_does`: `orient2d_rat` runs the `Rat`
        /// route first, so wherever it answers the `BigInt` arm is never consulted, and an error
        /// there would wait for the first overflowing input. `orient2d_big` is `pub(crate)` for
        /// exactly this call.
        #[test]
        fn the_big_orientation_answers_what_the_narrow_one_does(
            xs in prop::array::uniform6(-(1i64 << 20)..(1i64 << 20)),
            ds in prop::array::uniform6(1i64..(1i64 << 20)),
        ) {
            let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
            let (a, b, c) = ([r(0), r(1)], [r(2), r(3)], [r(4), r(5)]);
            // Small operands: the narrow route always answers, so this compares the arms.
            prop_assert_eq!(orient2d_rat(a, b, c), orient2d_big(a, b, c));
        }
    }

    /// ★★ **The population the narrow route cannot reach — and the reason the sign is total.**
    ///
    /// Two decimal-window denominators of `10^22` put the determinant's product denominator at
    /// `10^44`, past `i128` — precisely what `from_decimal`'d profile coordinates produce. The
    /// answers are hand-computable: `(0,0) → (1/d,0) → (0,1/d)` turns counter-clockwise, its
    /// mirror clockwise, and a doubled point on the same ray is collinear.
    #[test]
    fn an_orientation_past_i128_still_gets_its_sign() {
        let d = 10i128.pow(22);
        let r = |n: i128, den: i128| Rat::new(n, den).unwrap();
        let (o, x, y) = (
            [r(0, 1), r(0, 1)],
            [r(1, d), r(0, 1)],
            [r(0, 1), r(1, d + 1)],
        );
        // The narrow route must actually be dead here, or this pins nothing (self-qualification).
        assert!(
            r(1, d).checked_mul(r(1, d + 1)).is_none(),
            "the product denominator was expected to overflow i128"
        );
        assert_eq!(orient2d_rat(o, x, y), 1, "counter-clockwise");
        assert_eq!(orient2d_rat(o, y, x), -1, "clockwise");
        let far = [r(2, d), r(2, d)];
        let near = [r(1, d + 1), r(1, d + 1)];
        assert_eq!(orient2d_rat(o, near, far), 0, "one ray, three points");
    }

    proptest! {
        /// **The residual's two arms must be the same function** — `plane_residual_sign` runs the
        /// `Rat` substitution first, so wherever it answers the `BigInt` arm is never consulted
        /// (the `orient2d` arrangement, one dimension up).
        #[test]
        fn the_big_residual_answers_what_the_narrow_one_does(
            // Small enough that the canonical name always fits i128 (the lifted integers stay
            // near 2^40 and the derivation's peak near 2^122), so `narrow()` below cannot shrug.
            xs in prop::array::uniform9(-(1i64 << 10)..(1i64 << 10)),
            ds in prop::array::uniform9(1i64..(1i64 << 10)),
        ) {
            let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
            let (a, b, c) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)], [r(6), r(7), r(8)]);
            if let Some(name) = plane_name_exact(a, b, c) {
                let coeffs = name.narrow().expect("small operands stay narrow");
                let ci = coeffs.map(|x| num_bigint::BigInt::from(x.numer()));
                // The naming points themselves: both arms must call them on-plane...
                for p in [a, b, c] {
                    prop_assert_eq!(plane_residual_sign(&name, p), 0);
                    prop_assert_eq!(residual_sign_big(&ci, p), 0);
                }
                // ...and an off-plane probe (a naming point pushed along the normal) must get
                // the same nonzero sign from both.
                let n = [coeffs[0], coeffs[1], coeffs[2]];
                if let Some(q) = (|| -> Option<[Rat; 3]> {
                    Some([
                        a[0].checked_add(n[0])?,
                        a[1].checked_add(n[1])?,
                        a[2].checked_add(n[2])?,
                    ])
                })() {
                    let s = plane_residual_sign(&name, q);
                    prop_assert_eq!(s, 1, "n·n > 0: the push is to the positive side");
                    prop_assert_eq!(residual_sign_big(&ci, q), s);
                }
            }
        }
    }

    /// ★★ **A `Wide` name still judges its point** — the population `plane_residual_sign` exists
    /// for. Same fixture as `a_plane_too_wide_for_i128_is_named_wide`: canonical coefficients
    /// past `i128`, no `Rat` route to fall back on.
    #[test]
    fn a_residual_against_a_wide_name_still_gets_its_sign() {
        let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
        let big1 = (1i128 << 90) + 1;
        let big2 = (1i128 << 90) + 3;
        let a = [r(big1, 3), r(big2, 7), r(0, 1)];
        let b = [r(-big2, 5), r(big1, 11), r(0, 1)];
        let c = [r(1, 13), r(1, 17), r(1, 19)];
        let name = plane_name_exact(a, b, c).expect("a genuine plane");
        assert!(name.narrow().is_none(), "the fixture must actually be Wide");
        for p in [a, b, c] {
            assert_eq!(plane_residual_sign(&name, p), 0, "a naming point is on it");
        }
        // A point off the plane along ±z (the naming triangle is not vertical: a and b span
        // z = 0 and c leaves it, so the normal has a z-component): opposite pushes must get
        // opposite, nonzero signs.
        let one = Rat::from_int(1);
        let up = [c[0], c[1], c[2].checked_add(one).unwrap()];
        let down = [c[0], c[1], c[2].checked_sub(one).unwrap()];
        let (su, sd) = (
            plane_residual_sign(&name, up),
            plane_residual_sign(&name, down),
        );
        assert_ne!(su, 0);
        assert_eq!(su, -sd, "opposite sides, opposite signs");
    }

    /// ★ S7: `three_planes_rat` — the exact corner of three rational planes. Hand-checkable
    /// fixture, a degenerate (line-sharing) triple, and an i128-overflow decline.
    #[test]
    fn three_rational_planes_meet_where_they_should() {
        let r = Rat::from_int;
        // x = 2, y = 3, z = 5  (as a·x + d = 0 rows: [1,0,0,-2] etc.)
        let pt = three_planes_rat([
            [r(1), r(0), r(0), r(-2)],
            [r(0), r(1), r(0), r(-3)],
            [r(0), r(0), r(1), r(-5)],
        ])
        .expect("axis planes meet in one point");
        assert_eq!(pt, [r(2), r(3), r(5)]);
        // A tilted but exact triple: x+y=1, x−y=0, z=7 ⇒ (1/2, 1/2, 7).
        let pt = three_planes_rat([
            [r(1), r(1), r(0), r(-1)],
            [r(1), r(-1), r(0), r(0)],
            [r(0), r(0), r(1), r(-7)],
        ])
        .expect("a fair triple");
        assert_eq!(pt, [Rat::new(1, 2).unwrap(), Rat::new(1, 2).unwrap(), r(7)]);
        // Three planes through one line (z = 0, y = 0, y + z = 0): det = 0 — no unique point.
        assert_eq!(
            three_planes_rat([
                [r(0), r(0), r(1), r(0)],
                [r(0), r(1), r(0), r(0)],
                [r(0), r(1), r(1), r(0)],
            ]),
            None
        );
        // Coefficients near the i128 edge overflow the narrow Cramer — and the answer still
        // fits `Rat` (denominator 2¹²⁶ + 1), so since open item 0 the fallback answers it.
        // This clause used to assert the decline; the decline was the defect.
        let big = Rat::from_int(1 << 126);
        let edge = [
            [big, r(1), r(0), r(-1)],
            [r(1), big, r(0), r(-1)],
            [r(0), r(0), big, big],
        ];
        assert_eq!(
            three_planes_rat_narrow(edge),
            None,
            "the narrow route overflows here"
        );
        let inv = Rat::new(1, (1 << 126) + 1).unwrap();
        assert_eq!(three_planes_rat(edge), Some([inv, inv, r(-1)]));
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

    /// The angle value stays exact, and its cos/sin realization is f64 — but a *correctly rounded*
    /// one, so the irrational-realization boundary now costs at most half an ulp rather than
    /// whatever the platform's libm happened to do.
    ///
    /// **`cos 45°` is the case to look at**: `√2/2` cannot be an f64, so the realization is
    /// genuinely lossy, and it must land on the nearest f64 to the truth. Checked by asking a far
    /// deeper realization whether anything is closer.
    #[test]
    fn realization_is_f64_while_angle_stays_exact() {
        let a = Angle::from_deg(Rat::from_int(45)).unwrap();
        assert_eq!(a.deg(), Rat::from_int(45)); // angle exact
        let (c, _) = a.cos_sin_f64();
        assert!(a.try_exact_cos_sin().is_none()); // √2/2 is not rational

        // Nothing is nearer: both neighbours are further from the deep truth than `c` is.
        let (deep, _) = a.cos_sin_at(512);
        let dist = |f: f64| BigFloat::from_f64(f, 512).sub(&deep, 512, HP_RM).abs();
        let here = dist(c);
        for nb in [
            f64::from_bits(c.to_bits() - 1),
            f64::from_bits(c.to_bits() + 1),
        ] {
            assert!(
                here.cmp(&dist(nb)).is_some_and(|s| s < 0),
                "a neighbour of {c:e} is nearer the truth"
            );
        }
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

        // 45° has no exact rational realization → None (falls to the rounded high-precision path).
        let a45 = Angle::from_deg(Rat::from_int(45)).unwrap();
        assert!(a45.try_exact_cos_sin().is_none());
        // The high-precision realization of 90° is *not* zero either — it is ~2⁻¹²⁸ — which is
        // exactly why the branch above exists rather than being an optimization.
        assert!(!a90.cos_sin_at(128).0.is_zero());
    }

    /// `cos_sin_f64` snaps the 90°-family to exact `0.0`/`±1.0`, and the exactly-representable
    /// values Niven allows off that family come out exact too.
    ///
    /// ★★ **`cos 60° == 0.5` is the visible proof that libm left.** `1/2` is one of the three
    /// rational values a rational-degree cosine can take, and it *is* an f64 — but reaching it
    /// through `(60.0 * PI / 180.0).cos()` does not land on it. Rounding the high-precision value
    /// does.
    #[test]
    fn cos_sin_f64_is_exact_where_the_true_value_is_representable() {
        let deg = |d| Angle::from_deg(Rat::from_int(d)).unwrap();
        assert_eq!(deg(0).cos_sin_f64(), (1.0, 0.0));
        assert_eq!(deg(90).cos_sin_f64(), (0.0, 1.0));
        assert_eq!(deg(180).cos_sin_f64(), (-1.0, 0.0));
        assert_eq!(deg(270).cos_sin_f64(), (0.0, -1.0));
        for (d, want) in [(60, 0.5), (120, -0.5), (240, -0.5), (300, 0.5)] {
            assert_eq!(deg(d).cos_sin_f64().0, want, "cos {d}°");
        }
        for (d, want) in [(30, 0.5), (150, 0.5), (210, -0.5), (330, -0.5)] {
            assert_eq!(deg(d).cos_sin_f64().1, want, "sin {d}°");
        }
    }

    /// **The 90°-family realizes with no error at all, and that zero is load-bearing.**
    ///
    /// `WitnessPoint::rotate_about` reads [`Angle::realization_error_of`] to decide whether a rotation
    /// contributes any tolerance. For `cos`/`sin` in `{0, ±1}` the f64 values *are* the true ones,
    /// and every product and difference downstream is exact too — which is why a quadrantal origin
    /// rotation stays at tol 0 and an axis-aligned model never leaves the exact predicate path.
    /// A nonzero here would not fail loudly; it would quietly move those models.
    ///
    /// Everything else must report *something* — but only where the true value is *not*
    /// representable. Since the realization became correctly rounded, `cos 60°` really is `0.5`
    /// exactly, so the corpus below has to avoid the four angles where that happens or it would be
    /// asserting a nonzero error that does not exist.
    #[test]
    fn only_the_quadrantal_family_realizes_exactly() {
        let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
        let err = |a: Angle| {
            let (c, s) = a.cos_sin_f64();
            a.realization_error_of(c, s)
        };
        for d in [0, 90, 180, 270] {
            assert_eq!(err(deg(d, 1)), (0.0, 0.0), "{d} deg");
        }
        // ★ And the zero is of the *pair*, not of the angle: hand a 90°-family angle some other
        // realization and it must be measured like anything else. Otherwise a caller that got its
        // cos/sin from somewhere other than `cos_sin_f64` would be handed an exactness claim that
        // does not hold of what it is holding.
        let a90 = deg(90, 1);
        assert!(a90.realization_error_of(1e-17, 1.0).0 > 0.0);
        for (n, d) in [(1, 1), (37, 1), (45, 1), (337, 1), (1, 3), (359999, 1000)] {
            let (dc, ds) = err(deg(n, d));
            assert!(
                dc > 0.0 && ds > 0.0,
                "{n}/{d} deg claimed an exact realization"
            );
            // ε-scale: a bound this large would mean the measurement, not the platform, is wrong.
            assert!(
                dc < 64.0 * f64::EPSILON && ds < 64.0 * f64::EPSILON,
                "{n}/{d} deg: {dc:e}"
            );
        }
    }

    /// **The f64 read-out is a re-encoding, not a computation — checked by round trip.**
    ///
    /// `set_precision(53, ToEven)` is where a value loses bits; everything after it is supposed to
    /// be pure bookkeeping over the mantissa words. So any `f64` put in must come back out
    /// unchanged. **The 32-bit-`Word` path is the one this is really for** — the mantissa is
    /// assembled across two words there, wasm is a 32-bit target, and no amount of reading the
    /// crate source substitutes for running it on the target that ships.
    #[test]
    fn the_f64_readout_round_trips() {
        let mut cases = vec![
            1.0,
            0.5,
            -0.5,
            0.9999999999999999,
            1e-300,
            -3.7e17,
            f64::MIN_POSITIVE,
        ];
        let mut st = 0x2545_F491_4F6C_DD1Du64;
        for _ in 0..2000 {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            // Any finite normal double; the exponent is squeezed into the normal range.
            let bits = (st & !(0x7ffu64 << 52)) | ((1 + (st >> 53) % 2045) << 52);
            let x = f64::from_bits(bits);
            if x.is_normal() {
                cases.push(x);
            }
        }
        for x in cases {
            let back = to_f64_exact(&BigFloat::from_f64(x, 128));
            assert_eq!(back, Some(x), "{x:e} did not survive the round trip");
        }
        assert_eq!(to_f64_exact(&BigFloat::from_f64(0.0, 128)), Some(0.0));
    }

    /// **The rounding check has to be able to say no.**
    ///
    /// [`round_to_f64`] answers `None` when the interval straddles a rounding boundary, and that
    /// branch is the whole reason the function is not just "round the midpoint". A check that can
    /// only say yes is indistinguishable from no check — and there is a specific way to build one
    /// here, by forming `mid ± rad` at the realization's own precision so the radius rounds away.
    /// So: a radius wide enough to be undecidable must be refused, and a tight one accepted.
    #[test]
    fn the_rounding_check_refuses_an_undecidable_interval() {
        let a = Angle::from_deg(Rat::new(37, 1).unwrap()).unwrap();
        let (c, _, rc, _) = a.cos_sin_bounded(128);
        assert!(
            round_to_f64(&c, rc, 128).is_some(),
            "a 2^-128 radius is decidable"
        );
        // An ulp-wide radius cannot be: it reaches both neighbours.
        assert!(round_to_f64(&c, Mag::pow2(-52), 128).is_none());
        // And zero is the case no precision resolves — cos 90 is exactly 0, so its interval
        // straddles zero forever. This is why `cos_sin_f64` resolves the family first.
        let a90 = Angle::from_deg(Rat::from_int(90)).unwrap();
        for prec in [128usize, 256, 512] {
            let (c90, _, r90, _) = a90.cos_sin_bounded(prec);
            assert!(
                round_to_f64(&c90, r90, prec).is_none(),
                "cos 90 became decidable at {prec}, which would make the quadrantal branch optional"
            );
        }
    }

    /// The escalation and fallback rungs, counted — over a corpus that reaches for them.
    ///
    /// Near an axis `|cos|` is tiny while its error bound is absolute, so the relative radius grows
    /// and 128 bits stops being enough; that is the only place the second rung is reachable. **The
    /// count is reported rather than asserted nonzero**: with `Rat` bounded by `i128` an angle
    /// cannot get closer to 90° than ~6e-39 degrees, so it is entirely possible that nothing in a
    /// finite corpus needs it — but a silent zero and an unreachable branch look identical, and
    /// this at least says which corpus produced the zero.
    #[test]
    fn the_escalation_rung_is_reachable() {
        let before = ROUND_ESCALATED.with_borrow(|c| *c);
        let mut asked = 0usize;
        for k in 1..400i128 {
            // Just off 90 degrees, by ever smaller amounts.
            for d in [
                10i128.pow(9),
                10i128.pow(18),
                10i128.pow(30),
                i128::MAX / 91,
            ] {
                if let Some(a) = Rat::new(90 * d + k, d).and_then(Angle::from_deg) {
                    let (c, s) = a.cos_sin_f64();
                    assert!(
                        c.is_finite() && s.is_finite(),
                        "{a:?} realized to a non-number"
                    );
                    // Whatever rung answered, the answer must still bound its own error.
                    let (dc, ds) = a.realization_error_of(c, s);
                    assert!(dc >= 0.0 && ds >= 0.0);
                    asked += 1;
                }
            }
        }
        let after = ROUND_ESCALATED.with_borrow(|c| *c);
        eprintln!(
            "[round_to_f64] {} escalated to 256, {} fell back, over {asked} near-axis angles",
            after.0 - before.0,
            after.1 - before.1
        );
        assert!(asked > 500, "corpus shrank to {asked}");
        assert_eq!(
            after.1, before.1,
            "the undecidable fallback fired, which is derived not to"
        );
    }

    /// The memo answers the second ask, keys on the angle's *value* rather than its spelling — same
    /// failure mode as [`TRIG`]'s, showing up as work done twice rather than a wrong answer — **and
    /// keys on the realized pair**, which is the part that is about correctness rather than cost.
    #[test]
    fn the_realization_error_memo_keys_on_the_angle_and_the_pair() {
        let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
        let entries = || F64_ERR.with_borrow(|m| m.len());
        // An angle no other test in this thread asks for.
        let a = deg(1234567, 9973);
        let (c, s) = a.cos_sin_f64();
        let before = entries();
        let first = a.realization_error_of(c, s);
        assert_eq!(entries(), before + 1, "the first ask must insert");
        assert_eq!(a.realization_error_of(c, s), first);
        assert_eq!(
            entries(),
            before + 1,
            "the second ask must be answered, not recomputed"
        );
        // The same angle, unreduced, with the same pair, is the same entry.
        assert_eq!(deg(2469134, 19946).realization_error_of(c, s), first);
        assert_eq!(
            entries(),
            before + 1,
            "the key is the angle, not its spelling"
        );
        // ★ A neighbouring realization of that same angle is a *different question* and must get
        // its own answer — an entry keyed on the angle alone would hand back `first`, an error
        // measured against a pair this caller is not holding.
        let nudged = f64::from_bits(c.to_bits() + 1);
        assert_ne!(a.realization_error_of(nudged, s), first);
        assert_eq!(entries(), before + 2, "the pair is part of the key");
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

        // Non-quadrantal: the same arithmetic on the same realized pair, so this pins the *route*
        // (`px + u·c − v·s`, in that order) rather than the values.
        let a = Angle::from_deg(Rat::from_int(37)).unwrap();
        let (c, s) = a.cos_sin_f64();
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
            let ec = hp_err_exp(&a.cos_sin_at(GT_PREC).0, exact_cos);
            assert!(ec < GATE_EXP, "cos {deg}° error 2^{ec} exceeds gate");
            worst = worst.max(ec);
            if let Some(s) = exact_sin {
                let es = hp_err_exp(&a.cos_sin_at(GT_PREC).1, s);
                assert!(es < GATE_EXP, "sin {deg}° error 2^{es} exceeds gate");
                worst = worst.max(es);
            }
        }
        eprintln!(
            "[H1.5] astro-float {GT_PREC}-bit worst cos/sin error at rational angles: 2^{worst}"
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

    /// **The memo memoizes, and its key is the *value* of the angle rather than its spelling.**
    ///
    /// Two failure modes, and neither shows up as a wrong answer — [`Angle::cos_sin_bounded`] is a
    /// pure function either way, so a broken memo is only slow. That is exactly why it needs a
    /// test: a rotated boolean spent 22% of its time re-running Ziv's loop for angles it had
    /// already realized, and nothing but a measurement would say so.
    ///
    /// - **It does not memoize** (a rewrite drops the lookup): asking twice would insert twice.
    /// - ★ **The key fragments**: `90/1` and `180/2` are the same angle. If they hashed apart the
    ///   memo would still be *correct* and still show a high hit rate, while paying twice for every
    ///   angle a caller happened to spell in two ways. `Angle::from_deg` reduces through
    ///   `Rat::new`, so they are one entry — this is the guard on that.
    ///
    /// Deltas, not absolute counts: the memo is a `thread_local`, and the harness gives each test
    /// its own thread, but nothing here should depend on which tests ran first.
    #[test]
    fn the_trig_memo_keys_on_the_angle_not_its_spelling() {
        let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
        // A precision no other test asks for, so this thread's map cannot be pre-warmed for it.
        let prec = 704;
        let before = trig_entries();
        let first = deg(37, 1).cos_sin_bounded(prec);
        assert_eq!(trig_entries(), before + 1, "the first ask must insert");

        let again = deg(37, 1).cos_sin_bounded(prec);
        assert_eq!(
            trig_entries(),
            before + 1,
            "the second ask must be answered, not recomputed"
        );
        assert_eq!(
            (first.0.clone(), first.2),
            (again.0.clone(), again.2),
            "and answered with the same value"
        );

        // ★ The same angle, spelled as an unreduced ratio, is the same entry.
        let spelled = deg(74, 2).cos_sin_bounded(prec);
        assert_eq!(
            trig_entries(),
            before + 1,
            "74/2 is 37/1: a second entry means the key is the spelling, not the angle"
        );
        assert_eq!((first.0, first.2), (spelled.0, spelled.2));

        // …and precision *is* part of the key: a different depth is a different answer, so reusing
        // an entry across depths would hand back coordinates realized at the wrong one.
        let deeper = deg(37, 1).cos_sin_bounded(prec + 64);
        assert_eq!(trig_entries(), before + 2, "precision must key the memo");
        assert_ne!(
            deeper.2, again.2,
            "a deeper realization has a smaller bound"
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
                    let observed = Mag::pow2(de as i64);
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

    /// Squared lengths a frame normal actually produces, plus awkward ones.
    ///
    /// `(0,0,1)` and `(1,1,0)` give `1` and `2`; a profile edge `(1,2)` gives a wall normal
    /// `(2,−1,0)` and so `5`; the fractions are what a normal reduced by its own content leaves;
    /// and the last two exercise the `i128 → f64` term at and past 2⁵³.
    const INV_SQRT_CASES: [(i128, i128); 12] = [
        (1, 1),
        (2, 1),
        (3, 1),
        (5, 1),
        (4, 1),
        (9, 25),
        (1, 2),
        (13, 7),
        (1_000_003, 3),
        (1, 1_000_000),
        (1, 1i128 << 60),
        (1, (1i128 << 60) + 1),
    ];

    #[test]
    fn the_inverse_sqrt_bound_holds_against_a_far_deeper_realization() {
        let mut worst = std::collections::BTreeMap::<usize, (i64, String)>::new();
        for prec in [128usize, 256, 512] {
            for (num, den) in INV_SQRT_CASES {
                let v = Rat::new(num, den).unwrap();
                let (z, bound) = inv_sqrt_bounded(v, prec).unwrap();
                let deep = prec + 512;
                let (rz, _) = inv_sqrt_bounded(v, deep).unwrap();
                let diff = z.sub(&rz, deep, HP_RM);
                let Some(de) = (if diff.is_zero() {
                    None
                } else {
                    diff.exponent()
                }) else {
                    continue; // exactly equal — nothing to bound
                };
                // `|diff| < 2^de`; the bound must be at least that.
                assert!(
                    !bound.lt(Mag::pow2(de as i64)),
                    "1/sqrt({num}/{den}) at {prec} bits: error 2^{de} exceeds its bound 2^{:?}",
                    bound.exp2()
                );
                // ★★ Slack is only meaningful where the observation actually *measures* this
                // bound. A `prec`-bit realization of a value of magnitude `2^zexp` must carry an
                // error near `2^(zexp − prec)`; when the observed gap is far below that, the two
                // realizations agreed better than either one's own accuracy warrants — measured,
                // `1/√(1/(2⁶⁰+1))` at 128 bits agrees to 154 bits where 98 is all that is earned,
                // and `1/√(1/10⁶)` is exactly `1000` at every precision. Those samples understate
                // the true error, so counting them would read a *lucky observation* as a loose
                // bound. Same shape as the trig test's reference-noise filter.
                let earned = z.exponent().unwrap_or(0) as i64 - prec as i64 - 8;
                if let (Some(be), true) = (bound.exp2(), (de as i64) >= earned) {
                    let slack = be - de as i64;
                    let e = worst.entry(prec).or_insert((i64::MIN, String::new()));
                    if slack > e.0 {
                        *e = (slack, format!("{num}/{den}, bound 2^{be} vs error 2^{de}"));
                    }
                }
            }
        }
        for (prec, (slack, at)) in &worst {
            eprintln!("[cip] {prec}-bit inverse-sqrt bound slack: 2^{slack}  ({at})");
        }
        for (prec, (slack, at)) in &worst {
            assert!(
                *slack <= 8,
                "at {prec} bits the bound is 2^{slack} above the worst observed error ({at}) — \
                 that is a fudge factor wearing a derivation's clothes, not a tight bound"
            );
        }
    }

    /// ★ The claim the axis-aligned corpus rests on: a frame whose normal is a coordinate
    /// direction realizes its `1/|n|` **exactly**, so it never reaches arbitrary precision and
    /// carries no realization error at all.
    #[test]
    fn an_axis_aligned_frame_needs_no_arbitrary_precision() {
        // `n·n` for (0,0,1) and (2,0,0) — an axis-aligned face and a wall — then the Pythagorean
        // ones, where the answer is exactly rational but (`5/3`) need not be an f64.
        for (v, want) in [
            ((1, 1), 1.0),
            ((4, 1), 0.5),
            ((25, 1), 0.2),
            ((9, 25), 5.0 / 3.0),
        ] {
            let r = Rat::new(v.0, v.1).unwrap();
            assert_eq!(
                inv_sqrt_exact(r).map(Rat::to_f64),
                Some(want),
                "1/sqrt{v:?}"
            );
            assert_eq!(inv_sqrt_f64(r), Some(want));
        }
        // …and a tilted one is honestly irrational, so it falls to the ladder.
        for v in [(2, 1), (3, 1), (5, 1), (1, 2)] {
            assert_eq!(inv_sqrt_exact(Rat::new(v.0, v.1).unwrap()), None, "{v:?}");
        }
        // There is no frame with a non-positive squared length; that is a broken premise.
        for v in [(0, 1), (-1, 1)] {
            let r = Rat::new(v.0, v.1).unwrap();
            assert_eq!(inv_sqrt_exact(r), None);
            assert_eq!(inv_sqrt_f64(r), None);
            assert_eq!(inv_sqrt_bounded(r, 128), None);
        }
    }

    /// **Correctly rounded, which is what keeps debug and release identical** — the failure this
    /// crate already paid for once in the trig path. Checked against a realization deep enough
    /// that its own rounding cannot reach the 53rd bit.
    #[test]
    fn inv_sqrt_f64_lands_on_the_nearest_f64() {
        // Deltas, not absolute counts: the memo is a `thread_local` and the harness may have run
        // other tests on this thread first.
        let before = INV_SQRT_ESCALATED.with_borrow(|c| *c);
        for (num, den) in INV_SQRT_CASES {
            let v = Rat::new(num, den).unwrap();
            let (deep, _) = inv_sqrt_bounded(v, 1024).unwrap();
            let want = to_f64_exact(&deep).unwrap();
            assert_eq!(
                inv_sqrt_f64(v),
                Some(want),
                "1/sqrt({num}/{den}) is not the nearest f64"
            );
        }
        let after = INV_SQRT_ESCALATED.with_borrow(|c| *c);
        let (deeper, undecided) = (after.0 - before.0, after.1 - before.1);
        eprintln!("[cip] inverse-sqrt escalations: {deeper} deeper, {undecided} undecided");
        // ★ 128 bits sufficed for every case here, and the second rung is a proven cap rather
        // than a rung anything reaches: `1/√v` is never zero for a positive `v`, so its interval
        // cannot straddle zero the way `cos 90°`'s does and the rounding always decides.
        assert_eq!(
            (deeper, undecided),
            (0, 0),
            "a case escalated — the ladder's first rung no longer covers the population"
        );
    }

    /// **The reported realization error really does cover the error that is there** — checked
    /// against a realization far deeper than the one the measurement itself uses (128 bits), so
    /// the ruler and the thing being measured are not the same instrument.
    #[test]
    fn the_inverse_sqrt_realization_error_covers_the_error_that_is_there() {
        let mut worst = (i32::MIN, String::new());
        for (num, den) in INV_SQRT_CASES {
            let v = Rat::new(num, den).unwrap();
            let f = inv_sqrt_f64(v).unwrap();
            let reported = inv_sqrt_error_of(v, f).unwrap();
            let (deep, deep_rad) = inv_sqrt_bounded(v, 1024).unwrap();
            let true_err = hp_err_exp(&deep, f);
            // A `BigFloat` is `m · 2^e` with `m ∈ [0.5, 1)`, so the exponent gives
            // `2^(e−1) ≤ |f − deep| < 2^e` — the *lower* end is what a bound has to clear. Using
            // `2^e` would demand the reported value exceed an over-estimate.
            //
            // ★ **And the reference is not exact either**, so its own radius comes off:
            // `|f − true| ≥ |f − deep| − rad`. Without that, `1/√(1/10⁶)` — which is exactly
            // `1000.0`, correctly reported as a zero error — fails against the 1024-bit
            // realization's `2⁻¹⁰¹³` residue, and the test would be measuring its own ruler.
            let ref_rad = deep_rad.exp2().map_or(0.0, |e| 2f64.powi(e as i32));
            let floor = if true_err == i32::MIN {
                0.0
            } else {
                (2f64.powi(true_err - 1) - ref_rad).max(0.0)
            };
            assert!(
                reported >= floor,
                "1/sqrt({num}/{den}): reported {reported:e} is below the true error 2^{true_err}"
            );
            if true_err > worst.0 {
                worst = (true_err, format!("{num}/{den}"));
            }
        }
        // ★★★★ **Exact is not the same as dyadic, and only dyadic earns the zero.** An
        // axis-aligned frame lands on `1`/`½` and must report a literal `0.0` — that is what keeps
        // it on the exact predicate path. A Pythagorean one lands on `1/5` or `5/3`, which are
        // exactly rational and still not f64, so they must report the rounding that really
        // happened. Reporting `0` there (the shape the trig path can safely use) is the bug this
        // pair of assertions exists to keep out.
        // `1/9` is here rather than below because `1/√(1/9)` is `3` — an integer is dyadic too.
        for (num, den) in [(1i128, 1i128), (4, 1), (1, 4), (1, 64), (1, 9)] {
            let v = Rat::new(num, den).unwrap();
            let f = inv_sqrt_f64(v).unwrap();
            assert_eq!(inv_sqrt_error_of(v, f), Some(0.0), "1/sqrt({num}/{den})");
        }
        // `(3,4,0)` gives `n·n = 25` and `1/|n| = 1/5`; `(3,4,0)` scaled gives `9/25` and `5/3`.
        for (num, den) in [(25i128, 1i128), (9, 25), (49, 1)] {
            let v = Rat::new(num, den).unwrap();
            let f = inv_sqrt_f64(v).unwrap();
            assert!(
                inv_sqrt_exact(v).is_some() && inv_sqrt_error_of(v, f).unwrap() > 0.0,
                "1/sqrt({num}/{den}) is exactly rational but not an f64 — it must be charged"
            );
        }
        eprintln!(
            "[cip] worst inverse-sqrt f64 realization error: 2^{} ({})",
            worst.0, worst.1
        );
    }

    /// ★★ **The naive f64 route is not good enough, measured** — `1.0 / v.to_f64().sqrt()` is
    /// three roundings and lands on the wrong `f64` often enough to see. Without a case that
    /// actually differs, the exact realization above would be a cost with no effect, and this
    /// test would be indistinguishable from one that is not running.
    #[test]
    fn the_naive_f64_route_gets_the_last_bit_wrong() {
        let (mut n, mut differ) = (0usize, 0usize);
        for num in 1i128..400 {
            for den in 1i128..7 {
                let v = Rat::new(num, den).unwrap();
                let naive = 1.0 / v.to_f64().sqrt();
                n += 1;
                if inv_sqrt_f64(v) != Some(naive) {
                    differ += 1;
                }
            }
        }
        eprintln!("[cip] naive 1/sqrt disagrees with the correctly rounded one: {differ}/{n}");
        assert!(
            differ > 0,
            "the naive route agreed everywhere on {n} cases — either the correctly rounded path \
             is not running, or this population cannot tell them apart"
        );
    }

    /// Error of a high-precision realization `a` against an exact f64 `b`, as a
    /// power-of-two exponent (`i32::MIN` when exactly equal).
    fn hp_err_exp(a: &BigFloat, b: f64) -> i32 {
        let err = a.sub(&BigFloat::from_f64(b, GT_PREC), GT_PREC, HP_RM);
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

#[cfg(test)]
mod symprobe {
    use crate::{Angle, Rat};
    #[test]
    #[ignore = "probe"]
    fn mirror_pairs() {
        let a = |d: i128| Angle::from_deg(Rat::from_int(d)).unwrap();
        for (x, y) in [(72i128, 288i128), (9, 351), (117, 243), (45, 315)] {
            let (cx, sx) = a(x).cos_sin_f64();
            let (cy, sy) = a(y).cos_sin_f64();
            println!(
                "[sym] {x} vs {y}: cos equal {}  sin negated {}   ({cx:.20e} / {cy:.20e})",
                cx == cy,
                sx == -sy
            );
        }
    }
}
