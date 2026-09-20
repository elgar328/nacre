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
pub mod bounded;
pub mod mag;
pub mod quad;
pub use bounded::{Bounded, HpBounded, rat_to_big};
pub use mag::Mag;
pub use quad::{
    QuadVal, biquad_sign, cylinder_radial_side, cylinders_clear, cylinders_nested,
    exceeds_root_sum, segment_meets_cylinder,
};

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
///
/// Public because the judge (`nacre-cip`) rounds with it too — it used to carry an identical
/// private copy, and one mode in two places is one more thing that can drift.
pub const HP_RM: RoundingMode = RoundingMode::ToEven;

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
    static INV_SQRT: RefCell<HashMap<(Rat, usize), HpBounded>> =
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
type TrigAt = (HpBounded, HpBounded);

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
/// ★★★ **Built in `BigFloat`, not through `f64`.** An earlier spelling wrote
/// `BigFloat::from_f64(2f64.powi(e))` and refused `|e| > 1000` — safe, because flushing a radius
/// to zero would read as a *tighter* interval than the realization earned. But a realization at
/// 1024 bits carries a radius near `2⁻¹⁰¹⁹`, so that refusal turned **more** precision into
/// "undecided": climbing a ladder made the answer worse, measured. A power of two is exact in
/// `BigFloat` at any precision, so the bound is assembled there and the range guard is gone.
///
/// `None` only when the exponent is absurd enough that the multiply loop would not terminate
/// usefully — far outside anything a realization produces.
fn rad_upper_big(rad: Mag, p: usize) -> Option<BigFloat> {
    let Some(e) = rad.exp2() else {
        return Some(BigFloat::from_f64(0.0, p)); // an exact realization: a zero radius is honest
    };
    if !(-1_000_000..=1_000_000).contains(&e) {
        return None;
    }
    let mut out = BigFloat::from_f64(1.0, p);
    let mut rem = e;
    while rem != 0 {
        let step = rem.clamp(-1000, 1000);
        out = out.mul(&BigFloat::from_f64(2f64.powi(step as i32), p), p, HP_RM);
        rem -= step;
    }
    Some(out)
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

// ---------------------------------------------------------------------------------------------
// Decimal places — the twin of `round_to_f64`, and the one long division both arms share
// ---------------------------------------------------------------------------------------------

/// The precision `x` actually carries, in bits — the mantissa the words spell.
fn bits_of(x: &BigFloat) -> usize {
    x.as_raw_parts()
        .map_or(0, |(w, _, _, _, _)| w.len() * astro_float::WORD_BIT_SIZE)
}

/// A `BigFloat` as the exact `(numerator, positive denominator)` it *is*.
///
/// A binary float is a dyadic rational, so this loses nothing: astro-float stores `m · 2^e` with
/// the mantissa normalized into `[0.5, 1)`, which makes the integer the words spell equal to
/// `m · 2^bits` and the value `± words · 2^(e - bits)`.
///
/// ★ **The scale comes from the words being read, not from `as_raw_parts`' second field.**
/// Measured (6 precisions x 6 values, including 53 and 100): astro-float pads the mantissa to
/// whole words and reports that padded length, so the two agree everywhere. Deriving it from the
/// same words the integer is assembled from keeps them agreeing by construction rather than by
/// coincidence.
fn big_to_ratio(x: &BigFloat) -> Option<(num_bigint::BigInt, num_bigint::BigInt)> {
    use num_bigint::BigInt;
    if x.is_zero() {
        return Some((BigInt::from(0), BigInt::from(1)));
    }
    let (words, _bits, sign, e, _inexact) = x.as_raw_parts()?;
    const WB: usize = astro_float::WORD_BIT_SIZE;
    let mut m = BigInt::from(0);
    for w in words.iter().rev() {
        m = (m << WB) | BigInt::from(*w);
    }
    if sign == astro_float::Sign::Neg {
        m = -m;
    }
    let shift = i64::from(e) - (words.len() * WB) as i64;
    Some(match shift >= 0 {
        true => (m << shift as usize, BigInt::from(1)),
        false => (m, BigInt::from(1) << (-shift) as usize),
    })
}

/// Round `q + r/den` to an integer, **ties to even** — IEEE's rule, and the kernel's.
///
/// ⚠★★★ **Not "away from zero", which is what this crate shipped first.** At `2⁵² + ½` both
/// neighbours are representable, so the tie rule decides, and half-away answered
/// `4503599627370497` where [`Rat::to_f64`] — the road the point cache is built on — answers
/// `4503599627370496`. Two roads for one quantity disagreeing at a tie is exactly what the
/// vertex door exists to rule out. A 20,000-pair sweep missed it: random pairs are never ties,
/// and the boundary family it was checked against (`2ᵏ⁺¹ − 1` over 2) lands on the one tie both
/// rules resolve the same way.
///
/// Both roundings this crate hands out — binary ([`round_shifted`]) and decimal
/// ([`round_scaled`]) — go through here, so they cannot drift apart.
fn round_ties_even(
    q: num_bigint::BigInt,
    r: &num_bigint::BigInt,
    den: &num_bigint::BigInt,
) -> num_bigint::BigInt {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;
    let twice = (r * BigInt::from(2)).magnitude().clone();
    let step = match twice.cmp(den.magnitude()) {
        core::cmp::Ordering::Less => false,
        core::cmp::Ordering::Greater => true,
        // The exact tie: move only if it would land on an even integer.
        core::cmp::Ordering::Equal => q.is_odd(),
    };
    if !step {
        return q;
    }
    match r.sign() {
        Sign::Minus => q - 1,
        _ => q + 1,
    }
}

/// The nearest integer to `num · 10^places / den`, **ties to even**. `den` must be positive.
///
/// Exact: every step is integer arithmetic, so this is the *decision* about the last digit and
/// not an approximation of it.
fn round_scaled(
    num: &num_bigint::BigInt,
    den: &num_bigint::BigInt,
    places: usize,
) -> num_bigint::BigInt {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let scaled = num * BigInt::from(10).pow(places as u32);
    // `div_rem` truncates toward zero and `r` carries the dividend's sign, so the tie decision
    // works on magnitudes and steps away from zero.
    let (q, r) = scaled.div_rem(den);
    round_ties_even(q, &r, den)
}

/// The nearest integer to `num · 2^k / den` (**ties to even**), `den` positive, `k` any sign.
///
/// [`round_scaled`]'s binary twin — same decision through [`round_ties_even`], a different radix,
/// so the two roundings this crate hands out cannot drift apart in their tie rule.
fn round_shifted(num: &num_bigint::BigInt, den: &num_bigint::BigInt, k: i64) -> num_bigint::BigInt {
    use num_integer::Integer;
    let (n, d) = match k >= 0 {
        true => (num << k as usize, den.clone()),
        false => (num.clone(), den << (-k) as usize),
    };
    let (q, r) = n.div_rem(&d);
    round_ties_even(q, &r, &d)
}

/// **The nearest `f64` to the exact rational `num/den`** (`den` positive) — one rounding, from
/// integers.
///
/// [`nearest_f64`]'s unbounded twin: that one is the `Rat` road and caps at `u128`, this one takes
/// the `(numerator, common denominator)` pair [`MeetPoint::lift`] hands out, which is how a
/// coordinate too wide for `Rat` still gets a correctly rounded coordinate rather than a refusal.
///
/// `None` when the value is outside `f64`'s normal range.
pub fn nearest_f64_big(num: &num_bigint::BigInt, den: &num_bigint::BigInt) -> Option<f64> {
    nearest_f64_big_exact(num, den).map(|(v, _)| v)
}

/// [`nearest_f64_big`], **and whether the rounding lost anything**.
///
/// ★★★ **The second half is not optional information.** A three-plane meet is an exact rational,
/// which is a fact about the *realization*, not about the `f64` it is then read out as: a 59-bit
/// coordinate does not fit a 53-bit mantissa, so the readout rounds and the value handed over is
/// **not** the point. Reporting that as a zero error is the same lie the cache tells — measured on
/// the tilted-frame family, where the exact value is `0.130864196953086372` and its `f64` is
/// `0.13086419695308637578…`.
///
/// `true` means the `f64` **is** the rational, so a caller may honestly say its error is zero.
pub fn nearest_f64_big_exact(
    num: &num_bigint::BigInt,
    den: &num_bigint::BigInt,
) -> Option<(f64, bool)> {
    use num_bigint::BigInt;
    if num.sign() == num_bigint::Sign::NoSign {
        return Some((0.0, true));
    }
    // `value ∈ (2^(e-1), 2^(e+1))`, so `52 - e` aims the mantissa at 53 bits and lands one short
    // at worst — one nudge, never two.
    //
    // ★ **A carry out of the top bit needs no correction, measured.** Rounding can push `q` to
    // exactly `2^53` (from `2^53 - ½`, the only way), and `2^53 · 2^-k` is representable exactly,
    // so the wide `q` still spells the right `f64`. A branch for it was written, planted against a
    // 20,000-pair sweep plus a power-of-two boundary family, and never changed an answer.
    let e = num.magnitude().bits() as i64 - den.magnitude().bits() as i64;
    let mut k = 52 - e;
    let mut q = round_shifted(num, den, k);
    if q.magnitude().bits() < 53 {
        k += 1;
        q = round_shifted(num, den, k);
    }
    let m = i64::try_from(&q).ok()? as f64;
    if !(-1200..=1200).contains(&k) {
        return None;
    }
    // ⚠★★★ **Scaled in steps, because `2f64.powi` underflows before the product does.** `m` is
    // ~2⁵³, so `m · 2⁻ᵏ` can be an ordinary `f64` while `2⁻ᵏ` alone is zero: measured, `2⁻¹⁰⁰⁰` —
    // a normal `f64` at 9.33e-302 — came back **0.0**, and so did every subnormal. A silent wrong
    // answer, not a refusal. Powers of two are exact, so splitting the scale costs nothing and
    // the value underflows only where it genuinely should.
    //
    // ★ Third place in this cell where a `2f64.powi` of a realization-sized exponent was wrong
    // (`rad_upper_big`, `Realized::to_f64`'s error readout, here).
    let mut v = m;
    let mut rem = -k;
    while rem != 0 && v != 0.0 && v.is_finite() {
        let step = rem.clamp(-1000, 1000);
        v *= 2f64.powi(step as i32);
        rem -= step;
    }
    // Exact iff the scaled division left no remainder: `value = q · 2^-k` exactly means
    // `num · 2^k == q · den`. Both sides are integers, so this is a decision and not an estimate.
    let (lhs, rhs) = match k >= 0 {
        true => (num << k as usize, &q * den),
        false => (num.clone(), (&q * den) << (-k) as usize),
    };
    // ⚠★★★★ **A flushed value is not an exact one.** `exact` is decided on `q · 2^-k`, but `v` is
    // what the caller gets, and below the smallest subnormal the scaling loop flushes `v` to `0.0`
    // while the rational is nowhere near zero. That came back as `(0.0, true)` — and
    // [`Realized::to_f64`] turns the flag into `Mag::ZERO`, a *proven* bound — so the cache would
    // have published "this coordinate is exactly zero" for `2⁻¹¹⁰⁰`.
    //
    // ★ Only the claim goes, not the value: `0.0` **is** the nearest `f64` there, and
    // `a_tiny_rational_still_names_its_f64` locks that deliberately ("genuinely below the range").
    // Refusing would have lost a correct answer and broken that lock.
    let exact = lhs == rhs && !(v == 0.0 && num.sign() != num_bigint::Sign::NoSign);
    let _ = BigInt::from(0);
    v.is_finite().then_some((v, exact))
}

/// `places` decimal places of the **exact** rational `num/den` (`den` positive).
///
/// ★ For a value that really is rational — a three-plane meet, a `QuadVal` whose radical vanishes
/// — every digit this prints is a digit of the coordinate itself, not of an approximation to it.
/// That is the thing an exact kernel can say and a floating-point one cannot.
pub fn decimals_of_ratio(
    num: &num_bigint::BigInt,
    den: &num_bigint::BigInt,
    places: usize,
) -> String {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let n = round_scaled(num, den, places);
    let sign = if n.sign() == num_bigint::Sign::Minus {
        "-"
    } else {
        ""
    };
    let a = n.magnitude();
    if places == 0 {
        return format!("{sign}{a}");
    }
    let (int, frac) = BigInt::from(a.clone()).div_rem(&BigInt::from(10).pow(places as u32));
    format!("{sign}{int}.{:0>width$}", frac.magnitude(), width = places)
}

/// `places` decimal places of the realization `mid ± rad`, or `None` when `prec` bits **do not
/// decide them** — the twin of [`round_to_f64`], and `None` means the same thing there.
///
/// ★ The caller's move on `None` is to realize again at higher precision, exactly as
/// `realize_inv_sqrt_rounded` escalates. That is why this reports undecided rather than picking:
/// a digit invented here would be indistinguishable, to everything downstream, from one the
/// definition actually determines.
pub fn round_to_digits(mid: &BigFloat, rad: Mag, places: usize) -> Option<String> {
    if mid.is_nan() || mid.is_inf() {
        return None;
    }
    // ★★★ **The working precision is `mid`'s own, not a parameter.** An earlier spelling took one,
    // and a caller that climbed a ladder and then passed the bottom rung would widen `mid` back
    // down to it — printing digits of a 192-bit rounding as if they were the value's. With a zero
    // radius that is silent (the interval ends agree, because they are the same rounded number),
    // which makes it worse than a refusal. Reading the precision off the value it is about is the
    // one spelling that cannot disagree with itself.
    let p = bits_of(mid) + 64;
    let r = rad_upper_big(rad, p)?;
    let (lo, hi) = (mid.sub(&r, p, HP_RM), mid.add(&r, p, HP_RM));
    let (nlo, dlo) = big_to_ratio(&lo)?;
    let (nhi, dhi) = big_to_ratio(&hi)?;
    let (slo, shi) = (
        round_scaled(&nlo, &dlo, places),
        round_scaled(&nhi, &dhi, places),
    );
    // Both ends rounding to the same scaled integer is what "these digits are determined" means.
    (slo == shi).then(|| decimals_of_ratio(&nlo, &dlo, places))
}

// ---------------------------------------------------------------------------------------------
// Realizing an algebraic coordinate at a precision — the arithmetic behind a curved vertex
// ---------------------------------------------------------------------------------------------
//
// The arithmetic itself is [`HpBounded`]'s (`bounded.rs`) — what lives here are the realizations
// that drive it. A second spelling of that arithmetic used to sit here (a tuple alias and five
// free functions); its magnitude reader charged an exact zero a rounding and its rational entry
// truncated below 128 bits, both of which the one spelling does not.

/// **`√v` realized at `p` bits, with its error** — `v · (1/√v)`, so the one radical primitive this
/// crate already has ([`inv_sqrt_bounded`]) is the only place a square root is approached.
///
/// Crate-private: its consumer is [`realize_quad`] (measured); a radius stated as its square is
/// wide and takes [`sqrt_bounded_big`]. `None` for a negative `v`.
pub(crate) fn sqrt_bounded(v: Rat, p: usize) -> Option<HpBounded> {
    if v < Rat::from_int(0) {
        return None;
    }
    if v == Rat::from_int(0) {
        return Some(HpBounded::exact(BigFloat::from_f64(0.0, p)));
    }
    let inv = inv_sqrt_bounded(v, p)?;
    Some(HpBounded::of_rat(v, p).mul(&inv, p))
}

/// [`sqrt_bounded`] for a wide radicand — `√(n/d) = n · (1/√(n·d))`, so the one `BigInt` radical
/// primitive ([`inv_sqrt_bigint_bounded`]) is the only place a root is approached.
pub fn sqrt_bounded_big(v: &BigRat, p: usize) -> Option<HpBounded> {
    if v.is_negative() {
        return None;
    }
    if v.is_zero() {
        return Some(HpBounded::exact(BigFloat::from_f64(0.0, p)));
    }
    let nd = v.numer() * v.denom();
    let inv = inv_sqrt_bigint_bounded(&nd, p)?;
    Some(HpBounded::of_bigint(v.numer(), p).mul(&inv, p))
}

/// **A quadratic algebraic scalar `a + b√c` realized at `p` bits, with its error.**
///
/// ★ Exact when the radical vanishes or resolves ([`quad::QuadVal::as_rat`]) — the value is asked,
/// not its provenance, so a tangency's rational root takes the rational road even though it
/// arrived through the same variant as an irrational one.
pub fn realize_quad(q: &quad::QuadVal, p: usize) -> Option<HpBounded> {
    if let Some(r) = q.as_rat() {
        // The rational is the value; the only error is this realization's own rounding — charged
        // unconditionally rather than through `HpBounded::of_rat`'s exact branch, so that this
        // arm never hands an approached coordinate a zero radius (its consumer reads a zero as the
        // exact arm's answer).
        let v = rat_to_big(r, p);
        return Some(HpBounded::new(v.clone(), HpBounded::round_off(&v, p)));
    }
    let root = sqrt_bounded(q.c(), p)?;
    let term = HpBounded::of_rat(q.b(), p).mul(&root, p);
    Some(HpBounded::of_rat(q.a(), p).add(&term, p))
}

/// `base + dir·s` for exact rational `base`/`dir` and a realized `s` — the last step of a point
/// that lives at a parameter along an exactly-stated line.
pub fn affine_bounded(base: Rat, dir: Rat, s: &HpBounded, p: usize) -> Option<HpBounded> {
    let term = HpBounded::of_rat(dir, p).mul(s, p);
    Some(HpBounded::of_rat(base, p).add(&term, p))
}

/// **A point on a circle's `+ref` seam, realized at `p` bits.**
///
/// `centre + r · e₁/|e₁|`, where `e₁` is the reference direction's component perpendicular to the
/// axis and the radius is stated as its square `r2` (the cylinder truth's form). `centre`, `e₁`
/// and `r2` are exact rationals. When `r2` is a rational's square — every radius a user writes —
/// the single irrational step is `1/|e₁|`, one `inv_sqrt_bounded` with exact arithmetic around it,
/// the road this always took; otherwise `√r2` is realized beside it.
pub fn realize_seam_point(
    centre: [Rat; 3],
    e1: [Rat; 3],
    r2: &BigRat,
    p: usize,
) -> Option<[HpBounded; 3]> {
    let mut sq = Rat::from_int(0);
    for c in &e1 {
        sq = sq.checked_add(c.checked_mul(*c)?)?;
    }
    let inv = inv_sqrt_bounded(sq, p)?;
    let scale = match rat_sqrt_exact_big(r2) {
        Some(radius) => HpBounded::of_rat(radius, p).mul(&inv, p),
        None => sqrt_bounded_big(r2, p)?.mul(&inv, p),
    };
    let coord = |k: usize| {
        let radial = HpBounded::of_rat(e1[k], p).mul(&scale, p);
        HpBounded::of_rat(centre[k], p).add(&radial, p)
    };
    Some([coord(0), coord(1), coord(2)])
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
/// missing name before — a wide name does not open a frame.
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
    ///
    /// ⚠ **Primitive over *four* coefficients — the normal alone is not**.
    /// The content divided out is `gcd(a,b,c,d)`, so `(a,b,c)` keeps a factor of
    /// `gcd(a,b,c)/gcd(a,b,c,d)`; for an axis-aligned plane at an offset needing a long decimal
    /// that factor is the offset's **denominator**, arbitrarily large on the plainest of planes.
    /// A consumer that reads only the normal **and multiplies** must divide it by its own gcd
    /// first, or it pays that factor — squared, if it takes a cross product.
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
/// needed** — what lets a datum through [`MeetPoint::Wide`] vertices keep a name.
/// The same lift-and-join [`plane_name_big`] performs, with each point's
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

// ---------------------------------------------------------------------------
// Total rational-vector predicates.
//
// A question whose answer is a **sign** has no width: the answer is one of three values, and
// only the road to it can overflow. These clear the denominators once and run in `BigInt`, so
// they cannot decline — which is what lets their callers say `None`/reject for the geometry
// alone. The same move `three_planes_rat` makes for its value ("the caller resolves the
// conflation by asking the integer core") and `plane_name_exact` for its name.
//
// ★ **The lift keeps each vector's denominator.** Clearing denominators multiplies a vector by
// a positive factor, which is harmless for an expression that is *homogeneous* in that vector
// (a cross, a dot) and **wrong** for one that is not — an expression mixing a term in `p` with
// a constant would state a different proposition after scaling `p` alone. Every predicate here
// says which case it is in, so the ones that need the scale can multiply it back in.
// ---------------------------------------------------------------------------

/// A rational 3-vector as **(integer components, positive denominator)**: `v = out.0 / out.1`.
fn lift3(v: &[Rat; 3]) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let den: [BigInt; 3] = core::array::from_fn(|i| BigInt::from(v[i].denom()));
    let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
    let num = core::array::from_fn(|i| BigInt::from(v[i].numer()) * (&d / &den[i]));
    (num, d)
}

/// **Whether two rational directions are parallel** — `a × b = 0`, exactly and always.
///
/// The cross is homogeneous in each argument, so each vector's own denominator is a positive
/// factor of the result and clearing it cannot move the zero test.
///
/// **A zero vector is parallel to everything**, this predicate included: `0 × b = 0`. Callers
/// that mean "these two span a plane" therefore get the answer they want — a zero direction
/// spans nothing — without a separate zero check.
pub fn parallel_rat(a: &[Rat; 3], b: &[Rat; 3]) -> bool {
    use num_traits::Zero;
    let (x, _) = lift3(a);
    let (y, _) = lift3(b);
    let term = |i: usize, j: usize| &x[i] * &y[j] - &x[j] * &y[i];
    term(1, 2).is_zero() && term(2, 0).is_zero() && term(0, 1).is_zero()
}

/// **The sign of `a · b`** for rational vectors, exactly and always — [`Orient::Zero`] exactly
/// when they are perpendicular. Bilinear, so each vector's positive denominator factors out of
/// the answer and the lift's scales are dropped.
pub fn dot_sign_rat(a: &[Rat; 3], b: &[Rat; 3]) -> Orient {
    let (x, _) = lift3(a);
    let (y, _) = lift3(b);
    orient_of(big_sign(
        &(0..3).map(|i| &x[i] * &y[i]).sum::<num_bigint::BigInt>(),
    ))
}

/// **A cylinder's rational unit cross-section frame** — `û₁ = ref_dir/‖ref_dir‖`,
/// `û₂ = (dir × ref_dir)/(‖dir‖‖ref_dir‖)`, both exact rationals ⊥ the axis and to each other.
///
/// `None` — and the three causes are told apart on purpose:
/// * `ref_dir` is **parallel to the axis** (a zero vector included — [`parallel_rat`] says so):
///   the statement pins no cross-section at all;
/// * `ref_dir` is **not perpendicular** to the axis ([`dot_sign_rat`], exact and total). The
///   general recipe would project it, and that is deliberately **not** done here: a caller who
///   needs the frame of a `CylinderDef` gets a frame whose `û·dir = 0` is a *fact*, not a
///   derivation, because everything downstream (`quad::cylinder_radial_side` reduces to
///   `‖dir‖²·r²·(û·û − 1)`) is sound only when that holds exactly.
/// * a norm is **irrational**, or the arithmetic left `i128` ([`inv_sqrt_exact`]).
///
/// ★ Every cylinder the modelling road states satisfies the first two by construction — a sketch
/// frame is checked orthonormal exactly, and `ref_dir` is either its `x̂` or a rim chord divided
/// by its own radius — so `None` there means the third cause alone.
pub fn cyl_unit_frame(dir: &[Rat; 3], ref_dir: &[Rat; 3]) -> Option<([Rat; 3], [Rat; 3])> {
    if parallel_rat(ref_dir, dir) || dot_sign_rat(ref_dir, dir) != Orient::Zero {
        return None;
    }
    let dot = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<Rat> {
        let mut acc = Rat::from_int(0);
        for k in 0..3 {
            acc = acc.checked_add(a[k].checked_mul(b[k])?)?;
        }
        Some(acc)
    };
    let scale = |v: &[Rat; 3], k: Rat| -> Option<[Rat; 3]> {
        Some([
            v[0].checked_mul(k)?,
            v[1].checked_mul(k)?,
            v[2].checked_mul(k)?,
        ])
    };
    let cross = [
        dir[1]
            .checked_mul(ref_dir[2])?
            .checked_sub(dir[2].checked_mul(ref_dir[1])?)?,
        dir[2]
            .checked_mul(ref_dir[0])?
            .checked_sub(dir[0].checked_mul(ref_dir[2])?)?,
        dir[0]
            .checked_mul(ref_dir[1])?
            .checked_sub(dir[1].checked_mul(ref_dir[0])?)?,
    ];
    let inv_e = inv_sqrt_exact(dot(ref_dir, ref_dir)?)?;
    let inv_m = inv_sqrt_exact(dot(dir, dir)?)?;
    Some((
        scale(ref_dir, inv_e)?,
        scale(&cross, inv_e.checked_mul(inv_m)?)?,
    ))
}

/// **How a point's distance from a plane compares with `r`**, the radius stated as its square
/// `r2` — [`Orient::Negative`] inside the slab of half-width `r` about the plane, [`Orient::Zero`]
/// exactly at distance `r`, [`Orient::Positive`] clear of it. Exact and total.
///
/// `sign((n·p + d)² − r²|n|²)`, which is `sign(dist² − r²)` scaled by the positive `|n|²` — no
/// normalization and no square root. Stated in **plane** vocabulary on purpose: the cylinder
/// gate asks it about an axis point (a wall parallel to the axis is the same distance from every
/// point of it), but nothing here is about cylinders.
///
/// ★ **Not homogeneous in `p`** — the `d` term is why — so `p`'s denominator and the radius'
/// ride into the formula instead of dropping out: for `coeffs = C/Dc`, `p = P/Dp`, `r² = R/S`
/// the answer is `sign(S(C₀₋₂·P + C₃·Dp)² − R|C₀₋₂|²Dp²)` (`Dc²` *is* a positive common factor
/// and does cancel). Dropping `Dp` instead — the obvious spelling — states a different
/// proposition, and disagrees with the truth on 1.3% of mixed-denominator inputs (measured).
///
/// **Precondition:** `r2 ≥ 0`; a negative squared radius has no distance to compare with.
pub fn point_plane_clearance_rat(coeffs: &[Rat; 4], p: &[Rat; 3], r2: &BigRat) -> Orient {
    use num_bigint::BigInt;
    debug_assert!(
        !r2.is_negative(),
        "clearance compares against a non-negative squared radius"
    );
    let (c, _dc) = lift4(coeffs);
    let (pp, dp) = lift3(p);
    let (rn, rd): (BigInt, BigInt) = (r2.numer().clone(), r2.denom().clone());
    let dot: BigInt = (0..3).map(|i| &c[i] * &pp[i]).sum::<BigInt>() + &c[3] * &dp;
    let nn: BigInt = (0..3).map(|i| &c[i] * &c[i]).sum();
    orient_of(big_sign(&(&rd * (&dot * &dot) - &rn * nn * (&dp * &dp))))
}

/// **Does `p` satisfy these coefficients exactly?** — `n·p + c = 0`, at any width.
///
/// The consumer-side net for [`cylinder_strip_side`]'s on-plane precondition. A plane has two
/// exact descriptions — its coefficients and the triangle its points span — and they need not
/// agree; a face merged into a class by *rounded* coefficients can therefore have vertices that
/// do not satisfy the class root's exact name. Asking here turns that into a fact the caller can
/// refuse on, instead of a `debug_assert` that says nothing in a release build.
pub fn point_on_plane_exact(coeffs: &[Rat; 4], p: &MeetPoint) -> bool {
    use num_bigint::BigInt;
    let (c, _dc) = lift4(coeffs);
    let (pp, dp) = p.lift();
    // (C₀₋₂·P + C₃·Dp) / (Dc·Dp) — the denominators are positive, so only the numerator decides.
    (0..3).map(|i| &c[i] * &pp[i]).sum::<BigInt>() + &c[3] * &dp == BigInt::from(0)
}

/// Which side of the **strip** a cylinder cuts out of a plane parallel to its axis a point lies
/// on — see [`cylinder_strip_side`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StripSide {
    /// Clear of the strip, on the `+(n × m)` side.
    Plus,
    /// Clear of the strip, on the `−(n × m)` side.
    Minus,
    /// Inside the strip, boundary included — the point is within `r` of the axis.
    ///
    /// For a piece with **width** this is "reaches the strip but does not span a boundary":
    /// wholly between the two rulings, or touching one at a single point.
    Inside,
    /// **The piece spans a boundary** — its interior lies strictly on both sides of one ruling.
    ///
    /// ★ A point can never be this (it has no width), so every answer
    /// [`cylinder_strip_side`] and [`cylinder_strip_side_branch`] give is one of the three above,
    /// unchanged. Only [`cylinder_strip_side_margin`] with `rho > 0` produces it.
    Crosses,
}

/// **Where a point on a plane parallel to a cylinder's axis stands relative to the strip the
/// cylinder cuts out of that plane.** Exact and total, at any width.
///
/// A plane with `n · m = 0` meets the solid cylinder in a strip of half-width `h = √(r² − d²)`
/// (empty when the plane clears, `d ≥ r`), running along the axis direction. Decomposing
/// `p − o` into three **mutually orthogonal** parts — along `n` (magnitude `d`), along `m`, and
/// along `e = n × m` — gives `dist(p, axis)² = d² + t²` where `t` is the `e` component. So with
/// `U = (p − o) · e = t·|e|` and `|e| = |n||m|`:
///
/// ```text
/// clear of the strip  ⟺  U² > (r²|n|² − (n·o + c)²) · |m|²
/// ```
///
/// No normalization and **no square root** — `h` never has to be formed.
///
/// ★ **This is one of two axes, not a rule of its own.** What the cylinder occupies in this plane
/// is a *rectangle*: this strip across, and the lateral face's axis-parameter span along
/// ([`point_axis_side`]). Clearing either axis clears the rectangle, so a caller asking whether a
/// face misses the cylinder asks both and stops at the first that separates.
///
/// ★ **`U = 0` answers [`StripSide::Inside`] before the magnitude test.** With a non-empty strip
/// that is the truth (the point sits on the axis' own in-plane line). With an empty one it is
/// merely conservative — and a caller that cares has already passed the plane through
/// [`point_plane_clearance_rat`], which is the cheaper question and the one that decides
/// emptiness. That ordering makes this function total with no precondition to forget.
///
/// **Preconditions:** `r2 ≥ 0` (the radius stated as its square); `n · m = 0` (a plane that is not
/// parallel to the axis cuts a conic, not a strip); and **`p` lies on that plane** — the decomposition takes the perpendicular
/// distance from the *axis*, so an off-plane point would be judged against a distance that is not
/// its own. All three are `debug_assert`ed.
pub fn cylinder_strip_side(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    // ★ **A door, so the two can never drift.** A point is a disk of radius zero, and at that
    // radius the margin form's answers collapse onto this one exactly: `Crosses` needs a strictly
    // positive width to be possible at all, and the remaining comparison is term-for-term the one
    // this function used to spell for itself.
    cylinder_strip_side_margin(coeffs, p, &BigRat::zero(), o, m, r2)
}

/// The three squared quantities the strip questions compare, in **one common positive scale** —
/// the whole derivation, shared by every door below.
///
/// With `e = n × m` (so `|e|² = |n|²|m|²`, the two being perpendicular):
///
/// ```text
/// U  = (p − o)·e        W² = (r²|n|² − (n·o + c)²)|m|²        ρ'² = ρ²|e|²
/// ```
///
/// A disk of radius `ρ` about `p` sweeps `U ± ρ'`, and the cylinder's two rulings stand at
/// `U = ±W`; every question below is a comparison among those three. The radii arrive **squared**
/// (`r2 = r²`, `rho2 = ρ²`, the form every radius takes in this family), so cleared of
/// denominators the scale is `dp²·d_o²·rd·sd` times the `1/(dc²dm²)` the plane's and direction's
/// own denominators contribute — positive throughout, so only signs survive.
struct StripScale {
    /// `sign(U)` — which side of the axis plane the centre is on.
    u_sign: Orient,
    /// `U²`, `W²` and `ρ'²` in the common scale. `ww` is negative exactly when the plane clears
    /// the cylinder altogether, and then there is no strip.
    uu: num_bigint::BigInt,
    ww: num_bigint::BigInt,
    rr: num_bigint::BigInt,
}

fn strip_scale(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    rho2: &BigRat,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripScale {
    use num_bigint::BigInt;
    debug_assert!(!r2.is_negative(), "a squared radius is not negative");
    debug_assert!(!rho2.is_negative(), "a squared margin is not negative");
    debug_assert!(
        dot_sign_rat(&[coeffs[0], coeffs[1], coeffs[2]], m) == Orient::Zero,
        "the strip only exists on a plane parallel to the axis"
    );
    // ★ `_dc` and `_dm` go unused, and that is the derivation showing: the plane's and the
    // direction's own denominators enter only as positive squares in the common denominator, so
    // they cancel. The point's and the origin's do **not** — the `c` term breaks homogeneity in
    // `p`, the same way it does in `point_plane_clearance_rat`.
    let (c, _dc) = lift4(coeffs);
    let (oo, d_o) = lift3(o);
    let (mm, _dm) = lift3(m);
    let (pp, dp) = p.lift();
    let (rn, rd): (BigInt, BigInt) = (r2.numer().clone(), r2.denom().clone());
    let (sn, sd): (BigInt, BigInt) = (rho2.numer().clone(), rho2.denom().clone());
    let n = [c[0].clone(), c[1].clone(), c[2].clone()];
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    debug_assert!(
        (0..3).map(|i| &n[i] * &pp[i]).sum::<BigInt>() + &c[3] * &dp == BigInt::from(0),
        "the point must lie on the plane — the decomposition takes its distance from the axis's, \
         so an off-plane point would be measured against a distance that is not its own"
    );
    // e = n × m, and w = (p − o) scaled by dp·d_o — both integer, both carrying their own
    // positive factor, which is why only the sign of the combination below matters.
    let e: [BigInt; 3] = core::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &n[j] * &mm[k] - &n[k] * &mm[j]
    });
    let w: [BigInt; 3] = core::array::from_fn(|i| &pp[i] * &d_o - &oo[i] * &dp);
    let u = dot(&w, &e);
    // `n·o + c` over the common denominator `dc·d_o`, and the two squared magnitudes.
    let g = dot(&n, &oo) + &c[3] * &d_o;
    let nn = dot(&n, &n);
    let m2 = dot(&mm, &mm);
    StripScale {
        u_sign: orient_of(big_sign(&u)),
        uu: &u * &u * &rd * &sd,
        ww: (&rn * &nn * (&d_o * &d_o) - &g * &g * &rd) * &m2 * (&dp * &dp) * &sd,
        rr: &sn * &nn * &m2 * (&dp * &dp) * (&d_o * &d_o) * &rd,
    }
}

/// **Where a disk of squared radius `rho2` about `p`, lying on the plane, stands relative to the
/// strip.** Exact and total, at any width and any margin.
///
/// The four answers are what a piece with **extent** can say where a point could only say three:
/// clear on the `+` side, clear on the `−` side, reaching the strip without spanning a boundary,
/// or **spanning** one. Consumers fold them differently — the footprint rule treats the last two
/// alike ("did not clear"), while the tangency road needs `Crosses` apart, because a piece that
/// spans the line by itself is exactly the straddle two separate corners would otherwise have to
/// witness between them.
///
/// ```text
/// clear    |U| > W + ρ'        crosses   | |U| − W | < ρ'        else inside
/// ```
///
/// Each is a comparison of a rational against `2√(xy)` for rational `x, y`, so **one sign case
/// and one further squaring** closes it — no radical tower, no new number type. `Crosses` is
/// strict on purpose: a disk touching a ruling at one point spans nothing, which is the same
/// answer a corner sitting on the line gives.
pub fn cylinder_strip_side_margin(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    rho2: &BigRat,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    cylinder_strip_side_extent(
        coeffs,
        &StripReach {
            lo: (p, rho2),
            hi: None,
        },
        o,
        m,
        r2,
    )
}

/// **A piece's reach across the strip, as its two ends** — each end a point on the plane and the
/// **square** of the margin the doors add to it (`ρ²`, the form every radius takes in this family).
///
/// `hi = None` is the symmetric piece a point or a disk is: one end answers both. An arc's ends
/// differ, and stating them is the only way that shape can ask these questions at all.
pub struct StripReach<'a> {
    pub lo: (&'a MeetPoint, &'a BigRat),
    pub hi: Option<(&'a MeetPoint, &'a BigRat)>,
}

/// **Where a piece whose extent across the strip is stated by its two *ends* stands** — the
/// general form of [`cylinder_strip_side_margin`], and the only one an **arc** can use.
///
/// An end is a point on the plane and a squared margin: the low end is `U(p_lo) − ρ_lo'`, the high end
/// `U(p_hi) + ρ_hi'`. A point is both ends with no margin, a disk is one point with the same
/// margin twice — those pass `hi = None` — and an arc's two ends differ, because its angular
/// extent reaches further one way than the other.
///
/// ★★ **`Plus`/`Minus` are exact for every shape; `Crosses` is claimed only where it is proved.**
/// Spanning a boundary is a *positive* fact the tangency road acts on, so under-reporting it is
/// the safe direction — and for an asymmetric extent this door reports [`StripSide::Inside`]
/// rather than prove it. The symmetric case (`hi = None`) keeps the complete answer it always
/// had. An arc that genuinely spans a ruling is therefore "reached the strip, spanned nothing",
/// which every consumer folds as "did not clear"; the day a caller needs the stronger claim from
/// an arc, the proof is `A < W < B` on the two ends and it goes here.
pub fn cylinder_strip_side_extent(
    coeffs: &[Rat; 4],
    reach: &StripReach<'_>,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    use num_bigint::BigInt;
    let (lo, hi) = (reach.lo, reach.hi);
    let side = |o: Orient| match o {
        Orient::Positive => StripSide::Plus,
        Orient::Negative => StripSide::Minus,
        Orient::Zero => StripSide::Inside,
    };
    let s = strip_scale(coeffs, lo.0, lo.1, o, m, r2);
    // The plane clears the cylinder: no strip exists, so nothing can reach it.
    if s.ww.sign() == num_bigint::Sign::Minus {
        return side(s.u_sign);
    }
    let Some(hi) = hi else {
        // Symmetric: one scale answers both boundaries, and `Crosses` with it.
        // `|U| > W + ρ'` — the whole clearance, in the one comparison the scalar family shares.
        if crate::quad::sqrt_exceeds_root_sum(&s.uu, &s.ww, &s.rr) {
            return side(s.u_sign);
        }
        let four = BigInt::from(4);
        let span = &s.uu + &s.ww - &s.rr;
        if span.sign() == num_bigint::Sign::Minus || &span * &span < &four * &s.uu * &s.ww {
            return StripSide::Crosses;
        }
        return StripSide::Inside;
    };
    // Asymmetric: each end answers in its own scale, which is sound because each verdict is a
    // self-contained comparison of that end against `±W`.
    //
    // `U_lo − ρ_lo' > W` is `|U_lo| > W + ρ_lo'` with `U_lo` on the `+` side — the same line,
    // read once per end.
    if s.u_sign == Orient::Positive && crate::quad::sqrt_exceeds_root_sum(&s.uu, &s.ww, &s.rr) {
        return StripSide::Plus;
    }
    let t = strip_scale(coeffs, hi.0, hi.1, o, m, r2);
    if t.ww.sign() == num_bigint::Sign::Minus {
        return side(t.u_sign);
    }
    if t.u_sign == Orient::Negative && crate::quad::sqrt_exceeds_root_sum(&t.uu, &t.ww, &t.rr) {
        return StripSide::Minus;
    }
    StripSide::Inside
}

/// **Does that disk reach *one named* ruling?** — `side` names the ruling in **this** family's
/// vocabulary: the sign of `(p − o)·(n × m)`, the same one [`StripSide::Plus`] is written about.
///
/// ⚠ **The arrangement's `side` is the opposite sign.** `arrangement::ruling_side_signed` measures
/// against `m × n`, so a caller carrying a `MergedRuling` must negate before asking here. Stated
/// rather than absorbed: a predicate that silently accepted either convention would answer about
/// the wrong ruling for whichever caller it was not written for.
///
/// The strip form above answers about **both** boundaries at once, which is what a face-clearance
/// question wants. A caller holding the class's actual edges may carry only one of the two
/// rulings, and asking about the strip would then count a circle that reaches the ruling **not
/// present**. Same derivation, one boundary: `|σW − U| ≤ ρ'`, closed (a touch counts).
pub fn cylinder_ruling_reached(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    rho2: &BigRat,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
    side: i8,
) -> bool {
    cylinder_ruling_reached_extent(
        coeffs,
        &StripReach {
            lo: (p, rho2),
            hi: None,
        },
        o,
        m,
        r2,
        side,
        true,
    )
}

/// **Does a piece whose reach across the strip is stated by its two *ends* touch the named
/// ruling?** — the general form of [`cylinder_ruling_reached`], and the only one an **arc** can use.
///
/// The piece occupies `[U_lo − ρ_lo', U_hi + ρ_hi']` across the strip; the ruling sits at `σW`. It
/// is touched unless the whole reach lies on one side, which is two comparisons of the same shape
/// — each a signed sum of three roots, and each falling to [`quad::sqrt_root_sum_cmp`] once the
/// signs of `U` and `σ` say which term goes where.
///
/// `hi = None` is the symmetric piece (a point, a disk), which is what this door has always been
/// handed; then one scale answers both ends and the verdict is the one it always gave.
///
/// ★★ **`touch_counts` names the boundary rather than assuming one**. Its predecessor is
/// closed — a piece touching the ruling at one point has *reached* it — and that is the right
/// reading for a clearance. The arrangement's net wants the other one: what it cannot mint is a
/// **crossing**, and an edge tangent to another divides nothing. Two propositions, one door, and
/// the caller says which.
pub fn cylinder_ruling_reached_extent(
    coeffs: &[Rat; 4],
    reach: &StripReach<'_>,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
    side: i8,
    touch_counts: bool,
) -> bool {
    use num_bigint::Sign;
    let (lo, hi) = (reach.lo, reach.hi);
    let s_lo = strip_scale(coeffs, lo.0, lo.1, o, m, r2);
    if s_lo.ww.sign() == Sign::Minus {
        return false; // the plane clears the cylinder — it has no ruling here at all
    }
    let sigma = side.signum();
    // `U_lo − √rr > σ√ww`: the whole reach starts above the ruling. When a touch does **not**
    // count, "starts above" includes starting exactly on it, so the comparison relaxes.
    if strip_end_beyond(&s_lo, s_lo.u_sign, sigma, touch_counts) {
        return false;
    }
    let s_hi = match hi {
        None => s_lo,
        Some(h) => {
            let t = strip_scale(coeffs, h.0, h.1, o, m, r2);
            if t.ww.sign() == Sign::Minus {
                return false;
            }
            t
        }
    };
    // `U_hi + √rr < σ√ww` is the same question with `U` and `σ` both negated.
    !strip_end_beyond(&s_hi, orient_neg(s_hi.u_sign), -sigma, touch_counts)
}

/// Is `su·√uu − √rr > σ·√ww` (or `≥`, when `strict` is false)? — one end of a reach against the
/// named ruling, with the signs of `U` and `σ` deciding which side of the identity each root
/// belongs on, and the boundary named rather than assumed.
fn strip_end_beyond(s: &StripScale, su: Orient, sigma: i8, strict: bool) -> bool {
    use crate::quad::sqrt_root_sum_cmp as cmp;
    use num_bigint::Sign;
    let nil = |x: &num_bigint::BigInt| x.sign() == Sign::NoSign;
    let ord = |a: &num_bigint::BigInt, b: &num_bigint::BigInt| if strict { a > b } else { a >= b };
    match (su, sigma) {
        // `√uu ? √rr + √ww`
        (Orient::Positive, 1) => cmp(&s.uu, &s.rr, &s.ww, strict),
        // `√uu ? √rr`
        (Orient::Positive, 0) => ord(&s.uu, &s.rr),
        // `√uu + √ww ? √rr`, i.e. not `√rr ?̄ √uu + √ww` with the boundary flipped
        (Orient::Positive, _) => !cmp(&s.rr, &s.uu, &s.ww, !strict),
        // `√ww ? √rr`
        (Orient::Zero, -1) => ord(&s.ww, &s.rr),
        // `0 ? √rr`, and `0 ? √rr + √ww` — only an equality can hold, and only when not strict.
        (Orient::Zero, 0) => !strict && nil(&s.rr),
        (Orient::Zero, _) => !strict && nil(&s.rr) && nil(&s.ww),
        // `√ww ? √uu + √rr`
        (Orient::Negative, -1) => cmp(&s.ww, &s.uu, &s.rr, strict),
        // `−√uu ? √rr (+ √ww)` — likewise, only as an equality of zeros.
        (Orient::Negative, 0) => !strict && nil(&s.uu) && nil(&s.rr),
        (Orient::Negative, _) => !strict && nil(&s.uu) && nil(&s.rr) && nil(&s.ww),
    }
}

/// `Orient`'s negation — the sign of `−x` from the sign of `x`.
fn orient_neg(o: Orient) -> Orient {
    match o {
        Orient::Positive => Orient::Negative,
        Orient::Negative => Orient::Positive,
        Orient::Zero => Orient::Zero,
    }
}

/// **Which side of the plane at axis parameter `t` a point stands on.** Exact and total, at any
/// width.
///
/// ★ **This is the second separating axis of one rectangle, not a second rule.** A wall parallel
/// to a cylinder's axis meets that cylinder in a rectangle of the wall's own plane: the strip
/// across ([`cylinder_strip_side`], the first axis) and the lateral face's axis-parameter span
/// along (this one). A face misses the rectangle as soon as it clears *either* axis, because a
/// rectangle is the intersection of the two bands — so the two are read together and neither
/// stands as a rule of its own.
///
/// The parameter is written in the **raw** direction's scale — `axis(t) = o + t·m` with `m`
/// unnormalized, the scale `axis_param_of_plane` produces and a lateral face's span is stored in.
/// Re-scaling to a unit axis here would silently mismatch those spans.
///
/// With `s = (p − o)·m / (m·m)` the point's own parameter and `m·m > 0`,
///
/// ```text
/// sign(s − t) = sign((p − o)·m − t·(m·m))
/// ```
///
/// so no division is formed and no square root ever appears. [`Orient::Zero`] is the point sitting
/// exactly on that plane — which the footprint reading treats as *clear*, the uniform-slab theorem
/// speaking of the **open** slab.
///
/// **Precondition:** `m ≠ 0` (`debug_assert`ed) — a zero direction names no axis. Unlike the strip
/// test there is no on-plane precondition: a point's axis parameter is defined wherever it sits.
pub fn point_axis_side(p: &MeetPoint, o: &[Rat; 3], m: &[Rat; 3], t: Rat) -> Orient {
    use num_bigint::BigInt;
    let (oo, d_o) = lift3(o);
    let (mm, d_m) = lift3(m);
    let (pp, dp) = p.lift();
    let (tn, td) = (BigInt::from(t.numer()), BigInt::from(t.denom()));
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    let m2 = dot(&mm, &mm);
    debug_assert!(
        m2 != BigInt::from(0),
        "a zero direction names no axis, so no parameter along it"
    );
    // `w = p − o` over the common denominator `dp·d_o`, integer throughout.
    let w: [BigInt; 3] = core::array::from_fn(|i| &pp[i] * &d_o - &oo[i] * &dp);
    // sign((p−o)·m − t·|m|²) with the positive common denominator `dp·d_o·d_m²·td` cleared:
    //   (W·M)·d_m·td  vs  tn·|M|²·dp·d_o
    // ★ Every cleared factor is positive — `lift3`/`lift` build denominators from an lcm of
    // `Rat` denominators, and `Rat` keeps its sign in the numerator — so the comparison is the
    // sign it claims to be.
    let lhs = dot(&w, &mm) * &d_m * &td;
    let rhs = &tn * &m2 * &dp * &d_o;
    match big_sign(&(lhs - rhs)) {
        1 => Orient::Positive,
        -1 => Orient::Negative,
        _ => Orient::Zero,
    }
}

/// **A branch point's coordinates as `(A + B√C) / D`** — integer throughout, `D > 0`, `C ≥ 0`.
///
/// The two footprint predicates below are the same questions [`cylinder_strip_side`] and
/// [`point_axis_side`] ask; only the point's *description* differs. A `Vertex::Pierce` has no
/// rational coordinates at all — it is `line.base() + s·line.dir()` with `s` quadratic-irrational
/// — so its coordinates live in `ℚ(√c)`, and every quantity those predicates form from a point is
/// a polynomial in them, hence of the shape `X + Y√C` whose sign the tower already answers.
///
/// ★ **The radicand is integerised here, once.** `√(p/q) = √(p·q)/q`, so `s` is rewritten over a
/// single positive denominator with an **integer** radicand `C = p·q`; every sum downstream then
/// carries one `C` and never has to reconcile two spellings of the same surd.
fn lift_branch(
    line: &quad::MeetLine,
    s: &quad::QuadVal,
) -> (
    [num_bigint::BigInt; 3],
    [num_bigint::BigInt; 3],
    num_bigint::BigInt,
    num_bigint::BigInt,
) {
    use num_bigint::BigInt;
    // `s = (Sa + Sb√C) / Ds`, integer and `Ds > 0`: fold the radical's denominator into the
    // rational part (`√(cp/cq) = √(cp·cq)/cq`), then clear both coefficients' denominators.
    let (cp, cq) = (BigInt::from(s.c().numer()), BigInt::from(s.c().denom()));
    let big_c = &cp * &cq;
    let (san, sad) = (BigInt::from(s.a().numer()), BigInt::from(s.a().denom()));
    let (sbn, sbd) = (BigInt::from(s.b().numer()), BigInt::from(s.b().denom()));
    let s_a = &san * &cq * &sbd;
    let s_b = &sbn * &sad;
    let d_s = &cq * &sad * &sbd;
    // `p = base + s·dir` over `db·Ds·dd`, which is positive because every factor is.
    let (bb, db) = lift3(&line.base());
    let (dd, d_d) = lift3(&line.dir());
    let den = &db * &d_s * &d_d;
    let a: [BigInt; 3] = core::array::from_fn(|i| &bb[i] * &d_s * &d_d + &s_a * &dd[i] * &db);
    let b: [BigInt; 3] = core::array::from_fn(|i| &s_b * &dd[i] * &db);
    (a, b, den, big_c)
}

/// **[`cylinder_strip_side`] asked of a branch point** — the same three-way answer, the same
/// expression, and **total** for the same reason: a sign has no width, so the road to it runs in
/// `BigInt` rather than checked `Rat` (the discipline [`quad::cylinder_radial_side`] states).
///
/// **Preconditions** are that function's, and the on-plane one is the caller's to establish —
/// `quad::plane_side(coeffs, line, s) == Orient::Zero` is the branch spelling of it.
pub fn cylinder_strip_side_branch(
    coeffs: &[Rat; 4],
    line: &quad::MeetLine,
    s: &quad::QuadVal,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    use num_bigint::BigInt;
    debug_assert!(!r2.is_negative(), "a squared radius is not negative");
    debug_assert!(
        dot_sign_rat(&[coeffs[0], coeffs[1], coeffs[2]], m) == Orient::Zero,
        "the strip only exists on a plane parallel to the axis"
    );
    let (c, _dc) = lift4(coeffs);
    let (oo, d_o) = lift3(o);
    let (mm, _dm) = lift3(m);
    let (pa, pb, dp, big_c) = lift_branch(line, s);
    let (rn, rd): (BigInt, BigInt) = (r2.numer().clone(), r2.denom().clone());
    let n = [c[0].clone(), c[1].clone(), c[2].clone()];
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    // Exactly the rational road's `e` and `w`, with `w` now a **pair**: the point's rational part
    // and its `√C` part travel together through every linear step.
    let e: [BigInt; 3] = core::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &n[j] * &mm[k] - &n[k] * &mm[j]
    });
    let wa: [BigInt; 3] = core::array::from_fn(|i| &pa[i] * &d_o - &oo[i] * &dp);
    let wb: [BigInt; 3] = core::array::from_fn(|i| &pb[i] * &d_o);
    let (ua, ub) = (dot(&wa, &e), dot(&wb, &e));
    let u_sign = quad::sign1_int(&ua, &ub, &big_c);
    if u_sign == Orient::Zero {
        return StripSide::Inside;
    }
    let g = dot(&n, &oo) + &c[3] * &d_o;
    let nn = dot(&n, &n);
    let m2 = dot(&mm, &mm);
    // `U² = (Ua² + Ub²C) + 2·Ua·Ub·√C`, and the right-hand side is rational — so the difference
    // is one `X + Y√C` and the tower reads its sign. `rd` is the squared radius' denominator.
    let rhs = (&rn * &nn * (&d_o * &d_o) - &g * &g * &rd) * &m2 * (&dp * &dp);
    let x = (&ua * &ua + &ub * &ub * &big_c) * &rd - rhs;
    let y = BigInt::from(2) * &ua * &ub * &rd;
    match quad::sign1_int(&x, &y, &big_c) {
        Orient::Positive if u_sign == Orient::Positive => StripSide::Plus,
        Orient::Positive => StripSide::Minus,
        _ => StripSide::Inside,
    }
}

/// **[`point_axis_side`] asked of a branch point** — same expression, same total contract.
pub fn point_axis_side_branch(
    line: &quad::MeetLine,
    s: &quad::QuadVal,
    o: &[Rat; 3],
    m: &[Rat; 3],
    t: Rat,
) -> Orient {
    use num_bigint::BigInt;
    let (oo, d_o) = lift3(o);
    let (mm, d_m) = lift3(m);
    let (pa, pb, dp, big_c) = lift_branch(line, s);
    let (tn, td) = (BigInt::from(t.numer()), BigInt::from(t.denom()));
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    let m2 = dot(&mm, &mm);
    debug_assert!(
        m2 != BigInt::from(0),
        "a zero direction names no axis, so no parameter along it"
    );
    let wa: [BigInt; 3] = core::array::from_fn(|i| &pa[i] * &d_o - &oo[i] * &dp);
    let wb: [BigInt; 3] = core::array::from_fn(|i| &pb[i] * &d_o);
    let x = dot(&wa, &mm) * &d_m * &td - &tn * &m2 * &dp * &d_o;
    let y = dot(&wb, &mm) * &d_m * &td;
    quad::sign1_int(&x, &y, &big_c)
}

/// A rational 4-vector (plane coefficients) as **(integer components, positive denominator)**.
fn lift4(v: &[Rat; 4]) -> ([num_bigint::BigInt; 4], num_bigint::BigInt) {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let den: [BigInt; 4] = core::array::from_fn(|i| BigInt::from(v[i].denom()));
    let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
    let num = core::array::from_fn(|i| BigInt::from(v[i].numer()) * (&d / &den[i]));
    (num, d)
}

/// `big_sign`'s answer as the geometry's [`Orient`].
fn orient_of(s: i8) -> Orient {
    match s {
        1 => Orient::Positive,
        -1 => Orient::Negative,
        _ => Orient::Zero,
    }
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
/// corners solve to points that fit `Rat` (measured 8/8), but the road there overflowed.
/// The fallback runs only on the decline path, so the narrow
/// route's cost and answers are untouched.
///
/// This is what replaces a dissolved sketch-frame base vertex: the frame-shared triple of a
/// prism corner, solved in the frame the planes are stated in, is the corner's exact base —
/// measured bit-identical to the stored base-and-replay road (8/8).
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

/// [`three_planes_rat`]'s narrow route — Cramer over `Rat` with checked arithmetic, whose
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

    /// **(integer components, positive common denominator)** — `lift3`'s twin for a point whose
    /// width may exceed `Rat`.
    ///
    /// This is what lets a sign predicate take a meet *whatever* its width: the answer is a
    /// polynomial in these integers, so `Narrow` and `Wide` travel the same road and no
    /// consumer has to carry a width limit inside its refusal.
    pub fn lift(&self) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
        use num_bigint::BigInt;
        use num_integer::Integer;
        match self {
            MeetPoint::Narrow(p) => lift3(p),
            MeetPoint::Wide(p) => {
                let d = p.iter().fold(BigInt::from(1), |l, (_, den)| l.lcm(den));
                let num = core::array::from_fn(|i| {
                    let (n, den) = &p[i];
                    n * (&d / den)
                });
                (num, d)
            }
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
/// `√v` when it is rational — numerator and denominator both perfect squares — else `None`
/// (also for `v < 0`). The sketch's arcs ask this of `|start − centre|²`: a radius the kernel can
/// state is a rational one, and this is the only place that question is answered.
/// **The sign of `a + b·π`** for integers `a`, `b`, exactly. π is irrational, so the sum is zero
/// only when both are; otherwise `f(x) = a + b·x` is monotone and its sign at π is the sign it has
/// at both ends of a rational bracket around π — here π's decimal expansion to 37 places. If the
/// two ends disagree, π sits where `f` changes sign and the answer is `None`: undecidable at this
/// width, never a guess.
///
/// The consumer is the prism builder's winding question: a ring whose arcs are quarter turns has
/// twice its signed area equal to such a sum (`nacre-ops`, `Ring2d::winding_sign`).
pub fn sign_a_plus_b_pi_int(a: &num_bigint::BigInt, b: &num_bigint::BigInt) -> Option<Orient> {
    use num_bigint::BigInt;
    use num_traits::{Signed, Zero};
    let sign = |x: &BigInt| {
        if x.is_zero() {
            Orient::Zero
        } else if x.is_positive() {
            Orient::Positive
        } else {
            Orient::Negative
        }
    };
    if b.is_zero() {
        return Some(sign(a));
    }
    let den: BigInt = "10000000000000000000000000000000000000"
        .parse()
        .expect("10^37");
    let lo: BigInt = "31415926535897932384626433832795028841"
        .parse()
        .expect("π low");
    let hi: BigInt = "31415926535897932384626433832795028842"
        .parse()
        .expect("π high");
    // f(x) at x = num/den has the sign of a·den + b·num.
    let at = |num: &BigInt| sign(&(a * &den + b * num));
    let (s_lo, s_hi) = (at(&lo), at(&hi));
    (s_lo == s_hi && s_lo != Orient::Zero).then_some(s_lo)
}

/// A circular arc for [`winding_sign_quarter_arcs`]: `start → end` about `center`,
/// counter-clockwise when `ccw`; `start == end` is the whole circle. `r2` is the squared radius —
/// the form the sketch truth holds, and the one the area term wants.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QuarterArc {
    pub center: [Rat; 2],
    pub r2: Rat,
    pub start: [Rat; 2],
    pub end: [Rat; 2],
    pub ccw: bool,
}

/// **Which way a closed ring of straight steps and quarter-turn arcs runs**, exactly: `Positive`
/// is counter-clockwise (`+x → +y`). Twice the signed area is `∮ x dy − y dx`; a straight step
/// `p → q` contributes the shoelace term `p × q`, an arc `S → E` about `C` through the signed
/// angle `φ` contributes `C × (E − S) + r²·φ` (parametrize `p = C + r·u(t)`: `x dy − y dx =
/// r·(C × u′) + r²` dt). With every `φ` a multiple of `π/2` the total is `a + b·π`, whose sign
/// [`sign_a_plus_b_pi_int`] decides.
///
/// ★ **Integers, not `Rat`.** Every coordinate goes over one common denominator and the sums run
/// in `BigInt`: two 17-digit decimals multiplied already leave `i128` (the reason [`orient2d_rat`]
/// carries a `BigInt` arm), and a winding read that failed open on overflow would build a solid
/// inside out. The area is doubled so `r²·k/2` stays integral. `None` for an arc that is not a
/// quarter-turn multiple (no producer makes one today — the kit's fillet is axis-aligned, so its
/// `arc_rat` arcs are quarter-turns; the type allows more) or for a sum that lands inside the π
/// bracket.
pub fn winding_sign_quarter_arcs(lines: &[[[Rat; 2]; 2]], arcs: &[QuarterArc]) -> Option<Orient> {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;
    let mut den = BigInt::from(1);
    let mut fold = |r: &Rat| den = den.lcm(&BigInt::from(r.denom()));
    for [p, q] in lines {
        for c in p.iter().chain(q.iter()) {
            fold(c);
        }
    }
    for a in arcs {
        for c in a.center.iter().chain(a.start.iter()).chain(a.end.iter()) {
            fold(c);
        }
        fold(&a.r2);
    }
    let int = |r: &Rat| -> BigInt { BigInt::from(r.numer()) * (&den / BigInt::from(r.denom())) };
    let pt = |p: &[Rat; 2]| -> [BigInt; 2] { [int(&p[0]), int(&p[1])] };
    let cross2 = |a: &[BigInt; 2], b: &[BigInt; 2]| -> BigInt { &a[0] * &b[1] - &a[1] * &b[0] };
    let dot2 = |a: &[BigInt; 2], b: &[BigInt; 2]| -> BigInt { &a[0] * &b[0] + &a[1] * &b[1] };
    let sub2 = |a: &[BigInt; 2], b: &[BigInt; 2]| -> [BigInt; 2] { [&a[0] - &b[0], &a[1] - &b[1]] };
    // `a_int` collects 2·(shoelace + C × (E − S)), `b_int` collects r²·k — units of den².
    let (mut a_int, mut b_int) = (BigInt::from(0), BigInt::from(0));
    for [p, q] in lines {
        a_int += 2 * cross2(&pt(p), &pt(q));
    }
    for arc in arcs {
        let (c, s0, e0) = (pt(&arc.center), pt(&arc.start), pt(&arc.end));
        a_int += 2 * cross2(&c, &sub2(&e0, &s0));
        // The counter-clockwise quarter turns from S to E.
        let quarters_ccw: i128 = if s0 == e0 {
            4
        } else {
            let (v1, v2) = (sub2(&s0, &c), sub2(&e0, &c));
            match (dot2(&v1, &v2).sign(), cross2(&v1, &v2).sign()) {
                (Sign::NoSign, Sign::Plus) => 1,
                (Sign::Minus, Sign::NoSign) => 2,
                (Sign::NoSign, Sign::Minus) => 3,
                _ => return None,
            }
        };
        // φ = k·π/2 signed by the direction walked — a whole turn is ±4 either way, a partial arc
        // walked clockwise is the complement, negated; 2·r²·φ = (r²·k)·π.
        let k = match (arc.ccw, quarters_ccw) {
            (true, q) => q,
            (false, 4) => -4,
            (false, q) => q - 4,
        };
        // `r²` is stated, so `int` puts it in units of `den`; one more `den` brings it to the
        // `den²` the area terms carry.
        b_int += int(&arc.r2) * &den * BigInt::from(k);
    }
    sign_a_plus_b_pi_int(&a_int, &b_int)
}

/// [`sign_a_plus_b_pi_int`] for rationals: `a + b·π` and `(a·db + b·da·π)` share a sign since
/// `da·db > 0`.
pub fn sign_a_plus_b_pi(a: Rat, b: Rat) -> Option<Orient> {
    use num_bigint::BigInt;
    let (an, ad) = (BigInt::from(a.numer()), BigInt::from(a.denom()));
    let (bn, bd) = (BigInt::from(b.numer()), BigInt::from(b.denom()));
    sign_a_plus_b_pi_int(&(an * &bd), &(bn * &ad))
}

/// **A rational of any width** — the truth's spelling for a quantity that is the *square* of a
/// stated number: a cylinder's `r²`. Every `Rat` has a square, and this is where it fits when
/// `i128` does not (a 16-digit decimal below `1e-4` has a denominator whose square leaves `i128`),
/// so a statement is never refused for its square being wide — the reason `MeetPoint::Wide`
/// stands beside `Narrow`. Lowest terms, positive denominator (`Ratio`'s invariant).
///
/// A *storage and door* type: the exact predicates lift to `BigInt` anyway and take it directly;
/// the few `Rat` arithmetic sites on a squared radius ask [`BigRat::narrow`] and decline exactly
/// where they used to overflow on `r·r`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BigRat(Ratio<num_bigint::BigInt>);

impl BigRat {
    pub fn from_rat(r: Rat) -> BigRat {
        BigRat(Ratio::new(
            num_bigint::BigInt::from(r.numer()),
            num_bigint::BigInt::from(r.denom()),
        ))
    }

    /// `r²`, exactly — the one product this type exists for. Never overflows; `r` in lowest terms
    /// makes `r²` so too.
    pub fn square_of(r: Rat) -> BigRat {
        let (n, d) = (
            num_bigint::BigInt::from(r.numer()),
            num_bigint::BigInt::from(r.denom()),
        );
        BigRat(Ratio::new_raw(&n * &n, &d * &d))
    }

    pub fn zero() -> BigRat {
        BigRat(Ratio::from_integer(num_bigint::BigInt::from(0)))
    }

    pub fn numer(&self) -> &num_bigint::BigInt {
        self.0.numer()
    }

    pub fn denom(&self) -> &num_bigint::BigInt {
        self.0.denom()
    }

    /// The `Rat` this is, when both parts fit `i128` — `None` is width, not a value.
    pub fn narrow(&self) -> Option<Rat> {
        use num_traits::ToPrimitive;
        Rat::new(self.0.numer().to_i128()?, self.0.denom().to_i128()?)
    }

    pub fn is_positive(&self) -> bool {
        self.0.numer().sign() == num_bigint::Sign::Plus
    }

    pub fn is_negative(&self) -> bool {
        self.0.numer().sign() == num_bigint::Sign::Minus
    }

    pub fn is_zero(&self) -> bool {
        self.0.numer().sign() == num_bigint::Sign::NoSign
    }

    /// `self · r`, exactly.
    pub fn mul_rat(&self, r: Rat) -> BigRat {
        BigRat(&self.0 * BigRat::from_rat(r).0)
    }
}

impl From<Rat> for BigRat {
    fn from(r: Rat) -> BigRat {
        BigRat::from_rat(r)
    }
}

/// [`rat_sqrt_exact`] for a wide radicand: the root, when it is rational **and** fits `Rat` —
/// which every stated radius does, its square being what widened.
pub fn rat_sqrt_exact_big(v: &BigRat) -> Option<Rat> {
    if v.is_negative() {
        return None;
    }
    let (n, d) = (v.numer(), v.denom());
    let (rn, rd) = (n.sqrt(), d.sqrt());
    if &rn * &rn != *n || &rd * &rd != *d {
        return None;
    }
    BigRat(Ratio::new_raw(rn, rd)).narrow()
}

pub fn rat_sqrt_exact(v: Rat) -> Option<Rat> {
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
pub fn inv_sqrt_bounded(v: Rat, prec: usize) -> Option<HpBounded> {
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
/// The wide-frame realization feeds plane data wider than `i128` through this; keeping the
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

/// [`inv_sqrt_bounded`] for a squared length wider than `i128` — the wide-frame twin.
///
/// Same ladder, same derived bound: the integer enters **exactly** ([`bigint_to_bigfloat`] at
/// its own bit length), the square root rounds once (halving the incoming relative error, which
/// is zero here) and the reciprocal rounds once — under the narrow twin's `4u`, which is kept
/// for symmetry rather than tightened.
///
/// **Uncached** — the wide population is a fraction of a percent of pushes and the memo key
/// would be a `BigInt`; measured before optimizing, per the cache philosophy.
pub fn inv_sqrt_bigint_bounded(v: &num_bigint::BigInt, prec: usize) -> Option<HpBounded> {
    if v.sign() != num_bigint::Sign::Plus {
        return None;
    }
    let ip = prec.max(128);
    let x = bigint_to_bigfloat(v, ip);
    let one = BigFloat::from_i128(1, ip);
    let z = one.div(&x.sqrt(prec, HP_RM), prec, HP_RM);
    let u = Mag::pow2(-(prec as i64));
    let rel = u.times(Mag::of(4.0));
    let error = Mag::above(&z).times(rel);
    Some(HpBounded::new(z, error))
}

/// [`inv_sqrt_bounded`] without the memo — the evaluation itself, kept separate so no `INV_SQRT`
/// borrow is held across the arbitrary-precision work.
fn realize_inv_sqrt(v: Rat, prec: usize) -> HpBounded {
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
    let error = Mag::above(&z).times(rel); // `|x| < 2^exponent`
    HpBounded::new(z, error)
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
        let HpBounded {
            value: z,
            error: rad,
        } = realize_inv_sqrt_memoized(v, prec);
        if let Some(f) = round_to_f64(&z, rad, prec) {
            if i > 0 {
                INV_SQRT_ESCALATED.with_borrow_mut(|(e, _)| *e += 1);
            }
            return f;
        }
    }
    INV_SQRT_ESCALATED.with_borrow_mut(|(_, f)| *f += 1);
    let z = realize_inv_sqrt_memoized(v, 256).value;
    to_f64_exact(&z).unwrap_or(f64::NAN)
}

/// [`inv_sqrt_bounded`] for a `v` already known positive.
fn realize_inv_sqrt_memoized(v: Rat, prec: usize) -> HpBounded {
    inv_sqrt_bounded(v, prec).expect("v > 0 checked by the caller")
}

/// **`√v` as the `f64` nearest the true value**, or `None` when `v < 0` — the f64 realization of
/// a radius stated as its square, at any width.
///
/// The twin of [`inv_sqrt_f64`]: the exact branch runs first ([`rat_sqrt_exact`] — every radius a
/// user writes as a decimal lands here, and comes back as the `to_f64` of the rational it always
/// was, bit for bit), and only a genuinely irrational root reaches the 128 → 256 bit ladder, under
/// the same correct-rounding contract. `√v` is nonzero for `v > 0`, so the interval never
/// straddles zero and [`round_to_f64`] cannot decline forever; the documented fallback to the
/// midpoint is kept for the same reason the inverse has it.
pub fn sqrt_f64(v: &BigRat) -> Option<f64> {
    if v.is_negative() {
        return None;
    }
    if let Some(r) = rat_sqrt_exact_big(v) {
        return Some(r.to_f64());
    }
    for prec in [128usize, 256] {
        let HpBounded {
            value: z,
            error: rad,
        } = sqrt_bounded_big(v, prec)?;
        if let Some(f) = round_to_f64(&z, rad, prec) {
            return Some(f);
        }
    }
    let z = sqrt_bounded_big(v, 256)?.value;
    Some(to_f64_exact(&z).unwrap_or(f64::NAN))
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
    let HpBounded {
        value: h,
        error: rad,
    } = inv_sqrt_bounded(v, P)?;
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
/// the faces' coordinates instead.
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
        let (c, s) = self.cos_sin_bounded(prec);
        (c.value, s.value)
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
    /// Returns `(cos, sin)`, each with the radius it carries.
    pub fn cos_sin_bounded(self, prec: usize) -> (HpBounded, HpBounded) {
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
    fn realize_cos_sin(self, prec: usize) -> (HpBounded, HpBounded) {
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
            let d_theta = Mag::above(&rad).times(rel);
            let (c, s) = (rad.cos(prec, HP_RM, cc), rad.sin(prec, HP_RM, cc));
            // `|x| < 2^exponent` — the slope of the *other* function, and the scale of the
            // half-ulp of this one.
            let (uc, us) = (Mag::above(&c), Mag::above(&s));
            let err_cos = us.times(d_theta).plus(uc.times(u));
            let err_sin = uc.times(d_theta).plus(us.times(u));
            (HpBounded::new(c, err_cos), HpBounded::new(s, err_sin))
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
            let (c, s) = self.cos_sin_bounded(prec);
            if let (Some(cf), Some(sf)) = (
                round_to_f64(&c.value, c.error, prec),
                round_to_f64(&s.value, s.error, prec),
            ) {
                if i > 0 {
                    ROUND_ESCALATED.with_borrow_mut(|(e, _)| *e += 1);
                }
                return (cf, sf);
            }
        }
        ROUND_ESCALATED.with_borrow_mut(|(_, f)| *f += 1);
        let (c, s) = self.cos_sin_bounded(256);
        (
            to_f64_exact(&c.value).unwrap_or(f64::NAN),
            to_f64_exact(&s.value).unwrap_or(f64::NAN),
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
    /// The reading is at octave granularity (`2^exponent` of the residual), so it can sit up to 2×
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
        let (hc, hs) = self.cos_sin_bounded(P);
        let (rc, rs) = (hc.error, hs.error);
        let (hc, hs) = (hc.value, hs.value);
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
/// `pivot`) by the rational `angle`. Exact for the 90°-family (`try_exact_cos_sin`);
/// otherwise the realized coordinate is irrational (cos/sin) and carries tol.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotation {
    pub axis: Axis,
    pub pivot: [Rat; 3],
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
            let (px, py) = (r.pivot[i].to_f64(), r.pivot[j].to_f64());
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
                (m, r.pivot)
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
                let u = x[i].checked_sub(r.pivot[i])?;
                let v = x[j].checked_sub(r.pivot[j])?;
                let mut m = x;
                m[i] = r.pivot[i]
                    .checked_add(u.checked_mul(cos)?.checked_sub(v.checked_mul(sin)?)?)?;
                m[j] = r.pivot[j]
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

    /// **This isometry's rotation applied to a rational direction**, exactly — the rational twin
    /// of [`apply_dir`](Isometry::apply_dir). A direction is a difference of points, so the pivot
    /// and the translation cancel and only the turn remains. `None` under the same conditions as
    /// [`point_rat`](Isometry::point_rat) — a rotation the rationals cannot state, or `i128`
    /// overflow.
    pub fn dir_rat(&self, d: [Rat; 3]) -> Option<[Rat; 3]> {
        match self.rotate {
            None => Some(d),
            Some(r) => {
                let (cos, sin) = r.angle.try_exact_cos_sin()?;
                let (i, j) = r.axis.plane();
                let mut m = d;
                m[i] = d[i].checked_mul(cos)?.checked_sub(d[j].checked_mul(sin)?)?;
                m[j] = d[i].checked_mul(sin)?.checked_add(d[j].checked_mul(cos)?)?;
                Some(m)
            }
        }
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
#[path = "tests/lib.rs"]
mod tests;

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

#[cfg(test)]
#[path = "tests/decimal_realization.rs"]
mod decimal_realization;
