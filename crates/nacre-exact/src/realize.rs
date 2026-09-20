use super::*;
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
/// `compute_hp` already buys **only where a judgement actually needs it**. (`nacre-judge` depends on
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
pub(super) fn to_f64_exact(x: &BigFloat) -> Option<f64> {
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
pub(super) fn big_to_ratio(x: &BigFloat) -> Option<(num_bigint::BigInt, num_bigint::BigInt)> {
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

/// How many entries [`TRIG`] holds — for the tests that pin the memo actually memoizes.
///
/// **The count, not a hit tally.** A hit rate cannot tell "the memo works" from "the memo is
/// fragmenting": a caller that spelled one angle two ways would show high hits while paying twice
/// for every angle. Every miss inserts, so the entry count is the direct reading.
#[cfg(test)]
pub(super) fn trig_entries() -> usize {
    TRIG.with_borrow(|t| t.len())
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
pub(super) fn nearest_f64(n: u128, d: u128) -> f64 {
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
