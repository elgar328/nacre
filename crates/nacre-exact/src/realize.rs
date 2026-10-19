use super::*;
/// **The `f64` nearest the true value that `mid ± rad` encloses — or `None` when `mid ± rad` is
/// not narrow enough to say, or when that `f64` is not normal.**
///
/// An arbitrary-precision realization is an *interval*, and rounding its midpoint to 53 bits is
/// rounding an approximation: near a rounding boundary the answer would be the approximation's,
/// not the truth's. So both ends are rounded and compared. Round-to-nearest is **monotonic**, so
/// `lo ≤ v ≤ hi` gives `round(lo) ≤ round(v) ≤ round(hi)` — when the ends agree, the value between
/// them cannot round anywhere else, whatever it is. `None` says "realize deeper and ask again",
/// which is the same shape as `precision_for` escalating a judgement.
///
/// **The ends are formed at `mid`'s own precision plus 64 bits, and that line is load-bearing** —
/// the working precision is read off the value, never taken from the caller, and
/// [`round_to_digits`] forms its ends here too. Fewer bits than `mid` has round `mid` itself first:
/// with a zero radius both ends are then that rounding, and a value just past a rounding boundary
/// reads as the boundary and ties the other way (at 192 bits `1 + 2⁻⁵³ + 2⁻²⁵⁰` reads as `1.0`,
/// where the nearest is `1 + 2⁻⁵²`). No guard bits at all round the radius away instead — `rad` is
/// about `2⁻ᵖ` of a `p`-bit `mid` — so both ends collapse onto `mid` and the comparison passes for
/// every input, a check that verifies nothing. With the guard bits the sums are exact whenever `rad`
/// is zero, or no finer than `2⁻⁶⁴` of `mid`'s last bit (the rounding every `HpBounded` operation
/// charges its own result keeps it above that) and no wider than `mid` (a wider one reaches across
/// zero, where the ends never agree) — and then `Some` is a proof.
///
/// **Not for a value whose true magnitude may be zero.** `cos 90°` is exactly `0`, so its interval
/// straddles zero, the two ends have opposite signs, and no precision ever makes them agree — this
/// would return `None` forever. Callers must resolve the exactly-representable cases first;
/// [`Angle::cos_sin_f64`] does that with `cos_sin_rat`.
///
/// `pub` because **STEP export** will want the same rounding, and two implementations of "round
/// this interval to f64" is exactly the drift this crate keeps deleting. STEP carries `f64`
/// coordinates and nothing else — not the motion chain, not the three-plane definition — so export
/// is the one place where a point's exact truth must be realized as well as `f64` allows, and it
/// can be *verified* as it goes: this returning `Some` **is** the proof that the written coordinate
/// is the nearest `f64` to the exact one.
///
/// **Not for the coordinates the kernel works with.** Rounding `px + u·c − v·s` as a whole
/// would make every `WitnessPoint::coord_approx` literally `round(coord_hp_uncached)` and take
/// `tol` to a half-ulp, and it was **declined**: that is arbitrary precision *per vertex*, where
/// [`Angle::cos_sin_f64`]'s is per *angle* and memoised, so it would pay at construction for a
/// precision that `WitnessPoint`'s lazy `coord_hp_uncached` already buys **only where a judgement
/// actually needs it**. (`nacre-judge` depends on this crate, so that type cannot be named here as
/// a link.) Cheaper ways to shrink that arithmetic (an FMA, a compensated
/// evaluation) stay in `f64` and keep the laziness, so they are the candidates if the term ever
/// needs to move.
pub fn round_to_f64(mid: &BigFloat, rad: Mag) -> Option<f64> {
    if mid.is_nan() || mid.is_inf() {
        return None;
    }
    let (lo, hi) = interval_ends(mid, rad)?;
    let (rlo, rhi) = (bigfloat_to_f64(&lo)?, bigfloat_to_f64(&hi)?);
    (rlo == rhi && !rlo.is_nan()).then_some(rlo)
}

/// `mid ± rad`, formed at `mid`'s own precision plus 64 bits — see [`round_to_f64`] for why.
fn interval_ends(mid: &BigFloat, rad: Mag) -> Option<(BigFloat, BigFloat)> {
    let p = bits_of(mid) + 64;
    let r = rad_upper_bigfloat(rad, p)?;
    Some((mid.sub(&r, p, HP_RM), mid.add(&r, p, HP_RM)))
}

/// Is every value in `mid ± rad` one whose nearest `f64` is not normal, **on one side** — below the
/// range or above it, with one sign? Then no precision makes [`round_to_f64`] decide it, and
/// [`ladder_f64`] decides it as `NaN` instead.
///
/// Near the range the ends are formed as [`round_to_f64`] forms them and read the way
/// [`bigfloat_to_f64`] reads them: both refused, both on the same side. Far from it (`|exponent| >
/// 2048`) they are not formed — the radius there can sit past [`rad_upper_bigfloat`]'s exponent
/// guard — and a radius below a quarter of `|mid|` keeps the whole interval on `mid`'s side.
fn outside_normal_range(mid: &BigFloat, rad: Mag) -> bool {
    let Some(e) = mid.exponent().filter(|_| !mid.is_zero()) else {
        return false;
    };
    let e = i64::from(e);
    if e.abs() > 2048 {
        return rad.exp2().is_none_or(|r| r < e - 2);
    }
    let Some((lo, hi)) = interval_ends(mid, rad) else {
        return false;
    };
    let side = |x: &BigFloat| match (bigfloat_to_f64(x), x.exponent()) {
        (None, Some(xe)) if !x.is_zero() => Some((xe > 0, x.is_negative())),
        _ => None,
    };
    matches!((side(&lo), side(&hi)), (Some(a), Some(b)) if a == b)
}

/// The depth past which `Angle::cos_sin_f64`'s ladder does not climb — `2¹⁵` bits — **an alarm,
/// not a proof.** The cosine of a rational number of degrees is algebraic of degree up to the
/// angle's denominator, so no usable bound says how close it can come to a rounding boundary; a
/// legitimate angle reaches this only if its cosine or sine lies within about `2⁻⁽³²⁷⁶⁸⁻ᶜ⁾` of one
/// (`c` the constant its realization's radius carries near an axis). It is as high as cost allows,
/// measured (native release, one realization): `cos_sin_hp` takes 68 ms at 2¹⁵ bits and 281 ms at
/// 2¹⁶, so a climb to the cap costs about 92 ms over all its rungs and one to 2¹⁶ about 373 ms. A
/// broken premise should cost a third of a second, not stall the page — and the product runs as
/// wasm, slower than native — so the cap is the last power of two whose climb stays well under it.
pub(crate) const TRIG_LADDER_CAP: usize = 1 << 15;

/// **The depth that decides every value a `√`-family realization hands [`ladder_f64`]** — the
/// smallest rung at or above `max_bits + 128`, where `max_bits` is the wider of `p` and `q` for the
/// value's square `p/q`, in any representation (lowest terms are not needed).
///
/// Every point a rung can stick at — a rounding boundary, the bottom edge `(2⁵³ − 1)·2⁻¹⁰⁷⁵` or the
/// overflow edge of the range — is a dyadic `b = B/2ᵗ` with an odd `B < 2⁵⁴` (a midpoint's
/// mantissa `2M + 1`), and the value `x` is none of them (its caller's exact branch takes every
/// rational root). So `x² − b² = (p·2²ᵗ − q·B²)/(q·2²ᵗ)` is a nonzero rational and
/// `|x² − b²| ≥ 1/(q·2²ᵗ)`. For `x` in `[2ᴱ, 2ᴱ⁺¹)` and a boundary in its binade or the next one
/// either side, `2ᵗ ≤ 2⁵⁴⁻ᴱ` and the exponents cancel: `2²ᵗ·x·(x + b) < 2¹¹¹`, so
/// `|x − b|/x > 2^-(bits(q) + 111)`; a boundary farther off is farther. An integer boundary
/// (`t < 0`) gives `1/q` instead, and `bits(p)` takes the role with room to spare. Against that gap
/// the realizations (`inv_sqrt_hp`, `inv_sqrt_hp_big`, `sqrt_hp_big`, a unit vector's products)
/// carry a relative radius of at most `18 · 2⁻ᵖʳᵉᶜ`, which [`round_to_f64`] widens to a power of
/// two, so both ends lie within `72 · 2⁻ᵖʳᵉᶜ` of the truth — about ten bits short of the gap at
/// `max_bits + 128`. The bound is tight: `√(p/q)` with `p·2¹⁰⁸ − q·(2⁵⁴−1)² = 1`, both 108 bits,
/// sits a relative `2⁻²¹⁷` from `1 − 2⁻⁵⁴`, two bits inside it, and needs the 256 bits it earns
/// (`the_tightest_radicand_needs_the_whole_separation_margin`).
///
/// The proof holds while the depth is at most `2¹⁹`: deeper, a realization's radius falls below the
/// `2^-10⁶` that `rad_upper_bigfloat` assembles and no rung decides, so a wider input that near a
/// boundary reaches the cap and panics. Converting such an input alone takes minutes.
pub(crate) fn separation_cap(max_bits: u64) -> usize {
    (max_bits as usize + 128).next_power_of_two()
}

/// **The `f64`s nearest the true values a realization encloses**, climbing `128, 256, 512, …` bits
/// until each value is decided, up to `cap`. The one ladder for a value this crate can only
/// approach (`Angle::cos_sin_f64`, `inv_sqrt_f64`, `sqrt_f64`, `unit_vector_f64` — their exact
/// branches round a rational once, by `Rat::to_f64` or [`nearest_normal_f64`]).
///
/// `realize(prec)` gives the values at `prec` bits. A value is decided when [`round_to_f64`] names
/// its nearest `f64`, or when every value its interval holds has a nearest `f64` that is not
/// normal, on one side ([`outside_normal_range`]) — then it is `NaN`, a value an `f64` cannot say.
/// A rung answers when it decides every value. `None` when a realization is `None`.
///
/// **Why it terminates, and where.** A rung fails to decide only a value within its radius of a
/// rounding boundary (the midpoint between two neighbouring `f64`s), of an edge of the range — both
/// dyadic rationals — or of zero (an interval across zero has subnormal ends of two signs). Each
/// caller hands the ladder only values that are none of those, or an exact zero with a zero radius:
/// its exact branch takes every value that could be one, and what is left is irrational or a
/// rational that is neither (each caller says which). Every realization's radius shrinks as
/// `c · 2⁻ᵖʳᵉᶜ` — relative to the value, or for a cosine near an axis absolute (its error comes
/// through the derivative) — so some rung clears the gap. How deep is the caller's to say, as
/// `cap`: for the `√` family [`separation_cap`] proves a depth that decides every value, so
/// reaching it means the caller's premise broke; for a cosine no such bound exists, and
/// [`TRIG_LADDER_CAP`] is a cost alarm. Either way the ladder panics there rather than climbing
/// forever or handing back a guess.
///
/// **Doubling, not the computed jump judging makes.** The structural cases are gone before the
/// ladder, so an undecided rung is a genuine lack of bits, and the distance to the nearest boundary
/// is not known until it is computed — so each rung doubles, and a value that needs `n` bits is
/// realized at no more than about `2n`.
///
/// **128 is `nacre_judge`'s trial precision**, so a memoized realization (`cos_sin_hp`,
/// `inv_sqrt_hp`) a judgement also asks for is shared rather than paid twice — which is why the
/// ladder does not start lower: 64 bits measured no cheaper for the trigonometric realization
/// (31.9µs against 32.3µs — the cost is setup, not bit count) and would be a different memo key.
/// The wide realizations (`sqrt_hp_big`, `inv_sqrt_hp_big`) keep no memo; they climb the same
/// rungs so that all four answer under one rule.
///
/// A test build counts the answers decided at 256 bits and past them (`LADDER_CLIMBED`), so a climb
/// nothing in the corpus needs stays visible when something starts to.
pub(crate) fn ladder_f64<const N: usize>(
    cap: usize,
    mut realize: impl FnMut(usize) -> Option<[HpBounded; N]>,
) -> Option<[f64; N]> {
    let mut prec = 128;
    loop {
        let r = realize(prec)?;
        let mut out = [0.0; N];
        let decided = out.iter_mut().zip(&r).all(|(o, x)| {
            // Rounding first: a decided value has both ends normal (or zero), so it is never also
            // wholly outside the range, and the common case forms its ends once.
            let v = match round_to_f64(&x.value, x.error) {
                Some(v) => v,
                None if outside_normal_range(&x.value, x.error) => f64::NAN,
                None => return false,
            };
            *o = v;
            true
        });
        if decided {
            #[cfg(test)]
            match prec {
                128 => {}
                256 => LADDER_CLIMBED.with_borrow_mut(|(at, _)| *at += 1),
                _ => LADDER_CLIMBED.with_borrow_mut(|(_, past)| *past += 1),
            }
            return Some(out);
        }
        assert!(
            prec < cap,
            "the nearest-f64 ladder reached its cap ({cap} bits) undecided: a value on a \
             rounding boundary, a normal-range edge or zero reached it, or one closer to such a \
             point than its caller's bound allows"
        );
        prec *= 2;
    }
}

/// The nearest `f64` to the exact rational `num/den`, or `NaN` when it is nonzero and that `f64` is
/// not normal — the exact branches' reading, the rule [`ladder_f64`] decides by too, so a value
/// gets one answer whichever road decided it. [`nearest_f64_big`] rounds once, on the subnormal
/// grid below the range, so the `f64` it names is the nearest even there and `is_normal` on it is
/// the rule. Zero is read off the numerator, not the result: [`nearest_f64_big`] rounds a nonzero
/// value below the subnormals to `0.0`.
pub(crate) fn nearest_normal_f64(num: &num_bigint::BigInt, den: &num_bigint::BigInt) -> f64 {
    if num.sign() == num_bigint::Sign::NoSign {
        return 0.0;
    }
    match nearest_f64_big(num, den) {
        Some((v, _)) if v.is_normal() => v,
        _ => f64::NAN,
    }
}

/// `2^exp2()` as a `BigFloat` — an **upper bound** on the radius, exactly representable.
///
/// [`Mag`] is `m · 2^e` with `m ∈ [0.5, 1)`, so `2^e` is above it; the mantissa is left out
/// because widening the interval can only cost an escalation, never buy a wrong acceptance, and
/// `Mag` does not expose its mantissa. A power of two is exact in `BigFloat` at any precision.
///
/// **Built in `BigFloat`, not through `f64`.** An earlier spelling wrote
/// `BigFloat::from_f64(2f64.powi(e))` and refused `|e| > 1000` — safe, because flushing a radius
/// to zero would read as a *tighter* interval than the realization earned. But a realization at
/// 1024 bits carries a radius near `2⁻¹⁰¹⁹`, so that refusal turned **more** precision into
/// "undecided": climbing a ladder made the answer worse, measured. A power of two is exact in
/// `BigFloat` at any precision, so the bound is assembled there and the range guard is gone.
///
/// `None` only when the exponent is absurd enough that the multiply loop would not terminate
/// usefully — far outside anything a realization produces.
fn rad_upper_bigfloat(rad: Mag, prec: usize) -> Option<BigFloat> {
    let Some(e) = rad.exp2() else {
        return Some(BigFloat::from_f64(0.0, prec)); // an exact realization: a zero radius is honest
    };
    if !(-1_000_000..=1_000_000).contains(&e) {
        return None;
    }
    let mut out = BigFloat::from_f64(1.0, prec);
    let mut rem = e;
    while rem != 0 {
        let step = rem.clamp(-1000, 1000);
        out = out.mul(
            &BigFloat::from_f64(2f64.powi(step as i32), prec),
            prec,
            HP_RM,
        );
        rem -= step;
    }
    Some(out)
}

/// **The `f64` nearest `x`, when that `f64` is normal** — `None` when it is not.
///
/// Two steps that must not be confused: `set_precision(53, ToEven)` performs the *rounding* (this
/// is the only place a value loses bits), and the assembly below is a pure re-encoding of what
/// that produced. Splitting them is why the caller can round two interval ends and compare.
///
/// Restricted to a normal nearest on purpose: an `f64` cannot hold a subnormal at full precision
/// nor an overflow at all, and `None` says so rather than handing back a rounded-away answer.
/// [`round_to_f64`] then decides at no precision, and [`ladder_f64`] decides such a value as `NaN`.
///
/// **The bottom edge is IEEE's, not 53 bits'.** Below `2⁻¹⁰²²` the `f64` grid is the subnormals',
/// so a value rounds up to `MIN_POSITIVE` from as low as that grid's midpoint `(2⁵³ − 1)·2⁻¹⁰⁷⁵` (a
/// tie, which goes to the even `MIN_POSITIVE`) — lower than 53-bit rounding reaches it. That band
/// is read off `x` itself, so this agrees with `is_normal` on [`nearest_f64_big`]'s answer for
/// every value. The top edge needs no such reading: both roundings overflow from the same point.
pub(super) fn bigfloat_to_f64(x: &BigFloat) -> Option<f64> {
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
    // same value is `e - 1 + 1023`. Outside the normal range an `f64` cannot say the value.
    let biased = i64::from(e) + 1022;
    if biased == 0 {
        // Rounded to 53 bits, `x` sits just under `2⁻¹⁰²²`; on the subnormal grid it may still
        // reach `MIN_POSITIVE`.
        let edge = BigFloat::from_f64(f64::MIN_POSITIVE, 64).mul(
            &BigFloat::from_f64(1.0 - f64::EPSILON / 2.0, 64),
            64,
            HP_RM,
        );
        // The magnitude is taken by hand: astro-float 0.9's `BigFloat::abs_cmp` compares the signed
        // values, so it reads every negative `x` as below the edge.
        let neg = sign == astro_float::Sign::Neg;
        let mag = if neg { x.neg() } else { x.clone() };
        let up = mag.cmp(&edge).is_some_and(|c| c >= 0);
        let mp = f64::MIN_POSITIVE;
        return up.then_some(if neg { -mp } else { mp });
    }
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
/// **The scale comes from the words being read, not from `as_raw_parts`' second field.**
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
/// **Not "away from zero", which is what this crate shipped first.** At `2⁵² + ½` both
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
/// integers, on the subnormal grid below the range — **and whether the rounding lost anything**.
///
/// [`nearest_f64`]'s unbounded twin: that one is the `Rat` road and caps at `u128`, this one takes
/// the `(numerator, common denominator)` pair [`MeetPoint::lift`] hands out, which is how a
/// coordinate too wide for `Rat` still gets a correctly rounded coordinate rather than a refusal.
/// `None` only where the value overflows; below the subnormals the nearest `f64` is `±0.0`.
///
/// **The second half is not optional information.** A three-plane meet is an exact rational,
/// which is a fact about the *realization*, not about the `f64` it is then read out as: a 59-bit
/// coordinate does not fit a 53-bit mantissa, so the readout rounds and the value handed over is
/// **not** the point. Reporting that as a zero error is the same lie the cache tells — measured on
/// the tilted-frame family, where the exact value is `0.130864196953086372` and its `f64` is
/// `0.13086419695308637578…`.
///
/// `true` means the `f64` **is** the rational, so a caller may honestly say its error is zero.
pub fn nearest_f64_big(num: &num_bigint::BigInt, den: &num_bigint::BigInt) -> Option<(f64, bool)> {
    if num.sign() == num_bigint::Sign::NoSign {
        return Some((0.0, true));
    }
    // `value ∈ (2^(e-1), 2^(e+1))`, so `52 - e` aims the mantissa at 53 bits and lands one short
    // at worst — one nudge, never two. Below the normal range the scale stops at `2¹⁰⁷⁴`, the
    // subnormal grid's: `q` is then the grid index itself, shorter than 53 bits and rounded once.
    // Rounding to 53 bits there and again onto the grid would round twice, and two roundings do not
    // compose into one.
    //
    // A carry out of the top bit needs no correction, measured. Rounding can push `q` to exactly
    // `2^53` (from `2^53 - ½`, the only way), or on the grid to `2^52`, and `q · 2^-k` is
    // representable exactly either way, so the wide `q` still spells the right `f64`. A branch for
    // it was written, planted against a 20,000-pair sweep plus a power-of-two boundary family, and
    // never changed an answer.
    const SUBNORMAL_SCALE: i64 = 1074;
    let e = num.magnitude().bits() as i64 - den.magnitude().bits() as i64;
    let mut k = (52 - e).min(SUBNORMAL_SCALE);
    let mut q = round_shifted(num, den, k);
    if q.magnitude().bits() < 53 && k < SUBNORMAL_SCALE {
        k += 1;
        q = round_shifted(num, den, k);
    }
    let m = i64::try_from(&q).ok()? as f64;
    // A value that rounds to zero keeps its sign, as IEEE's does: `q` is then the integer `0`.
    let m = if m == 0.0 && num.sign() == num_bigint::Sign::Minus {
        -0.0
    } else {
        m
    };
    // Scaled in steps, because `2f64.powi` underflows before the product does. `m` is up to ~2⁵³,
    // so `m · 2⁻ᵏ` can be an ordinary `f64` while `2⁻ᵏ` alone is zero: measured, one `powi` reads
    // `2⁻¹⁰⁰⁰` — a normal `f64` at 9.33e-302 — as **0.0**, and every subnormal too. A silent wrong
    // answer, not a refusal. Powers of two are exact, and `q · 2^-k` is an `f64` value whenever it
    // is finite, so the steps only move the exponent. The same trap waits at every `2f64.powi` of a
    // realization-sized exponent (`rad_upper_bigfloat`, `Realization::to_f64`'s error readout).
    let mut v = m;
    let mut rem = -k;
    while rem != 0 && v != 0.0 && v.is_finite() {
        let step = rem.clamp(-1000, 1000);
        v *= 2f64.powi(step as i32);
        rem -= step;
    }
    // An overflow is refused before the exactness test pays for a shift as wide as the input.
    if !v.is_finite() {
        return None;
    }
    // Exact iff the scaled division left no remainder: `value = q · 2^-k` exactly means
    // `num · 2^k == q · den`. Both sides are integers, so this is a decision and not an estimate —
    // and `v` is `q · 2^-k` itself, so a value rounded to `0.0` (`q = 0`) is never called exact.
    // `Realization::to_f64` turns `true` into a proven zero bound, which is why that matters; the
    // value itself is kept, because `0.0` is the nearest `f64` there
    // (`a_tiny_rational_still_names_its_f64`).
    let (lhs, rhs) = match k >= 0 {
        true => (num << k as usize, &q * den),
        false => (num.clone(), (&q * den) << (-k) as usize),
    };
    Some((v, lhs == rhs))
}

/// `places` decimal places of the **exact** rational `num/den` (`den` positive).
///
/// For a value that really is rational — a three-plane meet, a `QuadVal` whose radical vanishes
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
/// The caller's move on `None` is to realize again at higher precision, exactly as
/// `ladder_f64` climbs. That is why this reports undecided rather than picking:
/// a digit invented here would be indistinguishable, to everything downstream, from one the
/// definition actually determines.
pub fn round_to_digits(mid: &BigFloat, rad: Mag, places: usize) -> Option<String> {
    if mid.is_nan() || mid.is_inf() {
        return None;
    }
    // The ends `round_to_f64` forms, at `mid`'s own precision: one taken from the caller would round
    // `mid` first and, with a zero radius, print that rounding's digits as the value's.
    let (lo, hi) = interval_ends(mid, rad)?;
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
// that drive it.

/// **`√v` realized at `prec` bits, with its error** — `v · (1/√v)`, so the one radical primitive
/// this crate already has ([`inv_sqrt_hp`]) is the only place a square root is approached.
///
/// Crate-private: its consumer is [`quad_hp`] (measured); a radius stated as its square is
/// wide and takes [`sqrt_hp_big`]. `None` for a negative `v`.
pub(crate) fn sqrt_hp(v: Rat, prec: usize) -> Option<HpBounded> {
    if v < Rat::from_int(0) {
        return None;
    }
    if v == Rat::from_int(0) {
        return Some(HpBounded::exact(BigFloat::from_f64(0.0, prec)));
    }
    let inv = inv_sqrt_hp(v, prec)?;
    Some(HpBounded::of_rat(v, prec).mul(&inv, prec))
}

/// [`sqrt_hp`] for a wide radicand — `√(n/d) = n · (1/√(n·d))`, so the one `BigInt` radical
/// primitive ([`inv_sqrt_hp_big`]) is the only place a root is approached.
pub fn sqrt_hp_big(v: &BigRat, prec: usize) -> Option<HpBounded> {
    if v.is_negative() {
        return None;
    }
    if v.is_zero() {
        return Some(HpBounded::exact(BigFloat::from_f64(0.0, prec)));
    }
    let nd = v.numer() * v.denom();
    let inv = inv_sqrt_hp_big(&nd, prec)?;
    Some(HpBounded::of_bigint(v.numer(), prec).mul(&inv, prec))
}

/// **A quadratic algebraic scalar `a + b√c` realized at `prec` bits, with its error.**
///
/// Exact when the radical vanishes or resolves ([`quad::QuadVal::as_rat`]) — the value is asked,
/// not its provenance, so a tangency's rational root takes the rational road even though it
/// arrived through the same variant as an irrational one.
pub fn quad_hp(q: &quad::QuadVal, prec: usize) -> Option<HpBounded> {
    if let Some(r) = q.as_rat() {
        // The rational is the value; the only error is this realization's own rounding — charged
        // unconditionally rather than through `HpBounded::of_rat`'s exact branch, so that this
        // arm never hands an approached coordinate a zero radius (its consumer reads a zero as the
        // exact arm's answer).
        let v = rat_to_bigfloat(r, prec);
        return Some(HpBounded::new(v.clone(), HpBounded::round_off(&v, prec)));
    }
    let root = sqrt_hp(q.c(), prec)?;
    let term = HpBounded::of_rat(q.b(), prec).mul(&root, prec);
    Some(HpBounded::of_rat(q.a(), prec).add(&term, prec))
}

/// `base + dir·s` for exact rational `base`/`dir` and a realized `s` — the last step of a point
/// that lives at a parameter along an exactly-stated line.
pub fn affine_hp(base: Rat, dir: Rat, s: &HpBounded, prec: usize) -> Option<HpBounded> {
    let term = HpBounded::of_rat(dir, prec).mul(s, prec);
    Some(HpBounded::of_rat(base, prec).add(&term, prec))
}

/// **A point on a circle's `+ref` seam, realized at `prec` bits.**
///
/// `centre + r · e₁/|e₁|`, where `e₁` is the reference direction's component perpendicular to the
/// axis and the radius is stated as its square `r2` (the cylinder truth's form). `centre`, `e₁`
/// and `r2` are exact rationals. When `r2` is a rational's square — every radius a user writes —
/// the single irrational step is `1/|e₁|`, one `inv_sqrt_hp` with exact arithmetic around it,
/// the road this always took; otherwise `√r2` is realized beside it.
pub fn seam_point_hp(
    centre: [Rat; 3],
    e1: [Rat; 3],
    r2: &BigRat,
    prec: usize,
) -> Option<[HpBounded; 3]> {
    let mut sq = Rat::from_int(0);
    for c in &e1 {
        sq = sq.checked_add(c.checked_mul(*c)?)?;
    }
    let inv = inv_sqrt_hp(sq, prec)?;
    let scale = match sqrt_rat_big(r2) {
        Some(radius) => HpBounded::of_rat(radius, prec).mul(&inv, prec),
        None => sqrt_hp_big(r2, prec)?.mul(&inv, prec),
    };
    let coord = |k: usize| {
        let radial = HpBounded::of_rat(e1[k], prec).mul(&scale, prec);
        HpBounded::of_rat(centre[k], prec).add(&radial, prec)
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

#[cfg(test)]
#[path = "tests/decimal_realization.rs"]
mod decimal_realization;
