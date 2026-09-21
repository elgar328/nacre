use super::*;
/// The result of an orientation judgment — the shared sign vocabulary used by both
/// the 2D and 3D toleranced-sign judges (now in `nacre-judge`) and their downstream consumers.
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
    /// judgement that ran out of precision all reach it looking alike. `nacre_judge::Decision`
    /// is what keeps them apart — and what turns the last of them into a named reject rather
    /// than a silent merge.
    Zero,
}

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
