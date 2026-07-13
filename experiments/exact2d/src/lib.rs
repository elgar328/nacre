//! Exact 2D sandbox — overhaul checkpoint 1, hypotheses H1 and H1.5.
//!
//! Isolated non-member crate (see this repo's `docs/overhaul.md`). It exists to
//! *measure* the overhaul's novel mechanisms in 2D before any 3D port, not to be
//! production code.
//!
//! - [`Rat`] — rational scalar for exact user-input dimensions/angles (tol 0,
//!   `§4`). Fixed-width `Ratio<i128>` with **checked** arithmetic: overflow is a
//!   *signal* (the §4 downgrade trigger), not a panic or silent wrap.
//! - [`Angle`] — rational degrees, mod-360, exact accumulation/closure. cos/sin
//!   realize in f64 (working) and arbitrary precision (judgment); the 90°-family
//!   realizes to exact rationals so those rotations stay tol 0.
//!
//! H1.5 finding: the double-double `twofloat` failed the trig accuracy gate
//! (~1e-16 at zero-crossings), so the high-precision judgment realization uses
//! `astro-float` (arbitrary precision, ~1e-58 at 160 bits) instead.

use astro_float::{BigFloat, Consts, RoundingMode};
use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedMul, CheckedSub};
use std::cell::RefCell;
use std::f64::consts::PI;

/// H2/H3 — explicit reference-based sharing (identity by `Handle`, not coord).
pub mod share;

/// H5 — operation log (`Document = Vec<Op>` + derived model; replay determinism).
pub mod oplog;

/// Precision (bits) and rounding for the high-precision realization layer
/// (astro-float). ~160 bits ≈ 48 decimal digits — far below CAD tolerance and
/// dial-able if H4-amplification ever needs more (the ceiling twofloat lacked).
const HP_PREC: usize = 160;
const HP_RM: RoundingMode = RoundingMode::ToEven;

thread_local! {
    /// Transcendental-constant cache (π, …) for the high-precision layer.
    static HP_CONSTS: RefCell<Consts> = RefCell::new(Consts::new().expect("astro-float consts"));
}

/// A rational scalar (exact, tol 0). Arithmetic returns `None` on i128 overflow
/// so the caller sees the §4 downgrade trigger explicitly; the real kernel would
/// switch that value to f64/double-double and tag its `Origin` with the tol.
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

    /// Bit-width of the larger of |numer|, |denom| — the "size" §4/H1 watches
    /// for bit growth under chained rational arithmetic.
    pub fn bit_width(self) -> u32 {
        let n = self.numer().unsigned_abs();
        let d = self.denom().unsigned_abs();
        (128 - n.max(d).leading_zeros()).max(1)
    }
}

/// A *direction* angle in degrees, kept normalized to `[0, 360)` exactly
/// (rational). Accumulation is exact: turning by a rational angle repeatedly and
/// completing a full turn lands back on exactly `0` — no f64 drift (H1). `cos`/
/// `sin` realization crosses into f64 (deg→rad via π): the irrational-realization
/// boundary (§4), where the real kernel would use double-double for *judgment*
/// (H1.5/H4). The angle stays exact; only its realized coordinate carries tol.
///
/// This is the direction type. A multi-turn *amount* (helix pitch × turns,
/// revolve sweep > 360°) must preserve the turn count, so it belongs to a
/// separate *unnormalized* `Sweep` type — flat 2D sketches (H1) never need it, so
/// it is left as a documented companion, not built here.
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

    /// `(cos, sin)` realized in arbitrary precision at `prec` bits — the judgment
    /// path (§4/§6, H4). astro-float replaces twofloat here (H1.5 finding). Higher
    /// `prec` gives ground truth; the default [`HP_PREC`] gives the working
    /// judgment realization. (numer/denom pass through f64, exact for the small
    /// values used here; a general large-rational path would build from a string.)
    fn cos_sin_at(self, prec: usize) -> (BigFloat, BigFloat) {
        HP_CONSTS.with_borrow_mut(|cc| {
            let pi = cc.pi(prec, HP_RM);
            let d180 = BigFloat::from_f64(180.0, prec);
            let n = BigFloat::from_f64(self.0.numer() as f64, prec);
            let d = BigFloat::from_f64(self.0.denom() as f64, prec);
            let rad = n
                .div(&d, prec, HP_RM)
                .mul(&pi, prec, HP_RM)
                .div(&d180, prec, HP_RM);
            (rad.cos(prec, HP_RM, cc), rad.sin(prec, HP_RM, cc))
        })
    }

    /// cos realized at the default judgment precision.
    pub fn cos_hp(self) -> BigFloat {
        self.cos_sin_at(HP_PREC).0
    }

    /// sin realized at the default judgment precision.
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
}

/// f64 trig-realization error bound used by the directional tol formula — a
/// conservative multiple of ulp covering cos/sin rounding, the deg→rad
/// conversion, and the combining arithmetic.
const DA_F64: f64 = 16.0 * f64::EPSILON;

/// A rational base point rotated about the origin by a rational `angle`, realized
/// in f64 with the design's *tangential* directional tol (§TIP: `tol = da × moment
/// arm`, tangential; the base is exact so the transported term `|R|·old_tol` is
/// 0). H4-soundness tests whether this tol upper-bounds the real realization error.
#[derive(Clone, Debug)]
pub struct Pt2 {
    pub base: (Rat, Rat),
    pub angle: Angle,
    pub coord: (f64, f64),
    pub tol: (f64, f64),
}

impl Pt2 {
    /// Rotate `base` about the origin by `angle`, realized in f64. Tol per the
    /// **H4-soundness-corrected** coordinate-mixing formula: `err_x = |bx·δc −
    /// by·δs|`, so each component is bounded by `(|bx|+|by|)·da` (the tangential
    /// `da·|coord|` was refuted — it under-predicts near an axis, see the
    /// soundness tests).
    pub fn rotated_about_origin(base: (Rat, Rat), angle: Angle) -> Self {
        let (c, s) = (angle.cos(), angle.sin());
        let (bx, by) = (base.0.to_f64(), base.1.to_f64());
        let coord = (bx * c - by * s, bx * s + by * c);
        let mix = (bx.abs() + by.abs()) * DA_F64;
        let tol = (mix, mix);
        Pt2 {
            base,
            angle,
            coord,
            tol,
        }
    }

    /// The coordinate realized at `prec` bits (astro-float) from the point's
    /// **definition** (base rotated by angle) — the escalation/ground-truth
    /// realization. Because it is a function of the definition, not the f64 cache,
    /// two points with the same definition realize identically (H4-consistency).
    pub fn hp_coord(&self, prec: usize) -> (BigFloat, BigFloat) {
        let (c, s) = self.angle.cos_sin_at(prec);
        let bx = rat_to_big(self.base.0, prec);
        let by = rat_to_big(self.base.1, prec);
        let x = bx
            .mul(&c, prec, HP_RM)
            .sub(&by.mul(&s, prec, HP_RM), prec, HP_RM);
        let y = bx
            .mul(&s, prec, HP_RM)
            .add(&by.mul(&c, prec, HP_RM), prec, HP_RM);
        (x, y)
    }
}

/// A rational base point carried through a *chain* of rotations about the origin,
/// realized **incrementally** in f64 (each rotation operates on the already-
/// realized coordinate, not a lazy single realization). It tracks the transported
/// directional tol `|R|·old + new` for H4-soundness piece 2, and — for comparison
/// — a `tol_no_transport` that omits the `|R|·old` term, to show that term is
/// necessary (an existing error rotates with the point; a tol that does not rotate
/// with it under-predicts).
#[derive(Clone, Debug)]
pub struct PtChain {
    pub base: (Rat, Rat),
    pub angles: Vec<Angle>,
    pub coord: (f64, f64),
    pub tol: (f64, f64),
    pub tol_no_transport: (f64, f64),
}

impl PtChain {
    /// Start at an exact rational base (tol 0).
    pub fn new(base: (Rat, Rat)) -> Self {
        PtChain {
            base,
            angles: Vec::new(),
            coord: (base.0.to_f64(), base.1.to_f64()),
            tol: (0.0, 0.0),
            tol_no_transport: (0.0, 0.0),
        }
    }

    /// Rotate the current (already-realized) point about the origin by `a`.
    /// Updates the transported tol `|R|·old + new` and the comparison
    /// `old + new` (no transport).
    pub fn rotate(&mut self, a: Angle) {
        let (c, s) = (a.cos(), a.sin());
        let (x, y) = self.coord;
        self.coord = (x * c - y * s, x * s + y * c);
        let new_err = (x.abs() + y.abs()) * DA_F64;

        // (b) transport old tol via |R|, then add this step's new error (a).
        let (tx, ty) = self.tol;
        self.tol = (
            c.abs() * tx + s.abs() * ty + new_err,
            s.abs() * tx + c.abs() * ty + new_err,
        );

        // comparison: carry old tol un-rotated (no |R|) + new error.
        let (nx, ny) = self.tol_no_transport;
        self.tol_no_transport = (nx + new_err, ny + new_err);

        self.angles.push(a);
    }

    /// Ground-truth coordinate: the same rotation chain applied at `prec` bits
    /// (test-only). At high `prec` this is effectively the exact composition.
    #[cfg(test)]
    fn truth(&self, prec: usize) -> (BigFloat, BigFloat) {
        let mut x = rat_to_big(self.base.0, prec);
        let mut y = rat_to_big(self.base.1, prec);
        for a in &self.angles {
            let (c, s) = a.cos_sin_at(prec);
            let nx = x
                .mul(&c, prec, HP_RM)
                .sub(&y.mul(&s, prec, HP_RM), prec, HP_RM);
            let ny = x
                .mul(&s, prec, HP_RM)
                .add(&y.mul(&c, prec, HP_RM), prec, HP_RM);
            x = nx;
            y = ny;
        }
        (x, y)
    }
}

/// A rational as an arbitrary-precision float (exact for the small values here).
fn rat_to_big(r: Rat, prec: usize) -> BigFloat {
    BigFloat::from_f64(r.numer() as f64, prec).div(
        &BigFloat::from_f64(r.denom() as f64, prec),
        prec,
        HP_RM,
    )
}

/// Magnitude of a `BigFloat` as an f64 power of two (0 when exactly zero).
fn bf_mag(bf: &BigFloat) -> f64 {
    if bf.is_zero() {
        0.0
    } else {
        2f64.powi(bf.exponent().unwrap_or(0))
    }
}

/// The result of an orientation judgment (§6 TIP).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orient {
    Positive,
    Negative,
    /// Declared 0 — collinear within the escalation precision cap. This is the
    /// `declare-0 조건` where the kernel stops and would ask the user (§6).
    Zero,
}

/// Precision (bits) the orient2d judge escalates to before declaring 0.
const JUDGE_PREC: usize = 200;

/// f64 orient2d determinant `(b−a) × (c−a)`.
fn det_f64(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

/// Sound error bound on the f64 orient2d determinant from each point's
/// directional tol (interval-style propagation through `d1x·d2y − d1y·d2x`;
/// the design prefers this explicit bound over CGAL-style intervals for speed).
fn det_bound(
    a: (f64, f64),
    ta: (f64, f64),
    b: (f64, f64),
    tb: (f64, f64),
    c: (f64, f64),
    tc: (f64, f64),
) -> f64 {
    let (d1x, d1y) = (b.0 - a.0, b.1 - a.1);
    let (d2x, d2y) = (c.0 - a.0, c.1 - a.1);
    let (t1x, t1y) = (tb.0 + ta.0, tb.1 + ta.1); // tol on b−a
    let (t2x, t2y) = (tc.0 + ta.0, tc.1 + ta.1); // tol on c−a
    let e = d1x.abs() * t2y
        + d2y.abs() * t1x
        + t1x * t2y
        + d1y.abs() * t2x
        + d2x.abs() * t1y
        + t1y * t2x;
    let mag = d1x.abs() * d2y.abs() + d1y.abs() * d2x.abs();
    e + 8.0 * f64::EPSILON * mag // + f64 arithmetic rounding of the determinant
}

/// orient2d determinant realized at `prec` bits (astro-float) from the points'
/// definitions — path-independent (a function of the definitions, not the cache).
fn det_hp(pa: &Pt2, pb: &Pt2, pc: &Pt2, prec: usize) -> BigFloat {
    let a = pa.hp_coord(prec);
    let b = pb.hp_coord(prec);
    let c = pc.hp_coord(prec);
    let d1x = b.0.sub(&a.0, prec, HP_RM);
    let d1y = b.1.sub(&a.1, prec, HP_RM);
    let d2x = c.0.sub(&a.0, prec, HP_RM);
    let d2y = c.1.sub(&a.1, prec, HP_RM);
    d1x.mul(&d2y, prec, HP_RM)
        .sub(&d1y.mul(&d2x, prec, HP_RM), prec, HP_RM)
}

/// TIP orient2d: f64 filter (`|det| > tol-bound` → trust the sign), else escalate
/// to astro-float at [`JUDGE_PREC`]; if the determinant there is below the
/// precision floor it is `Zero` (declare-0). The decision is a pure function of
/// the three point *definitions*, so it is path-independent (H4-consistency).
pub fn orient2d_judge(pa: &Pt2, pb: &Pt2, pc: &Pt2) -> Orient {
    let (a, b, c) = (pa.coord, pb.coord, pc.coord);
    let det = det_f64(a, b, c);
    let bound = det_bound(a, pa.tol, b, pb.tol, c, pc.tol);
    if det > bound {
        return Orient::Positive;
    }
    if det < -bound {
        return Orient::Negative;
    }
    // f64 filter ambiguous → escalate.
    let d = det_hp(pa, pb, pc, JUDGE_PREC);
    let scale =
        a.0.abs()
            .max(a.1.abs())
            .max(b.0.abs())
            .max(b.1.abs())
            .max(c.0.abs())
            .max(c.1.abs())
            .max(1.0);
    let floor = 16.0 * scale * scale * 2f64.powi(-(JUDGE_PREC as i32));
    // Sign via `is_positive`/`is_negative`, NOT `d > BigFloat::from_f64(0.0)`:
    // astro-float 0.9.5's `cmp` mishandles a zero right-hand operand — any
    // positive value below 0.5 (exponent < 0) is ordered `Less` than zero, since
    // the comparison falls through to an exponent compare without special-casing
    // zero/sign. Known upstream: issue stencillogic/astro-float#44 (open),
    // fix in PR #45 (open, not yet released as of 0.9.5). `is_positive`/
    // `is_negative` read the sign flag directly and sidestep it.
    if d.is_zero() || bf_mag(&d) <= floor {
        Orient::Zero
    } else if d.is_positive() {
        Orient::Positive
    } else {
        Orient::Negative
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The core H1 property in miniature: exact rational accumulation does not
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

    /// The §4 "thin film" example: `1.1 × 7` must be exactly `7.7`. In rationals
    /// `11/10 × 7 = 77/10`; in f64 `1.1 * 7.0` is not `7.7`.
    #[test]
    fn one_point_one_times_seven_is_exact() {
        let a = Rat::new(11, 10).unwrap();
        let seven = Rat::from_int(7);
        assert_eq!(a.checked_mul(seven).unwrap(), Rat::new(77, 10).unwrap());

        assert_ne!(1.1_f64 * 7.0, 7.7); // f64 cannot represent 7.7 exactly
    }

    /// H1 measurement: with fixed-width i128, chained coprime-denominator
    /// accumulation *does* overflow (the finite-precision cliff §4 handles by
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
        // H1 measurement record (run with `-- --nocapture`).
        eprintln!(
            "[H1] coprime 1/p accumulation overflows at prime index {idx} (p={}); \
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

    /// H1 closure: a rational angle accumulates exactly, so a full turn lands back
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

    // ---- H1.5: high-precision trig accuracy gate (before H4) ----
    //
    // FINDING (fail-fast): twofloat 0.7 (double-double) FAILED this gate. Its π
    // constant is exact, but its trig is "preliminary": intrinsic cos²+sin²−1
    // ≈ 5.3e-22, and absolute error at zero-crossings (cos 90°) degraded to
    // ≈ 1.8e-16 — f64-level, exactly where the judgment path needs accuracy most
    // (sign decisions happen near degeneracies = zero-crossings). Replaced by
    // astro-float (arbitrary precision), measured below, which passes by orders.

    /// Base-2 exponent of `|a − b|` — the error magnitude as a power of two
    /// (`i32::MIN` when exactly equal; astro-float reports a zero's exponent as
    /// `Some(0)`, so guard `is_zero` explicitly). More negative = more accurate.
    fn hp_err_exp(a: &BigFloat, b: f64) -> i32 {
        let err = a.sub(&BigFloat::from_f64(b, HP_PREC), HP_PREC, HP_RM);
        if err.is_zero() {
            i32::MIN
        } else {
            err.exponent().unwrap_or(i32::MIN)
        }
    }

    /// H1.5 gate: astro-float cos/sin at rational-exact angles must be accurate
    /// to well under 1e-30. As a power of two, `2^-100 ≈ 7.9e-31 < 1e-30`, so the
    /// error exponent must be below -100.
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

    /// H1.5 data point (not a gate — H4 owns speed): rough cost of one
    /// high-precision cos realization, recorded for later comparison.
    #[test]
    fn high_precision_trig_speed_sample() {
        use std::time::Instant;
        let a = Angle::from_deg(Rat::from_int(37)).unwrap();
        let n = 2000;
        let t = Instant::now();
        let mut acc = 0i32;
        for _ in 0..n {
            acc = acc.wrapping_add(a.cos_hp().exponent().unwrap_or(0));
        }
        let per = t.elapsed().as_nanos() as f64 / f64::from(n);
        std::hint::black_box(acc);
        eprintln!("[H1.5] astro-float {HP_PREC}-bit cos ~ {per:.0} ns/call");
    }

    // ---- H4-soundness (new-error term): is the tol formula an upper bound? ----

    /// H4-soundness. Does the design's *tangential* directional tol
    /// (`tol_x = da·|y|`, `tol_y = da·|x|`) upper-bound the real f64 realization
    /// error, component-wise, for a point rotated about the origin? Ground truth
    /// is astro-float at 256 bits.
    ///
    /// FINDING: it does **not**. Near an axis one coordinate is tiny, so
    /// `da·|small|` under-predicts the *other* coordinate's ~ulp rounding error —
    /// the tangential model assumes a pure rotation, but independent cos/sin
    /// rounding adds a radial error the tangential term misses. The moment-arm-
    /// independent bound `(|bx|+|by|)·da` per component IS sound (each coordinate
    /// error mixes both base components: `err_x = |bx·δc − by·δs|`).
    #[test]
    fn h4_soundness_new_error_term() {
        const GT: usize = 256;
        // (bx_num, bx_den, by_num, by_den)
        let bases = [(1, 1, 0, 1), (3, 1, 4, 1), (5, 2, 7, 3), (1, 1, 1, 1)];
        // degrees (num, den): include near-axis (89.9°, 0.1°) and generic.
        let angles = [
            (30, 1),
            (45, 1),
            (60, 1),
            (37, 1),
            (899, 10),
            (1, 10),
            (3, 1),
        ];

        let mut tangential_sound = true;
        let mut cross_sound = true;
        let mut first_break = None;
        for (bxn, bxd, byn, byd) in bases {
            for (an, ad) in angles {
                let base = (Rat::new(bxn, bxd).unwrap(), Rat::new(byn, byd).unwrap());
                let angle = Angle::from_deg(Rat::new(an, ad).unwrap()).unwrap();
                let p = Pt2::rotated_about_origin(base, angle);
                let (xt, yt) = p.hp_coord(GT);
                let ex = BigFloat::from_f64(p.coord.0, GT).sub(&xt, GT, HP_RM).abs();
                let ey = BigFloat::from_f64(p.coord.1, GT).sub(&yt, GT, HP_RM).abs();

                // (1) old tangential prediction (da·|y|, da·|x|), computed locally.
                let px = BigFloat::from_f64(DA_F64 * p.coord.1.abs(), GT);
                let py = BigFloat::from_f64(DA_F64 * p.coord.0.abs(), GT);
                if ex > px || ey > py {
                    tangential_sound = false;
                    first_break.get_or_insert((bxn, bxd, byn, byd, an, ad));
                }

                // (2) the adopted coordinate-mixing tol `(|bx|+|by|)·da` = p.tol.
                let a = BigFloat::from_f64(p.tol.0, GT);
                if ex > a || ey > a {
                    cross_sound = false;
                }
            }
        }
        eprintln!("[H4-soundness] design tangential (da·|y|,da·|x|) sound: {tangential_sound}");
        eprintln!("[H4-soundness] cross-term (|bx|+|by|)·da sound       : {cross_sound}");
        if let Some(b) = first_break {
            eprintln!(
                "[H4-soundness] tangential first breaks at (bxn,bxd,byn,byd,degn,degd) = {b:?}"
            );
        }
        assert!(
            !tangential_sound,
            "expected the tangential formula to break near an axis"
        );
        assert!(cross_sound, "(|bx|+|by|)·da should bound every sample");
    }

    /// Deterministic PRNG (splitmix64) so the random stress test is reproducible.
    fn splitmix64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn rng_i128(state: &mut u64, lo: i128, hi: i128) -> i128 {
        let span = (hi - lo + 1) as u128;
        lo + (u128::from(splitmix64(state)) % span) as i128
    }

    /// H4-soundness stress test: over many random (base, angle) — including
    /// **inexact rational degrees** (large denominators, so the f64 realization
    /// also rounds the angle) — check which directional-tol formula upper-bounds
    /// the real f64 realization error, and how tight it is. Ground truth is
    /// astro-float at 256 bits. Candidates: T (tangential/design), X (base-mixing
    /// `(|bx|+|by|)·da`), and X+T (component sum).
    #[test]
    fn h4_soundness_random_stress() {
        const GT: usize = 256;
        const N: usize = 10_000;
        let mut st = 0xC0FF_EE00_1234_5678u64;
        let (mut t_bad, mut x_bad, mut xt_bad) = (0usize, 0usize, 0usize);
        let mut worst_xt_tightness = 0.0_f64;

        for _ in 0..N {
            // wide magnitude range (µm-scale to km-scale) to check scale invariance.
            let bx = Rat::new(
                rng_i128(&mut st, -100_000, 100_000),
                rng_i128(&mut st, 1, 100),
            )
            .unwrap();
            let by = Rat::new(
                rng_i128(&mut st, -100_000, 100_000),
                rng_i128(&mut st, 1, 100),
            )
            .unwrap();
            let deg = Rat::new(rng_i128(&mut st, 0, 360_000), rng_i128(&mut st, 1, 9973)).unwrap();
            let angle = Angle::from_deg(deg).unwrap();
            let p = Pt2::rotated_about_origin((bx, by), angle);

            let (xt, yt) = p.hp_coord(GT);
            let ex = BigFloat::from_f64(p.coord.0, GT).sub(&xt, GT, HP_RM).abs();
            let ey = BigFloat::from_f64(p.coord.1, GT).sub(&yt, GT, HP_RM).abs();

            // T: old tangential (da·|y|, da·|x|), computed locally.
            let (tf_x, tf_y) = (DA_F64 * p.coord.1.abs(), DA_F64 * p.coord.0.abs());
            if ex > BigFloat::from_f64(tf_x, GT) || ey > BigFloat::from_f64(tf_y, GT) {
                t_bad += 1;
            }

            // X: the adopted coordinate-mixing tol = p.tol (= (|bx|+|by|)·da).
            let xval = p.tol.0;
            let xbf = BigFloat::from_f64(xval, GT);
            if ex > xbf || ey > xbf {
                x_bad += 1;
            }

            // X + T componentwise
            let (xtx, xty) = (xval + tf_x, xval + tf_y);
            if ex > BigFloat::from_f64(xtx, GT) || ey > BigFloat::from_f64(xty, GT) {
                xt_bad += 1;
            }
            worst_xt_tightness = worst_xt_tightness
                .max(bf_mag(&ex) / xtx)
                .max(bf_mag(&ey) / xty);
        }

        eprintln!(
            "[H4-rand N={N}] violations — tangential T: {t_bad}, base-mixing X: {x_bad}, X+T: {xt_bad}"
        );
        eprintln!("[H4-rand] X+T worst tightness (actual/pred): {worst_xt_tightness:.3}");
        assert_eq!(xt_bad, 0, "X+T should be sound over all random samples");
    }

    /// H4-soundness piece 2 (transport term `|R|·old`): carry a point through a
    /// chain of rotations, realized incrementally, and check the transported tol
    /// `|R|·old + new` upper-bounds the real error — and that OMITTING `|R|`
    /// (carrying old tol un-rotated) breaks it, proving the term is necessary.
    #[test]
    fn h4_soundness_transport_term() {
        const GT: usize = 256;
        const N: usize = 5000;
        let mut st = 0xBEEF_1234_5678_9ABCu64;
        let (mut with_bad, mut without_bad) = (0usize, 0usize);
        let mut worst_with_tightness = 0.0_f64;

        for _ in 0..N {
            let base = (
                Rat::new(rng_i128(&mut st, -1000, 1000), rng_i128(&mut st, 1, 50)).unwrap(),
                Rat::new(rng_i128(&mut st, -1000, 1000), rng_i128(&mut st, 1, 50)).unwrap(),
            );
            let mut p = PtChain::new(base);
            let steps = 2 + rng_i128(&mut st, 0, 2) as usize; // 2..4 rotations
            for _ in 0..steps {
                let deg =
                    Rat::new(rng_i128(&mut st, 0, 360_000), rng_i128(&mut st, 1, 997)).unwrap();
                p.rotate(Angle::from_deg(deg).unwrap());
            }

            let (xt, yt) = p.truth(GT);
            let ex = BigFloat::from_f64(p.coord.0, GT).sub(&xt, GT, HP_RM).abs();
            let ey = BigFloat::from_f64(p.coord.1, GT).sub(&yt, GT, HP_RM).abs();

            // with |R| transport
            if ex > BigFloat::from_f64(p.tol.0, GT) || ey > BigFloat::from_f64(p.tol.1, GT) {
                with_bad += 1;
            }
            // without transport (old tol carried un-rotated)
            if ex > BigFloat::from_f64(p.tol_no_transport.0, GT)
                || ey > BigFloat::from_f64(p.tol_no_transport.1, GT)
            {
                without_bad += 1;
            }
            if p.tol.0 > 0.0 && p.tol.1 > 0.0 {
                worst_with_tightness = worst_with_tightness
                    .max(bf_mag(&ex) / p.tol.0)
                    .max(bf_mag(&ey) / p.tol.1);
            }
        }

        eprintln!(
            "[H4-transport N={N}] |R| transport: {with_bad} violations; WITHOUT transport: {without_bad}"
        );
        eprintln!("[H4-transport] with-transport worst tightness: {worst_with_tightness:.3}");
        assert_eq!(
            with_bad, 0,
            "transported tol |R|·old + new should bound every chain"
        );
        // NOTE: `without_bad` is typically 0 too — v1's symmetric discretization tol
        // keeps transport's effect within the (~3×) margin, so random chains do not
        // expose the term's necessity. Asymmetric tol does (next test).
        eprintln!(
            "[H4-transport] (without-transport violations on symmetric-tol chains: {without_bad})"
        );
    }

    /// H4-soundness piece 2, necessity: with **asymmetric** tol the `|R|·old`
    /// transport term is required. A point whose error/tol lies entirely in x,
    /// rotated 90°, has its error move to y — transport follows it, omitting `|R|`
    /// leaves `tol_y = 0` and under-predicts. (Asymmetric tol arises v1-후: √
    /// distances, constraint-solver solutions.)
    #[test]
    fn h4_transport_is_necessary_for_asymmetric_tol() {
        // error e1 and tol1 entirely in x (T ≥ E > 0).
        let (e, t) = (1e-14_f64, 2e-14_f64);
        // 90° rotation matrix: cos = 0, sin = 1.
        let (c, s) = (0.0_f64, 1.0_f64);

        // actual error transports exactly: e2 = R·e1 = (0, e).
        let e2 = (c * e - s * 0.0, s * e + c * 0.0);
        // with transport: |R|·tol1 = (0, t).
        let tol_with = (c.abs() * t + s.abs() * 0.0, s.abs() * t + c.abs() * 0.0);
        // without transport: tol1 un-rotated = (t, 0).
        let tol_without = (t, 0.0);

        assert!(e2.1.abs() <= tol_with.1, "transport bounds the moved error");
        assert!(
            e2.1.abs() > tol_without.1,
            "without transport, tol_y = 0 under-predicts the error that moved to y"
        );
    }

    // ---- H4-consistency: path-independent orient2d decisions ----

    /// A rational base point rotated about the origin by `deg`.
    fn pt(bx: (i128, i128), by: (i128, i128), deg: (i128, i128)) -> Pt2 {
        Pt2::rotated_about_origin(
            (Rat::new(bx.0, bx.1).unwrap(), Rat::new(by.0, by.1).unwrap()),
            Angle::from_deg(Rat::new(deg.0, deg.1).unwrap()).unwrap(),
        )
    }

    /// orient2d decided purely at `prec` bits (a from-scratch high-precision
    /// reference, no f64 filter) — for checking the filter never disagrees.
    fn orient_at_prec(pa: &Pt2, pb: &Pt2, pc: &Pt2, prec: usize) -> Orient {
        let d = det_hp(pa, pb, pc, prec);
        let m = |p: &Pt2| p.coord.0.abs().max(p.coord.1.abs());
        let scale = m(pa).max(m(pb)).max(m(pc)).max(1.0);
        let floor = 16.0 * scale * scale * 2f64.powi(-(prec as i32));
        // sign via is_positive, not `> 0` — astro-float#44 (see orient2d_judge).
        if d.is_zero() || bf_mag(&d) <= floor {
            Orient::Zero
        } else if d.is_positive() {
            Orient::Positive
        } else {
            Orient::Negative
        }
    }

    /// H4-consistency: rotation is orientation-preserving, so the *same* triangle
    /// realized at *different* rotation angles must yield the *same* orient2d sign
    /// — the decision cannot depend on the realization path.
    #[test]
    fn orient2d_sign_is_preserved_under_shared_rotation() {
        // CCW triangle (0,0),(4,0),(0,3) → Positive at every angle.
        let bases = [((0, 1), (0, 1)), ((4, 1), (0, 1)), ((0, 1), (3, 1))];
        for (dn, dd) in [(0, 1), (30, 1), (37, 1), (123, 1), (1, 3), (359, 1)] {
            let p: Vec<Pt2> = bases.iter().map(|&(bx, by)| pt(bx, by, (dn, dd))).collect();
            assert_eq!(
                orient2d_judge(&p[0], &p[1], &p[2]),
                Orient::Positive,
                "CCW at {dn}/{dd} deg"
            );
            // reversed vertex order flips the sign, also at every angle.
            assert_eq!(
                orient2d_judge(&p[0], &p[2], &p[1]),
                Orient::Negative,
                "CW at {dn}/{dd} deg"
            );
        }
    }

    /// H4-consistency: exactly collinear points stay collinear under rotation, so
    /// the judge must declare `Zero` at every angle (never a spurious sign).
    #[test]
    fn collinear_is_zero_under_shared_rotation() {
        // (0,0),(1,1),(2,2) are collinear.
        let bases = [((0, 1), (0, 1)), ((1, 1), (1, 1)), ((2, 1), (2, 1))];
        for (dn, dd) in [(0, 1), (30, 1), (37, 1), (1, 7), (250, 1)] {
            let p: Vec<Pt2> = bases.iter().map(|&(bx, by)| pt(bx, by, (dn, dd))).collect();
            assert_eq!(
                orient2d_judge(&p[0], &p[1], &p[2]),
                Orient::Zero,
                "collinear at {dn}/{dd} deg"
            );
        }
    }

    /// H4-consistency: the f64 filter must never disagree with a from-scratch
    /// high-precision judgment (256 bits). Over random triples, the two paths
    /// agree on every one — the filter is sound, so the decision is precision-path
    /// independent.
    #[test]
    fn filter_agrees_with_full_precision_over_random() {
        const N: usize = 2000;
        let mut st = 0xF117_0B7E_C0DE_9999u64;
        let mut disagreements = 0usize;
        let mut escalated = 0usize;
        for _ in 0..N {
            let mk = |st: &mut u64| {
                let b = (
                    Rat::new(rng_i128(st, -100, 100), rng_i128(st, 1, 20)).unwrap(),
                    Rat::new(rng_i128(st, -100, 100), rng_i128(st, 1, 20)).unwrap(),
                );
                let a = Angle::from_deg(
                    Rat::new(rng_i128(st, 0, 360_000), rng_i128(st, 1, 997)).unwrap(),
                )
                .unwrap();
                Pt2::rotated_about_origin(b, a)
            };
            let (pa, pb, pc) = (mk(&mut st), mk(&mut st), mk(&mut st));
            let judged = orient2d_judge(&pa, &pb, &pc);
            let reference = orient_at_prec(&pa, &pb, &pc, 256);
            if judged != reference {
                disagreements += 1;
            }
            // count how often the f64 filter was insufficient (for context).
            let det = det_f64(pa.coord, pb.coord, pc.coord);
            let bound = det_bound(pa.coord, pa.tol, pb.coord, pb.tol, pc.coord, pc.tol);
            if det.abs() <= bound {
                escalated += 1;
            }
        }
        eprintln!(
            "[H4-consistency N={N}] filter↔full-precision disagreements: {disagreements}; escalated: {escalated}"
        );
        assert_eq!(
            disagreements, 0,
            "f64 filter must never disagree with high-precision judgment"
        );
    }

    // ---- H4-amplification: does tol approach the 1e-9 judgment tolerance? ----

    /// H4-amplification.
    ///
    /// (A) BUNDLED path (§4 — the path the kernel uses): a run of same-axis
    /// rational rotations accumulates into one angle, realized **once**. The tol
    /// is then `(|bx|+|by|)·da` — linear in coordinate scale, independent of how
    /// many rotations were accumulated. Even at 1 km the f64 tol is ~6e-9, and
    /// escalation to a modest astro-float precision (≪ our 160 bits) drives it far
    /// below 1e-9. So amplification never threatens correctness (dial-able
    /// precision); it only sets filter success rate.
    ///
    /// (B) INCREMENTAL (un-bundled) realization: the transport term `|R|·old`
    /// uses the component-wise **absolute** matrix, whose row sum `|cos|+|sin| ≥ 1`
    /// multiplies the tol every step → the bound grows **exponentially** in the
    /// step count, even though the *true* error (transported by the norm-
    /// preserving rotation `R`) grows only linearly. So the bound stays sound but
    /// becomes uselessly loose. ⟹ §4 bundling is **essential, not an optimization**.
    #[test]
    fn h4_amplification() {
        // (A) bundled: tol vs scale, N-independent.
        let target = 1e-15;
        let mut worst_prec = 0.0_f64;
        for &s in &[1.0_f64, 1e3, 1e6] {
            let si = s as i128;
            let p = Pt2::rotated_about_origin(
                (Rat::new(si, 1).unwrap(), Rat::new(si * 7, 10).unwrap()),
                Angle::from_deg(Rat::from_int(37)).unwrap(),
            );
            let tol = p.tol.0;
            let prec = if tol > target {
                52.0 + (tol / target).log2()
            } else {
                52.0
            };
            worst_prec = worst_prec.max(prec);
            eprintln!(
                "[H4-amp bundled] scale={s:>8.0}mm  f64 tol={tol:.2e}  → escalate ~{prec:.0}b for <{target:e}"
            );
        }
        eprintln!(
            "[H4-amp] worst bundled precision need: ~{worst_prec:.0} bits (we use {HP_PREC})"
        );

        // (B) incremental un-bundled: exponential bound growth vs linear true error.
        let mut p = PtChain::new((Rat::from_int(1000), Rat::from_int(700)));
        let base_step = 1700.0 * DA_F64; // ≈ per-step new error (linear rate)
        for n in 1..=30usize {
            p.rotate(Angle::from_deg(Rat::from_int(37)).unwrap());
            if [1, 10, 20, 30].contains(&n) {
                let (xt, yt) = p.truth(256);
                let ex = bf_mag(
                    &BigFloat::from_f64(p.coord.0, 256)
                        .sub(&xt, 256, HP_RM)
                        .abs(),
                );
                let ey = bf_mag(
                    &BigFloat::from_f64(p.coord.1, 256)
                        .sub(&yt, 256, HP_RM)
                        .abs(),
                );
                let actual = ex.max(ey);
                eprintln!(
                    "[H4-amp incr] N={n:>2}  bound tol={:.2e}  actual err={actual:.2e}  linear~{:.2e}",
                    p.tol.0,
                    f64::from(n as u32) * base_step
                );
            }
        }

        // (A) bundled worst combo is well within our escalation precision.
        assert!(
            worst_prec < HP_PREC as f64,
            "astro-float {HP_PREC}b covers the worst bundled amplification combo"
        );
        // (B) the un-bundled bound blew up far past a linear rate (≫), confirming
        // bundling is required. (30 steps of ×~1.4 ⇒ ~1e4× a linear bound.)
        assert!(
            p.tol.0 > 1000.0 * f64::from(30u32) * base_step,
            "incremental un-bundled tol should blow up super-linearly"
        );
    }

    // ---- H4-speed: per-judgment cost, escalation frequency, filter success ----

    /// H4-speed. Times the two orient2d paths — the f64 filter (common) and the
    /// astro-float escalation (rare) — and measures how often random geometry
    /// needs escalation. Escalation is far slower but only fires for genuinely
    /// near-degenerate configs, so the amortized cost tracks the fast filter.
    #[test]
    fn h4_speed() {
        use std::time::Instant;
        let mut st = 0x5EED_1234_ABCD_0001u64;
        let rnd = |st: &mut u64| {
            Pt2::rotated_about_origin(
                (
                    Rat::new(rng_i128(st, -100, 100), rng_i128(st, 1, 20)).unwrap(),
                    Rat::new(rng_i128(st, -100, 100), rng_i128(st, 1, 20)).unwrap(),
                ),
                Angle::from_deg(Rat::new(rng_i128(st, 0, 360_000), rng_i128(st, 1, 997)).unwrap())
                    .unwrap(),
            )
        };

        // (1) filter path: random well-separated triples.
        const N: usize = 5000;
        let tris: Vec<(Pt2, Pt2, Pt2)> = (0..N)
            .map(|_| (rnd(&mut st), rnd(&mut st), rnd(&mut st)))
            .collect();
        let mut escalations = 0usize;
        for (a, b, c) in &tris {
            let det = det_f64(a.coord, b.coord, c.coord);
            let bound = det_bound(a.coord, a.tol, b.coord, b.tol, c.coord, c.tol);
            if det.abs() <= bound {
                escalations += 1;
            }
        }
        let t = Instant::now();
        for (a, b, c) in &tris {
            std::hint::black_box(orient2d_judge(a, b, c));
        }
        let filter_ns = t.elapsed().as_nanos() as f64 / N as f64;

        // (2) escalation path: exactly collinear triples (force astro-float).
        const M: usize = 100;
        let col: Vec<(Pt2, Pt2, Pt2)> = (0..M)
            .map(|i| {
                let k = i as i128 + 1;
                // collinear points rotated 37° — still collinear, but the
                // escalation now realizes real cos/sin (not the cos 0° fast path).
                let p = |m: i128| {
                    Pt2::rotated_about_origin(
                        (Rat::new(m * k, 1).unwrap(), Rat::new(m * k, 1).unwrap()),
                        Angle::from_deg(Rat::from_int(37)).unwrap(),
                    )
                };
                (p(0), p(1), p(2))
            })
            .collect();
        assert_eq!(
            orient2d_judge(&col[0].0, &col[0].1, &col[0].2),
            Orient::Zero,
            "collinear must escalate to Zero"
        );
        let t2 = Instant::now();
        for (a, b, c) in &col {
            std::hint::black_box(orient2d_judge(a, b, c));
        }
        let esc_us = t2.elapsed().as_micros() as f64 / M as f64;

        eprintln!(
            "[H4-speed] filter path (random): {filter_ns:.0} ns/judgment; escalations {escalations}/{N}"
        );
        eprintln!(
            "[H4-speed] escalation path (collinear): {esc_us:.0} µs/judgment (un-cached astro-float; a real kernel caches per-angle realizations)"
        );
    }
}
