//! The 2D-frame exact-sign judgment (overhaul §3).
//!
//! Rotating a point by a rational angle produces irrational coordinates
//! (cos/sin), so a rotated point cannot be held exactly. But the *sign* of an
//! `orient2d` determinant over such points can still be decided **soundly**
//! (never silently wrong): an f64 filter with a sound error bound handles the
//! easy cases, and the ambiguous ones **escalate** to arbitrary precision
//! (astro-float) from the points' exact definitions; a determinant below the
//! precision floor is a **declare-0** ([`Orient::Zero`]) — the collinear-within-
//! judgment case where the kernel (later) stops and asks the user (§6).
//!
//! [`Pt2`] here is the verified *minimal* point-construction vehicle (rotation
//! about the origin only). The real deliverable is [`orient2d_judge`], which is
//! general over any coordinate + tol; real sketch-frame placement (translation,
//! arbitrary centre) is kernel wiring.

use crate::{Angle, HP_RM, Rat};
use astro_float::BigFloat;

/// f64 trig-realization error bound used by the directional tol formula — a
/// conservative multiple of ulp covering cos/sin rounding, the deg→rad
/// conversion, and the combining arithmetic.
const DA_F64: f64 = 16.0 * f64::EPSILON;

/// Precision (bits) the orient2d judge escalates to before declaring 0.
const JUDGE_PREC: usize = 200;

/// A rational base point rotated about the origin by a rational `angle`. The
/// `base` + `angle` are the exact **definition** (never lost); `coord` is the f64
/// realization (a cache), and `tol` bounds its error by the **H4-soundness-
/// corrected coordinate-mixing formula** `(|bx|+|by|)·da` per component. (The
/// design's earlier tangential `da·|coord|` was refuted — it under-predicts near
/// an axis, where independent cos/sin rounding produces a *radial* error the
/// tangential bound misses.)
#[derive(Clone, Debug)]
pub struct Pt2 {
    pub base: (Rat, Rat),
    pub angle: Angle,
    pub coord: (f64, f64),
    pub tol: (f64, f64),
}

impl Pt2 {
    /// Rotate `base` about the origin by `angle`, realized in f64, with the
    /// coordinate-mixing tol `err = (|bx|+|by|)·da` on each component.
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
    /// two points with the same definition realize identically (path-independent).
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
    /// `declare-0` case where the kernel (later) stops and asks the user (§6).
    Zero,
}

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
/// the three point *definitions*, so it is path-independent.
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
    // zero/sign. Known upstream: issue stencillogic/astro-float#44 (open), fix in
    // PR #45 (open, not yet released as of 0.9.5). `is_positive`/`is_negative`
    // read the sign flag directly and sidestep it.
    if d.is_zero() || bf_mag(&d) <= floor {
        Orient::Zero
    } else if d.is_positive() {
        Orient::Positive
    } else {
        Orient::Negative
    }
}
