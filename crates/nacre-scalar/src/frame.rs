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

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic PRNG (splitmix64) so the random stress tests are reproducible.
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

    /// H4-soundness. The design's *tangential* directional tol (`tol_x = da·|y|`,
    /// `tol_y = da·|x|`) does **not** upper-bound the real f64 realization error
    /// near an axis (independent cos/sin rounding adds a radial error the
    /// tangential model misses), but the adopted coordinate-mixing bound
    /// `(|bx|+|by|)·da` per component does. Ground truth is astro-float at 256 bits.
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

    /// H4-soundness stress: over many random (base, angle) — including inexact
    /// rational degrees (large denominators) and a wide magnitude range (scale
    /// invariance) — the adopted tol upper-bounds the real f64 realization error
    /// on every sample. Ground truth is astro-float at 256 bits.
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
    /// high-precision judgment (256 bits). Over random triples the two paths agree
    /// on every one — the filter is sound, so the decision is precision-path
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
}
