//! 3D exactness experiment (overhaul stage 0.5) — the go/no-go on the 3D TIP math.
//!
//! The 2D experiment (`experiments/exact2d`, ported to `nacre-scalar::frame`)
//! validated the *foundation*: a rotated point cannot be held exactly (cos/sin are
//! irrational), but the *sign* of an `orient2d` determinant over such points is
//! decided **soundly** by an f64 filter with a sound error bound, escalating the
//! ambiguous cases to arbitrary precision (astro-float). What 2D could not see is
//! **3D-specific** (overhaul.md §11-B): a 3×3 `orient3d` determinant, the implicit
//! point where three rotated planes meet, direction-wise (xyz) tol propagation, and
//! rotation-chain composition across different axes.
//!
//! This crate builds `Pt3` + `orient3d_judge` — the 3D analogue of `frame::Pt2` /
//! `orient2d_judge` — and stress-tests the new determinant error bound against an
//! astro-float ground truth. `Rat`/`Angle` are reused from the verified
//! `nacre-scalar`; the arbitrary-precision trig is reimplemented here (nacre-scalar's
//! `cos_sin_at` is `pub(crate)`), production unchanged.
//!
//! **Ground-truth trap.** The truth is the astro-float realization of the *exact
//! definition* (rational angle at high precision), never the exact sign of the
//! already-approximated f64 coefficients — that would confirm f64≈f64 and validate
//! nothing about rotation soundness.

use astro_float::{BigFloat, Consts, RoundingMode};
use nacre_scalar::{Angle, Rat};
use std::cell::RefCell;

/// Rounding mode for the high-precision layer (matches nacre-scalar).
const HP_RM: RoundingMode = RoundingMode::ToEven;

/// f64 trig-realization error bound (per the coordinate-mixing tol formula) — a
/// conservative ulp multiple covering cos/sin rounding, deg→rad, and combining.
const DA_F64: f64 = 16.0 * f64::EPSILON;

/// Precision (bits) the orient3d judge escalates to before declaring 0.
const JUDGE_PREC: usize = 200;

thread_local! {
    static HP_CONSTS: RefCell<Consts> = RefCell::new(Consts::new().expect("astro-float consts"));
}

/// A rational as an arbitrary-precision float (exact for the small values here).
fn rat_to_big(r: Rat, prec: usize) -> BigFloat {
    BigFloat::from_f64(r.numer() as f64, prec).div(
        &BigFloat::from_f64(r.denom() as f64, prec),
        prec,
        HP_RM,
    )
}

/// `(cos, sin)` of `angle` at `prec` bits (astro-float). Reimplements
/// `nacre-scalar::Angle::cos_sin_at` (which is `pub(crate)`) so the experiment can
/// dial the precision (ground truth vs judgment). `rad = deg·π/180`.
fn cos_sin_hp(angle: Angle, prec: usize) -> (BigFloat, BigFloat) {
    let deg = angle.deg();
    HP_CONSTS.with_borrow_mut(|cc| {
        let pi = cc.pi(prec, HP_RM);
        let d180 = BigFloat::from_f64(180.0, prec);
        let n = BigFloat::from_f64(deg.numer() as f64, prec);
        let d = BigFloat::from_f64(deg.denom() as f64, prec);
        let rad = n
            .div(&d, prec, HP_RM)
            .mul(&pi, prec, HP_RM)
            .div(&d180, prec, HP_RM);
        (rad.cos(prec, HP_RM, cc), rad.sin(prec, HP_RM, cc))
    })
}

/// Magnitude of a `BigFloat` as an f64 power of two (0 when exactly zero).
fn bf_mag(bf: &BigFloat) -> f64 {
    if bf.is_zero() {
        0.0
    } else {
        2f64.powi(bf.exponent().unwrap_or(0))
    }
}

/// A coordinate axis — a rotation node turns a point about one of these. (Arbitrary
/// rational axes via Rodrigues are a later extension; axis-aligned rotations already
/// exercise per-axis heterogeneous tol and, chained across different axes, the
/// axis-change composition — enough for the 3D determinant bound.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

/// One rotation in a point's definition: turn about `axis` by the rational `angle`.
#[derive(Clone, Copy, Debug)]
pub struct RotNode {
    pub axis: Axis,
    pub angle: Angle,
}

/// A rational base point carried through a chain of axis rotations (§TIP ⑦ rotation
/// history). `base` + `chain` are the exact **definition** (never lost); `coord` is
/// the f64 realization (a cache), and `tol` bounds its error as a **direction-wise
/// xyz vector** (§TIP ⑤). `hp_coord` realizes the chain at arbitrary precision from
/// the definition, so two points with the same definition realize identically
/// (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct Pt3 {
    pub base: [Rat; 3],
    pub chain: Vec<RotNode>,
    pub coord: [f64; 3],
    pub tol: [f64; 3],
}

/// The two in-plane axis indices for a rotation about `axis` (the third is fixed).
fn plane_of(axis: Axis) -> (usize, usize) {
    match axis {
        Axis::X => (1, 2), // rotate y,z
        Axis::Y => (2, 0), // rotate z,x
        Axis::Z => (0, 1), // rotate x,y
    }
}

impl Pt3 {
    /// A point at `base` with no rotation (coord = base, tol 0).
    pub fn at(base: [Rat; 3]) -> Self {
        Pt3 {
            coord: [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
            tol: [0.0; 3],
            base,
            chain: Vec::new(),
        }
    }

    /// Extend the definition by a rotation about `axis` by `angle`, updating the f64
    /// cache and its tol. Same-axis 90°-family angles rotate exactly (tol 0, Niven);
    /// otherwise the two in-plane coords gain the coordinate-mixing realization error
    /// `(|u|+|v|)·da`, and every existing tol is transported by the component-wise
    /// absolute rotation `|R|` (§TIP: `new tol = |R|·old + mix`).
    pub fn rotate(mut self, axis: Axis, angle: Angle) -> Self {
        let (i, j) = plane_of(axis);
        let (u, v) = (self.coord[i], self.coord[j]);
        // cos/sin — exact (rational) for the 90°-family, else f64 (with realization tol).
        let (c, s, exact) = match angle.try_exact_cos_sin() {
            Some((cr, sr)) => (cr.to_f64(), sr.to_f64(), true),
            None => (angle.cos(), angle.sin(), false),
        };
        self.coord[i] = u * c - v * s;
        self.coord[j] = u * s + v * c;
        let mix = if exact {
            0.0
        } else {
            (u.abs() + v.abs()) * DA_F64
        };
        let (ti, tj) = (self.tol[i], self.tol[j]);
        self.tol[i] = c.abs() * ti + s.abs() * tj + mix;
        self.tol[j] = s.abs() * ti + c.abs() * tj + mix;
        self.chain.push(RotNode { axis, angle });
        self
    }

    /// The coordinate realized at `prec` bits from the **definition** (base rotated
    /// through the chain) — path-independent ground truth / escalation realization.
    pub fn hp_coord(&self, prec: usize) -> [BigFloat; 3] {
        let mut p = [
            rat_to_big(self.base[0], prec),
            rat_to_big(self.base[1], prec),
            rat_to_big(self.base[2], prec),
        ];
        for node in &self.chain {
            let (i, j) = plane_of(node.axis);
            let (c, s) = cos_sin_hp(node.angle, prec);
            let u = p[i].clone();
            let v = p[j].clone();
            p[i] = u
                .mul(&c, prec, HP_RM)
                .sub(&v.mul(&s, prec, HP_RM), prec, HP_RM);
            p[j] = u
                .mul(&s, prec, HP_RM)
                .add(&v.mul(&c, prec, HP_RM), prec, HP_RM);
        }
        p
    }
}

/// The result of an orientation judgment (§6 TIP).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orient {
    Positive,
    Negative,
    /// Declared 0 — coplanar within the escalation precision cap; the §6 "ask the
    /// user" case.
    Zero,
}

/// The three edge rows of `orient3d(a,b,c,d) = det[a−d, b−d, c−d]`.
fn rows(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] - d[0], a[1] - d[1], a[2] - d[2]],
        [b[0] - d[0], b[1] - d[1], b[2] - d[2]],
        [c[0] - d[0], c[1] - d[1], c[2] - d[2]],
    ]
}

/// f64 `orient3d` determinant `det[a−d, b−d, c−d]`.
fn det3_f64(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let r = rows(a, b, c, d);
    r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
}

/// Error radius of a product `a·b·c` given each factor's value and tol (≥0):
/// `Π(|v|+τ) − Π|v|` — the interval product radius (worst-case, sound).
fn prod_err(va: f64, ta: f64, vb: f64, tb: f64, vc: f64, tc: f64) -> f64 {
    (va.abs() + ta) * (vb.abs() + tb) * (vc.abs() + tc) - va.abs() * vb.abs() * vc.abs()
}

/// Sound error bound on the f64 `orient3d` determinant from each point's directional
/// tol. The determinant is six signed triple-products of the edge entries; each
/// entry `(a−d)[k]` carries tol `tol_a[k] + tol_d[k]`. The bound sums the six product
/// radii (input-tol propagation, triangle-inequality worst case) plus a term for the
/// f64 rounding of the determinant's own arithmetic. This is the 3D analogue of
/// `nacre-scalar::frame::det_bound` and the experiment's new-math (§TIP ②).
fn det3_bound(p: [[f64; 3]; 4], t: [[f64; 3]; 4]) -> f64 {
    // p/t are points a,b,c,d with their tols; d = p[3] is the apex.
    let r = rows(p[0], p[1], p[2], p[3]);
    let td = t[3]; // apex tol adds to every edge on subtraction
    let t = [
        [t[0][0] + td[0], t[0][1] + td[1], t[0][2] + td[2]],
        [t[1][0] + td[0], t[1][1] + td[1], t[1][2] + td[2]],
        [t[2][0] + td[0], t[2][1] + td[1], t[2][2] + td[2]],
    ];
    // The six signed triple products of the cofactor expansion (indices into r/t).
    // (row, col) triples with sign folded into the sum below.
    let terms = [
        [(0, 0), (1, 1), (2, 2)],
        [(0, 0), (1, 2), (2, 1)],
        [(0, 1), (1, 0), (2, 2)],
        [(0, 1), (1, 2), (2, 0)],
        [(0, 2), (1, 0), (2, 1)],
        [(0, 2), (1, 1), (2, 0)],
    ];
    let mut input_tol = 0.0;
    let mut mag = 0.0;
    for tri in terms {
        let [(r0, c0), (r1, c1), (r2, c2)] = tri;
        input_tol += prod_err(
            r[r0][c0], t[r0][c0], r[r1][c1], t[r1][c1], r[r2][c2], t[r2][c2],
        );
        mag += r[r0][c0].abs() * r[r1][c1].abs() * r[r2][c2].abs();
    }
    // f64 rounding of the ~17-op determinant, bounded generously.
    input_tol + 16.0 * f64::EPSILON * mag
}

/// `orient3d` determinant realized at `prec` bits (astro-float) from the point
/// definitions — path-independent ground truth.
fn det3_hp(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3, prec: usize) -> BigFloat {
    let a = pa.hp_coord(prec);
    let b = pb.hp_coord(prec);
    let c = pc.hp_coord(prec);
    let d = pd.hp_coord(prec);
    let sub = |x: &BigFloat, y: &BigFloat| x.sub(y, prec, HP_RM);
    let mul = |x: &BigFloat, y: &BigFloat| x.mul(y, prec, HP_RM);
    let r = [
        [sub(&a[0], &d[0]), sub(&a[1], &d[1]), sub(&a[2], &d[2])],
        [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
        [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
    ];
    let m0 = mul(&r[1][1], &r[2][2]).sub(&mul(&r[1][2], &r[2][1]), prec, HP_RM);
    let m1 = mul(&r[1][0], &r[2][2]).sub(&mul(&r[1][2], &r[2][0]), prec, HP_RM);
    let m2 = mul(&r[1][0], &r[2][1]).sub(&mul(&r[1][1], &r[2][0]), prec, HP_RM);
    mul(&r[0][0], &m0)
        .sub(&mul(&r[0][1], &m1), prec, HP_RM)
        .add(&mul(&r[0][2], &m2), prec, HP_RM)
}

/// The magnitude scale of four points (for the declare-0 floor).
fn scale4(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let m = |p: [f64; 3]| p[0].abs().max(p[1].abs()).max(p[2].abs());
    m(a).max(m(b)).max(m(c)).max(m(d)).max(1.0)
}

/// TIP `orient3d`: f64 filter (`|det| > bound` → trust the sign), else escalate to
/// astro-float at [`JUDGE_PREC`]; a determinant below the precision floor (`~scale³`,
/// cubic for the 3×3 case) is `Zero` (declare-0). Path-independent (a function of the
/// four point definitions).
pub fn orient3d_judge(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3) -> Orient {
    let (a, b, c, d) = (pa.coord, pb.coord, pc.coord, pd.coord);
    let det = det3_f64(a, b, c, d);
    let bound = det3_bound([a, b, c, d], [pa.tol, pb.tol, pc.tol, pd.tol]);
    if det > bound {
        return Orient::Positive;
    }
    if det < -bound {
        return Orient::Negative;
    }
    let dh = det3_hp(pa, pb, pc, pd, JUDGE_PREC);
    let scale = scale4(a, b, c, d);
    let floor = 16.0 * scale * scale * scale * 2f64.powi(-(JUDGE_PREC as i32));
    // Sign via is_positive (astro-float#44 workaround — see nacre-scalar::frame).
    if dh.is_zero() || bf_mag(&dh) <= floor {
        Orient::Zero
    } else if dh.is_positive() {
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

    fn ri(n: i128, d: i128) -> Rat {
        Rat::new(n, d).unwrap()
    }

    fn axis_of(k: i128) -> Axis {
        match k.rem_euclid(3) {
            0 => Axis::X,
            1 => Axis::Y,
            _ => Axis::Z,
        }
    }

    /// A random rational base point (wide magnitude for scale invariance).
    fn rand_base(st: &mut u64) -> [Rat; 3] {
        [
            ri(rng_i128(st, -100_000, 100_000), rng_i128(st, 1, 100)),
            ri(rng_i128(st, -100_000, 100_000), rng_i128(st, 1, 100)),
            ri(rng_i128(st, -100_000, 100_000), rng_i128(st, 1, 100)),
        ]
    }

    /// A random point rotated about one random axis by an inexact rational angle
    /// (heterogeneous provenance: each point its own rotation).
    fn rand_point(st: &mut u64) -> Pt3 {
        let base = rand_base(st);
        let axis = axis_of(rng_i128(st, 0, 2));
        let ang = Angle::from_deg(ri(rng_i128(st, 0, 360_000), rng_i128(st, 1, 9973))).unwrap();
        Pt3::at(base).rotate(axis, ang)
    }

    /// Absolute value of a `BigFloat` det as an f64 magnitude bound (for comparison).
    fn abs_err(f: f64, truth: &BigFloat, gt: usize) -> f64 {
        let e = BigFloat::from_f64(f, gt).sub(truth, gt, HP_RM);
        bf_mag(&e.abs())
    }

    /// Sanity: a known positive tetrahedron stays `Positive` (and its mirror
    /// `Negative`) at several shared rotations; a coplanar quad is `Zero`.
    #[test]
    fn orient3d_sign_is_stable_under_rotation() {
        // orient3d(a,b,c,d) with d the apex below the CCW triangle a,b,c.
        let a = [ri(0, 1), ri(0, 1), ri(0, 1)];
        let b = [ri(4, 1), ri(0, 1), ri(0, 1)];
        let c = [ri(0, 1), ri(4, 1), ri(0, 1)];
        let d = [ri(0, 1), ri(0, 1), ri(3, 1)]; // above the z=0 plane
        for (an, ad, ax) in [(0, 1, 0), (30, 1, 2), (37, 1, 1), (123, 1, 0), (1, 3, 2)] {
            let ang = Angle::from_deg(ri(an, ad)).unwrap();
            let axis = axis_of(ax);
            let mk = |p: [Rat; 3]| Pt3::at(p).rotate(axis, ang);
            let (pa, pb, pc, pd) = (mk(a), mk(b), mk(c), mk(d));
            let s = orient3d_judge(&pa, &pb, &pc, &pd);
            let sr = orient3d_judge(&pa, &pc, &pb, &pd); // swapped → opposite sign
            assert!(s == Orient::Positive || s == Orient::Negative, "{an}/{ad}");
            assert_ne!(s, sr, "swap must flip sign at {an}/{ad}");
        }
        // Four coplanar points (all z=0) → Zero at every rotation.
        let q = [
            [ri(0, 1), ri(0, 1), ri(0, 1)],
            [ri(2, 1), ri(0, 1), ri(0, 1)],
            [ri(2, 1), ri(2, 1), ri(0, 1)],
            [ri(1, 1), ri(3, 1), ri(0, 1)],
        ];
        for (an, ad, ax) in [(0, 1, 0), (30, 1, 1), (1, 7, 2)] {
            let ang = Angle::from_deg(ri(an, ad)).unwrap();
            let axis = axis_of(ax);
            let p: Vec<Pt3> = q.iter().map(|&b| Pt3::at(b).rotate(axis, ang)).collect();
            assert_eq!(
                orient3d_judge(&p[0], &p[1], &p[2], &p[3]),
                Orient::Zero,
                "coplanar at {an}/{ad}"
            );
        }
    }

    /// H-a — explicit orient3d filter soundness (the crux new math). Over many
    /// random heterogeneous-rotation 4-point configs, `det3_bound` must upper-bound
    /// the real error of the f64 determinant relative to the astro-float truth of
    /// the exact definitions — never once exceeded. Records worst tightness
    /// (actual/predicted); a value < 1 means sound, close to 1 means tight.
    #[test]
    fn h_a_orient3d_bound_soundness() {
        const GT: usize = 512;
        const N: usize = 10_000;
        let mut st = 0x3D00_1234_ABCD_EF01u64;
        let mut bad = 0usize;
        let mut worst_tightness = 0.0_f64;

        for _ in 0..N {
            let (pa, pb, pc, pd) = (
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
            );
            let det = det3_f64(pa.coord, pb.coord, pc.coord, pd.coord);
            let truth = det3_hp(&pa, &pb, &pc, &pd, GT);
            let err = abs_err(det, &truth, GT);
            let bound = det3_bound(
                [pa.coord, pb.coord, pc.coord, pd.coord],
                [pa.tol, pb.tol, pc.tol, pd.tol],
            );
            if err > bound {
                bad += 1;
            }
            if bound > 0.0 {
                worst_tightness = worst_tightness.max(err / bound);
            }
        }
        eprintln!(
            "[H-a N={N}] det3_bound violations: {bad}; worst tightness (err/bound): {worst_tightness:.3}"
        );
        assert_eq!(
            bad, 0,
            "det3_bound must bound the determinant error on every sample"
        );
    }
}
