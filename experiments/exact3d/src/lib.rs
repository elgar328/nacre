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

// ---- indirect orient3d: three rotated planes meet at an implicit point ----
//
// A vertex is `∩` of three planes, each defined by three rotated points, so its
// coordinates are never exact. The kernel's `indirect_orient3d` decides the sign of
// `orient3d(V, q, r, s)` without materializing `V`: `sign = sign(D)·sign(M)` where
// `D = det(normals)` and `M = (Dvec − D·s)·((q−s)×(r−s))` (Cramer, no division). We
// mirror that here, but over **intervals** (value ± tol) for the f64 filter — the
// design's "dynamic filter" (§TIP ⑨): interval arithmetic is a sound worst-case
// bound by construction, so no per-predicate bound formula is hand-derived. Ambiguous
// (an interval straddling 0) escalates to astro-float, exactly as `orient3d_judge`.

/// A value with a symmetric error radius (`mid ± rad`, `rad ≥ 0`). Arithmetic keeps
/// `rad` a sound upper bound (worst case), plus a per-op f64-rounding inflation.
#[derive(Clone, Copy, Debug)]
struct Iv {
    mid: f64,
    rad: f64,
}

impl Iv {
    fn new(mid: f64, rad: f64) -> Self {
        Iv { mid, rad }
    }
    fn sub(self, o: Iv) -> Iv {
        let mid = self.mid - o.mid;
        Iv::new(mid, self.rad + o.rad + 2.0 * f64::EPSILON * mid.abs())
    }
    fn add(self, o: Iv) -> Iv {
        let mid = self.mid + o.mid;
        Iv::new(mid, self.rad + o.rad + 2.0 * f64::EPSILON * mid.abs())
    }
    fn mul(self, o: Iv) -> Iv {
        let mid = self.mid * o.mid;
        let rad = self.mid.abs() * o.rad + o.mid.abs() * self.rad + self.rad * o.rad;
        Iv::new(mid, rad + 2.0 * f64::EPSILON * mid.abs())
    }
    /// `Some(true)` if definitely positive, `Some(false)` if definitely negative,
    /// `None` if the interval straddles 0 (escalate).
    fn sign(self) -> Option<bool> {
        if self.mid > self.rad {
            Some(true)
        } else if self.mid < -self.rad {
            Some(false)
        } else {
            None
        }
    }
}

/// 3×3 determinant of interval rows.
fn det3_iv(r: [[Iv; 3]; 3]) -> Iv {
    let m0 = r[1][1].mul(r[2][2]).sub(r[1][2].mul(r[2][1]));
    let m1 = r[1][0].mul(r[2][2]).sub(r[1][2].mul(r[2][0]));
    let m2 = r[1][0].mul(r[2][1]).sub(r[1][1].mul(r[2][0]));
    r[0][0].mul(m0).sub(r[0][1].mul(m1)).add(r[0][2].mul(m2))
}

/// A point's coord+tol as an interval per component.
fn pt_iv(p: &Pt3) -> [Iv; 3] {
    [
        Iv::new(p.coord[0], p.tol[0]),
        Iv::new(p.coord[1], p.tol[1]),
        Iv::new(p.coord[2], p.tol[2]),
    ]
}

/// Plane `[a,b,c,d]` (n·X + d = 0) through three points, as intervals: `n =
/// (p1−p0)×(p2−p0)`, `d = −n·p0`. Coefficient tol propagates from the point tols
/// through the subtraction/cross/dot — the "coefficient tol is a corollary of point
/// tol" of design.md §539.
fn plane_iv(p0: &Pt3, p1: &Pt3, p2: &Pt3) -> [Iv; 4] {
    let (a, b, c) = (pt_iv(p0), pt_iv(p1), pt_iv(p2));
    let e1 = [b[0].sub(a[0]), b[1].sub(a[1]), b[2].sub(a[2])];
    let e2 = [c[0].sub(a[0]), c[1].sub(a[1]), c[2].sub(a[2])];
    let n = [
        e1[1].mul(e2[2]).sub(e1[2].mul(e2[1])),
        e1[2].mul(e2[0]).sub(e1[0].mul(e2[2])),
        e1[0].mul(e2[1]).sub(e1[1].mul(e2[0])),
    ];
    let d = Iv::new(0.0, 0.0)
        .sub(n[0].mul(a[0]))
        .sub(n[1].mul(a[1]))
        .sub(n[2].mul(a[2]));
    [n[0], n[1], n[2], d]
}

/// The plane's four coefficients realized at `prec` bits from the definitions
/// (ground truth for the coefficient tol).
fn plane_hp(p0: &Pt3, p1: &Pt3, p2: &Pt3, prec: usize) -> [BigFloat; 4] {
    let (a, b, c) = (p0.hp_coord(prec), p1.hp_coord(prec), p2.hp_coord(prec));
    let sub = |x: &BigFloat, y: &BigFloat| x.sub(y, prec, HP_RM);
    let mul = |x: &BigFloat, y: &BigFloat| x.mul(y, prec, HP_RM);
    let e1 = [sub(&b[0], &a[0]), sub(&b[1], &a[1]), sub(&b[2], &a[2])];
    let e2 = [sub(&c[0], &a[0]), sub(&c[1], &a[1]), sub(&c[2], &a[2])];
    let n = [
        mul(&e1[1], &e2[2]).sub(&mul(&e1[2], &e2[1]), prec, HP_RM),
        mul(&e1[2], &e2[0]).sub(&mul(&e1[0], &e2[2]), prec, HP_RM),
        mul(&e1[0], &e2[1]).sub(&mul(&e1[1], &e2[0]), prec, HP_RM),
    ];
    let d = BigFloat::from_f64(0.0, prec)
        .sub(&mul(&n[0], &a[0]), prec, HP_RM)
        .sub(&mul(&n[1], &a[1]), prec, HP_RM)
        .sub(&mul(&n[2], &a[2]), prec, HP_RM);
    [n[0].clone(), n[1].clone(), n[2].clone(), d]
}

/// `sign(D)·sign(M)` combined into an orientation (both must be definite).
fn combine(dsign: Option<bool>, msign: Option<bool>) -> Option<Orient> {
    match (dsign, msign) {
        (Some(dp), Some(mp)) => Some(if dp == mp {
            Orient::Positive
        } else {
            Orient::Negative
        }),
        _ => None,
    }
}

/// The interval f64 filter for `orient3d(V, q, r, s)`, `V = ∩(planes)`. `None` if
/// either `D` or `M` straddles 0 (escalate). Coefficient-direct (no division).
fn indirect_filter(planes: [[Iv; 4]; 3], q: [Iv; 3], r: [Iv; 3], s: [Iv; 3]) -> Option<Orient> {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    let h = |k: usize| Iv::new(0.0, 0.0).sub(planes[k][3]); // n·X = h, h = −d
    let (n0, n1, n2) = (n(0), n(1), n(2));
    let d = det3_iv([n0, n1, n2]);
    // Cramer numerators Dvec: replace column j with h.
    let col_h = [h(0), h(1), h(2)];
    let dvec = [
        det3_iv([
            [col_h[0], n0[1], n0[2]],
            [col_h[1], n1[1], n1[2]],
            [col_h[2], n2[1], n2[2]],
        ]),
        det3_iv([
            [n0[0], col_h[0], n0[2]],
            [n1[0], col_h[1], n1[2]],
            [n2[0], col_h[2], n2[2]],
        ]),
        det3_iv([
            [n0[0], n0[1], col_h[0]],
            [n1[0], n1[1], col_h[1]],
            [n2[0], n2[1], col_h[2]],
        ]),
    ];
    // row1 = Dvec − D·s ; cross = (q−s)×(r−s) ; M = row1·cross.
    let row1 = [
        dvec[0].sub(d.mul(s[0])),
        dvec[1].sub(d.mul(s[1])),
        dvec[2].sub(d.mul(s[2])),
    ];
    let dq = [q[0].sub(s[0]), q[1].sub(s[1]), q[2].sub(s[2])];
    let dr = [r[0].sub(s[0]), r[1].sub(s[1]), r[2].sub(s[2])];
    let cross = [
        dq[1].mul(dr[2]).sub(dq[2].mul(dr[1])),
        dq[2].mul(dr[0]).sub(dq[0].mul(dr[2])),
        dq[0].mul(dr[1]).sub(dq[1].mul(dr[0])),
    ];
    let m = row1[0]
        .mul(cross[0])
        .add(row1[1].mul(cross[1]))
        .add(row1[2].mul(cross[2]));
    combine(d.sign(), m.sign())
}

/// The same `sign(D)·sign(M)` realized at `prec` bits (astro-float) — ground truth /
/// escalation. Returns `(D, M, mag_d, mag_m)`: the two determinants and the f64
/// magnitudes of their term sums, so the caller can floor a below-precision result to
/// declare-0 (`sign_with_floor`).
fn indirect_hp(
    planes: [[BigFloat; 4]; 3],
    q: [BigFloat; 3],
    r: [BigFloat; 3],
    s: [BigFloat; 3],
    prec: usize,
) -> (BigFloat, BigFloat, f64, f64) {
    let sub = |x: &BigFloat, y: &BigFloat| x.sub(y, prec, HP_RM);
    let mul = |x: &BigFloat, y: &BigFloat| x.mul(y, prec, HP_RM);
    let add = |x: &BigFloat, y: &BigFloat| x.add(y, prec, HP_RM);
    let zero = BigFloat::from_f64(0.0, prec);
    let det3 = |r: &[[BigFloat; 3]; 3]| {
        let m0 = sub(&mul(&r[1][1], &r[2][2]), &mul(&r[1][2], &r[2][1]));
        let m1 = sub(&mul(&r[1][0], &r[2][2]), &mul(&r[1][2], &r[2][0]));
        let m2 = sub(&mul(&r[1][0], &r[2][1]), &mul(&r[1][1], &r[2][0]));
        add(
            &sub(&mul(&r[0][0], &m0), &mul(&r[0][1], &m1)),
            &mul(&r[0][2], &m2),
        )
    };
    // Sum of the six triple-product magnitudes of a 3×3 (bounds the term scale).
    let det3_mag = |r: &[[BigFloat; 3]; 3]| -> f64 {
        let t = |a: &BigFloat, b: &BigFloat, c: &BigFloat| bf_mag(a) * bf_mag(b) * bf_mag(c);
        t(&r[0][0], &r[1][1], &r[2][2])
            + t(&r[0][0], &r[1][2], &r[2][1])
            + t(&r[0][1], &r[1][0], &r[2][2])
            + t(&r[0][1], &r[1][2], &r[2][0])
            + t(&r[0][2], &r[1][0], &r[2][1])
            + t(&r[0][2], &r[1][1], &r[2][0])
    };
    let nrm = |k: usize| {
        [
            planes[k][0].clone(),
            planes[k][1].clone(),
            planes[k][2].clone(),
        ]
    };
    let hh = |k: usize| sub(&zero, &planes[k][3]);
    let (n0, n1, n2) = (nrm(0), nrm(1), nrm(2));
    let nrows = [n0.clone(), n1.clone(), n2.clone()];
    let d = det3(&nrows);
    let mag_d = det3_mag(&nrows);
    let hc = [hh(0), hh(1), hh(2)];
    let dvec = [
        det3(&[
            [hc[0].clone(), n0[1].clone(), n0[2].clone()],
            [hc[1].clone(), n1[1].clone(), n1[2].clone()],
            [hc[2].clone(), n2[1].clone(), n2[2].clone()],
        ]),
        det3(&[
            [n0[0].clone(), hc[0].clone(), n0[2].clone()],
            [n1[0].clone(), hc[1].clone(), n1[2].clone()],
            [n2[0].clone(), hc[2].clone(), n2[2].clone()],
        ]),
        det3(&[
            [n0[0].clone(), n0[1].clone(), hc[0].clone()],
            [n1[0].clone(), n1[1].clone(), hc[1].clone()],
            [n2[0].clone(), n2[1].clone(), hc[2].clone()],
        ]),
    ];
    let row1 = [
        sub(&dvec[0], &mul(&d, &s[0])),
        sub(&dvec[1], &mul(&d, &s[1])),
        sub(&dvec[2], &mul(&d, &s[2])),
    ];
    let dq = [sub(&q[0], &s[0]), sub(&q[1], &s[1]), sub(&q[2], &s[2])];
    let dr = [sub(&r[0], &s[0]), sub(&r[1], &s[1]), sub(&r[2], &s[2])];
    let cross = [
        sub(&mul(&dq[1], &dr[2]), &mul(&dq[2], &dr[1])),
        sub(&mul(&dq[2], &dr[0]), &mul(&dq[0], &dr[2])),
        sub(&mul(&dq[0], &dr[1]), &mul(&dq[1], &dr[0])),
    ];
    let m = add(
        &add(&mul(&row1[0], &cross[0]), &mul(&row1[1], &cross[1])),
        &mul(&row1[2], &cross[2]),
    );
    let mag_m = bf_mag(&row1[0]) * bf_mag(&cross[0])
        + bf_mag(&row1[1]) * bf_mag(&cross[1])
        + bf_mag(&row1[2]) * bf_mag(&cross[2]);
    (d, m, mag_d, mag_m)
}

/// Multiplier over `mag·2⁻ᵖʳᵉᶜ` for the declare-0 floor — generously above the
/// accumulated rounding of the ~100-op indirect computation (soundness first; a
/// tighter constant would only reduce the rare declare-0 rate, never soundness).
const FLOOR_K: f64 = 1.0e6;

/// A `BigFloat`'s sign, floored to declare-0 (`None`) when its magnitude is at or
/// below the rounding floor for a `prec`-bit computation whose terms scale to `mag`.
fn sign_with_floor(val: &BigFloat, mag: f64, prec: usize) -> Option<i8> {
    let floor = FLOOR_K * mag * 2f64.powi(-(prec as i32));
    if val.is_zero() || bf_mag(val) <= floor {
        None
    } else {
        Some(bf_sign(val))
    }
}

/// Sign of a `BigFloat` as `+1/−1/0` (astro-float#44 workaround via `is_positive`).
fn bf_sign(x: &BigFloat) -> i8 {
    if x.is_zero() {
        0
    } else if x.is_positive() {
        1
    } else {
        -1
    }
}

/// TIP indirect `orient3d(V, q, r, s)`, `V = ∩(3 planes)` — each plane through three
/// rotated points, the triangle three rotated points. Interval filter → astro-float
/// escalation. `Zero` when the high-precision `D` or `M` is exactly zero (a genuine
/// coincidence; the near-degenerate declare-0 floor is a later refinement).
#[allow(clippy::too_many_arguments)]
pub fn indirect_orient3d_judge(
    plane_a: (&Pt3, &Pt3, &Pt3),
    plane_b: (&Pt3, &Pt3, &Pt3),
    plane_c: (&Pt3, &Pt3, &Pt3),
    q: &Pt3,
    r: &Pt3,
    s: &Pt3,
) -> Orient {
    let planes = [
        plane_iv(plane_a.0, plane_a.1, plane_a.2),
        plane_iv(plane_b.0, plane_b.1, plane_b.2),
        plane_iv(plane_c.0, plane_c.1, plane_c.2),
    ];
    if let Some(o) = indirect_filter(planes, pt_iv(q), pt_iv(r), pt_iv(s)) {
        return o;
    }
    // Escalate. Trust each sign only if its magnitude clears the rounding floor for a
    // JUDGE_PREC computation (`sign_with_floor`); a below-floor `D` or `M` is declare-0
    // (§6 "ask the user"), never a sub-floor noisy sign.
    let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, JUDGE_PREC);
    let (d, m, mag_d, mag_m) = indirect_hp(
        [ph(plane_a), ph(plane_b), ph(plane_c)],
        q.hp_coord(JUDGE_PREC),
        r.hp_coord(JUDGE_PREC),
        s.hp_coord(JUDGE_PREC),
        JUDGE_PREC,
    );
    match (
        sign_with_floor(&d, mag_d, JUDGE_PREC),
        sign_with_floor(&m, mag_m, JUDGE_PREC),
    ) {
        (Some(ds), Some(ms)) => {
            if ds == ms {
                Orient::Positive
            } else {
                Orient::Negative
            }
        }
        _ => Orient::Zero,
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

    /// H-b — rotated plane coefficient tol. A plane's four coefficients derive from
    /// three rotated points by subtraction/cross/dot; the interval `rad` on each must
    /// upper-bound the real error of the f64 coefficient vs the astro-float truth.
    #[test]
    fn h_b_plane_coefficient_tol_soundness() {
        const GT: usize = 512;
        const N: usize = 10_000;
        let mut st = 0xB0B0_5555_1111_2222u64;
        let mut bad = 0usize;
        let mut worst = 0.0_f64;
        for _ in 0..N {
            let (p0, p1, p2) = (
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
            );
            let iv = plane_iv(&p0, &p1, &p2);
            let hp = plane_hp(&p0, &p1, &p2, GT);
            for k in 0..4 {
                let err = abs_err(iv[k].mid, &hp[k], GT);
                if err > iv[k].rad {
                    bad += 1;
                }
                if iv[k].rad > 0.0 {
                    worst = worst.max(err / iv[k].rad);
                }
            }
        }
        eprintln!("[H-b N={N}] coefficient tol violations: {bad}; worst tightness: {worst:.3}");
        assert_eq!(
            bad, 0,
            "plane coefficient tol must bound the error on every sample"
        );
    }

    /// Sign of the indirect predicate at `prec` bits with a stability flag: `None` if
    /// even the ground truth cannot resolve it (recomputing at `prec+128` disagrees,
    /// or a component is exactly 0 — a genuine degeneracy), else `Some(orient)`.
    #[allow(clippy::too_many_arguments)]
    fn indirect_truth(
        pa: (&Pt3, &Pt3, &Pt3),
        pb: (&Pt3, &Pt3, &Pt3),
        pc: (&Pt3, &Pt3, &Pt3),
        q: &Pt3,
        r: &Pt3,
        s: &Pt3,
        prec: usize,
    ) -> Option<Orient> {
        let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, prec);
        let (d, m, mag_d, mag_m) = indirect_hp(
            [ph(pa), ph(pb), ph(pc)],
            q.hp_coord(prec),
            r.hp_coord(prec),
            s.hp_coord(prec),
            prec,
        );
        match (
            sign_with_floor(&d, mag_d, prec),
            sign_with_floor(&m, mag_m, prec),
        ) {
            (Some(ds), Some(ms)) => Some(if ds == ms {
                Orient::Positive
            } else {
                Orient::Negative
            }),
            _ => None, // below the GT floor → genuine degeneracy
        }
    }

    /// H-c — indirect orient3d soundness with heterogeneous provenance and a
    /// near-degenerate corpus (the escalation path is otherwise never exercised). The
    /// interval filter → astro-float judge must never disagree with a GT-stable
    /// 512-bit truth. Also reports how often the filter escalated (must be > 0, else
    /// the test is vacuous).
    #[test]
    fn h_c_indirect_orient3d_soundness() {
        const GT: usize = 512;
        // wrong = judge claims a definite sign opposite to GT (SOUNDNESS violation).
        // declined = judge declares 0 where GT resolves a sign (sound: deferred to the
        // user, never wrong). skipped = even GT cannot resolve (genuine degeneracy).
        let (mut wrong, mut declined, mut escalated, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0);

        // Evaluate one config: accumulate the counters.
        let mut check = |a: (&Pt3, &Pt3, &Pt3),
                         b: (&Pt3, &Pt3, &Pt3),
                         c: (&Pt3, &Pt3, &Pt3),
                         q: &Pt3,
                         r: &Pt3,
                         s: &Pt3| {
            let judged = indirect_orient3d_judge(a, b, c, q, r, s);
            let planes = [
                plane_iv(a.0, a.1, a.2),
                plane_iv(b.0, b.1, b.2),
                plane_iv(c.0, c.1, c.2),
            ];
            if indirect_filter(planes, pt_iv(q), pt_iv(r), pt_iv(s)).is_none() {
                escalated += 1;
            }
            match indirect_truth(a, b, c, q, r, s, GT) {
                None => skipped += 1,
                Some(truth) => {
                    tested += 1;
                    if judged == Orient::Zero {
                        declined += 1;
                    } else if judged != truth {
                        wrong += 1;
                    }
                }
            }
        };

        // Corpus A — heterogeneous provenance, generic position (each of the twelve
        // points its own rotation). The transversal two-solid case: mixed tol in one
        // predicate. Generic → mostly filter-resolved; validates mixed-provenance
        // soundness (no disagreement).
        let mut st = 0xC0C0_9999_ABAB_CDCDu64;
        for _ in 0..2000 {
            let p: Vec<Pt3> = (0..12).map(|_| rand_point(&mut st)).collect();
            check(
                (&p[0], &p[1], &p[2]),
                (&p[3], &p[4], &p[5]),
                (&p[6], &p[7], &p[8]),
                &p[9],
                &p[10],
                &p[11],
            );
        }

        // Corpus B — shared rotation, near-coplanar (escalation-forcing). Rotation
        // preserves coplanarity, so a triangle built coplanar with the planes'
        // intersection V (with s pushed a tiny ε off the plane) stays near-degenerate
        // after rotation: the true orient is tiny, the interval filter is ambiguous,
        // and the astro-float escalation must still match the GT-stable sign.
        let add3 = |a: [Rat; 3], b: [Rat; 3]| {
            [
                a[0].checked_add(b[0]).unwrap(),
                a[1].checked_add(b[1]).unwrap(),
                a[2].checked_add(b[2]).unwrap(),
            ]
        };
        let smul = |k: Rat, a: [Rat; 3]| {
            [
                a[0].checked_mul(k).unwrap(),
                a[1].checked_mul(k).unwrap(),
                a[2].checked_mul(k).unwrap(),
            ]
        };
        let cross = |a: [Rat; 3], b: [Rat; 3]| {
            [
                a[1].checked_mul(b[2])
                    .unwrap()
                    .checked_sub(a[2].checked_mul(b[1]).unwrap())
                    .unwrap(),
                a[2].checked_mul(b[0])
                    .unwrap()
                    .checked_sub(a[0].checked_mul(b[2]).unwrap())
                    .unwrap(),
                a[0].checked_mul(b[1])
                    .unwrap()
                    .checked_sub(a[1].checked_mul(b[0]).unwrap())
                    .unwrap(),
            ]
        };
        for _ in 0..2000 {
            let off = |st: &mut u64| {
                [
                    ri(rng_i128(st, -20, 20), rng_i128(st, 1, 5)),
                    ri(rng_i128(st, -20, 20), rng_i128(st, 1, 5)),
                    ri(rng_i128(st, -20, 20), rng_i128(st, 1, 5)),
                ]
            };
            let v = [
                ri(rng_i128(&mut st, -200, 200), 1),
                ri(rng_i128(&mut st, -200, 200), 1),
                ri(rng_i128(&mut st, -200, 200), 1),
            ];
            let (u, w) = (off(&mut st), off(&mut st));
            let n = cross(u, w); // normal to the V-plane
            // ε off-plane push (tiny): 0 (exactly degenerate) or 1/big (near).
            // ε spans exactly-degenerate (0) through a wide off-plane range: large
            // denominators (down to ~1e-15) land in the window where the interval
            // filter is ambiguous but the true sign is definite, so the astro-float
            // escalation's sign-resolution is exercised (not just declare-0).
            let eps = match rng_i128(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng_i128(&mut st, 5_000, 200_000)),
                _ => ri(1, rng_i128(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            // Plane k through V + two random in-space offsets.
            let plane_pts = |st: &mut u64| [v, add3(v, off(st)), add3(v, off(st))];
            let a = plane_pts(&mut st);
            let b = plane_pts(&mut st);
            let c = plane_pts(&mut st);
            // Triangle q=v+u, r=v+w, s=v+u+w — coplanar with V=v (true orient 0),
            // non-degenerate (q−s=−w, r−s=−u independent), s pushed ε off-plane along n.
            let qb = add3(v, u);
            let rb = add3(v, w);
            let sb = add3(add3(add3(v, u), w), smul(eps, n));
            // One shared rotation for all twelve points.
            let axis = axis_of(rng_i128(&mut st, 0, 2));
            let ang = Angle::from_deg(ri(
                rng_i128(&mut st, 0, 360_000),
                rng_i128(&mut st, 1, 9973),
            ))
            .unwrap();
            let rp = |p: [Rat; 3]| Pt3::at(p).rotate(axis, ang);
            let (a0, a1, a2) = (rp(a[0]), rp(a[1]), rp(a[2]));
            let (b0, b1, b2) = (rp(b[0]), rp(b[1]), rp(b[2]));
            let (c0, c1, c2) = (rp(c[0]), rp(c[1]), rp(c[2]));
            let (q, r, s) = (rp(qb), rp(rb), rp(sb));
            check(
                (&a0, &a1, &a2),
                (&b0, &b1, &b2),
                (&c0, &c1, &c2),
                &q,
                &r,
                &s,
            );
        }

        eprintln!(
            "[H-c] wrong-sign: {wrong}/{tested} tested; declined (sound): {declined}; escalated: {escalated}; degenerate skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "indirect judge must never claim a definite sign opposite to GT (soundness)"
        );
        assert!(
            escalated > 0,
            "corpus must exercise the escalation path (else vacuous)"
        );
    }
}
