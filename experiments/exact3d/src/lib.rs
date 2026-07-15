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

/// One rotation in a point's definition: turn about `axis` (through the rational
/// pivot `point`) by the rational `angle`. `point = [0,0,0]` is the origin-pivot case
/// the 2D-style experiment used; the kernel's `Rotation` carries an arbitrary rational
/// pivot, so the port must track it (pivot-relative realization + its rounding).
#[derive(Clone, Copy, Debug)]
pub struct RotNode {
    pub axis: Axis,
    pub angle: Angle,
    pub point: [Rat; 3],
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
    /// A point at `base`. The coordinate is `base` realized in f64, which carries the
    /// rational→f64 rounding error (a division), so the initial tol is that error —
    /// exactly 0 for an f64-representable base (e.g. a small integer), positive
    /// otherwise. An axis a later chain never rotates keeps exactly this, which is why
    /// a per-axis tol check needs it (a determinant/coefficient check hides it in the
    /// rotated-plane mix).
    pub fn at(base: [Rat; 3]) -> Self {
        let coord = [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()];
        // Actual rounding = |coord − base| at high precision; 2× the power-of-two
        // magnitude is a sound upper bound (0 when the base is exact).
        let round_tol = |b: Rat, c: f64| {
            let e = BigFloat::from_f64(c, 120).sub(&rat_to_big(b, 120), 120, HP_RM);
            2.0 * bf_mag(&e.abs())
        };
        Pt3 {
            tol: [
                round_tol(base[0], coord[0]),
                round_tol(base[1], coord[1]),
                round_tol(base[2], coord[2]),
            ],
            coord,
            base,
            chain: Vec::new(),
        }
    }

    /// Extend the definition by a rotation about `axis` by `angle`, updating the f64
    /// cache and its tol. Same-axis 90°-family angles rotate exactly (tol 0, Niven);
    /// otherwise the two in-plane coords gain the coordinate-mixing realization error
    /// `(|u|+|v|)·da`, and every existing tol is transported by the component-wise
    /// absolute rotation `|R|` (§TIP: `new tol = |R|·old + mix`).
    pub fn rotate(self, axis: Axis, angle: Angle) -> Self {
        let z = Rat::from_int(0);
        self.rotate_about(axis, angle, [z; 3])
    }

    /// [`rotate`] about an axis through an arbitrary rational pivot `point`. The two
    /// in-plane coords are taken relative to the pivot (`u = coord − p`), rotated, and
    /// shifted back. Origin pivot (`point = 0`) reproduces [`rotate`] bit-for-bit
    /// (subtracting `0.0` is exact). A non-origin pivot adds its own f64 rounding
    /// (the `coord − p` and `p + …` arithmetic), which the tol must cover even for an
    /// exact (90°-family) angle — the one soundness question the origin-only corpus
    /// never asked.
    pub fn rotate_about(mut self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let (i, j) = plane_of(axis);
        let (px, py) = (point[i].to_f64(), point[j].to_f64());
        let (ci, cj) = (self.coord[i], self.coord[j]); // pre-rotation magnitudes for tol
        let (u, v) = (ci - px, cj - py);
        // cos/sin — exact (rational) for the 90°-family, else f64 (with realization tol).
        let (c, s, exact) = match angle.try_exact_cos_sin() {
            Some((cr, sr)) => (cr.to_f64(), sr.to_f64(), true),
            None => (angle.cos(), angle.sin(), false),
        };
        self.coord[i] = px + u * c - v * s;
        self.coord[j] = py + u * s + v * c;
        // Rotation-realization error (coordinate-mixing), 0 for an exact angle.
        let rot = if exact {
            0.0
        } else {
            (u.abs() + v.abs()) * DA_F64
        };
        // Pivot arithmetic — the `coord − p` subtraction, the `p + …` re-add, and the
        // pivot's own Rat→f64 rounding. Exactly 0 for an origin pivot (subtracting/
        // adding 0.0 is exact); otherwise ~ulp of every magnitude involved (present
        // even when the angle is exact, and even for a point sitting on the pivot).
        let piv = if px != 0.0 || py != 0.0 {
            (ci.abs() + cj.abs() + px.abs() + py.abs()) * DA_F64
        } else {
            0.0
        };
        let mix = rot + piv;
        let (ti, tj) = (self.tol[i], self.tol[j]);
        self.tol[i] = c.abs() * ti + s.abs() * tj + mix;
        self.tol[j] = s.abs() * ti + c.abs() * tj + mix;
        self.chain.push(RotNode { axis, angle, point });
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
            let (px, py) = (
                rat_to_big(node.point[i], prec),
                rat_to_big(node.point[j], prec),
            );
            // pivot-relative: u = p − pivot, rotate, shift back (exact at `prec` bits).
            let u = p[i].sub(&px, prec, HP_RM);
            let v = p[j].sub(&py, prec, HP_RM);
            p[i] = px.add(
                &u.mul(&c, prec, HP_RM)
                    .sub(&v.mul(&s, prec, HP_RM), prec, HP_RM),
                prec,
                HP_RM,
            );
            p[j] = py.add(
                &u.mul(&s, prec, HP_RM)
                    .add(&v.mul(&c, prec, HP_RM), prec, HP_RM),
                prec,
                HP_RM,
            );
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

    /// H-d — three-axis tol propagation across a rotation chain. A point carried
    /// through several rotations about mixed axes accumulates its tol by
    /// `new = |R|·old + mix` at each node (§TIP ⑦ traversal); the accumulated xyz tol
    /// must still bound the real coordinate error vs the astro-float truth — on all
    /// three axes, which a single-axis rotation never exercises.
    #[test]
    fn h_d_chain_tol_propagation() {
        const GT: usize = 512;
        const N: usize = 5_000;
        let mut st = 0xD00D_3333_7777_4444u64;
        let mut bad = 0usize;
        let mut worst = 0.0_f64;
        let mut axes_seen = [false; 3];
        for _ in 0..N {
            let mut p = Pt3::at(rand_base(&mut st));
            let k = rng_i128(&mut st, 2, 5);
            for _ in 0..k {
                let ax = axis_of(rng_i128(&mut st, 0, 2));
                match ax {
                    Axis::X => axes_seen[0] = true,
                    Axis::Y => axes_seen[1] = true,
                    Axis::Z => axes_seen[2] = true,
                }
                let ang = Angle::from_deg(ri(
                    rng_i128(&mut st, 0, 360_000),
                    rng_i128(&mut st, 1, 9973),
                ))
                .unwrap();
                p = p.rotate(ax, ang);
            }
            let hp = p.hp_coord(GT);
            for (axis, hp_a) in hp.iter().enumerate() {
                let err = abs_err(p.coord[axis], hp_a, GT);
                if err > p.tol[axis] {
                    bad += 1;
                }
                if p.tol[axis] > 0.0 {
                    worst = worst.max(err / p.tol[axis]);
                }
            }
        }
        eprintln!(
            "[H-d N={N}] chain-tol violations: {bad}; worst tightness: {worst:.3}; axes exercised: {axes_seen:?}"
        );
        assert_eq!(
            bad, 0,
            "chain-accumulated tol must bound the error on every axis"
        );
        assert!(
            axes_seen == [true; 3],
            "the corpus must rotate about all three axes"
        );
    }

    /// H-f — pivot soundness. The origin-only corpus (H-a..H-e) never exercised the
    /// kernel's arbitrary rational pivots. A rotation about an axis through a non-origin
    /// point adds f64 rounding (`coord − pivot`, the `pivot + …` re-add, and the pivot's
    /// own Rat→f64) that the tol must cover — even for an exact (90°-family) angle, and
    /// even for a point sitting on the pivot. Random chains with mixed pivots (origin,
    /// on-point, ordinary, far) and mixed angles (exact / inexact) must never violate.
    #[test]
    fn h_f_pivot_soundness() {
        const GT: usize = 512;
        const N: usize = 5_000;
        let mut st = 0x9111_2222_3333_4444u64;
        let mut bad = 0usize;
        let mut worst = 0.0_f64;
        let (mut exact_seen, mut far_seen) = (false, false);
        for _ in 0..N {
            let base = rand_base(&mut st);
            let mut p = Pt3::at(base);
            let k = rng_i128(&mut st, 2, 5);
            for _ in 0..k {
                let ax = axis_of(rng_i128(&mut st, 0, 2));
                let pivot = match rng_i128(&mut st, 0, 3) {
                    0 => [Rat::from_int(0); 3], // origin (must reproduce rotate())
                    1 => base,                  // point sits on the pivot (fixed point)
                    2 => rand_base(&mut st),    // ordinary
                    _ => {
                        far_seen = true;
                        let f = |st: &mut u64| ri(rng_i128(st, -10_000, 10_000) * 1_000, 1);
                        [f(&mut st), f(&mut st), f(&mut st)] // far pivot
                    }
                };
                let ang = if rng_i128(&mut st, 0, 3) == 0 {
                    exact_seen = true;
                    Angle::from_deg(ri(90 * rng_i128(&mut st, 0, 3), 1)).unwrap() // 0/90/180/270
                } else {
                    Angle::from_deg(ri(
                        rng_i128(&mut st, 0, 360_000),
                        rng_i128(&mut st, 1, 9973),
                    ))
                    .unwrap()
                };
                p = p.rotate_about(ax, ang, pivot);
            }
            let hp = p.hp_coord(GT);
            for (axis, hp_a) in hp.iter().enumerate() {
                let err = abs_err(p.coord[axis], hp_a, GT);
                // Ground-truth realization noise floor: the 512-bit astro-float truth
                // carries ~2^-GT·(chain magnitude) noise (a far pivot amplifies it to
                // ~1e-147). A *real* f64 error is never below ~1e-16·|coord|, so any
                // `err` under 1e-100 is pure GT noise — floor it. (Only matters for the
                // tol-0 exact cases the origin-only corpus never mixed with far pivots.)
                if err > p.tol[axis] && err > 1e-100 {
                    bad += 1;
                }
                if p.tol[axis] > 0.0 {
                    worst = worst.max(err / p.tol[axis]);
                }
            }
        }
        eprintln!(
            "[H-f N={N}] pivot-tol violations: {bad}; worst tightness: {worst:.3}; exact_seen={exact_seen} far_seen={far_seen}"
        );
        assert_eq!(bad, 0, "pivot-aware tol must bound the error on every axis");
        assert!(
            exact_seen && far_seen,
            "corpus must include exact and far-pivot rotations"
        );
    }

    /// H-e — axis-change and bundling. (a) A chain of 90°-family rotations stays
    /// exact (tol 0, Niven), whatever the axes. (b) Same-axis rotation bundles: K
    /// incremental steps of θ realize the same geometry as one step of Kθ, but the
    /// incremental tol *amplifies* (each node transports the tol by `|c|+|s| ≥ 1`),
    /// while the bundled tol stays flat — both sound, so bundling is mandatory (the
    /// design's H4-amplification, here across a shared axis). A different-axis node
    /// cannot be bundled (the composition is not one axis-angle), so its tol simply
    /// accumulates — already covered sound by H-d.
    #[test]
    fn h_e_axis_change_and_bundling() {
        const GT: usize = 512;
        // (a) 90°-family chain across mixed axes → tol exactly 0 (exact realization).
        let p = Pt3::at([ri(3, 1), ri(5, 1), ri(7, 1)])
            .rotate(Axis::Z, Angle::from_deg(ri(90, 1)).unwrap())
            .rotate(Axis::X, Angle::from_deg(ri(180, 1)).unwrap())
            .rotate(Axis::Y, Angle::from_deg(ri(270, 1)).unwrap());
        assert_eq!(p.tol, [0.0; 3], "90°-family chain must stay tol 0");

        // (b) same-axis bundling: K incremental θ vs one Kθ about Z.
        const K: i128 = 30;
        let theta = ri(1, 7); // 1/7 degree, inexact
        let base = [ri(11, 1), ri(-7, 1), ri(4, 1)];
        let mut incr = Pt3::at(base);
        for _ in 0..K {
            incr = incr.rotate(Axis::Z, Angle::from_deg(theta).unwrap());
        }
        let bundled = Pt3::at(base).rotate(
            Axis::Z,
            Angle::from_deg(theta.checked_mul(Rat::from_int(K)).unwrap()).unwrap(),
        );

        // Both realize the same geometry, so both must bound the (shared) true error.
        let hp = bundled.hp_coord(GT);
        for (axis, hp_a) in hp.iter().enumerate() {
            let err = abs_err(bundled.coord[axis], hp_a, GT);
            assert!(err <= bundled.tol[axis], "bundled unsound on axis {axis}");
            assert!(err <= incr.tol[axis], "incremental unsound on axis {axis}");
        }
        let incr_tol = incr.tol[0].max(incr.tol[1]);
        let bund_tol = bundled.tol[0].max(bundled.tol[1]);
        eprintln!(
            "[H-e] bundled tol {bund_tol:.3e} vs incremental tol {incr_tol:.3e} (amplification {:.1}x over K={K})",
            incr_tol / bund_tol
        );
        assert!(
            bund_tol < incr_tol,
            "bundling must be tighter than incremental (amplification)"
        );
    }

    /// Aux — the go/no-go performance signal: how often does the f64 filter escalate,
    /// and what does escalation cost? On generic (non-degenerate) configs the filter
    /// resolves almost everything (cheap ~ns); escalation fires only near degeneracy
    /// (~µs, rare). Also confirms the explicit judge never disagrees with the 512-bit
    /// truth end to end (soundness). If escalation were frequent on generic geometry,
    /// that would be the performance NO-GO — it is not.
    #[test]
    fn aux_escalation_frequency_and_speed() {
        use std::time::Instant;
        const GT: usize = 512;

        // (1) generic random 4-point configs: escalation must be rare, judge sound.
        let mut st = 0xAABB_1122_3344_5566u64;
        let (mut esc_generic, mut wrong, mut tested) = (0usize, 0usize, 0usize);
        const NG: usize = 4000;
        for _ in 0..NG {
            let p: Vec<Pt3> = (0..4).map(|_| rand_point(&mut st)).collect();
            let det = det3_f64(p[0].coord, p[1].coord, p[2].coord, p[3].coord);
            let bound = det3_bound(
                [p[0].coord, p[1].coord, p[2].coord, p[3].coord],
                [p[0].tol, p[1].tol, p[2].tol, p[3].tol],
            );
            if det.abs() <= bound {
                esc_generic += 1;
            }
            // soundness end-to-end vs GT (generic configs are never truly degenerate).
            let truth = det3_hp(&p[0], &p[1], &p[2], &p[3], GT);
            if !truth.is_zero() {
                let want = if truth.is_positive() {
                    Orient::Positive
                } else {
                    Orient::Negative
                };
                let judged = orient3d_judge(&p[0], &p[1], &p[2], &p[3]);
                if judged != Orient::Zero {
                    tested += 1;
                    if judged != want {
                        wrong += 1;
                    }
                }
            }
        }

        // (2) near-coplanar configs (shared rotation preserves coplanarity): escalation
        // is frequent — the filter correctly defers the hard cases.
        let mut esc_near = 0usize;
        const NN: usize = 2000;
        for _ in 0..NN {
            let b = |st: &mut u64| {
                [
                    ri(rng_i128(st, -50, 50), 1),
                    ri(rng_i128(st, -50, 50), 1),
                    ri(rng_i128(st, -50, 50), 1),
                ]
            };
            let p0 = b(&mut st);
            let p1 = b(&mut st);
            let p2 = b(&mut st);
            // p3 = p0 + a(p1−p0) + c(p2−p0): exactly coplanar with p0,p1,p2.
            let (aa, cc) = (rng_i128(&mut st, -3, 3), rng_i128(&mut st, -3, 3));
            let mix = |i: usize| {
                p0[i]
                    .checked_add(
                        (p1[i].checked_sub(p0[i]).unwrap())
                            .checked_mul(ri(aa, 1))
                            .unwrap(),
                    )
                    .unwrap()
                    .checked_add(
                        (p2[i].checked_sub(p0[i]).unwrap())
                            .checked_mul(ri(cc, 1))
                            .unwrap(),
                    )
                    .unwrap()
            };
            let p3 = [mix(0), mix(1), mix(2)];
            let axis = axis_of(rng_i128(&mut st, 0, 2));
            let ang = Angle::from_deg(ri(
                rng_i128(&mut st, 0, 360_000),
                rng_i128(&mut st, 1, 9973),
            ))
            .unwrap();
            let rp = |p: [Rat; 3]| Pt3::at(p).rotate(axis, ang);
            let (q0, q1, q2, q3) = (rp(p0), rp(p1), rp(p2), rp(p3));
            let det = det3_f64(q0.coord, q1.coord, q2.coord, q3.coord);
            let bound = det3_bound(
                [q0.coord, q1.coord, q2.coord, q3.coord],
                [q0.tol, q1.tol, q2.tol, q3.tol],
            );
            if det.abs() <= bound {
                esc_near += 1;
            }
        }

        // (3) timing: the f64 filter vs one astro-float escalation.
        let p: Vec<Pt3> = (0..4).map(|_| rand_point(&mut st)).collect();
        let t0 = Instant::now();
        let mut acc = 0.0;
        for _ in 0..100_000 {
            acc += det3_f64(p[0].coord, p[1].coord, p[2].coord, p[3].coord)
                - det3_bound(
                    [p[0].coord, p[1].coord, p[2].coord, p[3].coord],
                    [p[0].tol, p[1].tol, p[2].tol, p[3].tol],
                );
        }
        let filter_ns = t0.elapsed().as_nanos() as f64 / 100_000.0;
        let t1 = Instant::now();
        for _ in 0..2_000 {
            let _ = det3_hp(&p[0], &p[1], &p[2], &p[3], JUDGE_PREC);
        }
        let esc_us = t1.elapsed().as_nanos() as f64 / 2_000.0 / 1000.0;
        std::hint::black_box(acc);

        eprintln!(
            "[aux] escalation: generic {esc_generic}/{NG} ({:.2}%), near-coplanar {esc_near}/{NN} ({:.1}%); \
             wrong-sign {wrong}/{tested}; filter ~{filter_ns:.0}ns, escalation ~{esc_us:.0}µs",
            100.0 * esc_generic as f64 / NG as f64,
            100.0 * esc_near as f64 / NN as f64,
        );
        assert_eq!(wrong, 0, "explicit judge must match GT (soundness)");
        assert!(
            esc_generic * 20 < NG,
            "escalation on generic geometry must be rare (< 5%)"
        );
        assert!(
            esc_near > 0,
            "near-coplanar corpus must exercise escalation"
        );
    }
}
