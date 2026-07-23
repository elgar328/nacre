//! The 3D toleranced point + its `orient3d` judgment (design.md §9 CIP, stage 2). The
//! 3D analogue of [`crate::frame2`] (`Pt2` / `orient2d_judge`).
//!
//! A rotated point cannot be held exactly (cos/sin are irrational), but its f64
//! realization carries a **direction-wise xyz tol** (§CIP ⑤) that soundly bounds
//! the error, accumulated as the definition is turned through a chain of
//! axis-aligned rotations (§CIP ①: `new tol = |R|·old + mix`). [`orient3d_judge`]
//! consumes that tol: an f64 determinant filter with a sound error bound
//! ([`det3_bound`], §CIP ②) decides the easy cases, the ambiguous ones **escalate**
//! to astro-float from the point definitions, and a determinant below the precision
//! floor is a **declare-0** ([`Orient::Zero`]). This is the judgment **layer only**
//! ("층만") — not yet wired into boolean (that is stage 3), and it is the *tol > 0*
//! path: a tol-0 (`Constructed`) config is faster/exact via `nacre-predicates`
//! (Shewchuk), routed by a higher layer, not here.
//!
//! Validated in `experiments/exact3d` (FINDINGS = GO): H-a (`det3_bound` soundness),
//! H-d/H-f (chain + arbitrary-pivot tol) — the bound never under-estimates the true
//! error (astro-float ground truth) over random heterogeneous-rotation configs.

use crate::frame2::{DA_F64, JUDGE_PREC, bf_mag, rat_to_big};
use crate::{Angle, Axis, HP_RM, Orient, Rat};
use astro_float::BigFloat;
#[cfg(feature = "parallel")]
use std::sync::{Arc as HpRc, OnceLock as HpOnce};
#[cfg(not(feature = "parallel"))]
use std::{cell::OnceCell as HpOnce, rc::Rc as HpRc};

/// Shared, lazily-initialized cell for the memoized high-precision realization.
/// Under `parallel` it is `Arc<OnceLock>` — `Send + Sync`, so a `Pt3` crosses rayon
/// worker threads and the cache is shared (a hot definition-point's 200-bit realization
/// is computed once and reused across workers, not per thread). Otherwise it is
/// `Rc<OnceCell>` — single-threaded, no atomic overhead. `get_or_init` has the identical
/// signature on both, so the consumer ([`Pt3::hp_coord`]) is unchanged by the choice.
type HpCell = HpRc<HpOnce<[BigFloat; 3]>>;

/// One rotation in a point's definition: turn about `axis` (the line through the
/// rational pivot `point`) by the rational `angle`. `point = [0,0,0]` is the
/// origin-pivot case; the kernel's `Rotation` carries an arbitrary rational pivot.
#[derive(Clone, Copy, Debug)]
pub struct RotNode {
    pub axis: Axis,
    pub angle: Angle,
    pub point: [Rat; 3],
}

/// A rational base point carried through a chain of axis rotations (§CIP ⑦ rotation
/// history). `base` + `chain` are the exact **definition** (never lost); `coord` is
/// the f64 realization (a cache), and `tol` bounds its error as a **direction-wise
/// xyz vector** (§CIP ⑤). [`hp_coord`](Self::hp_coord) realizes the chain at
/// arbitrary precision from the definition, so two points with the same definition
/// realize identically (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct Pt3 {
    pub base: [Rat; 3],
    pub chain: Vec<RotNode>,
    pub coord: [f64; 3],
    pub tol: [f64; 3],
    /// Memoized `hp_coord(JUDGE_PREC)` — the astro-float realization of the chain, computed
    /// once per definition and shared across clones (`Rc`). A judge escalates the *same*
    /// definition-point dozens of times per boolean (`plane_def` clones `tri_pt3` per call);
    /// without this each escalation replays the rotation's cos/sin at 200 bits, which dominates
    /// the rotated-boolean cost. `base`/`chain` never change after construction except through
    /// [`rotate_about`], which resets this cell, so the cached value always matches the
    /// definition (a pure, path-independent function). See [`HpCell`] for the
    /// `Arc<OnceLock>` (parallel) vs `Rc<OnceCell>` (serial) choice.
    hp: HpCell,
}

impl Pt3 {
    /// A point at `base`, tol seeded with the base→f64 rounding (a division; exactly
    /// 0 for an f64-representable base, positive otherwise). An axis a later chain
    /// never rotates keeps exactly this, which a per-axis tol check needs.
    pub fn at(base: [Rat; 3]) -> Self {
        let coord = [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()];
        // Actual rounding = |coord − base| at high precision; 2× the power-of-two
        // magnitude is a sound upper bound (0 when the base is exact).
        let round_tol = |b: Rat, c: f64| {
            let e = BigFloat::from_f64(c, 120).sub(&rat_to_big(b, 120), 120, HP_RM);
            2.0 * bf_mag(&e.abs())
        };
        Self::at_with_tol(
            base,
            [
                round_tol(base[0], coord[0]),
                round_tol(base[1], coord[1]),
                round_tol(base[2], coord[2]),
            ],
        )
    }

    /// A point at `base` with an explicit initial `tol` — for a root that already
    /// carries tol (a `Discovered` boolean seam), whose tol the chain then transports
    /// (`|R|·old`). [`at`](Self::at) is the Constructed case (initial tol = base
    /// rounding).
    pub fn at_with_tol(base: [Rat; 3], tol: [f64; 3]) -> Self {
        Pt3 {
            coord: [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
            base,
            chain: Vec::new(),
            tol,
            hp: HpCell::default(),
        }
    }

    /// Extend the definition by a rotation about `axis` through the origin. See
    /// [`rotate_about`](Self::rotate_about).
    pub fn rotate(self, axis: Axis, angle: Angle) -> Self {
        self.rotate_about(axis, angle, [Rat::from_int(0); 3])
    }

    /// Extend the definition by a rotation about `axis` through the rational pivot
    /// `point`, updating the f64 cache and its direction-wise tol. The two in-plane
    /// coords are taken relative to the pivot, rotated, and shifted back. Same-axis
    /// 90°-family angles about the origin rotate exactly (tol 0, Niven). A non-origin
    /// pivot adds its own f64 rounding (`coord − p`, `p + …`, the pivot's Rat→f64) —
    /// present even for an exact angle — covered by the `piv` term (validated H-f).
    pub fn rotate_about(mut self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let (i, j) = axis.plane();
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
        // Rotation-realization error (coordinate-mixing), 0 for an exact angle; plus
        // the pivot arithmetic (subtract p, re-add p, the pivot's own rounding),
        // exactly 0 for an origin pivot and present even for an exact angle otherwise.
        let rot = if exact {
            0.0
        } else {
            (u.abs() + v.abs()) * DA_F64
        };
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
        // The definition changed — invalidate the memoized hp of the old definition. A fresh
        // (unshared) cell, so clones made before this rotation keep their own cached value.
        self.hp = HpCell::default();
        self
    }

    /// The coordinate realized at `prec` bits from the **definition** (base rotated
    /// through the chain, each node about its pivot) — path-independent ground truth /
    /// escalation realization. At [`JUDGE_PREC`] (every escalation) the result is memoized in
    /// [`Pt3::hp`] and shared across clones of the same definition, so a definition-point pays
    /// the astro-float cos/sin once per boolean rather than once per predicate.
    pub fn hp_coord(&self, prec: usize) -> [BigFloat; 3] {
        if prec == JUDGE_PREC {
            return self.hp.get_or_init(|| self.compute_hp(prec)).clone();
        }
        self.compute_hp(prec)
    }

    /// The uncached realization (the body of [`hp_coord`]).
    fn compute_hp(&self, prec: usize) -> [BigFloat; 3] {
        let mut p = [
            rat_to_big(self.base[0], prec),
            rat_to_big(self.base[1], prec),
            rat_to_big(self.base[2], prec),
        ];
        for node in &self.chain {
            let (i, j) = node.axis.plane();
            let (c, s) = node.angle.cos_sin_at(prec);
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
/// tol (§CIP ②). The determinant is six signed triple-products of the edge entries;
/// each entry `(a−d)[k]` carries tol `tol_a[k] + tol_d[k]`. The bound sums the six
/// product radii (input-tol propagation, triangle-inequality worst case) plus a term
/// for the f64 rounding of the determinant's own arithmetic. The 3D analogue of
/// [`crate::frame2`]'s 2D `det_bound`. Validated in exact3d (H-a).
fn det3_bound(p: [[f64; 3]; 4], t: [[f64; 3]; 4]) -> f64 {
    let r = rows(p[0], p[1], p[2], p[3]);
    let td = t[3]; // apex tol adds to every edge on subtraction
    let t = [
        [t[0][0] + td[0], t[0][1] + td[1], t[0][2] + td[2]],
        [t[1][0] + td[0], t[1][1] + td[1], t[1][2] + td[2]],
        [t[2][0] + td[0], t[2][1] + td[1], t[2][2] + td[2]],
    ];
    // The six signed triple products of the cofactor expansion (indices into r/t).
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
    // f64 rounding of the ~17-op determinant, bounded generously (2D uses 8ε).
    input_tol + 16.0 * f64::EPSILON * mag
}

/// 3×3 determinant of `BigFloat` rows at `prec` bits (astro-float). The row-space
/// analogue of [`det3_iv`]; both [`det3_hp`] (edge rows) and [`dir_orient3d_judge`]
/// (a direction row) build their rows and call this.
fn det3_big(r: [[BigFloat; 3]; 3], prec: usize) -> BigFloat {
    let mul = |x: &BigFloat, y: &BigFloat| x.mul(y, prec, HP_RM);
    let m0 = mul(&r[1][1], &r[2][2]).sub(&mul(&r[1][2], &r[2][1]), prec, HP_RM);
    let m1 = mul(&r[1][0], &r[2][2]).sub(&mul(&r[1][2], &r[2][0]), prec, HP_RM);
    let m2 = mul(&r[1][0], &r[2][1]).sub(&mul(&r[1][1], &r[2][0]), prec, HP_RM);
    mul(&r[0][0], &m0)
        .sub(&mul(&r[0][1], &m1), prec, HP_RM)
        .add(&mul(&r[0][2], &m2), prec, HP_RM)
}

/// `orient3d` determinant realized at `prec` bits (astro-float) from the point
/// definitions — path-independent ground truth / escalation realization.
fn det3_hp(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3, prec: usize) -> BigFloat {
    let (a, b, c, d) = (
        pa.hp_coord(prec),
        pb.hp_coord(prec),
        pc.hp_coord(prec),
        pd.hp_coord(prec),
    );
    let sub = |x: &BigFloat, y: &BigFloat| x.sub(y, prec, HP_RM);
    det3_big(
        [
            [sub(&a[0], &d[0]), sub(&a[1], &d[1]), sub(&a[2], &d[2])],
            [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
            [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
        ],
        prec,
    )
}

/// The magnitude scale of four points (for the declare-0 floor).
fn scale4(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let m = |p: [f64; 3]| p[0].abs().max(p[1].abs()).max(p[2].abs());
    m(a).max(m(b)).max(m(c)).max(m(d)).max(1.0)
}

/// CIP `orient3d`: f64 filter (`|det| > bound` → trust the sign), else escalate to
/// astro-float at [`JUDGE_PREC`]; a determinant below the precision floor (`~scale³`,
/// cubic for the 3×3 case) is [`Orient::Zero`] (declare-0). Path-independent (a
/// function of the four point definitions). This is the *tol > 0* path — a tol-0
/// config is exact/faster via `nacre-predicates` (Shewchuk), routed above this crate.
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
    // Sign via is_positive (astro-float#44 workaround — see crate::frame2).
    if dh.is_zero() || bf_mag(&dh) <= floor {
        Orient::Zero
    } else if dh.is_positive() {
        Orient::Positive
    } else {
        Orient::Negative
    }
}

/// CIP `dir_orient3d`: the sign of `det[d, x−base, y−base] = d·((x−base)×(y−base))` — the
/// orientation of the ray direction `d` against the edge fan `(base→x, base→y)`. The
/// **direction analogue** of [`orient3d_judge`]: `point_in_solid`'s ray-triangle test asks
/// `orient3d(p, p+d, ·, ·)`, but `p+d` (a rotated point plus a rational offset) has no exact
/// `base+chain` `Pt3` (`R⁻¹d` is irrational). Every such determinant reduces to this form,
/// where `d` enters as one **exact** (rad-0) column and only `x, y, base` carry rotation tol.
/// Interval filter → astro-float escalation, exactly like the indirect judges; a below-floor
/// determinant is [`Orient::Zero`] (declare-0, absorbed by the caller's ray retry).
///
/// `d` is realized through an unrotated [`Pt3`] purely to reuse `pt_iv`/`hp_coord` for its
/// coord/tol/hp — the row is `d` itself, never `d − base`.
pub fn dir_orient3d_judge(d: [Rat; 3], base: &Pt3, x: &Pt3, y: &Pt3) -> Orient {
    let dp = Pt3::at(d);
    let (bi, xi, yi, di) = (pt_iv(base), pt_iv(x), pt_iv(y), pt_iv(&dp));
    let sub_iv = |u: [Iv; 3], v: [Iv; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
    if let Some(pos) = det3_iv([di, sub_iv(xi, bi), sub_iv(yi, bi)]).sign() {
        return orient_of(pos);
    }
    // Escalate: the same determinant at JUDGE_PREC from the exact definitions.
    let (bh, xh, yh, dh) = (
        base.hp_coord(JUDGE_PREC),
        x.hp_coord(JUDGE_PREC),
        y.hp_coord(JUDGE_PREC),
        dp.hp_coord(JUDGE_PREC),
    );
    let sub_hp = |u: &BigFloat, v: &BigFloat| u.sub(v, JUDGE_PREC, HP_RM);
    let rows = [
        [dh[0].clone(), dh[1].clone(), dh[2].clone()],
        [
            sub_hp(&xh[0], &bh[0]),
            sub_hp(&xh[1], &bh[1]),
            sub_hp(&xh[2], &bh[2]),
        ],
        [
            sub_hp(&yh[0], &bh[0]),
            sub_hp(&yh[1], &bh[1]),
            sub_hp(&yh[2], &bh[2]),
        ],
    ];
    let det = det3_big(rows, JUDGE_PREC);
    // Declare-0 floor: the six |triple products| of the f64 rows (term-magnitude, as the
    // indirect judges — `d` exact, edge rows from the f64 coords).
    let sub_f = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    let mag = det3_mag([
        dp.coord,
        sub_f(x.coord, base.coord),
        sub_f(y.coord, base.coord),
    ]);
    match sign_with_floor(&det, mag, JUDGE_PREC) {
        Some(pos) => orient_of(pos),
        None => Orient::Zero,
    }
}

/// CIP `orient3d(base, base+dir, x, y)` — an `orient3d` whose 2nd point is the ideal point
/// in direction `dir` (a ray from `base`). This is exactly the shape `ray_triangle_cross`'s
/// edge tests take (`orient3d(p, p+d, ·, ·)`), so a caller mirrors the f64 predicate
/// argument-for-argument instead of hand-reducing the direction column.
///
/// The determinant `det[base−y, (base+dir)−y, x−y]` column-reduces (`R1−R0 = dir`) and, after
/// the two swaps that move `dir` to the front, is `dir·((x−y)×(base−y))` — i.e. a plain
/// [`dir_orient3d_judge`] with the points permuted, no sign fix needed.
pub fn orient3d_ray(base: &Pt3, dir: [Rat; 3], x: &Pt3, y: &Pt3) -> Orient {
    dir_orient3d_judge(dir, y, x, base)
}

/// Sum of the six `|triple products|` of a 3×3 f64 matrix's rows — the term-magnitude
/// scale for [`dir_orient3d_judge`]'s declare-0 floor (cf. [`det3_bound`]'s `mag`).
fn det3_mag(r: [[f64; 3]; 3]) -> f64 {
    [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ]
    .iter()
    .map(|c| (r[0][c[0]] * r[1][c[1]] * r[2][c[2]]).abs())
    .sum()
}

// ---- indirect orient3d: three rotated planes meet at an implicit point ----
//
// A `Discovered` seam vertex is `∩` of three planes, each a plane through three
// rotated points, so its coordinates are never exact — an **implicit point** `V`.
// The kernel decides `orient3d(V, q, r, s)` without materializing `V`:
// `sign = sign(D)·sign(M)` where `D = det(normals)` and
// `M = (Dvec − D·s)·((q−s)×(r−s))` (Cramer, no division). The f64 filter is over
// **intervals** (value ± tol) — the design's "dynamic filter" (§CIP ⑨): interval
// arithmetic is a sound worst-case bound by construction, so no per-predicate bound
// formula is hand-derived. An interval straddling 0 escalates to astro-float from the
// point definitions, exactly as [`orient3d_judge`]. This is the *tol > 0* path, the
// judgment **layer only** ("층만") — the boolean wiring (seam → plane `Pt3`s) is
// stage 3. Validated in `experiments/exact3d` (H-b coefficient tol, H-c indirect
// soundness); the constant `mag`-floor policy is indirect-only (distinct from the
// explicit `16·scale³` floor of [`orient3d_judge`]).

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
        // Worst-case product radius `|a|·rad_b + |b|·rad_a + rad_a·rad_b`, plus the
        // f64 rounding of the product itself.
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

/// A point's coord+tol as an interval per component (the pivot is already folded into
/// `coord`/`tol` by [`Pt3::rotate_about`]).
fn pt_iv(p: &Pt3) -> [Iv; 3] {
    [
        Iv::new(p.coord[0], p.tol[0]),
        Iv::new(p.coord[1], p.tol[1]),
        Iv::new(p.coord[2], p.tol[2]),
    ]
}

/// Plane `[a,b,c,d]` (`n·X + d = 0`) through three points, as intervals: `n =
/// (p1−p0)×(p2−p0)`, `d = −n·p0`. Coefficient tol propagates from the point tols
/// through the subtraction/cross/dot — "coefficient tol is a corollary of point tol"
/// (§CIP ②, design.md §539). Validated H-b.
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
/// (ground truth for the coefficient tol; no trig of its own — consumes point
/// `hp_coord`).
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

/// The interval Cramer parts of an implicit point `V = ∩(planes)`: `D = det(normals)`
/// and the numerator vector `Dvec` (column `j` replaced by `h = −d`), so `V[j] =
/// Dvec[j]/D` (no division taken here). Shared by [`indirect_filter`] (orient3d) and
/// [`cmp_filter`] (cmp_coord).
fn cramer_iv(planes: [[Iv; 4]; 3]) -> (Iv, [Iv; 3]) {
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
    (d, dvec)
}

/// The interval f64 filter for `orient3d(V, q, r, s)`, `V = ∩(planes)`. `None` if
/// either `D` or `M` straddles 0 (escalate). Coefficient-direct (no division).
fn indirect_filter(planes: [[Iv; 4]; 3], q: [Iv; 3], r: [Iv; 3], s: [Iv; 3]) -> Option<Orient> {
    let (d, dvec) = cramer_iv(planes);
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

/// The astro-float Cramer parts of `V = ∩(planes)` at `prec` bits: `(D, Dvec, mag_d,
/// mag_dvec)` — the determinant, its numerator vector, and the cancellation-free
/// magnitude bounds of each (for the declare-0 floors). Shared by [`indirect_hp`]
/// (orient3d, which ignores `mag_dvec`) and [`cmp_hp`] (cmp_coord).
fn cramer_hp(planes: [[BigFloat; 4]; 3], prec: usize) -> (BigFloat, [BigFloat; 3], f64, [f64; 3]) {
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
    let cols = [
        [
            [hc[0].clone(), n0[1].clone(), n0[2].clone()],
            [hc[1].clone(), n1[1].clone(), n1[2].clone()],
            [hc[2].clone(), n2[1].clone(), n2[2].clone()],
        ],
        [
            [n0[0].clone(), hc[0].clone(), n0[2].clone()],
            [n1[0].clone(), hc[1].clone(), n1[2].clone()],
            [n2[0].clone(), hc[2].clone(), n2[2].clone()],
        ],
        [
            [n0[0].clone(), n0[1].clone(), hc[0].clone()],
            [n1[0].clone(), n1[1].clone(), hc[1].clone()],
            [n2[0].clone(), n2[1].clone(), hc[2].clone()],
        ],
    ];
    let dvec = [det3(&cols[0]), det3(&cols[1]), det3(&cols[2])];
    let mag_dvec = [det3_mag(&cols[0]), det3_mag(&cols[1]), det3_mag(&cols[2])];
    (d, dvec, mag_d, mag_dvec)
}

/// The same `sign(D)·sign(M)` realized at `prec` bits (astro-float) — ground truth /
/// escalation. Returns `(D, M, mag_d, mag_m)`: the two determinants and the f64
/// magnitudes of their term sums, so the caller can floor a below-precision result to
/// declare-0 ([`sign_with_floor`]).
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
    let (d, dvec, mag_d, _mag_dvec) = cramer_hp(planes, prec);
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

/// Multiplier over `mag·2⁻ᵖʳᵉᶜ` for the indirect declare-0 floor — generously above
/// the accumulated rounding of the ~100-op indirect computation (soundness first; a
/// tighter constant would only reduce the rare declare-0 rate, never soundness). This
/// is the **indirect-only** floor policy (term-magnitude based), distinct from the
/// explicit `16·scale³` floor of [`orient3d_judge`].
const FLOOR_K: f64 = 1.0e6;

/// A `BigFloat`'s sign as `Some(is_positive)`, floored to declare-0 (`None`) when its
/// magnitude is at or below the rounding floor for a `prec`-bit computation whose terms
/// scale to `mag`. Reuses [`orient3d_judge`]'s inline sign idiom (`is_zero`/`bf_mag`/
/// `is_positive`, astro-float#44 workaround) — one sign-extraction path in frame3.
fn sign_with_floor(val: &BigFloat, mag: f64, prec: usize) -> Option<bool> {
    let floor = FLOOR_K * mag * 2f64.powi(-(prec as i32));
    if val.is_zero() || bf_mag(val) <= floor {
        None
    } else {
        Some(val.is_positive())
    }
}

/// CIP indirect `orient3d(V, q, r, s)`, `V = ∩(3 planes)` — each plane through three
/// rotated points, the triangle three rotated points. Interval filter → astro-float
/// escalation; a below-floor `D` or `M` is [`Orient::Zero`] (declare-0, §6 "ask the
/// user"). The `Discovered`-seam analogue of [`orient3d_judge`]; boolean wiring is
/// stage 3.
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
    // JUDGE_PREC computation; a below-floor `D` or `M` is declare-0.
    let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, JUDGE_PREC);
    let (d, m, mag_d, mag_m) = indirect_hp(
        [ph(plane_a), ph(plane_b), ph(plane_c)],
        q.hp_coord(JUDGE_PREC),
        r.hp_coord(JUDGE_PREC),
        s.hp_coord(JUDGE_PREC),
        JUDGE_PREC,
    );
    combine(
        sign_with_floor(&d, mag_d, JUDGE_PREC),
        sign_with_floor(&m, mag_m, JUDGE_PREC),
    )
    .unwrap_or(Orient::Zero)
}

// ---- indirect cmp_coord: order two implicit points along one axis ----
//
// Each implicit point's axis coordinate is `Dvec[axis]/D` (Cramer), so
// `sign(a[axis] − b[axis]) = sign(Dvec_a[axis]·D_b − Dvec_b[axis]·D_a)·sign(D_a)·
// sign(D_b)` — three exact signs, no division, matching `nacre_predicates::
// indirect_cmp_coord`. Interval filter → astro-float escalation, like orient3d. The
// result is invariant to each triple's plane-normal orientation: flipping a triple's
// normals negates both its `D` and `Dvec[axis]` (so the point is unchanged) and negates
// `M`, leaving `sign(M)·sign(D)` fixed — so a caller may define each plane by any three
// non-collinear points on it (the bridge relies on this).

/// Combine the three definite signs of `M·D_a·D_b` (`true` = positive) into an
/// ordering: `Positive` (`a[axis] > b[axis]`) iff an even number are negative; `None`
/// if any sign is indefinite.
fn cmp_combine(sm: Option<bool>, sda: Option<bool>, sdb: Option<bool>) -> Option<Orient> {
    match (sm, sda, sdb) {
        (Some(m), Some(a), Some(b)) => {
            let negatives = [m, a, b].iter().filter(|&&s| !s).count();
            Some(if negatives % 2 == 0 {
                Orient::Positive
            } else {
                Orient::Negative
            })
        }
        _ => None,
    }
}

/// The interval f64 filter for `cmp_coord(a, b, axis)` — the sign of `a[axis] −
/// b[axis]` between two implicit points. `None` if any of `M`, `D_a`, `D_b` straddles 0.
fn cmp_filter(a: [[Iv; 4]; 3], b: [[Iv; 4]; 3], axis: usize) -> Option<Orient> {
    let (da, dva) = cramer_iv(a);
    let (db, dvb) = cramer_iv(b);
    let m = dva[axis].mul(db).sub(dvb[axis].mul(da));
    cmp_combine(m.sign(), da.sign(), db.sign())
}

/// The astro-float escalation for `cmp_coord`: the same three signs at `prec` bits,
/// each floored to declare-0 ([`sign_with_floor`]); `Orient::Zero` if any is below its
/// floor (the two coordinates are equal or too close to separate).
fn cmp_hp(a: [[BigFloat; 4]; 3], b: [[BigFloat; 4]; 3], axis: usize, prec: usize) -> Orient {
    let (da, dva, mag_da, mag_dva) = cramer_hp(a, prec);
    let (db, dvb, mag_db, mag_dvb) = cramer_hp(b, prec);
    let m = dva[axis]
        .mul(&db, prec, HP_RM)
        .sub(&dvb[axis].mul(&da, prec, HP_RM), prec, HP_RM);
    // Sound bound on |M|'s two term magnitudes (cancellation-free Dvec bound × |D|).
    let mag_m = mag_dva[axis] * bf_mag(&db) + mag_dvb[axis] * bf_mag(&da);
    cmp_combine(
        sign_with_floor(&m, mag_m, prec),
        sign_with_floor(&da, mag_da, prec),
        sign_with_floor(&db, mag_db, prec),
    )
    .unwrap_or(Orient::Zero)
}

/// CIP indirect `cmp_coord`: the sign of `a[axis] − b[axis]` where `a`, `b` are the
/// implicit points at which each three-plane triple meets — each plane through three
/// rotated points. Interval filter → astro-float escalation. `Positive` = `a[axis] >
/// b[axis]`, `Negative` = `<`, `Zero` = equal **or** below the declare-0 floor (unlike
/// the exact `nacre_predicates::indirect_cmp_coord`, whose `0` means exactly equal). The
/// two-implicit companion of [`indirect_orient3d_judge`]; boolean wiring is stage 3.
pub fn indirect_cmp_coord_judge(
    a: [(&Pt3, &Pt3, &Pt3); 3],
    b: [(&Pt3, &Pt3, &Pt3); 3],
    axis: usize,
) -> Orient {
    let iv = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
        [
            plane_iv(t[0].0, t[0].1, t[0].2),
            plane_iv(t[1].0, t[1].1, t[1].2),
            plane_iv(t[2].0, t[2].1, t[2].2),
        ]
    };
    if let Some(o) = cmp_filter(iv(a), iv(b), axis) {
        return o;
    }
    let hp = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
        [
            plane_hp(t[0].0, t[0].1, t[0].2, JUDGE_PREC),
            plane_hp(t[1].0, t[1].1, t[1].2, JUDGE_PREC),
            plane_hp(t[2].0, t[2].1, t[2].2, JUDGE_PREC),
        ]
    };
    cmp_hp(hp(a), hp(b), axis, JUDGE_PREC)
}

/// A definite `bool` sign (`true` = positive) as an [`Orient`].
fn orient_of(pos: bool) -> Orient {
    if pos {
        Orient::Positive
    } else {
        Orient::Negative
    }
}

/// The sign of `det[n_a; n_b; n_c]`, the determinant of the three planes' normals — the
/// toleranced twin of `nacre_geom`'s `plane_pair_dir_sign` (`sign((n_a × n_b) · n_c)`),
/// which decides how the line `a ∩ b` runs relative to plane `c`. This is exactly the
/// Cramer `D` of the three planes ([`cramer_iv`] / [`cramer_hp`]); the interval filter
/// resolves the easy cases, an ambiguous one escalates, and a below-floor `D` is
/// [`Orient::Zero`] (the planes' normals are coincident — degenerate).
///
/// **Winding-dependent** (unlike [`orient3d_judge`] / [`indirect_cmp_coord_judge`], which
/// are normal-orientation invariant): each plane's normal is `(p1−p0)×(p2−p0)`, so the
/// result follows each triple's point order. The caller must pass the three points in a
/// consistent order (the ops wrapper passes each face's outward-oriented `tri`).
pub fn dir_sign_judge(
    a: (&Pt3, &Pt3, &Pt3),
    b: (&Pt3, &Pt3, &Pt3),
    c: (&Pt3, &Pt3, &Pt3),
) -> Orient {
    let (d, _) = cramer_iv([
        plane_iv(a.0, a.1, a.2),
        plane_iv(b.0, b.1, b.2),
        plane_iv(c.0, c.1, c.2),
    ]);
    if let Some(pos) = d.sign() {
        return orient_of(pos);
    }
    let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, JUDGE_PREC);
    let (dh, _, mag_d, _) = cramer_hp([ph(a), ph(b), ph(c)], JUDGE_PREC);
    match sign_with_floor(&dh, mag_d, JUDGE_PREC) {
        Some(pos) => orient_of(pos),
        None => Orient::Zero,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ri(n: i128, d: i128) -> Rat {
        Rat::new(n, d).unwrap()
    }

    fn deg(n: i128, d: i128) -> Angle {
        Angle::from_deg(ri(n, d)).unwrap()
    }

    /// Deterministic PRNG (splitmix64) for reproducible stress corpora.
    fn splitmix64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn rng(state: &mut u64, lo: i128, hi: i128) -> i128 {
        lo + (u128::from(splitmix64(state)) % (hi - lo + 1) as u128) as i128
    }
    fn axis_of(k: i128) -> Axis {
        match k.rem_euclid(3) {
            0 => Axis::X,
            1 => Axis::Y,
            _ => Axis::Z,
        }
    }
    fn rand_base(st: &mut u64) -> [Rat; 3] {
        [
            ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
            ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
            ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
        ]
    }
    fn abs_err(f: f64, truth: &BigFloat, gt: usize) -> f64 {
        bf_mag(&BigFloat::from_f64(f, gt).sub(truth, gt, HP_RM).abs())
    }

    /// Soundness: over random rotation chains (mixed axes, arbitrary pivots, exact and
    /// inexact angles), the direction-wise tol must bound the true f64 error on every
    /// axis (astro-float 512-bit ground truth) — the production mirror of exact3d
    /// H-d/H-f. (An `err` below 1e-100 is 512-bit GT noise, not a real f64 error.)
    /// `#[ignore]`: astro-float ground truth is slow; run with `--ignored` (+ CI). The
    /// full statistical validation lives in `experiments/exact3d`.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn tol_bounds_error_over_random_chains() {
        const GT: usize = 512;
        let mut st = 0x2A5C_1234_ABCD_9999u64;
        let (mut exact_seen, mut pivot_seen) = (false, false);
        for _ in 0..2000 {
            let base = rand_base(&mut st);
            let mut p = Pt3::at(base);
            for _ in 0..rng(&mut st, 2, 5) {
                let ax = axis_of(rng(&mut st, 0, 2));
                let pivot = match rng(&mut st, 0, 2) {
                    0 => [Rat::from_int(0); 3],
                    1 => {
                        pivot_seen = true;
                        rand_base(&mut st)
                    }
                    _ => {
                        pivot_seen = true;
                        base // point on the pivot
                    }
                };
                let ang = if rng(&mut st, 0, 3) == 0 {
                    exact_seen = true;
                    deg(90 * rng(&mut st, 0, 3), 1)
                } else {
                    deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973))
                };
                p = p.rotate_about(ax, ang, pivot);
            }
            let hp = p.hp_coord(GT);
            for (axis, hp_a) in hp.iter().enumerate() {
                let err = abs_err(p.coord[axis], hp_a, GT);
                assert!(
                    err <= p.tol[axis] || err < 1e-100,
                    "tol must bound the error: axis {axis}, err {err:e} > tol {:e}",
                    p.tol[axis]
                );
            }
        }
        assert!(
            exact_seen && pivot_seen,
            "corpus must mix exact angles and pivots"
        );
    }

    /// A 90°-family chain about the origin stays exactly tol 0 (Niven exact realization).
    #[test]
    fn quadrantal_origin_chain_is_tol_zero() {
        let p = Pt3::at([ri(3, 1), ri(5, 1), ri(7, 1)])
            .rotate(Axis::Z, deg(90, 1))
            .rotate(Axis::X, deg(180, 1))
            .rotate(Axis::Y, deg(270, 1));
        assert_eq!(p.tol, [0.0; 3]);
        // and the realized coords are the exact permutation/negation (no spurious term).
        assert_eq!(p.coord, [7.0, -3.0, -5.0]);
    }

    /// `at` seeds the base→f64 rounding: 0 for an integer base, positive for a base
    /// that is not f64-representable (e.g. 1/3).
    #[test]
    fn at_seeds_base_rounding_tol() {
        assert_eq!(Pt3::at([ri(2, 1), ri(3, 1), ri(4, 1)]).tol, [0.0; 3]);
        let third = Pt3::at([ri(1, 3), ri(0, 1), ri(0, 1)]);
        assert!(third.tol[0] > 0.0 && third.tol[1] == 0.0);
    }

    /// `at_with_tol` seeds a nonzero root tol (a Discovered seam) and the chain
    /// transports it (`|R|·old`): a 90° rotation about Z swaps the x/y tol components.
    #[test]
    fn seeded_tol_transports_through_rotation() {
        let p = Pt3::at_with_tol([ri(1, 1), ri(0, 1), ri(0, 1)], [1e-9, 2e-9, 3e-9])
            .rotate(Axis::Z, deg(90, 1));
        // 90° about Z: |R| swaps x,y → tol[0]=old tol[1], tol[1]=old tol[0]; z unchanged.
        assert_eq!(p.tol, [2e-9, 1e-9, 3e-9]);
    }

    /// Same-axis bundling is tighter than incremental (H-e): K steps of θ amplify the
    /// tol (each transported by |c|+|s| ≥ 1) vs one step of Kθ.
    #[test]
    fn bundling_is_tighter_than_incremental() {
        let base = [ri(11, 1), ri(-7, 1), ri(4, 1)];
        let (k, theta) = (30i128, ri(1, 7));
        let mut incr = Pt3::at(base);
        for _ in 0..k {
            incr = incr.rotate(Axis::Z, Angle::from_deg(theta).unwrap());
        }
        let k_theta = theta.checked_mul(Rat::from_int(k)).unwrap();
        let bundled = Pt3::at(base).rotate(Axis::Z, Angle::from_deg(k_theta).unwrap());
        assert!(
            bundled.tol[0] < incr.tol[0] && bundled.tol[1] < incr.tol[1],
            "bundled {:?} must be tighter than incremental {:?}",
            bundled.tol,
            incr.tol
        );
    }

    // ---- orient3d_judge (2a-ii) ----

    /// A point rotated about one random axis (through a random pivot) by one **inexact**
    /// rational angle — generic non-degenerate, heterogeneous provenance (each point its
    /// own rotation), the H-a corpus shape. Inexact angles only (an exact/near-coplanar
    /// mix would expose the tol-0 GT noise floor — that is the declare-0 test's job).
    fn rand_point(st: &mut u64) -> Pt3 {
        let base = rand_base(st);
        let axis = axis_of(rng(st, 0, 2));
        let ang = deg(rng(st, 0, 360_000), rng(st, 1, 9973)); // inexact
        let pivot = if rng(st, 0, 1) == 0 {
            [Rat::from_int(0); 3]
        } else {
            rand_base(st) // non-origin pivot: exercises the pivot tol → det3_bound path
        };
        Pt3::at(base).rotate_about(axis, ang, pivot)
    }

    /// H-a — `det3_bound` soundness: over many random heterogeneous-rotation 4-point
    /// configs (inexact angles, arbitrary pivots), the bound must upper-bound the real
    /// error of the f64 determinant vs the astro-float truth — never exceeded. The
    /// production mirror of exact3d H-a. `#[ignore]`: slow (4× astro-float per sample).
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_a_det3_bound_soundness() {
        const GT: usize = 512;
        let mut st = 0x3D00_1234_ABCD_EF01u64;
        let (mut bad, mut worst) = (0usize, 0.0_f64);
        for _ in 0..5000 {
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
                worst = worst.max(err / bound);
            }
        }
        eprintln!("[H-a N=5000] det3_bound violations: {bad}; worst tightness: {worst:.3}");
        assert_eq!(
            bad, 0,
            "det3_bound must bound the determinant error on every sample"
        );
    }

    /// Fast port check (default suite): a handful of random configs must not violate
    /// `det3_bound` (catches a transcription bug; the full statistical soundness is the
    /// `#[ignore]`d H-a + exact3d).
    #[test]
    fn det3_bound_port_check() {
        const GT: usize = 512;
        let mut st = 0x0BAD_F00D_1234_5678u64;
        for _ in 0..40 {
            let (pa, pb, pc, pd) = (
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
            );
            let det = det3_f64(pa.coord, pb.coord, pc.coord, pd.coord);
            let err = abs_err(det, &det3_hp(&pa, &pb, &pc, &pd, GT), GT);
            let bound = det3_bound(
                [pa.coord, pb.coord, pc.coord, pd.coord],
                [pa.tol, pb.tol, pc.tol, pd.tol],
            );
            assert!(err <= bound, "det3_bound {bound:e} < err {err:e}");
        }
    }

    /// A known tetrahedron judges `Positive`, its mirror `Negative`, and the sign is
    /// preserved under a shared rotation (rotating all four points does not flip
    /// orientation). Exact (tol 0) points take the filter's fast path.
    #[test]
    fn orient3d_sanity_and_rotation_invariance() {
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        // d at origin; a,b,c along +x,+y,+z → right-handed → Positive.
        let (a, b, c, d) = (pt(1, 0, 0), pt(0, 1, 0), pt(0, 0, 1), pt(0, 0, 0));
        assert_eq!(orient3d_judge(&a, &b, &c, &d), Orient::Positive);
        assert_eq!(
            orient3d_judge(&b, &a, &c, &d),
            Orient::Negative,
            "swap → mirror"
        );

        // Shared rotation of all four (37° about Z through a rational pivot) keeps the sign.
        let rot = |p: &Pt3| {
            p.clone()
                .rotate_about(Axis::Z, deg(37, 1), [ri(2, 1), ri(-3, 1), ri(0, 1)])
        };
        assert_eq!(
            orient3d_judge(&rot(&a), &rot(&b), &rot(&c), &rot(&d)),
            Orient::Positive,
            "orientation is rotation-invariant"
        );
    }

    /// Declare-0: four **exactly coplanar** rational points (tol 0) → the determinant is
    /// exactly 0 → `Orient::Zero` (the §6 ask-the-user case).
    #[test]
    fn orient3d_coplanar_is_zero() {
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        // all four in the plane z = 0.
        let (a, b, c, d) = (pt(0, 0, 0), pt(3, 0, 0), pt(0, 5, 0), pt(2, 7, 0));
        assert_eq!(orient3d_judge(&a, &b, &c, &d), Orient::Zero);
    }

    // ---- indirect orient3d (2c-i) ----

    /// `Iv` arithmetic keeps a sound radius (worst-case interval), and `sign` is
    /// definite exactly when the interval clears 0.
    #[test]
    fn iv_arithmetic_is_sound_and_sign_decides() {
        let (a, b) = (Iv::new(3.0, 0.1), Iv::new(-2.0, 0.2));
        // sub radius ≥ sum of input radii.
        assert!(a.sub(b).rad >= 0.1 + 0.2);
        // mul radius ≥ |mid_a|·rad_b + |mid_b|·rad_a (+ rad·rad + rounding).
        assert!(a.mul(b).rad >= 3.0 * 0.2 + 2.0 * 0.1);
        assert_eq!(Iv::new(1.0, 0.5).sign(), Some(true));
        assert_eq!(Iv::new(-1.0, 0.5).sign(), Some(false));
        assert_eq!(Iv::new(0.3, 0.5).sign(), None); // straddles 0 → escalate
    }

    /// Sanity: three coordinate planes meet at the origin `V=(0,0,0)`; `orient3d(V,q,r,s)`
    /// judges a known sign, flips on a q/r swap, is `Zero` when `s` is coplanar with
    /// `V,q,r`, and is invariant under a shared rotation (rotations preserve orientation).
    #[test]
    fn indirect_sanity_and_rotation_invariance() {
        fn tri(t: &(Pt3, Pt3, Pt3)) -> (&Pt3, &Pt3, &Pt3) {
            (&t.0, &t.1, &t.2)
        }
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        let pa = (pt(0, 0, 0), pt(1, 0, 0), pt(0, 1, 0)); // z = 0
        let pb = (pt(0, 0, 0), pt(1, 0, 0), pt(0, 0, 1)); // y = 0
        let pc = (pt(0, 0, 0), pt(0, 1, 0), pt(0, 0, 1)); // x = 0  → V = (0,0,0)
        let (q, r, s) = (pt(1, 0, 0), pt(0, 1, 0), pt(0, 0, 1));
        assert_eq!(
            indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &q, &r, &s),
            Orient::Negative
        );
        assert_eq!(
            indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &r, &q, &s),
            Orient::Positive,
            "q/r swap flips the sign"
        );
        // s coplanar with V,q,r (all z = 0) → orient exactly 0 → declare-0.
        let s0 = pt(1, 1, 0);
        assert_eq!(
            indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &q, &r, &s0),
            Orient::Zero
        );
        // Shared rotation of all twelve points keeps the definite sign.
        let rot = |p: &Pt3| {
            p.clone()
                .rotate_about(Axis::Z, deg(37, 1), [ri(2, 1), ri(-1, 1), ri(0, 1)])
        };
        let rp = |t: (&Pt3, &Pt3, &Pt3)| (rot(t.0), rot(t.1), rot(t.2));
        let (ra, rb, rc) = (rp(tri(&pa)), rp(tri(&pb)), rp(tri(&pc)));
        let (rq, rr, rs) = (rot(&q), rot(&r), rot(&s));
        assert_eq!(
            indirect_orient3d_judge(tri(&ra), tri(&rb), tri(&rc), &rq, &rr, &rs),
            Orient::Negative,
            "indirect orient is rotation-invariant"
        );
    }

    /// The high-precision indirect truth with a stability flag: `None` if even the
    /// ground truth cannot resolve `D` or `M` above its floor (a genuine degeneracy).
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
        combine(
            sign_with_floor(&d, mag_d, prec),
            sign_with_floor(&m, mag_m, prec),
        )
    }

    /// H-b — rotated plane coefficient tol soundness. A plane's four coefficients
    /// derive from three rotated points (subtraction/cross/dot); the interval `rad` on
    /// each must upper-bound the real f64 error vs the astro-float truth. Corpus mixes
    /// origin and arbitrary pivots (`rand_point`). `#[ignore]`: slow astro-float GT.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_b_plane_coefficient_tol_soundness() {
        const GT: usize = 512;
        const N: usize = 10_000;
        let mut st = 0xB0B0_5555_1111_2222u64;
        let (mut bad, mut worst) = (0usize, 0.0_f64);
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
        assert_eq!(bad, 0, "plane coefficient tol must bound the error");
    }

    /// H-c — indirect orient3d soundness over heterogeneous provenance (corpus A) and a
    /// near-coplanar escalation-forcing corpus (corpus B). The interval filter →
    /// astro-float judge must never claim a sign opposite to a GT-stable 512-bit truth.
    /// Asserts both paths are live: `escalated > 0` (filter defers) **and**
    /// `filter_resolved > 0` (fast path resolves — a filter stuck at `None` would pass
    /// wrong-sign 0 vacuously). `#[ignore]`: slow astro-float GT.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_c_indirect_orient3d_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
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
            } else {
                filter_resolved += 1;
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

        // Corpus A — heterogeneous provenance (each of the twelve points its own
        // rotation, origin or arbitrary pivot). Generic → mostly filter-resolved.
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

        // Corpus B — one shared rotation about a shared pivot (affine → coplanarity
        // preserved), near-coplanar so the escalation's sign-resolution is exercised.
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
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                ]
            };
            let v = [
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
            ];
            let (u, w) = (off(&mut st), off(&mut st));
            let n = cross(u, w); // normal to the V-plane
            // ε off-plane push: 0 (exactly degenerate), or a wide window down to ~1e-15
            // where the interval filter is ambiguous but the true sign is definite.
            let eps = match rng(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng(&mut st, 5_000, 200_000)),
                _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            let plane_pts = |st: &mut u64| [v, add3(v, off(st)), add3(v, off(st))];
            let a = plane_pts(&mut st);
            let b = plane_pts(&mut st);
            let c = plane_pts(&mut st);
            // Triangle q=v+u, r=v+w, s=v+u+w+ε·n (coplanar with V=v, s pushed ε off).
            let qb = add3(v, u);
            let rb = add3(v, w);
            let sb = add3(add3(add3(v, u), w), smul(eps, n));
            // One shared rotation (shared pivot) for all twelve points.
            let axis = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let rp = |p: [Rat; 3]| Pt3::at(p).rotate_about(axis, ang, piv);
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
            "[H-c] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "indirect judge must never disagree with GT (soundness)"
        );
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path (else a stuck filter passes vacuously)"
        );
    }

    /// Fast port check (default suite): a handful of generic configs — the indirect
    /// judge must match a moderate-precision truth (catches a transcription bug; full
    /// statistical soundness is the `#[ignore]`d H-b/H-c + exact3d).
    #[test]
    fn indirect_port_check() {
        const GT: usize = 384;
        let mut st = 0x1DEA_2C11_9F00_5A5Au64;
        for _ in 0..16 {
            let p: Vec<Pt3> = (0..12).map(|_| rand_point(&mut st)).collect();
            let a = (&p[0], &p[1], &p[2]);
            let b = (&p[3], &p[4], &p[5]);
            let c = (&p[6], &p[7], &p[8]);
            let judged = indirect_orient3d_judge(a, b, c, &p[9], &p[10], &p[11]);
            if let Some(truth) = indirect_truth(a, b, c, &p[9], &p[10], &p[11], GT) {
                assert!(
                    judged == truth || judged == Orient::Zero,
                    "indirect judge {judged:?} disagrees with truth {truth:?}"
                );
            }
        }
    }

    // ---- indirect cmp_coord (cmp-i) ----

    fn add3(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
        [
            a[0].checked_add(b[0]).unwrap(),
            a[1].checked_add(b[1]).unwrap(),
            a[2].checked_add(b[2]).unwrap(),
        ]
    }

    /// Borrow an owned plane-triple as the `&Pt3` tuples the judge takes.
    fn tr(t: &[[Pt3; 3]; 3]) -> [(&Pt3, &Pt3, &Pt3); 3] {
        [
            (&t[0][0], &t[0][1], &t[0][2]),
            (&t[1][0], &t[1][1], &t[1][2]),
            (&t[2][0], &t[2][1], &t[2][2]),
        ]
    }

    /// Three planes meeting at `v` (each through `v` + two small offsets), all rotated by
    /// `(ax, ang, piv)` — an implicit point at `rotate(v)` with heterogeneous provenance.
    fn triple_pts(v: [Rat; 3], st: &mut u64, ax: Axis, ang: Angle, piv: [Rat; 3]) -> [[Pt3; 3]; 3] {
        let plane = |st: &mut u64| {
            let off = |st: &mut u64| {
                [
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                ]
            };
            let (o1, o2) = (off(st), off(st));
            [
                Pt3::at(v).rotate_about(ax, ang, piv),
                Pt3::at(add3(v, o1)).rotate_about(ax, ang, piv),
                Pt3::at(add3(v, o2)).rotate_about(ax, ang, piv),
            ]
        };
        [plane(st), plane(st), plane(st)]
    }

    /// A random rotated three-plane triple (own random center, axis, inexact angle,
    /// pivot) — heterogeneous provenance. Sequences the RNG draws so each `&mut st`
    /// borrow ends before the next.
    fn rand_triple(st: &mut u64) -> [[Pt3; 3]; 3] {
        let v = rand_base(st);
        let ax = axis_of(rng(st, 0, 2));
        let angle = deg(rng(st, 0, 360_000), rng(st, 1, 9973));
        let piv = rand_base(st);
        triple_pts(v, st, ax, angle, piv)
    }

    /// The three axis-perpendicular planes through integer point `p` (meet exactly at `p`,
    /// tol 0) — an exact axis-aligned implicit point for the sanity oracle.
    fn axis_planes(p: [i128; 3]) -> [[Pt3; 3]; 3] {
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        let [x, y, z] = p;
        [
            [pt(x, y, z), pt(x, y + 1, z), pt(x, y, z + 1)], // ⊥ x
            [pt(x, y, z), pt(x + 1, y, z), pt(x, y, z + 1)], // ⊥ y
            [pt(x, y, z), pt(x + 1, y, z), pt(x, y + 1, z)], // ⊥ z
        ]
    }

    /// The high-precision cmp truth (`None` when the coordinates are equal or below the
    /// GT floor — a genuine tie).
    fn cmp_truth(
        a: [(&Pt3, &Pt3, &Pt3); 3],
        b: [(&Pt3, &Pt3, &Pt3); 3],
        axis: usize,
        prec: usize,
    ) -> Option<Orient> {
        let hp = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
            [
                plane_hp(t[0].0, t[0].1, t[0].2, prec),
                plane_hp(t[1].0, t[1].1, t[1].2, prec),
                plane_hp(t[2].0, t[2].1, t[2].2, prec),
            ]
        };
        match cmp_hp(hp(a), hp(b), axis, prec) {
            Orient::Zero => None,
            o => Some(o),
        }
    }

    /// `cmp_combine` maps the parity of negative signs to the ordering.
    #[test]
    fn cmp_combine_counts_negatives() {
        assert_eq!(
            cmp_combine(Some(true), Some(true), Some(true)),
            Some(Orient::Positive)
        );
        assert_eq!(
            cmp_combine(Some(false), Some(true), Some(true)),
            Some(Orient::Negative)
        );
        assert_eq!(
            cmp_combine(Some(false), Some(false), Some(true)),
            Some(Orient::Positive)
        );
        assert_eq!(cmp_combine(None, Some(true), Some(true)), None);
    }

    /// Sanity: two exact axis-aligned implicit points order by the compared axis; a swap
    /// flips it; an equal coordinate is `Zero`.
    #[test]
    fn cmp_sanity_axis_aligned() {
        let a = axis_planes([1, 2, 3]);
        let b = axis_planes([1, 5, 3]);
        // y: 2 < 5 → a below b → Negative; swap → Positive.
        assert_eq!(
            indirect_cmp_coord_judge(tr(&a), tr(&b), 1),
            Orient::Negative
        );
        assert_eq!(
            indirect_cmp_coord_judge(tr(&b), tr(&a), 1),
            Orient::Positive
        );
        // x and z equal → Zero.
        assert_eq!(indirect_cmp_coord_judge(tr(&a), tr(&b), 0), Orient::Zero);
        assert_eq!(indirect_cmp_coord_judge(tr(&a), tr(&b), 2), Orient::Zero);
    }

    /// Fast port check: the cmp judge matches a moderate-precision truth on generic
    /// rotated configs (catches a transcription bug; full soundness is `#[ignore]`d H-g).
    #[test]
    fn cmp_port_check() {
        const GT: usize = 384;
        let mut st = 0xC301_7A5E_2266_9911u64;
        for _ in 0..16 {
            let a = rand_triple(&mut st);
            let b = rand_triple(&mut st);
            let axis = rng(&mut st, 0, 2) as usize;
            let judged = indirect_cmp_coord_judge(tr(&a), tr(&b), axis);
            if let Some(truth) = cmp_truth(tr(&a), tr(&b), axis, GT) {
                assert!(
                    judged == truth || judged == Orient::Zero,
                    "cmp judge {judged:?} disagrees with truth {truth:?}"
                );
            }
        }
    }

    /// H-g — indirect cmp_coord soundness over heterogeneous provenance (corpus A) and a
    /// near-tie corpus (corpus B: two points sharing an axis coordinate up to ε, rotated
    /// about that same axis so the near-tie survives). The judge must never disagree with
    /// a GT-stable 512-bit truth; both the fast filter and the escalation are exercised.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_g_indirect_cmp_coord_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut check = |a: &[[Pt3; 3]; 3], b: &[[Pt3; 3]; 3], axis: usize| {
            let (ta, tb) = (tr(a), tr(b));
            let judged = indirect_cmp_coord_judge(ta, tb, axis);
            let iv = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
                [
                    plane_iv(t[0].0, t[0].1, t[0].2),
                    plane_iv(t[1].0, t[1].1, t[1].2),
                    plane_iv(t[2].0, t[2].1, t[2].2),
                ]
            };
            if cmp_filter(iv(ta), iv(tb), axis).is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match cmp_truth(ta, tb, axis, GT) {
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

        // Corpus A — heterogeneous provenance (each triple its own rotation and pivot).
        let mut st = 0x6A11_C0DE_5151_2323u64;
        for _ in 0..2000 {
            let a = rand_triple(&mut st);
            let b = rand_triple(&mut st);
            check(&a, &b, rng(&mut st, 0, 2) as usize);
        }

        // Corpus B — near-tie in z, shared Z-rotation (a Z-rotation leaves z unchanged, so
        // the pre-rotation z near-equality survives → M near 0 → escalation forced).
        for _ in 0..2000 {
            let z0 = rng(&mut st, -200, 200);
            let eps = match rng(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng(&mut st, 5_000, 200_000)),
                _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            let va = [
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
                ri(z0, 1),
            ];
            let vb = [
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
                ri(z0, 1).checked_add(eps).unwrap(),
            ];
            let angle = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let a = triple_pts(va, &mut st, Axis::Z, angle, piv);
            let b = triple_pts(vb, &mut st, Axis::Z, angle, piv);
            check(&a, &b, 2);
        }

        eprintln!(
            "[H-g] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "cmp judge must never disagree with GT (soundness)"
        );
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path"
        );
    }

    // ---- dir_sign_judge (3a-iii) ----

    fn rat_cross(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
        let m = |x: Rat, y: Rat| x.checked_mul(y).unwrap();
        let s = |x: Rat, y: Rat| x.checked_sub(y).unwrap();
        [
            s(m(a[1], b[2]), m(a[2], b[1])),
            s(m(a[2], b[0]), m(a[0], b[2])),
            s(m(a[0], b[1]), m(a[1], b[0])),
        ]
    }

    fn smul(k: Rat, a: [Rat; 3]) -> [Rat; 3] {
        [
            a[0].checked_mul(k).unwrap(),
            a[1].checked_mul(k).unwrap(),
            a[2].checked_mul(k).unwrap(),
        ]
    }

    /// Three base points defining a plane whose normal is parallel to `n` — `p0` and two
    /// in-plane edges `n × e1`, `n × e2` (single crosses, small magnitude).
    fn plane_norm(n: [Rat; 3], p0: [Rat; 3]) -> [[Rat; 3]; 3] {
        let e1 = [ri(1, 1), ri(2, 1), ri(3, 1)];
        let e2 = [ri(2, 1), ri(3, 1), ri(1, 1)];
        [p0, add3(p0, rat_cross(n, e1)), add3(p0, rat_cross(n, e2))]
    }

    fn rot_plane(b: [[Rat; 3]; 3], ax: Axis, ang: Angle, piv: [Rat; 3]) -> [Pt3; 3] {
        b.map(|p| Pt3::at(p).rotate_about(ax, ang, piv))
    }

    fn t3(p: &[Pt3; 3]) -> (&Pt3, &Pt3, &Pt3) {
        (&p[0], &p[1], &p[2])
    }

    /// The three points are far from collinear (`sin²` of the corner angle above a
    /// threshold) — a well-formed plane. A degenerate plane (tiny normal) is not the
    /// near-coplanar-*normals* regime under test, so the corpus skips it.
    fn well_conditioned(p: &[Pt3; 3]) -> bool {
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let (e1, e2) = (sub(p[1].coord, p[0].coord), sub(p[2].coord, p[0].coord));
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        dot(n, n) > 1e-6 * dot(e1, e1) * dot(e2, e2)
    }

    /// The sign of `D` (det of the three normals) at `prec` bits, floored to declare-0.
    fn dir_d_sign_at(
        a: (&Pt3, &Pt3, &Pt3),
        b: (&Pt3, &Pt3, &Pt3),
        c: (&Pt3, &Pt3, &Pt3),
        prec: usize,
    ) -> Option<bool> {
        let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, prec);
        let (dh, _, mag_d, _) = cramer_hp([ph(a), ph(b), ph(c)], prec);
        sign_with_floor(&dh, mag_d, prec)
    }

    /// The **GT-stable** `dir_sign` truth: `Some` only when `prec` and `prec + 128` agree
    /// on a definite sign — otherwise the config is degenerate beyond what the ground
    /// truth itself resolves, so it is `None` (skip), not a spurious "wrong". (Mirrors
    /// exact3d's `indirect_truth` stability check.)
    fn dir_sign_truth(
        a: (&Pt3, &Pt3, &Pt3),
        b: (&Pt3, &Pt3, &Pt3),
        c: (&Pt3, &Pt3, &Pt3),
        prec: usize,
    ) -> Option<Orient> {
        match (
            dir_d_sign_at(a, b, c, prec),
            dir_d_sign_at(a, b, c, prec + 128),
        ) {
            (Some(x), Some(y)) if x == y => Some(orient_of(x)),
            _ => None,
        }
    }

    /// Sanity: the three axis-perpendicular planes have normals `+x, −y, +z`, so their
    /// determinant is `−1`; a swap flips it, and three coplanar normals give `Zero`.
    #[test]
    fn dir_sign_judge_sanity() {
        let ap = axis_planes([0, 0, 0]);
        assert_eq!(
            dir_sign_judge(t3(&ap[0]), t3(&ap[1]), t3(&ap[2])),
            Orient::Negative,
            "det[+x, -y, +z] = -1"
        );
        assert_eq!(
            dir_sign_judge(t3(&ap[0]), t3(&ap[2]), t3(&ap[1])),
            Orient::Positive,
            "one swap flips the sign"
        );
        // Three normals in the plane z = 0 → coplanar → D = 0 → Zero.
        let mk = |n: [i128; 3]| {
            plane_norm([ri(n[0], 1), ri(n[1], 1), ri(n[2], 1)], [ri(0, 1); 3]).map(Pt3::at)
        };
        let (a, b, c) = (mk([1, 0, 0]), mk([0, 1, 0]), mk([1, 1, 0]));
        assert_eq!(
            dir_sign_judge(t3(&a), t3(&b), t3(&c)),
            Orient::Zero,
            "coplanar normals → D = 0"
        );
    }

    /// `dir_sign` is a determinant of normals, invariant under a shared rotation
    /// (`det(R·n) = det(R)·det(n) = det(n)`).
    #[test]
    fn dir_sign_rotation_invariant() {
        let mut st = 0x0D12_5157_ABCD_0007u64;
        for _ in 0..40 {
            let mk = |st: &mut u64| plane_norm(rand_base(st), rand_base(st));
            let (ba, bb, bc) = (mk(&mut st), mk(&mut st), mk(&mut st));
            let un = |b: [[Rat; 3]; 3]| b.map(Pt3::at);
            let (ua, ub, uc) = (un(ba), un(bb), un(bc));
            let s = dir_sign_judge(t3(&ua), t3(&ub), t3(&uc));
            if s == Orient::Zero {
                continue;
            }
            let ax = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let (ra, rb, rc) = (
                rot_plane(ba, ax, ang, piv),
                rot_plane(bb, ax, ang, piv),
                rot_plane(bc, ax, ang, piv),
            );
            assert_eq!(
                dir_sign_judge(t3(&ra), t3(&rb), t3(&rc)),
                s,
                "rotation-invariant"
            );
        }
    }

    /// H-i — dir_sign soundness over a **near-coplanar-normals** corpus (which H-c/H-g do
    /// not stress: they force `M ≈ 0`, not `D ≈ 0`). Three plane normals `n0, n1,
    /// n2 = α·n0 + β·n1 + ε·(n0×n1)` (ε tiny → `D ≈ ε` → escalation), shared rotation. The
    /// judge must never disagree with a GT-stable 512-bit `D` sign; both paths exercised.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_i_dir_sign_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut st = 0x0D18_5160_C0C0_2323u64;
        let rand_n = |st: &mut u64| {
            [
                ri(rng(st, -20, 20), 1),
                ri(rng(st, -20, 20), 1),
                ri(rng(st, -20, 20), 1),
            ]
        };
        for _ in 0..2000 {
            let n0 = rand_n(&mut st);
            let n1 = rand_n(&mut st);
            let (alpha, beta) = (ri(rng(&mut st, -5, 5), 1), ri(rng(&mut st, -5, 5), 1));
            let eps = match rng(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng(&mut st, 5_000, 200_000)),
                _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            // n2 = α·n0 + β·n1 + ε·(n0×n1) — near-coplanar with n0, n1.
            let n2 = add3(
                add3(smul(alpha, n0), smul(beta, n1)),
                smul(eps, rat_cross(n0, n1)),
            );
            let ax = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let mk =
                |n: [Rat; 3], st: &mut u64| rot_plane(plane_norm(n, rand_base(st)), ax, ang, piv);
            let pa = mk(n0, &mut st);
            let pb = mk(n1, &mut st);
            let pc = mk(n2, &mut st);
            if !(well_conditioned(&pa) && well_conditioned(&pb) && well_conditioned(&pc)) {
                continue; // a degenerate plane is not the regime under test
            }
            let judged = dir_sign_judge(t3(&pa), t3(&pb), t3(&pc));
            let (d, _) = cramer_iv([
                plane_iv(&pa[0], &pa[1], &pa[2]),
                plane_iv(&pb[0], &pb[1], &pb[2]),
                plane_iv(&pc[0], &pc[1], &pc[2]),
            ]);
            if d.sign().is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match dir_sign_truth(t3(&pa), t3(&pb), t3(&pc), GT) {
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
        }
        eprintln!(
            "[H-i] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(wrong, 0, "dir_sign must never disagree with GT (soundness)");
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path"
        );
    }

    /// The `dir_orient3d` determinant `det[d, x−base, y−base]` at `prec` bits, floored to
    /// declare-0 — the GT / escalation realization (mirrors the judge's hp path).
    fn dir_orient_at(d: [Rat; 3], base: &Pt3, x: &Pt3, y: &Pt3, prec: usize) -> Option<bool> {
        let dp = Pt3::at(d);
        let sub = |u: &BigFloat, v: &BigFloat| u.sub(v, prec, HP_RM);
        let (bh, xh, yh, dh) = (
            base.hp_coord(prec),
            x.hp_coord(prec),
            y.hp_coord(prec),
            dp.hp_coord(prec),
        );
        let det = det3_big(
            [
                [dh[0].clone(), dh[1].clone(), dh[2].clone()],
                [
                    sub(&xh[0], &bh[0]),
                    sub(&xh[1], &bh[1]),
                    sub(&xh[2], &bh[2]),
                ],
                [
                    sub(&yh[0], &bh[0]),
                    sub(&yh[1], &bh[1]),
                    sub(&yh[2], &bh[2]),
                ],
            ],
            prec,
        );
        let subf = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
        let mag = det3_mag([
            dp.coord,
            subf(x.coord, base.coord),
            subf(y.coord, base.coord),
        ]);
        sign_with_floor(&det, mag, prec)
    }

    /// GT-stable truth: `Some` only when `prec` and `prec + 128` agree (else too degenerate
    /// for the ground truth itself — skip, not a spurious "wrong"). Mirrors `dir_sign_truth`.
    fn dir_orient_truth(d: [Rat; 3], base: &Pt3, x: &Pt3, y: &Pt3, prec: usize) -> Option<Orient> {
        match (
            dir_orient_at(d, base, x, y, prec),
            dir_orient_at(d, base, x, y, prec + 128),
        ) {
            (Some(a), Some(b)) if a == b => Some(orient_of(a)),
            _ => None,
        }
    }

    /// Sanity: `det[d, x−base, y−base] = d·((x−base)×(y−base))`. With `base=0, x=e_x, y=e_y`
    /// the edge normal is `+e_z`, so `d=+e_z → Positive`, `−e_z → Negative`, in-plane →
    /// `Zero`; swapping the edge pair flips the sign.
    #[test]
    fn dir_orient3d_judge_sanity() {
        let o = Pt3::at([ri(0, 1); 3]);
        let x = Pt3::at([ri(1, 1), ri(0, 1), ri(0, 1)]);
        let y = Pt3::at([ri(0, 1), ri(1, 1), ri(0, 1)]);
        let e = |a: i128, b: i128, c: i128| [ri(a, 1), ri(b, 1), ri(c, 1)];
        assert_eq!(dir_orient3d_judge(e(0, 0, 1), &o, &x, &y), Orient::Positive);
        assert_eq!(
            dir_orient3d_judge(e(0, 0, -1), &o, &x, &y),
            Orient::Negative
        );
        assert_eq!(
            dir_orient3d_judge(e(1, 0, 0), &o, &x, &y),
            Orient::Zero,
            "d in the edge plane → det 0"
        );
        assert_eq!(
            dir_orient3d_judge(e(0, 0, 1), &o, &y, &x),
            Orient::Negative,
            "swapping x,y flips the sign"
        );
    }

    /// `orient3d_ray(base, dir, x, y)` is exactly `orient3d(base, base+dir, x, y)`: on
    /// unrotated points (where `base+dir` *is* an exact `Pt3`) it must equal `orient3d_judge`
    /// with the ideal point materialized — the argument-for-argument reduction the ops
    /// ray-triangle caller relies on.
    #[test]
    fn orient3d_ray_matches_materialized() {
        let mut st = 0x0D1B_7A44_0F0F_9001u64;
        let at = |b: [Rat; 3]| Pt3::at(b);
        for _ in 0..200 {
            let (bb, xb, yb) = (rand_base(&mut st), rand_base(&mut st), rand_base(&mut st));
            let dir = [
                ri(rng(&mut st, -9, 9), 1),
                ri(rng(&mut st, -9, 9), 1),
                ri(rng(&mut st, -9, 9), 1),
            ];
            let q = [
                bb[0].checked_add(dir[0]).unwrap(),
                bb[1].checked_add(dir[1]).unwrap(),
                bb[2].checked_add(dir[2]).unwrap(),
            ];
            let (base, x, y) = (at(bb), at(xb), at(yb));
            assert_eq!(
                orient3d_ray(&base, dir, &x, &y),
                orient3d_judge(&base, &at(q), &x, &y),
                "orient3d_ray == orient3d(base, base+dir, x, y)"
            );
        }
    }

    /// Soundness — `dir_orient3d` over a rotated corpus, oracle = a GT-stable 512-bit
    /// determinant (NOT rotation-invariance: `d` is fixed while the points rotate, so the
    /// sign is not preserved). Two regimes: random directions, and **near-grazing** (`d`
    /// almost in the edge plane, `det ≈ 0`) to force escalation. The judge must never claim
    /// a sign opposite to the GT; both the fast filter and astro-float paths are exercised.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_dir_orient3d_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut st = 0x0D1A_5170_BEEF_4242u64;
        for _ in 0..1500 {
            let (bb, xb, yb) = (rand_base(&mut st), rand_base(&mut st), rand_base(&mut st));
            let grazing = rng(&mut st, 0, 1) == 1;
            let ax = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let rp = |b: [Rat; 3]| Pt3::at(b).rotate_about(ax, ang, piv);
            let (base, x, y) = (rp(bb), rp(xb), rp(yb));
            let tri = [base.clone(), x.clone(), y.clone()];
            if !well_conditioned(&tri) {
                continue; // a degenerate edge fan is not the regime under test
            }
            // Two direction regimes. Random: any integer `d`. near-grazing: a large integer
            // multiple of the **rotated** in-plane edge `x−base` (built from the f64 coords),
            // so `d` is nearly ⊥ the rotated normal → `det ≈ 0` → the fast filter fails and
            // the astro-float path decides. Built after rotation, since the rotated normal is
            // what `d` must graze.
            let d = if !grazing {
                [
                    ri(rng(&mut st, -19, 19), 1),
                    ri(rng(&mut st, -19, 19), 1),
                    ri(rng(&mut st, -19, 19), 1),
                ]
            } else {
                let big = rng(&mut st, 100_000_000, 9_000_000_000) as f64;
                let comp = |k: usize| ri((big * (x.coord[k] - base.coord[k])).round() as i128, 1);
                [comp(0), comp(1), comp(2)]
            };
            let judged = dir_orient3d_judge(d, &base, &x, &y);
            // Recompute the Iv filter to tally which path resolved (mirrors the judge).
            let dp = Pt3::at(d);
            let (bi, xi, yi, di) = (pt_iv(&base), pt_iv(&x), pt_iv(&y), pt_iv(&dp));
            let subi = |u: [Iv; 3], v: [Iv; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
            if det3_iv([di, subi(xi, bi), subi(yi, bi)]).sign().is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match dir_orient_truth(d, &base, &x, &y, GT) {
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
        }
        eprintln!(
            "[dir_orient3d] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "dir_orient3d must never disagree with GT (soundness)"
        );
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path"
        );
    }
}
