//! The 3D toleranced point + its `orient3d` judgment (overhaul §TIP, stage 2). The
//! 3D analogue of [`crate::frame`] (`Pt2` / `orient2d_judge`).
//!
//! A rotated point cannot be held exactly (cos/sin are irrational), but its f64
//! realization carries a **direction-wise xyz tol** (§TIP ⑤) that soundly bounds
//! the error, accumulated as the definition is turned through a chain of
//! axis-aligned rotations (§TIP ①: `new tol = |R|·old + mix`). [`orient3d_judge`]
//! consumes that tol: an f64 determinant filter with a sound error bound
//! ([`det3_bound`], §TIP ②) decides the easy cases, the ambiguous ones **escalate**
//! to astro-float from the point definitions, and a determinant below the precision
//! floor is a **declare-0** ([`Orient::Zero`]). This is the judgment **layer only**
//! ("층만") — not yet wired into boolean (that is stage 3), and it is the *tol > 0*
//! path: a tol-0 (`Constructed`) config is faster/exact via `nacre-predicates`
//! (Shewchuk), routed by a higher layer, not here.
//!
//! Validated in `experiments/exact3d` (FINDINGS = GO): H-a (`det3_bound` soundness),
//! H-d/H-f (chain + arbitrary-pivot tol) — the bound never under-estimates the true
//! error (astro-float ground truth) over random heterogeneous-rotation configs.

use crate::frame::{DA_F64, JUDGE_PREC, bf_mag, rat_to_big};
use crate::{Angle, Axis, HP_RM, Orient, Rat};
use astro_float::BigFloat;

/// One rotation in a point's definition: turn about `axis` (the line through the
/// rational pivot `point`) by the rational `angle`. `point = [0,0,0]` is the
/// origin-pivot case; the kernel's `Rotation` carries an arbitrary rational pivot.
#[derive(Clone, Copy, Debug)]
pub struct RotNode {
    pub axis: Axis,
    pub angle: Angle,
    pub point: [Rat; 3],
}

/// A rational base point carried through a chain of axis rotations (§TIP ⑦ rotation
/// history). `base` + `chain` are the exact **definition** (never lost); `coord` is
/// the f64 realization (a cache), and `tol` bounds its error as a **direction-wise
/// xyz vector** (§TIP ⑤). [`hp_coord`](Self::hp_coord) realizes the chain at
/// arbitrary precision from the definition, so two points with the same definition
/// realize identically (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct Pt3 {
    pub base: [Rat; 3],
    pub chain: Vec<RotNode>,
    pub coord: [f64; 3],
    pub tol: [f64; 3],
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
        self
    }

    /// The coordinate realized at `prec` bits from the **definition** (base rotated
    /// through the chain, each node about its pivot) — path-independent ground truth /
    /// escalation realization.
    pub fn hp_coord(&self, prec: usize) -> [BigFloat; 3] {
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
/// tol (§TIP ②). The determinant is six signed triple-products of the edge entries;
/// each entry `(a−d)[k]` carries tol `tol_a[k] + tol_d[k]`. The bound sums the six
/// product radii (input-tol propagation, triangle-inequality worst case) plus a term
/// for the f64 rounding of the determinant's own arithmetic. The 3D analogue of
/// [`crate::frame`]'s 2D `det_bound`. Validated in exact3d (H-a).
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
    // Sign via is_positive (astro-float#44 workaround — see crate::frame).
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
}
