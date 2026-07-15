//! The 3D toleranced point — direction-wise tol propagation through a rotation
//! chain (overhaul §TIP, stage 2). The 3D analogue of [`crate::frame::Pt2`].
//!
//! A rotated point cannot be held exactly (cos/sin are irrational), but its f64
//! realization carries a **direction-wise xyz tol** (§TIP ⑤) that soundly bounds
//! the error, accumulated as the definition is turned through a chain of
//! axis-aligned rotations (§TIP ①: `new tol = |R|·old + mix`). The judgment
//! predicate that consumes this tol (`orient3d_judge`) is a later cell; this module
//! is the tol computation only ("층만" — no boolean wiring).
//!
//! Validated in `experiments/exact3d` (FINDINGS = GO, incl. H-f arbitrary pivots):
//! the tol soundly bounds the true error (astro-float ground truth) over random
//! rotation chains with arbitrary rational pivots, exact and inexact angles.

use crate::frame::{DA_F64, bf_mag, rat_to_big};
use crate::{Angle, Axis, HP_RM, Rat};
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
    #[test]
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
}
