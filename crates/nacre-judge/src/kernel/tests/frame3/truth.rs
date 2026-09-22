//! **Soundness oracles that do not borrow the judge's arithmetic.**
//!
//! A judge that escalates realizes its question at high precision through its own functions
//! (`plane_hp`, `cramer_hp`, `indirect_hp`, `cmp_hp_with_gap`, `det3_big_rows`). An oracle built
//! from those same functions can only compare a function with itself at another precision: a
//! mistake in them is made twice and passes. So the truth here shares only the **input** — each
//! point's own realization, `WitnessPoint::hp_coord` (cross-checked against the separately
//! written f64 walk by `tol_bounds_error_over_random_chains`) — and takes a different road from
//! there:
//!
//! * a plane's coefficients are formed here from its three points;
//! * a three-plane point `V` is **solved for** by Gaussian elimination with partial pivoting and
//!   division — the judge never materializes `V` (it reads `sign(D)·sign(M)` off Cramer's rule);
//! * a determinant's sign is the product of that elimination's pivots and its row-swap parity.
//!
//! Values are raw `BigFloat`s, no radius — confidence comes from two precisions (`prec` and
//! `2·prec`) agreeing on a nonzero sign, not from any bound the judge also ships. A question the
//! two precisions do not agree on is skipped, not scored.
//!
//! ★★ **Agreement between two precisions cannot tell a zero.** Where the true value is exactly
//! zero, each precision leaves a rounding residue, and two residues share a sign about half the
//! time — measured, this oracle and the one it replaced each read a confident sign on roughly half
//! of the constructed zeros (1,108 of 2,259 for `indirect_orient3d`), and on the other half they
//! disagreed with each other. So a zero is never asked of an oracle: the corpora **construct** it
//! (a pushed point with `ε = 0`, the meet placed on the query point, a normal in the span of the
//! other two, a tie in an axis a rotation keeps) and the soundness tests assert the judge answers
//! `Orient::Zero` there. That is the check that catches a confident sign on a true zero — the
//! shape of the cancellation bug `indirect_hp`'s doc records. Planted in the judge's
//! high-precision path, that defect passed the tests this replaced (`indirect_orient3d`,
//! `indirect_cmp_coord`) and fails these.
//!
//! What was *not* shown: that this oracle catches a defect on a **nonzero** answer the old one
//! missed. The judge jumps straight to the precision its condition number asks for, so a nonzero
//! case is decided with a margin no radius-sized planted error (even 256× the radius) reached — no
//! plant changed a nonzero answer. A formula error in the high-precision path shows anyway, in the
//! filter-resolved cases, where the judge does not use that path. The oracle is kept for what it
//! is: an answer that does not borrow what it checks.

use super::*;

type P3 = [BigFloat; 3];

fn add(a: &BigFloat, b: &BigFloat, p: usize) -> BigFloat {
    a.add(b, p, HP_RM)
}
fn sub(a: &BigFloat, b: &BigFloat, p: usize) -> BigFloat {
    a.sub(b, p, HP_RM)
}
fn mul(a: &BigFloat, b: &BigFloat, p: usize) -> BigFloat {
    a.mul(b, p, HP_RM)
}

/// `Some(true)` positive, `Some(false)` negative, `None` exactly zero.
fn sign(x: &BigFloat) -> Option<bool> {
    (!x.is_zero()).then(|| x.is_positive())
}

fn point(w: &WitnessPoint, p: usize) -> P3 {
    let h = w.hp_coord(p);
    [h[0].value.clone(), h[1].value.clone(), h[2].value.clone()]
}

fn sub3(a: &P3, b: &P3, p: usize) -> P3 {
    [
        sub(&a[0], &b[0], p),
        sub(&a[1], &b[1], p),
        sub(&a[2], &b[2], p),
    ]
}

fn cross(a: &P3, b: &P3, p: usize) -> P3 {
    [
        sub(&mul(&a[1], &b[2], p), &mul(&a[2], &b[1], p), p),
        sub(&mul(&a[2], &b[0], p), &mul(&a[0], &b[2], p), p),
        sub(&mul(&a[0], &b[1], p), &mul(&a[1], &b[0], p), p),
    ]
}

fn dot(a: &P3, b: &P3, p: usize) -> BigFloat {
    add(
        &add(&mul(&a[0], &b[0], p), &mul(&a[1], &b[1], p), p),
        &mul(&a[2], &b[2], p),
        p,
    )
}

/// The plane through three points as `(n, rhs)` with `n·X = rhs`: `n = (p1−p0)×(p2−p0)` — the
/// orientation the judge's planes carry — and `rhs = n·p0`.
fn plane(t: (&WitnessPoint, &WitnessPoint, &WitnessPoint), p: usize) -> (P3, BigFloat) {
    let (a, b, c) = (point(t.0, p), point(t.1, p), point(t.2, p));
    let n = cross(&sub3(&b, &a, p), &sub3(&c, &a, p), p);
    let rhs = dot(&n, &a, p);
    (n, rhs)
}

/// Is `|a| > |b|`? Zeros are settled before comparing: astro-float's `cmp` answers "less" in
/// both directions against a zero (see `nacre-exact`'s `bounded` tests).
fn larger(a: &BigFloat, b: &BigFloat) -> bool {
    match (a.is_zero(), b.is_zero()) {
        (true, _) => false,
        (false, true) => true,
        (false, false) => a.abs().cmp(&b.abs()).is_some_and(|s| s > 0),
    }
}

/// Gaussian elimination with partial pivoting on `m` (with an optional right-hand side). Returns
/// the sign of `det m` — `None` when a pivot column is all zero — and, when a right-hand side was
/// given and the matrix is not singular, the solution.
fn eliminate(mut m: [P3; 3], mut rhs: Option<P3>, p: usize) -> (Option<bool>, Option<P3>) {
    let mut negative = false;
    for col in 0..3 {
        let mut best = col;
        for r in col + 1..3 {
            if larger(&m[r][col], &m[best][col]) {
                best = r;
            }
        }
        if m[best][col].is_zero() {
            return (None, None);
        }
        if best != col {
            m.swap(best, col);
            if let Some(v) = rhs.as_mut() {
                v.swap(best, col);
            }
            negative = !negative;
        }
        if m[col][col].is_negative() {
            negative = !negative;
        }
        let pivot_row = m[col].clone();
        for r in col + 1..3 {
            let f = m[r][col].div(&pivot_row[col], p, HP_RM);
            for (x, pv) in m[r].iter_mut().zip(&pivot_row).skip(col) {
                *x = sub(x, &mul(&f, pv, p), p);
            }
            if let Some(v) = rhs.as_mut() {
                let t = mul(&f, &v[col], p);
                v[r] = sub(&v[r], &t, p);
            }
        }
    }
    let solution = rhs.map(|v| {
        let zero = BigFloat::from_f64(0.0, p);
        let mut x = [zero.clone(), zero.clone(), zero];
        for i in (0..3).rev() {
            let mut acc = v[i].clone();
            for j in i + 1..3 {
                acc = sub(&acc, &mul(&m[i][j], &x[j], p), p);
            }
            x[i] = acc.div(&m[i][i], p, HP_RM);
        }
        x
    });
    (Some(!negative), solution)
}

/// The point where three planes meet, solved for.
fn meet(planes: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3], p: usize) -> Option<P3> {
    let [(na, ra), (nb, rb), (nc, rc)] = planes.map(|t| plane(t, p));
    eliminate([na, nb, nc], Some([ra, rb, rc]), p).1
}

/// A sign read at `prec` and `2·prec`; the truth only when both are the same nonzero sign.
fn stable(at: impl Fn(usize) -> Option<bool>, prec: usize) -> Option<Orient> {
    match (at(prec), at(2 * prec)) {
        (Some(a), Some(b)) if a == b => Some(orient_of(a)),
        _ => None,
    }
}

/// `orient3d(V, q, r, s)` for `V = ∩(a, b, c)`: the sign of `((q−s)×(r−s))·(V−s)` — the
/// convention `indirect_orient3d_judge` reports (`sign(D)·sign(M)`, `M = D·(V−s)·cross`).
pub(super) fn orient(
    a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    q: &WitnessPoint,
    r: &WitnessPoint,
    s: &WitnessPoint,
    prec: usize,
) -> Option<Orient> {
    stable(
        |p| {
            let v = meet([a, b, c], p)?;
            let (q, r, s) = (point(q, p), point(r, p), point(s, p));
            let n = cross(&sub3(&q, &s, p), &sub3(&r, &s, p), p);
            sign(&dot(&n, &sub3(&v, &s, p), p))
        },
        prec,
    )
}

/// The sign of `A[axis] − B[axis]` for the two three-plane points — `Positive` when `A` is
/// larger, as `indirect_cmp_coord_judge` reports.
pub(super) fn cmp(
    a: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    b: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    axis: usize,
    prec: usize,
) -> Option<Orient> {
    stable(
        |p| {
            let (va, vb) = (meet(a, p)?, meet(b, p)?);
            sign(&sub(&va[axis], &vb[axis], p))
        },
        prec,
    )
}

/// The sign of `det[n_a; n_b; n_c]`, the three planes' normals in order — `dir_sign_judge`'s `D`.
pub(super) fn dir_sign(
    a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    prec: usize,
) -> Option<Orient> {
    stable(
        |p| {
            let [na, nb, nc] = [a, b, c].map(|t| plane(t, p).0);
            eliminate([na, nb, nc], None, p).0
        },
        prec,
    )
}

/// The sign of `det[d; x−base; y−base]` — `dir_orient3d_judge`'s question.
pub(super) fn dir_orient(
    d: [Rat; 3],
    base: &WitnessPoint,
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Option<Orient> {
    stable(
        |p| {
            let dv = point(&WitnessPoint::at(d), p);
            let (b, x, y) = (point(base, p), point(x, p), point(y, p));
            eliminate([dv, sub3(&x, &b, p), sub3(&y, &b, p)], None, p).0
        },
        prec,
    )
}
