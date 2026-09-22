use super::*;
/// The three edge rows of `orient3d(a,b,c,d) = det[a−d, b−d, c−d]`.
fn rows(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] - d[0], a[1] - d[1], a[2] - d[2]],
        [b[0] - d[0], b[1] - d[1], b[2] - d[2]],
        [c[0] - d[0], c[1] - d[1], c[2] - d[2]],
    ]
}

/// f64 `orient3d` determinant `det[a−d, b−d, c−d]`.
pub(super) fn det3_f64(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
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
/// tol. The determinant is six signed triple-products of the edge entries;
/// each entry `(a−d)[k]` carries tol `tol_a[k] + tol_d[k]`. The bound sums the six
/// product radii (input-tol propagation, triangle-inequality worst case) plus a term
/// for the f64 rounding of the determinant's own arithmetic. Checked by
/// `det3_bound_soundness`.
pub(super) fn det3_bound(p: [[f64; 3]; 4], t: [[f64; 3]; 4]) -> f64 {
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
    // f64 rounding of the determinant's own ~17 operations. Round-to-nearest costs `ε/2` each, so
    // the accumulation is about `8.5·ε·mag`; `16·ε` is that with a factor of two over it.
    input_tol + 16.0 * f64::EPSILON * mag
}

/// 3×3 determinant of high-precision interval rows at `prec` bits — the same expression as
/// [`det3_iv`], one precision up, so the filter and the escalation cannot drift apart.
///
/// The entries come in **by reference**, and that is not a micro-optimization: Cramer's four
/// matrices are the same nine or twelve values in different arrangements, so a signature that
/// owns its rows makes the caller copy each value up to three times over. An `HpBounded` copy is a
/// mantissa heap allocation, and this determinant is the hot path — [`cramer_hp`] alone made
/// forty-five of them per call for values it only ever read.
pub(super) fn det3_big(r: [[&HpBounded; 3]; 3], prec: usize) -> HpBounded {
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let m0 = mul(r[1][1], r[2][2]).sub(&mul(r[1][2], r[2][1]), prec);
    let m1 = mul(r[1][0], r[2][2]).sub(&mul(r[1][2], r[2][0]), prec);
    let m2 = mul(r[1][0], r[2][1]).sub(&mul(r[1][1], r[2][0]), prec);
    mul(r[0][0], &m0)
        .sub(&mul(r[0][1], &m1), prec)
        .add(&mul(r[0][2], &m2), prec)
}

/// [`det3_big`] over rows the caller already owns — the borrow, spelled once.
pub(super) fn det3_big_rows(r: &[[HpBounded; 3]; 3], prec: usize) -> HpBounded {
    det3_big(
        [
            [&r[0][0], &r[0][1], &r[0][2]],
            [&r[1][0], &r[1][1], &r[1][2]],
            [&r[2][0], &r[2][1], &r[2][2]],
        ],
        prec,
    )
}

/// `orient3d` determinant realized at `prec` bits (astro-float) from the point
/// definitions — path-independent ground truth / escalation realization.
pub(super) fn det3_hp(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
    prec: usize,
) -> HpBounded {
    let (a, b, c, d) = (
        pa.hp_coord(prec),
        pb.hp_coord(prec),
        pc.hp_coord(prec),
        pd.hp_coord(prec),
    );
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let rows = [
        [sub(&a[0], &d[0]), sub(&a[1], &d[1]), sub(&a[2], &d[2])],
        [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
        [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
    ];
    det3_big_rows(&rows, prec)
}

/// The points' pre-motion coordinates, **brought into the moved frame's handedness**, when they
/// all carry one and the same motion and every base is `f64`-representable.
///
/// A motion preserves the determinants these judges take *up to its own determinant*, so under
/// that condition the answer on the bases is the answer — exactly. Chains are compared
/// structurally, so two spellings of one motion read as different: a missed cancellation is
/// slower, never wrong. Pivots must match too, since a rotation about `c` is `Rx + (c − Rc)` and
/// two pivots leave two translations that do not cancel in a difference.
///
/// **A reflection is improper (`det = −1`), and the correction belongs here rather than at each
/// caller.** Every point carries the same chain, so an odd number of reflections flips the
/// determinant *uniformly* — the shortcut would return a confidently wrong sign, not a
/// conservative miss. Rather than hand each judge a parity to multiply in (four sites, each
/// taking a different quantity, one of them already carrying a second sign convention), the bases
/// are handed back **already reflected once** when the parity is odd: the base frame then has the
/// moved frame's handedness and every determinant question transfers unchanged.
///
/// The canonical reflection is a sign flip on x, which is exact for every finite `f64`.
pub(super) fn shared_base<const N: usize>(pts: &[&WitnessPoint; N]) -> Option<[[f64; 3]; N]> {
    let first = pts[0];
    if first.chain.is_empty() {
        return None; // nothing to cancel; the caller's own fast path already handled this
    }
    for p in pts {
        if p.chain.len() != first.chain.len() {
            return None;
        }
        // Structural equality over the **whole** node, variant included: two motions that compare
        // equal here are declared one motion and the exact predicate then answers in their shared
        // pre-motion frame. A comparison that ignored the variant would answer a different
        // question with full confidence — silently wrong, not slow.
        if p.chain.iter().zip(first.chain.iter()).any(|(a, b)| a != b) {
            return None;
        }
    }
    let exact = |r: Rat| Rat::try_from_f64(r.to_f64()) == Some(r);
    if !pts.iter().all(|p| p.base.iter().all(|&r| exact(r))) {
        return None; // the exact predicate takes f64; a base that does not round-trip cannot go
    }
    // Odd parity ⇒ hand back the mirror image, so the base frame is handed like the moved one.
    let flip = if chain_parity(&first.chain) < 0 {
        -1.0
    } else {
        1.0
    };
    Some(std::array::from_fn(|i| {
        [
            flip * pts[i].base[0].to_f64(),
            pts[i].base[1].to_f64(),
            pts[i].base[2].to_f64(),
        ]
    }))
}

/// `−1` when a chain contains an odd number of reflections, `+1` otherwise — the factor an
/// improper motion puts on any determinant of its images.
pub fn chain_parity(chain: &[MoveNode]) -> i8 {
    let mirrors = chain
        .iter()
        .filter(|n| matches!(n, MoveNode::Mirror { .. }))
        .count();
    if mirrors % 2 == 0 { 1 } else { -1 }
}
