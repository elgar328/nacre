use super::*;
/// **The `orient3d` determinant of four witness points as an `f64` interval** — the filter's first
/// rung. It is [`det3_iv`] over the rows `a−d, b−d, c−d`, each a [`Bounded`] difference: the same
/// interval machine the indirect judges' Cramer determinants run on, so the radius comes from each
/// point's proven tol and every operation's rounding, never a hand-derived formula or a chosen
/// constant. A sign is trusted only when the interval excludes zero ([`Bounded::sign`]).
pub(super) fn orient3d_iv(p: [&WitnessPoint; 4]) -> Bounded {
    let c = p.map(|w| w.coord_iv);
    let row = |a: [Bounded; 3]| [0, 1, 2].map(|k| a[k].sub(c[3][k]));
    det3_iv([row(c[0]), row(c[1]), row(c[2])])
}

/// 3×3 determinant of high-precision interval rows at `prec` bits — the same expression as
/// [`det3_iv`], one precision up, so the filter and the escalation cannot drift apart.
///
/// The entries come in **by reference**, and that is not a micro-optimization: Cramer's four
/// matrices are the same nine or twelve values in different arrangements, so a signature that
/// owns its rows makes the caller copy each value up to three times over. An `HpBounded` copy is a
/// mantissa heap allocation, and this determinant is the hot path — [`cramer_hp`] alone made
/// forty-five of them per call for values it only ever read.
pub(super) fn det3_refs_hp(r: [[&HpBounded; 3]; 3], prec: usize) -> HpBounded {
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let m0 = mul(r[1][1], r[2][2]).sub(&mul(r[1][2], r[2][1]), prec);
    let m1 = mul(r[1][0], r[2][2]).sub(&mul(r[1][2], r[2][0]), prec);
    let m2 = mul(r[1][0], r[2][1]).sub(&mul(r[1][1], r[2][0]), prec);
    mul(r[0][0], &m0)
        .sub(&mul(r[0][1], &m1), prec)
        .add(&mul(r[0][2], &m2), prec)
}

/// [`det3_refs_hp`] over rows the caller already owns — the borrow, spelled once.
pub(super) fn det3_rows_hp(r: &[[HpBounded; 3]; 3], prec: usize) -> HpBounded {
    det3_refs_hp(
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
        pa.coord_hp(prec),
        pb.coord_hp(prec),
        pc.coord_hp(prec),
        pd.coord_hp(prec),
    );
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let rows = [
        [sub(&a[0], &d[0]), sub(&a[1], &d[1]), sub(&a[2], &d[2])],
        [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
        [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
    ];
    det3_rows_hp(&rows, prec)
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
    if !pts
        .iter()
        .all(|p| p.base.iter().all(|r| r.as_f64().is_some()))
    {
        return None; // the exact predicate takes f64; a base that is not an `f64` cannot go
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
