use super::*;
/// **Carry a realized point through `nodes`, continuing from wherever it already is.**
///
/// [`WitnessPoint::compute_hp`] is this function seeded from the definition's base, and that is the
/// whole relationship: the fold is a left fold whose entire state is the three coordinates, so the
/// value after `k` nodes *is* the prefix, and resuming from it reaches the same bits as folding the
/// prefix again. That is what lets a caller who kept the previous generation's value pay for one
/// node instead of the whole history.
///
/// ⚠ **One fold, not two.** The per-node arithmetic lives in [`fold_one`] and nothing else spells
/// it, because two copies of a rounding rule are two rules that drift.
pub fn fold_suffix(prefix: [HpBounded; 3], nodes: &[MoveNode], prec: usize) -> [HpBounded; 3] {
    nodes.iter().fold(prefix, |p, n| fold_one(p, n, prec))
}

/// One node applied to a realized point — the body of the fold, spelled once.
///
/// ⚠ The degenerate arms **return `p` untouched**. They read as `continue` when this was a loop
/// body, and "skip this node" is what they meant; extracted, the same meaning is "this node moves
/// nothing". All three are unreachable (their producers refuse a non-positive squared length), and
/// they stay total rather than inventing a bound.
fn fold_one(mut p: [HpBounded; 3], node: &MoveNode, prec: usize) -> [HpBounded; 3] {
    match node {
        MoveNode::Rotate { axis, angle, pivot } => {
            let (i, j) = axis.plane();
            let (c, s) = angle.cos_sin_bounded(prec);
            let (px, py) = (
                HpBounded::of_rat(pivot[i], prec),
                HpBounded::of_rat(pivot[j], prec),
            );
            // pivot-relative: u = p − pivot, rotate, shift back.
            let u = p[i].sub(&px, prec);
            let v = p[j].sub(&py, prec);
            p[i] = px.add(&u.mul(&c, prec).sub(&v.mul(&s, prec), prec), prec);
            p[j] = py.add(&u.mul(&s, prec).add(&v.mul(&c, prec), prec), prec);
        }
        // A translation is exact input: the only error is `rat_to_hp`'s own division
        // rounding and the add's half-ulp, both of which the interval carries.
        MoveNode::Translate { offset } => {
            for k in 0..3 {
                p[k] = p[k].add(&HpBounded::of_rat(offset[k], prec), prec);
            }
        }
        // `2·offset − x` on one coordinate: exact input, so the interval carries only its
        // own rounding. The other two coordinates are untouched.
        MoveNode::Mirror { axis, offset } => {
            let k = axis.index();
            let c = HpBounded::of_rat(*offset, prec);
            p[k] = c.add(&c, prec).sub(&p[k], prec);
        }
        // The same derivation as [`WitnessPoint::frame`], in the domain that carries its own
        // error: two `1/√` intervals, two scaled basis vectors, their cross product, and
        // the combination. Nothing is charged by hand — every rounding is inside an
        // `HpBounded`, which is the point of realizing here rather than trusting the f64 tol.
        //
        // A degenerate frame cannot arise here: `WitnessPoint::frame` is the only producer of this
        // node and it refuses a non-positive squared length, so the chain never holds one.
        MoveNode::Frame { frame } => {
            let inv = |v: Rat| nacre_scalar::inv_sqrt_bounded(v, prec);
            let (Some(iu), Some(iw)) = (inv(frame.uu), inv(frame.nn)) else {
                return p; // unreachable — `plane_frame` refuses a non-positive length
            };
            let scaled = |v: [Rat; 3], s: &HpBounded| {
                [0, 1, 2].map(|k| HpBounded::of_rat(v[k], prec).mul(s, prec))
            };
            let (uh, wh) = (scaled(frame.u_raw, &iu), scaled(frame.n, &iw));
            // The same two routes `WitnessPoint::frame` takes, in the domain that carries its own
            // error interval.
            let vh = match frame.v.and_then(|(v_raw, vv)| Some((v_raw, inv(vv)?))) {
                Some((v_raw, iv)) => scaled(v_raw, &iv),
                None => [0, 1, 2].map(|k| {
                    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                    wh[i].mul(&uh[j], prec).sub(&wh[j].mul(&uh[i], prec), prec)
                }),
            };
            p = [0, 1, 2].map(|k| {
                HpBounded::of_rat(frame.origin[k], prec)
                    .add(&p[0].mul(&uh[k], prec), prec)
                    .add(&p[1].mul(&vh[k], prec), prec)
                    .add(&p[2].mul(&wh[k], prec), prec)
            });
        }
        // [`MoveNode::Frame`]'s wide twin — same shape, arbitrary-precision inputs.
        // `v_raw`/`vv` always exist (nothing overflows a `BigInt`), so there is no
        // cross-product fallback branch here.
        MoveNode::FrameWide(f) => {
            let inv = |v: &num_bigint::BigInt| nacre_scalar::inv_sqrt_bigint_bounded(v, prec);
            let (Some(iu), Some(iv2), Some(iw)) = (inv(&f.uu), inv(&f.vv), inv(&f.nn)) else {
                return p; // unreachable — the builder derives positive squared lengths
            };
            let scaled = |v: &[num_bigint::BigInt; 3], s: &HpBounded| {
                [0, 1, 2].map(|k| HpBounded::of_bigint(&v[k], prec).mul(s, prec))
            };
            let (uh, vh, wh) = (
                scaled(&f.u_raw, &iu),
                scaled(&f.v_raw, &iv2),
                scaled(&f.n, &iw),
            );
            let od = HpBounded::of_bigint(&f.origin_den, prec);
            let origin = (|| {
                let o =
                    |k: usize| HpBounded::of_bigint(&f.origin_num[k], prec).div_exact(&od, prec);
                Some([o(0)?, o(1)?, o(2)?])
            })();
            let Some(o) = origin else {
                return p; // unreachable — the builder's denominator is n·n > 0
            };
            p = [0, 1, 2].map(|k| {
                o[k].add(&p[0].mul(&uh[k], prec), prec)
                    .add(&p[1].mul(&vh[k], prec), prec)
                    .add(&p[2].mul(&wh[k], prec), prec)
            });
        }
        // The judged frame — the same derivation the f64 rung took, at this precision.
        MoveNode::FrameThrough(f) => {
            // Construction proved the derivation at the fixed rung, and every input
            // radius shrinks as `prec` grows (they are the points' own realization
            // errors), so the primary call is expected to succeed; the rung fallback
            // keeps the arm total without inventing a bound if that expectation is ever
            // wrong — a 128-bit basis is a *sound*, merely wider, interval for the same
            // true frame.
            let Some(basis) = judged_basis(f, prec).or_else(|| judged_basis(f, FrameThrough::RUNG))
            else {
                return p; // unreachable — `FrameThrough::of` proved the rung derivation
            };
            let [o, u, v, w] = basis;
            p = [0, 1, 2].map(|k| {
                o[k].add(&p[0].mul(&u[k], prec), prec)
                    .add(&p[1].mul(&v[k], prec), prec)
                    .add(&p[2].mul(&w[k], prec), prec)
            });
        }
    }
    p
}
