use super::*;
/// **Everything [`orient3d_judge`] can answer without spending unbounded precision** — the f64
/// filter, then the shared-motion exact path. `None` means "not settled here".
///
/// **It never escalates, and that is the point.** A caller that only wants to *skip work when it
/// can prove the work is pointless* must not pay an astro-float climb to find out; and because
/// `None` sends it back to doing the work, a missed proof costs time and nothing else. So this
/// can only ever lose an optimization — never change an answer.
///
/// Split out of [`orient3d_judge`] rather than written beside it: two implementations of one
/// filter would be two error budgets to keep in step, and this kernel has already been bitten
/// twice by a question with two answering sites. The judge calls this, so they cannot drift, and
/// soundness here is not a property to test but a consequence of being the same code.
pub fn orient3d_filter(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
) -> Option<Orient> {
    let (a, b, c, d) = (pa.coord(), pb.coord(), pc.coord(), pd.coord());
    let det = det3_f64(a, b, c, d);
    let bound = det3_bound([a, b, c, d], [pa.tol(), pb.tol(), pc.tol(), pd.tol()]);
    if det > bound {
        return Some(Orient::Positive);
    }
    if det < -bound {
        return Some(Orient::Negative);
    }
    // One shared rigid motion cancels out of `det[a−d, b−d, c−d]` (`det(R) = 1`), so the same
    // question is answered exactly on the pre-rotation coordinates — no tolerance, no escalation.
    // Its zero is a **proved** zero, which is why it is a sign and not a coincidence.
    let base = shared_base(&[pa, pb, pc, pd])?;
    Some(
        match nacre_predicates::orient3d(base[0], base[1], base[2], base[3]) {
            x if x > 0.0 => Orient::Positive,
            x if x < 0.0 => Orient::Negative,
            _ => Orient::Zero,
        },
    )
}

/// CIP `orient3d`: f64 filter (`|det| > bound` → trust the sign), else escalate to
/// astro-float at the operation's precision; a determinant that still straddles zero there is
/// turned into a distance and answered by [`escalate`]. Path-independent (a
/// function of the four point definitions). This is the *tol > 0* path — a tol-0
/// config is exact/faster via `nacre-predicates` (Shewchuk), routed above this crate.
pub fn orient3d_judge(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
    j: Standard,
) -> Decision {
    if let Some(o) = orient3d_filter(pa, pb, pc, pd) {
        return Decision::Sign(o);
    }
    escalate(j, j.coincidence, |prec| {
        let det = det3_hp(pa, pb, pc, pd, prec);
        match det.sign() {
            Some(pos) => Ok(orient_of(pos)),
            // Only now is the area term worth forming: the sign decides on the first try in all
            // but the coincident cases, and the cross product is six high-precision products.
            None => Err(distance_bound(det.error, &cross_of(pb, pc, pd, prec), prec)),
        }
    })
}

/// `(b−d) × (c−d)` at `prec` bits — the area term that turns [`orient3d_judge`]'s determinant
/// into a height above the plane through `b, c, d`.
fn cross_of(
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
    prec: usize,
) -> [HpBounded; 3] {
    let (b, c, d) = (pb.hp_coord(prec), pc.hp_coord(prec), pd.hp_coord(prec));
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let (u, v) = (
        [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
        [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
    );
    [
        u[1].mul(&v[2], prec).sub(&u[2].mul(&v[1], prec), prec),
        u[2].mul(&v[0], prec).sub(&u[0].mul(&v[2], prec), prec),
        u[0].mul(&v[1], prec).sub(&u[1].mul(&v[0], prec), prec),
    ]
}

/// CIP `dir_orient3d`: the sign of `det[d, x−base, y−base] = d·((x−base)×(y−base))` — the
/// orientation of the ray direction `d` against the edge fan `(base→x, base→y)`. The
/// **direction analogue** of [`orient3d_judge`]: `point_in_solid`'s ray-triangle test asks
/// `orient3d(p, p+d, ·, ·)`, but `p+d` (a rotated point plus a rational offset) has no exact
/// `base+chain` `WitnessPoint` (`R⁻¹d` is irrational). Every such determinant reduces to this form,
/// where `d` enters as one **exact** (error-0) column and only `x, y, base` carry rotation tol.
/// Interval filter → astro-float escalation, exactly like the indirect judges; a below-floor
/// determinant that stays undecided is [`Orient::Zero`], absorbed by the caller's ray retry.
///
/// **The one judge with no metric normalization**, and the reason is that its caller does not
/// need one: a zero here means the ray runs along the face's plane, and `boolean` answers that by
/// casting a different ray rather than by asking how close it came.
///
/// `d` is realized through an unrotated [`WitnessPoint`] purely to reuse `pt_iv`/`hp_coord` for its
/// coord/tol/hp — the row is `d` itself, never `d − base`.
pub fn dir_orient3d_judge(
    d: [Rat; 3],
    base: &WitnessPoint,
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Orient {
    let dp = WitnessPoint::at(d);
    let (bi, xi, yi, di) = (base.realized, x.realized, y.realized, dp.realized);
    let sub_iv =
        |u: [Bounded; 3], v: [Bounded; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
    if let Some(pos) = det3_iv([di, sub_iv(xi, bi), sub_iv(yi, bi)]).sign() {
        return orient_of(pos);
    }
    // Escalate: the same determinant at prec from the exact definitions.
    let (bh, xh, yh, dh) = (
        base.hp_coord(prec),
        x.hp_coord(prec),
        y.hp_coord(prec),
        dp.hp_coord(prec),
    );
    let sub_hp = |u: &HpBounded, v: &HpBounded| u.sub(v, prec);
    let rows = [
        dh,
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
    match det3_big_rows(&rows, prec).sign() {
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
pub fn orient3d_ray(
    base: &WitnessPoint,
    dir: [Rat; 3],
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Orient {
    dir_orient3d_judge(dir, y, x, base, prec)
}
