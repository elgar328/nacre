use super::*;
/// **Everything [`orient3d_judge`] can answer without spending unbounded precision** — the f64
/// filter, then the shared-motion exact path. `None` means "not settled here".
///
/// **It never escalates, and that is the point.** A caller that only wants to *skip work when it
/// can prove the work is pointless* must not pay an astro-float climb to find out; and because
/// `None` sends it back to doing the work, a missed proof costs time and nothing else. So this
/// can only ever lose an optimization — never change an answer.
///
/// Split out of [`orient3d_judge`] rather than written beside it: two implementations of one filter
/// would be two error budgets to keep in step — two answering sites for one question, which drift.
/// The judge calls this, so they cannot drift, and soundness here is not a property to test but a
/// consequence of being the same code.
pub fn orient3d_filter(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
) -> Option<Orient> {
    match orient3d_iv([pa, pb, pc, pd]).sign() {
        Some(true) => return Some(Orient::Positive),
        Some(false) => return Some(Orient::Negative),
        None => {}
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
    j: ProofStandard,
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
    let (b, c, d) = (pb.coord_hp(prec), pc.coord_hp(prec), pd.coord_hp(prec));
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
