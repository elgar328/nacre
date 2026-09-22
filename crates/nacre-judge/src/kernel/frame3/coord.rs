use super::*;
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
fn cmp_filter(a: [[Bounded; 4]; 3], b: [[Bounded; 4]; 3], axis: usize) -> Option<Orient> {
    let (da, dva) = cramer_iv(a);
    let (db, dvb) = cramer_iv(b);
    let m = dva[axis].mul(db).sub(dvb[axis].mul(da));
    cmp_combine(m.sign(), da.sign(), db.sign())
}

/// The astro-float escalation for `cmp_coord`: the three signs at `prec` bits, and — when they do
/// not combine into an ordering — **how far apart the two coordinates may be** in the model's own
/// units.
///
/// Each implicit point's coordinate is `Dvec[axis]/D` (Cramer), so their difference is
/// `M / (D_a · D_b)` — a length, once divided. `M` alone is not: it carries both denominators, so
/// a threshold applied to it would move with the planes' scaling. The gap is `None` when either
/// denominator cannot be bounded away from zero, which is the near-parallel-planes case where
/// the implicit point itself is not well defined.
fn cmp_hp_with_gap(
    a: [[HpBounded; 4]; 3],
    b: [[HpBounded; 4]; 3],
    axis: usize,
    prec: usize,
) -> Result<Orient, Gap> {
    let (da, dva) = cramer_hp(&a, prec);
    let (db, dvb) = cramer_hp(&b, prec);
    let m = dva[axis]
        .mul(&db, prec)
        .sub(&dvb[axis].mul(&da, prec), prec);
    match cmp_combine(m.sign(), da.sign(), db.sign()) {
        Some(o) => Ok(o),
        None => Err(coord_gap(&m, &da, &db)),
    }
}

/// `|M| / |D_a · D_b|`, the separation `M` stands for — bounded above from `M`'s own radius, and
/// below by the denominators' lower bounds (dividing by an upper bound would understate the gap).
fn coord_gap(m: &HpBounded, da: &HpBounded, db: &HpBounded) -> Gap {
    let (lo_a, lo_b) = match (denom_lo(da), denom_lo(db)) {
        (Ok(a), Ok(b)) => (a, b),
        // Both may be short; the deeper shortfall is the one that has to be covered.
        (Err(Gap::Unresolved { short: x }), Err(Gap::Unresolved { short: y })) => {
            return Gap::Unresolved { short: x.max(y) };
        }
        (Err(g), _) | (_, Err(g)) => return g,
    };
    match m.error.over(lo_a.times(lo_b)) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}

/// CIP indirect `cmp_coord`: the sign of `a[axis] − b[axis]` where `a`, `b` are the
/// implicit points at which each three-plane triple meets — each plane through three
/// rotated points. Interval filter → astro-float escalation. `Positive` = `a[axis] >
/// b[axis]`, `Negative` = `<`, and a zero that is either **proved** or **proved within the
/// coincidence limit** (unlike the exact `nacre_predicates::indirect_cmp_coord`, whose `0` means
/// exactly equal — see [`Decision`] for which of the two this was). The
/// two-implicit companion of [`indirect_orient3d_judge`]; the boolean reaches it through
/// [`crate::predicate::Judge::cmp_coord`].
pub fn indirect_cmp_coord_judge(
    a: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    b: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    axis: usize,
    j: Standard,
) -> Decision {
    let iv = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
        [
            plane_iv(t[0].0, t[0].1, t[0].2),
            plane_iv(t[1].0, t[1].1, t[1].2),
            plane_iv(t[2].0, t[2].1, t[2].2),
        ]
    };
    if let Some(o) = cmp_filter(iv(a), iv(b), axis) {
        return Decision::Sign(o);
    }
    escalate(j, j.coincidence, |prec| {
        let hp = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
            [
                plane_hp(t[0].0, t[0].1, t[0].2, prec),
                plane_hp(t[1].0, t[1].1, t[1].2, prec),
                plane_hp(t[2].0, t[2].1, t[2].2, prec),
            ]
        };
        cmp_hp_with_gap(hp(a), hp(b), axis, prec)
    })
}

/// A definite `bool` sign (`true` = positive) as an [`Orient`].
pub(super) fn orient_of(pos: bool) -> Orient {
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
    a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    j: Standard,
) -> Decision {
    // ★ **Not fed from `Judge`'s interval-plane cache, and that is measured.** Wiring it here
    // bought 1.007x on the crossing collector: this judgement's cost is almost entirely the
    // escalation below, because a *true* zero — the walls meeting `wc` in no point, 8.9% of the
    // pairs the collector tests — can never be settled by an interval filter.
    //
    // ★★★ **And this is where the kernel's escalations come from.** Measured on `fin_fold(80)`:
    // this judgement is 5.4% of its own calls but **78% of every escalation in the boolean**, and
    // **92.6% of those escalations have a midpoint of exactly `0.0`** — genuine degeneracies the
    // arrangement builds by construction. That is why tightening a coordinate's tol does not buy
    // speed: setting the rotation's whole tol contribution to zero left the escalation count here
    // **bit-identical**. The lever, if one is ever wanted, is recognising a true zero *before*
    // climbing — not more accurate coordinates.
    let d = normals_det_iv([
        plane_iv(a.0, a.1, a.2),
        plane_iv(b.0, b.1, b.2),
        plane_iv(c.0, c.1, c.2),
    ]);
    if let Some(pos) = d.sign() {
        return Decision::Sign(orient_of(pos));
    }
    // **The one judgement whose limit is an angle.** Its question is about directions, so the
    // coincidence length has to be converted: a deviation of `coincidence` across a model of size
    // `scale` subtends `coincidence / scale`. A model too small to divide by leaves no angle to
    // compare against, so the judgement is degenerate rather than guessed.
    let Some(limit) = j.coincidence.over(j.scale) else {
        return Decision::Degenerate;
    };
    escalate(j, limit, |prec| {
        let ph = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| plane_hp(t.0, t.1, t.2, prec);
        let planes = [ph(a), ph(b), ph(c)];
        let dh = normals_det_hp(&planes, prec);
        match dh.sign() {
            Some(pos) => Ok(orient_of(pos)),
            None => Err(dir_gap(&dh, &planes, prec)),
        }
    })
}

/// **How far three plane normals are from being coplanar — as an angle, not a length.**
///
/// The other judges reduce to a distance; this one cannot, because it asks a question about
/// *directions*: whether three planes share a common line direction. Its determinant is the
/// triple product of the normals, which carries each normal's own magnitude — and those are
/// arbitrary, since a plane's coefficients may be scaled freely. Dividing by `|n_a||n_b||n_c|`
/// leaves the triple product of the **unit** normals, which is dimensionless and, near zero, is
/// the sine of the angle by which the third normal misses the other two's plane.
///
/// A caller compares it against `coincidence_precision / scale`: the angle a deviation of the
/// coincidence limit subtends across the model. `None` when a normal cannot be bounded away from
/// zero — a degenerate plane has no direction to be off by.
fn dir_gap(d: &HpBounded, planes: &[[HpBounded; 4]; 3], prec: usize) -> Gap {
    let mut denom = Mag::of(1.0);
    for p in planes {
        // `|n| ≥ max|n_k|`, which is enough and needs no square root. One component clearing zero
        // is all a normal needs, so a shortfall only counts when **every** component fell short —
        // and then the smallest of them is the cheapest way out.
        let mut lo = Mag::ZERO;
        let mut short: Option<usize> = None;
        let mut vanished = 0;
        for c in p.iter().take(3) {
            match denom_lo(c) {
                Ok(b) if lo.lt(b) => lo = b,
                Ok(_) => {}
                Err(Gap::Unresolved { short: s }) => {
                    short = Some(short.map_or(s, |t: usize| t.min(s)));
                }
                Err(_) => vanished += 1,
            }
        }
        if lo.is_zero() {
            return match short {
                Some(short) => Gap::Unresolved { short },
                // Every component of this normal is exactly zero: the plane has no direction.
                None => {
                    debug_assert_eq!(vanished, 3);
                    Gap::Vanished
                }
            };
        }
        denom = denom.times(lo);
    }
    let _ = prec;
    match d.error.over(denom) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}

#[cfg(test)]
#[path = "../tests/frame3/coord.rs"]
mod tests;
