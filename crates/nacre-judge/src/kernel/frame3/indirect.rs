use super::*;
// ---- indirect orient3d: three rotated planes meet at an implicit point ----
//
// A boolean's seam vertex is `∩` of three planes, each a plane through three
// rotated points, so its coordinates are never exact — an **implicit point** `V`.
// The kernel decides `orient3d(V, q, r, s)` without materializing `V`:
// `sign = sign(D)·sign(M)` where `D = det(normals)` and
// `M = (Dvec − D·s)·((q−s)×(r−s))` (Cramer, no division). The f64 filter is over
// **intervals** (value ± tol) — a dynamic filter: interval
// arithmetic is a sound worst-case bound by construction, so no per-predicate bound
// formula is hand-derived. An interval straddling 0 escalates to astro-float from the
// point definitions, exactly as [`orient3d_judge`]. This is the *tol > 0* path; the boolean
// reaches it through `Judge::orient3d`, whose rotated arm hands the seam's planes over as
// `WitnessPoint`s. Soundness is tested by `plane_coefficient_tol_soundness` and
// `indirect_orient3d_soundness`; the constant `mag`-floor policy is indirect-only (distinct from the
// explicit `16·scale³` floor of [`orient3d_judge`]).

/// 3×3 determinant of interval rows.
pub(super) fn det3_iv(r: [[Bounded; 3]; 3]) -> Bounded {
    let m0 = r[1][1].mul(r[2][2]).sub(r[1][2].mul(r[2][1]));
    let m1 = r[1][0].mul(r[2][2]).sub(r[1][2].mul(r[2][0]));
    let m2 = r[1][0].mul(r[2][1]).sub(r[1][1].mul(r[2][0]));
    r[0][0].mul(m0).sub(r[0][1].mul(m1)).add(r[0][2].mul(m2))
}

/// Plane `[a,b,c,d]` (`n·X + d = 0`) through three points, as intervals: `n =
/// (p1−p0)×(p2−p0)`, `d = −n·p0`. Coefficient tol propagates from the point tols
/// through the subtraction/cross/dot — "coefficient tol is a corollary of point tol"
/// (checked by `plane_coefficient_tol_soundness`).
pub(crate) fn plane_iv(p0: &WitnessPoint, p1: &WitnessPoint, p2: &WitnessPoint) -> [Bounded; 4] {
    let (a, b, c) = (p0.realized, p1.realized, p2.realized);
    let e1 = [b[0].sub(a[0]), b[1].sub(a[1]), b[2].sub(a[2])];
    let e2 = [c[0].sub(a[0]), c[1].sub(a[1]), c[2].sub(a[2])];
    let n = [
        e1[1].mul(e2[2]).sub(e1[2].mul(e2[1])),
        e1[2].mul(e2[0]).sub(e1[0].mul(e2[2])),
        e1[0].mul(e2[1]).sub(e1[1].mul(e2[0])),
    ];
    let d = Bounded::new(0.0, 0.0)
        .sub(n[0].mul(a[0]))
        .sub(n[1].mul(a[1]))
        .sub(n[2].mul(a[2]));
    [n[0], n[1], n[2], d]
}

/// The plane's four coefficients realized at `prec` bits from the definitions
/// (ground truth for the coefficient tol; no trig of its own — consumes point
/// `hp_coord`). Public for a caller that meets one plane in many points and realizes it once
/// ([`meet_hp`]).
pub fn plane_hp(
    p0: &WitnessPoint,
    p1: &WitnessPoint,
    p2: &WitnessPoint,
    prec: usize,
) -> [HpBounded; 4] {
    let (a, b, c) = (p0.hp_coord(prec), p1.hp_coord(prec), p2.hp_coord(prec));
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let e1 = [sub(&b[0], &a[0]), sub(&b[1], &a[1]), sub(&b[2], &a[2])];
    let e2 = [sub(&c[0], &a[0]), sub(&c[1], &a[1]), sub(&c[2], &a[2])];
    let n = [
        mul(&e1[1], &e2[2]).sub(&mul(&e1[2], &e2[1]), prec),
        mul(&e1[2], &e2[0]).sub(&mul(&e1[0], &e2[2]), prec),
        mul(&e1[0], &e2[1]).sub(&mul(&e1[1], &e2[0]), prec),
    ];
    let d = HpBounded::exact(BigFloat::from_f64(0.0, prec))
        .sub(&mul(&n[0], &a[0]), prec)
        .sub(&mul(&n[1], &a[1]), prec)
        .sub(&mul(&n[2], &a[2]), prec);
    let [n0, n1, n2] = n; // `d` is already built from them, so the normal moves out rather than copying
    [n0, n1, n2, d]
}

/// `true` when an odd number of the three scale factors is negative — the case where the join
/// comes out pointing the other way. `sign()` reports `true` for positive.
fn negatives_odd(signs: [bool; 3]) -> bool {
    signs.iter().filter(|positive| !**positive).count() % 2 == 1
}

/// **The plane through three points that have no coordinates** — the projective join, the dual
/// of [`cramer_hp`]: that one meets three planes into a homogeneous point `[Dvec : D]`, this
/// joins three such points back into a plane. Both are the same 3×3 determinant expansion read
/// the other way round, which is why the whole thing is four [`det3_big`] calls and no new
/// arithmetic. Consumed by [`judged_coeffs`]' homogeneous route (a judged frame whose defining
/// point is a [`JudgedPoint::Meet`]).
///
/// ★★★ **Nothing is divided, and that is the point.** A point named as an intersection has no
/// rational coordinates to hand out — dividing `Dvec/D` would manufacture an irrational and bake
/// its rounding into everything downstream, which is exactly what the homogeneous form exists to
/// avoid. The plane is built in the same currency: multiplies and subtracts.
///
/// ★★★★ **The scale's sign is removed here, not handed to the caller.** The join of
/// `P_i = D_i·[p_i : 1]` is multilinear in its rows, so the answer is `D0·D1·D2` times the plane
/// through the affine points. A positive multiple is harmless — every question asked of a plane
/// is a sign question — but a **negative** one silently reverses the plane's direction, and
/// direction is what `frame_sign`, outward normals and the whole label frame are built on; the
/// function negates its own result instead of reporting (the kernel's existing idiom:
/// `canonical_plane_coeffs` also settles sign *inside* the value).
///
/// `None` when any `D_i`'s sign is undecided — then it is not even known whether those three
/// planes meet in a point, so there is nothing to normalize against and the caller must climb.
/// Degree 9 in the nine input coefficients.
///
/// (It has no f64/`Bounded` filter twin: the judging table needs no interval-coefficient route
/// at all. The lock lives in the Pure-vs-Meet differential, which pins this
/// function's basis against the affine shortcut end to end.)
pub(super) fn plane_hp_through(
    pts: [(&HpBounded, &[HpBounded; 3]); 3],
    prec: usize,
) -> Option<[HpBounded; 4]> {
    let zero = HpBounded::exact(BigFloat::from_f64(0.0, prec));
    let row = |i: usize| {
        let (d, v) = pts[i];
        [&v[0], &v[1], &v[2], d]
    };
    let (r0, r1, r2) = (row(0), row(1), row(2));
    let minor = |a: usize, b: usize, c: usize| {
        det3_big(
            [
                [r0[a], r0[b], r0[c]],
                [r1[a], r1[b], r1[c]],
                [r2[a], r2[b], r2[c]],
            ],
            prec,
        )
    };
    let mut pi = [
        minor(1, 2, 3),
        zero.sub(&minor(0, 2, 3), prec),
        minor(0, 1, 3),
        zero.sub(&minor(0, 1, 2), prec),
    ];
    if negatives_odd([pts[0].0.sign()?, pts[1].0.sign()?, pts[2].0.sign()?]) {
        for c in &mut pi {
            *c = zero.sub(c, prec);
        }
    }
    Some(pi)
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
/// Dvec[j]/D` (no division taken here). Shared by `indirect_filter` (orient3d) and
/// [`cmp_filter`] (cmp_coord).
pub(crate) fn cramer_iv(planes: [[Bounded; 4]; 3]) -> (Bounded, [Bounded; 3]) {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    let h = |k: usize| Bounded::new(0.0, 0.0).sub(planes[k][3]); // n·X = h, h = −d
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

/// The Cramer `D` **alone** — the determinant of the three planes' normals.
///
/// ★ **A caller that wants only `D` must not go through [`cramer_iv`]**, which builds four
/// determinants and the `h` column for the three it does not need. [`dir_sign_judge`] asks exactly
/// this question ("how does the line `a ∩ b` run relative to `c`"), and going through
/// [`cramer_iv`] would throw away three-quarters of the work on every call.
pub(super) fn normals_det_iv(planes: [[Bounded; 4]; 3]) -> Bounded {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    det3_iv([n(0), n(1), n(2)])
}

/// [`normals_det_iv`] at `prec` bits — the `d` half of [`cramer_hp`] without its `h` column or its
/// three `Dvec` determinants. The escalation is where a wasted determinant costs the most.
pub(super) fn normals_det_hp(planes: &[[HpBounded; 4]; 3], prec: usize) -> HpBounded {
    let n = |k: usize, j: usize| &planes[k][j];
    det3_big(
        [
            [n(0, 0), n(0, 1), n(0, 2)],
            [n(1, 0), n(1, 1), n(1, 2)],
            [n(2, 0), n(2, 1), n(2, 2)],
        ],
        prec,
    )
}

/// The interval f64 filter for `orient3d(V, q, r, s)`, `V = ∩(planes)`. `None` if
/// either `D` or `M` straddles 0 (escalate). Coefficient-direct (no division).
#[cfg(test)]
pub(super) fn indirect_filter(
    planes: [[Bounded; 4]; 3],
    q: [Bounded; 3],
    r: [Bounded; 3],
    s: [Bounded; 3],
) -> Option<Orient> {
    filter_from_cramer(cramer_iv(planes), q, r, s)
}

/// `indirect_filter` with the implicit point's Cramer parts **already in hand**.
///
/// ★ **The split exists because `(D, Dvec)` *is* the point** — it does not mention `q, r, s`, so a
/// caller asking about one point against several query triangles pays for it once. The arrangement's
/// crossing collector asks twice in a row (a segment's two endpoints), which is what
/// [`crate::predicate::ImplicitPoint`] exploits.
fn filter_from_cramer(
    (d, dvec): (Bounded, [Bounded; 3]),
    q: [Bounded; 3],
    r: [Bounded; 3],
    s: [Bounded; 3],
) -> Option<Orient> {
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

/// **The point three planes meet in**, realized at `prec` bits — [`cramer_hp`] divided out, each
/// coordinate with the radius the division carries. `None` when `D` may be zero at this
/// precision: the planes are not proven to meet in a point, and more bits may prove it.
///
/// ★ This is how a meet of carriers whose motion histories differ gets a coordinate: no shared
/// frame states it exactly, and each carrier's witness triangle does ([`plane_hp`]).
pub fn meet_hp(planes: &[[HpBounded; 4]; 3], prec: usize) -> Option<[HpBounded; 3]> {
    let (d, n) = cramer_hp(planes, prec);
    Some([
        n[0].div(&d, prec)?,
        n[1].div(&d, prec)?,
        n[2].div(&d, prec)?,
    ])
}

/// The Cramer parts of `V = ∩(planes)` at `prec` bits: `(D, Dvec)` — the determinant and its
/// numerator vector, each carrying the error radius accumulated along the way. Shared by
/// `indirect_hp` (orient3d) and `cmp_hp_with_gap` (cmp_coord).
///
/// No magnitude bound comes back alongside: the radius rides *with* the value, so there is
/// nothing for a caller to forget to use — a separate bound is how a scale that has already
/// cancelled ends up bounding `M`.
pub(super) fn cramer_hp(planes: &[[HpBounded; 4]; 3], prec: usize) -> (HpBounded, [HpBounded; 3]) {
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let zero = HpBounded::exact(BigFloat::from_f64(0.0, prec));
    // The four matrices are the same twelve values rearranged, so they are built as *views* of
    // the coefficients rather than copies of them — see [`det3_big`].
    let n = |k: usize, j: usize| &planes[k][j];
    let hc = [
        sub(&zero, &planes[0][3]),
        sub(&zero, &planes[1][3]),
        sub(&zero, &planes[2][3]),
    ];
    let d = det3_big(
        [
            [n(0, 0), n(0, 1), n(0, 2)],
            [n(1, 0), n(1, 1), n(1, 2)],
            [n(2, 0), n(2, 1), n(2, 2)],
        ],
        prec,
    );
    let dvec = [
        det3_big(
            [
                [&hc[0], n(0, 1), n(0, 2)],
                [&hc[1], n(1, 1), n(1, 2)],
                [&hc[2], n(2, 1), n(2, 2)],
            ],
            prec,
        ),
        det3_big(
            [
                [n(0, 0), &hc[0], n(0, 2)],
                [n(1, 0), &hc[1], n(1, 2)],
                [n(2, 0), &hc[2], n(2, 2)],
            ],
            prec,
        ),
        det3_big(
            [
                [n(0, 0), n(0, 1), &hc[0]],
                [n(1, 0), n(1, 1), &hc[1]],
                [n(2, 0), n(2, 1), &hc[2]],
            ],
            prec,
        ),
    ];
    (d, dvec)
}

/// The same `sign(D)·sign(M)` realized at `prec` bits (astro-float) — ground truth /
/// escalation. Returns the two determinants as intervals; their radii are what decides whether
/// either sign may be used.
///
/// **Where the cancellation happens.** `row1 = Dvec − D·s` collapses to nothing when the query
/// point coincides with the implicit point, and a declare-0 floor sized off that collapsed value
/// lets a rounding residue through as a confident sign (measured: a floor of `1e-131` passed an
/// `8.6e-78` residue while two other ways of asking the same question answered zero). An interval
/// cannot make that mistake: the radius of `row1` is the sum of what went into it, and a
/// subtraction that cancels leaves the radius behind.
fn indirect_hp(
    planes: [[HpBounded; 4]; 3],
    q: [HpBounded; 3],
    r: [HpBounded; 3],
    s: [HpBounded; 3],
    prec: usize,
) -> (HpBounded, HpBounded, [HpBounded; 3]) {
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let add = |x: &HpBounded, y: &HpBounded| x.add(y, prec);
    let (d, dvec) = cramer_hp(&planes, prec);
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
    // `cross` rides out with the two determinants rather than the finished gap: forming the gap
    // costs a division and an exponent walk, and the path that decides a sign never needs it.
    (d, m, cross)
}

/// **How far the implicit point may be from the triangle's plane**, in the model's own units.
///
/// `row1 = Dvec − D·s = D·(V − s)` and `M = row1 · cross`, so `M = D · (V−s)·cross`. The signed
/// distance from `V` to the plane through `q, r, s` is `(V−s)·n̂ = (V−s)·cross / |cross|`, hence
///
/// ```text
///     distance = M / (|D| · |cross|)
/// ```
///
/// Both denominators are needed: `|D|` because `V` is `Dvec/D` and the point's coordinates carry
/// that division, `|cross|` because the dot product carries the triangle's area. Bounded above
/// from `M`'s radius and below by the denominators' own lower bounds. `None` when either cannot
/// be kept away from zero — three near-parallel planes have no well-defined meeting point, and a
/// degenerate triangle no plane.
fn plane_gap(m: &HpBounded, d: &HpBounded, cross: &[HpBounded; 3], prec: usize) -> Gap {
    let d_lo = match denom_lo(d) {
        Ok(b) => b,
        Err(g) => return g,
    };
    match distance_bound(m.error, cross, prec) {
        Gap::Of(b) => match b.over(d_lo) {
            Some(b) => Gap::Of(b),
            None => Gap::Vanished,
        },
        g => g,
    }
}

/// CIP indirect `orient3d(V, q, r, s)`, `V = ∩(3 planes)` — each plane through three
/// rotated points, the triangle three rotated points. Interval filter → astro-float
/// escalation; an undecided `D` or `M` becomes the distance from the implicit point to the
/// triangle's plane and is answered by [`escalate`]. The seam-vertex analogue of
/// [`orient3d_judge`]. The boolean reaches its body (`orient3d_from_cramer`) through
/// [`crate::predicate::Judge::orient3d`].
#[allow(clippy::too_many_arguments)]
pub fn indirect_orient3d_judge(
    plane_a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    q: &WitnessPoint,
    r: &WitnessPoint,
    s: &WitnessPoint,
    j: Standard,
) -> Decision {
    let iv = [
        plane_iv(plane_a.0, plane_a.1, plane_a.2),
        plane_iv(plane_b.0, plane_b.1, plane_b.2),
        plane_iv(plane_c.0, plane_c.1, plane_c.2),
    ];
    orient3d_from_cramer(cramer_iv(iv), plane_a, plane_b, plane_c, q, r, s, j)
}

/// [`indirect_orient3d_judge`] with the implicit point's Cramer parts **already in hand** — the one
/// body of the certified `orient3d`, reached by both [`crate::predicate::Judge::orient3d`] and
/// [`crate::predicate::ImplicitPoint::orient3d`].
///
/// The `WitnessPoint` definitions stay in the signature because the escalation still needs them: `plane_hp`
/// realizes the coefficients afresh at each precision, and an interval cannot be sharpened after
/// the fact.
#[allow(clippy::too_many_arguments)]
pub(crate) fn orient3d_from_cramer(
    cr: (Bounded, [Bounded; 3]),
    plane_a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    q: &WitnessPoint,
    r: &WitnessPoint,
    s: &WitnessPoint,
    j: Standard,
) -> Decision {
    if let Some(o) = filter_from_cramer(cr, q.realized, r.realized, s.realized) {
        return Decision::Sign(o);
    }
    // Escalate: the same two determinants at prec, each carrying the radius accumulated
    // along its own computation. A determinant whose interval straddles zero is undecided,
    // and what it *did* establish is how far the implicit point may be from the triangle's plane.
    escalate(j, j.coincidence, |prec| {
        let ph = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| plane_hp(t.0, t.1, t.2, prec);
        let (d, m, cross) = indirect_hp(
            [ph(plane_a), ph(plane_b), ph(plane_c)],
            q.hp_coord(prec),
            r.hp_coord(prec),
            s.hp_coord(prec),
            prec,
        );
        match combine(d.sign(), m.sign()) {
            Some(o) => Ok(o),
            None => Err(plane_gap(&m, &d, &cross, prec)),
        }
    })
}

#[cfg(test)]
#[path = "../tests/frame3/indirect.rs"]
mod tests;
