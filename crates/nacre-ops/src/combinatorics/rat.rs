use super::*;
/// **Do these two exact 4-vectors describe one plane, and does `b`'s normal point the same way?**
///
/// `Some(+1)` same plane same sense, `Some(-1)` same plane opposite sense, `None` different planes
/// (or overflow). The comparison is cross-multiplication against a nonzero component — the idiom
/// [`nacre_exact::quad::plane_plane_cylinder`]'s parallel arm already uses ("coincident iff the
/// full 4-vectors are proportional"), spelled once here because a *second* caller now needs it.
pub(super) fn plane_sense(a: &[nacre_exact::Rat; 4], b: &[nacre_exact::Rat; 4]) -> Option<i8> {
    let zero = nacre_exact::Rat::from_int(0);
    let i = (0..3).find(|&k| a[k] != zero)?;
    if b[i] == zero {
        return None;
    }
    for j in 0..4 {
        if b[j].checked_mul(a[i])? != a[j].checked_mul(b[i])? {
            return None;
        }
    }
    Some(if (a[i] > zero) == (b[i] > zero) {
        1
    } else {
        -1
    })
}

/// **An operand vertex's own pierce name, restated in this arrangement's class space.**
///
/// ★★★★★ **The second half of a correspondence the forward direction gets for free.** When a
/// boolean *mints* a [`nacre_topo::Vertex::Pierce`] it writes the two planes as **its own
/// classes' representative surfaces**, so restating class order as handle order is the only
/// correction it needs ([`nacre_topo::QuadRoot::canonical`], which the assembly calls). Coming back
/// the other way the handles are **given**, and they may be a surface that merged into a class
/// under a different representative — and, because a class holds faces whose normals oppose, under
/// the **opposite sign**. `Vertex::Pierce`'s own doc says this correspondence "has to be
/// established a second time"; this is that time.
///
/// ★★ **Both corrections are the same rule.** `Lo`/`Hi` are the order along `ℓ = n₁ × n₂`, and
/// `plane_plane_cylinder` fixes the base by `{n₁·x = −d₁, n₂·x = −d₂, ℓ·x = 0}` — a condition
/// `−ℓ` satisfies identically. So negating **either** normal leaves the two planes, and the base,
/// exactly where they were and only reverses `ℓ`: the two roots trade places, which is what
/// `flipped` says. Swapping the pair reverses `ℓ` too (the derivation `canonical` already carries).
/// ⇒ **flip once per reversal, and an even number of reversals is no flip at all.** The swap is
/// [`NodeId::pierce`]'s to count; the two signs are this function's.
/// ☑ The cylinder needs no correction: negating its axis direction does not move the surface, so
/// the two roots are the same two points in the same order.
pub(crate) fn pierce_name_from_def(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    v: Handle<Vertex>,
    cyl: usize,
    candidates: &[usize],
) -> Option<NodeId> {
    let nacre_topo::Vertex::Pierce { planes, root, .. } = *model.vertex(v) else {
        return None;
    };
    // Which candidate class each stored handle *is*, and with which sense. The match decides the
    // correspondence and the sign in one comparison — asking them separately would be two chances
    // to disagree.
    let mut seen: [Option<(usize, i8)>; 2] = [None, None];
    for (i, &h) in planes.iter().enumerate() {
        let name = model.world_plane_name(h)?;
        let c = name.narrow()?;
        for &k in candidates {
            let Some(sense) = plane_sense(c, &class_coeffs_rat(jd, k)?) else {
                continue;
            };
            if seen[i].is_some() {
                return None; // one handle answering to both classes is not a correspondence
            }
            seen[i] = Some((k, sense));
        }
    }
    let ((k0, s0), (k1, s1)) = (seen[0]?, seen[1]?);
    if k0 == k1 {
        return None;
    }
    let root = if s0 == s1 { root } else { root.flipped() };
    Some(NodeId::pierce(k0, k1, cyl, root))
}

/// **Does this class carry that cylinder's *circle*?** — its normal is parallel to the axis, so
/// the section is a circle and not an ellipse. Exact and **total** ([`nacre_exact::parallel_rat`]
/// clears denominators into `BigInt`), so a caller's `false` means the geometry, never the width.
///
/// ★★★★★ **The one place this rule is named, and it is load-bearing far past its callers.**
/// Nine sites reason from "a circle's class is ⊥ to its axis" — [`crate::nesting`]'s rim witnesses
/// and `disk_in_disk`, the segment-vs-arc turn sign below, `segment_meets_cylinder`'s
/// precondition, `circle_crosses_ruling`'s extent derivation, an `unreachable!` in the
/// arrangement's circle crossings, and three of the merge road's f64 arguments — and the
/// proposition they share is held by the gate.
///
/// **What actually holds it** ([`crate::planes::cylinder_gate`]): the gate runs whenever the input
/// has a cylinder at all and sweeps **every (plane class × cylinder)** pair; an oblique pair whose
/// lateral faces cannot be *proved* to miss the plane is refused
/// ([`crate::RejectReason::ObliqueCylinderCut`] — the section would be an ellipse). So a
/// class that carries a circle is a class whose plane **meets** that cylinder, and had it been
/// oblique the gate would already have refused. ★ The seated producer has a second, stronger
/// reason that survives a gate change: `LoopRing::Circle` is minted only from a **single closed
/// edge**, so the cell's boundary is a whole circle lying in the class — and a circle determines
/// its own plane, while a cylinder's *circular* section is ⊥ to the axis.
///
/// ⚠ **The argument is easy to misread**,
/// which is why the two consumers that would answer *silently wrong* now ask instead of assume.
/// ⚠ The relation calls a **zero** vector parallel to everything (`Across`), so a zero normal or
/// axis answers `true` here — "carries" is the positive reading and a guard spelled `!` fails **open**.
/// Unreachable ([`nacre_topo::CylinderDef::new`] refuses a zero `dir`; the gate refuses a class
/// with no rational description), and the eight inline spellings this will replace already inherit
/// that convention — but the name reads the convention backwards, so it is written down here.
pub(crate) fn class_carries_circle(n: &[nacre_exact::Rat; 3], dir: &[nacre_exact::Rat; 3]) -> bool {
    nacre_exact::axis_relation(n, dir) == nacre_exact::AxisRelation::Across
}

/// **Where a cylinder's axis meets a plane class** — the centre of the circle that cylinder traces
/// on the plane, exact.
///
/// ★ It is rational **whatever way the axis points**: the class has rational coefficients or this
/// says nothing, and the meet is one division. That is why a circle can always name a witness of
/// its own where a *ring* cannot — a ring's corners are pierce points and carry radicals.
///
/// ★ It takes the cylinder's **statement**, not an arrangement element: the coplanar merge asks
/// the same question of a `Bound::Circle` it is carrying into a merged region, and one spelling
/// serves both.
///
/// ★★ It lives here, beside [`class_coeffs_rat`] which it reads, because it has **two** consumers
/// and belongs to neither: `nesting`'s witness supply asks it for a circle's own point, and the
/// arrangement's mixed-class net asks it for the circle to measure against a ruling.
/// Keeping it inside `nesting` would have meant either a second spelling or opening that module's
/// witness atoms, and both are the shape those atoms were made private to prevent.
pub(crate) fn circle_centre_rat(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 3]> {
    use nacre_exact::Rat;
    let coeffs = class_coeffs_rat(jd, wc)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let (o, m) = (def.origin(), def.dir());
    let nm = nacre_exact::dot3_rat(&n, &m)?;
    let no_d = nacre_exact::dot3_rat(&n, &o)?.checked_add(coeffs[3])?;
    let t = Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)?;
    let mut p = o;
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(m[k])?)?;
    }
    Some(p)
}

pub(crate) fn class_coeffs_rat(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
) -> Option<[nacre_exact::Rat; 4]> {
    jd.planes[c].world_rat()
}

/// A node's exact coordinates: the rational meet of its three classes' descriptions. `None` when
/// any class lacks a narrow rational description — the caller declines rather than guessing.
///
/// ★ **This reads the identity directly rather than going through [`three_plane_name`]**, because
/// it is one of the two places whose answer for a second variant is *its own*: a pierce point's
/// coordinate is `a + b√c`, not a rational meet, so this `match` is where that decision belongs —
/// not behind a door whose one answer is "no three-plane name, nothing to give".
///
/// ★★ **A pierce node's `None` is a type fact, not a width decline** — the vessel is rational and
/// the coordinate is not, so no amount of precision reaches it. Resist answering `Some` for a
/// tangency because *its* coordinate happens to be rational: a function right for one root and not
/// the other is the "sometimes right" trap.
///
/// ★★★ **So a caller must not read that `None` as a width decline** (`WitnessNotRational`, "a
/// wider rational would lift this", is **false** for a pierce point).
/// `arrangement::split_circles` asks it of both ends of every segment and uses the answer only as
/// a filter; the order along the segment is asked without it.
pub(crate) fn node_coords_rat(
    jd: &Judge<'_, WorkingPlane>,
    n: NodeId,
) -> Option<[nacre_exact::Rat; 3]> {
    match n.kind() {
        NodeKind::Pierce { .. } => None,
        NodeKind::ThreePlane(t) => nacre_exact::three_planes_rat([
            class_coeffs_rat(jd, t[0])?,
            class_coeffs_rat(jd, t[1])?,
            class_coeffs_rat(jd, t[2])?,
        ]),
    }
}

/// Whether a ring's corners name a class the gate proved silent about some cylinder — the rings
/// whose corners have no rational coordinates for the gate's own reason
/// ([`crate::planes::WorkingCyl::silent`]).
pub(crate) fn names_a_silent_class(cyls: &[crate::planes::WorkingCyl], rb: &[RingEdge]) -> bool {
    rb.iter().any(|e| {
        three_plane_name(e.node)
            .is_some_and(|t| t.iter().any(|&k| cyls.iter().any(|cy| cy.is_silent(k))))
    })
}

/// **Points whose box holds a ring**, in the judge's interval arithmetic — every corner, whatever
/// names it ([`node_point_hp`]), and for every arc the corners of its circle's box (an arc lies in
/// its circle, and the circle in its centre ± `r` on every axis, `r` rounded up to a rational by
/// [`crate::planes::sqrt_at_least`]). A straight edge lies in its two corners' hull and an arc in
/// its circle's box, so the ring lies in the hull of these points, and a point outside their box
/// on some axis is outside the ring. `wc` is the ring's class — an arc's circle centre is where its
/// axis meets that plane. `None` for a point the precision could not separate, or a class or
/// circle with no exact description.
pub(crate) fn ring_box_points_hp(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    rb: &[RingEdge],
) -> Option<Vec<[nacre_exact::HpBounded; 3]>> {
    let prec = jd.standard.prec;
    let mut out = Vec::with_capacity(rb.len());
    for e in rb {
        out.push(node_point_hp(jd, cyls, e.node)?);
        if let Carrier::Arc(ac) = &e.carrier {
            let c = circle_centre_rat(jd, wc, &ac.def)?;
            let r = crate::planes::sqrt_at_least(ac.def.r2().narrow()?)?;
            let lo = [
                c[0].checked_sub(r)?,
                c[1].checked_sub(r)?,
                c[2].checked_sub(r)?,
            ];
            let hi = [
                c[0].checked_add(r)?,
                c[1].checked_add(r)?,
                c[2].checked_add(r)?,
            ];
            for bits in 0..8usize {
                out.push(std::array::from_fn(|k| {
                    let v = if bits >> k & 1 == 0 { lo[k] } else { hi[k] };
                    nacre_exact::HpBounded::of_rat(v, prec)
                }));
            }
        }
    }
    Some(out)
}

/// **A node as an interval, whichever kind of name it is** — [`three_point_hp`] for three planes,
/// [`PierceMeet::point_hp`] for a pierce point — exact or, with a plane of no world statement
/// among its two, an interval box ([`PierceMeet::Bounded`]). `None` for a meet the precision could
/// not separate, a pierce point no box states (a tangency on the interval road, a named meet that
/// fails its exact solve), or a cylinder missing from the table.
pub(crate) fn node_point_hp(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
) -> Option<[nacre_exact::HpBounded; 3]> {
    match n.kind() {
        NodeKind::ThreePlane(t) => three_point_hp(jd, t),
        NodeKind::Pierce { cyl, .. } => pierce_meet(jd, cyl, &cyls.get(cyl)?.def, n)?.point_hp(jd),
    }
}

/// **A three-plane point as an interval** — its three classes' witness triangles
/// ([`nacre_judge::plane_hp`]) met at the model's precision ([`nacre_judge::meet_hp`]): a box that
/// holds the true point whether or not the classes have a world statement. `None` for a meet the
/// precision could not separate.
pub(crate) fn three_point_hp(
    jd: &Judge<'_, WorkingPlane>,
    t: [usize; 3],
) -> Option<[nacre_exact::HpBounded; 3]> {
    let prec = jd.standard.prec;
    let plane = |k: usize| {
        let w = &jd.planes[k].tri_pt3;
        nacre_judge::plane_hp(&w[0], &w[1], &w[2], prec)
    };
    nacre_judge::meet_hp(&[plane(t[0]), plane(t[1]), plane(t[2])], prec)
}

/// **What a question that could not be stated exactly is refused as**, given the classes it read:
/// one the gate passed as silent about the cylinders has no world statement, so the cause is the
/// gate's ([`RejectReason::CylinderGateUndecided`]); otherwise the world holds the classes and the
/// rational vessel does not ([`RejectReason::WitnessNotRational`] — a `Wide` name, an overflow).
pub(crate) fn undecided_reason(
    cyls: &[crate::planes::WorkingCyl],
    classes: &[usize],
) -> RejectReason {
    if classes
        .iter()
        .any(|&c| cyls.iter().any(|cy| cy.is_silent(c)))
    {
        RejectReason::CylinderGateUndecided
    } else {
        RejectReason::WitnessNotRational
    }
}

/// **The sign of a rational plane `a·x + b·y + c·z + d` at an interval point** — `Some(true)` when
/// positive, `None` when the interval straddles zero.
pub(crate) fn plane_sign_hp(
    plane: &[nacre_exact::Rat; 4],
    x: &[nacre_exact::HpBounded; 3],
    prec: usize,
) -> Option<bool> {
    use nacre_exact::HpBounded;
    let mut v = HpBounded::of_rat(plane[3], prec);
    for k in 0..3 {
        v = v.add(&HpBounded::of_rat(plane[k], prec).mul(&x[k], prec), prec);
    }
    v.sign()
}

/// **A node's realized coordinate, whichever kind of name it is** — the `f64` sibling of
/// [`node_coords_rat`] and [`pierce_point`], which each answer for one variant only.
///
/// ★★ **It lives here because the `match` does.** Reaching into a [`NodeId`] variant outside this
/// file (and `reuse`) is what [`three_plane_name`]'s gate forbids, and the first draft of this
/// function sat in `arrangement` and broke it — with the whole suite green, exactly as that gate's
/// doc predicts. The two roads it dispatches between are already both here, so this is where the
/// third question about the same name belongs.
///
/// ★ `cfg(test)` only while the audit is its one consumer. The seam table's `Vertex` minting
/// asks for the same pair of roads and will want it in production.
#[cfg(test)]
pub(crate) fn node_point_f64(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
) -> Option<[f64; 3]> {
    match n.kind() {
        NodeKind::ThreePlane(t) => nacre_geom::intersect::three_planes(
            &jd.planes[t[0]].plane,
            &jd.planes[t[1]].plane,
            &jd.planes[t[2]].plane,
        )
        .map(|p| p.as_array()),
        NodeKind::Pierce { cyl, .. } => pierce_point(jd, cyl, &cyls[cyl].def, n),
    }
}

/// **A pierce node's realized coordinate, derived from the name.**
///
/// The sibling of [`node_coords_rat`] for the other variant, and the two return types *are* the
/// distinction: a pierce coordinate is `a + b√c`, so the rational vessel next door cannot hold it
/// and only the cache can.
///
/// ★★ **It re-solves from the name rather than taking the producer's `(line, s)`** — the truth is
/// the definition and the coordinate is its cache. That is also what makes canonicalization
/// load-bearing: a name whose root failed to follow its pair through the sort designates the
/// *other* crossing, and the point moves where a test can see it.
///
/// ★ **Strict about a tangency.** `Double` is answered only by `Tangent` and `Lo`/`Hi` only by
/// `Pair`, and the mismatches are `None` rather than a nearest guess — so a second name for a
/// tangency's one point is unrepresentable here, not merely discouraged.
///
/// `def` must be the cylinder class `n` names; the caller holds the class table's row and this has
/// no way to look one up, so `cyl` comes with it and the two are checked against the name rather
/// than promised — a mismatched pair would otherwise realize a real point of the *wrong* cylinder.
/// `None` is a three-plane node, a class with no rational description, a root the meet does not
/// have, or checked-`Rat` overflow.
pub(crate) fn pierce_point(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<[f64; 3]> {
    let (line, s) = pierce_meet(jd, cyl, def, n)?.into_exact()?;
    Some(nacre_exact::quad::branch_point_f64(&line, &s))
}

/// **The exact half of [`pierce_point`]** — the `(line, s)` the name designates, before it is
/// realized.
///
/// ★ «the definition is the truth, the coordinate a cache»: the pair *is* the point and the `[f64;
/// 3]` beside it is its
/// realization, so the two are one function split in the middle rather than two solves. Every
/// exact question about a pierce point — its order along the line, its side of a plane, its θ about
/// the seam — takes this and never the realization.
/// **A pierce corner as a rational point, when it is one** — the meet line's point at its root,
/// for a root [`nacre_exact::quad::QuadVal::as_rat`] can state (a wall through or perpendicular
/// to the axis, a tangent wall's double root); `None` for any other corner or name.
///
/// ★ A **witness supply**, not a coordinate vessel: [`node_coords_rat`]'s `None` for a pierce node
/// is a type fact ("the coordinate is `a + b√c`") that its callers route on, and it must stay so.
/// This answers a different question — "is there a rational point *here* to cast from?" — and is
/// total over its input: the corners it cannot state simply do not join the probe list, the way
/// [`edge_interior_points`] abstains per edge. A
/// half-cylinder prism's cap has two corners, both pierce, both rational, and no other point.
pub(crate) fn pierce_coords_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
) -> Option<[nacre_exact::Rat; 3]> {
    let (_, cyl, _) = pierce_name(n)?;
    let (line, s) = pierce_meet(jd, cyl, &cyls.get(cyl)?.def, n)?.into_exact()?;
    let sv = s.as_rat()?;
    let (b, d) = (line.base(), line.dir());
    let mut p = b;
    for k in 0..3 {
        p[k] = p[k].checked_add(sv.checked_mul(d[k])?)?;
    }
    Some(p)
}
/// **Every rational point that names the interior of this *straight* ring edge** — the one rule,
/// with one arm per way the two ends can be described.
///
/// ★★★★★ **One sentence, which must not be four spellings.** "A chord's midpoint" and "an edge's
/// interior point" are one rule — a chord rule refuses `Carrier::Arc` in as many words, so both
/// are *a point inside a straight edge*, differing only in how the ends are named. Spelled per
/// supply, it leaves a road with no edge witness at all, and a planar component whose every
/// corner grazes then has nothing left to say and refuses `NoClearRay` where the shape's truth is
/// `SelfTouchingResult`.
///
/// | arm | the ends | why it is inside |
/// |---|---|---|
/// | [`conjugate_midpoint`] | one solve's two roots (`Lo`/`Hi`) | the shared `mid`, `disc > 0` |
/// | [`pierce_ends_between`] | two solves on one line | a rational verified strictly between |
/// | [`rational_ends_midpoint`] | both rational | the midpoint of two rationals |
///
/// ★★ **An iterator, not an `Option`, and that is load-bearing.** The first two arms **both** match
/// a conjugate-rooted edge — the two roots share the plane pair and the cylinder, so `pierce_meet`
/// hands each end the same line, which is all the second arm asks — and they name **different**
/// points (a shared `mid` against a realized-then-verified midpoint). Folding them into one answer
/// would delete a witness silently.
///
/// ☑ **Measured over the lib + census corpus**, because the argument above is read off the
/// guards and a reader should not have to re-derive it: **284** edges where the first two arms both
/// answer and **0** where the first answers alone — so a fold to "the first arm that matches" would
/// drop 284 points. Those two are **stable across runs**; the other two are not, because proptest
/// fixtures reach this function, so they are given as orders: ~2·10³ edges where only the second arm
/// answers, and **~10⁴** where only the third — by a wide margin the
/// largest supply. (Three runs of the same tree: 10,963 / 13,027 / 11,075 for the third.)
///
/// ★ **The order**: conjugate, then between, then the rational midpoint (which is disjoint
/// from both — a pierce end is not rational).
///
/// **Every arm lands *on* the edge**, which is what lets the component road wrap these in
/// [`Probe::Coord`]: that type's invariant is a point **on** the boundary, and an interior witness
/// "could be separated from the boundary by another component's wall".
pub(crate) fn edge_interior_points(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &RingEdge,
) -> impl Iterator<Item = [nacre_exact::Rat; 3]> {
    [
        conjugate_midpoint(jd, cyls, e),
        pierce_ends_between(jd, cyls, e),
        rational_ends_midpoint(jd, e),
    ]
    .into_iter()
    .flatten()
}

/// **Arm b of [`edge_interior_points`] — a rational point strictly inside a straight edge whose
/// two ends are pierce corners on one line.** The chord witness needs the two ends to be one
/// solve's two roots; a cap's
/// section between the rulings of two *coaxial* cylinders — a bore inside a fillet, cut by a wall
/// within both radii — has its ends on two solves, one radical each, and every corner of that cell
/// irrational. Both ends still lie on one rational line (the pair `{wc, wall}`'s meet, the same
/// parametrization from either solve), so a **rational parameter between the two** names a point
/// of the edge's interior exactly: chosen by the realized midpoint, then **verified** against each
/// end in its own radical ([`rational_between`]). `None` for any other edge shape, or when the
/// two solves do not parametrize one line.
pub(crate) fn pierce_ends_between(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &RingEdge,
) -> Option<[nacre_exact::Rat; 3]> {
    if !matches!(e.carrier, Carrier::Plane { .. }) {
        return None;
    }
    let (pa, ca, _) = pierce_name(e.node)?;
    let (pb, cb, _) = pierce_name(e.to)?;
    if pa != pb {
        return None;
    }
    let (la, sa) = pierce_meet(jd, ca, &cyls.get(ca)?.def, e.node)?.into_exact()?;
    let (lb, sb) = pierce_meet(jd, cb, &cyls.get(cb)?.def, e.to)?.into_exact()?;
    if la.base() != lb.base() || la.dir() != lb.dir() {
        return None;
    }
    let t = rational_between(&sa, &sb)?;
    let (b, d) = (la.base(), la.dir());
    let mut p = b;
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(d[k])?)?;
    }
    Some(p)
}

/// A rational strictly between two quadratic values that need not share a radical: the realized
/// midpoint, taken exactly as the `f64` it is (`Rat::try_from_f64`), then **verified** against
/// each end in that end's own radical — a comparison with a rational is always formable. If the
/// midpoint lands outside (ends closer than the realization resolves), a few bisections toward
/// the realized interval's middle are tried; `None` when none is inside.
fn rational_between(
    a: &nacre_exact::quad::QuadVal,
    b: &nacre_exact::quad::QuadVal,
) -> Option<nacre_exact::Rat> {
    use nacre_exact::{Orient, Rat, quad::QuadVal};
    let (mut lo, mut hi) = (a.to_f64(), b.to_f64());
    if lo > hi {
        std::mem::swap(&mut lo, &mut hi);
    }
    let inside = |t: Rat| -> Option<bool> {
        let q = QuadVal::from_rat(t);
        let da = q.checked_sub(a)?.sign();
        let db = q.checked_sub(b)?.sign();
        // strictly between: on opposite sides of the two ends
        Some(matches!(
            (da, db),
            (Orient::Positive, Orient::Negative) | (Orient::Negative, Orient::Positive)
        ))
    };
    let mut mid = (lo + hi) / 2.0;
    for _ in 0..8 {
        let t = Rat::try_from_f64(mid)?;
        if inside(t)? {
            return Some(t);
        }
        // the realization put it outside: pull toward the interval's middle
        mid = (mid + (lo + hi) / 2.0) / 2.0;
    }
    None
}

/// **Every rational point a ring edge offers a witness supply** — its start corner, then the points
/// its interior names ([`edge_interior_points`]).
///
/// ★ **The whole per-edge chain, in one place**, so that `nesting`'s witness supply and its
/// diagnostic twin do not each hold a copy — the interior half being one rule would still leave
/// the *corner* half written twice. The order: the corner first — rational if the node has
/// a three-plane name, else the pierce root's coordinates — then the interior arms.
///
/// ⚠ **The component road one dimension up does not call this**: it already offers every corner as
/// a [`Probe::Named`], which is exact without coordinates at all, so it takes
/// [`edge_interior_points`] alone rather than minting a second description of a point it has.
pub(crate) fn edge_witness_points(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &RingEdge,
) -> impl Iterator<Item = [nacre_exact::Rat; 3]> {
    node_coords_rat(jd, e.node)
        .or_else(|| pierce_coords_rat(jd, cyls, e.node))
        .into_iter()
        .chain(edge_interior_points(jd, cyls, e))
}

/// **Arm a of [`edge_interior_points`] — the midpoint of a ring edge whose two ends are one
/// solve's two roots** — rational, exactly, and strictly between them.
///
/// ★★★★★ **A chord names its own middle.** `plane_plane_cylinder` builds the pair as
/// `lo = (mid, −half, disc)` and `hi = (mid, +half, disc)` — **one `mid`, shared** — so a segment
/// whose ends are that pair has `base + s.a()·dir` for its midpoint whichever end is asked, with
/// no second solve and no approximation. `disc > 0` for a `Pair`, so it is strictly inside.
///
/// **Conjugacy is a question about names, not values**: the two ends must carry the same canonical
/// plane pair, the same cylinder, and the two roots. That is also what keeps a *piece* of a chord
/// out — an edge cut short by another feature has a different node at one end, and
/// `split_at_crossings` states that "whether a crossing is on this segment is the caller's
/// question", which this answers by refusing to guess.
pub(crate) fn conjugate_midpoint(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &RingEdge,
) -> Option<[nacre_exact::Rat; 3]> {
    use nacre_topo::QuadRoot::{Hi, Lo};
    let (pa, ca, ra) = pierce_name(e.node)?;
    let (pb, cb, rb) = pierce_name(e.to)?;
    if pa != pb || ca != cb || !matches!((ra, rb), (Lo, Hi) | (Hi, Lo)) {
        return None;
    }
    // ★★★★★ **The edge must *be* the segment between its ends, and an arc is not.** Two ends can
    // be one solve's two roots and still be joined by a **curve**: a plane cutting a circle names
    // both crossings, and *either* arc between them carries that same pair of names. The chord's
    // midpoint is then a point strictly inside the circle and **not on this ring at all** — and a
    // point off the ring is not a witness for it, since containment is read from a point *of* `a`
    // and a point in `a`'s interior answers a different question wherever `b` nests inside it.
    // ☑ Measured over the whole lib suite: 120 acceptances, **not one** curved carrier — an
    // all-arc ring is answered by `ring_own_circle` one arm up, and every ring that reaches here
    // offered two straight chords. The guard states the precondition; it does not describe a
    // population.
    if matches!(e.carrier, Carrier::Arc(_)) {
        return None;
    }
    let def = &cyls.get(ca)?.def;
    let (line, s) = pierce_meet(jd, ca, def, e.node)?.into_exact()?;
    let (b, d) = (line.base(), line.dir());
    let t = s.a();
    let mut p = [b[0], b[1], b[2]];
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(d[k])?)?;
    }
    // ★★★★★ **Asserted where the fact is made, not where it is consumed.** Both claims this
    // function rests on are checkable here and nowhere cheaper: that `MeetLine`'s `base`/`dir`
    // really do parameterize the two planes' meet (so `base + t·dir` is on both), and that the
    // shared `a()` lands **strictly between** the two roots (so it is strictly inside the
    // cylinder, which is what `disc > 0` buys). A producer change that broke either would
    // otherwise surface as a wrong containment answer two layers up.
    debug_assert!(
        [pa[0], pa[1]].iter().all(|&k| {
            class_coeffs_rat(jd, k).is_none_or(|c| {
                let n = [c[0], c[1], c[2]];
                nacre_exact::dot3_rat(&n, &p)
                    .and_then(|v| v.checked_add(c[3]))
                    .is_none_or(|v| v == nacre_exact::Rat::from_int(0))
            })
        }),
        "a chord midpoint is on both of its planes"
    );
    debug_assert_eq!(
        nacre_exact::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.r2()),
        nacre_exact::Orient::Negative,
        "a chord midpoint is strictly inside the cylinder"
    );
    Some(p)
}

/// **The midpoint of a straight ring edge whose two ends are rational** — the arm the
/// other two never covered, because both of them start by asking for a pierce name.
///
/// The plainest case there is, and the one a planar component is made of: two three-plane corners
/// joined by a straight step. Its midpoint is the average of two rationals, exact, and strictly
/// between the ends, so it is a point of the edge's interior — on the ring, which is what a
/// witness for the ring has to be.
fn rational_ends_midpoint(
    jd: &Judge<'_, WorkingPlane>,
    e: &RingEdge,
) -> Option<[nacre_exact::Rat; 3]> {
    if !matches!(e.carrier, Carrier::Plane { .. }) {
        return None;
    }
    let (a, b) = (node_coords_rat(jd, e.node)?, node_coords_rat(jd, e.to)?);
    let half = nacre_exact::Rat::new(1, 2)?;
    let mut p = [nacre_exact::Rat::from_int(0); 3];
    for k in 0..3 {
        p[k] = a[k].checked_add(b[k])?.checked_mul(half)?;
    }
    // ★ Asserted where the fact is made: the two ends share the two planes this edge rides, so the
    // midpoint is on both of them. The same postcondition `conjugate_midpoint` states below, in the
    // vocabulary this arm's ends come in.
    debug_assert!(
        {
            let shared = three_plane_name(e.node)
                .zip(three_plane_name(e.to))
                .map(|(x, y)| x.into_iter().filter(|k| y.contains(k)).collect::<Vec<_>>());
            shared.is_none_or(|ks| {
                ks.iter().all(|&k| {
                    class_coeffs_rat(jd, k).is_none_or(|c| {
                        let n = [c[0], c[1], c[2]];
                        nacre_exact::dot3_rat(&n, &p)
                            .and_then(|v| v.checked_add(c[3]))
                            .is_none_or(|v| v == nacre_exact::Rat::from_int(0))
                    })
                })
            })
        },
        "an edge midpoint is on the planes its two ends share"
    );
    Some(p)
}

/// **A pierce node's meet, re-solved from its name** ([`PierceMeet`]) — the two classes' world
/// coefficients and the cylinder, never a producer's coordinate.
///
/// `None` folds four causes, and only one of them reaches here on the product road:
/// * **a three-plane node** — a point on a cylinder is represented on the cylinder (`Aliases`'
///   `rep_rank`: a pierce name outranks a three-plane one), so a ruling's or an arc's end is a
///   pierce name; only a hand-built fixture hands one in. `arc_at` and `coord_key` ask the kind
///   first regardless;
/// * **a named class with no narrow world coefficients** — a `Wide` name. The cylinder gate asks
///   every class of a model with a cylinder for them and refuses before any arrangement. A class
///   with **no** world name is not this cause: its pierce points are solved as interval boxes
///   ([`PierceMeet::Bounded`]), and `None` there is a box the precision cannot form or a root an
///   interval cannot state (a tangency);
/// * **an overflow** in `plane_plane_cylinder` — reachable (long decimals at a small scale), and
///   the one cause `WitnessNotRational`, the name most callers give this `None`, is true of;
/// * **a meet without the name's root** — the name was minted from this same meet, of the same
///   coefficients by the same function, with the pair and the root carried together
///   (`QuadRoot::canonical`); only a naming defect produces it.
///
/// Measured over the suite, the census and the ignored sweep: no cause fires on the product road —
/// the only hits are the first cause, from `a_circle_meeting_a_ruling_is_refused_by_name`'s
/// hand-built ruling.
///
/// ⚠ A new node kind — a point two curved surfaces meet at — falls into the first arm and reaches
/// those callers as `WitnessNotRational`, a width sentence for a naming fact. Give it its own cause
/// here before it ships.
pub(crate) fn pierce_meet(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<PierceMeet> {
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    let (planes, root) = match n.kind() {
        NodeKind::Pierce {
            planes,
            cyl: named,
            root,
        } => {
            debug_assert_eq!(
                named, cyl,
                "the def handed in is not the cylinder the name says"
            );
            (planes, root)
        }
        NodeKind::ThreePlane(_) => return None,
    };
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    // ★ **The interval road is for a plane with no world statement, and only for it.** A named pair
    // whose exact solve fails (a wide name, an overflow) keeps declining here under its own name:
    // answering it by intervals would put a second number on a population the exact road owns.
    if planes.iter().any(|&c| jd.planes[c].world.is_none()) {
        let prec = jd.standard.prec;
        let line = nacre_exact::quad::plane_plane_cylinder_hp(
            &jd.planes[planes[0]].canon_hp(prec),
            &jd.planes[planes[1]].canon_hp(prec),
            &o,
            &m,
            r2,
            prec,
        )?;
        // No tangency on this road: an interval proves two roots apart, never one double root.
        let at = match root {
            QuadRoot::Lo => line.root(false, prec)?,
            QuadRoot::Hi => line.root(true, prec)?,
            QuadRoot::Double => return None,
        };
        return Some(PierceMeet::Bounded {
            at: Box::new(at),
            planes,
            cyl,
        });
    }
    let (p1, p2) = (
        class_coeffs_rat(jd, planes[0])?,
        class_coeffs_rat(jd, planes[1])?,
    );
    let (line, s) = match (
        nacre_exact::quad::plane_plane_cylinder(&p1, &p2, &o, &m, r2)?,
        root,
    ) {
        (CylinderMeet::Pair { line, s }, QuadRoot::Lo) => (line, s[0]),
        (CylinderMeet::Pair { line, s }, QuadRoot::Hi) => (line, s[1]),
        (CylinderMeet::Tangent { line, s }, QuadRoot::Double) => (line, QuadVal::from_rat(s)),
        _ => return None,
    };
    #[cfg(debug_assertions)]
    debug_assert_bounded_holds_exact(jd, planes, &o, &m, r2, root, &line, &s);
    Some(PierceMeet::Exact(line, s))
}

/// ★ **The interval road answers for the named population too, where it can be checked.** Every
/// exact pierce point is also solved the interval way ([`nacre_exact::quad::plane_plane_cylinder_hp`]
/// on the same coefficients, the same root by name), and the two boxes must overlap — both hold
/// the true point. A convention that differed between the roads (a root order, a line direction)
/// would put the interval point on the other ruling and fail here on the first named crossing.
#[cfg(debug_assertions)]
#[allow(clippy::too_many_arguments)]
fn debug_assert_bounded_holds_exact(
    jd: &Judge<'_, WorkingPlane>,
    planes: [usize; 2],
    o: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    r2: &nacre_exact::BigRat,
    root: nacre_topo::QuadRoot,
    line: &nacre_exact::quad::MeetLine,
    s: &nacre_exact::quad::QuadVal,
) {
    let prec = jd.standard.prec;
    let hi = match root {
        nacre_topo::QuadRoot::Lo => false,
        nacre_topo::QuadRoot::Hi => true,
        nacre_topo::QuadRoot::Double => return,
    };
    let Some(at) = nacre_exact::quad::plane_plane_cylinder_hp(
        &jd.planes[planes[0]].canon_hp(prec),
        &jd.planes[planes[1]].canon_hp(prec),
        o,
        m,
        r2,
        prec,
    )
    .and_then(|l| l.root(hi, prec)) else {
        return;
    };
    let Some(ex) = PierceMeet::Exact(line.clone(), *s).point_hp(jd) else {
        return;
    };
    debug_assert!(
        (0..3).all(|k| ex[k].sub(&at[k], prec).sign().is_none()),
        "the interval road put a named pierce point somewhere the exact road does not"
    );
}

/// **Where a pierce point is** — what [`pierce_meet`] solves a pierce name to, and the one vessel
/// every question about the point reads.
///
/// ★ The questions are this type's methods, so a caller asks *what it wants to know* (a side, an
/// interval) and never which number answers it. A caller that needs the exact pair itself — an
/// `a + b√c` comparison, a rational root — takes [`PierceMeet::into_exact`] or
/// [`PierceMeet::exact`] and declines where there is none.
// The exact arm is the wide one and the common one; its holders (`CoordKey::Pierce`, `PierceAt`)
// already box the whole meet, so boxing the arm again would buy nothing.
#[allow(clippy::large_enum_variant)]
pub(crate) enum PierceMeet {
    /// The meet line of the name's two planes and the root along it: `base + s·dir` with
    /// `s = a + b√c`, both planes stated in the world.
    Exact(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    /// A point one of whose planes has no world statement: an interval box that holds it, with
    /// the name's two planes and its cylinder — what the questions below answer **by name** where
    /// an interval cannot (a point on one of its own planes; two points of one ruling).
    Bounded {
        at: Box<[nacre_exact::HpBounded; 3]>,
        planes: [usize; 2],
        cyl: usize,
    },
}

impl PierceMeet {
    /// The exact pair, borrowed.
    pub(crate) fn exact(
        &self,
    ) -> Option<(&nacre_exact::quad::MeetLine, &nacre_exact::quad::QuadVal)> {
        match self {
            Self::Exact(line, s) => Some((line, s)),
            Self::Bounded { .. } => None,
        }
    }

    /// The exact pair, moved out.
    pub(crate) fn into_exact(
        self,
    ) -> Option<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)> {
        match self {
            Self::Exact(line, s) => Some((line, s)),
            Self::Bounded { .. } => None,
        }
    }

    /// **The side of plane `q` this point lies on, from `q`'s statement** — one `a + b√c` sign
    /// ([`nacre_exact::quad::plane_side`]) against `q`'s world coefficients, bridged from the name's
    /// facing to the outward one ([`outward_fix`]), so it answers in [`side_of`]'s convention.
    /// `None` for a named class with no narrow rational description. Where the point or `q` has no
    /// world statement the sign is read off interval boxes — `q`'s outward plane
    /// ([`WorkingPlane::outward_hp`]) at the point's box, `None` where the box straddles `q` — and a
    /// bounded point is on its own two planes by name.
    pub(crate) fn side(&self, jd: &Judge<'_, WorkingPlane>, q: usize) -> Option<i8> {
        let by_box = |at: &[nacre_exact::HpBounded; 3]| {
            let prec = jd.standard.prec;
            let pl = jd.planes[q].outward_hp(prec);
            pl[0]
                .mul(&at[0], prec)
                .add(&pl[1].mul(&at[1], prec), prec)
                .add(&pl[2].mul(&at[2], prec), prec)
                .add(&pl[3], prec)
                .sign()
                .map(|pos| if pos { 1 } else { -1 })
        };
        match self {
            Self::Bounded { at, planes, .. } => {
                if planes.contains(&q) {
                    return Some(0);
                }
                by_box(at)
            }
            // A plane with no world statement has no coefficients for the exact sign; the point's
            // own box is read against its interval plane instead.
            Self::Exact(..) if jd.planes[q].world.is_none() => by_box(&self.point_hp(jd)?),
            Self::Exact(line, s) => {
                let co = class_coeffs_rat(jd, q)?;
                let fix = outward_fix(jd, q)?;
                Some(
                    fix * match nacre_exact::quad::plane_side(&co, line, s) {
                        nacre_exact::Orient::Positive => 1,
                        nacre_exact::Orient::Negative => -1,
                        nacre_exact::Orient::Zero => 0,
                    },
                )
            }
        }
    }

    /// **This point as an interval** — `base + s·dir` along its exact meet line, the root realized
    /// at the model's precision ([`nacre_exact::realize_quad`]).
    pub(crate) fn point_hp(
        &self,
        jd: &Judge<'_, WorkingPlane>,
    ) -> Option<[nacre_exact::HpBounded; 3]> {
        match self {
            Self::Bounded { at, .. } => Some((**at).clone()),
            Self::Exact(line, s) => {
                let prec = jd.standard.prec;
                let sb = nacre_exact::realize_quad(s, prec)?;
                let (b, d) = (line.base(), line.dir());
                let at = |k: usize| nacre_exact::affine_bounded(b[k], d[k], &sb, prec);
                Some([at(0)?, at(1)?, at(2)?])
            }
        }
    }

    /// The name's two planes and its cylinder, for a point with no exact meet — what the
    /// by-name rules read.
    pub(crate) fn name_planes_and_cyl(&self) -> Option<([usize; 2], usize)> {
        match self {
            Self::Exact(..) => None,
            Self::Bounded { planes, cyl, .. } => Some((*planes, *cyl)),
        }
    }

    /// The name's two planes, for a point with no exact meet.
    pub(crate) fn name_planes(&self) -> Option<[usize; 2]> {
        self.name_planes_and_cyl().map(|(p, _)| p)
    }

    /// **Which of wall `wall`'s two rulings of `def` this point lies on** — the sign of
    /// `(x − o) · (m × n_wall)`, the convention [`ruling_side_signed`] states, with `n_wall` the
    /// wall's coefficients as a pierce name orders them ([`WorkingPlane::canon_hp`]). Exact for an
    /// exact point against a named wall, `None` there for a wall with no narrow rational description
    /// (a `Wide` name — the exact road's population, which the interval road does not answer for, as
    /// [`pierce_meet`] and [`Self::side`] rule); otherwise read off the interval box, `None` where it
    /// straddles the plane through the axis.
    pub(crate) fn ruling_side(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        wall: usize,
        def: &nacre_topo::CylinderDef,
    ) -> Option<i8> {
        if let Some((line, s)) = self.exact()
            && jd.planes[wall].world.is_some()
        {
            return ruling_side_signed(&class_coeffs_rat(jd, wall)?, def, (line, s));
        }
        let prec = jd.standard.prec;
        let at = self.point_hp(jd)?;
        let n = jd.planes[wall].canon_hp(prec);
        let m = def.dir().map(|x| nacre_exact::HpBounded::of_rat(x, prec));
        let o = def
            .origin()
            .map(|x| nacre_exact::HpBounded::of_rat(x, prec));
        let c = [
            m[1].mul(&n[2], prec).sub(&m[2].mul(&n[1], prec), prec),
            m[2].mul(&n[0], prec).sub(&m[0].mul(&n[2], prec), prec),
            m[0].mul(&n[1], prec).sub(&m[1].mul(&n[0], prec), prec),
        ];
        let d = [
            at[0].sub(&o[0], prec),
            at[1].sub(&o[1], prec),
            at[2].sub(&o[2], prec),
        ];
        d[0].mul(&c[0], prec)
            .add(&d[1].mul(&c[1], prec), prec)
            .add(&d[2].mul(&c[2], prec), prec)
            .sign()
            .map(|pos| if pos { 1 } else { -1 })
    }
}

/// **A node as an interval, for a question about cylinder `cyl`** — a three-plane corner
/// ([`three_point_hp`]) or a pierce point of that cylinder ([`PierceMeet::point_hp`]); `None` for a
/// pierce point of another cylinder or a box the precision cannot form.
pub(crate) fn node_hp_on(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<[nacre_exact::HpBounded; 3]> {
    match n.kind() {
        NodeKind::ThreePlane(t) => three_point_hp(jd, t),
        NodeKind::Pierce { cyl: c, .. } if c == cyl => pierce_meet(jd, cyl, def, n)?.point_hp(jd),
        NodeKind::Pierce { .. } => None,
    }
}

/// Whether a node is a pierce point one of whose planes has no world statement — the population
/// the interval readings below are for, and the only one.
pub(crate) fn names_a_nameless_plane(jd: &Judge<'_, WorkingPlane>, n: NodeId) -> bool {
    pierce_name(n)
        .is_some_and(|(planes, _, _)| planes.iter().any(|&c| jd.planes[c].world.is_none()))
}

/// **Is `n` strictly inside the counter-clockwise arc `lo → hi`**, three points of one circle of
/// cylinder `cyl`, read from their interval boxes: travelling counter-clockwise about the axis
/// `m`, `(lo, n, hi)` is in that order exactly when `m · ((n − lo) × (hi − lo)) > 0` — the sign of
/// the triangle's turn about the axis. `None` where the boxes cannot say (two names for one point
/// among them, or too little precision); the interval twin of the name road's
/// `circular_order`, and asked only where a point has no exact meet.
///
/// ★ **`lo == hi` is the closed arc** — the circle a single node cuts, its whole turn from that
/// node back to it — and every other point of the circle is inside it: `n` is inside exactly when
/// it is not the node, which the boxes say by being apart (`None` when they cannot prove it).
pub(crate) fn between_on_circle_hp(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    lo: NodeId,
    hi: NodeId,
    n: NodeId,
) -> Option<bool> {
    let prec = jd.standard.prec;
    let (a, b) = (node_hp_on(jd, cyl, def, lo)?, node_hp_on(jd, cyl, def, n)?);
    if lo == hi {
        if n == lo {
            return Some(false);
        }
        return (0..3)
            .any(|k| b[k].sub(&a[k], prec).sign().is_some())
            .then_some(true);
    }
    let c = node_hp_on(jd, cyl, def, hi)?;
    let m = def.dir().map(|x| nacre_exact::HpBounded::of_rat(x, prec));
    let e1 = [
        b[0].sub(&a[0], prec),
        b[1].sub(&a[1], prec),
        b[2].sub(&a[2], prec),
    ];
    let e2 = [
        c[0].sub(&a[0], prec),
        c[1].sub(&a[1], prec),
        c[2].sub(&a[2], prec),
    ];
    let x = [
        e1[1].mul(&e2[2], prec).sub(&e1[2].mul(&e2[1], prec), prec),
        e1[2].mul(&e2[0], prec).sub(&e1[0].mul(&e2[2], prec), prec),
        e1[0].mul(&e2[1], prec).sub(&e1[1].mul(&e2[0], prec), prec),
    ];
    m[0].mul(&x[0], prec)
        .add(&m[1].mul(&x[1], prec), prec)
        .add(&m[2].mul(&x[2], prec), prec)
        .sign()
}

/// **Where the line `p ∩ q` crosses the arc `a → b` strictly inside it**, for a pair one of whose
/// planes has no world statement — the interval twin of the arrangement's exact `arc_crossings`,
/// with its contract: the roots `{p, q, cyl}` the travelled arc contains, in the order its travel
/// meets them, and a root that *is* the arc's end on the line (`on_end`) not counted. The meet's
/// shape is read off the discriminant's sign; a tangency or an axis-parallel line is a statement an
/// interval cannot make, so those are `None`, as is a root the boxes cannot place against the arc.
#[allow(clippy::too_many_arguments)]
pub(crate) fn arc_crossings_hp(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    p: usize,
    q: usize,
    cyl: usize,
    ccw: bool,
    [a, b]: [NodeId; 2],
    on_end: Option<NodeId>,
) -> Option<Vec<NodeId>> {
    let prec = jd.standard.prec;
    let line = nacre_exact::quad::plane_plane_cylinder_hp(
        &jd.planes[p].canon_hp(prec),
        &jd.planes[q].canon_hp(prec),
        &def.origin(),
        &def.dir(),
        def.r2(),
        prec,
    )?;
    if line.a.sign() != Some(true) {
        return None;
    }
    if !line.disc(prec).sign()? {
        return Some(Vec::new()); // the line misses the circle
    }
    let (lo, hi) = if ccw { (a, b) } else { (b, a) };
    let mut inside: Vec<NodeId> = Vec::with_capacity(2);
    for root in [nacre_topo::QuadRoot::Lo, nacre_topo::QuadRoot::Hi] {
        let node = NodeId::pierce(p, q, cyl, root);
        if let Some(e) = on_end {
            if node == e {
                continue;
            }
            // Two names, one point: an interval cannot say so, and the exact road's `Coincident`
            // answer has no twin here — the boxes must prove the two apart.
            let (x, y) = (
                node_hp_on(jd, cyl, def, e)?,
                node_hp_on(jd, cyl, def, node)?,
            );
            if !(0..3).any(|k| x[k].sub(&y[k], prec).sign().is_some()) {
                return None;
            }
        }
        if between_on_circle_hp(jd, cyl, def, lo, hi, node)? {
            inside.push(node);
        }
    }
    if let [r, other] = inside[..] {
        let r_first = if ccw {
            between_on_circle_hp(jd, cyl, def, a, other, r)?
        } else {
            between_on_circle_hp(jd, cyl, def, other, a, r)?
        };
        if !r_first {
            inside.swap(0, 1);
        }
    }
    Some(inside)
}

/// **A plane's normal in primitive form** — divided by the gcd of its own three components.
///
/// ★★★★★ **The canonicalisation that made the name narrow was over *four* coefficients, and a
/// reader of three does not inherit it.** A [`nacre_exact::PlaneName`] is normalised by
/// clearing denominators, dividing out the **content of all four**, and fixing a sign; so the
/// normal `(a, b, c)` keeps a factor of `gcd(a,b,c) / gcd(a,b,c,d)`. For an axis-aligned class at
/// an offset that needs a long decimal — `z = s`, coefficients `(0, 0, D, −N)` with `gcd(D,N) = 1`
/// — that leftover factor is exactly **the offset's denominator**, and it can be arbitrarily
/// large while the plane itself is the plainest one there is.
///
/// The rescale is by a **positive** rational, so every direction derived from the normal is
/// unchanged and only the width moves: over the distinct class normals sampled from this corpus,
/// the widest component fell from a median of 51 bits to **1**.
///
/// **Total.** Class coefficients arrive as a canonical *integer* 4-vector (`class_coeffs_rat` is a
/// read of `world_rat`, which is `PlaneName::narrow()`), so the gcd is an integer one and the
/// division is exact. A component that is somehow not an integer is returned untouched rather than
/// guessed at — the structural argument can rot without this quietly changing an answer.
fn primitive_normal(n: &[nacre_exact::Rat; 3]) -> [nacre_exact::Rat; 3] {
    // `unsigned_abs` rather than `abs`: the latter panics on `i128::MIN`, and a panic is a worse
    // answer than the wide arithmetic this exists to avoid.
    fn gcd(a: u128, b: u128) -> u128 {
        let (mut a, mut b) = (a, b);
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a
    }
    let mut g = 0u128;
    for c in n.iter() {
        if c.denom() != 1 {
            return *n;
        }
        g = gcd(g, c.numer().unsigned_abs());
    }
    // `g == 0` is the zero normal (no plane, and `of_normal` says so); `g == 1` is already
    // primitive. Both leave the statement alone, and so does the one `u128` that has no `i128`
    // (a lone `i128::MIN` component) — declining to divide is never wrong here.
    let Ok(g) = i128::try_from(g).map(|g| g.max(1)) else {
        return *n;
    };
    if g == 1 {
        return *n;
    }
    core::array::from_fn(|k| nacre_exact::Rat::from_int(n[k].numer() / g))
}

/// **A rational 2D chart of a plane**, for running a parity test in it.
///
/// `e₁ = ê_k × n` for the first basis axis giving a nonzero cross, `e₂ = n × e₁`. The chart is
/// deliberately **not** orthonormal: crossing parity is invariant under any affine isomorphism of
/// the plane, and demanding unit vectors would need square roots that leave the rationals. One
/// copy of this rule, because a second spelling of it is how two consumers start disagreeing
/// about which side of a ring a point is on.
///
/// ★★★★★ **"Not orthonormal" is not "not orthogonal", and three cheaper charts die on the
/// difference.** `e₁·e₂ = e₁·(n × e₁) = 0`, so this is an **orthogonal, non-unit frame
/// of the plane** — and that is what makes [`Self::axes`]'s sentence ("the parity walks its ray
/// along `e₁`") a true statement about the *world*: in a skew frame the direction of "y fixed, x
/// increasing" is not `e₁`. Measured refutations, so the next reader does not re-derive them:
/// - **Drop a coordinate** (`e₁ = ê_i`, `e₂ = ê_j`, zero arithmetic): census **398 → 389 rows**,
///   `arcwalls rrect-box` and `roundplate` falling to `NoClearRay`. Parity is affine-invariant but
///   the **degeneracy pattern is not**, and a corner on the ray is what the census is made of.
/// - **`e₂ = ê_k`** (zero arithmetic, `e₁` untouched, and on the plane `p·e₂_old = |n|²p_k + n_k·d`
///   is a *positive affine* image of `p_k`, so every comparison and `orient2d` sign is identical):
///   refuted because [`ring_interior_candidates`] walks these axes as **3-D directions in the
///   plane** to mint cap witnesses, and `ê_k·n = n_k ≠ 0` leaves it. An axis here is a direction,
///   not only a coordinate functional.
/// - **`e₂ = ê_j × n`** (degree 1 in `n`, so no squaring, and it *is* in the plane): not
///   orthogonal to `e₁`, so it rotates the ray's level set — the same axis the first one died on.
///
/// ⇒ the only change that provably moves nothing is a **positive rescale** of `n`, which is what
/// [`primitive_normal`] does.
pub(crate) struct Chart2dRat {
    e1: [nacre_exact::Rat; 3],
    e2: [nacre_exact::Rat; 3],
}

impl Chart2dRat {
    /// The chart of the plane with normal `n`. `None` on a zero normal, or when `e₂ = n × e₁`
    /// leaves `i128` — which after [`primitive_normal`] means the **primitive** normal is itself
    /// past ~2⁶³, not that a spurious factor rode in on it.
    ///
    /// ★ **The rescale is first, and it is why this is not a behaviour change**:
    /// `ê_k × n` cannot overflow (its factors are 0 and 1) and parallelism is scale-invariant, so
    /// `k` is the same index either way, and both axes come out along the same directions — only
    /// narrower. What the corpus met before was never geometry: an axis-aligned class at an offset
    /// with a long decimal carries that offset's denominator in its normal, and `e₂` squares it.
    pub(crate) fn of_normal(n: &[nacre_exact::Rat; 3]) -> Option<Self> {
        let zero = nacre_exact::Rat::from_int(0);
        let basis = |k: usize| -> [nacre_exact::Rat; 3] {
            let mut e = [zero; 3];
            e[k] = nacre_exact::Rat::from_int(1);
            e
        };
        let n = &primitive_normal(n);
        let e1 = (0..3)
            .filter_map(|k| nacre_exact::cross3_rat(&basis(k), n))
            .find(|e| e.iter().any(|c| *c != zero))?;
        let e2 = nacre_exact::cross3_rat(n, &e1)?;
        Some(Chart2dRat { e1, e2 })
    }

    /// A point's chart coordinates.
    pub(crate) fn project(&self, p: &[nacre_exact::Rat; 3]) -> Option<[nacre_exact::Rat; 2]> {
        Some([
            nacre_exact::dot3_rat(p, &self.e1)?,
            nacre_exact::dot3_rat(p, &self.e2)?,
        ])
    }

    /// The chart's two in-plane axes — the mixed-ring parity walks its ray along `e1` (one
    /// decision rule with this chart, not a second spelling of a basis).
    pub(crate) fn axes(&self) -> (&[nacre_exact::Rat; 3], &[nacre_exact::Rat; 3]) {
        (&self.e1, &self.e2)
    }

    /// A ring of nodes in chart coordinates.
    pub(crate) fn ring(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        ring: &[NodeId],
    ) -> Option<Vec<[nacre_exact::Rat; 2]>> {
        ring.iter()
            .map(|&n| self.project(&node_coords_rat(jd, n)?))
            .collect()
    }
}
