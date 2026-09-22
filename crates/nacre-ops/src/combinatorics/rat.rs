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
    candidates: [usize; 2],
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
        for &k in &candidates {
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
/// ⚠ `parallel_rat` calls a **zero** vector parallel to everything, so a zero normal or axis
/// answers `true` here — "carries" is the positive reading and a guard spelled `!` fails **open**.
/// Unreachable ([`nacre_topo::CylinderDef::new`] refuses a zero `dir`; the gate refuses a class
/// with no rational description), and the eight inline spellings this will replace already inherit
/// that convention — but the name reads the convention backwards, so it is written down here.
pub(crate) fn class_carries_circle(n: &[nacre_exact::Rat; 3], dir: &[nacre_exact::Rat; 3]) -> bool {
    nacre_exact::parallel_rat(n, dir)
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
    jd.planes[c].world_rat
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
    let (line, s) = pierce_meet(jd, cyl, def, n)?;
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
    let (line, s) = pierce_meet(jd, cyl, &cyls.get(cyl)?.def, n)?;
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
    let (la, sa) = pierce_meet(jd, ca, &cyls.get(ca)?.def, e.node)?;
    let (lb, sb) = pierce_meet(jd, cb, &cyls.get(cb)?.def, e.to)?;
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
    let (line, s) = pierce_meet(jd, ca, def, e.node)?;
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

pub(crate) fn pierce_meet(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)> {
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
    let (p1, p2) = (
        class_coeffs_rat(jd, planes[0])?,
        class_coeffs_rat(jd, planes[1])?,
    );
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let (line, s) = match (
        nacre_exact::quad::plane_plane_cylinder(&p1, &p2, &o, &m, r2)?,
        root,
    ) {
        (CylinderMeet::Pair { line, s }, QuadRoot::Lo) => (line, s[0]),
        (CylinderMeet::Pair { line, s }, QuadRoot::Hi) => (line, s[1]),
        (CylinderMeet::Tangent { line, s }, QuadRoot::Double) => (line, QuadVal::from_rat(s)),
        _ => return None,
    };
    Some((line, s))
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
