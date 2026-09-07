//! **Is this cell of a plane class inside that one?** — the one engine, and the one place a
//! witness for it is named.
//!
//! ★★★★★ **Why a module of its own** (cell 13). The question is answered from two roads — the
//! arrangement's nesting (`arrangement::nest_cells`, `innermost_host`) and the coplanar merge's
//! ownership (`boolean::unify_coplanar_faces`) — and each had grown its own copy: five spellings
//! of «what point may I use as a witness», four of «what does that point say about the target»,
//! three of «try them until one decides, and name the failure honestly», and two of the four-arm
//! dispatch itself. The copies drifted, and the shortest of them — the arm that asks a **ring**
//! against a **disk** — was one witness kind wide. A plate with every corner filleted has a cap
//! ring of eight tangencies and no three-plane corner at all, so that arm found nothing and
//! refused, and the plate could not enter any boolean at all.
//!
//! Nothing here is new geometry: every predicate below is one the kernel already had. What is new
//! is that there is **one** of each, behind a module wall, so the sixth copy cannot be written
//! outside it.

use super::*;
use crate::arrangement::MergedCircle;
use crate::combinatorics::{self, NodeId};
use crate::planes::*;
use crate::tolerant::Judge;

/// Whether a circle's **center** lies inside a polygon ring of the class — the containment
/// witness `nest_cells` uses for a circle contour (the loops cannot **cross** — the gate proved
/// clearance or recorded the crossing, and a tangency touches at most at a point — so one point
/// decides). Exact: the center is `axis ∩ W` (rational), the
/// ring corners are rational meets, and the parity runs in a rational 2D basis of `W`
/// (`point_in_ring_2d_rat` — parity is invariant under the affine projection).
///
/// ★ It takes the cylinder's **statement**, not an arrangement element: the coplanar merge asks
/// the same question of a `Bound::Circle` it is carrying into a merged region, and one spelling
/// serves both.
/// **Where a cylinder's axis meets this class** — the circle's centre on the plane, exact.
///
/// ★ It is rational **whatever way the axis points**: the class has rational coefficients or this
/// says nothing, and the meet is one division. That is why a circle can always name a witness of
/// its own where a *ring* cannot — a ring's corners are branch points and carry radicals.
fn circle_centre_rat(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_scalar::Rat; 3]> {
    use nacre_scalar::Rat;
    let coeffs = combinatorics::class_coeffs_rat(jd, wc)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let (o, m) = (def.origin(), def.dir());
    let dot3 = crate::planes::dot3;
    let nm = dot3(&n, &m)?;
    let no_d = dot3(&n, &o)?.checked_add(coeffs[3])?;
    let t = Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)?;
    let mut p = o;
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(m[k])?)?;
    }
    Some(p)
}

/// **The circle a ring *is***, when every one of its edges is an arc of one cylinder and they
/// chain the whole way round — `None` otherwise.
///
/// ☑ **Its clauses beyond "the first edge is an arc" are guards, and none of them fired** over the
/// suite (88 acceptances, 0 rejections past that first test). They are kept because each states a
/// proposition proved somewhere else — two circles on one class cannot meet (the gate), a ring is
/// a chain (measured over 153,798 rings), a mixed sense would retrace one arc — and a guard that
/// stops holding is how a producer change is meant to surface here rather than two layers down.
///
/// ★★★★★ **A cut circle's cell is still a disk, and must be asked a disk's question.** The
/// dispatch below reaches [`circle_center_in_ring`] only for a cell whose half-edge is literally a
/// `Circle`; a circle a wall has **split into arcs** takes the polygon road instead, where the
/// probes are plane-triple names and an arc's corners are branch points — so the list comes out
/// **empty** and the road refuses with a name about rays it never cast. The shape is the same
/// circle either way, and this is what says so.
///
/// **Why the chain is the whole circle.** The gate proves every pair of cylinders clear — their
/// surfaces by more than the radius sum, or their faces along an axis
/// ([`crate::planes::cylinder_gate`], else `CylinderPairContact`) — so two
/// circles on one class **cannot meet** — arcs that chain head-to-tail therefore all ride the same
/// circle, and running one way round (`ccw` all equal — a mixed pair would retrace one arc) closes
/// it exactly once. ☑ Measured before this was written: **every** ring the road refused for an
/// empty probe list is either such a chain — and each of those the circle's own centre answers —
/// or a wall panel with no circle at all, which still refuses and now by its own name
/// ([`crate::RejectReason::RingHasNoWitness`]). In cells: the crossing census's `NoClearRay` 44
/// became **8 built** and **20 renamed**, the rest being the other road's (`boolean.rs`).
/// ★ Counted as *cells*, not as raises — a traced boolean runs twice in `debug`, and
/// `reject_census`'s own note forbids reading raise counts as populations.
fn ring_own_circle<'a>(ring: &'a [combinatorics::RingEdge]) -> Option<&'a nacre_topo::CylinderDef> {
    let n = ring.len();
    if n < 2 {
        return None;
    }
    let arc = |e: &'a combinatorics::RingEdge| match &e.carrier {
        combinatorics::Carrier::Arc(ac) => Some(&**ac),
        _ => None,
    };
    let first = arc(&ring[0])?;
    for (i, e) in ring.iter().enumerate() {
        let a = arc(e)?;
        if a.cyl != first.cyl || a.ccw != first.ccw || e.to != ring[(i + 1) % n].node {
            return None;
        }
    }
    Some(&first.def)
}

pub(crate) fn circle_center_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    def: &nacre_topo::CylinderDef,
    ring: &[combinatorics::RingEdge],
) -> Result<bool, BoolError> {
    // ★ Every refusal below is a **value** that could not be formed exactly — a class with no
    // narrow description, a coordinate past `Rat` — which is the road's name, not the gate's.
    // (The gate's own questions are signs and were made total; borrowing its name here pointed
    // at a layer that had already answered.)
    let undecided = || reject(RejectReason::WitnessNotRational);
    let center = circle_centre_rat(jd, wc, def).ok_or_else(undecided)?;
    rational_point_in_ring(jd, cyls, wc, &center, ring)?.ok_or_else(undecided)
}

/// **A rational point strictly inside a straight edge whose two ends are branch corners on one
/// line** (cell ⑩). The chord witness needs the two ends to be one solve's two roots; a cap's
/// section between the rulings of two *coaxial* cylinders — a bore inside a fillet, cut by a wall
/// within both radii — has its ends on two solves, one radical each, and every corner of that cell
/// irrational. Both ends still lie on one rational line (the pair `{wc, wall}`'s meet, the same
/// parametrization from either solve), so a **rational parameter between the two** names a point
/// of the edge's interior exactly: chosen by the realized midpoint, then **verified** against each
/// end in its own radical ([`rational_between`]). `None` for any other edge shape, or when the
/// two solves do not parametrize one line.
fn edge_interior_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &combinatorics::RingEdge,
) -> Option<[nacre_scalar::Rat; 3]> {
    if !matches!(e.carrier, combinatorics::Carrier::Plane { .. }) {
        return None;
    }
    let (pa, ca, _) = combinatorics::branch_name(e.node)?;
    let (pb, cb, _) = combinatorics::branch_name(e.to)?;
    if pa != pb {
        return None;
    }
    let (la, sa) = combinatorics::branch_meet(jd, ca, &cyls.get(ca)?.def, e.node)?;
    let (lb, sb) = combinatorics::branch_meet(jd, cb, &cyls.get(cb)?.def, e.to)?;
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
    a: &nacre_scalar::quad::QuadVal,
    b: &nacre_scalar::quad::QuadVal,
) -> Option<nacre_scalar::Rat> {
    use nacre_scalar::{Orient, Rat, quad::QuadVal};
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

/// **The midpoint of a ring edge whose two ends are one solve's two roots** — rational, exactly,
/// and strictly between them.
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
fn chord_midpoint_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &combinatorics::RingEdge,
) -> Option<[nacre_scalar::Rat; 3]> {
    use nacre_topo::QuadRoot::{Hi, Lo};
    let (pa, ca, ra) = combinatorics::branch_name(e.node)?;
    let (pb, cb, rb) = combinatorics::branch_name(e.to)?;
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
    if matches!(e.carrier, combinatorics::Carrier::Arc(_)) {
        return None;
    }
    let def = &cyls.get(ca)?.def;
    let (line, s) = combinatorics::branch_meet(jd, ca, def, e.node)?;
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
            combinatorics::class_coeffs_rat(jd, k).is_none_or(|c| {
                let n = [c[0], c[1], c[2]];
                combinatorics::dot3_rat(&n, &p)
                    .and_then(|v| v.checked_add(c[3]))
                    .is_none_or(|v| v == nacre_scalar::Rat::from_int(0))
            })
        }),
        "a chord midpoint is on both of its planes"
    );
    debug_assert_eq!(
        nacre_scalar::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.radius()),
        nacre_scalar::Orient::Negative,
        "a chord midpoint is strictly inside the cylinder"
    );
    Some(p)
}

/// **Is a rational point inside a ring?** — `Ok(None)` when *this point* cannot answer, `Err` when
/// a **value** could not be formed exactly.
///
/// ★★★★★ **The two are different facts and the split is the point of this function.** The body
/// below used to sit inside [`circle_center_in_ring`], where one witness is all there is, so a
/// point that landed *on* the ring could be folded into the same refusal as a class with no
/// rational description. A caller with **several** witnesses must not read them the same way: a
/// point on the ring has the **next witness as its remedy**, a value that cannot be formed does
/// not — the lesson `point_in_component`'s doc states for the road one dimension up.
///
/// ★ **And it dispatches to the strongest road it can.** A ring the rational chart can name takes
/// [`nacre_geom::intersect::point_in_ring_2d_rat`], which is **half-open in y** and so decides
/// even where a ring corner sits on the ray; only a ring with branch corners or arc steps falls
/// to `point_in_mixed_ring`, which abstains there. Anything that hands this a witness gets the
/// better answer for free.
fn rational_point_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    p: &[nacre_scalar::Rat; 3],
    ring: &[combinatorics::RingEdge],
) -> Result<Option<bool>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let coeffs = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let center = *p;
    // The class's rational chart — the one copy of that rule ([`combinatorics::Chart2dRat`]);
    // parity is affine-invariant, so the basis need not be orthonormal.
    // ★ **A ring the chart road cannot name takes the mixed road** (M6-2b chaining ladder,
    // wall 3): branch corners have no rational coordinates and arc steps no straight chart
    // image, so the parity walks the ring step by step in ℚ(√c) instead. Rings the old road
    // could always name still take it — the mixed arm activates on exactly the population the
    // old road refused, which is what keeps every green census row bit-identical.
    if combinatorics::ring_is_mixed(ring) {
        return Ok(combinatorics::point_in_mixed_ring(
            jd, cyls, &coeffs, &center, ring,
        ));
    }
    let chart = combinatorics::Chart2dRat::of_normal(&n).ok_or_else(undecided)?;
    let p2 = chart.project(&center).ok_or_else(undecided)?;
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let ring2 = chart.ring(jd, &nodes).ok_or_else(undecided)?;
    Ok(
        match nacre_geom::intersect::point_in_ring_2d_rat(p2, &ring2) {
            nacre_geom::intersect::RingSide::Inside => Some(true),
            nacre_geom::intersect::RingSide::Outside => Some(false),
            // On the boundary this witness says nothing — the caller's next one may.
            nacre_geom::intersect::RingSide::OnBoundary => None,
        },
    )
}

/// Whether a polygon contour lies inside a disk — **all but impossible** (its edges ride wall
/// faces, and a wall face that meets the lateral is either recorded as a crossing, which puts its
/// edges on the ruling road, or *tangent*, which since cell ⑥ passes: a corner of such a face can
/// sit exactly on the tangent line and land a ring node exactly on the circle). Computed honestly
/// from one node's radial side rather than assumed, and the `Zero` arm below is what that leftover
/// reaches.
pub(crate) fn node_in_circle(
    jd: &Judge<'_, WorkingPlane>,
    ring: &[combinatorics::RingEdge],
    def: &nacre_topo::CylinderDef,
) -> Result<bool, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    // Any node decides (disjoint loops put every node on one side), so take the first
    // *rational* one — a bitten ring's branch corners have none, but its wall-meet corners do.
    //
    // An all-branch ring (no rational corner at all) still refuses.
    //
    // ★ **A fallback for that was written here and then measured away.** The plan for the road
    // above predicted it would revive contours that arrive here with no rational corner, and the
    // remedy looked free — a ring that *is* a circle can hand over its centre
    // ([`ring_own_circle`] + [`circle_centre_rat`]). ☑ It fired **0** times over the suite: the
    // contours that road revives are asked against *polygons*, which have rational corners. So the
    // note above stands as it was, and the machinery is not here waiting for a population that
    // does not exist.
    let p = ring
        .iter()
        .find_map(|e| combinatorics::node_coords_rat(jd, e.node))
        .ok_or_else(undecided)?;
    match nacre_scalar::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.radius()) {
        nacre_scalar::Orient::Negative => Ok(true),
        nacre_scalar::Orient::Positive => Ok(false),
        // On the surface: a contact the gate admits but this road cannot rank — a *geometric*
        // degeneracy, so it keeps the gate's name ("the gate could not decide exactly") while the
        // width causes above take the road's. ★ It used to say "gate-impossible"; the tangent arm
        // opened in cell ⑥ and a face corner on the tangent line reaches here.
        nacre_scalar::Orient::Zero => Err(reject(RejectReason::CylinderGateUndecided)),
    }
}

/// **Is polygon ring `ra` inside polygon ring `rb`?** — the two-polygon arm of [`cell_in_cell`],
/// written once because the coplanar merge asks the same question of a hole and an outer
/// (cell ⑩: it used to mirror this arm by hand, minus the witness supplies, and named a hole with
/// no three-plane corner `NoClearRay` for a ray never cast). `Ok(None)` is adjacency (a shared
/// node): not nested, not comparable.
pub(crate) fn ring_in_ring_by_witness(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    ra: &[combinatorics::RingEdge],
    rb: &[combinatorics::RingEdge],
) -> Result<Option<bool>, BoolError> {
    if ra.iter().any(|e| rb.iter().any(|f| f.node == e.node)) {
        return Ok(None);
    }
    if combinatorics::ring_is_mixed(rb) {
        let undecided = || reject(RejectReason::WitnessNotRational);
        let coeffs = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
        // From each rational witness of `a` until one ray is clear (`None` is the probe
        // on the ring, a tangent ray, or a seam-incident root — a corner on the ray is
        // decided, cell ②); an exhausted ring is the same degeneracy `ring_in_ring`
        // names. ★ The witnesses, in order: the three-plane corners, the **rational
        // branch corners** ([`combinatorics::branch_coords_rat`] — a half-cylinder
        // prism's cap has no other kind, cell ⑩), and the chords' midpoints
        // ([`chord_midpoint_rat`]) — the supplies the named road below has, so the two
        // roads offer the same points.
        let mut asked = false;
        let witnesses = ra.iter().flat_map(|e| {
            combinatorics::node_coords_rat(jd, e.node)
                .or_else(|| combinatorics::branch_coords_rat(jd, cyls, e.node))
                .into_iter()
                .chain(chord_midpoint_rat(jd, cyls, e))
                .chain(edge_interior_rat(jd, cyls, e))
        });
        for p in witnesses {
            asked = true;
            if let Some(hit) = combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, rb) {
                return Ok(Some(hit));
            }
        }
        // ★ Same split as the named road below: a list that was **empty** is not a list
        // that ran out. ☑ Measured (cell ⑩): the half-cylinder pair reached this arm with
        // an empty list before the branch corners joined it.
        return Err(reject(if asked {
            RejectReason::NoClearRay
        } else {
            RejectReason::RingHasNoWitness
        }));
    }
    // ★ `ring_in_ring` casts from each of `a`'s nodes until one gives a clear ray; an
    // exhausted ring is the genuine degeneracy it rejects for.
    //
    // ★★★★★ **This used to end «— which is also why dropping a branch node from the probe
    // list here is honest», and that was the defect wearing the word.** Dropping them left
    // *nothing*: measured, every refusal this road raised came from a list that was empty
    // before the first cast. What is honest is to say so, which the two arms above now do.
    let probes = combinatorics::three_plane_probes(ra.iter().map(|e| e.node));
    // ★★★★★ **A ring that *is* a circle is asked the circle's question.** The probe list
    // above is plane-triple **names**, and a circle a wall has split into arcs has none —
    // every corner is a branch point — so it comes out empty and `ring_in_ring` refuses
    // with a name about rays it never cast. The cell is a disk either way, and the arm one
    // match-arm up already answers disks exactly ([`circle_center_in_ring`], whose witness
    // is the centre and so is rational whatever way the axis points).
    //
    // ★ **Only where the names run out.** Where they do not, today's road and today's
    // order are untouched — a coordinate is not a name, and a rotated class has names but
    // no rational coefficients (`WorkingPlane::world_rat` is `None` there), so keying this
    // on "no rational witness" instead would divert the rotation sweep's own population.
    // Asking the circle question *always* is the tidier end state and should be measured
    // as an agreement first; this is the strict extension.
    if probes.is_empty() {
        if let Some(def) = ring_own_circle(ra) {
            return Ok(Some(
                circle_center_in_ring(jd, cyls, wc, def, rb)? && !node_in_circle(jd, rb, def)?,
            ));
        }
        // ★★★★★ **And when the corners cannot name a witness, an edge can.** A ring edge
        // whose two ends are one solve's two roots has a **rational midpoint**
        // ([`chord_midpoint_rat`]) — the pair is built from a shared `mid`, so it costs one
        // `branch_meet` and no approximation. That is the wall panel's case: four branch
        // corners, no circle to take a centre from, and two perpendicular traces that are
        // each a whole chord.
        //
        // ★ **`Ok(None)` is this witness's abstention, not the ring's** — the next edge's
        // midpoint may still answer, which is why [`rational_point_in_ring`] hands the two
        // apart. ☑ Measured before this was written: every ring here offers **two**
        // midpoints and the two always agree, and none of them abstains.
        //
        // ★ It sits after [`ring_own_circle`] for the reader's sake only: the two supplies
        // are **disjoint by construction** — that one needs every edge to be an arc, this
        // one needs an edge that is not.
        for e in ra {
            let Some(mid) = chord_midpoint_rat(jd, cyls, e) else {
                continue;
            };
            if let Some(hit) = rational_point_in_ring(jd, cyls, wc, &mid, rb)? {
                return Ok(Some(hit));
            }
        }
        // ★ And a **rational branch corner** is a witness too (cell ⑩): a slot's stadium
        // has four, all tangent points, and neither a circle nor a whole chord — the
        // supply the mixed road above offers ([`combinatorics::branch_coords_rat`]).
        for e in ra {
            let Some(p) = combinatorics::branch_coords_rat(jd, cyls, e.node) else {
                continue;
            };
            if let Some(hit) = rational_point_in_ring(jd, cyls, wc, &p, rb)? {
                return Ok(Some(hit));
            }
        }
        // ★ And a rational point **inside an edge** whose ends are two solves' roots
        // ([`edge_interior_rat`]) — the cell between a bore's and a fillet's rulings.
        for e in ra {
            let Some(p) = edge_interior_rat(jd, cyls, e) else {
                continue;
            };
            if let Some(hit) = rational_point_in_ring(jd, cyls, wc, &p, rb)? {
                return Ok(Some(hit));
            }
        }
        // ★ And when neither a circle nor a chord names one, say **that** —
        // `ring_in_ring` below would report an exhausted probe list, which is a different
        // fact and one that never happened here.
        return Err(reject(RejectReason::RingHasNoWitness));
    }
    combinatorics::ring_in_ring(jd, wc, &probes, rb).map(Some)
}

/// **Whether the disk cylinder `a` cuts on class `wc` lies inside the disk `b` cuts there** — one
/// rational inequality on the plane: `r_a < r_b` and `(r_b − r_a)² > |c_b − c_a|²`. A disk lies
/// inside another iff its rim does, and the rims of two classes never meet (the gate proves the
/// faces clear, or the solid is valid), so the two centres and radii decide — concentric or not.
///
/// ★ **Written once because it is asked from two directions** (cell ⑩): the arrangement's nesting
/// ([`cell_in_cell`], a pin's trace under a boss's cap) and the coplanar merge (a region whose
/// outer bound is a circle asking which circle holes it owns).
pub(crate) fn disk_in_disk(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    a: &nacre_topo::CylinderDef,
    b: &nacre_topo::CylinderDef,
) -> Result<bool, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let pa = circle_centre_rat(jd, wc, a).ok_or_else(undecided)?;
    let pb = circle_centre_rat(jd, wc, b).ok_or_else(undecided)?;
    let (ra, rb) = (a.radius(), b.radius());
    if ra >= rb {
        return Ok(false);
    }
    let dr = rb.checked_sub(ra).ok_or_else(undecided)?;
    let mut dist2 = nacre_scalar::Rat::from_int(0);
    for k in 0..3 {
        let d = pb[k].checked_sub(pa[k]).ok_or_else(undecided)?;
        dist2 = dist2
            .checked_add(d.checked_mul(d).ok_or_else(undecided)?)
            .ok_or_else(undecided)?;
    }
    let dr2 = dr.checked_mul(dr).ok_or_else(undecided)?;
    Ok(dr2 > dist2)
}

/// **Which cell a question is about** — a loop of the class's arrangement, or a whole disk.
///
/// The two roads that ask this question hold their cells differently (the arrangement by index
/// into its walk's rings and circles, the coplanar merge by `Bound`), so the engine takes neither:
/// it takes the shape itself.
#[derive(Clone, Copy)]
pub(crate) enum Cell<'a> {
    Ring(&'a [combinatorics::RingEdge]),
    Disk(&'a nacre_topo::CylinderDef),
}

/// **Where a witness is, and in what vocabulary** — the 2D sibling of [`combinatorics::Probe`],
/// which does the same job one dimension up (a component's witness carries a ray direction too,
/// which is why the two are not one type).
///
/// A `Named` witness needs no coordinates, and that is not a convenience: a **rotated** class has
/// no rational description, so its points have no rational coordinates at all — the ray road is
/// the only road there, and it runs on names.
#[derive(Clone, Copy, Debug)]
enum Where {
    Named([usize; 3]),
    Coord([nacre_scalar::Rat; 3]),
}

/// **A point that can decide a nesting question — and how much it decides.**
///
/// ★★★★★ **The distinction is the whole reason this is a type** (cell 13). The two loops do not
/// cross, so a point of `a`'s **boundary** settles the question outright: `∂a` is connected and
/// misses `∂b`, so all of `∂a` is on the side that point is on. A point of `a`'s **interior**
/// settles something weaker — «this point is in `b`» — which is *not* the same claim, because
/// `b` may sit inside `a` and then `a`'s interior points are in `b` while `a` is not. The old
/// dispatch carried that as a hand-written second clause on one arm (`circle_center_in_ring` and
/// then «and the ring is not inside the circle»), and the merge road, which had copied the arm,
/// had left the second clause out.
enum Witness {
    On(Where),
    In(Where),
    /// A witness this cell *has* but whose value could not be formed exactly — checked-`Rat`
    /// overflow, a class with no rational description. Not the same as having none.
    Unformed,
}

/// What one witness said about the target.
enum Said {
    In,
    Out,
    /// This witness cannot rank itself against the target — it sits on the target's boundary, or
    /// every ray from it grazes. **The next witness is the remedy**, which is why this is not an
    /// error.
    Abstain,
    /// The value needed to ask could not be formed. The next witness may still answer, but if
    /// none does, the honest name is the arithmetic one and not "no clear ray".
    Unformed,
}

/// **Every witness `cell` can offer, in one order, once** (cell 13).
///
/// Boundary witnesses first and in this order: the corners' three-plane **names** (the ray road,
/// and the only road a rotated class has), then the corners' rational coordinates, then the
/// rational **branch** corners (a fillet tangency — cell 10), then a whole chord's **midpoint**,
/// then a rational point **inside** an edge whose ends are two solves' roots. An interior witness
/// last, and only where there is one: a disk's centre, and a ring that *is* a circle.
///
/// ★ **Lazy on purpose.** The corpus asks this question 33,781 times in one census pass, and every
/// kind after the first is a `branch_meet` solve. The name that decides is almost always the first
/// one offered (measured: 23,581 of 24,449 ring questions), so nothing after it should be built.
fn witnesses<'a>(
    jd: &'a Judge<'a, WorkingPlane>,
    cyls: &'a [crate::planes::WorkingCyl],
    wc: usize,
    cell: Cell<'a>,
) -> Box<dyn Iterator<Item = Witness> + 'a> {
    match cell {
        Cell::Disk(def) => Box::new(std::iter::once(match circle_centre_rat(jd, wc, def) {
            Some(c) => Witness::In(Where::Coord(c)),
            None => Witness::Unformed,
        })),
        Cell::Ring(r) => {
            let named = r.iter().filter_map(|e| {
                combinatorics::three_plane_name(e.node).map(|t| Witness::On(Where::Named(t)))
            });
            let coords = r.iter().flat_map(move |e| {
                combinatorics::node_coords_rat(jd, e.node)
                    .or_else(|| combinatorics::branch_coords_rat(jd, cyls, e.node))
                    .into_iter()
                    .chain(chord_midpoint_rat(jd, cyls, e))
                    .chain(edge_interior_rat(jd, cyls, e))
                    .map(|p| Witness::On(Where::Coord(p)))
            });
            // A ring that is a whole circle has no boundary point anyone can name exactly unless
            // its corners give one; its centre always answers, as an interior witness.
            let centre = ring_own_circle(r)
                .map(|def| match circle_centre_rat(jd, wc, def) {
                    Some(c) => Witness::In(Where::Coord(c)),
                    None => Witness::Unformed,
                })
                .into_iter();
            Box::new(named.chain(coords).chain(centre))
        }
    }
}

/// **What this witness says about that cell** — the one place a point is ranked against a target.
fn ask(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    w: &Where,
    target: Cell<'_>,
) -> Said {
    let radial = |p: &[nacre_scalar::Rat; 3], def: &nacre_topo::CylinderDef| {
        match nacre_scalar::cylinder_radial_side(p, &def.origin(), &def.dir(), def.radius()) {
            nacre_scalar::Orient::Negative => Said::In,
            nacre_scalar::Orient::Positive => Said::Out,
            // On the rim. The ring road has always read this as "this witness says nothing, the
            // next may"; the disk road used to reject here, and the two are one road now.
            nacre_scalar::Orient::Zero => Said::Abstain,
        }
    };
    match (w, target) {
        (Where::Named(t), Cell::Ring(rb)) => match combinatorics::point_in_ring(jd, wc, *t, rb) {
            Ok(true) => Said::In,
            Ok(false) => Said::Out,
            // ★ Every error here is this probe's abstention, and the retry over the next probe is
            // the remedy — which is exactly what `ring_in_ring` did with the same errors before
            // this loop replaced it (its swallowed-error ledger moved here with it).
            Err(e) => {
                #[cfg(test)]
                if !matches!(
                    e,
                    BoolError::Rejected {
                        reason: RejectReason::NoClearRay,
                        ..
                    }
                ) {
                    *combinatorics::swallowed_probe::COUNT
                        .lock()
                        .expect("the probe's lock is never held across a panic") += 1;
                }
                #[cfg(not(test))]
                let _ = e;
                Said::Abstain
            }
        },
        (Where::Coord(p), Cell::Ring(rb)) => match rational_point_in_ring(jd, cyls, wc, p, rb) {
            Ok(Some(true)) => Said::In,
            Ok(Some(false)) => Said::Out,
            Ok(None) => Said::Abstain,
            // A **value** that could not be formed — the class has no rational description. Not
            // an abstention: no coordinate witness will do better on this class.
            Err(_) => Said::Unformed,
        },
        (Where::Named(t), Cell::Disk(def)) => {
            match combinatorics::node_coords_rat(
                jd,
                NodeId::three_planes(combinatorics::Canon3::three(*t)),
            ) {
                Some(p) => radial(&p, def),
                None => Said::Unformed,
            }
        }
        (Where::Coord(p), Cell::Disk(def)) => radial(p, def),
    }
}

/// **Is cell `a` inside cell `b`?** — the one engine: `a`'s witnesses, in one order, against `b`,
/// until one decides.
///
/// The precondition is that the two loops do not **cross**; each road carries its own argument for
/// why (the arrangement's cells are faces of one class's DCEL; the merge's are faces of one valid
/// solid on one plane). Under it, a boundary witness of `a` decides outright and an interior one
/// needs the converse — see [`Witness`].
///
/// The three refusals are one rule, not three sites: a value that could not be formed is
/// [`RejectReason::WitnessNotRational`], an empty offer is [`RejectReason::RingHasNoWitness`], and
/// an offer every one of whose members abstained is [`RejectReason::NoClearRay`].
pub(crate) fn cell_inside(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    a: Cell<'_>,
    b: Cell<'_>,
) -> Result<bool, BoolError> {
    inside(jd, cyls, wc, a, b, true)
}

fn inside(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    a: Cell<'_>,
    b: Cell<'_>,
    // ★ **The converse is asked one level deep and no further.** An interior witness of `a` that
    // lands in `b` leaves the question open, and `b`'s **boundary** witnesses settle it; if `b`
    // has none, the honest answer is the refusal below — not a second descent, which two cells
    // that each offer only an interior point (two rings that are both whole circles) would turn
    // into an infinite one.
    may_ask_interior: bool,
) -> Result<bool, BoolError> {
    // ★ S1a preserves today's two road-choices exactly, so that this rung changes nothing but the
    // **supply** of the arm that was one witness wide. Both go in S1b, where the list is one list:
    //   · a **mixed** target has no ray road at all today (every ray answers `Unnameable` at the
    //     first branch corner), so its names are not offered;
    //   · a plain ring target that named any probe answers on the ray road **alone**.
    let target_mixed = matches!(b, Cell::Ring(rb) if combinatorics::ring_is_mixed(rb));
    let stop_after_names = matches!(b, Cell::Ring(_)) && !target_mixed;
    let (mut asked, mut unformed, mut saw_named) = (false, false, false);
    for w in witnesses(jd, cyls, wc, a) {
        let (p, interior) = match w {
            Witness::Unformed => {
                unformed = true;
                continue;
            }
            Witness::On(p) => (p, false),
            Witness::In(p) => {
                if !may_ask_interior {
                    continue;
                }
                (p, true)
            }
        };
        if target_mixed && matches!(p, Where::Named(_)) {
            continue;
        }
        if stop_after_names && saw_named && !matches!(p, Where::Named(_)) {
            break;
        }
        saw_named |= matches!(p, Where::Named(_));
        match ask(jd, cyls, wc, &p, b) {
            Said::Unformed => unformed = true,
            Said::Abstain => asked = true,
            // A point of `a` outside `b`: whether it is a boundary point or an interior one, `a`
            // is not inside `b`.
            Said::Out => return Ok(false),
            Said::In if !interior => return Ok(true),
            Said::In => return Ok(!inside(jd, cyls, wc, b, a, false)?),
        }
    }
    Err(reject(if unformed {
        RejectReason::WitnessNotRational
    } else if asked {
        RejectReason::NoClearRay
    } else {
        RejectReason::RingHasNoWitness
    }))
}

/// **Is cell `a`'s loop inside cell `b`'s?** — the one dispatch, four arms by carrier kind.
///
/// `Ok(None)` is *not comparable*, and it covers two shapes that both mean "these are neighbours,
/// not nested": two circles (see below), and two polygons that **share a node** — a shared node is
/// a split point, so the loops touch rather than one wrapping the other (that also excludes a
/// contour's own `+1` partner, which carries the same ring).
///
/// ★★ **Written once because it is asked from two directions.** [`nest_cells`] asks it of
/// (contour, `+1` cell) to find hosts, and [`innermost_host`] asks it of (host, host) to order
/// them. The four arms are the same question either way; two spellings of it would be free to
/// drift, which is this repository's most-repeated defect.
///
/// The arms:
/// - **two circles** — [`disk_in_disk`] on the radii and the centre distance (cell ⑩ S1: two disks
///   of one class *do* nest — a pin fused onto a boss — the gate keeps only different classes apart);
/// - **circle in polygon** — the centre (`axis ∩ wc`, rational) inside the ring, **and the ring not
///   inside the circle**. ★★ That second clause is not belt-and-braces: with disjoint loops the
///   centre test alone says "inside" for *both* nestings when the polygon happens to straddle the
///   centre — a boss standing over a bore, footprint `[7,9]×[9,11]` around the axis `(8,10)`, put
///   the **circle** inside the **square**. `circle_center_in_ring`'s own doc says one point decides
///   *"the loops cannot cross"*, and cross they do not; what it does not settle is **which way
///   round**. The other direction does, so the pair is the predicate;
/// - **polygon in circle** — one node's radial side ([`node_in_circle`]), decisive on its own:
///   disjoint loops put every node on one side;
/// - **polygon in polygon** — a ray from each of `a`'s nodes until one is clear
///   ([`combinatorics::ring_in_ring`]); when `b` is a **mixed** ring (branch corners, arc steps —
///   a plate's section bitten by a boss, once the scan names crossings on rulings), the ray is
///   [`combinatorics::point_in_mixed_ring`]'s from each of `a`'s rational nodes, the predicate the circle arm
///   already asks of a centre. ★ The old road *always* exhausts on a mixed `b` (every ray answers
///   `Unnameable` at the first branch corner), so this arm activates on exactly the population it
///   refused — every other row stays bit-identical.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cell_in_cell(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    circles: &[MergedCircle],
    circle_ix: &[Option<usize>],
    a: usize,
    b: usize,
) -> Result<Option<bool>, BoolError> {
    #[cfg(test)]
    if nesting_probe::on() {
        let inv = |i: usize| -> (usize, usize, usize, usize, usize, bool) {
            if circle_ix[i].is_some() {
                return (0, 0, 0, 0, 0, false);
            }
            let r = &rings[i];
            (
                r.iter()
                    .filter(|e| combinatorics::three_plane_name(e.node).is_some())
                    .count(),
                r.iter()
                    .filter(|e| combinatorics::node_coords_rat(jd, e.node).is_some())
                    .count(),
                r.iter()
                    .filter(|e| combinatorics::branch_coords_rat(jd, cyls, e.node).is_some())
                    .count(),
                r.iter()
                    .filter(|e| chord_midpoint_rat(jd, cyls, e).is_some())
                    .count(),
                r.iter()
                    .filter(|e| edge_interior_rat(jd, cyls, e).is_some())
                    .count(),
                ring_own_circle(r).is_some(),
            )
        };
        let (named, coords, branch, chord, edge, circle) = inv(a);
        // Where both roads can answer, ask both — the premise under test.
        let roads = (|| {
            if circle_ix[a].is_some() || circle_ix[b].is_some() {
                return None;
            }
            let (ra, rb) = (&rings[a], &rings[b]);
            // Through the shared door, as production asks it — the retry over the probes lives
            // in `ring_in_ring`, and comparing anything else would compare a road nobody walks.
            let probes = combinatorics::three_plane_probes(ra.iter().map(|e| e.node));
            let by_name = combinatorics::ring_in_ring(jd, wc, &probes, rb).ok()?;
            let by_coord = ra
                .iter()
                .find_map(|e| {
                    combinatorics::node_coords_rat(jd, e.node)
                        .or_else(|| combinatorics::branch_coords_rat(jd, cyls, e.node))
                        .or_else(|| chord_midpoint_rat(jd, cyls, e))
                        .or_else(|| edge_interior_rat(jd, cyls, e))
                })
                .and_then(|p| rational_point_in_ring(jd, cyls, wc, &p, rb).ok().flatten())?;
            Some((by_name, by_coord))
        })();
        nesting_probe::push(nesting_probe::Row {
            a_disk: circle_ix[a].is_some(),
            b_disk: circle_ix[b].is_some(),
            b_mixed: circle_ix[b].is_none() && combinatorics::ring_is_mixed(&rings[b]),
            named,
            coords,
            branch,
            chord,
            edge,
            circle,
            roads,
        });
    }
    // ★ **Adjacency stays a ring–ring rule, where it has always been.** A shared node is a split
    // point, so the two loops touch rather than one wrapping the other (that also excludes a
    // contour's own `+1` partner, which carries the same ring). It is *not* asked of the disk arms
    // today, and hoisting it into the engine would answer «not comparable» where they answer — and
    // would change the merge road, whose rings do share nodes. The asymmetry is left where it is
    // and written down; whether it is right is a question for its own cell.
    if circle_ix[a].is_none() && circle_ix[b].is_none() {
        let (ra, rb) = (&rings[a], &rings[b]);
        if ra.iter().any(|e| rb.iter().any(|f| f.node == e.node)) {
            return Ok(None);
        }
    }
    // ★ **Two disks are not a witness question at all** (cell 10): a disk lies inside another iff
    // its rim does, and the rims of two classes never meet, so the radii and the centre distance
    // decide it exactly. ★ This arm used to answer `None`, which left two disk cells unnested and
    // *silently* kept both operands' caps: a pin stacked on a boss fused into **two** untouched
    // bodies. The old refusal was masking a gap.
    if let (Some(ca), Some(cb)) = (circle_ix[a], circle_ix[b]) {
        return disk_in_disk(jd, wc, &circles[ca].def, &circles[cb].def).map(Some);
    }
    let cell = |i: usize| match circle_ix[i] {
        Some(c) => Cell::Disk(&circles[c].def),
        None => Cell::Ring(&rings[i]),
    };
    cell_inside(jd, cyls, wc, cell(a), cell(b)).map(Some)
}

/// **What every nesting question was asked with** (cell 13, test-only).
///
/// The defect this cell fixes was invisible for one reason: the corpus never put a **ring with no
/// three-plane corner** against a **disk**, so the arm whose witness supply was truncated to that
/// one kind never fired. A reject count could not have seen it. So the instrument measures the
/// **population**: for each question, what the source cell had to offer and which road could answer.
///
/// Off by default and switched on by the audit test — the road-agreement column runs *both* roads
/// where both are available, which is work no ordinary test should pay for.
#[cfg(test)]
pub(crate) mod nesting_probe {
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{Mutex, OnceLock};

    /// One question, as it arrived.
    #[derive(Clone, Copy, Debug)]
    pub(crate) struct Row {
        pub a_disk: bool,
        pub b_disk: bool,
        pub b_mixed: bool,
        /// `a`'s witnesses by kind — the five supplies today's four spellings draw from.
        pub named: usize,
        pub coords: usize,
        pub branch: usize,
        pub chord: usize,
        pub edge: usize,
        pub circle: bool,
        /// Where **both** roads could answer, what each said — the direct measurement of this
        /// cell's premise, that any witness gives the same answer.
        pub roads: Option<(bool, bool)>,
    }

    static ENABLED: AtomicBool = AtomicBool::new(false);

    pub(crate) fn on() -> bool {
        ENABLED.load(Ordering::Relaxed)
    }
    pub(crate) fn enable() {
        ENABLED.store(true, Ordering::Relaxed);
    }
    pub(crate) fn disable() {
        ENABLED.store(false, Ordering::Relaxed);
    }

    fn rows() -> &'static Mutex<Vec<Row>> {
        static ROWS: OnceLock<Mutex<Vec<Row>>> = OnceLock::new();
        ROWS.get_or_init(|| Mutex::new(Vec::new()))
    }
    pub(crate) fn push(r: Row) {
        rows()
            .lock()
            .expect("the probe's lock is never held across a panic")
            .push(r);
    }
    /// Take everything recorded so far, leaving the list empty.
    pub(crate) fn take() -> Vec<Row> {
        std::mem::take(
            &mut *rows()
                .lock()
                .expect("the probe's lock is never held across a panic"),
        )
    }
}
