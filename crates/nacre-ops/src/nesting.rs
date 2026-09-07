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
/// dispatch calls a cell a `Disk` only where its half-edge is literally a `Circle`; a circle a wall
/// has **split into arcs** arrives as a ring instead, whose corners are branch points — so the
/// three-plane names come out
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
/// below used to sit in the disk-against-a-ring arm, where one witness is all there is, so a point
/// that landed *on* the ring could be folded into the same refusal as a class with no rational
/// description. A caller with **several** witnesses must not read them the same way: a
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
    // A disk's interior witness, wherever the disk came from — the cell itself, or a ring that
    // turned out to be a whole circle. Written once because this cell is about supplies that were
    // written twice.
    let centre_of = |def: &nacre_topo::CylinderDef| match circle_centre_rat(jd, wc, def) {
        Some(c) => Witness::In(Where::Coord(c)),
        None => Witness::Unformed,
    };
    match cell {
        Cell::Disk(def) => Box::new(std::iter::once(centre_of(def))),
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
            let centre = ring_own_circle(r).map(centre_of).into_iter();
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
    // ★★★★★ **The list is one list, and it is walked to the end.** Two road-choices used to cut it
    // short, and both were statements about a *road* rather than about the question:
    //   · a **mixed** target was never offered the ray road, because every ray answers
    //     `Unnameable` at the first branch corner — true, and now said by the ray itself, which
    //     abstains and hands on to the next witness at no cost but its own;
    //   · a plain ring target that named any probe answered on names **alone**, so a ring whose
    //     every ray grazed was refused with coordinates still in hand.
    // Neither is a fact about *whether `a` is inside `b`*, and the premise that lets them go is
    // measured, not assumed: where both roads can answer they agree
    // (`the_two_roads_never_disagree`).
    let (mut asked, mut unformed) = (false, false);
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

/// **The arrangement's adapter onto [`cell_inside`]** — its cells arrive as indices into the
/// walk's rings and circles, and two questions are settled before a witness is asked.
///
/// `Ok(None)` is *not comparable*: two polygons that **share a node** touch rather than nest (a
/// shared node is a split point, and that also excludes a contour's own `+1` partner, which
/// carries the same ring). ★ It is a **ring-ring** rule and stays one: the disk arms have never
/// asked it, and the coplanar merge — the engine's other road — works with rings that do share
/// nodes. Whether the asymmetry is right is a question for its own cell, not a thing to change
/// while unifying.
///
/// ★ **The engine's precondition holds here structurally**: these cells are faces of one class's
/// DCEL, and the split gave every crossing a vertex — so two of their loops cannot cross, they can
/// only nest or share nodes.
///
/// ★★ **Asked from two directions here and a third elsewhere.** [`crate::arrangement::nest_cells`]
/// asks it of (contour, `+1` cell) to find hosts and `innermost_host` of (host, host) to order
/// them; the coplanar merge asks the same question of its own bounds. That third one used to hold
/// a copy of this dispatch, and the copy is what cell 13 removed.
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
            // ★ Through the engine's own per-witness door, so the comparison is between the two
            // **roads** and not between two hand-written loops: the first name that decides
            // against the first coordinate that decides.
            let decided = |w: Where| match ask(jd, cyls, wc, &w, Cell::Ring(rb)) {
                Said::In => Some(true),
                Said::Out => Some(false),
                Said::Abstain | Said::Unformed => None,
            };
            let by_name = ra
                .iter()
                .filter_map(|e| combinatorics::three_plane_name(e.node))
                .find_map(|t| decided(Where::Named(t)))?;
            let by_coord = ra
                .iter()
                .flat_map(|e| {
                    combinatorics::node_coords_rat(jd, e.node)
                        .or_else(|| combinatorics::branch_coords_rat(jd, cyls, e.node))
                        .into_iter()
                        .chain(chord_midpoint_rat(jd, cyls, e))
                        .chain(edge_interior_rat(jd, cyls, e))
                })
                .find_map(|p| decided(Where::Coord(p)))?;
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
