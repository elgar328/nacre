use super::*;
/// Why an order around a circle could not be formed — two causes, because two callers turn them
/// into different words (the split into a reject, the tracer into a decline), the same split
/// the tracer's `RingFail`/`decline_of` makes on the segment side.
pub(crate) enum CircleOrderFail {
    /// A θ comparison could not be formed exactly (checked-`Rat` overflow, a missing description).
    Undecided,
    /// Two of the nodes are one point wearing two names.
    Coincident,
}

impl CircleOrderFail {
    /// The refusal a boolean gives — the one mapping every road that refuses on an order reads:
    /// an order past `Rat` is width, one point under two names is the four-plane limit
    /// ([`RejectReason::CoincidentNodes`]).
    pub(crate) fn reason(self) -> RejectReason {
        match self {
            Self::Undecided => RejectReason::WitnessNotRational,
            Self::Coincident => RejectReason::CoincidentNodes,
        }
    }
}

/// **Order nodes around one circle by θ, the seam first** — the nodes' reading of the cyclic order
/// [`nacre_exact::quad::SeamOrder::seam_first`] states once (the arcs' spans read it too).
///
/// ★★★ **A crossing on the seam is ordered, not refused — it is the cut point.** (Measured: the
/// very first fixture puts a crossing there — a boss on a plate's edge cuts its own rim exactly on
/// the seam generator, so this is the common case, not an exotic one.)
///
/// ★ At most one node can be seam-incident: two would be the same point, which the adjacency check
/// below refuses as the two-names-for-one-point it is.
///
/// ★★ **Two names for one point.** Two walls crossing the circle at one point are two *different*
/// pierce names — a `dedup` cannot see it, and the θ sort puts them adjacent — so an arc of zero
/// length would follow. The segment side asks this question already; asking it here too is what
/// keeps the two sides from disagreeing about what "one point" means.
///
/// Returns the permutation of `nodes` in θ order. `nodes` must already be deduped by name.
pub(crate) fn circular_order(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    nodes: &[NodeId],
) -> Result<Vec<usize>, CircleOrderFail> {
    use nacre_exact::quad::{SeamOrder, circular_order_about_seam};
    let meets = nodes
        .iter()
        .map(|&n| {
            pierce_meet(jd, cyl, def, n)
                .and_then(PierceMeet::into_exact)
                .ok_or(CircleOrderFail::Undecided)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let cmp = |i: usize, j: usize| {
        circular_order_about_seam(
            &def.origin(),
            &def.dir(),
            &def.ref_dir(),
            (&meets[i].0, &meets[i].1),
            (&meets[j].0, &meets[j].1),
        )
    };
    // ★ Counted before the sort, so two seam nodes are `Coincident` ahead of any order the sort
    // could not form (`Undecided`).
    let mut seam_nodes = 0usize;
    for i in 0..meets.len() {
        match cmp(i, i) {
            Some(SeamOrder::SeamIncident { first: true, .. }) => seam_nodes += 1,
            Some(_) => {}
            None => return Err(CircleOrderFail::Undecided),
        }
    }
    if seam_nodes > 1 {
        return Err(CircleOrderFail::Coincident);
    }
    let mut order: Vec<usize> = (0..meets.len()).collect();
    try_sort_by(&mut order, |&i, &j| {
        cmp(i, j)
            .map(SeamOrder::seam_first)
            .ok_or(CircleOrderFail::Undecided)
    })?;
    for w in order.windows(2) {
        match cmp(w[0], w[1]) {
            Some(SeamOrder::Ordered(core::cmp::Ordering::Equal)) => {
                return Err(CircleOrderFail::Coincident);
            }
            Some(_) => {}
            // The sort may have asked this pair the other way round; a failure is no answer.
            None => return Err(CircleOrderFail::Undecided),
        }
    }
    Ok(order)
}

/// **Is `node` strictly inside the counter-clockwise arc from `lo` to `hi`?** — cyclic order about
/// the axis, by name ([`circular_order`]). For two ends that differ, `None` when the order cannot
/// be formed or two of the three are one point wearing two names — the caller's population, not
/// a tie to shrug at.
///
/// ★ **`lo == hi` is the closed arc**: one node cuts the circle and the arc runs its whole turn
/// from that node back to it, so every other point of the circle is inside — and the node itself,
/// under its own name or another (`Coincident`), is not. Read as the three-point order it would be
/// two points, and refuse.
pub(crate) fn theta_between(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    lo: NodeId,
    hi: NodeId,
    node: NodeId,
) -> Option<bool> {
    let nodes = [lo, hi, node];
    // A point with no exact meet is placed by its interval box — and only such a point: a named
    // trio keeps the name road, whose «one point, two names» answer an interval has no twin for.
    if nodes.iter().any(|&n| names_a_nameless_plane(jd, n)) {
        return between_on_circle_hp(jd, k, def, lo, hi, node);
    }
    if lo == hi {
        let inside = node != lo
            && match circular_order(jd, k, def, &[lo, node]) {
                Ok(_) => true,
                Err(CircleOrderFail::Coincident) => false,
                Err(CircleOrderFail::Undecided) => return None,
            };
        debug_assert!(
            between_on_circle_hp(jd, k, def, lo, hi, node).is_none_or(|b| b == inside),
            "the interval reading of the closed arc disagrees with the circular order"
        );
        return Some(inside);
    }
    let order = circular_order(jd, k, def, &nodes).ok()?;
    let pos = |x: usize| {
        order
            .iter()
            .position(|&i| i == x)
            .expect("a permutation of the three")
    };
    let (p_lo, p_hi, p_n) = (pos(0), pos(1), pos(2));
    let inside = (p_n + 3 - p_lo) % 3 < (p_hi + 3 - p_lo) % 3;
    // ★ The interval reading's turn convention, checked against the name road wherever both answer.
    debug_assert!(
        between_on_circle_hp(jd, k, def, lo, hi, node).is_none_or(|b| b == inside),
        "the interval turn about the axis disagrees with the circular order"
    );
    Some(inside)
}

/// **Where the class line `L = wc ∩ fc` crosses an arc strictly inside it** — the pierce nodes
/// `wc ∩ fc ∩ cyl` the travelled arc `a → b` contains, in the order its travel meets them (the
/// arc runs counter-clockwise about the axis when `ccw`). One answer serves both questions the
/// line roads ask of an arc ([`ClassLine`]): the walk's *how many* (its length, where the ends
/// cannot say) and a crossing's *which* (its `nth` entry, [`Feature::Crossing`]).
///
/// The arc lies in the face's own plane `fc` (a cap's ⊥ plane), so the meet decides first:
/// - two parallel planes (`fc` against another ⊥ class — every cap meets its twin this way), a
///   line missing the circle, or a **tangent** one cross nothing. The tangent case is common: a
///   fillet's straight edge lies on `L` and its arc leaves the end tangentially, a touch with no
///   crossing beyond it;
/// - a secant line meets the circle in two roots, and each crosses the arc iff the arc strictly
///   contains it ([`theta_between`], the ruling sweep's containment);
/// - overflow, or a meet this shape cannot produce, is `None`.
///
/// ★★ **A root at the arc's end on `L` is that end, not a crossing.** With `on_end` — the end
/// the walk found on `L` — one root *is* that end, under the same name or another (a fillet's
/// tangent corner is named by its own two planes, the root by `{wc, fc}`). That is asked of the two
/// points alone: equal names, or the circle's order answering "one point" ([`circular_order`]'s
/// `Coincident`). The containment is then asked of three **distinct** points, so its `None` is a
/// real failure — never a coincidence mistaken for one.
///
/// ★ The same canonical name the lateral's ruling sweep gives this point as a station
/// (`crossing_on_ruling(fc, wc, …)` — [`NodeId::pierce`] folds the pair and root together), so
/// the arrangement's `merge_coincident` reads one point, not two.
#[allow(clippy::too_many_arguments)]
pub(crate) fn arc_crossings(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    fc: usize,
    cyl: usize,
    ccw: bool,
    [a, b]: [NodeId; 2],
    on_end: Option<NodeId>,
) -> Option<Vec<NodeId>> {
    use nacre_exact::quad::CylinderMeet;
    // A plane with no world statement has no exact meet to read; its crossings are read off
    // interval boxes by the same contract ([`arc_crossings_hp`]).
    if jd.planes[wc].world.is_none() || jd.planes[fc].world.is_none() {
        return arc_crossings_hp(jd, def, wc, fc, cyl, ccw, [a, b], on_end);
    }
    let w = class_coeffs_rat(jd, wc)?;
    let v = class_coeffs_rat(jd, fc)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    match nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2)? {
        CylinderMeet::ParallelPlanes | CylinderMeet::Miss(_) | CylinderMeet::Tangent { .. } => {
            return Some(Vec::new());
        }
        CylinderMeet::Pair { .. } => {}
        CylinderMeet::CoincidentPlanes
        | CylinderMeet::OnRuling(_)
        | CylinderMeet::AxisParallelMiss(_) => return None,
    }
    let (lo, hi) = if ccw { (a, b) } else { (b, a) };
    let mut inside: Vec<NodeId> = Vec::with_capacity(2);
    for root in [nacre_topo::QuadRoot::Lo, nacre_topo::QuadRoot::Hi] {
        let node = NodeId::pierce(wc, fc, cyl, root);
        if let Some(e) = on_end {
            if node == e {
                continue;
            }
            match circular_order(jd, cyl, def, &[e, node]) {
                Err(CircleOrderFail::Coincident) => continue,
                Err(CircleOrderFail::Undecided) => return None,
                Ok(_) => {}
            }
        }
        if theta_between(jd, cyl, def, lo, hi, node)? {
            inside.push(node);
        }
    }
    if let [r, q] = inside[..] {
        // Travelled from `a`: counter-clockwise meets first the root counter-clockwise between
        // `a` and the other, clockwise the one counter-clockwise between the other and `a`.
        let r_first = if ccw {
            theta_between(jd, cyl, def, a, q, r)
        } else {
            theta_between(jd, cyl, def, q, a, r)
        }?;
        if !r_first {
            inside.swap(0, 1);
        }
    }
    // ★ The interval road's contract (which roots, in what travel order), checked against the
    // exact one wherever both answer — the road a class with no world statement takes.
    debug_assert!(
        arc_crossings_hp(jd, def, wc, fc, cyl, ccw, [a, b], on_end).is_none_or(|hp| hp == inside),
        "the interval arc crossings disagree with the exact ones"
    );
    Some(inside)
}

/// [`arc_crossings`] behind the one shortcut the gate's record gives: an arc of a cylinder the gate
/// proved clear of the class `wc` lies wholly on one side of it, so it crosses nothing — and a class
/// with no world statement is not asked. What [`ClassLine`] counts and names an arc's crossings by.
#[allow(clippy::too_many_arguments)]
pub(crate) fn arc_crossings_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    fc: usize,
    cyl: usize,
    ccw: bool,
    ends: [NodeId; 2],
    on_end: Option<NodeId>,
) -> Option<Vec<NodeId>> {
    let wcy = cyls.get(cyl)?;
    if wcy.beside(wc).is_some() {
        return Some(Vec::new());
    }
    arc_crossings(jd, &wcy.def, wc, fc, cyl, ccw, ends, on_end)
}

/// **The planar scan's crossing on a ruling** — the point where the class line `L = wc ∩ fc`
/// leaves the face across an edge riding a cylinder's ruling, named as the pierce node
/// `wc ∩ fc ∩ cyl` at the root that *is* this ruling.
///
/// A ruling edge lies in the face's own plane `fc` (a plane holding a ruling runs through the
/// axis — or is tangent, which the gate passes but does **not** record, so no ruling of it ever
/// reaches here), so the pair `{wc, fc}` cuts the cylinder in two
/// points, one on each of `fc`'s two rulings, and `(cyl, side)` — the identity
/// [`Wall::Ruling`] carries, measured by
/// [`ruling_side`] against its `plane` when the ring was named (on a plane
/// face's ring, `fc` itself; on a lateral's, the plane across the edge, which the caller passes
/// as `fc`) — says which. The same predicate asked of each root picks it. `wc` must be ⊥
/// to the axis for the class to cross a ruling in a point at all (∥ contains it; a tilt is
/// refused at the gate).
///
/// ★ Solved in the caller's order `(wc, fc)`, as `rulings_on_class` and `chord_stations`
/// spell it — [`NodeId::pierce`] canonicalizes the pair and the root together, so this is the one
/// name every road gives the point. The lateral face names the same point when its hole ring
/// crosses the class (`cycle_on_class`'s crossing arm restates the hole corner's own root to
/// `wc` by the axis senses of the two ⊥ classes); the two spellings agree because the roots of
/// `{⊥, wall}` are ordered along `ε·k·(m̂ × n_wall)`, so «which root» and «which side of `wall`»
/// are the same question — a derivation the exact-volume rows of the crossing census check.
///
/// `Err` says why there is no such point, in two causes ([`RulingNameFail`]): a description past
/// `Rat`, or classes that name no crossing of that wall's ruling — a wall not parallel to the
/// axis, no pair of roots (or no double root for a tangent wall), or both roots on one side (the
/// two rulings of a wall within the radius — through the axis or offset from it — are symmetric
/// about the plane through the axis with normal `m × n̂`, so `ruling_side` tells them apart).
pub(crate) fn crossing_on_ruling(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    fc: usize,
    cyl: usize,
    side: i8,
) -> Result<NodeId, RulingNameFail> {
    use RulingNameFail::{NoCrossing, Undecided, Width};
    use nacre_exact::quad::CylinderMeet;
    // ★ **A wall with no world statement picks its root by the same predicate, read off boxes.**
    // The class `wc` crossing the ruling is ⊥ the axis and named; the wall sections the cylinder,
    // so its two rulings are two distinct lines and each root is on one of them
    // ([`PierceMeet::ruling_side`]).
    if jd.planes[fc].world.is_none() {
        let w = class_coeffs_rat(jd, wc).ok_or(NoCrossing)?;
        if side == 0 || !nacre_exact::parallel_rat(&[w[0], w[1], w[2]], &def.dir()) {
            return Err(NoCrossing);
        }
        let mut found = None;
        for root in [nacre_topo::QuadRoot::Lo, nacre_topo::QuadRoot::Hi] {
            let node = NodeId::pierce(wc, fc, cyl, root);
            let meet = pierce_meet(jd, cyl, def, node).ok_or(Undecided)?;
            if meet.ruling_side(jd, fc, def).ok_or(Undecided)? == side {
                if found.is_some() {
                    return Err(NoCrossing);
                }
                found = Some(node);
            }
        }
        return found.ok_or(NoCrossing);
    }
    let w = class_coeffs_rat(jd, wc).ok_or(NoCrossing)?;
    let v = class_coeffs_rat(jd, fc).ok_or(NoCrossing)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    if !nacre_exact::parallel_rat(&[w[0], w[1], w[2]], &m) {
        return Err(NoCrossing);
    }
    // ★ A **tangent** wall (`side == 0`) has one ruling and one root — `Double`.
    let meet = nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2).ok_or(Width)?;
    if side == 0 {
        return match meet {
            CylinderMeet::Tangent { .. } => {
                Ok(NodeId::pierce(wc, fc, cyl, nacre_topo::QuadRoot::Double))
            }
            _ => Err(NoCrossing),
        };
    }
    let CylinderMeet::Pair { line, s } = meet else {
        return Err(NoCrossing);
    };
    // A root whose side overflowed is not on `side` as far as this can tell; the other root may
    // still be, so the overflow is the answer only when no root is.
    let (mut found, mut overflow) = (None, false);
    for (root, sv) in [
        (nacre_topo::QuadRoot::Lo, &s[0]),
        (nacre_topo::QuadRoot::Hi, &s[1]),
    ] {
        match ruling_side_signed(&v, def, (&line, sv)) {
            None => overflow = true,
            Some(sd) if sd == side => {
                if found.is_some() {
                    return Err(NoCrossing); // both roots on one side: not a pair of rulings
                }
                found = Some(root);
            }
            Some(_) => {}
        }
    }
    match found {
        Some(root) => Ok(NodeId::pierce(wc, fc, cyl, root)),
        None if overflow => Err(Width),
        None => Err(NoCrossing),
    }
}

/// Why [`crossing_on_ruling`] names no point — two causes, because two readers turn them into
/// different words (the tracer into its decline, the chart into a refusal), the same split
/// [`CircleOrderFail`] makes for an order around a circle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RulingNameFail {
    /// A description past `Rat`: the meet of the two planes and the cylinder, or a root's side
    /// of the wall, overflowed.
    Width,
    /// The classes name no crossing of that wall's ruling — or a class has no narrow world
    /// statement, which the cylinder gate demands of every class it does not section.
    NoCrossing,
    /// A wall with no world statement (one that sections the cylinder): the roots' boxes could not
    /// say which ruling each is on.
    Undecided,
}

impl RulingNameFail {
    /// The refusal where the caller knows the crossing exists — a station the split or the
    /// chart already placed: past `Rat` is width, no crossing is the stages disagreeing.
    pub(crate) fn reason(self) -> crate::RejectReason {
        match self {
            Self::Width => crate::RejectReason::WitnessNotRational,
            Self::NoCrossing => crate::RejectReason::CylinderStagesDisagree,
            Self::Undecided => crate::RejectReason::CylinderGateUndecided,
        }
    }
}
