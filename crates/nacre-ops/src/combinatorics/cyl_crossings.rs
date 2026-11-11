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

/// A pierce node's exact point on cylinder `k` — `WitnessNotRational` where it has none.
pub(crate) fn exact_point(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Result<ExactPoint, RejectReason> {
    pierce_meet(jd, k, def, n)
        .and_then(PiercePoint::into_exact)
        .ok_or(RejectReason::WitnessNotRational)
}

/// A point of a meet line: the line and the parameter along it.
pub(crate) type ExactPoint = (nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal);

/// **Which of two lines of a cylinder's chart is higher** — the sections of plane classes `p` and
/// `q`, compared on the generator through point `at` ([`Toward::At`](nacre_exact::quad::Toward)),
/// just after or just before it, or with no generator where the two do not cross anywhere
/// ([`nacre_exact::quad::section_height_order`], the one spelling). Refused as
/// `WitnessNotRational` where a plane has no narrow world statement, and as
/// `CylinderStagesDisagree` where the question has no answer the chart should have asked: a plane
/// along the axis is no line, and two lines that cross were asked at no generator.
pub(crate) fn line_order(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    p: usize,
    q: usize,
    at: Option<&ExactPoint>,
    toward: nacre_exact::quad::Toward,
) -> Result<core::cmp::Ordering, RejectReason> {
    use nacre_exact::quad::SectionOrder;
    let (Some(wp), Some(wq)) = (jd.world_rat(p), jd.world_rat(q)) else {
        return Err(RejectReason::WitnessNotRational);
    };
    match nacre_exact::quad::section_height_order(&wp, &wq, def, at.map(|(l, s)| (l, s)), toward) {
        Some(SectionOrder::Ordered(o)) => Ok(o),
        Some(SectionOrder::Crosses) | None => Err(RejectReason::CylinderStagesDisagree),
    }
}

/// **Order points of the cylinder by θ, the seam first, points at one θ together** — the nodes'
/// reading of the cyclic order [`nacre_exact::quad::SeamOrder::seam_first`] states once (the arcs'
/// spans read it too). Returns the groups of indices into `nodes` in θ order; a group holds the
/// nodes on one generator — on one circle that is one point, on the lateral's chart one station
/// several lines pass.
///
/// **A crossing on the seam is ordered, not refused — it is the cut point.** (Measured: the
/// very first fixture puts a crossing there — a boss on a plate's edge cuts its own rim exactly on
/// the seam generator, so this is the common case, not an exotic one.) Points on the seam are one
/// θ and share the first group.
pub(crate) fn theta_groups(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    nodes: &[NodeId],
) -> Result<Vec<Vec<usize>>, CircleOrderFail> {
    use nacre_exact::quad::{SeamOrder, circular_order_about_seam};
    let meets = nodes
        .iter()
        .map(|&n| {
            pierce_meet(jd, cyl, def, n)
                .and_then(PiercePoint::into_exact)
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
        .map(SeamOrder::seam_first)
        .ok_or(CircleOrderFail::Undecided)
    };
    // Every point's own θ is asked, a point alone included.
    for i in 0..meets.len() {
        cmp(i, i)?;
    }
    let mut order: Vec<usize> = (0..meets.len()).collect();
    try_sort_by(&mut order, |&i, &j| cmp(i, j))?;
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for &i in &order {
        match groups.last_mut() {
            // The sort may have asked this pair the other way round; a failure is no answer.
            Some(g) if cmp(g[0], i)? == core::cmp::Ordering::Equal => g.push(i),
            _ => groups.push(vec![i]),
        }
    }
    Ok(groups)
}

/// **Order nodes around one circle by θ, the seam first** — [`theta_groups`] on points of one
/// circle, where one θ is one point.
///
/// **Two names for one point.** Two walls crossing the circle at one point are two *different*
/// pierce names — a `dedup` cannot see it, and the θ sort puts them together — so an arc of zero
/// length would follow. A group of two or more is that point, refused as `Coincident`, and so are
/// two nodes on the seam. This is the trace's order, where two cuts at one θ are one boundary
/// cutting one point twice (`cyl_trace`'s `CylHoleFeature`); the curve splits read the same groups
/// ([`theta_groups`]) inside the alias rounds and record a tie as one point.
///
/// Returns the permutation of `nodes` in θ order. `nodes` must already be deduped by name.
pub(crate) fn circular_order(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    nodes: &[NodeId],
) -> Result<Vec<usize>, CircleOrderFail> {
    let groups = theta_groups(jd, cyl, def, nodes)?;
    if groups.iter().any(|g| g.len() > 1) {
        return Err(CircleOrderFail::Coincident);
    }
    Ok(groups.into_iter().flatten().collect())
}

/// **Where `node` lies against the counter-clockwise arc from `lo` to `hi`** — one cyclic order of
/// the three about the axis ([`theta_groups`]), which says at once whether `node` is strictly
/// inside, outside, or one point with an end (`At(0)` for `lo`, `At(1)` for `hi`; `lo` when it is
/// both). Points with no exact meet are placed by their interval boxes, which never tie; the exact
/// points among them are still grouped exactly, so a `node` at an exact end is `At` that end.
/// `None` when the order cannot be formed, and when the two ends of an arc that is not closed are
/// one point.
///
/// **`lo == hi` is the closed arc**: one node cuts the circle and the arc runs its whole turn
/// from that node back to it, so every other point of the circle is inside.
pub(crate) fn theta_place(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    lo: NodeId,
    hi: NodeId,
    node: NodeId,
) -> Option<ThetaPlace> {
    let closed = lo == hi;
    if node == lo {
        return Some(ThetaPlace::At(0));
    }
    if node == hi {
        return Some(ThetaPlace::At(1));
    }
    let nodes: &[NodeId] = if closed { &[lo, node] } else { &[lo, hi, node] };
    let ix = nodes.len() - 1;
    let exact: Vec<usize> = (0..nodes.len())
        .filter(|&i| !names_a_nameless_plane(jd, nodes[i]))
        .collect();
    let between = |inside: bool| {
        if inside {
            ThetaPlace::Inside
        } else {
            ThetaPlace::Outside
        }
    };
    // A point with no exact meet is placed by its interval box — and only such a point: a named
    // trio keeps the name road, whose «one point, two names» answer an interval has no twin for.
    if exact.len() < nodes.len() {
        if exact.len() > 1 && exact.contains(&ix) {
            let sub: Vec<NodeId> = exact.iter().map(|&i| nodes[i]).collect();
            if let Ok(groups) = theta_groups(jd, k, def, &sub) {
                let g = groups
                    .iter()
                    .find(|g| g.contains(&(sub.len() - 1)))
                    .expect("every point is in a group");
                if let Some(&e) = g.iter().map(|&j| &exact[j]).find(|&&i| i != ix) {
                    return Some(ThetaPlace::At(e));
                }
            }
        }
        return between_on_circle_hp(jd, k, def, lo, hi, node).map(between);
    }
    let groups = theta_groups(jd, k, def, nodes).ok()?;
    let group_of = |i: usize| {
        groups
            .iter()
            .position(|g| g.contains(&i))
            .expect("every point is in a group")
    };
    let (g_lo, g_n) = (group_of(0), group_of(ix));
    let place = if g_n == g_lo {
        ThetaPlace::At(0)
    } else if closed {
        ThetaPlace::Inside
    } else {
        let g_hi = group_of(1);
        if g_n == g_hi {
            ThetaPlace::At(1)
        } else if g_hi == g_lo {
            return None;
        } else {
            between((g_n + 3 - g_lo) % 3 < (g_hi + 3 - g_lo) % 3)
        }
    };
    // The interval reading's turn convention, checked against the name road wherever both answer.
    debug_assert!(
        match place {
            ThetaPlace::Inside => between_on_circle_hp(jd, k, def, lo, hi, node) != Some(false),
            ThetaPlace::Outside => between_on_circle_hp(jd, k, def, lo, hi, node) != Some(true),
            ThetaPlace::At(_) => true,
        },
        "the interval turn about the axis disagrees with the circular order"
    );
    Some(place)
}

/// Where a point lies against an arc ([`theta_place`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ThetaPlace {
    /// Strictly inside the arc.
    Inside,
    /// Strictly outside it.
    Outside,
    /// One point with the arc's end `0` (`lo`) or `1` (`hi`).
    At(usize),
}

/// **Is `node` strictly inside the counter-clockwise arc from `lo` to `hi`?** — [`theta_place`]
/// read as a yes or no. For two ends that differ, `None` when the order cannot be formed or two of
/// the three are one point wearing two names — the caller's population, not a tie to shrug at. On
/// the closed arc (`lo == hi`) the node itself, under its own name or another, is not inside.
pub(crate) fn theta_between(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    lo: NodeId,
    hi: NodeId,
    node: NodeId,
) -> Option<bool> {
    match theta_place(jd, k, def, lo, hi, node)? {
        ThetaPlace::Inside => Some(true),
        ThetaPlace::Outside => Some(false),
        ThetaPlace::At(_) => (lo == hi).then_some(false),
    }
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
///   contains it ([`theta_between`], the containment the rulings road reads too);
/// - overflow, a `Wide` name, a meet this shape cannot produce, or — with a plane of no world
///   statement — a shape the intervals cannot say ([`Judge::line_meet`] answers both roads), is
///   `None`.
///
/// **A root at the arc's end on `L` is that end, not a crossing.** With `on_end` — the end
/// the walk found on `L` — one root *is* that end, under the same name or another (a fillet's
/// tangent corner is named by its own two planes, the root by `{wc, fc}`). That is asked of the two
/// points alone: equal names, or the closed arc from the end not holding the root
/// ([`theta_between`] — the circle's order answering "one point" on the exact road, the boxes not
/// proved apart on the interval road, which is `None`). The containment is then asked of three **distinct** points, so its `None` is a
/// real failure — never a coincidence mistaken for one.
///
/// The same canonical name the lateral's ruling sweep gives this point as a station
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
    use nacre_judge::cylinder::LineMeet;
    match jd.line_meet(wc, fc, def, false) {
        LineMeet::ParallelPlanes | LineMeet::Miss | LineMeet::Tangent(..) => {
            return Some(Vec::new());
        }
        LineMeet::Pair(_) => {}
        LineMeet::CoincidentPlanes
        | LineMeet::OnRuling(_)
        | LineMeet::AxisParallelMiss
        | LineMeet::Undecided
        | LineMeet::Wide
        | LineMeet::Overflow => return None,
    }
    let (lo, hi) = if ccw { (a, b) } else { (b, a) };
    let mut inside: Vec<NodeId> = Vec::with_capacity(2);
    for root in [QuadRoot::Lo, QuadRoot::Hi] {
        let node = NodeId::pierce(wc, fc, cyl, root);
        if let Some(e) = on_end {
            if node == e {
                continue;
            }
            // The end under another name: the closed arc from `e` holds every point of the circle
            // but `e` itself ([`theta_between`]).
            if !theta_between(jd, cyl, def, e, e, node)? {
                continue;
            }
        }
        if theta_between(jd, cyl, def, lo, hi, node)? {
            inside.push(node);
        }
    }
    if let [r, q] = inside[..]
        && !met_first(jd, cyl, def, ccw, a, r, q)?
    {
        inside.swap(0, 1);
    }
    Some(inside)
}

/// **Whether an arc travelled from its end `a` meets `r` before `q`**, two points strictly inside
/// it — counter-clockwise when `ccw`. Counter-clockwise travel meets first the point
/// counter-clockwise between `a` and the other, clockwise travel the one counter-clockwise between
/// the other and `a` ([`theta_between`]). `None` where that order cannot be formed.
pub(crate) fn met_first(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    ccw: bool,
    a: NodeId,
    r: NodeId,
    q: NodeId,
) -> Option<bool> {
    if ccw {
        theta_between(jd, cyl, def, a, q, r)
    } else {
        theta_between(jd, cyl, def, q, a, r)
    }
}

/// **Where a class oblique to a cylinder's axis crosses one of its rim circles** — the circle on
/// class `plane`, across the axis: the line `wc ∩ plane` meets it in the two roots named on the
/// cylinder (`NodeId::pierce(wc, plane, cyl, Lo | Hi)`), the names a cap's chord and the section's
/// arcs both end at. `Some(Some(_))` for the two crossings, `Some(None)` for a line that misses
/// the circle, `None` for any other meet — a touch (`planes::boundary_stand` refuses it at the
/// gate), or one that cannot be stated. The two readers take a miss differently: a chord is then
/// absent, while the section's tracer asks only of a rim it has read as crossed.
pub(crate) fn rim_crossings(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    plane: usize,
    cyl: usize,
) -> Option<Option<[NodeId; 2]>> {
    use nacre_judge::cylinder::LineMeet;
    match jd.line_meet(wc, plane, def, false) {
        LineMeet::Pair(_) => Some(Some([
            NodeId::pierce(wc, plane, cyl, QuadRoot::Lo),
            NodeId::pierce(wc, plane, cyl, QuadRoot::Hi),
        ])),
        LineMeet::Miss | LineMeet::ParallelPlanes => Some(None),
        _ => None,
    }
}

/// [`arc_crossings`] behind the gate's record: an arc of a cylinder the gate proved clear of the class
/// `wc` lies wholly on one side of it, and an arc of a cylinder `wc` touches along one ruling
/// ([`crate::planes::Silence::Grazes`]) lies on one side but for that line — so neither crosses, and
/// a class with no world statement is not asked. What [`ClassLine`] counts and names an arc's
/// crossings by.
///
/// The record answers where the intervals cannot: a nameless plane's line grazing an arc is a
/// double root, which the interval road never counts.
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
    if wcy.beside(wc).is_some() || wcy.grazes(wc) {
        return Some(Vec::new());
    }
    arc_crossings(jd, &wcy.def, wc, fc, cyl, ccw, ends, on_end)
}

/// **Where the class line `L = wc ∩ fc` touches the arc of cylinder `cyl` from `ends[0]` to
/// `ends[1]` without crossing it** — the tangent point, when `L` is tangent to the circle there and
/// the point lies strictly inside the arc; [`arc_crossings`]'s other half, which reads the same
/// meet and drops the tangent.
///
/// The point is named as [`arc_crossings`] names a root, with the one root a tangent has
/// (`QuadRoot::Double`). A tangent point at an arc's end is not this answer: that end is a node of
/// the ring, named by a pierce name. A class the gate passed silent about the cylinder is not
/// asked: where the whole cylinder lies beside its plane the line misses the arc, and on the other
/// silences every face on it stands apart from every lateral face's box — the class's plane may
/// cut the cylinder, but the caller's segment lies on such a face, outside the box that holds the
/// arc. `None` where the meet or the angular order cannot be read.
pub(crate) fn arc_touch(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    fc: usize,
    cyl: usize,
    ccw: bool,
    [a, b]: [NodeId; 2],
) -> Option<Option<NodeId>> {
    use nacre_judge::cylinder::LineMeet;
    let wcy = cyls.get(cyl)?;
    if wcy.is_silent(wc) {
        return Some(None);
    }
    match jd.line_meet(wc, fc, &wcy.def, false) {
        LineMeet::Tangent(..) => {}
        LineMeet::Pair(_) | LineMeet::ParallelPlanes | LineMeet::Miss => return Some(None),
        // `fc` carries an arc, so its plane meets the axis: a line on it is no ruling and not
        // parallel to the axis, and `wc` is another class. These three do not arise.
        LineMeet::CoincidentPlanes
        | LineMeet::OnRuling(_)
        | LineMeet::AxisParallelMiss
        | LineMeet::Undecided
        | LineMeet::Wide
        | LineMeet::Overflow => return None,
    }
    let (lo, hi) = if ccw { (a, b) } else { (b, a) };
    let node = NodeId::pierce(wc, fc, cyl, QuadRoot::Double);
    match theta_place(jd, cyl, &wcy.def, lo, hi, node)? {
        ThetaPlace::Inside => Some(Some(node)),
        ThetaPlace::Outside | ThetaPlace::At(_) => Some(None),
    }
}

/// **The planar scan's crossing on a ruling** — the point where the class line `L = wc ∩ fc`
/// leaves the face across an edge riding a cylinder's ruling, named as the pierce node
/// `wc ∩ fc ∩ cyl` at the root that *is* this ruling.
///
/// A ruling edge lies in the face's own plane `fc` (a plane holding a ruling runs through the
/// axis — or is tangent, which the gate passes without a cut, a class with no world statement as
/// grazing, and neither draws a ruling, so no ruling of it ever reaches here), so the pair
/// `{wc, fc}` cuts the cylinder in two
/// points, one on each of `fc`'s two rulings, and `(cyl, side)` — the identity
/// [`Wall::Ruling`] carries, measured by
/// [`ruling_side_signed`] against its `plane` when the ring was named (on a plane
/// face's ring, `fc` itself; on a lateral's, the plane across the edge, which the caller passes
/// as `fc`) — says which. The same predicate asked of each root picks it. `wc` must cross the
/// axis for the class to cross a ruling in a point at all (∥ contains it) — across it or oblique to
/// it; this asks nothing more.
///
/// Solved in the caller's order `(wc, fc)`, as `rulings_on_class` and `chord_stations`
/// spell it — [`NodeId::pierce`] canonicalizes the pair and the root together, so this is the one
/// name every road gives the point. The lateral face names the same point when its hole ring
/// crosses the class — `cycle_on_class`'s crossing arm asks this very door with the ruling's own
/// side — so the planar scan and the lateral's walk cannot name it apart.
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
    use nacre_judge::cylinder::LineMeet;
    let w = jd.world_rat(wc).ok_or(NoCrossing)?;
    if nacre_exact::dot_sign_rat(&[w[0], w[1], w[2]], &def.dir()) == nacre_exact::Orient::Zero {
        return Err(NoCrossing);
    }
    let meet = jd.line_meet(wc, fc, def, false);
    // A **tangent** wall (`side == 0`) has one ruling and one root — `Double`, which only the exact
    // road states.
    if side == 0 {
        return match meet {
            LineMeet::Tangent(..) => Ok(NodeId::pierce(wc, fc, cyl, QuadRoot::Double)),
            LineMeet::Overflow => Err(Width),
            _ => Err(NoCrossing),
        };
    }
    let roots = match meet {
        LineMeet::Pair(roots) => roots,
        LineMeet::Undecided => return Err(Undecided),
        LineMeet::Overflow => return Err(Width),
        _ => return Err(NoCrossing),
    };
    // A root whose side overflowed is not on `side` as far as this can tell; the other root may
    // still be, so the overflow is the answer only when no root is. A box that cannot say which
    // ruling it is on is the answer at once: the wall sections the cylinder, so its two rulings are
    // two distinct lines and each root is on one of them.
    let (mut found, mut overflow) = (None, false);
    for root in [QuadRoot::Lo, QuadRoot::Hi] {
        let at = roots.point(jd, cyl, def, root).ok_or(Undecided)?;
        match at.ruling_side(jd, fc, def) {
            None if at.exact().is_none() => return Err(Undecided),
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

/// Why [`crossing_on_ruling`] names no point — three causes, which its readers turn into their own
/// words (the tracer into its decline).
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
