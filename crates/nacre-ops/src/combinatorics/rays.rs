use super::*;
/// The parity every clear ray reports. The ring is simple, so they must all agree; a golden
/// says so, which is a second machine for free.
pub(crate) fn every_ray(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    v: [usize; 3],
    ring: &[RingEdge],
) -> Result<Vec<bool>, BoolError> {
    // The vertex name is a plane triple, so it obeys the same rule as a ring's: class roots only.
    // A caller holding face indices (a hand-built table, a test) is normalized here rather than
    // silently comparing a face against a class.
    let mut v = [v[0], v[1], v[2]];
    v.sort_unstable();
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    // **The unnameable check is the walk's, so it is asked per `qa` rather than once up
    // front**: a ring carrying a node the walk cannot read is refused by the first ray that
    // actually asks, not before any ray is cast. The two readings differ only when *every* `qa`
    // answers `AllOn` — a ring lying in both cut planes — and that pairs with a pierce node nothing
    // produces here, so the difference is unreachable twice over.
    //
    // A **ring**, not a probe list: `ring_against_plane` reads it as a cyclic sign sequence, so
    // a dropped member would be a different polygon answered about confidently — which is why the
    // walk is handed the ring **whole** and answers `Unnameable` for the ring rather than letting
    // a caller drop a node (see [`three_plane_probes`], where dropping *is* honest).
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let mut out = Vec::new();
    for &qa in v.iter().filter(|&&x| x != p) {
        // Where the ring meets the line — the reading every road that walks a plane face's ring
        // along a line shares ([`ClassLine`]).
        //
        // **A straight edge's `Crossing` needs no per-edge derivation.** "Is `X = {P, Q_a, R}` strictly inside
        // the edge", asked with two `order_along`s — `a` the sign of `X − From` along `P ∩ R` and
        // `b` that of `X − To`, both normalized to the same direction on the same line — holds iff
        // `a·b < 0` iff `From` and `To` lie on opposite sides of `Q_a`, which the walk already
        // knows from their sides. The parallel guard goes with it: an edge
        // whose line is parallel to `P ∩ Q_a` has both endpoints on one side and is not a crossing.
        // No cylinder table on this road for a class with a world statement ([`silent_table`]):
        // a pierce node in a **result** cell's ring declines here (`nesting` reads it as the road
        // declining, and the coordinate road answers), and so does an arc whose crossings or
        // departure need the table. A class the gate passed as silent reads the gate's record
        // instead: against it an arc crosses nowhere, or — where the class sections the arc's
        // cylinder — as often as the boxes count.
        let table = silent_table(cyls, qa);
        let line = ClassLine::new(jd, table, p, qa, &nodes, |i| (&ring[i].carrier).into());
        let features = match line.walk() {
            RingWalk::Met(f) => f,
            RingWalk::Unnameable => return Err(reject(RejectReason::PierceVertexUnnamed)),
            RingWalk::AllOn => continue, // the whole ring lies on `Q_a`
        };
        // The plane `v` is ordered against along the line, by the rule every handle on a line
        // takes: one of `v`'s own planes that cuts `p ∩ qa`. When `p` is not among `v`'s planes
        // one of the other two can hold the line (a pencil), and read there it orders nothing.
        let qb = pin_on_line(jd, p, qa, &v).ok_or_else(|| reject(RejectReason::RingNaming))?;
        // Where a point on the line lies against `v` along `n_P × n_Qa` ([`order_on`]). Runs and
        // crossings both ask it.
        let against_v = |at: PointOn| -> Result<i8, BoolError> {
            order_on(jd, table, p, qa, at, PointOn::Class(qb))
                .ok_or_else(|| reject(RejectReason::WitnessNotRational))
        };
        // Which side of `v` each run sits on, and whether the ring crossed the line there at all.
        // A run is one interval of the line and the ring is simple, so it cannot double back
        // inside itself: its two ends bracket it, and ends that disagree mean `v` is *between*
        // them — on the ring.
        let mut run_hits: Vec<i8> = Vec::new();
        let mut unnameable = false;
        for f in &features {
            let Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } = *f
            else {
                continue;
            };
            let ends = [first, (first + len - 1) % nodes.len()];
            // A node on the line is a crossing point in its own right, pinned there by the one
            // rule ([`pin_for`]): a plane of its own that cuts the line, or — for a pierce corner of
            // the line's own pair — its cylinder. A name nothing of its own pins here (every plane
            // of it holds the line) leaves this ray uncast.
            let pin = |i: usize| end_on_line(jd, p, qa, nodes[i]);
            let (Some(lo), Some(hi)) = (pin(ends[0]), pin(ends[1])) else {
                unnameable = true;
                break;
            };
            let o = [against_v(lo)?, against_v(hi)?];
            if o[0] != o[1] || o[0] == 0 {
                return Err(reject(RejectReason::PointOnRing)); // `v` inside the run, or one of it
            }
            run_hits.push(if flanks_differ { o[0] } else { 0 });
        }
        if unnameable {
            continue;
        }
        // Where each crossing lies against `v`, once, in ring order — so the first refusal is the
        // one met first along the ring, as a pass per direction met it. A plane wall's crossing is
        // the point `{P, Q_a, wall}`; an arc's is the pierce point the walk counted.
        // Crossings strictly ahead of `v` along `+d` and along `−d`.
        let mut ahead = [0usize; 2];
        for f in &features {
            let Feature::Crossing { edge, nth, .. } = *f else {
                continue; // runs are counted above
            };
            let at = line
                .crossing(edge, nth)
                .ok_or_else(|| reject(RejectReason::RingNaming))?;
            let o = against_v(at)?;
            if o == 0 {
                return Err(reject(RejectReason::PointOnRing)); // `X == v`, inside an edge
            }
            ahead[usize::from(o < 0)] += 1;
        }
        for (k, dir) in [1i8, -1].into_iter().enumerate() {
            let crossings = run_hits.iter().filter(|&&o| o == dir).count() + ahead[k];
            out.push(crossings % 2 == 1);
        }
    }
    Ok(out)
}

/// **Does the open segment between two named points on `P ∩ w` meet this face's material?**
///
/// [`every_ray`]'s sibling, and not its copy. That one asks about a **point** and may bail with
/// `PointOnRing` when the point lands on the boundary — here the two endpoints are *expected* to,
/// because the defect this answers is an edge whose ends sit on a face's ring while its middle
/// crosses the interior. An interval query has no candidate to fall back to, so every case that one
/// declines has to become a value.
///
/// **One algorithm, no special cases.** The rings meet the line `P ∩ w` at a set of places; the
/// two unbounded ends of the line are outside the face and each genuine crossing flips that, so the
/// line reads **outside / inside / outside / …**. The answer is whether any *inside* stretch
/// overlaps the open `(u, v)`. Written as branches — "is there a crossing between them", "is it in
/// a hole" — each branch re-derives a piece of that structure and loses another. Read as one
/// alternation it also subsumes the endpoint test: an endpoint strictly inside puts `u` in an
/// inside stretch.
///
/// **Holes come along for free.** All bounds go into one bag: the bounds of a face are disjoint, so
/// even-odd over the union *is* the material region (a point inside a hole has crossed twice) — the
/// rule `design.md` states one dimension down for 2D sketches. A circle hole meets the line at the
/// two points where it crosses the hole's cylinder ([`class_line_meets_cylinder`] — a tangent line
/// touches without crossing), and they toggle like any crossing. The circles are read only when
/// the polygon bounds alone leave an inside stretch overlapping the segment: a hole only takes
/// material away, so a «no» needs no circle, and asking one can fail where the answer needed none.
/// A face bounded outside by a circle (a cap disk) is not this road's: its polygon bounds are holes.
///
/// **A `Run` is a stretch, not a place.** `Run { len >= 2 }` means the boundary *lies along* the
/// line, so it occupies an interval where the segment would be **on** the face rather than inside
/// it — two faces sharing an edge, which is ordinary adjacency. Those stretches are boundary and
/// are not counted as inside; `flanks_differ` still says whether crossing the run flips the side,
/// which is the rule the tracer and the ray caster already read a ring with.
///
/// **A ring a cylinder trimmed is read as the tracer reads it** — through [`ClassLine`] with the
/// whole cylinder table, so a pierce corner is sided by [`side_of`]'s pierce arm and an arc's and
/// a ruling's crossings are counted and named; a run's ends are pinned by [`pin_for`], a pierce
/// corner by the plane of its own pair that cuts the line. No part of that reading sits
/// behind the gate: one of the edge's own faces may lie on `P` facing the other way (a box's
/// bottom over a bore, beside the plate's top), the other may then stand parallel to the
/// cylinder's axis, which the gate passes, and the line crosses the ring's arc.
///
/// **The segment's own two ends are pinned by the same rule** ([`pin_for`]), so an end a cylinder
/// made is asked like a three-plane one: a lateral's own ruling edge, its ends pierce points of the
/// flat face it bounds and a cap, is pinned by the cap. An end the rule does not pin is
/// `RingNaming`. The risk is a pencil: where planes beyond `P` and `w` hold the edge's line, a pierce
/// end named by a pair other than `P` and `w` can have both its planes holding the line, and then
/// neither cuts it (named by `P` and `w` themselves, its cylinder pins it) — which pair the name
/// carries, against which of the edge's faces gave `w`, comes down to the order of the faces, and
/// so does whether the end is pinned. Pinning it regardless needs its root carried from one pair to the other,
/// which nothing does.
pub(crate) fn segment_meets_face(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    w: usize,
    u: NodeId,
    v: NodeId,
    face: &CompFace,
) -> Result<bool, BoolError> {
    let (Some(hu), Some(hv)) = (end_on_line(jd, p, w, u), end_on_line(jd, p, w, v)) else {
        return Err(reject(RejectReason::RingNaming));
    };
    let BoundEdges::Ring(_) = &face.outer else {
        return Err(reject(RejectReason::RingNaming));
    };
    // Two places on the line, ordered ([`order_on`]).
    let order = |a: PointOn, b: PointOn| -> Result<i8, BoolError> {
        order_on(jd, cyls, p, w, a, b).ok_or_else(|| reject(RejectReason::WitnessNotRational))
    };
    // Every place a bound meets the line: `[lo, hi]` (equal for a crossing at a point) and whether
    // passing it flips inside/outside.
    let mut events: Vec<([PointOn; 2], bool)> = Vec::new();
    let mut circles: Vec<&nacre_topo::CylinderDef> = Vec::new();
    for bound in core::iter::once(&face.outer).chain(&face.inner) {
        let ring = match bound {
            BoundEdges::Ring(ring) => ring,
            BoundEdges::Circle(def) => {
                circles.push(def);
                continue;
            }
            // Loops bound a cylinder, never a plane — a producer error, not an input.
            BoundEdges::Lateral(_) => return Err(reject(RejectReason::RingNaming)),
        };
        // A ring, whole: the walk reads it as a cyclic sign sequence — see [`three_plane_probes`]
        // for where dropping a node *is* honest.
        let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
        let line = ClassLine::new(jd, cyls, p, w, &nodes, |i| (&ring[i].carrier).into());
        let features = match line.walk() {
            RingWalk::Met(f) => f,
            RingWalk::Unnameable => return Err(reject(RejectReason::PierceVertexUnnamed)),
            // The whole ring lies on `w`: this face's boundary is the line itself, and the
            // alternation has no crossings to read. Refusing to guess.
            RingWalk::AllOn => return Err(reject(RejectReason::PointOnRing)),
        };
        for f in &features {
            match *f {
                Feature::Crossing { edge, nth, .. } => {
                    let at = line
                        .crossing(edge, nth)
                        .ok_or_else(|| reject(RejectReason::RingNaming))?;
                    events.push(([at; 2], true));
                }
                Feature::Run {
                    first,
                    len,
                    flanks_differ,
                    ..
                } => {
                    let ends = [first, (first + len - 1) % nodes.len()];
                    // A run's end is pinned on the line by the one rule ([`pin_for`]): a pierce
                    // corner — the end of a ruling the line runs along, when the line is parallel
                    // to the axis — by the plane of its own pair that cuts the line.
                    let pin = |i: usize| end_on_line(jd, p, w, nodes[i]);
                    let (Some(a), Some(b)) = (pin(ends[0]), pin(ends[1])) else {
                        return Err(reject(RejectReason::RingNaming));
                    };
                    let lo_first = order(a, b)? <= 0;
                    events.push((if lo_first { [a, b] } else { [b, a] }, flanks_differ));
                }
            }
        }
    }
    let (lo, hi) = if order(hu, hv)? <= 0 {
        (hu, hv)
    } else {
        (hv, hu)
    };
    let by_place = |x: &([PointOn; 2], bool), y: &([PointOn; 2], bool)| {
        order(x.0[0], y.0[0]).map(|o| o.cmp(&0))
    };
    try_sort_by(&mut events, by_place)?;
    let polygons = inside_overlaps(&events, lo, hi, &order)?;
    if !polygons || circles.is_empty() {
        return Ok(polygons);
    }
    for def in circles {
        let cyl = cyls
            .iter()
            .position(|c| c.def == *def)
            .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
        // A class with no world statement the gate proved this cylinder beside: the line lies in
        // that plane and the whole cylinder off it, so the line misses the circle. Proved apart
        // only, it may cross the circle at points no world statement names — not read.
        if cyls[cyl].beside(w).is_some() {
            continue;
        }
        use nacre_exact::quad::CylinderMeet;
        match class_line_meets_cylinder(jd, p, w, def) {
            Some(CylinderMeet::Pair { .. }) => {
                for root in [QuadRoot::Lo, QuadRoot::Hi] {
                    let at = PointOn::Pierce(NodeId::pierce(p, w, cyl, root));
                    events.push(([at; 2], true));
                }
            }
            // A tangent line touches the circle without crossing it; a missing one meets nothing.
            Some(
                CylinderMeet::Tangent { .. }
                | CylinderMeet::Miss(_)
                | CylinderMeet::AxisParallelMiss(_),
            ) => {}
            _ => return Err(reject(RejectReason::WitnessNotRational)),
        }
    }
    try_sort_by(&mut events, by_place)?;
    inside_overlaps(&events, lo, hi, &order)
}

/// **Walk the line** over `events` sorted along it: outside before the first event, flipping as
/// each genuine crossing is passed. The stretch between two events is a cell; an inside cell that
/// overlaps the open `(lo, hi)` is the surface meeting itself.
fn inside_overlaps(
    events: &[([PointOn; 2], bool)],
    lo: PointOn,
    hi: PointOn,
    order: &impl Fn(PointOn, PointOn) -> Result<i8, BoolError>,
) -> Result<bool, BoolError> {
    let mut inside = false;
    for i in 0..events.len() {
        if events[i].1 {
            inside = !inside;
        }
        if !inside {
            continue;
        }
        // The cell runs from this event's far end to the next event's near end.
        let cell_start = events[i].0[1];
        let Some(next) = events.get(i + 1) else {
            // Reaching the unbounded tail while *inside* means the bounds crossed the line an
            // odd number of times, which a closed curve cannot do. Reading it as "outside" would
            // let a real contact past in silence, so it is named instead — the same rule the rest
            // of this engine follows for an invariant it cannot verify.
            return Err(reject(RejectReason::RingParity));
        };
        let cell_end = next.0[0];
        // Overlap with the **open** interval: strictly, so touching at `lo` or `hi` is not inside.
        if order(cell_start, hi)? < 0 && order(cell_end, lo)? > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Every **face** incident to `vh`, as `planes`-table slots. (It returns `inc`'s pairs verbatim,
/// and those are faces. Its one caller maps them through `plane_ix`.)
pub(crate) fn vertex_face_indices(vh: Handle<Vertex>, inc: &EdgeFaces) -> Vec<usize> {
    // A lookup — the incidence table carries the vertex → faces map (built once per
    // operand), rather than a scan of every edge for every vertex asked.
    inc.faces_at(vh).to_vec()
}

/// Whether the point named by plane triple `v` lies on any edge of `ring` (a ring on plane `p`).
/// Used by [`point_in_component`] to abandon a non-generic ray rather than guess on a boundary.
pub(super) fn point_on_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    mut v: [usize; 3],
    ring: &[RingEdge],
) -> Result<bool, BoolError> {
    v.sort_unstable();
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    for e in ring {
        let (Some(si), Some(sj)) = (e.from_h.class(), e.to_h.class()) else {
            return Err(reject(RejectReason::RingNaming));
        };
        let r = e
            .carrier
            .wall()
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        // A literal triple, so `side_of`'s pierce arm is unreachable here and the empty
        // cylinder table is never consulted — this asks about a *point*, not a ring.
        #[expect(
            clippy::disallowed_methods,
            reason = "a point question: is the literal vertex `v` on this edge's line"
        )]
        let on_line = side_of(jd, &[], NodeId::three_planes(Canon3::three(v)), r) == Some(0);
        if !on_line {
            continue; // `v` is not even on the edge's line
        }
        // Name `v` as a point of that line: `{p, r, s}` for one of its own planes `s` off the line.
        let s = pin_on_line(jd, p, r, &v).ok_or_else(|| reject(RejectReason::RingNaming))?;
        let (a, b) = (order_along(jd, p, r, s, si), order_along(jd, p, r, s, sj));
        if a * b <= 0 {
            return Ok(true); // between the endpoints (or on one)
        }
    }
    Ok(false)
}

/// **A named point on the line `a ∩ b`, as the place a comparison along it takes** — pinned by
/// [`pin_for`]; `None` where nothing of the name pins it there.
fn end_on_line(jd: &Judge<'_, WorkingPlane>, a: usize, b: usize, name: NodeId) -> Option<PointOn> {
    Some(PointOn::of((name, pin_for(jd, a, b, name)?)))
}

/// **The cylinder table [`every_ray`] reads a ring against plane `q` with** ([`ClassLine`]). This road
/// carries none for a class with a world statement: a pierce corner there has no side here and the
/// walk declines, as an arc's crossings and departure do — the coordinate road answers instead. A
/// class the gate passed as silent about the cylinders is the exception the gate's record makes:
/// every pierce corner and every arc of a cylinder that lies beside it is on the recorded side
/// ([`crate::planes::WorkingCyl::beside`]), so the table is handed over for [`side_of`] to read,
/// and such an arc crosses `q` nowhere ([`arc_crossings_of`]). A class that **sections** a cylinder
/// (its plane cuts the lateral, its faces stand apart — `Silence::Sectioned`) has no recorded side
/// for that cylinder's arcs: how often one crosses the line `p ∩ q` is read off interval boxes
/// ([`arc_crossings`] on the judge's interval road), the count the arrangement's own arc road gives the same arc. A class
/// silent by `Silence::Apart` meets no arc of the cylinder at all.
fn silent_table(cyls: &[crate::planes::WorkingCyl], q: usize) -> &[crate::planes::WorkingCyl] {
    if crate::planes::silent_among(cyls, &[q]) {
        cyls
    } else {
        &[]
    }
}
