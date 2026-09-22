use super::*;
/// The parity every clear ray reports. The ring is simple, so they must all agree; a golden
/// says so, which is a second machine for free.
pub(crate) fn every_ray(
    jd: &Judge<'_, WorkingPlane>,
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
    // ★★ **The unnameable check is the walk's, so it is asked per `qa` rather than once up
    // front**: a ring carrying a node the walk cannot read is refused by the first ray that
    // actually asks, not before any ray is cast. The two readings differ only when *every* `qa`
    // answers `AllOn` — a ring lying in both cut planes — and that pairs with a pierce node nothing
    // produces here, so the difference is unreachable twice over.
    //
    // ★ A **ring**, not a probe list: `ring_against_plane` reads it as a cyclic sign sequence, so
    // a dropped member would be a different polygon answered about confidently — which is why the
    // walk is handed the ring **whole** and answers `Unnameable` for the ring rather than letting
    // a caller drop a node (see [`three_plane_probes`], where dropping *is* honest).
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let mut out = Vec::new();
    for &qa in v.iter().filter(|&&x| x != p) {
        // Where the ring meets the line — the walk `trace_transversal_face` reads too.
        //
        // ★ **A `Crossing` needs no per-edge derivation.** "Is `X = {P, Q_a, R}` strictly inside
        // the edge", asked with two `order_along`s — `a` the sign of `X − From` along `P ∩ R` and
        // `b` that of `X − To`, both normalized to the same direction on the same line — holds iff
        // `a·b < 0` iff `From` and `To` lie on opposite sides of `Q_a`, which the walk already
        // knows from their sides. The parallel guard goes with it: an edge
        // whose line is parallel to `P ∩ Q_a` has both endpoints on one side and is not a crossing.
        // ★ No cylinder table on this road: a pierce node in a **result** cell's ring
        // declines here, and threading one is the arc road's business.
        // ★ The meet here is a **line** (`p ∩ Q_a`), and two points fix a line — so a straight
        // edge between two on-line nodes is on it and a curved one is not. Same reading as
        // `arrangement::trace_transversal_face`'s, which walks the same kind of ring.
        // ★ An arc between two on-line nodes would need a cylinder table to side
        // (`arc_departure_side`), and this road carries none — `None` is the honest answer,
        // and the walk's `Unnameable` is the same refusal a pierce node already meets here.
        let on_meet =
            |i: usize| (!matches!(ring[i].carrier, Carrier::Arc(_))).then_some(EdgeMeet::On);
        let features = match ring_against_plane(jd, &[], &nodes, qa, on_meet) {
            RingWalk::Met(f) => f,
            RingWalk::Unnameable => return Err(reject(RejectReason::PierceVertexUnnamed)),
            RingWalk::AllOn => continue, // the whole ring lies on `Q_a`
        };
        let qb = *v
            .iter()
            .find(|&&x| x != p && x != qa)
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        // A node on the line is a crossing point in its own right, so it must be nameable as one:
        // some plane of its own, off the line, pins it there. (`third_on_l` picks a handle the
        // same way, and for the same reason — one parallel to the line names no point on it.)
        // ★ The projection is safe *here* and nowhere earlier: the walk answered `Met`, which it
        // only does when every node is three planes.
        let namer = |i: usize| pin_on_line(jd, p, qa, three_plane_name(nodes[i])?);
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
            let (Some(lo), Some(hi)) = (namer(ends[0]), namer(ends[1])) else {
                unnameable = true;
                break;
            };
            let o = [
                order_along(jd, p, qa, lo, qb),
                order_along(jd, p, qa, hi, qb),
            ];
            if o[0] != o[1] || o[0] == 0 {
                return Err(reject(RejectReason::PointOnRing)); // `v` inside the run, or one of it
            }
            run_hits.push(if flanks_differ { o[0] } else { 0 });
        }
        if unnameable {
            continue;
        }
        for dir in [1i8, -1] {
            let mut crossings = run_hits.iter().filter(|&&o| o == dir).count();
            for f in &features {
                let Feature::Crossing { edge, .. } = *f else {
                    continue; // runs are counted above
                };
                // Strictly ahead of `v` along `dir · (n_P × n_Qa)`?
                let carrier = ring[edge]
                    .carrier
                    .wall()
                    .ok_or_else(|| reject(RejectReason::RingNaming))?;
                match order_along(jd, p, qa, carrier, qb) {
                    0 => return Err(reject(RejectReason::PointOnRing)), // `X == v`, inside an edge
                    o if o == dir => crossings += 1,
                    _ => {}
                }
            }
            out.push(crossings % 2 == 1);
        }
    }
    Ok(out)
}

/// **Does the open segment between two named points on one line meet this face's material?**
///
/// [`every_ray`]'s sibling, and not its copy. That one asks about a **point** and may bail with
/// `PointOnRing` when the point lands on the boundary — here the two endpoints are *expected* to,
/// because the defect this answers is an edge whose ends sit on a face's ring while its middle
/// crosses the interior. An interval query has no candidate to fall back to, so every case that one
/// declines has to become a value.
///
/// ★ **One algorithm, no special cases.** The rings meet the line `P ∩ w` at a set of places; the
/// two unbounded ends of the line are outside the face and each genuine crossing flips that, so the
/// line reads **outside / inside / outside / …**. The answer is whether any *inside* stretch
/// overlaps the open `(u, v)`. Written as branches — "is there a crossing between them", "is it in a
/// hole" — it was twice wrong, because each branch re-derived a piece of that structure and lost
/// another. Read as one alternation it also subsumes the endpoint test: an endpoint strictly inside
/// puts `u` in an inside stretch.
///
/// **Holes come along for free.** All rings go into one bag: the rings of a face are disjoint, so
/// even-odd over the union *is* the material region (a point inside a hole has crossed twice) — the
/// rule `design.md` states one dimension down for 2D sketches.
///
/// ★★ **A `Run` is a stretch, not a place.** `Run { len >= 2 }` means the boundary *lies along* the
/// line, so it occupies an interval where the segment would be **on** the face rather than inside
/// it — two faces sharing an edge, which is ordinary adjacency. Those stretches are boundary and
/// are not counted as inside; `flanks_differ` still says whether crossing the run flips the side,
/// which is the rule the tracer and the ray caster already read a ring with.
pub(crate) fn segment_meets_face(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    w: usize,
    u: [usize; 3],
    v: [usize; 3],
    rings: &[Vec<RingEdge>],
) -> Result<bool, BoolError> {
    // A point on `P ∩ w` is named there by a third plane of its own that **cuts** the line; one
    // parallel to it names nothing (the same duty `every_ray`'s `namer` states).
    let handle = |t: [usize; 3]| -> Option<usize> {
        t.iter()
            .copied()
            .find(|&x| x != p && x != w && jd.plane_pair_dir_sign(p, w, x) != 0)
    };
    let (Some(hu), Some(hv)) = (handle(u), handle(v)) else {
        return Err(reject(RejectReason::RingNaming));
    };
    // Every place a ring meets the line: `[lo, hi]` handles (equal for a crossing at a point) and
    // whether passing it flips inside/outside.
    let mut events: Vec<([usize; 2], bool)> = Vec::new();
    for ring in rings {
        // ★ A ring, whole: the walk reads it as a cyclic sign sequence — see [`three_plane_probes`]
        // for where dropping a node *is* honest.
        let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
        // No cylinder table on this road either (see `every_ray`): an arc answers `None`.
        let on_meet =
            |i: usize| (!matches!(ring[i].carrier, Carrier::Arc(_))).then_some(EdgeMeet::On);
        let features = match ring_against_plane(jd, &[], &nodes, w, on_meet) {
            RingWalk::Met(f) => f,
            RingWalk::Unnameable => return Err(reject(RejectReason::PierceVertexUnnamed)),
            // The whole ring lies on `w`: this face's boundary is the line itself, and the
            // alternation has no crossings to read. Refusing to guess.
            RingWalk::AllOn => return Err(reject(RejectReason::PointOnRing)),
        };
        for f in &features {
            match *f {
                Feature::Crossing { edge, .. } => {
                    let h = ring[edge]
                        .carrier
                        .wall()
                        .ok_or_else(|| reject(RejectReason::RingNaming))?;
                    events.push(([h, h], true));
                }
                Feature::Run {
                    first,
                    len,
                    flanks_differ,
                    ..
                } => {
                    let ends = [first, (first + len - 1) % nodes.len()];
                    // ★ The projection is safe *here*: the walk answered `Met`, which it only does
                    // when every node is three planes.
                    let name = |i: usize| -> Option<usize> { handle(three_plane_name(nodes[i])?) };
                    let (Some(a), Some(b)) = (name(ends[0]), name(ends[1])) else {
                        return Err(reject(RejectReason::RingNaming));
                    };
                    let lo_first = order_along(jd, p, w, a, b) <= 0;
                    events.push((if lo_first { [a, b] } else { [b, a] }, flanks_differ));
                }
            }
        }
    }
    events.sort_by(|x, y| match order_along(jd, p, w, x.0[0], y.0[0]) {
        -1 => std::cmp::Ordering::Less,
        1 => std::cmp::Ordering::Greater,
        _ => std::cmp::Ordering::Equal,
    });
    // Walk the line: outside before the first event, flipping as each genuine crossing is passed.
    // The stretch between two events is a cell; an inside cell that overlaps the open `(u, v)` is
    // the surface meeting itself.
    let (lo, hi) = if order_along(jd, p, w, hu, hv) <= 0 {
        (hu, hv)
    } else {
        (hv, hu)
    };
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
            // ★ Reaching the unbounded tail while *inside* means the rings crossed the line an odd
            // number of times, which a closed curve cannot do. Reading it as "outside" would let a
            // real contact past in silence, so it is named instead — the same rule the rest of this
            // engine follows for an invariant it cannot verify.
            return Err(reject(RejectReason::RingParity));
        };
        let cell_end = next.0[0];
        // Overlap with the **open** interval: strictly, so touching at `u` or `v` is not inside.
        if order_along(jd, p, w, cell_start, hi) < 0 && order_along(jd, p, w, cell_end, lo) > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Every **face** incident to `vh`, as `planes`-table slots. (It returns `inc`'s pairs verbatim,
/// and those are faces. Its one caller maps them through `plane_ix`.)
pub(crate) fn vertex_face_indices(vh: Handle<Vertex>, inc: &EdgeFaces) -> Vec<usize> {
    // ★ A lookup — the incidence table carries the vertex → faces map (built once per
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
        // ☑ A literal triple, so `side_of`'s pierce arm is unreachable here and the empty
        // cylinder table is never consulted — this asks about a *point*, not a ring.
        if side_of(jd, &[], NodeId::three_planes(Canon3::three(v)), r) != Some(0) {
            continue; // `v` is not even on the edge's line
        }
        // Name `v` as a point of that line: `{p, r, s}` for one of its own planes `s` off the line.
        let s = *v
            .iter()
            .find(|&&x| x != p && jd.plane_pair_dir_sign(p, r, x) != 0)
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        let (a, b) = (order_along(jd, p, r, s, si), order_along(jd, p, r, s, sj));
        if a * b <= 0 {
            return Ok(true); // between the endpoints (or on one)
        }
    }
    Ok(false)
}
