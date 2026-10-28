use super::*;
/// `split_circles`' product: `None` when nothing crossed (the caller keeps its slices), else the
/// split segments, the still-whole circles, the arcs, and each cut circle's [`CutRim`].
type SplitCircles = Option<(
    Vec<MergedSeg>,
    Vec<MergedCircle>,
    Vec<Piece>,
    Vec<(usize, CutRim)>,
)>;

/// **Cut every circle a segment crosses into arcs, and the segments with it.**
///
/// This is the arc split. What it does *not* do is find the crossings — [`circle_crossings`] names
/// them. The work here is ordering: around the circle (θ,
/// [`nacre_exact::quad::circular_order_about_seam`]) to make arcs, and along each segment
/// ([`combinatorics::order_located`]) to make sub-segments.
///
/// A circle nothing crosses is returned **whole**, on its own road. The parallel
/// road shrinks to the population it is actually about.
pub(super) fn split_circles(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    circles: &[MergedCircle],
    aliases: &Aliases,
) -> Result<SplitCircles, BoolError> {
    // Nothing to cross. The ends below are built for the whole slice at once, so the empty case has
    // to be answered before them or an empty class would start declining. `split_rulings` opens
    // with the same line. Measured: **61993** times over the suite — this is the common case, not
    // an edge one.
    if circles.is_empty() {
        return Ok(None);
    }
    // Which nodes land on each circle, and which on each segment. Collected together because one
    // crossing is a point of both — splitting only one of them would leave the other's edge running
    // through a vertex it does not have.
    let canon = |x: NodeId| aliases.canon_point(x);
    let mut on_circle: Vec<Vec<NodeId>> = vec![Vec::new(); circles.len()];
    // **A contribution's own ends are cut points, and that is a rule rather than a hope.**
    // Where a lateral face's hole starts is where the material a class sees changes, and the
    // arrangement's own principle is that such a place is a vertex. Feeding both ends of every
    // extent in here is what makes it one — after which no arc straddles a boundary and the
    // distribution below is a containment test with no partial case to arbitrate.
    for (ci, circ) in circles.iter().enumerate() {
        for &(_, _, arc) in &circ.merged {
            if let Some(ends) = arc {
                on_circle[ci].extend(ends);
            }
        }
    }
    // **Where a segment's ends are, asked once per segment — and a missing one is not
    // fatal.** A pierce end has no rational coordinate at any width, and nothing below needs a
    // coordinate to decide anything; the one thing that reads it is a **filter**.
    let seg_coords: Vec<[Option<[Rat; 3]>; 2]> = segs
        .iter()
        .map(|s| s.end.map(|t| combinatorics::node_coords_rat(jd, t)))
        .collect();
    // The circle of cylinder `cyl` — one per cylinder class, as `merge_circles` groups them.
    let index = |cyl: usize| {
        circles
            .iter()
            .position(|c| c.cyl == cyl)
            .expect("a cylinder this split listed")
    };
    let here: Vec<usize> = circles.iter().map(|c| c.cyl).collect();
    let (mut on_seg, cuts) = collect(
        jd,
        cyls,
        wc,
        segs,
        &here,
        |si, sg, cyl| {
            let circ = &circles[index(cyl)];
            let def = &cyls[cyl].def;
            let (o, m, r2) = (def.origin(), def.dir(), def.r2());
            // **The cheap rejection runs where it can be formed, and skipping it asks
            // *more*, not less.** `segment_meets_cylinder` needs both endpoints' coordinates; where
            // one is a pierce point the question goes straight to `circle_crossings`, which answers
            // `Miss` exactly when the line misses. Losing the filter costs a solve, never an answer.
            if let [Some(p0), Some(p1)] = &seg_coords[si]
                && (p0 == p1 || !nacre_exact::segment_meets_cylinder(p0, p1, &o, &m, r2))
            {
                return Ok(Vec::new());
            }
            let Some(xs) = circle_crossings(jd, wc, circ, def, sg)? else {
                return Err(reject(RejectReason::WitnessNotRational));
            };
            Ok(xs
                .into_iter()
                // A tangency touches without separating — `segment_meets_cylinder` above already
                // let that shape through, and cutting there would make a zero-length arc.
                //
                // **The skip is right, and a tangency mints no vertex** — the circle it grazes
                // stays whole, so a circle the split leaves with one node is always one a vertex
                // lies on (`combinatorics::ring_floor`). Three reasons: measured, the link at such
                // a touch is a *single* circle and OCCT returns the same body, so minting a vertex
                // would invent structure that is not there; the `UnorderedEdges` that lifting this
                // skip would reach is the **rulings** vocabulary (`side` ±1), which the tangent
                // wall never enters because the gate records no crossing for it; and this arm is
                // not gated on the gate's cut record, so a ⊥ class takes the same road whatever the gate
                // recorded. The 15 measured cells with no third plane on the tangent line assemble
                // with it in place.
                .filter(|&n| {
                    !matches!(
                        combinatorics::pierce_name(n),
                        Some((_, _, crate::combinatorics::QuadRoot::Double))
                    )
                })
                .map(|n| (n, None))
                .collect())
        },
        // **A crossing shown off every arc the class carries is no cut point** — asked against each
        // arc's two ends alone ([`combinatorics::theta_between`]), under the table's name, so the
        // points of the circle's continuation are never ordered against each other: two names
        // there for one point would refuse a point no edge of this class has. An arc's own ends
        // are cut points by construction (fed in above), and a crossing that is another name of
        // one is that end. Three circles drop nothing: one with a contribution that is the whole
        // circle (every crossing is on it), one with no contribution (where «off every arc» would
        // be an empty truth), and one with a closed arc (`p == q`), which `theta_between` reads as
        // the whole turn less its end and so answers `false` for another name of the end itself.
        // What cannot be shown off stays a cut point: off an arc it is a vertex where nothing
        // crosses the segment, which costs a vertex, where dropping a crossing on an arc would
        // lose one. Cutting every crossing of the whole circle instead puts a vertex on the segment
        // where nothing crosses it (measured: 4 per op on the `tangentline slab 1.6` controls, 16 on
        // `rrect-box`, 0 on the other 353 census rows), and where that point *is* a vertex of the
        // segment already (a face through the top of a fillet's circle: the user's gusset foot,
        // `y = 0 = −2 + r`) the segment split refuses it as two names for one point.
        //
        // The crossing travels under the table's name: a class through a corner names the corner
        // again here, and the tables key by name.
        |cyl, n, _meet| {
            let ci = index(cyl);
            let x = canon(n);
            let on = Landing::On {
                cut: Some(ci),
                node: x,
            };
            let Some(extents) = circles[ci]
                .merged
                .iter()
                .map(|&(_, _, arc)| arc)
                .collect::<Option<Vec<[NodeId; 2]>>>()
            else {
                return on;
            };
            let off = !extents.is_empty()
                && extents.iter().all(|&[p, q]| canon(p) != canon(q))
                && !extents.iter().flatten().any(|&e| canon(e) == x)
                && extents.iter().all(|&[p, q]| {
                    combinatorics::theta_between(jd, cyl, &cyls[cyl].def, p, q, x) == Some(false)
                });
            if off { Landing::Off } else { on }
        },
    )?;
    for (ci, n) in cuts {
        on_circle[ci].push(n);
    }

    // Nothing crossed: the caller keeps its own slices and no edge is copied. The common case
    // pays for the question and not for an answer it does not need.
    if on_circle.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    let mut out_circles = Vec::new();
    let mut arcs = Vec::new();
    let mut cut_rims: Vec<(usize, CutRim)> = Vec::new();

    for (ci, circ) in circles.iter().enumerate() {
        let raw = std::mem::take(&mut on_circle[ci]);
        if raw.is_empty() {
            out_circles.push(circ.clone());
            continue;
        }
        let is_end = |n: NodeId| {
            circ.merged
                .iter()
                .any(|&(_, _, arc)| arc.is_some_and(|e| e.contains(&n)))
        };
        // One node per point, asked of the table; an arc's own end keeps its name, because the
        // arc is looked up by it below.
        let mut nodes: Vec<NodeId> = Vec::with_capacity(raw.len());
        for n in raw {
            match nodes.iter().position(|&m| canon(m) == canon(n)) {
                Some(i) => {
                    if is_end(n) && !is_end(nodes[i]) {
                        nodes[i] = n;
                    }
                }
                None => nodes.push(n),
            }
        }
        nodes.sort_unstable();
        // ---- the circle → arcs, in θ order about the seam ----
        let order = combinatorics::circular_order(jd, circ.cyl, &cyls[circ.cyl].def, &nodes)
            .map_err(|e| reject(e.reason()))?;
        // θ position of each node, so an extent can be tested against an arc by rank alone.
        let mut at = vec![0usize; order.len()];
        for (k, &i) in order.iter().enumerate() {
            at[i] = k;
        }
        let place = |x: NodeId| {
            nodes
                .iter()
                .position(|&m| canon(m) == canon(x))
                .map(|i| at[i])
                .expect("an extent's ends were fed into this very split")
        };
        for k in 0..order.len() {
            // **Each arc takes only the contributions that cover it.** A lateral face with a
            // hole marks its class over an *arc*, and the arc inside the hole is a piece of circle
            // that face does not bound at all — copying the circle's whole contribution list onto
            // it would flip the bits of a face that is not there. Cyclic containment by rank:
            // `[p, q)` runs counter-clockwise from `p`, so arc `k` is inside it exactly when `k`
            // sits between their ranks the same way round.
            let merged = circ
                .merged
                .iter()
                .filter(|&&(_, _, arc)| match arc {
                    None => true,
                    Some([p, q]) => {
                        let (a, b) = (place(p), place(q));
                        if a <= b {
                            a <= k && k < b
                        } else {
                            a <= k || k < b
                        }
                    }
                })
                .map(|&(s, kind, _)| (s, kind))
                .collect::<Vec<_>>();
            // **A piece no contribution covers is not an edge**. A partial rim — a
            // fillet's quarter, a slot's half — states its arc alone, and the rest of the circle
            // bounds no face on this class; emitted, it stood as a phantom edge that tied with
            // the tangent wall's line at the fillet's corner (`UnorderedEdges`, measured). The
            // straight road drops its newsless pieces (`drop_newsless`); this is the arc's twin.
            // The rim's cut nodes stay in `cut_rims` either way.
            if merged.is_empty() {
                continue;
            }
            arcs.push(Piece {
                on: OnCurve::Arc { cyl: circ.cyl },
                end: [nodes[order[k]], nodes[order[(k + 1) % order.len()]]],
                merged,
            });
        }
        cut_rims.push((
            circ.cyl,
            CutRim {
                nodes: order.iter().map(|&k| nodes[k]).collect(),
            },
        ));
    }
    // ---- segments → sub-segments, in line order ----
    let out_segs = split_segments_at(jd, cyls, wc, segs, &mut on_seg, aliases)?;

    Ok(Some((out_segs, out_circles, arcs, cut_rims)))
}

/// **Cut every segment at the points collected on it** — the one emitter the arc split and the
/// ruling split share.
///
/// **One emitter for both splits.** Both roads collect crossings onto their segments and then have
/// to put them in order, refuse the ones that coincide, and hand each piece the sense the whole
/// segment ran in — one piece of work, so one copy. One consequence of sharing it: an end a
/// cylinder pins passes both splits or neither, rather than one road being taught and the other
/// forgotten.
///
/// `on_seg[si]` is drained: a segment nobody crossed is passed through untouched, and the common
/// case copies nothing.
pub(super) fn split_segments_at(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    on_seg: &mut [Vec<NodeId>],
    aliases: &Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let mut out_segs = Vec::with_capacity(segs.len());
    for (si, sg) in segs.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_seg[si]);
        if nodes.is_empty() {
            out_segs.push(sg);
            continue;
        }
        nodes.sort_unstable();
        // Identity is the table's — two crossings the table knows as one point are one.
        nodes.dedup_by(|a, b| aliases.canon_point(*a) == aliases.canon_point(*b));
        // **No ruler: the points are ordered by the rule, not by a parameter.** Solving every
        // crossing's `(line, s)`, putting both endpoints on that same line with `along`, and
        // sorting the parameters would stop at an end a cylinder pinned: it has no rational
        // coordinate for `along` to take. [`combinatorics::order_located`] answers the same
        // question from the pair `{wc, sg.wall}` and the points' own names, so there is no ruler to
        // lay and nothing to convert.
        //
        // **The pair is the sorted one, and that is not cosmetic.** `pierce_meet`'s canonical
        // line runs along `n_min × n_max`; asking in call order would flip every comparison on a
        // class where `wc > wall`, taking the emitted pieces and the sense with it.
        //
        // **The guard is a check on the crossings' *names*.** Every crossing on this segment must
        // be a point of `wc ∩ sg.wall`, and the name says so more cheaply than two `MeetLine`
        // comparisons. Measured **unexercised** over the suite: it is here because the sort has no
        // way to notice, not because a fixture is red.
        let pair = {
            let mut p = [wc, sg.wall];
            p.sort_unstable();
            p
        };
        let mut keyed: Vec<(combinatorics::PointOn, NodeId, combinatorics::EndPin)> = Vec::new();
        for &n in &nodes {
            if combinatorics::pierce_name(n).is_none() {
                return Err(reject(RejectReason::RingNaming));
            }
            // The crossing arrives under the table's name, which may carry another
            // pair — a corner's own, or a secant's on a line two classes share on the cylinder —
            // so its pin on *this* line is derived, not assumed (`pin_for`: the cylinder when the
            // pair is the line's, else the plane of its pair that cuts the line), and the point
            // is located by that pin, the ends' own rule: a pierce name is a locator only on
            // the line its pair spells (`on_line`'s `names_the_line`).
            let Some(pin) = combinatorics::pin_for(jd, pair[0], pair[1], n) else {
                return Err(reject(RejectReason::RingNaming));
            };
            keyed.push((Split::of(n, pin).on(wc, sg.wall), n, pin));
        }
        // **A crossing that *is* an endpoint is one point with one name, so it is deduped rather
        // than refused.** "One point wearing two names, a three-plane one and a pierce one" is true
        // and refused, but only for the case it describes: a crossing at a *three-plane* end really
        // does carry a second name, and the equality check below catches it. Where the end was
        // pinned by a cylinder, its name **is** the crossing's `Pierce{planes, cyl, root}` — the
        // same vertex, arrived at twice — and shipping it twice would put a zero-length piece
        // between a point and itself.
        //
        // **"the same name" is asked of the alias table.** A tangent corner is a point
        // the cylinder crossing names `Pierce{[cap, wc], root}` and the plane road names
        // `ThreePlane([cap, t, wc])`; `Aliases::record_on_cylinder` has joined them, so the
        // crossing folds onto the end here as it would had the names been equal. The end's own
        // name and pin stay on the slot — they are the canonical ones (`merge_coincident`) and
        // the pieces emitted below are keyed by them. What the table does *not* know still
        // meets the equality check below, and that refusal is the honest floor.
        //
        // This is why the ends' indices are remembered instead of read off the tail of the
        // vector: after a dedup the last two entries are no longer the two ends.
        // Measured over the suite: the dedup fires **400** times and the equality check below
        // still refuses **1** — the two cases are separated by `NodeId`, and both happen.
        let mut ends = (0usize, 0usize);
        for k in 0..2 {
            let same = |e: &(combinatorics::PointOn, NodeId, combinatorics::EndPin)| {
                aliases.canon_point(e.1) == aliases.canon_point(sg.end[k])
            };
            let slot = match keyed.iter().position(same) {
                Some(ix) => {
                    if keyed[ix].1 != sg.end[k] {
                        keyed[ix] = (
                            Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
                            sg.end[k],
                            sg.end_h[k],
                        );
                    }
                    ix
                }
                None => {
                    keyed.push((
                        Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
                        sg.end[k],
                        sg.end_h[k],
                    ));
                    keyed.len() - 1
                }
            };
            if k == 0 { ends.0 = slot } else { ends.1 = slot }
        }
        // **Nothing between them: every crossing was an end, so the segment is not cut.** It goes
        // back whole rather than being re-emitted as a single piece — which would recompute a sense
        // it already carries, and could decline where the untouched segment does not.
        if keyed.len() == 2 {
            out_segs.push(sg);
            continue;
        }
        // The comparison's **second**-argument form, once per point — the hoist `Located`'s own doc
        // measured: a first argument is built a few times, a second is stored and asked repeatedly.
        let on: Vec<combinatorics::OnLine> = keyed
            .iter()
            .map(|k| combinatorics::on_line(jd, cyls, pair[0], pair[1], k.0).ok_or_else(undecided))
            .collect::<Result<_, BoolError>>()?;
        let cmp = |i: usize, j: usize| -> Option<core::cmp::Ordering> {
            let a = combinatorics::locate(jd, cyls, pair[0], pair[1], keyed[i].0)?;
            Some(
                match combinatorics::order_located(jd, pair[0], pair[1], &a, &on[j])? {
                    -1 => core::cmp::Ordering::Less,
                    1 => core::cmp::Ordering::Greater,
                    _ => core::cmp::Ordering::Equal,
                },
            )
        };
        // A crossing that lands **on** an endpoint is one point wearing two names — a three-plane
        // one and a pierce one. Names the alias table knows as one point were folded above; what
        // reaches the order is a coincidence no discovery event recorded ([`order_distinct`]).
        let order = order_distinct(keyed.len(), cmp, undecided)?;
        // **The sub-segments' travel sense, taken from the whole segment's own two named
        // ends — exactly, with no frame arithmetic at all.** Every sub-segment runs the way the
        // segment ran, and the pieces below are emitted in **ascending** order; so the only
        // question is whether ascending *is* the segment's direction, and the two ends answer it.
        // (The alternative — turning a meet line's `dir` into an `EdgeDir` sense — needs the
        // canonical-vs-stored turn on two classes and an `f64` dot to decide it. This needs
        // neither: `edge_dir` is still the one place a sense is made.)
        let forward = match cmp(ends.0, ends.1) {
            Some(core::cmp::Ordering::Less) => true,
            Some(core::cmp::Ordering::Greater) => false,
            // Equal is the coincident-endpoint case the windows check above already refused, and
            // `None` is a width decline.
            _ => return Err(undecided()),
        };
        let whole = combinatorics::edge_dir(
            jd,
            cyls,
            wc,
            sg.wall,
            (sg.end[0], sg.end_h[0]),
            (sg.end[1], sg.end_h[1]),
        )?
        .sense()
        .ok_or_else(undecided)?;
        let sense = if forward { whole } else { -whole };
        for w in order.windows(2) {
            let (a, b) = (&keyed[w[0]], &keyed[w[1]]);
            out_segs.push(MergedSeg {
                wall: sg.wall,
                end: [a.1, b.1],
                end_h: [a.2, b.2],
                merged: sg.merged.clone(),
                sense: Some(sense),
            });
        }
    }
    Ok(out_segs)
}

/// **Where a segment crosses a circle**, exactly — the points a `Vertex::Pierce` names.
///
/// The segment rides `wc ∩ sg.wall` and the circle is `cylinder ∩ wc`, so a crossing is
/// `plane ∩ plane ∩ cylinder` — the very shape [`nacre_exact::quad::plane_plane_cylinder`]
/// answers and
/// [`nacre_topo::Vertex::Pierce`] names. Solving along the segment instead would be shorter and
/// would yield a point with **no name**, which the next stage (splitting the circle into arcs)
/// would have to re-derive.
///
/// **It asks the judge's one door for how a plane pair's line meets a cylinder** ([`Judge::line_meet`]),
/// and only where both planes have a world statement: the arc split orders its points about the
/// circle exactly and the emitter reads the class's world statement, so a crossing the interval
/// road names has nowhere to go. On that road the door's answer is the exact meet's.
///
/// **Only a crossing pair, a tangency or a miss of [`nacre_judge::cylinder::LineMeet`]'s shapes can
/// arrive here (or an overflow), and that is provable.** A plane cuts a cylinder in a *circle* only
/// when it is ⊥ to the axis, so `wc ∩ sg.wall` is always ⊥ to the axis too, and the other shapes
/// cannot occur:
/// - `CoincidentPlanes`/`ParallelPlanes` — the segment lies on that line, so the line exists;
/// - `OnRuling` — a line on the cylinder and on `wc` would have to lie inside `cylinder ∩ wc`,
///   which is a circle, and a circle contains no line;
/// - `AxisParallelMiss` — a line ⊥ to the axis is not ∥ to it.
///
/// The last two are the gate's promise that no class oblique to the axis carries a circle
/// ([`class_carries_circle`]); if it breaks, the answer would be a silent «no
/// crossing», so it is refused by that promise's name (`ObliqueCircleClass`) where a debug build
/// stops. The quadratic cannot degenerate either: its leading coefficient is `|d|² > 0` because
/// `d ⊥ m`.
///
/// `Ok(None)` is a plane with no narrow world statement, or overflow (the caller refuses by width).
///
/// **It returns the crossings' *names*, and the coordinate is derived from the name** (see
/// [`combinatorics::pierce_point`]). The pair is solved in this function's own call order —
/// `(wc, sg.wall)` — and `NodeId::pierce` puts it in canonical order, restating the root with it.
/// Solving in ascending order instead would make the correspondence true by construction and leave
/// the canonicalization unexercised, which is where a wrong rule hides; the arc split, walking
/// segments in DCEL order, does not have that luxury either.
/// **It names where the *line* crosses, and the caller says which of those are on the
/// segment.** Both roots come back: dropping the ones outside takes a fence through each
/// endpoint, and a fence is a *plane*, which an end a cylinder pins does not have. Extent is an
/// ordering question, so it belongs with the rule that answers ordering for both kinds of end
/// ([`combinatorics::closed_contains`]) — and the caller already post-filters here anyway, for the
/// tangency.
fn circle_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    circ: &MergedCircle,
    def: &nacre_topo::CylinderDef,
    sg: &MergedSeg,
) -> Result<Option<Vec<combinatorics::NodeId>>, BoolError> {
    use crate::combinatorics::QuadRoot;
    use nacre_judge::cylinder::LineMeet;
    if jd.world_rat(wc).is_none() || jd.world_rat(sg.wall).is_none() {
        return Ok(None);
    }
    let roots = match jd.line_meet(wc, sg.wall, def, false) {
        LineMeet::Pair(_) => vec![QuadRoot::Lo, QuadRoot::Hi],
        // `Double`, not `Lo`: the two roots coincide, so a re-sort must leave the name alone.
        LineMeet::Tangent(..) => vec![QuadRoot::Double],
        // **Reachable, and the honest answer.** The caller's `segment_meets_cylinder` filter (a
        // segment that meets the solid cylinder has a line that meets its surface) needs both
        // endpoints' coordinates, so it is skipped where one end is a pierce point, and there the
        // line genuinely can miss. "No crossings" is what a miss means.
        // Measured: **6** times over the suite.
        LineMeet::Miss => return Ok(Some(Vec::new())),
        LineMeet::Wide | LineMeet::Overflow => return Ok(None),
        shape @ (LineMeet::OnRuling(_) | LineMeet::AxisParallelMiss) => {
            debug_assert!(
                false,
                "a circle's class is ⊥ to the axis, so its meet with any wall is ⊥ to the axis and \
                 can only miss, touch or cross the cylinder — got {shape:?}"
            );
            return Err(reject(RejectReason::ObliqueCircleClass));
        }
        other => unreachable!(
            "a segment lies on the line of its two planes, and both are named — got {other:?}"
        ),
    };
    Ok(Some(
        roots
            .into_iter()
            .map(|root| combinatorics::NodeId::pierce(wc, sg.wall, circ.cyl, root))
            .collect(),
    ))
}
/// **Where a segment crosses a cylinder's rulings on a ∥ class** — the straight sibling of
/// [`circle_crossings`], same machinery ([`nacre_exact::quad::plane_plane_cylinder`],
/// the same fences, the same names), different reachable arms: on a class **parallel** to the axis
/// the meet can be tangent to the lateral — **skipped**, because a touch divides nothing and
/// cutting there would make a zero-length piece — or lie on it: a segment at an angle no lateral
/// face holds crosses nothing, and one a face holds the ruling split refuses by name rather than
/// reasoning it unreachable.
///
/// Returns each crossing with the meet it was solved as ([`combinatorics::PiercePoint`] — exact, or
/// an interval box where a class has no world statement; the side and extent tests below read it
/// without re-solving). `None` where the exact road cannot say: the wall has no narrow world
/// statement, a `Wide` name, an overflow, or a ruling whose angle has no rational spelling.
fn lateral_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    sg: &MergedSeg,
    on_a_shared_line: bool,
    pieces: impl FnOnce() -> Vec<crate::planes::Footprint>,
) -> Result<Option<Vec<LateralCrossing>>, BoolError> {
    use crate::combinatorics::QuadRoot;
    use nacre_judge::cylinder::LineMeet;
    let Some(v) = jd.world_rat(sg.wall) else {
        return Ok(None);
    };
    let m = def.dir();
    // **The line runs along the axis when both normals are perpendicular to it.** The wall's is
    // stated exactly; so is `wc`'s where it has a narrow world name. A class with no world statement
    // here sections the cylinder — it carries rulings, so the gate found it along the axis
    // (`Silence::Sectioned`). A `Wide` name says neither, and the exact road refuses it for its
    // width. Along the axis the line crosses no ruling: it misses the surface, or lies on it, which
    // no record states here. The interval road cannot see that parallelism from the planes.
    let perpendicular = |w: &[Rat; 4]| {
        nacre_exact::dot_sign_rat(&[w[0], w[1], w[2]], &m) == nacre_exact::Orient::Zero
    };
    let along_axis = perpendicular(&v)
        && match jd.planes[wc].world {
            None => true,
            Some(_) => jd.world_rat(wc).is_some_and(|w| perpendicular(&w)),
        };
    let undecided = || reject(RejectReason::CylinderGateUndecided);
    let crossing = |root: QuadRoot, at: combinatorics::PiercePoint| {
        (combinatorics::NodeId::pierce(wc, sg.wall, cyl, root), at)
    };
    // A point of the segment the interval road can box: a three-plane end.
    let end_point = sg.end.iter().find_map(|&e| match e.kind() {
        NodeKind::ThreePlane(t) => Some(nacre_judge::cylinder::AxisLine::Through(
            nacre_judge::cylinder::CylPoint::Three(t),
        )),
        NodeKind::Pierce { .. } => None,
    });
    match jd.line_meet(wc, sg.wall, def, along_axis) {
        LineMeet::Pair(roots) => {
            let mut out = Vec::with_capacity(2);
            for root in [QuadRoot::Lo, QuadRoot::Hi] {
                let at = roots.point(jd, cyl, def, root).ok_or_else(undecided)?;
                out.push(crossing(root, at));
            }
            Ok(Some(out))
        }
        // A **tangent** class touches the cylinder along one ruling: the segment's line meets it in
        // one point, root `Double` — the crossing the seated face's own tangent piece (`side == 0`)
        // is cut at, named as the scan names it.
        LineMeet::Tangent(line, s) => Ok(Some(vec![crossing(
            QuadRoot::Double,
            combinatorics::PiercePoint::Exact(line, nacre_exact::quad::QuadVal::from_rat(s)),
        )])),
        LineMeet::Miss | LineMeet::AxisParallelMiss => Ok(Some(Vec::new())),
        // The segment's line lies on the cylinder: a line the gate recorded
        // ([`crate::planes::SharedRuling`]) crosses no ruling — it *is* one, and the fold in
        // [`split_rulings`] states it there.
        LineMeet::OnRuling(_) if on_a_shared_line => Ok(Some(Vec::new())),
        LineMeet::OnRuling(line) => {
            match crate::planes::line_off_faces(jd, def, Some(&line), end_point, &pieces()) {
                Some(true) => Ok(Some(Vec::new())),
                // On a face: the ruling split does not arrange the shape yet — nor where the planes
                // degenerate.
                Some(false) => Err(reject(RejectReason::RulingBoundNotYet)),
                None => Ok(None),
            }
        }
        LineMeet::CoincidentPlanes | LineMeet::ParallelPlanes => {
            Err(reject(RejectReason::RulingBoundNotYet))
        }
        // Along the axis, the interval road cannot tell a line on the surface from one beside it;
        // either way it crosses no ruling, and it touches a face only if it lies at an angle one
        // holds — the question [`crate::planes::line_off_faces`] answers from a point of the segment.
        LineMeet::Undecided
            if along_axis
                && crate::planes::line_off_faces(jd, def, None, end_point, &pieces())
                    == Some(true) =>
        {
            Ok(Some(Vec::new()))
        }
        LineMeet::Undecided => Err(undecided()),
        LineMeet::Wide | LineMeet::Overflow => Ok(None),
    }
}

/// One crossing [`lateral_crossings`] found: its Pierce name, and the meet it was solved as (the
/// side and extent tests read it without re-solving).
type LateralCrossing = (combinatorics::NodeId, combinatorics::PiercePoint);

/// `x · y`, checked.
fn dot(x: &[Rat; 3], y: &[Rat; 3]) -> Option<Rat> {
    x[0].checked_mul(y[0])?
        .checked_add(x[1].checked_mul(y[1])?)?
        .checked_add(x[2].checked_mul(y[2])?)
}

/// **A point's coordinate along a cylinder's axis** (`m · x`), for either kind of name — the one
/// reader the ruling split and the fold both order by. Exact where the point states it: three
/// planes with exact coordinates, a plane of the name across the axis, a pierce point's exact meet
/// (`a + b√c`); else the point itself, whose interval box the judge reads ([`Judge::axis_order`]).
/// **Exact comes first** because an interval cannot prove two coordinates equal — two ends on one
/// cap must compare as `Equal`, which only a stated value can say.
#[derive(Clone)]
enum AxisCoord {
    Exact(nacre_exact::quad::QuadVal),
    Pierce(Box<combinatorics::PiercePoint>),
    Three([usize; 3]),
}

impl AxisCoord {
    /// The coordinate of point `n` along cylinder `cyl`'s axis. A pierce point is solved on its own
    /// name's cylinder, which need not be `cyl` (a segment's end where another cylinder cut it).
    fn of(
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        cyl: usize,
        n: combinatorics::NodeId,
    ) -> Option<Self> {
        let def = &cyls[cyl].def;
        match n.kind() {
            NodeKind::ThreePlane(t) => {
                if let Some(p) = combinatorics::node_coords_rat(jd, n) {
                    return Some(Self::Exact(nacre_exact::quad::QuadVal::from_rat(dot(
                        &p,
                        &def.dir(),
                    )?)));
                }
                Some(Self::across_axis(jd, def, &t)?.unwrap_or(Self::Three(t)))
            }
            NodeKind::Pierce { cyl: own, .. } => Self::of_meet(
                jd,
                def,
                n,
                &combinatorics::pierce_meet(jd, own, &cyls[own].def, n)?,
            ),
        }
    }

    /// The coordinate of pierce point `n` already solved as `meet`.
    fn of_meet(
        jd: &Judge<'_, WorkingPlane>,
        def: &nacre_topo::CylinderDef,
        n: combinatorics::NodeId,
        meet: &combinatorics::PiercePoint,
    ) -> Option<Self> {
        if meet.exact().is_some() {
            return Some(Self::Exact(meet.stated_axis(def)?));
        }
        let (planes, _, _) = combinatorics::pierce_name(n)?;
        Some(
            Self::across_axis(jd, def, &planes)?
                .unwrap_or_else(|| Self::Pierce(Box::new(meet.clone()))),
        )
    }

    /// A plane of the name across the axis holds every point of it at one axis coordinate,
    /// `m·o + t·|m|²` with `t` where the axis meets it — stated exactly, so two points on one cap
    /// compare as equal. `Some(None)` where no plane of the name is across the axis; `None` on
    /// overflow.
    fn across_axis(
        jd: &Judge<'_, WorkingPlane>,
        def: &nacre_topo::CylinderDef,
        planes: &[usize],
    ) -> Option<Option<Self>> {
        let m = def.dir();
        for &c in planes {
            if let Some(w) = jd.world_rat(c)
                && class_carries_circle(&[w[0], w[1], w[2]], &m)
            {
                let t = crate::planes::axis_param_of_plane(&w, def)?;
                return Some(Some(Self::Exact(nacre_exact::quad::QuadVal::from_rat(
                    dot(&def.origin(), &m)?.checked_add(t.checked_mul(dot(&m, &m)?)?)?,
                ))));
            }
        }
        Some(None)
    }

    fn value(&self) -> nacre_judge::cylinder::AxisValue<'_> {
        use nacre_judge::cylinder::{AxisValue, CylPoint};
        match self {
            Self::Exact(q) => AxisValue::Stated(q),
            Self::Pierce(p) => AxisValue::Point(CylPoint::Pierce(p)),
            Self::Three(t) => AxisValue::Point(CylPoint::Three(*t)),
        }
    }

    fn cmp(
        &self,
        other: &Self,
        jd: &Judge<'_, WorkingPlane>,
        def: &nacre_topo::CylinderDef,
    ) -> Option<core::cmp::Ordering> {
        jd.axis_order(def, self.value(), other.value())
    }
}

/// `split_rulings`' product — `None` when nothing crossed.
type SplitRulings = Option<(Vec<MergedSeg>, Vec<MergedRuling>)>;

/// **Cut every ruling a segment crosses into pieces, and the segments with it** — the straight
/// twin of [`split_circles`](fn@split_circles), run after it.
///
/// **It cuts rulings against *segments*, and deliberately not against arcs** — as its twin cuts
/// circles against segments and not against rulings. One class *can* carry both populations
/// (there is no blanket refusal that says otherwise), so what keeps the pair of splits
/// complete is not «they never share a class» but the sharper fact that **a class's edge is
/// some face's boundary**: two edges of different cylinders meeting would put two lateral faces on
/// one point, which `crate::planes::lateral_faces_clear` denies. That one sentence covers
/// every pair no split cuts — circle×circle and ruling×ruling as well as this one — and
/// `ClassEdges::of` carries a shipped backstop over the pair that has an exact predicate
/// ([`crate::RejectReason::CircleCrossesRuling`]).
///
/// The ordering vocabulary is the established one: along a segment,
/// [`nacre_exact::QuadVal::cmp_value`] on the canonical meet line; along a ruling, the same
/// comparator on the **axis coordinate** (`m · x`, evaluated exactly from each node's re-solved
/// name).
pub(super) fn split_rulings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    wc: usize,
    segs: &[MergedSeg],
    rulings: &[MergedRuling],
    aliases: &Aliases,
) -> Result<SplitRulings, BoolError> {
    // A class with no world statement (it sections the cylinder) refuses by the gate's name.
    let undecided = || reject(combinatorics::undecided_reason(cyls, &[wc]));
    if rulings.is_empty() {
        return Ok(None);
    }
    // Each ruling's end coordinates, re-solved from the names once.
    let ruling_ends: Vec<[AxisCoord; 2]> = rulings
        .iter()
        .map(|r| {
            let coord = |n| AxisCoord::of(jd, cyls, r.cyl, n);
            match (coord(r.end[0]), coord(r.end[1])) {
                (Some(a), Some(b)) => Ok([a, b]),
                _ => Err(undecided()),
            }
        })
        .collect::<Result<_, BoolError>>()?;
    // Distinct cylinders carrying rulings here, each solved once per segment.
    let mut cyl_list: Vec<usize> = Vec::new();
    for r in rulings {
        if !cyl_list.contains(&r.cyl) {
            cyl_list.push(r.cyl);
        }
    }
    // The lines on this class that lie on a cylinder and that this class states as a ruling —
    // `(the other class, cylinder, which of this class's two rulings)`, the pair read off the
    // gate's record, never decided here, and the side asked of the line once
    // ([`combinatorics::shared_line`], which the straight-angle cleaning reads too): it is the
    // identity a ruling piece carries (`MergedRuling::side`).
    let mut shared_here: Vec<(usize, usize, i8)> = Vec::new();
    for (ci, cy) in cyls.iter().enumerate() {
        for sr in cy.shared.iter().filter(|sr| sr.stated_by.contains(&wc)) {
            let l = combinatorics::shared_line(jd, cyls, ci, sr, wc).ok_or_else(undecided)?;
            shared_here.push((l.other, ci, l.side));
        }
    }
    let mut on_ruling: Vec<Vec<NodeId>> = vec![Vec::new(); rulings.len()];
    let (mut on_seg, cuts) = collect(
        jd,
        cyls,
        wc,
        segs,
        &cyl_list,
        |_si, sg, cyl| {
            let shared = shared_here
                .iter()
                .any(|&(t, c, _)| t == sg.wall && c == cyl);
            let pieces = || crate::planes::lateral_pieces(faces, cyls[cyl].surf);
            let Some(xs) = lateral_crossings(jd, wc, cyl, &cyls[cyl].def, sg, shared, pieces)?
            else {
                return Err(reject(combinatorics::undecided_reason(cyls, &[wc])));
            };
            Ok(xs.into_iter().map(|(n, meet)| (n, Some(meet))).collect())
        },
        // Which of `wc`'s rulings the crossing is on, `0` being `wc`'s tangent one — asked of the
        // point ([`ruling_side_signed`]), not read off the crossing's own root (`Double`, i.e. the
        // *line* `wc ∩ wall` touches the cylinder), which is the same statement only while `wc` is
        // the wall that touches: a line tangent to the cylinder can cross a `wc` that cuts it, and
        // then the root's `0` matches no ruling of `wc` at all. Then which piece of that ruling
        // (several faces of one solid can put several pieces on one side), by the closed extent of
        // the axis coordinate; beyond every piece, or no piece on that side, is off.
        |cyl, n, meet| {
            let def = &cyls[cyl].def;
            let undecided = || Landing::Fails(reject(combinatorics::undecided_reason(cyls, &[wc])));
            let meet = meet.expect("a ruling crossing carries the meet it was solved as");
            let Some(side) = meet.ruling_side(jd, wc, def) else {
                return undecided();
            };
            let Some(c) = AxisCoord::of_meet(jd, def, n, meet) else {
                return undecided();
            };
            let (mut cut, mut landed) = (None, false);
            for (ri, r) in rulings.iter().enumerate() {
                if r.cyl != cyl || r.side != side {
                    continue;
                }
                let [lo, hi] = &ruling_ends[ri];
                let (Some(a), Some(b)) = (c.cmp(lo, jd, def), c.cmp(hi, jd, def)) else {
                    return undecided();
                };
                use core::cmp::Ordering::*;
                match (a, b) {
                    // Strictly inside: the ruling splits here, and the segment with it. Pieces on
                    // one side do not overlap, so at most one holds the point inside.
                    (Greater, Less) => {
                        debug_assert!(cut.is_none(), "two ruling pieces on one side overlap");
                        cut = Some(ri);
                        landed = true;
                    }
                    // On an end: the T-junction. One point, and it must wear **one** name — the
                    // crossing's `{wc, wall}` pair must be the very pair that named the ruling's
                    // end, or the table must know the two names as one point — else two names
                    // share the point.
                    (Equal, _) | (_, Equal) => {
                        let end_n = if a == Equal { r.end[0] } else { r.end[1] };
                        if aliases.canon_point(end_n) != aliases.canon_point(n) {
                            return Landing::Fails(reject(RejectReason::CoincidentNodes));
                        }
                        landed = true;
                    }
                    _ => {} // beyond this piece: the surface continues, the face does not
                }
            }
            if landed {
                Landing::On {
                    cut,
                    node: aliases.canon_point(n),
                }
            } else {
                Landing::Off
            }
        },
    )?;
    for (ri, n) in cuts {
        on_ruling[ri].push(n);
    }
    let folds = segs
        .iter()
        .any(|sg| shared_here.iter().any(|&(t, _, _)| t == sg.wall));
    if on_ruling.iter().all(|v| v.is_empty()) && on_seg.iter().all(|v| v.is_empty()) && !folds {
        return Ok(None);
    }
    // ---- segments → sub-segments, the arc split's own idiom ----
    let out_segs = split_segments_at(jd, cyls, wc, segs, &mut on_seg, aliases)?;
    // ---- rulings → pieces, in axis order ----
    let mut out_rulings = Vec::with_capacity(rulings.len());
    for (ri, r) in rulings.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_ruling[ri]);
        if nodes.is_empty() {
            out_rulings.push(r);
            continue;
        }
        nodes.sort_unstable();
        nodes.dedup();
        let mut keyed: Vec<(AxisCoord, combinatorics::NodeId)> = Vec::new();
        for &n in &nodes {
            keyed.push((AxisCoord::of(jd, cyls, r.cyl, n).ok_or_else(undecided)?, n));
        }
        let [lo, hi] = ruling_ends[ri].clone();
        keyed.push((lo, r.end[0]));
        keyed.push((hi, r.end[1]));
        let order = order_distinct(
            keyed.len(),
            |i, j| keyed[i].0.cmp(&keyed[j].0, jd, &cyls[r.cyl].def),
            undecided,
        )?;
        // Ascending axis coordinate **is** the stored convention (`end[0] → end[1]` ascends
        // `+m`), so each window is already a well-formed piece.
        for pair in order.windows(2) {
            out_rulings.push(MergedRuling {
                cyl: r.cyl,
                side: r.side,
                end: [keyed[pair[0]].1, keyed[pair[1]].1],
                merged: r.merged.clone(),
                #[cfg(test)]
                orient: r.orient,
            });
        }
    }
    let out_segs = if folds {
        fold_onto_rulings(jd, cyls, out_segs, &mut out_rulings, &shared_here, aliases)?
    } else {
        out_segs
    };
    Ok(Some((out_segs, out_rulings)))
}

/// **A segment on a line the gate recorded is that ruling's edge** — so it is stated once, in
/// the ruling vocabulary: its `(solid, kind)` contributions join the ruling piece with the same
/// two ends, and the segment goes. The chart's vertical lines come from rulings alone, and the
/// lateral continues across such a line (the prism's corner is the *other* solid's edge), so the
/// station has to be there; stated as a plane segment as well, one line would reach the walk as
/// two edges.
///
/// The splits above cut both at every point the other crosses, so a piece either matches a
/// ruling piece end for end or shares no stretch with one. No stretch shared is a segment off the
/// lateral — past its cap, or where the ruling is silent because the line is the lateral face's
/// own edge (`ruling_sweep`'s `quiet_nodes`: there the plane vocabulary states it, the other half
/// of the one rule) — and it stays a segment. A stretch shared without a match is a shape the
/// splits did not cut, and is refused by name rather than stated twice.
///
/// `Seated`, `Graze` and `Tangent` say which side of the class plane a face fills, so they mean
/// the same on a ruling. `Transversal` says the face crosses the class plane along the line — a
/// fact about the face and the plane, whichever way the line runs — and that is what
/// `cells::edge_mask` reads of it: both of its solid's bits flip. Its `mat` stays the contributing
/// face's, a side of the directed line `wc × wall` and not of the ruling's `+m`; read on a ruling
/// piece as the latter it would be wrong. A second mark of the same solid on the piece needs a
/// plane face whose interior holds its own lateral's ruling where the lateral is — two faces of
/// one valid solid meeting inside, which no valid operand has (`edge_mask` refuses two transversals
/// or a `Graze` beside one, `EdgeOccupancyConflict`, and takes a `Seated` or `Tangent` beside one as
/// seated). Which side of
/// the station holds the other solid the chart reads off the cells, not off these marks.
///
/// A plane face carried here may still touch the lateral's own edge from inside the same result
/// solid — the lateral ending on a flat face where another part of the body crosses that line. That
/// is the self-touch verdict's to refuse (`assembly::self_touch_reject`, which asks a lateral's
/// ruling edge against the solid's plane faces).
fn fold_onto_rulings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    segs: Vec<MergedSeg>,
    rulings: &mut [MergedRuling],
    shared_here: &[(usize, usize, i8)],
    aliases: &Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let canon = |e: [combinatorics::NodeId; 2]| {
        let mut x = [aliases.canon_point(e[0]), aliases.canon_point(e[1])];
        x.sort_unstable();
        x
    };
    let mut kept = Vec::with_capacity(segs.len());
    for sg in segs {
        let Some(&(_, cyl, side)) = shared_here.iter().find(|&&(t, _, _)| t == sg.wall) else {
            kept.push(sg);
            continue;
        };
        let key = canon(sg.end);
        if let Some(r) = rulings
            .iter_mut()
            .find(|r| r.cyl == cyl && canon(r.end) == key)
        {
            r.merged.extend(sg.merged.iter().copied());
            continue;
        }
        // No piece end for end: off the lateral, or a stretch the splits left uncut.
        let def = &cyls[cyl].def;
        let coord = |n| AxisCoord::of(jd, cyls, cyl, n);
        let (Some(a), Some(b)) = (coord(sg.end[0]), coord(sg.end[1])) else {
            return Err(undecided());
        };
        let (a, b) = match a.cmp(&b, jd, def).ok_or_else(undecided)? {
            core::cmp::Ordering::Greater => (b, a),
            _ => (a, b),
        };
        for r in rulings.iter().filter(|r| r.cyl == cyl && r.side == side) {
            let (Some(lo), Some(hi)) = (coord(r.end[0]), coord(r.end[1])) else {
                return Err(undecided());
            };
            use core::cmp::Ordering::Less;
            if a.cmp(&hi, jd, def).ok_or_else(undecided)? == Less
                && lo.cmp(&b, jd, def).ok_or_else(undecided)? == Less
            {
                return Err(reject(RejectReason::RulingBoundNotYet));
            }
        }
        kept.push(sg);
    }
    Ok(kept)
}

/// [`collect`]'s product: each segment's cut points, and each `(piece, node)` a crossing cuts a piece
/// of the kind at, in the order the walk found them.
type Collected = (Vec<Vec<NodeId>>, Vec<(usize, NodeId)>);

/// A crossing, by name, with the meet it was solved as where its kind reads it again.
type Crossing = (NodeId, Option<combinatorics::PiercePoint>);

/// **Every crossing a segment holds, for one kind of curve** — the one walk pass B makes over a
/// class's segments, whatever it cuts them against. The kind supplies only its own answers: the
/// cylinders it has curves of here (`here`), where a segment's line crosses one of them
/// (`crossings`, each crossing's name and the meet it was solved as where the kind reads it again),
/// and where a crossing lands among its pieces (`land`, a [`Landing`]). The cut points left on each
/// segment come back (`on_seg[si]`), with each `(piece, node)` a crossing cuts a piece of the kind at,
/// in the order the walk found them.
///
/// A segment whose two ends are one point cuts nothing — asked by **name**, which is the identity:
/// two spellings of one point are refused upstream, not tolerated here (measured unexercised). A
/// segment on a class the gate proved clear of a cylinder cannot reach its curves — the gate's
/// record, asked before anything the segment's class would have to describe.
///
/// **Whether a crossing is on the segment is asked here, in the one vocabulary that can answer it
/// for both kinds of end** ([`combinatorics::closed_contains`]). Two plane fences would need a
/// *plane* through each endpoint, and so could only be built where every end is three-plane named.
/// No segment is skipped for an end a cylinder pins: «such a segment shares only endpoints with a
/// curve» is an argument about one narrow population — the sort of claim that goes quietly false
/// when the population widens — so it is asked instead of assumed. The segment's two ends are
/// located on its own line when a crossing first asks, not for every segment: a segment no
/// crossing reaches is never asked, and its class may have no description to locate it by.
///
/// `land` is asked before the segment's question, so a crossing **shown off every piece** of the kind
/// on this class leaves first. Only a proof counts: a landing the pieces cannot judge waits for the
/// segment's answer and is refused only where the segment holds the crossing.
fn collect(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    here: &[usize],
    crossings: impl Fn(usize, &MergedSeg, usize) -> Result<Vec<Crossing>, BoolError>,
    land: impl Fn(usize, NodeId, Option<&combinatorics::PiercePoint>) -> Landing,
) -> Result<Collected, BoolError> {
    // A class with no world statement (it sections the cylinder) refuses by the gate's name.
    let undecided = || reject(combinatorics::undecided_reason(cyls, &[wc]));
    let mut on_seg: Vec<Vec<NodeId>> = vec![Vec::new(); segs.len()];
    let mut cuts: Vec<(usize, NodeId)> = Vec::new();
    for (si, sg) in segs.iter().enumerate() {
        if sg.end[0] == sg.end[1] {
            continue;
        }
        let mut ends: Option<[combinatorics::OnLine; 2]> = None;
        for &cyl in here {
            if cyls.get(cyl).is_some_and(|cy| cy.is_silent(sg.wall)) {
                continue;
            }
            for (n, meet) in crossings(si, sg, cyl)? {
                // **A proof of absence outranks an undecided question.** A crossing shown off every
                // piece of its kind cuts nothing, whatever the segment says; asked after, the
                // segment's undecided answer would refuse a point no edge of this class has.
                let landing = land(cyl, n, meet.as_ref());
                if matches!(landing, Landing::Off) {
                    continue;
                }
                let at =
                    combinatorics::locate(jd, cyls, wc, sg.wall, combinatorics::PointOn::Pierce(n))
                        .ok_or_else(undecided)?;
                if ends.is_none() {
                    let mut it = (0..2).map(|k| {
                        combinatorics::on_line(
                            jd,
                            cyls,
                            wc,
                            sg.wall,
                            Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
                        )
                        .ok_or_else(undecided)
                    });
                    ends = Some([it.next().unwrap()?, it.next().unwrap()?]);
                }
                let located = ends.as_ref().expect("located just above");
                if !combinatorics::closed_contains(jd, wc, sg.wall, &at, located)
                    .ok_or_else(undecided)?
                {
                    continue;
                }
                match landing {
                    Landing::Off => {}
                    Landing::Fails(e) => return Err(e),
                    Landing::On { cut, node } => {
                        on_seg[si].push(node);
                        cuts.extend(cut.map(|c| (c, node)));
                    }
                }
            }
        }
    }
    Ok((on_seg, cuts))
}

/// **Where a crossing lands among a kind's pieces**, judged once, before the segment is asked
/// ([`collect`]).
enum Landing {
    /// Shown off every piece of the kind on this class: it cuts nothing.
    Off,
    /// On a piece: `node` (the crossing under the table's name) cuts the segment, and the piece
    /// `cut` — an index into the kind's own list — where it lands strictly inside one.
    On { cut: Option<usize>, node: NodeId },
    /// The pieces cannot place it — refused, but only once the segment is shown to hold it: a
    /// crossing off the segment is no point of this class, and refusing it would refuse nothing
    /// the class has.
    Fails(BoolError),
}

/// **The order of distinct points along one curve**, by `cmp` over their indices `0..n`.
///
/// The arrangement keys vertices by name, so two names at one place would make two vertices where
/// there is one: an `Equal` between neighbours is refused (`CoincidentNodes`) — names the alias
/// table knows as one point were folded before the order is asked, so what reaches here is a
/// coincidence no discovery event recorded, and there is no ground to invent the identity on. A
/// pair that cannot be told apart is `undecided`'s, and so is a failure the sort met: it may have
/// asked the pair the other way round, and a failure is no answer.
fn order_distinct(
    n: usize,
    cmp: impl Fn(usize, usize) -> Option<core::cmp::Ordering>,
    undecided: impl Fn() -> BoolError,
) -> Result<Vec<usize>, BoolError> {
    let mut order: Vec<usize> = (0..n).collect();
    combinatorics::try_sort_by(&mut order, |&i, &j| cmp(i, j).ok_or(()))
        .map_err(|()| undecided())?;
    for w in order.windows(2) {
        match cmp(w[0], w[1]) {
            Some(core::cmp::Ordering::Equal) => return Err(reject(RejectReason::CoincidentNodes)),
            Some(_) => {}
            None => return Err(undecided()),
        }
    }
    Ok(order)
}
