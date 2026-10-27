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
    let undecided = || reject(RejectReason::WitnessNotRational);
    // Nothing to cross. The ends below are built for the whole slice at once, so the empty case has
    // to be answered before them or an empty class would start declining. `split_rulings` opens
    // with the same line. Measured: **61993** times over the suite — this is the common case, not
    // an edge one.
    if circles.is_empty() {
        return Ok(None);
    }
    // **Both ends of a segment, located on that segment's own line, once** — the hoist
    // `split_at_crossings` makes, for the same reason and with the same obligation: each pair is
    // located for `(wc, sg.wall)`, which is the pair [`combinatorics::closed_contains`] must be
    // called with below. Located when a crossing first asks, not for every segment: a segment on a
    // class the gate proved silent about the circle's cylinder is never asked, and its class may
    // have no description to locate it by.
    let locate_ends = |s: &MergedSeg| -> Result<[combinatorics::OnLine; 2], BoolError> {
        let mut it = (0..2).map(|k| {
            combinatorics::on_line(
                jd,
                cyls,
                wc,
                s.wall,
                Split::of(s.end[k], s.end_h[k]).on(wc, s.wall),
            )
            .ok_or_else(undecided)
        });
        Ok([it.next().unwrap()?, it.next().unwrap()?])
    };
    let mut seg_ends: Vec<Option<[combinatorics::OnLine; 2]>> =
        (0..segs.len()).map(|_| None).collect();
    // Which nodes land on each circle, and which on each segment. Collected together because one
    // crossing is a point of both — splitting only one of them would leave the other's edge running
    // through a vertex it does not have.
    let canon = |x: NodeId| aliases.canon_point(x);
    let mut on_circle: Vec<Vec<NodeId>> = vec![Vec::new(); circles.len()];
    let mut on_seg: Vec<Vec<NodeId>> = vec![Vec::new(); segs.len()];
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
    for (ci, circ) in circles.iter().enumerate() {
        let def = &cyls[circ.cyl].def;
        let (o, m, r2) = (def.origin(), def.dir(), def.r2());
        for (si, sg) in segs.iter().enumerate() {
            // A segment whose two ends are one point cuts nothing. Asked by **name**, which is the
            // identity: two spellings of one point are refused upstream, not tolerated here.
            // Measured unexercised.
            if sg.end[0] == sg.end[1] {
                continue;
            }
            // A segment on a class the gate proved clear of this cylinder misses its circle — the
            // gate's record, asked before anything the segment's class would have to describe.
            if cyls.get(circ.cyl).is_some_and(|cy| cy.is_silent(sg.wall)) {
                continue;
            }
            // **The cheap rejection runs where it can be formed, and skipping it asks
            // *more*, not less.** `segment_meets_cylinder` needs both endpoints' coordinates; where
            // one is a pierce point the question goes straight to `circle_crossings`, which answers
            // `Miss` exactly when the line misses. Losing the filter costs a solve, never an answer.
            if let [Some(p0), Some(p1)] = &seg_coords[si]
                && (p0 == p1 || !nacre_exact::segment_meets_cylinder(p0, p1, &o, &m, r2))
            {
                continue;
            }
            let Some(xs) = circle_crossings(jd, wc, circ, def, sg) else {
                return Err(undecided());
            };
            for n in xs {
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
                // not gated on `crossings`, so a ⊥ class takes the same road whatever the gate
                // recorded. 15 of 21 measured cells assemble with it in place (the other 6 hold a
                // third plane on the tangent line and are refused earlier, by `CoincidentNodes`).
                if matches!(
                    combinatorics::pierce_name(n),
                    Some((_, _, crate::combinatorics::QuadRoot::Double))
                ) {
                    continue;
                }
                // **Whether the crossing is on this segment is asked here, in the one
                // vocabulary that can answer it for both kinds of end.** Two plane fences would
                // need a *plane* through each endpoint, and so could only be built where every
                // end is three-plane named.
                let at =
                    combinatorics::locate(jd, cyls, wc, sg.wall, combinatorics::PointOn::Pierce(n))
                        .ok_or_else(undecided)?;
                if seg_ends[si].is_none() {
                    seg_ends[si] = Some(locate_ends(sg)?);
                }
                let ends = seg_ends[si].as_ref().expect("located just above");
                if !combinatorics::closed_contains(jd, wc, sg.wall, &at, ends)
                    .ok_or_else(undecided)?
                {
                    continue;
                }
                // The crossing's name, as the table knows it: a class through a
                // corner names the corner again here, and the tables key by name.
                let n = canon(n);
                on_circle[ci].push(n);
                on_seg[si].push(n);
            }
        }
    }

    // Nothing crossed: the caller keeps its own slices and no edge is copied. The common case
    // pays for the question and not for an answer it does not need.
    if on_circle.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    let mut out_circles = Vec::new();
    let mut arcs = Vec::new();
    let mut cut_rims: Vec<(usize, CutRim)> = Vec::new();

    // **The circle side runs first, because it decides which crossings are cut points
    // at all**. `circle_crossings` names where the segment's *line* meets the whole
    // circle; but the class's circle is only where its arcs are (`MergedCircle::merged` — a
    // fillet's quarter, a slot's half), and a crossing on the circle's **continuation** lies on
    // no edge of this class. Cutting it anyway puts a vertex on the segment where nothing
    // crosses it (measured: 4 per op on the `tangentline slab 1.6` controls, 16 on
    // `rrect-box`, 0 on the other 353 census rows), and where that point *is* a vertex of the
    // segment already (a face through the top of a fillet's circle: the user's gusset foot,
    // `y = 0 = −2 + r`) the segment split refuses it as two names for one point.
    // **The rule: a crossing shown off every arc the class carries is no cut point.** An arc's own
    // ends are cut points by construction (fed in above); every other crossing is asked against
    // each arc's two ends, and one shown off all of them leaves before the order is formed. One
    // that cannot be shown off stays a cut point: off an arc it is a vertex where nothing crosses
    // the segment, which costs a vertex, where dropping a crossing on an arc would lose one.
    //
    // Segments first ("the sharper question before the vaguer one") would be an unmeasured order,
    // and the order changes no refusal: both halves refuse a point wearing two names with the same
    // reject, and identity on both halves is the alias table's (`canon`), so what reaches either
    // refusal is a coincidence no discovery event recorded — the honest floor.
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
        // **A crossing shown off every arc leaves before the order is formed** — asked against
        // each arc's two ends alone ([`combinatorics::theta_between`]), so the points of the
        // circle's continuation are never ordered against each other: two names there for one
        // point would refuse a point no edge of this class has. What this cannot show stays, and
        // the order below refuses it if it is another name of a point already there. Three
        // circles drop nothing here: one with a contribution that is the whole circle (every
        // crossing is on it), one with no contribution (where «off every arc» would be an empty
        // truth), and one with a closed arc (`p == q`), which `theta_between` reads as the whole
        // turn less its end and so answers `false` for another name of the end itself.
        let extents: Option<Vec<[NodeId; 2]>> =
            circ.merged.iter().map(|&(_, _, arc)| arc).collect();
        let mut off_arcs: Vec<NodeId> = Vec::new();
        if let Some(extents) = extents
            && !extents.is_empty()
            && extents.iter().all(|&[p, q]| canon(p) != canon(q))
        {
            nodes.retain(|&x| {
                let off = !is_end(x)
                    && extents.iter().all(|&[p, q]| {
                        combinatorics::theta_between(jd, circ.cyl, &cyls[circ.cyl].def, p, q, x)
                            == Some(false)
                    });
                if off {
                    off_arcs.push(x);
                }
                !off
            });
            for v in on_seg.iter_mut() {
                v.retain(|&n| !off_arcs.iter().any(|&d| canon(d) == canon(n)));
            }
        }
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
        let mut order: Vec<usize> = (0..keyed.len()).collect();
        combinatorics::try_sort_by(&mut order, |&i, &j| cmp(i, j).ok_or(()))
            .map_err(|()| undecided())?;
        // A crossing that lands **on** an endpoint is one point wearing two names — a three-plane
        // one and a pierce one — and the DCEL keys vertices by name, so shipping both would make
        // two vertices where there is one. Names the alias table knows as one point were folded
        // above; what reaches here is a coincidence no discovery event recorded, and
        // refusing it is honest — there is no ground to invent the identity on.
        for w in order.windows(2) {
            match cmp(w[0], w[1]) {
                Some(core::cmp::Ordering::Equal) => {
                    return Err(reject(RejectReason::CoincidentNodes));
                }
                Some(_) => {}
                // The sort may have asked this pair the other way round; a failure is no answer.
                None => return Err(undecided()),
            }
        }
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
/// **Only three of `CylinderMeet`'s variants can arrive here, and that is provable.** A plane
/// cuts a cylinder in a *circle* only when it is ⊥ to the axis, so `wc ∩ sg.wall` is always ⊥ to
/// the axis too:
/// - `CoincidentPlanes`/`ParallelPlanes` — the segment lies on that line, so the line exists;
/// - `OnRuling` — a line on the cylinder and on `wc` would have to lie inside `cylinder ∩ wc`,
///   which is a circle, and a circle contains no line;
/// - `AxisParallelMiss` — a line ⊥ to the axis is not ∥ to it.
///
/// The quadratic cannot degenerate either: its leading coefficient is `|d|² > 0` because `d ⊥ m`.
///
/// `None` is overflow (the caller drops the witness, never the rejection).
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
) -> Option<Vec<combinatorics::NodeId>> {
    use crate::combinatorics::QuadRoot;
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    let w = jd.world_rat(wc)?;
    let v = jd.world_rat(sg.wall)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let roots = match nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2)? {
        CylinderMeet::Pair { s, .. } => vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])],
        // `Double`, not `Lo`: the two roots coincide, so a re-sort must leave the name alone.
        CylinderMeet::Tangent { s, .. } => vec![(QuadRoot::Double, QuadVal::from_rat(s))],
        // **Reachable, and the honest answer.** The caller's `segment_meets_cylinder` filter (a
        // segment that meets the solid cylinder has a line that meets its surface) needs both
        // endpoints' coordinates, so it is skipped where one end is a pierce point, and there the
        // line genuinely can miss. "No crossings" is what a miss means.
        // Measured: **6** times over the suite.
        CylinderMeet::Miss(_) => return Some(Vec::new()),
        other => unreachable!(
            "a circle's class is ⊥ to the axis, so its meet with any wall is ⊥ to the axis and \
             can only miss, touch or cross the cylinder — got {other:?}"
        ),
    };
    Some(
        roots
            .into_iter()
            .map(|(root, _)| combinatorics::NodeId::pierce(wc, sg.wall, circ.cyl, root))
            .collect(),
    )
}
