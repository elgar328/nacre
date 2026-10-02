use super::*;
/// `split_circles`' product: `None` when nothing crossed (the caller keeps its slices), else the
/// split segments, the still-whole circles, the arcs, and each cut circle's [`CutRim`].
type SplitCircles = Option<(
    Vec<MergedSeg>,
    Vec<MergedCircle>,
    Vec<MergedArc>,
    Vec<(usize, CutRim)>,
)>;

/// **Cut every circle a segment crosses into arcs, and the segments with it.**
///
/// ★★★ This is the arc split. What it does *not* do is find the crossings — [`circle_crossings`]
/// names them.
/// The work here is ordering: around the circle (θ, [`nacre_exact::quad::circular_order_about_seam`])
/// to make arcs, and along each segment (the line parameter, [`cmp_along`]) to make sub-segments.
///
/// ★ A circle nothing crosses is returned **whole**, on the road it has always taken. The parallel
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
    // ★ Nothing to cross. This used to be answered further down, after the loop had asked every
    // segment for nothing; the ends below are built for the whole slice at once, so the empty case
    // has to be answered before them or an empty class would start declining. `split_rulings`
    // opens with the same line. ☑ Measured: **61993** times over the suite — this is the common
    // case, not an edge one.
    if circles.is_empty() {
        return Ok(None);
    }
    // ★★ **Both ends of every segment, located on that segment's own line, once** — the hoist
    // `split_at_crossings` makes, for the same reason and with the same obligation: each pair is
    // located for `(wc, sg.wall)`, which is the pair [`combinatorics::closed_contains`] must be
    // called with below.
    let seg_ends: Vec<[combinatorics::OnLine; 2]> = segs
        .iter()
        .map(|s| {
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
        })
        .collect::<Result<_, BoolError>>()?;
    // Which nodes land on each circle, and which on each segment. Collected together because one
    // crossing is a point of both — splitting only one of them would leave the other's edge running
    // through a vertex it does not have.
    let canon = |x: NodeId| aliases.canon_point(x);
    let mut on_circle: Vec<Vec<NodeId>> = vec![Vec::new(); circles.len()];
    let mut on_seg: Vec<Vec<NodeId>> = vec![Vec::new(); segs.len()];
    // ★★★★★ **A contribution's own ends are cut points, and that is a rule rather than a hope.**
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
    // ★★★★★ **Where a segment's ends are, asked once per segment — and a missing one is not
    // fatal.** A pierce end has no rational coordinate at any width, and nothing below needs a
    // coordinate to decide anything; the one thing that reads it is a **filter**.
    let seg_coords: Vec<[Option<[Rat; 3]>; 2]> = segs
        .iter()
        .map(|s| s.end.map(|t| combinatorics::node_coords_rat(jd, t)))
        .collect();
    for (ci, circ) in circles.iter().enumerate() {
        let (o, m, r2) = (circ.def.origin(), circ.def.dir(), circ.def.r2());
        for (si, sg) in segs.iter().enumerate() {
            // A segment whose two ends are one point cuts nothing. Asked by **name**, which is the
            // identity: two spellings of one point are refused upstream, not tolerated here.
            // ☑ Measured unexercised.
            if sg.end[0] == sg.end[1] {
                continue;
            }
            // ★★★★★ **The cheap rejection runs where it can be formed, and skipping it asks
            // *more*, not less.** `segment_meets_cylinder` needs both endpoints' coordinates; where
            // one is a pierce point the question goes straight to `circle_crossings`, which answers
            // `Miss` exactly when the line misses. Losing the filter costs a solve, never an answer.
            if let [Some(p0), Some(p1)] = &seg_coords[si]
                && (p0 == p1 || !nacre_exact::segment_meets_cylinder(p0, p1, &o, &m, r2))
            {
                continue;
            }
            let Some(xs) = circle_crossings(jd, wc, circ, sg) else {
                return Err(undecided());
            };
            for n in xs {
                // A tangency touches without separating — `segment_meets_cylinder` above already
                // let that shape through, and cutting there would make a zero-length arc.
                //
                // ★★★★★ **That question is answered: the skip is right, and a tangency mints no
                // vertex**. The note 60 lines down reads
                // *"a circle cut at exactly one point is **slit**, not divided"* and the two
                // agree. Three reasons: measured,
                // the link at such a touch is a *single* circle and that OCCT returns the same
                // body, so minting a vertex would invent structure that is not there; the
                // `UnorderedEdges` this skip's lifting once reached is the **rulings** vocabulary
                // (`side` ±1), which the tangent wall never enters because the gate records no
                // crossing for it; and this arm is not gated on `crossings`, so ⊥ classes behave
                // exactly as they always did once the gate opens. ☑ 15 of 21 measured cells
                // assemble with it in place (the other 6 hold a third plane on the tangent line
                // and are refused earlier, by `CoincidentNodes`).
                if matches!(
                    combinatorics::pierce_name(n),
                    Some((_, _, nacre_topo::QuadRoot::Double))
                ) {
                    continue;
                }
                // ★★★★★ **Whether the crossing is on this segment is asked here, in the one
                // vocabulary that can answer it for both kinds of end.** Two plane fences would
                // need a *plane* through each endpoint, and so could only be built where every
                // end is three-plane named.
                let at =
                    combinatorics::locate(jd, cyls, wc, sg.wall, combinatorics::PointOn::Pierce(n))
                        .ok_or_else(undecided)?;
                #[cfg(test)]
                if let [Some(p0), Some(p1)] = &seg_coords[si] {
                    extent_probe::against_the_fences(
                        jd,
                        sg,
                        [p0, p1],
                        n,
                        &circ.def,
                        &at,
                        &seg_ends[si],
                        wc,
                    );
                }
                if !combinatorics::closed_contains(jd, wc, sg.wall, &at, &seg_ends[si])
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

    // ★ Nothing crossed: the caller keeps its own slices and no edge is copied. The common case
    // pays for the question and not for an answer it does not need.
    if on_circle.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    let mut out_circles = Vec::new();
    let mut arcs = Vec::new();
    let mut cut_rims: Vec<(usize, CutRim)> = Vec::new();

    // ★★★★★ **The circle's θ order runs first, because it decides which crossings are cut points
    // at all**. `circle_crossings` names where the segment's *line* meets the whole
    // circle; but the class's circle is only where its arcs are (`MergedCircle::merged` — a
    // fillet's quarter, a slot's half), and a crossing on the circle's **continuation** lies on
    // no edge of this class. Cutting it anyway puts a vertex on the segment where nothing
    // crosses it (☑ measured: 4 per op on the `tangentline slab 1.6` controls, 16 on
    // `rrect-box`, 0 on the other 353 census rows), and where that point *is* a vertex of the
    // segment already (a face through the top of a fillet's circle: the user's gusset foot,
    // `y = 0 = −2 + r`) the segment split refuses it as two names for one point.
    // **The rule: a crossing is a cut point iff it lies on an arc the class carries.** An arc's
    // own ends are cut points by construction (fed in above); a crossing is on an arc iff the
    // piece of circle before or after it is covered — being on an arc's boundary is being its
    // end, and ends are kept by name.
    //
    // ★ The old order — segments first, "the sharper question before the vaguer one" — was
    // unmeasured and said so; both halves refuse a point wearing two names with the same reject,
    // and identity on both halves is now the alias table's (`canon`), so what reaches either
    // refusal is a coincidence no discovery event recorded — the honest floor.
    struct Ordered {
        nodes: Vec<NodeId>,
        order: Vec<usize>,
        seam_is_node: bool,
    }
    let mut ordered: Vec<Option<Ordered>> = Vec::with_capacity(circles.len());
    for (ci, circ) in circles.iter().enumerate() {
        let raw = std::mem::take(&mut on_circle[ci]);
        if raw.is_empty() {
            ordered.push(None);
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
        let (order, seam_is_node) =
            circular_order(jd, circ.cyl, &circ.def, &nodes).map_err(|e| reject(e.reason()))?;
        let len = order.len();
        let mut at = vec![0usize; len];
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
        let covers = |k: usize| {
            circ.merged.iter().any(|&(_, _, arc)| match arc {
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
        };
        let keep: Vec<bool> = (0..len)
            .map(|k| is_end(nodes[order[k]]) || covers(k) || covers((k + len - 1) % len))
            .collect();
        if keep.iter().all(|&b| b) {
            ordered.push(Some(Ordered {
                nodes,
                order,
                seam_is_node,
            }));
            continue;
        }
        // The crossings on no arc leave — the circle and every segment they were put on.
        let dropped: Vec<NodeId> = (0..len)
            .filter(|&k| !keep[k])
            .map(|k| nodes[order[k]])
            .collect();
        for v in on_seg.iter_mut() {
            v.retain(|&n| !dropped.iter().any(|&d| canon(d) == canon(n)));
        }
        let kept: Vec<NodeId> = nodes
            .iter()
            .copied()
            .filter(|n| !dropped.contains(n))
            .collect();
        // ★ Nothing left to cut at: the circle is whole, and saying so here is what keeps it in
        // the class at all — an `Ordered` with no nodes would emit no arc *and* drop the circle,
        // losing an edge in silence. (No corpus row reaches it: every measured drop leaves at
        // least one arc end. It is spelled because the alternative fails quietly.)
        if kept.is_empty() {
            ordered.push(None);
            continue;
        }
        // The reduced order is the full one filtered: the relative θ order is unchanged, and it
        // still runs from the seam — which is a node only if it was one and stayed.
        let seam_is_node = seam_is_node && keep[0];
        let order: Vec<usize> = order
            .iter()
            .filter(|&&i| keep[at[i]])
            .map(|&i| kept.iter().position(|&m| m == nodes[i]).expect("kept"))
            .collect();
        ordered.push(Some(Ordered {
            nodes: kept,
            order,
            seam_is_node,
        }));
    }
    // ---- segments → sub-segments, in line order ----
    let out_segs = split_segments_at(jd, cyls, wc, segs, &mut on_seg, aliases)?;
    // ---- circles → arcs, in θ order about the seam ----
    for (ci, circ) in circles.iter().cloned().enumerate() {
        let Some(Ordered {
            nodes,
            order,
            seam_is_node,
        }) = ordered[ci].take()
        else {
            out_circles.push(circ);
            continue;
        };
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
            // ★★★★★ **Each arc takes only the contributions that cover it.** A lateral face with a
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
            // ★ **A piece no contribution covers is not an edge**. A partial rim — a
            // fillet's quarter, a slot's half — states its arc alone, and the rest of the circle
            // bounds no face on this class; emitted, it stood as a phantom edge that tied with
            // the tangent wall's line at the fillet's corner (`UnorderedEdges`, measured). The
            // straight road drops its newsless pieces (`drop_newsless`); this is the arc's twin.
            // The rim's cut nodes stay in `cut_rims` either way.
            if merged.is_empty() {
                continue;
            }
            arcs.push(MergedArc {
                cyl: circ.cyl,
                def: circ.def.clone(),
                end: [nodes[order[k]], nodes[order[(k + 1) % order.len()]]],
                merged,
            });
        }
        cut_rims.push((
            circ.cyl,
            CutRim {
                nodes: order.iter().map(|&k| nodes[k]).collect(),
                seam_is_node,
            },
        ));
    }

    Ok(Some((out_segs, out_circles, arcs, cut_rims)))
}

/// Why an order around a circle could not be formed — two causes, because two callers turn them
/// into different words (the split into a reject, the tracer into a decline), the same split
/// [`RingFail`]/`decline_of` makes on the segment side.
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
/// the seam generator, so this is the common case, not an exotic one.) Callers read `order[0]` as
/// that point when the flag says it is a node.
///
/// ★ At most one node can be seam-incident: two would be the same point, which the adjacency check
/// below refuses as the two-names-for-one-point it is.
///
/// ★★ **Two names for one point.** Two walls crossing the circle at one point are two *different*
/// pierce names — a `dedup` cannot see it, and the θ sort puts them adjacent — so an arc of zero
/// length would follow. The segment side asks this question already; asking it here too is what
/// keeps the two sides from disagreeing about what "one point" means.
///
/// Returns the permutation of `nodes` in θ order and whether the first of them is the seam point.
/// `nodes` must already be deduped by name.
pub(crate) fn circular_order(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    nodes: &[NodeId],
) -> Result<(Vec<usize>, bool), CircleOrderFail> {
    use nacre_exact::quad::{SeamOrder, circular_order_about_seam};
    let meets = nodes
        .iter()
        .map(|&n| combinatorics::pierce_meet(jd, cyl, def, n).ok_or(CircleOrderFail::Undecided))
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
    let mut bad_theta = false;
    let mut order: Vec<usize> = (0..meets.len()).collect();
    order.sort_by(|&i, &j| match cmp(i, j) {
        Some(o) => o.seam_first(),
        None => {
            bad_theta = true;
            core::cmp::Ordering::Equal
        }
    });
    if bad_theta {
        return Err(CircleOrderFail::Undecided);
    }
    let seam_is_node = seam_nodes == 1;
    for w in order.windows(2) {
        if matches!(
            cmp(w[0], w[1]),
            Some(SeamOrder::Ordered(core::cmp::Ordering::Equal))
        ) {
            return Err(CircleOrderFail::Coincident);
        }
    }
    Ok((order, seam_is_node))
}

/// **Cut every segment at the points collected on it** — the one emitter the arc split and the
/// ruling split share.
///
/// ★★★★★ **It was written twice, almost to the character.** Both roads collect crossings onto their
/// segments and then have to put them in order, refuse the ones that coincide, and hand each piece
/// the sense the whole segment ran in; the two copies drifted only in the names of their locals. One
/// consequence of joining them is that the arc split's rung is the ruling split's rung too — an end
/// a cylinder pins passes both or neither, rather than one road being taught and the other
/// forgotten, which is the shape of defect this file has already had.
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
        // ★ Identity is the table's — two crossings the table knows as one point are one.
        nodes.dedup_by(|a, b| aliases.canon_point(*a) == aliases.canon_point(*b));
        // ★★★★★ **The ruler is gone; the points are ordered by the rule, not by a parameter.**
        // This used to solve every crossing's `(line, s)`, put both endpoints on that same line with
        // `along`, and sort the parameters — which is why an end a cylinder pinned stopped it: it
        // has no rational coordinate for `along` to take. [`combinatorics::order_located`] answers
        // the same question from the pair `{wc, sg.wall}` and the points' own names, so there is no
        // ruler to lay and nothing to convert.
        //
        // ★★ **The pair is the sorted one, and that is not cosmetic.** `pierce_meet`'s canonical
        // line runs along `n_min × n_max`; asking in call order would flip every comparison on a
        // class where `wc > wall`, taking the emitted pieces and the sense with it. ☑ Differenced
        // against that ruler (`order_probe`).
        //
        // ★★ **The guard that checked all crossings share one line becomes a check on their
        // *names*.** It was there because the sort silently mixed two rulers if it stopped holding;
        // there is no ruler now, but the proposition it stood for — every crossing on this segment
        // is a point of `wc ∩ sg.wall` — still has to hold, and the name says so more cheaply than
        // two `MeetLine` comparisons did. ☑ Measured **unexercised** over the suite, like the ruler
        // check before it: it is here because the sort has no way to notice, not because a fixture
        // is red.
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
            // ★ The crossing arrives under the table's name, which may carry another
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
        // ★★★★★ **A crossing that *is* an endpoint is one point with one name, so it is deduped
        // rather than refused.** The old sentence here — "one point wearing two names, a three-plane
        // one and a pierce one" — is still true and still refused, but only for the case it
        // describes: a crossing at a *three-plane* end really does carry a second name, and the
        // equality check below catches it. Where the end was pinned by a cylinder, its name **is**
        // the crossing's `Pierce{planes, cyl, root}` — the same vertex, arrived at twice — and
        // shipping it twice would put a zero-length piece between a point and itself.
        //
        // ★ **"the same name" is asked of the alias table.** A tangent corner is a point
        // the cylinder crossing names `Pierce{[cap, wc], root}` and the plane road names
        // `ThreePlane([cap, t, wc])`; `Aliases::record_on_cylinder` has joined them, so the
        // crossing folds onto the end here as it would had the names been equal. The end's own
        // name and pin stay on the slot — they are the canonical ones (`merge_coincident`) and
        // the pieces emitted below are keyed by them. What the table does *not* know still
        // meets the equality check below, and that refusal is the honest floor.
        //
        // ★ This is why the ends' indices are remembered instead of read off the tail of the
        // vector: after a dedup the last two entries are no longer the two ends.
        // ☑ Measured over the suite: the dedup fires **400** times and the equality check below
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
        // ★★ **Nothing between them: every crossing was an end, so the segment is not cut.** It goes
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
        let mut bad = false;
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
        order.sort_by(|&i, &j| match cmp(i, j) {
            Some(o) => o,
            None => {
                bad = true;
                core::cmp::Ordering::Equal
            }
        });
        if bad {
            return Err(undecided());
        }
        #[cfg(test)]
        order_probe::against_the_ruler(jd, cyls, wc, &sg, &keyed, pair);
        // ★ A crossing that lands **on** an endpoint is one point wearing two names — a three-plane
        // one and a pierce one — and the DCEL keys vertices by name, so shipping both would make
        // two vertices where there is one. Names the alias table knows as one point were folded
        // above; what reaches here is a coincidence no discovery event recorded, and
        // refusing it is honest — there is no ground to invent the identity on.
        for w in order.windows(2) {
            if cmp(w[0], w[1]) == Some(core::cmp::Ordering::Equal) {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        }
        // ★★★ **The sub-segments' travel sense, taken from the whole segment's own two named
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

/// **Where a rational point sits along a line** — a ruler's half for a three-plane end (a pierce
/// end brings its own parameter from [`combinatorics::pierce_meet`]).
///
/// ★★ **Production lays no ruler**: it needs a *rational* point, and an end a cylinder pins has
/// none. This lives only inside [`order_probe`], which differences the order rule against it.
#[cfg(test)]
pub(super) fn along(
    line: &nacre_exact::quad::MeetLine,
    p: &[Rat; 3],
) -> Option<nacre_exact::quad::QuadVal> {
    use nacre_exact::quad::QuadVal;
    let (base, dir) = (line.base(), line.dir());
    let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    let mut rel = [Rat::from_int(0); 3];
    for k in 0..3 {
        rel[k] = p[k].checked_sub(base[k])?;
    }
    let dd = dot(&dir, &dir)?;
    let t = dot(&rel, &dir)?.checked_mul(Rat::new(dd.denom(), dd.numer())?)?;
    Some(QuadVal::from_rat(t))
}

/// The order of two line parameters. Same radicand is a subtraction and a sign; different ones
/// (two cylinders cutting one segment) land on `biquad_sign`, which answers there without
/// declining. `None` is checked-`Rat` overflow — the road's own name, not a shape answer.
pub(super) fn cmp_along(
    a: &nacre_exact::quad::QuadVal,
    b: &nacre_exact::quad::QuadVal,
) -> Option<core::cmp::Ordering> {
    use core::cmp::Ordering;
    use nacre_exact::Orient;
    let orient = if a.c() == b.c() {
        a.checked_sub(b)?.sign()
    } else {
        nacre_exact::quad::biquad_sign(
            a.a().checked_sub(b.a())?,
            a.b(),
            Rat::from_int(0).checked_sub(b.b())?,
            Rat::from_int(0),
            a.c(),
            b.c(),
        )?
    };
    Some(match orient {
        Orient::Negative => Ordering::Less,
        Orient::Positive => Ordering::Greater,
        Orient::Zero => Ordering::Equal,
    })
}

/// **Where a segment crosses a circle**, exactly — the points a `Vertex::Pierce` names.
///
/// The segment rides `wc ∩ sg.wall` and the circle is `cylinder ∩ wc`, so a crossing is
/// `plane ∩ plane ∩ cylinder` — the very shape [`nacre_exact::quad::plane_plane_cylinder`]
/// answers and
/// [`nacre_topo::Vertex::Pierce`] names. Solving along the segment instead would be shorter and
/// would yield a point with **no name**, which the next rung (splitting the circle into arcs)
/// would have to re-derive.
///
/// ★ **Only three of `CylinderMeet`'s variants can arrive here, and that is provable.** A plane
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
/// ★★ **It returns the crossings' *names*, and the coordinate is derived from the name** (see
/// [`combinatorics::pierce_point`]). The pair is solved in this function's own call order —
/// `(wc, sg.wall)` — and `NodeId::pierce` puts it in canonical order, restating the root with it.
/// Solving in ascending order instead would make the correspondence true by construction and leave
/// the canonicalization unexercised, which is where a wrong rule hides; the arc split, walking
/// segments in DCEL order, does not have that luxury either.
/// ★★★★★ **It names where the *line* crosses, and the caller says which of those are on the
/// segment.** Both roots come back: dropping the ones outside takes a fence through each
/// endpoint, and a fence is a *plane*, which an end a cylinder pins does not have. Extent is an
/// ordering question, so it belongs with the rule that answers ordering for both kinds of end
/// ([`combinatorics::closed_contains`]) — and the caller already post-filters here anyway, for the
/// tangency.
fn circle_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    circ: &MergedCircle,
    sg: &MergedSeg,
) -> Option<Vec<combinatorics::NodeId>> {
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    let w = combinatorics::class_coeffs_rat(jd, wc)?;
    let v = combinatorics::class_coeffs_rat(jd, sg.wall)?;
    let (o, m, r2) = (circ.def.origin(), circ.def.dir(), circ.def.r2());
    let roots = match nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2)? {
        CylinderMeet::Pair { s, .. } => vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])],
        // ★ `Double`, not `Lo`: the two roots coincide, so a re-sort must leave the name alone.
        CylinderMeet::Tangent { s, .. } => vec![(QuadRoot::Double, QuadVal::from_rat(s))],
        // ★★ **Reachable, and the honest answer.** The caller's `segment_meets_cylinder` filter (a
        // segment that meets the solid cylinder has a line that meets its surface) needs both
        // endpoints' coordinates, so it is skipped where one end is a pierce point, and there the
        // line genuinely can miss. "No crossings" is what a miss means.
        // ☑ Measured: **6** times over the suite.
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
