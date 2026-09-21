use super::*;
/// Merge every group of coplanar, same-facing result faces into one face per connected piece, then
/// dissolve straight-angle vertices — the mandatory post-boolean defeature.
///
/// **Erase the interior, do not stitch the exterior.** When two faces become one region the boundary
/// between them stops being a boundary, so the merge is: take every directed ring edge in the group,
/// drop the ones that appear as an opposed pair (`a→b` together with `b→a`), and re-thread what is
/// left. Nothing has to be spliced, which is what lets one rule cover every way the pieces can meet
/// — sharing an edge, a hole filled exactly by a neighbour, a hole filled by *several* neighbours,
/// and any chain of those (one erase settles them all at once). The earlier version stitched loops
/// with `splice_along` and so had to special-case "exactly two hole-free faces across one edge",
/// leaving `// holed — deferred` for the rest; a flush tool cap then stayed two faces forever.
///
/// A group is one plane class, one outward direction, one `flip` — mixing any of those would fold
/// material the wrong way — split further into **edge-connected components**, because faces that
/// merely lie on the same plane without touching must each survive on their own.
///
/// Reused, not reinvented: `canon` for coplanarity, `n_out.dot > 0` for facing,
/// [`combinatorics::loop_winding`] to tell an outer ring from a hole, [`combinatorics::point_in_ring`] to give
/// each hole its owner — the same two exact predicates `arrangement::nest_cells` uses for the same
/// question, neither of which reads a coordinate.
///
/// Rejects rather than guesses: a directed edge appearing twice the same way (two faces claiming the
/// same side), an undirected edge on three or more rings (non-manifold in the plane), or a node with
/// two outgoing edges after erasure (pieces meeting at a single point, where the cycle is not
/// unique).
pub(crate) fn unify_coplanar_faces(
    faces: Vec<LocalFace>,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Vec<LocalFace>, BoolError> {
    let n = faces.len();
    // One plane class, one flip.
    //
    // There used to be an outward-direction component here, and a `canon` lookup beside it. Neither
    // separated anything: `plane_idx` names a plane, so its class is itself, and the component's dot
    // product was `|n|² > 0` — identically true. It was redundant with `flip` besides, which is
    // already in the key: `assemble_fuse_cut` derives a result face's `Orientation` from exactly
    // `(the plane's frame, flip)`, so two faces in one group orient the same way by construction.
    let group_key = |lf: &LocalFace| -> (usize, bool) { (lf.surf.plane(), lf.flip) };

    // Edge-connected components within a group.
    let mut comp: Vec<usize> = (0..n).collect();
    let mut carriers: HashMap<(NodeId, NodeId), Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        // Holes count here too: a tool cap sitting flush inside another face touches it only
        // along that hole, so leaving `inner` out would put the two in different components and
        // nothing would merge at all.
        for ring in lf.poly_rings() {
            for (a, b) in ring_edges(ring) {
                carriers.entry(norm_edge(a, b)).or_default().push(fi);
            }
        }
    }
    for fs in carriers.values() {
        for w in fs.windows(2) {
            if group_key(&faces[w[0]]) == group_key(&faces[w[1]]) {
                let (ri, rj) = (uf_find(&mut comp, w[0]), uf_find(&mut comp, w[1]));
                if ri != rj {
                    comp[ri] = rj;
                }
            }
        }
    }
    // ★★ **A circle joins two faces too, and it is the only fact that says so.** A cap disk's
    // boundary is a `Bound::Circle` with **no nodes at all** (the arrangement's `bound_of`: "a
    // circle cell has no nodes — its boundary is the cylinder class itself"), so the edge table
    // above cannot see it and the disk sits alone in its own component. What connects it to the
    // face it lies in is that the other one holds *the same cylinder class* as a hole — and one
    // plane class with one cylinder class names **one circle**, so that coincidence is adjacency.
    // No position test is needed, which is why this is a table lookup and not a predicate.
    let mut circle_outer: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut circle_hole: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        if let Bound::Circle { cyl } = lf.outer {
            circle_outer.entry(cyl).or_default().push(fi);
        }
        for b in &lf.inner {
            if let Bound::Circle { cyl } = b {
                circle_hole.entry(*cyl).or_default().push(fi);
            }
        }
    }
    for (cyl, outs) in &circle_outer {
        let Some(ins) = circle_hole.get(cyl) else {
            continue; // a disk whose circle no member holds as a hole: nothing to join
        };
        for (&o, &i) in outs.iter().flat_map(|o| ins.iter().map(move |i| (o, i))) {
            if group_key(&faces[o]) == group_key(&faces[i]) {
                let (ri, rj) = (uf_find(&mut comp, o), uf_find(&mut comp, i));
                if ri != rj {
                    comp[ri] = rj;
                }
            }
        }
    }

    // Group members by component, in face order so the result is replay-stable.
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); n];
    for fi in 0..n {
        let r = uf_find(&mut comp, fi);
        members[r].push(fi);
    }

    let mut merged: Vec<LocalFace> = Vec::new();
    let mut kept: Vec<Option<LocalFace>> = faces.into_iter().map(Some).collect();
    for mem in &members {
        if mem.len() < 2 {
            continue; // nothing to merge; the face (if any) is emitted as-is below
        }
        // ★ A face whose **outer** bound is curved contributes no node edge, so the re-threading
        // below has nothing of it to thread. A **disk whose circle another member holds as a
        // hole** is an interior seam `merge_component` erases; a **disk whose circle no member
        // holds** is the merged region's own outer bound, which `merge_component` states as a
        // circle (a boss's cap that a pin's disk fills at its hole — left unmerged, the cut of a
        // stacked pin has the pin's seam vertex naming a cylinder
        // the result has no face on, `VertexNamesAbsentSurface`, measured). What still cannot be
        // threaded is a band (no producer today); spelling it keeps "a band never merges" an
        // *invariant* rather than an accident of which producer exists.
        //
        // Merging is a **correctness** matter, not the tidiness this pass once claimed: unmerged,
        // the corners on the erased boundary stay corners, keep naming a surface the result drops,
        // and the assembly refuses the whole boolean (`VertexNamesAbsentSurface`, measured on both
        // the contact fuse and the contact cut).
        if mem.iter().any(|&fi| {
            kept[fi]
                .as_ref()
                .is_some_and(|lf| matches!(lf.outer, Bound::Band { .. }))
        }) {
            continue;
        }
        let group: Vec<&LocalFace> = mem
            .iter()
            .map(|&fi| kept[fi].as_ref().expect("member present"))
            .collect();
        // Abstention: the group pinches at a point, so its faces are emitted as-is — the
        // whole-result judgement (`self_touch_reject`, run before any cell is minted) owns
        // what this shape is about to be named for.
        let Some(rings) = merge_component(&group, jd, cyls)? else {
            continue;
        };
        let (plane_idx, flip) = (group[0].surf.plane(), group[0].flip);
        merged.extend(rings.into_iter().map(|(outer, inner)| LocalFace {
            surf: crate::planes::ClassIx::Plane(plane_idx),
            outer,
            inner,
            flip,
        }));
        for &fi in mem {
            kept[fi] = None;
        }
    }
    let mut out: Vec<LocalFace> = kept.into_iter().flatten().collect();
    out.extend(merged);
    let all: Vec<usize> = (0..out.len()).collect();
    dissolve_straight_angles(&mut out, &all);
    Ok(out)
}

/// A ring's directed edges, `i → i+1` around.
fn ring_edges(ring: &[NodeId]) -> impl Iterator<Item = (NodeId, NodeId)> + '_ {
    (0..ring.len()).map(move |i| (ring[i], ring[(i + 1) % ring.len()]))
}

/// A ring's directed edges **with the carrier each rides**.
fn ring_edges_walled(ring: &Ring) -> impl Iterator<Item = ((NodeId, NodeId), Wall)> + '_ {
    (0..ring.len()).map(move |i| {
        (
            (ring.nodes[i], ring.nodes[(i + 1) % ring.len()]),
            ring.walls[i],
        )
    })
}

/// An outer bound with the bounds that belong to it — what one merged region looks like before it
/// becomes a `LocalFace`. The inner side is [`Bound`] rather than [`Ring`] because a member's
/// **circle** hole rides through the merge (a bore under a plate the merge is unifying); it has
/// no nodes, so it is carried rather than re-threaded. ★ The outer is a [`Bound`] too since cell
/// ⑩: a region whose boundary is a member's whole circle (a boss cap whose pin hole the pin's disk
/// filled) has no ring to thread either, and is stated as the circle it is.
type RegionRings = (Bound, Vec<Bound>);

/// One edge-connected group → its faces after erasing the interior boundary: each outer ring with
/// the bounds that belong to it — ring holes re-threaded, circle holes carried.
///
/// `Ok(None)` is **abstention**, and two shapes take it: the group's pieces meet at a point (the
/// merged contour would be a figure-8, so there is no cycle set to re-thread), or two of them
/// claim the same directed edge (they overlap rather than tile). Either way the caller emits the
/// group as-is.
///
/// ★ **Why abstaining is safe is not "merging is only tidiness"** — that was this pass's old
/// charter and the contact fuse refuted it: unmerged coplanar faces leave the corners on a
/// dropped wall naming a plane the result has no face on, and the assembly refuses the whole
/// boolean. What makes abstention safe is narrower: emitting the group as-is puts the pipeline
/// exactly where it stood *before this pass ran*, so nothing can come out of it that would not
/// have come out without the pass at all. The whole-result judgements ([`self_touch_reject`],
/// the closed-shell guard, the every-vertex-re-solves check) then name the shape with their own
/// sentences — and with a witness, which a capability name raised from here has none of.
/// (Rejecting `CoplanarPinch` here would be a capability-limit name for what is, on
/// every input that reaches it, a self-touching result.)
fn merge_component(
    group: &[&LocalFace],
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Option<Vec<RegionRings>>, BoolError> {
    // 1. Collect directed edges **with their walls in the key**: between one pair of pierce
    //    vertices a chord and an arc — or two complementary arcs — are *different edges*, and
    //    a node-pair key made them collide (measured: the flush half-disk pair abstained here
    //    and shipped a stated plane-self edge for the whole-result guard to refuse).
    let mut cnt: HashMap<((NodeId, NodeId), Wall), usize> = HashMap::new();
    for lf in group {
        for ring in lf.poly_rings() {
            for (e, wall) in ring_edges_walled(ring) {
                *cnt.entry((e, wall)).or_insert(0) += 1;
            }
        }
    }
    // 2a. First pass: erase exact interior pairs — the same edge walked back with the
    //     **reversed carrier** (`Wall::reversed`: the same arc back flips `ccw`, the
    //     complementary arc keeps it, which is what tells the curved twins apart; a plane
    //     wall is direction-blind). Erased 1:1 by `min` of the two multiplicities, never by
    //     fiat.
    let snapshot: Vec<((NodeId, NodeId), Wall)> = cnt.keys().copied().collect();
    let mut visited: HashSet<((NodeId, NodeId), Wall)> = HashSet::new();
    for k in snapshot {
        if !visited.insert(k) {
            continue;
        }
        let ((a, b), w) = k;
        let rk = ((b, a), w.reversed());
        visited.insert(rk);
        let c1 = cnt.get(&k).copied().unwrap_or(0);
        let c2 = cnt.get(&rk).copied().unwrap_or(0);
        let m = c1.min(c2);
        if m > 0 {
            if c1 == m {
                cnt.remove(&k);
            } else {
                cnt.insert(k, c1 - m);
            }
            if c2 == m {
                cnt.remove(&rk);
            } else {
                cnt.insert(rk, c2 - m);
            }
        }
    }
    // ★★ **Two faces claiming one side is an abstention, not a refusal** — the same shape the
    // figure-8 case takes, and for the same reason: this guard protects *the merge*, not the
    // result. Emitting the group as-is puts the pipeline back exactly where it stood before this
    // cleaning pass ran, and the whole-result judgements (`self_touch_reject`, the closed-shell
    // guard, `check_result_topology`'s plane-self edge guard) name the shape with their own
    // sentences — and with a witness. Asked **carrier-aware and after the exact erase**, so a
    // chord and an arc sharing a direction (the flush pair — the population that used to abstain
    // here into an invalid ship) are no longer a collision.
    if cnt.values().any(|&c| c > 1) {
        return Ok(None);
    }
    // (The old "carried three or more times is non-manifold" check is gone: it counted
    // *distinct directed keys* per node pair, of which there are at most two, so it could
    // never fire — a vacuous guard, measured vacuous by its own construction.)
    //
    // 2b. Second pass, wall-agnostic and **plane-only**: a pencil alias — one line under two
    //     of its planes' names — arrives opposed with two different plane walls and stays
    //     interior exactly as it always was. Only a 1:1 singleton pair may erase; an opposed
    //     *curved* leftover is not the same edge (complementary arcs) and stays a boundary;
    //     any multiplicity across an opposed pair is ambiguous and the merge abstains.
    let mut by_dir: HashMap<(NodeId, NodeId), Vec<Wall>> = HashMap::new();
    for k in cnt.keys() {
        by_dir.entry(k.0).or_default().push(k.1);
    }
    let mut erased: HashSet<(NodeId, NodeId)> = HashSet::new();
    for (&e, ws) in &by_dir {
        let re = (e.1, e.0);
        if erased.contains(&e) || erased.contains(&re) {
            continue;
        }
        let Some(rws) = by_dir.get(&re) else {
            continue;
        };
        if ws.len() == 1
            && rws.len() == 1
            && matches!(ws[0], Wall::Plane(_))
            && matches!(rws[0], Wall::Plane(_))
        {
            erased.insert(e);
            erased.insert(re);
        } else if ws.len() != 1 || rws.len() != 1 {
            return Ok(None); // multiplicity across an opposed pair: ambiguous
        }
        // An opposed 1:1 with a curved wall on either side is two different edges — both stay
        // boundaries and thread below.
    }
    // 3. Re-thread what survives. Two outgoing edges at one node means the pieces meet at a
    //    point and the cycles are not determined — the merged contour would be a figure-8:
    //    nothing to re-thread, so the merge abstains (see the function doc). The wall rides
    //    with the step, so the threaded cycle never re-derives a carrier from node names.
    let mut next: HashMap<NodeId, (NodeId, Wall)> = HashMap::new();
    for k in cnt.keys() {
        let (e, w) = (k.0, k.1);
        if erased.contains(&e) {
            continue;
        }
        if next.insert(e.0, (e.1, w)).is_some() {
            return Ok(None);
        }
    }
    let mut starts: Vec<NodeId> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut seen: HashSet<NodeId> = HashSet::new();
    let mut cycles: Vec<Ring> = Vec::new();
    for start in starts {
        if seen.contains(&start) {
            continue;
        }
        let (mut cur, w0) = next[&start];
        let mut nodes = vec![start];
        let mut walls = vec![w0];
        seen.insert(start);
        while cur != start {
            if !seen.insert(cur) {
                return Err(reject(RejectReason::CoplanarMerge)); // walk re-entered another cycle
            }
            let &(nx, w) = next
                .get(&cur)
                .ok_or_else(|| reject(RejectReason::CoplanarMerge))?;
            nodes.push(cur);
            walls.push(w);
            cur = nx;
        }
        // A 2-gon closed by a curved wall is a legal ring (`loop_winding`'s own floor); a
        // two-node polygon is not.
        let floor = if walls.iter().any(|w| !matches!(w, Wall::Plane(_))) {
            2
        } else {
            3
        };
        if nodes.len() < floor {
            return Err(reject(RejectReason::CoplanarMerge));
        }
        cycles.push(Ring::new(nodes, walls));
    }
    // 4. Winding tells an outer ring from a hole; the plane is the frame both are read in. (This
    // used to canon the index first — `plane_idx` names a plane now, so there is nothing to fold.)
    let wc = group[0].surf.plane();
    let mut outers: Vec<Ring> = Vec::new();
    let mut holes: Vec<Ring> = Vec::new();
    for cyc in cycles {
        // ★ Built from the walls the cycle carries, not from the node names. This is where a
        // four-plane concurrency used to break the merge: a canonical name need not mention the
        // plane its edge rides, and two names can share nothing but `wc`.
        let ring = cyc.edges(jd, cyls, Some(wc))?;
        match combinatorics::loop_winding(jd, cyls, wc, &ring)? {
            1 => outers.push(cyc),
            -1 => holes.push(cyc),
            _ => return Err(reject(RejectReason::CoplanarMerge)),
        }
    }
    // ★★ **A circle another member holds as its own outer bound is not a hole — it is the seam
    // between the two, and merging erases it.** That is the disk-into-face case: the face around
    // it holds the circle as a hole, the disk *is* the circle, and the merged region is simply the
    // face without it. Erasing here rather than in the caller keeps one rule in one place, and
    // costs the owner search below one question fewer.
    //
    // ★ **A circle no member holds as a hole is a region's outer bound**: the disk it
    // bounds joined the group through its *own* holes (a boss cap whose pin hole the pin's disk
    // fills) or through node edges inside it, and what is left after the seams are erased is a
    // region bounded by that circle — stated as the circle, with no ring to thread.
    //
    // Exactly one outer per circle, and at most one hole, or this abstains: two disks on one
    // circle (coincident faces) or a circle held as a hole twice is a shape this pass does not
    // arrange, and the whole-result judgements name it with their own sentences.
    let mut erased: Vec<usize> = Vec::new();
    let mut disk_outers: Vec<usize> = Vec::new();
    {
        let mut outer_of: HashMap<usize, usize> = HashMap::new();
        for lf in group {
            if let Bound::Circle { cyl } = lf.outer {
                *outer_of.entry(cyl).or_insert(0) += 1;
            }
        }
        for (&cyl, &outers) in &outer_of {
            let holes = group
                .iter()
                .flat_map(|lf| lf.inner.iter())
                .filter(|b| matches!(b, Bound::Circle { cyl: c } if *c == cyl))
                .count();
            match (outers, holes) {
                (1, 1) => erased.push(cyl),
                (1, 0) => disk_outers.push(cyl),
                _ => return Ok(None),
            }
        }
    }
    // In class order, so the result is replay-stable whatever the table's iteration order.
    disk_outers.sort_unstable();
    if outers.is_empty() && disk_outers.is_empty() {
        return Err(reject(RejectReason::CoplanarMerge)); // everything erased: not a region
    }
    // 5. Each hole belongs to the outer bound that contains it — the same question `nest_cells`
    //    asks of the arrangement's cells, answered by the same predicates.
    let mut faces: Vec<RegionRings> = outers
        .into_iter()
        .map(|o| (Bound::Ring(o), Vec::new()))
        .chain(
            disk_outers
                .into_iter()
                .map(|cyl| (Bound::Circle { cyl }, Vec::new())),
        )
        .collect();
    for hole in holes {
        // Every node is a probe, not just the first: which vertex can cast a clear ray is a fact
        // about that vertex, and settling for `nodes[0]` is what lost whole bands of rotation
        // angles here. `nesting::cell_inside` holds that retry now — one witness list walked to
        // the end — for this caller and the arrangement's alike.
        let mut owner = None;
        for (i, (outer, _)) in faces.iter().enumerate() {
            // ★ A region bounded by a **circle** asks the disk's question of the hole's nodes —
            // the one `nest_cells` asks of a ring under a disk.
            let outer = match outer {
                Bound::Ring(o) => o,
                Bound::Circle { cyl } => {
                    let hole_edges = hole.edges(jd, cyls, Some(wc))?;
                    // ★ The same engine the arrangement's nesting asks, through the same
                    // door.
                    if crate::nesting::cell_inside(
                        jd,
                        cyls,
                        wc,
                        crate::nesting::Cell::Ring(&hole_edges),
                        crate::nesting::Cell::Disk(&cyls[*cyl].def),
                    )? {
                        if owner.is_some() {
                            return Err(reject(RejectReason::CoplanarMerge)); // nested deeper than this brick names
                        }
                        owner = Some(i);
                    }
                    continue;
                }
                Bound::Band { .. } => return Ok(None),
            };
            let ring = outer.edges(jd, cyls, Some(wc))?;
            // A mixed outer takes the rational road - the (None, None) arm of the
            // arrangement's `cell_in_cell`, mirrored: each rational node of the hole is
            // a probe against the mixed walk, and the vocabulary is that arm's — a list
            // that ran out is `NoClearRay`, a list that was **empty** (a hole with no
            // rational corner) is `RingHasNoWitness`. ☑ Neither has a population today
            // (0 raises at this site); named so the first arrives under
            // the right word.
            // ★ The rational road serves a mixed outer **and a hole with no three-plane corner**
            // (a slot's stadium, four tangent corners and nothing else): its witnesses
            // are the hole's rational corners, pierce ones included
            // (`combinatorics::pierce_coords_rat`). Handing an empty probe list to the ray road
            // said `NoClearRay` for a ray never cast; the engine names an empty offer
            // `RingHasNoWitness` and an exhausted one `NoClearRay`, which are different facts.
            // ★ The same engine the arrangement's nesting asks — one rule, one set of
            // witnesses.
            let hole_edges = hole.edges(jd, cyls, Some(wc))?;
            // A hole and an outer of one coplanar group can share a node — they are pieces of one
            // boundary — so «share a node ⇒ not comparable» is not this road's rule, and the engine
            // does not carry it (it stays in the arrangement's own adapter, where it has always
            // been). Here the question is only whether the hole sits in this outer.
            let hit = crate::nesting::cell_inside(
                jd,
                cyls,
                wc,
                crate::nesting::Cell::Ring(&hole_edges),
                crate::nesting::Cell::Ring(&ring),
            )?;
            if hit {
                if owner.is_some() {
                    return Err(reject(RejectReason::CoplanarMerge)); // nested deeper than this brick names
                }
                owner = Some(i);
            }
        }
        faces[owner.ok_or_else(|| reject(RejectReason::CoplanarMerge))?]
            .1
            .push(Bound::Ring(hole));
    }
    // ★★ **The members' circle holes ride through.** A circle has no nodes, so the re-threading
    // above cannot see it — which is why this pass used to skip such a component altogether, and
    // why lifting that skip without this loop drops the hole and opens the shell (measured).
    // ★ **True of nodes, with a consequence that does not follow**: "no node" is
    // not "no boundary point that can be named" — a circle
    // has four exact ones — its rim over the cylinder's own unit frame. The kernel already mints
    // that very point as a seam vertex; what the arrangement drops is the node, not the point.
    // ★ **The engine's precondition holds here too, for a different reason**: these bounds are
    // faces of one *valid* solid lying on one plane, so no two of their loops cross — a crossing
    // would be a self-intersection the solid does not have.
    //
    // Which merged region owns a circle is a **question, not a fact about counts**: the group can
    // come apart into several outer rings, so the disk is tested against each through the engine
    // (`nesting::cell_inside` — the same door the arrangement's own nesting uses, which is also
    // where the converse a centre alone cannot give lives). None, or more than one, is refused
    // rather than guessed.
    let mut circles: Vec<usize> = group
        .iter()
        .flat_map(|lf| lf.inner.iter())
        .filter_map(|b| match b {
            Bound::Circle { cyl } => Some(*cyl),
            _ => None,
        })
        .collect();
    circles.sort_unstable();
    // One rim belongs to one face; a duplicate would make the mint count one face twice.
    circles.dedup();
    circles.retain(|c| !erased.contains(c));
    for cyl in circles {
        let def = &cyls[cyl].def;
        let mut owner = None;
        for (i, (outer, _)) in faces.iter().enumerate() {
            // A ring outer is asked the centre's parity; a circle outer the two-disk inequality
            // ([`crate::arrangement::disk_in_disk`], the nesting's own spelling).
            let inside = match outer {
                Bound::Ring(o) => {
                    let ring = o.edges(jd, cyls, Some(wc))?;
                    // ★ A disk's only witness is its **centre**, which is an interior
                    // point — it settles «the centre is in that ring», not «the disk is inside
                    // it», because a ring inside the disk contains the centre too. So
                    // the engine also asks the converse the rule calls for.
                    crate::nesting::cell_inside(
                        jd,
                        cyls,
                        wc,
                        crate::nesting::Cell::Disk(def),
                        crate::nesting::Cell::Ring(&ring),
                    )?
                }
                Bound::Circle { cyl: oc } => {
                    crate::nesting::disk_in_disk(jd, wc, def, &cyls[*oc].def)?
                }
                Bound::Band { .. } => return Ok(None),
            };
            if inside {
                if owner.is_some() {
                    return Err(reject(RejectReason::CoplanarMerge));
                }
                owner = Some(i);
            }
        }
        let owner = owner.ok_or_else(|| reject(RejectReason::CoplanarMerge))?;
        // ★ …and it must be in the region's **material**, not in one of its holes. No producer
        // makes a hole inside a hole today; asking is one call, and assuming is the shape that
        // answers about geometry that is not there.
        for hole in faces[owner].1.iter() {
            if let Some(h) = hole.ring() {
                let ring = h.edges(jd, cyls, Some(wc))?;
                if crate::nesting::cell_inside(
                    jd,
                    cyls,
                    wc,
                    crate::nesting::Cell::Disk(def),
                    crate::nesting::Cell::Ring(&ring),
                )? {
                    return Err(reject(RejectReason::CoplanarMerge));
                }
            }
        }
        faces[owner].1.push(Bound::Circle { cyl });
    }
    Ok(Some(faces))
}

/// Drop straight-angle vertices: a node whose only two neighbours across **all** rings continue
/// along the same line.
///
/// ★ **The test is the two incident edges' walls**, which the rings carry. It used to be
/// combinatorial on the names — "the node is on a line iff some pair of its planes is shared by
/// both neighbours" — which is sound only while every vertex lies on exactly three planes, the same
/// assumption that broke the merge under a four-plane concurrency. The walls say it outright.
///
/// Dropping from every incident ring at once keeps a vertex that is a real corner somewhere
/// (degree > 2), which is what stops a T-junction from opening.
///
/// ★★ **Arc-bearing rings flow through here, and the carrier keeps equality honest**.
/// While the wall was the `usize::MAX` sentinel, two *consecutive arcs* compared as "same wall"
/// and the pierce vertex between them would have dissolved — a recorded hazard, unfired only
/// because in both arc fixtures a chord or a segment sits between any two arcs. `Wall`'s derived
/// equality killed it structurally: arcs of different circles, or of one circle in different
/// directions, now compare unequal. The one pair still equal — two *same-direction* arcs of one
/// circle — is the pair for which "no turn" is geometrically true on this face, so equality
/// answers right there too; a crossing at such a vertex is kept by the other faces' rings
/// (degree > 2), the function's own rule. No population produces that consecutive pair today.
///
/// ★★ **Lateral faces participate, keyed by their own class** (the label's miss-first cell).
/// The per-solid caller hands over every face of a group — laterals included — and the wall
/// key is `(node, ClassIx)` rather than a plane index, so a curved face compares its walls
/// in its own key space. The drop proposition is the same sentence there: two ruling edges
/// on one wall are collinear exactly as two plane edges are, two same-direction arcs of one
/// circle are smooth continuation, and a cusp (the `ccw` flip) lands in `bent`. Two risks
/// are *recorded*, both with no population today: a rim circle divided only by another
/// solid's T-nodes could dissolve to fewer than two nodes (`Ring::new` does not check —
/// a comment rather than a debug_assert, because an assertion no population can fire is
/// vacuous); and a merged arc can exceed a half circle, which the edge key's CCW twin rule,
/// the seam predicate and `mass_props` all read — the crossing census's exact-volume
/// assertions are the standing measurement of that contract.
pub(super) fn dissolve_straight_angles(out: &mut [LocalFace], which: &[usize]) {
    // ★ The walls are per `(node, face plane)`. Globally they cannot be: the two result faces
    // that share a 3D edge each ride *the other's* plane as their wall, so a node in the middle of
    // that edge always sees two walls overall and would never dissolve. On each face separately it
    // sees one, which is exactly "the ring runs straight through here".
    // Per node: its neighbours, and whether any face bends there. `first_wall` remembers one wall
    // per `(node, face)` and `bent` records the first disagreement — flat maps, because a nested
    // one per node costs more than the whole pass is worth (measured: 4% of a fold).
    let mut nbrs: HashMap<NodeId, HashSet<NodeId>> = HashMap::new();
    let mut first_wall: HashMap<(NodeId, ClassIx), Wall> = HashMap::new();
    let mut bent: HashSet<NodeId> = HashSet::new();
    for &fi in which {
        let lf = &out[fi];
        for ring in lf.poly_rings() {
            let k = ring.len();
            for i in 0..k {
                let (a, b) = (ring.nodes[i], ring.nodes[(i + 1) % k]);
                nbrs.entry(a).or_default().insert(b);
                nbrs.entry(b).or_default().insert(a);
                for nd in [a, b] {
                    match first_wall.entry((nd, lf.surf)) {
                        std::collections::hash_map::Entry::Occupied(e) => {
                            if *e.get() != ring.walls[i] {
                                bent.insert(nd);
                            }
                        }
                        std::collections::hash_map::Entry::Vacant(e) => {
                            e.insert(ring.walls[i]);
                        }
                    }
                }
            }
        }
    }
    let mut drop: HashSet<NodeId> = HashSet::new();
    for (&node, ns) in &nbrs {
        // Exactly two neighbours, and on **every** face it appears in the two edges ride one wall.
        if ns.len() == 2 && !bent.contains(&node) {
            drop.insert(node);
        }
    }
    if drop.is_empty() {
        return;
    }
    for &fi in which {
        let lf = &mut out[fi];
        for ring in lf.poly_rings_mut() {
            if !ring.nodes.iter().any(|nd| drop.contains(nd)) {
                continue; // untouched rings keep their allocation
            }
            let k = ring.len();
            let mut nodes = Vec::with_capacity(k);
            let mut walls = Vec::with_capacity(k);
            for i in 0..k {
                if drop.contains(&ring.nodes[i]) {
                    continue; // its two edges are one; keep the wall of the one that survives
                }
                nodes.push(ring.nodes[i]);
                walls.push(ring.walls[i]);
            }
            *ring = Ring::new(nodes, walls);
        }
    }
}
