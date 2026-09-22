use super::*;
/// The DCEL walk itself: half-edges into cells, over the **three** ranges an arrangement can hold.
///
/// - segments `[0, 2·ns)` — an ordinary two-ended edge on a plane's meet with `P`;
/// - arcs `[2·ns, 2·(ns+na))` — a piece of a circle a segment cut, equally two-ended;
/// - uncut circles, appended **after** the walk on pseudo-half-edges `2·(ns+na) + 2i`, because a
///   closed curve with no vertex is not an orbit the walk can express.
pub(super) fn walk_cells(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    edges: &ClassEdges<'_>,
) -> Result<(Vec<Cell>, HashMap<usize, usize>), BoolError> {
    let (segs, arcs, circles) = (&edges.segs, &edges.arcs, &edges.circles);
    let n = segs.len();
    let he_count = edges.he_count();
    let is_arc = |he: usize| matches!(edges.kind(he), HalfEdgeKind::Arc(_));
    let origin = |he: usize| edges.origin(he);
    let edge_of = |he: usize| edges.edge_at(he);

    // Outgoing half-edges per vertex.
    let mut outgoing: HashMap<NodeId, Vec<usize>> = HashMap::new();
    for he in 0..he_count {
        outgoing.entry(origin(he)).or_default().push(he);
    }

    // For each vertex, the cyclic order of its outgoing half-edges (indices into its `outs` list).
    let mut cyclic: HashMap<NodeId, (Vec<usize>, Vec<usize>)> = HashMap::new();
    for (&v, outs) in &outgoing {
        let mut edges = Vec::with_capacity(outs.len());
        {
            watch!(E_ORDER);
            for &he in outs {
                // Direction sign away from `v` toward the far end — one spelling, `edge_dir`'s (a
                // swapped-argument `order_along(target, origin)` agrees with it only by the
                // antisymmetry of a difference cancelling an inversion).
                // ★ The maker pairs the carrier with the sense, so this site cannot pair them
                // wrong.
                // ★ The half-edge **leaves** `v`, so its direction is read at `v` — and `dir_at`
                // is where that is checked. Built once per half-edge, which is the hoist the
                // sort below rests on.
                edges.push(combinatorics::dir_at(jd, cyls, wc, &edge_of(he), v)?);
            }
        }
        let ord = timed!(E_ANGULAR, angular_order(jd, wc, &edges))?;
        cyclic.insert(v, (outs.clone(), ord));
    }

    watch!(E_WALK);
    let components = component_count(segs, arcs, &edges.rulings);

    // Try both step directions (predecessor / successor); keep the one whose bounded/outer split
    // is right. Handedness is fixed by `orient_sign(w)` and unknown up front.
    for predecessor in [true, false] {
        let mut next = vec![usize::MAX; he_count];
        for (outs, ord) in cyclic.values() {
            let len = ord.len();
            let step = if predecessor { len - 1 } else { 1 };
            for (local, &he_out) in outs.iter().enumerate() {
                let incoming = he_out ^ 1; // twin arrives at this vertex
                let pos = ord.iter().position(|&x| x == local).unwrap();
                next[incoming] = outs[ord[(pos + step) % len]];
            }
        }

        // Walk next-orbits into faces.
        let mut face_of: HashMap<usize, usize> = HashMap::new();
        let mut cells: Vec<Cell> = Vec::new();
        let mut ok = true;
        for start in 0..he_count {
            if face_of.contains_key(&start) {
                continue;
            }
            let mut cyc = Vec::new();
            let mut he = start;
            loop {
                if face_of.contains_key(&he) {
                    ok = false;
                    break;
                }
                face_of.insert(he, cells.len());
                cyc.push(he);
                he = next[he];
                if he == usize::MAX || cyc.len() > he_count {
                    ok = false;
                    break;
                }
                if he == start {
                    break;
                }
            }
            // ★★ **Three half-edges is the *straight* floor, and it is a fact about lines.** Two
            // straight edges between two points are one edge traced twice; two *arcs* between two
            // points are a lens, and an arc and a chord are a circular segment — both are honest
            // cells. So the floor reads the carriers, and a one-half-edge orbit (a circle a
            // segment only grazed, slit but not divided) is refused in either.
            let floor = if cyc.iter().any(|&h| is_arc(h)) { 2 } else { 3 };
            if !ok || cyc.len() < floor {
                ok = false;
                break;
            }
            // The cell's edges come from the walk, which knows each one's carrier and both
            // handles. The endpoint names cannot supply them: a canonical name need not mention
            // this line's planes at all, and an arc's carrier is not a plane class it could name
            // anyway.
            let ring: Vec<combinatorics::RingEdge> = cyc.iter().map(|&h| edge_of(h)).collect();
            let w = combinatorics::loop_winding(jd, cyls, wc, &ring)?;
            cells.push(Cell {
                half_edges: cyc,
                winding: w,
            });
        }
        // ★★★★ **Euler, and it is not decoration — it is the half of "this is a subdivision" the
        // contour count cannot see.** The outer-contour rule counts one number and a walk can
        // satisfy it while attaching the wrong edges: measured, dropping the canonical-to-stored
        // turn in `combinatorics::arc_side` makes this class trace **one 8-half-edge orbit where
        // there are four cells**, and the contour count passes it — a silently wrong subdivision,
        // which is the thing this file exists to refuse. `V − E + F = 2C` sees it (`0 ≠ 2`).
        //
        // `F` is the walk's cell count, which is not the topological face count: each extra
        // component contributes a second contour bounding the same outer region, so `F = f + C − 1`
        // and `V − E + f = 1 + C` becomes this. Measured across the whole suite before it was
        // asserted: **111,750 accepted walks, slack 0 in every one** — so it decides nothing that
        // was passing, and the census is bit-identical.
        let subdivides = |cells: &[Cell]| -> bool {
            let (v, e, f) = (
                outgoing.len() as i64,
                (n + arcs.len() + edges.rulings.len()) as i64,
                cells.len() as i64,
            );
            v - e + f == 2 * components as i64
        };
        if ok
            && cells.iter().filter(|c| c.winding == -1).count() == components
            && subdivides(&cells)
        {
            // ★ **Circle cells, appended after the segment orbits**. Each circle is a
            // 1-edge component the DCEL walk cannot express (orbits need ≥3): a `+1` disk cell
            // and a `−1` contour, with pseudo-half-edges numbered past the segment range —
            // `2n + 2i` (disk side) and its `^1` twin (outside), so the label propagation's
            // twin arithmetic works unmodified (`2n` is even). No crossing machinery is owed: a
            // circle and a segment cannot **cross** on this class — the gate cleared the pair or
            // recorded it — and the one shape that touches, a tangency, is skipped by
            // `split_circles`' `Double` arm rather than cut (a touch divides nothing).
            for (i, _) in circles.iter().enumerate() {
                let he_in = he_count + 2 * i;
                face_of.insert(he_in, cells.len());
                cells.push(Cell {
                    half_edges: vec![he_in],
                    winding: 1,
                });
                face_of.insert(he_in + 1, cells.len());
                cells.push(Cell {
                    half_edges: vec![he_in + 1],
                    winding: -1,
                });
            }
            return Ok((cells, face_of));
        }
    }
    Err(reject(RejectReason::RingOrientation))
}

/// Classify every cell as root / hole / plain `+1` (see [`Nesting`]). A face may carry **any number
/// of holes** (`holes[host]` is a list — e.g. a slab pierced by a U's two prongs), nested **any
/// number of levels** deep (a pocket sealed by a slab, a boss cut after fusing), and the plane may
/// carry **any number of separate bodies** (each contributes its own unbounded contour).
///
/// A root is a `-1` contour inside no `+1` ring at all, so the region it bounds lies outside every
/// cell — that is the one unbounded region, however many cycles bound it. Having none of them is
/// the only impossibility (a closed figure always has an outside), and that is `HoleRoots`. A hole
/// whose owner is not uniquely determined is `HoleDepth` (see [`innermost_host`]).
pub(super) fn nest_cells(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    cells: &[Cell],
    edges: &ClassEdges<'_>,
) -> Result<Nesting, BoolError> {
    let circles = &edges.circles;
    let n = cells.len();
    // Which circle a cell is, `None` for a cell the walk built. ★ The range test is
    // `ClassEdges::kind`'s now: an arc's half-edge is past the segments too, and reading one as a
    // circle would ask the *disk's* exact statement about a shape that is not a disk.
    let circle_of = |c: &Cell| -> Option<usize> {
        match edges.kind(*c.half_edges.first()?) {
            HalfEdgeKind::Circle(i) => Some(i),
            HalfEdgeKind::Seg(_) | HalfEdgeKind::Arc(_) | HalfEdgeKind::Ruling(_) => None,
        }
    };
    // The walk knows each edge's carrier and both handles, so the ring carries them instead of
    // leaving them to be re-derived from the endpoint names. A circle cell has no ring — its
    // containment questions run on the cylinder's exact statement instead.
    let ring_of = |c: &Cell| -> Vec<combinatorics::RingEdge> {
        if circle_of(c).is_some() {
            return Vec::new();
        }
        c.half_edges.iter().map(|&he| edges.edge_at(he)).collect()
    };
    let rings: Vec<Vec<combinatorics::RingEdge>> = cells.iter().map(ring_of).collect();
    let pos: Vec<usize> = (0..n).filter(|&i| cells[i].winding == 1).collect();

    // Union-find over cells (the pattern of `component_count`, but joining cells, not vertices).
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    let circle_ix: Vec<Option<usize>> = cells.iter().map(circle_of).collect();
    // The nesting engine asks a circle for its exact statement and nothing else, so it is
    // handed the defs rather than this pass's merged circles.
    let circle_defs: Vec<&nacre_topo::CylinderDef> = circles.iter().map(|c| &c.def).collect();
    let mut holes: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut roots: Vec<usize> = Vec::new();
    for c in (0..n).filter(|&i| cells[i].winding == -1) {
        let mut hosts: Vec<usize> = Vec::new();
        for &r in &pos {
            if crate::nesting::cell_in_cell(jd, cyls, wc, &rings, &circle_defs, &circle_ix, c, r)?
                == Some(true)
            {
                hosts.push(r);
            }
        }
        if hosts.is_empty() {
            roots.push(c);
        } else {
            // ★ **A disk may host.** Dropping the host would send the contour to `roots`, where
            // `label_cells` seeds it **void** — the boss over a bore would lose its base face and
            // the shell would open.
            let host = innermost_host(jd, cyls, wc, &rings, &circle_defs, &circle_ix, &hosts)?;
            let (rc, rr) = (find(&mut parent, c), find(&mut parent, host));
            parent[rc] = rr;
            holes.entry(host).or_default().push(c);
        }
    }
    if roots.is_empty() {
        return Err(reject(RejectReason::HoleRoots)); // no unbounded contour: not a closed arrangement
    }
    let group_of: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
    let mut root_groups: Vec<usize> = roots.iter().map(|&r| group_of[r]).collect();
    root_groups.sort_unstable();
    root_groups.dedup();
    Ok(Nesting {
        group_of,
        root_groups,
        holes,
    })
}

/// Which of `hosts` owns the hole: the **innermost** one — the candidate contained in all the
/// others.
///
/// Being inside several `+1` rings at once is ordinary two-level nesting, not a degeneracy: for
/// `A⁺ ⊃ D⁻ ⊃ B⁺ ⊃ c⁻`, `c` lies inside `B`'s ring *and* `A`'s, since `A`'s ring encloses its own
/// hole. Only `B` actually wraps `c` in material, and `emit_faces` must hang `c` off `B` — hanging
/// it off `A` would punch a hole through a face the hole is not even on.
///
/// **The innermost candidate always exists.** Cell rings are simple closed curves that do not cross
/// (a crossing would have been split into a node), so containment among them is a *total* order;
/// the candidates are a chain and its minimum is unique. Measured across the OCCT corpus:
/// every one of the 8 multi-host cases was a chain of exactly two. So the reject below is a net for
/// a broken invariant — if it ever fires, rings are crossing and the fault is upstream in
/// `split_at_crossings`, not here.
///
/// Containment is read the same way [`nest_cells`] reads it: one representative vertex through
/// [`combinatorics::point_in_ring`], and candidates sharing a node are adjacent rather than nested, so
/// they cannot be ordered and the honest answer is to reject.
fn innermost_host(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    circle_defs: &[&nacre_topo::CylinderDef],
    circle_ix: &[Option<usize>],
    hosts: &[usize],
) -> Result<usize, BoolError> {
    // ★ **Two answers, two names.** "Adjacent, so not comparable" is `None` here; "every ray was
    // spoiled" is `NoClearRay`, its cause — not `HoleDepth`, whose sentence is about nesting
    // depth and is wrong for a spoiled ray.
    //
    // ★★ **A disk can be a host, so this asks the same four-way question `nest_cells` does** —
    // [`cell_in_cell`], not a polygon-only copy of one arm (sound only while a disk is filtered
    // out before it gets here).
    let inside = |a: usize, b: usize| {
        crate::nesting::cell_in_cell(jd, cyls, wc, rings, circle_defs, circle_ix, a, b)
    };
    let mut found = None;
    for &h in hosts {
        let mut wraps_all = true;
        for &o in hosts {
            if o != h && inside(h, o)? != Some(true) {
                wraps_all = false;
                break;
            }
        }
        if wraps_all {
            if found.is_some() {
                return Err(reject(RejectReason::HoleDepth)); // two minima: not a chain
            }
            found = Some(h);
        }
    }
    found.ok_or_else(|| reject(RejectReason::HoleDepth)) // no minimum: not a chain
}

/// A per-solid, per-side material label of one cell: `[A_above, A_below, B_above, B_below]`.
///
/// ★★ **"Above" is the side of the class's *stored* plane normal** — `trace_one`'s `w_normal`,
/// which is what every `body_above` in this file is measured against. Not the class's canonical
/// rational name (that points the other way on half the classes) and not the root face's outward
/// (`frame_sign` relates the two, and `emit`'s `flip` is what folds it back in). The absence of
/// this sentence is what let a lateral face's rim contribution be written in the wrong frame.
pub(crate) type Label = [bool; 4];

/// The flip mask an edge applies when crossed, grouping its `merged` contributions **per solid**.
/// Crossing the edge XORs this into the cell label.
///
/// **The circle around the edge.** An arrangement edge is the intersection of exactly two planes —
/// `W` and the wall — so the little circle around it has only two arcs to fill, above `W` and below.
/// A solid's material fills an arc or it does not, and crossing the edge inside `W` flips a label
/// exactly for the arcs the solid's boundary separates there:
///
/// ```text
///        above (n_s(W))          [T,F] one arc  → flip the above bit
///     ────────┼────────  W       [F,T] one arc  → flip the below bit
///        below                   [T,T] both     → flip both
/// ```
///
/// This is the BRep form of a Nef local pyramid (Hachenberger & Kettner, CGAL `Nef_3`) collapsed to
/// two planes, and the XOR itself is binary winding-number propagation (Zhou et al. 2016, "Mesh
/// Arrangements for Solid Geometry"; libigl `propagate_winding_numbers`).
///
/// Precedence per solid is **Graze > Seated > Transversal**:
/// - A `Graze` is a wall touching W along this edge with its body on one side — the solid's true
///   material boundary here. It overrides a coincident `Seated` cap-rim: the two agree on a convex
///   cap (same side) but disagree at a **reflex dihedral** in W, where the graze side is correct and
///   seated-wins would flip the wrong bit.
/// - A `Seated` face is ground truth where no graze coincides: it directly says which side the body
///   is on, so it **wins over a coincident `Transversal`** (the transversal's "flip both" assumes the
///   solid straddles W, which is false exactly where a cap crosses it).
/// - A lone `Transversal` genuinely straddles → flips both sides.
///
/// `> 1 Transversal`, disagreeing graze/seated sides, or a graze coincident with a same-solid true
/// crossing is a coincident-wall degeneracy outside the corpus → honest reject.
pub(super) fn edge_mask(merged: &[(SolidSide, SegKind)]) -> Result<Label, BoolError> {
    let mut mask = [false; 4];
    for (solid, base) in [(SolidSide::A, 0usize), (SolidSide::B, 2)] {
        let kinds: Vec<SegKind> = merged
            .iter()
            .filter(|(s, _)| *s == solid)
            .map(|(_, k)| *k)
            .collect();
        let side_of_kind = |want_graze: bool| -> Vec<bool> {
            kinds
                .iter()
                .filter_map(|k| match k {
                    SegKind::Graze { body_above } if want_graze => Some(*body_above),
                    SegKind::Seated { body_above } if !want_graze => Some(*body_above),
                    // A tangent ruling flips its axis side — the seated rule with that side.
                    SegKind::Tangent { axis_above } if !want_graze => Some(*axis_above),
                    _ => None,
                })
                .collect()
        };
        let grazes = side_of_kind(true);
        let seated = side_of_kind(false);
        let transversals = kinds
            .iter()
            .filter(|k| matches!(k, SegKind::Transversal { .. }))
            .count();
        if !grazes.is_empty() {
            // Graze wins: it is the real boundary. A same-solid true crossing must not coincide.
            if transversals > 0 {
                // A graze coincident with a same-solid true crossing: the two disagree about
                // whether the solid straddles here, and nothing says which to believe.
                return Err(reject(RejectReason::EdgeOccupancyConflict));
            }
            // A valid solid brings exactly two faces to an edge, so more than two grazes on one
            // merged edge is degenerate input, not a case to arbitrate.
            if grazes.len() > 2 {
                return Err(reject(RejectReason::EdgeOccupancyConflict));
            }
            // ★ Grazes combine by **parity per side**, not by union. One graze is a step (the
            // face fills that arc — flip it). An opposite pair is a solid whose edge lies in `W`
            // with material above on one in-plane side and below on the other — both flip, the
            // `[T,T]` of the diagram, which is what a four-plane concurrency is made of. And a
            // **same-side pair flips nothing**: whether it is a knife edge (the wedge between the
            // two faces is the material, which pinches to measure zero at `W`) or its reflex
            // complement (the wedge is the void), what is immediately above `W` is the same on
            // both sides of this line. The union rule read that pair as one fill and flipped a
            // bit that changes nothing — a mask the propagation check then caught two layers
            // later as LabelConflict.
            for side in [true, false] {
                if grazes.iter().filter(|&&a| a == side).count() % 2 == 1 {
                    mask[base + usize::from(!side)] = true;
                }
            }
        } else if !seated.is_empty() {
            if seated.iter().any(|&b| b != seated[0]) {
                return Err(reject(RejectReason::EdgeOccupancyConflict)); // disagreeing seated sides
            }
            // seated wins over a coincident transversal: flip above if body_above, else below.
            mask[base + usize::from(!seated[0])] ^= true;
        } else if transversals == 1 {
            // pure transversal: the solid straddles W, flip both.
            mask[base] ^= true;
            mask[base + 1] ^= true;
        } else if transversals > 1 {
            return Err(reject(RejectReason::EdgeOccupancyConflict)); // >1 transversal, same solid
        }
        // no contributions ⇒ solid absent from this edge ⇒ no flip.
    }
    Ok(mask)
}

/// **The disk side of a cylinder's trace carries that cylinder's own solid** — `false` when it
/// does not, which is a labelling that cannot be right whatever the frame conventions are.
///
/// ★★★★★ **The absolute claim the arrangement had none of.** `label_cells` seeds the unbounded
/// contours void and flips per edge, then verifies **every** edge — and that verification is blind
/// to a global flip of one solid's two bits, because `l ^ mask == l'` survives flipping both sides
/// together. So a wrong seed inverts a whole component's material for one solid and every check
/// passes. Nothing in the crate stated an *absolute* fact about a label until this.
///
/// The fact: where the mask sets **both** of a solid's bits from a lone transversal
/// (`edge_mask`'s *"the solid straddles W, flip both"*), that solid's material is on the inside of
/// the cylinder exactly when the lateral face's outward is radially outward —
/// [`crate::planes::CylFaceInfo::orient_sign`] `= +1`, which is the `mat` the trace carries. So
/// the disk-side cell has the solid on **both** sides for a boss and on **neither** for a bore.
///
/// ★ It reads `edge_mask` rather than re-deriving the branch, so graze and seated keep their
/// priority; and it says nothing where both bits come from an **opposite graze pair** instead,
/// which carries no `mat` to compare against.
#[cfg(test)]
pub(super) fn disk_side_agrees(
    merged: &[(SolidSide, SegKind)],
    label: &Label,
) -> Result<bool, BoolError> {
    let mask = edge_mask(merged)?;
    for (solid, base) in [(SolidSide::A, 0usize), (SolidSide::B, 2)] {
        if !(mask[base] && mask[base + 1]) {
            continue;
        }
        let Some(mat) = merged.iter().find_map(|(s, k)| match k {
            SegKind::Transversal { mat } if *s == solid => Some(*mat),
            _ => None,
        }) else {
            continue;
        };
        let want = mat > 0;
        if label[base] != want || label[base + 1] != want {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Label every cell by propagating from the unbounded contours across edges, flipping per
/// `edge_mask`. The cell interior cannot be point-queried (no plane-triple name), so propagation is
/// the only route. After propagating, **every** edge's flip
/// relation is verified (`label[c] XOR mask == label[neighbour]`); a violation means the trace was
/// incomplete and is an honest reject.
///
/// ★ **`seed` is the one place global information enters the arrangement.** Everything else here
/// is a flip relation between neighbours — purely local — so a plane's labelling is determined by
/// its own segments *plus* one known classification to start from. Whole-model arrangements pass
/// `[false; 4]`, which is the statement "the region outside every contour is void" and is true
/// exactly because the arrangement covers the whole plane: `Nesting`'s unbounded region reaches
/// infinity, where neither solid is. An arrangement restricted to a region of space cannot say
/// that — its unbounded region is an artifact of the restriction — and must be told instead. That
/// is the parameter's whole reason to exist; the only caller passes `[false; 4]`.
pub(super) fn label_cells(
    cells: &[Cell],
    face_of: &HashMap<usize, usize>,
    edges: &ClassEdges<'_>,
    nesting: &Nesting,
    seed: Label,
) -> Result<Vec<Label>, BoolError> {
    // Crossing a circle flips like crossing any edge — its occupancy list is the mask's input;
    // the pseudo-half-edge numbering (`≥ 2·segs.len()`) picks the table.
    let mask_of = |he: usize| -> Result<Label, BoolError> {
        match edges.kind(he) {
            HalfEdgeKind::Seg(i) => edge_mask(&edges.segs[i].merged),
            // An arc is a piece of its circle's trace, so it carries the same contributions —
            // and a ruling piece is a piece of the lateral's trace, likewise.
            HalfEdgeKind::Arc(i) => edge_mask(&edges.arcs[i].merged),
            HalfEdgeKind::Ruling(i) => edge_mask(&edges.rulings[i].merged),
            HalfEdgeKind::Circle(i) => edge_mask(&edges.circles[i].whole_marks()?),
        }
    };
    // A face-with-holes is one region: label its group as a unit. Group representative → members.
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, &g) in nesting.group_of.iter().enumerate() {
        members.entry(g).or_default().push(i);
    }
    let mut label = vec![None; cells.len()];
    let mut queue = std::collections::VecDeque::new();
    // Seed every unbounded contour's group, enqueuing every member. The outside region may be
    // bounded by more than one cycle, and they all describe the *same* region — so one seed serves
    // them all, whatever it says.
    for g in &nesting.root_groups {
        for &i in &members[g] {
            label[i] = Some(seed);
            queue.push_back(i);
        }
    }
    while let Some(c) = queue.pop_front() {
        let lc = label[c].unwrap();
        for &he in &cells[c].half_edges {
            let nb = face_of[&(he ^ 1)];
            if label[nb].is_none() {
                let mask = mask_of(he)?;
                let lab: Label = std::array::from_fn(|i| lc[i] ^ mask[i]);
                // A hole and its host bound the same region: label the whole group at once and
                // enqueue every member, so the hole cell's edges bridge to the cell inside it —
                // the two components share no edge, so nothing else reaches the interior one.
                for &mem in &members[&nesting.group_of[nb]] {
                    if label[mem].is_none() {
                        label[mem] = Some(lab);
                        queue.push_back(mem);
                    }
                }
            }
        }
    }
    let out: Vec<Label> = label
        .into_iter()
        .collect::<Option<_>>()
        .ok_or_else(|| reject(RejectReason::UnreachedCell))?; // a cell never reached
    // Verify every edge (tree and non-tree): the flip relation must hold everywhere.
    for (he, &c) in face_of {
        let nb = face_of[&(he ^ 1)];
        let mask = mask_of(*he)?;
        if std::array::from_fn::<bool, 4, _>(|i| out[c][i] ^ mask[i]) != out[nb] {
            return Err(reject(RejectReason::LabelConflict)); // inconsistent propagation
        }
    }
    Ok(out)
}
