use super::*;
/// Face `f`'s outer-loop vertices as three-plane triples: `f`'s own plane, and the
/// neighbouring planes of the two edges meeting there.
///
/// An original vertex is as implicit a point as a seam node, so a cycle's ring — which mixes
/// them — is one uniform list and [`point_in_ring`] need not know the difference. Two
/// adjacent edges on one neighbour plane would be a straight angle, and it rejects.
#[allow(clippy::too_many_arguments)]
pub(crate) fn face_vertex_triples(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<LoopRing, BoolError> {
    loop_triples(
        model,
        &model.face(f).outer.half_edges,
        p,
        inc,
        jd,
        plane_ix,
        cyls,
    )
}

/// **A lateral face's boundary cycles**: the outer loop cut at its slit edges — the
/// self-adjacent `[lateral, lateral]` edges the assembly's outer walk climbs and descends the seam
/// on — into the pieces that walk was made of, then the pieces read back as cycles: a piece that
/// closes on itself is a rim circle or a wrapping chain; two open pieces that end where the other
/// starts are a hole the walk had spliced in; a loop with no slit is a panel. The inner loops
/// follow. Each cycle is named by [`loop_triples`] like every other loop.
///
/// The pairing is unique because a slit is never of zero length (`ZeroLengthEdge`) and at most
/// one hole is spliced (`band_loop`'s own bound), so an unpaired or odd set of pieces is a shape
/// the producer never makes — refused by the producer's own name for it.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lateral_cycles(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Vec<(CycleKind, LoopRing)>, BoolError> {
    use nacre_topo::HalfEdge;
    let face = model.face(f);
    let hes = &face.outer.half_edges;
    let n = hes.len();
    let is_slit = |he: &HalfEdge| -> Result<bool, BoolError> {
        let (_, pair) = inc
            .get(&he.edge)
            .copied()
            .ok_or_else(|| reject(RejectReason::MissingSeam))?;
        Ok(pair == [p, p])
    };
    let start_of = |piece: &[HalfEdge]| crate::he_start(model, piece[0]);
    let end_of = |piece: &[HalfEdge]| {
        let he = piece[piece.len() - 1];
        let e = model.edge(he.edge);
        if he.forward {
            e.vertices[1]
        } else {
            e.vertices[0]
        }
    };
    // Pieces between slits, in loop order from the first slit (from index 0 when there is none).
    let mut slits = 0usize;
    let mut first = 0usize;
    for (i, he) in hes.iter().enumerate() {
        if is_slit(he)? {
            if slits == 0 {
                first = i;
            }
            slits += 1;
        }
    }
    let mut pieces: Vec<Vec<HalfEdge>> = Vec::new();
    let mut cur: Vec<HalfEdge> = Vec::new();
    for k in 0..n {
        let he = hes[(first + k) % n];
        if is_slit(&he)? {
            if !cur.is_empty() {
                pieces.push(std::mem::take(&mut cur));
            }
        } else {
            cur.push(he);
        }
    }
    if !cur.is_empty() {
        pieces.push(cur);
    }
    let name = |hes: &[HalfEdge]| loop_triples(model, hes, p, inc, jd, plane_ix, cyls);
    let mut out: Vec<(CycleKind, LoopRing)> = Vec::new();
    let mut open: Vec<Vec<HalfEdge>> = Vec::new();
    for piece in pieces {
        if start_of(&piece) == end_of(&piece) {
            let ring = name(&piece)?;
            let kind = match (&ring, slits) {
                (LoopRing::Rim { .. }, _) => CycleKind::Rim,
                (_, 0) => CycleKind::Panel,
                _ => CycleKind::Chain,
            };
            out.push((kind, ring));
        } else {
            open.push(piece);
        }
    }
    while let Some(a) = open.pop() {
        let Some(j) = open
            .iter()
            .position(|b| start_of(b) == end_of(&a) && end_of(b) == start_of(&a))
        else {
            return Err(reject(RejectReason::ArcBoundNotYet));
        };
        let b = open.remove(j);
        let joined: Vec<HalfEdge> = a.into_iter().chain(b).collect();
        out.push((CycleKind::Hole, name(&joined)?));
    }
    for l in &face.inner {
        out.push((CycleKind::Hole, name(&l.half_edges)?));
    }
    Ok(out)
}

/// One loop in class form with each edge's **carried wall** beside it: `walls[i]` is the plane
/// class the edge `triples[i] → triples[i+1]` rides — the far face's class, read off `inc` where
/// the triples are produced ([`loop_triples`]), never re-derived from the endpoint names. The
/// same trust model as the merge's `Ring { nodes, walls }`: deriving a wall from two names is
/// sound only while every vertex lies on exactly three planes, and the carried value is total
/// even where the *names* degenerate (the fallback-named vertices still know their edges).
///
/// ★★★★★ **A wall is a *carrier*, not a plane index** — the same [`crate::combinatorics::Wall`] the
/// result side's `Ring` uses, and deliberately not a second vocabulary. An operand can be a
/// previous boolean's result, and then a face's ring runs along a cylinder: a boss on a wall bites
/// an arc out of the plate's caps and splits the wall with its rulings. `usize` had nowhere to
/// write that, which is why the ring was declined there rather than described.
#[derive(Clone, Debug)]
pub(crate) struct NamedRing {
    pub triples: Vec<NodeId>,
    pub walls: Vec<crate::combinatorics::Wall>,
    /// For an arc edge of a **lateral** face's ring, which way it runs about the axis — the
    /// producer's own convention (`derive_edge_curve`: a circle carrier's `[A, B]` is A to B
    /// counter-clockwise, so walking the edge `forward` is walking it CCW), read off the
    /// half-edge as `curved_wall` reads it for a plane face's arc. `None` on a plane face's ring
    /// and on a ruling. ★ Carried because the flank of an on-class run says which side the ring's
    /// *interior* is, which is the arc's direction only for a convex hole; a wrapping rim has no
    /// interior side. `cycle_on_class`'s Run arm reads this (measured equal to the flank's
    /// reading on every on-class arc of today's holes).
    pub arc_ccw: Vec<Option<bool>>,
    /// ★ Every **concurrency** this loop's corners revealed — a corner with four or more
    /// plane classes incident, as the full sorted class set. Its name in `triples` is
    /// [`canonical_triple`] of that set; the set itself goes to the arrangement's alias table
    /// (`arrangement::Aliases`) before any class is traced, so the arrangement's own discoveries
    /// fold onto the same representative.
    pub concurrencies: Vec<Vec<usize>>,
}

/// One loop of a face, in the vocabulary the tracer speaks: a polygon of three-plane
/// triples, a **full circle** — one rim edge whose far face is a cylinder, named by that
/// cylinder's class — or, on a lateral face, a **rim** — one closed rim edge whose far face is
/// a cap plane, named by that plane's class. Neither closed loop has triples, walls or
/// endpoints; forcing them through `NamedRing` was `RingNaming`'s job before the vocabulary
/// existed.
#[derive(Clone, Debug)]
pub(crate) enum LoopRing {
    Poly(NamedRing),
    Circle { cyl: usize },
    Rim { plane: usize },
}

impl LoopRing {
    /// The polygon ring, `None` for a closed circle or rim — the poly-only consumers' filter.
    pub(crate) fn poly(&self) -> Option<&NamedRing> {
        match self {
            LoopRing::Poly(nr) => Some(nr),
            LoopRing::Circle { .. } | LoopRing::Rim { .. } => None,
        }
    }
}

/// What a boundary cycle of a lateral face **is**, read off the structure of its outer loop
/// rather than any angle: the loop is the assembly's own spelling — a `lo` walk, a slit, a `hi`
/// walk, a slit, with a hole's two runs spliced in between further slits — and cutting it at the
/// slits inverts that spelling exactly (`band_loop`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CycleKind {
    /// A whole rim circle: one closed edge whose far face is a cap plane.
    Rim,
    /// A closed piece between two slits that is not a whole circle: a wrapping chain of arcs and
    /// rulings (the chain rim).
    Chain,
    /// The whole outer loop with no slit at all: a face that does not wrap the cylinder.
    Panel,
    /// An inner loop, or two open pieces between slits joined end to start — a hole the
    /// assembly spliced into the outer walk.
    Hole,
}

/// Every loop of one face in class form — **the only thing the tracer needs from `Model`**.
///
/// A loop that could not be named is `None` rather than an error, because the two failures have
/// *different names at the call site* (`OuterRing` vs `HoleRing`) and which one applies is the
/// tracer's to say, not this table's.
#[derive(Clone, Debug, Default)]
pub(crate) struct FaceLoops {
    /// The outer loop, or `None` if [`face_vertex_triples`] declined — and always `None` for a
    /// **lateral** face, whose outer loop is rims joined by self-adjacent slit edges that no
    /// triple names: `cycles` carries it cut into its rims and holes instead. Its `holes` beside
    /// it are named like any other face's.
    pub outer: Option<LoopRing>,
    /// One entry per hole ring, or `None` if [`hole_rings`] declined for **any** of them — a hole
    /// that cannot be named is not "no hole". Not filled for a **lateral** face (always `None`
    /// there): its holes are among its `cycles`, and the lateral roads read only those.
    pub holes: Option<Vec<LoopRing>>,
    /// A **lateral** face's every boundary cycle ([`lateral_cycles`]): its rims, chains and
    /// panel, with the holes after them (the spliced ones recovered) — or `None` when the outer
    /// loop could not be cut into cycles or a cycle could not be named. `None` on a plane row.
    /// The two lateral roads read these (`arrangement::lateral_shape`) and nothing else of a lateral's loops.
    pub cycles: Option<Vec<(CycleKind, LoopRing)>>,
}

/// What one boolean's tracer reads instead of the `Model`: every face's loops, plus which slots
/// belong to which operand.
///
/// ★ **Both fields are independent of the plane being traced onto.** They were nevertheless
/// re-derived once per plane class — measured on an 80-fin fold at **2,280,285** calls (one per
/// face per class) for **13,492** distinct answers, 3.3% of the boolean. Hoisting them is what
/// makes the tracer a function of a face table rather than of a topology store, and the 169×
/// reduction comes along for free.
pub(crate) struct TraceInput {
    /// Each operand's faces: the `planes`-table slot and that face's loops, in the order the
    /// shells list them.
    ///
    /// **Compact on purpose.** This was a full-length `Vec<FaceLoops>` beside a list of slots — one
    /// row per table slot whether or not that slot's face was in the input. Pairing the slot with
    /// its loops makes the length the number of faces actually traced, which is what lets a caller
    /// hand the tracer a *subset* without the table's size leaking into the cost.
    pub faces: [Vec<(usize, FaceLoops)>; 2],
    /// **The population gate's own answer, carried — never re-derived**:
    /// the `(plane class, cylinder class)` pairs the gate let through **without** proving the
    /// class's faces clear of the lateral. The tracer's ruling and chord arms fire only on pairs
    /// listed here.
    ///
    /// ★★ **It is not always empty in production**: the record-and-pass arm fills it for a
    /// wall whose plane holds the axis exactly
    /// (`planes`, `crossings.insert`). ☑ Measured: listed for the wall/boss pair 2 times in
    /// a plate-and-wall-boss fuse and 16 in the operation after it.
    pub crossings: std::collections::HashSet<(usize, usize)>,
}

/// Derive [`TraceInput`] for one boolean, once.
///
/// Walked exactly as the tracer walked: shell by shell, `surf_ix` naming each face's slot. That is
/// what keeps the table's index space the `planes` one — a face missing from `surf_ix` cannot
/// happen, since `surf_ix` was built from the same two solids.
#[allow(clippy::too_many_arguments)]
pub(crate) fn trace_input(
    model: &Model,
    operands: [(Handle<Solid>, &EdgeFaces); 2],
    surf_ix: &HashMap<Handle<Face>, usize>,
    n_faces: usize,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
    crossings: std::collections::HashSet<(usize, usize)>,
) -> TraceInput {
    let _ = n_faces;
    let mut faces = [Vec::new(), Vec::new()];
    for (side, (solid, inc)) in operands.into_iter().enumerate() {
        for sh in crate::planes::solid_shell_handles(model, solid) {
            for &fh in &model.shell(sh).faces {
                let fp = surf_ix[&fh];
                // ★ A **lateral face's outer loop** is not named as one loop: it is rims joined
                // by the chart's seam, the slit edges are self-adjacent, and no triple describes
                // their corners. `cycles` carries it cut at the slits into its rims and holes
                // and the tracer's cylinder roads read those. So `outer` stays `None`
                // here; it is a skip, not a decline.
                //
                // ★★★★★ **Its holes are named like every other face's.** A fuse can burn a hole
                // into a band (a boss straddling a plate's wall), and that loop is an ordinary
                // closed ring whose corners are `plane ∩ plane ∩ cylinder` — the shape
                // [`pierce_name_from_def`] restates. The tracer needs it to answer *per angle*
                // instead of claiming the whole circle, and naming it here is what puts the lateral
                // on the same road as everything else: one walk ([`ring_against_plane`]), not a
                // second description of the same loop.
                let loops = if matches!(plane_ix[fp], ClassIx::Cyl(_)) {
                    // A lateral's holes are among its cycles; naming them twice would ring
                    // every reject twice and read the same loop by two spellings.
                    FaceLoops {
                        outer: None,
                        holes: None,
                        cycles: lateral_cycles(model, fh, fp, inc, jd, plane_ix, cyls).ok(),
                    }
                } else {
                    FaceLoops {
                        outer: face_vertex_triples(model, fh, fp, inc, jd, plane_ix, cyls).ok(),
                        holes: hole_rings(model, fh, fp, inc, jd, plane_ix, cyls).ok(),
                        cycles: None,
                    }
                };
                faces[side].push((fp, loops));
            }
        }
    }
    TraceInput { faces, crossings }
}

/// Each hole ring of face `f`, as three-plane triples.
///
/// A rim edge's incidence is `[p, wall]`, so the neighbour plane is the wall on the other
/// side of the rim — the same construction as an outer vertex, and the same rejection of a
/// straight angle. The ring keeps its stored direction: clockwise about `f`'s outward
/// normal, which is what makes it a hole.
#[allow(clippy::too_many_arguments)]
pub(crate) fn hole_rings(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Vec<LoopRing>, BoolError> {
    model
        .face(f)
        .inner
        .iter()
        .map(|l| loop_triples(model, &l.half_edges, p, inc, jd, plane_ix, cyls))
        .collect()
}

/// A loop's vertices as names — a three-plane triple, or a pierce point where a cylinder is one of
/// the three surfaces (the face's own, or a neighbour's).
///
/// ★ **A plane vertex names itself, not the loop.** Its name is [`canonical_triple`] of the
/// plane classes of its incident faces — the set the incidence table carries ([`EdgeFaces`]) — so
/// every loop that visits the vertex hands the tracer the same name. Three classes is the ordinary
/// corner and its name is the face's own plane with the two neighbours' (no judgement asked).
/// Four or more is a concurrency: one name, and the
/// class set travels in [`NamedRing::concurrencies`] to the arrangement's alias table.
///
/// Naming from the loop instead — the face's own plane and the two neighbours the
/// meeting edges carry, with a fallback to every plane touching the vertex only when both
/// neighbours lie on one plane (a loop running straight through a shared line) — gives a
/// four-plane vertex
/// whose neighbours differ a name per face, and from the face whose two edges ride
/// the planes that share a line with its own, a triple that names no point: the judge reads that
/// «point» as lying on every class and the alias table folds everything onto a corner elsewhere.
/// The straight-through case is the one place a three-class corner is asked about
/// independence ([`Judge::plane_pair_dir_sign`]).
#[allow(clippy::too_many_arguments)]
fn loop_triples(
    model: &Model,
    hes: &[nacre_topo::HalfEdge],
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<LoopRing, BoolError> {
    let edge = |he: &nacre_topo::HalfEdge| -> Result<([Handle<Vertex>; 2], [usize; 2]), BoolError> {
        inc.get(&he.edge)
            .copied()
            .ok_or_else(|| reject(RejectReason::MissingSeam))
    };
    let other = |pair: [usize; 2]| if pair[0] == p { pair[1] } else { pair[0] };
    let n = hes.len();
    // A **full-circle loop**: one rim edge whose far face is a cylinder — no triples,
    // no walls, no endpoints; it is named by the cylinder's class. Detected structurally
    // (`[v, v]` rims are the only single-edge loops a producer makes), so no curve is read.
    if n == 1 {
        let (_, pair) = edge(&hes[0])?;
        match (plane_ix[p], plane_ix[other(pair)]) {
            (_, ClassIx::Cyl(k)) => return Ok(LoopRing::Circle { cyl: k }),
            // A lateral face's own whole rim: one closed edge whose far face is a cap plane —
            // a cycle of its outer loop, named by the plane it rides.
            (ClassIx::Cyl(_), ClassIx::Plane(w)) => return Ok(LoopRing::Rim { plane: w }),
            _ => {}
        }
    }
    let mut out = Vec::with_capacity(n);
    let mut walls = Vec::with_capacity(n);
    let mut arc_ccw = Vec::with_capacity(n);
    let mut concurrencies: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        // Vertex `i` starts edge `i` and ends edge `i - 1`.
        let (in_bounds, in_pair) = edge(&hes[(i + n - 1) % n])?;
        let (_, out_pair) = edge(&hes[i])?;
        let (a, b) = (other(in_pair), other(out_pair));
        // ★★★★★ **The filter [`crate::planes::ClassIx::plane`] asks its callers for.** That
        // accessor panics on a cylinder — "a loud panic beats a silently wrong plane" — and it is
        // right to, for the forty-odd callers whose upstream really does filter. This is the one
        // that has to *do* the filtering, because its input is an **operand**, and an operand can
        // be a previous boolean's result: a boss on a wall leaves the plate's caps bitten by an
        // arc and the wall split by two rulings, so those rings run along a cylinder. Both things
        // read below are plane-only — the carried wall and the vertex's three-plane name — so the
        // honest answer is to decline the ring rather than to name it wrongly or to abort.
        //
        // The one curved loop this road *does* speak is the full circle handled above.
        // ★★★★★ **A curved neighbour is described now, not declined.** The two questions are
        // **independent**: edge `i`'s carrier is `b`'s business, and the corner's name is both
        // neighbours' — a ruling can arrive at this corner and a plane leave it. What this road
        // lacked was the vertex's restatement into class space ([`pierce_name_from_def`] — the
        // vertex already carries its own name) and somewhere to write a curved carrier
        // ([`NamedRing`]'s `Wall`). The one curved loop answered before this point is the circle.
        //
        // ★★★★★ **And the face itself may be the cylinder.** A lateral face's *hole* is a loop
        // like any other — a rectangle of two arcs and two rulings in the chart — and its corners
        // are `plane ∩ plane ∩ cylinder`, the very shape [`pierce_name_from_def`] restates. The
        // only loop of a lateral face this road cannot walk *whole* is its **outer** one, whose
        // slit edges are self-adjacent (`other` gives back `p`) and whose corners are therefore
        // not three-surface points; `lateral_cycles` cuts it at the slits first.
        // ★★ **The corner is where half-edge `i` starts** — not "the one vertex the two edges
        // share", which is the same vertex wherever that is unique (a loop's edge `i` starts where
        // edge `i − 1` ends, validate's `OpenLoop`) and no vertex at all for a two-gon (an arc and
        // its chord share both ends). The half-edge already knows. `boolean` does not validate its
        // inputs, so the invariant this leans on is restated here, where it is leaned on.
        let corner = crate::he_start(model, hes[i]);
        debug_assert_eq!(
            if hes[(i + n - 1) % n].forward {
                in_bounds[1]
            } else {
                in_bounds[0]
            },
            corner,
            "a loop's edge starts where the previous one ends (OpenLoop)"
        );
        // ★★ **A seam joint is not a corner.** An `OnSeam` vertex lies on two surfaces only — a
        // rim's θ = 0 point, where the loop builder split a wrap arc in two — so no third
        // surface names it and the ring does not turn there. The two legs meeting at it are one
        // step of the ring: no triple, and the step's wall was pushed with its first leg. Total
        // over faces: a cap's bitten arc (legs on the cylinder) and a lateral hole's rim (legs on
        // the cap plane) read the same way.
        if matches!(*model.vertex(corner), nacre_topo::Vertex::OnSeam(_)) {
            let prev = &hes[(i + n - 1) % n];
            debug_assert!(
                matches!(model.edge_curve(prev.edge), nacre_geom::Curve::Circle(_))
                    && matches!(model.edge_curve(hes[i].edge), nacre_geom::Curve::Circle(_)),
                "a seam joint joins two arc legs"
            );
            debug_assert_eq!(a, b, "one far face on both legs of a seam joint");
            debug_assert_eq!(prev.forward, hes[i].forward, "one sense about the axis");
            continue;
        }
        let pierce = match plane_ix[p] {
            ClassIx::Plane(near) => match (plane_ix[a], plane_ix[b]) {
                (ClassIx::Cyl(k), ClassIx::Plane(far)) | (ClassIx::Plane(far), ClassIx::Cyl(k)) => {
                    Some(
                        pierce_name_from_def(model, jd, corner, k, [near, far])
                            .ok_or_else(|| reject(RejectReason::CurvedOperandBoundary))?,
                    )
                }
                // Two laterals meeting at one corner is M6b's cylinder pair, not this road's
                // (a seam joint between two legs of one arc was taken out above).
                (ClassIx::Cyl(_), ClassIx::Cyl(_)) => {
                    return Err(reject(RejectReason::CurvedOperandBoundary));
                }
                (ClassIx::Plane(_), ClassIx::Plane(_)) => None,
            },
            ClassIx::Cyl(k) => match (plane_ix[a], plane_ix[b]) {
                (ClassIx::Plane(x), ClassIx::Plane(y)) => Some(
                    pierce_name_from_def(model, jd, corner, k, [x, y])
                        .ok_or_else(|| reject(RejectReason::CurvedOperandBoundary))?,
                ),
                // A lateral's loop running along a second lateral is M6b's cylinder pair too.
                _ => return Err(reject(RejectReason::CurvedOperandBoundary)),
            },
        };
        // Edge `i`'s carried wall: the far face's class, read off `inc` — total even where the
        // vertex *names* below have to fall back or decline (see [`NamedRing`]).
        walls.push(match (plane_ix[b], plane_ix[p]) {
            (ClassIx::Plane(w), _) => crate::combinatorics::Wall::Plane(w),
            (ClassIx::Cyl(k), ClassIx::Plane(near)) => {
                let end = pierce.expect("a curved edge's corner is a pierce point");
                curved_wall(model, jd, cyls, &hes[i], k, near, end)?
            }
            (ClassIx::Cyl(_), ClassIx::Cyl(_)) => {
                unreachable!("a lateral face beside a lateral neighbour was rejected above")
            }
        });
        // A lateral face's own arc: its direction about the axis is the producer's (see
        // [`NamedRing::arc_ccw`]).
        arc_ccw.push(match (plane_ix[p], plane_ix[b]) {
            (ClassIx::Cyl(_), ClassIx::Plane(_))
                if matches!(model.edge_curve(hes[i].edge), nacre_geom::Curve::Circle(_)) =>
            {
                Some(hes[i].forward)
            }
            _ => None,
        });
        if let Some(n) = pierce {
            out.push(n);
            continue;
        }
        let (ClassIx::Plane(near), ClassIx::Plane(wall), ClassIx::Plane(far)) =
            (plane_ix[p], plane_ix[b], plane_ix[a])
        else {
            unreachable!("the curved arms are handled above")
        };
        // ★ **The vertex names itself, not the face loop.** The classes through this
        // corner are the classes of its incident faces, which the incidence table knows, and the
        // name is [`canonical_triple`] of that set. Three classes is the ordinary corner, and its
        // name is the face's own plane and both neighbours' — no judgement is asked.
        // Four or more is a concurrency: **one** name for every loop that visits the vertex, and
        // the set is handed on (`concurrencies`) so the arrangement's alias table starts from it.
        //
        // Building `[near, far, wall]` per face — falling back to the incident set only when the
        // two neighbours are one plane — gives a four-plane vertex whose neighbours differ a name
        // per face: four names, and from the face whose two edges ride the planes sharing a line
        // with its own, a triple that names no point. The judge reads that «point» as lying on
        // every class, and the alias table folds everything onto a corner elsewhere.
        let mut classes: Vec<usize> = vertex_face_indices(corner, inc)
            .into_iter()
            .filter_map(|k| match plane_ix[k] {
                ClassIx::Plane(c) => Some(c),
                ClassIx::Cyl(_) => None,
            })
            .collect();
        classes.sort_unstable();
        classes.dedup();
        let Some(t) = canonical_triple(jd, &classes) else {
            // Fewer than three classes: the vertex lies on fewer than three planes, so no triple
            // names it — a genuine straight angle (or a coplanar seam, which no producer
            // makes). Three or more with no
            // independent triple cannot happen: two distinct planes through a point meet in a
            // line, and a third off that line completes the point.
            return Err(reject(if classes.len() >= 3 {
                RejectReason::ThreePlanes
            } else {
                RejectReason::RingNaming
            }));
        };
        if classes.len() == 3 {
            // Both neighbours on one plane: the loop runs straight through and the three classes
            // may share a line — the one place this road has always asked (and still asks) the
            // judge about a three-class corner.
            let tp = t.planes();
            if far == wall && jd.plane_pair_dir_sign(tp[0], tp[1], tp[2]) == 0 {
                return Err(reject(RejectReason::ThreePlanes)); // three planes through one line, not one point
            }
            debug_assert!(
                far == wall || {
                    let mut u = [near, far, wall];
                    u.sort_unstable();
                    u == t.planes()
                },
                "a three-class corner keeps the loop's own name: {t:?} vs {near} {far} {wall}"
            );
        } else {
            concurrencies.push(classes);
        }
        out.push(NodeId::three_planes(t));
    }
    // Every joint was a seam joint: a one-edge cap rim seen from a lateral face's hole, which is
    // the cylinder pair's shape — today's answer for it, kept as a backstop.
    if out.is_empty() {
        return Err(reject(RejectReason::CurvedOperandBoundary));
    }
    Ok(LoopRing::Poly(NamedRing {
        triples: out,
        walls,
        arc_ccw,
        concurrencies,
    }))
}

/// **Canonical → outward**: the sign that carries a class's `world_rat` name to the frame
/// `orient3d` answers in, for the predicates that read the name ([`side_of`]'s pierce arm,
/// [`arc_departure_side`]).
///
/// ★★★★★ **`world_rat` is the plane's *name*, not an oriented normal.** `orient3d` answers
/// against the class's **outward** normal, and `world_rat` may be any nonzero multiple of the
/// stored one — including a negative. Both describe one plane, so they are proportional; the sign
/// of that constant is read off the first component `world_rat` makes nonzero, and `frame_sign`
/// carries stored → outward.
/// ★ Both must be nonzero, not just the rational one: they are proportional so their zero sets
/// agree *exactly*, but `raw` is `f64` and a component it rounds to zero would make `raw[i] > 0.0`
/// false and hand back a sign with nothing behind it. Requiring both turns that into a refusal.
pub(super) fn outward_fix(jd: &Judge<'_, WorkingPlane>, q: usize) -> Option<i8> {
    let co = class_coeffs_rat(jd, q)?;
    let raw = jd.planes[q].plane.coefficients();
    let zero = nacre_exact::Rat::from_int(0);
    let i = (0..4).find(|&i| co[i] != zero && raw[i] != 0.0)?;
    let k = if (co[i] > zero) == (raw[i] > 0.0) {
        1
    } else {
        -1
    };
    Some(k * jd.planes[q].frame_sign)
}

/// **Which side of `q` an arc leaves to**: the arc starts at `node` — a pierce point of
/// cylinder `cyl` on the line `q` cuts — and travels counter-clockwise about the axis when `ccw`.
/// Its tangent there is `±m̂ × (a − c)`, and the side of `q` that points to is the side the whole
/// excursion lies on (a circle meets a plane in two points).
///
/// Coordinate-free: `(a − o) · (m × n_q) = −n_q · (m̂ × (a − c))` (the axial part of `a − o` drops
/// against `m`), so the counter-clockwise tangent's side is **minus** [`ruling_side`]
/// at `a` — the very predicate that names which ruling a lateral point lies on — and a clockwise
/// arc's is plus. Then the same canonical → outward bridge as [`side_of`] ([`outward_fix`]), so
/// the answer sits in the walk's frame. `None` when the arc is tangent to `q` at `a`
/// (`ruling_side` reads zero) or a description is missing — the walk answers `Unnameable`.
///
/// ★ Not [`arc_side`], which is the *turn* of an arc against a segment in the face's own plane
/// (the winding walk's question, in the stored frame); this is a half-space of `q`.
pub(crate) fn arc_departure_side(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    node: NodeId,
    q: usize,
    cyl: usize,
    ccw: bool,
) -> Option<i8> {
    let def = &cyls.get(cyl)?.def;
    let w = class_coeffs_rat(jd, q)?;
    let (line, sv) = pierce_meet(jd, cyl, def, node)?;
    let rs = ruling_side(&w, def, (&line, &sv))?;
    Some(outward_fix(jd, q)? * if ccw { -rs } else { rs })
}

// ---------------------------------------------------------------------------
// The **rational road**: the same containment questions as above, asked about a
// point that has coordinates instead of a name.
//
// Everything else in this module names a point by three planes, which is what makes it exact.
// A cylinder's band has no such name to offer — its witness is a rational point on the axis —
// so the uniform-slab theorem needs a road that starts from coordinates and stays exact anyway.
// It does: a rational point, a rational direction, plane classes with narrow rational
// descriptions, and `point_in_ring_2d_rat` for the in-face parity. No `f64` decides anything
// here either.
// ---------------------------------------------------------------------------

/// A plane class's exact description, or `None` where it has none to give (a rotated class'
/// realized coefficients are not its truth, so they are refused rather than read).
/// **The carrier of a ring edge that rides a cylinder** — an arc or a ruling, filled or refused.
///
/// ★★ **The direction bits do not come from `edge_at`.** That convention ("even half-edge = CCW /
/// up") is the *arrangement's* edge indexing, and an operand's face loop has no such index. What is
/// true here is the model's own:
///
/// * an **arc** has a stated convention — `derive_edge_curve`'s (Plane, Cylinder) arm: "on a circle
///   carrier the vertex *order* says which arc; `[A, B]` is A to B **counter-clockwise about the
///   axis**". So walking the edge `forward` is walking it CCW.
/// * a **ruling** has none — the same arm says a plane parallel to the axis meets the lateral along
///   rulings and "the endpoints decide". `MergedRuling::end` ascending the axis is the
///   *arrangement's* convention, so `up` is **derived** here from the two endpoints' axial
///   coordinate rather than read off a rule.
///
/// `side` is [`ruling_side`]'s one spelling, and it needs a point on the ruling
/// *exactly* — which is why the pierce name comes in: [`pierce_meet`] realizes it as the `(line, s)`
/// that function takes. A ruling whose end is not a pierce point (a seam end) has no such point and
/// is refused rather than guessed.
fn curved_wall(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    he: &nacre_topo::HalfEdge,
    cyl: usize,
    near: usize,
    end: NodeId,
) -> Result<crate::combinatorics::Wall, BoolError> {
    let curved = || reject(RejectReason::CurvedOperandBoundary);
    match model.edge_curve(he.edge) {
        nacre_geom::Curve::Circle(_) => Ok(crate::combinatorics::Wall::Arc {
            cyl,
            ccw: he.forward,
        }),
        nacre_geom::Curve::Line(_) => {
            let def = cyls.get(cyl).ok_or_else(curved)?.def.clone();
            let at = pierce_meet(jd, cyl, &def, end).ok_or_else(curved)?;
            let w = class_coeffs_rat(jd, near).ok_or_else(curved)?;
            // ★ **A tangent wall has one ruling, and its side is `0`.** A fillet's
            // or a slot's own walls are tangent to their cylinder, so without this every such
            // operand falls here as `CurvedOperandBoundary`.
            // ★ The reading is the **point**'s, not the corner's **root**'s
            // ([`ruling_side_signed`], which answers the axis plane instead
            // of abstaining as `ruling_side` does for the ray caster's sake): the root says `0`
            // only when `near` is the very wall the name pairs, and `near` may be a plane through
            // the axis carrying that same corner on one of its two rulings.
            let side = ruling_side_signed(&w, &def, (&at.0, &at.1)).ok_or_else(curved)?;
            // Which way travel runs along the axis: the stored edge ascends when its second
            // endpoint does, and `forward` says whether this half-edge walks it that way.
            //
            // ★★★ **There is nothing to read, so it is derived — from the two ends' definitions,
            // never from their realizations.** A ruling edge's stored pair carries no order:
            // `edge_for` keys it `unordered(va, vb)` ("a ruling edge is straight, so the unordered
            // pair orders it"), unlike an arc, whose `[A, B]` *is* the CCW convention. So the two
            // ends have to be compared — and each end is a pierce point of `near`, the cylinder,
            // and **one other plane**. `near` *holds* the ruling, so it is that other plane that
            // **cuts** it, and where it crosses the axis is a rational question
            // ([`crate::planes::axis_param_of_plane`]).
            //
            // ☑ **The two parameters cannot tie**: a plane that *meets* this cylinder's faces is
            // parallel to the axis or perpendicular to it — the gate admits an oblique class only
            // after proving it misses every lateral face, so no pierce vertex names one
            // — and a parallel plane cannot cut a ruling. So both cutting planes are caps, and
            // distinct caps cross the axis at distinct parameters. The strict `>` therefore
            // restates the comparison it replaces exactly, rather than growing a decline for a
            // case that has none.
            let other_param = |v: Handle<Vertex>| -> Option<nacre_exact::Rat> {
                let nacre_topo::Vertex::Pierce { planes, .. } = *model.vertex(v) else {
                    return None;
                };
                let mut cut = None;
                for &h in &planes {
                    let c = *model.world_plane_name(h)?.narrow()?;
                    // `near`'s own class, in whichever of the two spellings this vertex carries.
                    if plane_sense(&c, &w).is_some() {
                        continue;
                    }
                    if cut.replace(c).is_some() {
                        return None;
                    }
                }
                crate::planes::axis_param_of_plane(&cut?, &def)
            };
            let [v0, v1] = model.edge(he.edge).vertices;
            let (t0, t1) = (
                other_param(v0).ok_or_else(curved)?,
                other_param(v1).ok_or_else(curved)?,
            );
            let ascends = t1 > t0;
            Ok(crate::combinatorics::Wall::Ruling {
                cyl,
                side,
                up: ascends == he.forward,
            })
        }
    }
}
