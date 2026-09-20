use super::*;
/// Rebuild the result solids from the arrangement's faces. Pushes into the arena; it does not take
/// the operands at all, which is the point — the live set is its caller's to move.
#[allow(clippy::too_many_arguments)]
pub(crate) fn reconstruct(
    model: &mut Model,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::draft::CutRims,
    deferred: Option<BoolError>,
    tangencies: Tangencies<'_>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes = jd.planes;
    // No faces means no result — `Common` of two solids that miss each other, `Cut` of a box that
    // is wholly inside what cuts it. That is an answer, not a failure: a solid is bounded by faces,
    // so a non-empty result cannot have none. The inputs are still consumed, exactly as they are on
    // any other successful boolean — this returns `Ok`, so the caller's retire runs.
    if faces.is_empty() {
        return Ok(Vec::new());
    }
    let named = name_result_vertices(jd, seam, faces, cyls, cut_rims);
    // The naming's failure yields to the deferred stopper like every stage before it; the raise
    // itself now stands at the assembly's very end below.
    let Named {
        grouping,
        group_of,
        per_solid,
        defs: def_triple,
    } = match named {
        Ok(n) => n,
        Err(e) => return Err(deferred.unwrap_or(e)),
    };
    let faces: &[LocalFace] = per_solid.as_deref().unwrap_or(faces);
    // ★ **The tangency verdict stands beside the self-touch sieve, and for the same reason** —
    // both are whole-result judgements that need the grouping and must speak *before* a handle is
    // minted, so a refusal leaves the arena as it found it. A held grouping error stays held: the
    // question "do these two lumps share a solid" is only askable of a grouping that answered.
    if let Ok(g) = &grouping {
        tangency_reject(&tangencies, faces, g)?;
    }

    // Vertices (deterministic: first appearance across faces in order).
    let mut vh: HashMap<(usize, NodeId), Handle<Vertex>> = HashMap::new();
    let mut node_handle =
        |model: &mut Model, g: usize, node: NodeId| -> Result<Handle<Vertex>, BoolError> {
            if let Some(&h) = vh.get(&(g, node)) {
                return Ok(h);
            }
            // A face references a seam node that was not welded into `seam` — a reconstruction
            // dropped a crossing. Reject (never panic): an unmodeled flush topology must decline
            // honestly, not abort the kernel (DNA).
            //
            // ★ This used to `match` the node to compare its triple against `SeamVertex.triple`.
            // Both are identities now, so the comparison is the identity's own `==` and the
            // destructure that existed only to reach the payload is gone.
            let handle = {
                let sv = seam
                    .iter()
                    .find(|s| s.triple == node)
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                // A vertex that is a corner of no face at all has no name in the result's own
                // planes — a degeneracy, and the honest answer is the one this reason already
                // carries ("a corner with no turn").
                let tri = match def_triple.get(&(g, node)).copied() {
                    Some(Def::Three(t)) => t.planes(),
                    // ★★★ **A pierce vertex is minted from its declaration** — the def already
                    // names two result plane classes and the cylinder, so this is the class→handle
                    // mapping and nothing else. That mapping is where `QuadRoot::canonical`
                    // answers a **second** time: `NodeId::Pierce` is canonical in *class* order,
                    // `Vertex::Pierce` in *handle* order, and the class→handle map is not
                    // monotone in general — a re-sort must carry the root through
                    // (`transform`'s remap already locks the same rule on the way back out).
                    // ★ The flip is unexercised **at this call**: dropping it leaves every fence
                    // green, because both fixtures' class order happens to match their handle
                    // order. ★★★★★ **The rule is not unexercised, though — its inverse
                    // is red.** The population that exercises it is a
                    // boolean's *result used as the next operand*, where a second boolean builds
                    // its classes afresh and their order is not the handles': dropping the same
                    // restatement in `combinatorics::pierce_name_from_def` turns
                    // `an_operand_bounded_by_a_cylinder_is_named_in_class_space` red. So what is
                    // still owed here is only a fixture that reaches *this* line, not the rule.
                    // The tolerance is the Three arm's rule below, one surface swapped: measured
                    // against the very `model.surface` objects `validate` reads, maxed with
                    // `sv.tol` (which `pierce_vertex_tol` built, meet-line term included).
                    Some(Def::Pierce {
                        planes: p2,
                        cyl,
                        root,
                    }) => {
                        let (pair, root) = nacre_topo::QuadRoot::canonical(
                            [planes[p2[0]].surf, planes[p2[1]].surf],
                            root,
                        );
                        let cylinder = cyls[cyl].surf;
                        let def = Vertex::Pierce {
                            planes: pair,
                            cylinder,
                            root,
                        };
                        let h = crate::realize::push_vertex_realized(
                            model,
                            def,
                            PointCache::Unrealized { coord: sv.point },
                            crate::realize::ChainLink::Fresh,
                        );
                        vh.insert((g, node), h);
                        return Ok(h);
                    }
                    None => return Err(reject(RejectReason::StraightAngle)),
                };
                let def = Vertex::ThreePlane([
                    planes[tri[0]].surf,
                    planes[tri[1]].surf,
                    planes[tri[2]].surf,
                ]);
                // ★★ **The arrangement's figure is a fallback coordinate now, and nothing more.**
                //
                // This site used to measure a residual here — the distance from `sv.point` to the
                // planes the vertex is *re-named* in, `max`ed with `sv.tol` — and store it in the
                // cache. It measured the right thing (the re-naming swaps the triple wherever four
                // planes concur, and the old code's carried tolerance was short on 454 of 91,394
                // corpus vertices), but a cache that stores a residual is storing the wrong *kind*
                // of knowledge: a residual is one distance to the carriers and says nothing about
                // how far the coordinate is from the truth. The realization answers that, and it
                // runs first (`push_vertex_realized`). What the arrangement still owes the kernel
                // is the seam table's `tol`, which the self-touch sieve reads directly.
                crate::realize::push_vertex_realized(
                    model,
                    def,
                    PointCache::Unrealized { coord: sv.point },
                    crate::realize::ChainLink::Fresh,
                )
            };
            vh.insert((g, node), handle);
            Ok(handle)
        };
    // Materialize all vertex handles first. Outer then inner, rings in order: `vh`'s
    // first-appearance order fixes the vertex handles, and replay depends on it.
    let mut materialized = Ok(());
    'mat: for (fi, lf) in faces.iter().enumerate() {
        for &node in lf.poly_rings().flat_map(|r| r.iter()) {
            if let Err(e) = node_handle(model, group_of[fi], node) {
                materialized = Err(e);
                break 'mat;
            }
        }
    }
    // The materialization's failure yields to the deferred stopper like the naming's above; the
    // raise itself now stands at the assembly's very end below.
    if let Err(e) = materialized {
        return Err(deferred.unwrap_or(e));
    }

    // ★★ **The rim table**: a curved bound has no node and no triple, so it is not
    // welded through the seam table at all — it is minted here, once per `(group, cylinder class,
    // plane class)`, and every face that meets that rim asks for the same pair back. That sharing
    // is what makes the rim edge's use count two (the cap face once, the band once) rather than
    // two separate edges the closed-shell guard would call dangling.
    //
    // The group is in the key for the reason it is in `vh`'s: two result solids must not share a
    // handle.
    let mut rim: RimTable = HashMap::new();
    let rims_built = (|| -> Result<(), BoolError> {
        for (fi, lf) in faces.iter().enumerate() {
            let g = group_of[fi];
            let mut keys: Vec<(usize, usize)> = std::iter::once(&lf.outer)
                .chain(lf.inner.iter())
                .flat_map(|b| match (b, lf.surf) {
                    (Bound::Circle { cyl }, ClassIx::Plane(c)) => vec![(*cyl, c)],
                    (Bound::Band { lo, hi }, ClassIx::Cyl(k)) => [lo, hi]
                        .into_iter()
                        .filter_map(Rim::circle)
                        .map(|c| (k, c))
                        .collect(),
                    _ => Vec::new(),
                })
                .collect();
            // ★★ **And every rim some ring's *wrap arc* rides.** The loop builder splits such an
            // arc at that rim's seam vertex, so the vertex has to exist — and until now it existed
            // only because a **band** happened to claim the same rim. That is a coupling, not a
            // rule: a cap's arc is split for the cap's own reason, and the merged lateral face
            // this ladder is heading for keeps only its outermost rims, so nothing would register
            // the cut ones at all (measured: the notch's wrap arc dies `MissingSeam`).
            //
            // ★ **Only a wrap arc**, never every arc: a rim nothing splits would gain an orphan
            // seam vertex. And **appended**, never prepended, so a key a band also names is still
            // minted in the old order — the values are identical either way, but the vertex
            // handles are not, and a renumbering is exactly what the census cannot see.
            //
            // The ambiguous-naming decline is swallowed here on purpose: this table *offers*
            // keys, it does not judge. `ring` runs the same derivation and raises the same name
            // when it reaches that face, which keeps the honest reject where its witness is.
            for b in std::iter::once(&lf.outer).chain(lf.inner.iter()) {
                // ★ `rings()`, so a band's chain rims offer their wrap arcs too — stated for
                // totality rather than need: a chain's arc is where the lateral meets a plane
                // face, so that face's ring offers the same rim (☑ dropping this arm changed
                // nothing on the chain population; the control is vacuous by construction).
                for r in b.rings() {
                    let n = r.nodes.len();
                    for t in 0..n {
                        if let Ok(Some(key)) = wrapping_rim(
                            lf.surf,
                            r.nodes[t],
                            r.nodes[(t + 1) % n],
                            r.walls[t],
                            cut_rims,
                        ) {
                            keys.push(key);
                        }
                    }
                }
            }
            for (k, c) in keys {
                if rim.contains_key(&(g, k, c)) {
                    continue;
                }
                let cut = cut_rims.get(&(k, c));
                let (lat, plane) = (cyls[k].surf, planes[c].surf);
                // The seam point of this rim, spelled as `add_cylinder` spells one: the axis meets the
                // plane at the circle's centre, and `θ = 0` is the `+ref_dir` side of it.
                let realized = cyls[k].realized;
                let centre = nacre_geom::intersect::line_plane(
                    &realized.axis(),
                    match model.surface_cache(plane) {
                        nacre_geom::Surface::Plane(p) => p,
                        // ★ **Spelled out rather than `_`, so a third surface kind is a compile
                        // error here.**
                        //
                        // The arm is unreachable today: `plane` comes from `jd.planes`, the
                        // *plane* class table (`WorkingPlane` holds a `geom::Plane`), while a
                        // cylinder's rows live in `cyls` (`WorkingCyl`) — and `c` indexes the
                        // cutting plane of a `cut_rims` key. The reject is kept rather than
                        // `unreachable!` because that is a **table** invariant, not a type one:
                        // `planes`'s `unreachable!("a cylinder truth carries a cylinder
                        // cache")` is the type-guaranteed shape, and this is not that.
                        //
                        // ⚠ The reason's wording ("three planes that should meet do not") does
                        // not describe a non-planar cap; it is harmless only because nothing
                        // reaches it. An author who makes it reachable owes it a real name.
                        nacre_geom::Surface::Cylinder(_) => {
                            return Err(reject(RejectReason::ThreePlanes));
                        }
                    },
                )
                .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                let point = centre + realized.ref_dir() * realized.radius();
                // ★★ **A cut circle mints no closed rim edge — but its seam vertex stands.**
                // The `[v, v]` spelling is false for a cut circle (measured: both arc fixtures
                // minted one that only the reject then discarded), while θ = 0 is still where the
                // band's joint must sit (a `[lat, lat]` seam edge derives a line from its
                // endpoints, so both must share one θ — and the seam is model geometry, fixed at
                // `+ref_dir`). Degenerate case first: when a pierce vertex lies **on** the seam
                // generator — the split's own `SeamIncident` classification, carried in
                // `CutRim::seam_is_node`, never re-derived from coordinates — that vertex *is*
                // the seam point and minting another would stand a second handle on the same
                // point (a zero-length arc piece nothing downstream could see).
                let v = match cut {
                    Some(cr) if cr.seam_is_node => vh[&(g, cr.nodes[0])],
                    _ => crate::realize::push_vertex_realized(
                        model,
                        Vertex::OnSeam([lat, plane]),
                        PointCache::Unrealized { coord: point },
                        crate::realize::ChainLink::Fresh,
                    ),
                };
                let e = match cut {
                    Some(_) => None,
                    None => Some(
                        model
                            .push_edge(Edge::carrier_pair(lat, plane), [v, v])
                            .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?,
                    ),
                };
                rim.insert((g, k, c), (v, e));
            }
        }
        Ok(())
    })();
    // A rim's failure yields to the deferred stopper like every stage before it.
    if let Err(e) = rims_built {
        return Err(deferred.unwrap_or(e));
    }

    // ★ **An edge's carriers are the two faces that use it — read off the whole result, not
    // guessed from one side.** The obvious per-face answer ("my plane, plus the wall my
    // arrangement says the edge rides") is wrong exactly where the resolved four-plane
    // concurrency lives: three planes through one LINE make each face's wall a *third* plane
    // that legitimately contains the edge but does not bound it here — measured, the two faces
    // sharing such an edge named two different walls. The faces themselves are the ground
    // truth, and every ring is already in hand, so one pre-scan reads it.
    let mut pair_surfs: HashMap<(usize, usize), Vec<Handle<Surface>>> = HashMap::new();
    // ★ **The same reading for a ruling edge's plane carrier.** Deriving it
    // from the two ends' *names* (the plane they share) states the line only while a
    // point's name is the pair that minted it: once the alias table represents a tangent corner
    // by another pair — a class through the axis crossing the cylinder there — the shared plane
    // is that class, which contains the line but does not bound the edge. The face that bounds
    // it does, and the scan already walks every face. Keyed like the edge itself (`EdgeKey`),
    // because a ruling and a plane-pair line can share one vertex pair.
    /// A ruling edge's key in that scan — `EdgeKey::Ruling`'s fields, which is the point.
    type RulingKey = (usize, i8, (usize, usize));
    let mut ruling_surfs: HashMap<RulingKey, Vec<Handle<Surface>>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        // A band contributes no node edge — its rims and seam are minted with their carriers
        // stated (`[lateral, plane]`, `[lateral, lateral]`), so there is nothing for this scan
        // to read off it.
        let ClassIx::Plane(fc) = lf.surf else {
            continue;
        };
        let fsurf = planes[fc].surf;
        for r in lf.poly_rings() {
            let k = r.nodes.len();
            for t in 0..k {
                // ★ Only plane-carried edges feed the scan. An arc edge shares its vertex pair
                // with the chord between the same pierce vertices, and pushing this face's plane
                // here would pollute the chord's line-key entry into the fallback arm — the arc
                // states its carriers directly instead (`edge_for`'s arc arm).
                let (va, vb) = (vh[&(g, r.nodes[t])], vh[&(g, r.nodes[(t + 1) % k])]);
                let pair = unordered(va.index() as usize, vb.index() as usize);
                match r.walls[t] {
                    Wall::Plane(_) => pair_surfs.entry(pair).or_default().push(fsurf),
                    Wall::Ruling { cyl, side, .. } => ruling_surfs
                        .entry((cyl, side, pair))
                        .or_default()
                        .push(fsurf),
                    // An arc edge shares its vertex pair with the chord between the same pierce
                    // vertices, and pushing this face's plane here would pollute the chord's
                    // line-key entry into the fallback arm — the arc states its carriers directly
                    // instead (`edge_for`'s arc arm).
                    Wall::Arc { .. } => continue,
                }
            }
        }
    }

    // Edges keyed by `EdgeKey` (lookup only).
    let mut edge_of: HashMap<EdgeKey, Handle<Edge>> = HashMap::new();
    // A ring's two consecutive nodes are distinct arrangement vertices, so their points differ and
    // the line through them exists. Reject rather than panic if it does not: an aborting kernel is
    // below the floor (`overview.md`: out-of-coverage input declines honestly). `SEAM_ALIAS`
    // catches the known way this happens — two triples on one point — at the seam table, where the
    // names are still in hand, so this is a backstop with no firing test (cf. `NON_MANIFOLD_EDGE`).
    let mut edge_for = |model: &mut Model,
                        va: Handle<Vertex>,
                        vb: Handle<Vertex>,
                        // The two ends' names — `Some` wherever the caller walks a ring of
                        // named nodes; the band chain passes `None` (its seam vertex has no
                        // name), and only the ruling arm requires them.
                        ends: Option<(NodeId, NodeId)>,
                        wall: Wall,
                        face_surf: Handle<Surface>|
     -> Result<Handle<Edge>, BoolError> {
        match wall {
            Wall::Plane(w) => {
                let pair = unordered(va.index() as usize, vb.index() as usize);
                let key = EdgeKey::Line(pair);
                if let Some(&e) = edge_of.get(&key) {
                    return Ok(e);
                }
                // Manifold edges (everything a green result contains) have exactly two uses. Any
                // other count is on its way to the existing non-manifold reject — the fallback
                // (this face's wall + plane) keeps construction deterministic until that reject
                // fires, deciding nothing new.
                let surfaces = match pair_surfs[&pair][..] {
                    [a, b] => Edge::carrier_pair(a, b),
                    _ => Edge::carrier_pair(planes[w].surf, face_surf),
                };
                let e = model
                    .push_edge(surfaces, [va, vb])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
                edge_of.insert(key, e);
                Ok(e)
            }
            Wall::Ruling { cyl, side, .. } => {
                // ★ A ruling edge is straight, so the unordered pair orders it (no complementary
                // pieces — the arc's problem does not arise); `(cyl, side)` keys it apart from
                // plane edges. The carriers are **the edge's own fact**, stated from its end
                // names (the shared wall plane) and the cylinder — never from `face_surf`, whose
                // value depends on which face minted first (the wall face or the panel), and a
                // carrier that depends on mint order is exactly the kind of drift
                // `EdgeCarrierMismatch` exists to catch.
                let pair = unordered(va.index() as usize, vb.index() as usize);
                let key = EdgeKey::Ruling { cyl, side, pair };
                if let Some(&e) = edge_of.get(&key) {
                    return Ok(e);
                }
                // The plane that **bounds** it, read off the face that does (the scan above);
                // the ends' shared plane is the fallback for an edge no plane face carries.
                let bounding = match ruling_surfs.get(&(cyl, side, pair)).map(|v| {
                    let mut v = v.clone();
                    v.sort_unstable();
                    v.dedup();
                    v
                }) {
                    Some(v) if v.len() == 1 => Some(v[0]),
                    _ => None,
                };
                let wall_surf = match bounding {
                    Some(s) => s,
                    None => {
                        planes[ends
                            .and_then(|(a, b)| shared_pierce_plane(a, b))
                            .ok_or_else(|| {
                                // Two Pierce ends that share no single wall plane: a naming this
                                // ladder does not arrange yet.
                                reject(RejectReason::RulingBoundNotYet)
                            })?]
                        .surf
                    }
                };
                let e = model
                    .push_edge(Edge::carrier_pair(cyls[cyl].surf, wall_surf), [va, vb])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
                edge_of.insert(key, e);
                Ok(e)
            }
            Wall::Arc { cyl, ccw } => {
                // ★★★ **An arc edge is minted in CCW order** — `[A, B]` is the piece from A to
                // B counter-clockwise about the axis, so the two complementary arcs between one
                // pierce pair are `[A, B]` and `[B, A]`: the vertex order is the last bit the
                // endpoints alone cannot give (`EdgeKey`'s note). The carriers are stated
                // directly — the two faces using an arc edge are this cap and the cylinder's
                // side, so there is nothing for the scan to read — which also keeps the chord's
                // line key clean (`pair_surfs` skips arc edges for the same reason).
                let (from, to) = if ccw { (va, vb) } else { (vb, va) };
                let key = EdgeKey::Arc {
                    cyl,
                    from: from.index() as usize,
                    to: to.index() as usize,
                };
                if let Some(&e) = edge_of.get(&key) {
                    return Ok(e);
                }
                let e = model
                    .push_edge(Edge::carrier_pair(cyls[cyl].surf, face_surf), [from, to])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
                edge_of.insert(key, e);
                Ok(e)
            }
        }
    };

    // ★★ **A cut rim's boundary chain, minted once per `(group, cylinder, plane)` — before the
    // face loop, in face order.** The chain walks the circle CCW from the seam vertex through the
    // pierce nodes and back (the wrap arc's two seam-split pieces included), so a cut rim's
    // pieces exist under their `EdgeKey`s before either consumer asks: the cap faces' ring steps
    // weld to these very handles by key, and `band_loop` reads the chain directly — which is
    // also why this is not minted inside `band_loop`: two closures cannot both own `edge_for`.
    // Iterated over the faces' curved bounds (not the map) so the minting order is
    // deterministic, the same discipline as the rim table above.
    let mut band_chains: HashMap<(usize, usize, usize), Vec<HalfEdge>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let keys: Vec<(usize, usize)> = std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .flat_map(|b| match (b, lf.surf) {
                (Bound::Circle { cyl }, ClassIx::Plane(c)) => vec![(*cyl, c)],
                (Bound::Band { lo, hi }, ClassIx::Cyl(k)) => [lo, hi]
                    .into_iter()
                    .filter_map(Rim::circle)
                    .map(|c| (k, c))
                    .collect(),
                _ => Vec::new(),
            })
            .collect();
        for (k, c) in keys {
            if band_chains.contains_key(&(g, k, c)) {
                continue;
            }
            let Some(cr) = cut_rims.get(&(k, c)) else {
                continue; // an uncut rim's boundary is its closed edge, no chain to build
            };
            // The rim table ran this very enumeration, so the entry exists whenever this loop
            // reaches the key — the arm is a backstop, not a population (cf. the ring closure's
            // `MissingSeam`, which *is* reachable: a cap ring's arc steps do not put the key in
            // either enumeration, only a curved bound does).
            let Some(&(sv, _)) = rim.get(&(g, k, c)) else {
                continue;
            };
            let mut vs: Vec<Handle<Vertex>> = Vec::with_capacity(cr.nodes.len() + 1);
            if !cr.seam_is_node {
                vs.push(sv);
            }
            vs.extend(cr.nodes.iter().map(|n| vh[&(g, *n)]));
            debug_assert_eq!(vs[0], sv, "the chain starts at the seam vertex");
            let m = vs.len();
            let mut chain = Vec::with_capacity(m);
            for i in 0..m {
                let (u, v) = (vs[i], vs[(i + 1) % m]);
                let e = edge_for(
                    model,
                    u,
                    v,
                    None,
                    Wall::Arc { cyl: k, ccw: true },
                    planes[c].surf,
                )?;
                let forward = model.edge(e).vertices[0] == u;
                chain.push(HalfEdge { edge: e, forward });
            }
            band_chains.insert((g, k, c), chain);
        }
    }

    let mut face_handles = Vec::new();
    let mut assembled: Result<(), BoolError> = Ok(());
    'faces: for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let face_surf = match lf.surf {
            ClassIx::Plane(c) => planes[c].surf,
            ClassIx::Cyl(k) => cyls[k].surf,
        };
        let mut ring = |model: &mut Model, r: &Ring| -> Result<Loop, BoolError> {
            let handles: Vec<Handle<Vertex>> = r.nodes.iter().map(|nd| vh[&(g, *nd)]).collect();
            let k = handles.len();
            let mut half_edges: Vec<HalfEdge> = Vec::with_capacity(k);
            for t in 0..k {
                let (va, vb) = (handles[t], handles[(t + 1) % k]);
                // ★★ **The wrap arc is minted as two pieces, split at the seam vertex.** θ = 0
                // lies inside exactly one arc of a cut circle (unless a pierce vertex sits on the
                // seam — `CutRim::seam_is_node`, in which case nothing splits), and the band's
                // seam edge must end there; splitting on the *cap* side is what hands the band
                // the same two edges and keeps the shell guard's use count at two.
                //
                // ★ The wrap test is **directed**: with two pierce nodes the two complementary
                // arcs share one unordered endpoint pair, so the match orients the ring step by
                // the `ccw` bit the wall carries and compares against the split's own θ order
                // (`nodes.last() → nodes[0]` is the piece that wraps past θ = 0).
                let split_at = match wrapping_rim(
                    lf.surf,
                    r.nodes[t],
                    r.nodes[(t + 1) % k],
                    r.walls[t],
                    cut_rims,
                )? {
                    Some(key) => {
                        let &(v, _) = rim
                            .get(&(g, key.0, key.1))
                            .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                        Some(v)
                    }
                    None => None,
                };
                let legs = match split_at {
                    Some(s) => vec![[va, s], [s, vb]],
                    None => vec![[va, vb]],
                };
                for [u, v] in legs {
                    let e = edge_for(
                        model,
                        u,
                        v,
                        Some((r.nodes[t], r.nodes[(t + 1) % k])),
                        r.walls[t],
                        face_surf,
                    )?;
                    // For an arc edge the stored order is CCW, so this reads back exactly the
                    // `ccw` bit the wall carried in.
                    let forward = model.edge(e).vertices[0] == u;
                    half_edges.push(HalfEdge { edge: e, forward });
                }
            }
            Ok(Loop { half_edges })
        };
        // ★★ **Which way a rim circle is walked.** `derive_edge_curve` builds it as
        // `Circle::from_center_normal(centre, axis, ref_dir, r)`, so its parameter runs **CCW
        // about the axis direction `d`**. A loop must run CCW about its own face's outward
        // normal, so a bound on a face whose normal agrees with `d` is walked forward and one
        // whose normal opposes it backward — exactly the convention `add_cylinder` writes down
        // (top cap `forward: true`, bottom cap `false`). A *hole* runs the other way again,
        // because an inner loop keeps the material on its left by winding against the outer.
        //
        // The face's outward normal is the stored surface normal when the face is `Forward` and
        // its negation when `Reversed`, which is the `flip` decision made below — so the sign is
        // read off `frame_sign` and `flip` here rather than carried in from anywhere.
        // The **unflipped** face normal against the axis: `flip` is applied once, to every kind of
        // bound, right below — reading it here too would toggle the winding twice.
        //
        // ★★★★★ **This `f64` stays, and `world_rat` cannot replace it.** Every other axial
        // decision on this road was moved onto exact descriptions, because they all ask *where* a
        // plane crosses the axis — a question about the plane, whose answer does not depend on
        // which way its coefficients are written. This one asks something else: which way this
        // class's **frame** points relative to the axis. A class has no outward normal at all
        // ([`crate::planes::WorkingPlane`]'s own words — it holds faces from both operands and two
        // of them can oppose); what it has is the frame, and `plane.normal()` with `frame_sign`
        // *is* that frame. `world_rat` is the plane's **name**, canonicalised, and carries no
        // frame direction: measured over the suite, `world_rat · m` came out positive on all 436
        // calls while this sign varied, so substituting it would have inverted the winding on 203
        // of them. There is no precision to gain either — a plane bounding a circle is
        // perpendicular to the axis, so `|n · m|` is maximal, as far from a close call as the
        // quantity gets.
        let axis_sign = |k: usize| -> f64 {
            let c = lf.surf.plane();
            let s = if planes[c].frame_sign > 0 { 1.0 } else { -1.0 };
            planes[c]
                .plane
                .normal()
                .dot(cyls[k].realized.axis().direction())
                * s
        };
        let circle_loop =
            |model: &mut Model, cyl: usize, cls: usize, hole: bool| -> Result<Loop, BoolError> {
                // ★★ **The cut check comes before the rim lookup, and answers with the
                // population's own name.** A cut circle cannot bound a whole disk — the trace
                // subdivides that disk into cells — so reaching here with one is a producer
                // inconsistency this backstop names honestly rather than as a dropped crossing
                // (`MissingSeam` would misdiagnose a `SuspectedDefect`).
                if cut_rims.contains_key(&(cyl, cls)) {
                    return Err(reject(RejectReason::ArcBoundNotYet));
                }
                let (_, e) = *rim
                    .get(&(g, cyl, cls))
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                // An uncut circle's rim always carries its closed edge; the cut arm was refused
                // one line above, so `None` here is the same dropped-crossing shape.
                let e = e.ok_or_else(|| reject(RejectReason::MissingSeam))?;
                let _ = model;
                let mut forward = axis_sign(cyl) > 0.0;
                if hole {
                    forward = !forward;
                }
                Ok(Loop {
                    half_edges: vec![HalfEdge { edge: e, forward }],
                })
            };
        // **A periodic face's outer boundary: its loops, joined along the seam generator.**
        //
        // The four-half-edge spelling `add_cylinder` builds — two rims and the seam used twice in
        // opposite senses — is the case where **nothing is in the seam's way**. When a hole meets
        // the seam, that same generator is *interrupted* by it: the slit runs from the `lo` rim up
        // to the hole, the hole's own boundary carries the walk across, and a second slit finishes
        // the climb to the `hi` rim; coming back down uses each piece the other way. One rule,
        // whose no-hole case is exactly today's four half-edges.
        //
        // ★★★★★ **Which of the hole's two runs the climb takes is forced, not chosen.** Arriving
        // at the lower contact the walk must leave along an edge that keeps material on its left,
        // and the hole is *already* wound that way (it is this face's hole) — so the climb takes
        // the run leaving that contact **in the hole's own direction**, and the descent takes the
        // other. No side test, no sign. It is also what keeps the `seam_is_node` shape honest:
        // there the seam segment *is* the hole's own ruling, and it lands in the climb as a single
        // hole edge used **once**, rather than being minted again as a slit and used a third time.
        //
        // The hole is **consumed** — it stops being an inner loop, because it is now part of the
        // outer walk. A hole that does not meet the seam is left alone.
        // ★★ **Where the seam generator meets a cylinder's cut circles** — the only points a
        // hole or a chain rim can touch it at. One rule, two spellings, both carried rather
        // than re-derived: a pierce vertex sitting **on** the seam *is* the contact
        // (`CutRim::seam_is_node` — the split's own classification), and otherwise it is the
        // `OnSeam` vertex the rim table minted for that circle.
        let contacts_of = |k: usize| -> Vec<(Handle<Vertex>, usize)> {
            cut_rims
                .iter()
                .filter(|((kk, _), _)| *kk == k)
                .filter_map(|(&(_, c), cr)| {
                    let v = if cr.seam_is_node {
                        cr.nodes.first().and_then(|n| vh.get(&(g, *n)).copied())
                    } else {
                        rim.get(&(g, k, c)).map(|&(v, _)| v)
                    };
                    v.map(|v| (v, c))
                })
                .collect()
        };
        // ★★ A contact's **station**: the axial parameter of the cut circle it sits on. A cut
        // circle's plane is perpendicular to the axis, so equal stations would be the same
        // plane, hence one class — two contacts on distinct circles never tie. A rational
        // question, asked of the class the contact table carries beside the vertex.
        // [`crate::planes::WorkingCyl::cache`] states the rule this follows: the f64 twin is
        // what a measurement reads, and every decision reads `def`.
        let station_of = |k: usize,
                          contacts: &[(Handle<Vertex>, usize)],
                          v: Handle<Vertex>|
         -> Option<nacre_exact::Rat> {
            let c = contacts.iter().find(|&&(x, _)| x == v).map(|&(_, c)| c)?;
            crate::planes::axis_param_of_plane(
                &crate::combinatorics::class_coeffs_rat(jd, c)?,
                &cyls[k].def,
            )
        };
        // ★★ **A whole rim's traversal: the closed edge, or the pre-minted chain.** Both walk
        // the circle CCW from the seam vertex, which is where the slit attaches.
        let rim_walk = |k: usize, c: usize| -> Result<(Vec<HalfEdge>, Handle<Vertex>), BoolError> {
            let (v, closed) = *rim
                .get(&(g, k, c))
                .ok_or_else(|| reject(RejectReason::MissingSeam))?;
            let hes = match closed {
                Some(e) => vec![HalfEdge {
                    edge: e,
                    forward: true,
                }],
                None => band_chains
                    .get(&(g, k, c))
                    .cloned()
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?,
            };
            Ok((hes, v))
        };
        // **A periodic face's outer boundary from its two rim walks**, each already in walk
        // order and starting at the vertex the slit attaches to: `lo` walked forward, `hi`
        // backward — a whole circle's closed edge or chain, or a wrapping chain's ring
        // rotated to its contact (`build` makes both).
        let band_loop = |model: &mut Model,
                         k: usize,
                         (lo_hes, v_lo): (Vec<HalfEdge>, Handle<Vertex>),
                         (hi_hes, v_hi): (Vec<HalfEdge>, Handle<Vertex>),
                         holes: &mut Vec<Loop>|
         -> Result<Loop, BoolError> {
            let lat = cyls[k].surf;
            let contacts = contacts_of(k);
            let is_contact = |v: Handle<Vertex>| contacts.iter().any(|&(x, _)| x == v);
            let taken = holes.iter().position(|lp| {
                lp.half_edges
                    .iter()
                    .any(|&he| is_contact(model.he_start(he)))
            });
            let cut = match taken {
                None => None,
                Some(idx) => {
                    let lp = holes.remove(idx);
                    let hits: Vec<usize> = (0..lp.half_edges.len())
                        .filter(|&i| is_contact(model.he_start(lp.half_edges[i])))
                        .collect();
                    // Two contacts, or this is a shape the walk does not arrange. The merge
                    // that produces such a hole abstains on the same count, so a face reaching
                    // here with any other number is a producer inconsistency, not an input.
                    let [i, j] = hits[..] else {
                        return Err(reject(RejectReason::ArcBoundNotYet));
                    };
                    let (pi, pj) = (
                        model.he_start(lp.half_edges[i]),
                        model.he_start(lp.half_edges[j]),
                    );
                    let run_i = lp.half_edges[i..j].to_vec();
                    let mut run_j = lp.half_edges[j..].to_vec();
                    run_j.extend_from_slice(&lp.half_edges[..i]);
                    // ★★ Which contact the **lower** slit reaches: the one nearer the `lo` rim
                    // along the axis. A contact sits on a cut circle, and a cut circle's station
                    // *is* the axial parameter of the plane that cut it — a rational question,
                    // asked of the class the contact table is already carrying beside it.
                    // [`crate::planes::WorkingCyl::cache`] states the rule this follows: the f64
                    // twin is what a measurement reads, and every decision reads `def`.
                    //
                    // ☑ **The two stations cannot tie**, so the `<=` restates the comparison it
                    // replaces rather than growing a case: `cut_rims` is keyed by `(k, c)`, so a
                    // circle offers at most one contact and two contacts are two distinct
                    // circles — and a cut circle's plane is perpendicular to the axis, so equal
                    // stations would be the same plane, hence one class and one `c`.
                    let station = |v: Handle<Vertex>| station_of(k, &contacts, v);
                    let (si, sj) = match (station(pi), station(pj)) {
                        (Some(a), Some(b)) => (a, b),
                        // The gate already required a world description of every plane class
                        // before a band could be emitted, so this is spelled rather than assumed
                        // away — the same discipline the both-rims-cut arm above follows.
                        _ => return Err(reject(RejectReason::ArcBoundNotYet)),
                    };
                    Some(if si <= sj {
                        (run_i, run_j, pi, pj)
                    } else {
                        (run_j, run_i, pj, pi)
                    })
                }
            };
            let slit = |model: &mut Model, a, b| {
                model
                    .push_edge(Edge::carrier_pair(lat, lat), [a, b])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))
            };
            let (seam_edge, seam_hi, up, down) = match cut {
                None => {
                    let e = slit(model, v_lo, v_hi)?;
                    (e, e, Vec::new(), Vec::new())
                }
                Some((up, down, lower, upper)) => {
                    let a = slit(model, v_lo, lower)?;
                    let b = slit(model, upper, v_hi)?;
                    (a, b, up, down)
                }
            };
            let bridged = seam_edge != seam_hi;
            let mut half_edges = lo_hes;
            half_edges.push(HalfEdge {
                edge: seam_edge,
                forward: true,
            });
            half_edges.extend(up);
            if bridged {
                half_edges.push(HalfEdge {
                    edge: seam_hi,
                    forward: true,
                });
            }
            half_edges.extend(hi_hes);
            if bridged {
                half_edges.push(HalfEdge {
                    edge: seam_hi,
                    forward: false,
                });
            }
            half_edges.extend(down);
            half_edges.push(HalfEdge {
                edge: seam_edge,
                forward: false,
            });
            Ok(Loop { half_edges })
        };
        // ★★★★★ **The inner loops are built first, and `flip` is applied to none of them yet.**
        // A band whose hole meets the seam splices that hole's *own half-edges* into its outer
        // walk, so it has to be handed them in the sense every bound is written in — the
        // unflipped one. Reversing per bound as they were built (the old shape) would splice a
        // reversed run into an unreversed walk, and the result is a closed loop that is quietly
        // wound wrong: watertight, right triangle count, wrong solid.
        let mut build = |model: &mut Model,
                         holes: &mut Vec<Loop>,
                         b: &Bound,
                         hole: bool|
         -> Result<Loop, BoolError> {
            match b {
                Bound::Ring(r) => ring(model, r),
                Bound::Circle { cyl } => circle_loop(model, *cyl, lf.surf.plane(), hole),
                Bound::Band { lo, hi } => {
                    let k = lf
                        .surf
                        .cyl()
                        .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                    // ★ A cut rim never arrives as `Rim::Circle`: the emitter spells it
                    // once, as the chain of its arcs (`cyl_chart::regions`), so a whole-circle
                    // rim the split cut is a producer inconsistency — named, not walked.
                    for rim in [lo, hi] {
                        if let Some(c) = rim.circle()
                            && cut_rims.contains_key(&(k, c))
                        {
                            return Err(reject(RejectReason::ArcBoundNotYet));
                        }
                    }
                    let contacts = contacts_of(k);
                    // ★★ **A rim's walk, and where the slit attaches to it.** A whole circle:
                    // the closed edge or the pre-minted chain, `lo` forward and `hi` backward
                    // (reverse the pieces and flip each sense), from its seam vertex. A
                    // wrapping chain: its own ring as the loop builder spells it (wrap arcs
                    // split at the seam vertex), **rotated to start at its contact** and walked
                    // as stored — the cleaning pass kept it in walk order. With several
                    // contacts the slit takes the one nearest the other rim along the axis:
                    // a `lo` chain's **highest** station, a `hi` chain's **lowest** — the
                    // stretch of θ = 0 between them is then interior to the face, as a
                    // spliced hole's contacts are ordered. ☑ Both senses measured (a boss
                    // whose cap sits inside the other body, and its mirror).
                    let mut walk_of = |model: &mut Model,
                                       rim: &Rim,
                                       forward: bool|
                     -> Result<(Vec<HalfEdge>, Handle<Vertex>), BoolError> {
                        match rim {
                            Rim::Circle(c) => {
                                let (mut hes, v) = rim_walk(k, *c)?;
                                if !forward {
                                    hes.reverse();
                                    for he in &mut hes {
                                        he.forward = !he.forward;
                                    }
                                }
                                Ok((hes, v))
                            }
                            Rim::Chain(r) => {
                                let mut hes = ring(model, r)?.half_edges;
                                let mut best: Option<(usize, nacre_exact::Rat)> = None;
                                for (i, &he) in hes.iter().enumerate() {
                                    let v = model.he_start(he);
                                    if !contacts.iter().any(|&(x, _)| x == v) {
                                        continue;
                                    }
                                    let s = station_of(k, &contacts, v)
                                        .ok_or_else(|| reject(RejectReason::ArcBoundNotYet))?;
                                    let better = match &best {
                                        None => true,
                                        Some((_, b)) => {
                                            if forward {
                                                s > *b
                                            } else {
                                                s < *b
                                            }
                                        }
                                    };
                                    if better {
                                        best = Some((i, s));
                                    }
                                }
                                // A chain winds once, so it passes the seam somewhere: a chain
                                // with no contact is a producer inconsistency, named.
                                let (i, _) =
                                    best.ok_or_else(|| reject(RejectReason::ArcBoundNotYet))?;
                                hes.rotate_left(i);
                                let v = model.he_start(hes[0]);
                                Ok((hes, v))
                            }
                        }
                    };
                    let lo_walk = walk_of(model, lo, true)?;
                    let hi_walk = walk_of(model, hi, false)?;
                    band_loop(model, k, lo_walk, hi_walk, holes)
                }
            }
        };
        let mut none: Vec<Loop> = Vec::new();
        let mut inner: Vec<Loop> = match lf
            .inner
            .iter()
            .map(|b| build(model, &mut none, b, true))
            .collect::<Result<Vec<_>, BoolError>>()
        {
            Ok(v) => v,
            Err(e) => {
                assembled = Err(e);
                break 'faces;
            }
        };
        let mut outer = match build(model, &mut inner, &lf.outer, false) {
            Ok(lp) => lp,
            Err(e) => {
                assembled = Err(e);
                break 'faces;
            }
        };
        if lf.flip {
            // Cut's inside-A B-pieces, and a bore's wall: reverse every loop and toggle the
            // orientation below, so the outward normal points into the removed region and
            // each loop still keeps material on its left. ★ One place for every bound kind —
            // a curved loop reversed only by its own builder was how the first drilled box
            // came out with its rim used twice in the same sense (`NonOpposedEdge`).
            for lp in std::iter::once(&mut outer).chain(inner.iter_mut()) {
                lp.half_edges.reverse();
                for he in &mut lp.half_edges {
                    he.forward = !he.forward;
                }
            }
        }
        debug_assert!(
            none.is_empty(),
            "only the outer bound bridges, and only a band does"
        );
        // The plane's frame *is* the root face's orientation: `frame_sign` carries that face's
        // `Forward`/`Reversed` as a sign (read off the stored flag since the cutover). Reading it
        // here is what used to be `planes[plane_idx].orient` — a face field indexed by a plane,
        // the shape of every bug this split exists to prevent.
        // A cylinder's stored normal points away from its axis, and that *is* the face's outward
        // normal for a boss — so its "framed" sense is `Forward`, and `flip` (material outside
        // the wall: a hole) is what reverses it. Planes read their class's frame instead.
        let framed = match lf.surf {
            ClassIx::Cyl(_) => Orientation::Forward,
            ClassIx::Plane(c) if planes[c].frame_sign > 0 => Orientation::Forward,
            ClassIx::Plane(_) => Orientation::Reversed,
        };
        let orientation = if lf.flip {
            match framed {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
            }
        } else {
            framed
        };
        face_handles.push(model.push_face(Face {
            surface: face_surf,
            outer,
            inner,
            orientation,
        }));
    }
    // The face loop's failure yields to the deferred stopper like every stage before it; the
    // raise itself now stands at the assembly's very end below (an incomplete face set has no
    // shell to count, so an assembly failure still stops here).
    if let Err(e) = assembled {
        return Err(deferred.unwrap_or(e));
    }
    // Closed-shell guard: a 2-manifold b-rep uses every edge exactly twice (once from each of the
    // two faces that share it). A reconstruction that emits a face set with a dangling edge (use
    // count 1) or a pinched one (>2) is not a solid — `validate` would call it `NonManifoldEdge`,
    // but `boolean` never runs `validate` on its own output, so without this the caller receives a
    // silently invalid solid. Honest-reject instead ("honest-reject > silent-wrong").
    // Counted on the welded `Handle<Edge>`s, so it is exact and coordinate-free.
    {
        let mut uses: HashMap<Handle<Edge>, usize> = HashMap::new();
        for &fh in &face_handles {
            let f = model.face(fh);
            for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &l.half_edges {
                    *uses.entry(he.edge).or_default() += 1;
                }
            }
        }
        // A use count above two is a *pinch*: the two bodies meet exactly along that edge, and no
        // 2-manifold solid contains it (the edge twin of `NonManifoldVertex`). A count of one is a
        // *dangling* edge — the assembly dropped a face, which is ours to fix, not the input's.
        // The witness is the minimum offending edge handle, not the map's first hit: `HashMap`
        // iteration order would make the reported location differ between runs.
        if let Some(eh) = uses
            .iter()
            .filter(|&(_, &n)| n > 2)
            .map(|(&e, _)| e)
            .min_by_key(|e| e.index())
        {
            let [va, vb] = model.edge(eh).vertices;
            return Err(deferred.unwrap_or(crate::reject_at(
                RejectReason::NonManifoldResultEdge,
                crate::RejectWhere::Segment([model.vertex_point(va), model.vertex_point(vb)]),
            )));
        }
        if uses.values().any(|&n| n != 2) {
            return Err(deferred.unwrap_or(reject(RejectReason::OpenResultShell)));
        }
    }
    // ★ **The grouping decided at the top of this function, raised here** — where the old code
    // decided it, so a boolean that declines leaves the arena cells it always left. What it says:
    // one component is the whole result; several mean either an enclosed void (a cavity, an
    // inward-oriented shell) or a severed operand (two or more material-enclosing shells), and a
    // surviving cavity belongs to the piece whose outer shell nests it. See [`group_faces`] — and
    // for why the *group*, not the component, is the unit the handles above were minted per.
    // The grouping's failure yields to the deferred stopper like every stage before it; the
    // raise itself now stands at the very end, past the solid assembly below.
    let Grouping {
        labels,
        n,
        positives,
        mut comps_of,
        ..
    } = match grouping {
        Ok(g) => g,
        Err(e) => return Err(deferred.unwrap_or(e)),
    };
    let mut by_comp: Vec<Vec<Handle<Face>>> = vec![Vec::new(); n];
    for (i, &fh) in face_handles.iter().enumerate() {
        by_comp[labels[i]].push(fh);
    }
    let shells: Vec<Handle<Shell>> = by_comp
        .iter()
        .map(|faces| {
            model.push_shell(Shell {
                faces: faces.clone(),
            })
        })
        .collect();
    // Emit each material solid (with its cavities) in a canonical, replay-stable order keyed on
    // geometry, so a downstream op can index the returned Vec deterministically. ★ One material
    // needs no key and must not pay for one — `comp_key` walks every face of the piece.
    let order: Vec<usize> = if positives.len() == 1 {
        vec![0]
    } else {
        let keys: Vec<Vec<[f64; 3]>> = positives
            .iter()
            .map(|&c| comp_key(model, &by_comp[c]))
            .collect();
        let mut order: Vec<usize> = (0..positives.len()).collect();
        order.sort_by(|&x, &y| {
            keys[x]
                .partial_cmp(&keys[y])
                .expect("finite vertex coordinates")
        });
        order
    };
    let out: Result<Vec<Handle<Solid>>, BoolError> = Ok(order
        .into_iter()
        .map(|oi| {
            let c = positives[oi];
            let comps = comps_of
                .remove(&c)
                .expect("every material component owns a component list");
            // A cavity shell's faces already point into the void (the material is outside it, so
            // the material-on-correct-side reconstruction winds them inward) — measured, no flip.
            let cavities = comps
                .iter()
                .copied()
                .filter(|&x| x != c)
                .map(|x| shells[x])
                .collect();
            model.push_solid(Solid {
                outer: shells[c],
                cavities,
            })
        })
        .collect());
    // ★★★ **The deferred stopper's raise — the very end of the assembly.** The socket is empty
    // today, but the interception ladder it crowns is architecture: seven
    // layers (per-class → seam stretch → naming → vertex materialization → face loop → shell
    // guard → grouping and the solid assembly) all yield to `deferred`, so the next
    // out-of-coverage class that plugs a stopper into `arrange`'s socket carries one name out
    // however far this pipeline gets.
    //
    // ★★ **An arc reject therefore leaves a complete garbage solid in the store —
    // deliberately.** Vertices, edges, faces, shells and the solid itself: cells outside the
    // live set (`assemble_fuse_cut` retires operands only on `Ok`, and nothing pushed here is
    // reachable from a live solid), the same class of residue a late reject's arena cells have
    // always been. The live-set fences stay green, and a session that keeps recording after a
    // reject rebuilds from the log (`replay`'s discipline, stated in `docs/design.md`). Holding
    // the raise any earlier would put the solid assembly behind an interception nothing can see
    // past — the unreachable-machinery trap this ladder keeps refusing.
    if let Some(d) = deferred {
        return Err(d);
    }
    // ★★ **Every result vertex must re-solve from the result's own faces** — the property
    // transform and replay stand on. It is a shipped check, not a debug_assert, because the
    // chained contact-cut reaches it: a bored plate cut by a boss that only touches its top kept
    // the
    // top ring's pierce vertices, defs still naming the boss's cylinder with every boss face
    // gone — right volume, wrong names, and release builds shipped it silently. Refusing here
    // is the floor until the assembly learns to shed the stale corners; the garbage-solid
    // residue is the same class every late reject leaves (see the raise above).
    if let Ok(solids) = &out {
        for &s in solids {
            if let Some(vh) = crate::transform::foreign_named_vertex(model, s) {
                return Err(crate::reject_at(
                    RejectReason::VertexNamesAbsentSurface,
                    crate::RejectWhere::Point(model.vertex_point(vh)),
                ));
            }
        }
    }
    out
}

/// The one plane class two **Pierce** end names share — the fact a curved-face ring cannot
/// read off its own surface (a panel's `surf` is the cylinder): an arc's two ends share the
/// class of the circle's plane, a ruling's two ends share the wall's. `None` is the honest
/// answer for a pair that shares none or both (a degenerate naming this ladder does not
/// arrange) — callers refuse by the ladder's name rather than unwrap.
fn shared_pierce_plane(a: NodeId, b: NodeId) -> Option<usize> {
    let (pa, _, _) = combinatorics::pierce_name(a)?;
    let (pb, _, _) = combinatorics::pierce_name(b)?;
    let mut shared = pa.iter().filter(|x| pb.contains(x));
    let c = *shared.next()?;
    shared.next().is_none().then_some(c)
}

/// An unordered edge key: the two nodes in a fixed order, so `{a,b}` and `{b,a}` collide.
/// (`pub(crate)` for the twin-match fence, which counts exactly these keys.)
pub(crate) fn norm_edge(a: NodeId, b: NodeId) -> (NodeId, NodeId) {
    if a <= b { (a, b) } else { (b, a) }
}
