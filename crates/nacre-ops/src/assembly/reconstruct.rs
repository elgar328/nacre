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
    rims: &crate::draft::HeldRims,
    deferred: Option<BoolError>,
    tangencies: Tangencies<'_>,
    memo: &mut crate::realize::PlaneMemo,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes = jd.planes;
    // No faces means no result — `Common` of two solids that miss each other, `Cut` of a box that
    // is wholly inside what cuts it. That is an answer, not a failure: a solid is bounded by faces,
    // so a non-empty result cannot have none. The inputs are still consumed, exactly as they are on
    // any other successful boolean — this returns `Ok`, so the caller's retire runs.
    if faces.is_empty() {
        return Ok(Vec::new());
    }
    let named = name_result_vertices(jd, seam, faces, cyls);
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
    // The per-solid straight-angle pass may have dropped nodes of its own — per solid, so a rim
    // can be whole for one body and cut for another it touches; everything below asks the rims
    // as these faces hold them, with the group.
    let cut = rims.cut_per_group(faces, &group_of);
    // ★ **The tangency verdict stands beside the self-touch sieve, and for the same reason** —
    // both are whole-result judgements that need the grouping and must speak *before* a handle is
    // minted, so a refusal leaves the arena as it found it. A held grouping error stays held: the
    // question "do these two lumps share a solid" is only askable of a grouping that answered.
    if let Ok(g) = &grouping {
        tangency_reject(&tangencies, g)?;
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
            let handle = {
                let sv = seam
                    .iter()
                    .find(|s| s.triple == node)
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                // A vertex that is a corner of no face at all has no name in the result's own
                // planes — a degeneracy, and the honest answer is the one this reason already
                // carries ("a corner with no turn").
                let Some(def) = def_triple.get(&(g, node)) else {
                    return Err(reject(RejectReason::StraightAngle));
                };
                // ★★ **The seam table already asked this vertex's question.** It realized the
                // node's own definition through the push funnel's rule (`realize::point_cache`),
                // so where the naming defines the vertex by that same definition its answer is
                // pushed as it stands; asking twice paid a second realization on every vertex
                // (fold 40, rotated: 211 → 237 ms). Where four planes concur the naming may
                // pick another triple through the same point, and that definition is realized
                // afresh, with the seam's coordinate as its figure.
                let v = def.vertex(planes, cyls);
                if Def::of_name(node) == Some(*def) {
                    crate::realize::push_vertex_asked(model, v, sv.cache)
                } else {
                    crate::realize::push_vertex_realized(
                        model,
                        v,
                        PointCache::Unrealized {
                            coord: sv.cache.coord(),
                        },
                        crate::realize::ChainLink::Fresh,
                    )
                }
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

    // ★★ **The rim table**: a whole circle has no node and no triple, so it is not welded through
    // the seam table at all — it is minted here, once per `(group, cylinder class, plane class)`,
    // with its `OnSeam` vertex, and every face that meets that rim asks for the same pair back.
    // That sharing is what makes the rim edge's use count two (the cap face once, the lateral
    // once) rather than two separate edges the closed-shell guard would call dangling. A cut
    // circle has no entry: its boundary is its arcs between nodes, minted by the rings that walk
    // them.
    //
    // The group is in the key for the reason it is in `vh`'s: two result solids must not share a
    // handle.
    let mut rim: RimTable = HashMap::new();
    let rims_built = (|| -> Result<(), BoolError> {
        for (fi, lf) in faces.iter().enumerate() {
            let g = group_of[fi];
            let keys: Vec<(usize, usize)> = std::iter::once(&lf.outer)
                .chain(lf.inner.iter())
                .filter_map(|b| match (b, lf.surf) {
                    (Bound::Circle { cyl }, ClassIx::Plane(c)) => Some((*cyl, c)),
                    (Bound::Rim { plane, .. }, ClassIx::Cyl(k)) => Some((k, *plane)),
                    _ => None,
                })
                .collect();
            for (k, c) in keys {
                if rim.contains_key(&(g, k, c)) {
                    continue;
                }
                // ★ A circle this solid still cuts cannot bound it as a whole — one of its faces
                // on that plane holds one of its nodes ([`HeldRims::cut_per_group`]) — so reaching
                // here with one is the split and the faces disagreeing, named as that.
                if cut.contains(&(g, k, c)) {
                    return Err(reject(RejectReason::CylinderStagesDisagree));
                }
                let (lat, plane) = (cyls[k].surf, planes[c].surf);
                // The seam point of this rim, spelled as a circle prism spells one: the axis meets
                // the plane at the circle's centre, and `θ = 0` is the `+ref_dir` side of it.
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
                        // cutting plane of a rim key. The reject is kept rather than
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
                let v = crate::realize::push_vertex_realized(
                    model,
                    Vertex::OnSeam([lat, plane]),
                    PointCache::Unrealized { coord: point },
                    crate::realize::ChainLink::Fresh,
                );
                let e = crate::realize::push_edge_realized(model, [lat, plane], [v, v], memo)
                    .map_err(edge_refused)?;
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
    // ★ **A ruling edge reads its carriers the same way, from the same table.** A straight
    // edge is its two ends (`EdgeKey`), whichever vocabulary a face states it in, so a ruling and
    // a plane-pair line between one vertex pair are one edge and one entry here: a panel pushes
    // its cylinder, a plane face its plane. Deriving a ruling's plane from the two ends' *names*
    // (the plane they share) states the line only while a point's name is the pair that minted
    // it: once the alias table represents a tangent corner by another pair — a class through the
    // axis crossing the cylinder there — the shared plane is that class, which contains the line
    // but does not bound the edge. And a ruling no panel uses at all — two planes meeting on the
    // cylinder's line, the lateral gone from both sides — is carried by those two planes, which
    // no rule starting from the cylinder can state.
    let mut pair_surfs: HashMap<(usize, usize), Vec<Handle<Surface>>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        // A whole circle contributes no node edge — its rim is minted with its carriers stated
        // (`[lateral, plane]`); what a curved face does walk as a ring (a chain, a panel, a hole)
        // is read here like any ring.
        let fsurf = match lf.surf {
            ClassIx::Plane(fc) => planes[fc].surf,
            ClassIx::Cyl(k) => cyls[k].surf,
        };
        for r in lf.poly_rings() {
            let k = r.nodes.len();
            for t in 0..k {
                let (va, vb) = (vh[&(g, r.nodes[t])], vh[&(g, r.nodes[(t + 1) % k])]);
                let pair = unordered(va.index() as usize, vb.index() as usize);
                match r.walls[t] {
                    Wall::Plane(_) | Wall::Ruling { .. } => {
                        pair_surfs.entry(pair).or_default().push(fsurf);
                    }
                    // An arc edge shares its vertex pair with the chord between the same pierce
                    // vertices, and pushing this face's surface here would pollute the chord's
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
    // below the floor (`overview.md`: out-of-coverage input declines honestly). `SeamAlias`
    // catches the known way this happens — two triples on one point — at the seam table, where the
    // names are still in hand, so this is a backstop with no firing test (cf. `NonManifoldEdge`).
    let mut edge_for = |model: &mut Model,
                        memo: &mut crate::realize::PlaneMemo,
                        va: Handle<Vertex>,
                        vb: Handle<Vertex>,
                        wall: Wall,
                        face_surf: Handle<Surface>|
     -> Result<Handle<Edge>, BoolError> {
        let wall_surf = match wall {
            Wall::Plane(w) => planes[w].surf,
            Wall::Ruling { cyl, .. } => cyls[cyl].surf,
            Wall::Arc { cyl, ccw, plane } => {
                // ★★★ **An arc edge is minted in CCW order** — `[A, B]` is the piece from A to
                // B counter-clockwise about the axis, so the two complementary arcs between one
                // pierce pair are `[A, B]` and `[B, A]`: the vertex order is the last bit the
                // endpoints alone cannot give (`EdgeKey`'s note). The carriers are the wall's own
                // — the cylinder and the plane of the rim the arc lies on — so there is nothing for
                // the scan to read, which also keeps the chord's line key clean (`pair_surfs`
                // skips arc edges for the same reason).
                //
                // ★★ **Not the face's surface**: a lateral's is its cylinder, and an arc stated off
                // it would be `[cyl, cyl]`. The wall names the rim's plane on a cap and a lateral
                // alike, so whichever face mints the arc states one pair — and the second face to
                // walk it must name the plane the edge already carries, or the planar arrangement
                // and the cylinder chart disagree about which rim the arc lies on.
                let (from, to) = if ccw { (va, vb) } else { (vb, va) };
                let key = EdgeKey::Arc {
                    cyl,
                    from: from.index() as usize,
                    to: to.index() as usize,
                };
                let rim_plane = planes[plane].surf;
                if let Some(&e) = edge_of.get(&key) {
                    if !model.edge(e).surfaces.contains(&rim_plane) {
                        return Err(reject(RejectReason::CylinderStagesDisagree));
                    }
                    return Ok(e);
                }
                let e = crate::realize::push_edge_realized(
                    model,
                    [cyls[cyl].surf, rim_plane],
                    [from, to],
                    memo,
                )
                .map_err(edge_refused)?;
                edge_of.insert(key, e);
                return Ok(e);
            }
        };
        // ★ **A straight edge is one edge, whichever vocabulary states it**: keyed by its two
        // ends, carried by the surfaces of the two faces that use it (the scan above) — one
        // expression for a plane-pair line and a ruling alike, so which face mints it first
        // cannot change what it says. Two faces on one plane state `(P, P)`, which the result
        // check refuses as the merge that did not happen (`CoplanarMerge`); two on one cylinder
        // state `(C, C)`, which `push_edge` refuses. Manifold edges (everything a green result
        // contains) have exactly two uses; any other count is on its way to the shell guard's
        // reject, and the fallback — this face's wall and its own surface — only keeps
        // construction deterministic until then.
        let pair = unordered(va.index() as usize, vb.index() as usize);
        let key = EdgeKey::Line(pair);
        if let Some(&e) = edge_of.get(&key) {
            return Ok(e);
        }
        let surfaces = match pair_surfs[&pair][..] {
            [a, b] => Edge::carrier_pair(a, b),
            _ => Edge::carrier_pair(wall_surf, face_surf),
        };
        let e = crate::realize::push_edge_realized(model, surfaces, [va, vb], memo)
            .map_err(edge_refused)?;
        edge_of.insert(key, e);
        Ok(e)
    };

    let mut face_handles = Vec::new();
    let mut assembled: Result<(), BoolError> = Ok(());
    'faces: for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let face_surf = match lf.surf {
            ClassIx::Plane(c) => planes[c].surf,
            ClassIx::Cyl(k) => cyls[k].surf,
        };
        let mut ring = |model: &mut Model,
                        memo: &mut crate::realize::PlaneMemo,
                        r: &Ring|
         -> Result<Loop, BoolError> {
            let handles: Vec<Handle<Vertex>> = r.nodes.iter().map(|nd| vh[&(g, *nd)]).collect();
            let k = handles.len();
            let mut half_edges: Vec<HalfEdge> = Vec::with_capacity(k);
            for t in 0..k {
                let (va, vb) = (handles[t], handles[(t + 1) % k]);
                let e = edge_for(model, memo, va, vb, r.walls[t], face_surf)?;
                // ★ An arc runs the way its wall says — `ccw`, against the edge's CCW-stored
                // order. For two distinct ends that is the same as asking which end the edge starts
                // at (`edge_for` minted it `[from, to]` by that bit); for a rim cut at one point,
                // the arc from the node back to itself, the ends are one vertex and say neither.
                let forward = match r.walls[t] {
                    Wall::Arc { ccw, .. } => ccw,
                    Wall::Plane(_) | Wall::Ruling { .. } => model.edge(e).vertices[0] == va,
                };
                half_edges.push(HalfEdge { edge: e, forward });
            }
            Ok(Loop { half_edges })
        };
        // ★★ **Which way a rim circle is walked.** `derive_edge_curve` builds it on the
        // cylinder cache's frame (`Circle::from_unit_frame(centre, axis, ref_dir, r)`), so its
        // parameter runs **CCW about the cylinder's axis** — the cache's, which runs the same way
        // as the statement's `def.dir()` (the cylinder door asserts it). A loop must run CCW about
        // its own face's outward normal, so a bound on a face whose normal agrees with the axis is
        // walked forward and one whose normal opposes it backward — exactly how a circle prism's
        // caps walk their rims (top cap forward, base cap backward). A *hole* runs the other way
        // again, because an inner loop keeps the material on its left by winding against the outer.
        //
        // The face's outward normal is the way its plane faces when the face is `Forward` and the
        // reverse when `Reversed` — the `flip` decision made below — so the sign is read off
        // `frame_sign` and `flip` here rather than carried in from anywhere, and the facing read
        // against the axis is the **unflipped** one: `flip` is applied once to every kind of bound
        // right below, and reading it here too would toggle the winding twice.
        //
        // ★★★★★ **`world_rat` alone cannot answer this; its sense can.** Every other axial
        // decision on this road asks *where* a plane crosses the axis — a question about the plane,
        // whose answer does not depend on which way its coefficients are written. This one asks
        // which way the class's **frame** points relative to the axis. `world_rat` is the plane's
        // **name**, canonicalised, and carries no direction: measured over the suite, `world_rat ·
        // m` came out positive on all 436 calls while this sign varied, so reading the name alone
        // would have inverted the winding on 203 of them. The name **with its sense**
        // ([`crate::planes::plus_t_is_above`]) is the frame, exactly — read off the truth, where the
        // `f64` plane and axis caches it replaced were two rounded images of it.
        let axis_up = |k: usize| -> Result<bool, BoolError> {
            let c = lf.surf.plane();
            let world = planes[c]
                .world
                .as_ref()
                .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
            Ok(crate::planes::plus_t_is_above(world, &cyls[k].def) == (planes[c].frame_sign > 0))
        };
        let circle_loop =
            |model: &mut Model, cyl: usize, cls: usize, hole: bool| -> Result<Loop, BoolError> {
                let (_, e) = *rim
                    .get(&(g, cyl, cls))
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                let _ = model;
                let mut forward = axis_up(cyl)?;
                if hole {
                    forward = !forward;
                }
                Ok(Loop {
                    half_edges: vec![HalfEdge { edge: e, forward }],
                })
            };
        // One bound's loop, unflipped (`flip` is applied once to every kind of bound below). A
        // lateral's whole rim is its closed edge walked the way the chart walked it: `ccw` — the
        // face above it, the lower rim — runs with the edge's own CCW parametrization.
        let mut build = |model: &mut Model, b: &Bound, hole: bool| -> Result<Loop, BoolError> {
            match b {
                Bound::Ring(r) => ring(model, memo, r),
                // ★ **Each whole circle is spelled by its face** — the cap's on a plane, the rim
                // on a lateral (`dissolve_straight_angles`'s `closed`). On the other face kind it
                // is an earlier stage stating the wrong fact, and is named as that.
                Bound::Circle { cyl } => {
                    let c = match lf.surf {
                        ClassIx::Plane(c) => c,
                        ClassIx::Cyl(_) => {
                            return Err(reject(RejectReason::CylinderStagesDisagree));
                        }
                    };
                    circle_loop(model, *cyl, c, hole)
                }
                Bound::Rim { plane, ccw } => {
                    let k = match lf.surf {
                        ClassIx::Cyl(k) => k,
                        ClassIx::Plane(_) => {
                            return Err(reject(RejectReason::CylinderStagesDisagree));
                        }
                    };
                    let &(_, edge) = rim
                        .get(&(g, k, *plane))
                        .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                    Ok(Loop {
                        half_edges: vec![HalfEdge {
                            edge,
                            forward: *ccw,
                        }],
                    })
                }
            }
        };
        let mut inner: Vec<Loop> = match lf
            .inner
            .iter()
            .map(|b| build(model, b, true))
            .collect::<Result<Vec<_>, BoolError>>()
        {
            Ok(v) => v,
            Err(e) => {
                assembled = Err(e);
                break 'faces;
            }
        };
        let mut outer = match build(model, &lf.outer, false) {
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
        // The plane's frame *is* the root face's orientation: `frame_sign` carries that face's
        // `Forward`/`Reversed` as a sign (read off the stored flag) — not a face field indexed by
        // a plane (`planes[plane_idx].orient`), the shape of every bug this split exists to
        // prevent.
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
    // ★ **The grouping decided at the top of this function, raised here.** What it says:
    // one component is the whole result; several mean either an enclosed void (a cavity, an
    // inward-oriented shell) or a severed operand (two or more material-enclosing shells), and a
    // surviving cavity belongs to the piece whose outer shell nests it. See [`group_faces`] — and
    // for why the *group*, not the component, is the unit the handles above were minted per.
    // The grouping's failure yields to the deferred stopper like every stage before it; the
    // raise itself stands at the very end, past the solid assembly below.
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
    // transform and replay stand on. It is a shipped check, not a debug_assert: a vertex whose
    // definition names a surface with no face in the result is a wrong name with the right
    // volume, which a release build would otherwise ship silently (no input in the reject-trace
    // sweeps reaches it — the coplanar merge dissolves such corners). The garbage-solid residue
    // is the same class every late reject leaves (see the raise above).
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

/// An unordered edge key: the two nodes in a fixed order, so `{a,b}` and `{b,a}` collide.
/// (`pub(crate)` for the twin-match fence, which counts exactly these keys.)
pub(crate) fn norm_edge(a: NodeId, b: NodeId) -> (NodeId, NodeId) {
    if a <= b { (a, b) } else { (b, a) }
}

/// **What an edge the assembly asked for and could not have means here** — one reject per
/// [`nacre_topo::EdgeDecline`] cause. Every edge the assembly pushes joins two distinct nodes
/// on classes the gate admitted, so each cause is a broken promise, named for the promise.
fn edge_refused(d: nacre_topo::EdgeDecline) -> crate::BoolError {
    use nacre_topo::EdgeDecline::*;
    reject(match d {
        Coincident => RejectReason::ZeroLengthEdge,
        Oblique | Unstated => RejectReason::ObliqueCircleClass,
        TwoCylinders | Degenerate => RejectReason::EdgeCurveUnderived,
    })
}
