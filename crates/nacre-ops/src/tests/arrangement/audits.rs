use super::*;
/// One arrangement vertex a boolean named, with **every** plane through it.
///
/// `planes` is ground truth: it is found by asking every plane in the table whether it passes
/// through the point, not by the rule the engine uses to notice concurrencies. That is the whole
/// value of it — the discovery rule can then be measured against something other than itself.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct Concurrency {
    /// The plane class whose arrangement used the name.
    pub wc: usize,
    /// The name that class used for the point.
    pub triple: NodeId,
    /// Every plane through that point, sorted. Longer than 3 exactly when the point is concurrent.
    pub planes: Vec<usize>,
    /// The sub-triples of `planes` that share a **line** rather than meeting at the point — the
    /// second face of the same degeneracy, and what a line-identity rule would have to fold. Read
    /// off `planes`, which is why one discovery answers both questions.
    pub lines: Vec<[usize; 3]>,
}

/// **Trace every class of the two solids' arrangement, past any decline** — for the ledgers.
///
/// The production driver returns at the first declined class, and in serial mode (no default
/// features) it never traces the classes after it; a lock that reads a road's ledger for a face
/// whose boolean is still refused one road later would then see different records per traversal
/// mode. This runs the tracer on every class, like [`concurrency_audit`], and keeps nothing but
/// what the roads recorded on the way.
#[cfg(test)]
pub(crate) fn trace_every_class(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(), BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings,
    );
    // The table production starts from the operands' own concurrencies and corners.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    for wc in 0..geom.len() {
        let _ = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
    }
    Ok(())
}

/// Every **concurrent** vertex (four or more planes) the two solids' arrangement would name.
///
/// Runs the per-class front half — trace, merge, split — over **all** classes and keeps going past
/// a decline, which the production driver cannot do: it returns at the first declined class, so on
/// exactly the models that have concurrencies it would see one and stop. Both halves are called,
/// not reimplemented, so what is observed is what the engine does.
#[cfg(test)]
pub(crate) fn concurrency_audit(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Concurrency>, BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    let mut out = Vec::new();
    // The table production starts from the operands' own concurrencies and corners.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        // Names this class used: segment endpoints, single-point touches, and — since a crossing
        // the arrangement mints is a vertex too — the split's endpoints where it got that far.
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let mut names: Vec<NodeId> = tr.touches.clone();
        for s in merged.iter().chain(
            split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default())
                .as_deref()
                .unwrap_or(&[])
                .iter(),
        ) {
            names.extend(s.end);
        }
        // ★ And the input the *trace's* rule actually sees: the operand faces' ring vertices.
        // Emitted segment endpoints are named `[wc, wall, r]`, so they always mention `wc` and can
        // never exercise the `wc ∉ t` branch — collecting only those measured nothing, which is
        // what the "was the rule exercised?" counter in the driver test caught.
        for (solid, inc) in [(a, &inc_a), (b, &inc_b)] {
            let face_handles: Vec<Handle<Face>> = solid_shell_handles(model, solid)
                .into_iter()
                .flat_map(|sh| model.shell(sh).faces.clone())
                .collect();
            for fh in face_handles {
                let Some(&fp) = surf_ix.get(&fh) else {
                    continue;
                };
                let mut rings =
                    combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix, &cyls)
                        .map(|lr| lr.poly().map(|nr| nr.triples.clone()).unwrap_or_default())
                        .unwrap_or_default();
                rings.extend(
                    combinatorics::hole_rings(model, fh, fp, inc, &jd, &plane_ix, &cyls)
                        .unwrap_or_default()
                        .into_iter()
                        .filter_map(|lr| lr.poly().map(|nr| nr.triples.clone()))
                        .flatten(),
                );
                // Only vertices the class would actually name: those lying on it (`side == 0`),
                // which is exactly the run condition the trace's rule fires under.
                //
                // ★ No re-canonicalization around this filter: a ring's nodes are `NodeId`s, and
                // the only constructor sorts.
                // ★★ **No `expect` on the projection.** A corner a cylinder made is *excluded* —
                // the honest answer for a four-plane concurrency hunt — instead of panicking.
                // ★★★★★ **And the exclusion is said, not left to an argument value.** Left to an
                // empty cylinder table handed to `side_of` (which declines a pierce node without
                // one), giving that call a real table "for consistency" would turn the `expect`
                // below into a panic. A concurrency is a fact about plane triples,
                // so scope is the reason and it belongs in the filter — after which `&[]` is
                // provably never read, which is what the other two sites say by wrapping their own
                // triple.
                // ☑ Measured: it drops **0** nodes across the workspace suite — this audit's own
                // corpus is plane-only, so the filter is a precondition made explicit.
                // ★ The call is one line because the source-text meta-test
                // `no_production_code_walks_a_ring_past_the_shared_walk` allow-lists this site by
                // its argument text, and it reads line by line.
                names.extend(
                    rings
                        .into_iter()
                        .filter(|&n| three_plane_name(n).is_some())
                        .filter(|&n| combinatorics::side_of(&jd, &[], n, wc) == Some(0)),
                );
            }
        }
        names.sort_unstable();
        names.dedup();
        for n in names {
            let t = three_plane_name(n).expect("a three-plane node");
            if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
                continue; // names no point, so "the planes through it" is not a question
            }
            let mut planes: Vec<usize> = t.to_vec();
            planes.extend(
                (0..geom.len())
                    .filter(|q| !t.contains(q) && jd.orient3d(t[0], t[1], t[2], *q) == 0),
            );
            if planes.len() > 3 {
                planes.sort_unstable();
                let mut lines = Vec::new();
                for i in 0..planes.len() {
                    for j in (i + 1)..planes.len() {
                        for k in (j + 1)..planes.len() {
                            let tri = [planes[i], planes[j], planes[k]];
                            if jd.plane_pair_dir_sign(tri[0], tri[1], tri[2]) == 0 {
                                lines.push(tri);
                            }
                        }
                    }
                }
                out.push(Concurrency {
                    wc,
                    triple: n,
                    planes,
                    lines,
                });
            }
        }
    }
    Ok(out)
}

/// **One operand vertex, as every producer names it** (an instrument).
///
/// A point's identity is meant to be a function of the set of planes through it. This report
/// holds, for one vertex of one operand, what the three sources say that set is and what names
/// come out: the **topology** (the plane classes of the faces incident to the vertex — what the
/// operand itself knows), the **geometry** (every class of the arrangement whose plane passes
/// through the point, asked plane by plane), and the **names** the ring road hands the tracer
/// for that vertex — one per face loop that visits it, with whether each is an independent
/// triple (three planes sharing a line name no point) and what the alias table folds it to after
/// every class has been traced.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct OperandVertexReport {
    /// Which operand (`0` = A, `1` = B) and the vertex.
    pub side: usize,
    pub vertex: Handle<Vertex>,
    pub point: [f64; 3],
    /// Plane classes of the incident faces, sorted; cylinder classes are not counted.
    pub topo: Vec<usize>,
    /// Every class through the point, from the first independent name; empty when no name is
    /// independent (then no point exists to ask about).
    pub geom: Vec<usize>,
    /// `(face class, name)` per loop visit; `dependent[i]` says the name's planes share a line.
    pub names: Vec<(usize, NodeId)>,
    pub dependent: Vec<bool>,
    /// `names[i]` after the alias fold of a full trace over every class.
    pub folded: Vec<NodeId>,
    /// Whether `folded[i]` names a point on every plane of `topo` — a fold that lands on another
    /// point is the silent shape of a naming defect (a dependent fold counts as off).
    pub folded_on_vertex: Vec<bool>,
}

/// **A pierce corner of an operand seen from a class plane that passes through it** (an
/// instrument). The point has more names than the corner's own: the three-plane name of its two
/// planes with the class, and the class's ruling crossing with the corner's cap plane at whichever
/// root is this point. Whether the alias table knows they are one point is what this reports.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct PierceCornerReport {
    pub side: usize,
    pub vertex: Handle<Vertex>,
    pub point: [f64; 3],
    /// The corner's own name (its two planes, its cylinder, its root).
    pub corner: NodeId,
    /// The class plane through the point that is not one of the corner's own.
    pub class: usize,
    /// The candidate names for the same point: the corner, the three-plane name, and the two
    /// ruling-crossing names (`side ±1`) of `class` with the corner's planes — one of the two is
    /// this point, the other the class's other ruling on the same cap.
    pub candidates: Vec<NodeId>,
    /// `candidates` after the alias fold of a full trace over every class.
    pub folded: Vec<NodeId>,
}

/// One declined (class, face, kind) of a full trace — [`trace_declines`]' row.
#[cfg(test)]
pub(crate) type ClassDecline = (usize, Option<Handle<Face>>, DeclineKind);

/// Every declined (class, face) of a full trace over every class — the audits' companion, since
/// the production driver stops at the first.
#[cfg(test)]
pub(crate) fn trace_declines(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<ClassDecline>, BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    let mut out = Vec::new();
    // The table production starts from the operands' own concurrencies and corners.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        for &(fp, kind) in &tr.declined {
            out.push((wc, faces_tab[fp].face(), kind));
        }
    }
    Ok(out)
}

/// Every pierce corner of `a` and `b` that lies on a class plane not its own, as
/// [`PierceCornerReport`]s — the cylinder twin of [`operand_vertex_audit`].
#[cfg(test)]
pub(crate) fn pierce_corner_audit(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<PierceCornerReport>, BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        aliases.absorb(&tr.aliases);
        let merged = merge_coincident(&jd, &tr.segs, wc, &aliases);
        let _ = split_at_crossings(&jd, &cyls, wc, &merged, &mut aliases);
    }
    // The corners: every pierce name a plane face's ring carries, once per vertex.
    let mut out = Vec::new();
    for (side, (solid, inc)) in [(a, &inc_a), (b, &inc_b)].into_iter().enumerate() {
        let mut seen: Vec<Handle<Vertex>> = Vec::new();
        let face_handles: Vec<Handle<Face>> = solid_shell_handles(model, solid)
            .into_iter()
            .flat_map(|sh| model.shell(sh).faces.clone())
            .collect();
        for fh in face_handles {
            let Some(&fp) = surf_ix.get(&fh) else {
                continue;
            };
            let Ok(lr) =
                combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix, &cyls)
            else {
                continue;
            };
            let Some(nr) = lr.poly() else {
                continue;
            };
            let face = model.face(fh);
            if nr.triples.len() != face.outer.half_edges.len() {
                continue;
            }
            for (he, &n) in face.outer.half_edges.iter().zip(nr.triples.iter()) {
                let Some((planes, cyl, _root)) = combinatorics::pierce_name(n) else {
                    continue;
                };
                let vh = crate::he_start(model, *he);
                if seen.contains(&vh) {
                    continue;
                }
                seen.push(vh);
                let def = &cyls[cyl].def;
                for class in 0..geom.len() {
                    if planes.contains(&class) {
                        continue;
                    }
                    if combinatorics::side_of(&jd, &cyls, n, class) != Some(0) {
                        continue;
                    }
                    let mut cands = vec![n];
                    let mut s = vec![planes[0], planes[1], class];
                    s.sort_unstable();
                    if let Some(t) = combinatorics::canonical_triple(&jd, &s) {
                        cands.push(NodeId::three_planes(t));
                    }
                    for fc in planes {
                        for side_r in [1i8, -1i8] {
                            // `crossing_on_ruling` takes the cap (normal ∥ axis) first and the
                            // class holding the axis second; the corner's wall plane errs out.
                            if let Ok(id) = crossing_on_ruling(&jd, def, fc, class, cyl, side_r) {
                                if !cands.contains(&id) {
                                    cands.push(id);
                                }
                            }
                        }
                    }
                    let folded = cands.iter().map(|&c| aliases.canon_point(c)).collect();
                    let p = model.vertex_point(vh);
                    out.push(PierceCornerReport {
                        side,
                        vertex: vh,
                        point: [p[0], p[1], p[2]],
                        corner: n,
                        class,
                        candidates: cands,
                        folded,
                    });
                }
            }
        }
    }
    Ok(out)
}

/// Every plane-only operand vertex of `a` and `b`, as [`OperandVertexReport`]s — the
/// audit of the ring road's naming rule against the topology and the geometry.
///
/// Like [`concurrency_audit`] it drives the real front half (trace, merge, split) over **every**
/// class and keeps going past a decline, accumulating the alias table so the fold it reports is
/// the one production would reach. Vertices a cylinder touches are skipped: their names are
/// pierce points, a different vocabulary.
#[cfg(test)]
pub(crate) fn operand_vertex_audit(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<OperandVertexReport>, BoolError> {
    use crate::planes::ClassIx;
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    // The alias fold production would reach: every class traced, merged and split, discoveries
    // accumulated across classes.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        aliases.absorb(&tr.aliases);
        let merged = merge_coincident(&jd, &tr.segs, wc, &aliases);
        let _ = split_at_crossings(&jd, &cyls, wc, &merged, &mut aliases);
    }
    let plane_class = |k: usize| match plane_ix[k] {
        ClassIx::Plane(c) => Some(c),
        ClassIx::Cyl(_) => None,
    };
    let mut out = Vec::new();
    for (side, (solid, inc)) in [(a, &inc_a), (b, &inc_b)].into_iter().enumerate() {
        // Topology: vertex → plane classes of its incident faces; a vertex on any cylinder face
        // is left out (its name is a pierce point).
        let mut topo: HashMap<Handle<Vertex>, Vec<usize>> = HashMap::new();
        let mut curved: std::collections::HashSet<Handle<Vertex>> = Default::default();
        for (bounds, pair) in inc.edges() {
            for &vh in bounds {
                for &k in pair {
                    match plane_class(k) {
                        Some(c) => topo.entry(vh).or_default().push(c),
                        None => {
                            curved.insert(vh);
                        }
                    }
                }
            }
        }
        // Names: per face loop, the ring road's triples, position `i` being the start of
        // half-edge `i` (a seam joint pushes no name, so a loop whose lengths differ is skipped).
        let mut names: HashMap<Handle<Vertex>, Vec<(usize, NodeId)>> = HashMap::new();
        let face_handles: Vec<Handle<Face>> = solid_shell_handles(model, solid)
            .into_iter()
            .flat_map(|sh| model.shell(sh).faces.clone())
            .collect();
        for fh in face_handles {
            let Some(&fp) = surf_ix.get(&fh) else {
                continue;
            };
            let Some(fc) = plane_class(fp) else {
                continue;
            };
            let face = model.face(fh);
            let mut loops: Vec<(Vec<nacre_topo::HalfEdge>, Option<Vec<NodeId>>)> = Vec::new();
            loops.push((
                face.outer.half_edges.clone(),
                combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix, &cyls)
                    .ok()
                    .and_then(|lr| lr.poly().map(|nr| nr.triples.clone())),
            ));
            if let Ok(holes) = combinatorics::hole_rings(model, fh, fp, inc, &jd, &plane_ix, &cyls)
            {
                for (lp, lr) in face.inner.iter().zip(holes) {
                    loops.push((
                        lp.half_edges.clone(),
                        lr.poly().map(|nr| nr.triples.clone()),
                    ));
                }
            }
            for (hes, triples) in loops {
                let Some(triples) = triples else {
                    continue;
                };
                if triples.len() != hes.len() {
                    continue;
                }
                for (he, &n) in hes.iter().zip(triples.iter()) {
                    let vh = crate::he_start(model, *he);
                    names.entry(vh).or_default().push((fc, n));
                }
            }
        }
        let mut vertices: Vec<Handle<Vertex>> = topo.keys().copied().collect();
        vertices.sort_by_key(|v| v.index());
        for vh in vertices {
            if curved.contains(&vh) {
                continue;
            }
            let mut t = topo[&vh].clone();
            t.sort_unstable();
            t.dedup();
            let vnames = names.get(&vh).cloned().unwrap_or_default();
            let dependent: Vec<bool> = vnames
                .iter()
                .map(|(_, n)| match three_plane_name(*n) {
                    Some(tr) => jd.plane_pair_dir_sign(tr[0], tr[1], tr[2]) == 0,
                    None => false,
                })
                .collect();
            let geom_set =
                vnames
                    .iter()
                    .zip(&dependent)
                    .find(|(_, d)| !**d)
                    .and_then(|((_, n), _)| three_plane_name(*n))
                    .map(|tr| {
                        let mut g: Vec<usize> = tr.to_vec();
                        g.extend((0..geom.len()).filter(|q| {
                            !tr.contains(q) && jd.orient3d(tr[0], tr[1], tr[2], *q) == 0
                        }));
                        g.sort_unstable();
                        g
                    })
                    .unwrap_or_default();
            let folded: Vec<NodeId> = vnames
                .iter()
                .map(|(_, n)| aliases.canon_point(*n))
                .collect();
            let folded_on_vertex = folded
                .iter()
                .map(|n| match three_plane_name(*n) {
                    Some(tr) => {
                        jd.plane_pair_dir_sign(tr[0], tr[1], tr[2]) != 0
                            && t.iter().all(|&q| {
                                tr.contains(&q) || jd.orient3d(tr[0], tr[1], tr[2], q) == 0
                            })
                    }
                    None => false,
                })
                .collect();
            let p = model.vertex_point(vh);
            out.push(OperandVertexReport {
                side,
                vertex: vh,
                point: [p[0], p[1], p[2]],
                topo: t,
                geom: geom_set,
                names: vnames,
                dependent,
                folded,
                folded_on_vertex,
            });
        }
    }
    Ok(out)
}

/// One plane class's **label-frame audit** (family #2 diagnostic): what each producer says about
/// "above" on this class, plus how far the per-class pipeline gets. `#[cfg(test)]`, `pub(crate)`
/// so the driver test can live in `crate::tests` where the two-solid fixtures are (the same reason
/// [`boolean`] is `pub(crate)`).
#[cfg(test)]
#[derive(Debug)]
pub(crate) struct ClassAudit {
    pub wc: usize,
    /// The root plane resolved to real coordinates — a point on it and its stored normal. Plane
    /// indices alone have twice carried a wrong geometric story into a write-up.
    pub root_point: [f64; 3],
    pub root_normal: [f64; 3],
    /// `sign(stored normal · n_out)` — `-1` exactly when the class root face is `Reversed`, which
    /// is when the stored-normal label frame and `side_of`'s outward frame disagree.
    pub orient_sign: i8,
    /// `body_above` of every seated segment, and `body_above` of every graze segment.
    pub seated: Vec<bool>,
    pub grazes: Vec<bool>,
    pub transversals: usize,
    pub declined: Vec<(usize, DeclineKind)>,
    /// The reason the per-class pipeline raised, if any (`None` = the class went through).
    pub failed_at: Option<RejectReason>,
    /// **What the arrangement produced, before anything refused it** — cells walked, arcs the
    /// split cut, and the nesting's two counts.
    ///
    /// ★★★ Born while the arc stopper swallowed every arrangement stage's answer (a probe was
    /// the only way to see the population, and a probe is deleted before the commit); the
    /// population is green now, and this stays as the audit's direct read of what a class
    /// produced. `None` when the class stopped before those stages ran.
    pub produced: Option<Produced>,
    /// **Every emitted face's outer ring, as coordinates** — one `Vec` per emitted face whose outer
    /// bound is a polygon ring (a face bounded by an uncut circle has no nodes and contributes
    /// nothing), sorted by their rotation-normalized first coordinate. `None` when the class
    /// stopped before `emit_faces`.
    ///
    /// ★★★★★ **Coordinates, not node identities, and the reason is that the lock must not be
    /// circular.** A node name is `ThreePlane([0, 2, 5])` — plane **class** indices, which no
    /// fixture derives; writing the expected value means running the engine and copying what it
    /// said, and a test whose oracle is its subject measures nothing. Coordinates come from the
    /// fixture's own numbers.
    ///
    /// ★★★★★ **Rotation is normalized, reversal deliberately is not.** The walk tries both
    /// handednesses and pushes cells in its own order, so *which* vertex a ring starts at is not a
    /// fact about the geometry — the minimum coordinate goes first. Reversal *is* a fact: it is
    /// what a wrong `sense` on a split segment produces, and it is invisible in every
    /// order-independent summary the previous rungs measured (cell counts, nesting, sorted
    /// labels). Folding it here would delete the one thing this field exists to see.
    ///
    /// ★★★ **And rotation-normalization is itself reversal-blind at `n <= 2`**: reversing a
    /// two-cycle *is* a rotation of it. Such a ring cannot carry this lock at all — which is why
    /// the fence's red probe is read on the turned boss (5- and 3-rings), not the straddling one
    /// (6- and **2**-rings).
    ///
    /// ★ It is a sibling of [`Produced`] rather than a field in it because `[f64; 3]` is not `Eq`
    /// (the derive would break) and one of these coordinates is irrational, so it needs `near()`
    /// rather than `==` — the two could not ride the same `assert_eq!` regardless.
    ///
    /// ★★★ **Two things are left out on purpose, and saying so is the point** — a silent omission
    /// reads as an oversight to whoever comes next:
    /// - **`flip`.** It is `keep_above == (frame_sign > 0)`, so it turns on the *class's stored
    ///   normal direction* — an implementation fact, not one of the fixture's numbers. Carrying it
    ///   would put one bit in here that a fence could only fill by copying a run, and reversal —
    ///   the thing this field exists for — is already caught by the sequence itself.
    /// - **Inner rings.** A face's holes are not carried. Today's arc classes have none
    ///   (`nesting.holes` is empty, measured), so there is nothing to lose *yet*; an arc class with
    ///   a hole would go unchecked here, and that is the day this grows.
    pub outer_rings: Option<Vec<Vec<[f64; 3]>>>,
}

/// The per-class arrangement's output, for [`ClassAudit`].
#[cfg(test)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Produced {
    pub cells: usize,
    pub arcs: usize,
    pub roots: usize,
    pub holes: usize,
    /// Every `+1` cell's label, **sorted** — the labelling stage's answer, in a form the walk's
    /// handedness cannot reorder.
    ///
    /// ★★ `Copy` goes with this, and that was the choice rather than the discovery: a sibling field
    /// on [`ClassAudit`] would keep it, but then one class's two facts live in two places and the
    /// fence has to join them by index — a seam where they can drift. `PartialEq` stays, so the
    /// fence remains a single `assert_eq!` against a value derived from the geometry.
    pub pos_labels: Vec<Label>,
}

/// Lexicographic order on realized coordinates. `total_cmp` rather than `partial_cmp`: a node that
/// failed to realize comes through as `NaN` below, and this must still be a total order.
#[cfg(test)]
fn cmp_pt(a: &[f64; 3], b: &[f64; 3]) -> std::cmp::Ordering {
    a.iter()
        .zip(b)
        .map(|(x, y)| x.total_cmp(y))
        .find(|o| o.is_ne())
        .unwrap_or(std::cmp::Ordering::Equal)
}

/// Every emitted face's outer ring as coordinates, for [`ClassAudit::outer_rings`] — that field
/// carries the argument for coordinates, for normalizing rotation, and for leaving reversal alone.
///
/// ★★ **A ring here can need both roads to a coordinate** — three of its nodes are plane triples
/// and two are pierce points — and the dispatch between them lives in
/// [`combinatorics::node_point_f64`], not here: spelling a [`combinatorics::NodeId`] variant
/// outside that file is what `three_plane_name`'s gate forbids, and the first draft of this
/// function broke it with the suite green.
///
/// ★ A node that cannot be realized becomes `NaN`, not a dropped element: the ring keeps its
/// length, so the fence reads "this vertex had no coordinate" instead of "the ring is short".
#[cfg(test)]
fn face_ring_coords(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
) -> Vec<Vec<[f64; 3]>> {
    let point = |n: NodeId| combinatorics::node_point_f64(jd, cyls, n).unwrap_or([f64::NAN; 3]);
    let mut out: Vec<Vec<[f64; 3]>> = faces
        .iter()
        .filter_map(|f| f.outer.ring())
        .map(|r| {
            let pts: Vec<[f64; 3]> = r.nodes.iter().map(|&n| point(n)).collect();
            // Rotation only. Where the walk started is not a fact about the geometry; which way it
            // went is, and it is the one this instrument was added to see.
            match pts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| cmp_pt(a, b))
                .map(|(i, _)| i)
            {
                Some(i) => pts[i..].iter().chain(&pts[..i]).copied().collect(),
                None => pts,
            }
        })
        .collect();
    // Between faces the order *is* the walk's, so it is normalized — lexicographically over the
    // whole (already rotated) sequence.
    //
    // ★ **Total on purpose, not just "by the first vertex".** Two faces of one class can share
    // their minimum vertex — they meet there — and ordering on that alone would leave such a pair
    // in walk order, which is a handedness coin flip inside an instrument whose whole job is to be
    // stable. Today's two fixtures have distinct minima, so this costs nothing and removes the
    // class of flake rather than relying on the fixture to avoid it.
    out.sort_by(|a, b| {
        a.iter()
            .zip(b)
            .map(|(x, y)| cmp_pt(x, y))
            .find(|o| o.is_ne())
            .unwrap_or_else(|| a.len().cmp(&b.len()))
    });
    out
}

/// Audit every plane class of one boolean input pair. Runs the same per-class pipeline as
/// [`trace_result_faces`] but never aborts, so one failing class does not hide the rest.
#[cfg(test)]
pub(crate) fn frame_audit(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<ClassAudit>, BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        cyls,
        crossings,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    // The alias fixpoint, exactly as the boolean runs it (sequentially): classes discover names
    // for one another's features, so auditing each class against an empty table is auditing a
    // *different* pipeline — it diverged from the boolean's answer the day the rounds arrived,
    // and the divergence surfaced when this file's test moved to a fixture that needs them.
    let mut aliases = Aliases::default();
    loop {
        let before = aliases.len();
        for wc in 0..geom.len() {
            let mut tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
            aliases.absorb(&std::mem::take(&mut tr.aliases));
            if tr.declined.is_empty() {
                let merged = merge_coincident(&jd, &tr.segs, wc, &aliases);
                let _ = split_at_crossings(&jd, &cyls, wc, &merged, &mut aliases);
            }
        }
        if aliases.len() == before {
            break;
        }
    }
    let mut out = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        let side = |f: fn(&SegKind) -> Option<bool>| -> Vec<bool> {
            tr.segs.iter().filter_map(|s| f(&s.kind)).collect()
        };
        let mut audit = ClassAudit {
            wc,
            root_point: geom[wc].tri[0].as_array(),
            root_normal: geom[wc].plane.normal().as_array(),
            orient_sign: geom[wc].frame_sign,
            seated: side(|k| match k {
                SegKind::Seated { body_above } => Some(*body_above),
                _ => None,
            }),
            grazes: side(|k| match k {
                SegKind::Graze { body_above } => Some(*body_above),
                _ => None,
            }),
            transversals: tr
                .segs
                .iter()
                .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
                .count(),
            declined: tr.declined.clone(),
            failed_at: None,
            produced: None,
            outer_rings: None,
        };
        // Run the rest of the per-class pipeline, recording where it stops.
        audit.failed_at = if let Some(&(fp, kind)) = audit.declined.first() {
            Some(decline_to_reject(kind, faces_tab[fp].face()))
        } else {
            let mut produced = None;
            let mut outer_rings = None;
            let mut run = || -> Result<(), BoolError> {
                // A copy, so one class's split discoveries cannot leak into the next class's
                // audit — the fixpoint above already holds everything the boolean would know.
                let mut local = aliases.clone();
                let merged = merge_coincident(&jd, &tr.segs, wc, &local);
                let split = split_at_crossings(&jd, &cyls, wc, &merged, &mut local)?;
                let split = drop_newsless(split)?;
                let circles = merge_circles(&tr.circles, &cyls, &local)?;
                let rulings = merge_rulings(&tr.rulings, &cyls, &local);
                // ★★★★ **The same edges the boolean uses.** This copy of the pipeline is what
                // makes the audit an instrument; feeding it un-split edges would let it run to
                // the end and report `failed_at: None` for an input the boolean refuses —
                // *"the worst possible time to be lying"* (`decline_to_reject`). The arc fence
                // in `bands.rs` locks the two together as «both
                // succeed» (a stopper plugged into the socket would raise here per class while
                // the boolean defers it to the assembly's end; same classes, same name).
                let edges = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &rulings, &local)?;
                let staged = per_class(&jd, &cyls, kind, wc, &edges);
                if let Ok(s) = &staged {
                    let mut pos_labels: Vec<Label> = s
                        .cells
                        .iter()
                        .zip(&s.labels)
                        .filter(|(c, _)| c.winding == 1)
                        .map(|(_, &l)| l)
                        .collect();
                    pos_labels.sort_unstable();
                    produced = Some(Produced {
                        cells: s.cells.len(),
                        arcs: edges.arcs.len(),
                        roots: s.nesting.root_groups.len(),
                        holes: s.nesting.holes.len(),
                        pos_labels,
                    });
                    outer_rings = Some(face_ring_coords(&jd, &cyls, &s.faces));
                }
                let _ = staged?;
                Ok(())
            };
            let stopped = run().err().and_then(|e| match e {
                BoolError::Rejected { reason, .. } => Some(reason),
                // No live-set check runs inside the pipeline, so this arm is unreachable.
                BoolError::InputNotLive => None,
            });
            audit.produced = produced;
            audit.outer_rings = outer_rings;
            stopped
        };
        out.push(audit);
    }
    Ok(out)
}
