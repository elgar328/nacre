use super::*;
/// **A result face as the ray reads it** — its surface's truth and every boundary in the
/// engine's own vocabulary.
///
/// ★ **Every boundary travels** — not `poly_rings()` only, which would silently drop a face's
/// circular and banded bounds, and a ray that counts an incomplete component answers confidently
/// and wrongly.
///
/// A **lateral** face's bounds become one list of loops on its cylinder's chart
/// ([`combinatorics::LateralLoop`]): a band's rims (a whole circle by its class, a chain by its
/// edges), a panel's ring, and the holes — outer and holes together, because on the chart's
/// annulus a wrapping loop has no inside and the face is a parity over all of them.
/// A lateral ring's corners are pierce names (the region walk's stations and rim nodes), so its
/// edges need no plane class to pin a three-plane end and refuse one by name (`RingNaming`).
pub(crate) fn comp_face(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    lf: &LocalFace,
) -> Result<combinatorics::CompFace, BoolError> {
    let circle = |cyl: usize| combinatorics::BoundEdges::Circle(Box::new(cyls[cyl].def.clone()));
    match lf.surf {
        ClassIx::Plane(c) => {
            let bound = |b: &Bound| -> Result<combinatorics::BoundEdges, BoolError> {
                Ok(match b {
                    Bound::Ring(r) => {
                        combinatorics::BoundEdges::Ring(r.edges(jd, cyls, Some(c))?)
                    }
                    Bound::Circle { cyl } => circle(*cyl),
                    // Rims bound a cylinder, never a plane — a producer error the ray
                    // abstains on by name.
                    Bound::Band { .. } => combinatorics::BoundEdges::Lateral(Vec::new()),
                })
            };
            Ok(combinatorics::CompFace {
                surf: combinatorics::CompSurf::Plane(c),
                outer: bound(&lf.outer)?,
                inner: lf
                    .inner
                    .iter()
                    .map(bound)
                    .collect::<Result<Vec<_>, BoolError>>()?,
            })
        }
        ClassIx::Cyl(k) => {
            let ring_loop = |r: &Ring| -> Result<combinatorics::LateralLoop, BoolError> {
                Ok(combinatorics::LateralLoop::Ring(r.edges(jd, cyls, None)?))
            };
            let rim_loop = |rim: &Rim| -> Result<combinatorics::LateralLoop, BoolError> {
                Ok(match rim {
                    Rim::Circle(c) => combinatorics::LateralLoop::Circle(*c),
                    Rim::Chain(r) => ring_loop(r)?,
                })
            };
            let mut loops = Vec::new();
            for b in core::iter::once(&lf.outer).chain(lf.inner.iter()) {
                match b {
                    Bound::Ring(r) => loops.push(ring_loop(r)?),
                    Bound::Band { lo, hi } => {
                        loops.push(rim_loop(lo)?);
                        loops.push(rim_loop(hi)?);
                    }
                    // A circle bound on a cylinder face has no producer; the ray abstains on
                    // it by name rather than guess.
                    Bound::Circle { cyl } => {
                        return Ok(combinatorics::CompFace {
                            surf: combinatorics::CompSurf::Cylinder(Box::new(cyls[k].def.clone())),
                            outer: circle(*cyl),
                            inner: Vec::new(),
                        });
                    }
                }
            }
            Ok(combinatorics::CompFace {
                // ★ The truth travels, not the index — the probe should not have to hold the
                // class table to solve a crossing.
                surf: combinatorics::CompSurf::Cylinder(Box::new(cyls[k].def.clone())),
                outer: combinatorics::BoundEdges::Lateral(loops),
                inner: Vec::new(),
            })
        }
    }
}

/// Partition the faces into connected components, then into output solids. One component is the
/// whole result; several mean either an enclosed void (a cavity — an inward-oriented shell) or a
/// severed operand (two or more material-enclosing shells).
///
/// ★ **Deciding early is not the same as *rejecting* early.** [`reconstruct`] holds this `Result`
/// and raises it late, past the minting; the arena a declining boolean leaves is the quantity
/// `replay::a_late_reject_is_not_index_neutral` measures.
pub(super) fn group_faces(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::draft::CutRims,
) -> Result<Grouping, BoolError> {
    let (labels, n) = face_components(faces, cut_rims);
    let mut by_comp_lf: Vec<Vec<&LocalFace>> = vec![Vec::new(); n];
    for (i, lf) in faces.iter().enumerate() {
        by_comp_lf[labels[i]].push(lf);
    }
    // A component as `(plane, rings)` per face, and its nodes — what the containment predicate
    // takes. Both the material/void label below and the cavity-owner search further down ask the
    // same question of them.
    // ★ **Only when there is something to compare.** One component is the whole result — depth 0,
    // material, nothing to ask — and that is nearly every boolean. Building these rings costs an
    // exact predicate per node, so doing it unconditionally put **12x** on the 80-fold star
    // (measured). Built once per component rather than once per query for the same reason: the
    // label below retries with another node when one grazes.
    let comp_faces: Vec<combinatorics::ComponentFaces> = (0..if n > 1 { n } else { 0 })
        .map(|c| {
            by_comp_lf[c]
                .iter()
                .map(|lf| comp_face(jd, cyls, lf))
                .collect::<Result<_, BoolError>>()
        })
        .collect::<Result<_, BoolError>>()?;
    // A **probe list**, so a pierce node may be dropped: `first_deciding` below tries each in
    // turn, and an exhausted list is already `Ok(None)` — the caller's own rejection.
    //
    // ★★★ **Coordinates only when there is no vertex to name.** A named probe is exact without
    // coordinates at all, so it keeps working where no rational one exists (a rotated class), and
    // it is the answer for every component that has a polygon anywhere on it. What has none is a
    // lone cylinder, whose boundary is two disks and a band and whose vertex count is therefore
    // **zero**. `coord_probes` names a
    // point on that boundary instead of a vertex of it.
    let probes_of = |c: usize| -> Vec<combinatorics::Probe> {
        // ★ Each name once, in first-seen order: a vertex sits on three faces' rings and used
        // to be offered three times, so an exhausted list cast the same rays three times over
        // (the corner Commons: six probes for two points). The first deciding probe is the same,
        // so the answers are; only the ledger's «offered» count falls.
        let mut seen = std::collections::HashSet::new();
        let named: Vec<combinatorics::Probe> = combinatorics::three_plane_probes(
            by_comp_lf[c]
                .iter()
                .flat_map(|lf| lf.poly_rings().flat_map(|r| r.iter().copied())),
        )
        .into_iter()
        .filter(|x| seen.insert(*x))
        .map(combinatorics::Probe::Named)
        .collect();
        if named.is_empty() {
            // ★ Rational pierce corners first (a slot or a fillet-cornered prism has
            // no other vertex), then the caps' derived points.
            let mut out = combinatorics::corner_probes(
                jd,
                cyls,
                by_comp_lf[c]
                    .iter()
                    .flat_map(|lf| lf.poly_rings().flat_map(|r| r.iter().copied())),
            );
            out.extend(combinatorics::coord_probes(jd, cyls, &comp_faces[c]));
            out
        } else {
            named
        }
    };
    // ★★★★★ **The supply a component's *edges* name, tried only when its corners are exhausted**.
    // `probes_of` above offers a polyhedral component nothing but its **vertices** — the
    // two fallbacks beside it both need a cylinder (`corner_probes` wants pierce corners,
    // `coord_probes` a cap's circle), so a planar component has none. A void whose every corner
    // sits on its host's wall would therefore run out of witnesses and be refused `NoClearRay`,
    // while the *shape's* truth — its surface meets itself along the corner's edge — is what the
    // self-touch test says once the depth is decided. Its one-grazing-corner sibling has always
    // said exactly that, with the same witness segment; the two differed only in **how many**
    // corners were degenerate, which is a fact about the supply and not about the geometry.
    //
    // Every point comes from [`combinatorics::edge_interior_points`] — the same rule the 2-D
    // witness supply reads, not a second spelling of it — so each is **on** this component's
    // boundary, which is [`combinatorics::Probe`]'s invariant.
    //
    // ★ **Built only on exhaustion**, because `probes_of` is eager and asked at three sites: a
    // component with 8 vertices has some 30 edge points behind it, and paying for them where a
    // vertex answers first is waste the ledger already warns about two doc comments above.
    let edge_probes_of = |c: usize| -> Vec<combinatorics::Probe> {
        let mut out = Vec::new();
        for f in comp_faces[c].iter() {
            let normal = match &f.surf {
                combinatorics::CompSurf::Plane(q) => {
                    match combinatorics::class_coeffs_rat(jd, *q) {
                        Some(k) => [k[0], k[1], k[2]],
                        None => continue,
                    }
                }
                // A curved face's own edges are still edges; its class has no plane to point along,
                // so the ray directions come from the cylinder the `Carrier` names — which is what
                // `corner_probes` already does for that population. Nothing here yet.
                _ => continue,
            };
            // One direction set per face: it is a function of the normal alone, and rebuilding it
            // per point would allocate once for every edge of every ring.
            let dirs = combinatorics::probe_dirs(&normal);
            let rings = std::iter::once(&f.outer).chain(f.inner.iter());
            for b in rings {
                let combinatorics::BoundEdges::Ring(r) = b else {
                    continue;
                };
                for e in r {
                    for p in combinatorics::edge_interior_points(jd, cyls, e) {
                        out.extend(
                            dirs.iter()
                                .map(|&dir| combinatorics::Probe::Coord { p, dir }),
                        );
                    }
                }
            }
        }
        out
    };
    // Try `f` at each node in turn: the first node that **decides** wins, a node that abstains
    // (`Ok(None)` — it grazed a boundary) is passed over for the next, and a failed judgement
    // (`Err`) propagates immediately. The last part is the point of the shape: an abstention has
    // other nodes as its remedy, a failed judgement does not — retrying on `Err(_)` too would let
    // a real cause masquerade as `NoClearRay` once every node hit it. `Ok(None)` here means every
    // node abstained; that being a reject is the *caller's*
    // proposition to raise.
    /// ★★ **Two stages, the second built only on exhaustion**: `more` is the supply that
    /// costs something to derive, and a component whose first corner decides never pays for it.
    /// Trying it **after** the primary list means every question the primary list decides is
    /// decided on the same probe, in the same order; the second stage answers only questions the
    /// first ran out on.
    ///
    /// The ledger's `offered` is now **what was actually built**: the primary list when it decided,
    /// the sum when it did not. One row per question either way.
    fn first_deciding<T>(
        probes: &[combinatorics::Probe],
        more: impl FnOnce() -> Vec<combinatorics::Probe>,
        mut f: impl FnMut(&combinatorics::Probe) -> Result<Option<T>, BoolError>,
    ) -> Result<Option<T>, BoolError> {
        #[cfg(test)]
        let tie0 = combinatorics::tie_probe::len();
        for (tried, x) in probes.iter().enumerate() {
            #[cfg(not(test))]
            let _ = tried;
            if let Some(v) = f(x)? {
                #[cfg(test)]
                probe::deciding::record(tried + 1, probes.len(), true, Vec::new());
                return Ok(Some(v));
            }
        }
        let extra = more();
        let offered = probes.len() + extra.len();
        for (tried, x) in extra.iter().enumerate() {
            #[cfg(not(test))]
            let _ = tried;
            if let Some(v) = f(x)? {
                #[cfg(test)]
                probe::deciding::record(probes.len() + tried + 1, offered, true, Vec::new());
                return Ok(Some(v));
            }
        }
        #[cfg(not(test))]
        let _ = offered;
        #[cfg(test)]
        probe::deciding::record(
            offered,
            offered,
            false,
            combinatorics::tie_probe::since(tie0),
        );
        Ok(None)
    }
    // ★ A list that **started empty** is not a list that ran out: a component whose
    // every corner is a pierce name and whose caps offer no candidate their ring says is inside
    // — a thin segment of a disk — has no witness at all, and says so by the name the rings use
    // for the same proposition.
    ///
    /// ⚠ **Asked of both stages**: a component whose corners are all pierce names *and*
    /// whose edges name no interior point has started empty; one whose edges did offer points has
    /// run out.
    fn no_witness(primary: &[combinatorics::Probe], extra: usize) -> RejectReason {
        if primary.is_empty() && extra == 0 {
            RejectReason::RingHasNoWitness
        } else {
            RejectReason::NoClearRay
        }
    }
    // ★★ **Material or void is a question about nesting, not about normals.**
    //
    // Not at the component's lexicographically-minimal vertex `v*` ("outward iff *some* face there
    // has an outward normal with `n_x < 0`"): that existential is a shortcut for "the face with the
    // largest `|n_x|` faces −x", equivalent only while the normals are **axis-aligned**. A slanted
    // sketch breaks it with no rotation in sight: a wedge void cut inside a box comes back as its
    // own *material* solid of **negative volume**, with the box unchanged beside it — a wrong
    // model `validate` has nothing to say about (measured).
    //
    // So ask the question the kernel already answers one dimension down. `sketch::from_rings`
    // decides a ring by **containment depth — even is material, odd is a hole** (design.md: there
    // is no
    // fill-rule parameter — even-odd is the only rule), and the cavity-owner search
    // below already picks the *innermost* container, the other half of that same rule. This is the
    // 3D reading of it, with `point_in_component` where the 2D one uses `point_in_ring`.
    //
    // Nothing here reads an orientation, so a shell that winds either way is labelled the same —
    // which is the point: the label does not try to recover the winding from orientations.
    // `positives` stays in ascending `c` order (a downstream contract).
    let mut positives: Vec<usize> = Vec::new();
    for c in 0..n {
        if n == 1 {
            positives.push(0); // the whole result: depth 0, and no other component to be inside
            break;
        }
        // One origin for the whole row: a node of `c` that classifies against *every* other
        // component. Trying them in turn is what the cavity search does, and for the same reason —
        // a node that grazes one component's boundary is a fact about that node, not about the
        // components.
        let probes = probes_of(c);
        let mut extra_len = 0usize;
        let depth = first_deciding(
            &probes,
            || {
                let v = edge_probes_of(c);
                extra_len = v.len();
                v
            },
            |x| {
                let mut d = 0usize;
                for other in (0..n).filter(|&o| o != c) {
                    match combinatorics::probe_in_component(jd, cyls, x, &comp_faces[other])? {
                        Some(true) => d += 1,
                        Some(false) => {}
                        None => return Ok(None), // grazed — this node abstains
                    }
                }
                Ok(Some(d))
            },
        )?
        .ok_or_else(|| reject(no_witness(&probes, extra_len)))?;
        if depth % 2 == 0 {
            positives.push(c);
        }
    }
    // A result solid is a material component and the cavity components it owns, and that grouping —
    // not the component — is the unit a self-contact question is asked about, and the unit handles
    // are minted per.
    let mut comps_of: HashMap<usize, Vec<usize>> =
        positives.iter().map(|&m| (m, vec![m])).collect();
    match positives.len() {
        0 => return Err(reject(RejectReason::NoOutwardShell)),
        // One material: every other component is a cavity of it, and no search is needed to say so.
        1 => {
            comps_of.insert(positives[0], (0..n).collect());
        }
        _ => {
            // Several material solids. Assign each surviving cavity (an inward component) to the
            // material whose outer shell nests it — the innermost, if materials themselves nest —
            // by the exact point-in-solid `point_in_component`. With no cavity this loop is empty
            // and every material emits cavity-free (the plain sever, unchanged).
            // ★ The rings hand over their walls; nothing here derives one from a name —
            // `comp_faces`/`probes_of` are the same pair the material/void label above uses.
            for d in (0..n).filter(|c| !positives.contains(c)) {
                // A cavity node that classifies cleanly against *every* material (one shared origin
                // keeps the nesting consistent); its `true` materials nest, so take the innermost.
                let probes = probes_of(d);
                let mut extra_len = 0usize;
                let containers = first_deciding(
                    &probes,
                    || {
                        let v = edge_probes_of(d);
                        extra_len = v.len();
                        v
                    },
                    |x| {
                        let mut cs = Vec::new();
                        for &m in &positives {
                            match combinatorics::probe_in_component(jd, cyls, x, &comp_faces[m])? {
                                Some(true) => cs.push(m),
                                Some(false) => {}
                                None => return Ok(None), // grazed against a material — abstain
                            }
                        }
                        Ok(Some(cs))
                    },
                )?
                .ok_or_else(|| reject(no_witness(&probes, extra_len)))?;
                let owner = match containers.as_slice() {
                    [] => return Err(reject(RejectReason::CavityNoOwner)),
                    [only] => *only,
                    _ => {
                        // Innermost: inside every other container. The pairwise facts are
                        // decided **before** the search, because the search closures speak
                        // bool and a failed judgement must propagate, not vanish into "not
                        // inside". A pair where every node abstains stays `false` — no
                        // evidence places `c` inside `o`. (Nested containers are a rare
                        // population; deciding
                        // the few pairs up front costs nothing measurable.)
                        let mut inside = HashMap::new();
                        for &c in &containers {
                            for &o in &containers {
                                if o == c {
                                    continue;
                                }
                                let v = first_deciding(
                                    &probes_of(c),
                                    || edge_probes_of(c),
                                    |x| {
                                        combinatorics::probe_in_component(
                                            jd,
                                            cyls,
                                            x,
                                            &comp_faces[o],
                                        )
                                    },
                                )?
                                .unwrap_or(false);
                                inside.insert((c, o), v);
                            }
                        }
                        *containers
                            .iter()
                            .find(|&&c| containers.iter().all(|&o| o == c || inside[&(c, o)]))
                            .ok_or_else(|| reject(RejectReason::CavityNoOwner))?
                    }
                };
                comps_of.get_mut(&owner).unwrap().push(d);
            }
        }
    }
    let mut comp_group = vec![0usize; n];
    for (gi, m) in positives.iter().enumerate() {
        for &c in &comps_of[m] {
            comp_group[c] = gi;
        }
    }
    let group_of = labels.iter().map(|&l| comp_group[l]).collect();
    Ok(Grouping {
        labels,
        n,
        positives,
        comps_of,
        group_of,
    })
}

/// A canonical, replay-stable sort key for a severed component: its outer-loop vertex
/// coordinates, sorted lexicographically. Total order for disjoint components — distinct pieces
/// occupy different space, so their coordinate multisets differ, and the lex-min vertex alone can
/// tie (identity is by `Handle`, not coordinates, so two vertices may coincide). Coordinates are a
/// derived cache used only to order multi-solid output; no judgment reads this (tol-irrelevant).
pub(super) fn comp_key(model: &Model, faces: &[Handle<Face>]) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = faces
        .iter()
        .flat_map(|&fh| model.face(fh).outer.half_edges.iter().copied())
        .map(|he| model.vertex_point(he_start(model, he)).as_array())
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite vertex coordinates"));
    pts
}
