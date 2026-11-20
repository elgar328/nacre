use super::*;
pub(crate) struct Named {
    /// The grouping — always one that answered: a grouping that cannot answer is the naming's
    /// refusal. (`pub(crate)`: the grouping fence reads it directly.)
    pub(crate) grouping: Grouping,
    pub(crate) group_of: Vec<usize>,
    /// The working face list where it differs from the caller's — the split-twin subdivision
    /// and/or the per-solid dissolve rewrote rings. `None` = the caller's list stands.
    pub(crate) per_solid: Option<Vec<LocalFace>>,
    pub(crate) defs: HashMap<(usize, NodeId), Def>,
}

/// **Everything `reconstruct` decides before a single handle is minted** — the grouping, the
/// per-solid straight-angle dissolve, the whole-result self-touch judgement, and every ring node's
/// defining triple.
///
/// A named function rather than the top of `reconstruct`, for the same reason `seam_table` is
/// one: a fence calls it directly on the faces production feeds it. Model-immutable by signature:
/// nothing here takes `&mut Model`.
pub(crate) fn name_result_vertices(
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    tg: &super::Tangencies<'_>,
) -> Result<Named, BoolError> {
    // **The split-twin subdivision — every pierce node is cut into every edge it lies on.**
    // The arrangement cannot do this: a pierce point needs the cylinder, and the neighbouring
    // plane class has no circle (the cylinder is parallel to it), so no vocabulary for the point
    // — measured, the T-junction finding. Only here, where every class's rings are in one hand,
    // can the plate-top's whole edge learn that the arc class subdivided its twin. Without this,
    // `norm_edge` keys never match across such a pair, the corner's incident-face set misses the
    // twin's plane (a result vertex is defined from the planes of the faces whose rings
    // visit it), and the future edge welding has no twins to weld.
    //
    // An edge's candidates are the pierce nodes **on its line** — on both its planes, which
    // `combinatorics::pierce_on_line` asks of the planes rather than of the node's name: the name a
    // ring carries is one pair's, picked by the class order, and a point on two lines has a pair per
    // line (a tangent point a notch's plane also passes through; a prism's edge along a line two
    // classes share on a cylinder, past the lateral's cap, the keyhole's box taller than the bore,
    // whose rim point is named by the cap and one of the two). Only the pierce nodes whose proved
    // box meets the edge's are asked — a point on the segment is in the hull of its ends — since
    // the exact question is the cost (measured over the census: every node against every line
    // its faces ride put the serial census +1.7% in instructions). Betweenness is
    // `between_on_line`'s.
    // An edge whose order cannot be formed is left unsplit (a corner short of a plane,
    // `StraightAngle`): the conservative arm degrades to the unsplit edge, never to something new.
    // Identity for every non-arc input: no pierce nodes, no candidates, no rewrite.
    let reach = node_reach(seam);
    let mut pierced: Vec<NodeId> = Vec::new();
    for lf in faces {
        for ring in lf.poly_rings() {
            for &n in ring.iter() {
                if combinatorics::pierce_name(n).is_some() && !pierced.contains(&n) {
                    pierced.push(n);
                }
            }
        }
    }
    let subdivided: Option<Vec<LocalFace>> = if pierced.is_empty() {
        None
    } else {
        // Per line and node, whether it holds the node; per node, its solved point — one point is
        // asked against every line its faces' edges ride.
        let mut on_line: HashMap<([usize; 2], NodeId), bool> = HashMap::new();
        let mut points: HashMap<NodeId, Option<Box<combinatorics::PierceAt>>> = HashMap::new();
        let mut v = faces.to_vec();
        for lf in &mut v {
            let ClassIx::Plane(own) = lf.surf else {
                continue;
            };
            for ring in lf.poly_rings_mut() {
                let k = ring.nodes.len();
                let mut nodes = Vec::with_capacity(k);
                let mut walls = Vec::with_capacity(k);
                for t in 0..k {
                    let (a, b, w) = (ring.nodes[t], ring.nodes[(t + 1) % k], ring.walls[t]);
                    nodes.push(a);
                    walls.push(w);
                    let Wall::Plane(w) = w else {
                        // An arc or ruling edge: the arrangement's own split made it, pierce
                        // points and all — there is no whole twin to subdivide.
                        continue;
                    };
                    let mut pair = [own, w];
                    pair.sort_unstable();
                    let (el, eh) = reach.edge_reach([a, b]);
                    let cands: Vec<NodeId> = pierced
                        .iter()
                        .copied()
                        .filter(|&n| {
                            let (nl, nh) = reach.reach(n);
                            (0..3).all(|i| nl[i] <= eh[i] && nh[i] >= el[i])
                                && *on_line.entry((pair, n)).or_insert_with(|| {
                                    combinatorics::pierce_on_line(jd, cyls, &mut points, n, own, w)
                                        == Some(true)
                                })
                        })
                        .collect();
                    if cands.is_empty() {
                        continue;
                    }
                    if let Some(bet) =
                        combinatorics::between_on_line(jd, cyls, own, w, a, b, &cands)
                    {
                        for x in bet {
                            nodes.push(x);
                            walls.push(Wall::Plane(w));
                        }
                    }
                }
                ring.nodes = nodes;
                ring.walls = walls;
            }
        }
        Some(v)
    };
    let faces: &[LocalFace] = subdivided.as_deref().unwrap_or(faces);
    // **Which faces make one solid, decided before a single handle exists** — see [`Grouping`].
    // Everything derived below (an edge's far plane, a vertex's defining triple, the handles
    // themselves) is scoped to one group, so no result solid can be named by — or share a handle
    // with — a solid it merely touches.
    //
    // **A grouping that cannot answer is the answer**, raised here: every step below reads the
    // groups, so none of them may speak over a stand-in. Putting every face in one group instead
    // welds a contact line's four uses onto one handle, and the edge-use guard then names
    // `NonManifoldResultEdge` for a result whose faces the component join could not pair.
    let grouping = group_faces(jd, faces, cyls)?;
    let group_of: Vec<usize> = grouping.group_of.clone();
    // **A straight angle is a per-solid question, and only a result in several pieces can make
    // the two answers differ.** The cleaning pass runs `dissolve_straight_angles` over the whole
    // result, so it keeps a node that is a real corner *somewhere* — right while the result is one
    // body. The moment it is two, the tip of one body's knife edge can land in the middle of
    // another body's wall: a corner of the first, a straight run of the second. Left in the
    // second's ring it is a vertex with no name there — only two of the three planes through it
    // bound that solid — and a whole-result derivation fills the gap by borrowing the *other*
    // body's plane. So each solid drops the nodes that are straight runs **for it**.
    let per_solid: Option<Vec<LocalFace>> = if grouping.n > 1 {
        let mut v = faces.to_vec();
        let shared = combinatorics::SharedLines::of(jd, cyls);
        for gi in 0..grouping.positives.len() {
            let which: Vec<usize> = (0..v.len()).filter(|&i| group_of[i] == gi).collect();
            dissolve_straight_angles(&mut v, &which, &shared);
        }
        Some(v)
    } else {
        None
    };
    let faces: &[LocalFace] = per_solid.as_deref().unwrap_or(faces);
    // **The whole-result judgement runs before a single cell is minted.** Everything
    // `self_touch_reject` reads — the seam table, the per-body component lists, the faces'
    // rings — exists right here, and an impossible result must be named by its truth, not by
    // whichever local derivation happens to fail first on its unnameable corners: a
    // self-touching body's pinch line cannot be honestly named by any face-local rule (its
    // in-plane edges are lobe-to-lobe, its touch edge carries four faces — measured).
    {
        let g = &grouping;
        let body_comps: Vec<Vec<usize>> =
            g.positives.iter().map(|m| g.comps_of[m].clone()).collect();
        let mut by_comp_lf: Vec<Vec<&LocalFace>> = vec![Vec::new(); g.n];
        for (i, lf) in faces.iter().enumerate() {
            by_comp_lf[g.labels[i]].push(lf);
        }
        self_touch_reject(jd, cyls, &reach, &body_comps, &by_comp_lf, tg)?;
    }
    face_bounds_meet(jd, cyls, &reach, faces)?;
    // **A result vertex is named by the faces that meet it.**
    //
    // The arrangement names a point by a canonical plane triple — lexicographically first among
    // the planes through it — and where four planes concur that choice can land on a plane the
    // result keeps **no face on**: a rotated copy's wall, say, buried inside the union it was
    // fused into. The definition *is* the vertex's identity, so such a result cannot describe
    // itself: `transform`'s remap (which walks this solid's own face surfaces) would refuse it two
    // operations later, and the assembly's closing fence refuses it here.
    //
    // **A result vertex is defined by the canonical triple of its incident faces'
    // planes.** Every face whose ring visits the node passes through the point, so the planes of
    // those faces are the result planes through it, and `canonical_triple` picks the name — the
    // one rule the operand road and the alias table use. Over *incident* planes only, on purpose:
    // the arrangement's alias representative ranges over every class, buried faces included, and
    // a definition naming a surface no face at the vertex lies on is a defect (the assembly's
    // closing fence, `VertexNamesAbsentSurface`, stands guard). A node fewer than three faces
    // name is a straight corner and keeps `StraightAngle` below.
    //
    // **A pierce vertex is named by them too — its cylinder included.** Its name carries two
    // plane classes and the cylinder, and while those are surfaces of faces the vertex's rings
    // put there — both planes, and a lateral of that cylinder — the name is its definition as it
    // stands. Where a line two classes share lies on the cylinder, one point has a pierce name
    // per pair (`Aliases::record_shared_ruling`), and the representative can name a plane no face
    // meets there. The cylinder can be that surface too: a corner of a solid lies on the class's
    // unbounded surface wherever the arrangement put it, and the lateral need not reach it — a
    // block standing on another solid's rim, its corner a straight run of that rim and the two
    // solids separate; a box whose edge runs along a ruling of the cylinder and on past its cap,
    // the edge's top corner on the cylinder's surface above the cap. So a pierce vertex keeps its
    // name when a lateral of its cylinder and both its planes meet it here; otherwise it takes the
    // canonical triple of the planes that do — three of them only when they meet in a point, not
    // when they hold one line ([`combinatorics::canonical_triple`]) — else the pierce name of the
    // first pair of them that restates it ([`combinatorics::restate_pierce`]). Where no pair does —
    // fewer than two planes, the pair meets the cylinder elsewhere, or the arithmetic cannot say,
    // which `restate_pierce`'s `None` does not tell apart — the name stands. A restatement stays on
    // the same cylinder, so with no lateral of it here and fewer than three planes there is no name
    // the vertex's own faces give: the fence refuses the vertex by that symptom, not the cause.
    let mut planes_at: HashMap<(usize, NodeId), Vec<usize>> = HashMap::new();
    let mut def_triple: HashMap<(usize, NodeId), Def> = HashMap::new();
    // `(group, node, cylinder class)` for every node a lateral's ring visits — the cylinders that
    // meet the vertex, as `planes_at` holds the planes.
    let mut cyls_at: HashSet<(usize, NodeId, usize)> = HashSet::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        for ring in lf.poly_rings() {
            for &node in &ring.nodes {
                match lf.surf {
                    ClassIx::Plane(c) => {
                        let at = planes_at.entry((g, node)).or_default();
                        if !at.contains(&c) {
                            at.push(c);
                        }
                    }
                    ClassIx::Cyl(k) => {
                        planes_at.entry((g, node)).or_default();
                        cyls_at.insert((g, node, k));
                    }
                }
            }
        }
    }
    let mut keys: Vec<(usize, NodeId)> = planes_at.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        let mut at = planes_at.remove(&key).unwrap_or_default();
        at.sort_unstable();
        let Some((own, cyl, _)) = combinatorics::pierce_name(key.1) else {
            if let Some(t) = combinatorics::canonical_triple(jd, &at) {
                def_triple.insert(key, Def::Three(t));
            }
            continue;
        };
        let def = if cyls_at.contains(&(key.0, key.1, cyl)) && own.iter().all(|p| at.contains(p)) {
            Def::of_pierce_name(key.1)
        } else if let Some(t) = combinatorics::canonical_triple(jd, &at) {
            Some(Def::Three(t))
        } else {
            at.iter()
                .enumerate()
                .flat_map(|(i, &a)| at[i + 1..].iter().map(move |&b| (a, b)))
                .find_map(|(a, b)| combinatorics::restate_pierce(jd, cyls, key.1, a, b))
                .and_then(Def::of_pierce_name)
                .or_else(|| Def::of_pierce_name(key.1))
        };
        if let Some(def) = def {
            def_triple.insert(key, def);
        }
    }

    Ok(Named {
        grouping,
        group_of,
        // The dissolve's product already derives from the subdivided list (it cloned `faces`
        // after the rebinding above), so the later layer wins and the earlier one backs it up.
        per_solid: per_solid.or(subdivided),
        defs: def_triple,
    })
}

/// [`NodeReach`] off the seam table.
fn node_reach(seam: &[SeamVertex]) -> NodeReach {
    NodeReach(
        seam.iter()
            .map(|sv| {
                (
                    sv.triple,
                    (sv.cache.coord(), crate::planes::cache_reach(&sv.cache)),
                )
            })
            .collect(),
    )
}
