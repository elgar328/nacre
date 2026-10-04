use super::*;
pub(crate) struct Named {
    /// The grouping, **held** — raised deep in the minting.
    /// (`pub(crate)`: the grouping fence reads the held result directly.)
    pub(crate) grouping: Result<Grouping, BoolError>,
    pub(crate) group_of: Vec<usize>,
    /// The working face list where it differs from the caller's — the split-twin subdivision
    /// and/or the per-solid dissolve rewrote rings. `None` = the caller's list stands.
    pub(crate) per_solid: Option<Vec<LocalFace>>,
    pub(crate) defs: HashMap<(usize, NodeId), Def>,
}

/// **Everything `reconstruct` decides before a single handle is minted** — the grouping (held,
/// not raised), the per-solid straight-angle dissolve, the whole-result self-touch judgement, and
/// every ring node's defining triple.
///
/// ★ A named function rather than the top of `reconstruct`, for the same reason `seam_table` is
/// one: the deferred-stopper socket stands behind it (at the assembly's very end) and
/// intercepts everything a plugged stopper would, so no
/// reject name can testify that the naming completed — only a fence that calls it directly on the
/// faces production feeds it can. Model-immutable by signature: nothing here takes `&mut Model`.
pub(crate) fn name_result_vertices(
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Named, BoolError> {
    // ★★★ **The split-twin subdivision — every pierce node is cut into every edge it lies on.**
    // The arrangement cannot do this: a pierce point needs the cylinder, and the neighbouring
    // plane class has no circle (the cylinder is parallel to it), so no vocabulary for the point
    // — measured, the T-junction finding. Only here, where every class's rings are in one hand,
    // can the plate-top's whole edge learn that the arc class subdivided its twin. Without this,
    // `norm_edge` keys never match across such a pair, the corner's incident-face set misses the
    // twin's plane (a result vertex is defined from the planes of the faces whose rings
    // visit it), and the future edge welding has no twins to weld.
    //
    // A pierce node lies on an edge's carrier line exactly when its plane pair *is* the edge's
    // `{own, wall}` — a name fact, no geometry — and betweenness is `pierce_between`'s exact
    // half. An edge whose order cannot be formed is left unsplit, which is precisely the
    // behaviour before this pass (a corner short of a plane, `StraightAngle`); the conservative
    // arm degrades to the state this pass improves, never to something new. Identity for every non-arc input:
    // no pierce nodes, no pairs, no rewrite.
    let mut by_pair: HashMap<[usize; 2], Vec<NodeId>> = HashMap::new();
    for lf in faces {
        for ring in lf.poly_rings() {
            for &n in ring.iter() {
                if let Some((pair, _, _)) = combinatorics::pierce_name(n) {
                    let v = by_pair.entry(pair).or_default();
                    if !v.contains(&n) {
                        v.push(n);
                    }
                }
            }
        }
    }
    // ★ **And a line two classes share on a cylinder** (`SharedRuling`): its pierce points are
    // named by a cap and *one* of the two classes, never by the pair itself, so the name fact
    // above cannot find them there. A prism's edge along that line runs past the lateral's cap
    // (the keyhole's box taller than the bore), and the lateral's ruling ends at the cap's rim
    // point while the prism's wall — tangent, so no class split it — carries the edge whole. The
    // point's being on the other class is asked of the point (`side_of`), once per line.
    let pierced: Vec<NodeId> = by_pair.values().flatten().copied().collect();
    for cy in cyls {
        for sr in &cy.shared {
            let (a, b) = sr.line;
            let on_line: Vec<NodeId> = pierced
                .iter()
                .copied()
                .filter(|&n| {
                    let Some((pl, _, _)) = combinatorics::pierce_name(n) else {
                        return false;
                    };
                    (pl.contains(&a) && combinatorics::side_of(jd, cyls, n, b) == Some(0))
                        || (pl.contains(&b) && combinatorics::side_of(jd, cyls, n, a) == Some(0))
                })
                .collect();
            if on_line.is_empty() {
                continue;
            }
            let v = by_pair.entry([a, b]).or_default();
            for n in on_line {
                if !v.contains(&n) {
                    v.push(n);
                }
            }
        }
    }
    let subdivided: Option<Vec<LocalFace>> = if by_pair.is_empty() {
        None
    } else {
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
                    let Some(cands) = by_pair.get(&pair) else {
                        continue;
                    };
                    if let Some(bet) = combinatorics::pierce_between(jd, cyls, a, b, cands) {
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
    // ★ **Which faces make one solid, decided before a single handle exists** — see [`Grouping`].
    // Everything derived below (an edge's far plane, a vertex's defining triple, the handles
    // themselves) is scoped to one group, so no result solid can be named by — or share a handle
    // with — a solid it merely touches.
    //
    // ★★ **Held, not raised.** A failed grouping is reported further down, and until then every
    // face is one group (`replay::a_late_reject_is_not_index_neutral` measures the arena a
    // declining boolean leaves).
    let grouping = group_faces(jd, faces, cyls);
    let group_of: Vec<usize> = match &grouping {
        Ok(g) => g.group_of.clone(),
        Err(_) => vec![0; faces.len()],
    };
    // ★★ **A straight angle is a per-solid question, and only a result in several pieces can make
    // the two answers differ.** The cleaning pass runs `dissolve_straight_angles` over the whole
    // result, so it keeps a node that is a real corner *somewhere* — right while the result is one
    // body. The moment it is two, the tip of one body's knife edge can land in the middle of
    // another body's wall: a corner of the first, a straight run of the second. Left in the
    // second's ring it is a vertex with no name there — only two of the three planes through it
    // bound that solid — and a whole-result derivation fills the gap by borrowing the *other*
    // body's plane. So each solid drops the nodes that are straight runs **for it**.
    let per_solid: Option<Vec<LocalFace>> = match &grouping {
        Ok(g) if g.n > 1 => {
            let mut v = faces.to_vec();
            for gi in 0..g.positives.len() {
                let which: Vec<usize> = (0..v.len()).filter(|&i| group_of[i] == gi).collect();
                dissolve_straight_angles(&mut v, &which);
            }
            Some(v)
        }
        _ => None,
    };
    let faces: &[LocalFace] = per_solid.as_deref().unwrap_or(faces);
    // ★ **The whole-result judgement runs before a single cell is minted.** Everything
    // `self_touch_reject` reads — the seam table, the per-body component lists, the faces'
    // rings — exists right here, and an impossible result must be named by its truth, not by
    // whichever local derivation happens to fail first on its unnameable corners: a
    // self-touching body's pinch line cannot be honestly named by any face-local rule (its
    // in-plane edges are lobe-to-lobe, its touch edge carries four faces — measured).
    // A held grouping error stays held (raised further down); the self-touch question is only
    // askable of a
    // grouping that answered.
    if let Ok(g) = &grouping {
        let body_comps: Vec<Vec<usize>> =
            g.positives.iter().map(|m| g.comps_of[m].clone()).collect();
        let mut by_comp_lf: Vec<Vec<&LocalFace>> = vec![Vec::new(); g.n];
        for (i, lf) in faces.iter().enumerate() {
            by_comp_lf[g.labels[i]].push(lf);
        }
        self_touch_reject(jd, cyls, seam, &body_comps, &by_comp_lf)?;
    }
    // ★★ **A result vertex is named by the faces that meet it.**
    //
    // The arrangement names a point by a canonical plane triple — lexicographically first among
    // the planes through it — and where four planes concur that choice can land on a plane the
    // result keeps **no face on**: a rotated copy's wall, say, buried inside the union it was
    // fused into. The definition *is* the vertex's identity, so such a result cannot describe
    // itself, and `transform`'s remap (which walks this solid's own face surfaces) refuses it two
    // operations later — a reject whose cause is here.
    //
    // ★ **A result vertex is defined by the canonical triple of its incident faces'
    // planes.** Every face whose ring visits the node passes through the point, so the planes of
    // those faces are the result planes through it, and `canonical_triple` picks the name — the
    // one rule the operand road and the alias table use. Over *incident* planes only, on purpose:
    // the arrangement's alias representative ranges over every class, buried faces included, and
    // a definition naming a surface the result has no face on is a defect
    // (`defs_are_remappable` still stands guard, true by construction). A node fewer than
    // three faces name is a straight corner and keeps `StraightAngle` below.
    //
    // ★ **A pierce vertex is named by them too.** Its name carries two plane classes and the
    // cylinder, and while those are the planes the result's faces meet it on, the name is its
    // definition as it stands. Where a line two classes share lies on the cylinder, one point
    // has a pierce name per pair (`Aliases::record_shared_ruling`), and the representative can
    // name a plane the result keeps no face on — the vertex then says where it is by a surface
    // the solid does not have. So a pierce vertex keeps its name when both its planes meet it
    // here, or when fewer than two planes do (nothing to restate it on); otherwise it takes the
    // canonical triple of the planes that do, else the pierce name of the first pair of them that
    // restates it ([`combinatorics::restate_pierce`]). Where no pair does — the pair meets the
    // cylinder elsewhere, or the arithmetic cannot say, which `restate_pierce`'s `None` does not
    // tell apart — the name stands; where it names a surface the solid has no face on, the
    // shipped fence refuses the vertex by that symptom (`VertexNamesAbsentSurface`), not the
    // cause. Lateral faces carry no plane and are skipped —
    // the cylinder is in every pierce definition already.
    let mut planes_at: HashMap<(usize, NodeId), Vec<usize>> = HashMap::new();
    let mut def_triple: HashMap<(usize, NodeId), Def> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        for ring in lf.poly_rings() {
            for &node in &ring.nodes {
                let at = planes_at.entry((g, node)).or_default();
                if let ClassIx::Plane(c) = lf.surf
                    && !at.contains(&c)
                {
                    at.push(c);
                }
            }
        }
    }
    let mut keys: Vec<(usize, NodeId)> = planes_at.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        let mut at = planes_at.remove(&key).unwrap_or_default();
        at.sort_unstable();
        let Some((own, ..)) = combinatorics::pierce_name(key.1) else {
            if let Some(t) = combinatorics::canonical_triple(jd, &at) {
                def_triple.insert(key, Def::Three(t));
            }
            continue;
        };
        let def = if own.iter().all(|p| at.contains(p)) || at.len() < 2 {
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
