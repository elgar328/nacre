use super::*;
use nacre_geom::intersect::RingSide;
/// **What the tangency verdict needs, carried as one parameter** — the operation (the only thing
/// that decides `keep`) and the rows the gate stated. Empty for every boolean without a tangent
/// wall, which is almost all of them.
#[derive(Clone, Copy)]
pub(crate) struct Tangencies<'a> {
    pub(crate) kind: BoolKind,
    pub(crate) rows: &'a [crate::planes::Tangency],
}

impl Tangencies<'_> {
    /// For the roads that assemble without a gate run (the band-road fixtures). The `kind` here
    /// is a placeholder, not a choice: with no rows the verdict returns before reading it, which is
    /// why naming an operation at this one site is not the casework the rule forbids.
    #[cfg(test)]
    pub(crate) fn none() -> Self {
        Self {
            kind: BoolKind::Fuse,
            rows: &[],
        }
    }
}

/// **The shape the material near a tangent line falls into, when it falls apart** — [`pinch`]'s
/// answer, one variant per adjacency the kept regions lose.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pinch {
    /// `W2 ∧ ¬L ∧ ¬F`: the two wedges between the parabola and the plane, which cannot reach each
    /// other.
    Wedges,
    /// `L ∧ F ∧ ¬W2`: the lens inside the cylinder and the far side, which meet only on the line.
    LensAndFar,
}

/// **Does the material near a tangent line fall into more than one piece, and which?** — the three
/// regions and their two adjacencies, asked through [`crate::draft::keep`] and nothing else.
///
/// **The operation enters here and only here.** Near the tangent line the material is three
/// regions — the lens `L` inside the cylinder, the **two** wedges `W2` between the parabola and the
/// plane, and the far half-space `F` — and the adjacencies are `L–W2` and `F–W2` only: `L` and `F`
/// meet along the line itself and nowhere else, because at `v = 0` the wedges are *inside* the
/// cylinder. So the kept regions fall apart exactly as the two [`Pinch`] variants say.
///
/// Each region's `(in_A, in_B)` is local: `lens_in_wall_solid` says whether the cylinder's side of
/// the wall plane is the wall **face**'s material side, and `cyl_orient` is `+1` for a boss
/// (material inside) and `-1` for a bore. Those two bits fix all six memberships, because near the
/// line the only boundaries are those two surfaces — which is why this is exact without touching
/// the arrangement. (It is exact only where no face on another plane class reaches the tangent
/// line; such a plane is a secant, so the picture there is six regions, and the population gate
/// refuses the row (`TangentLineInAnotherPlane`) — or, where that line is an edge both faces end
/// on, hands the question to the structure. A plane that holds the line with no face reaching it
/// divides nothing there.)
pub(crate) fn pinch(
    kind: BoolKind,
    wall_solid: crate::planes::SolidSide,
    lens_in_wall_solid: bool,
    cyl_orient: i8,
) -> Option<Pinch> {
    use crate::planes::SolidSide;
    // `(in the wall face's solid, in the lateral face's solid)`.
    let member = |in_lens_side: bool, inside_cyl: bool| {
        (
            in_lens_side == lens_in_wall_solid,
            inside_cyl == (cyl_orient > 0),
        )
    };
    let kept = |(w, c): (bool, bool)| match wall_solid {
        SolidSide::A => crate::draft::keep(kind, w, c),
        SolidSide::B => crate::draft::keep(kind, c, w),
    };
    let (ml, mf) = (member(true, true), member(false, false));
    let (l, w2, f) = (kept(ml), kept(member(true, false)), kept(mf));
    if w2 && !l && !f {
        Some(Pinch::Wedges)
    } else if l && f && !w2 {
        // What the verdict's reading of this shape stands on: the lens is the cylinder's
        // operand's material alone and the far side the wall's alone.
        debug_assert!(
            ml == (false, true) && mf == (true, false),
            "the lens and the far side are the two operands' own"
        );
        Some(Pinch::LensAndFar)
    } else {
        None
    }
}

/// **A tangency pinches the result when the material it leaves is two lumps that meet only along
/// the line — and both lumps end up in one solid.**
///
/// **Two lumps is not yet a defect** — and that was measured, not assumed. A boss tangent to a
/// wall *from outside* leaves the lens and the far side, and the honest answer is **two solids
/// touching along a line**, which the engine produces (`Ok(2)`, `validate` clean, and it meshes).
/// What makes it a defect is the two lumps landing in **one** body — and the [`Pinch`] says how to
/// read that:
///
/// - **The wedges** are both bounded by the one wall face the line runs through (`runs_through`,
///   this verdict's premise), and a tangency splits no face — it mints no vertex and no edge — so that
///   face ties them into one component of the result. A statement about the engine's result, not
///   the geometry: two pieces touching along two lines still share the face.
/// - **The lens and the far side** arise only where `keep` keeps both, which is the union of a boss
///   and a wall it touches from outside (the truth table in `bands`). The lens is then the
///   cylinder's operand's material alone and the far side the wall's alone — two different solids
///   (the gate keeps no row of a solid tangent to itself, `cylinder_gate`'s `retain`). An operand is
///   one lump (a solid is one outer shell and its cavities, and a boolean hands back one solid per
///   material component), so the union is one body or two, and the lumps are in one body exactly
///   when the result is **one** body: two is the valid pair touching along the line; one is the two
///   lumps joined elsewhere — a boss standing in a notch, tangent to the notch's wall and
///   overlapping the block beside it — whose boundary touches itself, which nothing else in the
///   kernel sees. (Bodies are face components joined across edges two faces use — where that
///   differs from the material, at an edge four faces use, another guard refuses first.)
///
/// **The shape is the whole verdict — no body is asked which classes it holds.** A row states a
/// real contact: the gate writes one only for a wall face that touches the line where the lateral
/// face lies on it, and `runs_through` says the line crosses that face's interior there
/// (holes and gaps read, `planes::cylinder_gate`'s `line_runs_through_face`). The wedges then lie
/// against that face and the lateral in one body, and the lens and the far side are in one body
/// exactly when the result is one. Asking whether a body holds a face on the wall's class and one
/// on the cylinder's is not the rule: read alone, it convicts the lens and the far side of two
/// bodies whenever one of them has a face anywhere on the wall's plane.
///
/// **Every row here is one this verdict can read.** A row that could not be stated, or whose line
/// a third class's face reaches other than as an edge (six regions, not three), was refused by the
/// gate before anything was arranged (`planes::cylinder_gate`); what reaches here is stated, and
/// no other class's face reaches its line, or the line is such an edge.
pub(super) fn tangency_reject(tg: &Tangencies<'_>, g: &Grouping) -> Result<(), BoolError> {
    for t in tg.rows {
        // The line is an edge here (`Tangency::line_is_an_edge`): two lumps meeting on it meet
        // on an edge four faces use, and the shell guard names that — this row has nothing to
        // add, and the six regions around a line a third class's face reaches are not its to judge. Said here
        // rather than left to the `!runs_through` acquittal below, which such a row also meets: that
        // is the gate's premise, and this verdict does not lean on it.
        if t.line_is_an_edge {
            continue;
        }
        // Two lumps, and the contact is a *segment* — a corner grazing the line at a point is a
        // valid tangency (measured), and this must not convict it.
        let Some(shape) = pinch(tg.kind, t.wall_solid, t.lens_in_wall_solid, t.cyl_orient) else {
            continue;
        };
        if !t.runs_through {
            continue;
        }
        let together = match shape {
            Pinch::Wedges => true,
            Pinch::LensAndFar => g.positives.len() == 1,
        };
        if together {
            return Err(crate::reject_at(
                RejectReason::SelfTouchingResult,
                crate::RejectWhere::Point(t.witness),
            ));
        }
    }
    Ok(())
}

/// Push the reconstructed result and supersede the inputs.
///
/// **The operands retire in exactly one place, and it is after the result is accepted.** Copied
/// into each arm that produces solids, every reject would have to know whether it stands before or
/// after its arm's copy; in one place, "a reject leaves the model alone" is true of the *shape* of
/// this function rather than of where the rejects happen to sit. [`reconstruct`](fn@reconstruct)
/// does the work and never touches the live set.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    rims: &crate::draft::HeldRims,
    deferred: Option<BoolError>,
    tangencies: Tangencies<'_>,
    memo: &mut crate::realize::PlaneMemo,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let out = reconstruct(
        model, jd, seam, faces, cyls, rims, deferred, tangencies, memo,
    )?;
    model.supersede_live(&[a, b]);
    Ok(out)
}

/// **A solid whose surface touches itself is not a solid.**
///
/// The proposition: *an edge of this solid lies in the interior of one of this same solid's faces,
/// or a vertex of it on the closed region of one it does not bound.* An embedded boundary cannot do
/// that — the surface would occupy the same points twice — so a result that does is refused rather
/// than returned. The unit is the **result solid**, outer shell
/// and cavities together: a wedge cut whose tip lands on the far wall touches its own cavity, and
/// asking the question per component would miss exactly that.
///
/// **Three sieves, cheapest first**, because the exact question is the expensive one:
///
/// 1. **Which faces are near** — a box per face against the edge's box (its two ends' reach),
///    the faces filed in a grid over their boxes ([`FaceGrid`]). Over the perf folds, 666,653
///    pairs.
/// 2. **Which of those lie under the edge** — both ends exactly on the face's plane
///    (`side_of`), 1,583 of those pairs. Not the planes the two ends' **names** share: a name
///    carries three of the planes through its point, and where three planes hold the edge's line (a
///    pencil) it carries two of them, which two decided by class order — that is, by operand
///    order, and an edge resting on the other operand's face went unasked in one of the two.
/// 3. The exact test, on what is left: [`combinatorics::segment_meets_face`] — **the open segment**,
///    not its endpoints. The endpoints alone miss the shape that motivated this: a cut taken all the
///    way *through* leaves a contact line whose ends sit on the face's ring while its middle crosses
///    the interior, and a point test reads both ends as "on the boundary, not inside" and passes a
///    body that cannot exist (measured: `Ok(n=1)`, volume exact, `validate` silent). An end a
///    cylinder made is asked like a three-plane one — that test pins both kinds on the line — so a
///    lateral's own ruling edge, where the lateral ends on a flat face of the solid, is asked against
///    the solid's other plane faces like any straight edge.
///
/// **One judge, not two.** The endpoint test is *subsumed* — an endpoint strictly inside puts that
/// end in an inside stretch of the line, so the segment test answers it too. Measured across the
/// corpus on edges with three-plane ends: 1,696 candidate pairs, the two verdicts agreed on every
/// one, and no pair was undecidable. Keeping both would answer the same question in two places.
///
/// **The boxes hold where the truth can be, so they can only over-keep.** They are built from
/// realized `f64` coordinates, which are rounded; a box used to *reject* a candidate before an
/// exact test must therefore be conservative, or a real self-contact is dropped in silence. Each
/// vertex contributes its reach ([`crate::planes::cache_reach`]) — the interval its cache proves
/// the truth lies in — so the box holds the true face, and a query asks whether the endpoint's
/// reach meets it. The case this protects is not hypothetical: the wedge that motivated this check
/// touches at exactly `x = 1`, which is a face of its own box.
///
/// A face's bounds — its rings and its circles (`comp_face`) — are read **lazily and once per
/// face**: `Ring::edges` spends an exact predicate per node, and building them for every result
/// face is the same unconditional cost that put 12x on the star when the void label did it.
///
/// **The curved arm** asks the same proposition of a **lateral** face, for the one edge family
/// that can lie inside it: an edge two of the solid's plane faces share along a ruling of that
/// cylinder. The sieve is the gate's record of such lines (`SharedLines::on_cylinder` — a line on
/// a lateral *face* of an operand, which every result lateral lies within), so a boolean whose
/// cylinders record none skips the arm, and the census asks it of no edge at all. The exact test is
/// [`combinatorics::ruling_meets_lateral`]. The shapes it refuses are ones no other guard can
/// see — each one body whose counts, `validate` and volume are all right: a plate's bore with a
/// wedge whose apex edge lies along the bore wall and that crosses the hole into the plate beyond,
/// and a cylinder with a C-shaped prism touching it from outside and running through it
/// elsewhere. A test whose exact question cannot be formed is held, not raised — a self-touch
/// found elsewhere is still the answer. Two bodies touching on that line
/// are two groups and are not asked; the line a lateral's own slit cut is its boundary, which the
/// test abstains on, and the edge-use guard names that shape.
pub(super) fn self_touch_reject(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    seam: &[SeamVertex],
    groups: &[Vec<usize>],
    by_comp_lf: &[Vec<&LocalFace>],
) -> Result<(), BoolError> {
    // Per node: the coordinate a reject names as its witness, and the reach the sieve reads.
    let pt: HashMap<NodeId, (Point3, crate::planes::CacheReach)> = seam
        .iter()
        .map(|sv| {
            (
                sv.triple,
                (sv.cache.coord(), crate::planes::cache_reach(&sv.cache)),
            )
        })
        .collect();
    // **The witness is the smallest offending edge, else the smallest offending vertex**, not
    // the first a hash map hands over — the rejection's `at` is a census key's neighbour and must
    // not move between runs. Only the rejecting path keeps scanning, so a result that builds pays
    // nothing for it.
    let mut found: Option<Touch> = None;
    let mut found_at = |e: Touch| {
        found = Some(found.map_or(e, |f| f.min(e)));
    };
    // **A proof beats a «cannot tell».** A candidate whose exact test cannot be formed (a ring
    // the walk cannot name, a ray it cannot cast) is held, not raised: the scan goes on, and a
    // self-touch found anywhere is the answer. Only a result with none raises the held error —
    // the one met at the smallest edge or vertex, so its kind does not follow the hash map's order
    // either.
    let mut held: Option<(Touch, BoolError)> = None;
    let mut hold = |e: Touch, err: BoolError| {
        if held.as_ref().is_none_or(|(k, _)| e < *k) {
            held = Some((e, err));
        }
    };
    // The curved arm's sieve, built the first time a solid with a lateral asks it.
    let mut shared: Option<combinatorics::SharedLines> = None;
    let any_shared = cyls.iter().any(|c| !c.shared.is_empty());
    for g in groups {
        // **Two arms, one proposition** — "an edge of this solid lies in the interior of one of
        // its own faces". The planar arm asks it of plane faces, through boxes and exact plane
        // membership. The curved arm asks it of a lateral face, for the edges two of the
        // solid's plane faces share along a ruling of that cylinder (the gate's record,
        // `SharedLines::on_cylinder`) — a wedge's apex edge lying along the solid's own bore wall,
        // say ([`combinatorics::ruling_meets_lateral`]). The point arm below asks a vertex
        // against a plane face. What none asks abstains rather than pretends, each for a reason
        // written down: an arc edge lies in a plane face only on its own circle's plane, where
        // the faces overlap in area (the coplanar merge's); and a lateral face touched at a single
        // point is not asked — where a plane face through that point is slanted to the axis the gate
        // refuses it first (`ObliqueCylinderCut`), and on planes along or across the axis the case
        // is not measured.
        let all: Vec<&LocalFace> = g
            .iter()
            .flat_map(|&c| by_comp_lf[c].iter().copied())
            .collect();
        let faces: Vec<&LocalFace> = all
            .iter()
            .copied()
            .filter(|lf| matches!(lf.surf, ClassIx::Plane(_)))
            .collect();
        let laterals: Vec<(usize, &LocalFace)> = all
            .iter()
            .filter_map(|lf| match lf.surf {
                ClassIx::Cyl(k) => Some((k, *lf)),
                ClassIx::Plane(_) => None,
            })
            .collect();
        // One pass: which two faces own each edge, which edges an arc carries, and a box per face
        // holding its vertices' reach and its arcs' circles.
        let mut owners: HashMap<[NodeId; 2], Vec<usize>> = HashMap::new();
        let mut visits: Vec<(NodeId, usize)> = Vec::new();
        let mut arc_edges: std::collections::HashSet<[NodeId; 2]> = Default::default();
        let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(faces.len());
        for (j, lf) in faces.iter().enumerate() {
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            let mut grow = |rl: [f64; 3], rh: [f64; 3]| {
                for k in 0..3 {
                    lo[k] = lo[k].min(rl[k]);
                    hi[k] = hi[k].max(rh[k]);
                }
            };
            for r in rings_of(lf) {
                let k = r.nodes.len();
                for i in 0..k {
                    let (u, v) = (r.nodes[i], r.nodes[(i + 1) % k]);
                    let key = if u < v { [u, v] } else { [v, u] };
                    owners.entry(key).or_default().push(j);
                    visits.push((u, j));
                    if let Some(&(_, (rl, rh))) = pt.get(&u) {
                        grow(rl, rh);
                    }
                    // An arc bulges past its two ends: the box takes its whole circle's, and with
                    // no exact circle it holds everything — a box with no proved bound must not
                    // drop a candidate.
                    if let Wall::Arc { cyl, .. } = r.walls[i] {
                        arc_edges.insert(key);
                        match cyls.get(cyl).and_then(|c| {
                            combinatorics::arc_box_rat(
                                jd,
                                lf.surf.plane(),
                                &c.def,
                                crate::planes::Theta::Whole,
                            )
                        }) {
                            Some([cl, ch]) => grow(
                                cl.map(|x| outward(x.to_f64(), -1.0)),
                                ch.map(|x| outward(x.to_f64(), 1.0)),
                            ),
                            None => grow([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]),
                        }
                    }
                }
            }
            boxes.push((lo, hi));
        }
        let grid = FaceGrid::new(&boxes);
        let mut near_buf = NearFaces::new(faces.len());
        // Does node `n` lie on plane `q`? Exact, asked afresh: neither the answers nor the
        // node's implicit point pay for a table (measured over the perf folds: 932k distinct
        // (node, plane) pairs of 1.33M asks; a table of implicit points, to share a turned
        // class's Cramer parts, cost more in hashing than it saved).
        #[expect(
            clippy::disallowed_methods,
            reason = "a point question: does result node `n` lie on face plane `q`"
        )]
        let on = |n: NodeId, q: usize| combinatorics::side_of(jd, cyls, n, q) == Some(0);
        let reach_of = |n: NodeId| {
            pt.get(&n)
                .map_or(([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]), |&(_, r)| r)
        };
        let mut read: HashMap<usize, combinatorics::CompFace> = HashMap::new();
        'edges: for (&[u, v], own) in &owners {
            // An arc lies in a plane face only on its own circle's plane, where the face would
            // overlap the arc's own face in area — the coplanar merge's case, not a line contact.
            if arc_edges.contains(&[u, v]) {
                continue;
            }
            let ((ul, uh), (vl, vh)) = (reach_of(u), reach_of(v));
            let el: [f64; 3] = core::array::from_fn(|k| ul[k].min(vl[k]));
            let eh: [f64; 3] = core::array::from_fn(|k| uh[k].max(vh[k]));
            for &j in near_buf.of(&grid, &boxes, el, eh) {
                if own.contains(&j) {
                    continue;
                }
                // **Which faces: the ones on whose plane both ends lie.** Not the planes the two
                // ends' names share — a name carries three of the planes through its point, and
                // where three planes hold the edge's line (a pencil) it carries two of them, which
                // two decided by class order, that is by operand order.
                let q = faces[j].surf.plane();
                if !(on(u, q) && on(v, q)) {
                    continue;
                }
                // A face bounded outside by a circle — a cap disk — is not asked. Its plane is
                // ⊥ to the cylinder's axis, and an edge of the same body lying on it from
                // outside needs a plane slanted to that axis through a point of the disk,
                // which cannot clear the cylinder's side: `ObliqueCylinderCut` refuses first.
                if matches!(faces[j].outer, Bound::Circle { .. }) {
                    continue;
                }
                if let std::collections::hash_map::Entry::Vacant(slot) = read.entry(j) {
                    match super::grouping::comp_face(jd, cyls, faces[j]) {
                        Ok(face) => {
                            slot.insert(face);
                        }
                        Err(err) => {
                            hold(Touch::Edge([u, v]), err);
                            continue;
                        }
                    }
                }
                // The line is `q ∩ w`, and `w` comes from **the edge's own two faces**, not
                // from the endpoint names and not from the ring's recorded wall. Both of those
                // fail, measured: an edge can lie on a pencil of three planes and the ring is
                // free to call it by any of them — including `q` itself — and where four planes
                // concur the endpoints' canonical triples need not share a second plane at all
                // (31 pairs in the rotation sweep). An edge's own plane faces give one, because
                // that is what the edge *is* — a lateral's ruling edge has one, the flat face it
                // ends on, and both its ends lie on `q` and on that plane.
                let Some(w2) = own.iter().map(|&o| faces[o].surf.plane()).find(|&x| x != q) else {
                    // Every plane face of the edge lies on `q`: face `j` holds the edge inside only
                    // where it overlaps that face in area — the coplanar merge's case.
                    continue;
                };
                match combinatorics::segment_meets_face(jd, cyls, q, w2, u, v, &read[&j]) {
                    Ok(true) => {
                        found_at(Touch::Edge([u, v]));
                        continue 'edges;
                    }
                    Ok(false) => {}
                    Err(err) => hold(Touch::Edge([u, v]), err),
                }
            }
        }
        // **The point arm** — a vertex of this solid on the closed region of one of its own plane
        // faces it does not bound: a corner of one part resting on a face of another part the body
        // joins elsewhere. No edge of it lies on that face, so the edges above never ask; inside the
        // face or on its boundary, the surface touches itself at that point. The faces are the
        // edges' sieve's: near the vertex's reach, its plane exactly through the vertex.
        //
        // `OnBoundary` counts: a vertex on another face's ring where it is no node of it touches
        // there whether the arrangement split one of the two faces at it (the edge-use guard would
        // name the open shell further on) or neither (the counts close and nothing else sees it).
        // The cost: a producer's T-junction — a node one face lists and its neighbour across the
        // edge does not — would read here as a self-touch, `Impossible` where `SuspectedDefect` is
        // the truth; the suite holds no such vertex (measured).
        //
        // A pierce vertex abstains. A plane face through it is slanted to the cylinder's axis, along
        // it, or across it: the first the gate refuses (`ObliqueCylinderCut`); the other two are not
        // measured here.
        visits.sort_unstable();
        for at in visits.chunk_by(|x, y| x.0 == y.0) {
            let n = at[0].0;
            let Some(name) = combinatorics::three_plane_name(n) else {
                continue;
            };
            let mine = |j: usize| at.binary_search(&(n, j)).is_ok();
            let (nl, nh) = reach_of(n);
            for &j in near_buf.of(&grid, &boxes, nl, nh) {
                if mine(j) || matches!(faces[j].outer, Bound::Circle { .. }) {
                    continue;
                }
                let q = faces[j].surf.plane();
                if !on(n, q) {
                    continue;
                }
                if let std::collections::hash_map::Entry::Vacant(slot) = read.entry(j) {
                    match super::grouping::comp_face(jd, cyls, faces[j]) {
                        Ok(face) => {
                            slot.insert(face);
                        }
                        Err(err) => {
                            hold(Touch::Point(n), err);
                            continue;
                        }
                    }
                }
                match combinatorics::named_point_in_face(jd, cyls, q, name, &read[&j]) {
                    Ok(Some(RingSide::Inside | RingSide::OnBoundary)) => {
                        found_at(Touch::Point(n));
                        break;
                    }
                    Ok(_) => {}
                    Err(err) => hold(Touch::Point(n), err),
                }
            }
        }
        if laterals.is_empty() || !any_shared {
            continue;
        }
        let shared = shared.get_or_insert_with(|| combinatorics::SharedLines::of(jd, cyls));
        // The lateral's loops, read once per face and only for a face some edge is asked against;
        // one the loops road cannot read (`RingNaming`) abstains like every other unreadable case.
        let mut loops_of: HashMap<usize, Option<Vec<combinatorics::LateralLoop>>> = HashMap::new();
        'ruling_edges: for (&[u, v], own) in &owners {
            // Two faces on one plane (or one face passing the edge twice) ask `on_cylinder(p, p, k)`,
            // which is false: the gate records a line of two classes, never of one with itself.
            let [f, h] = own[..] else { continue };
            let w = [faces[f].surf.plane(), faces[h].surf.plane()];
            for (li, &(k, lf)) in laterals.iter().enumerate() {
                if !shared.on_cylinder(w[0], w[1], k) {
                    continue;
                }
                let loops =
                    loops_of.entry(li).or_insert_with(|| {
                        match super::grouping::comp_face(jd, cyls, lf).ok()?.outer {
                            combinatorics::BoundEdges::Lateral(loops) => Some(loops),
                            _ => None,
                        }
                    });
                let Some(loops) = loops else { continue };
                if combinatorics::ruling_meets_lateral(jd, &cyls[k].def, loops, w, [u, v])
                    == Some(true)
                {
                    found_at(Touch::Edge([u, v]));
                    continue 'ruling_edges;
                }
            }
        }
    }
    if let Some(touch) = found {
        // The offending edge itself, as the witness — or the end the seam table realized — or the
        // offending vertex.
        let at = match touch {
            Touch::Edge([u, v]) => match (pt.get(&u), pt.get(&v)) {
                (Some(a), Some(b)) => Some(crate::RejectWhere::Segment([a.0, b.0])),
                (Some(a), None) | (None, Some(a)) => Some(crate::RejectWhere::Point(a.0)),
                (None, None) => None,
            },
            Touch::Point(n) => pt.get(&n).map(|a| crate::RejectWhere::Point(a.0)),
        };
        return Err(match at {
            Some(at) => crate::reject_at(RejectReason::SelfTouchingResult, at),
            None => reject(RejectReason::SelfTouchingResult),
        });
    }
    if let Some((_, err)) = held {
        return Err(err);
    }
    Ok(())
}

/// A face's rings, outer first.
pub(super) fn rings_of(lf: &LocalFace) -> impl Iterator<Item = &Ring> {
    lf.poly_rings()
}

/// `x` moved one rounding step further `away` (`−1` down, `+1` up) — a rounded rational bound
/// made to hold the rational.
fn outward(x: f64, away: f64) -> f64 {
    x + away * crate::planes::round_slack(x)
}

/// **Where a self-touch was found** — on an edge of the body, or at a vertex of it with no edge on
/// the face it touches. Ordered edges first, each by its names, so the witness and the held error
/// are the same run after run.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Touch {
    Edge([NodeId; 2]),
    Point(NodeId),
}

/// **The faces filed by their boxes** — a uniform grid over the group's finite boxes, about `∛F`
/// cells a side, so a box reads the faces of the cells it covers rather than every face. A sweep
/// on one axis is not enough: the fins of a turned fold overlap on every axis, and a window on the
/// faces sorted by low `x` read about 120 faces per edge there (`y` and `z` read more) — over the
/// perf folds both arms took 242 ms through the window and 196 ms through this grid (measured).
/// A query collects into one buffer and drops repeats by a stamp ([`NearFaces`]); collecting a
/// fresh vector per query and sorting it to drop them was slower than the window. A face whose box
/// has no finite bound (an arc whose circle has no exact description) is filed everywhere; an empty
/// box (a ring-less cap disk) nowhere.
struct FaceGrid {
    lo: [f64; 3],
    step: [f64; 3],
    n: [usize; 3],
    cells: Vec<Vec<usize>>,
    everywhere: Vec<usize>,
}

impl FaceGrid {
    fn new(boxes: &[([f64; 3], [f64; 3])]) -> FaceGrid {
        let finite = |(lo, hi): &([f64; 3], [f64; 3])| {
            (0..3).all(|k| lo[k].is_finite() && hi[k].is_finite() && lo[k] <= hi[k])
        };
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for b in boxes.iter().filter(|b| finite(b)) {
            for k in 0..3 {
                lo[k] = lo[k].min(b.0[k]);
                hi[k] = hi[k].max(b.1[k]);
            }
        }
        let side = ((boxes.len() as f64).cbrt().ceil() as usize).clamp(1, 32);
        let (mut n, mut step) = ([1; 3], [1.0; 3]);
        for k in 0..3 {
            if hi[k] > lo[k] {
                n[k] = side;
                step[k] = (hi[k] - lo[k]) / side as f64;
            }
        }
        let mut grid = FaceGrid {
            lo,
            step,
            n,
            cells: vec![Vec::new(); n[0] * n[1] * n[2]],
            everywhere: Vec::new(),
        };
        for (j, b) in boxes.iter().enumerate() {
            if finite(b) {
                for c in grid.cells_of(b.0, b.1) {
                    grid.cells[c].push(j);
                }
            } else if (0..3).all(|k| b.0[k] <= b.1[k]) {
                grid.everywhere.push(j);
            }
        }
        grid
    }

    /// The cells a box covers, clamped to the grid.
    fn cells_of(&self, lo: [f64; 3], hi: [f64; 3]) -> impl Iterator<Item = usize> + use<> {
        let r: [[usize; 2]; 3] = core::array::from_fn(|k| {
            let at = |x: f64| -> usize {
                let i = ((x - self.lo[k]) / self.step[k]).floor();
                if i.is_nan() || i < 0.0 {
                    0
                } else {
                    (i as usize).min(self.n[k] - 1)
                }
            };
            [at(lo[k]), at(hi[k])]
        });
        let n = self.n;
        (r[0][0]..=r[0][1]).flat_map(move |x| {
            (r[1][0]..=r[1][1])
                .flat_map(move |y| (r[2][0]..=r[2][1]).map(move |z| (x * n[1] + y) * n[2] + z))
        })
    }
}

/// A query's buffer over a [`FaceGrid`]: the faces whose boxes meet a box, each once, in
/// ascending order — so which face a caller asks first does not follow the grid. A face filed in
/// several cells is kept once by the stamp of the query that met it.
struct NearFaces {
    stamp: Vec<u32>,
    query: u32,
    out: Vec<usize>,
}

impl NearFaces {
    fn new(faces: usize) -> NearFaces {
        NearFaces {
            stamp: vec![0; faces],
            query: 0,
            out: Vec::new(),
        }
    }

    fn of(
        &mut self,
        grid: &FaceGrid,
        boxes: &[([f64; 3], [f64; 3])],
        lo: [f64; 3],
        hi: [f64; 3],
    ) -> &[usize] {
        self.query += 1;
        self.out.clear();
        let meets = |j: usize| {
            let (fl, fh) = boxes[j];
            (0..3).all(|k| hi[k] >= fl[k] && lo[k] <= fh[k])
        };
        for c in grid.cells_of(lo, hi) {
            for &j in &grid.cells[c] {
                if self.stamp[j] != self.query {
                    self.stamp[j] = self.query;
                    if meets(j) {
                        self.out.push(j);
                    }
                }
            }
        }
        for &j in &grid.everywhere {
            if self.stamp[j] != self.query && meets(j) {
                self.stamp[j] = self.query;
                self.out.push(j);
            }
        }
        self.out.sort_unstable();
        &self.out
    }
}
