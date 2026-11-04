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
/// The proposition: *an edge of this solid meets the closed region of one of this same solid's
/// faces it does not bound at a point that is no vertex — inside the face, or where it touches a
/// curved bound — or a vertex of it lies on the closed region of one it does not bound.* An
/// embedded boundary cannot do that — the surface would occupy the same points twice — so a result
/// that does is refused rather than returned. The unit is the **result solid**, outer shell
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
///    the solid's other plane faces like any straight edge. A face bounded by a circle — a cap
///    disk — is asked like any other. The test reads the face's **closed** region where no vertex
///    stands: a line tangent to a circle or an arc of the face's bounds touches it there, the shape
///    of an edge lying on a cap's plane tangent to its rim (measured before the test read it:
///    `Ok(n=1)`, volume exact, `validate` silent). A polygon bound is touched at no such point — a
///    line through a point inside one of its sides runs along it or enters the face, and one through
///    a corner meets a vertex, the point arm's to ask.
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
/// **The lateral arms** ask the same proposition of a **lateral** face: a vertex of the solid on
/// its closed region ([`combinatorics::named_point_in_face`]'s lateral reading) and a straight
/// edge meeting it away from the edge's ends ([`combinatorics::segment_meets_lateral`]) — along a
/// ruling the gate recorded (`WorkingCyl::shared_of`, asked by
/// [`combinatorics::ruling_meets_lateral`]) or at the roots of the edge's line on the cylinder. The
/// shapes they refuse are ones no other guard can see — each one body whose counts, `validate` and
/// volume are all right: a plate's bore with a wedge whose apex edge lies along the bore wall and
/// that crosses the hole into the plate beyond, a cylinder with a C-shaped prism touching it from
/// outside and running through it elsewhere, a corner resting on a boss's side, an edge tangent to
/// it. A point on the stretch of a tangent line a tangency row speaks for, which a result face on
/// the row's wall covers, is a point of that line contact, the tangency verdict's (`tangent_hand`);
/// an edge meeting the lateral's boundary where a face it bounds shares that curve with the lateral
/// is that face's two bounds meeting, its pinch (`shares_a_curve`). Two bodies touching there are
/// two groups and are not asked; the line a lateral's own slit cut is its boundary, which the
/// reading abstains on, and the edge-use guard names that shape.
pub(super) fn self_touch_reject(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    seam: &[SeamVertex],
    groups: &[Vec<usize>],
    by_comp_lf: &[Vec<&LocalFace>],
    tg: &Tangencies<'_>,
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
    for g in groups {
        // **The arms, one proposition** — "an edge or a vertex of this solid meets the closed
        // region of one of its own faces it does not bound". The planar arm asks an edge of every
        // plane face, through boxes and exact plane membership; the point arm below asks a vertex
        // against every plane face; the lateral arms after it ask both of every lateral face. What
        // none asks abstains rather than pretends, each for a reason written down. An arc edge:
        // where it lies in a plane face the faces overlap in area (the coplanar merge's), and where
        // it touches one at a point no vertex names, a straight edge of that face touches the arc's
        // own cap there (the planar arm's touch), or the face lies on a wall tangent along the
        // ruling (the tangency verdict), or it is slanted to the axis (`ObliqueCylinderCut`), or the
        // face's bound is another cylinder's arc (`CylinderPairContact`) — the last two lean on the
        // gate refusing first. A vertex resting on a rim: the arrangement cuts the rim's two faces
        // at it, so it bounds both, and the result check names the pinch (`NonManifoldVertex`). A
        // point of a lateral on the stretch of a tangent line a tangency row speaks for, a result
        // face of the row's wall covering it: a point of that line contact, which the tangency
        // verdict judges (`tangent_hand`). An edge meeting a lateral's boundary where a face it
        // bounds shares that curve with the lateral: that face's two bounds meeting, its pinch
        // (`shares_a_curve`). And what a reading cannot place, dropped — the point arm's (a point
        // with no rational coordinates against a mixed ring, and a point the mixed ring's
        // coordinate walk abstains on: over the suite, the census, the ignored sweep and nacre-kit,
        // 722 such pairs, every one off the face when read approximately) and the lateral arms'
        // (over the same corpus, 44 vertices and 73 edges at a root, beside the 100 edges
        // the ruling reading drops; these are not read approximately).
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
            // A circle bulges past every node, so a curved bound takes its whole circle's box — and
            // with no exact circle the box holds everything: a box with no proved bound must not
            // drop a candidate.
            let circle_box = |cyl: usize| match cyls.get(cyl).and_then(|c| {
                combinatorics::arc_box_rat(jd, lf.surf.plane(), &c.def, crate::planes::Theta::Whole)
            }) {
                Some([cl, ch]) => (
                    cl.map(|x| outward(x.to_f64(), -1.0)),
                    ch.map(|x| outward(x.to_f64(), 1.0)),
                ),
                None => ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]),
            };
            // A face bounded outside by a circle — a cap disk — has no ring to grow from; its
            // box is its circle's.
            if let Bound::Circle { cyl } = lf.outer {
                let (cl, ch) = circle_box(cyl);
                grow(cl, ch);
            }
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
                    // An arc bulges past its two ends: the box takes its whole circle's.
                    if let Wall::Arc { cyl, .. } = r.walls[i] {
                        arc_edges.insert(key);
                        let (cl, ch) = circle_box(cyl);
                        grow(cl, ch);
                    }
                }
            }
            boxes.push((lo, hi));
        }
        let grid = FaceGrid::new(&boxes);
        let mut near_buf = NearFaces::new(faces.len());
        // Which side of plane `q` node `n` lies on, `0` on it, `None` unread. Exact, asked afresh:
        // neither the answers nor the node's implicit point pay for a table (measured over the
        // perf folds: 932k distinct (node, plane) pairs of 1.33M asks; a table of implicit points,
        // to share a turned class's Cramer parts, cost more in hashing than it saved).
        let side = |n: NodeId, q: usize| node_side(jd, cyls, n, q);
        let on = |n: NodeId, q: usize| node_side(jd, cyls, n, q) == Some(0);
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
        // the truth; the suite holds no such vertex, three-plane or pierce-named (measured).
        //
        // Every vertex is asked, whatever kind its name is: a point the alias table knows on a
        // cylinder's surface — a three-plane corner there included — is represented by its pierce
        // name (`Aliases::rep_rank`), and filtering on the name's kind let such a corner pressing a
        // face of its own body through. The reading is `named_point_in_face`'s; what it cannot read
        // it abstains on, and an abstention here is dropped.
        visits.sort_unstable();
        for at in visits.chunk_by(|x, y| x.0 == y.0) {
            let n = at[0].0;
            let mine = |j: usize| at.binary_search(&(n, j)).is_ok();
            let (nl, nh) = reach_of(n);
            for &j in near_buf.of(&grid, &boxes, nl, nh) {
                if mine(j) {
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
                match combinatorics::named_point_in_face(jd, cyls, n, &read[&j]) {
                    Ok(Some(RingSide::Inside | RingSide::OnBoundary)) => {
                        found_at(Touch::Point(n));
                        break;
                    }
                    Ok(_) => {}
                    Err(err) => hold(Touch::Point(n), err),
                }
            }
        }
        // **The lateral arms** — the same proposition asked of each lateral face of this solid:
        // a vertex it does not bound on its closed region, and a straight edge (two plane faces,
        // so it bounds no lateral) meeting that region away from its ends. The readings are
        // `named_point_in_face`'s lateral arm and `segment_meets_lateral` — the second reads an
        // edge along a ruling the gate recorded as the ruling question
        // ([`combinatorics::ruling_meets_lateral`]), any other line at its roots. A point the
        // lateral is met at is then asked whose question it is (`tangent_hand`): on the stretch of
        // a tangent line where a tangency row says the wall face meets this lateral, and covered
        // by a result face on that wall's class, it is a point of a line contact, which the
        // tangency verdict judges; anywhere else it is a touch at a point. An edge meeting the
        // lateral's boundary where a face it bounds shares that curve is asked nothing here: the
        // point is that face's two bounds meeting (`shares_a_curve`). A lateral face whose
        // bounds the road cannot read holds its failure for the candidates those bounds decide —
        // a point on the cylinder's surface, an edge whose line meets it — as the planar arms hold
        // theirs past their plane test; and a reading that cannot answer is dropped, as the point
        // arm drops its abstentions. The bounds are read at the first candidate past the sieve,
        // which is the lateral's box: the hull of the whole circles its bounds lie on, which holds
        // the face; with no exact circle it holds everything.
        if laterals.is_empty() {
            continue;
        }
        let mut nodes: Vec<NodeId> = visits.iter().map(|&(n, _)| n).collect();
        nodes.sort_unstable();
        nodes.dedup();
        for &(k, lf) in &laterals {
            let mut cf: Option<Result<combinatorics::CompFace, BoolError>> = None;
            let (mut bl, mut bh) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            let mut grow = |plane: usize, cyl: usize| {
                let (cl, ch) = match cyls.get(cyl).and_then(|c| {
                    combinatorics::arc_box_rat(jd, plane, &c.def, crate::planes::Theta::Whole)
                }) {
                    Some([cl, ch]) => (
                        cl.map(|x| outward(x.to_f64(), -1.0)),
                        ch.map(|x| outward(x.to_f64(), 1.0)),
                    ),
                    None => ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]),
                };
                for i in 0..3 {
                    bl[i] = bl[i].min(cl[i]);
                    bh[i] = bh[i].max(ch[i]);
                }
            };
            let mut own: std::collections::HashSet<NodeId> = Default::default();
            // A lateral is bounded by rims and rings; a circle bound is a cap's, which no lateral
            // carries — read as a box that holds everything rather than one that drops a candidate.
            let mut unbounded = false;
            for b in std::iter::once(&lf.outer).chain(&lf.inner) {
                match b {
                    Bound::Rim { plane, .. } => grow(*plane, k),
                    Bound::Ring(r) => {
                        own.extend(r.nodes.iter().copied());
                        for w in &r.walls {
                            if let Wall::Arc { cyl, plane, .. } = *w {
                                grow(plane, cyl);
                            }
                        }
                    }
                    Bound::Circle { .. } => unbounded = true,
                }
            }
            if unbounded {
                (bl, bh) = ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]);
            }
            let near =
                |(l, h): ([f64; 3], [f64; 3])| (0..3).all(|i| l[i] <= bh[i] && h[i] >= bl[i]);
            for &n in &nodes {
                if own.contains(&n) || !near(reach_of(n)) {
                    continue;
                }
                let cf = match cf.get_or_insert_with(|| super::grouping::comp_face(jd, cyls, lf)) {
                    Ok(cf) => &*cf,
                    // Held only where the bounds decide: a point off the cylinder's surface is
                    // off the face whatever they are.
                    Err(err) => {
                        if !matches!(
                            combinatorics::named_point_in_disk(jd, cyls, n, &cyls[k].def),
                            Some(RingSide::Inside | RingSide::Outside)
                        ) {
                            hold(Touch::Point(n), *err);
                        }
                        continue;
                    }
                };
                if !matches!(
                    combinatorics::named_point_in_face(jd, cyls, n, cf),
                    Ok(Some(RingSide::Inside | RingSide::OnBoundary))
                ) {
                    continue;
                }
                match tangent_hand(jd, cyls, &faces, &mut read, tg.rows, &side, k, n, &[]) {
                    Hand::Touch => found_at(Touch::Point(n)),
                    Hand::Tangency => {}
                    Hand::Held(err) => hold(Touch::Point(n), err),
                }
            }
            for (&[u, v], owned) in &owners {
                // An edge of two faces on one plane (or one face passing it twice) has no line of
                // two planes; an arc is no straight edge.
                let [f, h] = owned[..] else { continue };
                let w = [faces[f].surf.plane(), faces[h].surf.plane()];
                if w[0] == w[1] || arc_edges.contains(&[u, v]) {
                    continue;
                }
                let ((ul, uh), (vl, vh)) = (reach_of(u), reach_of(v));
                let el: [f64; 3] = core::array::from_fn(|i| ul[i].min(vl[i]));
                let eh: [f64; 3] = core::array::from_fn(|i| uh[i].max(vh[i]));
                if !near((el, eh)) {
                    continue;
                }
                let loops = match cf.get_or_insert_with(|| super::grouping::comp_face(jd, cyls, lf))
                {
                    Ok(combinatorics::CompFace {
                        outer: combinatorics::BoundEdges::Lateral(loops),
                        ..
                    }) => &*loops,
                    Ok(_) => continue,
                    // Held only where the bounds decide: a line that misses the cylinder meets
                    // no face of it.
                    Err(err) => {
                        if !matches!(
                            jd.line_meet(w[0], w[1], &cyls[k].def, false),
                            nacre_judge::cylinder::LineMeet::Miss
                                | nacre_judge::cylinder::LineMeet::AxisParallelMiss
                        ) {
                            hold(Touch::Edge([u, v]), *err);
                        }
                        continue;
                    }
                };
                match combinatorics::segment_meets_lateral(jd, cyls, k, loops, w, [u, v]) {
                    Some(combinatorics::EdgeOnLateral::AlongRuling) => {
                        found_at(Touch::Edge([u, v]));
                    }
                    Some(combinatorics::EdgeOnLateral::AtRoots(roots)) => {
                        for (at, on_lateral) in roots {
                            // On the lateral's boundary, where a face this edge bounds shares a
                            // curve with the lateral: that face's two bounds meeting, its pinch.
                            if on_lateral == RingSide::OnBoundary
                                && owned.iter().any(|&o| shares_a_curve(faces[o], lf, k))
                            {
                                continue;
                            }
                            match tangent_hand(
                                jd, cyls, &faces, &mut read, tg.rows, &side, k, at, owned,
                            ) {
                                Hand::Touch => found_at(Touch::Edge([u, v])),
                                Hand::Tangency => {}
                                Hand::Held(err) => hold(Touch::Edge([u, v]), err),
                            }
                        }
                    }
                    Some(combinatorics::EdgeOnLateral::Clear) | None => {}
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

/// Whose question a point the lateral arms found on a lateral is — [`tangent_hand`]'s answer.
enum Hand {
    /// A touch at a point: the self-touch check's.
    Touch,
    /// A point of a line contact a tangency row speaks for: the tangency verdict's.
    Tangency,
    /// The point touches the lateral, and whose question it is could not be read.
    Held(BoolError),
}

/// **Is the point `n`, found on a lateral of cylinder `k`, a point of a line contact the tangency
/// verdict judges?** — it lies on a tangent line on the stretch where a row says the wall face
/// meets that lateral (`planes::Tangency::contact`, closed), and a result face on the wall's class
/// covers it (a node of the face's rings, or read inside it or on its boundary). The row's
/// stretch is the verdict's ground, so the two ask one population: a corner meeting the line at a
/// point, a face reaching the line past the lateral's span, a wall tangent to its own solid's
/// cylinder — no row speaks for those, and the point is a touch. The row is the operands'; the
/// result face is asked too, because where the boolean cut the wall face away around `n` the
/// lateral is met there by something else.
///
/// Whatever cannot be read — which side of the wall's plane `n` lies on, its height, a face's
/// reading — holds its failure: the touch is proved, only its owner is unknown, and a proof
/// elsewhere still wins. A reading with no error of its own holds `WitnessNotRational`, the name
/// for a point those readings have no statement of (none over the suite, the census, the ignored
/// sweep and nacre-kit).
#[allow(clippy::too_many_arguments)]
fn tangent_hand(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[&LocalFace],
    read: &mut HashMap<usize, combinatorics::CompFace>,
    rows: &[crate::planes::Tangency],
    side: &impl Fn(NodeId, usize) -> Option<i8>,
    k: usize,
    n: NodeId,
    own: &[usize],
) -> Hand {
    use core::cmp::Ordering;
    let unread = || reject(RejectReason::WitnessNotRational);
    let mut held = None;
    let mut height = None;
    for t in rows.iter().filter(|t| t.cyl == k) {
        match side(n, t.wall) {
            Some(0) => {}
            Some(_) => continue,
            None => {
                held = Some(unread());
                continue;
            }
        }
        let Some(v) = *height.get_or_insert_with(|| {
            combinatorics::axis_param_on_cylinder(jd, cyls, n, &cyls[k].def)
        }) else {
            held = Some(unread());
            continue;
        };
        let mut within = Some(false);
        for [lo, hi] in &t.contact {
            match (lo.cmp_value(&v), v.cmp_value(hi)) {
                (Some(a), Some(b)) if a != Ordering::Greater && b != Ordering::Greater => {
                    within = Some(true);
                    break;
                }
                (Some(_), Some(_)) => {}
                _ => within = None,
            }
        }
        match within {
            Some(true) => {}
            Some(false) => continue,
            None => {
                held = Some(unread());
                continue;
            }
        }
        // A point inside an edge of a face is on that face; a node of a face's rings is too.
        if own.iter().any(|&o| faces[o].surf.plane() == t.wall) {
            return Hand::Tangency;
        }
        for (j, lf) in faces.iter().enumerate() {
            if lf.surf.plane() != t.wall {
                continue;
            }
            if rings_of(lf).any(|r| r.nodes.contains(&n)) {
                return Hand::Tangency;
            }
            if let std::collections::hash_map::Entry::Vacant(slot) = read.entry(j) {
                match super::grouping::comp_face(jd, cyls, lf) {
                    Ok(face) => {
                        slot.insert(face);
                    }
                    Err(err) => {
                        held = Some(err);
                        continue;
                    }
                }
            }
            match combinatorics::named_point_in_face(jd, cyls, n, &read[&j]) {
                Ok(Some(RingSide::Inside | RingSide::OnBoundary)) => return Hand::Tangency,
                Ok(Some(RingSide::Outside)) => {}
                Ok(None) => held = Some(unread()),
                Err(err) => held = Some(err),
            }
        }
    }
    held.map_or(Hand::Touch, Hand::Held)
}

/// Which side of plane class `q` the result node `n` lies on — `0` on it, `None` unread.
#[expect(
    clippy::disallowed_methods,
    reason = "a point question: where a result node lies against a face's plane"
)]
fn node_side(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
    q: usize,
) -> Option<i8> {
    combinatorics::side_of(jd, cyls, n, q)
}

/// **Do the plane face `g` and the lateral `lat` of cylinder `k` share a curve** — a bound of the
/// lateral on `g`'s plane (a rim, or an arc on it) and a bound of `g` on the cylinder (its circle,
/// or an arc of it)? Every point of a lateral's boundary is a point of the face across it, so an
/// element of that face meeting the lateral's boundary meets another bound of its own face there —
/// the face's pinch, which the surface takes (a pinched face is a valid solid), and which the
/// planar arms do not ask of a face's own bounds either.
fn shares_a_curve(g: &LocalFace, lat: &LocalFace, k: usize) -> bool {
    let p = g.surf.plane();
    fn bounds(lf: &LocalFace) -> impl Iterator<Item = &Bound> {
        std::iter::once(&lf.outer).chain(&lf.inner)
    }
    let lateral_on_p = bounds(lat).any(|b| match b {
        Bound::Rim { plane, .. } => *plane == p,
        Bound::Ring(r) => r
            .walls
            .iter()
            .any(|w| matches!(*w, Wall::Arc { cyl, plane, .. } if cyl == k && plane == p)),
        Bound::Circle { .. } => false,
    });
    let g_on_k = bounds(g).any(|b| match b {
        Bound::Circle { cyl } => *cyl == k,
        Bound::Ring(r) => r
            .walls
            .iter()
            .any(|w| matches!(*w, Wall::Arc { cyl, .. } if cyl == k)),
        Bound::Rim { .. } => false,
    });
    lateral_on_p && g_on_k
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
/// has no finite bound (a circle with no exact description) is filed everywhere; an empty box
/// nowhere.
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
