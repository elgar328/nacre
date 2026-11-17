use super::*;
use std::collections::BTreeMap;
/// The first topological defect in a boolean's output, checked **per solid** so one bad piece is
/// never masked by another (a sum-of-all Euler parity misses two odd-χ solids; a global count
/// misses much), or `None` if every result solid is a valid closed 2-manifold. In one pass per
/// solid it builds each edge's using faces (`edge_uses`), the vertex set and the `V/E/F/L/S` counts,
/// then, most-specific first:
///  0. an edge whose uses lie in **two shells** — a shell is one closed boundary, so one split
///     between two shells leaves each open though the solid's edge counts are whole and its Euler
///     count reads the extra shell as genus ([`RejectReason::OpenResultShell`]);
///  1. a **non-manifold vertex** (pinch) — [`nacre_topo::nonmanifold_vertices`]; sound for any
///     number of pinches (unlike Euler parity, which a second pinch flips back to even);
///  2. an edge stated on **one plane twice** — an interior boundary the coplanar merge had to
///     erase ([`RejectReason::CoplanarMerge`]);
///  3. an edge whose stated carriers are **not the surfaces its two faces lie on**
///     ([`RejectReason::EdgeCarrierMismatch`]);
///  4. an **odd Euler characteristic** `χ = V − E + F − L_i` (a valid closed solid's is even);
///  5. a **negative genus** `S − χ/2 < 0` (an even χ that still cannot be a solid).
///
/// A boolean output is freshly built, so its cells never alias another live solid's — the scoped
/// maps are exact.
pub(crate) fn check_result_topology(
    model: &Model,
    solids: &[Handle<Solid>],
) -> Option<(RejectReason, Option<crate::RejectWhere>)> {
    for &sh in solids {
        let solid = model.solid(sh);
        let mut shells: HashSet<Handle<Shell>> = HashSet::new();
        let mut faces: HashSet<Handle<Face>> = HashSet::new();
        let mut edge_uses: HashMap<Handle<Edge>, Vec<Handle<Face>>> = HashMap::new();
        let mut vertices: HashSet<Handle<Vertex>> = HashSet::new();
        let mut inner_loops: i64 = 0;
        let mut edge_shell: HashMap<Handle<Edge>, Handle<Shell>> = HashMap::new();
        let mut across: Vec<Handle<Edge>> = Vec::new();
        for &shell_h in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
            if !shells.insert(shell_h) {
                continue;
            }
            for &fh in &model.shell(shell_h).faces {
                if !faces.insert(fh) {
                    continue;
                }
                let face = model.face(fh);
                inner_loops += face.inner.len() as i64;
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        edge_uses.entry(he.edge).or_default().push(fh);
                        vertices.extend(model.edge(he.edge).vertices);
                        if *edge_shell.entry(he.edge).or_insert(shell_h) != shell_h {
                            across.push(he.edge);
                        }
                    }
                }
            }
        }
        // **Each shell closes on its own.** Witness: the smallest such edge's first vertex, for
        // replay determinism.
        if let Some(&eh) = across.iter().min() {
            let [va, _] = model.edge(eh).vertices;
            return Some((
                RejectReason::OpenResultShell,
                Some(crate::RejectWhere::Point(model.vertex_point(va))),
            ));
        }
        // The returned list is sorted by handle index, so the first entry is a deterministic
        // witness (one of possibly several pinch vertices).
        if let Some(&vh) = nacre_topo::nonmanifold_vertices(model, faces.iter().copied()).first() {
            return Some((
                RejectReason::NonManifoldVertex,
                Some(crate::RejectWhere::Point(model.vertex_point(vh))),
            ));
        }
        // **No edge may say it separates one plane from itself.** An edge whose *stated*
        // surface pair is one plane twice is an interior boundary the coplanar cleaning was
        // obliged to erase. When the merge abstains instead — a 2-node cap's chord and its
        // arcs collide in the node-pair key, so the group ships as-is — the pair would sail
        // through without this, and `validate` names the spelling a producer bug
        // (`EdgeCarrierMismatch`: no edge separates a surface from itself). Measured on the
        // straddling flush boss (`rul flush` Fuse). An *arc*-carried boundary between two
        // coplanar faces states a cylinder in its pair and ships by design, so only a stated
        // plane-self pair is asked here — a cylinder-self pair is refused one door earlier, by
        // the edge push (`TwoCylinders`). The reject is the cleaning's own name: an abstention that
        // leaves this shape is a merge that was mandatory and did not happen. Witness: the
        // offending edge's first vertex, smallest edge handle for replay determinism.
        let mut same_plane: Vec<Handle<Edge>> = Vec::new();
        for &eh in edge_uses.keys() {
            let e = model.edge(eh);
            let [s0, s1] = e.surfaces;
            if s0 == s1 && matches!(model.surface(s0), nacre_topo::Surface::Plane { .. }) {
                same_plane.push(eh);
            }
        }
        if let Some(&eh) = same_plane.iter().min() {
            let [va, _] = model.edge(eh).vertices;
            return Some((
                RejectReason::CoplanarMerge,
                Some(crate::RejectWhere::Point(model.vertex_point(va))),
            ));
        }
        // **And the other half of that proposition: the stated pair is the pair of surfaces
        // the two using faces lie on** — `validate`'s rule, read the same way: two faces on two
        // surfaces must *be* the stated pair; two faces on one surface (two laterals of one
        // cylinder sharing an arc) must have that surface among the two. The
        // edge's curve derives from the stated pair, so a pair the faces do not keep is a wrong
        // curve on a shell every count above calls closed, and only `validate` would say so — the
        // shape is a straight edge on a cylinder's ruling stated `(cylinder, wall)` between two
        // plane faces (`an_edge_stated_on_a_pair_its_faces_do_not_keep_is_refused` plants one).
        // An edge with other than two uses is the shell guard's to name.
        let mut disagreeing: Vec<Handle<Edge>> = Vec::new();
        for (&eh, uses) in &edge_uses {
            let [f0, f1] = uses[..] else {
                continue;
            };
            let stated = model.edge(eh).surfaces;
            let mut observed = [model.face(f0).surface, model.face(f1).surface];
            if observed[1].index() < observed[0].index() {
                observed.swap(0, 1);
            }
            let agrees = if observed[0] == observed[1] {
                stated.contains(&observed[0])
            } else {
                stated == observed
            };
            if !agrees {
                disagreeing.push(eh);
            }
        }
        if let Some(&eh) = disagreeing.iter().min() {
            let [va, _] = model.edge(eh).vertices;
            return Some((
                RejectReason::EdgeCarrierMismatch,
                Some(crate::RejectWhere::Point(model.vertex_point(va))),
            ));
        }
        let (v, e, f) = (
            vertices.len() as i64,
            edge_uses.len() as i64,
            faces.len() as i64,
        );
        let chi = v - e + f - inner_loops;
        if chi % 2 != 0 {
            return Some((RejectReason::EulerParity, None));
        }
        if shells.len() as i64 - chi / 2 < 0 {
            return Some((RejectReason::NegativeGenus, None));
        }
    }
    None
}

/// Connected components of the reconstructed faces, joined **where two of them bound one material
/// across a boundary curve**: round each curve, a face joins the next face round on its material
/// side — the two faces bounding one kept sector of the curve's star
/// ([`LineStar::bounding_pairs`](crate::draft::LineStar::bounding_pairs)). Returns a component label
/// per face (dense `0..n` in order of first appearance, for replay determinism) and the component
/// count.
///
/// **A curve two faces use joins them without reading anything** — two spokes have one cyclic
/// order, and each is the other's next face either way round. That is nearly every curve. Where
/// more faces use one, the curve's kind supplies the spokes' order and each face's material side,
/// and the one rule pairs them: [`circle_star`] for a whole rim or a cut arc. A star that cannot be
/// read is refused by its cause, never left unjoined — an unjoined face is what turns two bodies
/// into a false `NonManifoldResultEdge`.
///
/// **Touching is not joining.** Two bodies that meet along a line or a circle share it — four
/// faces use it, two of each — and each face's next on its material side is its own body's, so
/// they land in two components and come back as the two bodies they are, not as one pinched
/// pseudo-solid. Where a body touches **itself** the material still runs around the contact, and
/// the component stays one: the guard fires exactly where no pair of solids exists.
///
/// This is the 3D reading of a rule the assembly stage already applies one dimension down:
/// [`unify_coplanar_faces`] splits a coplanar group "into **edge-connected components**, because
/// faces that merely lie on the same plane without touching must each survive on their own".
///
/// An enclosed void is its own component — its boundary shares no edge with the outer.
pub(super) fn face_components(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
) -> Result<(Vec<usize>, usize), BoolError> {
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    // The faces round each ring edge, with where each meets it. One use is a dangling edge — the
    // edge-use guard's to name. Both tables are walked in key order, so a refusal names the same
    // curve run after run.
    let mut users: BTreeMap<NodeEdge, Vec<(usize, Wall)>> = BTreeMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for ring in rings_of(lf) {
            for t in 0..ring.len() {
                users
                    .entry(NodeEdge::of(ring, t))
                    .or_default()
                    .push((i, ring.walls[t]));
            }
        }
    }
    // **A whole rim** has no node at all — a circle has none — so a cap face and the lateral that
    // meet on it are keyed by the rim, exactly as `reconstruct` keys it, and where each meets it is
    // read off the kind of bound.
    let mut rims: BTreeMap<(usize, usize), Vec<(usize, RimSpoke)>> = BTreeMap::new();
    for (i, lf) in faces.iter().enumerate() {
        let bounds = std::iter::once((&lf.outer, true)).chain(lf.inner.iter().map(|b| (b, false)));
        for (b, outer) in bounds {
            match (b, lf.surf) {
                (Bound::Circle { cyl }, ClassIx::Plane(c)) => {
                    let spoke = if outer {
                        RimSpoke::CapInside
                    } else {
                        RimSpoke::CapOutside
                    };
                    rims.entry((*cyl, c)).or_default().push((i, spoke));
                }
                (Bound::Rim { plane, ccw }, ClassIx::Cyl(k)) => {
                    let spoke = if *ccw {
                        RimSpoke::LateralAbove
                    } else {
                        RimSpoke::LateralBelow
                    };
                    rims.entry((k, *plane)).or_default().push((i, spoke));
                }
                _ => {}
            }
        }
    }
    let mut join = |a: usize, b: usize| {
        let (ra, rb) = (uf_find(&mut parent, a), uf_find(&mut parent, b));
        parent[ra] = rb;
    };
    for (key, us) in &users {
        match (us.len(), key) {
            (0 | 1, _) => {}
            (2, _) => join(us[0].0, us[1].0),
            (_, NodeEdge::Arc { cyl, plane, .. }) => {
                let spokes = us
                    .iter()
                    .map(|&(i, wall)| Ok((i, arc_spoke(jd, cyls, faces, i, wall, *cyl)?)))
                    .collect::<Result<Vec<_>, BoolError>>()?;
                for (a, b) in circle_star(jd, cyls, faces, *cyl, *plane, &spokes)? {
                    join(a, b);
                }
            }
            // A straight edge more faces use: not paired yet.
            (_, NodeEdge::Line(_)) => {}
        }
    }
    for (&(k, c), us) in &rims {
        match us.len() {
            0 | 1 => {}
            2 => join(us[0].0, us[1].0),
            _ => {
                for (a, b) in circle_star(jd, cyls, faces, k, c, us)? {
                    join(a, b);
                }
            }
        }
    }
    let mut label: HashMap<usize, usize> = HashMap::new();
    let mut labels = Vec::with_capacity(faces.len());
    let mut n = 0;
    for i in 0..faces.len() {
        let root = uf_find(&mut parent, i);
        let next = label.len();
        let l = *label.entry(root).or_insert(next);
        labels.push(l);
        n = n.max(l + 1);
    }
    Ok((labels, n))
}

/// **Where a face meets a circle, round it** — the four directions, in the cross-section through
/// the axis, a face can leave a cylinder's circle on a plane along: the cap inside the circle, the
/// lateral above the plane, the cap outside it, the lateral below. In that cyclic order whatever
/// the plane's tilt, and these are all there are — a circle lies on its plane and its cylinder
/// alone (an arc's key names both).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RimSpoke {
    CapInside = 0,
    LateralAbove = 1,
    CapOutside = 2,
    LateralBelow = 3,
}

impl RimSpoke {
    /// The four, in their cyclic order.
    const ROUND: [RimSpoke; 4] = [
        RimSpoke::CapInside,
        RimSpoke::LateralAbove,
        RimSpoke::CapOutside,
        RimSpoke::LateralBelow,
    ];
}

/// **Where face `i` meets the cut arc of cylinder class `k` it walks along `wall`** — the arc's
/// travel against the face's frame says which side of it the face lies on, the way a whole circle's
/// bound kind does. Its ring keeps the face on its left about the frame's normal before `flip`
/// (`arrangement::emit`), so travelling counter-clockwise about the axis a lateral lies above the
/// arc (the convention `Bound::Rim` states), and a cap lies inside it exactly when its frame faces
/// up the axis (the even half-edge borders the disk exactly then — `emit_faces`' arc labels).
///
/// A lateral on the arc is one of cylinder `k`'s, the cylinder the arc's key names: an arc on two
/// cylinders' laterals is a contact the gate refuses (`CylinderPairContact`), so another reads as
/// the stages disagreeing.
fn arc_spoke(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    i: usize,
    wall: Wall,
    k: usize,
) -> Result<RimSpoke, BoolError> {
    let Wall::Arc { ccw, plane, .. } = wall else {
        return Err(reject(RejectReason::CylinderStagesDisagree));
    };
    Ok(match faces[i].surf {
        ClassIx::Cyl(c) if c != k => return Err(reject(RejectReason::CylinderStagesDisagree)),
        ClassIx::Cyl(_) if ccw => RimSpoke::LateralAbove,
        ClassIx::Cyl(_) => RimSpoke::LateralBelow,
        ClassIx::Plane(_) => {
            let up = crate::planes::frame_faces_up_the_axis(&jd.planes[plane], &cyls[k].def)
                .ok_or_else(|| reject(jd.planes[plane].world_absent()))?;
            if up == ccw {
                RimSpoke::CapInside
            } else {
                RimSpoke::CapOutside
            }
        }
    })
}

/// **The faces that bound one kept sector round the circle of cylinder class `k` on plane class
/// `c`** — `users` are the faces meeting it and where. Each face's material side is read off its
/// truth:
///
/// - a cap's material lies below the plane along the axis exactly when its outward normal points up
///   the axis — the class's frame facing up it ([`crate::planes::frame_faces_up_the_axis`]) and the
///   face not flipped, or the frame facing down and the face flipped;
/// - a lateral's lies inside the cylinder exactly when it is not flipped (its outward normal looks
///   away from the axis).
///
/// Two faces leaving along one direction have no order (`UnorderedEdges`); a class whose facing
/// cannot be read (by the cause of its missing world statement,
/// [`crate::planes::WorkingPlane::world_absent`], as `reconstruct` reads the same facing — the gate
/// keeps a class with no world statement off a cylinder's circles) and neighbours whose sides
/// disagree (`LabelConflict`) are refused.
fn circle_star(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    k: usize,
    c: usize,
    users: &[(usize, RimSpoke)],
) -> Result<Vec<(usize, usize)>, BoolError> {
    let mut at: [Option<usize>; 4] = [None; 4];
    for &(i, spoke) in users {
        if at[spoke as usize].replace(i).is_some() {
            return Err(reject(RejectReason::UnorderedEdges));
        }
    }
    // Whether the face's material lies the next spoke on in `RimSpoke` order.
    let forward = |i: usize, spoke: RimSpoke| -> Result<bool, BoolError> {
        let flip = faces[i].flip;
        Ok(match spoke {
            RimSpoke::CapInside | RimSpoke::CapOutside => {
                let up = crate::planes::frame_faces_up_the_axis(&jd.planes[c], &cyls[k].def)
                    .ok_or_else(|| reject(jd.planes[c].world_absent()))?;
                // Material above the plane is toward `LateralAbove`: on from inside, back from
                // outside.
                let above = up == flip;
                above == (spoke == RimSpoke::CapInside)
            }
            RimSpoke::LateralAbove | RimSpoke::LateralBelow => {
                // Material inside is toward `CapInside`: back from above, on from below.
                let inside = !flip;
                inside == (spoke == RimSpoke::LateralBelow)
            }
        })
    };
    let mut star = Vec::with_capacity(users.len());
    for (spoke, slot) in RimSpoke::ROUND.into_iter().zip(at) {
        if let Some(i) = slot {
            star.push((i, forward(i, spoke)?));
        }
    }
    crate::draft::LineStar::new(star)
        .bounding_pairs()
        .ok_or_else(|| reject(RejectReason::LabelConflict))
}
