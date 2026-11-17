use super::*;
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

/// Connected components of the reconstructed faces, joined **only where two of them bound one
/// material across a boundary** — a ring edge that exactly two of them use, or a whole rim where a
/// face meets the next one round the rim on its material side. Returns a component label per face
/// (dense `0..n` in order of first appearance, for replay determinism) and the component count.
///
/// **Touching is not joining.** Two bodies that meet along a line share that line's nodes but
/// no *manifold* edge — four faces use the contact line, not two — so they land in two components
/// and come back as the two bodies they are, not as one pinched pseudo-solid. Where a body touches
/// **itself** the material still runs around the contact, every edge of that path is used twice,
/// and the component stays one: the guard fires exactly where no pair of solids exists. A face
/// whose every boundary is such a contact line has no other way to its body, so it lands in a
/// component of its own. Where the grouping places it by containment, one closed boundary comes
/// out split between two shells and the shell guard refuses it (`OpenResultShell` — a D pin under
/// its D bore); where the grouping abstains, the edge-use guard names the line
/// (`NonManifoldResultEdge` — a square pin under its square bore).
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
) -> (Vec<usize>, usize) {
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    // Which faces use each ring edge. A count other than two is not a contact between neighbours:
    // one is a dangling edge and more is a pinch, and both are the edge-use guard's to name.
    let mut users: HashMap<NodeEdge, Vec<usize>> = HashMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for ring in rings_of(lf) {
            for t in 0..ring.len() {
                users.entry(NodeEdge::of(ring, t)).or_default().push(i);
            }
        }
    }
    for us in users.values() {
        if let [a, b] = us[..] {
            let (ra, rb) = (uf_find(&mut parent, a), uf_find(&mut parent, b));
            parent[ra] = rb;
        }
    }
    // **The second joining rule: a whole rim.** A cap face and the lateral that meets it on a
    // whole rim share **no node at all** — a circle has none — so the node rule above would leave
    // a drilled box's wall in its own component and send the result down the cavity/containment
    // branch. They do share the rim, keyed exactly as `reconstruct` keys it.
    //
    // Round the rim, in the cross-section through the axis, the faces that can meet it stand in
    // one fixed order: the cap inside the circle, the lateral above the plane, the cap outside,
    // the lateral below ([`RimSpoke`]). **A face is joined with the next face round the rim on
    // its material side** — the two bound the same material there. Two faces are each other's
    // next face whichever way round, so a rim two faces use always joins them (a drilled box, a
    // blind bore's floor, a boss on a plate); four faces — two bodies meeting on the circle, a pin
    // under the bore of a block it touches — alternate material and empty round it and join in
    // the two pairs that bound each body.
    let mut rim_users: HashMap<(usize, usize), Vec<(usize, RimSpoke)>> = HashMap::new();
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
                    rim_users.entry((*cyl, c)).or_default().push((i, spoke));
                }
                (Bound::Rim { plane, ccw }, ClassIx::Cyl(k)) => {
                    let spoke = if *ccw {
                        RimSpoke::LateralAbove
                    } else {
                        RimSpoke::LateralBelow
                    };
                    rim_users.entry((k, *plane)).or_default().push((i, spoke));
                }
                _ => {}
            }
        }
    }
    for (&(k, c), us) in &rim_users {
        for (a, b) in rim_pairs(jd, cyls, faces, k, c, us) {
            let (ra, rb) = (uf_find(&mut parent, a), uf_find(&mut parent, b));
            parent[ra] = rb;
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
    (labels, n)
}

/// **Where a face meets a whole rim, round it** — the four directions, in the cross-section
/// through the axis, a face can leave the rim along: the cap inside the circle, the lateral above
/// the plane, the cap outside it, the lateral below. In that cyclic order whatever the plane's tilt.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RimSpoke {
    CapInside = 0,
    LateralAbove = 1,
    CapOutside = 2,
    LateralBelow = 3,
}

impl RimSpoke {
    fn at(i: usize) -> Self {
        [
            RimSpoke::CapInside,
            RimSpoke::LateralAbove,
            RimSpoke::CapOutside,
            RimSpoke::LateralBelow,
        ][i % 4]
    }
}

/// **The faces that join across the whole rim of cylinder class `k` on plane class `c`** — each
/// face with the next face round the rim on its material side, where the two point at each other
/// ([`face_components`]' second rule). `users` are the faces meeting the rim and where; each face's
/// material side is read off its truth:
///
/// - a cap's material lies below the plane along the axis exactly when its outward normal points up
///   the axis — the class's frame facing up it ([`crate::planes::frame_faces_up_the_axis`]) and the
///   face not flipped, or the frame facing down and the face flipped;
/// - a lateral's lies inside the cylinder exactly when it is not flipped (its outward normal looks
///   away from the axis).
///
/// A face whose side cannot be read (a class with no world statement — the gate keeps such a
/// class off a cylinder's circles), two faces at one spoke, or a next face that points elsewhere
/// joins nothing: the faces stay apart and the shell guard names what that leaves.
fn rim_pairs(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    k: usize,
    c: usize,
    users: &[(usize, RimSpoke)],
) -> Vec<(usize, usize)> {
    let mut at: [Option<usize>; 4] = [None; 4];
    for &(i, spoke) in users {
        if at[spoke as usize].replace(i).is_some() {
            return Vec::new();
        }
    }
    // The step round the rim, `+1` or `−1` in `RimSpoke` order, toward the face's material.
    let toward = |i: usize, spoke: RimSpoke| -> Option<usize> {
        let flip = faces[i].flip;
        let forward = match spoke {
            RimSpoke::CapInside | RimSpoke::CapOutside => {
                let up = crate::planes::frame_faces_up_the_axis(&jd.planes[c], &cyls[k].def)?;
                // Material above the plane is toward `LateralAbove`: `+1` from inside, `−1` from
                // outside.
                let above = up == flip;
                above == (spoke == RimSpoke::CapInside)
            }
            RimSpoke::LateralAbove | RimSpoke::LateralBelow => {
                // Material inside is toward `CapInside`: `−1` from above, `+1` from below.
                let inside = !flip;
                inside == (spoke == RimSpoke::LateralBelow)
            }
        };
        let step = if forward { 1 } else { 3 };
        (1..4)
            .map(|n| RimSpoke::at(spoke as usize + n * step))
            .find_map(|s| at[s as usize])
    };
    let mut out = Vec::new();
    for &(i, spoke) in users {
        let Some(j) = toward(i, spoke) else { continue };
        let back = users
            .iter()
            .find(|&&(x, _)| x == j)
            .and_then(|&(x, s)| toward(x, s));
        if i < j && back == Some(i) {
            out.push((i, j));
        }
    }
    out
}
