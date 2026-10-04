use super::*;
/// The first topological defect in a boolean's output, checked **per solid** so one bad piece is
/// never masked by another (a sum-of-all Euler parity misses two odd-χ solids; a global count
/// misses much), or `None` if every result solid is a valid closed 2-manifold. In one pass per
/// solid it builds the reverse-index maps (`edge_uses`, `vertex_edges`) and the `V/E/F/L/S` counts,
/// then, most-specific first:
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
        let mut edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>> = HashMap::new();
        let mut vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>> = HashMap::new();
        let mut edges_seen: HashSet<Handle<Edge>> = HashSet::new();
        let mut inner_loops: i64 = 0;
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
                        edge_uses.entry(he.edge).or_default().push((fh, he.forward));
                        if edges_seen.insert(he.edge) {
                            let [va, vb] = model.edge(he.edge).vertices;
                            vertex_edges.entry(va).or_default().push(he.edge);
                            vertex_edges.entry(vb).or_default().push(he.edge);
                        }
                    }
                }
            }
        }
        // The returned list is sorted by handle index, so the first entry is a deterministic
        // witness (one of possibly several pinch vertices).
        if let Some(&vh) = nacre_topo::nonmanifold_vertices(&vertex_edges, &edge_uses).first() {
            return Some((
                RejectReason::NonManifoldVertex,
                Some(crate::RejectWhere::Point(model.vertex_point(vh))),
            ));
        }
        // ★ **No edge may say it separates one plane from itself.** An edge whose *stated*
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
        // ★ **And the other half of that proposition: the stated pair is the pair of surfaces
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
            let [(f0, _), (f1, _)] = uses[..] else {
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
            vertex_edges.len() as i64,
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

/// Connected components of the reconstructed faces, joined **only across a manifold contact** — a
/// ring edge that exactly two of them use. Returns a component label per face (dense `0..n` in
/// order of first appearance, for replay determinism) and the component count.
///
/// ★★ **Touching is not joining.** Two bodies that meet along a line share that line's nodes but
/// no *manifold* edge — four faces use the contact line, not two — so they land in two components
/// and come back as the two bodies they are, not as one pinched pseudo-solid. Where a body touches
/// **itself** the material still runs around the contact, every edge of that path is used twice,
/// and the component stays one: the guard fires exactly where no pair of solids exists.
///
/// This is the 3D reading of a rule this file already applies one dimension down:
/// [`unify_coplanar_faces`] splits a coplanar group "into **edge-connected components**, because
/// faces that merely lie on the same plane without touching must each survive on their own".
///
/// An enclosed void is its own component, as before — its boundary shares no edge with the outer.
pub(super) fn face_components(faces: &[LocalFace]) -> (Vec<usize>, usize) {
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    /// The joining key — the third appearance of "a line is unordered, a circle is ordered"
    /// (`Wall` is the type half, `EdgeKey` the handle-space half; this is the node-space half).
    /// A 2-node circle folds its chord and both complementary arcs into one `norm_edge` pair, so
    /// an unordered key counts six users where each piece really has two — the CCW order is the
    /// bit that keeps them apart, exactly as it is in the edge welding.
    #[derive(Clone, Copy, PartialEq, Eq, Hash)]
    enum JoinKey {
        Line((NodeId, NodeId)),
        Arc {
            cyl: usize,
            from: NodeId,
            to: NodeId,
        },
    }
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    // Which faces use each ring edge. A count other than two is not a contact between neighbours:
    // one is a dangling edge and more is a pinch, and both are the edge-use guard's to name.
    let mut users: HashMap<JoinKey, Vec<usize>> = HashMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for ring in rings_of(lf) {
            let k = ring.len();
            for t in 0..k {
                let (a, b) = (ring[t], ring[(t + 1) % k]);
                // ★ A ruling step is a `Line` here: straight, so its unordered pair names it —
                // what the order bit keeps apart is a circle's complementary pieces, and two
                // points bound one straight segment whichever vocabulary a face states it in.
                let key = match ring.walls[t] {
                    Wall::Plane(_) | Wall::Ruling { .. } => JoinKey::Line(norm_edge(a, b)),
                    Wall::Arc { cyl, ccw } => {
                        let (from, to) = if ccw { (a, b) } else { (b, a) };
                        JoinKey::Arc { cyl, from, to }
                    }
                };
                users.entry(key).or_default().push(i);
            }
        }
    }
    // ★ **The second joining rule**: a cap face and the lateral that meets it on a whole rim
    // share **no node at all** — a circle has none — so the node rule above would leave a drilled
    // box's wall in its own component and send the result down the cavity/containment branch.
    // They do share the rim, keyed exactly as `reconstruct` keys it.
    //
    // The condition is the node rule's, deliberately: join only where **exactly two** faces use
    // the rim. More than two is a contact, not a join, and the closed-shell guard is what names
    // it.
    let mut rim_users: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for b in std::iter::once(&lf.outer).chain(lf.inner.iter()) {
            match (b, lf.surf) {
                (Bound::Circle { cyl }, ClassIx::Plane(c)) => {
                    rim_users.entry((*cyl, c)).or_default().push(i)
                }
                (Bound::Rim { plane, .. }, ClassIx::Cyl(k)) => {
                    rim_users.entry((k, *plane)).or_default().push(i)
                }
                _ => {}
            }
        }
    }
    for us in users.values().chain(rim_users.values()) {
        if let [a, b] = us[..] {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            parent[ra] = rb;
        }
    }
    let mut label: HashMap<usize, usize> = HashMap::new();
    let mut labels = Vec::with_capacity(faces.len());
    let mut n = 0;
    for i in 0..faces.len() {
        let root = find(&mut parent, i);
        let next = label.len();
        let l = *label.entry(root).or_insert(next);
        labels.push(l);
        n = n.max(l + 1);
    }
    (labels, n)
}
