//! Boolean assembly: turn the arrangement engine's per-face output (`LocalFace`, named by
//! `Node` seam triples) into result solids, and the mandatory coplanar-face cleaning pass.
//!
//! The public [`boolean`] entry lives here and delegates the cell-complex work to
//! [`crate::arrangement`]; the engine calls back into [`assemble_fuse_cut`] and
//! [`unify_coplanar_faces`] to build and clean the shells (a legal module cycle).

use crate::combinatorics;
use crate::planes::{WorkingPlane, uf_find};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, he_start, reject, unordered};
use nacre_cip::predicate::{Evidence, Notes, Site};
use nacre_cip::{Decision, dir_orient3d_judge};
use nacre_geom::Surface;
use nacre_math::Point3;
use nacre_store::Handle;
use nacre_topo::{Edge, Face, HalfEdge, Loop, Model, Orientation, Shell, Solid, Vertex, VertexDef};
use std::collections::{HashMap, HashSet};

/// Boolean of two live solids (design §8 M5, overview 불리언 전략 — 정직하게 거절).
///
/// **Coverage:** planar solids. All three kinds go through the single per-plane-class arrangement
/// engine ([`crate::arrangement::boolean`]), which handles transverse, coplanar-contact,
/// coincident, contained and disjoint cases in one path and cleans its own output (coplanar-face
/// merge) so results are chainable. Anything it cannot resolve is rejected with [`BoolError`] —
/// never a silent wrong answer (DNA). Transactional: it computes the result in local structures
/// and pushes only after every degeneracy check passes, so a rejected boolean leaves the model
/// untouched.
pub fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    boolean_with_report(model, kind, a, b).map(|(solids, _)| solids)
}

/// [`boolean`], and **what the kernel had to assume to get there** — see [`BoolReport`].
///
/// A parallel entry point rather than a wider return type: the report is wanted by roughly one
/// caller in thirty, and changing `boolean`'s signature would rewrite every other one for nothing.
pub fn boolean_with_report(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(Vec<Handle<Solid>>, BoolReport), BoolError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids.contains(&a) || !model.live_solids.contains(&b) {
        return Err(BoolError::InputNotLive);
    }
    // The arrangement engine (`arrangement.rs`) is the sole boolean path: one per-plane-class 2D
    // arrangement handles transverse, coplanar-contact, coincident, contained and disjoint cases,
    // and cleans its own output (coplanar-face merge) so results are chainable.
    let snapshot = model.live_solids.clone();
    let (result, notes, _class_of) = crate::arrangement::boolean(model, kind, a, b)?;
    // Topological self-check on the assembled result (DNA: never return a malformed solid). Reject
    // rather than return, restoring the pre-op live set so the reject leaves the *live* model
    // untouched (orphaned result cells stay in the append-only arena, unreachable, as any superseded
    // solid's do). Valid results always pass, so this never false-rejects; the traversal reads the
    // topology stores directly (no adjacency rebuild, no coordinates).
    if let Some(t) = check_result_topology(model, &result) {
        model.live_solids = snapshot;
        return Err(reject(t));
    }
    Ok((result, BoolReport::of(&notes)))
}

/// [`boolean`], and **what each input face's plane became** — its plane class's representative
/// surface, which is what every result face on that plane carries.
///
/// ★★★ A caller that needs to find *its own* face in the result must ask this rather than compare
/// handles or coordinates afterwards. The engine settled "are these one plane?" with evidence —
/// including for two faces turned by **different motion chains**, where only a composed-rotation
/// proof or a coincidence within the limit can answer — and none of that is visible to a later
/// geometric test. `ops::find_face_coplanar_with` is the one caller today.
///
/// A parallel entry point for the same reason [`boolean_with_report`] is one: the map is wanted by
/// one caller, and widening `boolean`'s return would rewrite every other one for nothing.
pub(crate) fn boolean_with_classes(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(Vec<Handle<Solid>>, crate::arrangement::ClassOf), BoolError> {
    let snapshot = model.live_solids.clone();
    let (result, notes, class_of) = crate::arrangement::boolean(model, kind, a, b)?;
    if let Some(t) = check_result_topology(model, &result) {
        model.live_solids = snapshot;
        return Err(reject(t));
    }
    let _ = notes;
    Ok((result, class_of))
}

/// **What the boolean had to take on faith**, so a user can see it and act on it.
///
/// Rotated geometry has no exact zero: a plane through turned points meets another at coordinates
/// no finite precision writes down, so "these two faces are the same plane" is *proved to within a
/// distance*, never proved outright. The kernel decides such a question by proving the separation
/// is below the coincidence limit — and then this says which questions those were and how close
/// the closest call came.
///
/// **It is a diagnosis, not a prompt.** Nothing here asks the user to choose; the choice was made
/// on evidence and the evidence is here. What it is *for* is noticing an unintended coincidence —
/// two features that met because a dimension made them meet — and going back to fix the design,
/// which is a thing only the author of the model can do.
#[derive(Clone, Debug, Default)]
pub struct BoolReport {
    /// **Faces merged into one plane class on toleranced evidence — first, because everything
    /// else follows from them.** A merge decided here changes which planes exist before a single
    /// vertex is computed, so a surprise in this list explains surprises everywhere else.
    pub merges: Vec<Evidence>,
    /// How many judgements were answered by a proved coincidence rather than a proved sign.
    pub coincidences: usize,
    /// The closest call: the coincidence with the **widest** bound, the one nearest to having
    /// been wrong. `None` when nothing was assumed at all — an axis-aligned model, typically,
    /// where every question has an exact answer.
    pub loosest: Option<Evidence>,
}

impl BoolReport {
    fn of(notes: &Notes) -> BoolReport {
        let all = notes.sorted();
        let widest = |e: &Evidence| match e.outcome {
            Decision::Coincident { within } => within.exp2(),
            _ => None,
        };
        BoolReport {
            merges: all
                .iter()
                .filter(|e| matches!(e.site, Site::PlanesCoplanar { .. }))
                .copied()
                .collect(),
            coincidences: all
                .iter()
                .filter(|e| matches!(e.outcome, Decision::Coincident { .. }))
                .count(),
            // Ties keep the first in sorted order, so the answer does not depend on the schedule.
            loosest: all
                .iter()
                .filter(|e| matches!(e.outcome, Decision::Coincident { .. }))
                .max_by_key(|e| widest(e))
                .copied(),
        }
    }
}

/// The first topological defect in a boolean's output, checked **per solid** so one bad piece is
/// never masked by another (a sum-of-all Euler parity misses two odd-χ solids; a global count
/// misses much), or `None` if every result solid is a valid closed 2-manifold. In one pass per
/// solid it builds the reverse-index maps (`edge_uses`, `vertex_edges`) and the `V/E/F/L/S` counts,
/// then, most-specific first:
///  1. a **non-manifold vertex** (pinch) — [`nacre_topo::nonmanifold_vertices`]; sound for any
///     number of pinches (unlike Euler parity, which a second pinch flips back to even);
///  2. an **odd Euler characteristic** `χ = V − E + F − L_i` (a valid closed solid's is even);
///  3. a **negative genus** `S − χ/2 < 0` (an even χ that still cannot be a solid).
///
/// A boolean output is freshly built, so its cells never alias another live solid's — the scoped
/// maps are exact.
fn check_result_topology(model: &Model, solids: &[Handle<Solid>]) -> Option<RejectReason> {
    for &sh in solids {
        let solid = model.solids.get(sh);
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
            for &fh in &model.shells.get(shell_h).faces {
                if !faces.insert(fh) {
                    continue;
                }
                let face = model.faces.get(fh);
                inner_loops += face.inner.len() as i64;
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        edge_uses.entry(he.edge).or_default().push((fh, he.forward));
                        if edges_seen.insert(he.edge) {
                            let [va, vb] = model.edges.get(he.edge).vertices;
                            vertex_edges.entry(va).or_default().push(he.edge);
                            vertex_edges.entry(vb).or_default().push(he.edge);
                        }
                    }
                }
            }
        }
        if !nacre_topo::nonmanifold_vertices(&vertex_edges, &edge_uses).is_empty() {
            return Some(RejectReason::NonManifoldVertex);
        }
        let (v, e, f) = (
            vertex_edges.len() as i64,
            edge_uses.len() as i64,
            faces.len() as i64,
        );
        let chi = v - e + f - inner_loops;
        if chi % 2 != 0 {
            return Some(RejectReason::EulerParity);
        }
        if shells.len() as i64 - chi / 2 < 0 {
            return Some(RejectReason::NegativeGenus);
        }
    }
    None
}

/// Connected components of the reconstructed faces by shared `Node` — the same identity
/// `assemble_fuse_cut` welds result vertices by. Returns a component label per face (dense
/// `0..n` in order of first appearance, for replay determinism) and the component count. An
/// enclosed void is its own component: its boundary shares no vertex with the outer.
fn face_components(faces: &[LocalFace]) -> (Vec<usize>, usize) {
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    let mut owner: HashMap<Node, usize> = HashMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for &nd in lf
            .loop_nodes
            .iter()
            .chain(lf.inner.iter().flat_map(|r| r.iter()))
        {
            match owner.entry(nd) {
                std::collections::hash_map::Entry::Occupied(e) => {
                    let (ra, rb) = (find(&mut parent, *e.get()), find(&mut parent, i));
                    parent[ra] = rb;
                }
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(i);
                }
            }
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

/// A canonical, replay-stable sort key for a severed component: its outer-loop vertex
/// coordinates, sorted lexicographically. Total order for disjoint components — distinct pieces
/// occupy different space, so their coordinate multisets differ, and the lex-min vertex alone can
/// tie (identity is by `Handle`, not coordinates, so two vertices may coincide). Coordinates are a
/// derived cache used only to order multi-solid output; no judgment reads this (tol-irrelevant).
fn comp_key(model: &Model, faces: &[Handle<Face>]) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = faces
        .iter()
        .flat_map(|&fh| model.faces.get(fh).outer.half_edges.iter().copied())
        .map(|he| model.vertex_point(he_start(model, he)).as_array())
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite vertex coordinates"));
    pts
}

/// Whether a closed component shell (a set of oriented faces) is **outward**
/// (material-enclosing, positive signed volume — an outer shell) versus
/// **inward** (a void/cavity shell). The `assemble_fuse_cut` cavity-vs-outer
/// label ((5d)#5 retired the f64 signed-volume flux this replaces), reading no
/// coordinate arithmetic: only a lexicographic vertex ordering and one
/// axis-aligned plane-coefficient sign.
///
/// At the component's lexicographically-minimal vertex `v*` (min x, then y, then
/// z) the shell is a convex corner, and the material lies toward increasing
/// coordinates. So an outward shell has a `−x`-facing boundary face at `v*` (its
/// materialized outward normal `n_x < 0`), while a void's three walls all face
/// into the void (`n_x ≥ 0` at its own `v*`). Hence: **outward iff some face
/// incident to `v*` has materialized outward normal with `n_x < 0`**. Two
/// antiparallel `x`-perpendicular faces cannot share a vertex, so this `∃`-test
/// is equivalent to (and simpler than) picking the max-`|n_x|` face.
///
/// Exact for the axis-aligned M5 corpus: face normals are exactly `±eₓ/±e_y/±e_z`
/// so `sign(n_x)` is the exact sign of the plane's `x`-coefficient (times the
/// face orientation), and the `v*` search is exact coordinate ordering — both
/// hold even for non-representable coordinates (e.g. a `0.3`-offset void face).
/// Rotated shells break the "axis-aligned normal / unique x-perpendicular face"
/// premises and are CIP's job (design §9 (5d)#5, honest scope).
pub(crate) fn is_shell_outward(model: &Model, faces: &[Handle<Face>]) -> bool {
    // Lexicographically-minimal vertex over the component's outer loops.
    let mut vstar: Option<Handle<Vertex>> = None;
    let mut pstar = [f64::INFINITY; 3];
    for &fh in faces {
        for &he in &model.faces.get(fh).outer.half_edges {
            let vh = he_start(model, he);
            let p = model.vertex_point(vh).as_array();
            if p < pstar {
                pstar = p;
                vstar = Some(vh);
            }
        }
    }
    let Some(vstar) = vstar else { return false };
    // Outward iff some face at v* faces −x (materialized outward normal n_x < 0).
    // n_x's sign is the plane x-coefficient's sign times the orientation sign
    // (no normalization — exact for axis-aligned faces).
    for &fh in faces {
        let face = model.faces.get(fh);
        if !face
            .outer
            .half_edges
            .iter()
            .any(|&he| he_start(model, he) == vstar)
        {
            continue;
        }
        let Surface::Plane(plane) = model.surface(face.surface) else {
            continue;
        };
        let sign = match face.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        if plane.coefficients()[0] * sign < 0.0 {
            return true;
        }
    }
    false
}

/// Rotation-sound twin of [`is_shell_outward`]: whether a result component's shell is
/// **outward** (material, an outer shell) versus **inward** (a void/cavity shell), decided on
/// the faces' exact `WitnessPoint` definitions through `frame3`, so it is sound when the coordinates
/// are rounded irrationals (rotation). Same algorithm as the f64 [`is_shell_outward`] — the
/// lexicographically-minimal vertex `v*` is a convex extreme corner and the shell is outward
/// iff some face there has an outward normal with `n_x < 0` — but both numeric steps become
/// exact CIP predicates:
///
/// - **`v*`** by [`Judge::cmp_coord`](crate::tolerant) over each node's three-plane triple, the same
///   lex-min scan as [`loop_winding`](crate::combinatorics::loop_winding). A node's triple is its
///   own face plane plus the neighbour planes of its two loop edges (the
///   [`loop_triples`](crate::combinatorics) construction), whose meet *is* that vertex — so an
///   original corner is as implicit a point as a seam node, no mixed compare needed.
/// - **`sign(n_x)`** by [`dir_orient3d_judge`]`([1,0,0], tri…)` on the face's exact plane
///   definition ([`plane_def`](crate::tolerant), mixed-rotation safe): the x-component of the
///   RH normal, flipped by `lf.flip` to the result face's materialized outward normal.
///
/// Operates on the **pre-assembly** `LocalFace`s (not the result faces), so it never reads a
/// result vertex whose exact rotation provenance `assemble_fuse_cut` drops — every exact
/// definition it needs lives in `planes` and in the loop adjacency. Reads no `Model`. `Err`
/// on a non-simple/degenerate component (coincident nodes, a straight angle, or a non-manifold
/// edge) — an honest reject, never a silent wrong label. Routed from `assemble_fuse_cut`'s
/// per-component outward test by [`any_rotated`](crate::tolerant) (cell 3c-vi-b).
fn component_is_outward_tol(
    jd: &Judge<'_, WorkingPlane>,
    comp: &[&LocalFace],
) -> Result<bool, BoolError> {
    let planes = jd.planes;
    use nacre_scalar::{Orient, Rat};

    // The unordered edge key: `Node` is `Ord`, so order the pair canonically.
    let ekey = |a: Node, b: Node| if a <= b { (a, b) } else { (b, a) };

    // Edge -> the planes carrying it, over every loop (outer + inner): a hole-rim edge is the
    // outer edge of its wall and an inner edge of the holed face, so building over both loops
    // gives it both planes. A manifold edge yields exactly two.
    type EKey = (Node, Node);
    let mut edge_faces: HashMap<EKey, Vec<usize>> = HashMap::new();
    for lf in comp {
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            let k = ring.len();
            for t in 0..k {
                edge_faces
                    .entry(ekey(ring[t], ring[(t + 1) % k]))
                    .or_default()
                    .push(lf.plane_idx);
            }
        }
    }
    let other_plane = |a: Node, b: Node, own: usize| -> Result<usize, BoolError> {
        let ps = edge_faces
            .get(&ekey(a, b))
            .ok_or_else(|| reject(RejectReason::MissingSeam))?;
        let mut others = ps.iter().copied().filter(|&x| x != own);
        let o = others
            .next()
            .ok_or_else(|| reject(RejectReason::UnpairedSeamEdge))?;
        if others.any(|x| x != o) {
            return Err(reject(RejectReason::NonManifoldEdge)); // edge shared by >2 distinct planes
        }
        Ok(o)
    };

    // Each unique outer node -> its three-plane triple (meet = that vertex). Any incident face
    // yields a valid triple (all its planes pass through the vertex); first occurrence wins.
    let mut triple_of: HashMap<Node, [usize; 3]> = HashMap::new();
    for lf in comp {
        let ring = &lf.loop_nodes;
        let k = ring.len();
        for t in 0..k {
            let node = ring[t];
            if triple_of.contains_key(&node) {
                continue;
            }
            let prev = other_plane(ring[(t + k - 1) % k], node, lf.plane_idx)?;
            let next = other_plane(node, ring[(t + 1) % k], lf.plane_idx)?;
            if prev == next {
                return Err(reject(RejectReason::StraightAngle)); // a straight angle
            }
            let mut tri = [lf.plane_idx, prev, next];
            tri.sort_unstable();
            triple_of.insert(node, tri);
        }
    }

    // Lexicographically-minimal vertex over the unique outer nodes (`loop_winding`'s scan;
    // `Judge::cmp_coord` is exact for the rotated triples). Sort candidates for replay determinism.
    let mut nodes: Vec<Node> = triple_of.keys().copied().collect();
    nodes.sort_unstable();
    let Some((&first, rest)) = nodes.split_first() else {
        return Ok(false); // empty component
    };
    let mut lo = first;
    for &node in rest {
        let ord = (0..3)
            .map(|axis| jd.cmp_coord(triple_of[&node], triple_of[&lo], axis))
            .find(|&c| c != 0);
        match ord {
            Some(c) if c < 0 => lo = node,
            Some(_) => {}
            None => return Err(reject(RejectReason::CoincidentNodes)), // two distinct nodes coincide
        }
    }
    // ★★ **Same postcondition as `loop_winding`'s scan, and here the stake is higher**: `v*`
    // decides whether this shell is the outside or the inside, so a node that is not extreme
    // turns the solid inside out. The forward scan finds a minimum only if the lexicographic
    // relation is an order, which is a property of the predicates rather than of this loop.
    debug_assert!(
        !nodes.iter().any(|n| {
            *n != lo
                && (0..3)
                    .map(|axis| jd.cmp_coord(triple_of[n], triple_of[&lo], axis))
                    .find(|&c| c != 0)
                    == Some(-1)
        }),
        "the lexicographic scan did not find a minimum — the comparison is not an order here"
    );

    // Outward iff some outer face at v* has a result outward normal with n_x < 0. n_x's sign is
    // the RH-normal x-component (`dir_orient3d_judge` on the exact plane def), flipped by `flip`.
    for lf in comp {
        if !lf.loop_nodes.contains(&lo) {
            continue;
        }
        let tri = crate::tolerant::plane_def(planes, lf.plane_idx);
        let ex = [Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)];
        let nx = dir_orient3d_judge(ex, &tri[0], &tri[1], &tri[2], jd.standard.prec);
        let nx = if lf.flip {
            match nx {
                Orient::Positive => Orient::Negative,
                Orient::Negative => Orient::Positive,
                Orient::Zero => Orient::Zero,
            }
        } else {
            nx
        };
        if nx == Orient::Negative {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A seam vertex — a three-plane point on both `∂A` and `∂B` (2 A-planes + 1
/// B-plane, or 1 A + 2 B). Shared (one `Handle`) by every incident result piece.
pub(crate) struct SeamVertex {
    pub(crate) point: Point3,
    pub(crate) triple: [usize; 3], // sorted combined-plane indices
    pub(crate) tol: f64,
}

/// A node in a reconstructed face loop: the sorted plane triple naming a seam vertex.
/// `Eq`/`Hash` give identity dedup so an A-piece and a B-piece that meet at a seam node
/// share one result vertex/edge; `Ord` gives the deterministic node order replay needs.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub(crate) enum Node {
    Seam([usize; 3]), // sorted triple (key into the seam map)
}

/// One ring of a result face: its nodes, and **the plane each edge rides**.
///
/// ★ **The walls are carried, not derived.** Reading an edge's supporting plane back out of its two
/// endpoint names is sound only while every vertex lies on exactly three planes — see
/// [`combinatorics::ring_edges_with_walls`]. Every producer here knows the wall (the arrangement's
/// half-edge was told it; a pass-through face reads it off the edge's other face), so it hands it
/// over instead of leaving it to be guessed.
///
/// Derefs to its nodes, so the many places that only walk the ring read unchanged.
#[derive(Clone, Debug, Default)]
pub(crate) struct Ring {
    pub(crate) nodes: Vec<Node>,
    /// `walls[i]` is the plane class the edge `nodes[i] -> nodes[i+1]` rides.
    pub(crate) walls: Vec<usize>,
}

impl Ring {
    /// A ring whose walls are **derived from the node names** — for hand-built fixtures, where
    /// every vertex is a clean three-plane point so "the class the two names share besides `p`" is
    /// well defined. Production never derives; see the type's note.
    #[cfg(test)]
    pub(crate) fn from_clean_names(p: usize, nodes: Vec<Node>) -> Ring {
        let k = nodes.len();
        let walls = (0..k)
            .map(|i| {
                let (Node::Seam(a), Node::Seam(b)) = (nodes[i], nodes[(i + 1) % k]);
                a.iter()
                    .copied()
                    .find(|&c| c != p && b.contains(&c))
                    .expect("a clean fixture ring edge rides one wall")
            })
            .collect();
        Ring { nodes, walls }
    }
}

impl std::ops::Deref for Ring {
    type Target = [Node];
    fn deref(&self) -> &[Node] {
        &self.nodes
    }
}

impl Ring {
    pub(crate) fn new(nodes: Vec<Node>, walls: Vec<usize>) -> Ring {
        debug_assert_eq!(nodes.len(), walls.len(), "one wall per edge");
        Ring { nodes, walls }
    }

    /// This ring's edges, ready for the exact predicates.
    pub(crate) fn edges(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        p: usize,
    ) -> Result<Vec<combinatorics::RingEdge>, BoolError> {
        combinatorics::ring_edges_with_walls(jd, p, &seam_ring(&self.nodes), &self.walls)
    }
}

/// A reconstructed result face: which combined plane it is on, its loop as
/// nodes, and whether to flip it (cut's inside-A B-pieces).
#[derive(Clone)]
pub(crate) struct LocalFace {
    pub(crate) plane_idx: usize,
    pub(crate) loop_nodes: Ring,
    /// Hole rings, each already wound so the kept material stays on its left
    /// about the face's outward normal. Only the non-convex path ever fills this.
    pub(crate) inner: Vec<Ring>,
    pub(crate) flip: bool,
}

/// Push the reconstructed result and supersede the inputs (mirrors `assemble`).
pub(crate) fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes = jd.planes;
    // No faces means no result — `Common` of two solids that miss each other, `Cut` of a box that
    // is wholly inside what cuts it. That is an answer, not a failure: a solid is bounded by faces,
    // so a non-empty result cannot have none. The inputs are still consumed, exactly as they are on
    // any other successful boolean — the retire below sits inside the `positives` match, which this
    // early return skips, so it has to happen here too or the operands stay live.
    if faces.is_empty() {
        model.live_solids.retain(|&s| s != a && s != b);
        return Ok(Vec::new());
    }
    // Vertices (deterministic: first appearance across faces in order).
    let mut vh: HashMap<Node, Handle<Vertex>> = HashMap::new();
    let mut node_handle = |model: &mut Model, node: Node| -> Result<Handle<Vertex>, BoolError> {
        if let Some(&h) = vh.get(&node) {
            return Ok(h);
        }
        let handle = match node {
            Node::Seam(triple) => {
                // A face references a seam node whose triple was not welded into `seam` — a
                // reconstruction dropped a crossing. Reject (never panic): an unmodeled flush
                // topology must decline honestly, not abort the kernel (DNA).
                let sv = seam
                    .iter()
                    .find(|s| s.triple == triple)
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                let def = VertexDef::ThreePlane([
                    planes[triple[0]].surf,
                    planes[triple[1]].surf,
                    planes[triple[2]].surf,
                ]);
                // The coordinate and its measured tolerance travel together into the cache
                // (the arrangement made them as a pair — `three_planes` + `vertex_tol`).
                model.push_vertex(def, sv.point, Some(sv.tol))
            }
        };
        vh.insert(node, handle);
        Ok(handle)
    };
    // Materialize all vertex handles first. Outer then inner, rings in order: `vh`'s
    // first-appearance order fixes the vertex handles, and replay depends on it.
    for lf in faces {
        for &node in lf
            .loop_nodes
            .iter()
            .chain(lf.inner.iter().flat_map(|r| r.iter()))
        {
            node_handle(model, node)?;
        }
    }

    // ★ **An edge's carriers are the two faces that use it — read off the whole result, not
    // guessed from one side.** The obvious per-face answer ("my plane, plus the wall my
    // arrangement says the edge rides") is wrong exactly where the resolved four-plane
    // concurrency lives: three planes through one LINE make each face's wall a *third* plane
    // that legitimately contains the edge but does not bound it here — measured, the two faces
    // sharing such an edge named two different walls. The faces themselves are the ground
    // truth, and every ring is already in hand, so one pre-scan reads it.
    let mut pair_surfs: HashMap<(usize, usize), Vec<Handle<Surface>>> = HashMap::new();
    for lf in faces {
        let fsurf = planes[lf.plane_idx].surf;
        for r in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            let k = r.nodes.len();
            for t in 0..k {
                let (va, vb) = (vh[&r.nodes[t]], vh[&r.nodes[(t + 1) % k]]);
                pair_surfs
                    .entry(unordered(va.index() as usize, vb.index() as usize))
                    .or_default()
                    .push(fsurf);
            }
        }
    }

    // Edges keyed by unordered handle-index pair (lookup only).
    let mut edge_of: HashMap<(usize, usize), Handle<Edge>> = HashMap::new();
    // A ring's two consecutive nodes are distinct arrangement vertices, so their points differ and
    // the line through them exists. Reject rather than panic if it does not: an aborting kernel is
    // below the floor (`overview.md`: out-of-coverage input declines honestly). `SEAM_ALIAS`
    // catches the known way this happens — two triples on one point — at the seam table, where the
    // names are still in hand, so this is a backstop with no firing test (cf. `NON_MANIFOLD_EDGE`).
    let mut edge_for = |model: &mut Model,
                        va: Handle<Vertex>,
                        vb: Handle<Vertex>,
                        fallback: [Handle<Surface>; 2]|
     -> Result<Handle<Edge>, BoolError> {
        let key = unordered(va.index() as usize, vb.index() as usize);
        if let Some(&e) = edge_of.get(&key) {
            return Ok(e);
        }
        // Manifold edges (everything a green result contains) have exactly two uses. Any other
        // count is on its way to the existing non-manifold reject — the fallback (this face's
        // wall + plane) keeps construction deterministic until that reject fires, deciding
        // nothing new.
        let surfaces = match pair_surfs[&key][..] {
            [a, b] => Edge::carrier_pair(a, b),
            _ => Edge::carrier_pair(fallback[0], fallback[1]),
        };
        let e = model
            .push_edge(surfaces, [va, vb])
            .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
        edge_of.insert(key, e);
        Ok(e)
    };

    let mut face_handles = Vec::new();
    for lf in faces {
        let face_surf = planes[lf.plane_idx].surf;
        let mut ring = |model: &mut Model, r: &Ring| -> Result<Loop, BoolError> {
            let handles: Vec<Handle<Vertex>> = r.nodes.iter().map(|nd| vh[nd]).collect();
            let k = handles.len();
            let mut half_edges: Vec<HalfEdge> = (0..k)
                .map(|t| {
                    let (va, vb) = (handles[t], handles[(t + 1) % k]);
                    let e = edge_for(model, va, vb, [planes[r.walls[t]].surf, face_surf])?;
                    let forward = model.edges.get(e).vertices[0] == va;
                    Ok(HalfEdge { edge: e, forward })
                })
                .collect::<Result<Vec<_>, BoolError>>()?;
            if lf.flip {
                // Cut's inside-A B-pieces: reverse every loop and toggle the
                // orientation below, so the outward normal points into the removed
                // region and each loop still keeps material on its left.
                half_edges.reverse();
                for he in &mut half_edges {
                    he.forward = !he.forward;
                }
            }
            Ok(Loop { half_edges })
        };
        let outer = ring(model, &lf.loop_nodes)?;
        let inner: Vec<Loop> = lf
            .inner
            .iter()
            .map(|h| ring(model, h))
            .collect::<Result<Vec<_>, BoolError>>()?;
        // The plane's frame *is* the root face's orientation: `frame_sign` is
        // `sign(stored normal · that face's n_out)`, and `collect_planes` asserts that sign equals
        // `Forward`/`Reversed`. Reading it here is what used to be `planes[plane_idx].orient` — a
        // face field indexed by a plane, the shape of every bug this split exists to prevent.
        let framed = if planes[lf.plane_idx].frame_sign > 0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        let orientation = if lf.flip {
            match framed {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
            }
        } else {
            framed
        };
        face_handles.push(model.faces.push(Face {
            surface: planes[lf.plane_idx].surf,
            outer,
            inner,
            orientation,
        }));
    }
    // Closed-shell guard: a 2-manifold b-rep uses every edge exactly twice (once from each of the
    // two faces that share it). A reconstruction that emits a face set with a dangling edge (use
    // count 1) or a pinched one (>2) is not a solid — `validate` would call it `NonManifoldEdge`,
    // but `boolean` never runs `validate` on its own output, so without this the caller receives a
    // silently invalid solid. Honest-reject instead ("honest-reject > silent-wrong", overview §1).
    // Counted on the welded `Handle<Edge>`s, so it is exact and coordinate-free.
    {
        let mut uses: HashMap<Handle<Edge>, usize> = HashMap::new();
        for &fh in &face_handles {
            let f = model.faces.get(fh);
            for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &l.half_edges {
                    *uses.entry(he.edge).or_default() += 1;
                }
            }
        }
        // A use count above two is a *pinch*: the two bodies meet exactly along that edge, and no
        // 2-manifold solid contains it (the edge twin of `NonManifoldVertex`). A count of one is a
        // *dangling* edge — the assembly dropped a face, which is ours to fix, not the input's.
        if uses.values().any(|&n| n > 2) {
            return Err(reject(RejectReason::NonManifoldResultEdge));
        }
        if uses.values().any(|&n| n != 2) {
            return Err(reject(RejectReason::OpenResultShell));
        }
    }
    // Partition the faces into connected components (by shared node). One component is the
    // whole result; several mean either an enclosed void (a cavity — an inward-oriented shell)
    // or a severed operand (two or more outward, material-enclosing shells). `is_shell_outward`
    // (exact extreme-vertex sign) tells them apart. One outward component ⇒ outer shell + the
    // rest as its cavities. Several outward components ⇒ the result severed into that many
    // solids (cell 0.4); a cavity that also survives is assigned to the piece whose outer shell
    // nests it (`combinatorics::point_in_component`). No outward component is impossible.
    let (labels, n) = face_components(faces);
    let mut by_comp: Vec<Vec<Handle<Face>>> = vec![Vec::new(); n];
    for (i, &fh) in face_handles.iter().enumerate() {
        by_comp[labels[i]].push(fh);
    }
    // Outward/void label per component, routed by rotation (overhaul 3c-vi): a component with
    // any rotated plane is decided on exact `WitnessPoint` definitions (`component_is_outward_tol` over
    // its pre-assembly `LocalFace`s), else the axis-aligned f64 `is_shell_outward` — unchanged,
    // so an unrotated result is bit-identical. `positives` stays in ascending `c` order.
    let mut by_comp_lf: Vec<Vec<&LocalFace>> = vec![Vec::new(); n];
    for (i, lf) in faces.iter().enumerate() {
        by_comp_lf[labels[i]].push(lf);
    }
    let mut positives: Vec<usize> = Vec::new();
    for c in 0..n {
        let idxs: Vec<usize> = by_comp_lf[c].iter().map(|lf| lf.plane_idx).collect();
        let outward = if crate::tolerant::any_rotated(planes, &idxs) {
            component_is_outward_tol(jd, &by_comp_lf[c])?
        } else {
            is_shell_outward(model, &by_comp[c])
        };
        if outward {
            positives.push(c);
        }
    }
    let shells: Vec<Handle<Shell>> = by_comp
        .iter()
        .map(|faces| {
            model.shells.push(Shell {
                faces: faces.clone(),
            })
        })
        .collect();
    match positives.len() {
        0 => Err(reject(RejectReason::NoOutwardShell)),
        1 => {
            let outer_c = positives[0];
            // A cavity shell's faces already point into the void (the material is outside it, so
            // the material-on-correct-side reconstruction winds them inward) — measured, no flip.
            let cavities = (0..n)
                .filter(|&c| c != outer_c)
                .map(|c| shells[c])
                .collect();
            let solid = model.push_solid(Solid {
                outer: shells[outer_c],
                cavities,
            });
            model.live_solids.retain(|&s| s != a && s != b);
            Ok(vec![solid])
        }
        _ => {
            // Several material solids. Assign each surviving cavity (an inward component) to the
            // material whose outer shell nests it — the innermost, if materials themselves nest —
            // by the exact point-in-solid `point_in_component`. With no cavity this loop is empty
            // and every material emits cavity-free (the plain sever, unchanged).
            // ★ The rings hand over their walls; nothing here derives one from a name.
            let comp_faces = |c: usize| -> Result<combinatorics::ComponentFaces, BoolError> {
                by_comp_lf[c]
                    .iter()
                    .map(|lf| {
                        let mut rings = vec![lf.loop_nodes.edges(jd, lf.plane_idx)?];
                        for h in &lf.inner {
                            rings.push(h.edges(jd, lf.plane_idx)?);
                        }
                        Ok((lf.plane_idx, rings))
                    })
                    .collect()
            };
            let nodes_of = |c: usize| -> Vec<[usize; 3]> {
                by_comp_lf[c]
                    .iter()
                    .flat_map(|lf| lf.loop_nodes.iter().map(|Node::Seam(t)| *t))
                    .collect()
            };
            let mut cavities_of: std::collections::HashMap<usize, Vec<Handle<Shell>>> =
                positives.iter().map(|&m| (m, Vec::new())).collect();
            for d in (0..n).filter(|c| !positives.contains(c)) {
                // A cavity node that classifies cleanly against *every* material (one shared origin
                // keeps the nesting consistent); its `true` materials nest, so take the innermost.
                let containers = nodes_of(d).iter().find_map(|&x| {
                    let mut cs = Vec::new();
                    for &m in &positives {
                        match comp_faces(m)
                            .and_then(|f| combinatorics::point_in_component(jd, x, &f))
                        {
                            Ok(true) => cs.push(m),
                            Ok(false) => {}
                            Err(_) => return None, // grazed against a material — try next node
                        }
                    }
                    Some(cs)
                });
                let containers = containers.ok_or_else(|| reject(RejectReason::NoClearRay))?;
                let owner = match containers.as_slice() {
                    [] => return Err(reject(RejectReason::CavityNoOwner)),
                    [only] => *only,
                    _ => *containers
                        .iter()
                        .find(|&&c| {
                            // Innermost: inside every other container.
                            containers.iter().all(|&o| {
                                o == c
                                    || nodes_of(c)
                                        .iter()
                                        .find_map(|&x| {
                                            comp_faces(o)
                                                .and_then(|f| {
                                                    combinatorics::point_in_component(jd, x, &f)
                                                })
                                                .ok()
                                        })
                                        .unwrap_or(false)
                            })
                        })
                        .ok_or_else(|| reject(RejectReason::CavityNoOwner))?,
                };
                cavities_of.get_mut(&owner).unwrap().push(shells[d]);
            }
            // Emit each material solid (with its cavities) in a canonical, replay-stable order
            // keyed on geometry, so a downstream op can index the returned Vec deterministically.
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
            let solids: Vec<Handle<Solid>> = order
                .into_iter()
                .map(|oi| {
                    let c = positives[oi];
                    model.push_solid(Solid {
                        outer: shells[c],
                        cavities: cavities_of.remove(&c).unwrap(),
                    })
                })
                .collect();
            model.live_solids.retain(|&s| s != a && s != b);
            Ok(solids)
        }
    }
}

/// An unordered edge key: the two nodes in a fixed order, so `{a,b}` and `{b,a}` collide.
fn norm_edge(a: Node, b: Node) -> (Node, Node) {
    if a <= b { (a, b) } else { (b, a) }
}

/// Merge every group of coplanar, same-facing result faces into one face per connected piece, then
/// dissolve straight-angle vertices — the mandatory post-boolean defeature.
///
/// **Erase the interior, do not stitch the exterior.** When two faces become one region the boundary
/// between them stops being a boundary, so the merge is: take every directed ring edge in the group,
/// drop the ones that appear as an opposed pair (`a→b` together with `b→a`), and re-thread what is
/// left. Nothing has to be spliced, which is what lets one rule cover every way the pieces can meet
/// — sharing an edge, a hole filled exactly by a neighbour, a hole filled by *several* neighbours,
/// and any chain of those (one erase settles them all at once). The earlier version stitched loops
/// with `splice_along` and so had to special-case "exactly two hole-free faces across one edge",
/// leaving `// holed — deferred` for the rest; a flush tool cap then stayed two faces forever.
///
/// A group is one plane class, one outward direction, one `flip` — mixing any of those would fold
/// material the wrong way — split further into **edge-connected components**, because faces that
/// merely lie on the same plane without touching must each survive on their own.
///
/// Reused, not reinvented: `canon` for coplanarity, `n_out.dot > 0` for facing,
/// [`combinatorics::loop_winding`] to tell an outer ring from a hole, [`combinatorics::point_in_ring`] to give
/// each hole its owner — the same two exact predicates `arrangement::nest_cells` uses for the same
/// question, neither of which reads a coordinate.
///
/// Rejects rather than guesses: a directed edge appearing twice the same way (two faces claiming the
/// same side), an undirected edge on three or more rings (non-manifold in the plane), or a node with
/// two outgoing edges after erasure (pieces meeting at a single point, where the cycle is not
/// unique).
pub(crate) fn unify_coplanar_faces(
    faces: Vec<LocalFace>,
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Vec<LocalFace>, BoolError> {
    let n = faces.len();
    // One plane class, one flip.
    //
    // There used to be an outward-direction component here, and a `canon` lookup beside it. Neither
    // separated anything: `plane_idx` names a plane, so its class is itself, and the component's dot
    // product was `|n|² > 0` — identically true. It was redundant with `flip` besides, which is
    // already in the key: `assemble_fuse_cut` derives a result face's `Orientation` from exactly
    // `(the plane's frame, flip)`, so two faces in one group orient the same way by construction.
    let group_key = |lf: &LocalFace| -> (usize, bool) { (lf.plane_idx, lf.flip) };

    // Edge-connected components within a group.
    let mut comp: Vec<usize> = (0..n).collect();
    let mut carriers: HashMap<(Node, Node), Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        // Holes count here too: a tool cap sitting flush inside another face touches it only
        // along that hole, so leaving `inner` out would put the two in different components and
        // nothing would merge at all.
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            for (a, b) in ring_edges(ring) {
                carriers.entry(norm_edge(a, b)).or_default().push(fi);
            }
        }
    }
    for fs in carriers.values() {
        for w in fs.windows(2) {
            if group_key(&faces[w[0]]) == group_key(&faces[w[1]]) {
                let (ri, rj) = (uf_find(&mut comp, w[0]), uf_find(&mut comp, w[1]));
                if ri != rj {
                    comp[ri] = rj;
                }
            }
        }
    }

    // Group members by component, in face order so the result is replay-stable.
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); n];
    for fi in 0..n {
        let r = uf_find(&mut comp, fi);
        members[r].push(fi);
    }

    let mut merged: Vec<LocalFace> = Vec::new();
    let mut kept: Vec<Option<LocalFace>> = faces.into_iter().map(Some).collect();
    for mem in &members {
        if mem.len() < 2 {
            continue; // nothing to merge; the face (if any) is emitted as-is below
        }
        let group: Vec<&LocalFace> = mem
            .iter()
            .map(|&fi| kept[fi].as_ref().expect("member present"))
            .collect();
        let rings = merge_component(&group, jd)?;
        let (plane_idx, flip) = (group[0].plane_idx, group[0].flip);
        merged.extend(rings.into_iter().map(|(outer, inner)| LocalFace {
            plane_idx,
            loop_nodes: outer,
            inner,
            flip,
        }));
        for &fi in mem {
            kept[fi] = None;
        }
    }
    let mut out: Vec<LocalFace> = kept.into_iter().flatten().collect();
    out.extend(merged);
    dissolve_straight_angles(&mut out);
    Ok(out)
}

/// A ring's directed edges, `i → i+1` around.
fn ring_edges(ring: &[Node]) -> impl Iterator<Item = (Node, Node)> + '_ {
    (0..ring.len()).map(move |i| (ring[i], ring[(i + 1) % ring.len()]))
}

/// A ring's directed edges **with the wall each rides**.
fn ring_edges_walled(ring: &Ring) -> impl Iterator<Item = ((Node, Node), usize)> + '_ {
    (0..ring.len()).map(move |i| {
        (
            (ring.nodes[i], ring.nodes[(i + 1) % ring.len()]),
            ring.walls[i],
        )
    })
}

/// The `[usize; 3]` form a ring's nodes carry, for the exact predicates.
fn seam_ring(ring: &[Node]) -> Vec<[usize; 3]> {
    ring.iter().map(|Node::Seam(t)| *t).collect()
}

/// An outer ring with the holes that belong to it — what one merged region looks like before it
/// becomes a `LocalFace`.
type RegionRings = (Ring, Vec<Ring>);

/// One edge-connected group → its faces after erasing the interior boundary: each outer ring with
/// the holes that belong to it.
fn merge_component(
    group: &[&LocalFace],
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Vec<RegionRings>, BoolError> {
    // 1. Collect directed edges **with their walls**. A repeat in the same direction means two
    //    faces claim the same side.
    // One map, `(count, wall)` — the wall rides along rather than in a second table.
    let mut dirs: HashMap<(Node, Node), (usize, usize)> = HashMap::new();
    for lf in group {
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            for (e, wall) in ring_edges_walled(ring) {
                let slot = dirs.entry(e).or_insert((0, wall));
                slot.0 += 1;
            }
        }
    }
    if dirs.values().any(|&(c, _)| c > 1) {
        return Err(reject(RejectReason::CoplanarMerge));
    }
    // 2. An edge carried in both directions is interior — it separates nothing. Anything carried
    //    three or more times (either direction) is non-manifold in the plane.
    let mut undirected: HashMap<(Node, Node), usize> = HashMap::new();
    for &(a, b) in dirs.keys() {
        *undirected.entry(norm_edge(a, b)).or_insert(0) += 1;
    }
    if undirected.values().any(|&c| c > 2) {
        return Err(reject(RejectReason::CoplanarMerge));
    }
    // 3. Re-thread what survives. Two outgoing edges at one node means the pieces meet at a point
    //    and the cycles are not determined.
    let mut next: HashMap<Node, Node> = HashMap::new();
    for &(a, b) in dirs.keys() {
        if dirs.contains_key(&(b, a)) {
            continue; // interior
        }
        if next.insert(a, b).is_some() {
            return Err(reject(RejectReason::CoplanarMerge));
        }
    }
    let mut starts: Vec<Node> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut seen: HashSet<Node> = HashSet::new();
    let mut cycles: Vec<Ring> = Vec::new();
    for start in starts {
        if seen.contains(&start) {
            continue;
        }
        // ★ The threaded cycle keeps each surviving edge's own wall, so the merged ring never has
        // to have it read back out of the node names.
        let mut nodes = vec![start];
        let mut walls = vec![dirs[&(start, next[&start])].1];
        seen.insert(start);
        let mut cur = next[&start];
        while cur != start {
            if !seen.insert(cur) {
                return Err(reject(RejectReason::CoplanarMerge)); // walk re-entered another cycle
            }
            let nx = *next
                .get(&cur)
                .ok_or_else(|| reject(RejectReason::CoplanarMerge))?;
            nodes.push(cur);
            walls.push(dirs[&(cur, nx)].1);
            cur = nx;
        }
        if nodes.len() < 3 {
            return Err(reject(RejectReason::CoplanarMerge));
        }
        cycles.push(Ring::new(nodes, walls));
    }
    if cycles.is_empty() {
        return Err(reject(RejectReason::CoplanarMerge)); // everything erased: not a region
    }
    // 4. Winding tells an outer ring from a hole; the plane is the frame both are read in. (This
    // used to canon the index first — `plane_idx` names a plane now, so there is nothing to fold.)
    let wc = group[0].plane_idx;
    let mut outers: Vec<Ring> = Vec::new();
    let mut holes: Vec<Ring> = Vec::new();
    for cyc in cycles {
        // ★ Built from the walls the cycle carries, not from the node names. This is where a
        // four-plane concurrency used to break the merge: a canonical name need not mention the
        // plane its edge rides, and two names can share nothing but `wc`.
        let ring = cyc.edges(jd, wc)?;
        match combinatorics::loop_winding(jd, wc, &ring)? {
            1 => outers.push(cyc),
            -1 => holes.push(cyc),
            _ => return Err(reject(RejectReason::CoplanarMerge)),
        }
    }
    // 5. Each hole belongs to the outer ring that contains it — the same question `nest_cells` asks
    //    of the arrangement's cells, answered by the same predicate.
    let mut faces: Vec<RegionRings> = outers.into_iter().map(|o| (o, Vec::new())).collect();
    for hole in holes {
        let Node::Seam(probe) = hole.nodes[0];
        let mut owner = None;
        for (i, (outer, _)) in faces.iter().enumerate() {
            let ring = outer.edges(jd, wc)?;
            if combinatorics::point_in_ring(jd, wc, probe, &ring)? {
                if owner.is_some() {
                    return Err(reject(RejectReason::CoplanarMerge)); // nested deeper than this brick names
                }
                owner = Some(i);
            }
        }
        faces[owner.ok_or_else(|| reject(RejectReason::CoplanarMerge))?]
            .1
            .push(hole);
    }
    Ok(faces)
}

/// Drop straight-angle vertices: a node whose only two neighbours across **all** rings continue
/// along the same line.
///
/// ★ **The test is the two incident edges' walls**, which the rings carry. It used to be
/// combinatorial on the names — "the node is on a line iff some pair of its planes is shared by
/// both neighbours" — which is sound only while every vertex lies on exactly three planes, the same
/// assumption that broke the merge under a four-plane concurrency. The walls say it outright.
///
/// Dropping from every incident ring at once keeps a vertex that is a real corner somewhere
/// (degree > 2), which is what stops a T-junction from opening.
fn dissolve_straight_angles(out: &mut [LocalFace]) {
    // ★ The walls are per `(node, face plane)`. Globally they cannot be: the two result faces
    // that share a 3D edge each ride *the other's* plane as their wall, so a node in the middle of
    // that edge always sees two walls overall and would never dissolve. On each face separately it
    // sees one, which is exactly "the ring runs straight through here".
    // Per node: its neighbours, and whether any face bends there. `first_wall` remembers one wall
    // per `(node, face)` and `bent` records the first disagreement — flat maps, because a nested
    // one per node costs more than the whole pass is worth (measured: 4% of a fold).
    let mut nbrs: HashMap<Node, HashSet<Node>> = HashMap::new();
    let mut first_wall: HashMap<(Node, usize), usize> = HashMap::new();
    let mut bent: HashSet<Node> = HashSet::new();
    for lf in out.iter() {
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            let k = ring.len();
            for i in 0..k {
                let (a, b) = (ring.nodes[i], ring.nodes[(i + 1) % k]);
                nbrs.entry(a).or_default().insert(b);
                nbrs.entry(b).or_default().insert(a);
                for nd in [a, b] {
                    match first_wall.entry((nd, lf.plane_idx)) {
                        std::collections::hash_map::Entry::Occupied(e) => {
                            if *e.get() != ring.walls[i] {
                                bent.insert(nd);
                            }
                        }
                        std::collections::hash_map::Entry::Vacant(e) => {
                            e.insert(ring.walls[i]);
                        }
                    }
                }
            }
        }
    }
    let mut drop: HashSet<Node> = HashSet::new();
    for (&node, ns) in &nbrs {
        // Exactly two neighbours, and on **every** face it appears in the two edges ride one wall.
        if ns.len() == 2 && !bent.contains(&node) {
            drop.insert(node);
        }
    }
    if drop.is_empty() {
        return;
    }
    for lf in out.iter_mut() {
        for ring in std::iter::once(&mut lf.loop_nodes).chain(lf.inner.iter_mut()) {
            if !ring.nodes.iter().any(|nd| drop.contains(nd)) {
                continue; // untouched rings keep their allocation
            }
            let k = ring.len();
            let mut nodes = Vec::with_capacity(k);
            let mut walls = Vec::with_capacity(k);
            for i in 0..k {
                if drop.contains(&ring.nodes[i]) {
                    continue; // its two edges are one; keep the wall of the one that survives
                }
                nodes.push(ring.nodes[i]);
                walls.push(ring.walls[i]);
            }
            *ring = Ring::new(nodes, walls);
        }
    }
}
