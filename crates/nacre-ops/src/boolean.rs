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
use nacre_cip::Decision;
use nacre_cip::predicate::{Evidence, Notes, Site};
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
    // ★ **A rejected boolean leaves the live set as it found it** — every reject, not just the
    // topology one below. The engine retires the operands once the result is accepted, so a reject
    // raised after that point would otherwise hand back a model whose inputs had vanished. Restoring
    // here makes the rule hold no matter where inside the engine a reject is raised or added later.
    // (Orphaned result cells stay in the append-only arena, unreachable, as any superseded solid's
    // do — the arena is not restored and is not meant to be.)
    let snapshot = model.live_solids.clone();
    let (result, notes, _class_of) = match crate::arrangement::boolean(model, kind, a, b) {
        Ok(v) => v,
        Err(e) => {
            model.live_solids = snapshot;
            return Err(surfacing(e));
        }
    };
    // Topological self-check on the assembled result (DNA: never return a malformed solid). Reject
    // rather than return. Valid results always pass, so this never false-rejects; the traversal
    // reads the topology stores directly (no adjacency rebuild, no coordinates).
    if let Some((t, at)) = check_result_topology(model, &result) {
        model.live_solids = snapshot;
        return Err(surfacing(match at {
            Some(w) => crate::reject_at(t, w),
            None => reject(t),
        }));
    }
    Ok((result, BoolReport::of(&notes)))
}

/// Record that `e` is leaving the kernel — the *surfaced* column of [`crate::reject_census`].
///
/// Sits at the public entry points rather than at [`reject`] because those are two different
/// populations: guards are raised and swallowed (a retry tries the ring's next node), and only
/// what comes back here is something a caller ever sees. `InputNotLive` is not a
/// [`RejectReason`] — it is a caller mistake, not a guard — so it stays out of the census.
fn surfacing(e: BoolError) -> BoolError {
    if let BoolError::Rejected { reason, .. } = e {
        crate::reject_census::surfaced(reason);
    }
    e
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
    // The live-set restore is [`boolean`]'s rule, held here too: a reject leaves the live model
    // alone whichever entry point raised it.
    let snapshot = model.live_solids.clone();
    let (result, notes, class_of) = match crate::arrangement::boolean(model, kind, a, b) {
        Ok(v) => v,
        Err(e) => {
            model.live_solids = snapshot;
            return Err(surfacing(e));
        }
    };
    if let Some((t, at)) = check_result_topology(model, &result) {
        model.live_solids = snapshot;
        return Err(surfacing(match at {
            Some(w) => crate::reject_at(t, w),
            None => reject(t),
        }));
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
fn check_result_topology(
    model: &Model,
    solids: &[Handle<Solid>],
) -> Option<(RejectReason, Option<crate::RejectWhere>)> {
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
        // The returned list is sorted by handle index, so the first entry is a deterministic
        // witness (one of possibly several pinch vertices).
        if let Some(&vh) = nacre_topo::nonmanifold_vertices(&vertex_edges, &edge_uses).first() {
            return Some((
                RejectReason::NonManifoldVertex,
                Some(crate::RejectWhere::Point(model.vertex_point(vh))),
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
/// ★★ **Touching is not joining.** Two bodies that meet along a line share that line's nodes, and
/// used to come back as one component and so one pinched pseudo-solid, which the edge-use guard
/// then had to reject. They share no *manifold* edge — four faces use the contact line, not two —
/// so they land in two components and come back as the two bodies they are. Where a body touches
/// **itself** the material still runs around the contact, every edge of that path is used twice,
/// and the component stays one: the guard fires exactly where no pair of solids exists.
///
/// This is the 3D reading of a rule this file already applies one dimension down:
/// [`unify_coplanar_faces`] splits a coplanar group "into **edge-connected components**, because
/// faces that merely lie on the same plane without touching must each survive on their own".
///
/// An enclosed void is its own component, as before — its boundary shares no edge with the outer.
fn face_components(faces: &[LocalFace]) -> (Vec<usize>, usize) {
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    // Which faces use each ring edge. A count other than two is not a contact between neighbours:
    // one is a dangling edge and more is a pinch, and both are the edge-use guard's to name.
    let mut users: HashMap<(Node, Node), Vec<usize>> = HashMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for ring in rings_of(lf) {
            let k = ring.len();
            for t in 0..k {
                users
                    .entry(norm_edge(ring[t], ring[(t + 1) % k]))
                    .or_default()
                    .push(i);
            }
        }
    }
    for us in users.values() {
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

/// **Which faces make one output solid — decided before any topology exists.**
///
/// A result solid is a *material* component plus the cavity components nested in it, and that
/// grouping is the unit [`reconstruct`] mints vertex and edge handles per: two faces in one group
/// may share a handle, two faces in different groups never do. Every question asked here is asked
/// of `LocalFace`/`Node`/`jd` alone — none of it reads the arena — which is what lets it run first.
///
/// ★★ **The group, not the component, is the right unit.** A cavity that touches its host's outer
/// shell is one solid with a pinch, and a pinch is counted on *handles* ([`check_result_topology`],
/// `validate`): split the handles there and the defect stops being visible while the model keeps
/// its zero-thickness material. Grouping keeps a cavity with its host, so that case still reaches
/// the reject it deserves.
struct Grouping {
    /// Connected-component label per face, dense `0..n`.
    labels: Vec<usize>,
    /// Component count.
    n: usize,
    /// The material components (even nesting depth), ascending.
    positives: Vec<usize>,
    /// Per material component, the components its solid is made of: itself first, then its
    /// cavities in ascending order.
    comps_of: HashMap<usize, Vec<usize>>,
    /// Group index per face — the unit handles are minted per. Numbered by position in
    /// `positives`; the canonical *output* order is decided later, from geometry (`comp_key`).
    group_of: Vec<usize>,
}

/// Partition the faces into connected components, then into output solids. One component is the
/// whole result; several mean either an enclosed void (a cavity — an inward-oriented shell) or a
/// severed operand (two or more material-enclosing shells).
///
/// ★ **Deciding early is not the same as *rejecting* early.** [`reconstruct`] holds this `Result`
/// and raises it exactly where the old code did, so a boolean that declines here leaves the arena
/// it left before — the quantity `replay::a_late_reject_is_not_index_neutral` measures.
fn group_faces(jd: &Judge<'_, WorkingPlane>, faces: &[LocalFace]) -> Result<Grouping, BoolError> {
    let (labels, n) = face_components(faces);
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
                .map(|lf| {
                    // Circle bounds carry no ring edges — their component adjacency is the
                    // rim rule's (C4b); a node-level containment probe reads polygons only.
                    let rings = lf
                        .poly_rings()
                        .map(|r| r.edges(jd, lf.plane_idx))
                        .collect::<Result<Vec<_>, BoolError>>()?;
                    Ok((lf.plane_idx, rings))
                })
                .collect::<Result<_, BoolError>>()
        })
        .collect::<Result<_, BoolError>>()?;
    let nodes_of = |c: usize| -> Vec<[usize; 3]> {
        by_comp_lf[c]
            .iter()
            .flat_map(|lf| {
                lf.poly_rings()
                    .flat_map(|r| r.iter().map(|Node::Seam(t)| *t))
            })
            .collect()
    };
    // Try `f` at each node in turn: the first node that **decides** wins, a node that abstains
    // (`Ok(None)` — it grazed a boundary) is passed over for the next, and a failed judgement
    // (`Err`) propagates immediately. The last part is the point of the shape: an abstention has
    // other nodes as its remedy, a failed judgement does not — the old `Err(_) => try the next
    // node` arms retried both, so a real cause could masquerade as `NoClearRay` once every node
    // hit it. `Ok(None)` here means every node abstained; that being a reject is the *caller's*
    // proposition to raise.
    fn first_deciding<T>(
        nodes: &[[usize; 3]],
        mut f: impl FnMut([usize; 3]) -> Result<Option<T>, BoolError>,
    ) -> Result<Option<T>, BoolError> {
        for &x in nodes {
            if let Some(v) = f(x)? {
                return Ok(Some(v));
            }
        }
        Ok(None)
    }
    // ★★ **Material or void is a question about nesting, not about normals.**
    //
    // It used to be answered at the component's lexicographically-minimal vertex `v*`: outward iff
    // *some* face there has an outward normal with `n_x < 0`. That existential is a shortcut for
    // "the face with the largest `|n_x|` faces −x", and the shortcut is only equivalent while the
    // normals are **axis-aligned** — the old `is_shell_outward`'s own doc said so. A slanted sketch
    // breaks it with no rotation in sight: a wedge void cut inside a box came back as its own
    // *material* solid of **negative volume**, with the box unchanged beside it. A wrong model,
    // and `validate` had nothing to say about it (measured, `docs/dev-log.md`).
    //
    // So ask the question the kernel already answers one dimension down. `sketch::from_rings`
    // decides a ring by **containment depth — even is material, odd is a hole** (design.md: "채우기
    // 규칙 파라미터는 두지 않는다 — 짝수-홀수가 유일한 규칙이다"), and the cavity-owner search
    // below already picks the *innermost* container, the other half of that same rule. This is the
    // 3D reading of it, with `point_in_component` where the 2D one uses `point_in_ring`.
    //
    // Nothing here reads an orientation, so a shell that winds either way is labelled the same —
    // which is the point: the winding is what the old test was trying, and failing, to recover.
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
        let depth = first_deciding(&nodes_of(c), |x| {
            let mut d = 0usize;
            for other in (0..n).filter(|&o| o != c) {
                match combinatorics::point_in_component(jd, x, &comp_faces[other])? {
                    Some(true) => d += 1,
                    Some(false) => {}
                    None => return Ok(None), // grazed — this node abstains
                }
            }
            Ok(Some(d))
        })?
        .ok_or_else(|| reject(RejectReason::NoClearRay))?;
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
            // `comp_faces`/`nodes_of` are the same pair the material/void label above uses.
            for d in (0..n).filter(|c| !positives.contains(c)) {
                // A cavity node that classifies cleanly against *every* material (one shared origin
                // keeps the nesting consistent); its `true` materials nest, so take the innermost.
                let containers = first_deciding(&nodes_of(d), |x| {
                    let mut cs = Vec::new();
                    for &m in &positives {
                        match combinatorics::point_in_component(jd, x, &comp_faces[m])? {
                            Some(true) => cs.push(m),
                            Some(false) => {}
                            None => return Ok(None), // grazed against a material — abstain
                        }
                    }
                    Ok(Some(cs))
                })?
                .ok_or_else(|| reject(RejectReason::NoClearRay))?;
                let owner = match containers.as_slice() {
                    [] => return Err(reject(RejectReason::CavityNoOwner)),
                    [only] => *only,
                    _ => {
                        // Innermost: inside every other container. The pairwise facts are
                        // decided **before** the search, because the search closures speak
                        // bool and a failed judgement must propagate, not vanish into "not
                        // inside". A pair where every node abstains stays `false` — no
                        // evidence places `c` inside `o`, the same answer the old `.ok()`
                        // retry produced. (Nested containers are a rare population; deciding
                        // the few pairs up front costs nothing measurable.)
                        let mut inside = HashMap::new();
                        for &c in &containers {
                            for &o in &containers {
                                if o == c {
                                    continue;
                                }
                                let v = first_deciding(&nodes_of(c), |x| {
                                    combinatorics::point_in_component(jd, x, &comp_faces[o])
                                })?
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
fn comp_key(model: &Model, faces: &[Handle<Face>]) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = faces
        .iter()
        .flat_map(|&fh| model.faces.get(fh).outer.half_edges.iter().copied())
        .map(|he| model.vertex_point(he_start(model, he)).as_array())
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite vertex coordinates"));
    pts
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

/// One boundary of a result face (M6-2a): a polygon of seam nodes, or a **full circle** of a
/// cylinder class. A circle has no nodes and no walls; it is assembled through the rim
/// machinery (`push_edge([lateral, plane], [v, v])` + an `OnSeam` vertex), never through the
/// seam-vertex table.
#[derive(Clone, Debug)]
pub(crate) enum Bound {
    Ring(Ring),
    Circle { cyl: usize },
}

impl Bound {
    /// The polygon ring, `None` for a circle — the node-walking consumers' filter.
    pub(crate) fn ring(&self) -> Option<&Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } => None,
        }
    }

    pub(crate) fn ring_mut(&mut self) -> Option<&mut Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } => None,
        }
    }

    /// The polygon ring, asserted — for consumers whose population cannot carry circles (and
    /// tests). Panics on a circle with the caller's location.
    #[cfg_attr(not(test), allow(dead_code))]
    #[track_caller]
    pub(crate) fn expect_ring(&self) -> &Ring {
        match self {
            Bound::Ring(r) => r,
            Bound::Circle { cyl } => panic!("a polygon-only path got a circle bound (cyl {cyl})"),
        }
    }
}

/// A reconstructed result face: which combined plane it is on, its boundaries, and whether to
/// flip it (cut's inside-A B-pieces).
#[derive(Clone, Debug)]
pub(crate) struct LocalFace {
    pub(crate) plane_idx: usize,
    pub(crate) outer: Bound,
    /// Hole boundaries, each already wound so the kept material stays on its left
    /// about the face's outward normal.
    pub(crate) inner: Vec<Bound>,
    pub(crate) flip: bool,
}

impl LocalFace {
    /// Every **polygon** ring of the face, outer first. Circle bounds carry no nodes and are
    /// deliberately absent — node-level consumers (the seam table, edge scans, self-touch)
    /// have nothing to read off them.
    pub(crate) fn poly_rings(&self) -> impl Iterator<Item = &Ring> {
        std::iter::once(&self.outer)
            .chain(self.inner.iter())
            .filter_map(Bound::ring)
    }

    pub(crate) fn poly_rings_mut(&mut self) -> impl Iterator<Item = &mut Ring> {
        std::iter::once(&mut self.outer)
            .chain(self.inner.iter_mut())
            .filter_map(Bound::ring_mut)
    }

    /// Whether any boundary is a circle — the unify partition reads this (circle-bounded
    /// faces are never merge candidates).
    pub(crate) fn has_circle(&self) -> bool {
        std::iter::once(&self.outer)
            .chain(self.inner.iter())
            .any(|b| matches!(b, Bound::Circle { .. }))
    }
}

/// Push the reconstructed result and supersede the inputs.
///
/// ★ **The operands retire in exactly one place, and it is after the result is accepted.** The
/// retire used to be copied into each arm that produced solids, which meant every future reject had
/// to know whether it stood before or after its arm's copy. Hoisting it here makes "a reject leaves
/// the model alone" true of the *shape* of this function rather than of a fact about where the
/// rejects happen to sit today. [`reconstruct`] does the work and never touches the live set.
pub(crate) fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let out = reconstruct(model, jd, seam, faces)?;
    model.live_solids.retain(|&s| s != a && s != b);
    Ok(out)
}

/// **A solid whose surface touches itself is not a solid.**
///
/// The proposition: *an edge of this solid lies in the interior of one of this same solid's faces.*
/// An embedded boundary cannot do that — the surface would occupy the same points twice — so a
/// result that does is refused rather than returned. The unit is the **result solid**, outer shell
/// and cavities together: a wedge cut whose tip lands on the far wall touches its own cavity, and
/// asking the question per component would miss exactly that.
///
/// **Three sieves, cheapest first**, because the exact question is the expensive one:
///
/// 1. ★ **Which planes** — both endpoints must lie on the face's plane, and a vertex's plane triple
///    *is* the set of planes it lies on, so the candidates for an edge are `triple(u) ∩ triple(v)`,
///    usually two. Exact, no coordinates, a hash lookup per edge. Scanning every face for every
///    edge instead costs ~40% on the 80-fold star (measured).
/// 2. ★ **Which faces on that plane** — one plane can carry *eighty* faces in a folded star, so the
///    lookup above is not the end of it. A box per face rejects 99.99% of what survives step 1
///    (measured: 2,538,703 candidate pairs down to 203), and on the booleans that cost the most it
///    rejects all of them.
/// 3. The exact test, on what is left: [`combinatorics::segment_meets_face`] — **the open segment**,
///    not its endpoints. The endpoints alone miss the shape that motivated this: a cut taken all the
///    way *through* leaves a contact line whose ends sit on the face's ring while its middle crosses
///    the interior, and a point test reads both ends as "on the boundary, not inside" and passes a
///    body that cannot exist (measured: `Ok(n=1)`, volume exact, `validate` silent).
///
/// ★ **One judge, not two.** The endpoint test was here first and is *subsumed* — an endpoint
/// strictly inside puts that end in an inside stretch of the line, so the segment test answers it
/// too. Measured across the corpus: 1,696 candidate pairs, the two verdicts agreed on every one,
/// and no pair was undecidable. Keeping both would have left the same question answered in two
/// places.
///
/// ★ **The boxes are inflated by the vertices' own tolerance and so can only over-keep.** They are
/// built from realized `f64` coordinates, which are rounded; a box used to *reject* a candidate
/// before an exact test must therefore be conservative, or a real self-contact is dropped in
/// silence. `SeamVertex::tol` is the measured bound on that vertex's realization — both the box and
/// the query point are widened by it. The case this protects is not hypothetical: the wedge that
/// motivated this check touches at exactly `x = 1`, which is a face of its own box.
///
/// The rings a face is trimmed by are built **lazily and once per face**: `Ring::edges` spends an
/// exact predicate per node, and building them for every result face is the same unconditional cost
/// that put 12x on the star when the void label did it.
fn self_touch_reject(
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    groups: &[Vec<usize>],
    by_comp_lf: &[Vec<&LocalFace>],
) -> Result<(), BoolError> {
    let pt: HashMap<[usize; 3], (Point3, f64)> = seam
        .iter()
        .map(|sv| (sv.triple, (sv.point, sv.tol)))
        .collect();
    for g in groups {
        let faces: Vec<&LocalFace> = g
            .iter()
            .flat_map(|&c| by_comp_lf[c].iter().copied())
            .collect();
        // One pass: which faces sit on each plane, which two faces own each edge, and a box per
        // face already widened by its own vertices' tolerances.
        let mut by_plane: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut owners: HashMap<[[usize; 3]; 2], Vec<usize>> = HashMap::new();
        let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(faces.len());
        for (j, lf) in faces.iter().enumerate() {
            by_plane.entry(lf.plane_idx).or_default().push(j);
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for r in rings_of(lf) {
                let k = r.nodes.len();
                for i in 0..k {
                    let (Node::Seam(u), Node::Seam(v)) = (r.nodes[i], r.nodes[(i + 1) % k]);
                    owners
                        .entry(if u < v { [u, v] } else { [v, u] })
                        .or_default()
                        .push(j);
                    if let Some(&(p, tol)) = pt.get(&u) {
                        let a = p.as_array();
                        for k in 0..3 {
                            lo[k] = lo[k].min(a[k] - tol);
                            hi[k] = hi[k].max(a[k] + tol);
                        }
                    }
                }
            }
            boxes.push((lo, hi));
        }
        let mut rings_of_face: HashMap<usize, Vec<Vec<combinatorics::RingEdge>>> = HashMap::new();
        for (&[u, v], own) in &owners {
            for q in u.iter().copied().filter(|x| v.contains(x)) {
                let Some(js) = by_plane.get(&q) else { continue };
                for &j in js.iter().filter(|j| !own.contains(j)) {
                    let (lo, hi) = boxes[j];
                    let in_box = [u, v].iter().all(|t| {
                        pt.get(t).is_some_and(|&(p, tol)| {
                            let a = p.as_array();
                            (0..3).all(|k| a[k] >= lo[k] - tol && a[k] <= hi[k] + tol)
                        })
                    });
                    if !in_box {
                        continue;
                    }
                    if let std::collections::hash_map::Entry::Vacant(slot) = rings_of_face.entry(j)
                    {
                        slot.insert(
                            rings_of(faces[j])
                                .map(|r| r.edges(jd, q))
                                .collect::<Result<_, BoolError>>()?,
                        );
                    }
                    // ★ The line is `q ∩ w`, and `w` comes from **the edge's own two faces**, not
                    // from the endpoint names and not from the ring's recorded wall. Both of those
                    // fail, measured: an edge can lie on a pencil of three planes and the ring is
                    // free to call it by any of them — including `q` itself — and where four planes
                    // concur the endpoints' canonical triples need not share a second plane at all
                    // (31 pairs in the rotation sweep). An edge's two faces always give one, because
                    // that is what the edge *is*.
                    let Some(w2) = own.iter().map(|&o| faces[o].plane_idx).find(|&x| x != q) else {
                        continue; // an edge whose faces are both on `q` is not an edge
                    };
                    if combinatorics::segment_meets_face(jd, q, w2, u, v, &rings_of_face[&j])? {
                        // The offending edge itself, as the witness: both endpoints are in `pt`
                        // (the `in_box` guard above already looked them up).
                        return Err(crate::reject_at(
                            RejectReason::SelfTouchingResult,
                            crate::RejectWhere::Segment([pt[&u].0, pt[&v].0]),
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

/// A face's rings, outer first.
fn rings_of(lf: &LocalFace) -> impl Iterator<Item = &Ring> {
    lf.poly_rings()
}

/// Rebuild the result solids from the arrangement's faces. Pushes into the arena; it does not take
/// the operands at all, which is the point — the live set is its caller's to move.
fn reconstruct(
    model: &mut Model,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes = jd.planes;
    // No faces means no result — `Common` of two solids that miss each other, `Cut` of a box that
    // is wholly inside what cuts it. That is an answer, not a failure: a solid is bounded by faces,
    // so a non-empty result cannot have none. The inputs are still consumed, exactly as they are on
    // any other successful boolean — this returns `Ok`, so the caller's retire runs.
    if faces.is_empty() {
        return Ok(Vec::new());
    }
    // ★ **Which faces make one solid, decided before a single handle exists** — see [`Grouping`].
    // Everything derived below (an edge's far plane, a vertex's defining triple, the handles
    // themselves) is scoped to one group, so no result solid can be named by — or share a handle
    // with — a solid it merely touches.
    //
    // ★★ **Held, not raised.** A failed grouping is reported further down, where the old code
    // reported it, and until then every face is one group — which is exactly the keying this
    // function used before groups existed. So a boolean that declines pushes the arena cells it
    // always did (`replay::a_late_reject_is_not_index_neutral` measures that).
    let grouping = group_faces(jd, faces);
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
    // bound that solid — and the whole-result derivation used to fill the gap by borrowing the
    // *other* body's plane, which is exactly the defect scoping the derivation closes. So each
    // solid now drops the nodes that are straight runs **for it**.
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
    // in-plane edges are lobe-to-lobe, its touch edge carries four faces — measured, the
    // Stage-A probe, dev-log 2026-08-17). A held grouping error stays held (raised further
    // down, where the old code raised it); the self-touch question is only askable of a
    // grouping that answered.
    if let Ok(g) = &grouping {
        let body_comps: Vec<Vec<usize>> =
            g.positives.iter().map(|m| g.comps_of[m].clone()).collect();
        let mut by_comp_lf: Vec<Vec<&LocalFace>> = vec![Vec::new(); g.n];
        for (i, lf) in faces.iter().enumerate() {
            by_comp_lf[g.labels[i]].push(lf);
        }
        self_touch_reject(jd, seam, &body_comps, &by_comp_lf)?;
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
    // So the triple is derived the way this function already derives an edge's carriers
    // (`pair_surfs` below: "read off the whole result, not guessed from one side") and the way
    // `component_is_outward_tol` derives a corner's: **one incident face, plus the far planes of
    // its two edges at that corner.** Those three are result faces by construction, so
    // `defs_are_remappable` holds by construction — and their meet is exactly this vertex, since
    // the two edge lines through it are distinct (checked, not assumed — see the corner guards).
    let mut edge_faces: HashMap<(usize, (Node, Node)), Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        for ring in lf.poly_rings() {
            let k = ring.len();
            for t in 0..k {
                edge_faces
                    .entry((group_of[fi], norm_edge(ring[t], ring[(t + 1) % k])))
                    .or_default()
                    .push(lf.plane_idx);
            }
        }
    }
    // The plane on the other side of an edge, **within this solid**. `None` rather than an error:
    // a corner this cannot
    // resolve is one to skip, and the edge-use guard further down is what judges the face set.
    let far_plane = |g: usize, a: Node, b: Node, own: usize| -> Option<usize> {
        let mut others = edge_faces
            .get(&(g, norm_edge(a, b)))?
            .iter()
            .copied()
            .filter(|&x| x != own);
        let o = others.next()?;
        others.all(|x| x == o).then_some(o)
    };
    let mut def_triple: HashMap<(usize, Node), [usize; 3]> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        for ring in lf.poly_rings() {
            let k = ring.len();
            for t in 0..k {
                let node = ring[t];
                if def_triple.contains_key(&(g, node)) {
                    continue;
                }
                let (Some(prev), Some(next)) = (
                    far_plane(g, ring[(t + k - 1) % k], node, lf.plane_idx),
                    far_plane(g, node, ring[(t + 1) % k], lf.plane_idx),
                ) else {
                    continue;
                };
                // ★ A straight corner is **skipped, not refused**. At a four-plane vertex a face's
                // ring can run straight through the point — both its edges on one line — and then
                // this face has no triple to give, while another face meeting the same vertex
                // does. (`component_is_outward_tol` rejects here instead, and is right to: it
                // asks about one lex-minimal corner, not about every vertex.)
                if prev == next {
                    continue;
                }
                // ★ …and neither is a corner whose two edges ride **different planes that carry
                // one line**. `prev != next` does not rule that out, and three planes through a
                // line name no point — the def would then be a triple that defines nothing. The
                // same determinant the arrangement uses to spot it (`Aliases::record`: "shares a
                // line: names no point") answers here. Measured: it fires nowhere in the suite,
                // which is why the check is here rather than trusted — the argument above claims
                // the meet *is* this vertex, and this is what makes that true by construction.
                if jd.plane_pair_dir_sign(lf.plane_idx, prev, next) == 0 {
                    continue;
                }
                let mut tri = [lf.plane_idx, prev, next];
                tri.sort_unstable();
                def_triple.insert((g, node), tri);
            }
        }
    }

    // Vertices (deterministic: first appearance across faces in order).
    let mut vh: HashMap<(usize, Node), Handle<Vertex>> = HashMap::new();
    let mut node_handle =
        |model: &mut Model, g: usize, node: Node| -> Result<Handle<Vertex>, BoolError> {
            if let Some(&h) = vh.get(&(g, node)) {
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
                    // A vertex that is a corner of no face at all has no name in the result's own
                    // planes — a degeneracy, and the honest answer is the one this reason already
                    // carries ("a corner with no turn").
                    let tri = def_triple
                        .get(&(g, node))
                        .copied()
                        .ok_or_else(|| reject(RejectReason::StraightAngle))?;
                    let def = VertexDef::ThreePlane([
                        planes[tri[0]].surf,
                        planes[tri[1]].surf,
                        planes[tri[2]].surf,
                    ]);
                    // ★★ **The tolerance measures the planes the vertex is *defined* by.**
                    //
                    // The arrangement made the coordinate and `sv.tol` as a pair — `three_planes` on
                    // one triple, `vertex_tol` on the same one. The lines above then **re-name** the
                    // vertex in the result's own surfaces, which is a *different* triple wherever four
                    // planes concur (that is why the derivation exists). The old code carried the
                    // tolerance through that swap on the grounds that "every plane through the point
                    // contains it exactly" — **true only while nothing is rotated**. A rotated plane
                    // passes a realized point within an ulp or two, not through it, so the carried
                    // tolerance bounded a distance nobody was going to measure while `validate`
                    // measured a different one (found by `replay`'s proptest; 454 of 91,394 result
                    // vertices in the corpus were short, by at most 2.1e-14 — rounding, as it should
                    // be, but rounding the record has to admit to).
                    //
                    // ★ Measured on `model.surface(..)`, the very object `validate` reads — not on the
                    // class's own `plane` copy, or this would compare two descriptions of one plane
                    // again. And `max`ed with `sv.tol` rather than replacing it: the arrangement's
                    // figure also covers the pairwise meet lines, which is a real part of what this
                    // number means.
                    let tol = [tri[0], tri[1], tri[2]]
                        .iter()
                        .map(|&i| model.surface(planes[i].surf).distance(sv.point))
                        .fold(sv.tol, f64::max);
                    model.push_vertex(def, sv.point, Some(tol))
                }
            };
            vh.insert((g, node), handle);
            Ok(handle)
        };
    // Materialize all vertex handles first. Outer then inner, rings in order: `vh`'s
    // first-appearance order fixes the vertex handles, and replay depends on it.
    for (fi, lf) in faces.iter().enumerate() {
        for &node in lf.poly_rings().flat_map(|r| r.iter()) {
            node_handle(model, group_of[fi], node)?;
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
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let fsurf = planes[lf.plane_idx].surf;
        for r in lf.poly_rings() {
            let k = r.nodes.len();
            for t in 0..k {
                let (va, vb) = (vh[&(g, r.nodes[t])], vh[&(g, r.nodes[(t + 1) % k])]);
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
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let face_surf = planes[lf.plane_idx].surf;
        let mut ring = |model: &mut Model, r: &Ring| -> Result<Loop, BoolError> {
            let handles: Vec<Handle<Vertex>> = r.nodes.iter().map(|nd| vh[&(g, *nd)]).collect();
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
        // ★ Circle bounds do not assemble through the seam table — their loops are the rim
        // machinery's (`push_edge([lateral, plane], [v, v])` + `OnSeam`), which lands with the
        // cylinder arrangement (C4b). Until that commit no production path emits one (the C2
        // population stopper holds), so reaching here with a circle is a wiring bug, not an
        // input.
        let mut ring_of = |b: &Bound| -> Result<Loop, BoolError> {
            match b {
                Bound::Ring(r) => ring(model, r),
                Bound::Circle { cyl } => {
                    unreachable!("circle bound (cyl {cyl}) reached assembly before C4b")
                }
            }
        };
        let outer = ring_of(&lf.outer)?;
        let inner: Vec<Loop> = lf
            .inner
            .iter()
            .map(ring_of)
            .collect::<Result<Vec<_>, BoolError>>()?;
        // The plane's frame *is* the root face's orientation: `frame_sign` carries that face's
        // `Forward`/`Reversed` as a sign (read off the stored flag since the cutover). Reading it
        // here is what used to be `planes[plane_idx].orient` — a face field indexed by a plane,
        // the shape of every bug this split exists to prevent.
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
        // The witness is the minimum offending edge handle, not the map's first hit: `HashMap`
        // iteration order would make the reported location differ between runs.
        if let Some(eh) = uses
            .iter()
            .filter(|&(_, &n)| n > 2)
            .map(|(&e, _)| e)
            .min_by_key(|e| e.index())
        {
            let [va, vb] = model.edges.get(eh).vertices;
            return Err(crate::reject_at(
                RejectReason::NonManifoldResultEdge,
                crate::RejectWhere::Segment([model.vertex_point(va), model.vertex_point(vb)]),
            ));
        }
        if uses.values().any(|&n| n != 2) {
            return Err(reject(RejectReason::OpenResultShell));
        }
    }
    // ★ **The grouping decided at the top of this function, raised here** — where the old code
    // decided it, so a boolean that declines leaves the arena cells it always left. What it says:
    // one component is the whole result; several mean either an enclosed void (a cavity, an
    // inward-oriented shell) or a severed operand (two or more material-enclosing shells), and a
    // surviving cavity belongs to the piece whose outer shell nests it. See [`group_faces`] — and
    // for why the *group*, not the component, is the unit the handles above were minted per.
    let Grouping {
        labels,
        n,
        positives,
        mut comps_of,
        ..
    } = grouping?;
    let mut by_comp: Vec<Vec<Handle<Face>>> = vec![Vec::new(); n];
    for (i, &fh) in face_handles.iter().enumerate() {
        by_comp[labels[i]].push(fh);
    }
    let shells: Vec<Handle<Shell>> = by_comp
        .iter()
        .map(|faces| {
            model.shells.push(Shell {
                faces: faces.clone(),
            })
        })
        .collect();
    // Emit each material solid (with its cavities) in a canonical, replay-stable order keyed on
    // geometry, so a downstream op can index the returned Vec deterministically. ★ One material
    // needs no key and must not pay for one — `comp_key` walks every face of the piece.
    let order: Vec<usize> = if positives.len() == 1 {
        vec![0]
    } else {
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
        order
    };
    let out: Result<Vec<Handle<Solid>>, BoolError> = Ok(order
        .into_iter()
        .map(|oi| {
            let c = positives[oi];
            let comps = comps_of
                .remove(&c)
                .expect("every material component owns a component list");
            // A cavity shell's faces already point into the void (the material is outside it, so
            // the material-on-correct-side reconstruction winds them inward) — measured, no flip.
            let cavities = comps
                .iter()
                .copied()
                .filter(|&x| x != c)
                .map(|x| shells[x])
                .collect();
            model.push_solid(Solid {
                outer: shells[c],
                cavities,
            })
        })
        .collect());
    debug_assert!(
        out.as_ref().is_ok_and(|solids| solids
            .iter()
            .all(|&s| crate::transform::defs_are_remappable(model, s)))
            || out.is_err(),
        "a result vertex names a surface this solid has no face on — see the def derivation above"
    );
    out
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
        for ring in lf.poly_rings() {
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
        // ★ A circle-bounded face never merges (M6-2a): `merge_component` rebuilds a
        // component's boundary from node rings alone, so a member's circle hole would
        // silently vanish from the rebuilt face — the exact "quiet pass-through" the
        // integration survey flagged. Skipping keeps a drilled cap's coplanar neighbours as
        // separate faces: a tidiness loss, never a correctness one (this pass's own charter).
        if mem
            .iter()
            .any(|&fi| kept[fi].as_ref().is_some_and(LocalFace::has_circle))
        {
            continue;
        }
        let group: Vec<&LocalFace> = mem
            .iter()
            .map(|&fi| kept[fi].as_ref().expect("member present"))
            .collect();
        // Abstention: the group pinches at a point, so its faces are emitted as-is — the
        // whole-result judgement (`self_touch_reject`, run before any cell is minted) owns
        // what this shape is about to be named for.
        let Some(rings) = merge_component(&group, jd)? else {
            continue;
        };
        let (plane_idx, flip) = (group[0].plane_idx, group[0].flip);
        merged.extend(rings.into_iter().map(|(outer, inner)| LocalFace {
            plane_idx,
            outer: Bound::Ring(outer),
            inner: inner.into_iter().map(Bound::Ring).collect(),
            flip,
        }));
        for &fi in mem {
            kept[fi] = None;
        }
    }
    let mut out: Vec<LocalFace> = kept.into_iter().flatten().collect();
    out.extend(merged);
    let all: Vec<usize> = (0..out.len()).collect();
    dissolve_straight_angles(&mut out, &all);
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
/// the holes that belong to it. `Ok(None)` is **abstention**: the group's pieces meet at a point,
/// the merged contour would be a figure-8, and there is no cycle set to re-thread — the caller
/// emits the group as-is. Merging is a tidiness pass, not a correctness one, so a merge that
/// cannot merge passes through; the whole-result judgement such a shape is headed for
/// ([`self_touch_reject`], which runs before any cell is minted) is not this function's to
/// pre-empt. (It used to reject `CoplanarPinch` here — a capability-limit name for what is, on
/// every input that reaches it, a self-touching result; the detour of naming it from here was
/// measured and declined twice before the abstention landed, dev-log 2026-08-16/17.)
fn merge_component(
    group: &[&LocalFace],
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Option<Vec<RegionRings>>, BoolError> {
    // 1. Collect directed edges **with their walls**. A repeat in the same direction means two
    //    faces claim the same side.
    // One map, `(count, wall)` — the wall rides along rather than in a second table.
    let mut dirs: HashMap<(Node, Node), (usize, usize)> = HashMap::new();
    for lf in group {
        for ring in lf.poly_rings() {
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
    // 3. Re-thread what survives. Two outgoing edges at one node means the pieces meet at a
    //    point and the cycles are not determined — the merged contour would be a figure-8:
    //    nothing to re-thread, so the merge abstains (see the function doc). Which node
    //    collided is irrelevant to the outcome, so the first collision answers.
    let mut next: HashMap<Node, Node> = HashMap::new();
    for &(a, b) in dirs.keys() {
        if dirs.contains_key(&(b, a)) {
            continue; // interior
        }
        if next.insert(a, b).is_some() {
            return Ok(None);
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
        // Every node is a probe, not just the first: which vertex can cast a clear ray is a fact
        // about that vertex, and settling for `nodes[0]` is what lost whole bands of rotation
        // angles here. `ring_in_ring` holds that retry now, for this caller and the two in the
        // arrangement alike.
        let probes: Vec<[usize; 3]> = hole.nodes.iter().map(|Node::Seam(t)| *t).collect();
        let mut owner = None;
        for (i, (outer, _)) in faces.iter().enumerate() {
            let ring = outer.edges(jd, wc)?;
            if combinatorics::ring_in_ring(jd, wc, &probes, &ring)? {
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
    Ok(Some(faces))
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
fn dissolve_straight_angles(out: &mut [LocalFace], which: &[usize]) {
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
    for &fi in which {
        let lf = &out[fi];
        for ring in lf.poly_rings() {
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
    for &fi in which {
        let lf = &mut out[fi];
        for ring in lf.poly_rings_mut() {
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
