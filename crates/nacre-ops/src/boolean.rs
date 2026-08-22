//! Boolean assembly: turn the arrangement engine's per-face output (`LocalFace`, named by
//! `NodeId` seam triples) into result solids, and the mandatory coplanar-face cleaning pass.
//!
//! The public [`boolean`] entry lives here and delegates the cell-complex work to
//! [`crate::arrangement`]; the engine calls back into [`assemble_fuse_cut`] and
//! [`unify_coplanar_faces`] to build and clean the shells (a legal module cycle).

use crate::combinatorics::{self, NodeId};
use crate::planes::{ClassIx, WorkingPlane, uf_find};
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
fn face_components(
    faces: &[LocalFace],
    cut_rims: &crate::arrangement::CutRims,
) -> (Vec<usize>, usize) {
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
                let key = match ring.walls[t] {
                    Wall::Plane(_) => JoinKey::Line(norm_edge(a, b)),
                    Wall::Arc { cyl, ccw } => {
                        let (from, to) = if ccw { (a, b) } else { (b, a) };
                        JoinKey::Arc { cyl, from, to }
                    }
                };
                users.entry(key).or_default().push(i);
            }
        }
    }
    // ★★ **A band whose rim is cut joins by its arcs, not by a rim key.** The cut circle bounds
    // no whole disk, so the second rule below cannot see it — and the cap side of each arc is a
    // *ring step*, so the join lands in the node rule instead: the band registers the same CCW
    // pairs its chain is assembled from (`CutRim.nodes`, cyclic — carried from the split, the one
    // source), and every arc meets exactly its cap face there. The wrap arc is one node pair
    // here even where the seam vertex splits it into two edges — S is a handle, not a node.
    for (i, lf) in faces.iter().enumerate() {
        let ClassIx::Cyl(k) = lf.surf else { continue };
        for b in std::iter::once(&lf.outer).chain(lf.inner.iter()) {
            let Bound::Band { lo, hi } = b else { continue };
            for c in [*lo, *hi] {
                let Some(cr) = cut_rims.get(&(k, c)) else {
                    continue;
                };
                let m = cr.nodes.len();
                for j in 0..m {
                    users
                        .entry(JoinKey::Arc {
                            cyl: k,
                            from: cr.nodes[j],
                            to: cr.nodes[(j + 1) % m],
                        })
                        .or_default()
                        .push(i);
                }
            }
        }
    }
    // ★ **The second joining rule** (M6-2a C4b): a cap face and the band that meets it on a rim
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
                (Bound::Band { lo, hi }, ClassIx::Cyl(k)) => {
                    rim_users.entry((k, *lo)).or_default().push(i);
                    rim_users.entry((k, *hi)).or_default().push(i);
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

/// **Which faces make one output solid — decided before any topology exists.**
///
/// A result solid is a *material* component plus the cavity components nested in it, and that
/// grouping is the unit [`reconstruct`] mints vertex and edge handles per: two faces in one group
/// may share a handle, two faces in different groups never do. Every question asked here is asked
/// of `LocalFace`/`NodeId`/`jd` alone — none of it reads the arena — which is what lets it run first.
///
/// ★★ **The group, not the component, is the right unit.** A cavity that touches its host's outer
/// shell is one solid with a pinch, and a pinch is counted on *handles* ([`check_result_topology`],
/// `validate`): split the handles there and the defect stops being visible while the model keeps
/// its zero-thickness material. Grouping keeps a cavity with its host, so that case still reaches
/// the reject it deserves.
pub(crate) struct Grouping {
    /// Connected-component label per face, dense `0..n`.
    labels: Vec<usize>,
    /// Component count. (`pub(crate)`: the grouping fence in `bands` asserts a cut-rim result
    /// joins into one component through the door production uses.)
    pub(crate) n: usize,
    /// The material components (even nesting depth), ascending.
    pub(crate) positives: Vec<usize>,
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
fn group_faces(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::arrangement::CutRims,
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
                .map(|lf| {
                    // ★ **Every boundary travels, in the engine's own vocabulary.** This used to
                    // keep `poly_rings()` only, which silently dropped a face's circular and
                    // banded bounds — and a ray that counts an incomplete component answers
                    // confidently and wrongly. That drop is what the old `curved_component_depth`
                    // refusal stood in for; carrying the bounds is the first half of what retired
                    // it, name and all.
                    let bound = |b: &Bound| -> Result<combinatorics::BoundEdges, BoolError> {
                        Ok(match b {
                            // `plane()` is the plane-only projection whose panic is the
                            // upstream-filter-bug detector: a polygon bound on a cylinder face
                            // would be a producer error, not an input.
                            Bound::Ring(r) => {
                                combinatorics::BoundEdges::Ring(r.edges(jd, lf.surf.plane())?)
                            }
                            Bound::Circle { cyl } => {
                                combinatorics::BoundEdges::Circle(Box::new(cyls[*cyl].def.clone()))
                            }
                            Bound::Band { lo, hi } => {
                                combinatorics::BoundEdges::Band { lo: *lo, hi: *hi }
                            }
                        })
                    };
                    let surf = match lf.surf {
                        ClassIx::Plane(c) => combinatorics::CompSurf::Plane(c),
                        // ★ The truth travels, not the index — the probe should not have to hold
                        // the class table to solve a crossing.
                        ClassIx::Cyl(k) => {
                            combinatorics::CompSurf::Cylinder(Box::new(cyls[k].def.clone()))
                        }
                    };
                    Ok(combinatorics::CompFace {
                        surf,
                        outer: bound(&lf.outer)?,
                        inner: lf
                            .inner
                            .iter()
                            .map(bound)
                            .collect::<Result<Vec<_>, BoolError>>()?,
                    })
                })
                .collect::<Result<_, BoolError>>()
        })
        .collect::<Result<_, BoolError>>()?;
    // A **probe list**, so a branch node may be dropped: `first_deciding` below tries each in
    // turn, and an exhausted list is already `Ok(None)` — the caller's own rejection.
    //
    // ★★★ **Coordinates only when there is no vertex to name.** A named probe is exact without
    // coordinates at all, so it keeps working where no rational one exists (a rotated class), and
    // it is the answer for every component that has a polygon anywhere on it. What has none is the
    // shape the retired `curved_component_depth` refused outright: a lone cylinder, whose boundary
    // is two disks and a band and whose vertex count is therefore **zero**. `coord_probes` names a
    // point on that boundary instead of a vertex of it.
    let probes_of = |c: usize| -> Vec<combinatorics::Probe> {
        let named: Vec<combinatorics::Probe> = combinatorics::three_plane_probes(
            by_comp_lf[c]
                .iter()
                .flat_map(|lf| lf.poly_rings().flat_map(|r| r.iter().copied())),
        )
        .into_iter()
        .map(combinatorics::Probe::Named)
        .collect();
        if named.is_empty() {
            combinatorics::coord_probes(jd, &comp_faces[c])
        } else {
            named
        }
    };
    // Try `f` at each node in turn: the first node that **decides** wins, a node that abstains
    // (`Ok(None)` — it grazed a boundary) is passed over for the next, and a failed judgement
    // (`Err`) propagates immediately. The last part is the point of the shape: an abstention has
    // other nodes as its remedy, a failed judgement does not — the old `Err(_) => try the next
    // node` arms retried both, so a real cause could masquerade as `NoClearRay` once every node
    // hit it. `Ok(None)` here means every node abstained; that being a reject is the *caller's*
    // proposition to raise.
    fn first_deciding<T>(
        probes: &[combinatorics::Probe],
        mut f: impl FnMut(&combinatorics::Probe) -> Result<Option<T>, BoolError>,
    ) -> Result<Option<T>, BoolError> {
        for x in probes {
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
        let depth = first_deciding(&probes_of(c), |x| {
            let mut d = 0usize;
            for other in (0..n).filter(|&o| o != c) {
                match combinatorics::probe_in_component(jd, x, &comp_faces[other])? {
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
            // `comp_faces`/`probes_of` are the same pair the material/void label above uses.
            for d in (0..n).filter(|c| !positives.contains(c)) {
                // A cavity node that classifies cleanly against *every* material (one shared origin
                // keeps the nesting consistent); its `true` materials nest, so take the innermost.
                let containers = first_deciding(&probes_of(d), |x| {
                    let mut cs = Vec::new();
                    for &m in &positives {
                        match combinatorics::probe_in_component(jd, x, &comp_faces[m])? {
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
                                let v = first_deciding(&probes_of(c), |x| {
                                    combinatorics::probe_in_component(jd, x, &comp_faces[o])
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
    pub(crate) triple: NodeId,
    pub(crate) tol: f64,
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
    pub(crate) nodes: Vec<NodeId>,
    /// `walls[i]` is the carrier of the edge `nodes[i] -> nodes[i+1]`.
    pub(crate) walls: Vec<Wall>,
}

/// **A ring edge's carrier** — the type half of "a line is unordered, a circle is ordered".
///
/// A plane-carried edge rides one wall class, as `Ring.walls` always said. An arc rides a
/// cylinder, and for it the ring additionally remembers **which way around the axis this edge
/// runs**: `ccw` restates `ClassEdges::edge_at`'s own convention (*"`MergedArc::end` runs
/// counter-clockwise about the axis"* — the even half-edge travels that way, its twin the other),
/// carried rather than re-derived. That bit is what will let `edge_for` tell the two
/// complementary arcs between one pair of branch vertices apart.
///
/// ★ Replacing the `usize::MAX` sentinel with a variant also kills a recorded hazard for free:
/// `dissolve_straight_angles` folds on wall *equality*, and two arcs of different circles — or
/// of one circle in different directions — now compare unequal instead of `MAX == MAX`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Wall {
    Plane(usize),
    Arc { cyl: usize, ccw: bool },
}

/// **The edge-welding key** — the key half of "a line is unordered, a circle is ordered"
/// (`Wall` is the type half).
///
/// Keys live in *handle* space, the space `edge_of` always keyed. A line is its unordered
/// endpoint pair, as before. Between one pair of branch vertices a circle offers **two**
/// complementary pieces, so the endpoints alone cannot name an arc — the key carries them **in
/// CCW order about the axis** (`from → to`), and the two complementary arcs get the two orders.
/// The minted edge stores its vertices in that same order, which is what the `[A, B]`-CCW
/// convention means downstream (`Model::derive_edge_curve`'s circle arm).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum EdgeKey {
    Line((usize, usize)),
    Arc { cyl: usize, from: usize, to: usize },
}

impl Ring {
    /// A ring whose walls are **derived from the node names** — for hand-built fixtures, where
    /// every vertex is a clean three-plane point so "the class the two names share besides `p`" is
    /// well defined. Production never derives; see the type's note.
    #[cfg(test)]
    pub(crate) fn from_clean_names(p: usize, nodes: Vec<NodeId>) -> Ring {
        let k = nodes.len();
        let walls = (0..k)
            .map(|i| {
                // Its own contract is "clean fixture rings only", so a branch node here is a
                // fixture bug, not an input the kernel must survive.
                let name = |n| {
                    combinatorics::three_plane_name(n).expect("a clean fixture ring names triples")
                };
                let (a, b) = (name(nodes[i]), name(nodes[(i + 1) % k]));
                Wall::Plane(
                    a.iter()
                        .copied()
                        .find(|&c| c != p && b.contains(&c))
                        .expect("a clean fixture ring edge rides one wall"),
                )
            })
            .collect();
        Ring { nodes, walls }
    }
}

impl std::ops::Deref for Ring {
    type Target = [NodeId];
    fn deref(&self) -> &[NodeId] {
        &self.nodes
    }
}

impl Ring {
    pub(crate) fn new(nodes: Vec<NodeId>, walls: Vec<Wall>) -> Ring {
        debug_assert_eq!(nodes.len(), walls.len(), "one wall per edge");
        Ring { nodes, walls }
    }

    /// This ring's edges, ready for the exact predicates.
    pub(crate) fn edges(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        p: usize,
    ) -> Result<Vec<combinatorics::RingEdge>, BoolError> {
        // ★ **A legacy shim, on purpose.** This is the names-road (the grouping's component
        // machinery and `self_touch` build their edges here), and its carrier-ization is the
        // grouping-arm cell's own item — today it must behave exactly as it always has, so an arc
        // carrier maps back to the sentinel it used to be and the road's own branch-node guard
        // stays the thing that answers.
        let legacy: Vec<usize> = self
            .walls
            .iter()
            .map(|w| match w {
                Wall::Plane(c) => *c,
                Wall::Arc { .. } => usize::MAX,
            })
            .collect();
        combinatorics::ring_edges_with_walls(jd, p, &self.nodes, &legacy)
    }
}

/// One boundary of a result face (M6-2a): a polygon of seam nodes, a **full circle** of a
/// cylinder class, or a lateral **band** between two of them. Neither curved bound has nodes or
/// walls; both are assembled through the rim machinery (`push_edge([lateral, plane], [v, v])` +
/// an `OnSeam` vertex), never through the seam-vertex table.
#[derive(Clone, Debug)]
pub(crate) enum Bound {
    Ring(Ring),
    Circle {
        cyl: usize,
    },
    /// A lateral band's whole boundary: the two plane classes its rims sit on, `lo` the one with
    /// the smaller axis parameter. The face it bounds is the cylinder itself, so the class is on
    /// [`LocalFace::surf`] rather than repeated here.
    #[allow(dead_code)] // the producer is C4b-2's band pass
    Band {
        lo: usize,
        hi: usize,
    },
}

impl Bound {
    /// The polygon ring, `None` for a curved bound — the node-walking consumers' filter.
    pub(crate) fn ring(&self) -> Option<&Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } | Bound::Band { .. } => None,
        }
    }

    pub(crate) fn ring_mut(&mut self) -> Option<&mut Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } | Bound::Band { .. } => None,
        }
    }

    /// The polygon ring, asserted — for consumers whose population cannot carry curved bounds
    /// (and tests). Panics with the caller's location.
    #[cfg_attr(not(test), allow(dead_code))]
    #[track_caller]
    pub(crate) fn expect_ring(&self) -> &Ring {
        match self {
            Bound::Ring(r) => r,
            Bound::Circle { cyl } => panic!("a polygon-only path got a circle bound (cyl {cyl})"),
            Bound::Band { lo, hi } => panic!("a polygon-only path got a band bound ({lo}..{hi})"),
        }
    }
}

/// The rim a curved bound is assembled on: `(group, cylinder class, plane class) → (seam vertex,
/// rim edge)`. The group is in the key for the reason it is in the seam vertices' — two result
/// solids never share a handle — and the rest is what makes a cap face and the band that meets it
/// pick up the *same* edge. A **cut** circle's entry carries its seam vertex and `None`: the
/// closed `[v, v]` edge spelling is false for it, and its boundary is assembled from arc pieces
/// instead.
type RimTable = HashMap<(usize, usize, usize), (Handle<Vertex>, Option<Handle<Edge>>)>;

/// A reconstructed result face: which combined plane it is on, its boundaries, and whether to
/// flip it (cut's inside-A B-pieces).
#[derive(Clone, Debug)]
pub(crate) struct LocalFace {
    /// Which class table this face's surface lives in — a plane class, or (M6-2a C4b) a cylinder
    /// class whose band this face is. Plane-only consumers project with [`ClassIx::plane`], whose
    /// panic is the upstream-filter-bug detector the type was introduced with.
    pub(crate) surf: crate::planes::ClassIx,
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

    /// Whether any boundary is **curved** — the unify partition reads this (a face with a
    /// curved bound is never a merge candidate: `merge_component` rebuilds a boundary from node
    /// rings, and a bound with no nodes would vanish from the rebuilt face).
    ///
    /// ★ Bands are covered as well as circles, even though a band reaches the merge path only if
    /// something else changes: it has no node ring, so it joins no edge-connected component and
    /// is emitted as a singleton. Stating the property is what keeps that an *invariant* rather
    /// than an accident of today's component rule.
    pub(crate) fn has_curved_bound(&self) -> bool {
        std::iter::once(&self.outer)
            .chain(self.inner.iter())
            .any(|b| matches!(b, Bound::Circle { .. } | Bound::Band { .. }))
    }
}

/// Push the reconstructed result and supersede the inputs.
///
/// ★ **The operands retire in exactly one place, and it is after the result is accepted.** The
/// retire used to be copied into each arm that produced solids, which meant every future reject had
/// to know whether it stood before or after its arm's copy. Hoisting it here makes "a reject leaves
/// the model alone" true of the *shape* of this function rather than of a fact about where the
/// rejects happen to sit today. [`reconstruct`] does the work and never touches the live set.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::arrangement::CutRims,
    deferred: Option<BoolError>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let out = reconstruct(model, jd, seam, faces, cyls, cut_rims, deferred)?;
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
    let pt: HashMap<NodeId, (Point3, f64)> = seam
        .iter()
        .map(|sv| (sv.triple, (sv.point, sv.tol)))
        .collect();
    for g in groups {
        // ★ **Planar faces only** (M6-2a C4b). The proposition this sieve tests — "an edge of
        // this solid lies in the interior of one of its own faces" — is asked through plane
        // membership and in-plane boxes, and a lateral band has neither a plane nor a node to
        // offer. It abstains rather than pretending: the edge-use guard, the non-manifold vertex
        // check and `validate` all still stand behind it, so a curved self-touch surfaces there
        // instead of being missed silently. A band that pinches against a *plane* face is
        // therefore the case this does not see yet, and it is written down rather than assumed
        // away. ★ **An edge with a branch endpoint joins that population** (M6-2b): it has no
        // three-plane name for the plane-membership question below, so it is skipped by the same
        // rule and for the same reason.
        let faces: Vec<&LocalFace> = g
            .iter()
            .flat_map(|&c| by_comp_lf[c].iter().copied())
            .filter(|lf| matches!(lf.surf, ClassIx::Plane(_)))
            .collect();
        // One pass: which faces sit on each plane, which two faces own each edge, and a box per
        // face already widened by its own vertices' tolerances.
        let mut by_plane: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut owners: HashMap<[NodeId; 2], Vec<usize>> = HashMap::new();
        let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(faces.len());
        for (j, lf) in faces.iter().enumerate() {
            by_plane.entry(lf.surf.plane()).or_default().push(j);
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for r in rings_of(lf) {
                let k = r.nodes.len();
                for i in 0..k {
                    let (u, v) = (r.nodes[i], r.nodes[(i + 1) % k]);
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
            // ★ **Abstain, exactly as the paragraph above already does for a lateral band.** An
            // edge with a branch endpoint has no three-plane name to intersect, and this is a
            // *rejection guard*: declining here would turn "I cannot check this edge" into "this
            // model is invalid" — a reject on a possibly-sound solid. The nets named above (the
            // edge-use guard, the non-manifold vertex check, `validate`) still stand behind it,
            // and this population joins the one already written down there.
            let (Some(up), Some(vp)) = (
                combinatorics::three_plane_name(u),
                combinatorics::three_plane_name(v),
            ) else {
                continue;
            };
            for q in up.iter().copied().filter(|x| vp.contains(x)) {
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
                    let Some(w2) = own.iter().map(|&o| faces[o].surf.plane()).find(|&x| x != q)
                    else {
                        continue; // an edge whose faces are both on `q` is not an edge
                    };
                    if combinatorics::segment_meets_face(jd, q, w2, up, vp, &rings_of_face[&j])? {
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

/// **Everything `reconstruct` decides before a single handle is minted** — the grouping (held,
/// not raised), the per-solid straight-angle dissolve, the whole-result self-touch judgement, and
/// every ring node's defining triple.
///
/// ★ A named function rather than the top of `reconstruct`, for the same reason `seam_table` is
/// one: the deferred arc stopper stands behind it (at the assembly's very end now) and
/// intercepts everything, so no
/// reject name can testify that the naming completed — only a fence that calls it directly on the
/// faces production feeds it can. Model-immutable by signature: nothing here takes `&mut Model`.
/// **A result vertex's definition, in class space** — what the minting turns into a `VertexDef`.
///
/// ★ `Three` is the derived triple the pre-pass has always built. `Branch` is a **declaration,
/// not a derivation**: `NodeId::Branch` already names two result plane classes and the cylinder,
/// so its def is the name's own payload (the class→handle mapping and `QuadRoot::canonical`'s
/// second answer belong to the minting, which the deferred stopper still stands in front of).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Def {
    Three([usize; 3]),
    Branch {
        planes: [usize; 2],
        cyl: usize,
        root: nacre_topo::QuadRoot,
    },
}

pub(crate) struct Named {
    /// The grouping, **held** — raised where the old code raised it, deep in the minting.
    /// (`pub(crate)`: the grouping fence reads the held result directly.)
    pub(crate) grouping: Result<Grouping, BoolError>,
    pub(crate) group_of: Vec<usize>,
    /// The working face list where it differs from the caller's — the split-twin subdivision
    /// and/or the per-solid dissolve rewrote rings. `None` = the caller's list stands.
    pub(crate) per_solid: Option<Vec<LocalFace>>,
    pub(crate) defs: HashMap<(usize, NodeId), Def>,
}

pub(crate) fn name_result_vertices(
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::arrangement::CutRims,
) -> Result<Named, BoolError> {
    // ★★★ **The split-twin subdivision — every branch node is cut into every edge it lies on.**
    // The arrangement cannot do this: a branch point needs the cylinder, and the neighbouring
    // plane class has no circle (the cylinder is parallel to it), so no vocabulary for the point
    // — measured, the T-junction finding. Only here, where every class's rings are in one hand,
    // can the plate-top's whole edge learn that the arc class subdivided its twin. Without this,
    // `norm_edge` keys never match across such a pair, `far_plane` starves, and the future edge
    // welding has no twins to weld.
    //
    // A branch node lies on an edge's carrier line exactly when its plane pair *is* the edge's
    // `{own, wall}` — a name fact, no geometry — and betweenness is `branch_between`'s exact
    // half. An edge whose order cannot be formed is left unsplit, which is precisely today's
    // behaviour (starved `far_plane`, the walls-fallback net); the conservative arm degrades to
    // the state this pass improves, never to something new. Identity for every non-arc input:
    // no branch nodes, no pairs, no rewrite.
    let mut by_pair: HashMap<[usize; 2], Vec<NodeId>> = HashMap::new();
    for lf in faces {
        for ring in lf.poly_rings() {
            for &n in ring.iter() {
                if let Some((pair, _, _)) = combinatorics::branch_name(n) {
                    let v = by_pair.entry(pair).or_default();
                    if !v.contains(&n) {
                        v.push(n);
                    }
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
                        continue; // an arc edge: the arrangement's own split made it
                    };
                    let mut pair = [own, w];
                    pair.sort_unstable();
                    let Some(cands) = by_pair.get(&pair) else {
                        continue;
                    };
                    if let Some(bet) = combinatorics::branch_between(jd, cyls, a, b, cands) {
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
    // ★★ **Held, not raised.** A failed grouping is reported further down, where the old code
    // reported it, and until then every face is one group — which is exactly the keying this
    // function used before groups existed. So a boolean that declines pushes the arena cells it
    // always did (`replay::a_late_reject_is_not_index_neutral` measures that).
    let grouping = group_faces(jd, faces, cyls, cut_rims);
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
    let mut edge_faces: HashMap<(usize, (NodeId, NodeId)), Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        for ring in lf.poly_rings() {
            let k = ring.len();
            for t in 0..k {
                edge_faces
                    .entry((group_of[fi], norm_edge(ring[t], ring[(t + 1) % k])))
                    .or_default()
                    .push(lf.surf.plane());
            }
        }
    }
    // The plane on the other side of an edge, **within this solid**. `None` rather than an error:
    // a corner this cannot
    // resolve is one to skip, and the edge-use guard further down is what judges the face set.
    let far_plane = |g: usize, a: NodeId, b: NodeId, own: usize| -> Option<usize> {
        let mut others = edge_faces
            .get(&(g, norm_edge(a, b)))?
            .iter()
            .copied()
            .filter(|&x| x != own);
        let o = others.next()?;
        others.all(|x| x == o).then_some(o)
    };
    let mut def_triple: HashMap<(usize, NodeId), Def> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        for ring in lf.poly_rings() {
            let k = ring.len();
            for t in 0..k {
                let node = ring[t];
                if def_triple.contains_key(&(g, node)) {
                    continue;
                }
                // ★ A branch vertex's def is a **declaration, not a derivation** — the name
                // already carries its two result plane classes and its cylinder, so the far-plane
                // road (which can only ever answer in planes) is not asked.
                if let Some((planes2, cyl, root)) = combinatorics::branch_name(node) {
                    def_triple.insert(
                        (g, node),
                        Def::Branch {
                            planes: planes2,
                            cyl,
                            root,
                        },
                    );
                    continue;
                }
                let (Some(prev), Some(next)) = (
                    far_plane(g, ring[(t + k - 1) % k], node, lf.surf.plane()),
                    far_plane(g, node, ring[(t + 1) % k], lf.surf.plane()),
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
                if jd.plane_pair_dir_sign(lf.surf.plane(), prev, next) == 0 {
                    continue;
                }
                let mut tri = [lf.surf.plane(), prev, next];
                tri.sort_unstable();
                def_triple.insert((g, node), Def::Three(tri));
            }
        }
    }
    // ★★ **The walls-fallback — a zero-population net since the split-twin subdivision.** It was
    // built for the arc population's bitten corner (a disk eating a plate corner starved every
    // incident face's far-plane road — measured, exactly one such node); the subdivision above
    // now matches those twins, so the corner derives on the main road and nothing reaches here
    // today. Kept as the net for any future input whose far-plane road starves in a way the
    // subdivision does not repair (an unsplittable edge takes exactly that path), under the same
    // two guards as the derivation above plus one more: **every deriving face must agree** — the
    // four-plane hazard (a wall that carries the edge's line without bounding the solid) is
    // exactly where faces could disagree, so a split vote stays def-less rather than picking a
    // winner. A genuinely straight corner is still refused by the guards and keeps its
    // `StraightAngle`.
    let mut fallback: HashMap<(usize, NodeId), Option<[usize; 3]>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        for ring in lf.poly_rings() {
            let k = ring.len();
            for t in 0..k {
                let node = ring[t];
                if def_triple.contains_key(&(g, node)) {
                    continue;
                }
                // An arc edge carries a cylinder, not a plane — no wall to name a triple with.
                let (Wall::Plane(prev), Wall::Plane(next)) =
                    (ring.walls[(t + k - 1) % k], ring.walls[t])
                else {
                    continue;
                };
                if prev == next {
                    continue;
                }
                if jd.plane_pair_dir_sign(lf.surf.plane(), prev, next) == 0 {
                    continue;
                }
                let mut tri = [lf.surf.plane(), prev, next];
                tri.sort_unstable();
                let vote = fallback.entry((g, node)).or_insert(Some(tri));
                if matches!(vote, Some(seen) if *seen != tri) {
                    *vote = None; // a split vote stays def-less
                }
            }
        }
    }
    for ((g, node), tri) in fallback {
        if let Some(tri) = tri {
            def_triple.insert((g, node), Def::Three(tri));
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

/// Rebuild the result solids from the arrangement's faces. Pushes into the arena; it does not take
/// the operands at all, which is the point — the live set is its caller's to move.
fn reconstruct(
    model: &mut Model,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::arrangement::CutRims,
    deferred: Option<BoolError>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes = jd.planes;
    // No faces means no result — `Common` of two solids that miss each other, `Cut` of a box that
    // is wholly inside what cuts it. That is an answer, not a failure: a solid is bounded by faces,
    // so a non-empty result cannot have none. The inputs are still consumed, exactly as they are on
    // any other successful boolean — this returns `Ok`, so the caller's retire runs.
    if faces.is_empty() {
        return Ok(Vec::new());
    }
    let named = name_result_vertices(jd, seam, faces, cyls, cut_rims);
    // The naming's failure yields to the deferred stopper like every stage before it; the raise
    // itself now stands at the assembly's very end below.
    let Named {
        grouping,
        group_of,
        per_solid,
        defs: def_triple,
    } = match named {
        Ok(n) => n,
        Err(e) => return Err(deferred.unwrap_or(e)),
    };
    let faces: &[LocalFace] = per_solid.as_deref().unwrap_or(faces);

    // Vertices (deterministic: first appearance across faces in order).
    let mut vh: HashMap<(usize, NodeId), Handle<Vertex>> = HashMap::new();
    let mut node_handle =
        |model: &mut Model, g: usize, node: NodeId| -> Result<Handle<Vertex>, BoolError> {
            if let Some(&h) = vh.get(&(g, node)) {
                return Ok(h);
            }
            // A face references a seam node that was not welded into `seam` — a reconstruction
            // dropped a crossing. Reject (never panic): an unmodeled flush topology must decline
            // honestly, not abort the kernel (DNA).
            //
            // ★ This used to `match` the node to compare its triple against `SeamVertex.triple`.
            // Both are identities now, so the comparison is the identity's own `==` and the
            // destructure that existed only to reach the payload is gone.
            let handle = {
                let sv = seam
                    .iter()
                    .find(|s| s.triple == node)
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                // A vertex that is a corner of no face at all has no name in the result's own
                // planes — a degeneracy, and the honest answer is the one this reason already
                // carries ("a corner with no turn").
                let tri = match def_triple.get(&(g, node)).copied() {
                    Some(Def::Three(t)) => t,
                    // ★★★ **A branch vertex is minted from its declaration** — the def already
                    // names two result plane classes and the cylinder, so this is the class→handle
                    // mapping and nothing else. That mapping is where `QuadRoot::canonical`
                    // answers a **second** time: `NodeId::Branch` is canonical in *class* order,
                    // `VertexDef::Branch` in *handle* order, and the class→handle map is not
                    // monotone in general — a re-sort must carry the root through
                    // (`transform`'s remap already locks the same rule on the way back out).
                    // ★ The flip is **unexercised in today's corpus, measured** — dropping the
                    // canonical call leaves every fence green, because both fixtures' class
                    // order happens to match their handle order. The rule itself is pinned by
                    // `QuadRoot::canonical`'s own locks in `nacre-topo` (the +X-axis fixtures
                    // that forced `Lo|Hi|Double`); the population that turns this call red is a
                    // boss whose cylinder classes arrive in the other order, and it gets its
                    // fixture when it arrives rather than a contorted one now.
                    // The tolerance is the Three arm's rule below, one surface swapped: measured
                    // against the very `model.surface` objects `validate` reads, maxed with
                    // `sv.tol` (which `branch_vertex_tol` built, meet-line term included).
                    Some(Def::Branch {
                        planes: p2,
                        cyl,
                        root,
                    }) => {
                        let (pair, root) = nacre_topo::QuadRoot::canonical(
                            [planes[p2[0]].surf, planes[p2[1]].surf],
                            root,
                        );
                        let cylinder = cyls[cyl].surf;
                        let tol = pair
                            .iter()
                            .chain(std::iter::once(&cylinder))
                            .map(|&s| model.surface(s).distance(sv.point))
                            .fold(sv.tol, f64::max);
                        let def = VertexDef::Branch {
                            planes: pair,
                            cylinder,
                            root,
                        };
                        let h = model.push_vertex(def, sv.point, Some(tol));
                        vh.insert((g, node), h);
                        return Ok(h);
                    }
                    None => return Err(reject(RejectReason::StraightAngle)),
                };
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
            };
            vh.insert((g, node), handle);
            Ok(handle)
        };
    // Materialize all vertex handles first. Outer then inner, rings in order: `vh`'s
    // first-appearance order fixes the vertex handles, and replay depends on it.
    let mut materialized = Ok(());
    'mat: for (fi, lf) in faces.iter().enumerate() {
        for &node in lf.poly_rings().flat_map(|r| r.iter()) {
            if let Err(e) = node_handle(model, group_of[fi], node) {
                materialized = Err(e);
                break 'mat;
            }
        }
    }
    // The materialization's failure yields to the deferred stopper like the naming's above; the
    // raise itself now stands at the assembly's very end below.
    if let Err(e) = materialized {
        return Err(deferred.unwrap_or(e));
    }

    // ★★ **The rim table** (M6-2a C4b): a curved bound has no node and no triple, so it is not
    // welded through the seam table at all — it is minted here, once per `(group, cylinder class,
    // plane class)`, and every face that meets that rim asks for the same pair back. That sharing
    // is what makes the rim edge's use count two (the cap face once, the band once) rather than
    // two separate edges the closed-shell guard would call dangling.
    //
    // The group is in the key for the reason it is in `vh`'s: two result solids must not share a
    // handle.
    let mut rim: RimTable = HashMap::new();
    let rims_built = (|| -> Result<(), BoolError> {
        for (fi, lf) in faces.iter().enumerate() {
            let g = group_of[fi];
            let keys: Vec<(usize, usize)> = std::iter::once(&lf.outer)
                .chain(lf.inner.iter())
                .flat_map(|b| match (b, lf.surf) {
                    (Bound::Circle { cyl }, ClassIx::Plane(c)) => vec![(*cyl, c)],
                    (Bound::Band { lo, hi }, ClassIx::Cyl(k)) => vec![(k, *lo), (k, *hi)],
                    _ => Vec::new(),
                })
                .collect();
            for (k, c) in keys {
                if rim.contains_key(&(g, k, c)) {
                    continue;
                }
                let cut = cut_rims.get(&(k, c));
                let (lat, plane) = (cyls[k].surf, planes[c].surf);
                // The seam point of this rim, spelled as `add_cylinder` spells one: the axis meets the
                // plane at the circle's centre, and `θ = 0` is the `+ref_dir` side of it.
                let cache = cyls[k].cache;
                let centre = nacre_geom::intersect::line_plane(
                    &cache.axis(),
                    match model.surface(plane) {
                        nacre_geom::Surface::Plane(p) => p,
                        _ => return Err(reject(RejectReason::ThreePlanes)),
                    },
                )
                .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                let point = centre + cache.ref_dir() * cache.radius();
                // The tolerance is measured, not assumed — the same rule the three-plane vertices
                // above follow: how far the realized point sits from each surface that defines it.
                let tol = model
                    .surface(lat)
                    .distance(point)
                    .max(model.surface(plane).distance(point));
                // ★★ **A cut circle mints no closed rim edge — but its seam vertex stands.**
                // The `[v, v]` spelling is false for a cut circle (measured: both arc fixtures
                // minted one that only the reject then discarded), while θ = 0 is still where the
                // band's joint must sit (a `[lat, lat]` seam edge derives a line from its
                // endpoints, so both must share one θ — and the seam is model geometry, fixed at
                // `+ref_dir`). Degenerate case first: when a branch vertex lies **on** the seam
                // generator — the split's own `SeamIncident` classification, carried in
                // `CutRim::seam_is_node`, never re-derived from coordinates — that vertex *is*
                // the seam point and minting another would stand a second handle on the same
                // point (a zero-length arc piece nothing downstream could see).
                let v = match cut {
                    Some(cr) if cr.seam_is_node => vh[&(g, cr.nodes[0])],
                    _ => model.push_vertex(VertexDef::OnSeam([lat, plane]), point, Some(tol)),
                };
                let e = match cut {
                    Some(_) => None,
                    None => Some(
                        model
                            .push_edge(Edge::carrier_pair(lat, plane), [v, v])
                            .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?,
                    ),
                };
                rim.insert((g, k, c), (v, e));
            }
        }
        Ok(())
    })();
    // A rim's failure yields to the deferred stopper like every stage before it.
    if let Err(e) = rims_built {
        return Err(deferred.unwrap_or(e));
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
        // A band contributes no node edge — its rims and seam are minted with their carriers
        // stated (`[lateral, plane]`, `[lateral, lateral]`), so there is nothing for this scan
        // to read off it.
        let ClassIx::Plane(fc) = lf.surf else {
            continue;
        };
        let fsurf = planes[fc].surf;
        for r in lf.poly_rings() {
            let k = r.nodes.len();
            for t in 0..k {
                // ★ Only plane-carried edges feed the scan. An arc edge shares its vertex pair
                // with the chord between the same branch vertices, and pushing this face's plane
                // here would pollute the chord's line-key entry into the fallback arm — the arc
                // states its carriers directly instead (`edge_for`'s arc arm).
                if !matches!(r.walls[t], Wall::Plane(_)) {
                    continue;
                }
                let (va, vb) = (vh[&(g, r.nodes[t])], vh[&(g, r.nodes[(t + 1) % k])]);
                pair_surfs
                    .entry(unordered(va.index() as usize, vb.index() as usize))
                    .or_default()
                    .push(fsurf);
            }
        }
    }

    // Edges keyed by `EdgeKey` (lookup only).
    let mut edge_of: HashMap<EdgeKey, Handle<Edge>> = HashMap::new();
    // A ring's two consecutive nodes are distinct arrangement vertices, so their points differ and
    // the line through them exists. Reject rather than panic if it does not: an aborting kernel is
    // below the floor (`overview.md`: out-of-coverage input declines honestly). `SEAM_ALIAS`
    // catches the known way this happens — two triples on one point — at the seam table, where the
    // names are still in hand, so this is a backstop with no firing test (cf. `NON_MANIFOLD_EDGE`).
    let mut edge_for = |model: &mut Model,
                        va: Handle<Vertex>,
                        vb: Handle<Vertex>,
                        wall: Wall,
                        face_surf: Handle<Surface>|
     -> Result<Handle<Edge>, BoolError> {
        match wall {
            Wall::Plane(w) => {
                let pair = unordered(va.index() as usize, vb.index() as usize);
                let key = EdgeKey::Line(pair);
                if let Some(&e) = edge_of.get(&key) {
                    return Ok(e);
                }
                // Manifold edges (everything a green result contains) have exactly two uses. Any
                // other count is on its way to the existing non-manifold reject — the fallback
                // (this face's wall + plane) keeps construction deterministic until that reject
                // fires, deciding nothing new.
                let surfaces = match pair_surfs[&pair][..] {
                    [a, b] => Edge::carrier_pair(a, b),
                    _ => Edge::carrier_pair(planes[w].surf, face_surf),
                };
                let e = model
                    .push_edge(surfaces, [va, vb])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
                edge_of.insert(key, e);
                Ok(e)
            }
            Wall::Arc { cyl, ccw } => {
                // ★★★ **An arc edge is minted in CCW order** — `[A, B]` is the piece from A to
                // B counter-clockwise about the axis, so the two complementary arcs between one
                // branch pair are `[A, B]` and `[B, A]`: the vertex order is the last bit the
                // endpoints alone cannot give (`EdgeKey`'s note). The carriers are stated
                // directly — the two faces using an arc edge are this cap and the cylinder's
                // side, so there is nothing for the scan to read — which also keeps the chord's
                // line key clean (`pair_surfs` skips arc edges for the same reason).
                let (from, to) = if ccw { (va, vb) } else { (vb, va) };
                let key = EdgeKey::Arc {
                    cyl,
                    from: from.index() as usize,
                    to: to.index() as usize,
                };
                if let Some(&e) = edge_of.get(&key) {
                    return Ok(e);
                }
                let e = model
                    .push_edge(Edge::carrier_pair(cyls[cyl].surf, face_surf), [from, to])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
                edge_of.insert(key, e);
                Ok(e)
            }
        }
    };

    // ★★ **A cut rim's boundary chain, minted once per `(group, cylinder, plane)` — before the
    // face loop, in face order.** The chain walks the circle CCW from the seam vertex through the
    // branch nodes and back (the wrap arc's two seam-split pieces included), so a cut rim's
    // pieces exist under their `EdgeKey`s before either consumer asks: the cap faces' ring steps
    // weld to these very handles by key, and `band_loop` reads the chain directly — which is
    // also why this is not minted inside `band_loop`: two closures cannot both own `edge_for`.
    // Iterated over the faces' curved bounds (not the map) so the minting order is
    // deterministic, the same discipline as the rim table above.
    let mut band_chains: HashMap<(usize, usize, usize), Vec<HalfEdge>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let keys: Vec<(usize, usize)> = std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .flat_map(|b| match (b, lf.surf) {
                (Bound::Circle { cyl }, ClassIx::Plane(c)) => vec![(*cyl, c)],
                (Bound::Band { lo, hi }, ClassIx::Cyl(k)) => vec![(k, *lo), (k, *hi)],
                _ => Vec::new(),
            })
            .collect();
        for (k, c) in keys {
            if band_chains.contains_key(&(g, k, c)) {
                continue;
            }
            let Some(cr) = cut_rims.get(&(k, c)) else {
                continue; // an uncut rim's boundary is its closed edge, no chain to build
            };
            // The rim table ran this very enumeration, so the entry exists whenever this loop
            // reaches the key — the arm is a backstop, not a population (cf. the ring closure's
            // `MissingSeam`, which *is* reachable: a cap ring's arc steps do not put the key in
            // either enumeration, only a curved bound does).
            let Some(&(sv, _)) = rim.get(&(g, k, c)) else {
                continue;
            };
            let mut vs: Vec<Handle<Vertex>> = Vec::with_capacity(cr.nodes.len() + 1);
            if !cr.seam_is_node {
                vs.push(sv);
            }
            vs.extend(cr.nodes.iter().map(|n| vh[&(g, *n)]));
            debug_assert_eq!(vs[0], sv, "the chain starts at the seam vertex");
            let m = vs.len();
            let mut chain = Vec::with_capacity(m);
            for i in 0..m {
                let (u, v) = (vs[i], vs[(i + 1) % m]);
                let e = edge_for(model, u, v, Wall::Arc { cyl: k, ccw: true }, planes[c].surf)?;
                let forward = model.edges.get(e).vertices[0] == u;
                chain.push(HalfEdge { edge: e, forward });
            }
            band_chains.insert((g, k, c), chain);
        }
    }

    let mut face_handles = Vec::new();
    let mut assembled: Result<(), BoolError> = Ok(());
    'faces: for (fi, lf) in faces.iter().enumerate() {
        let g = group_of[fi];
        let face_surf = match lf.surf {
            ClassIx::Plane(c) => planes[c].surf,
            ClassIx::Cyl(k) => cyls[k].surf,
        };
        let mut ring = |model: &mut Model, r: &Ring| -> Result<Loop, BoolError> {
            let handles: Vec<Handle<Vertex>> = r.nodes.iter().map(|nd| vh[&(g, *nd)]).collect();
            let k = handles.len();
            let mut half_edges: Vec<HalfEdge> = Vec::with_capacity(k);
            for t in 0..k {
                let (va, vb) = (handles[t], handles[(t + 1) % k]);
                // ★★ **The wrap arc is minted as two pieces, split at the seam vertex.** θ = 0
                // lies inside exactly one arc of a cut circle (unless a branch vertex sits on the
                // seam — `CutRim::seam_is_node`, in which case nothing splits), and the band's
                // seam edge must end there; splitting on the *cap* side is what hands the band
                // the same two edges and keeps the shell guard's use count at two.
                //
                // ★ The wrap test is **directed**: with two branch nodes the two complementary
                // arcs share one unordered endpoint pair, so the match orients the ring step by
                // the `ccw` bit the wall carries and compares against the split's own θ order
                // (`nodes.last() → nodes[0]` is the piece that wraps past θ = 0).
                let split_at = if let Wall::Arc { cyl, ccw } = r.walls[t] {
                    let c = lf.surf.plane();
                    match cut_rims.get(&(cyl, c)) {
                        Some(cr) if !cr.seam_is_node => {
                            let (a, b) = (r.nodes[t], r.nodes[(t + 1) % k]);
                            let ccw_pair = if ccw { (a, b) } else { (b, a) };
                            let last = *cr.nodes.last().expect("a cut circle has branch nodes");
                            if ccw_pair == (last, cr.nodes[0]) {
                                let &(v, _) = rim
                                    .get(&(g, cyl, c))
                                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                                Some(v)
                            } else {
                                None
                            }
                        }
                        _ => None,
                    }
                } else {
                    None
                };
                let legs = match split_at {
                    Some(s) => vec![[va, s], [s, vb]],
                    None => vec![[va, vb]],
                };
                for [u, v] in legs {
                    let e = edge_for(model, u, v, r.walls[t], face_surf)?;
                    // For an arc edge the stored order is CCW, so this reads back exactly the
                    // `ccw` bit the wall carried in.
                    let forward = model.edges.get(e).vertices[0] == u;
                    half_edges.push(HalfEdge { edge: e, forward });
                }
            }
            Ok(Loop { half_edges })
        };
        // ★★ **Which way a rim circle is walked.** `derive_edge_curve` builds it as
        // `Circle::from_center_normal(centre, axis, ref_dir, r)`, so its parameter runs **CCW
        // about the axis direction `d`**. A loop must run CCW about its own face's outward
        // normal, so a bound on a face whose normal agrees with `d` is walked forward and one
        // whose normal opposes it backward — exactly the convention `add_cylinder` writes down
        // (top cap `forward: true`, bottom cap `false`). A *hole* runs the other way again,
        // because an inner loop keeps the material on its left by winding against the outer.
        //
        // The face's outward normal is the stored surface normal when the face is `Forward` and
        // its negation when `Reversed`, which is the `flip` decision made below — so the sign is
        // read off `frame_sign` and `flip` here rather than carried in from anywhere.
        // The **unflipped** face normal against the axis: `flip` is applied once, to every kind of
        // bound, right below — reading it here too would toggle the winding twice.
        let axis_sign = |k: usize| -> f64 {
            let c = lf.surf.plane();
            let s = if planes[c].frame_sign > 0 { 1.0 } else { -1.0 };
            planes[c]
                .plane
                .normal()
                .dot(cyls[k].cache.axis().direction())
                * s
        };
        let circle_loop =
            |model: &mut Model, cyl: usize, cls: usize, hole: bool| -> Result<Loop, BoolError> {
                // ★★ **The cut check comes before the rim lookup, and answers with the
                // population's own name.** A cut circle cannot bound a whole disk — the trace
                // subdivides that disk into cells — so reaching here with one is a producer
                // inconsistency this backstop names honestly rather than as a dropped crossing
                // (`MissingSeam` would misdiagnose a `SuspectedDefect`).
                if cut_rims.contains_key(&(cyl, cls)) {
                    return Err(reject(RejectReason::ArcBoundNotYet));
                }
                let (_, e) = *rim
                    .get(&(g, cyl, cls))
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                // An uncut circle's rim always carries its closed edge; the cut arm was refused
                // one line above, so `None` here is the same dropped-crossing shape.
                let e = e.ok_or_else(|| reject(RejectReason::MissingSeam))?;
                let _ = model;
                let mut forward = axis_sign(cyl) > 0.0;
                if hole {
                    forward = !forward;
                }
                Ok(Loop {
                    half_edges: vec![HalfEdge { edge: e, forward }],
                })
            };
        // A band's boundary is the lateral face `add_cylinder` builds: two rims joined by the
        // seam, which the one face uses twice in opposite senses.
        let band_loop =
            |model: &mut Model, k: usize, lo: usize, hi: usize| -> Result<Loop, BoolError> {
                // ★ Both rims cut is unreachable today — `chamber` finds no disk label at either
                // end and refuses the band before it is emitted — so no chain code is written
                // for a population nothing can reach; the honest name stands in its place.
                if cut_rims.contains_key(&(k, lo)) && cut_rims.contains_key(&(k, hi)) {
                    return Err(reject(RejectReason::ArcBoundNotYet));
                }
                let (v_lo, e_lo) = *rim
                    .get(&(g, k, lo))
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                let (v_hi, e_hi) = *rim
                    .get(&(g, k, hi))
                    .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                let lat = cyls[k].surf;
                let seam_edge = model
                    .push_edge(Edge::carrier_pair(lat, lat), [v_lo, v_hi])
                    .ok_or_else(|| reject(RejectReason::ZeroLengthEdge))?;
                // ★★ **A rim's traversal: the closed edge, or the pre-minted chain.** Both walk
                // the circle CCW; the `lo` rim is walked forward and the `hi` rim backward
                // (reverse the pieces and flip each sense — one rule for the closed edge and the
                // chain alike), exactly the senses the four-half-edge spelling always had.
                let walk =
                    |c: usize, closed: Option<Handle<Edge>>| -> Result<Vec<HalfEdge>, BoolError> {
                        match closed {
                            Some(e) => Ok(vec![HalfEdge {
                                edge: e,
                                forward: true,
                            }]),
                            None => band_chains
                                .get(&(g, k, c))
                                .cloned()
                                .ok_or_else(|| reject(RejectReason::MissingSeam)),
                        }
                    };
                let mut half_edges = walk(lo, e_lo)?;
                half_edges.push(HalfEdge {
                    edge: seam_edge,
                    forward: true,
                });
                let mut hi_hes = walk(hi, e_hi)?;
                hi_hes.reverse();
                for he in &mut hi_hes {
                    he.forward = !he.forward;
                }
                half_edges.extend(hi_hes);
                half_edges.push(HalfEdge {
                    edge: seam_edge,
                    forward: false,
                });
                Ok(Loop { half_edges })
            };
        let mut ring_of = |b: &Bound, hole: bool| -> Result<Loop, BoolError> {
            let mut lp = match b {
                Bound::Ring(r) => ring(model, r)?,
                Bound::Circle { cyl } => circle_loop(model, *cyl, lf.surf.plane(), hole)?,
                Bound::Band { lo, hi } => {
                    let k = lf
                        .surf
                        .cyl()
                        .ok_or_else(|| reject(RejectReason::MissingSeam))?;
                    band_loop(model, k, *lo, *hi)?
                }
            };
            if lf.flip {
                // Cut's inside-A B-pieces, and a bore's wall: reverse every loop and toggle the
                // orientation below, so the outward normal points into the removed region and
                // each loop still keeps material on its left. ★ One place for every bound kind —
                // a curved loop reversed only by its own builder was how the first drilled box
                // came out with its rim used twice in the same sense (`NonOpposedEdge`).
                lp.half_edges.reverse();
                for he in &mut lp.half_edges {
                    he.forward = !he.forward;
                }
            }
            Ok(lp)
        };
        let outer = match ring_of(&lf.outer, false) {
            Ok(lp) => lp,
            Err(e) => {
                assembled = Err(e);
                break 'faces;
            }
        };
        let inner: Vec<Loop> = match lf
            .inner
            .iter()
            .map(|b| ring_of(b, true))
            .collect::<Result<Vec<_>, BoolError>>()
        {
            Ok(v) => v,
            Err(e) => {
                assembled = Err(e);
                break 'faces;
            }
        };
        // The plane's frame *is* the root face's orientation: `frame_sign` carries that face's
        // `Forward`/`Reversed` as a sign (read off the stored flag since the cutover). Reading it
        // here is what used to be `planes[plane_idx].orient` — a face field indexed by a plane,
        // the shape of every bug this split exists to prevent.
        // A cylinder's stored normal points away from its axis, and that *is* the face's outward
        // normal for a boss — so its "framed" sense is `Forward`, and `flip` (material outside
        // the wall: a hole) is what reverses it. Planes read their class's frame instead.
        let framed = match lf.surf {
            ClassIx::Cyl(_) => Orientation::Forward,
            ClassIx::Plane(c) if planes[c].frame_sign > 0 => Orientation::Forward,
            ClassIx::Plane(_) => Orientation::Reversed,
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
            surface: face_surf,
            outer,
            inner,
            orientation,
        }));
    }
    // The face loop's failure yields to the deferred stopper like every stage before it; the
    // raise itself now stands at the assembly's very end below (an incomplete face set has no
    // shell to count, so an assembly failure still stops here).
    if let Err(e) = assembled {
        return Err(deferred.unwrap_or(e));
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
            return Err(deferred.unwrap_or(crate::reject_at(
                RejectReason::NonManifoldResultEdge,
                crate::RejectWhere::Segment([model.vertex_point(va), model.vertex_point(vb)]),
            )));
        }
        if uses.values().any(|&n| n != 2) {
            return Err(deferred.unwrap_or(reject(RejectReason::OpenResultShell)));
        }
    }
    // ★ **The grouping decided at the top of this function, raised here** — where the old code
    // decided it, so a boolean that declines leaves the arena cells it always left. What it says:
    // one component is the whole result; several mean either an enclosed void (a cavity, an
    // inward-oriented shell) or a severed operand (two or more material-enclosing shells), and a
    // surviving cavity belongs to the piece whose outer shell nests it. See [`group_faces`] — and
    // for why the *group*, not the component, is the unit the handles above were minted per.
    // The grouping's failure yields to the deferred stopper like every stage before it; the
    // raise itself now stands at the very end, past the solid assembly below.
    let Grouping {
        labels,
        n,
        positives,
        mut comps_of,
        ..
    } = match grouping {
        Ok(g) => g,
        Err(e) => return Err(deferred.unwrap_or(e)),
    };
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
    // ★★★ **The deferred arc stopper's raise — the very end of the assembly.** The interception
    // is the same shape at its seventh layer (per-class → seam stretch → naming → vertex
    // materialization → face loop → shell guard → grouping and the solid assembly): everything
    // ran — the shells and solids stand in the store — and the stopper's reject wins over
    // whatever any stage said, so an arc population's name never depends on how far the
    // pipeline got.
    //
    // ★★ **An arc reject therefore leaves a complete garbage solid in the store —
    // deliberately.** Vertices, edges, faces, shells and the solid itself: cells outside the
    // live set (`assemble_fuse_cut` retires operands only on `Ok`, and nothing pushed here is
    // reachable from a live solid), the same class of residue a late reject's arena cells have
    // always been. The live-set fences stay green, and a session that keeps recording after a
    // reject rebuilds from the log (`replay`'s discipline, stated in `docs/design.md`). Holding
    // the raise any earlier would put the solid assembly behind an interception nothing can see
    // past — the unreachable-machinery trap this ladder keeps refusing.
    if let Some(d) = deferred {
        return Err(d);
    }
    out
}

/// An unordered edge key: the two nodes in a fixed order, so `{a,b}` and `{b,a}` collide.
/// (`pub(crate)` for the twin-match fence, which counts exactly these keys.)
pub(crate) fn norm_edge(a: NodeId, b: NodeId) -> (NodeId, NodeId) {
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
    let group_key = |lf: &LocalFace| -> (usize, bool) { (lf.surf.plane(), lf.flip) };

    // Edge-connected components within a group.
    let mut comp: Vec<usize> = (0..n).collect();
    let mut carriers: HashMap<(NodeId, NodeId), Vec<usize>> = HashMap::new();
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
            .any(|&fi| kept[fi].as_ref().is_some_and(LocalFace::has_curved_bound))
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
        let (plane_idx, flip) = (group[0].surf.plane(), group[0].flip);
        merged.extend(rings.into_iter().map(|(outer, inner)| LocalFace {
            surf: crate::planes::ClassIx::Plane(plane_idx),
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
fn ring_edges(ring: &[NodeId]) -> impl Iterator<Item = (NodeId, NodeId)> + '_ {
    (0..ring.len()).map(move |i| (ring[i], ring[(i + 1) % ring.len()]))
}

/// A ring's directed edges **with the carrier each rides**.
fn ring_edges_walled(ring: &Ring) -> impl Iterator<Item = ((NodeId, NodeId), Wall)> + '_ {
    (0..ring.len()).map(move |i| {
        (
            (ring.nodes[i], ring.nodes[(i + 1) % ring.len()]),
            ring.walls[i],
        )
    })
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
    let mut dirs: HashMap<(NodeId, NodeId), (usize, Wall)> = HashMap::new();
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
    let mut undirected: HashMap<(NodeId, NodeId), usize> = HashMap::new();
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
    let mut next: HashMap<NodeId, NodeId> = HashMap::new();
    for &(a, b) in dirs.keys() {
        if dirs.contains_key(&(b, a)) {
            continue; // interior
        }
        if next.insert(a, b).is_some() {
            return Ok(None);
        }
    }
    let mut starts: Vec<NodeId> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut seen: HashSet<NodeId> = HashSet::new();
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
    let wc = group[0].surf.plane();
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
        let probes = combinatorics::three_plane_probes(hole.nodes.iter().copied());
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
///
/// ★★ **Arc-bearing rings flow through here, and the carrier keeps equality honest** (2026-08-23).
/// While the wall was the `usize::MAX` sentinel, two *consecutive arcs* compared as "same wall"
/// and the branch vertex between them would have dissolved — a recorded hazard, unfired only
/// because in both arc fixtures a chord or a segment sits between any two arcs. `Wall`'s derived
/// equality killed it structurally: arcs of different circles, or of one circle in different
/// directions, now compare unequal. The one pair still equal — two *same-direction* arcs of one
/// circle — is the pair for which "no turn" is geometrically true on this face, so equality
/// answers right there too; a crossing at such a vertex is kept by the other faces' rings
/// (degree > 2), the function's own rule. No population produces that consecutive pair today.
fn dissolve_straight_angles(out: &mut [LocalFace], which: &[usize]) {
    // ★ The walls are per `(node, face plane)`. Globally they cannot be: the two result faces
    // that share a 3D edge each ride *the other's* plane as their wall, so a node in the middle of
    // that edge always sees two walls overall and would never dissolve. On each face separately it
    // sees one, which is exactly "the ring runs straight through here".
    // Per node: its neighbours, and whether any face bends there. `first_wall` remembers one wall
    // per `(node, face)` and `bent` records the first disagreement — flat maps, because a nested
    // one per node costs more than the whole pass is worth (measured: 4% of a fold).
    let mut nbrs: HashMap<NodeId, HashSet<NodeId>> = HashMap::new();
    let mut first_wall: HashMap<(NodeId, usize), Wall> = HashMap::new();
    let mut bent: HashSet<NodeId> = HashSet::new();
    for &fi in which {
        let lf = &out[fi];
        for ring in lf.poly_rings() {
            let k = ring.len();
            for i in 0..k {
                let (a, b) = (ring.nodes[i], ring.nodes[(i + 1) % k]);
                nbrs.entry(a).or_default().insert(b);
                nbrs.entry(b).or_default().insert(a);
                for nd in [a, b] {
                    match first_wall.entry((nd, lf.surf.plane())) {
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
    let mut drop: HashSet<NodeId> = HashSet::new();
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
