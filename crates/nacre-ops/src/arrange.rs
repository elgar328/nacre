//! Per-face planar arrangement for the polyhedral boolean (design §8 M5, sub-unit 3).
//!
//! This module answers combinatorial questions about how one solid's boundary cuts
//! a face of the other. It lives in `nacre-ops` and not in `nacre-geom` because it
//! needs `Model`/`Face`/`Edge`, and geom sits below topo (design §1: dependencies
//! flow upward only).
//!
//! Everything decided here is decided by an exact predicate. Coordinates that
//! appear (`three_planes`' cache) are never the basis of a decision — the truth of
//! a seam point is its plane triple, as it is for `Origin::Discovered` (design §4).

// Nothing in the boolean calls this yet — cell 3c replaces `reconstruct_face` with
// the arrangement and wires it in. Until then only tests exercise it, so a non-test
// build rightly sees it as unreachable.
#![cfg_attr(not(test), allow(dead_code))]

use crate::{
    BoolError, PlaneInfo, edge_incidence, face_half_edges, reject, solid_shell_handles, tag,
};
use nacre_geom::Surface;
use nacre_geom::intersect::{
    plane_pair_dir_sign, plane_side, three_plane_cmp_coord, three_plane_orient3d, three_planes,
};
use nacre_math::Point3;
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Loop, Model, Orientation, Solid, Vertex};
use std::cmp::Ordering;
use std::collections::HashMap;

/// A solid's edges, each with its endpoints and the planes of its two faces.
pub(crate) type EdgePlanes = HashMap<Handle<Edge>, ([Handle<Vertex>; 2], [usize; 2])>;

/// Where the seam line `P ∩ Q` enters and leaves the material common to face `f`
/// (plane `P`) and face `g` (plane `Q`).
///
/// `ends` are the truth: each is the sorted triple of planes whose exact meet is
/// that endpoint — the same identity `SeamVertex` uses. `points` are f64 caches of
/// those meets, kept for tolerance bookkeeping and debugging; no decision in this
/// module reads them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SeamSegment {
    pub(crate) plane_pair: [usize; 2],
    pub(crate) ends: [[usize; 3]; 2],
    pub(crate) points: [Point3; 2],
    /// The edge of `∂f` each end lies on, recorded where the crossing is produced.
    /// `Some` for a crossing of `∂f` with `g` (a boundary node); `None` for a
    /// crossing of `f` with `∂g` (an interior node, where a `g`-edge pierces `f`).
    ///
    /// Never back-derive this from the endpoint's triple `{P, Q, R}`: a non-convex
    /// face can carry two collinear edges on one neighbour plane `R`, and the seam
    /// would silently attach to the wrong one.
    ///
    /// An interior node loses no information by holding `None`. It is a boundary
    /// node of the *other* solid's face `g`, and running this on `g` records the
    /// `g`-edge there. Each face records only the edges of its own boundary.
    pub(crate) on_edge: [Option<Handle<Edge>>; 2],
}

/// One node of an assembled seam path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SeamEnd {
    pub(crate) triple: [usize; 3],
    /// f64 cache of the triple's exact meet. No decision reads it.
    pub(crate) point: Point3,
    pub(crate) on_edge: Option<Handle<Edge>>,
}

/// The seam of one solid on a face `f` of the other, decomposed into its components.
///
/// `Open` runs between two nodes of `∂f`; `Closed` never touches `∂f` — that is an
/// inner loop, a hole punched through `f`. `Closed`'s node list does **not** repeat
/// its first node; the cycle closes implicitly.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum SeamPath {
    Open(Vec<SeamEnd>),
    Closed(Vec<SeamEnd>),
}

impl SeamPath {
    pub(crate) fn nodes(&self) -> &[SeamEnd] {
        match self {
            SeamPath::Open(n) | SeamPath::Closed(n) => n,
        }
    }
}

/// Index a solid's edges by handle — [`edge_incidence`] keyed for lookup.
///
/// Hole-ring edges are in here too: `edge_incidence` walks every loop of every
/// face, and rejects an edge whose incidence is not a pair.
pub(crate) fn edge_planes(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Surface>, usize>,
) -> Result<EdgePlanes, BoolError> {
    let mut out = EdgePlanes::new();
    for (eh, bounds, inc) in edge_incidence(model, solid, surf_ix)? {
        out.insert(eh, (bounds, inc));
    }
    Ok(out)
}

/// The seam segments that solid `other` carves on face `f`.
///
/// For each face `g` of `other` the seam line is `L = P ∩ Q` (`P` = `f`'s plane,
/// `Q` = `g`'s). **`L` is never built.** Its direction `d = n_P × n_Q` only ever
/// enters as a sign, and two parallel planes simply contribute no crossings.
///
/// A crossing of `L` with `∂f` lies on an edge of `f`, which lies on `P ∩ R` for the
/// neighbouring face's plane `R` — so the crossing is exactly `P ∩ Q ∩ R`, the same
/// species of three-plane point the seam machinery already builds. Whether it lies
/// *within* that edge and *inside* `g` is a single question, and [`edge_crosses_face`]
/// answers it exactly.
pub(crate) fn seam_segments_on(
    model: &Model,
    f: Handle<Face>,
    other: Handle<Solid>,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Surface>, usize>,
    inc_x: &EdgePlanes,
    inc_y: &EdgePlanes,
) -> Result<Vec<SeamSegment>, BoolError> {
    let p = surf_ix[&model.faces.get(f).surface];
    let f_rings = face_rings(model, f, p, inc_x)?;

    let mut out = Vec::new();
    for shell in solid_shell_handles(model, other) {
        for &g in &model.shells.get(shell).faces {
            let q = surf_ix[&model.faces.get(g).surface];
            let g_rings = face_rings(model, g, q, inc_y)?;

            // `∂(f ∩ g) ⊆ (∂f ∩ g) ∪ (f ∩ ∂g)`, and these two sweeps collect exactly
            // those parts: a boundary crossing that lands outside the other face never
            // enters. That is what lets the sorted crossings be paired off directly.
            //
            // Sweep 1 walks `∂f`, so its crossings are boundary nodes and carry the edge
            // they lie on. Sweep 2 walks `∂g`: interior nodes, no edge of `f`.
            //
            // `∂f` and `∂g` are every loop, holes included. A hole rim of `g` crossing `P`
            // inside `f` is a seam node like any other; skipping it left an odd crossing
            // count, and `ARRANGEMENT_DEGENERATE` cried in place of the honest guard.
            let mut third: Vec<(usize, Option<Handle<Edge>>)> = Vec::new();
            for (face, rings, inc, own, into, on_f) in [
                (f, &g_rings, inc_x, p, q, true),
                (g, &f_rings, inc_y, q, p, false),
            ] {
                for he in face_half_edges(model.faces.get(face)) {
                    let (bounds, pair) = inc[&he.edge];
                    let [pa, pb] = pair;
                    let r = if pa == own { pb } else { pa };
                    let (v0, v1) = (
                        model.vertices.get(bounds[0]).point,
                        model.vertices.get(bounds[1]).point,
                    );
                    if edge_crosses_face(planes, pair, v0, v1, into, rings)? {
                        third.push((r, on_f.then_some(he.edge)));
                    }
                }
            }
            if third.is_empty() {
                continue;
            }

            third.sort_by(|&(i, _), &(j, _)| match order_along(planes, p, q, i, j) {
                -1 => Ordering::Less,
                1 => Ordering::Greater,
                _ => Ordering::Equal,
            });
            for w in third.windows(2) {
                if order_along(planes, p, q, w[0].0, w[1].0) != -1 {
                    // Coincident crossings, or `L` through a vertex of a face: four planes
                    // meet at one point. A duplicate `R` (two collinear edges of a
                    // non-convex face on one neighbour plane) lands here too — `order` of a
                    // plane against itself is 0 — which is what keeps two distinct nodes
                    // from ever collapsing onto one triple.
                    return Err(reject(tag::FOURPLANE));
                }
            }
            if third.len() % 2 != 0 {
                // Crossings alternate enter/exit and both ends of `L` lie outside, so an
                // odd count means one was missed: `L` grazed a vertex or ran along an edge.
                // Never guess where the material is.
                return Err(reject(tag::ARRANGEMENT_DEGENERATE));
            }

            for w in third.chunks_exact(2) {
                let ((r0, e0), (r1, e1)) = (w[0], w[1]);
                let point_of = |r: usize| {
                    three_planes(&planes[p].plane, &planes[q].plane, &planes[r].plane)
                        .ok_or_else(|| reject(tag::THREE_PLANES))
                };
                let seg = SeamSegment {
                    plane_pair: [p, q],
                    ends: [triple(p, q, r0), triple(p, q, r1)],
                    points: [point_of(r0)?, point_of(r1)?],
                    on_edge: [e0, e1],
                };
                debug_assert_eq!(
                    shared_planes(&seg.ends),
                    sorted_pair(p, q),
                    "a segment's endpoint triples must share exactly P and Q"
                );
                out.push(seg);
            }
        }
    }
    Ok(out)
}

/// Order two crossings of the line `P ∩ Q` along that line: `-1` if `V_i` precedes
/// `V_j`, `+1` if it follows, `0` if they coincide.
///
/// With `V_k = P ∩ Q ∩ R_k` and `d = n_P × n_Q`, we want `sign((V_i − V_j)·d)`. Since
/// `V_j ∈ R_j`, [`three_plane_orient3d`] gives `sign((V_i − V_j)·N_j)` for `N_j` the
/// right-hand normal of `R_j.tri`; multiplying by `sign(d·N_j)` recovers the order.
/// Both factors are exact predicates, so the comparator is a true total order.
///
/// The pair `(P, Q)` is a parameter, not the seam pair: sub-unit 3d orders two seam
/// crossings along an *edge* of `f` by calling this with `(P, R)`, the edge's own
/// two planes. No new predicate is needed for that.
fn order_along(planes: &[PlaneInfo], p: usize, q: usize, i: usize, j: usize) -> i8 {
    three_plane_orient3d(
        &planes[p].plane,
        &planes[q].plane,
        &planes[i].plane,
        planes[j].tri[0],
        planes[j].tri[1],
        planes[j].tri[2],
    ) * dir_sign(planes, p, q, j)
}

/// `+1` when a plane's stored normal already points out of its solid, `-1` when the
/// face is `Reversed` and the two oppose.
fn orient_sign(planes: &[PlaneInfo], i: usize) -> i8 {
    let dot = planes[i].plane.normal().dot(planes[i].n_out);
    debug_assert!(
        dot.abs() > 0.5,
        "a plane's normal must be parallel to n_out"
    );
    debug_assert_eq!(
        dot > 0.0,
        planes[i].orient == Orientation::Forward,
        "n_out's sign against the surface normal is the face's orientation"
    );
    if dot > 0.0 { 1 } else { -1 }
}

/// Group a face's seam arcs into rings, by index alone.
///
/// Each arc has one end on a kept→dropped transition (`kd`) and one on a dropped→kept
/// transition (`dk`) — the arc separates `f` in two, and the two pieces of `∂f` it lands
/// between belong to opposite sides, so `∂f` flips class at one end and back at the other.
/// Walking `∂f` forward from `dk[a] + 1` while the positions are kept therefore arrives at
/// some `kd[b]`, and `b` is the next arc of `a`'s ring. That successor map is a
/// permutation of the arcs; its cycles are the face's kept regions, one `LocalFace` each.
///
/// With one arc this is the old splice verbatim: the single cycle is `[0]`.
///
/// **Only integers go in and out**, and it does not know what they count. Cell 3e-2 fed it
/// `∂f`'s vertices, one per edge. Cell 3e-3 feeds it [`BoundaryRuns`] instead, because an
/// edge may carry two crossings — and there the classes alternate, so the walk above is
/// always exactly one step. The function is the more general of the two; its goldens say
/// so, and its callers no longer need that generality.
///
/// `next` is the successor over the crossing/run index space: `next[i]` is the position one
/// step forward along `∂f`. For a single ring it is `(i + 1) % n`; when `∂f` is several rings
/// each occupies a contiguous block and `next` cycles within it, so a run's walk never leaves
/// its ring. An arc bridges two rings by having its ends on both, and `succ` following that
/// arc is what threads a cycle across them — the index space alone carries the topology.
///
/// The caller turns a cycle into a ring by emitting, for each arc `b` in it, the boundary
/// positions `next[dk[prev]] ..= kd[b]` (walked with `next`) and then `b`'s nodes.
///
/// Three things are checked rather than assumed, because they are exactly where
/// `classof`'s ray casting and the arrangement's exact crossings would disagree — the same
/// two machines `SEAM_COUNT_MISMATCH` already arbitrates, so it is the same rejection:
/// as many `kd` transitions as arcs, no two arcs claiming one `kd`, and every arc used
/// exactly once. The walk is bounded, so a broken successor map rejects rather than hangs.
pub(crate) fn stitch_cycles(
    kept: &[bool],
    kd: &[usize],
    dk: &[usize],
    next: &[usize],
) -> Result<Vec<Vec<usize>>, BoolError> {
    let n = kept.len();
    let m = kd.len();
    let bad = || reject(tag::SEAM_COUNT_MISMATCH);
    if m == 0 || dk.len() != m || next.len() != n {
        return Err(bad());
    }
    // `kd`s and `dk`s are the two halves of the transition set, one per arc.
    if (0..n).filter(|&i| kept[i] && !kept[next[i]]).count() != m {
        return Err(bad());
    }
    let mut arc_at_kd: HashMap<usize, usize> = HashMap::new();
    for (a, &t) in kd.iter().enumerate() {
        if t >= n || !kept[t] || kept[next[t]] || arc_at_kd.insert(t, a).is_some() {
            return Err(bad());
        }
    }
    for &t in dk {
        if t >= n || kept[t] || !kept[next[t]] {
            return Err(bad());
        }
    }

    // The successor of arc `a`: walk the kept run that starts just past `dk[a]`. `next` is the
    // per-ring cyclic successor over the crossing/run index space, so the walk stays on the
    // ring `dk[a]` lies on, and an arc whose ends sit on two different rings is what carries
    // the walk across — no ring index is named, the successor map already encodes it.
    let succ = |a: usize| -> Result<usize, BoolError> {
        let mut i = next[dk[a]];
        for _ in 0..n {
            if kept[i] && !kept[next[i]] {
                return arc_at_kd.get(&i).copied().ok_or_else(bad);
            }
            if !kept[i] {
                return Err(bad()); // the run died before reaching a `kd`
            }
            i = next[i];
        }
        Err(bad())
    };

    let mut used = vec![false; m];
    let mut cycles = Vec::new();
    // Start at the lowest unused arc: `opens` is already in a deterministic order and
    // replay rests on the face order this fixes.
    for start in 0..m {
        if used[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut a = start;
        for _ in 0..=m {
            if used[a] {
                return Err(bad()); // re-entered a used arc: not a permutation
            }
            used[a] = true;
            cycle.push(a);
            a = succ(a)?;
            if a == start {
                break;
            }
        }
        if a != start {
            return Err(bad());
        }
        cycles.push(cycle);
    }
    Ok(cycles)
}

/// The two planes an ordered ring's edge `i → i+1` shares beyond `P`, plus the two nodes'
/// third planes — everything [`order_along`] needs to say which way that edge runs.
fn ring_edge(p: usize, ring: &[[usize; 3]], i: usize) -> Result<(usize, usize, usize), BoolError> {
    let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
    let shared: Vec<usize> = a.iter().copied().filter(|x| b.contains(x)).collect();
    if shared.len() != 2 || !shared.contains(&p) {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    let q = shared[usize::from(shared[0] == p)];
    let third = |t: [usize; 3]| t.iter().copied().find(|&x| x != p && x != q);
    let (Some(ri), Some(rj)) = (third(a), third(b)) else {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    };
    Ok((q, ri, rj))
}

/// `+1` when the ring's edge `i → i+1` runs along `d = n_P × n_Q`, `-1` against it.
fn edge_sign(
    planes: &[PlaneInfo],
    p: usize,
    ring: &[[usize; 3]],
    i: usize,
) -> Result<i8, BoolError> {
    let (q, ri, rj) = ring_edge(p, ring, i)?;
    // `order_along` is `sign((V_i − V_j)·d)`, so `-1` — `V_i` precedes `V_j` — is the edge
    // running along `+d`. Recomputed rather than remembered from the assembly walk.
    match order_along(planes, p, q, ri, rj) {
        -1 => Ok(1),
        1 => Ok(-1),
        _ => Err(reject(tag::LOOP_ORIENT_MISMATCH)), // two nodes coincide
    }
}

/// `∂f` cut into runs by the seam's boundary crossings.
///
/// A crossing always flips the kept/dropped class, so the runs alternate — which is what
/// lets the crossings, rather than the edges, index the splice. An edge may carry more
/// than one crossing (the other solid's boundary enters and leaves through it), and then
/// the run between them holds **no vertex at all**. That is the case an edge-indexed
/// model cannot express, and the whole reason this exists.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BoundaryRuns {
    /// Crossing `j`'s plane triple, in walk order around `∂f`.
    pub(crate) crossings: Vec<[usize; 3]>,
    /// The original vertices strictly inside run `j` (crossing `j` → crossing `j+1`), as
    /// indices into `f.outer.half_edges`. Possibly empty.
    pub(crate) runs: Vec<Vec<usize>>,
}

/// Assemble [`BoundaryRuns`] from the crossings grouped by the edge each rides.
///
/// `bnd` is `f`'s outer-vertex triples ([`face_vertex_triples`]); it is read only for
/// edges carrying more than one crossing, so the caller may pass an empty slice when no
/// edge does — building it rejects a straight angle, and there is no reason to run that
/// on a face this cell does not need it for.
pub(crate) fn boundary_runs(
    planes: &[PlaneInfo],
    p: usize,
    bnd: &[[usize; 3]],
    by_edge: &[Vec<[usize; 3]>],
) -> Result<BoundaryRuns, BoolError> {
    let mut crossings: Vec<[usize; 3]> = Vec::new();
    let mut runs: Vec<Vec<usize>> = Vec::new();
    // Vertices seen before the first crossing belong to the last run, which wraps.
    let mut pre: Vec<usize> = Vec::new();
    for (i, list) in by_edge.iter().enumerate() {
        match runs.last_mut() {
            Some(r) => r.push(i),
            None => pre.push(i),
        }
        for t in ordered_on_edge(planes, p, bnd, list, i)? {
            crossings.push(t);
            runs.push(Vec::new());
        }
    }
    let Some(last) = runs.last_mut() else {
        return Err(reject(tag::SEAM_COUNT_MISMATCH)); // an arc, but no boundary crossing
    };
    last.extend(pre);
    Ok(BoundaryRuns { crossings, runs })
}

/// The crossings on edge `i`, in the order the boundary walk meets them.
///
/// They all lie on `P ∩ R_i`, so [`order_along`] with `(P, R_i)` sorts them — the use its
/// doc reserved. The walk's direction along that line is [`edge_sign`], because the edge's
/// own endpoints are three-plane points on it too.
fn ordered_on_edge(
    planes: &[PlaneInfo],
    p: usize,
    bnd: &[[usize; 3]],
    list: &[[usize; 3]],
    i: usize,
) -> Result<Vec<[usize; 3]>, BoolError> {
    if list.len() < 2 {
        return Ok(list.to_vec());
    }
    let (r, _, _) = ring_edge(p, bnd, i)?;
    let q_of = |t: &[usize; 3]| -> Result<usize, BoolError> {
        let mut it = t.iter().copied().filter(|&x| x != p && x != r);
        match (it.next(), it.next()) {
            (Some(q), None) => Ok(q),
            // The node does not ride this edge's line: the two machines disagree.
            _ => Err(reject(tag::SEAM_COUNT_MISMATCH)),
        }
    };
    let qs: Vec<usize> = list.iter().map(q_of).collect::<Result<_, _>>()?;
    let mut ord: Vec<usize> = (0..list.len()).collect();
    // Insertion sort: the lists are two long in practice, and `?` stays available.
    for a in 1..ord.len() {
        let mut j = a;
        while j > 0 && order_along(planes, p, r, qs[ord[j]], qs[ord[j - 1]]) == -1 {
            ord.swap(j, j - 1);
            j -= 1;
        }
    }
    for w in ord.windows(2) {
        if order_along(planes, p, r, qs[w[0]], qs[w[1]]) != -1 {
            return Err(reject(tag::FOURPLANE)); // two crossings at one point on the edge
        }
    }
    if edge_sign(planes, p, bnd, i)? < 0 {
        ord.reverse();
    }
    Ok(ord.into_iter().map(|k| list[k]).collect())
}

/// The kept/dropped class of every run, from two machines that must agree.
///
/// Alternation gives the shape: a crossing flips the class, so run `j+1` opposes run `j`.
/// `kept_vert` gives the anchor: a run holding an original vertex takes that vertex's
/// class. Seed from any such run, propagate, then check **every** vertex against the
/// propagation — `classof`'s ray casting against the arrangement's exact crossings, which
/// is what `SEAM_COUNT_MISMATCH` has always arbitrated.
///
/// A run with no vertex is decided by alternation alone. That is exactly the run between
/// two crossings of one edge, and it has no other source of truth.
///
/// `ring_lens` is the run count of each ring of `∂f`, in the concatenated order `runs` are
/// laid out. Alternation only wraps within a ring — a crossing on the outer loop says nothing
/// about a run on a hole rim — so each ring is seeded and propagated on its own. A single ring
/// passes `&[runs.len()]`. Each ring must have evenly many runs (its crossings alternate
/// around a closed loop) and at least one vertex to anchor it; a ring of vertex-free runs has
/// no source of truth and is rejected. The vertex check is global: every run, every ring.
pub(crate) fn run_classes(
    runs: &[Vec<usize>],
    kept_vert: &[bool],
    ring_lens: &[usize],
) -> Result<Vec<bool>, BoolError> {
    if ring_lens.iter().sum::<usize>() != runs.len() {
        return Err(reject(tag::SEAM_COUNT_MISMATCH));
    }
    let mut kept = vec![false; runs.len()];
    let mut base = 0;
    for &len in ring_lens {
        // Crossings alternate enter/exit around a closed ring, so there are evenly many.
        if len == 0 || len % 2 != 0 {
            return Err(reject(tag::SEAM_COUNT_MISMATCH));
        }
        let ring = &runs[base..base + len];
        let seed = (0..len)
            .find(|&j| !ring[j].is_empty())
            .ok_or_else(|| reject(tag::SEAM_COUNT_MISMATCH))?;
        let s = kept_vert[ring[seed][0]];
        for j in 0..len {
            kept[base + (seed + j) % len] = if j % 2 == 0 { s } else { !s };
        }
        base += len;
    }
    for (j, run) in runs.iter().enumerate() {
        if run.iter().any(|&v| kept_vert[v] != kept[j]) {
            return Err(reject(tag::SEAM_COUNT_MISMATCH));
        }
    }
    Ok(kept)
}

/// The turn at ring node `i`, about the face's **outward** normal: `+1` left, `-1` right.
///
/// **No point is materialized, and no coordinate is read.** The incoming edge runs along
/// `s_a·(n_P × n_A)` and the outgoing along `s_b·(n_P × n_B)`, where `A` and `B` are node
/// `i`'s two non-`P` planes. Their cross product, dotted with `n_P`:
///
/// ```text
///   (n_P × n_A) × (n_P × n_B) = n_P · det[n_P, n_A, n_B]      (a×b)×(a×c) = a·det(a,b,c)
///     ⇒  turn = s_a · s_b · sign(det[n_P, n_A, n_B]) · orient_sign(P)
/// ```
///
/// and `sign(det[…])` is [`plane_pair_dir_sign`], already exact. It is never `0`: node `i`
/// lies on all three planes, and a point exists there only if their normals are independent.
///
/// A ring is not convex, so this is **not** the winding — at a reflex node it is its
/// opposite. [`loop_winding`] asks it at a hull vertex, where the two agree.
pub(crate) fn turn_at(
    planes: &[PlaneInfo],
    p: usize,
    ring: &[[usize; 3]],
    i: usize,
) -> Result<i8, BoolError> {
    let n = ring.len();
    let prev = (i + n - 1) % n;
    let (a, _, _) = ring_edge(p, ring, prev)?; // plane of the edge arriving at `i`
    let (b, _, _) = ring_edge(p, ring, i)?; // plane of the edge leaving `i`
    let sa = edge_sign(planes, p, ring, prev)?;
    let sb = edge_sign(planes, p, ring, i)?;
    let det = plane_pair_dir_sign(&planes[p].plane, &planes[a].plane, &planes[b].plane);
    if det == 0 {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    Ok(sa * sb * det * orient_sign(planes, p))
}

/// Face `f`'s outer-loop vertices as three-plane triples: `f`'s own plane, and the
/// neighbouring planes of the two edges meeting there.
///
/// An original vertex is as implicit a point as a seam node, so a cycle's ring — which mixes
/// them — is one uniform list and [`point_in_ring`] need not know the difference. Two
/// adjacent edges on one neighbour plane would be a straight angle, and it rejects.
pub(crate) fn face_vertex_triples(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgePlanes,
) -> Result<Vec<[usize; 3]>, BoolError> {
    loop_triples(&model.faces.get(f).outer, p, inc)
}

/// Each hole ring of face `f`, as three-plane triples.
///
/// A rim edge's incidence is `[p, wall]`, so the neighbour plane is the wall on the other
/// side of the rim — the same construction as an outer vertex, and the same rejection of a
/// straight angle. The ring keeps its stored direction: clockwise about `f`'s outward
/// normal, which is what makes it a hole.
pub(crate) fn hole_rings(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgePlanes,
) -> Result<Vec<Vec<[usize; 3]>>, BoolError> {
    model
        .faces
        .get(f)
        .inner
        .iter()
        .map(|l| loop_triples(l, p, inc))
        .collect()
}

fn loop_triples(l: &Loop, p: usize, inc: &EdgePlanes) -> Result<Vec<[usize; 3]>, BoolError> {
    let hes = &l.half_edges;
    let other = |he: &nacre_topo::HalfEdge| -> Result<usize, BoolError> {
        let (_, [pa, pb]) = *inc.get(&he.edge).ok_or_else(|| reject(tag::MISSING_SEAM))?;
        Ok(if pa == p { pb } else { pa })
    };
    let n = hes.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        // Vertex `i` starts edge `i` and ends edge `i - 1`.
        let (a, b) = (other(&hes[(i + n - 1) % n])?, other(&hes[i])?);
        if a == b {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH)); // a straight angle
        }
        let mut t = [p, a, b];
        t.sort_unstable();
        out.push(t);
    }
    Ok(out)
}

/// The exact side of plane `q` that the implicit point `t` lies on: `0` means *on* it.
fn side_of(planes: &[PlaneInfo], t: [usize; 3], q: usize) -> i8 {
    three_plane_orient3d(
        &planes[t[0]].plane,
        &planes[t[1]].plane,
        &planes[t[2]].plane,
        planes[q].tri[0],
        planes[q].tri[1],
        planes[q].tri[2],
    )
}

/// Is the implicit point `v` inside the simple ring `ring`, both on face plane `p`?
///
/// **A ray, cast along a line we already have.** Every ring edge lies on `P ∩ R`, and `v`
/// lies on `P ∩ Q_a` for either of its own two planes. Those two lines meet at
/// `X = {P, Q_a, R}`, which is *itself* a three-plane point — so "is `X` inside the edge"
/// and "is `X` ahead of `v`" are both [`order_along`], the comparator cell 3d already built
/// for two three-plane points on one line. **No coordinate is read and no point is built.**
///
/// **Choosing the ray first deletes the special cases.** A ring node on the ray's *line*
/// would make the parity ambiguous; `side_of` decides that exactly. Once no node lies on the
/// line, an edge's line cannot *be* the ray's line (its endpoints would be on it), so a zero
/// determinant always means "parallel and distinct" — no crossing, no collinearity to handle
/// — and every crossing is transversal, so parity is containment.
///
/// Candidates are each node's two non-`P` planes, in ring order, `+d` before `-d`; the first
/// clear one wins, which keeps the answer deterministic. `no_clear_ray` if none is clear.
/// The answer must not depend on which was chosen, and a golden says so.
///
/// `v` must not lie *on* `ring`. That is `SeamPath::Closed`'s standing claim — a closed seam
/// loop never touches `∂f` — and this is where it is finally checked: an intersection at
/// `X == v` strictly inside an edge is `point_on_ring`.
pub(crate) fn point_in_ring(
    planes: &[PlaneInfo],
    p: usize,
    v: [usize; 3],
    ring: &[[usize; 3]],
) -> Result<bool, BoolError> {
    every_ray(planes, p, v, ring)?
        .first()
        .copied()
        .ok_or_else(|| reject(tag::NO_CLEAR_RAY))
}

/// Is the implicit point `v` **on** the simple ring `ring` — on an edge, endpoints included?
///
/// [`point_in_ring`] cannot answer this. It casts a ray along `P ∩ Q_a` for each of `v`'s own
/// planes, and skips a line that carries a ring node; when `v` *is* a ring node, both of its
/// lines carry it and every candidate is skipped, so the honest `point_on_ring` comes back as
/// `no_clear_ray`. Asking first, and separately, is what keeps the tag truthful.
///
/// The question is two exact signs per edge. `v` lies on edge `i`'s line `P ∩ R` iff `v` lies
/// on `R` ([`side_of`]); it lies within the edge iff it does not fall on the same side of both
/// endpoints along that line ([`order_along`]). No new predicate, and no coordinate.
pub(crate) fn point_on_ring(
    planes: &[PlaneInfo],
    p: usize,
    v: [usize; 3],
    ring: &[[usize; 3]],
) -> Result<bool, BoolError> {
    if ring.len() < 3 {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    for i in 0..ring.len() {
        let (r, si, sj) = ring_edge(p, ring, i)?;
        if side_of(planes, v, r) != 0 {
            continue; // `v` is not even on the edge's line
        }
        // Name `v` as a point of that line: `{P, R, S}` for one of its own planes `S`. Both
        // cannot fail — if neither `E0` nor `E1` is independent of `P, R`, then `E0` and `E1`
        // lie in `span(n_P, n_R)` and `v = {P, E0, E1}` would not have been a point.
        let s = *v
            .iter()
            .find(|&&x| {
                x != p
                    && plane_pair_dir_sign(&planes[p].plane, &planes[r].plane, &planes[x].plane)
                        != 0
            })
            .ok_or_else(|| reject(tag::LOOP_ORIENT_MISMATCH))?;
        let (a, b) = (
            order_along(planes, p, r, s, si),
            order_along(planes, p, r, s, sj),
        );
        if a * b <= 0 {
            return Ok(true); // between the endpoints, or on one of them
        }
    }
    Ok(false)
}

/// Face `f`'s rings as three-plane triples: its outer loop, then each hole.
///
/// The triple form of [`face_loops`], and what [`edge_crosses_face`] needs. Building it costs
/// the straight-angle rejection of [`face_vertex_triples`] — a face with two adjacent edges on
/// one neighbour plane has no triple for the vertex between them.
pub(crate) fn face_rings(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgePlanes,
) -> Result<Vec<Vec<[usize; 3]>>, BoolError> {
    let mut out = vec![face_vertex_triples(model, f, p, inc)?];
    out.extend(hole_rings(model, f, p, inc)?);
    Ok(out)
}

/// Does the edge on planes `[e0, e1]`, running `p0 → p1`, pierce the **material** of face `g`
/// (plane `q`, `rings` outer-first)?
///
/// **The piercing point is a three-plane point.** The edge lies on `E0` and `E1`, the face on
/// `Q`, so the point is `{E0, E1, Q}` — the same species the seam machinery already builds,
/// and [`point_in_ring`] already answers containment for. The fan that stood here instead cut
/// the face into triangles and asked which one the point fell in; its own diagonals then
/// grazed the point whenever the geometry was symmetric, and a square's centre lies on the
/// diagonal of **every** apex. Design §9 guessed the fix would be a real triangulation, or a
/// choice of apex. It is neither: no triangulation is free of diagonals.
///
/// **One coordinate is read, and it is exact.** Whether `p0` and `p1` straddle `Q` is
/// [`plane_side`], an exact `orient3d` on the stored points. Everything after is combinatorial:
/// straddling means the segment meets `Q` exactly once, so the crossing is inside the edge by
/// construction, and containment is a ray cast along a line the planes already give.
///
/// **Contacts are named, not lumped.** An endpoint on `Q` — or the whole edge lying in it — is
/// `vertex_on_face_plane`. A crossing exactly on `∂g` is `point_on_ring`. Neither is a graze:
/// they are the two ways an edge can touch a face without properly piercing it.
pub(crate) fn edge_crosses_face(
    planes: &[PlaneInfo],
    [e0, e1]: [usize; 2],
    p0: Point3,
    p1: Point3,
    q: usize,
    rings: &[Vec<[usize; 3]>],
) -> Result<bool, BoolError> {
    let (s0, s1) = (plane_side(planes[q].tri, p0), plane_side(planes[q].tri, p1));
    if s0 == 0 || s1 == 0 {
        return Err(reject(tag::VERTEX_ON_FACE_PLANE));
    }
    if s0 == s1 {
        return Ok(false); // the segment never reaches `Q`
    }
    let x = triple(e0, e1, q);
    for ring in rings {
        if point_on_ring(planes, q, x, ring)? {
            return Err(reject(tag::POINT_ON_RING));
        }
    }
    if !point_in_ring(planes, q, x, &rings[0])? {
        return Ok(false);
    }
    for hole in &rings[1..] {
        if point_in_ring(planes, q, x, hole)? {
            return Ok(false); // through the hole, not the material
        }
    }
    Ok(true)
}

/// The parity every clear ray reports. The ring is simple, so they must all agree; a golden
/// says so, which is a second machine for free.
pub(crate) fn every_ray(
    planes: &[PlaneInfo],
    p: usize,
    v: [usize; 3],
    ring: &[[usize; 3]],
) -> Result<Vec<bool>, BoolError> {
    if ring.len() < 3 {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    let mut out = Vec::new();
    for &qa in v.iter().filter(|&&x| x != p) {
        // Clear iff no ring node sits on `Q_a`, hence none on the line `P ∩ Q_a`.
        if ring.iter().any(|&m| side_of(planes, m, qa) == 0) {
            continue;
        }
        let qb = *v
            .iter()
            .find(|&&x| x != p && x != qa)
            .ok_or_else(|| reject(tag::LOOP_ORIENT_MISMATCH))?;
        for dir in [1i8, -1] {
            let mut crossings = 0usize;
            for i in 0..ring.len() {
                let (r, si, sj) = ring_edge(p, ring, i)?;
                // Parallel: distinct lines, because no ring node lies on `P ∩ Q_a`.
                if plane_pair_dir_sign(&planes[p].plane, &planes[qa].plane, &planes[r].plane) == 0 {
                    continue;
                }
                // `X = {P, Q_a, R}` strictly inside the edge?
                let (a, b) = (
                    order_along(planes, p, r, qa, si),
                    order_along(planes, p, r, qa, sj),
                );
                if a * b >= 0 {
                    continue; // outside the edge, or on an endpoint (excluded above)
                }
                // Strictly ahead of `v` along `dir · (n_P × n_Qa)`?
                match order_along(planes, p, qa, r, qb) {
                    0 => return Err(reject(tag::POINT_ON_RING)), // `X == v`, inside an edge
                    o if o == dir => crossings += 1,
                    _ => {}
                }
            }
            out.push(crossings % 2 == 1);
        }
    }
    Ok(out)
}

/// An ordered ring's winding about the face's outward normal: `-1` clockwise — the material
/// is *outside* the ring, so it bounds a hole — and `+1` counter-clockwise, an island.
///
/// The turn at a convex-hull vertex is the winding, and the lexicographically smallest node
/// is one: it is an extreme point of the node set, which is planar, so it is a vertex of the
/// ring's hull. Finding it is the **only** thing here that needs two implicit points in one
/// decision, and [`three_plane_cmp_coord`] is that predicate.
///
/// A shortcut dies here, and is recorded so it is not walked twice: a *supporting edge* —
/// one whose plane `Q_j` has every other node on one side — would give a hull vertex from
/// the one-implicit `three_plane_orient3d` alone. But a simple polygon need not have an edge
/// on its hull (fold each side of a pentagon slightly inward), so no such edge is guaranteed.
/// A hull *vertex* always exists.
///
/// Two distinct triples that compare equal on all three axes are the same point, and the
/// loop is not simple. That is `LOOP_ORIENT_MISMATCH`, decided by exact equality rather than
/// by a tolerance — which is why design.md §9 could not check it before.
pub(crate) fn loop_winding(
    planes: &[PlaneInfo],
    p: usize,
    ring: &[[usize; 3]],
) -> Result<i8, BoolError> {
    if ring.len() < 3 {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    let tri = |t: [usize; 3]| {
        [
            &planes[t[0]].plane,
            &planes[t[1]].plane,
            &planes[t[2]].plane,
        ]
    };
    let mut lo = 0usize;
    for i in 1..ring.len() {
        let ord = (0..3)
            .map(|axis| three_plane_cmp_coord(tri(ring[i]), tri(ring[lo]), axis))
            .find(|&c| c != 0);
        match ord {
            Some(-1) => lo = i,
            Some(_) => {}
            None => return Err(reject(tag::LOOP_ORIENT_MISMATCH)), // two nodes coincide
        }
    }
    turn_at(planes, p, ring, lo)
}

/// A closed seam loop's triples, ordered so that every directed edge keeps the kept
/// material on its left — the convention **every** b-rep loop obeys, outer or inner.
/// So the caller wraps them in `Node::Seam` and needs no further knowledge; returning
/// a "should I reverse?" flag would leave it an obligation it could forget.
///
/// **No coordinate is read.** Every loop of a face winds CCW about that face's outward
/// normal, and for a planar ring the interior of a CCW loop lies along `n × t`. So with
/// `f` on plane `P` and the loop's edge on `f ∩ g` (`g` a face of the other solid),
/// wanting the material — the `+n_out_g` side when `keep == Outside` — on the left:
///
/// ```text
///   n_out_f × t ∝ +n_out_g   ⇒   t = n_out_g × n_out_f = −(keep)·s_f·s_g·d
/// ```
///
/// where `n_out_f = s_f·n_P`, `n_out_g = s_g·n_Q`, `d = n_P × n_Q`, and `keep = ±1` for
/// `Outside`/`Inside`. (Expand: `n_f × (n_g × n_f) = n_g − (n_f·n_g) n_f`, the in-plane
/// part of `n_out_g`, a positive multiple.) `d` is the direction [`seam_segments_on`]
/// already sorted its crossings along, so [`order_along`] answers "does this edge run
/// along `+d`?" exactly. Three signs, no area.
///
/// Every edge must give the same verdict, the wrap-around one included — a loop that
/// keeps material on its left does so everywhere. `order_along` is *recomputed* rather
/// than remembered from the assembly walk: the walk's memory and the exact predicate
/// are the two machines, and a check that reuses one machine's memory checks nothing.
pub(crate) fn orient_seam_loop(
    planes: &[PlaneInfo],
    p: usize,
    nodes: &[SeamEnd],
    material_outside: bool,
) -> Result<Vec<[usize; 3]>, BoolError> {
    let n = nodes.len();
    if n < 3 {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    let keep: i8 = if material_outside { 1 } else { -1 };
    let mut verdict: Option<bool> = None;

    for i in 0..n {
        let (a, b) = (nodes[i].triple, nodes[(i + 1) % n].triple);
        let shared: Vec<usize> = a.iter().copied().filter(|x| b.contains(x)).collect();
        // `P` plus exactly one `Q`: adjacent loop nodes bound one seam segment.
        if shared.len() != 2 || !shared.contains(&p) {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH));
        }
        let q = shared[usize::from(shared[0] == p)];
        let third = |t: [usize; 3]| t.iter().copied().find(|&x| x != p && x != q);
        let (Some(ri), Some(rj)) = (third(a), third(b)) else {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH));
        };

        let ord = order_along(planes, p, q, ri, rj);
        if ord == 0 {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH)); // two nodes coincide
        }
        let along_d = ord == -1;
        let want_along_d = -keep * orient_sign(planes, p) * orient_sign(planes, q) > 0;
        let reverse = along_d != want_along_d;
        match verdict {
            None => verdict = Some(reverse),
            Some(v) if v == reverse => {}
            Some(_) => return Err(reject(tag::LOOP_ORIENT_MISMATCH)),
        }
    }

    let it = nodes.iter().map(|nd| nd.triple);
    Ok(if verdict == Some(true) {
        it.rev().collect()
    } else {
        it.collect()
    })
}

/// The seam of `other` on face `f`, assembled into open arcs and closed loops.
///
/// `∂other ∩ P ∩ f` is a 1-manifold: the segments of [`seam_segments_on`] never cross,
/// they only meet at endpoints. So the assembly is pure combinatorics on the endpoint
/// triples — **no coordinate is ever read here**, and no convexity is assumed. That is
/// what makes the arc order *true* rather than a projection's guess.
///
/// Two nodes never collapse onto one triple, so the adjacency map is faithful:
///
/// * A boundary node's triple holds two `f`-solid planes, an interior node's two
///   `other`-solid planes — the two species cannot coincide.
/// * Crossings against different faces `g` name different `Q`.
/// * Within one `g`, a duplicate third plane compares `0` against itself and
///   [`seam_segments_on`] has already rejected it as `FOURPLANE`.
///
/// Degrees follow: a boundary node is produced by exactly one `g`, so degree 1; an
/// interior node sits on a pierced `other`-edge, whose two faces each produce it, so
/// degree 2. Anything else is `SEAM_BRANCH`.
///
/// **Determinism comes from the walk, not the map.** `HashMap` iteration order varies
/// per instance, so the map is only ever *queried*. Segments are produced in a
/// deterministic order (`Shell::faces` is a `Vec`, crossings are sorted by an exact
/// predicate), and paths start at the lowest-index unused segment — so the same model
/// yields the same paths, in the same order, walked in the same direction. Sub-unit
/// 3d's operation-log replay rests on this.
pub(crate) fn seam_paths_on(
    model: &Model,
    f: Handle<Face>,
    other: Handle<Solid>,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Surface>, usize>,
    inc_x: &EdgePlanes,
    inc_y: &EdgePlanes,
) -> Result<Vec<SeamPath>, BoolError> {
    let segs = seam_segments_on(model, f, other, planes, surf_ix, inc_x, inc_y)?;

    let mut adj: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        for e in s.ends {
            adj.entry(e).or_default().push(i);
        }
    }
    if segs.iter().any(|s| s.ends.iter().any(|e| adj[e].len() > 2)) {
        return Err(reject(tag::SEAM_BRANCH));
    }

    let node = |i: usize, k: usize| SeamEnd {
        triple: segs[i].ends[k],
        point: segs[i].points[k],
        on_edge: segs[i].on_edge[k],
    };
    // Walk from segment `i`, leaving by its end `k`; mark every segment consumed.
    let walk = |used: &mut Vec<bool>, i: usize, k: usize| -> Vec<SeamEnd> {
        let mut nodes = vec![node(i, k)];
        let (mut si, mut sk) = (i, k);
        loop {
            used[si] = true;
            let far = 1 - sk;
            nodes.push(node(si, far));
            let t = segs[si].ends[far];
            let Some(&next) = adj[&t].iter().find(|&&j| !used[j]) else {
                return nodes;
            };
            si = next;
            sk = usize::from(segs[next].ends[1] == t && segs[next].ends[0] != t);
        }
    };

    let mut used = vec![false; segs.len()];
    let mut out: Vec<SeamPath> = Vec::new();

    // Open arcs first: every one has exactly two degree-1 ends, and the lower-indexed
    // of its segments carries one of them — so index order fixes start and direction.
    for i in 0..segs.len() {
        for k in 0..2 {
            if !used[i] && adj[&segs[i].ends[k]].len() == 1 {
                out.push(SeamPath::Open(walk(&mut used, i, k)));
            }
        }
    }
    // What survives is all-degree-2: a disjoint union of cycles.
    for i in 0..segs.len() {
        if !used[i] {
            let mut nodes = walk(&mut used, i, 0);
            // `walk` re-emits the start triple when the cycle closes; drop it.
            debug_assert_eq!(nodes.last().map(|n| n.triple), Some(nodes[0].triple));
            nodes.pop();
            out.push(SeamPath::Closed(nodes));
        }
    }

    debug_assert!(out.iter().all(|p| match p {
        // `on_edge.is_some()` ⇔ boundary node ⇔ degree 1: an arc's two ends, and
        // nothing else. A closed loop touches `∂f` nowhere.
        SeamPath::Open(n) =>
            n.len() >= 2
                && n[1..n.len() - 1].iter().all(|e| e.on_edge.is_none())
                && n[0].on_edge.is_some()
                && n[n.len() - 1].on_edge.is_some(),
        SeamPath::Closed(n) => n.iter().all(|e| e.on_edge.is_none()),
    }));
    Ok(out)
}

/// `sign((n_P × n_Q) · N_R)`, where `N_R` is the right-hand normal of `R.tri`.
///
/// [`plane_pair_dir_sign`] gives the sign against `R`'s *stored* normal, exactly.
/// That normal is parallel to `N_R` but may oppose it on a `Reversed` face, so we
/// correct with their dot — two parallel unit vectors, `|·| ≈ 1`, nowhere near the
/// sign boundary.
///
/// It is tempting to read the correction off `PlaneInfo::orient` instead. Don't:
/// "`Reversed` ⇔ `n_out = −plane.normal()`" is an invariant nothing enforces, while
/// the predicate's convention is tied to `tri`'s RH normal by construction. Were the
/// invariant to break, an `orient`-based order would reverse silently. Assert the
/// agreement; do not depend on it.
fn dir_sign(planes: &[PlaneInfo], p: usize, q: usize, r: usize) -> i8 {
    plane_pair_dir_sign(&planes[p].plane, &planes[q].plane, &planes[r].plane)
        * orient_sign(planes, r)
}

fn triple(a: usize, b: usize, c: usize) -> [usize; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

fn sorted_pair(a: usize, b: usize) -> [usize; 2] {
    let mut t = [a, b];
    t.sort_unstable();
    t
}

/// The planes both endpoint triples name — `P` and `Q`, if the wiring is right.
fn shared_planes(ends: &[[usize; 3]; 2]) -> [usize; 2] {
    let both: Vec<usize> = ends[0]
        .iter()
        .filter(|i| ends[1].contains(i))
        .copied()
        .collect();
    debug_assert_eq!(both.len(), 2, "endpoints differ in exactly one plane");
    sorted_pair(both[0], both[1])
}
