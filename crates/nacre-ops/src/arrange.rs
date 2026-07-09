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

use crate::{BoolError, PlaneInfo, edge_incidence, face_loops, reject, segment_crosses_face, tag};
use nacre_geom::Surface;
use nacre_geom::intersect::{plane_pair_dir_sign, three_plane_orient3d, three_planes};
use nacre_math::Point3;
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Model, Orientation, Solid, Vertex};
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

/// Index a solid's edges by handle.
///
/// `edge_incidence` walks outer loops only, so a hole-ring edge is seen once — its
/// other use is on the lid's inner loop. Callers must have rejected holed operands
/// (`INNER_LOOP_OPERAND`); we check rather than index blindly, because indexing
/// `inc[1]` on such an edge is exactly what used to panic.
pub(crate) fn edge_planes(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Surface>, usize>,
) -> Result<EdgePlanes, BoolError> {
    let mut out = EdgePlanes::new();
    for (eh, bounds, inc) in edge_incidence(model, solid, surf_ix) {
        let [a, b] = match inc[..] {
            [a, b] => [a, b],
            _ => return Err(reject(tag::INNER_LOOP_OPERAND)),
        };
        out.insert(eh, (bounds, [a, b]));
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
/// *within* that edge and *inside* `g` is a single question, and `segment_crosses_face`
/// already answers it exactly.
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
    let f_rings = face_loops(model, f);
    let shell = model.solids.get(other).outer;

    let mut out = Vec::new();
    for &g in &model.shells.get(shell).faces {
        let q = surf_ix[&model.faces.get(g).surface];
        let g_rings = face_loops(model, g);

        // `∂(f ∩ g) ⊆ (∂f ∩ g) ∪ (f ∩ ∂g)`, and these two sweeps collect exactly
        // those parts: a boundary crossing that lands outside the other face never
        // enters. That is what lets the sorted crossings be paired off directly.
        //
        // Sweep 1 walks `∂f`, so its crossings are boundary nodes and carry the edge
        // they lie on. Sweep 2 walks `∂g`: interior nodes, no edge of `f`.
        let mut third: Vec<(usize, Option<Handle<Edge>>)> = Vec::new();
        for (face, rings, inc, own, on_f) in [
            (f, &g_rings, inc_x, p, true),
            (g, &f_rings, inc_y, q, false),
        ] {
            for he in &model.faces.get(face).outer.half_edges {
                let (bounds, [pa, pb]) = inc[&he.edge];
                let r = if pa == own { pb } else { pa };
                let (v0, v1) = (
                    model.vertices.get(bounds[0]).point,
                    model.vertices.get(bounds[1]).point,
                );
                if segment_crosses_face(v0, v1, rings)? {
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
pub(crate) fn orient_hole_loop(
    planes: &[PlaneInfo],
    p: usize,
    nodes: &[SeamEnd],
    material_outside: bool,
) -> Result<Vec<[usize; 3]>, BoolError> {
    let n = nodes.len();
    if n < 3 {
        return Err(reject(tag::HOLE_ORIENT_MISMATCH));
    }
    let keep: i8 = if material_outside { 1 } else { -1 };
    let mut verdict: Option<bool> = None;

    for i in 0..n {
        let (a, b) = (nodes[i].triple, nodes[(i + 1) % n].triple);
        let shared: Vec<usize> = a.iter().copied().filter(|x| b.contains(x)).collect();
        // `P` plus exactly one `Q`: adjacent loop nodes bound one seam segment.
        if shared.len() != 2 || !shared.contains(&p) {
            return Err(reject(tag::HOLE_ORIENT_MISMATCH));
        }
        let q = shared[usize::from(shared[0] == p)];
        let third = |t: [usize; 3]| t.iter().copied().find(|&x| x != p && x != q);
        let (Some(ri), Some(rj)) = (third(a), third(b)) else {
            return Err(reject(tag::HOLE_ORIENT_MISMATCH));
        };

        let ord = order_along(planes, p, q, ri, rj);
        if ord == 0 {
            return Err(reject(tag::HOLE_ORIENT_MISMATCH)); // two nodes coincide
        }
        let along_d = ord == -1;
        let want_along_d = -keep * orient_sign(planes, p) * orient_sign(planes, q) > 0;
        let reverse = along_d != want_along_d;
        match verdict {
            None => verdict = Some(reverse),
            Some(v) if v == reverse => {}
            Some(_) => return Err(reject(tag::HOLE_ORIENT_MISMATCH)),
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
