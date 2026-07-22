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
//!
//! # One `usize`, two meanings — the rule that keeps them apart
//!
//! `planes` is indexed **per face**, but an arrangement reasons **per plane**: an earlier boolean
//! can split one geometric plane between two faces (a base's exposed top and the cantilever
//! underside above it), and those are two `PlaneInfo` for one plane. When the two meanings meet in
//! one comparison the result is wrong *silently* — measured four times on this branch, most
//! recently as 117748 predicate calls that read "different plane" for two faces of one plane and
//! answered from rounding noise (2026-07-22). So:
//!
//! - **A plane triple's elements, any index compared for plane identity, and every exact-predicate
//!   argument are class roots** ([`class_of`]). Producers emit class form; consumers normalize on
//!   entry; [`crate::tolerant`]'s predicates assert it.
//! - **An index that reads *this face's* geometry is a face** — `n_out`, `orient`, `inc`
//!   (`EdgePlanes` incidences are faces, so `other()`/`pa == own` must match raw), and the
//!   `declined` log. Reading `n_out` off a class representative can flip a seated face inside out,
//!   because two faces of one plane may have **opposite** normals.
//! - **An index that states a *plane class's* label frame is a class root** — the `[*_above,
//!   *_below]` frame and `orient_sign(wc)` are defined about the class root's stored normal, by
//!   construction (see [`side_of`] for the trap that hides in reading it raw).
//! - **Exception:** the code that *defines* the classes (`crate::fill_classes` →
//!   `shares_or_coplanar`) runs before they exist and takes face indices.
//!
//! Name the two apart wherever both are in scope: `fp` for the face, `fc` for its class.

// Nothing in the boolean calls this yet — cell 3c replaces `reconstruct_face` with
// the arrangement and wires it in. Until then only tests exercise it, so a non-test
// build rightly sees it as unreachable.
#![cfg_attr(not(test), allow(dead_code))]

use crate::tolerant::{t_cmp_coord, t_orient3d, t_plane_pair_dir_sign};
use crate::{BoolError, PlaneInfo, edge_incidence, reject, tag};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Loop, Model, Orientation, Solid, Vertex};
use std::collections::HashMap;

/// A solid's edges, each with its endpoints and the planes of its two faces.
pub(crate) type EdgePlanes = HashMap<Handle<Edge>, ([Handle<Vertex>; 2], [usize; 2])>;

/// Index a solid's edges by handle — [`edge_incidence`] keyed for lookup.
///
/// Hole-ring edges are in here too: `edge_incidence` walks every loop of every
/// face, and rejects an edge whose incidence is not a pair.
pub(crate) fn edge_planes(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<EdgePlanes, BoolError> {
    let mut out = EdgePlanes::new();
    for (eh, bounds, inc) in edge_incidence(model, solid, surf_ix)? {
        out.insert(eh, (bounds, inc));
    }
    Ok(out)
}

/// Order two crossings of the line `P ∩ Q` along that line: `-1` if `V_i` precedes
/// `V_j`, `+1` if it follows, `0` if they coincide.
///
/// With `V_k = P ∩ Q ∩ R_k` and `d = n_P × n_Q`, we want `sign((V_i − V_j)·d)`. Since
/// `V_j ∈ R_j`, [`three_plane_orient3d`](nacre_geom::intersect::three_plane_orient3d) gives `sign((V_i − V_j)·N_j)` for `N_j` the
/// right-hand normal of `R_j.tri`; multiplying by `sign(d·N_j)` recovers the order.
/// Both factors are exact predicates, so the comparator is a true total order.
///
/// The pair `(P, Q)` is a parameter, not the seam pair: sub-unit 3d orders two seam
/// crossings along an *edge* of `f` by calling this with `(P, R)`, the edge's own
/// two planes. No new predicate is needed for that.
pub(crate) fn order_along(planes: &[PlaneInfo], p: usize, q: usize, i: usize, j: usize) -> i8 {
    t_orient3d(planes, p, q, i, j) * dir_sign(planes, p, q, j)
}

/// `+1` when a plane's stored normal already points out of its solid, `-1` when the
/// face is `Reversed` and the two oppose. `pub(crate)` so the toleranced `dir_sign`
/// wrapper ([`crate::tolerant`]) can bridge the frame3 `D` (over each face's outward
/// `tri`) to `plane_pair_dir_sign`'s stored-normal convention.
pub(crate) fn orient_sign(planes: &[PlaneInfo], i: usize) -> i8 {
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

/// The two planes an ordered ring's edge `i → i+1` shares beyond `P`, plus the two nodes'
/// third planes — everything [`order_along`] needs to say which way that edge runs.
pub(crate) fn ring_edge(
    p: usize,
    ring: &[[usize; 3]],
    i: usize,
) -> Result<(usize, usize, usize), BoolError> {
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
/// and `sign(det[…])` is [`plane_pair_dir_sign`](nacre_geom::intersect::plane_pair_dir_sign), already exact. It is never `0`: node `i`
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
    let det = t_plane_pair_dir_sign(planes, p, a, b);
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
    planes: &[PlaneInfo],
    canon: &[usize],
) -> Result<Vec<[usize; 3]>, BoolError> {
    loop_triples(&model.faces.get(f).outer, p, inc, planes, canon)
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
    planes: &[PlaneInfo],
    canon: &[usize],
) -> Result<Vec<Vec<[usize; 3]>>, BoolError> {
    model
        .faces
        .get(f)
        .inner
        .iter()
        .map(|l| loop_triples(l, p, inc, planes, canon))
        .collect()
}

/// A loop's vertices as three-plane triples.
///
/// The name normally comes from the loop itself — the face's own plane and the two neighbours the
/// meeting edges carry. **That fails when both neighbours lie on one plane**: an earlier boolean can
/// split a plane between two faces with opposite normals (a base's exposed top and the cantilever
/// underside above it), and a vertex where the loop runs straight through their shared line then
/// names one plane twice. Such a triple defines no point — `three_planes` answers `None`, and the
/// exact predicates, whose precondition is `D ≠ 0`, abort on it. Measured 2026-07-22: this is the
/// *only* path by which a degenerate triple reaches them.
///
/// So when the two neighbours are one class, the name is taken from **every plane touching the
/// vertex** ([`vertex_plane_indices`]) instead: exactly three classes ⇒ that is the name, and it is
/// the same set whichever face's loop asks, so welding stays consistent (a face whose loop does not
/// degenerate here derives the same three). More than three is a real four-plane concurrency and
/// fewer is a genuine straight angle — both decline. Three *dependent* planes share a line rather
/// than a point, so independence is checked too ([`t_plane_pair_dir_sign`], which reads the same
/// un-normalized coefficients the consumer does).
///
/// The common case is untouched, so no existing name moves.
fn loop_triples(
    l: &Loop,
    p: usize,
    inc: &EdgePlanes,
    planes: &[PlaneInfo],
    canon: &[usize],
) -> Result<Vec<[usize; 3]>, BoolError> {
    let hes = &l.half_edges;
    let edge = |he: &nacre_topo::HalfEdge| -> Result<([Handle<Vertex>; 2], [usize; 2]), BoolError> {
        inc.get(&he.edge)
            .copied()
            .ok_or_else(|| reject(tag::MISSING_SEAM))
    };
    let other = |pair: [usize; 2]| if pair[0] == p { pair[1] } else { pair[0] };
    let n = hes.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        // Vertex `i` starts edge `i` and ends edge `i - 1`.
        let (in_bounds, in_pair) = edge(&hes[(i + n - 1) % n])?;
        let (out_bounds, out_pair) = edge(&hes[i])?;
        let (a, b) = (other(in_pair), other(out_pair));
        // `inc` names faces, so `other` matches by face — but the triple names *planes*, and a
        // consumer's `==` on it must mean "same plane". Canonize here, once, at the source.
        let mut t = [canon[p], canon[a], canon[b]];
        t.sort_unstable();
        if t[0] != t[1] && t[1] != t[2] {
            out.push(t);
            continue;
        }
        // Both neighbours are one plane. The vertex is the two edges' shared endpoint — and only
        // if that is unambiguous (a two-gon or a self-bounded rim would give two, or none).
        let shared: Vec<Handle<Vertex>> = in_bounds
            .iter()
            .copied()
            .filter(|v| out_bounds.contains(v))
            .collect();
        let [vh] = shared[..] else {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH));
        };
        let mut classes: Vec<usize> = vertex_plane_indices(vh, inc)
            .into_iter()
            .map(|k| canon[k])
            .collect();
        classes.sort_unstable();
        classes.dedup();
        if classes.len() != 3 {
            // >3: a genuine four-plane concurrency, which this substrate cannot name.
            // <3: the vertex lies on fewer than three planes, so no triple names it — a genuine
            // straight angle, where dropping the vertex would preserve the polygon but this brick
            // does not drop vertices.
            //
            // A vertex can also lack a triple by sitting on a **coplanar seam** (an edge between
            // two coplanar, same-facing faces of one solid): its only edges are the seam's, so it
            // touches one plane class. That is *not* a straight angle — the vertex is a corner and
            // dropping it would change the shape — but no producer makes such an edge since
            // `ImprintSketch` was retired, so it cannot reach here. See the 2026-07-22 dev-log
            // cells if one ever does: the answer is that a whole-seam ring is not a boundary.
            return Err(reject(if classes.len() > 3 {
                tag::FOURPLANE
            } else {
                tag::LOOP_ORIENT_MISMATCH
            }));
        }
        // `canon[p]`, not `p`. An earlier revision kept `p` raw because consumers still matched the
        // face's own plane by raw index; they now compare classes (2026-07-22), and a triple that
        // mixed one face index with two class indices was exactly the ambiguity this brick exists
        // to remove.
        let mut t = [canon[p], 0, 0];
        let mut k = 1;
        for &c in &classes {
            if c != canon[p] {
                t[k] = c;
                k += 1;
            }
        }
        if t_plane_pair_dir_sign(planes, t[0], t[1], t[2]) == 0 {
            return Err(reject(tag::THREE_PLANES)); // three planes through one line, not one point
        }
        t.sort_unstable();
        out.push(t);
    }
    Ok(out)
}

/// The exact side of plane `q` that the implicit point `t` lies on: `0` means *on* it.
///
/// `+1` is the side [`PlaneInfo::tri`]'s right-hand normal points to — that is `n_out(q)`, the face's
/// **outward** side, since `outer_tri` winds the triangle outward.
///
/// ★ **That is not the frame the arrangement's labels are stated in.** A plane class's
/// `[*_above, *_below]` labels are about the class root's **stored surface normal** — the convention
/// `SegKind::Seated{body_above}` and `emit_faces`' `flip` are written against — and the two frames
/// differ by [`orient_sign`], which is `-1` exactly when the root face is `Reversed`. No
/// `add_cuboid` face ever is, but a face an earlier boolean re-emitted flipped is (a pocket wall),
/// so **a producer that turns raw `side_of` into an above/below *label* silently flips its bit on
/// such a class**; multiply by `orient_sign(q)` if that is what you are computing. Reading a sign
/// *difference* (does this edge cross `W`?) is frame-free and needs no correction.
pub(crate) fn side_of(planes: &[PlaneInfo], t: [usize; 3], q: usize) -> i8 {
    t_orient3d(planes, t[0], t[1], t[2], q)
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
/// `v` must not lie *on* `ring` — a hole ring never touches the outer ring it sits in, and a seam
/// loop never touches `∂f` — and this is where it is finally checked: an intersection at
/// `X == v` strictly inside an edge is the `POINT_ON_RING` reject.
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

/// The parity every clear ray reports. The ring is simple, so they must all agree; a golden
/// says so, which is a second machine for free.
pub(crate) fn every_ray(
    planes: &[PlaneInfo],
    p: usize,
    v: [usize; 3],
    ring: &[[usize; 3]],
) -> Result<Vec<bool>, BoolError> {
    let p = class_of(planes, p);
    // The vertex name is a plane triple, so it obeys the same rule as a ring's: class roots only.
    // A caller holding face indices (a hand-built table, a test) is normalized here rather than
    // silently comparing a face against a class.
    let mut v = [
        class_of(planes, v[0]),
        class_of(planes, v[1]),
        class_of(planes, v[2]),
    ];
    v.sort_unstable();
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
                if t_plane_pair_dir_sign(planes, p, qa, r) == 0 {
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

pub(crate) fn vertex_plane_indices(vh: Handle<Vertex>, inc: &EdgePlanes) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::new();
    for (bounds, pair) in inc.values() {
        if bounds.contains(&vh) {
            for &p in pair {
                if !out.contains(&p) {
                    out.push(p);
                }
            }
        }
    }
    out.sort_unstable();
    out
}

/// Classify vertex `vh` in/out of `other`, entirely on the exact index-plane substrate — the
/// 3D lift of [`point_in_ring`]. Non-convex containment is a winding-parity ray, but the ray
/// **line `L = A ∩ B` is built from two of the query vertex's own planes** (never an arbitrary
/// direction), so each crossing with a face on plane `q` is the three-plane point `{A,B,q}`,
/// judged inside/ahead by [`every_ray`]/[`order_along`] — no coordinate read, no f64.
///
/// `inc_v` is the query vertex's solid's [`EdgePlanes`] (to read `V`'s planes); `inc_o` is
/// `other`'s (its face rings). Both index into the shared `planes`/`surf_ix`. Faces of **every**
/// shell (outer + cavities) are summed, so a point in a void reads `Outside`.
///
/// A ray whose line grazes a node/edge, or leaves an in-face containment undecidable, is
/// abandoned for the vertex's next plane pair; if every pair is blocked the answer is honestly
/// `NO_CLEAR_RAY` — the degeneracy a later SoS layer resolves, and the same honest reject the
/// f64 ray gives via `RAY_DEGENERATE`. `V` lying on a face plane of `other` (a boundary-ish
/// query, its crossing at the ray origin) inside that face is one such abandon; off the face
/// (a disjoint-coplanar query) it is simply not a crossing.
// The per-op plane table (`planes`/`surf_ix`/`canon`) plus both edge-plane maps are all genuine
// inputs a classification needs; they travel together from one `plane_index_setup`.
/// An ordered ring's winding about the face's outward normal: `-1` clockwise — the material
/// is *outside* the ring, so it bounds a hole — and `+1` counter-clockwise, an island.
///
/// The turn at a convex-hull vertex is the winding, and the lexicographically smallest node
/// is one: it is an extreme point of the node set, which is planar, so it is a vertex of the
/// ring's hull. Finding it is the **only** thing here that needs two implicit points in one
/// decision, and [`three_plane_cmp_coord`](nacre_geom::intersect::three_plane_cmp_coord) is that predicate.
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
    let p = class_of(planes, p);
    if ring.len() < 3 {
        return Err(reject(tag::LOOP_ORIENT_MISMATCH));
    }
    let mut lo = 0usize;
    for i in 1..ring.len() {
        let ord = (0..3)
            .map(|axis| t_cmp_coord(planes, ring[i], ring[lo], axis))
            .find(|&c| c != 0);
        match ord {
            Some(-1) => lo = i,
            Some(_) => {}
            None => return Err(reject(tag::LOOP_ORIENT_MISMATCH)), // two nodes coincide
        }
    }
    turn_at(planes, p, ring, lo)
}

/// `sign((n_P × n_Q) · N_R)`, where `N_R` is the right-hand normal of `R.tri`.
///
/// [`plane_pair_dir_sign`](nacre_geom::intersect::plane_pair_dir_sign) gives the sign against `R`'s *stored* normal, exactly.
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
    t_plane_pair_dir_sign(planes, p, q, r) * orient_sign(planes, r)
}

/// The plane class a face index lies on, read off the table rather than a passed-around `canon`.
///
/// `PlaneInfo::class` is filled by `crate::fill_classes` for every table the boolean builds. A
/// hand-built table in a unit test may leave it unset, and then a face is its own class — which is
/// what the code did before classes existed, so such a table keeps its old behaviour instead of
/// indexing out of a `canon` it never had.
pub(crate) fn class_of(planes: &[PlaneInfo], k: usize) -> usize {
    if planes[k].class == usize::MAX {
        k
    } else {
        planes[k].class
    }
}
