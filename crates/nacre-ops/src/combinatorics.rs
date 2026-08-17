//! Combinatorial queries over the per-face planar arrangement (design §8 M5, sub-unit 3).
//!
//! This module answers combinatorial questions about how one solid's boundary cuts
//! a face of the other. It lives in `nacre-ops` and not in `nacre-geom` because it
//! needs `Model`/`Face`/`Edge`, and geom sits below topo (design §1: dependencies
//! flow upward only).
//!
//! Everything decided here is decided by an exact predicate. Coordinates that
//! appear (`three_planes`' cache) are never the basis of a decision — the truth of
//! a seam point is its plane triple, as it is for a measured vertex (design §4).
//!
//! # One `usize`, two meanings — now two tables
//!
//! An arrangement reasons **per plane**, but a solid gives you **faces**: an earlier boolean can
//! split one geometric plane between two faces (a base's exposed top and the cantilever underside
//! above it) whose outward normals **oppose**. Both were rows of one `planes` table, so the same
//! `usize` meant "face" here and "plane" there, and when the two meanings met in one comparison the
//! answer was wrong *silently* — four times on this branch, most recently as 117748 predicate calls
//! that read "different plane" for two faces of one plane and answered from rounding noise.
//!
//! That used to be held by a naming convention (`fp` / `fc`) and a debug-time net. It is now the
//! type: the predicates here take [`crate::WorkingPlane`], which has no face geometry to offer, and
//! `plane_ix` is the one place a face index becomes a plane index (in [`loop_triples`]).
//!
//! **Exception:** the code that *defines* the classes (`crate::fill_classes` →
//! `crate::shares_or_coplanar` → `Judge::planes_coplanar`) necessarily runs before a
//! plane table exists, so it takes face indices — hence that predicate's generic `Witness` bound.

use crate::planes::{ClassIx, WorkingPlane, edge_incidence};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, reject};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Loop, Model, Solid, Vertex};
use std::collections::HashMap;

/// A solid's edges, each with its endpoints and **the indices of its two faces**.
///
/// Named for what it holds: the pair is `planes`-table slots, i.e. *faces*, not plane classes —
/// `edge_faces` builds it from `surf_ix: HashMap<Handle<Face>, usize>`. It was `EdgePlanes`, and
/// that name is how a face index gets read as a plane one. The face→plane step is `plane_ix`, and
/// it happens in `loop_triples`, nowhere else.
pub(crate) type EdgeFaces = HashMap<Handle<Edge>, ([Handle<Vertex>; 2], [usize; 2])>;

/// Index a solid's edges by handle — [`edge_incidence`] keyed for lookup.
///
/// Hole-ring edges are in here too: `edge_incidence` walks every loop of every
/// face, and rejects an edge whose incidence is not a pair.
pub(crate) fn edge_faces(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<EdgeFaces, BoolError> {
    let mut out = EdgeFaces::new();
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
pub(crate) fn order_along(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    i: usize,
    j: usize,
) -> i8 {
    jd.orient3d(p, q, i, j) * dir_sign(jd, p, q, j)
}

/// One edge of a ring on plane `P`, carrying **its own geometry** rather than leaving it to be
/// recovered from the two endpoint names.
///
/// ★ **Why this type exists.** A vertex name is a plane triple, and for a long time the engine read
/// an edge's supporting plane back out of its endpoints — "the class the two names share besides
/// `P`". That works only while every vertex lies on exactly three planes. It is an accident of the
/// corpus, not an invariant: let four planes meet at a point, give the point one canonical name, and
/// the shared class is **some other plane than the one the edge rides**, silently. So the walker
/// that knows the edge — the DCEL half-edge, which was told its wall — hands the geometry over
/// instead, and only rings whose provenance is *names alone* go through [`ring_from_names`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RingEdge {
    /// Identity of the vertex this edge leaves.
    pub node: [usize; 3],
    /// The plane whose meet with `P` carries this edge.
    pub wall: usize,
    /// The endpoints as handles on `P ∩ wall` — the third plane pinning each there. Not a name:
    /// see [`RingEdge`]'s note.
    pub from_h: usize,
    pub to_h: usize,
}

/// Recover a ring's edges from its vertex names — the classic derivation, now in **one** place.
///
/// Sound exactly while each name lists all three of its planes and no more (see [`RingEdge`]).
/// ★ **Test-only since 2026-08-17**: the last production consumer (the tracer's crossed-edge
/// wall) reads the wall the producer carries (`NamedRing`) instead, so no production path
/// derives ring geometry from names any more. Kept for hand-built test rings, whose vertices
/// are clean three-plane points by construction.
#[cfg(test)]
pub(crate) fn ring_from_names(p: usize, ring: &[[usize; 3]]) -> Result<Vec<RingEdge>, BoolError> {
    (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let shared: Vec<usize> = a.iter().copied().filter(|x| b.contains(x)).collect();
            if shared.len() != 2 || !shared.contains(&p) {
                return Err(reject(RejectReason::RingNaming));
            }
            let wall = shared[usize::from(shared[0] == p)];
            let third = |t: [usize; 3]| t.iter().copied().find(|&x| x != p && x != wall);
            let (Some(from_h), Some(to_h)) = (third(a), third(b)) else {
                return Err(reject(RejectReason::RingNaming));
            };
            Ok(RingEdge {
                node: a,
                wall,
                from_h,
                to_h,
            })
        })
        .collect()
}

/// A ring's edges from its nodes and the **carried** wall of each edge.
///
/// ★ **This is the sound twin of [`ring_from_names`].** That one reads an edge's supporting plane
/// back out of its two endpoint names — "the class they share besides `P`" — which works only while
/// every vertex lies on exactly three planes. Let four meet at a point, give it one canonical name,
/// and the name need not mention the plane the edge rides at all; the two names can even share
/// nothing but `P`. Measured (2026-07, `docs/dev-log.md`): a boolean whose fourth plane fell on an
/// arrangement vertex produced exactly that, and `ring_from_names` declined `RingNaming`.
///
/// The **handle** is still derived, and soundly: it only has to be *some* plane through the node
/// that cuts the line `P ∩ wall`, which is what `order_along` asks of it. Under a concurrency there
/// are several and any will do, so the smallest is taken and the answer stays replay-stable. This is
/// the same rule `trace_transversal_face`'s `third_on_l` uses.
pub(crate) fn ring_edges_with_walls(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    nodes: &[[usize; 3]],
    walls: &[usize],
) -> Result<Vec<RingEdge>, BoolError> {
    if nodes.len() != walls.len() {
        return Err(reject(RejectReason::RingNaming));
    }
    let handle = |t: [usize; 3], wall: usize| -> Option<usize> {
        let mut cs: Vec<usize> = t
            .iter()
            .copied()
            .filter(|&c| c != p && c != wall && jd.plane_pair_dir_sign(p, wall, c) != 0)
            .collect();
        cs.sort_unstable();
        cs.first().copied()
    };
    (0..nodes.len())
        .map(|i| {
            let (a, b) = (nodes[i], nodes[(i + 1) % nodes.len()]);
            let wall = walls[i];
            let (Some(from_h), Some(to_h)) = (handle(a, wall), handle(b, wall)) else {
                return Err(reject(RejectReason::RingNaming));
            };
            Ok(RingEdge {
                node: a,
                wall,
                from_h,
                to_h,
            })
        })
        .collect()
}

/// `+1` when the ring's edge `i → i+1` runs along `d = n_P × n_Q`, `-1` against it.
fn edge_sign(jd: &Judge<'_, WorkingPlane>, p: usize, e: &RingEdge) -> Result<i8, BoolError> {
    // `order_along` is `sign((V_i − V_j)·d)`, so `-1` — `V_i` precedes `V_j` — is the edge
    // running along `+d`. Recomputed rather than remembered from the assembly walk.
    match order_along(jd, p, e.wall, e.from_h, e.to_h) {
        -1 => Ok(1),
        1 => Ok(-1),
        _ => Err(reject(RejectReason::CoincidentNodes)), // two nodes coincide
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
// Used by `loop_winding`'s tests and by the winding goldens in `lib.rs`; production reads the
// turn through `turn_between`, which lets the caller skip a straight stretch.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn turn_at(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    ring: &[RingEdge],
    i: usize,
) -> Result<i8, BoolError> {
    let n = ring.len();
    turn_between(jd, p, &ring[(i + n - 1) % n], &ring[i])
}

/// The turn from one edge to another, both on `P` — [`turn_at`] with the two edges named, so a
/// caller that had to look past a straight stretch can say which pair it means.
fn turn_between(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    arriving: &RingEdge,
    leaving: &RingEdge,
) -> Result<i8, BoolError> {
    let (a, b) = (arriving.wall, leaving.wall);
    let sa = edge_sign(jd, p, arriving)?;
    let sb = edge_sign(jd, p, leaving)?;
    let det = jd.plane_pair_dir_sign(p, a, b);
    if det == 0 {
        return Err(reject(RejectReason::StraightAngle));
    }
    Ok(sa * sb * det * jd.planes[p].frame_sign)
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
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
) -> Result<NamedRing, BoolError> {
    loop_triples(&model.faces.get(f).outer, p, inc, jd, plane_ix)
}

/// One loop in class form with each edge's **carried wall** beside it: `walls[i]` is the plane
/// class the edge `triples[i] → triples[i+1]` rides — the far face's class, read off `inc` where
/// the triples are produced ([`loop_triples`]), never re-derived from the endpoint names. The
/// same trust model as the merge's `Ring { nodes, walls }`: deriving a wall from two names is
/// sound only while every vertex lies on exactly three planes, and the carried value is total
/// even where the *names* degenerate (the fallback-named vertices still know their edges).
#[derive(Clone, Debug)]
pub(crate) struct NamedRing {
    pub triples: Vec<[usize; 3]>,
    pub walls: Vec<usize>,
}

/// Every loop of one face in class form — **the only thing the tracer needs from `Model`**.
///
/// A loop that could not be named is `None` rather than an error, because the two failures have
/// *different names at the call site* (`OuterRing` vs `HoleRing`) and which one applies is the
/// tracer's to say, not this table's.
#[derive(Clone, Debug, Default)]
pub(crate) struct FaceLoops {
    /// The outer loop, or `None` if [`face_vertex_triples`] declined.
    pub outer: Option<NamedRing>,
    /// One entry per hole ring, or `None` if [`hole_rings`] declined for **any** of them — a hole
    /// that cannot be named is not "no hole".
    pub holes: Option<Vec<NamedRing>>,
}

/// What one boolean's tracer reads instead of the `Model`: every face's loops, plus which slots
/// belong to which operand.
///
/// ★ **Both fields are independent of the plane being traced onto.** They were nevertheless
/// re-derived once per plane class — measured on an 80-fin fold at **2,280,285** calls (one per
/// face per class) for **13,492** distinct answers, 3.3% of the boolean. Hoisting them is what
/// makes the tracer a function of a face table rather than of a topology store, and the 169×
/// reduction comes along for free.
pub(crate) struct TraceInput {
    /// Each operand's faces: the `planes`-table slot and that face's loops, in the order the
    /// shells list them — the walk `trace_one` used to do over `Model`.
    ///
    /// **Compact on purpose.** This was a full-length `Vec<FaceLoops>` beside a list of slots — one
    /// row per table slot whether or not that slot's face was in the input. Pairing the slot with
    /// its loops makes the length the number of faces actually traced, which is what lets a caller
    /// hand the tracer a *subset* without the table's size leaking into the cost.
    pub faces: [Vec<(usize, FaceLoops)>; 2],
}

/// Derive [`TraceInput`] for one boolean, once.
///
/// Walked exactly as the tracer walked: shell by shell, `surf_ix` naming each face's slot. That is
/// what keeps the table's index space the `planes` one — a face missing from `surf_ix` cannot
/// happen, since `surf_ix` was built from the same two solids.
pub(crate) fn trace_input(
    model: &Model,
    operands: [(Handle<Solid>, &EdgeFaces); 2],
    surf_ix: &HashMap<Handle<Face>, usize>,
    n_faces: usize,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
) -> TraceInput {
    let _ = n_faces;
    let mut faces = [Vec::new(), Vec::new()];
    for (side, (solid, inc)) in operands.into_iter().enumerate() {
        for sh in crate::planes::solid_shell_handles(model, solid) {
            for &fh in &model.shells.get(sh).faces {
                let fp = surf_ix[&fh];
                faces[side].push((
                    fp,
                    FaceLoops {
                        outer: face_vertex_triples(model, fh, fp, inc, jd, plane_ix).ok(),
                        holes: hole_rings(model, fh, fp, inc, jd, plane_ix).ok(),
                    },
                ));
            }
        }
    }
    TraceInput { faces }
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
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
) -> Result<Vec<NamedRing>, BoolError> {
    model
        .faces
        .get(f)
        .inner
        .iter()
        .map(|l| loop_triples(l, p, inc, jd, plane_ix))
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
/// vertex** ([`vertex_face_indices`]) instead: exactly three classes ⇒ that is the name, and it is
/// the same set whichever face's loop asks, so welding stays consistent (a face whose loop does not
/// degenerate here derives the same three). More than three is a real four-plane concurrency and
/// fewer is a genuine straight angle — both decline. Three *dependent* planes share a line rather
/// than a point, so independence is checked too ([`Judge::plane_pair_dir_sign`], which reads the same
/// un-normalized coefficients the consumer does).
///
/// The common case is untouched, so no existing name moves.
fn loop_triples(
    l: &Loop,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
) -> Result<NamedRing, BoolError> {
    let hes = &l.half_edges;
    let edge = |he: &nacre_topo::HalfEdge| -> Result<([Handle<Vertex>; 2], [usize; 2]), BoolError> {
        inc.get(&he.edge)
            .copied()
            .ok_or_else(|| reject(RejectReason::MissingSeam))
    };
    let other = |pair: [usize; 2]| if pair[0] == p { pair[1] } else { pair[0] };
    let n = hes.len();
    let mut out = Vec::with_capacity(n);
    let mut walls = Vec::with_capacity(n);
    for i in 0..n {
        // Vertex `i` starts edge `i` and ends edge `i - 1`.
        let (in_bounds, in_pair) = edge(&hes[(i + n - 1) % n])?;
        let (out_bounds, out_pair) = edge(&hes[i])?;
        let (a, b) = (other(in_pair), other(out_pair));
        // Edge `i`'s carried wall: the far face's class, read off `inc` — total even where the
        // vertex *names* below have to fall back or decline (see [`NamedRing`]).
        walls.push(plane_ix[b].plane());
        // `inc` names faces, so `other` matches by face — but the triple names *planes*, and a
        // consumer's `==` on it must mean "same plane". Canonize here, once, at the source.
        let mut t = [
            plane_ix[p].plane(),
            plane_ix[a].plane(),
            plane_ix[b].plane(),
        ];
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
            return Err(reject(RejectReason::AmbiguousCorner));
        };
        let mut classes: Vec<usize> = vertex_face_indices(vh, inc)
            .into_iter()
            .map(|k| plane_ix[k].plane())
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
                RejectReason::FourPlane
            } else {
                // Fewer than three classes: there is no triple to name this vertex by.
                RejectReason::RingNaming
            }));
        }
        // `plane_ix[p].plane()`, not `p`. An earlier revision kept `p` raw because consumers still matched the
        // face's own plane by raw index; they now compare classes (2026-07-22), and a triple that
        // mixed one face index with two class indices was exactly the ambiguity this brick exists
        // to remove.
        let mut t = [plane_ix[p].plane(), 0, 0];
        let mut k = 1;
        for &c in &classes {
            if c != plane_ix[p].plane() {
                t[k] = c;
                k += 1;
            }
        }
        if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
            return Err(reject(RejectReason::ThreePlanes)); // three planes through one line, not one point
        }
        t.sort_unstable();
        out.push(t);
    }
    Ok(NamedRing {
        triples: out,
        walls,
    })
}

/// The exact side of plane `q` that the implicit point `t` lies on: `0` means *on* it.
///
/// `+1` is the side the witness triangle's right-hand normal points to — that is `n_out(q)`, the face's
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
pub(crate) fn side_of(jd: &Judge<'_, WorkingPlane>, t: [usize; 3], q: usize) -> i8 {
    jd.orient3d(t[0], t[1], t[2], q)
}

/// Where a ring meets the line that `q` cuts its plane along — see [`ring_against_plane`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Feature {
    /// Edge `edge` (node `edge` → node `edge + 1`) crosses the line strictly inside: both its
    /// endpoints are off `q` and on opposite sides of it.
    Crossing { edge: usize },
    /// `len` consecutive nodes from `first` lie *on* `q`. `len >= 2` means the edges between them
    /// lie on the line too — an on-line interval rather than a point.
    ///
    /// `flanks_differ` is the whole decision: the two off-line neighbours bracketing the run sit on
    /// **opposite** sides, so the ring genuinely crosses the line here; equal sides mean it touched
    /// and turned back, and nothing crossed.
    Run {
        first: usize,
        len: usize,
        flanks_differ: bool,
    },
}

/// Read a ring against one plane: where it meets the line, and whether it crosses or only touches.
///
/// ★ **One walk, two consumers.** [`trace_transversal_face`] clips a face's ring against a cut
/// plane and [`every_ray`] casts a parity ray along `P ∩ Q_a`; both must answer the same question
/// first — *does the boundary cross this line here?* — and a node sitting **on** the line is the
/// only hard part of it. The tracer had the rule (look at the node's two off-line neighbours:
/// opposite sides is a crossing, equal sides a touch) inlined in its scan, entangled with naming,
/// alias recording and decline kinds; the ray caster had no rule at all and threw such a candidate
/// away. Resilience that lives in one consumer is resilience the other does not have — the same
/// shape [`ring_in_ring`]'s retry was in.
///
/// What each consumer does *with* a feature stays its own: the tracer turns it into a **named
/// point** (four-plane aliases, `DeclineKind`, occupancy), the ray caster into one bit ("ahead of
/// `v`?"). That is where the sharing stops.
///
/// `None` when every node lies on `q` — a ring in the plane has no flanks to be decided by.
///
/// **Features come out in ring order from the first off-`q` node.** That is the order the tracer's
/// scan produced them in, and its naming step records aliases into a union-find as it goes, so the
/// order is contract, not incident.
pub(crate) fn ring_against_plane(
    jd: &Judge<'_, WorkingPlane>,
    nodes: &[[usize; 3]],
    q: usize,
) -> Option<Vec<Feature>> {
    let n = nodes.len();
    let side: Vec<i8> = (0..n).map(|i| side_of(jd, nodes[i], q)).collect();
    let start = side.iter().position(|&s| s != 0)?;
    let mut out = Vec::new();
    let mut j = 0;
    while j < n {
        let i = (start + j) % n;
        if side[i] != 0 {
            let ni = (i + 1) % n;
            if side[ni] != 0 && side[ni] != side[i] {
                out.push(Feature::Crossing { edge: i });
            }
            j += 1;
        } else {
            // A maximal run of on-line vertices. Two is the common case, but a vertex whose name
            // had to be taken from its touching planes (`loop_triples`) stays in the ring even
            // when the loop runs straight through it, so a run can be longer.
            let first = i;
            let mut len = 0;
            while j < n && side[(start + j) % n] == 0 {
                len += 1;
                j += 1;
            }
            let before = side[(first + n - 1) % n];
            let after = side[(start + j) % n];
            out.push(Feature::Run {
                first,
                len,
                flanks_differ: before != after,
            });
        }
    }
    Some(out)
}

/// Is the implicit point `v` inside the simple ring `ring`, both on face plane `p`?
///
/// **A ray, cast along a line we already have.** Every ring edge lies on `P ∩ R`, and `v`
/// lies on `P ∩ Q_a` for either of its own two planes. Those two lines meet at
/// `X = {P, Q_a, R}`, which is *itself* a three-plane point — so "is `X` inside the edge"
/// and "is `X` ahead of `v`" are both [`order_along`], the comparator cell 3d already built
/// for two three-plane points on one line. **No coordinate is read and no point is built.**
///
/// **The flanks delete the special case, not the choice of ray.** A ring node *on* the ray's line
/// leaves no room for "is `X` inside the edge" to decide anything, and this used to abandon the
/// candidate; with both of `v`'s candidates abandoned the question came back `no_clear_ray`, and a
/// band of rotation angles died of it. But the node is not ambiguous at all — its two off-line
/// **neighbours** settle it: opposite sides and the boundary crossed here, equal sides and it
/// touched and turned back. That is the rule `trace_transversal_face` has always read a ring with,
/// and [`ring_against_plane`] is now where both get it.
///
/// So a `Feature::Run` — one node, or a whole edge of the ring lying on the line — contributes one
/// crossing iff its flanks differ, and a `Feature::Crossing` contributes one where it always did.
/// Nothing is counted twice: a crossing's endpoints are both off the line by construction.
///
/// Candidates are each node's two non-`P` planes, in ring order, `+d` before `-d`; the first
/// usable one wins, which keeps the answer deterministic. `no_clear_ray` survives for the two
/// cases nothing can name: a ring lying wholly on `Q_a` (no flanks), and an on-line node whose own
/// two walls are both parallel to the line (nothing pins it there). The answer must not depend on
/// which candidate was chosen, and a golden says so.
///
/// `v` must not lie *on* `ring` — a hole ring never touches the outer ring it sits in, and a seam
/// loop never touches `∂f` — and this is where it is finally checked: an intersection at
/// `X == v` strictly inside an edge is the `POINT_ON_RING` reject.
pub(crate) fn point_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    v: [usize; 3],
    ring: &[RingEdge],
) -> Result<bool, BoolError> {
    every_ray(jd, p, v, ring)?
        .first()
        .copied()
        .ok_or_else(|| reject(RejectReason::NoClearRay))
}

/// Is the ring whose nodes are `probes` inside `outer` — **asked of the ring, not of one point**.
///
/// [`point_in_ring`] answers for a single vertex and can honestly fail there: every ray it can
/// cast from that vertex may have a ring node on its line. That is a fact about the *probe*, not
/// about the two rings, so the question is retried from the next node and only an exhausted ring
/// is a real degeneracy.
///
/// ★ **The retry belongs here, not in the callers.** It used to live in two of them, spelled two
/// different ways, and the third — the coplanar-merge cleaning pass — never got it: it probed
/// `nodes[0]` alone and rejected the whole boolean when that one vertex happened to be spoiled.
/// Measured, that lost a *band* of rotation angles at a stroke, because a small change of angle
/// leaves the topology (and therefore the unlucky first node) exactly where it was. Resilience
/// that lives in a caller is resilience the next caller does not have.
///
/// The probes are tried in order, so the answer is deterministic. Callers decide what *adjacency*
/// means for them (see [`arrangement::nest_cells`] and `innermost_host`, which disagree) and ask
/// this only about rings they have already established are disjoint.
pub(crate) fn ring_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    probes: &[[usize; 3]],
    outer: &[RingEdge],
) -> Result<bool, BoolError> {
    for &v in probes {
        if let Ok(hit) = point_in_ring(jd, p, v, outer) {
            return Ok(hit);
        }
    }
    Err(reject(RejectReason::NoClearRay))
}

/// The parity every clear ray reports. The ring is simple, so they must all agree; a golden
/// says so, which is a second machine for free.
pub(crate) fn every_ray(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    v: [usize; 3],
    ring: &[RingEdge],
) -> Result<Vec<bool>, BoolError> {
    // The vertex name is a plane triple, so it obeys the same rule as a ring's: class roots only.
    // A caller holding face indices (a hand-built table, a test) is normalized here rather than
    // silently comparing a face against a class.
    let mut v = [v[0], v[1], v[2]];
    v.sort_unstable();
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    let nodes: Vec<[usize; 3]> = ring.iter().map(|e| e.node).collect();
    let mut out = Vec::new();
    for &qa in v.iter().filter(|&&x| x != p) {
        // Where the ring meets the line — the walk `trace_transversal_face` reads too.
        //
        // ★ **A `Crossing` is exactly what this loop used to derive per edge.** It said "is
        // `X = {P, Q_a, R}` strictly inside the edge" with two `order_along`s: `a` is the sign of
        // `X − From` along `P ∩ R` and `b` that of `X − To`, both normalized to the same direction
        // on the same line, so `a·b < 0` iff `From` and `To` lie on opposite sides of `Q_a` — which
        // is what the walk already knows from their sides. The parallel guard goes with it: an edge
        // whose line is parallel to `P ∩ Q_a` has both endpoints on one side and is not a crossing.
        let Some(features) = ring_against_plane(jd, &nodes, qa) else {
            continue; // the whole ring lies on `Q_a`
        };
        let qb = *v
            .iter()
            .find(|&&x| x != p && x != qa)
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        // A node on the line is a crossing point in its own right, so it must be nameable as one:
        // some plane of its own, off the line, pins it there. (`third_on_l` picks a handle the
        // same way, and for the same reason — one parallel to the line names no point on it.)
        let namer = |i: usize| {
            nodes[i]
                .iter()
                .copied()
                .find(|&x| x != p && x != qa && jd.plane_pair_dir_sign(p, qa, x) != 0)
        };
        // Which side of `v` each run sits on, and whether the ring crossed the line there at all.
        // A run is one interval of the line and the ring is simple, so it cannot double back
        // inside itself: its two ends bracket it, and ends that disagree mean `v` is *between*
        // them — on the ring.
        let mut run_hits: Vec<i8> = Vec::new();
        let mut unnameable = false;
        for f in &features {
            let Feature::Run {
                first,
                len,
                flanks_differ,
            } = *f
            else {
                continue;
            };
            let ends = [first, (first + len - 1) % nodes.len()];
            let (Some(lo), Some(hi)) = (namer(ends[0]), namer(ends[1])) else {
                unnameable = true;
                break;
            };
            let o = [
                order_along(jd, p, qa, lo, qb),
                order_along(jd, p, qa, hi, qb),
            ];
            if o[0] != o[1] || o[0] == 0 {
                return Err(reject(RejectReason::PointOnRing)); // `v` inside the run, or one of it
            }
            run_hits.push(if flanks_differ { o[0] } else { 0 });
        }
        if unnameable {
            continue;
        }
        for dir in [1i8, -1] {
            let mut crossings = run_hits.iter().filter(|&&o| o == dir).count();
            for f in &features {
                let Feature::Crossing { edge } = *f else {
                    continue; // runs are counted above
                };
                // Strictly ahead of `v` along `dir · (n_P × n_Qa)`?
                match order_along(jd, p, qa, ring[edge].wall, qb) {
                    0 => return Err(reject(RejectReason::PointOnRing)), // `X == v`, inside an edge
                    o if o == dir => crossings += 1,
                    _ => {}
                }
            }
            out.push(crossings % 2 == 1);
        }
    }
    Ok(out)
}

/// **Does the open segment between two named points on one line meet this face's material?**
///
/// [`every_ray`]'s sibling, and not its copy. That one asks about a **point** and may bail with
/// `PointOnRing` when the point lands on the boundary — here the two endpoints are *expected* to,
/// because the defect this answers is an edge whose ends sit on a face's ring while its middle
/// crosses the interior. An interval query has no candidate to fall back to, so every case that one
/// declines has to become a value.
///
/// ★ **One algorithm, no special cases.** The rings meet the line `P ∩ w` at a set of places; the
/// two unbounded ends of the line are outside the face and each genuine crossing flips that, so the
/// line reads **outside / inside / outside / …**. The answer is whether any *inside* stretch
/// overlaps the open `(u, v)`. Written as branches — "is there a crossing between them", "is it in a
/// hole" — it was twice wrong, because each branch re-derived a piece of that structure and lost
/// another. Read as one alternation it also subsumes the endpoint test: an endpoint strictly inside
/// puts `u` in an inside stretch.
///
/// **Holes come along for free.** All rings go into one bag: the rings of a face are disjoint, so
/// even-odd over the union *is* the material region (a point inside a hole has crossed twice) — the
/// rule `design.md` states one dimension down for 2D sketches.
///
/// ★★ **A `Run` is a stretch, not a place.** `Run { len >= 2 }` means the boundary *lies along* the
/// line, so it occupies an interval where the segment would be **on** the face rather than inside
/// it — two faces sharing an edge, which is ordinary adjacency. Those stretches are boundary and
/// are not counted as inside; `flanks_differ` still says whether crossing the run flips the side,
/// which is the rule the tracer and the ray caster already read a ring with.
pub(crate) fn segment_meets_face(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    w: usize,
    u: [usize; 3],
    v: [usize; 3],
    rings: &[Vec<RingEdge>],
) -> Result<bool, BoolError> {
    // A point on `P ∩ w` is named there by a third plane of its own that **cuts** the line; one
    // parallel to it names nothing (the same duty `every_ray`'s `namer` states).
    let handle = |t: [usize; 3]| -> Option<usize> {
        t.iter()
            .copied()
            .find(|&x| x != p && x != w && jd.plane_pair_dir_sign(p, w, x) != 0)
    };
    let (Some(hu), Some(hv)) = (handle(u), handle(v)) else {
        return Err(reject(RejectReason::RingNaming));
    };
    // Every place a ring meets the line: `[lo, hi]` handles (equal for a crossing at a point) and
    // whether passing it flips inside/outside.
    let mut events: Vec<([usize; 2], bool)> = Vec::new();
    for ring in rings {
        let nodes: Vec<[usize; 3]> = ring.iter().map(|e| e.node).collect();
        let Some(features) = ring_against_plane(jd, &nodes, w) else {
            // The whole ring lies on `w`: this face's boundary is the line itself, and the
            // alternation has no crossings to read. Refusing to guess.
            return Err(reject(RejectReason::PointOnRing));
        };
        for f in &features {
            match *f {
                Feature::Crossing { edge } => {
                    let h = ring[edge].wall;
                    events.push(([h, h], true));
                }
                Feature::Run {
                    first,
                    len,
                    flanks_differ,
                } => {
                    let ends = [first, (first + len - 1) % nodes.len()];
                    let (Some(a), Some(b)) = (handle(nodes[ends[0]]), handle(nodes[ends[1]]))
                    else {
                        return Err(reject(RejectReason::RingNaming));
                    };
                    let lo_first = order_along(jd, p, w, a, b) <= 0;
                    events.push((if lo_first { [a, b] } else { [b, a] }, flanks_differ));
                }
            }
        }
    }
    events.sort_by(|x, y| match order_along(jd, p, w, x.0[0], y.0[0]) {
        -1 => std::cmp::Ordering::Less,
        1 => std::cmp::Ordering::Greater,
        _ => std::cmp::Ordering::Equal,
    });
    // Walk the line: outside before the first event, flipping as each genuine crossing is passed.
    // The stretch between two events is a cell; an inside cell that overlaps the open `(u, v)` is
    // the surface meeting itself.
    let (lo, hi) = if order_along(jd, p, w, hu, hv) <= 0 {
        (hu, hv)
    } else {
        (hv, hu)
    };
    let mut inside = false;
    for i in 0..events.len() {
        if events[i].1 {
            inside = !inside;
        }
        if !inside {
            continue;
        }
        // The cell runs from this event's far end to the next event's near end.
        let cell_start = events[i].0[1];
        let Some(next) = events.get(i + 1) else {
            // ★ Reaching the unbounded tail while *inside* means the rings crossed the line an odd
            // number of times, which a closed curve cannot do. Reading it as "outside" would let a
            // real contact past in silence, so it is named instead — the same rule the rest of this
            // engine follows for an invariant it cannot verify.
            return Err(reject(RejectReason::RingParity));
        };
        let cell_end = next.0[0];
        // Overlap with the **open** interval: strictly, so touching at `u` or `v` is not inside.
        if order_along(jd, p, w, cell_start, hi) < 0 && order_along(jd, p, w, cell_end, lo) > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Every **face** incident to `vh`, as `planes`-table slots. (Was `vertex_plane_indices`; it
/// returns `inc`'s pairs verbatim, and those are faces. Its one caller maps them through
/// `plane_ix`.)
pub(crate) fn vertex_face_indices(vh: Handle<Vertex>, inc: &EdgeFaces) -> Vec<usize> {
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

/// Whether the point named by plane triple `v` lies on any edge of `ring` (a ring on plane `p`).
/// Used by [`point_in_component`] to abandon a non-generic ray rather than guess on a boundary.
fn point_on_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    mut v: [usize; 3],
    ring: &[RingEdge],
) -> Result<bool, BoolError> {
    v.sort_unstable();
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    for e in ring {
        let (r, si, sj) = (e.wall, e.from_h, e.to_h);
        if side_of(jd, v, r) != 0 {
            continue; // `v` is not even on the edge's line
        }
        // Name `v` as a point of that line: `{p, r, s}` for one of its own planes `s` off the line.
        let s = *v
            .iter()
            .find(|&&x| x != p && jd.plane_pair_dir_sign(p, r, x) != 0)
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        let (a, b) = (order_along(jd, p, r, s, si), order_along(jd, p, r, s, sj));
        if a * b <= 0 {
            return Ok(true); // between the endpoints (or on one)
        }
    }
    Ok(false)
}

/// Whether the point named by plane triple `query` lies inside a connected component — the exact
/// 3D lift of [`point_in_ring`]. A winding-parity ray whose line `L = a ∩ b` is built from two of
/// the query's own planes (never an arbitrary direction): each crossing with a face on plane `q`
/// is the exact three-plane point `{a,b,q}`, so the whole test is on the index-plane substrate —
/// no coordinate read, no f64. Non-convex is native (parity, not a convex test).
///
/// `faces` is the component as `(plane, rings)` per face — `rings[0]` = outer, `rings[1..]` = holes.
/// Only this component's faces are summed, so testing a void's vertex against a material component
/// reads `true` iff that material's outer shell nests the void.
///
/// A ray grazing a face boundary, or an undecidable in-face containment, is abandoned for the
/// query's next plane pair. **`Ok(None)` means every pair of this query's planes was blocked** —
/// this *node* cannot decide the question, which is a fact about the node, not an error: the
/// callers hold other nodes to try, and "every node abstained" is *their* proposition to reject
/// (`NoClearRay`, raised where the retries actually run out). `Err` is reserved for the
/// judgement itself failing (`JudgeExhausted`, a ring that cannot be named, …) — those must
/// propagate, never be traded for the next node: an abstention has other nodes as its remedy,
/// a failed judgement does not, and retrying it would let a real cause masquerade as
/// "no clear ray" once every node hit it.
/// A component as `(plane, rings)` per face, each ring already carrying its edges' walls.
pub(crate) type ComponentFaces = Vec<(usize, Vec<Vec<RingEdge>>)>;

pub(crate) fn point_in_component(
    jd: &Judge<'_, WorkingPlane>,
    query: [usize; 3],
    faces: &[(usize, Vec<Vec<RingEdge>>)],
) -> Result<Option<bool>, BoolError> {
    let mut vplanes = query.to_vec();
    vplanes.sort_unstable();
    vplanes.dedup();

    // One ray attempt along `L = a ∩ b`, located by the query point `{a,b,c}`. `Ok(None)` = the
    // line grazed and the caller should try the next plane pair.
    let attempt = |a: usize, b: usize, c: usize| -> Result<Option<bool>, BoolError> {
        let mut count = 0usize;
        for (q, rings) in faces {
            let q = *q;
            if jd.plane_pair_dir_sign(a, b, q) == 0 {
                continue; // `L` parallel to (or in) plane `q` — no transversal crossing
            }
            let mut x = [a, b, q];
            x.sort_unstable();
            for ring in rings {
                if point_on_ring(jd, q, x, ring)? {
                    return Ok(None); // `x` on `q`'s boundary — non-generic, abandon
                }
            }
            // Inside `q`'s material: inside the outer ring, outside every hole.
            let inside_outer = match every_ray(jd, q, x, &rings[0])?.first().copied() {
                Some(v) => v,
                None => return Ok(None),
            };
            let mut in_g = inside_outer;
            if in_g {
                for hole in &rings[1..] {
                    match every_ray(jd, q, x, hole)?.first().copied() {
                        Some(true) => {
                            in_g = false;
                            break;
                        }
                        Some(false) => {}
                        None => return Ok(None),
                    }
                }
            }
            let fwd = order_along(jd, a, b, c, q);
            if fwd == 0 {
                // The crossing is the ray origin itself (query on plane `q`). Inside `q`'s
                // material ⇒ query on the component's surface ⇒ undecidable → abandon.
                if in_g {
                    return Ok(None);
                }
                continue;
            }
            if fwd == 1 && in_g {
                count += 1; // one side of the line — parity is the same on either half
            }
        }
        Ok(Some(count % 2 == 1))
    };

    for i in 0..vplanes.len() {
        for j in (i + 1)..vplanes.len() {
            let (a, b) = (vplanes[i], vplanes[j]);
            // Locator `c`: a plane of the query off the line `a ∩ b`, so `{a,b,c}` is the query.
            let Some(&c) = vplanes
                .iter()
                .find(|&&x| x != a && x != b && jd.plane_pair_dir_sign(a, b, x) != 0)
            else {
                continue;
            };
            if let Some(inside) = attempt(a, b, c)? {
                return Ok(Some(inside));
            }
        }
    }
    // Every plane pair of this query was blocked: the node abstains. The inner `attempt`
    // already speaks this language per pair (`Ok(None)`); the boundary now keeps it instead of
    // dressing the abstention up as an error for the caller to catch and swallow.
    Ok(None)
}
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
/// A ring may be *non-simple* — visiting one node twice — and still be a legitimate face: the
/// unbounded contour of two cells that meet at a single point pinches through that point, tracing
/// a figure-8. The winding is read from the turn at the lexicographically smallest node, and a
/// coincidence elsewhere in the ring does not affect that turn, so a repeated node is not by itself
/// an error. (The old check rejected on the first coincidence with the running minimum, which made
/// the verdict depend on the ring's arbitrary start index — one operand order rejected a pinch the
/// other accepted.) Only a pinch *at* the extreme node itself leaves the turn ambiguous; that stays
/// a `LOOP_ORIENT_MISMATCH`, decided by exact equality rather than by a tolerance.
pub(crate) fn loop_winding(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    ring: &[RingEdge],
) -> Result<i8, BoolError> {
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    // Lexicographically smallest node — a hull vertex, hence a valid turn site. A coincidence with
    // the running minimum just means "not strictly smaller", so keep it; do not reject.
    let mut lo = 0usize;
    for i in 1..ring.len() {
        let strictly_less = (0..3)
            .map(|axis| jd.cmp_coord(ring[i].node, ring[lo].node, axis))
            .find(|&c| c != 0)
            == Some(-1);
        if strictly_less {
            lo = i;
        }
    }
    // ★★★ **The scan above is only a minimum if the relation is an order, and that is an
    // assumption about the predicates, not about this loop.** It composes per-axis comparisons
    // lexicographically; a per-axis answer that is not a fact about the geometry — one plane
    // described two ways, say — makes the composition intransitive, and then a forward scan can
    // stop at a node with something smaller behind it. That node is not extreme, the turn read
    // there is not the winding, and the wrong sign comes back **confident**: the engine noticed
    // only two layers later, as "no outer contour", and named the symptom.
    //
    // ★ This is the postcondition the algorithm actually needs — cheaper than asking whether the
    // relation is transitive (`O(n)` against `O(n³)`, and rings here reach 95 nodes) and closer to
    // the point. **It holds however the predicates behave; it is the net under them.**
    debug_assert!(
        !ring.iter().enumerate().any(|(i, _)| {
            i != lo
                && (0..3)
                    .map(|axis| jd.cmp_coord(ring[i].node, ring[lo].node, axis))
                    .find(|&c| c != 0)
                    == Some(-1)
        }),
        "the lexicographic scan did not find a minimum — the comparison is not an order here"
    );
    // The turn is read at `lo`; if that exact point recurs the corner is a pinch and its turn is
    // ambiguous — honest-reject rather than guess.
    let pinched_extreme = ring.iter().enumerate().any(|(i, _)| {
        i != lo && (0..3).all(|axis| jd.cmp_coord(ring[i].node, ring[lo].node, axis) == 0)
    });
    if pinched_extreme {
        return Err(reject(RejectReason::CoincidentNodes));
    }
    // **A ring node need not be a corner.** The arrangement names a point wherever another feature
    // crosses an edge, and `loop_triples` keeps such a vertex even when the loop runs straight
    // through it — so one edge of the polygon can arrive as several collinear ring edges. Reading
    // the turn at `lo` against its immediate predecessor then asks about two halves of one
    // straight edge, which has no turn to give.
    //
    // The turn to read is the one between the directions the loop **actually** arrives and leaves
    // on: walk back past the edges collinear with the leaving one. `lo` stays a hull vertex — the
    // stretch lies on a line through it, so the polygon is still on one side of that line.
    //
    // **Only while the stretch keeps going the same way.** A collinear edge traversed the *other*
    // way means the ring doubles back along the line it came in on — an antenna, whose tip has no
    // turn and whose neighbours' turn belongs to a different vertex. Skipping past that would
    // read a turn from somewhere else and call it this vertex's: a wrong winding, silently.
    let n = ring.len();
    let dir = edge_sign(jd, p, &ring[lo])?;
    let mut back = (lo + n - 1) % n;
    while jd.plane_pair_dir_sign(p, ring[back].wall, ring[lo].wall) == 0 {
        if edge_sign(jd, p, &ring[back])? != dir {
            return Err(reject(RejectReason::StraightAngle)); // the ring doubles back here
        }
        back = (back + n - 1) % n;
        if back == lo {
            // Every edge of the ring lies on one line: it bounds nothing.
            return Err(reject(RejectReason::DegenerateRing));
        }
    }
    turn_between(jd, p, &ring[back], &ring[lo])
}

/// `sign((n_P × n_Q) · N_R)`, where `N_R` is the right-hand normal of `R.tri`.
///
/// [`plane_pair_dir_sign`](nacre_geom::intersect::plane_pair_dir_sign) gives the sign against `R`'s *stored* normal, exactly.
/// That normal is parallel to `N_R` but may oppose it on a `Reversed` face, so we
/// correct with their dot — two parallel unit vectors, `|·| ≈ 1`, nowhere near the
/// sign boundary.
///
/// The correction *is* the face's stated flag — since the stored-orientation
/// cutover, `frame_sign` is `Forward`/`Reversed` as a sign, and
/// "`Reversed` ⇔ `n_out = −plane.normal()`" holds by construction rather than by
/// hope. What keeps it honest is the winding: `collect_planes` debug_asserts the
/// witness triangle against `n_out`, and `validate` pins the loop itself as
/// `FaceMisoriented`.
pub(crate) fn dir_sign(jd: &Judge<'_, WorkingPlane>, p: usize, q: usize, r: usize) -> i8 {
    let planes = jd.planes;
    jd.plane_pair_dir_sign(p, q, r) * planes[r].frame_sign
}
