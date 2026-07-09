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
        let mut third: Vec<usize> = Vec::new();
        for (face, rings, inc, own) in [(f, &g_rings, inc_x, p), (g, &f_rings, inc_y, q)] {
            for he in &model.faces.get(face).outer.half_edges {
                let (bounds, [pa, pb]) = inc[&he.edge];
                let r = if pa == own { pb } else { pa };
                let (v0, v1) = (
                    model.vertices.get(bounds[0]).point,
                    model.vertices.get(bounds[1]).point,
                );
                if segment_crosses_face(v0, v1, rings)? {
                    third.push(r);
                }
            }
        }
        if third.is_empty() {
            continue;
        }

        // Order along `L`. With `V_i = P ∩ Q ∩ R_i`, we want `sign((V_i − V_j)·d)`.
        // Since `V_j ∈ R_j`, `three_plane_orient3d(P, Q, R_i, R_j.tri)` is
        // `sign((V_i − V_j)·N_j)` for `N_j` the RH normal of `R_j.tri`; multiplying by
        // `sign(d·N_j)` recovers the order. Both factors are exact predicates.
        let order = |i: usize, j: usize| -> i8 {
            three_plane_orient3d(
                &planes[p].plane,
                &planes[q].plane,
                &planes[i].plane,
                planes[j].tri[0],
                planes[j].tri[1],
                planes[j].tri[2],
            ) * dir_sign(planes, p, q, j)
        };
        third.sort_by(|&i, &j| match order(i, j) {
            -1 => Ordering::Less,
            1 => Ordering::Greater,
            _ => Ordering::Equal,
        });
        for w in third.windows(2) {
            if order(w[0], w[1]) != -1 {
                // Coincident crossings, or `L` through a vertex of a face: four planes
                // meet at one point.
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
            let (r0, r1) = (w[0], w[1]);
            let point_of = |r: usize| {
                three_planes(&planes[p].plane, &planes[q].plane, &planes[r].plane)
                    .ok_or_else(|| reject(tag::THREE_PLANES))
            };
            let seg = SeamSegment {
                plane_pair: [p, q],
                ends: [triple(p, q, r0), triple(p, q, r1)],
                points: [point_of(r0)?, point_of(r1)?],
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
    let dot = planes[r].plane.normal().dot(planes[r].n_out);
    debug_assert!(
        dot.abs() > 0.5,
        "a plane's normal must be parallel to n_out"
    );
    debug_assert_eq!(
        dot > 0.0,
        planes[r].orient == Orientation::Forward,
        "n_out's sign against the surface normal is the face's orientation"
    );
    let s = plane_pair_dir_sign(&planes[p].plane, &planes[q].plane, &planes[r].plane);
    if dot > 0.0 { s } else { -s }
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
