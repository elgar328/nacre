//! Boolean assembly: turn the arrangement engine's per-face output (`LocalFace`, named by
//! `NodeId` seam triples) into result solids, and the mandatory coplanar-face cleaning pass.
//!
//! The public [`boolean`](fn@crate::boolean) entry runs the arrangement engine
//! ([`crate::arrangement`]), and the engine calls [`assemble_fuse_cut`], [`unify_coplanar_faces`]
//! and [`dissolve_straight_angles`] here to build and clean the shells; nothing here calls the
//! engine.

use crate::combinatorics::{self, NodeId, Wall};
use crate::draft::*;
use crate::planes::{ClassIx, WorkingPlane, uf_find};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, he_start, reject, unordered};
use nacre_math::Point3;
use nacre_store::Handle;
use nacre_topo::PointCache;
use nacre_topo::{Edge, Face, HalfEdge, Loop, Model, Orientation, Shell, Solid, Surface, Vertex};
use std::collections::{HashMap, HashSet};

mod coplanar;
mod face_bounds;
mod grouping;
mod naming;
mod reconstruct;
mod self_touch;
#[cfg(test)]
#[path = "../tests/probes/tess_census.rs"]
pub(crate) mod tess_census;
mod topology;

pub(crate) use coplanar::*;
pub(crate) use face_bounds::*;
pub(crate) use grouping::*;
pub(crate) use naming::*;
pub(crate) use reconstruct::*;
pub(crate) use self_touch::*;
pub(crate) use topology::*;

/// **Per result node, by name: the coordinate a refusal names as its witness and the reach its cache
/// proves**, which the sieves read — built once for the result naming: its split-twin subdivision
/// ([`name_result_vertices`]) and the result contract's two judgements ([`self_touch_reject`],
/// [`face_bounds_meet`]) read it through the methods below.
#[derive(Default)]
pub(crate) struct NodeReach(HashMap<NodeId, (Point3, crate::planes::CacheReach)>);

impl NodeReach {
    /// The box node `n`'s cache proves it in — the box that holds everything for a node the table
    /// does not hold, which keeps every candidate a sieve reads it for.
    pub(super) fn reach(&self, n: NodeId) -> crate::planes::CacheReach {
        self.0
            .get(&n)
            .map_or(([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]), |&(_, r)| r)
    }

    /// The box a straight edge `[u, v]` lies in: the hull of its ends' [`Self::reach`].
    pub(super) fn edge_reach(&self, [u, v]: [NodeId; 2]) -> crate::planes::CacheReach {
        let ((ul, uh), (vl, vh)) = (self.reach(u), self.reach(v));
        (
            core::array::from_fn(|k| ul[k].min(vl[k])),
            core::array::from_fn(|k| uh[k].max(vh[k])),
        )
    }

    /// The witness a refusal names at node `n`: its coordinate, where the table holds it.
    pub(super) fn witness_point(&self, n: NodeId) -> Option<crate::RejectWhere> {
        self.0.get(&n).map(|a| crate::RejectWhere::Point(a.0))
    }

    /// The witness a refusal names at the edge `[u, v]`: the edge itself, or the end the table
    /// holds.
    pub(super) fn witness_segment(&self, [u, v]: [NodeId; 2]) -> Option<crate::RejectWhere> {
        match (self.0.get(&u), self.0.get(&v)) {
            (Some(a), Some(b)) => Some(crate::RejectWhere::Segment([a.0, b.0])),
            (Some(a), None) | (None, Some(a)) => Some(crate::RejectWhere::Point(a.0)),
            (None, None) => None,
        }
    }
}

/// **A ring edge's identity in node space** — "a line is unordered, a circle is ordered" (`Wall`
/// is the type half, [`EdgeKey`] the handle-space half; this is the node-space half), and the one
/// spelling of it every stage that asks «which faces use this edge» reads: the component join
/// ([`face_components`]), the coplanar merge's join ([`unify_coplanar_faces`]) and the self-touch
/// verdict's owners ([`self_touch_reject`]).
///
/// A straight piece is its unordered endpoint pair — two points bound one straight segment,
/// whichever vocabulary (plane-pair line or ruling) a face states it in. An arc is its two carriers
/// and its ends **in CCW order about the axis**: a 2-node circle folds its chord and both
/// complementary arcs onto one node pair, so an unordered pair counts six users where each piece
/// really has two, and the CCW order is the bit that keeps the two arcs apart, exactly as it is in
/// the edge welding. The plane is in it too: a cylinder's cap cut by a plane oblique to it leaves,
/// between the two points they share, the cap's circular arc and the section's arc, both
/// counter-clockwise from one to the other — the lens of lateral below a tilted cylinder's cap that
/// a box's floor cuts. Without the plane the two are one key with four users. A key of the pair
/// alone folds a straight edge onto an arc between the same two nodes, and the self-touch verdict,
/// which asks straight edges only, skips both — 4,580 straight edges over the workspace suite.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum NodeEdge {
    Line((NodeId, NodeId)),
    Arc {
        cyl: usize,
        plane: usize,
        from: NodeId,
        to: NodeId,
    },
}

impl NodeEdge {
    /// The identity of `ring`'s edge `t` — from `ring.nodes[t]` to the next node, on
    /// `ring.walls[t]`.
    pub(super) fn of(ring: &Ring, t: usize) -> Self {
        let (a, b) = (ring.nodes[t], ring.nodes[(t + 1) % ring.nodes.len()]);
        match ring.walls[t] {
            Wall::Plane(_) | Wall::Ruling { .. } => NodeEdge::Line(norm_edge(a, b)),
            Wall::Arc { cyl, ccw, plane } => {
                let (from, to) = if ccw { (a, b) } else { (b, a) };
                NodeEdge::Arc {
                    cyl,
                    plane,
                    from,
                    to,
                }
            }
        }
    }
}

/// **Which faces make one output solid — decided before any topology exists.**
///
/// A result solid is a *material* component plus the cavity components nested in it, and that
/// grouping is the unit [`reconstruct`](fn@reconstruct) mints vertex and edge handles per: two
/// faces in one group may share a handle, two faces in different groups never do. Every question
/// asked here is asked of `LocalFace`/`NodeId`/`jd` alone — none of it reads the arena — which is
/// what lets it run first.
///
/// **The group, not the component, is the right unit.** A cavity that touches its host's outer
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

#[cfg(test)]
#[path = "../tests/assembly.rs"]
mod tests;
