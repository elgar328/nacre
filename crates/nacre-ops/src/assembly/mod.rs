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
mod grouping;
mod naming;
mod reconstruct;
mod self_touch;
#[cfg(test)]
#[path = "../tests/probes/tess_census.rs"]
pub(crate) mod tess_census;
mod topology;

pub(crate) use coplanar::*;
pub(crate) use grouping::*;
pub(crate) use naming::*;
pub(crate) use reconstruct::*;
pub(crate) use self_touch::*;
pub(crate) use topology::*;

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
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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
