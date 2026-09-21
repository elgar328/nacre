//! Boolean assembly: turn the arrangement engine's per-face output (`LocalFace`, named by
//! `NodeId` seam triples) into result solids, and the mandatory coplanar-face cleaning pass.
//!
//! The public [`boolean`] entry lives here and delegates the cell-complex work to
//! [`crate::arrangement`]; the engine calls back into [`assemble_fuse_cut`] and
//! [`unify_coplanar_faces`] to build and clean the shells (a legal module cycle).

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
mod cycles;
mod grouping;
mod naming;
#[cfg(test)]
#[path = "../tests/probes/assembly_probe.rs"]
pub(crate) mod probe;
mod reconstruct;
mod self_touch;
#[cfg(test)]
#[path = "../tests/probes/tess_census.rs"]
pub(crate) mod tess_census;
mod topology;

pub(crate) use coplanar::*;
pub(crate) use cycles::*;
pub(crate) use grouping::*;
pub(crate) use naming::*;
pub(crate) use reconstruct::*;
pub(crate) use self_touch::*;
pub(crate) use topology::*;

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

#[cfg(test)]
#[path = "../tests/assembly.rs"]
mod tests;
