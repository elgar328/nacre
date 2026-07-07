//! `Adjacency` — the reverse index over the topology stores.

use crate::topology::{Edge, Face, Vertex};
use nacre_store::Handle;
use std::collections::HashMap;

/// Reverse index over the topology stores — a **cache, not truth** (design §4,
/// overview 절대원칙 5). Rebuilt from scratch by scanning the stores; discard
/// and regenerate any time. Validation ("every edge used by exactly two faces
/// in opposite directions") and adjacency queries read this.
///
/// Uses plain `Vec` rather than the design's `SmallVec<[_; 2]>` / `[_; 4]` for
/// now — one fewer dependency; the inline-storage win is a future perf tweak.
#[derive(Clone, Debug, Default)]
pub struct Adjacency {
    /// Edge → every `(face, half-edge.forward)` that uses it. A manifold edge
    /// has exactly two entries with opposite `bool`s.
    pub edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>>,
    /// Vertex → the edges bounded by it (either endpoint).
    pub vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>>,
}
