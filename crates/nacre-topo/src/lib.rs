//! b-rep topology and the truth-only `Model` aggregate (design.md §2, §4).
//!
//! The topology (`Vertex`/`Edge`/`Face`/`Loop`/`HalfEdge`/`Shell`/`Solid`)
//! references exact geometry only by `Handle` — geometry never knows about
//! topology, topology never inspects coordinates (overview 절대원칙 2).
//!
//! [`Model`] is the **truth**: exact geometry stores + topology stores + the
//! derived [`Adjacency`] cache. It holds no tessellation and no operation log —
//! a mesh cache and the op log are companions owned at higher layers (§0: "the
//! model is the replay result"; putting them here would make topo depend on
//! tess/ops and break the truth/cache split).

mod adjacency;
mod topology;

pub use adjacency::Adjacency;
pub use topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

use nacre_geom::{Curve, Surface};
use nacre_store::Store;

/// Provenance of a vertex or edge (design §4, overview 절대원칙 4).
///
/// `Constructed` elements know their identity by construction and carry no
/// tolerance; `Discovered` elements come from an intersection and hold the
/// *measured* accuracy relaxation achieved (M3+). In M1 every element is
/// `Constructed`, so the tolerance path is never exercised.
///
/// Holds an `f64`, so `PartialEq` only — no `Eq`/`Hash` (identity is by
/// `Handle`, never by value).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Origin {
    Constructed,
    Discovered { tol: f64 },
}

/// Whether a face uses its surface normal as-is (`Forward`) or flipped
/// (`Reversed`). A pure tag — full derives.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reversed,
}

/// The truth-only aggregate: exact geometry + topology stores + the derived
/// adjacency cache. No `tess`, no `ops` (see the crate docs).
#[derive(Debug, Default)]
pub struct Model {
    // exact geometry (truth)
    pub surfaces: Store<Surface>,
    pub curves: Store<Curve>,
    // topology (references geometry by Handle only)
    pub vertices: Store<Vertex>,
    pub edges: Store<Edge>,
    pub faces: Store<Face>,
    pub shells: Store<Shell>,
    pub solids: Store<Solid>,
    // derived cache (rebuilt on demand)
    pub adj: Adjacency,
}

impl Model {
    /// An empty model.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Recompute the [`Adjacency`] cache from the current topology stores.
    /// Call once after a batch of additions (the cache is otherwise stale).
    pub fn rebuild_adjacency(&mut self) {
        // Build against an immutable borrow, then move into place — avoids
        // borrowing `self.adj` mutably while iterating the other stores.
        let adj = Adjacency::rebuild(&*self);
        self.adj = adj;
    }
}
