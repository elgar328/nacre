//! b-rep topology cells (design.md §4).
//!
//! Every cell references exact geometry only by `Handle` — the topology never
//! looks at coordinates. Identity is by `Handle`, so a cell derives `Eq`/`Hash` only
//! when every field is a handle or a flag; anything that could transitively hold an
//! `f64` gets `PartialEq` for tests only. (No cell holds one today — the coordinate
//! became a cache in S7 — but the rule is about what a cell is *allowed* to hold.)

use crate::{Orientation, Surface, VertexDef};
use nacre_store::Handle;

/// A 0-cell: **its definition is all it is** (S7). The realized coordinate and its measured
/// tolerance live in the index-parallel point cache (`Model::vertex_point` /
/// `Model::vertex_tol`), filled by `Model::push_vertex` — definition first, coordinate second,
/// the inversion `docs/truth-and-cache.md` builds toward. `Origin` (Constructed / Discovered /
/// Moved) died here: the tag never meant exactness, the measured tolerance moved to the cache,
/// and a moved vertex's motion was always the *faces'* motion (they record it themselves).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vertex {
    pub def: VertexDef,
}

/// A 1-cell: a segment of the carriers' intersection, trimmed by its two endpoint vertices.
/// The realized curve is a cache beside the store (`Model::edge_cache`), not a field — the
/// carriers and endpoints decide it (S8).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    /// **The two surfaces whose faces this edge bounds** — the carriers (S8).
    ///
    /// ★ Not derivable from the endpoints' surface triples: where four planes pass through one
    /// point, the triples' intersection can name a plane the edge does *not* ride — silently
    /// (the four-plane concurrency trap, `nacre-ops`' `RingEdge`/`Ring.walls` doc). Every
    /// producer states the pair from the faces it is building, never derives it.
    ///
    /// ★ Nor is it "every plane containing the line": where three planes share a LINE (the
    /// resolved four-plane concurrency), a third plane legitimately contains the edge without
    /// bounding it here — measured, each side's arrangement named that third plane as its wall.
    /// The pair is the *adjacency* answer: the two faces that use the edge.
    ///
    /// Stored in ascending handle-index order — the pair is a set, not a sequence.
    /// A cylinder seam is self-adjacent: both entries are the lateral surface (its loop already
    /// uses the seam edge twice) — the **confirmed** spelling (M6-0): a seam is a
    /// parameterization joint of one surface, and the self-pair is that sentence's honest
    /// carrier form, guarded by validate's "self-adjacent ⇔ cylinder" rule.
    pub surfaces: [Handle<Surface>; 2],
    /// Endpoint vertices — the boundary (S8: no longer `Option`).
    ///
    /// A closed **solid**'s circular rim carries a seam vertex and is `[v, v]` (start == end),
    /// so the b-rep stays a valid CW-complex (`add_cylinder`; design §4). The old `None` was
    /// reserved for a standalone full circle with no seam — a wireframe/open-shell element,
    /// which is a v1 non-goal (§9) and had exactly one occupant: a validate fixture built to
    /// test the reject. The type now says what the kernel always required.
    pub vertices: [Handle<Vertex>; 2],
}

impl Edge {
    /// The carrier pair in its stored (canonical, ascending handle-index) order — the pair is
    /// a set, and one spelling keeps `==` on edges meaningful.
    pub fn carrier_pair(a: Handle<Surface>, b: Handle<Surface>) -> [Handle<Surface>; 2] {
        if b.index() < a.index() {
            [b, a]
        } else {
            [a, b]
        }
    }
}

/// A directed use of a shared [`Edge`] inside a [`Loop`].
///
/// `forward == true` traverses the edge start→end (`vertices[0]` → `vertices[1]`);
/// `false` traverses it in reverse.
///
/// ★ **A closed rim starts where it ends**, so those vertices settle nothing for it:
/// there, `forward` reads as *along the curve's own parameterization* — CCW about the
/// circle's normal, which is the cylinder's axis direction. That is the only reading
/// available and the one every producer writes (a bottom cap's rim is `false`, and its
/// face states `−axis` as its outward). `validate`'s face-orientation check is what
/// holds producers to it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge {
    pub edge: Handle<Edge>,
    pub forward: bool,
}

/// A closed cycle of half-edges bounding (part of) a face.
#[derive(Clone, Debug, PartialEq)]
pub struct Loop {
    pub half_edges: Vec<HalfEdge>,
}

/// A 2-cell: a surface patch bounded by one outer loop and zero+ inner loops
/// (holes). `orientation` says whether the face uses the surface normal as-is
/// (`Forward`) or flipped (`Reversed`).
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    pub surface: Handle<Surface>,
    pub outer: Loop,
    pub inner: Vec<Loop>,
    pub orientation: Orientation,
}

/// A closed set of faces forming one connected boundary surface.
#[derive(Clone, Debug, PartialEq)]
pub struct Shell {
    pub faces: Vec<Handle<Face>>,
}

/// A volume = one outer shell plus zero+ inner cavity shells.
#[derive(Clone, Debug, PartialEq)]
pub struct Solid {
    pub outer: Handle<Shell>,
    pub cavities: Vec<Handle<Shell>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Model;
    use nacre_math::Point3;

    #[test]
    fn cells_round_trip_their_fields() {
        let mut m = Model::new();
        let seeds = crate::VertexDef::ThreePlane([
            m.world_plane(nacre_scalar::Axis::Z),
            m.world_plane(nacre_scalar::Axis::X),
            m.world_plane(nacre_scalar::Axis::Y),
        ]);
        let v0 = m.push_vertex(seeds, Point3::origin(), None);
        let v1 = m.push_vertex(seeds, Point3::from_array([1.0, 0.0, 0.0]), None);
        let r = nacre_scalar::Rat::from_int;
        let sa = m.push_plane_unregistered(
            nacre_geom::Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
        );
        let sb = m.push_plane_unregistered(
            nacre_geom::Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
            )
            .unwrap(),
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(0), r(1)]],
        );
        let e = m.push_edge([sb, sa], [v0, v1]).expect("distinct endpoints");

        let stored = *m.edges.get(e);
        assert_eq!(stored.surfaces, [sa, sb], "already ascending — kept as-is");
        assert_eq!(stored.vertices, [v0, v1]);

        // HalfEdge is Copy + Eq + Hash (Handle + bool only).
        let he = HalfEdge {
            edge: e,
            forward: true,
        };
        let he_same = HalfEdge {
            edge: e,
            forward: true,
        };
        assert_eq!(he, he_same);
        let lp = Loop {
            half_edges: vec![he],
        };
        assert!(lp.half_edges[0].forward);
    }
}
