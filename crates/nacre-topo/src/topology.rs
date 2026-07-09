//! b-rep topology cells (design.md §4).
//!
//! Every cell references exact geometry only by `Handle` — the topology never
//! looks at coordinates. Identity is by `Handle`, so nothing here derives
//! `Eq`/`Hash` if it transitively holds an `f64` (a `Point3` or a
//! `Origin::Discovered { tol, definition }`); those get `PartialEq` for tests only.

use crate::{Orientation, Origin};
use nacre_geom::{Curve, Surface};
use nacre_math::Point3;
use nacre_store::Handle;

/// A 0-cell. The point is inlined — there is no point store (design §2: two
/// vertices never share a point by design, and `Point3` is small).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vertex {
    pub point: Point3,
    pub origin: Origin,
}

/// A 1-cell: an unbounded [`Curve`] trimmed by its two endpoint vertices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    pub curve: Handle<Curve>,
    /// Endpoint vertices, or `None` for a truly closed edge (no endpoints).
    ///
    /// A closed **solid**'s circular rim is *not* `None`: it carries a seam
    /// vertex and is `Some([v, v])` (start == end), so the b-rep stays a valid
    /// CW-complex (`add_cylinder`; design §4). `None` is reserved for a
    /// standalone full circle with no seam — a wireframe/open-shell element,
    /// which is a v1 non-goal (§9), so it is currently unused.
    pub bounds: Option<[Handle<Vertex>; 2]>,
    pub origin: Origin,
}

/// A directed use of a shared [`Edge`] inside a [`Loop`].
///
/// `forward == true` traverses the edge start→end (`bounds[0]` → `bounds[1]`);
/// `false` traverses it in reverse.
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
    use nacre_geom::{Curve, Line};

    #[test]
    fn cells_round_trip_their_fields() {
        let mut m = Model::new();
        let v0 = m.vertices.push(Vertex {
            point: Point3::origin(),
            origin: Origin::Constructed,
        });
        let v1 = m.vertices.push(Vertex {
            point: Point3::from_array([1.0, 0.0, 0.0]),
            origin: Origin::Constructed,
        });
        let curve = m.curves.push(Curve::Line(
            Line::through_points(Point3::origin(), Point3::from_array([1.0, 0.0, 0.0])).unwrap(),
        ));
        let e = m.edges.push(Edge {
            curve,
            bounds: Some([v0, v1]),
            origin: Origin::Constructed,
        });

        let stored = *m.edges.get(e);
        assert_eq!(stored.curve, curve);
        assert_eq!(stored.bounds, Some([v0, v1]));
        assert_eq!(stored.origin, Origin::Constructed);

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
