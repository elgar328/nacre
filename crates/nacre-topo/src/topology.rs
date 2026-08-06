//! b-rep topology cells (design.md §4).
//!
//! Every cell references exact geometry only by `Handle` — the topology never
//! looks at coordinates. Identity is by `Handle`, so nothing here derives
//! `Eq`/`Hash` if it transitively holds an `f64` (a `Point3` or a
//! `Origin::Discovered { tol, definition }`); those get `PartialEq` for tests only.

use crate::{Orientation, Origin, VertexDef};
use nacre_geom::{Curve, Surface};
use nacre_math::Point3;
use nacre_store::Handle;

/// A 0-cell. The point is inlined — there is no point store (design §2: two
/// vertices never share a point by design, and `Point3` is small).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vertex {
    pub point: Point3,
    pub origin: Origin,
    /// **The three surfaces that meet here** — what actually decides where this point is, as
    /// opposed to `point`, which is the rounded answer to that question.
    ///
    /// ★ Three planes meet in at most one point, so naming them names the vertex exactly, with no
    /// tolerance anywhere; `point` is then a cache that could in principle be recomputed to any
    /// precision. That inversion — definition first, coordinate second — is what
    /// `docs/truth-and-cache.md` builds toward, and this field is the first half of it. Nothing
    /// reads it yet.
    ///
    /// `None` where a producer could not name three planes: a curved surface (`ThreePlane` cannot
    /// speak about a cylinder), or a transform whose re-pointing found a surface it could not
    /// map. (A profile with a collinear vertex used to be a third cause — its two walls are one
    /// plane, so the triple names a line and not a point — but since S3 the profile constructor
    /// dissolves such corners, so that population no longer reaches the topology.) The rest are
    /// counted rather than rejected — a missing definition costs nothing today, and the count is
    /// what says whether `point` can ever be dropped.
    pub definition: Option<VertexDef>,
}

/// A 1-cell: an unbounded [`Curve`] trimmed by its two endpoint vertices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    pub curve: Handle<Curve>,
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
    /// uses the seam edge twice) — a provisional spelling until M6 decides the seam's carrier
    /// representation together with the cylinder's truth (`docs/truth-and-cache.md` open item 5).
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
            // A hand-built cell; nothing here names three planes.
            definition: None,
        });
        let v1 = m.vertices.push(Vertex {
            point: Point3::from_array([1.0, 0.0, 0.0]),
            origin: Origin::Constructed,
            // A hand-built cell; nothing here names three planes.
            definition: None,
        });
        let curve = m.curves.push(Curve::Line(
            Line::through_points(Point3::origin(), Point3::from_array([1.0, 0.0, 0.0])).unwrap(),
        ));
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
        let e = m.edges.push(Edge {
            curve,
            surfaces: Edge::carrier_pair(sa, sb),
            vertices: [v0, v1],
        });

        let stored = *m.edges.get(e);
        assert_eq!(stored.curve, curve);
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
