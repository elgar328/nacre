//! b-rep topology cells.
//!
//! Every cell references exact geometry only by `Handle` — the topology never
//! looks at coordinates. Identity is by `Handle`, so a cell derives `Eq`/`Hash` only
//! when every field is a handle or a flag; anything that could transitively hold an
//! `f64` gets `PartialEq` for tests only. (No cell holds one today — the coordinate
//! is a cache — but the rule is about what a cell is *allowed* to hold.)

use crate::{Orientation, Surface, Vertex};
use nacre_store::Handle;

/// A 1-cell: a segment of the carriers' intersection, trimmed by its two endpoint vertices.
/// The realized curve is a cache beside the store (`Model::edge_cache`), not a field — the
/// carriers and endpoints decide it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    /// **The two surfaces whose faces this edge bounds** — the carriers.
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
    /// Stored in ascending handle-index order — the pair is a set, not a sequence. The two entries
    /// are two surfaces: an edge separates two faces' surfaces, and no edge separates a surface from
    /// itself — a cylinder's lateral is bounded by its rims, with no seam edge
    /// ([`crate::Model::push_edge`] refuses a pair that is one cylinder, `validate` flags any
    /// self-pair).
    pub surfaces: [Handle<Surface>; 2],
    /// Endpoint vertices — the boundary (never optional).
    ///
    /// A closed **solid**'s circular rim carries a seam vertex and is `[v, v]` (start == end),
    /// so the b-rep stays a valid CW-complex. A standalone full circle
    /// with no seam would be a wireframe/open-shell element, which is a v1 non-goal, so the
    /// type has no spelling for it: it says what the kernel requires.
    pub vertices: [Handle<Vertex>; 2],
}

/// **Why an edge's curve does not derive** ([`crate::Model::derive_edge_curve`]) — one name per
/// cause, so a caller can say which promise broke.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeDecline {
    /// A straight edge (two planes, a seam, a ruling) whose endpoints are one point: zero length.
    Coincident,
    /// A plane oblique to a cylinder's axis — the section is an ellipse, which no edge carries.
    Oblique,
    /// A plane and a cylinder whose truths cannot be placed in one frame, so their relation is
    /// not stated.
    Unstated,
    /// Two cylinders — two distinct ones meet in a quartic, and one stated twice is a surface
    /// separated from itself; no edge carries either.
    TwoCylinders,
    /// The truth states a circle, but the cache could not build it (no centre, no radius).
    Degenerate,
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
    use crate::{Model, PointCache};
    use nacre_math::Point3;

    #[test]
    fn cells_round_trip_their_fields() {
        let mut m = Model::new();
        let seeds = crate::Vertex::ThreePlane([
            m.world_plane(nacre_exact::Axis::Z),
            m.world_plane(nacre_exact::Axis::X),
            m.world_plane(nacre_exact::Axis::Y),
        ]);
        let v0 = m.push_vertex(
            seeds,
            PointCache::Unrealized {
                coord: Point3::origin(),
            },
        );
        let v1 = m.push_vertex(
            seeds,
            PointCache::Unrealized {
                coord: Point3::from_array([1.0, 0.0, 0.0]),
            },
        );
        let r = nacre_exact::Rat::from_int;
        let sa = m.push_plane_unregistered(
            nacre_geom::Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
            Orientation::Forward,
        );
        let sb = m.push_plane_unregistered(
            nacre_geom::Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
            )
            .unwrap(),
            // `x̂ × ẑ = −ŷ`, against the `+y` cache.
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(0), r(1)]],
            Orientation::Reversed,
        );
        let e = m
            .push_edge([sb, sa], [v0, v1], |_| crate::EdgeGiven::NONE)
            .expect("distinct endpoints");

        let stored = *m.edge(e);
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
