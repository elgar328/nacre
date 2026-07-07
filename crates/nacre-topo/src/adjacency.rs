//! `Adjacency` — the reverse index over the topology stores.

use crate::Model;
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

impl Adjacency {
    /// Build a fresh index from a model's topology stores. Scans every face's
    /// outer + inner loops → half-edges → edges (recording `forward`), and every
    /// bounded edge → its two endpoint vertices. From scratch — cheap, always
    /// correct; incremental maintenance lands with the operation log (M2).
    pub fn rebuild(model: &Model) -> Adjacency {
        let mut edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>> = HashMap::new();
        for (fh, face) in model.faces.iter() {
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    edge_uses.entry(he.edge).or_default().push((fh, he.forward));
                }
            }
        }
        let mut vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>> = HashMap::new();
        for (eh, edge) in model.edges.iter() {
            if let Some([a, b]) = edge.bounds {
                vertex_edges.entry(a).or_default().push(eh);
                vertex_edges.entry(b).or_default().push(eh);
            }
        }
        Adjacency {
            edge_uses,
            vertex_edges,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::topology::{Edge, Face, HalfEdge, Loop};
    use crate::{Model, Orientation, Origin, Vertex};
    use nacre_geom::{Curve, Line, Plane, Surface};
    use nacre_math::Point3;

    /// Two triangles sharing edge v0–v1, wound so the shared edge is used once
    /// forward and once reversed. Exercises `rebuild` independent of the cube.
    #[test]
    fn rebuild_indexes_a_shared_edge() {
        let mut m = Model::new();
        let mk_v = |m: &mut Model, p: [f64; 3]| {
            m.vertices.push(Vertex {
                point: Point3::from_array(p),
                origin: Origin::Constructed,
            })
        };
        let v0 = mk_v(&mut m, [0.0, 0.0, 0.0]);
        let v1 = mk_v(&mut m, [1.0, 0.0, 0.0]);
        let v2 = mk_v(&mut m, [0.0, 1.0, 0.0]);
        let v3 = mk_v(&mut m, [0.0, -1.0, 0.0]);

        // rebuild never dereferences curve/surface handles, so one dummy of each
        // is enough for every edge/face.
        let curve = m.curves.push(Curve::Line(
            Line::through_points(Point3::origin(), Point3::from_array([1.0, 0.0, 0.0])).unwrap(),
        ));
        let surface = m.surfaces.push(Surface::Plane(
            Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
        ));
        let mk_e = |m: &mut Model, a, b| {
            m.edges.push(Edge {
                curve,
                bounds: Some([a, b]),
                origin: Origin::Constructed,
            })
        };
        let e_shared = mk_e(&mut m, v0, v1);
        let ea1 = mk_e(&mut m, v1, v2);
        let ea2 = mk_e(&mut m, v2, v0);
        let eb1 = mk_e(&mut m, v1, v3);
        let eb2 = mk_e(&mut m, v3, v0);

        let mk_face = |m: &mut Model, hes: Vec<HalfEdge>| {
            m.faces.push(Face {
                surface,
                outer: Loop { half_edges: hes },
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Triangle A: v0→v1→v2→v0, shared edge forward.
        mk_face(
            &mut m,
            vec![
                HalfEdge {
                    edge: e_shared,
                    forward: true,
                },
                HalfEdge {
                    edge: ea1,
                    forward: true,
                },
                HalfEdge {
                    edge: ea2,
                    forward: true,
                },
            ],
        );
        // Triangle B: v1→v0→v3→v1, shared edge reversed.
        mk_face(
            &mut m,
            vec![
                HalfEdge {
                    edge: e_shared,
                    forward: false,
                },
                HalfEdge {
                    edge: eb2,
                    forward: false,
                },
                HalfEdge {
                    edge: eb1,
                    forward: false,
                },
            ],
        );

        m.rebuild_adjacency();

        // shared edge: two uses, opposite flags.
        let shared = &m.adj.edge_uses[&e_shared];
        assert_eq!(shared.len(), 2);
        assert_ne!(shared[0].1, shared[1].1);
        // a private edge: one use.
        assert_eq!(m.adj.edge_uses[&ea1].len(), 1);
        // v0 is incident to e_shared, ea2, eb2.
        assert_eq!(m.adj.vertex_edges[&v0].len(), 3);
        assert_eq!(m.adj.vertex_edges[&v2].len(), 2);
    }
}
