use super::nonmanifold_vertices;
use crate::topology::{Edge, Face, HalfEdge, Loop};
use crate::{Handle, Model, Orientation, PointCache, Shell, Solid, Vertex};
use nacre_geom::Plane;
use nacre_math::Point3;
use std::collections::HashMap;

/// The pinch detector is per-vertex: a single face-fan at a vertex is manifold, two+ fans is a
/// pinch — and **every** pinch is flagged independently, so an *even* number of them is caught
/// (which an Euler-parity count cannot: two pinches flip χ back to even). Builds the reverse-index
/// maps directly with minted handles (the detector reads only the maps, never the stores).
#[test]
fn nonmanifold_vertices_flags_each_pinch_even_count() {
    let mut m = Model::new();
    let r = nacre_scalar::Rat::from_int;
    let surface = m.push_plane_unregistered(
        Plane::from_point_normal(
            Point3::origin(),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap(),
        [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
    );
    // Hand-built cells: the def names the three world seeds — real, distinct planes the
    // detector never dereferences (it reads only the maps).
    let mk_v = |m: &mut Model| {
        let def = crate::Vertex::ThreePlane([
            m.world_plane(nacre_scalar::Axis::Z),
            m.world_plane(nacre_scalar::Axis::X),
            m.world_plane(nacre_scalar::Axis::Y),
        ]);
        m.push_vertex(
            def,
            PointCache::Unrealized {
                coord: Point3::origin(),
            },
        )
    };
    let mut vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>> = HashMap::new();
    let mut edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>> = HashMap::new();
    // Add one triangular face-fan (3 edges from `apex`, 3 faces cyclically joining them).
    let add_fan =
        |m: &mut Model,
         apex: Handle<Vertex>,
         vertex_edges: &mut HashMap<Handle<Vertex>, Vec<Handle<Edge>>>,
         edge_uses: &mut HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>>| {
            let e: Vec<Handle<Edge>> = (0..3)
                .map(|_| {
                    let other = mk_v(m);
                    // The detector reads only the maps — a self-pair of the one dummy
                    // surface is enough (never dereferenced as adjacency here; a raw
                    // `edges.push` skips the curve cache, which nothing here reads).
                    m.edges.push(Edge {
                        surfaces: [surface, surface],
                        vertices: [apex, other],
                    })
                })
                .collect();
            vertex_edges
                .entry(apex)
                .or_default()
                .extend(e.iter().copied());
            for i in 0..3 {
                let f = m.faces.push(Face {
                    surface,
                    outer: Loop { half_edges: vec![] },
                    inner: vec![],
                    orientation: Orientation::Forward,
                });
                edge_uses.entry(e[i]).or_default().push((f, true));
                edge_uses.entry(e[(i + 1) % 3]).or_default().push((f, true));
            }
        };
    // Manifold apex: one fan. Two pinched apexes: two fans each.
    let ok = mk_v(&mut m);
    add_fan(&mut m, ok, &mut vertex_edges, &mut edge_uses);
    let pinch1 = mk_v(&mut m);
    add_fan(&mut m, pinch1, &mut vertex_edges, &mut edge_uses);
    add_fan(&mut m, pinch1, &mut vertex_edges, &mut edge_uses);
    let pinch2 = mk_v(&mut m);
    add_fan(&mut m, pinch2, &mut vertex_edges, &mut edge_uses);
    add_fan(&mut m, pinch2, &mut vertex_edges, &mut edge_uses);

    let mut got = nonmanifold_vertices(&vertex_edges, &edge_uses);
    got.sort_unstable_by_key(|h| h.index());
    let mut want = vec![pinch1, pinch2];
    want.sort_unstable_by_key(|h| h.index());
    assert_eq!(got, want, "both pinches flagged, the manifold fan is not");
}

/// Two triangles sharing edge v0–v1, wound so the shared edge is used once
/// forward and once reversed. Exercises `rebuild` independent of the cube.
#[test]
fn rebuild_indexes_a_shared_edge() {
    let mut m = Model::new();
    let mk_v = |m: &mut Model, p: [f64; 3]| {
        let def = crate::Vertex::ThreePlane([
            m.world_plane(nacre_scalar::Axis::Z),
            m.world_plane(nacre_scalar::Axis::X),
            m.world_plane(nacre_scalar::Axis::Y),
        ]);
        m.push_vertex(
            def,
            PointCache::Unrealized {
                coord: Point3::from_array(p),
            },
        )
    };
    let v0 = mk_v(&mut m, [0.0, 0.0, 0.0]);
    let v1 = mk_v(&mut m, [1.0, 0.0, 0.0]);
    let v2 = mk_v(&mut m, [0.0, 1.0, 0.0]);
    let v3 = mk_v(&mut m, [0.0, -1.0, 0.0]);

    // rebuild never dereferences surface handles, so one dummy is enough for
    // every edge/face (raw `edges.push` — nothing here reads the curve cache).
    let r = nacre_scalar::Rat::from_int;
    let surface = m.push_plane_unregistered(
        Plane::from_point_normal(
            Point3::origin(),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap(),
        [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
    );
    let mk_e = |m: &mut Model, a, b| {
        m.edges.push(Edge {
            // rebuild never dereferences carriers either — the dummy pair mirrors the
            // dummy surface above.
            surfaces: [surface, surface],
            vertices: [a, b],
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
    let fa = mk_face(
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
    let fb = mk_face(
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

    // Wrap the two faces in a live solid — adjacency now indexes only the
    // reachable closure, so loose faces would be invisible. This
    // pair is an open surface (validate would flag it), but the test only
    // inspects the rebuilt adjacency, which needs the faces reachable.
    let shell = m.shells.push(Shell {
        faces: vec![fa, fb],
    });
    m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
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
