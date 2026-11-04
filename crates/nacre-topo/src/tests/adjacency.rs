use super::{Adjacency, nonmanifold_vertices};
use crate::topology::{Edge, Face, HalfEdge, Loop};
use crate::{Handle, Model, Orientation, PointCache, Shell, Solid, Vertex};
use nacre_geom::Plane;
use nacre_math::Point3;

/// Hand-built cells for the pinch detector, which reads only the faces' loops and their edges'
/// endpoints — never a surface or a coordinate, so one dummy plane carries every face and edge.
struct Cells {
    m: Model,
    surface: Handle<crate::Surface>,
}

impl Cells {
    fn new() -> Cells {
        let mut m = Model::new();
        let r = nacre_exact::Rat::from_int;
        let surface = m.push_plane_unregistered(
            Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
            Orientation::Forward,
        );
        Cells { m, surface }
    }

    fn vertex(&mut self) -> Handle<Vertex> {
        let def = crate::Vertex::ThreePlane([
            self.m.seeded_plane(nacre_exact::Axis::Z),
            self.m.seeded_plane(nacre_exact::Axis::X),
            self.m.seeded_plane(nacre_exact::Axis::Y),
        ]);
        self.m.push_vertex(
            def,
            PointCache::Unrealized {
                coord: Point3::origin(),
            },
        )
    }

    fn edge(&mut self, a: Handle<Vertex>, b: Handle<Vertex>) -> Handle<Edge> {
        self.m.edges.push(Edge {
            surfaces: [self.surface, self.surface],
            vertices: [a, b],
        })
    }

    /// A face with one loop of `(edge, forward)` half-edges.
    fn face(&mut self, hes: &[(Handle<Edge>, bool)]) -> Handle<Face> {
        self.m.faces.push(Face {
            surface: self.surface,
            outer: Loop {
                half_edges: hes
                    .iter()
                    .map(|&(edge, forward)| HalfEdge { edge, forward })
                    .collect(),
            },
            inner: vec![],
            orientation: Orientation::Forward,
        })
    }

    /// A closed fan of three triangles about `apex` — `apex → o_i → o_{i+1} → apex`.
    fn fan(&mut self, apex: Handle<Vertex>, faces: &mut Vec<Handle<Face>>) {
        let o: Vec<Handle<Vertex>> = (0..3).map(|_| self.vertex()).collect();
        let spoke: Vec<Handle<Edge>> = o.iter().map(|&v| self.edge(apex, v)).collect();
        let rim: Vec<Handle<Edge>> = (0..3).map(|i| self.edge(o[i], o[(i + 1) % 3])).collect();
        for i in 0..3 {
            faces.push(self.face(&[
                (spoke[i], true),
                (rim[i], true),
                (spoke[(i + 1) % 3], false),
            ]));
        }
    }

    fn pinches(&self, faces: &[Handle<Face>]) -> Vec<Handle<Vertex>> {
        nonmanifold_vertices(&self.m, faces.iter().copied())
    }
}

/// **The detector is per vertex**: one fan of faces at a vertex is manifold, two fans a pinch — and
/// every pinch is flagged on its own, so an *even* number of them is caught (an Euler-parity count
/// is not: two pinches flip χ back to even).
#[test]
fn nonmanifold_vertices_flags_each_pinch_even_count() {
    let mut c = Cells::new();
    let mut faces = Vec::new();
    let ok = c.vertex();
    c.fan(ok, &mut faces);
    let pinch1 = c.vertex();
    c.fan(pinch1, &mut faces);
    c.fan(pinch1, &mut faces);
    let pinch2 = c.vertex();
    c.fan(pinch2, &mut faces);
    c.fan(pinch2, &mut faces);
    assert_eq!(
        c.pinches(&faces),
        vec![pinch1, pinch2],
        "both pinches flagged, the manifold fan is not"
    );
}

/// **A closed edge is two ends of the link.** A whole rim `[v, v]` bounds a disk cap on one side and
/// a lateral on the other: the two faces join the rim's two ends at `v`, one component, manifold.
/// A fan of other faces meeting `v` as well — a box whose corner rests on the rim — is a second
/// component: a pinch, where a link read in *edges* sees a face using one edge at `v` and cannot
/// say.
#[test]
fn a_fan_meeting_a_closed_rim_at_its_vertex_is_a_pinch() {
    let mut c = Cells::new();
    let v = c.vertex();
    let rim = c.edge(v, v);
    let mut faces = vec![c.face(&[(rim, true)]), c.face(&[(rim, false)])];
    assert_eq!(
        c.pinches(&faces),
        vec![],
        "a cap and its lateral are one fan"
    );
    c.fan(v, &mut faces);
    assert_eq!(
        c.pinches(&faces),
        vec![v],
        "a second fan at the rim's vertex"
    );
}

/// **A face passing a vertex twice in one loop reads as its two corners.** The loop
/// `w →a→ v →rim→ v →b→ w` takes a closed edge between two others — a plate's top whose ring runs
/// round a hole's rim at the point a boss's corner touches it. With the lateral on the rim's other
/// side and a face closing `a` and `b`, the link at `v` is one cycle through the four ends: manifold.
/// Reading the ends' slots the wrong way round splits `v` into two components — the rim's two
/// ends, and `a`'s and `b`'s — a false pinch.
#[test]
fn a_closed_edge_inside_a_longer_loop_is_read_by_its_ends() {
    let mut c = Cells::new();
    let (v, w) = (c.vertex(), c.vertex());
    let (a, rim, b) = (c.edge(w, v), c.edge(v, v), c.edge(v, w));
    let mut faces = vec![
        c.face(&[(a, true), (rim, true), (b, true)]),
        c.face(&[(rim, false)]),
        c.face(&[(b, false), (a, false)]),
    ];
    assert_eq!(c.pinches(&faces), vec![], "one fan at v, one at w");
    c.fan(v, &mut faces);
    assert_eq!(c.pinches(&faces), vec![v], "a second fan at v");
}

/// Two triangles sharing edge v0–v1, wound so the shared edge is used once
/// forward and once reversed. Exercises `rebuild` independent of the cube.
#[test]
fn rebuild_indexes_a_shared_edge() {
    let mut m = Model::new();
    let mk_v = |m: &mut Model, p: [f64; 3]| {
        let def = crate::Vertex::ThreePlane([
            m.seeded_plane(nacre_exact::Axis::Z),
            m.seeded_plane(nacre_exact::Axis::X),
            m.seeded_plane(nacre_exact::Axis::Y),
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
    let r = nacre_exact::Rat::from_int;
    let surface = m.push_plane_unregistered(
        Plane::from_point_normal(
            Point3::origin(),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap(),
        [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
        Orientation::Forward,
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

    // Wrap the two faces in a live solid — adjacency indexes only the reachable closure, so loose
    // faces are invisible. This
    // pair is an open surface (validate would flag it), but the test only
    // inspects the rebuilt adjacency, which needs the faces reachable.
    let shell = m.shells.push(Shell {
        faces: vec![fa, fb],
    });
    m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    let adj = Adjacency::rebuild(&m);

    // shared edge: two uses, opposite flags.
    let shared = &adj.edge_uses[&e_shared];
    assert_eq!(shared.len(), 2);
    assert_ne!(shared[0].1, shared[1].1);
    // a private edge: one use.
    assert_eq!(adj.edge_uses[&ea1].len(), 1);
    // v0 is incident to e_shared, ea2, eb2.
    assert_eq!(adj.vertex_edges[&v0].len(), 3);
    assert_eq!(adj.vertex_edges[&v2].len(), 2);
}
