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

use nacre_geom::{Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_store::{Handle, Store};

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

    /// Add an axis-aligned box `min`..`max` to this model and return its solid.
    ///
    /// Requires `max[i] > min[i]` on every axis (a degenerate box is a caller
    /// bug → panic). Every element is [`Origin::Constructed`] (design §8); each
    /// face is wound so its plane normal points outward, so all faces are
    /// [`Orientation::Forward`]. Does **not** rebuild adjacency — call
    /// [`Model::rebuild_adjacency`] once after all additions.
    pub fn add_cuboid(&mut self, min: Point3, max: Point3) -> Handle<Solid> {
        let [x0, y0, z0] = min.as_array();
        let [x1, y1, z1] = max.as_array();
        debug_assert!(
            x1 > x0 && y1 > y0 && z1 > z0,
            "cuboid needs max > min on every axis"
        );

        // 8 corners, numbered by bits (i·x, j·y, k·z).
        let corners = [
            Point3::from_array([x0, y0, z0]), // V0
            Point3::from_array([x1, y0, z0]), // V1
            Point3::from_array([x1, y1, z0]), // V2
            Point3::from_array([x0, y1, z0]), // V3
            Point3::from_array([x0, y0, z1]), // V4
            Point3::from_array([x1, y0, z1]), // V5
            Point3::from_array([x1, y1, z1]), // V6
            Point3::from_array([x0, y1, z1]), // V7
        ];
        let vh: [Handle<Vertex>; 8] = core::array::from_fn(|i| {
            self.vertices.push(Vertex {
                point: corners[i],
                origin: Origin::Constructed,
            })
        });

        // 12 edges as (start, end) vertex indices: bottom ring, top ring, verticals.
        const EDGES: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let eh: [Handle<Edge>; 12] = core::array::from_fn(|i| {
            let (a, b) = EDGES[i];
            let curve = self.curves.push(Curve::Line(
                Line::through_points(corners[a], corners[b]).expect("non-degenerate box"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([vh[a], vh[b]]),
                origin: Origin::Constructed,
            })
        });

        // 6 faces: (first 3 loop vertices for the outward plane, half-edges as
        // (edge index, forward)). Loops wound CCW seen from outside → outward
        // normal. See the design plan's winding table (hand-verified).
        type FaceDef = ([usize; 3], [(usize, bool); 4]);
        let faces_def: [FaceDef; 6] = [
            ([0, 3, 2], [(3, false), (2, false), (1, false), (0, false)]), // Bottom −Z
            ([4, 5, 6], [(4, true), (5, true), (6, true), (7, true)]),     // Top +Z
            ([0, 1, 5], [(0, true), (9, true), (4, false), (8, false)]),   // Front −Y
            ([2, 3, 7], [(2, true), (11, true), (6, false), (10, false)]), // Back +Y
            ([0, 4, 7], [(8, true), (7, false), (11, false), (3, true)]),  // Left −X
            ([1, 2, 6], [(1, true), (10, true), (5, false), (9, false)]),  // Right +X
        ];
        let fh: [Handle<Face>; 6] = core::array::from_fn(|i| {
            let (tri, hes) = &faces_def[i];
            let surface = self.surfaces.push(Surface::Plane(
                Plane::through_points(corners[tri[0]], corners[tri[1]], corners[tri[2]])
                    .expect("non-degenerate box"),
            ));
            let outer = Loop {
                half_edges: hes
                    .iter()
                    .map(|&(e, forward)| HalfEdge {
                        edge: eh[e],
                        forward,
                    })
                    .collect(),
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        });

        let shell = self.shells.push(Shell { faces: fh.to_vec() });
        self.solids.push(Solid {
            outer: shell,
            cavities: vec![],
        })
    }

    /// Add a closed cylinder solid and return it: the axis runs from `base` along
    /// `axis` for `height`, with the given `radius`.
    ///
    /// Built as the standard **seam b-rep** (V2/E3/F3): two seam vertices, two
    /// full-circle rim edges (`bounds: Some([seam, seam])`, start == end), one
    /// straight seam edge, a cylindrical lateral face whose loop uses the seam
    /// edge twice (opposite orientation), and two planar caps (each a
    /// single-half-edge rim loop). This forms a valid CW-complex (Euler χ = 2)
    /// that [`validate`](../nacre_validate/fn.validate.html) accepts — a closed
    /// periodic surface needs a seam vertex, so the rims are `Some([v, v])`, not
    /// `bounds: None` (that form is for a standalone full circle; §4).
    ///
    /// `radius`/`height` must be positive and `axis` nonzero (caller bug →
    /// panic). Every element is [`Origin::Constructed`]. Does **not** rebuild
    /// adjacency — call [`Model::rebuild_adjacency`] once after all additions.
    pub fn add_cylinder(
        &mut self,
        base: Point3,
        axis: Vector3,
        radius: f64,
        height: f64,
    ) -> Handle<Solid> {
        debug_assert!(
            radius > 0.0 && height > 0.0,
            "cylinder needs positive radius and height"
        );
        let d = axis.normalize().expect("cylinder axis must be nonzero");
        let u = d
            .any_perpendicular()
            .expect("a unit axis has a perpendicular");
        let c0 = base;
        let c1 = base + d * height;
        let p_bot = c0 + u * radius; // seam point on the bottom rim (angle 0)
        let p_top = c1 + u * radius; // seam point on the top rim

        let v_bot = self.vertices.push(Vertex {
            point: p_bot,
            origin: Origin::Constructed,
        });
        let v_top = self.vertices.push(Vertex {
            point: p_top,
            origin: Origin::Constructed,
        });

        // Rims are full circles seamed at their vertex (start == end); the seam is
        // a straight edge joining the two rim seam points.
        let bottom = {
            let curve = self.curves.push(Curve::Circle(
                Circle::from_center_normal(c0, d, u, radius).expect("non-degenerate rim"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([v_bot, v_bot]),
                origin: Origin::Constructed,
            })
        };
        let top = {
            let curve = self.curves.push(Curve::Circle(
                Circle::from_center_normal(c1, d, u, radius).expect("non-degenerate rim"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([v_top, v_top]),
                origin: Origin::Constructed,
            })
        };
        let seam = {
            let curve = self.curves.push(Curve::Line(
                Line::through_points(p_bot, p_top).expect("positive height"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([v_bot, v_top]),
                origin: Origin::Constructed,
            })
        };

        // Lateral cylindrical face: one loop wrapping the seam twice (opposite).
        let lateral = {
            let surface = self.surfaces.push(Surface::Cylinder(
                Cylinder::from_axis(c0, d, u, radius).expect("non-degenerate cylinder"),
            ));
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: bottom,
                        forward: true,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: true,
                    },
                    HalfEdge {
                        edge: top,
                        forward: false,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: false,
                    },
                ],
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Bottom cap: outward normal −d, the bottom rim reversed.
        let bottom_cap = {
            let surface = self.surfaces.push(Surface::Plane(
                Plane::from_point_normal(c0, -d).expect("nonzero axis"),
            ));
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: bottom,
                    forward: false,
                }],
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Top cap: outward normal +d, the top rim forward.
        let top_cap = {
            let surface = self.surfaces.push(Surface::Plane(
                Plane::from_point_normal(c1, d).expect("nonzero axis"),
            ));
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: top,
                    forward: true,
                }],
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };

        let shell = self.shells.push(Shell {
            faces: vec![lateral, bottom_cap, top_cap],
        });
        self.solids.push(Solid {
            outer: shell,
            cavities: vec![],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;
    use proptest::prelude::*;

    fn build(min: [f64; 3], max: [f64; 3]) -> Model {
        let mut m = Model::new();
        m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
        m.rebuild_adjacency();
        m
    }

    fn he_start(m: &Model, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = m.edges.get(he.edge).bounds.unwrap();
        if he.forward { a } else { b }
    }
    fn he_end(m: &Model, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = m.edges.get(he.edge).bounds.unwrap();
        if he.forward { b } else { a }
    }
    fn face_plane_normal(m: &Model, f: &Face) -> Vector3 {
        match m.surfaces.get(f.surface) {
            Surface::Plane(p) => p.normal(),
            // Planar-only helper: callers filter to plane faces (caps), never cylinders.
            Surface::Cylinder(_) => unreachable!("face_plane_normal called on a curved face"),
        }
    }
    fn face_centroid(m: &Model, f: &Face) -> Point3 {
        let pts: Vec<Point3> = f
            .outer
            .half_edges
            .iter()
            .map(|he| m.vertices.get(he_start(m, *he)).point)
            .collect();
        Point3::centroid(&pts).unwrap()
    }

    // --- golden ---

    #[test]
    fn cuboid_counts() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.edges.len(), 12);
        assert_eq!(m.faces.len(), 6);
        assert_eq!(m.shells.len(), 1);
        assert_eq!(m.solids.len(), 1);
        assert_eq!(m.surfaces.len(), 6);
        assert_eq!(m.curves.len(), 12);
    }

    #[test]
    fn cuboid_corner_points() {
        let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        let expected = vec![
            [-2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 4.0, 0.0],
            [-2.0, 4.0, 0.0],
            [-2.0, 1.0, 10.0],
            [3.0, 1.0, 10.0],
            [3.0, 4.0, 10.0],
            [-2.0, 4.0, 10.0],
        ];
        let got: Vec<[f64; 3]> = m.vertices.iter().map(|(_, v)| v.point.as_array()).collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn face_normals_point_outward() {
        let m = build([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let center = Point3::origin().lerp(Point3::from_array([2.0, 3.0, 4.0]), 0.5);
        for (_, f) in m.faces.iter() {
            let outward = face_plane_normal(&m, f).dot(face_centroid(&m, f) - center);
            assert!(outward > 0.0);
        }
    }

    #[test]
    fn every_edge_used_twice_opposite() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.adj.edge_uses.len(), 12);
        for uses in m.adj.edge_uses.values() {
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1);
        }
    }

    #[test]
    fn every_vertex_incident_to_three_edges() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.adj.vertex_edges.len(), 8);
        for edges in m.adj.vertex_edges.values() {
            assert_eq!(edges.len(), 3);
        }
    }

    #[test]
    fn outer_loops_are_closed() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        for (_, f) in m.faces.iter() {
            let hes = &f.outer.half_edges;
            assert_eq!(hes.len(), 4);
            for i in 0..hes.len() {
                assert_eq!(he_end(&m, hes[i]), he_start(&m, hes[(i + 1) % hes.len()]));
            }
        }
    }

    #[test]
    fn edge_endpoints_lie_on_their_curve() {
        let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        for (_, e) in m.edges.iter() {
            let curve = m.curves.get(e.curve);
            let [a, b] = e.bounds.unwrap();
            assert!(curve.contains(m.vertices.get(a).point, 1e-9));
            assert!(curve.contains(m.vertices.get(b).point, 1e-9));
        }
    }

    #[test]
    fn euler_poincare_holds() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (v, e, f) = (
            m.vertices.len() as i64,
            m.edges.len() as i64,
            m.faces.len() as i64,
        );
        // V − E + F = 2(S − G) + L_i, with S=1, G=0, L_i=0. Formal validate:
        // nacre-validate (next unit).
        assert_eq!(v - e + f, 2);
    }

    // --- cylinder (seam b-rep) --- (validate-clean lives in nacre-validate)

    fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
        let mut m = Model::new();
        m.add_cylinder(Point3::from_array(base), Vector3::from_array(axis), r, h);
        m.rebuild_adjacency();
        m
    }

    #[test]
    fn cylinder_counts_and_euler() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        assert_eq!(m.vertices.len(), 2);
        assert_eq!(m.edges.len(), 3);
        assert_eq!(m.faces.len(), 3);
        assert_eq!(m.shells.len(), 1);
        assert_eq!(m.solids.len(), 1);
        // Euler χ = V − E + F = 2 (one shell, genus 0, no inner loops).
        assert_eq!(
            m.vertices.len() as i64 - m.edges.len() as i64 + m.faces.len() as i64,
            2
        );
    }

    #[test]
    fn cylinder_geometry() {
        // +Z axis, r=2, h=5. Seam direction is X (least-aligned axis of +Z).
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let pts: Vec<[f64; 3]> = m.vertices.iter().map(|(_, v)| v.point.as_array()).collect();
        // Seam direction for +Z is any_perpendicular([0,0,1]) = X×Z = [0,-1,0], so
        // the seam vertices sit at radius 2 along −Y, at z=0 and z=5.
        assert_eq!(pts, vec![[0.0, -2.0, 0.0], [0.0, -2.0, 5.0]]);
        // Two rim circles carry a Circle; the straight seam a Line.
        let mut circles = 0;
        for (_, e) in m.edges.iter() {
            if let Curve::Circle(c) = m.curves.get(e.curve) {
                assert_eq!(c.radius(), 2.0);
                circles += 1;
            }
        }
        assert_eq!(circles, 2);
    }

    #[test]
    fn cylinder_caps_point_outward() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        // Planar caps only (the lateral cylindrical face has no single normal).
        let mut caps = 0;
        for (_, f) in m.faces.iter() {
            if matches!(m.surfaces.get(f.surface), Surface::Plane(_)) {
                let n = face_plane_normal(&m, f).as_array();
                // Bottom cap → −Z, top cap → +Z (outward along the axis).
                assert!(n == [0.0, 0.0, -1.0] || n == [0.0, 0.0, 1.0]);
                caps += 1;
            }
        }
        assert_eq!(caps, 2);
    }

    #[test]
    fn cylinder_seam_edge_is_self_adjacent() {
        // The novel topology: the seam edge is used twice by the SAME lateral
        // face with opposite orientation (a valid non-manifold-looking but
        // manifold seam). Every edge is still used exactly twice, opposite.
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let mut seam_uses = None;
        for (eh, e) in m.edges.iter() {
            let uses = &m.adj.edge_uses[&eh];
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1); // opposite orientation
            if matches!(m.curves.get(e.curve), Curve::Line(_)) {
                seam_uses = Some(uses.clone());
            }
        }
        let uses = seam_uses.expect("a seam line edge exists");
        assert_eq!(uses[0].0, uses[1].0); // both uses are the same (lateral) face
    }

    // --- proptest ---

    fn box_strategy() -> impl Strategy<Value = ([f64; 3], [f64; 3])> {
        (
            prop::array::uniform3(-1e3f64..1e3),
            prop::array::uniform3(1e-2f64..1e3),
        )
            .prop_map(|(min, ext)| {
                let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
                (min, max)
            })
    }

    proptest! {
        #[test]
        fn prop_cuboid_structural((min, max) in box_strategy()) {
            let m = build(min, max);
            prop_assert_eq!(m.vertices.len(), 8);
            prop_assert_eq!(m.edges.len(), 12);
            prop_assert_eq!(m.faces.len(), 6);
            prop_assert_eq!(m.adj.edge_uses.len(), 12);
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
            prop_assert_eq!(m.adj.vertex_edges.len(), 8);
            for edges in m.adj.vertex_edges.values() {
                prop_assert_eq!(edges.len(), 3);
            }
        }

        #[test]
        fn prop_cuboid_outward_normals((min, max) in box_strategy()) {
            let m = build(min, max);
            let center = Point3::from_array(min).lerp(Point3::from_array(max), 0.5);
            for (_, f) in m.faces.iter() {
                prop_assert!(face_plane_normal(&m, f).dot(face_centroid(&m, f) - center) > 0.0);
            }
        }

        #[test]
        fn prop_cuboid_endpoints_on_curves((min, max) in box_strategy()) {
            let m = build(min, max);
            let scale = 1e-6 * (max.iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
            for (_, e) in m.edges.iter() {
                let curve = m.curves.get(e.curve);
                let [a, b] = e.bounds.unwrap();
                prop_assert!(curve.contains(m.vertices.get(a).point, scale));
                prop_assert!(curve.contains(m.vertices.get(b).point, scale));
            }
        }

        #[test]
        fn prop_cylinder_structural(
            base in prop::array::uniform3(-1e3f64..1e3),
            axis in prop::array::uniform3(-1.0f64..1.0),
            r in 0.5f64..10.0,
            h in 0.1f64..10.0,
        ) {
            let axis = Vector3::from_array(axis);
            prop_assume!(axis.norm() > 0.1); // skip near-zero axes
            let mut m = Model::new();
            m.add_cylinder(Point3::from_array(base), axis, r, h);
            m.rebuild_adjacency();
            prop_assert_eq!(m.vertices.len(), 2);
            prop_assert_eq!(m.edges.len(), 3);
            prop_assert_eq!(m.faces.len(), 3);
            // Every edge used exactly twice, opposite orientation (seam included).
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
        }
    }
}
