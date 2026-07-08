//! Operations for the nacre kernel, plus a replayable operation log (design §6).
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::ImprintSketch`] (M4) splits an existing planar face along a
//! closed profile — the first op that consumes a prior op's face by `Handle`
//! (exposed via [`OpOutput`]) and supersedes a solid (design §2 live-solid
//! semantics). Ops are applied by [`apply`] and folded by [`replay`]; every
//! result is a **closed** solid, so `nacre-validate` applies fully.

use nacre_geom::{Curve, Line, Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, HalfEdge, Loop, Model, Orientation, Origin, Shell, Solid, Vertex};

/// A sketch-plane frame: a 2-D point `(u, v)` maps to `origin + u·x + v·y`.
/// `x_axis`/`y_axis` are assumed unit and orthogonal (the constructors ensure
/// it); the plane normal is `x × y`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SketchPlane {
    pub origin: Point3,
    pub x_axis: Vector3,
    pub y_axis: Vector3,
}

impl SketchPlane {
    /// The world XY plane (normal +Z).
    pub fn world_xy() -> Self {
        Self {
            origin: Point3::origin(),
            x_axis: Vector3::from_array([1.0, 0.0, 0.0]),
            y_axis: Vector3::from_array([0.0, 1.0, 0.0]),
        }
    }

    /// A plane through `origin` with the given `normal`; `x`/`y` axes are
    /// synthesized (`x = n.any_perpendicular()`, `y = n × x`). `None` if `normal`
    /// is zero.
    pub fn from_origin_normal(origin: Point3, normal: Vector3) -> Option<Self> {
        let n = normal.normalize()?;
        let x = n.any_perpendicular()?;
        Some(Self {
            origin,
            x_axis: x,
            y_axis: n.cross(x),
        })
    }

    /// The 3-D point for sketch coordinates `p = (u, v)`.
    #[inline]
    pub fn point(&self, p: Point2) -> Point3 {
        self.origin + self.x_axis * p[0] + self.y_axis * p[1]
    }

    /// The plane normal `x × y` (unit when the axes are unit and orthogonal).
    #[inline]
    pub fn normal(&self) -> Vector3 {
        self.x_axis.cross(self.y_axis)
    }
}

/// A closed, simple polygon profile (straight segments only), at least 3 points.
/// Extrusion normalizes the winding to counter-clockwise, so input orientation
/// does not matter; self-intersection is assumed absent (not checked in M2).
#[derive(Clone, Debug, PartialEq)]
pub struct Profile2d {
    pub points: Vec<Point2>,
}

/// A modelling operation.
#[derive(Clone, Debug, PartialEq)]
pub enum Operation {
    /// Extrude `profile` (on `plane`) by `dist` along the plane normal.
    Extrude {
        plane: SketchPlane,
        profile: Profile2d,
        dist: f64,
    },
    /// Imprint a closed `profile` onto an existing planar `face`, splitting it
    /// into an outer face with the profile as a hole plus a coplanar region
    /// face. The profile is expressed in a frame derived from the face (M4).
    ImprintSketch {
        face: Handle<Face>,
        profile: Profile2d,
    },
}

/// A failure while applying an operation.
#[derive(Debug, PartialEq)]
pub enum OpError {
    /// A profile with fewer than 3 points.
    DegenerateProfile,
    /// A non-positive extrusion distance.
    NonPositiveDistance,
    /// A curve/surface construction collapsed (collinear/coincident points, a
    /// zero-length profile edge).
    DegenerateGeometry,
    /// An imprint target face is not planar (only planar faces support imprint
    /// in M4; curved-face imprint arrives with the quadric milestones).
    NonPlanarFace,
    /// An imprint target face belongs to no live solid's outer shell (a stale
    /// or non-live handle).
    FaceNotInLiveSolid,
}

/// The handles an operation produced. Not `Copy`: `Extrude` carries a `Vec`.
#[derive(Clone, Debug, PartialEq)]
pub enum OpOutput {
    /// The created solid and its faces in push order: `faces[0]` base cap,
    /// `faces[1]` top cap, then one side face per profile edge.
    Extrude {
        solid: Handle<Solid>,
        faces: Vec<Handle<Face>>,
    },
    /// The superseding solid and the new coplanar region face (inside the hole).
    ImprintSketch {
        solid: Handle<Solid>,
        region_face: Handle<Face>,
    },
}

/// Apply one operation to `model`, returning the handles it created. Does not
/// rebuild the adjacency cache (do that once after a batch — see [`replay`]).
pub fn apply(model: &mut Model, op: &Operation) -> Result<OpOutput, OpError> {
    match op {
        Operation::Extrude {
            plane,
            profile,
            dist,
        } => {
            let (solid, faces) = extrude(model, plane, profile, *dist)?;
            Ok(OpOutput::Extrude { solid, faces })
        }
        Operation::ImprintSketch { face, profile } => {
            let (solid, region_face) = imprint(model, *face, profile)?;
            Ok(OpOutput::ImprintSketch { solid, region_face })
        }
    }
}

/// Replay an operation log into a fresh model. Deterministic: the same log
/// reproduces the same model down to handle indices.
pub fn replay(ops: &[Operation]) -> Result<Model, OpError> {
    let mut model = Model::new();
    for op in ops {
        apply(&mut model, op)?;
    }
    model.rebuild_adjacency();
    Ok(model)
}

/// Twice the signed area of the polygon (sign only is used): positive = CCW.
fn signed_area(pts: &[Point2]) -> f64 {
    let n = pts.len();
    (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
        })
        .sum()
}

fn push_line_edge(
    model: &mut Model,
    a: Handle<Vertex>,
    ap: Point3,
    b: Handle<Vertex>,
    bp: Point3,
) -> Result<Handle<Edge>, OpError> {
    let curve = model.curves.push(Curve::Line(
        Line::through_points(ap, bp).ok_or(OpError::DegenerateGeometry)?,
    ));
    Ok(model.edges.push(Edge {
        curve,
        bounds: Some([a, b]),
        origin: Origin::Constructed,
    }))
}

/// Build a prism: the profile forms the base and (translated by `normal·dist`)
/// the top; each profile edge grows a side quad. Winding generalizes the M1
/// cuboid — base loop reversed (normal −N, outward), top forward (+N), side
/// `(B_i, B_{i+1}, T_{i+1}, T_i)`; every edge is used twice with opposite flags.
fn extrude(
    model: &mut Model,
    plane: &SketchPlane,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    if profile.points.len() < 3 {
        return Err(OpError::DegenerateProfile);
    }

    // Normalize to CCW so the synthesized normals point outward.
    let mut pts = profile.points.clone();
    if signed_area(&pts) < 0.0 {
        pts.reverse();
    }
    let n = pts.len();
    let normal = plane.normal();
    let base_pts: Vec<Point3> = pts.iter().map(|p| plane.point(*p)).collect();
    let top_pts: Vec<Point3> = base_pts.iter().map(|b| *b + normal * dist).collect();

    let bv: Vec<Handle<Vertex>> = base_pts
        .iter()
        .map(|p| {
            model.vertices.push(Vertex {
                point: *p,
                origin: Origin::Constructed,
            })
        })
        .collect();
    let tv: Vec<Handle<Vertex>> = top_pts
        .iter()
        .map(|p| {
            model.vertices.push(Vertex {
                point: *p,
                origin: Origin::Constructed,
            })
        })
        .collect();

    let mut be = Vec::with_capacity(n); // base edges B_i -> B_{i+1}
    let mut te = Vec::with_capacity(n); // top edges  T_i -> T_{i+1}
    let mut ve = Vec::with_capacity(n); // vertical   B_i -> T_i
    for i in 0..n {
        let j = (i + 1) % n;
        be.push(push_line_edge(
            model,
            bv[i],
            base_pts[i],
            bv[j],
            base_pts[j],
        )?);
        te.push(push_line_edge(model, tv[i], top_pts[i], tv[j], top_pts[j])?);
        ve.push(push_line_edge(
            model,
            bv[i],
            base_pts[i],
            tv[i],
            top_pts[i],
        )?);
    }

    let mut faces = Vec::with_capacity(n + 2);

    // Base cap: outward normal −N, loop reversed (B_0 -> B_{n-1} -> ... -> B_1).
    let base_surface = model.surfaces.push(Surface::Plane(
        Plane::from_point_normal(base_pts[0], -normal).ok_or(OpError::DegenerateGeometry)?,
    ));
    let base_loop = Loop {
        half_edges: (0..n)
            .rev()
            .map(|i| HalfEdge {
                edge: be[i],
                forward: false,
            })
            .collect(),
    };
    faces.push(model.faces.push(Face {
        surface: base_surface,
        outer: base_loop,
        inner: vec![],
        orientation: Orientation::Forward,
    }));

    // Top cap: outward normal +N.
    let top_surface = model.surfaces.push(Surface::Plane(
        Plane::from_point_normal(top_pts[0], normal).ok_or(OpError::DegenerateGeometry)?,
    ));
    let top_loop = Loop {
        half_edges: (0..n)
            .map(|i| HalfEdge {
                edge: te[i],
                forward: true,
            })
            .collect(),
    };
    faces.push(model.faces.push(Face {
        surface: top_surface,
        outer: top_loop,
        inner: vec![],
        orientation: Orientation::Forward,
    }));

    // Side quads.
    for i in 0..n {
        let j = (i + 1) % n;
        let surface = model.surfaces.push(Surface::Plane(
            Plane::through_points(base_pts[i], base_pts[j], top_pts[i])
                .ok_or(OpError::DegenerateGeometry)?,
        ));
        let outer = Loop {
            half_edges: vec![
                HalfEdge {
                    edge: be[i],
                    forward: true,
                },
                HalfEdge {
                    edge: ve[j],
                    forward: true,
                },
                HalfEdge {
                    edge: te[i],
                    forward: false,
                },
                HalfEdge {
                    edge: ve[i],
                    forward: false,
                },
            ],
        };
        faces.push(model.faces.push(Face {
            surface,
            outer,
            inner: vec![],
            orientation: Orientation::Forward,
        }));
    }

    let shell = model.shells.push(Shell {
        faces: faces.clone(),
    });
    let solid = model.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    Ok((solid, faces))
}

/// The start vertex of a half-edge (`bounds[0]` if forward, else `bounds[1]`).
/// The target face of an imprint is part of a valid solid, so its edges are
/// bounded.
fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
    let [a, b] = model
        .edges
        .get(he.edge)
        .bounds
        .expect("a solid's loop edge is bounded");
    if he.forward { a } else { b }
}

/// Imprint a closed `profile` onto a planar `face`: split it into an outer face
/// carrying the profile as an inner-loop hole plus a coplanar region face inside
/// it. Supersedes the owning solid (design §2) — reuses the untouched cells,
/// pushes a new shell/solid, and swaps `live_solids`. Returns `(new solid,
/// region face)`. The model shape is unchanged (a coplanar subdivision).
fn imprint(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    if profile.points.len() < 3 {
        return Err(OpError::DegenerateProfile);
    }

    // Locate the live solid whose outer shell holds this face.
    let (solid_h, shell_h) = model
        .live_solids
        .iter()
        .map(|&s| (s, model.solids.get(s).outer))
        .find(|&(_, sh)| model.shells.get(sh).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;

    // Read the target face, then release the borrow before mutating.
    let f = model.faces.get(face);
    let surface_h = f.surface;
    let orientation = f.orientation;
    let outer_loop = f.outer.clone();
    let plane = match model.surfaces.get(surface_h) {
        Surface::Plane(p) => *p,
        Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    };

    // Outward normal and an in-plane right-handed frame (x × y = n), centred on
    // the face so the profile's (0, 0) lands at the face centre.
    let sign = match orientation {
        Orientation::Forward => 1.0,
        Orientation::Reversed => -1.0,
    };
    let n = plane.normal() * sign;
    let x = n.any_perpendicular().ok_or(OpError::DegenerateGeometry)?;
    let y = n.cross(x);
    let outer_pts: Vec<Point3> = outer_loop
        .half_edges
        .iter()
        .map(|he| model.vertices.get(he_start(model, *he)).point)
        .collect();
    let origin = Point3::centroid(&outer_pts).ok_or(OpError::DegenerateGeometry)?;

    // Profile → CCW in the frame (positive signed area) so its RH normal is +n.
    let mut pts = profile.points.clone();
    if signed_area(&pts) < 0.0 {
        pts.reverse();
    }
    let pts3: Vec<Point3> = pts.iter().map(|p| origin + x * p[0] + y * p[1]).collect();

    // Profile vertices and segment edges.
    let m = pts3.len();
    let pv: Vec<Handle<Vertex>> = pts3
        .iter()
        .map(|p| {
            model.vertices.push(Vertex {
                point: *p,
                origin: Origin::Constructed,
            })
        })
        .collect();
    let pe: Vec<Handle<Edge>> = (0..m)
        .map(|i| push_line_edge(model, pv[i], pts3[i], pv[(i + 1) % m], pts3[(i + 1) % m]))
        .collect::<Result<_, _>>()?;

    // Region face: profile forward (CCW), outward +n.
    let region_loop = Loop {
        half_edges: pe
            .iter()
            .map(|&edge| HalfEdge {
                edge,
                forward: true,
            })
            .collect(),
    };
    let region_face = model.faces.push(Face {
        surface: surface_h,
        outer: region_loop,
        inner: vec![],
        orientation,
    });

    // Outer face: the original boundary with the profile as a hole — same edges
    // reversed + `forward = false` (CW), so each profile edge pairs oppositely.
    let hole_loop = Loop {
        half_edges: pe
            .iter()
            .rev()
            .map(|&edge| HalfEdge {
                edge,
                forward: false,
            })
            .collect(),
    };
    let outer_face = model.faces.push(Face {
        surface: surface_h,
        outer: outer_loop,
        inner: vec![hole_loop],
        orientation,
    });

    // New shell: the old faces with `face` replaced by the two new faces.
    let old_faces = model.shells.get(shell_h).faces.clone();
    let mut new_faces = Vec::with_capacity(old_faces.len() + 1);
    for &fh in &old_faces {
        if fh == face {
            new_faces.push(outer_face);
            new_faces.push(region_face);
        } else {
            new_faces.push(fh);
        }
    }
    let cavities = model.solids.get(solid_h).cavities.clone();
    let new_shell = model.shells.push(Shell { faces: new_faces });
    let new_solid = model.push_solid(Solid {
        outer: new_shell,
        cavities,
    });

    // Supersede: the old solid is no longer live (its old face lingers as arena).
    model.live_solids.retain(|&s| s != solid_h);

    Ok((new_solid, region_face))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn p2(x: f64, y: f64) -> Point2 {
        Point2::from_array([x, y])
    }

    fn square() -> Profile2d {
        Profile2d {
            points: vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)],
        }
    }

    fn extrude_op(profile: Profile2d, dist: f64) -> Operation {
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist,
        }
    }

    fn regular_ngon(n: usize, r: f64) -> Profile2d {
        let points = (0..n)
            .map(|i| {
                let a = std::f64::consts::TAU * (i as f64) / (n as f64);
                p2(r * a.cos(), r * a.sin())
            })
            .collect();
        Profile2d { points }
    }

    #[test]
    fn square_extrudes_to_a_cube() {
        let m = replay(&[extrude_op(square(), 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.edges.len(), 12);
        assert_eq!(m.faces.len(), 6);
        assert_eq!(m.solids.len(), 1);

        let mut got: Vec<[f64; 3]> = m.vertices.iter().map(|(_, v)| v.point.as_array()).collect();
        let mut want = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let key = |p: &[f64; 3]| (p[0].to_bits(), p[1].to_bits(), p[2].to_bits());
        got.sort_by_key(key);
        want.sort_by_key(key);
        assert_eq!(got, want);
    }

    #[test]
    fn triangle_extrudes_to_a_prism() {
        let tri = Profile2d {
            points: vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(1.0, 1.5)],
        };
        let m = replay(&[extrude_op(tri, 3.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 6);
        assert_eq!(m.edges.len(), 9);
        assert_eq!(m.faces.len(), 5);
    }

    #[test]
    fn pentagon_extrudes_clean() {
        let m = replay(&[extrude_op(regular_ngon(5, 2.0), 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 10);
        assert_eq!(m.faces.len(), 7);
    }

    #[test]
    fn concave_l_profile_is_valid() {
        // An L-shape (a reflex vertex) — a simple concave hexagon.
        let l = Profile2d {
            points: vec![
                p2(0.0, 0.0),
                p2(2.0, 0.0),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
            ],
        };
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 12);
        assert_eq!(m.faces.len(), 8);
    }

    #[test]
    fn clockwise_input_is_auto_corrected() {
        // The square wound CW; auto-CCW makes it a valid cube anyway.
        let cw = Profile2d {
            points: vec![p2(0.0, 1.0), p2(1.0, 1.0), p2(1.0, 0.0), p2(0.0, 0.0)],
        };
        let m = replay(&[extrude_op(cw, 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.faces.len(), 6);
    }

    #[test]
    fn replay_is_deterministic() {
        let log = vec![extrude_op(square(), 1.0)];
        let m1 = replay(&log).unwrap();
        let m2 = replay(&log).unwrap();
        let pts = |m: &Model| {
            m.vertices
                .iter()
                .map(|(_, v)| v.point.as_array())
                .collect::<Vec<_>>()
        };
        assert_eq!(pts(&m1), pts(&m2));
        assert_eq!(m1.edges.len(), m2.edges.len());
        assert_eq!(m1.faces.len(), m2.faces.len());
    }

    #[test]
    fn two_extrudes_make_two_solids() {
        let far = SketchPlane {
            origin: Point3::from_array([5.0, 0.0, 0.0]),
            ..SketchPlane::world_xy()
        };
        let log = vec![
            extrude_op(square(), 1.0),
            Operation::Extrude {
                plane: far,
                profile: square(),
                dist: 1.0,
            },
        ];
        let m = replay(&log).unwrap();
        assert_eq!(m.solids.len(), 2);
        assert!(nacre_validate::validate(&m).is_empty());
    }

    #[test]
    fn degenerate_inputs_are_rejected() {
        let plane = SketchPlane::world_xy();
        let two = Profile2d {
            points: vec![p2(0.0, 0.0), p2(1.0, 0.0)],
        };
        assert_eq!(
            apply(
                &mut Model::new(),
                &Operation::Extrude {
                    plane,
                    profile: two,
                    dist: 1.0
                }
            ),
            Err(OpError::DegenerateProfile)
        );
        assert_eq!(
            apply(&mut Model::new(), &extrude_op(square(), 0.0)),
            Err(OpError::NonPositiveDistance)
        );
        let dup = Profile2d {
            points: vec![p2(0.0, 0.0), p2(0.0, 0.0), p2(1.0, 1.0)],
        };
        assert_eq!(
            apply(
                &mut Model::new(),
                &Operation::Extrude {
                    plane,
                    profile: dup,
                    dist: 1.0
                }
            ),
            Err(OpError::DegenerateGeometry)
        );
    }

    /// Extrude a unit cube and return `(model, top face handle)`.
    fn cube_with_top() -> (Model, Handle<Face>) {
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude_op(square(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        let top = faces[1]; // base, top, sides…
        (m, top)
    }

    fn small_square() -> Profile2d {
        Profile2d {
            points: vec![p2(-0.2, -0.2), p2(0.2, -0.2), p2(0.2, 0.2), p2(-0.2, 0.2)],
        }
    }

    #[test]
    fn imprint_square_hole_in_cube_top() {
        let (mut m, top) = cube_with_top();
        let out = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: top,
                profile: small_square(),
            },
        )
        .unwrap();
        let OpOutput::ImprintSketch { region_face, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();

        let v = nacre_validate::validate(&m);
        assert!(v.is_empty(), "{v:?}");

        let reach = m.reachable();
        assert_eq!(reach.faces.len(), 7); // 6 − top + (outer' + region)
        let inner: usize = reach
            .faces
            .iter()
            .map(|fh| m.faces.get(*fh).inner.len())
            .sum();
        assert_eq!(inner, 1); // one hole, in the outer face
        assert!(reach.faces.contains(&region_face));
        // Euler: V 12, E 16, F 7, L_i 1 → χ = 2.
        assert_eq!(reach.vertices.len(), 12);
        assert_eq!(reach.edges.len(), 16);
    }

    #[test]
    fn imprint_rejects_nonplanar_face() {
        let mut m = Model::new();
        m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        let shell = m.solids.get(m.live_solids[0]).outer;
        let lateral = *m
            .shells
            .get(shell)
            .faces
            .iter()
            .find(|&&fh| {
                matches!(
                    m.surfaces.get(m.faces.get(fh).surface),
                    Surface::Cylinder(_)
                )
            })
            .unwrap();
        assert!(matches!(
            apply(
                &mut m,
                &Operation::ImprintSketch {
                    face: lateral,
                    profile: small_square(),
                },
            ),
            Err(OpError::NonPlanarFace)
        ));
    }

    #[test]
    fn imprint_rejects_degenerate_profile() {
        let (mut m, top) = cube_with_top();
        let two = Profile2d {
            points: vec![p2(0.0, 0.0), p2(0.1, 0.0)],
        };
        assert!(matches!(
            apply(
                &mut m,
                &Operation::ImprintSketch {
                    face: top,
                    profile: two,
                },
            ),
            Err(OpError::DegenerateProfile)
        ));
    }

    #[test]
    fn imprint_step_roundtrips() {
        let (mut m, top) = cube_with_top();
        apply(
            &mut m,
            &Operation::ImprintSketch {
                face: top,
                profile: small_square(),
            },
        )
        .unwrap();
        // The imprinted solid exports (nacre-step handles the inner loop), and
        // the hole is emitted as a FACE_BOUND (distinct from FACE_OUTER_BOUND).
        let step = nacre_step::to_step(&m).expect("imprinted solid exports");
        assert!(
            step.contains("FACE_BOUND("),
            "hole should emit a FACE_BOUND"
        );
    }

    proptest! {
        #[test]
        fn prop_regular_ngon_on_xy_is_clean(
            n in 3usize..8,
            r in 0.5f64..10.0,
            dist in 0.1f64..10.0,
        ) {
            let m = replay(&[extrude_op(regular_ngon(n, r), dist)]).unwrap();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            prop_assert_eq!(m.vertices.len(), 2 * n);
            prop_assert_eq!(m.faces.len(), n + 2);
        }

        #[test]
        fn prop_ngon_on_arbitrary_plane_is_clean(
            n in 3usize..8,
            nx in -1.0f64..1.0,
            ny in -1.0f64..1.0,
            nz in -1.0f64..1.0,
            dist in 0.1f64..10.0,
        ) {
            let normal = Vector3::from_array([nx, ny, nz]);
            prop_assume!(normal.norm() > 0.1); // skip near-zero normals
            let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
            let m = replay(&[Operation::Extrude {
                plane,
                profile: regular_ngon(n, 2.0),
                dist,
            }])
            .unwrap();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            prop_assert_eq!(m.faces.len(), n + 2);
        }

        /// A random box, then a small centred square imprinted on its top face,
        /// stays a valid b-rep. The hole half-size `h` keeps its circumradius
        /// `h√2 < 0.15·√2 ≈ 0.21` below the top face's inradius `min(sx,sy)/2 ≥
        /// 0.25`, so the profile is interior regardless of the derived frame.
        #[test]
        fn prop_imprint_stays_valid(
            sx in 0.5f64..5.0,
            sy in 0.5f64..5.0,
            sz in 0.5f64..5.0,
            h in 0.05f64..0.15,
        ) {
            let rect = Profile2d {
                points: vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)],
            };
            let mut m = Model::new();
            let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: rect,
                dist: sz,
            }).unwrap() else { unreachable!() };
            let top = faces[1];
            let hole = Profile2d {
                points: vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)],
            };
            apply(&mut m, &Operation::ImprintSketch { face: top, profile: hole }).unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
        }
    }
}
