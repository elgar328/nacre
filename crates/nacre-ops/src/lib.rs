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
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Orientation, Origin, Shell, Solid, Vertex, VertexDef,
};

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
    /// Pad a boss: imprint `profile` on a planar `face`, then raise the region
    /// outward by `dist` (walls + a top cap). Adds material (design §6, M4).
    PadOnFace {
        face: Handle<Face>,
        profile: Profile2d,
        dist: f64,
    },
    /// Carve a blind pocket: imprint `profile` on a planar `face`, then sink the
    /// region inward by `dist` (walls + a floor). Removes material (design §6, M4).
    PocketOnFace {
        face: Handle<Face>,
        profile: Profile2d,
        dist: f64,
    },
    /// Boolean of two live solids (design §8 M5). M5-c3 implements only
    /// `Common` (intersection) of convex planar solids; other kinds/inputs are
    /// rejected with [`BoolError`].
    Boolean {
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    },
}

/// Which boolean to compute.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolKind {
    /// A ∪ B (not yet implemented — `Unsupported`).
    Fuse,
    /// A − B (not yet implemented — `Unsupported`).
    Cut,
    /// A ∩ B.
    Common,
}

/// Why a boolean could not be computed. The engine rejects out-of-coverage
/// input honestly rather than returning a plausibly-wrong solid (overview
/// 불리언 전략).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolError {
    /// Outside the current coverage: a non-`Common` kind, a non-planar face, a
    /// non-convex input, coplanar faces (across or within an input), a 4-plane
    /// concurrency, a tangential contact, or any other general-position
    /// violation.
    Unsupported,
    /// An input solid handle is not in `model.live_solids`.
    InputNotLive,
    /// The intersection is empty (or too degenerate to be a closed solid).
    EmptyResult,
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
    /// A boolean operation failed (design §8 M5).
    Boolean(BoolError),
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
    /// The superseding solid and the boss's top cap face.
    PadOnFace {
        solid: Handle<Solid>,
        top_face: Handle<Face>,
    },
    /// The superseding solid and the pocket's floor face.
    PocketOnFace {
        solid: Handle<Solid>,
        bottom_face: Handle<Face>,
    },
    /// The boolean result solid (supersedes both inputs).
    Boolean { solid: Handle<Solid> },
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
        Operation::PadOnFace {
            face,
            profile,
            dist,
        } => {
            let (solid, top_face) = pad(model, *face, profile, *dist)?;
            Ok(OpOutput::PadOnFace { solid, top_face })
        }
        Operation::PocketOnFace {
            face,
            profile,
            dist,
        } => {
            let (solid, bottom_face) = pocket(model, *face, profile, *dist)?;
            Ok(OpOutput::PocketOnFace { solid, bottom_face })
        }
        Operation::Boolean { kind, a, b } => {
            let solid = boolean(model, *kind, *a, *b).map_err(OpError::Boolean)?;
            Ok(OpOutput::Boolean { solid })
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

/// The result of splitting a planar face along a profile: the base profile
/// geometry (on the face's plane) plus the rebuilt outer face carrying the
/// profile as a hole. Shared by [`imprint`] and [`pad`]; each finishes by adding
/// its own faces (a coplanar region, or boss walls + a cap) and calling
/// [`finish_split`].
struct Split {
    solid_h: Handle<Solid>,
    shell_h: Handle<Shell>,
    /// The target face's surface and orientation (for a coplanar region face).
    surface_h: Handle<Surface>,
    orientation: Orientation,
    /// Outward normal of the target face.
    n: Vector3,
    base_pv: Vec<Handle<Vertex>>,
    base_pts: Vec<Point3>,
    /// Profile segment edges `base_i → base_{i+1}` (CCW about `n`).
    base_pe: Vec<Handle<Edge>>,
    /// The original face rebuilt with the profile as an inner-loop hole.
    f_outer: Handle<Face>,
}

/// Split a planar `face` along a closed `profile`: push the profile's vertices
/// and edges onto the face's plane (centred on the face, CCW about the outward
/// normal) and build the outer face carrying the profile as a hole. The caller
/// adds the faces that fill/raise the profile region, then calls [`finish_split`].
fn prepare_face_split(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
) -> Result<Split, OpError> {
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
    let base_pts: Vec<Point3> = pts.iter().map(|p| origin + x * p[0] + y * p[1]).collect();

    // Profile vertices and segment edges.
    let m = base_pts.len();
    let base_pv: Vec<Handle<Vertex>> = base_pts
        .iter()
        .map(|p| {
            model.vertices.push(Vertex {
                point: *p,
                origin: Origin::Constructed,
            })
        })
        .collect();
    let base_pe: Vec<Handle<Edge>> = (0..m)
        .map(|i| {
            push_line_edge(
                model,
                base_pv[i],
                base_pts[i],
                base_pv[(i + 1) % m],
                base_pts[(i + 1) % m],
            )
        })
        .collect::<Result<_, _>>()?;

    // Outer face: the original boundary with the profile as a hole — the profile
    // edges reversed + `forward = false` (CW), so each pairs oppositely with the
    // caller's region/wall use.
    let hole_loop = Loop {
        half_edges: base_pe
            .iter()
            .rev()
            .map(|&edge| HalfEdge {
                edge,
                forward: false,
            })
            .collect(),
    };
    let f_outer = model.faces.push(Face {
        surface: surface_h,
        outer: outer_loop,
        inner: vec![hole_loop],
        orientation,
    });

    Ok(Split {
        solid_h,
        shell_h,
        surface_h,
        orientation,
        n,
        base_pv,
        base_pts,
        base_pe,
        f_outer,
    })
}

/// Supersede the split solid (design §2): rebuild its shell with `face` replaced
/// by `new_faces`, push a new solid, and drop the old one from `live_solids`.
fn finish_split(
    model: &mut Model,
    solid_h: Handle<Solid>,
    shell_h: Handle<Shell>,
    face: Handle<Face>,
    new_faces: &[Handle<Face>],
) -> Handle<Solid> {
    let old_faces = model.shells.get(shell_h).faces.clone();
    let mut faces = Vec::with_capacity(old_faces.len() + new_faces.len());
    for &fh in &old_faces {
        if fh == face {
            faces.extend_from_slice(new_faces);
        } else {
            faces.push(fh);
        }
    }
    let cavities = model.solids.get(solid_h).cavities.clone();
    let new_shell = model.shells.push(Shell { faces });
    let new_solid = model.push_solid(Solid {
        outer: new_shell,
        cavities,
    });
    // The old solid is no longer live (its old face lingers as arena).
    model.live_solids.retain(|&s| s != solid_h);
    new_solid
}

/// Imprint a closed `profile` onto a planar `face`: split it into an outer face
/// carrying the profile as an inner-loop hole plus a coplanar region face inside
/// it. Returns `(new solid, region face)`. The model shape is unchanged.
fn imprint(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    let s = prepare_face_split(model, face, profile)?;

    // Region face: profile forward (CCW), same surface/orientation, outward +n.
    let region_loop = Loop {
        half_edges: s
            .base_pe
            .iter()
            .map(|&edge| HalfEdge {
                edge,
                forward: true,
            })
            .collect(),
    };
    let region_face = model.faces.push(Face {
        surface: s.surface_h,
        outer: region_loop,
        inner: vec![],
        orientation: s.orientation,
    });

    let new_solid = finish_split(model, s.solid_h, s.shell_h, face, &[s.f_outer, region_face]);
    Ok((new_solid, region_face))
}

/// Pad a boss on a planar `face`: raise the imprinted region **outward** by
/// `dist` (walls + a cap), adding `profile_area · dist` of material. Returns
/// `(new solid, top cap face)`.
fn pad(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    raise_region(model, face, profile, dist)
}

/// Carve a blind pocket on a planar `face`: raise the imprinted region
/// **inward** by `dist` (walls + a floor), removing `profile_area · dist` of
/// material. Geometrically this is [`pad`] with a negative displacement — the
/// walls face inward and the cap becomes the pocket floor. Returns `(new solid,
/// floor face)`. Precondition (unchecked): `dist` is less than the solid's
/// thickness at the face, so the pocket does not punch through.
fn pocket(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    raise_region(model, face, profile, -dist)
}

/// Imprint `profile` on a planar `face`, then displace its region along the
/// outward normal by `signed_dist` into a prism (walls + a cap) whose open base
/// is the hole in the outer face. `signed_dist > 0` grows a boss (outward);
/// `< 0` carves a pocket (inward — the walls flip to face inward and the cap
/// becomes the floor, since it is built from `base + n·signed_dist`). Callers
/// guarantee `signed_dist != 0`. Returns `(new solid, cap/floor face)`.
fn raise_region(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    signed_dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    // `prepare` mutates the model; callers check `dist > 0` first so an early
    // error leaves live_solids untouched (the op is atomic w.r.t. live).
    let s = prepare_face_split(model, face, profile)?;
    let n = s.base_pts.len();

    // Displaced ring: the base profile moved along the normal by `n · signed_dist`.
    let top_pts: Vec<Point3> = s.base_pts.iter().map(|p| *p + s.n * signed_dist).collect();
    let top_pv: Vec<Handle<Vertex>> = top_pts
        .iter()
        .map(|p| {
            model.vertices.push(Vertex {
                point: *p,
                origin: Origin::Constructed,
            })
        })
        .collect();
    let top_pe: Vec<Handle<Edge>> = (0..n)
        .map(|i| {
            push_line_edge(
                model,
                top_pv[i],
                top_pts[i],
                top_pv[(i + 1) % n],
                top_pts[(i + 1) % n],
            )
        })
        .collect::<Result<_, _>>()?;
    let vert_e: Vec<Handle<Edge>> = (0..n)
        .map(|i| push_line_edge(model, s.base_pv[i], s.base_pts[i], top_pv[i], top_pts[i]))
        .collect::<Result<_, _>>()?;

    // Side walls: quads (base_i, base_j, top_j, top_i) with the extrude side-quad
    // winding. Their normal is `edge × (n·signed_dist)`, so it faces outward for a
    // boss and inward for a pocket automatically. Each base edge pairs oppositely
    // with the outer face's hole; each top edge with the cap; each vertical with a
    // neighbouring wall.
    let mut new_faces = Vec::with_capacity(n + 2);
    new_faces.push(s.f_outer);
    for i in 0..n {
        let j = (i + 1) % n;
        let surface = model.surfaces.push(Surface::Plane(
            Plane::through_points(s.base_pts[i], s.base_pts[j], top_pts[i])
                .ok_or(OpError::DegenerateGeometry)?,
        ));
        let outer = Loop {
            half_edges: vec![
                HalfEdge {
                    edge: s.base_pe[i],
                    forward: true,
                },
                HalfEdge {
                    edge: vert_e[j],
                    forward: true,
                },
                HalfEdge {
                    edge: top_pe[i],
                    forward: false,
                },
                HalfEdge {
                    edge: vert_e[i],
                    forward: false,
                },
            ],
        };
        new_faces.push(model.faces.push(Face {
            surface,
            outer,
            inner: vec![],
            orientation: Orientation::Forward,
        }));
    }

    // Cap: outward normal +n (the boss top, or the pocket floor seen from the
    // opening), the displaced ring forward.
    let cap_surface = model.surfaces.push(Surface::Plane(
        Plane::from_point_normal(top_pts[0], s.n).ok_or(OpError::DegenerateGeometry)?,
    ));
    let cap = model.faces.push(Face {
        surface: cap_surface,
        outer: Loop {
            half_edges: top_pe
                .iter()
                .map(|&edge| HalfEdge {
                    edge,
                    forward: true,
                })
                .collect(),
        },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    new_faces.push(cap);

    let new_solid = finish_split(model, s.solid_h, s.shell_h, face, &new_faces);
    Ok((new_solid, cap))
}

// ---- boolean (M5-c3) ----

use nacre_geom::intersect::{plane_plane, three_plane_orient3d, three_planes};
use std::collections::{HashMap, HashSet};

/// A face's supporting plane plus the data the half-space enumeration needs.
struct PlaneInfo {
    surf: Handle<Surface>,
    plane: Plane,
    /// Three non-collinear outer-CCW loop points; their RH normal is outward.
    tri: [Point3; 3],
    /// Outward normal, `(tri[1]−tri[0])×(tri[2]−tri[0])` normalized — the single
    /// source of "outward" for both the in/out sign test and face ordering.
    n_out: Vector3,
    orient: Orientation,
}

/// A vertex of the result: a point on exactly three planes, inside all others.
struct ResultVertex {
    point: Point3,
    /// Indices (into the combined plane list) of the three defining planes.
    triple: [usize; 3],
    /// Measured accuracy: max distance of `point` to its 3 planes and 3 pairwise
    /// lines (so a `Constructed` edge's `VertexOffCurve` bound holds too).
    tol: f64,
}

/// Boolean of two live solids (design §8 M5, overview 불리언 전략 — 정직하게 거절).
/// **M5-c coverage:** convex, all-planar solids in general position — `Common`
/// (M5-c3, half-space enumeration) and `Fuse`/`Cut` (M5-c4, face clipping,
/// clean-seam). Anything else is rejected with [`BoolError`]. Transactional:
/// each variant computes its result in local structures and pushes only after
/// every degeneracy check passes, so a rejected boolean leaves the model
/// untouched.
pub fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    if !model.live_solids.contains(&a) || !model.live_solids.contains(&b) {
        return Err(BoolError::InputNotLive);
    }
    match kind {
        BoolKind::Common => common(model, a, b),
        BoolKind::Fuse | BoolKind::Cut => fuse_cut(model, kind, a, b),
    }
}

/// `A ∩ B` by half-space vertex enumeration: the intersection is the set of
/// points inside every face half-space of both solids, so each result vertex is
/// the intersection of three of those planes lying inside all the others (an
/// exact indirect-predicate decision — the point is never materialized).
fn common(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    let planes_a = collect_planes(model, a)?;
    let planes_b = collect_planes(model, b)?;
    // Each input must equal the intersection of its face half-spaces (convex).
    if !is_convex(&planes_a, &solid_vertices(model, a))
        || !is_convex(&planes_b, &solid_vertices(model, b))
    {
        return Err(BoolError::Unsupported);
    }
    let mut planes = planes_a;
    planes.extend(planes_b);
    // Coplanar faces (within one input — e.g. an imprinted face's outer+region —
    // or shared across inputs) give the enumeration duplicate half-spaces.
    if has_coplanar_pair(&planes) {
        return Err(BoolError::Unsupported);
    }

    // --- all-local computation (nothing pushed to the model yet) ---
    let verts = enumerate_vertices(&planes)?;
    if verts.len() < 4 {
        return Err(BoolError::EmptyResult); // disjoint or degenerate-empty
    }
    let edges = build_edges(&planes, &verts)?;
    let edge_pairs: HashSet<(usize, usize)> = edges.iter().map(|e| unordered(e.va, e.vb)).collect();
    let faces = build_faces(&planes, &verts, &edge_pairs)?;
    if faces.len() < 4 {
        return Err(BoolError::EmptyResult);
    }

    // --- push (deterministic order → reproducible handles) ---
    Ok(assemble(model, a, b, &planes, &verts, &edges, &faces))
}

/// The supporting planes of a solid's outer shell. `Unsupported` if any face is
/// non-planar or lacks three non-collinear loop points.
fn collect_planes(model: &Model, solid: Handle<Solid>) -> Result<Vec<PlaneInfo>, BoolError> {
    let shell = model.solids.get(solid).outer;
    let mut out = Vec::new();
    for &fh in &model.shells.get(shell).faces {
        let face = model.faces.get(fh);
        let plane = match model.surfaces.get(face.surface) {
            Surface::Plane(p) => *p,
            Surface::Cylinder(_) => return Err(BoolError::Unsupported),
        };
        let tri = outer_tri(model, face).ok_or(BoolError::Unsupported)?;
        let n_out = (tri[1] - tri[0])
            .cross(tri[2] - tri[0])
            .normalize()
            .ok_or(BoolError::Unsupported)?;
        out.push(PlaneInfo {
            surf: face.surface,
            plane,
            tri,
            n_out,
            orient: face.orientation,
        });
    }
    Ok(out)
}

/// The first three non-collinear consecutive start points of a face's outer loop
/// (outer-CCW, so their RH normal is the outward normal).
fn outer_tri(model: &Model, face: &Face) -> Option<[Point3; 3]> {
    let pts: Vec<Point3> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| model.vertices.get(he_start(model, he)).point)
        .collect();
    let n = pts.len();
    (0..n).find_map(|i| {
        let (a, b, c) = (pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        ((b - a).cross(c - a).norm() > 0.0).then_some([a, b, c])
    })
}

/// The distinct outer-shell vertices of a solid.
fn solid_vertices(model: &Model, solid: Handle<Solid>) -> Vec<Point3> {
    let shell = model.solids.get(solid).outer;
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for &fh in &model.shells.get(shell).faces {
        for &he in &model.faces.get(fh).outer.half_edges {
            let vh = he_start(model, he);
            if seen.insert(vh) {
                out.push(model.vertices.get(vh).point);
            }
        }
    }
    out
}

/// Convexity: every vertex is on the inner side of (or on) every face plane.
/// The signed distance along the outward normal `n_out` is `≤ 0` on the inside;
/// a small scale-relative tolerance absorbs the f64 non-coplanarity of a
/// constructed solid's own vertices with its faces (an exact predicate would
/// read that ~1e-13 slop as "outside"). A genuine reflex vertex pokes out by a
/// macroscopic amount, far above the tolerance.
fn is_convex(planes: &[PlaneInfo], verts: &[Point3]) -> bool {
    planes.iter().all(|pi| {
        verts.iter().all(|&v| {
            let sd = (v - pi.tri[0]).dot(pi.n_out);
            let scale = 1e-9 * (v.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
            sd <= scale
        })
    })
}

/// Whether any two planes in the set are coplanar (parallel normals + each
/// origin on the other plane).
fn has_coplanar_pair(planes: &[PlaneInfo]) -> bool {
    (0..planes.len())
        .any(|i| (i + 1..planes.len()).any(|j| coplanar(&planes[i].plane, &planes[j].plane)))
}

fn coplanar(a: &Plane, b: &Plane) -> bool {
    const EPS: f64 = 1e-9;
    a.normal().cross(b.normal()).norm() <= EPS
        && a.distance(b.origin()) <= EPS
        && b.distance(a.origin()) <= EPS
}

/// Enumerate result vertices: every plane triple whose intersection point lies
/// inside all other half-spaces. `three_plane_orient3d == −1` means the implicit
/// point is on the inner side of a face (derived in the plan); `+1` outside
/// (reject the triple); `0` on a 4th plane (a concurrency degeneracy →
/// `Unsupported`).
fn enumerate_vertices(planes: &[PlaneInfo]) -> Result<Vec<ResultVertex>, BoolError> {
    let n = planes.len();
    let mut verts = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let (pi, pj, pk) = (&planes[i].plane, &planes[j].plane, &planes[k].plane);
                let Some(point) = three_planes(pi, pj, pk) else {
                    continue; // parallel/near-coplanar triple — no vertex
                };
                let mut outside = false;
                let mut on_extra = false;
                for (m, pm) in planes.iter().enumerate() {
                    if m == i || m == j || m == k {
                        continue;
                    }
                    match three_plane_orient3d(pi, pj, pk, pm.tri[0], pm.tri[1], pm.tri[2]) {
                        1 => {
                            outside = true;
                            break;
                        }
                        0 => on_extra = true,
                        _ => {} // −1: inside this half-space
                    }
                }
                if outside {
                    continue;
                }
                if on_extra {
                    return Err(BoolError::Unsupported); // 4-plane concurrency
                }
                verts.push(ResultVertex {
                    point,
                    triple: [i, j, k],
                    tol: vertex_tol(point, pi, pj, pk),
                });
            }
        }
    }
    Ok(verts)
}

/// Max distance of `p` to its 3 planes and 3 pairwise lines (the measured
/// `Origin::Discovered` tolerance).
fn vertex_tol(p: Point3, a: &Plane, b: &Plane, c: &Plane) -> f64 {
    let mut tol = a.distance(p).max(b.distance(p)).max(c.distance(p));
    for (x, y) in [(a, b), (a, c), (b, c)] {
        if let Some(line) = plane_plane(x, y) {
            tol = tol.max(line.distance(p));
        }
    }
    tol
}

/// A result edge: two vertices sharing two planes, on the line of that plane pair.
struct ResultEdge {
    va: usize,
    vb: usize,
    planes: [usize; 2],
}

/// One edge per plane pair whose line carries exactly two result vertices.
fn build_edges(planes: &[PlaneInfo], verts: &[ResultVertex]) -> Result<Vec<ResultEdge>, BoolError> {
    let n = planes.len();
    let mut edges = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let on_pair: Vec<usize> = verts
                .iter()
                .enumerate()
                .filter(|(_, v)| v.triple.contains(&i) && v.triple.contains(&j))
                .map(|(vi, _)| vi)
                .collect();
            match on_pair.len() {
                0 => {}
                2 => edges.push(ResultEdge {
                    va: on_pair[0],
                    vb: on_pair[1],
                    planes: [i, j],
                }),
                _ => return Err(BoolError::Unsupported), // tangent / degenerate edge
            }
        }
    }
    Ok(edges)
}

/// One face per plane carrying ≥3 result vertices, ordered CCW about the outward
/// normal, with each consecutive pair confirmed to be a result edge.
fn build_faces(
    planes: &[PlaneInfo],
    verts: &[ResultVertex],
    edge_pairs: &HashSet<(usize, usize)>,
) -> Result<Vec<(usize, Vec<usize>)>, BoolError> {
    let mut faces = Vec::new();
    for (m, pm) in planes.iter().enumerate() {
        let on_m: Vec<usize> = verts
            .iter()
            .enumerate()
            .filter(|(_, v)| v.triple.contains(&m))
            .map(|(vi, _)| vi)
            .collect();
        if on_m.len() < 3 {
            continue; // this plane is not a face of the result
        }
        let ordered = order_ccw(&on_m, verts, pm.n_out)?;
        // Every consecutive pair must be a result edge (else non-manifold).
        let k = ordered.len();
        for t in 0..k {
            if !edge_pairs.contains(&unordered(ordered[t], ordered[(t + 1) % k])) {
                return Err(BoolError::Unsupported);
            }
        }
        faces.push((m, ordered));
    }
    Ok(faces)
}

/// Order coplanar points CCW about the outward normal `n_out` (angle about the
/// centroid). `Unsupported` if they are collinear (a tangential contact).
fn order_ccw(
    idxs: &[usize],
    verts: &[ResultVertex],
    n_out: Vector3,
) -> Result<Vec<usize>, BoolError> {
    let pts: Vec<Point3> = idxs.iter().map(|&vi| verts[vi].point).collect();
    let c = Point3::centroid(&pts).ok_or(BoolError::Unsupported)?;
    let u = (pts[0] - c).normalize().ok_or(BoolError::Unsupported)?;
    let w = n_out.cross(u);
    // Collinear ⇒ every point lies on the `u` axis (no `w` spread) ⇒ tangent.
    let spread = idxs
        .iter()
        .map(|&vi| (verts[vi].point - c).dot(w).abs())
        .fold(0.0, f64::max);
    let scale = idxs
        .iter()
        .map(|&vi| (verts[vi].point - c).norm())
        .fold(0.0, f64::max);
    if spread <= 1e-9 * scale.max(1.0) {
        return Err(BoolError::Unsupported);
    }
    let mut keyed: Vec<(f64, usize)> = idxs
        .iter()
        .map(|&vi| {
            let d = verts[vi].point - c;
            (d.dot(w).atan2(d.dot(u)), vi)
        })
        .collect();
    keyed.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite angles"));
    Ok(keyed.iter().map(|&(_, vi)| vi).collect())
}

fn unordered(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Push the computed result and supersede the inputs (mirrors `finish_split`).
fn assemble(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    planes: &[PlaneInfo],
    verts: &[ResultVertex],
    edges: &[ResultEdge],
    faces: &[(usize, Vec<usize>)],
) -> Handle<Solid> {
    // Vertices (discovered three-plane points).
    let vh: Vec<Handle<Vertex>> = verts
        .iter()
        .map(|v| {
            let def = VertexDef::ThreePlane([
                planes[v.triple[0]].surf,
                planes[v.triple[1]].surf,
                planes[v.triple[2]].surf,
            ]);
            model.vertices.push(Vertex {
                point: v.point,
                origin: Origin::Discovered {
                    tol: v.tol,
                    definition: def,
                },
            })
        })
        .collect();
    // Edges (constructed; lines are the plane-pair intersections).
    let mut edge_of: HashMap<(usize, usize), Handle<Edge>> = HashMap::new();
    for e in edges {
        let line = plane_plane(&planes[e.planes[0]].plane, &planes[e.planes[1]].plane)
            .expect("survivor edge planes meet in a line");
        let curve = model.curves.push(Curve::Line(line));
        let eh = model.edges.push(Edge {
            curve,
            bounds: Some([vh[e.va], vh[e.vb]]),
            origin: Origin::Constructed,
        });
        edge_of.insert(
            unordered(vh[e.va].index() as usize, vh[e.vb].index() as usize),
            eh,
        );
    }
    // Faces (outward-CCW loops of the plane's vertices).
    let mut face_handles = Vec::new();
    for (m, loop_verts) in faces {
        let k = loop_verts.len();
        let half_edges = (0..k)
            .map(|t| {
                let va = vh[loop_verts[t]];
                let vb = vh[loop_verts[(t + 1) % k]];
                let eh = edge_of[&unordered(va.index() as usize, vb.index() as usize)];
                let forward = model.edges.get(eh).bounds.expect("bounded")[0] == va;
                HalfEdge { edge: eh, forward }
            })
            .collect();
        face_handles.push(model.faces.push(Face {
            surface: planes[*m].surf,
            outer: Loop { half_edges },
            inner: vec![],
            orientation: planes[*m].orient,
        }));
    }
    let shell = model.shells.push(Shell {
        faces: face_handles,
    });
    let solid = model.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    model.live_solids.retain(|&s| s != a && s != b);
    solid
}

// ---- fuse / cut (M5-c4): face clipping, clean-seam convex ----

/// A vertex's side relative to the *other* solid.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Side {
    Inside,
    Outside,
}

/// A seam vertex — a three-plane point on both `∂A` and `∂B` (2 A-planes + 1
/// B-plane, or 1 A + 2 B). Shared (one `Handle`) by every incident result piece.
struct SeamVertex {
    point: Point3,
    triple: [usize; 3], // sorted combined-plane indices
    tol: f64,
}

/// A node in a reconstructed face loop. `Eq`/`Hash` give identity dedup so an
/// A-piece and a B-piece that meet at a seam node share one result vertex/edge.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Node {
    Orig(Handle<Vertex>),
    Seam([usize; 3]), // sorted triple (key into the seam map)
}

/// A reconstructed result face: which combined plane it is on, its loop as
/// nodes, and whether to flip it (cut's inside-A B-pieces).
struct LocalFace {
    plane_idx: usize,
    loop_nodes: Vec<Node>,
    flip: bool,
}

/// `Fuse` (A∪B) / `Cut` (A−B) of two convex solids by face clipping (M5-c4).
/// Clean-seam general position only; else [`BoolError::Unsupported`].
fn fuse_cut(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    if kind == BoolKind::Cut {
        return Err(BoolError::Unsupported); // M5-c4 commit 2
    }
    let planes_a = collect_planes(model, a)?;
    let planes_b = collect_planes(model, b)?;
    if !is_convex(&planes_a, &solid_vertices(model, a))
        || !is_convex(&planes_b, &solid_vertices(model, b))
    {
        return Err(BoolError::Unsupported);
    }
    let na = planes_a.len();
    let mut planes = planes_a;
    planes.extend(planes_b);
    if has_coplanar_pair(&planes) {
        return Err(BoolError::Unsupported);
    }
    let a_range = 0..na;
    let b_range = na..planes.len();

    let mut surf_ix: HashMap<Handle<Surface>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        surf_ix.insert(pi.surf, i);
    }

    // Classify every original vertex vs the other solid (direct signed-distance).
    let mut classof: HashMap<Handle<Vertex>, Side> = HashMap::new();
    for &vh in &solid_vertex_handles(model, a) {
        let side = classify_vertex(model.vertices.get(vh).point, &planes[b_range.clone()])?;
        classof.insert(vh, side);
    }
    for &vh in &solid_vertex_handles(model, b) {
        let side = classify_vertex(model.vertices.get(vh).point, &planes[a_range.clone()])?;
        classof.insert(vh, side);
    }

    // Seam vertices: A-edges piercing B (2A+1B) and B-edges piercing A (1A+2B).
    let edges_a = edge_incidence(model, a, &surf_ix);
    let edges_b = edge_incidence(model, b, &surf_ix);
    let mut seam: Vec<SeamVertex> = Vec::new();
    let mut seam_ix: HashMap<[usize; 3], usize> = HashMap::new();
    let mut edge_seam: HashMap<Handle<Edge>, [usize; 3]> = HashMap::new();
    for (edges, other) in [(&edges_a, b_range.clone()), (&edges_b, a_range.clone())] {
        for &(eh, bounds, ref inc) in edges {
            let [v0, v1] = bounds;
            let (s0, s1) = (classof[&v0], classof[&v1]);
            if s0 == s1 {
                continue; // edge does not cross the other solid
            }
            let (p_out, p_in) = if s0 == Side::Outside {
                (model.vertices.get(v0).point, model.vertices.get(v1).point)
            } else {
                (model.vertices.get(v1).point, model.vertices.get(v0).point)
            };
            let entry =
                enter_face(p_out, p_in, other.clone(), &planes).ok_or(BoolError::Unsupported)?;
            let [e0, e1] = [inc[0], inc[1]];
            let point = three_planes(&planes[e0].plane, &planes[e1].plane, &planes[entry].plane)
                .ok_or(BoolError::Unsupported)?;
            // Exact: the point must be inside the entry face (on ∂ of the other solid).
            for m in other.clone() {
                if m == entry {
                    continue;
                }
                match three_plane_orient3d(
                    &planes[e0].plane,
                    &planes[e1].plane,
                    &planes[entry].plane,
                    planes[m].tri[0],
                    planes[m].tri[1],
                    planes[m].tri[2],
                ) {
                    -1 => {}
                    _ => return Err(BoolError::Unsupported), // outside or on a 4th plane
                }
            }
            let mut triple = [e0, e1, entry];
            triple.sort_unstable();
            edge_seam.insert(eh, triple);
            if let std::collections::hash_map::Entry::Vacant(slot) = seam_ix.entry(triple) {
                let tol = vertex_tol(
                    point,
                    &planes[e0].plane,
                    &planes[e1].plane,
                    &planes[entry].plane,
                );
                slot.insert(seam.len());
                seam.push(SeamVertex { point, triple, tol });
            }
        }
    }

    // Which side of each solid each op keeps.
    let (keep_a, keep_b, flip_b) = match kind {
        BoolKind::Fuse => (Side::Outside, Side::Outside, false),
        BoolKind::Cut => (Side::Outside, Side::Inside, true),
        BoolKind::Common => unreachable!("common has its own path"),
    };
    if seam.is_empty() {
        // No crossing: either disjoint (empty) or containment (needs a cavity).
        return Err(disjoint_or_contained(model, a, b, &planes, na));
    }

    // Reconstruct faces of A (keep_a side) and B (keep_b side, flip_b).
    let mut faces: Vec<LocalFace> = Vec::new();
    for (solid, keep, flip) in [(a, keep_a, false), (b, keep_b, flip_b)] {
        let shell = model.solids.get(solid).outer;
        for &fh in &model.shells.get(shell).faces {
            let face = model.faces.get(fh);
            let pidx = surf_ix[&face.surface];
            if let Some(lf) = reconstruct_face(
                model, face, pidx, keep, flip, &classof, &edge_seam, &seam, &seam_ix,
            )? {
                faces.push(lf);
            }
        }
    }
    if faces.len() < 4 {
        return Err(BoolError::EmptyResult);
    }

    Ok(assemble_fuse_cut(model, a, b, &planes, &seam, &faces))
}

/// Inside/outside/on a convex solid by signed distance along each face's outward
/// normal. On any face plane (`|sd| ≤ scale`) ⇒ not general position ⇒
/// `Unsupported`.
fn classify_vertex(v: Point3, other: &[PlaneInfo]) -> Result<Side, BoolError> {
    let scale = 1e-9 * (v.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
    let mut side = Side::Inside;
    for pm in other {
        let sd = (v - pm.tri[0]).dot(pm.n_out);
        if sd.abs() <= scale {
            return Err(BoolError::Unsupported); // on the other solid's boundary
        }
        if sd > scale {
            side = Side::Outside; // outside this half-space ⇒ outside the convex
        }
    }
    Ok(side)
}

/// The face of the convex solid (`range`) that segment `p_out → p_in` enters
/// through: the last inward crossing (Cyrus–Beck argmax of `t`).
fn enter_face(
    p_out: Point3,
    p_in: Point3,
    range: std::ops::Range<usize>,
    planes: &[PlaneInfo],
) -> Option<usize> {
    let mut best: Option<(f64, usize)> = None;
    for m in range {
        let a = (p_out - planes[m].tri[0]).dot(planes[m].n_out);
        let b = (p_in - planes[m].tri[0]).dot(planes[m].n_out);
        if a > 0.0 && b < 0.0 {
            let t = a / (a - b);
            if best.is_none_or(|(bt, _)| t > bt) {
                best = Some((t, m));
            }
        }
    }
    best.map(|(_, m)| m)
}

/// Distinct outer-shell vertex handles of a solid, in shell→face→loop order.
fn solid_vertex_handles(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Vertex>> {
    let shell = model.solids.get(solid).outer;
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for &fh in &model.shells.get(shell).faces {
        for &he in &model.faces.get(fh).outer.half_edges {
            let vh = he_start(model, he);
            if seen.insert(vh) {
                out.push(vh);
            }
        }
    }
    out
}

/// Each outer-shell edge with its bound vertices and the two combined-plane
/// indices of its adjacent faces, in first-seen (deterministic) order.
#[allow(clippy::type_complexity)]
fn edge_incidence(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Surface>, usize>,
) -> Vec<(Handle<Edge>, [Handle<Vertex>; 2], Vec<usize>)> {
    let shell = model.solids.get(solid).outer;
    let mut order: Vec<Handle<Edge>> = Vec::new();
    let mut map: HashMap<Handle<Edge>, ([Handle<Vertex>; 2], Vec<usize>)> = HashMap::new();
    for &fh in &model.shells.get(shell).faces {
        let face = model.faces.get(fh);
        let pidx = surf_ix[&face.surface];
        for he in &face.outer.half_edges {
            let bounds = model.edges.get(he.edge).bounds.expect("bounded");
            let entry = map.entry(he.edge).or_insert_with(|| {
                order.push(he.edge);
                (bounds, Vec::new())
            });
            entry.1.push(pidx);
        }
    }
    order
        .into_iter()
        .map(|e| {
            let (b, p) = map.remove(&e).unwrap();
            (e, b, p)
        })
        .collect()
}

/// Classify a no-seam overlap: disjoint ⇒ `EmptyResult`, otherwise containment
/// ⇒ `Unsupported` (would need a cavity shell). Overlap ⇔ any vertex of one
/// solid is inside the other.
fn disjoint_or_contained(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    planes: &[PlaneInfo],
    na: usize,
) -> BoolError {
    let a_inside = solid_vertex_handles(model, a).iter().any(|&vh| {
        classify_vertex(model.vertices.get(vh).point, &planes[na..]) == Ok(Side::Inside)
    });
    let b_inside = solid_vertex_handles(model, b).iter().any(|&vh| {
        classify_vertex(model.vertices.get(vh).point, &planes[..na]) == Ok(Side::Inside)
    });
    if a_inside || b_inside {
        BoolError::Unsupported // containment: cavity needed
    } else {
        BoolError::EmptyResult // disjoint
    }
}

/// Reconstruct a face's kept portion. `None` if the whole face is dropped.
#[allow(clippy::too_many_arguments)]
fn reconstruct_face(
    model: &Model,
    face: &Face,
    plane_idx: usize,
    keep: Side,
    flip: bool,
    classof: &HashMap<Handle<Vertex>, Side>,
    edge_seam: &HashMap<Handle<Edge>, [usize; 3]>,
    seam: &[SeamVertex],
    seam_ix: &HashMap<[usize; 3], usize>,
) -> Result<Option<LocalFace>, BoolError> {
    let hes = &face.outer.half_edges;
    let n = hes.len();
    let verts: Vec<Handle<Vertex>> = hes.iter().map(|&he| he_start(model, he)).collect();
    let kept: Vec<bool> = verts.iter().map(|v| classof[v] == keep).collect();

    // Count boundary crossings (kept↔dropped edge transitions).
    let transitions: Vec<usize> = (0..n).filter(|&i| kept[i] != kept[(i + 1) % n]).collect();
    match transitions.len() {
        0 => {
            if kept[0] {
                // Whole face kept.
                let loop_nodes = verts.iter().map(|&v| Node::Orig(v)).collect();
                return Ok(Some(LocalFace {
                    plane_idx,
                    loop_nodes,
                    flip,
                }));
            }
            // Whole face dropped — but a seam loop interior to it would be a hole.
            if seam.iter().any(|s| s.triple.contains(&plane_idx)) {
                return Err(BoolError::Unsupported); // interior seam ⇒ hole
            }
            Ok(None)
        }
        2 => {
            // One kept run; splice the seam sub-path across the dropped run.
            // Transition edge i (kept[i]!=kept[i+1]) carries seam vertex edge_seam[edge].
            let s_at = |i: usize| -> Result<[usize; 3], BoolError> {
                edge_seam
                    .get(&hes[i].edge)
                    .copied()
                    .ok_or(BoolError::Unsupported)
            };
            // Boundary crossing on each transition edge.
            let (t0, t1) = (transitions[0], transitions[1]);
            let (b0, b1) = (s_at(t0)?, s_at(t1)?);
            // Kept run: vertices with kept==true, starting right after a drop→keep edge.
            // Identify the keep→drop edge (kept[t]) and drop→keep edge.
            let (kd, dk) = if kept[t0] { (t0, t1) } else { (t1, t0) };
            let (s_kd, s_dk) = if kept[t0] { (b0, b1) } else { (b1, b0) };
            // Walk kept run from dk+1 .. kd (inclusive), CCW.
            let mut nodes: Vec<Node> = Vec::new();
            let mut i = (dk + 1) % n;
            loop {
                nodes.push(Node::Orig(verts[i]));
                if i == kd {
                    break;
                }
                i = (i + 1) % n;
            }
            // Seam sub-path S_kd → bends → S_dk (bends = seam on this plane, not the
            // two boundary crossings), ordered along S_kd→S_dk.
            let p_kd = seam[seam_ix[&s_kd]].point;
            let p_dk = seam[seam_ix[&s_dk]].point;
            let dir = p_dk - p_kd;
            let mut bends: Vec<[usize; 3]> = seam
                .iter()
                .filter(|s| s.triple.contains(&plane_idx) && s.triple != s_kd && s.triple != s_dk)
                .map(|s| s.triple)
                .collect();
            bends.sort_by(|x, y| {
                let px = seam[seam_ix[x]].point;
                let py = seam[seam_ix[y]].point;
                (px - p_kd)
                    .dot(dir)
                    .partial_cmp(&(py - p_kd).dot(dir))
                    .expect("finite")
            });
            nodes.push(Node::Seam(s_kd));
            nodes.extend(bends.into_iter().map(Node::Seam));
            nodes.push(Node::Seam(s_dk));
            Ok(Some(LocalFace {
                plane_idx,
                loop_nodes: nodes,
                flip,
            }))
        }
        _ => Err(BoolError::Unsupported), // ≥4 crossings ⇒ multiple chords
    }
}

/// Push the reconstructed result and supersede the inputs (mirrors `assemble`).
fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    planes: &[PlaneInfo],
    seam: &[SeamVertex],
    faces: &[LocalFace],
) -> Handle<Solid> {
    // Vertices (deterministic: first appearance across faces in order).
    let mut vh: HashMap<Node, Handle<Vertex>> = HashMap::new();
    let mut node_handle = |model: &mut Model, node: Node| -> Handle<Vertex> {
        if let Some(&h) = vh.get(&node) {
            return h;
        }
        let handle = match node {
            Node::Orig(orig) => {
                let point = model.vertices.get(orig).point;
                model.vertices.push(Vertex {
                    point,
                    origin: Origin::Constructed,
                })
            }
            Node::Seam(triple) => {
                let sv = seam.iter().find(|s| s.triple == triple).expect("seam node");
                let def = VertexDef::ThreePlane([
                    planes[triple[0]].surf,
                    planes[triple[1]].surf,
                    planes[triple[2]].surf,
                ]);
                model.vertices.push(Vertex {
                    point: sv.point,
                    origin: Origin::Discovered {
                        tol: sv.tol,
                        definition: def,
                    },
                })
            }
        };
        vh.insert(node, handle);
        handle
    };
    // Materialize all vertex handles first (deterministic order).
    for lf in faces {
        for &node in &lf.loop_nodes {
            node_handle(model, node);
        }
    }

    // Edges keyed by unordered handle-index pair (lookup only).
    let mut edge_of: HashMap<(usize, usize), Handle<Edge>> = HashMap::new();
    let mut edge_for =
        |model: &mut Model, va: Handle<Vertex>, vb: Handle<Vertex>| -> Handle<Edge> {
            let key = unordered(va.index() as usize, vb.index() as usize);
            if let Some(&e) = edge_of.get(&key) {
                return e;
            }
            let pa = model.vertices.get(va).point;
            let pb = model.vertices.get(vb).point;
            let curve = model
                .curves
                .push(Curve::Line(Line::through_points(pa, pb).expect("distinct")));
            let e = model.edges.push(Edge {
                curve,
                bounds: Some([va, vb]),
                origin: Origin::Constructed,
            });
            edge_of.insert(key, e);
            e
        };

    let mut face_handles = Vec::new();
    for lf in faces {
        let handles: Vec<Handle<Vertex>> = lf.loop_nodes.iter().map(|&nd| vh[&nd]).collect();
        let k = handles.len();
        let mut half_edges: Vec<HalfEdge> = (0..k)
            .map(|t| {
                let (va, vb) = (handles[t], handles[(t + 1) % k]);
                let e = edge_for(model, va, vb);
                let forward = model.edges.get(e).bounds.expect("bounded")[0] == va;
                HalfEdge { edge: e, forward }
            })
            .collect();
        let orientation = if lf.flip {
            // Cut's inside-A B-pieces: reverse the loop and toggle orientation so
            // the outward normal points into the removed region.
            half_edges.reverse();
            for he in &mut half_edges {
                he.forward = !he.forward;
            }
            match planes[lf.plane_idx].orient {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
            }
        } else {
            planes[lf.plane_idx].orient
        };
        face_handles.push(model.faces.push(Face {
            surface: planes[lf.plane_idx].surf,
            outer: Loop { half_edges },
            inner: vec![],
            orientation,
        }));
    }
    let shell = model.shells.push(Shell {
        faces: face_handles,
    });
    let solid = model.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    model.live_solids.retain(|&s| s != a && s != b);
    solid
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

    fn pad_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
        Operation::PadOnFace {
            face,
            profile,
            dist,
        }
    }

    #[test]
    fn pad_boss_on_cube_top() {
        let (mut m, top) = cube_with_top();
        let out = apply(&mut m, &pad_op(top, small_square(), 0.5)).unwrap();
        let OpOutput::PadOnFace { top_face, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();

        let v = nacre_validate::validate(&m);
        assert!(v.is_empty(), "{v:?}");

        let reach = m.reachable();
        // 6 cube faces − top + (outer' + 4 walls + cap) = 11.
        assert_eq!(reach.faces.len(), 11);
        assert_eq!(reach.vertices.len(), 16); // 8 cube + 4 base + 4 top
        assert_eq!(reach.edges.len(), 24); // 12 cube + 4 base + 4 top + 4 vertical
        let inner: usize = reach
            .faces
            .iter()
            .map(|fh| m.faces.get(*fh).inner.len())
            .sum();
        assert_eq!(inner, 1);
        assert!(reach.faces.contains(&top_face));
    }

    #[test]
    fn pad_rejects_nonplanar_face() {
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
            apply(&mut m, &pad_op(lateral, small_square(), 0.5)),
            Err(OpError::NonPlanarFace)
        ));
    }

    #[test]
    fn pad_rejects_nonpositive_dist() {
        let (mut m, top) = cube_with_top();
        assert!(matches!(
            apply(&mut m, &pad_op(top, small_square(), 0.0)),
            Err(OpError::NonPositiveDistance)
        ));
    }

    #[test]
    fn pad_rejects_degenerate_profile() {
        let (mut m, top) = cube_with_top();
        let two = Profile2d {
            points: vec![p2(0.0, 0.0), p2(0.1, 0.0)],
        };
        assert!(matches!(
            apply(&mut m, &pad_op(top, two, 0.5)),
            Err(OpError::DegenerateProfile)
        ));
    }

    #[test]
    fn pad_step_exports() {
        // The boss (holed outer face + walls + cap) exports without error.
        let (mut m, top) = cube_with_top();
        apply(&mut m, &pad_op(top, small_square(), 0.5)).unwrap();
        let step = nacre_step::to_step(&m).expect("boss exports");
        assert!(step.contains("FACE_BOUND("), "the hole emits a FACE_BOUND");
    }

    fn pocket_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
        Operation::PocketOnFace {
            face,
            profile,
            dist,
        }
    }

    #[test]
    fn pocket_on_cube_top() {
        let (mut m, top) = cube_with_top();
        let out = apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap();
        let OpOutput::PocketOnFace { bottom_face, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();

        let v = nacre_validate::validate(&m);
        assert!(v.is_empty(), "{v:?}");

        // Same topology as a boss (a downward prism instead of upward).
        let reach = m.reachable();
        assert_eq!(reach.faces.len(), 11);
        assert_eq!(reach.vertices.len(), 16);
        assert_eq!(reach.edges.len(), 24);
        let inner: usize = reach
            .faces
            .iter()
            .map(|fh| m.faces.get(*fh).inner.len())
            .sum();
        assert_eq!(inner, 1);
        assert!(reach.faces.contains(&bottom_face));
    }

    #[test]
    fn pocket_rejects_nonplanar_face() {
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
            apply(&mut m, &pocket_op(lateral, small_square(), 0.5)),
            Err(OpError::NonPlanarFace)
        ));
    }

    #[test]
    fn pocket_rejects_nonpositive_dist() {
        let (mut m, top) = cube_with_top();
        assert!(matches!(
            apply(&mut m, &pocket_op(top, small_square(), 0.0)),
            Err(OpError::NonPositiveDistance)
        ));
    }

    #[test]
    fn pocket_rejects_degenerate_profile() {
        let (mut m, top) = cube_with_top();
        let two = Profile2d {
            points: vec![p2(0.0, 0.0), p2(0.1, 0.0)],
        };
        assert!(matches!(
            apply(&mut m, &pocket_op(top, two, 0.5)),
            Err(OpError::DegenerateProfile)
        ));
    }

    #[test]
    fn pocket_step_exports() {
        let (mut m, top) = cube_with_top();
        apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap();
        let step = nacre_step::to_step(&m).expect("pocket exports");
        assert!(step.contains("FACE_BOUND("), "the hole emits a FACE_BOUND");
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

        /// A random boss on a random box stays a valid b-rep (any interior
        /// profile, any positive height).
        #[test]
        fn prop_pad_stays_valid(
            sx in 0.5f64..5.0,
            sy in 0.5f64..5.0,
            sz in 0.5f64..5.0,
            h in 0.05f64..0.15,
            dist in 0.1f64..5.0,
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
            let hole = Profile2d {
                points: vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)],
            };
            apply(&mut m, &Operation::PadOnFace { face: faces[1], profile: hole, dist }).unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
        }

        /// A random blind pocket on a random box stays valid. `dist ≤ 0.8 < sz`
        /// keeps the pocket from punching through the box (height `sz ≥ 1`).
        #[test]
        fn prop_pocket_stays_valid(
            sx in 0.5f64..5.0,
            sy in 0.5f64..5.0,
            sz in 1.0f64..5.0,
            h in 0.05f64..0.15,
            dist in 0.1f64..0.8,
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
            let hole = Profile2d {
                points: vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)],
            };
            apply(&mut m, &Operation::PocketOnFace { face: faces[1], profile: hole, dist }).unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
        }
    }

    // ---- boolean API (M5-c3 commit 1) ----

    fn two_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        (m, a, b)
    }

    #[test]
    fn boolean_cut_not_yet_supported() {
        // Cut lands in M5-c4 commit 2; for now it is Unsupported.
        let (mut m, a, b) = two_boxes();
        assert_eq!(
            boolean(&mut m, BoolKind::Cut, a, b),
            Err(BoolError::Unsupported)
        );
    }

    #[test]
    fn fuse_of_two_cubes() {
        // A = [0,1]³, B = [0.5,1.5]³ ⇒ A∪B volume 1+1−0.125 = 1.875.
        let (mut m, a, b) = two_boxes();
        let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.875).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn fuse_common_inclusion_exclusion() {
        // vol(A∪B) + vol(A∩B) == vol(A) + vol(B). boolean supersedes inputs, so
        // fuse and common run on independent copies.
        let corner = |o: f64| {
            (
                Point3::from_array([o; 3]),
                Point3::from_array([o + 1.0, o + 1.0, o + 1.0]),
            )
        };
        let (amin, amax) = corner(0.0);
        let (bmin, bmax) = corner(0.5);
        let vol_of = |mn, mx| {
            let mut m = Model::new();
            let s = m.add_cuboid(mn, mx);
            nacre_props::mass_props(&m, s).unwrap().volume
        };
        let (va, vb) = (vol_of(amin, amax), vol_of(bmin, bmax));

        let mut m1 = Model::new();
        let (a1, b1) = (m1.add_cuboid(amin, amax), m1.add_cuboid(bmin, bmax));
        let rf = boolean(&mut m1, BoolKind::Fuse, a1, b1).unwrap();
        let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;

        let mut m2 = Model::new();
        let (a2, b2) = (m2.add_cuboid(amin, amax), m2.add_cuboid(bmin, bmax));
        let rc = boolean(&mut m2, BoolKind::Common, a2, b2).unwrap();
        let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;

        assert!((vf + vc - va - vb).abs() < 1e-9, "{vf}+{vc} vs {va}+{vb}");
    }

    #[test]
    fn boolean_rejects_non_live_input() {
        let (mut m, a, b) = two_boxes();
        m.live_solids.retain(|&s| s != b); // as if superseded
        assert_eq!(
            boolean(&mut m, BoolKind::Common, a, b),
            Err(BoolError::InputNotLive)
        );
    }

    #[test]
    fn boolean_op_applies_and_wraps_error() {
        // A degenerate boolean's error is surfaced as OpError::Boolean.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        // Offset in all axes so no faces are coplanar with A (else the coplanar
        // gate fires first).
        let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
        // Disjoint ⇒ EmptyResult, wrapped.
        assert_eq!(
            apply(
                &mut m,
                &Operation::Boolean {
                    kind: BoolKind::Common,
                    a,
                    b
                }
            ),
            Err(OpError::Boolean(BoolError::EmptyResult))
        );
    }

    // ---- boolean Common algorithm (M5-c3 commit 2) ----

    #[test]
    fn common_of_two_cubes_is_their_overlap() {
        // A = [0,1]³, B = [0.5,1.5]³ ⇒ A∩B = [0.5,1]³ (volume 0.125). This also
        // pins the in/out sign convention end to end: a flipped parity would take
        // the complement and give a wrong (or non-closed) result.
        let (mut m, a, b) = two_boxes();
        let r = boolean(&mut m, BoolKind::Common, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let reach = m.reachable();
        assert_eq!(reach.vertices.len(), 8);
        assert_eq!(reach.edges.len(), 12);
        assert_eq!(reach.faces.len(), 6);
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.125).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]); // A and B superseded
    }

    #[test]
    fn common_with_enclosing_box_is_the_inner_solid() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let c = m.add_cuboid(Point3::from_array([-5.0; 3]), Point3::from_array([5.0; 3]));
        let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
        let r = boolean(&mut m, BoolKind::Common, a, c).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let reach = m.reachable();
        assert_eq!(reach.vertices.len(), 8);
        assert_eq!(reach.faces.len(), 6);
        let vol_r = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol_r - vol_a).abs() < 1e-9, "{vol_r} vs {vol_a}");
    }

    #[test]
    fn common_of_disjoint_boxes_is_empty() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        // Offset in all axes so no faces are coplanar with A.
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        assert_eq!(
            boolean(&mut m, BoolKind::Common, a, b),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn common_rejects_non_planar_input() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let cyl = m.add_cylinder(
            Point3::from_array([1.0, 1.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        assert_eq!(
            boolean(&mut m, BoolKind::Common, a, cyl),
            Err(BoolError::Unsupported)
        );
    }

    #[test]
    fn common_rejects_non_convex_input() {
        // An L-shaped prism (reflex edge) is not the intersection of its face
        // half-spaces.
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
        let mut m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let lsolid = *m.live_solids.first().unwrap();
        let b = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.5, 1.5, 0.5]),
        );
        assert_eq!(
            boolean(&mut m, BoolKind::Common, lsolid, b),
            Err(BoolError::Unsupported)
        );
    }

    #[test]
    fn common_rejects_coplanar_faces_from_imprint() {
        // Imprinting splits a face into two coplanar faces (outer + region), so
        // the imprinted-but-still-convex cube has coplanar half-spaces.
        let (mut m, top) = cube_with_top();
        let hole = Profile2d {
            points: vec![p2(-0.2, -0.2), p2(0.2, -0.2), p2(0.2, 0.2), p2(-0.2, 0.2)],
        };
        let OpOutput::ImprintSketch { solid, .. } = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: top,
                profile: hole,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let b = m.add_cuboid(
            Point3::from_array([0.2, 0.2, 0.2]),
            Point3::from_array([1.2, 1.2, 1.2]),
        );
        assert_eq!(
            boolean(&mut m, BoolKind::Common, solid, b),
            Err(BoolError::Unsupported)
        );
    }

    proptest! {
        /// Overlapping axis-aligned boxes: the intersection volume equals the
        /// independent AABB-overlap product (mixed A/B axis-aligned vertices).
        #[test]
        fn common_axis_boxes_volume_matches_aabb_overlap(
            amin in prop::array::uniform3(-5.0f64..5.0),
            aext in prop::array::uniform3(1.0f64..4.0),
            t in prop::array::uniform3(0.05f64..0.7),
            bext in prop::array::uniform3(1.0f64..4.0),
        ) {
            let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
            let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
            let bmax: [f64; 3] = std::array::from_fn(|i| bmin[i] + bext[i]);
            let expected: f64 = (0..3)
                .map(|i| (amax[i].min(bmax[i]) - bmin[i]).max(0.0))
                .product();
            prop_assume!(expected > 1e-3);

            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
            let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
            let res = boolean(&mut m, BoolKind::Common, a, b);
            prop_assume!(res.is_ok()); // skip rare coplanar/degenerate configs
            let r = res.unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            prop_assert!((vol - expected).abs() <= 1e-9 * expected.max(1.0), "{vol} vs {expected}");
        }

        /// A tilted square prism (oblique planes) intersected with a big enclosing
        /// box is the prism — exercises non-axis-aligned face normals (the in/out
        /// sign and CCW ordering) with an independent oracle (the prism's own mass).
        #[test]
        fn common_tilted_prism_with_enclosing_box_is_the_prism(
            nx in -0.5f64..0.5,
            ny in -0.5f64..0.5,
        ) {
            let plane = SketchPlane::from_origin_normal(
                Point3::origin(),
                Vector3::from_array([nx, ny, 1.0]),
            )
            .unwrap();
            let mut m = replay(&[Operation::Extrude {
                plane,
                profile: square(),
                dist: 1.0,
            }])
            .unwrap();
            let prism = *m.live_solids.first().unwrap();
            let vol_prism = nacre_props::mass_props(&m, prism).unwrap().volume;
            let c = m.add_cuboid(Point3::from_array([-10.0; 3]), Point3::from_array([10.0; 3]));
            let res = boolean(&mut m, BoolKind::Common, prism, c);
            prop_assume!(res.is_ok());
            let r = res.unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            let vol_r = nacre_props::mass_props(&m, r).unwrap().volume;
            prop_assert!(
                (vol_r - vol_prism).abs() <= 1e-9 * vol_prism.max(1.0),
                "{vol_r} vs {vol_prism}"
            );
        }
    }
}
