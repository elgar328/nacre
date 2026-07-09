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

mod arrange;

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

/// Names of the `Unsupported` reject sites, shared by the guard that raises one
/// and the test that asserts it — a renamed tag then cannot silently drift out
/// of a test's expectation. Not `#[cfg(test)]`: the guards name these in release
/// builds too. They are `const`, so they inline away where the tag is unused.
pub(crate) mod tag {
    pub const HOLLOW_OPERAND: &str = "hollow_operand";
    pub const INNER_LOOP_OPERAND: &str = "inner_loop_operand";
    /// A backstop with no firing test yet — the degeneracies that would produce an
    /// odd crossing count are expected to trip `CONTACT_DEGENERATE` first. Recorded
    /// as unverified in design.md §9, alongside `fourplane`.
    pub const ARRANGEMENT_DEGENERATE: &str = "arrangement_degenerate";
    pub const COMMON_OVERLAP: &str = "common_overlap";
    pub const PIERCED_MULTI: &str = "pierced_multi";
    pub const COPLANAR_PAIR: &str = "coplanar_pair";
    pub const TUNNEL: &str = "tunnel";
    pub const NO_ENTRY_FACE: &str = "no_entry_face";
    pub const THREE_PLANES: &str = "three_planes";
    pub const FOURPLANE: &str = "fourplane";
    pub const NONCONVEX_OPERAND: &str = "nonconvex_operand";
    pub const CYLINDER_FACE: &str = "cylinder_face";
    pub const DEGENERATE_FACE: &str = "degenerate_face";
    pub const DEGENERATE_NORMAL: &str = "degenerate_normal";
    pub const TANGENT_EDGE: &str = "tangent_edge";
    pub const NONMANIFOLD_FACE: &str = "nonmanifold_face";
    pub const DEGENERATE_CENTROID: &str = "degenerate_centroid";
    pub const DEGENERATE_RADIUS: &str = "degenerate_radius";
    pub const COLLINEAR_FACE: &str = "collinear_face";
    pub const POKE_THROUGH: &str = "poke_through";
    pub const OUTSIDE_OR_FOURPLANE: &str = "outside_or_fourplane";
    pub const ON_BOUNDARY: &str = "on_boundary";
    pub const RAY_DEGENERATE: &str = "ray_degenerate";
    pub const CONTACT_DEGENERATE: &str = "contact_degenerate";
    pub const POKEHOLE: &str = "pokehole";
    pub const MISSING_SEAM: &str = "missing_seam";
    pub const STRICTARC: &str = "strictarc";
    pub const MULTICHORD: &str = "multichord";
}

#[cfg(test)]
thread_local! {
    static LAST_REJECT: std::cell::Cell<Option<&'static str>> =
        const { std::cell::Cell::new(None) };
}

/// Build an `Unsupported`, recording *which* guard raised it. `Err(Unsupported)`
/// alone cannot distinguish the guards, so a reject test whose fixture drifts
/// onto a different guard would still pass — [`assert_rejects`] closes that hole.
/// Every `Unsupported` site in this crate goes through here: two call sites
/// discard a `collect_planes` error (`detect_coincident_interface`), so only
/// exhaustive tagging makes "last tag written == the site that returned" hold.
#[inline]
#[cfg_attr(not(test), allow(unused_variables))]
pub(crate) fn reject(tag: &'static str) -> BoolError {
    #[cfg(test)]
    LAST_REJECT.with(|c| c.set(Some(tag)));
    BoolError::Unsupported
}

/// Assert that `f` rejects *through the intended guard*. Clears any stale tag
/// first, so a prior call in the same test cannot be mistaken for this one.
#[cfg(test)]
fn assert_rejects<T: std::fmt::Debug + PartialEq>(
    f: impl FnOnce() -> Result<T, BoolError>,
    expect: &'static str,
) {
    LAST_REJECT.with(|c| c.take());
    assert_eq!(f(), Err(BoolError::Unsupported));
    assert_eq!(LAST_REJECT.with(|c| c.take()), Some(expect));
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
pub(crate) fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
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

use nacre_geom::intersect::{
    RayCross, SegCross, plane_plane, ray_face_cross, segment_face_cross, three_plane_orient3d,
    three_planes,
};
use std::collections::{HashMap, HashSet};

/// A face's supporting plane plus the data the half-space enumeration needs.
///
/// `three_plane_orient3d(.., tri[0], tri[1], tri[2])` returns `+1` when the
/// implicit point lies on **`tri`'s right-hand-normal side** — the convention is
/// tied to the triangle, never to `plane`. `n_out` happens to equal that RH normal
/// only because `tri` is taken outer-CCW; `plane.normal()` is the *surface's*
/// normal and may point inward on a `Reversed` face. Every sign test here reads
/// `n_out` (or `tri`), and none reads `plane.normal()`.
pub(crate) struct PlaneInfo {
    pub(crate) surf: Handle<Surface>,
    pub(crate) plane: Plane,
    /// Three non-collinear outer-CCW loop points; their RH normal is outward.
    pub(crate) tri: [Point3; 3],
    /// Outward normal, `(tri[1]−tri[0])×(tri[2]−tri[0])` normalized — the single
    /// source of "outward" for both the in/out sign test and face ordering.
    pub(crate) n_out: Vector3,
    pub(crate) orient: Orientation,
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
    // A cavitied operand is out of coverage, and every path below gets it wrong
    // *silently*: `collect_planes`/`solid_vertices` see the outer shell alone, so a
    // hollow box reads as convex; `assemble_fuse_cut` then emits `cavities: vec![]`
    // and the void vanishes with `validate` still clean. Reject rather than answer
    // wrong. Whether to support cavities or keep rejecting is a sub-unit 5 decision.
    // Inputs only — a *result* may be hollow (containment `Cut` builds one).
    if !model.solids.get(a).cavities.is_empty() || !model.solids.get(b).cavities.is_empty() {
        return Err(reject(tag::HOLLOW_OPERAND));
    }
    // Coincident-coplanar degeneracy (M5-c5): a clean matched-interface stack.
    if let Some(iface) = detect_coincident_interface(model, a, b) {
        return coincident_merge(model, kind, a, b, &iface);
    }
    // Convexity selector (M5-d1): the convex paths (`common`/`fuse_cut`) treat each
    // operand as the intersection of its face half-spaces. A non-convex operand
    // instead takes the general seam-free path (containment/disjoint); a
    // non-convex *overlap* (a real seam) is out of this sub-unit's coverage.
    let planes_a = collect_planes(model, a)?;
    let planes_b = collect_planes(model, b)?;
    if is_convex(&planes_a, &solid_vertices(model, a))
        && is_convex(&planes_b, &solid_vertices(model, b))
    {
        match kind {
            BoolKind::Common => common(model, a, b),
            BoolKind::Fuse | BoolKind::Cut => fuse_cut(model, kind, a, b),
        }
    } else {
        nonconvex_seamfree(model, kind, a, b)
    }
}

/// The boolean of two solids, at least one non-convex, when their boundaries do
/// not cross (M5-d1). A real seam ⇒ `Unsupported` (the general arrangement is a
/// later sub-unit); otherwise one solid contains the other or they are disjoint,
/// classified by exact [`point_in_solid`] and assembled by [`contained_result`].
fn nonconvex_seamfree(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    if boundaries_intersect(model, a, b)? {
        // A genuine seam. Single-chord `Fuse`/`Cut` (M5-d2) is handled here; the
        // non-convex `Common`, and multi-chord/multi-loop seams, are later sub-units.
        return match kind {
            BoolKind::Fuse | BoolKind::Cut => overlap_fuse_cut(model, kind, a, b),
            BoolKind::Common => Err(reject(tag::COMMON_OVERLAP)),
        };
    }
    // Seam-free: each solid's vertices all fall on one side of the other.
    let mut classof: HashMap<Handle<Vertex>, Side> = HashMap::new();
    for (verts_solid, other) in [(a, b), (b, a)] {
        for vh in solid_vertex_handles(model, verts_solid) {
            let side = point_in_solid(model, model.vertices.get(vh).point, other)?;
            classof.insert(vh, side);
        }
    }
    debug_assert!(
        {
            let one = |s: Handle<Solid>| {
                let hs = solid_vertex_handles(model, s);
                let first = hs.first().and_then(|vh| classof.get(vh)).copied();
                hs.iter().all(|vh| classof.get(vh).copied() == first)
            };
            one(a) && one(b)
        },
        "seam-free classification must be consistent per solid"
    );
    contained_result(model, kind, a, b, &classof)
}

/// The single outer-shell face of `face_solid` that segment `p0 → p1` pierces,
/// returned as a combined-plane index (`surf_ix`). The non-convex analogue of
/// [`enter_face`]: each face is fanned and tested by [`segment_crosses_face`]
/// (oriented winding, apex retry). `Ok(None)` if no face is pierced, `Ok(Some)`
/// for exactly one, `Unsupported` if the segment pierces more than one face (it
/// threads multiple chords — beyond this sub-unit) or any contact is `Degenerate`.
fn pierced_face(
    model: &Model,
    p0: Point3,
    p1: Point3,
    face_solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Surface>, usize>,
) -> Result<Option<usize>, BoolError> {
    let shell = model.solids.get(face_solid).outer;
    let mut hit: Option<usize> = None;
    for &fh in &model.shells.get(shell).faces {
        let rings = face_loops(model, fh);
        if segment_crosses_face(p0, p1, &rings)? {
            if hit.is_some() {
                return Err(reject(tag::PIERCED_MULTI)); // pierces >1 face — multi-chord edge
            }
            hit = Some(surf_ix[&model.faces.get(fh).surface]);
        }
    }
    Ok(hit)
}

/// `Fuse`/`Cut` of two solids, at least one non-convex, whose boundaries cross in
/// a single chord per face (M5-d2). Mirrors [`fuse_cut`] but (a) classifies each
/// original vertex with exact [`point_in_solid`] instead of the convex
/// half-space [`classify_vertex`], (b) finds seam entries with [`pierced_face`]
/// instead of the convex Cyrus–Beck [`enter_face`]/[`segment_enters`], and (c)
/// reconstructs faces in `strict` mode (rejecting a non-convex/self-intersecting
/// seam arc). Multiple chords, poke-through holes, and non-convex `Common` are
/// honestly `Unsupported` (later sub-units).
/// Reject an operand carrying a face with holes.
///
/// The seam machinery and `solid_local_faces` walk outer rings only, so a hole-ring
/// edge is seen once — its other use is on the lid's inner loop. Three things then
/// go wrong, and all three were measured on a pocketed cube:
///
/// * `reconstruct_face` rebuilds the lid from its outer ring and `assemble_fuse_cut`
///   emits `inner: vec![]`, so the hole is dropped: `Ok`, volume 0.96996 instead of
///   0.916625, `validate` reporting five violations.
/// * `edge_incidence` gives such an edge one incident plane, and `overlap_fuse_cut`
///   indexes `inc[1]` — a panic, whenever a rim edge straddles the seam.
/// * `coincident_merge` drops the hole through `solid_local_faces`: `Ok`, volume
///   2.0533 instead of 2.0.
///
/// Cell 3f teaches the arrangement to carry inner loops and retires this.
fn reject_holed_operands(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(), BoolError> {
    for solid in [a, b] {
        let shell = model.solids.get(solid).outer;
        if model
            .shells
            .get(shell)
            .faces
            .iter()
            .any(|&fh| !model.faces.get(fh).inner.is_empty())
        {
            return Err(reject(tag::INNER_LOOP_OPERAND));
        }
    }
    Ok(())
}

fn overlap_fuse_cut(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    reject_holed_operands(model, a, b)?;
    let mut planes = collect_planes(model, a)?;
    planes.extend(collect_planes(model, b)?);
    if has_coplanar_pair(&planes) {
        return Err(reject(tag::COPLANAR_PAIR));
    }

    let mut surf_ix: HashMap<Handle<Surface>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        surf_ix.insert(pi.surf, i);
    }

    // Classify every original vertex vs the other solid (exact forward-ray winding).
    let mut classof: HashMap<Handle<Vertex>, Side> = HashMap::new();
    for (verts_solid, other) in [(a, b), (b, a)] {
        for &vh in &solid_vertex_handles(model, verts_solid) {
            let side = point_in_solid(model, model.vertices.get(vh).point, other)?;
            classof.insert(vh, side);
        }
    }

    // Seam vertices: an edge straddling the other boundary pierces exactly one of
    // its faces (that face's plane is the seam vertex's third plane).
    let edges_a = edge_incidence(model, a, &surf_ix);
    let edges_b = edge_incidence(model, b, &surf_ix);
    let mut seam: Vec<SeamVertex> = Vec::new();
    let mut seam_ix: HashMap<[usize; 3], usize> = HashMap::new();
    let mut edge_seam: HashMap<Handle<Edge>, [usize; 3]> = HashMap::new();
    for (edges, other) in [(&edges_a, b), (&edges_b, a)] {
        for &(eh, bounds, ref inc) in edges {
            let [v0, v1] = bounds;
            let (p0, p1) = (model.vertices.get(v0).point, model.vertices.get(v1).point);
            let (s0, s1) = (classof[&v0], classof[&v1]);
            if s0 == s1 {
                // Same-side endpoints, yet the edge pierces exactly one outer face of
                // the other solid. This is *not* the out→in→out tunnel it looks like:
                // such an edge crosses two faces and `pierced_face` already rejected
                // it above as `pierced_multi`. Nor is it parity-impossible — that
                // argument only holds for a cavity-free solid.
                //
                // The one way in is a cavitied operand. `point_in_solid` counts the
                // cavity shells, so an edge running from outside the solid to a point
                // inside a void classifies `Outside` at both ends; `pierced_face`
                // scans the outer shell alone, so it reports one crossing.
                //
                // `boolean` now rejects cavitied operands at the door
                // (`tag::HOLLOW_OPERAND`), so nothing reaches here through the public
                // entry. This stays as defense-in-depth: for direct callers of
                // `overlap_fuse_cut` (the test that pins it does exactly that), and
                // for the day sub-unit 5 relaxes that door to admit cavities — the
                // asymmetry above is what would still stop a seam being built from a
                // boundary the classifier disagrees with.
                if pierced_face(model, p0, p1, other, &surf_ix)?.is_some() {
                    return Err(reject(tag::TUNNEL));
                }
                continue;
            }
            let entry = pierced_face(model, p0, p1, other, &surf_ix)?
                .ok_or_else(|| reject(tag::NO_ENTRY_FACE))?;
            let [e0, e1] = [inc[0], inc[1]];
            let point = three_planes(&planes[e0].plane, &planes[e1].plane, &planes[entry].plane)
                .ok_or_else(|| reject(tag::THREE_PLANES))?;
            // Exact 4-plane-concurrency guard: the seam vertex must not lie on any
            // *other* plane (a degenerate 4-plane meet). Unlike the convex path we do
            // NOT require it inside every half-space (`== -1`): a seam vertex of a
            // non-convex solid can be outside some face's half-space yet on the
            // boundary — `pierced_face` already certified it is on the entry face.
            for (m, pm) in planes.iter().enumerate() {
                if m == e0 || m == e1 || m == entry {
                    continue;
                }
                if three_plane_orient3d(
                    &planes[e0].plane,
                    &planes[e1].plane,
                    &planes[entry].plane,
                    pm.tri[0],
                    pm.tri[1],
                    pm.tri[2],
                ) == 0
                {
                    return Err(reject(tag::FOURPLANE)); // seam vertex on a 4th plane
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

    let (keep_a, keep_b, flip_b) = match kind {
        BoolKind::Fuse => (Side::Outside, Side::Outside, false),
        BoolKind::Cut => (Side::Outside, Side::Inside, true),
        BoolKind::Common => return Err(reject(tag::COMMON_OVERLAP)),
    };
    if seam.is_empty() {
        return contained_result(model, kind, a, b, &classof);
    }

    let mut faces: Vec<LocalFace> = Vec::new();
    for (solid, keep, flip) in [(a, keep_a, false), (b, keep_b, flip_b)] {
        let shell = model.solids.get(solid).outer;
        for &fh in &model.shells.get(shell).faces {
            let face = model.faces.get(fh);
            let pidx = surf_ix[&face.surface];
            if let Some(lf) = reconstruct_face(
                model, face, pidx, keep, flip, &classof, &edge_seam, &seam, &seam_ix, true,
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
        return Err(reject(tag::NONCONVEX_OPERAND));
    }
    let mut planes = planes_a;
    planes.extend(planes_b);
    // Coplanar faces (within one input — e.g. an imprinted face's outer+region —
    // or shared across inputs) give the enumeration duplicate half-spaces.
    if has_coplanar_pair(&planes) {
        return Err(reject(tag::COPLANAR_PAIR));
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
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<PlaneInfo>, BoolError> {
    let shell = model.solids.get(solid).outer;
    let mut out = Vec::new();
    for &fh in &model.shells.get(shell).faces {
        let face = model.faces.get(fh);
        let plane = match model.surfaces.get(face.surface) {
            Surface::Plane(p) => *p,
            Surface::Cylinder(_) => return Err(reject(tag::CYLINDER_FACE)),
        };
        let tri = outer_tri(model, face).ok_or_else(|| reject(tag::DEGENERATE_FACE))?;
        let n_out = (tri[1] - tri[0])
            .cross(tri[2] - tri[0])
            .normalize()
            .ok_or_else(|| reject(tag::DEGENERATE_NORMAL))?;
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
///
/// Outer loops only, and that is complete: manifoldness puts every hole-ring
/// vertex on the outer loop of an adjacent wall face too.
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
                    // `+1` = the vertex is on `pm.tri`'s RH-normal side, and an
                    // outer-CCW triple's RH normal is the outward normal — so `+1`
                    // is "outside this half-space" (see `PlaneInfo`).
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
                    return Err(reject(tag::FOURPLANE)); // 4-plane concurrency
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
pub(crate) fn vertex_tol(p: Point3, a: &Plane, b: &Plane, c: &Plane) -> f64 {
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
                _ => return Err(reject(tag::TANGENT_EDGE)), // tangent / degenerate edge
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
                return Err(reject(tag::NONMANIFOLD_FACE));
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
    let c = Point3::centroid(&pts).ok_or_else(|| reject(tag::DEGENERATE_CENTROID))?;
    let u = (pts[0] - c)
        .normalize()
        .ok_or_else(|| reject(tag::DEGENERATE_RADIUS))?;
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
        return Err(reject(tag::COLLINEAR_FACE));
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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Side {
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
    // Defense-in-depth: a convex operand cannot reach here with a holed face, since
    // a hole in an outer-shell face is bounded either by inward walls (making the
    // solid non-convex) or by a coplanar region face (imprint), and `has_coplanar_pair`
    // below catches the latter. That argument depends on the order of the two checks,
    // so the guard is explicit rather than implied.
    reject_holed_operands(model, a, b)?;
    let planes_a = collect_planes(model, a)?;
    let planes_b = collect_planes(model, b)?;
    if !is_convex(&planes_a, &solid_vertices(model, a))
        || !is_convex(&planes_b, &solid_vertices(model, b))
    {
        return Err(reject(tag::NONCONVEX_OPERAND));
    }
    let na = planes_a.len();
    let mut planes = planes_a;
    planes.extend(planes_b);
    if has_coplanar_pair(&planes) {
        return Err(reject(tag::COPLANAR_PAIR));
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
            let (p0, p1) = (model.vertices.get(v0).point, model.vertices.get(v1).point);
            let (s0, s1) = (classof[&v0], classof[&v1]);
            if s0 == s1 {
                // Both endpoints outside but the segment passes through the other
                // solid ⇒ this edge pierces a face mid-face (a poke-through) ⇒
                // out of clean-seam coverage.
                if s0 == Side::Outside && segment_enters(p0, p1, other.clone(), &planes) {
                    return Err(reject(tag::POKE_THROUGH));
                }
                continue;
            }
            let (p_out, p_in) = if s0 == Side::Outside {
                (p0, p1)
            } else {
                (p1, p0)
            };
            let entry = enter_face(p_out, p_in, other.clone(), &planes)
                .ok_or_else(|| reject(tag::NO_ENTRY_FACE))?;
            let [e0, e1] = [inc[0], inc[1]];
            let point = three_planes(&planes[e0].plane, &planes[e1].plane, &planes[entry].plane)
                .ok_or_else(|| reject(tag::THREE_PLANES))?;
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
                    _ => return Err(reject(tag::OUTSIDE_OR_FOURPLANE)), // outside or on a 4th plane
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
        // No boundary crossing: one solid contains the other (Cut ⇒ a cavity,
        // Fuse ⇒ the container) or they are disjoint (empty).
        return contained_result(model, kind, a, b, &classof);
    }

    // Reconstruct faces of A (keep_a side) and B (keep_b side, flip_b).
    let mut faces: Vec<LocalFace> = Vec::new();
    for (solid, keep, flip) in [(a, keep_a, false), (b, keep_b, flip_b)] {
        let shell = model.solids.get(solid).outer;
        for &fh in &model.shells.get(shell).faces {
            let face = model.faces.get(fh);
            let pidx = surf_ix[&face.surface];
            if let Some(lf) = reconstruct_face(
                model, face, pidx, keep, flip, &classof, &edge_seam, &seam, &seam_ix, false,
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
            return Err(reject(tag::ON_BOUNDARY)); // on the other solid's boundary
        }
        if sd > scale {
            side = Side::Outside; // outside this half-space ⇒ outside the convex
        }
    }
    Ok(side)
}

/// Deterministic generic ray directions (small coprime integers, none
/// axis-aligned) for the point-in-polyhedron cast. Axis-aligned rays would
/// systematically graze the axis-aligned faces/edges of typical inputs, so the
/// list starts off-axis; retrying down it finds a degeneracy-free ray.
const RAY_DIRECTIONS: [[f64; 3]; 6] = [
    [2.0, 3.0, 5.0],
    [3.0, 5.0, 7.0],
    [5.0, 7.0, 11.0],
    [7.0, 11.0, 13.0],
    [11.0, 13.0, 17.0],
    [13.0, 17.0, 19.0],
];

/// Exact point-in-polyhedron for a general (non-convex, possibly hollow) planar
/// solid, by **forward-ray winding** (design §8 M5-d): cast a ray from `p` and
/// sum the oriented crossings ([`ray_face_cross`]) over every face of every shell
/// (outer + cavities), each face fanned from its first vertex `(v0, vi, vi+1)`.
/// The sum is the winding number about `p` — nonzero ⇒ inside the material. A
/// concave face's spurious fan triangles cancel by orientation, so no ear-clip is
/// needed; every decision is an exact `orient3d` sign (no tolerance).
///
/// Degeneracies (the ray grazing an edge/vertex) retry the next
/// [`RAY_DIRECTIONS`]; exhausting the list yields `Unsupported` (astronomically
/// unlikely without adversarial alignment — a symbolic-perturbation upgrade is
/// deferred). **Precondition:** `p` is not on the solid's boundary — the caller's
/// edge-face gate guarantees `p` is not coplanar with any face plane, so the
/// per-face `s0 = 0` case never arises here.
fn point_in_solid(model: &Model, p: Point3, solid: Handle<Solid>) -> Result<Side, BoolError> {
    let faces: Vec<Vec<Vec<Point3>>> = solid_faces(model, solid)
        .into_iter()
        .map(|fh| face_loops(model, fh))
        .collect();
    'dirs: for dir in RAY_DIRECTIONS {
        let d = Vector3::from_array(dir);
        let mut winding = 0i32;
        for rings in &faces {
            for ring in rings {
                for tri in fan_triangles(ring, 0) {
                    match ray_face_cross(p, d, tri) {
                        RayCross::Cross(sign) => winding += sign as i32,
                        RayCross::Miss => {}
                        RayCross::Degenerate => continue 'dirs, // grazed — try another direction
                    }
                }
            }
        }
        return Ok(if winding != 0 {
            Side::Inside
        } else {
            Side::Outside
        });
    }
    Err(reject(tag::RAY_DEGENERATE)) // every direction grazed the boundary (adversarial)
}

/// Every face of a solid's boundary — outer shell then each cavity shell.
fn solid_faces(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Face>> {
    let s = model.solids.get(solid);
    std::iter::once(&s.outer)
        .chain(s.cavities.iter())
        .flat_map(|&sh| model.shells.get(sh).faces.iter().copied())
        .collect()
}

/// A planar face's rings of vertex points: the outer loop, then each hole.
///
/// Fanning only the outer ring would fill the holes in, and every winding-number
/// classifier below would call a point over the hole "material". A hole ring runs
/// CW about the face's outward normal — not by convention but by manifoldness,
/// since each of its edges is used once here and once, oppositely, on the outer
/// loop of the adjacent wall face (`validate` enforces it as `NonOpposedEdge`).
/// So the hole's fan triangles are oriented against the outer ring's and its
/// signed crossings subtract. No ear-clipping, no tolerance — the same
/// orientation-cancellation the concave-outer fan already relies on.
pub(crate) fn face_loops(model: &Model, fh: Handle<Face>) -> Vec<Vec<Point3>> {
    let face = model.faces.get(fh);
    std::iter::once(&face.outer)
        .chain(face.inner.iter())
        .map(|lp| {
            lp.half_edges
                .iter()
                .map(|&he| model.vertices.get(he_start(model, he)).point)
                .collect()
        })
        .collect()
}

/// The non-degenerate fan triangles `(pts[apex], pts[apex+s], pts[apex+s+1])` of a
/// planar loop, fanned from vertex `apex` (indices mod `k`). A concave loop's
/// spurious (reflex) triangles are kept — they cancel by orientation in the
/// oriented crossing sum — but zero-area (collinear) triangles are dropped.
/// Varying `apex` changes which internal diagonals appear, which the segment gate
/// exploits to sidestep a diagonal that happens to be coplanar with a query edge.
fn fan_triangles(pts: &[Point3], apex: usize) -> Vec<[Point3; 3]> {
    let k = pts.len();
    let mut tris = Vec::new();
    for s in 1..k.saturating_sub(1) {
        let (t0, t1, t2) = (pts[apex], pts[(apex + s) % k], pts[(apex + s + 1) % k]);
        let (e1, e2) = (t1 - t0, t2 - t0);
        if e1.cross(e2).norm() <= 1e-12 * e1.norm() * e2.norm() {
            continue; // degenerate (collinear) fan triangle
        }
        tris.push([t0, t1, t2]);
    }
    tris
}

/// Distinct boundary edges of a solid (outer + cavity shells) as endpoint-point
/// pairs. First-seen order → deterministic.
///
/// Outer loops only, and that is complete: manifoldness uses every edge exactly
/// twice, so a hole ring's edges also appear on the outer loop of the adjacent
/// wall face.
fn solid_edges(model: &Model, solid: Handle<Solid>) -> Vec<(Point3, Point3)> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for fh in solid_faces(model, solid) {
        for &he in &model.faces.get(fh).outer.half_edges {
            if seen.insert(he.edge) {
                let [v0, v1] = model.edges.get(he.edge).bounds.expect("bounded");
                out.push((model.vertices.get(v0).point, model.vertices.get(v1).point));
            }
        }
    }
    out
}

/// Whether the boundaries of `a` and `b` actually cross — an edge of one pierces
/// a face of the other (either direction). Each face is fanned and the oriented
/// [`segment_face_cross`] sum decides "pierces this face"; a concave face's
/// spurious triangles cancel. A `Degenerate` contact (a coplanar face, or an edge
/// grazing another's edge/vertex) is out of clean coverage ⇒ `Unsupported`.
///
/// `Ok(false)` means seam-free: the two solids are disjoint or one strictly
/// contains the other (their interiors do not partially overlap) — the only cases
/// this sub-unit's boolean handles. This edge-face test, not vertex
/// classification, is what makes non-convex containment sound: an edge can pierce
/// a reflex region with both endpoints inside, or the solids can interlock with
/// every vertex outside, and only a direct boundary-crossing test catches those.
fn boundaries_intersect(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<bool, BoolError> {
    for (edge_solid, face_solid) in [(a, b), (b, a)] {
        let edges = solid_edges(model, edge_solid);
        let faces: Vec<Vec<Vec<Point3>>> = solid_faces(model, face_solid)
            .into_iter()
            .map(|fh| face_loops(model, fh))
            .collect();
        for (p0, p1) in &edges {
            for rings in &faces {
                if segment_crosses_face(*p0, *p1, rings)? {
                    return Ok(true); // a genuine seam
                }
            }
        }
    }
    Ok(false)
}

/// Whether segment `p0→p1` pierces the *material* of planar face `rings` (outer
/// loop first, then holes), by the oriented [`segment_face_cross`] winding summed
/// over every ring's fan. A segment through a hole gets `+1` from the outer fan and
/// `−1` from the hole's, so it correctly reports no crossing.
///
/// A fan diagonal may be coplanar with an axis-aligned query edge (a spurious
/// `Degenerate`); re-fanning from another apex uses different diagonals, so we
/// retry — for the face as a whole, since one ring's verdict is meaningless without
/// the others. A genuine contact (the segment grazing a *real* edge/vertex) is
/// `Degenerate` from every apex ⇒ `Unsupported`.
pub(crate) fn segment_crosses_face(
    p0: Point3,
    p1: Point3,
    rings: &[Vec<Point3>],
) -> Result<bool, BoolError> {
    let apexes = rings.iter().map(|r| r.len()).max().unwrap_or(0);
    'apex: for apex in 0..apexes {
        let mut crossing = 0i32;
        for ring in rings {
            for tri in fan_triangles(ring, apex % ring.len()) {
                match segment_face_cross(p0, p1, tri) {
                    SegCross::Cross(sign) => crossing += sign as i32,
                    SegCross::Miss => {}
                    SegCross::Degenerate => continue 'apex, // fan diagonal grazed — re-fan the face
                }
            }
        }
        return Ok(crossing != 0);
    }
    Err(reject(tag::CONTACT_DEGENERATE)) // grazed a real edge/vertex from every apex
}

/// Whether a segment (both endpoints outside the convex solid `range`) passes
/// *through* it — Cyrus–Beck line clip yields a non-empty interior interval.
fn segment_enters(
    p0: Point3,
    p1: Point3,
    range: std::ops::Range<usize>,
    planes: &[PlaneInfo],
) -> bool {
    let (mut t_enter, mut t_exit) = (0.0f64, 1.0f64);
    for m in range {
        let d0 = (p0 - planes[m].tri[0]).dot(planes[m].n_out); // >0 outside plane m
        let d1 = (p1 - planes[m].tri[0]).dot(planes[m].n_out);
        if d0 > 0.0 && d1 > 0.0 {
            return false; // segment entirely outside this half-space
        }
        if d0 <= 0.0 && d1 <= 0.0 {
            continue; // no constraint from this plane
        }
        let t = d0 / (d0 - d1);
        if d0 > 0.0 {
            t_enter = t_enter.max(t);
        } else {
            t_exit = t_exit.min(t);
        }
    }
    t_enter < t_exit
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
pub(crate) fn edge_incidence(
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

/// The result when no edge crosses the other solid's boundary: one solid
/// contains the other, or they are disjoint. Direction is read from `classof`
/// (already `Err`-free — a boundary vertex would have failed `classify_vertex`
/// before this point): every B vertex inside A ⇒ B⊂A, every A vertex inside B ⇒
/// A⊂B. `Cut(A−B)` with B⊂A adds B as an inward cavity; A⊂B removes A entirely;
/// `Fuse` yields the container; `Common` yields the contained one. The convex
/// path routes only `Fuse`/`Cut` here (its `Common` enumerates directly); the
/// non-convex seam-free path routes all three kinds here.
fn contained_result(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    classof: &HashMap<Handle<Vertex>, Side>,
) -> Result<Handle<Solid>, BoolError> {
    let all_inside = |s: Handle<Solid>| {
        solid_vertex_handles(model, s)
            .iter()
            .all(|vh| classof.get(vh) == Some(&Side::Inside))
    };
    let b_in_a = all_inside(b);
    let a_in_b = all_inside(a);
    match (kind, b_in_a, a_in_b) {
        (BoolKind::Cut, true, false) => Ok(supersede_with_cavity(model, a, b)), // A − B = A + void(B)
        (BoolKind::Cut, false, true) => Err(BoolError::EmptyResult),            // A ⊂ B ⇒ A removed
        (BoolKind::Cut, false, false) => Ok(supersede_reuse(model, a, b, a)), // disjoint ⇒ A − B = A
        (BoolKind::Fuse, true, false) => Ok(supersede_reuse(model, a, b, a)), // A ∪ B = A
        (BoolKind::Fuse, false, true) => Ok(supersede_reuse(model, a, b, b)), // A ∪ B = B
        (BoolKind::Fuse, false, false) => Err(BoolError::EmptyResult), // disconnected union unrepresentable
        (BoolKind::Common, true, false) => Ok(supersede_reuse(model, a, b, b)), // A ∩ B = B (B ⊂ A)
        (BoolKind::Common, false, true) => Ok(supersede_reuse(model, a, b, a)), // A ∩ B = A (A ⊂ B)
        (BoolKind::Common, false, false) => Err(BoolError::EmptyResult), // disjoint ⇒ empty
        (_, true, true) => {
            unreachable!("mutual containment means a shared boundary — rejected before this point")
        }
    }
}

/// `A` with `B`'s boundary attached as an inward-oriented cavity — the result of
/// `Cut(A − B)` when B is strictly inside A. Reuses A's outer shell and B's
/// cells (via [`Model::reversed_shell`]); supersedes both inputs.
fn supersede_with_cavity(model: &mut Model, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
    let b_outer = model.solids.get(b).outer;
    let void = model.reversed_shell(b_outer);
    let outer = model.solids.get(a).outer;
    let solid = model.push_solid(Solid {
        outer,
        cavities: vec![void],
    });
    model.live_solids.retain(|&s| s != a && s != b);
    solid
}

/// The solid `keep` re-emitted as a fresh live solid (reusing its outer shell),
/// superseding both inputs — the result of a containment `Fuse` (the union is
/// just the container).
fn supersede_reuse(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    keep: Handle<Solid>,
) -> Handle<Solid> {
    let outer = model.solids.get(keep).outer;
    let solid = model.push_solid(Solid {
        outer,
        cavities: vec![],
    });
    model.live_solids.retain(|&s| s != a && s != b);
    solid
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
    strict: bool,
) -> Result<Option<LocalFace>, BoolError> {
    let hes = &face.outer.half_edges;
    let n = hes.len();
    let verts: Vec<Handle<Vertex>> = hes.iter().map(|&he| he_start(model, he)).collect();
    let kept: Vec<bool> = verts.iter().map(|v| classof[v] == keep).collect();

    // Count boundary crossings (kept↔dropped edge transitions).
    let transitions: Vec<usize> = (0..n).filter(|&i| kept[i] != kept[(i + 1) % n]).collect();
    match transitions.len() {
        0 => {
            // No boundary crossing: an interior seam on this plane means the other
            // solid pokes through this face (a hole) — out of clean-seam coverage.
            if seam.iter().any(|s| s.triple.contains(&plane_idx)) {
                return Err(reject(tag::POKEHOLE));
            }
            if kept[0] {
                let loop_nodes = verts.iter().map(|&v| Node::Orig(v)).collect();
                Ok(Some(LocalFace {
                    plane_idx,
                    loop_nodes,
                    flip,
                }))
            } else {
                Ok(None)
            }
        }
        2 => {
            // One kept run; splice the seam sub-path across the dropped run.
            // Transition edge i (kept[i]!=kept[i+1]) carries seam vertex edge_seam[edge].
            let s_at = |i: usize| -> Result<[usize; 3], BoolError> {
                edge_seam
                    .get(&hes[i].edge)
                    .copied()
                    .ok_or_else(|| reject(tag::MISSING_SEAM))
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
            // `strict` (non-convex overlap): the bend-sort projects the seam arc
            // onto the chord direction, which assumes a convex (monotone) arc — a
            // non-convex operand can fold the arc back, and the linear sort would
            // then build a self-intersecting face that `validate` cannot catch. Guard
            // it: the sorted sub-path `p_kd → bends → p_dk` must turn one way about
            // the face normal (a convex chain). Any sign flip ⇒ honest `Unsupported`
            // (the true arrangement is a later sub-unit). Convex `fuse_cut` passes
            // `strict = false` and is unaffected.
            if strict && !bends.is_empty() {
                let n = {
                    let t = outer_tri(model, face).ok_or_else(|| reject(tag::DEGENERATE_FACE))?;
                    (t[1] - t[0]).cross(t[2] - t[0])
                };
                let mut arc: Vec<Point3> = Vec::with_capacity(bends.len() + 2);
                arc.push(p_kd);
                arc.extend(bends.iter().map(|t| seam[seam_ix[t]].point));
                arc.push(p_dk);
                let mut sign = 0i32;
                for w in arc.windows(3) {
                    let turn = (w[1] - w[0]).cross(w[2] - w[1]).dot(n);
                    let s = if turn > 0.0 {
                        1
                    } else if turn < 0.0 {
                        -1
                    } else {
                        0
                    };
                    if s != 0 {
                        if sign != 0 && sign != s {
                            return Err(reject(tag::STRICTARC)); // non-convex seam arc
                        }
                        sign = s;
                    }
                }
            }
            nodes.push(Node::Seam(s_kd));
            nodes.extend(bends.into_iter().map(Node::Seam));
            nodes.push(Node::Seam(s_dk));
            Ok(Some(LocalFace {
                plane_idx,
                loop_nodes: nodes,
                flip,
            }))
        }
        _ => Err(reject(tag::MULTICHORD)), // ≥4 crossings ⇒ multiple chords
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

// ---- coincident-coplanar merge (M5-c5): matched-interface stack ----

/// A detected clean coincident-coplanar glue: the two faces meeting on the
/// shared plane (dropped by fuse) and the B→A correspondence for their shared
/// interface ring.
struct Interface {
    fa: Handle<Face>,
    fb: Handle<Face>,
    remap: HashMap<Handle<Vertex>, Handle<Vertex>>,
}

/// `Some` iff A and B share exactly one fully-coincident, opposite-normal face
/// pair (identical boundary) — the clean stack/glue case. `None` (fall through to
/// the coplanar-rejecting paths) for anything else.
fn detect_coincident_interface(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<Interface> {
    let planes_a = collect_planes(model, a).ok()?;
    let planes_b = collect_planes(model, b).ok()?;
    if !is_convex(&planes_a, &solid_vertices(model, a))
        || !is_convex(&planes_b, &solid_vertices(model, b))
    {
        return None;
    }
    // Opposite-normal coplanar face pairs (same-normal coplanar side faces are
    // the expected coplanar-adjacent result and are ignored).
    let mut opposite: Vec<(usize, usize)> = Vec::new();
    for (i, pa) in planes_a.iter().enumerate() {
        for (j, pb) in planes_b.iter().enumerate() {
            if coplanar(&pa.plane, &pb.plane) && pa.n_out.dot(pb.n_out) < 0.0 {
                opposite.push((i, j));
            }
        }
    }
    if opposite.len() != 1 {
        return None; // 0 ⇒ not a stack; ≥2 ⇒ multiple interfaces (out of scope)
    }
    let (i, j) = opposite[0];
    let fa = model.shells.get(model.solids.get(a).outer).faces[i];
    let fb = model.shells.get(model.solids.get(b).outer).faces[j];
    let ring_a = face_ring(model, model.faces.get(fa));
    let ring_b = face_ring(model, model.faces.get(fb));
    let remap = interface_correspondence(&ring_a, &ring_b)?;
    Some(Interface { fa, fb, remap })
}

/// A face's outer-loop (start vertex, point) ring, in loop order.
fn face_ring(model: &Model, face: &Face) -> Vec<(Handle<Vertex>, Point3)> {
    face.outer
        .half_edges
        .iter()
        .map(|&he| {
            let vh = he_start(model, he);
            (vh, model.vertices.get(vh).point)
        })
        .collect()
}

/// A bijective coordinate match B-ring → A-ring within a scale-relative tol, or
/// `None` if the boundaries are not identical (different length, an unmatched or
/// ambiguous vertex, or two B vertices sharing an A vertex). Input analysis only.
fn interface_correspondence(
    ring_a: &[(Handle<Vertex>, Point3)],
    ring_b: &[(Handle<Vertex>, Point3)],
) -> Option<HashMap<Handle<Vertex>, Handle<Vertex>>> {
    if ring_a.len() != ring_b.len() {
        return None;
    }
    let mut remap = HashMap::new();
    for &(bh, bp) in ring_b {
        let scale = 1e-9 * (bp.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        let found: Vec<Handle<Vertex>> = ring_a
            .iter()
            .filter(|&&(_, ap)| (ap - bp).norm() <= scale)
            .map(|&(ah, _)| ah)
            .collect();
        if found.len() != 1 {
            return None; // unmatched or ambiguous ⇒ boundaries differ
        }
        remap.insert(bh, found[0]);
    }
    let distinct: HashSet<Handle<Vertex>> = remap.values().copied().collect();
    if distinct.len() != remap.len() {
        return None; // two B vertices mapped to one A vertex
    }
    Some(remap)
}

/// A solid's faces as `LocalFace`s (`Node::Orig`, `flip:false`), with
/// `plane_idx` offset by `plane_offset`, optionally skipping one face and
/// remapping vertices (for B's interface vertices → A's handle).
fn solid_local_faces(
    model: &Model,
    solid: Handle<Solid>,
    plane_offset: usize,
    skip: Option<Handle<Face>>,
    remap: Option<&HashMap<Handle<Vertex>, Handle<Vertex>>>,
) -> Vec<LocalFace> {
    let shell = model.solids.get(solid).outer;
    let mut out = Vec::new();
    for (pos, &fh) in model.shells.get(shell).faces.iter().enumerate() {
        if Some(fh) == skip {
            continue;
        }
        let face = model.faces.get(fh);
        let loop_nodes = face
            .outer
            .half_edges
            .iter()
            .map(|&he| {
                let vh = he_start(model, he);
                let mapped = remap.and_then(|r| r.get(&vh)).copied().unwrap_or(vh);
                Node::Orig(mapped)
            })
            .collect();
        out.push(LocalFace {
            plane_idx: plane_offset + pos,
            loop_nodes,
            flip: false,
        });
    }
    out
}

/// The coincident-merge result (reuses `assemble_fuse_cut`; no seam, no flip).
fn coincident_merge(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    iface: &Interface,
) -> Result<Handle<Solid>, BoolError> {
    // `detect_coincident_interface` counts only *cross-solid, opposite-normal*
    // coplanar pairs, so an imprint on some other face — whose coplanar region face
    // is same-normal and within one solid — does not disqualify the stack. Such an
    // operand reaches here with a holed face; `solid_local_faces` would drop the hole.
    reject_holed_operands(model, a, b)?;
    match kind {
        // The intersection is the flat shared face (zero volume).
        BoolKind::Common => Err(BoolError::EmptyResult),
        // A and B are on opposite sides of the interface, so B removes nothing:
        // A − B is a fresh copy of A (its interface face survives).
        BoolKind::Cut => {
            let planes_a = collect_planes(model, a)?;
            let faces = solid_local_faces(model, a, 0, None, None);
            Ok(assemble_fuse_cut(model, a, b, &planes_a, &[], &faces))
        }
        // Drop both interface faces; keep every other face, sewing B's interface
        // ring to A's shared vertices. Coplanar-adjacent side faces stay separate.
        BoolKind::Fuse => {
            let planes_a = collect_planes(model, a)?;
            let na = planes_a.len();
            let planes_b = collect_planes(model, b)?;
            let mut planes = planes_a;
            planes.extend(planes_b);
            let mut faces = solid_local_faces(model, a, 0, Some(iface.fa), None);
            faces.extend(solid_local_faces(
                model,
                b,
                na,
                Some(iface.fb),
                Some(&iface.remap),
            ));
            Ok(assemble_fuse_cut(model, a, b, &planes, &[], &faces))
        }
    }
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

    /// The L-prism: profile `[(0,0),(2,0),(2,1),(1,1),(1,2),(0,2)]` extruded to
    /// z ∈ [0,1]. Material = bottom bar (x∈[0,2],y∈[0,1]) ∪ left bar (x∈[0,1],
    /// y∈[1,2]); the notch (x∈[1,2],y∈[1,2]) is empty.
    fn l_prism() -> (Model, Handle<Solid>) {
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
        let s = m.live_solids[0];
        (m, s)
    }

    #[test]
    fn point_in_solid_classifies_concave_prism() {
        let (m, s) = l_prism();
        let side = |x, y, z| point_in_solid(&m, Point3::from_array([x, y, z]), s).unwrap();
        // Interior of the arm, near the reflex corner, and the left bar.
        assert_eq!(side(0.5, 0.5, 0.5), Side::Inside);
        assert_eq!(side(0.9, 0.9, 0.5), Side::Inside);
        assert_eq!(side(0.5, 1.5, 0.5), Side::Inside);
        // The notch is OUTSIDE the L though inside its bounding box — the concave
        // case a convex all-half-spaces test gets wrong.
        assert_eq!(side(1.5, 1.5, 0.5), Side::Outside);
        // Clearly outside (beside, below, above).
        assert_eq!(side(3.0, 3.0, 0.5), Side::Outside);
        assert_eq!(side(0.5, 0.5, -1.0), Side::Outside);
        assert_eq!(side(0.5, 0.5, 2.0), Side::Outside);
    }

    #[test]
    fn point_in_solid_classifies_cube() {
        let mut m = Model::new();
        let s = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        assert_eq!(
            point_in_solid(&m, Point3::from_array([1.0, 1.0, 1.0]), s).unwrap(),
            Side::Inside
        );
        assert_eq!(
            point_in_solid(&m, Point3::from_array([3.0, 1.0, 1.0]), s).unwrap(),
            Side::Outside
        );
    }

    #[test]
    fn point_in_solid_handles_cavity() {
        // 4-cube with a concentric 2-cube void [1,3]³: a point in the material
        // wall is Inside, a point in the empty void is Outside (the winding sums
        // outer + cavity shells).
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let b_outer = m.solids.get(b).outer;
        let void = m.reversed_shell(b_outer);
        let a_outer = m.solids.get(a).outer;
        let hollow = m.push_solid(Solid {
            outer: a_outer,
            cavities: vec![void],
        });
        m.live_solids.retain(|&x| x == hollow);
        assert_eq!(
            point_in_solid(&m, Point3::from_array([0.5, 0.5, 0.5]), hollow).unwrap(),
            Side::Inside // in the wall
        );
        assert_eq!(
            point_in_solid(&m, Point3::from_array([2.0, 2.0, 2.0]), hollow).unwrap(),
            Side::Outside // in the void
        );
    }

    /// The L-prism with a `[0.1,0.9]³` box strictly inside its bottom bar
    /// (non-coplanar coordinates ⇒ no shared face planes). `V_L = 3`, `V_box =
    /// 0.512`.
    fn l_and_inner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
        (m, l, bx)
    }

    #[test]
    fn cut_non_convex_containment_makes_cavity() {
        // Cut(L − box) with box ⊂ L ⇒ a hollow L (outer L shell + box void).
        let (mut m, l, bx) = l_and_inner_box();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
    }

    #[test]
    fn fuse_non_convex_containment_is_container() {
        let (mut m, l, bx) = l_and_inner_box();
        let vol_l = nacre_props::mass_props(&m, l).unwrap().volume;
        let r = boolean(&mut m, BoolKind::Fuse, l, bx).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_l).abs() < 1e-9);
        assert!(m.solids.get(r).cavities.is_empty());
    }

    #[test]
    fn common_non_convex_containment_is_inner() {
        let (mut m, l, bx) = l_and_inner_box();
        let vol_bx = nacre_props::mass_props(&m, bx).unwrap().volume;
        let r = boolean(&mut m, BoolKind::Common, l, bx).unwrap();
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_bx).abs() < 1e-9);
    }

    #[test]
    fn cut_box_inside_non_convex_is_empty() {
        // Cut(box − L): the box is wholly inside L ⇒ nothing remains.
        let (mut m, l, bx) = l_and_inner_box();
        assert_eq!(
            boolean(&mut m, BoolKind::Cut, bx, l),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn disjoint_non_convex_operand() {
        // L and a far box (non-coplanar): Cut ⇒ L, Fuse ⇒ empty, Common ⇒ empty.
        let far = || Point3::from_array([10.0; 3]);
        let far_max = || Point3::from_array([11.0; 3]);
        let (mut m, l) = l_prism();
        let vol_l = nacre_props::mass_props(&m, l).unwrap().volume;
        let d = m.add_cuboid(far(), far_max());
        let r = boolean(&mut m, BoolKind::Cut, l, d).unwrap();
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_l).abs() < 1e-9);
        let (mut m, l) = l_prism();
        let d = m.add_cuboid(far(), far_max());
        assert_eq!(
            boolean(&mut m, BoolKind::Fuse, l, d),
            Err(BoolError::EmptyResult)
        );
        let (mut m, l) = l_prism();
        let d = m.add_cuboid(far(), far_max());
        assert_eq!(
            boolean(&mut m, BoolKind::Common, l, d),
            Err(BoolError::EmptyResult)
        );
    }

    /// The L-prism with a box biting its convex corner `(2, 0)` — the first
    /// non-convex *overlap* (a real single-chord seam), M5-d2. The box spans
    /// `x∈[1.3,2.4]`, `y∈[-0.3,0.4]`, `z∈[0.2,1.4]`: it straddles the corner in x
    /// and y, and its z-range pokes above the L (`z=1`) while its floor `z=0.2`
    /// sits inside — so every crossing edge is a clean straddle (no edge tunnels
    /// fully through the other) and no box face is coplanar with an L face. The
    /// span is deliberately asymmetric so no seam point lands on a face centre
    /// (where both fan diagonals cross and every apex would graze).
    /// Overlap = `x∈[1.3,2]·y∈[0,0.4]·z∈[0.2,1]` = `0.224`;
    /// `V_L=3`, `V_box=1.1·0.7·1.2=0.924`.
    fn l_and_corner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([1.3, -0.3, 0.2]),
            Point3::from_array([2.4, 0.4, 1.4]),
        );
        (m, l, bx)
    }

    #[test]
    fn cut_non_convex_overlap_corner_bite() {
        // Cut(L − box): the corner bite carves 0.175 off the L.
        let (mut m, l, bx) = l_and_corner_box();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.224)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn fuse_non_convex_overlap_corner_bite() {
        // Fuse(L ∪ box): the protruding box adds (0.924 − 0.224) to the L.
        let (mut m, l, bx) = l_and_corner_box();
        let r = boolean(&mut m, BoolKind::Fuse, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 0.924 - 0.224)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn common_non_convex_overlap_is_unsupported() {
        // Non-convex ∩ needs the full arrangement (a later sub-unit) — honest reject.
        let (mut m, l, bx) = l_and_corner_box();
        assert_rejects(
            || boolean(&mut m, BoolKind::Common, l, bx),
            tag::COMMON_OVERLAP,
        );
    }

    #[test]
    fn cut_across_reflex_corner_bite() {
        // A box straddling the reflex corner (1,1) leaves a *single* chord with one
        // reflex bend — still transitions==2, so it reconstructs into a correct
        // L-shaped face. This is the strongest classification test: the box vertex
        // (1.6,1.6,·) sits in the L's notch (inside the convex hull, outside the L),
        // exactly where a convex half-space test would misclassify it Inside; only
        // exact `point_in_solid` gets the volume right.
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.6, 0.6, 0.2]),
            Point3::from_array([1.6, 1.6, 1.4]),
        );
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        // Overlap = xy(1.0 − notch 0.36 = 0.64) · z(0.8) = 0.512.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn poke_through_hole_is_unsupported() {
        // A thin rod skewering the L's bottom bar in z (both ends outside) would
        // drill a through-hole (an inner loop) — not yet representable. Its vertical
        // edges tunnel through the bar, so `pierced_face` sees >1 face ⇒ reject.
        let (mut m, l) = l_prism();
        let rod = m.add_cuboid(
            Point3::from_array([0.3, 0.3, -0.5]),
            Point3::from_array([0.5, 0.6, 1.5]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, l, rod),
            tag::PIERCED_MULTI,
        );
    }

    /// A U-prism: a bottom bar `y∈[0,1]` with two prongs rising from it. The prong
    /// tops sit at *different* heights (y=2.3 and y=2.0) on purpose — level tops
    /// would be coplanar faces and `has_coplanar_pair` would reject before the seam
    /// machinery ran. Area 3 + 1 + 1.3, extruded 1.0 ⇒ volume 5.3.
    fn u_prism() -> (Model, Handle<Solid>) {
        let u = Profile2d {
            points: vec![
                p2(0.0, 0.0),
                p2(3.0, 0.0),
                p2(3.0, 2.3),
                p2(2.0, 2.3),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
            ],
        };
        let m = replay(&[extrude_op(u, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// The U with a slab shearing off both prong tops. The slab overhangs the U in
    /// x and z, so **every slab edge lies outside the U** (pierces nothing) and every
    /// U edge either straddles cleanly (one crossing of the slab's `y=1.5` face) or
    /// misses. Neither `pierced_multi` nor `tunnel` can fire.
    fn u_and_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, u) = u_prism();
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.5, -0.5]),
            Point3::from_array([3.5, 2.5, 1.5]),
        );
        (m, u, slab)
    }

    #[test]
    fn u_prism_is_valid() {
        // Pin the fixture itself: a mistyped profile could still trip `multichord`
        // below, for the wrong reason.
        let (m, u) = u_prism();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, u).unwrap().volume;
        assert!((vol - 5.3).abs() < 1e-9, "volume {vol}");
    }

    // ---- arrangement: seam segment gathering (M5-d3 cell 3b) ----

    /// Combined plane list (A then B) plus the surface→index map, as the boolean
    /// builds them.
    fn combined(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
    ) -> (Vec<PlaneInfo>, HashMap<Handle<Surface>, usize>) {
        let mut planes = collect_planes(m, a).unwrap();
        planes.extend(collect_planes(m, b).unwrap());
        let surf_ix = planes
            .iter()
            .enumerate()
            .map(|(i, pi)| (pi.surf, i))
            .collect();
        (planes, surf_ix)
    }

    /// The face of `solid` whose outward normal is `n` (there is exactly one, for
    /// the axis-aligned fixtures here).
    fn face_facing(
        m: &Model,
        solid: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Surface>, usize>,
        n: [f64; 3],
    ) -> Handle<Face> {
        let want = Vector3::from_array(n);
        let shell = m.solids.get(solid).outer;
        let mut hit = None;
        for &fh in &m.shells.get(shell).faces {
            let pi = &planes[surf_ix[&m.faces.get(fh).surface]];
            if (pi.n_out - want).norm() < 1e-9 {
                assert!(hit.is_none(), "two faces share an outward normal");
                hit = Some(fh);
            }
        }
        hit.expect("no face with that outward normal")
    }

    fn segs(
        m: &Model,
        f: Handle<Face>,
        x: Handle<Solid>,
        y: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Surface>, usize>,
    ) -> Vec<arrange::SeamSegment> {
        let inc_x = arrange::edge_planes(m, x, surf_ix).unwrap();
        let inc_y = arrange::edge_planes(m, y, surf_ix).unwrap();
        arrange::seam_segments_on(m, f, y, planes, surf_ix, &inc_x, &inc_y).unwrap()
    }

    fn near(a: Point3, b: [f64; 3]) -> bool {
        (a - Point3::from_array(b)).norm() < 1e-9
    }

    proptest! {
        /// Whatever the sort, the pairing and the parity sweep do, an endpoint is the
        /// meet of the three planes its triple names. This holds independently of all
        /// of them — and breaks the instant `three_planes` and the triple disagree.
        #[test]
        fn prop_seam_endpoints_lie_on_their_three_planes(
            aext in prop::array::uniform3(1.0f64..2.0),
            bmin in prop::array::uniform3(0.3f64..0.9),
            bext in prop::array::uniform3(1.1f64..2.0),
        ) {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array(aext));
            let b = m.add_cuboid(
                Point3::from_array(bmin),
                Point3::from_array([bmin[0] + bext[0], bmin[1] + bext[1], bmin[2] + bext[2]]),
            );
            let (planes, surf_ix) = combined(&m, a, b);
            let (Ok(inc_x), Ok(inc_y)) = (
                arrange::edge_planes(&m, a, &surf_ix),
                arrange::edge_planes(&m, b, &surf_ix),
            ) else {
                return Ok(());
            };
            let shell = m.solids.get(a).outer;
            for &f in &m.shells.get(shell).faces {
                // A degenerate touch is an honest reject, not a counterexample.
                let Ok(segs) =
                    arrange::seam_segments_on(&m, f, b, &planes, &surf_ix, &inc_x, &inc_y)
                else {
                    continue;
                };
                for s in segs {
                    for (t, p) in s.ends.iter().zip(s.points) {
                        for &i in t {
                            prop_assert!(
                                planes[i].plane.distance(p) < 1e-9,
                                "endpoint {p:?} off plane {i}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn seam_segments_on_two_overlapping_boxes() {
        // A = [0,2]×[0,2.2]×[0,2.4] and B = [1,3.5]×[1,3.2]×[1,3.4] share a corner.
        // Extents are deliberately unequal: with two cubes the seam lands on a face's
        // *centre*, where both fan diagonals cross, and `segment_crosses_face` grazes
        // from every apex and honestly reports `contact_degenerate` (the same trap the
        // `l_and_corner_box` fixture avoids).
        //
        // On A's x=2 face (y∈[0,2.2], z∈[0,2.4]) two of B's faces cut a seam:
        //   B's z=1 face ⇒ the line {x=2, z=1} clipped to y∈[1, 2.2]
        //   B's y=1 face ⇒ the line {x=2, y=1} clipped to z∈[1, 2.4]
        // They meet at (2,1,1), the overlap's corner on this face. Hand-computed.
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.2, 2.4]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0; 3]),
            Point3::from_array([3.5, 3.2, 3.4]),
        );
        let (planes, surf_ix) = combined(&m, a, b);
        let f = face_facing(&m, a, &planes, &surf_ix, [1.0, 0.0, 0.0]);

        let out = segs(&m, f, a, b, &planes, &surf_ix);
        assert_eq!(
            out.len(),
            2,
            "one seam segment per B-face that cuts A's face: {out:?}"
        );

        let mut ends: Vec<[[f64; 3]; 2]> = out
            .iter()
            .map(|s| {
                let mut e = [s.points[0].as_array(), s.points[1].as_array()];
                e.sort_by(|x, y| x.partial_cmp(y).unwrap());
                e
            })
            .collect();
        ends.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let want = {
            let mut w = [
                [[2.0, 1.0, 1.0], [2.0, 2.2, 1.0]],
                [[2.0, 1.0, 1.0], [2.0, 1.0, 2.4]],
            ];
            for e in &mut w {
                e.sort_by(|x, y| x.partial_cmp(y).unwrap());
            }
            w.sort_by(|x, y| x.partial_cmp(y).unwrap());
            w
        };
        for (got, want) in ends.iter().zip(want.iter()) {
            for (g, w) in got.iter().zip(want.iter()) {
                assert!(near(Point3::from_array(*g), *w), "{ends:?} vs {want:?}");
            }
        }

        // Each segment's two endpoint triples differ in exactly one plane: the two they
        // share are the pair that defines the seam line.
        for s in &out {
            let shared: Vec<usize> = s.ends[0]
                .iter()
                .filter(|i| s.ends[1].contains(i))
                .copied()
                .collect();
            assert_eq!(shared.len(), 2);
            assert!(shared.contains(&s.plane_pair[0]) && shared.contains(&s.plane_pair[1]));
        }
    }

    #[test]
    fn seam_segments_split_a_face_into_two_chords() {
        // What `reconstruct_face` cannot do. The slab meets the U-prism's base cap
        // (z=0) along {z=0, y=1.5}; the U has material there only for x∈[0,1] and
        // x∈[2,3], so the seam is *two* segments. Today the boolean rejects this
        // fixture with `multichord`.
        //
        // All four crossings come from the base cap's own edges (neighbour planes
        // x=0,1,2,3); the slab's edges miss the cap, lying outside the U in x or off
        // the z=0 plane.
        let (m, u, slab) = u_and_slab();
        let (planes, surf_ix) = combined(&m, u, slab);
        let cap = face_facing(&m, u, &planes, &surf_ix, [0.0, 0.0, -1.0]);

        let out = segs(&m, cap, u, slab, &planes, &surf_ix);
        assert_eq!(out.len(), 2, "two chords on one face: {out:?}");
        let mut spans: Vec<(f64, f64)> = out
            .iter()
            .map(|s| {
                for e in s.points {
                    assert!(
                        (e[1] - 1.5).abs() < 1e-9 && e[2].abs() < 1e-9,
                        "off the seam line"
                    );
                }
                let (a, b) = (s.points[0][0], s.points[1][0]);
                (a.min(b), a.max(b))
            })
            .collect();
        spans.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        assert!(
            (spans[0].0 - 0.0).abs() < 1e-9 && (spans[0].1 - 1.0).abs() < 1e-9,
            "{spans:?}"
        );
        assert!(
            (spans[1].0 - 2.0).abs() < 1e-9 && (spans[1].1 - 3.0).abs() < 1e-9,
            "{spans:?}"
        );
    }

    /// Every seam-segment endpoint, as a sorted triple of surface handles.
    fn seam_endpoint_triples(
        m: &Model,
        x: Handle<Solid>,
        y: Handle<Solid>,
    ) -> std::collections::BTreeSet<[Handle<Surface>; 3]> {
        let (planes, surf_ix) = combined(m, x, y);
        let inc_x = arrange::edge_planes(m, x, &surf_ix).unwrap();
        let inc_y = arrange::edge_planes(m, y, &surf_ix).unwrap();
        let shell = m.solids.get(x).outer;
        let mut out = std::collections::BTreeSet::new();
        for &f in &m.shells.get(shell).faces {
            let segs =
                arrange::seam_segments_on(m, f, y, &planes, &surf_ix, &inc_x, &inc_y).unwrap();
            for s in segs {
                for t in s.ends {
                    let mut surfs = [planes[t[0]].surf, planes[t[1]].surf, planes[t[2]].surf];
                    surfs.sort_unstable();
                    out.insert(surfs);
                }
            }
        }
        out
    }

    /// The `ThreePlane` definitions of a solid's `Discovered` vertices.
    fn discovered_triples(m: &Model) -> std::collections::BTreeSet<[Handle<Surface>; 3]> {
        let reach = m.reachable();
        let mut out = std::collections::BTreeSet::new();
        for vh in &reach.vertices {
            if let Origin::Discovered {
                definition: VertexDef::ThreePlane(mut t),
                ..
            } = m.vertices.get(*vh).origin
            {
                t.sort_unstable();
                out.insert(t);
            }
        }
        out
    }

    #[test]
    fn seam_endpoints_match_the_boolean_result_vertices() {
        // An independent oracle. `assemble_fuse_cut` stamps every result vertex it
        // discovers with `VertexDef::ThreePlane` — the same identity a seam-segment
        // endpoint carries. The two are computed along completely different routes:
        // one walks edges through faces and reconstructs, the other clips a plane-pair
        // line to a face pair. On a clean single-chord overlap they must agree exactly.
        //
        // (They would not on `poke_through` or `tunnel` fixtures, where the reconstruction
        // rejects part-way through enumeration while the face-local gathering completes.)
        for kind in [BoolKind::Cut, BoolKind::Fuse] {
            let (mut m, l, bx) = l_and_corner_box();
            let from_arrange = seam_endpoint_triples(&m, l, bx);
            assert_eq!(from_arrange.len(), 6);
            boolean(&mut m, kind, l, bx).unwrap();
            m.rebuild_adjacency();
            assert_eq!(from_arrange, discovered_triples(&m), "{kind:?}");
        }

        // The reflex-corner bite too: eight seam vertices, one of them where the box
        // straddles the L's notch.
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.6, 0.6, 0.2]),
            Point3::from_array([1.6, 1.6, 1.4]),
        );
        let from_arrange = seam_endpoint_triples(&m, l, bx);
        assert_eq!(from_arrange.len(), 8);
        boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));
    }

    #[test]
    fn overlap_with_a_holed_face_is_unsupported() {
        // The box bites a corner far from the pocket, so no rim edge straddles the
        // seam and nothing stops `reconstruct_face` from rebuilding the lid without
        // its hole. Measured before the guard: `Ok`, volume 0.96996 against a correct
        // 0.916625, and `validate` reporting five violations.
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.85, 0.85, 0.85]),
            Point3::from_array([1.15, 1.15, 1.15]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, pc, bx),
            tag::INNER_LOOP_OPERAND,
        );
    }

    #[test]
    fn overlap_across_a_hole_rim_is_unsupported() {
        // Here the box crosses the pocket rim, so a rim edge straddles. That edge is
        // used once by the lid's inner loop and once by a wall's outer loop, but
        // `edge_incidence` walks outer loops only and hands back a single incident
        // plane. Measured before the guard: a panic indexing `inc[1]`.
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.55, 0.55, 0.85]),
            Point3::from_array([1.15, 1.15, 1.15]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, pc, bx),
            tag::INNER_LOOP_OPERAND,
        );
    }

    #[test]
    fn coincident_merge_with_a_holed_face_is_unsupported() {
        // `detect_coincident_interface` counts only cross-solid opposite-normal
        // coplanar pairs, so imprinting a *different* face leaves the stack looking
        // clean and routes into `coincident_merge`, whose `solid_local_faces` drops
        // the hole. Measured before the guard: `Ok`, volume 2.0533 against 2.0.
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude_op(square(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        let side = faces[3]; // the x=1 face
        let OpOutput::ImprintSketch { solid, .. } = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: side,
                profile: small_square(),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let bx = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, solid, bx),
            tag::INNER_LOOP_OPERAND,
        );
    }

    #[test]
    fn cut_multi_chord_slab_is_unsupported() {
        // The slab cuts both prongs, so the U's base cap alternates kept/dropped four
        // times around its loop — two chords on one face. Each crossing is contributed
        // by a *different* straddle edge piercing exactly one face, which is why the
        // per-edge guards let it through to `reconstruct_face`.
        let (mut m, u, slab) = u_and_slab();
        assert_rejects(|| boolean(&mut m, BoolKind::Cut, u, slab), tag::MULTICHORD);
    }

    /// The L with a box biting its reflex corner and poking out the top. The box top
    /// (z=1.2) clears the L's z=1 **deliberately**: sunk inside the L's slab, the L's
    /// vertical edges at (2,1) and (1,1) would pierce the box's bottom *and* top face
    /// and `pierced_multi` would reject first.
    fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.2]),
            Point3::from_array([2.5, 1.5, 1.2]),
        );
        (m, l, bx)
    }

    #[test]
    fn cut_staircase_seam_arc_is_unsupported() {
        // On the box's bottom face the seam runs (2,0.5) → (2,1) → (1,1) → (1,1.5):
        // a staircase whose two bends turn opposite ways, so `strict` rejects it.
        //
        // This is an **over-rejection**, and deliberately so. The arc is simple, the
        // projection-sort orders it correctly, and the reconstructed face would have
        // been valid. `strict` cannot tell a reflex turn from an arc folded back on
        // itself, and a folded arc builds a self-intersecting face that `validate`
        // accepts (still manifold, Euler holds). Rejecting both is the sound trade
        // until the real arrangement lands (sub-unit 3).
        let (mut m, l, bx) = l_and_popup_box();
        assert_rejects(|| boolean(&mut m, BoolKind::Cut, l, bx), tag::STRICTARC);
    }

    /// The L with a stub rising out of its top face, footprint strictly inside that
    /// face. Unlike the rod of `poke_through_hole_is_unsupported`, the stub enters
    /// the L from within, so each of its vertical edges crosses exactly one face
    /// (`pierced_multi` cannot fire) and its bottom ring stays inside (`tunnel`
    /// cannot fire either).
    fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let stub = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 0.5]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        (m, l, stub)
    }

    #[test]
    fn cut_blind_dimple_is_unsupported() {
        // The stub's footprint never reaches the top face's boundary, so that face's
        // loop has zero kept/dropped transitions — yet four seam vertices sit on its
        // plane. The seam is a closed ring in the face interior: an inner loop, which
        // `LocalFace` cannot carry. Honest reject until sub-unit 3 emits inner loops.
        let (mut m, l, stub) = l_and_dimple();
        assert_rejects(|| boolean(&mut m, BoolKind::Cut, l, stub), tag::POKEHOLE);
    }

    /// The unit cube with a 0.4-square pocket, 0.5 deep, in its top face: the void
    /// is `[0.3,0.7]² × [0.5,1]` and the solid measures `1 − 0.16·0.5 = 0.92`. Its
    /// lid is the only face in the suite that carries an inner loop.
    fn pocketed_cube() -> (Model, Handle<Solid>) {
        let (mut m, top) = cube_with_top();
        let OpOutput::PocketOnFace { solid, .. } =
            apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        (m, solid)
    }

    #[test]
    fn point_in_solid_sees_through_a_pocket() {
        // Fanning the lid's outer ring alone fills the pocket mouth in, and a ray
        // leaving the void through it counts a crossing that is not there. Worse, the
        // error is not even uniform: rays that exit sideways through a pocket wall
        // miss the lid entirely, so neighbouring points disagree.
        let (m, pc) = pocketed_cube();
        assert!((nacre_props::mass_props(&m, pc).unwrap().volume - 0.92).abs() < 1e-9);
        let at = |p: [f64; 3]| point_in_solid(&m, Point3::from_array(p), pc).unwrap();
        assert_eq!(at([0.5, 0.5, 0.75]), Side::Outside); // in the void
        assert_eq!(at([0.1, 0.1, 0.75]), Side::Inside); // in the wall around it
        assert_eq!(at([0.5, 0.5, 0.25]), Side::Inside); // under the pocket floor
    }

    #[test]
    fn cut_by_a_box_inside_the_pocket_is_a_no_op() {
        // The box sits wholly in the void, so the solids are disjoint and `A − B = A`.
        // Reading the lid as filled instead classified the box's eight corners five
        // Inside and three Outside, and the seam-free path's own `debug_assert`
        // ("classification must be consistent per solid") caught it.
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.6]),
            Point3::from_array([0.6, 0.6, 0.9]),
        );
        let r = boolean(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!(m.solids.get(r).cavities.is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.92).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn cut_with_a_hollow_operand_is_unsupported() {
        // The convex path would read this hollow box as convex (`solid_vertices` sees
        // the outer shell alone) and drop the void, returning a wrong volume that
        // `validate` accepts. Reject at the door instead.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        let cutter = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, hollow, cutter),
            tag::HOLLOW_OPERAND,
        );
    }

    #[test]
    fn cut_into_a_cavity_hits_the_tunnel_guard() {
        // `boolean` now rejects cavitied operands at the door, so this calls the seam
        // path directly — the `tunnel` guard is still load-bearing for direct callers
        // and for the day sub-unit 5 admits cavities, and an unverified safety guard
        // is its own kind of silent failure.
        //
        // The guard's only reachable path: `point_in_solid` counts cavity shells, so
        // the stub's vertical edges — running from below the hollow L up into its void
        // — classify Outside at both ends, while `pierced_face` (outer shell only)
        // sees a single crossing of the z=0 face.
        let (mut m, l, inner) = l_and_inner_box();
        let hollow = boolean(&mut m, BoolKind::Cut, l, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        let stub = m.add_cuboid(
            Point3::from_array([0.4, 0.4, -0.2]),
            Point3::from_array([0.6, 0.6, 0.5]),
        );
        assert_rejects(
            || overlap_fuse_cut(&mut m, BoolKind::Cut, hollow, stub),
            tag::TUNNEL,
        );
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

    /// Outer A = [0,3]³ (volume 27) with inner B = [1,2]³ (volume 1) strictly
    /// inside it — the containment fixture (returns `(model, outer, inner)`).
    fn nested_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        (m, a, b)
    }

    #[test]
    fn cut_of_two_cubes() {
        // A − B where A = [0,1]³, B = [0.5,1.5]³ ⇒ 1 − 0.125 = 0.875.
        let (mut m, a, b) = two_boxes();
        let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.875).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn fuse_of_disjoint_boxes_is_empty() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        assert_eq!(
            boolean(&mut m, BoolKind::Fuse, a, b),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn cut_containment_makes_a_cavity() {
        // A = [0,3]³ (27) with B = [1,2]³ (1) strictly inside ⇒ A − B is a
        // hollow solid: volume 26, an outer + one void shell (V16/E24/F12/S2).
        let (mut m, a, b) = nested_boxes();
        let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 26.0).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
        let reach = m.reachable();
        assert_eq!(reach.shells.len(), 2);
        assert_eq!(reach.faces.len(), 12);
        assert_eq!(reach.vertices.len(), 16);
        assert_eq!(reach.edges.len(), 24);
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn cut_containment_off_center_cavity() {
        // The inner box need not be concentric — any strictly-interior B works.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 3.5]),
        );
        let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (64.0 - 6.0)).abs() < 1e-12, "volume {vol}"); // 4³ − 1·2·3
        assert_eq!(m.solids.get(r).cavities.len(), 1);
    }

    #[test]
    fn fuse_containment_is_the_container() {
        // A ∪ B with B ⊂ A is just A (no cavity).
        let (mut m, a, b) = nested_boxes();
        let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
        let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - vol_a).abs() < 1e-12, "volume {vol}");
        assert!(m.solids.get(r).cavities.is_empty());
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn containment_symmetric_when_a_inside_b() {
        // Arguments swapped: A = inner ⊂ B = outer.
        let (mut m, outer, inner) = nested_boxes();
        let vol_outer = nacre_props::mass_props(&m, outer).unwrap().volume;
        // Cut(inner − outer): inner is wholly removed ⇒ empty.
        assert_eq!(
            boolean(&mut m, BoolKind::Cut, inner, outer),
            Err(BoolError::EmptyResult)
        );
        // Fuse(inner ∪ outer) = outer.
        let r = boolean(&mut m, BoolKind::Fuse, inner, outer).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - vol_outer).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn common_containment_is_the_inner_solid() {
        // Either argument order ⇒ the intersection is the inner solid (handled
        // by the existing half-space enumeration path, no cavity code).
        for swap in [false, true] {
            let (mut m, outer, inner) = nested_boxes();
            let vol_inner = nacre_props::mass_props(&m, inner).unwrap().volume;
            let (x, y) = if swap { (inner, outer) } else { (outer, inner) };
            let r = boolean(&mut m, BoolKind::Common, x, y).unwrap();
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            assert!((vol - vol_inner).abs() < 1e-9, "swap={swap} volume {vol}");
        }
    }

    #[test]
    fn cut_of_disjoint_is_a() {
        // A − B with B disjoint from A removes nothing ⇒ the result is A.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
        let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - vol_a).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]);
    }

    // ---- coincident-coplanar merge (M5-c5) ----

    fn stacked_cubes() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        (m, a, b)
    }

    #[test]
    fn fuse_stacked_cubes() {
        // A=[0,1]³ and B=[0,1]²×[1,2] share the z=1 face ⇒ merge into a 1×1×2 box
        // (10 faces: the 4 side pairs stay coplanar-adjacent, unmerged).
        let (mut m, a, b) = stacked_cubes();
        let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let reach = m.reachable();
        assert_eq!(reach.faces.len(), 10);
        assert_eq!(reach.vertices.len(), 12);
        assert_eq!(reach.edges.len(), 20);
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 2.0).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn common_stacked_cubes_is_empty() {
        let (mut m, a, b) = stacked_cubes();
        assert_eq!(
            boolean(&mut m, BoolKind::Common, a, b),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn cut_stacked_cubes_is_a() {
        let (mut m, a, b) = stacked_cubes();
        let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn coincident_merge_rejects_offset_footprint() {
        // Coplanar z=1 interface but B's footprint is offset ⇒ boundaries differ ⇒
        // partial 2D overlap, out of scope.
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        // Rejected by the combined-plane coplanar check, not by a merge-specific guard.
        assert_rejects(|| boolean(&mut m, BoolKind::Fuse, a, b), tag::COPLANAR_PAIR);
    }

    #[test]
    fn same_ground_overlap_is_unsupported() {
        // Two boxes sharing the z=0 ground with overlapping footprints: only
        // same-normal coplanar contact + 3D overlap ⇒ needs the general 2D
        // coplanar path (next unit) ⇒ Unsupported.
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.0]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        assert_rejects(|| boolean(&mut m, BoolKind::Fuse, a, b), tag::COPLANAR_PAIR);
    }

    #[test]
    fn fuse_rejects_poke_through_hole() {
        // A bar through A's interior pierces two A-faces mid-face (an interior
        // seam loop = a hole) — out of clean-seam coverage.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, a, bar),
            tag::POKE_THROUGH,
        );
    }

    proptest! {
        /// Diagonal corner overlaps (clean seam): fuse/cut volumes match the
        /// independent AABB formula (not nacre's own common).
        #[test]
        fn fuse_cut_diagonal_boxes_match_aabb(
            amin in prop::array::uniform3(-3.0f64..3.0),
            aext in prop::array::uniform3(1.0f64..3.0),
            t in prop::array::uniform3(0.15f64..0.6),
            s in prop::array::uniform3(0.3f64..2.0),
        ) {
            let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
            let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
            let bmax: [f64; 3] = std::array::from_fn(|i| amax[i] + s[i]);
            let ov: f64 = (0..3).map(|i| amax[i] - bmin[i]).product();
            let va: f64 = aext.iter().product();
            let vb: f64 = (0..3).map(|i| bmax[i] - bmin[i]).product();

            let build = || {
                let mut m = Model::new();
                let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
                let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
                (m, a, b)
            };

            let (mut m1, a1, b1) = build();
            let rf = boolean(&mut m1, BoolKind::Fuse, a1, b1);
            prop_assume!(rf.is_ok()); // skip rare coplanar/degenerate configs
            let rf = rf.unwrap();
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb - ov)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            let rc = boolean(&mut m2, BoolKind::Cut, a2, b2).unwrap();
            m2.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m2).is_empty());
            let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;
            prop_assert!((vc - (va - ov)).abs() <= 1e-9 * va, "cut {vc}");
        }

        /// Matched-footprint stacked boxes (coincident z-interface): fuse volume
        /// is the sum, cut is A, common is empty.
        #[test]
        fn stacked_boxes_merge_volumes(
            x0 in -3.0f64..3.0,
            y0 in -3.0f64..3.0,
            dx in 0.5f64..3.0,
            dy in 0.5f64..3.0,
            z0 in -3.0f64..3.0,
            h1 in 0.5f64..3.0,
            h2 in 0.5f64..3.0,
        ) {
            let (x1, y1) = (x0 + dx, y0 + dy);
            let (zm, z1) = (z0 + h1, z0 + h1 + h2);
            let build = || {
                let mut m = Model::new();
                let a = m.add_cuboid(Point3::from_array([x0, y0, z0]), Point3::from_array([x1, y1, zm]));
                let b = m.add_cuboid(Point3::from_array([x0, y0, zm]), Point3::from_array([x1, y1, z1]));
                (m, a, b)
            };
            let (va, vb) = (dx * dy * h1, dx * dy * h2);

            let (mut m1, a1, b1) = build();
            let rf = boolean(&mut m1, BoolKind::Fuse, a1, b1);
            prop_assume!(rf.is_ok());
            let rf = rf.unwrap();
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            prop_assert_eq!(
                boolean(&mut m2, BoolKind::Common, a2, b2),
                Err(BoolError::EmptyResult)
            );

            let (mut m3, a3, b3) = build();
            let rc = boolean(&mut m3, BoolKind::Cut, a3, b3).unwrap();
            let vc = nacre_props::mass_props(&m3, rc).unwrap().volume;
            prop_assert!((vc - va).abs() <= 1e-9 * va, "cut {vc}");
        }
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
        assert_rejects(
            || boolean(&mut m, BoolKind::Common, a, cyl),
            tag::CYLINDER_FACE,
        );
    }

    #[test]
    fn common_rejects_non_convex_input() {
        // An L-shaped prism (reflex edge) is not the intersection of its face
        // half-spaces. It never reaches `common`'s `is_convex` guard, though: the
        // box's corner grazes the L's boundary, so `boundaries_intersect` bails on a
        // degenerate contact first. Measured, not assumed.
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
        assert_rejects(
            || boolean(&mut m, BoolKind::Common, lsolid, b),
            tag::CONTACT_DEGENERATE,
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
        assert_rejects(
            || boolean(&mut m, BoolKind::Common, solid, b),
            tag::COPLANAR_PAIR,
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
