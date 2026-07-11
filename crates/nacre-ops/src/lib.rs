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
    /// A seam-free face whose rings do not agree about which side of the other solid
    /// they are on. A backstop with no firing test: were a rim vertex classified against
    /// the outer ring, the other boundary would separate the two rings and so would cut
    /// `f`, which is a seam. A cross-check between `point_in_solid` and the arrangement,
    /// not a defensive assert.
    pub const HOLE_CLASS_SPLIT: &str = "hole_class_split";
    /// An outer-shell edge used by other than two face loops. A backstop with no
    /// firing test: `validate` calls this `NonOpposedEdge` and every shell the
    /// operations build is manifold — but `boolean` never runs `validate` on its
    /// inputs, so a direct caller could still hand one in.
    pub const NON_MANIFOLD_EDGE: &str = "non_manifold_edge";
    /// A backstop with no firing test yet — the degeneracies that would produce an
    /// odd crossing count are expected to trip `VERTEX_ON_FACE_PLANE` or `POINT_ON_RING`
    /// first, since an edge that neither straddles a face's plane nor pierces it cleanly
    /// is named there. Recorded as unverified in design.md §9, alongside `fourplane`.
    pub const ARRANGEMENT_DEGENERATE: &str = "arrangement_degenerate";
    /// The arrangement and the vertex classification disagree about where the seam
    /// meets `∂f`, or the `seam` list and the arrangement disagree about whether it
    /// meets `f` at all.
    ///
    /// Not a defensive assert but a cross-check between two machineries, neither of
    /// which is allowed to be assumed right: `classof`'s ray casting (`point_in_solid`)
    /// and the arrangement's exact segment/face crossings. Cell 3c pinned that
    /// agreement in a test; here it is production.
    ///
    /// Cell 3e-3 gave it its real work. `∂f` is cut into runs at the seam's crossings,
    /// and the runs alternate kept/dropped because a crossing flips the class. That shape
    /// is one machine; `classof` at each original vertex is the other. A run with no vertex
    /// — the piece of an edge the other solid enters and leaves through — has only the
    /// alternation, so the check is what makes the rest of it trustworthy. It also holds
    /// the arc bookkeeping: every arc has two boundary ends whose classes oppose.
    pub const SEAM_COUNT_MISMATCH: &str = "seam_count_mismatch";
    /// A closed seam loop's edges disagree about which side its material lies on, or
    /// two of its nodes coincide, or two neighbours share no plane pair.
    ///
    /// Like `SEAM_COUNT_MISMATCH`, a cross-check rather than a defence: the exact
    /// predicate (`order_along`) and the orientation bookkeeping (`PlaneInfo::n_out`
    /// vs its plane's normal) must agree, edge by edge, and neither is assumed right.
    /// It survives release on purpose. Downstream, `validate` catches a wrong loop and
    /// `tessellate` catches a wrong *hole* — but a caller may run neither.
    ///
    /// It cannot catch a *globally* flipped loop: every edge would be wrong together.
    /// That is pinned before wiring, by a golden on `orient_seam_loop` itself.
    ///
    /// Unreachable today: "material on the left" is a global property, so a consistent
    /// loop makes every edge agree. Not dead code — relax `FOURPLANE` or cell 3c's
    /// node-identity argument and this is what speaks first.
    pub const LOOP_ORIENT_MISMATCH: &str = "loop_orient_mismatch";
    /// Three seam segments meeting at one node.
    ///
    /// Unreachable today, and *not* dead code. A node is a boundary node (two
    /// X-planes in its triple, degree 1) or an interior node (two Y-planes, degree
    /// 2 — the pierced Y-edge has exactly two faces); the only way to collapse two
    /// distinct nodes onto one triple is a duplicate third plane within one plane
    /// pair, which `FOURPLANE` rejects. Relax that guard (sub-unit 3e/3h) and a
    /// node can gain a third segment. This is the backstop for that day.
    pub const SEAM_BRANCH: &str = "seam_branch";
    pub const COPLANAR_PAIR: &str = "coplanar_pair";
    /// The result would be two solids, and `boolean` returns one handle. Reachable since
    /// cell 3e-3 let an edge thread the other solid: `Cut(rod, L)` leaves the rod's two
    /// ends on either side of the bar.
    pub const DISCONNECTED_RESULT: &str = "disconnected_result";
    /// The two machines that decide where the seam is disagree about *parity*: a segment
    /// crosses a closed surface an odd number of times exactly when its endpoints lie on
    /// opposite sides, and `point_in_solid`'s winding says one thing while
    /// `pierced_faces`' exact crossings say the other.
    ///
    /// The one way that happens is a cavitied operand, which is where the name comes from:
    /// the classifier counts the cavity shells, the crossing count scans the outer shell
    /// alone, so an edge reaching into a void reads `Outside` at both ends while piercing
    /// once. `boolean` rejects those at the door (`HOLLOW_OPERAND`); this stands for the
    /// tests that call `overlap_fuse_cut` directly, and for the day that door opens.
    ///
    /// It is *not* the out→in→out tunnel the name suggests. That edge crosses twice, the
    /// parity agrees, and since cell 3e-3 it earns two seam vertices and drills a hole.
    pub const TUNNEL: &str = "tunnel";
    pub const NO_ENTRY_FACE: &str = "no_entry_face";
    pub const THREE_PLANES: &str = "three_planes";
    pub const FOURPLANE: &str = "fourplane";
    pub const CYLINDER_FACE: &str = "cylinder_face";
    pub const DEGENERATE_FACE: &str = "degenerate_face";
    pub const DEGENERATE_NORMAL: &str = "degenerate_normal";
    pub const RAY_DEGENERATE: &str = "ray_degenerate";
    /// An endpoint of an edge lies exactly on a face's plane — and if both do, the edge lies
    /// in that plane. A tangential contact, not a crossing, and out of clean-seam coverage.
    pub const VERTEX_ON_FACE_PLANE: &str = "vertex_on_face_plane";
    /// A closed seam loop's winding and the face boundary's class disagree.
    ///
    /// A hole keeps its material outside, so walked material-left it runs clockwise about
    /// the face's outward normal; an island runs counter-clockwise. Whether the boundary is
    /// kept says the same thing, by ray casting. Three sources meet in that one equation and
    /// none is assumed right: the local sign rule (`n_out` bookkeeping and `order_along`),
    /// the ring's global winding (`det3` at a hull vertex, found by exact comparison), and
    /// `classof`. `keep` cancels — it flips the ring and `kept[0]` together — so this checks
    /// the geometry, not the operation.
    ///
    /// Unreachable today; the three agree on every fixture. Cell 3f-3 promotes the same
    /// equation into a *nesting* detector, where the winding classifies each loop and
    /// `kept[0]` only cross-checks depth zero.
    pub const LOOP_CLASS_MISMATCH: &str = "loop_class_mismatch";
    /// Every candidate ray from a loop's nodes has a ring node on its line.
    ///
    /// `point_in_ring` casts along `P ∩ Q_a` for a node's own plane `Q_a`; a ring node on
    /// that line makes the crossing parity ambiguous. Candidates are `2 · |loop|` lines and
    /// two directions, and half of them can be spoiled at once — `l_and_staple`'s loop and
    /// arc share both `y` planes, so only the `x` lines are clear there. Unfired today.
    pub const NO_CLEAR_RAY: &str = "no_clear_ray";
    /// A loop's node lies *on* the ring it is being tested against.
    ///
    /// `SeamPath::Closed` has claimed since cell 3f-1 that a closed seam loop never touches
    /// `∂f`. Nothing checked it. `point_in_ring` does, exactly: the ray's line meets an edge
    /// at `X`, and `X == v` strictly inside that edge means `v` is on the ring. Unfired.
    pub const POINT_ON_RING: &str = "point_on_ring";
    /// A face whose boundary never crosses the seam, yet the seam lies on its plane — the
    /// convex path only.
    ///
    pub const MISSING_SEAM: &str = "missing_seam";
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
    /// The imprint/pad/pocket profile is not strictly inside the target face's region
    /// (it crosses, sits outside, or overhangs the boundary, or touches a hole). A
    /// face-local imprint needs a clean inner-loop hole; profiles that reach past the face
    /// belong to a boolean pad/pocket (extrude + fuse/cut), not this path (cell
    /// imprint-containment).
    ProfileNotContainedInFace,
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
/// Drop the axis of the face normal's largest component to project a planar point to 2D.
/// An exact projection — a coordinate is discarded, never recomputed — so `orient2d` on the
/// result stays exact (cell imprint-containment).
fn planar_drop_axes(n: Vector3) -> (usize, usize) {
    let a = [n[0].abs(), n[1].abs(), n[2].abs()];
    if a[0] >= a[1] && a[0] >= a[2] {
        (1, 2)
    } else if a[1] >= a[2] {
        (0, 2)
    } else {
        (0, 1)
    }
}

fn proj2(p: Point3, (i, j): (usize, usize)) -> [f64; 2] {
    [p[i], p[j]]
}

/// `c` is within the axis-aligned bounding box of segment `ab` — paired with an exact
/// `orient2d == 0` collinearity test to decide on-segment.
fn in_bbox(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    c[0] >= a[0].min(b[0])
        && c[0] <= a[0].max(b[0])
        && c[1] >= a[1].min(b[1])
        && c[1] <= a[1].max(b[1])
}

/// Exact strict point-in-ring: `Some(true)` strictly inside, `Some(false)` strictly
/// outside, `None` on the boundary. Even-odd ray cast, `orient2d` for the crossing side.
fn point_in_ring2(p: [f64; 2], ring: &[[f64; 2]]) -> Option<bool> {
    let n = ring.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        if orient2d(a, b, p) == 0.0 && in_bbox(a, b, p) {
            return None;
        }
        if (a[1] > p[1]) != (b[1] > p[1]) && (orient2d(a, b, p) > 0.0) == (b[1] > a[1]) {
            inside = !inside;
        }
    }
    Some(inside)
}

/// Exact: do closed segments `p0p1` and `q0q1` meet at all (cross or merely touch)?
fn segments_meet(p0: [f64; 2], p1: [f64; 2], q0: [f64; 2], q1: [f64; 2]) -> bool {
    let (o1, o2) = (orient2d(p0, p1, q0), orient2d(p0, p1, q1));
    let (o3, o4) = (orient2d(q0, q1, p0), orient2d(q0, q1, p1));
    if o1 != 0.0 && o2 != 0.0 && o3 != 0.0 && o4 != 0.0 {
        return (o1 > 0.0) != (o2 > 0.0) && (o3 > 0.0) != (o4 > 0.0);
    }
    (o1 == 0.0 && in_bbox(p0, p1, q0))
        || (o2 == 0.0 && in_bbox(p0, p1, q1))
        || (o3 == 0.0 && in_bbox(q0, q1, p0))
        || (o4 == 0.0 && in_bbox(q0, q1, p1))
}

/// Whether `profile` is strictly inside the face region: inside `outer`, outside every
/// hole, and touching no boundary edge. The imprint contract — a profile that reaches past
/// the face makes an invalid inner loop (cell imprint-containment). All rings are exact 2D
/// projections onto the face plane (dominant axis dropped).
fn profile_strictly_in_region(
    profile: &[[f64; 2]],
    outer: &[[f64; 2]],
    holes: &[Vec<[f64; 2]>],
) -> bool {
    for &v in profile {
        if point_in_ring2(v, outer) != Some(true) {
            return false;
        }
        if holes.iter().any(|h| point_in_ring2(v, h) != Some(false)) {
            return false;
        }
    }
    let m = profile.len();
    let rings = std::iter::once(outer).chain(holes.iter().map(|h| h.as_slice()));
    for ring in rings {
        let rn = ring.len();
        for i in 0..m {
            let (p0, p1) = (profile[i], profile[(i + 1) % m]);
            if (0..rn).any(|j| segments_meet(p0, p1, ring[j], ring[(j + 1) % rn])) {
                return false;
            }
        }
    }
    true
}

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
    let inner_loops = f.inner.clone();
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

    // Strict containment: the profile must lie inside the face region (inside the outer
    // ring, outside every hole, touching no boundary). Otherwise the "hole" is not a clean
    // inner loop and the result is silently invalid — validate sees only topology, props
    // integrates the ring, and only tessellate's `NoEar` catches it (cell imprint-containment;
    // measured in design.md §10). A profile reaching past the face is a boolean pad/pocket.
    let drop = planar_drop_axes(n);
    let profile2: Vec<[f64; 2]> = base_pts.iter().map(|&p| proj2(p, drop)).collect();
    let outer2: Vec<[f64; 2]> = outer_pts.iter().map(|&p| proj2(p, drop)).collect();
    let holes2: Vec<Vec<[f64; 2]>> = inner_loops
        .iter()
        .map(|l| {
            l.half_edges
                .iter()
                .map(|he| proj2(model.vertices.get(he_start(model, *he)).point, drop))
                .collect()
        })
        .collect();
    if !profile_strictly_in_region(&profile2, &outer2, &holes2) {
        return Err(OpError::ProfileNotContainedInFace);
    }

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
    RayCross, orient2d, plane_plane, plane_side, planes_coplanar, ray_face_cross,
    three_plane_orient3d, three_planes,
};
#[cfg(test)]
use std::collections::BTreeSet;
use std::collections::{HashMap, HashSet};

/// A face's supporting plane plus the exact in/out data the seam path needs.
///
/// `three_plane_orient3d(.., tri[0], tri[1], tri[2])` returns `+1` when the
/// implicit point lies on **`tri`'s right-hand-normal side** — the convention is
/// tied to the triangle, never to `plane`. `n_out` happens to equal that RH normal
/// only because `tri` is taken outer-CCW; `plane.normal()` is the *surface's*
/// normal and may point inward on a `Reversed` face. Every sign test here reads
/// `n_out` (or `tri`), and none reads `plane.normal()`.
pub(crate) struct PlaneInfo {
    pub(crate) surf: Handle<Surface>,
    /// The face this plane came from. Distinguishes two coplanar faces that share one
    /// `Surface` (a Cut splits one face into disjoint pieces reusing its surface —
    /// cell coplanar-narrow), which `surf` alone collapses. `surf_ix` keys on this.
    pub(crate) face: Handle<Face>,
    pub(crate) plane: Plane,
    /// Three non-collinear outer-loop points, **ordered so their RH normal is outward**.
    /// The order need not follow the loop: at a reflex corner it is reversed.
    pub(crate) tri: [Point3; 3],
    /// Outward normal, `(tri[1]−tri[0])×(tri[2]−tri[0])` normalized — the single
    /// source of "outward" for both the in/out sign test and face ordering.
    pub(crate) n_out: Vector3,
    pub(crate) orient: Orientation,
}

/// Boolean of two live solids (design §8 M5, overview 불리언 전략 — 정직하게 거절).
/// **Coverage:** planar solids in general position — all three kinds go through one
/// exact seam path (`general_boolean`), the convex fast paths (`fuse_cut` cell (5b),
/// `common` cell 3g) retired. Anything else is rejected with [`BoolError`]. Transactional:
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
    // Cavitied operands are supported (cells (5c-in), (5c-in-2)): the seam front-end and
    // reconstruction walk all shells (outer + cavities) via `solid_shell_handles`, so a
    // void the cut misses is carried through, and a cut reaching into a void reconstructs
    // (the opened void merges with the outer shell).
    // Coincident-coplanar degeneracy (M5-c5): a clean matched-interface stack.
    if let Some(iface) = detect_coincident_interface(model, a, b) {
        return coincident_merge(model, kind, a, b, &iface);
    }
    // A boss sitting on a face inside its boundary — coplanar contact with a contained
    // footprint (cell coplanar-contact-boss). Fuse only for now; Cut/Common fall through.
    if kind == BoolKind::Fuse {
        if let Some(cc) = detect_contained_contact(model, a, b) {
            return contained_boss_fuse(model, &cc);
        }
    }
    // One general exact path for every kind. Cell (5b) retired the convex `fuse_cut`
    // (it decided in/out with a 1e-9 tolerance) and cell 3g the convex `common` (its
    // `order_ccw` was the last `atan2` enumeration) — the seam path answers all three
    // exactly, so there is no fast path left to gate.
    general_boolean(model, kind, a, b)
}

/// The general boolean entry point for all three kinds since cells (5b) and 3g deleted the
/// convex `fuse_cut` and `common`. If the boundaries cross it hands off to
/// [`overlap_fuse_cut`]; otherwise one solid contains the other or
/// they are disjoint, classified by exact [`point_in_solid`] and assembled by
/// [`contained_result`].
fn general_boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    // Exact containment reads a face's rings as three-plane triples, and an imprinted face
    // cannot give them: its hole rim's two neighbours are *coplanar* (the holed face and the
    // region face cut from it), so the rim's vertices have no triple and the line `P ∩ R` no
    // direction. That is a coplanar pair sharing an *edge*; `solid_has_coplanar_neighbour_edge`
    // rejects it while passing disjoint coplanar faces (a slot's two top strips), so
    // multi-feature parts chain (cell coplanar-narrow). A coplanar pair *across* the two
    // solids is harmless here, since a face's rings only ever name its own solid's planes.
    let mut planes = collect_planes(model, a)?;
    let na = planes.len();
    planes.extend(collect_planes(model, b)?);
    if solid_has_coplanar_neighbour_edge(model, a, &planes[..na])
        || solid_has_coplanar_neighbour_edge(model, b, &planes[na..])
    {
        return Err(reject(tag::COPLANAR_PAIR));
    }
    let surf_ix: HashMap<Handle<Face>, usize> = planes
        .iter()
        .enumerate()
        .map(|(i, pi)| (pi.face, i))
        .collect();

    if boundaries_intersect(model, a, b, &planes, &surf_ix)? {
        // A genuine seam — all three kinds arrange the same faces, differing only in the
        // keep/flip table (cell 3g opened `Common`).
        return overlap_fuse_cut(model, kind, a, b);
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

/// A solid's outer-shell faces, each as its combined-plane index and its rings in triple
/// form — the input [`pierced_faces`] tests against, built once per solid.
type FaceRings = Vec<(usize, Vec<Vec<[usize; 3]>>)>;

fn solid_face_rings(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &arrange::EdgePlanes,
) -> Result<FaceRings, BoolError> {
    // All shells, outer + cavities: `pierced_faces` must see the void walls too, so
    // its crossing parity matches `point_in_solid`'s all-shell winding — the old
    // outer-only asymmetry was what `HOLLOW_OPERAND` stood in for (cell (5c-in)).
    let mut out = FaceRings::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let q = surf_ix[&fh];
            out.push((q, arrange::face_rings(model, fh, q, inc)?));
        }
    }
    Ok(out)
}

/// Every face of `rings` that the edge on planes `pair`, running `p0 → p1`, pierces — as
/// combined-plane indices, in shell order. The exact analogue of the retired convex
/// `enter_face` (cell (5b)).
///
/// A segment meets a plane at most once, so no plane appears twice — and two faces on one
/// plane would already have been rejected as a `coplanar_pair`. So each hit is a distinct
/// seam vertex.
fn pierced_faces(
    planes: &[PlaneInfo],
    pair: [usize; 2],
    p0: Point3,
    p1: Point3,
    rings: &FaceRings,
) -> Result<Vec<usize>, BoolError> {
    let mut hits = Vec::new();
    for (q, r) in rings {
        if arrange::edge_crosses_face(planes, pair, p0, p1, *q, r)? {
            hits.push(*q);
        }
    }
    Ok(hits)
}

/// `Fuse`/`Cut` of two solids whose boundaries cross. Since cell (5b) this is the only
/// `Fuse`/`Cut` path — the convex `fuse_cut` is gone. It (a) classifies each original
/// vertex with exact [`point_in_solid`], (b) finds seam entries with [`pierced_faces`],
/// and (c) reconstructs faces from the arrangement. Multiple chords (3e-2), poke-through
/// holes and through-drilling (3e-3), and holed operands (3f-5) are handled; non-convex
/// `Common` is `COMMON_OVERLAP` (cell 3g). A result that falls into disconnected pieces is
/// `DISCONNECTED_RESULT`.
fn overlap_fuse_cut(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    let mut planes = collect_planes(model, a)?;
    let na = planes.len();
    planes.extend(collect_planes(model, b)?);
    // Reject a coplanar pair that shares an edge within either operand (imprint rim, Fuse
    // flat-edge) — disjoint coplanar faces (a slot) pass. A coplanar pair *across* the two
    // operands is a tangential face-to-face contact the seam machinery does not handle, so
    // it stays rejected (cell coplanar-narrow).
    if solid_has_coplanar_neighbour_edge(model, a, &planes[..na])
        || solid_has_coplanar_neighbour_edge(model, b, &planes[na..])
        || cross_coplanar(&planes[..na], &planes[na..])
    {
        return Err(reject(tag::COPLANAR_PAIR));
    }

    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        surf_ix.insert(pi.face, i);
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
    let edges_a = edge_incidence(model, a, &surf_ix)?;
    let edges_b = edge_incidence(model, b, &surf_ix)?;
    let inc_a = arrange::edge_planes(model, a, &surf_ix)?;
    let inc_b = arrange::edge_planes(model, b, &surf_ix)?;
    // Each solid's faces as triple rings, once. `edge_crosses_face` reads them per edge.
    let rings_a = solid_face_rings(model, a, &surf_ix, &inc_a)?;
    let rings_b = solid_face_rings(model, b, &surf_ix, &inc_b)?;
    let mut seam: Vec<SeamVertex> = Vec::new();
    let mut seam_ix: HashMap<[usize; 3], usize> = HashMap::new();
    for (edges, other) in [(&edges_a, &rings_b), (&edges_b, &rings_a)] {
        for &(_, bounds, inc) in edges {
            let [v0, v1] = bounds;
            let (p0, p1) = (model.vertices.get(v0).point, model.vertices.get(v1).point);
            let (s0, s1) = (classof[&v0], classof[&v1]);
            let hits = pierced_faces(&planes, inc, p0, p1, other)?;
            let straddles = s0 != s1;
            if straddles && hits.is_empty() {
                return Err(reject(tag::NO_ENTRY_FACE)); // an unfired backstop, now
            }
            // A segment crosses a closed surface an odd number of times exactly when its
            // endpoints lie on opposite sides of it. Two machines say those two things —
            // `point_in_solid`'s winding and `pierced_faces`' exact crossings — and they
            // must agree. The one way they can differ is a cavitied operand: the classifier
            // counts the cavity shells, the crossing count scans the outer shell alone, so
            // an edge reaching into a void reads `Outside` at both ends while piercing once.
            // `boolean` rejects those at the door (`tag::HOLLOW_OPERAND`); this stays for the
            // tests that call `overlap_fuse_cut` directly, and for the day that door opens.
            //
            // An edge threading the other solid — out, in, out — is no longer a tunnel: it
            // crosses twice, the parity agrees, and it earns two seam vertices. That is what
            // cell 3e-3 opened, and it is why the guard is a parity check and not a count.
            if (hits.len() % 2 == 1) != straddles {
                return Err(reject(tag::TUNNEL));
            }
            let [e0, e1] = inc;
            for entry in hits {
                let point =
                    three_planes(&planes[e0].plane, &planes[e1].plane, &planes[entry].plane)
                        .ok_or_else(|| reject(tag::THREE_PLANES))?;
                // Exact 4-plane-concurrency guard: the seam vertex must not lie on any
                // *other* plane (a degenerate 4-plane meet). Unlike the convex path we do
                // NOT require it inside every half-space (`== -1`): a seam vertex of a
                // non-convex solid can be outside some face's half-space yet on the
                // boundary — `pierced_faces` already certified it is on the entry face.
                for (m, pm) in planes.iter().enumerate() {
                    if m == e0 || m == e1 || m == entry {
                        continue;
                    }
                    // A coplanar twin of one of the triple's planes (two disjoint faces
                    // sharing a Surface — cell coplanar-narrow) is not a genuine 4th plane:
                    // the vertex lies on it only because it *is* one of the triple's planes.
                    if planes_coplanar(&pm.plane, &planes[e0].plane)
                        || planes_coplanar(&pm.plane, &planes[e1].plane)
                        || planes_coplanar(&pm.plane, &planes[entry].plane)
                    {
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
    }

    let (keep_a, keep_b, flip_b) = match kind {
        BoolKind::Fuse => (Side::Outside, Side::Outside, false),
        BoolKind::Cut => (Side::Outside, Side::Inside, true),
        // A∩B: keep each solid's material *inside* the other, neither shell flipped —
        // the De Morgan dual of Fuse. Cell 3g.
        BoolKind::Common => (Side::Inside, Side::Inside, false),
    };
    if seam.is_empty() {
        return contained_result(model, kind, a, b, &classof);
    }

    // `edge_seam` is gone from this path: the arrangement carries each boundary
    // crossing on the edge it was produced on (`SeamSegment::on_edge`), so nothing has
    // to key a splice by `Handle<Edge>` — the map that could hold only one crossing per
    // edge no longer exists here.
    let mut faces: Vec<LocalFace> = Vec::new();
    for (solid, other, inc_f, inc_o, keep, flip) in [
        (a, b, &inc_a, &inc_b, keep_a, false),
        (b, a, &inc_b, &inc_a, keep_b, flip_b),
    ] {
        for sh in solid_shell_handles(model, solid) {
            for &fh in &model.shells.get(sh).faces {
                let pidx = surf_ix[&fh];
                faces.extend(reconstruct_face_paths(
                    model, fh, other, pidx, keep, flip, &classof, &seam, &seam_ix, &planes,
                    &surf_ix, inc_f, inc_o,
                )?);
            }
        }
    }
    // Output faces, not surviving input faces: one input face may split into several.
    if faces.len() < 4 {
        return Err(BoolError::EmptyResult);
    }
    assemble_fuse_cut(model, a, b, &planes, &seam, &faces)
}

/// Connected components of the reconstructed faces by shared `Node` — the same identity
/// `assemble_fuse_cut` welds result vertices by. Returns a component label per face (dense
/// `0..n` in order of first appearance, for replay determinism) and the component count. An
/// enclosed void is its own component: its boundary shares no vertex with the outer.
fn face_components(faces: &[LocalFace]) -> (Vec<usize>, usize) {
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    let mut owner: HashMap<Node, usize> = HashMap::new();
    for (i, lf) in faces.iter().enumerate() {
        for &nd in lf.loop_nodes.iter().chain(lf.inner.iter().flatten()) {
            match owner.entry(nd) {
                std::collections::hash_map::Entry::Occupied(e) => {
                    let (ra, rb) = (find(&mut parent, *e.get()), find(&mut parent, i));
                    parent[ra] = rb;
                }
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(i);
                }
            }
        }
    }
    let mut label: HashMap<usize, usize> = HashMap::new();
    let mut labels = Vec::with_capacity(faces.len());
    let mut n = 0;
    for i in 0..faces.len() {
        let root = find(&mut parent, i);
        let next = label.len();
        let l = *label.entry(root).or_insert(next);
        labels.push(l);
        n = n.max(l + 1);
    }
    (labels, n)
}

/// Area and area-weighted centroid of a planar loop (straight polygon), mirroring
/// `nacre-props::polygon_area_centroid`: a signed triangle fan from the first vertex, exact
/// for concave rings.
fn loop_area_centroid(model: &Model, l: &Loop) -> (f64, Point3) {
    let pts: Vec<Point3> = l
        .half_edges
        .iter()
        .map(|&he| model.vertices.get(he_start(model, he)).point)
        .collect();
    let base = pts[0];
    let mut area_vec = Vector3::zero();
    for w in pts[1..].windows(2) {
        area_vec += (w[0] - base).cross(w[1] - base);
    }
    let unit = area_vec.normalize().unwrap_or(Vector3::zero());
    let (mut weighted, mut weight) = (Vector3::zero(), 0.0);
    for w in pts[1..].windows(2) {
        let signed = (w[0] - base).cross(w[1] - base).dot(unit);
        let c_rel = ((w[0] - base) + (w[1] - base)) * (1.0 / 3.0);
        weighted += c_rel * signed;
        weight += signed;
    }
    let area = 0.5 * area_vec.norm();
    let centroid = if weight != 0.0 {
        base + weighted * (1.0 / weight)
    } else {
        base
    };
    (area, centroid)
}

/// The signed volume flux `∮ (r − R)·n̂ dA` (three times the signed volume) over a set of
/// faces forming one closed shell — **positive** for an outward, material-enclosing shell,
/// **negative** for an inward void shell. That sign is exactly "solid region vs enclosed
/// void": a genuinely-separate solid piece holds material inside (faces outward, positive),
/// a void holds material outside (faces inward, negative). Planar faces only (M5). Mirrors
/// `nacre-props::face_contribution`, reading the materialized `Face` so orientation/flip is
/// already baked in.
fn shell_signed_flux(model: &Model, faces: &[Handle<Face>]) -> f64 {
    let f0 = model.faces.get(faces[0]);
    let reference = model
        .vertices
        .get(he_start(model, f0.outer.half_edges[0]))
        .point;
    let mut flux = 0.0;
    for &fh in faces {
        let face = model.faces.get(fh);
        let sign = match face.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        let Surface::Plane(plane) = model.surfaces.get(face.surface) else {
            continue;
        };
        let normal = plane.normal() * sign;
        let (area, centroid) = loop_area_centroid(model, &face.outer);
        flux += normal.dot(centroid - reference) * area;
        for hole in &face.inner {
            let (a, c) = loop_area_centroid(model, hole);
            flux -= normal.dot(c - reference) * a;
        }
    }
    flux
}

/// The supporting planes of a solid's outer shell. `Unsupported` if any face is
/// non-planar or lacks three non-collinear loop points.
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<PlaneInfo>, BoolError> {
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
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
                face: fh,
                plane,
                tri,
                n_out,
                orient: face.orientation,
            });
        }
    }
    Ok(out)
}

/// All shells of a solid — outer first, then cavities. The boolean seam
/// front-end walks these so a cavitied operand's void walls are seen (cell
/// (5c-in)); a non-hollow solid yields just its outer shell, unchanged.
pub(crate) fn solid_shell_handles(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Shell>> {
    let s = model.solids.get(solid);
    std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .collect()
}

/// Three non-collinear points of a face's outer loop, ordered so their right-hand
/// normal points **out** of the solid.
fn outer_tri(model: &Model, face: &Face) -> Option<[Point3; 3]> {
    let pts: Vec<Point3> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| model.vertices.get(he_start(model, he)).point)
        .collect();
    let n = pts.len();
    // The turn at one corner does not know which way the ring winds. Every b-rep loop is
    // CCW about its face's outward normal, but at a *reflex* corner the local turn
    // opposes the global winding, so three consecutive points can hand back an inward
    // normal. The Newell sum has no single corner to be fooled by.
    let newell = (0..n).fold(Vector3::zero(), |acc, i| {
        acc + (pts[i] - pts[0]).cross(pts[(i + 1) % n] - pts[0])
    });
    let [a, b, c] = (0..n).find_map(|i| {
        let (a, b, c) = (pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        ((b - a).cross(c - a).norm() > 0.0).then_some([a, b, c])
    })?;
    Some(if (b - a).cross(c - a).dot(newell) < 0.0 {
        [a, c, b]
    } else {
        [a, b, c]
    })
}

/// The distinct outer-shell vertices of a solid.
///
/// Convexity: every vertex is on the inner side of (or on) every face plane,
/// judged by the vertex's **definition**, not its f64 coordinate cache (cell
/// 5d-3). `pi.tri` is oriented so its right-hand normal is `n_out`, so `side ≤ 0`
/// is exactly "inner side or on the plane", and no tolerance enters:
/// - a `Constructed` vertex's coordinate *is* the truth → exact `plane_side`;
/// - a `Discovered` vertex is the meet of three planes → judge that implicit
///   point against the face directly with `three_plane_orient3d`, never reading
///   its rounded `point` cache.
///
/// A boolean *result* fed back as an operand carries `Discovered` corners (the
/// first such operand is cell 5d-3's fixture), so this branch is live. In M5 the
/// two paths agree (axis-aligned meets are f64-exact); the point is to keep the
/// convexity decision off the coordinate cache, as the sweep requires.
fn is_convex(model: &Model, planes: &[PlaneInfo], vhs: &[Handle<Vertex>]) -> bool {
    let plane_of = |s: Handle<Surface>| match model.surfaces.get(s) {
        Surface::Plane(p) => Some(p),
        _ => None,
    };
    planes.iter().all(|pi| {
        vhs.iter().all(|&vh| {
            let v = model.vertices.get(vh);
            let side = match v.origin {
                Origin::Discovered {
                    definition: VertexDef::ThreePlane([s0, s1, s2]),
                    ..
                } => match (plane_of(s0), plane_of(s1), plane_of(s2)) {
                    (Some(p0), Some(p1), Some(p2)) => {
                        three_plane_orient3d(p0, p1, p2, pi.tri[0], pi.tri[1], pi.tri[2])
                    }
                    // M5-impossible (a three-plane def whose surfaces aren't all
                    // planes); fall back to the cached point defensively.
                    _ => plane_side(pi.tri, v.point),
                },
                Origin::Constructed => plane_side(pi.tri, v.point),
            };
            side <= 0
        })
    })
}

/// Whether any two planes in the set are the same plane, by the exact rank-1
/// [`planes_coplanar`] test (design §3 (5d)-2). Superseded in production by
/// [`solid_has_coplanar_neighbour_edge`] (cell coplanar-narrow); kept as a test predicate
/// that characterises the coplanar content of a fixture.
#[cfg(test)]
fn has_coplanar_pair(planes: &[PlaneInfo]) -> bool {
    (0..planes.len())
        .any(|i| (i + 1..planes.len()).any(|j| planes_coplanar(&planes[i].plane, &planes[j].plane)))
}

/// Whether two faces of `solid` that meet at an edge are coplanar — an imprint rim (the
/// holed face and the region cut from it) or a Fuse flat-edge, both defeature artifacts the
/// seam machinery cannot arrange. Unlike [`has_coplanar_pair`] this passes *disjoint*
/// coplanar faces: a slot's two top strips share no edge, so multi-feature parts chain
/// (cell coplanar-narrow). `planes` are the faces of `solid` alone (per-operand).
fn solid_has_coplanar_neighbour_edge(
    model: &Model,
    solid: Handle<Solid>,
    planes: &[PlaneInfo],
) -> bool {
    let plane_of: HashMap<Handle<Face>, &Plane> =
        planes.iter().map(|pi| (pi.face, &pi.plane)).collect();
    let mut edge_faces: HashMap<Handle<Edge>, Vec<Handle<Face>>> = HashMap::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            for he in face_half_edges(model.faces.get(fh)) {
                edge_faces.entry(he.edge).or_default().push(fh);
            }
        }
    }
    edge_faces
        .values()
        .any(|fs| fs.len() == 2 && planes_coplanar(plane_of[&fs[0]], plane_of[&fs[1]]))
}

/// Whether a face of `a` is coplanar with a face of `b` — a tangential face-to-face contact
/// across the two operands, which the seam machinery does not handle (cell coplanar-narrow).
fn cross_coplanar(a: &[PlaneInfo], b: &[PlaneInfo]) -> bool {
    a.iter()
        .any(|pa| b.iter().any(|pb| planes_coplanar(&pa.plane, &pb.plane)))
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

fn unordered(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
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
    /// Hole rings, each already wound so the kept material stays on its left
    /// (`arrange::orient_seam_loop`). Only the non-convex path ever fills this.
    inner: Vec<Vec<Node>>,
    flip: bool,
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

/// Whether the boundaries of `a` and `b` actually cross — an edge of one pierces a face of
/// the other (either direction). The same question [`pierced_faces`] answers, so it is the
/// same function, and one degeneracy policy covers both.
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
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<bool, BoolError> {
    for (edge_solid, face_solid) in [(a, b), (b, a)] {
        let inc_f = arrange::edge_planes(model, face_solid, surf_ix)?;
        let rings = solid_face_rings(model, face_solid, surf_ix, &inc_f)?;
        for (_, bounds, pair) in edge_incidence(model, edge_solid, surf_ix)? {
            let (p0, p1) = (
                model.vertices.get(bounds[0]).point,
                model.vertices.get(bounds[1]).point,
            );
            if !pierced_faces(planes, pair, p0, p1, &rings)?.is_empty() {
                return Ok(true); // a genuine seam
            }
        }
    }
    Ok(false)
}

/// Distinct outer-shell vertex handles of a solid, in shell→face→loop order.
fn solid_vertex_handles(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Vertex>> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            for &he in &model.faces.get(fh).outer.half_edges {
                let vh = he_start(model, he);
                if seen.insert(vh) {
                    out.push(vh);
                }
            }
        }
    }
    out
}

/// Each outer-shell edge with its bound vertices and the two combined-plane
/// indices of its adjacent faces, in first-seen (deterministic) order.
///
/// Every loop of every face is walked, holes included: a hole-ring edge is used
/// once by the holed face's inner loop and once by the neighbouring wall's outer
/// loop, so it too has exactly two incident faces. Walking `outer` before `inner`
/// on each face leaves the order of a hole-free solid untouched.
///
/// The pair is returned as `[usize; 2]`, so no caller can index a third slot: an
/// edge with any other incidence count is a non-manifold shell and rejects here.
#[allow(clippy::type_complexity)]
pub(crate) fn edge_incidence(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<Vec<(Handle<Edge>, [Handle<Vertex>; 2], [usize; 2])>, BoolError> {
    let mut order: Vec<Handle<Edge>> = Vec::new();
    let mut map: HashMap<Handle<Edge>, ([Handle<Vertex>; 2], Vec<usize>)> = HashMap::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let pidx = surf_ix[&fh];
            for he in face_half_edges(face) {
                let bounds = model.edges.get(he.edge).bounds.expect("bounded");
                let entry = map.entry(he.edge).or_insert_with(|| {
                    order.push(he.edge);
                    (bounds, Vec::new())
                });
                entry.1.push(pidx);
            }
        }
    }
    order
        .into_iter()
        .map(|e| {
            let (b, p) = map.remove(&e).unwrap();
            match p[..] {
                [x, y] => Ok((e, b, [x, y])),
                // `validate` would call this `NonOpposedEdge`, but `boolean` never runs
                // `validate` on its inputs, so the guard stays. No firing test.
                _ => Err(reject(tag::NON_MANIFOLD_EDGE)),
            }
        })
        .collect()
}

/// Every half-edge of a face: its outer loop first, then each hole ring in order.
pub(crate) fn face_half_edges(face: &Face) -> impl Iterator<Item = &HalfEdge> {
    face.outer
        .half_edges
        .iter()
        .chain(face.inner.iter().flat_map(|l| l.half_edges.iter()))
}

/// The result when no edge crosses the other solid's boundary: one solid
/// contains the other, or they are disjoint. Direction is read from `classof`
/// (already `Err`-free — a boundary vertex would have failed `point_in_solid`
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

/// The same question for a face's whole set of loops: all holes, or all islands?
///
/// **Precondition: the face carries no open arc.** Then `∂f` is one class throughout, the
/// region holding it has depth zero, and a loop's winding is fixed by the parity of its
/// depth. Containment shifts depth by one, so loops that all wind the same way cannot
/// contain one another — and being flat, they all take the class `kept0` dictates.
///
/// With an arc present the argument dies: `∂f` is no longer one class, and a hole in the
/// kept region beside an island in the dropped one wind oppositely without nesting. Calling
/// this there would raise a false `nested_loops`. Cell 3f-4 pays for that with a real
/// containment test.
///
/// Nothing more can be checked in the mixed case. The outermost loops sit at depth zero and
/// wind the way `kept0` dictates — but "mixed" already means both signs are present, so
/// that sign is there by definition and the cross-check would be vacuous.
/// The winding a loop must have, given where it turned out to be.
///
/// A hole keeps its material outside itself and so runs clockwise; an island runs the other
/// way. `is_hole` comes from containment, the winding from an exact turn at a hull vertex,
/// and the ring's direction from the local sign rule. Three sources, none assumed right.
fn check_loop_class(is_hole: bool, winding: i8) -> Result<(), BoolError> {
    if winding == if is_hole { -1 } else { 1 } {
        Ok(())
    } else {
        Err(reject(tag::LOOP_CLASS_MISMATCH))
    }
}

/// Each loop's place in the containment forest, as pure geometry.
///
/// `containers` are the loops whose interior holds this loop's representative node; their count
/// is the loop's nesting depth. `region` is the one kept region holding it, if any — two would
/// mean the regions overlap, which `stitch_cycles` forbids, so `seam_count_mismatch`. Neither
/// class nor winding is read here; `classify_nesting` turns this into hole/island placement.
struct Nesting {
    containers: Vec<usize>,
    region: Option<usize>,
}

fn nest_loops(
    planes: &[PlaneInfo],
    p: usize,
    regions: &[Vec<[usize; 3]>],
    loops: &[Vec<[usize; 3]>],
) -> Result<Vec<Nesting>, BoolError> {
    let mut out = Vec::with_capacity(loops.len());
    for (i, l) in loops.iter().enumerate() {
        let mut containers = Vec::new();
        for (j, other) in loops.iter().enumerate() {
            if i != j && arrange::point_in_ring(planes, p, l[0], other)? {
                containers.push(j);
            }
        }
        let mut region = None;
        for (k, r) in regions.iter().enumerate() {
            if arrange::point_in_ring(planes, p, l[0], r)? {
                if region.is_some() {
                    return Err(reject(tag::SEAM_COUNT_MISMATCH));
                }
                region = Some(k);
            }
        }
        out.push(Nesting { containers, region });
    }
    Ok(out)
}

/// The material face a loop belongs to: a kept `∂f` region, or a seam island (an index into
/// the loop list, the island's own place).
#[derive(Clone, Copy, PartialEq, Debug)]
enum Owner {
    Region(usize),
    Island(usize),
}

/// What becomes of a loop: it is a material island in its own right, a hole of some owner, or
/// void that is discarded.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Placement {
    Island,
    Hole(Owner),
    Dropped,
}

/// Turn the containment forest into hole/island placement — pure combinatorics, no geometry.
///
/// A loop's interior is material iff its region's class flipped once per containing loop, since
/// material keeps to the left and crossing a loop reverses it: `is_island = base_kept XOR
/// (depth even)`, `depth` the number of containers. So a seam loop alternates with nesting
/// depth — a hole in a kept region, an island in that hole, a hole in that island. `f`'s own
/// holes (`i >= n_rings`) are never reclassified: they are stored clockwise and can only be
/// holes, whatever their depth (feeding them to the parity would call an unowned one an island
/// and reject a shape that is merely discarded). `check_loop_class` is the third, independent
/// source — the exact winding must agree with the parity.
///
/// A hole belongs to the innermost material face containing it: the deepest island among its
/// containers, else its region, else nothing — void in void, discarded as a hole in a dropped
/// `∂f` region always was. An island is a material face of its own.
fn classify_nesting(
    nest: &[Nesting],
    n_rings: usize,
    windings: &[i8],
) -> Result<Vec<Placement>, BoolError> {
    let is_island: Vec<bool> = nest
        .iter()
        .enumerate()
        .map(|(i, n)| i < n_rings && (n.region.is_some() ^ (n.containers.len() % 2 == 0)))
        .collect();
    for (i, &island) in is_island.iter().enumerate() {
        check_loop_class(!island, windings[i])?;
    }
    let place = is_island
        .iter()
        .enumerate()
        .map(|(i, &island)| {
            if island {
                return Placement::Island;
            }
            let inner_island = nest[i]
                .containers
                .iter()
                .copied()
                .filter(|&c| is_island[c])
                .max_by_key(|&c| nest[c].containers.len());
            match inner_island {
                Some(c) => Placement::Hole(Owner::Island(c)),
                None => match nest[i].region {
                    Some(k) => Placement::Hole(Owner::Region(k)),
                    None => Placement::Dropped,
                },
            }
        })
        .collect();
    Ok(place)
}

/// Reconstruct a face's kept portion. `None` if the face is dropped.
///
/// The seam sub-path comes from the arrangement (`arrange::seam_paths_on`), which walks
/// the true arc order, so
/// nothing here assumes the arc is monotone along its chord — on a reflex bend it is
/// not, and a bend can project outside the chord's endpoints entirely.
///
/// The arrangement is a source of *combinatorics*, never of geometry: the points that
/// reach the model stay `seam[..]`'s, with their tolerance. `SeamEnd::point` is a cache.
///
/// **This function reads no coordinate.** `touches_seam` compares plane triples, `bnd`
/// compares handles, the splice emits handles and triples. A folded arc needs no
/// turn test because the arc is a component of `∂other ∩ P ∩ f`, a 1-manifold: it is a
/// simple curve meeting `∂f` only at its two ends, so `kept run + arc` is a simple
/// polygon whichever way it turns. Checking that would take an axis projection and
/// `orient2d` — the 2D machinery this sub-unit has done without since cell 3b.
#[allow(clippy::too_many_arguments)]
fn reconstruct_face_paths(
    model: &Model,
    fh: Handle<Face>,
    other: Handle<Solid>,
    plane_idx: usize,
    keep: Side,
    flip: bool,
    classof: &HashMap<Handle<Vertex>, Side>,
    seam: &[SeamVertex],
    seam_ix: &HashMap<[usize; 3], usize>,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_f: &arrange::EdgePlanes,
    inc_o: &arrange::EdgePlanes,
) -> Result<Vec<LocalFace>, BoolError> {
    let face = model.faces.get(fh);
    let hes = &face.outer.half_edges;
    let n = hes.len();
    let verts: Vec<Handle<Vertex>> = hes.iter().map(|&he| he_start(model, he)).collect();
    let kept: Vec<bool> = verts.iter().map(|v| classof[v] == keep).collect();
    let transitions: Vec<usize> = (0..n).filter(|&i| kept[i] != kept[(i + 1) % n]).collect();

    // A dropped face contributes nothing, so `whole` returns zero faces or one — the
    // `Vec` a split face will need (cell 3e-2), not yet used for more than that. A
    // face the seam misses keeps its holes exactly as it had them; `flip` reverses
    // every ring alike, so a hole of a `Cut`'s inside-B piece stays a hole.
    let orig_ring = |l: &Loop| -> Vec<Node> {
        l.half_edges
            .iter()
            .map(|&he| Node::Orig(he_start(model, he)))
            .collect()
    };
    let whole = || {
        Vec::from_iter(kept[0].then(|| LocalFace {
            plane_idx,
            loop_nodes: verts.iter().map(|&v| Node::Orig(v)).collect(),
            inner: face.inner.iter().map(&orig_ring).collect(),
            flip,
        }))
    };

    // The seam misses this face: today's `transitions == 0` branch, unchanged. Running
    // the arrangement here would be worse than useless — it calls
    // `segment_crosses_face`, which can honestly reject a grazing contact on a face
    // that has nothing to do with the seam, and that would rewrite a reject tag.
    let touches_seam = seam.iter().any(|s| s.triple.contains(&plane_idx));
    if transitions.is_empty() && !touches_seam {
        // The whole face is kept or dropped, holes with it, so every ring must classify
        // alike. Were a rim vertex to disagree, the other boundary would separate the
        // rings and hence cut `f` — a seam, which we just established is not there.
        if face
            .inner
            .iter()
            .flat_map(|l| &l.half_edges)
            .any(|&he| (classof[&he_start(model, he)] == keep) != kept[0])
        {
            return Err(reject(tag::HOLE_CLASS_SPLIT));
        }
        return Ok(whole());
    }

    let paths = arrange::seam_paths_on(model, fh, other, planes, surf_ix, inc_f, inc_o)?;
    let opens: Vec<&arrange::SeamPath> = paths
        .iter()
        .filter(|p| matches!(p, arrange::SeamPath::Open(_)))
        .collect();
    let n_closed = paths.len() - opens.len();

    // `∂f` is the outer loop plus every hole rim the seam crosses. A boundary node riding a
    // rim edge is what takes `∂f` beyond one ring; `crossed` names those holes and is the one
    // source of truth for it — the same set joins the boundary here and is withheld from
    // `place_loops` below (a crossed hole is absorbed into the outer boundary, not placed).
    let ring_of = |e: Handle<Edge>| -> Option<usize> {
        // `None` is the outer loop; `Some(k)` is hole `k`.
        if face.outer.half_edges.iter().any(|he| he.edge == e) {
            None
        } else {
            face.inner
                .iter()
                .position(|l| l.half_edges.iter().any(|he| he.edge == e))
        }
    };
    let mut crossed: Vec<usize> = Vec::new();
    for nd in paths.iter().flat_map(|p| p.nodes()) {
        if let Some(k) = nd.on_edge.and_then(ring_of) {
            if !crossed.contains(&k) {
                crossed.push(k);
            }
        }
    }
    crossed.sort_unstable(); // `f.inner` order, so the boundary layout is deterministic

    // The boundary loops — outer, then each crossed rim — and the per-vertex arrays laid out
    // over them. A rim vertex is `classof`'d like any other. `bases[r]` starts ring `r` in the
    // concatenated edge/vertex space.
    let boundary: Vec<&Loop> = std::iter::once(&face.outer)
        .chain(crossed.iter().map(|&k| &face.inner[k]))
        .collect();
    let ring_edges: Vec<usize> = boundary.iter().map(|l| l.half_edges.len()).collect();
    let bases: Vec<usize> = ring_edges
        .iter()
        .scan(0, |s, &len| {
            let b = *s;
            *s += len;
            Some(b)
        })
        .collect();
    let verts: Vec<Handle<Vertex>> = boundary
        .iter()
        .flat_map(|l| l.half_edges.iter().map(|&he| he_start(model, he)))
        .collect();
    let kept: Vec<bool> = verts.iter().map(|v| classof[v] == keep).collect();
    let n = verts.len();

    // The concatenated-space position of the edge a boundary node rides. A crossed rim edge
    // now resolves to its slot — this is where `SEAM_ACROSS_HOLE_RIM` used to reject. An edge
    // on no boundary loop (an uncrossed rim, or foreign) is the two machines disagreeing.
    let edge_pos: HashMap<Handle<Edge>, usize> = boundary
        .iter()
        .flat_map(|l| l.half_edges.iter())
        .enumerate()
        .map(|(i, he)| (he.edge, i))
        .collect();
    let edge_ix = |e: Handle<Edge>| {
        edge_pos
            .get(&e)
            .copied()
            .ok_or_else(|| reject(tag::SEAM_COUNT_MISMATCH))
    };
    // The boundary crossings, grouped by the edge each rides. An edge may carry more than
    // one: that is an edge the other solid's boundary enters and leaves through, and it is
    // what an edge-indexed model cannot express.
    let mut by_edge: Vec<Vec<[usize; 3]>> = vec![Vec::new(); n];
    for nd in paths.iter().flat_map(|p| p.nodes()) {
        if let Some(e) = nd.on_edge {
            by_edge[edge_ix(e)?].push(nd.triple);
        }
    }
    let n_bnd: usize = by_edge.iter().map(|v| v.len()).sum();

    // The `seam` list (built by `pierced_faces`) and the arrangement's paths must agree that
    // this face is met at all; and every arc has exactly two boundary ends, so a boundary
    // node anywhere else — an arc's interior touching `∂f` — is caught by the count.
    let meets_face = !paths.is_empty();
    if touches_seam != meets_face || n_bnd != 2 * opens.len() {
        return Err(reject(tag::SEAM_COUNT_MISMATCH));
    }

    if opens.is_empty() && n_closed == 0 {
        // Defensive, and believed unreachable: `paths.is_empty()` forces
        // `touches_seam == false`, and then `transitions` must be empty too — the pair
        // that already returned above.
        return Ok(whole()); // no seam on this face after all
    }

    // `∂f`'s rings as triples, outer then each crossed rim, co-indexed with `by_edge`/`verts`.
    // Built whenever a run is emitted (`region_rings` reads `bnd[i]`) or an edge carries two —
    // its construction rejects a straight angle. `all_holes` also feeds `place_loops` below.
    let all_holes = arrange::hole_rings(model, fh, plane_idx, inc_f)?;
    let need_bnd = by_edge.iter().any(|v| v.len() > 1) || n_closed > 0 || !face.inner.is_empty();
    let bnd: Vec<[usize; 3]> = if need_bnd {
        let mut b = arrange::face_vertex_triples(model, fh, plane_idx, inc_f)?;
        for &k in &crossed {
            b.extend(all_holes[k].iter().copied());
        }
        b
    } else {
        Vec::new()
    };
    // Each ring is cut into runs on its own — a rim wraps to itself, not to the outer loop —
    // and the runs are concatenated. `boundary_runs` returns slice-relative vertex indices, so
    // `+base` lifts them into the concatenated `verts`/`bnd`. `ring_runs[r]` is ring `r`'s run
    // count, the size of its block in the crossing/run index space.
    let (crossings, runs, run_kept, ring_runs) = if opens.is_empty() {
        (
            Vec::new(),
            vec![(0..n).collect::<Vec<_>>()],
            Vec::new(),
            Vec::new(),
        )
    } else {
        let mut crossings = Vec::new();
        let mut runs: Vec<Vec<usize>> = Vec::new();
        let mut ring_runs = Vec::with_capacity(boundary.len());
        for (r, &base) in bases.iter().enumerate() {
            let span = base..base + ring_edges[r];
            let ring_by_edge = &by_edge[span.clone()];
            // `boundary_runs` reads `bnd` only to order two crossings sharing an edge; a ring
            // with at most one crossing per edge needs none, and `bnd` may be empty for it (a
            // single-crossing hole-free face never built it — a straight angle it carries has
            // no business rejecting).
            let ring_bnd: &[[usize; 3]] = if ring_by_edge.iter().any(|v| v.len() > 1) {
                &bnd[span]
            } else {
                &[]
            };
            let mut br = arrange::boundary_runs(planes, plane_idx, ring_bnd, ring_by_edge)?;
            for run in &mut br.runs {
                for v in run.iter_mut() {
                    *v += base;
                }
            }
            ring_runs.push(br.runs.len());
            crossings.extend(br.crossings);
            runs.extend(br.runs);
        }
        // Alternation gives the run classes their shape, `classof` anchors them, and every
        // vertex is checked against the propagation. Each ring is classed on its own, so a rim
        // and the outer loop may carry opposite classes.
        let rk = arrange::run_classes(&runs, &kept, &ring_runs)?;
        (crossings, runs, rk, ring_runs)
    };
    let n_x = crossings.len();
    debug_assert!(
        opens.is_empty() || ring_runs.iter().sum::<usize>() == n_x,
        "the per-ring run blocks tile the crossing/run index space"
    );
    // The per-ring cyclic successor over the crossing/run index space, and its inverse. `next`
    // steps forward one position along `∂f` without leaving its ring; `prev[c]` names the run
    // flowing into crossing `c`. A single ring makes these `(i ± 1) % n_x`.
    let mut next = vec![0usize; n_x];
    let mut prev = vec![0usize; n_x];
    {
        let mut base = 0;
        for &rl in &ring_runs {
            for j in 0..rl {
                let (cur, nxt) = (base + j, base + (j + 1) % rl);
                next[cur] = nxt;
                prev[nxt] = cur;
            }
            base += rl;
        }
    }
    let cross_ix: HashMap<[usize; 3], usize> =
        crossings.iter().enumerate().map(|(j, &t)| (t, j)).collect();
    // Crossing `c` is kept→dropped exactly when the run flowing into it is kept.
    let is_kd = |c: usize| run_kept[prev[c]];

    // Each arc, oriented `s_kd → bends → s_dk`: the direction its kept run is spliced in.
    // Its ends are found by *triple*, not by edge — both may ride the same one.
    //
    // `stitch_cycles` reads a transition as "position `t` kept, `t+1` dropped", which in this
    // index space names crossing `next[t]`. So an arc's slot is `prev` of its crossing — the
    // run flowing into it, the multi-ring form of "crossing minus one".
    let mut kd = Vec::with_capacity(opens.len());
    let mut dk = Vec::with_capacity(opens.len());
    let mut arc_nodes: Vec<Vec<[usize; 3]>> = Vec::with_capacity(opens.len());
    for arc in &opens {
        let ends = arc.nodes();
        let at = |nd: &arrange::SeamEnd| {
            cross_ix
                .get(&nd.triple)
                .copied()
                .ok_or_else(|| reject(tag::SEAM_COUNT_MISMATCH))
        };
        let (h, t) = (at(&ends[0])?, at(&ends[ends.len() - 1])?);
        let forward = is_kd(h);
        if forward == is_kd(t) {
            return Err(reject(tag::SEAM_COUNT_MISMATCH)); // an arc's ends must oppose
        }
        let ordered: Vec<[usize; 3]> = if forward {
            ends.iter().map(|nd| nd.triple).collect()
        } else {
            ends.iter().rev().map(|nd| nd.triple).collect()
        };
        for tr in &ordered {
            if !seam_ix.contains_key(tr) {
                return Err(reject(tag::MISSING_SEAM));
            }
        }
        let (kd_c, dk_c) = if forward { (h, t) } else { (t, h) };
        kd.push(prev[kd_c]);
        dk.push(prev[dk_c]);
        arc_nodes.push(ordered);
    }

    // The kept regions of `f`, as walks over indices — a boundary run of original vertices,
    // or an arc spliced whole. The arcs' successor map decomposes into cycles, one region
    // each; with no arc the only region is `f` itself, kept iff `∂f` is.
    enum Step {
        Run(usize),
        Arc(usize),
    }
    let regions: Vec<Vec<Step>> = if opens.is_empty() {
        Vec::from_iter(kept[0].then(|| vec![Step::Run(0)]))
    } else {
        arrange::stitch_cycles(&run_kept, &kd, &dk, &next)?
            .into_iter()
            .map(|cycle| {
                let mut steps = Vec::new();
                for w in 0..cycle.len() {
                    let (pre, this) = (cycle[(w + cycle.len() - 1) % cycle.len()], cycle[w]);
                    let mut i = next[dk[pre]];
                    loop {
                        steps.push(Step::Run(i));
                        if i == kd[this] {
                            break;
                        }
                        i = next[i];
                    }
                    steps.push(Step::Arc(this));
                }
                steps
            })
            .collect()
    };

    let region_face = |steps: &[Step], inner: Vec<Vec<Node>>| LocalFace {
        plane_idx,
        loop_nodes: steps
            .iter()
            .flat_map(|s| match s {
                Step::Run(j) => runs[*j]
                    .iter()
                    .map(|&i| Node::Orig(verts[i]))
                    .collect::<Vec<_>>(),
                Step::Arc(a) => arc_nodes[*a].iter().copied().map(Node::Seam).collect(),
            })
            .collect(),
        inner,
        flip,
    };

    if n_closed == 0 && face.inner.is_empty() {
        return Ok(regions.iter().map(|r| region_face(r, vec![])).collect());
    }

    // A loop is a hole of the kept region containing it, and an island when no kept region
    // does — its interior is then what survives. Only containment can say which, and with an
    // arc present nothing else can: a hole in the kept region beside an island in the dropped
    // one wind oppositely without nesting (cell 3f-3's counterexample, realized by
    // `l_and_staple`). The winding is kept as the cross-check, not the decision.
    //
    // `∂f`'s vertices are three-plane points too, so a region's ring is one uniform list.
    let region_rings: Vec<Vec<[usize; 3]>> = regions
        .iter()
        .map(|steps| {
            steps
                .iter()
                .flat_map(|s| match s {
                    Step::Run(j) => runs[*j].iter().map(|&i| bnd[i]).collect(),
                    Step::Arc(a) => arc_nodes[*a].clone(),
                })
                .collect()
        })
        .collect();

    let mut rings = Vec::with_capacity(n_closed);
    let mut windings = Vec::with_capacity(n_closed);
    for path in &paths {
        let arrange::SeamPath::Closed(nodes) = path else {
            continue;
        };
        let ring = arrange::orient_seam_loop(planes, plane_idx, nodes, keep == Side::Outside)?;
        for t in &ring {
            if !seam_ix.contains_key(t) {
                return Err(reject(tag::MISSING_SEAM));
            }
        }
        windings.push(arrange::loop_winding(planes, plane_idx, &ring)?);
        rings.push(ring);
    }

    // `f`'s own holes join the loop list, and `classify_nesting` places them too, with the
    // seam loops, in one containment forest: which material face — a region or a seam island —
    // owns each, and which loops are islands in their own right. A hole tangled inside a seam
    // loop, which the old pairwise guard over-rejected as `nested_loops`, now finds its owner.
    //
    // Only the *uncrossed* holes: a crossed one was absorbed into the outer boundary above, so
    // placing it again would count it twice. `all` is the seam rings then the holes; indices
    // `>= rings.len()` are holes, remapped through `uncrossed` back to `f.inner`.
    //
    // A hole's winding joins `windings` and is checked by `classify_nesting` like any loop's —
    // it must run clockwise, whatever region or island it falls in.
    let uncrossed: Vec<usize> = (0..face.inner.len())
        .filter(|k| !crossed.contains(k))
        .collect();
    let holes: Vec<Vec<[usize; 3]>> = uncrossed.iter().map(|&k| all_holes[k].clone()).collect();
    let n_rings = rings.len();
    let mut all: Vec<Vec<[usize; 3]>> = rings;
    all.extend(holes);
    for h in &all[n_rings..] {
        windings.push(arrange::loop_winding(planes, plane_idx, h)?);
    }
    let nest = nest_loops(planes, plane_idx, &region_rings, &all)?;
    let place = classify_nesting(&nest, n_rings, &windings)?;

    let node_ring = |r: &[[usize; 3]]| r.iter().copied().map(Node::Seam).collect::<Vec<Node>>();
    // Loop `j`'s nodes: a seam ring as `Node::Seam`, an `f` hole as `Node::Orig`. Emitting a
    // hole's own vertices as `Node::Seam` would mint a fresh `Discovered` vertex on each — the
    // arrangement is a source of combinatorics, never of geometry, and here that is a node kind.
    let inner_of = |j: usize| -> Vec<Node> {
        if j < n_rings {
            node_ring(&all[j])
        } else {
            orig_ring(&face.inner[uncrossed[j - n_rings]])
        }
    };
    let holes_of = |owner: Owner| -> Vec<Vec<Node>> {
        (0..all.len())
            .filter(|&j| place[j] == Placement::Hole(owner))
            .map(&inner_of)
            .collect()
    };

    // Regions first, each with its holes; then islands, each with its holes. `assemble_fuse_cut`
    // fixes vertex handles by first appearance across `faces`, and replay rests on that order.
    let mut out: Vec<LocalFace> = regions
        .iter()
        .enumerate()
        .map(|(k, steps)| region_face(steps, holes_of(Owner::Region(k))))
        .collect();
    for j in 0..n_rings {
        if place[j] == Placement::Island {
            out.push(LocalFace {
                plane_idx,
                loop_nodes: node_ring(&all[j]),
                inner: holes_of(Owner::Island(j)),
                flip,
            });
        }
    }
    Ok(out)
}

/// Push the reconstructed result and supersede the inputs (mirrors `assemble`).
fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    planes: &[PlaneInfo],
    seam: &[SeamVertex],
    faces: &[LocalFace],
) -> Result<Handle<Solid>, BoolError> {
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
    // Materialize all vertex handles first. Outer then inner, rings in order: `vh`'s
    // first-appearance order fixes the vertex handles, and replay depends on it.
    for lf in faces {
        for &node in lf.loop_nodes.iter().chain(lf.inner.iter().flatten()) {
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
        let mut ring = |model: &mut Model, nodes: &[Node]| -> Loop {
            let handles: Vec<Handle<Vertex>> = nodes.iter().map(|nd| vh[nd]).collect();
            let k = handles.len();
            let mut half_edges: Vec<HalfEdge> = (0..k)
                .map(|t| {
                    let (va, vb) = (handles[t], handles[(t + 1) % k]);
                    let e = edge_for(model, va, vb);
                    let forward = model.edges.get(e).bounds.expect("bounded")[0] == va;
                    HalfEdge { edge: e, forward }
                })
                .collect();
            if lf.flip {
                // Cut's inside-A B-pieces: reverse every loop and toggle the
                // orientation below, so the outward normal points into the removed
                // region and each loop still keeps material on its left.
                half_edges.reverse();
                for he in &mut half_edges {
                    he.forward = !he.forward;
                }
            }
            Loop { half_edges }
        };
        let outer = ring(model, &lf.loop_nodes);
        let inner: Vec<Loop> = lf.inner.iter().map(|h| ring(model, h)).collect();
        let orientation = if lf.flip {
            match planes[lf.plane_idx].orient {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
            }
        } else {
            planes[lf.plane_idx].orient
        };
        face_handles.push(model.faces.push(Face {
            surface: planes[lf.plane_idx].surf,
            outer,
            inner,
            orientation,
        }));
    }
    // Partition the faces into connected components (by shared node). One component is the
    // whole result; several mean either an enclosed void (a cavity — one component holds
    // material outside it, signed flux negative) or a severed operand (two solids — two
    // components each holding material inside, flux positive). The sign tells them apart:
    // exactly one positive component is the outer shell, the rest are cavities; anything else
    // is two solids one handle cannot answer (`Cut(rod, L)`), which stays `DISCONNECTED_RESULT`.
    let (labels, n) = face_components(faces);
    let mut by_comp: Vec<Vec<Handle<Face>>> = vec![Vec::new(); n];
    for (i, &fh) in face_handles.iter().enumerate() {
        by_comp[labels[i]].push(fh);
    }
    let positives: Vec<usize> = (0..n)
        .filter(|&c| shell_signed_flux(model, &by_comp[c]) > 0.0)
        .collect();
    let [outer_c] = positives[..] else {
        return Err(reject(tag::DISCONNECTED_RESULT));
    };
    let shells: Vec<Handle<Shell>> = by_comp
        .into_iter()
        .map(|faces| model.shells.push(Shell { faces }))
        .collect();
    // A cavity shell's faces already point into the void (the material is outside it, so the
    // material-on-correct-side reconstruction winds them inward) — measured, so no reversal.
    let cavities = (0..n)
        .filter(|&c| c != outer_c)
        .map(|c| shells[c])
        .collect();
    let solid = model.push_solid(Solid {
        outer: shells[outer_c],
        cavities,
    });
    model.live_solids.retain(|&s| s != a && s != b);
    Ok(solid)
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
    if !is_convex(model, &planes_a, &solid_vertex_handles(model, a))
        || !is_convex(model, &planes_b, &solid_vertex_handles(model, b))
    {
        return None;
    }
    // Opposite-normal coplanar face pairs (same-normal coplanar side faces are
    // the expected coplanar-adjacent result and are ignored).
    let mut opposite: Vec<(usize, usize)> = Vec::new();
    for (i, pa) in planes_a.iter().enumerate() {
        for (j, pb) in planes_b.iter().enumerate() {
            if planes_coplanar(&pa.plane, &pb.plane) && pa.n_out.dot(pb.n_out) < 0.0 {
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

/// A bijective **exact-coordinate** match B-ring → A-ring, or `None` if the
/// boundaries are not identical (different length, an unmatched or ambiguous
/// vertex, or two B vertices sharing an A vertex). A coincident interface's
/// corresponding vertices are literally the same point, so exact `Point3`
/// equality decides the match — no tolerance (design §3 (5d)-3). Input analysis
/// only.
fn interface_correspondence(
    ring_a: &[(Handle<Vertex>, Point3)],
    ring_b: &[(Handle<Vertex>, Point3)],
) -> Option<HashMap<Handle<Vertex>, Handle<Vertex>>> {
    if ring_a.len() != ring_b.len() {
        return None;
    }
    let mut remap = HashMap::new();
    for &(bh, bp) in ring_b {
        let found: Vec<Handle<Vertex>> = ring_a
            .iter()
            .filter(|&&(_, ap)| ap == bp)
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

/// A coplanar face-on-face contact where one face is strictly contained in the other — a
/// boss touching a face inside its boundary (cell coplanar-contact-boss). Unlike a coincident
/// interface the footprints differ, so the containing face keeps its boundary and gains the
/// contained footprint as a hole.
struct ContainedContact {
    big_solid: Handle<Solid>,
    big_face: Handle<Face>,
    small_solid: Handle<Solid>,
    small_face: Handle<Face>,
}

/// Whether `small`'s outer boundary is strictly inside `big`'s region (inside `big`'s outer
/// loop, outside its holes, touching no boundary), both on the shared plane with normal `n`.
/// Exact 2D containment on the plane (drop the dominant axis) — reuses the imprint-containment
/// predicate.
fn face_contains_face(model: &Model, big: Handle<Face>, small: Handle<Face>, n: Vector3) -> bool {
    let drop = planar_drop_axes(n);
    let proj_loop = |l: &Loop| -> Vec<[f64; 2]> {
        l.half_edges
            .iter()
            .map(|&he| proj2(model.vertices.get(he_start(model, he)).point, drop))
            .collect()
    };
    let bf = model.faces.get(big);
    let small_2d = proj_loop(&model.faces.get(small).outer);
    let outer_2d = proj_loop(&bf.outer);
    let holes_2d: Vec<Vec<[f64; 2]>> = bf.inner.iter().map(&proj_loop).collect();
    profile_strictly_in_region(&small_2d, &outer_2d, &holes_2d)
}

/// `Some` iff A and B meet at exactly one opposite-normal coplanar face pair whose footprints
/// are strictly nested (one contained in the other) — the boss-on-a-face case that
/// `detect_coincident_interface` (which needs identical boundaries) does not cover.
fn detect_contained_contact(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<ContainedContact> {
    let planes_a = collect_planes(model, a).ok()?;
    let planes_b = collect_planes(model, b).ok()?;
    if !is_convex(model, &planes_a, &solid_vertex_handles(model, a))
        || !is_convex(model, &planes_b, &solid_vertex_handles(model, b))
    {
        return None;
    }
    let mut opposite: Vec<(usize, usize)> = Vec::new();
    for (i, pa) in planes_a.iter().enumerate() {
        for (j, pb) in planes_b.iter().enumerate() {
            if planes_coplanar(&pa.plane, &pb.plane) && pa.n_out.dot(pb.n_out) < 0.0 {
                opposite.push((i, j));
            }
        }
    }
    if opposite.len() != 1 {
        return None;
    }
    let (i, j) = opposite[0];
    let fa = model.shells.get(model.solids.get(a).outer).faces[i];
    let fb = model.shells.get(model.solids.get(b).outer).faces[j];
    let n = planes_a[i].plane.normal();
    if face_contains_face(model, fa, fb, n) {
        Some(ContainedContact {
            big_solid: a,
            big_face: fa,
            small_solid: b,
            small_face: fb,
        })
    } else if face_contains_face(model, fb, fa, n) {
        Some(ContainedContact {
            big_solid: b,
            big_face: fb,
            small_solid: a,
            small_face: fa,
        })
    } else {
        None
    }
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
        let ring = |l: &Loop| -> Vec<Node> {
            l.half_edges
                .iter()
                .map(|&he| {
                    let vh = he_start(model, he);
                    let mapped = remap.and_then(|r| r.get(&vh)).copied().unwrap_or(vh);
                    Node::Orig(mapped)
                })
                .collect()
        };
        out.push(LocalFace {
            plane_idx: plane_offset + pos,
            loop_nodes: ring(&face.outer),
            inner: face.inner.iter().map(ring).collect(),
            flip: false,
        });
    }
    out
}

/// The vertex handle of an original-face node (coincident-merge loops carry no seam nodes).
fn node_vh(n: Node) -> Handle<Vertex> {
    match n {
        Node::Orig(v) => v,
        Node::Seam(_) => unreachable!("coincident-merge faces are all Node::Orig"),
    }
}

/// Splice A's and B's coplanar side faces, which share the interface edge at `lf_a`'s
/// position `ia` (`a[ia] → a[ia+1]`), into one loop: walk A from the edge's far end all the
/// way round to its near end, then insert B's complementary path. The shared interface edge
/// is dropped; its two endpoints remain as (collinear, for a right prism) boundary points.
fn splice_side_faces(
    lf_a: &LocalFace,
    lf_b: &LocalFace,
    ia: usize,
    dissolve: &HashSet<Handle<Vertex>>,
) -> LocalFace {
    let (a, b) = (&lf_a.loop_nodes, &lf_b.loop_nodes);
    let (m, n) = (a.len(), b.len());
    let x = node_vh(a[ia]);
    let y = node_vh(a[(ia + 1) % m]);
    // B traverses the shared edge the opposite way: `b[ib] = Y`, `b[ib+1] = X`.
    let ib = (0..n)
        .find(|&j| node_vh(b[j]) == y && node_vh(b[(j + 1) % n]) == x)
        .expect("B's side face shares the interface edge, opposite orientation");
    let mut nodes = Vec::with_capacity(m + n - 2);
    // A from Y (ia+1) forward all the way to X (ia): the whole A loop, less the dropped edge.
    for t in 0..m {
        nodes.push(a[(ia + 1 + t) % m]);
    }
    // B's interior, strictly between X and Y (skip the shared edge's two endpoints).
    for t in 2..n {
        nodes.push(b[(ib + t) % n]);
    }
    // Drop the straight-angle interface vertices: each is flanked by its A and B far
    // corners, so removing it fuses the split vertical edge into one (cell fuse-coplanar-merge).
    nodes.retain(|nd| !dissolve.contains(&node_vh(*nd)));
    LocalFace {
        plane_idx: lf_a.plane_idx,
        loop_nodes: nodes,
        inner: Vec::new(),
        flip: false,
    }
}

/// For a coincident Fuse, merge each interface edge's two incident side faces (one from A,
/// one from B) into a single face **when they are coplanar** — the flat-edge defeature that
/// lets the fused solid chain (cell fuse-coplanar-merge). A non-coplanar pair is a genuine
/// dihedral (a slanted operand) and stays separate; a holed side face stays separate; cap
/// faces (no interface edge) pass through. `faces_a`/`faces_b` already skip the interface
/// faces, with B's interface vertices remapped to A's handles.
fn merge_coincident_fuse_faces(
    model: &Model,
    iface: &Interface,
    faces_a: Vec<LocalFace>,
    faces_b: Vec<LocalFace>,
    planes: &[PlaneInfo],
) -> Vec<LocalFace> {
    let ring: Vec<Handle<Vertex>> = model
        .faces
        .get(iface.fa)
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let iface_edges: HashSet<(usize, usize)> = (0..ring.len())
        .map(|i| {
            unordered(
                ring[i].index() as usize,
                ring[(i + 1) % ring.len()].index() as usize,
            )
        })
        .collect();
    // A side face has exactly one interface edge; return its key and the position of the
    // edge's start node in the loop. A cap face has none.
    let side_edge = |lf: &LocalFace| -> Option<((usize, usize), usize)> {
        let k = lf.loop_nodes.len();
        (0..k).find_map(|i| {
            let u = node_vh(lf.loop_nodes[i]).index() as usize;
            let w = node_vh(lf.loop_nodes[(i + 1) % k]).index() as usize;
            let key = unordered(u, w);
            iface_edges.contains(&key).then_some((key, i))
        })
    };
    let edge_to_face = |faces: &[LocalFace]| -> HashMap<(usize, usize), usize> {
        faces
            .iter()
            .enumerate()
            .filter_map(|(idx, lf)| side_edge(lf).map(|(key, _)| (key, idx)))
            .collect()
    };
    let a_side = edge_to_face(&faces_a);
    let b_side = edge_to_face(&faces_b);
    // An interface edge merges when both incident side faces are coplanar and hole-free.
    let merged: HashSet<(usize, usize)> = iface_edges
        .iter()
        .filter(|key| match (a_side.get(key), b_side.get(key)) {
            (Some(&ia), Some(&ib)) => {
                let (la, lb) = (&faces_a[ia], &faces_b[ib]);
                la.inner.is_empty()
                    && lb.inner.is_empty()
                    && planes_coplanar(&planes[la.plane_idx].plane, &planes[lb.plane_idx].plane)
            }
            _ => false,
        })
        .copied()
        .collect();
    // An interface vertex whose *both* incident interface edges merge is a straight angle:
    // its A and B vertical edges are the same line `P_i ∩ P_{i-1}`. Dissolve it (drop from
    // the merged loops) so the split vertical edge fuses into one — a clean box.
    let nr = ring.len();
    let dissolve: HashSet<Handle<Vertex>> = (0..nr)
        .filter(|&i| {
            let prev = unordered(
                ring[(i + nr - 1) % nr].index() as usize,
                ring[i].index() as usize,
            );
            let next = unordered(
                ring[i].index() as usize,
                ring[(i + 1) % nr].index() as usize,
            );
            merged.contains(&prev) && merged.contains(&next)
        })
        .map(|i| ring[i])
        .collect();

    let mut b_taken = vec![false; faces_b.len()];
    let mut out = Vec::new();
    for lf_a in faces_a {
        let Some((key, ia)) = side_edge(&lf_a) else {
            out.push(lf_a); // cap
            continue;
        };
        if merged.contains(&key) {
            let ib = b_side[&key];
            b_taken[ib] = true;
            out.push(splice_side_faces(&lf_a, &faces_b[ib], ia, &dissolve));
        } else {
            out.push(lf_a); // real dihedral, or holed: keep separate
        }
    }
    for (ib, lf_b) in faces_b.into_iter().enumerate() {
        if !b_taken[ib] {
            out.push(lf_b);
        }
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
    // `detect_coincident_interface` counts only *cross-solid, opposite-normal* coplanar
    // pairs, so an imprint on some other face — whose coplanar region face is
    // same-normal and within one solid — does not disqualify the stack. Such an operand
    // reaches here with a holed face, and `solid_local_faces` carries the hole through.
    // This is the one path that admits an imprint: the seam paths never see one, because
    // `has_coplanar_pair` stops it at their door.
    match kind {
        // The intersection is the flat shared face (zero volume).
        BoolKind::Common => Err(BoolError::EmptyResult),
        // A and B are on opposite sides of the interface, so B removes nothing:
        // A − B is a fresh copy of A (its interface face survives).
        BoolKind::Cut => {
            let planes_a = collect_planes(model, a)?;
            let faces = solid_local_faces(model, a, 0, None, None);
            assemble_fuse_cut(model, a, b, &planes_a, &[], &faces)
        }
        // Drop both interface faces; keep every other face, sewing B's interface ring to A's
        // shared vertices. Each interface edge's two coplanar side faces (one from A, one
        // from B) are spliced into one, dropping the flat edge so the result chains
        // (cell fuse-coplanar-merge).
        BoolKind::Fuse => {
            let planes_a = collect_planes(model, a)?;
            let na = planes_a.len();
            let planes_b = collect_planes(model, b)?;
            let mut planes = planes_a;
            planes.extend(planes_b);
            let faces_a = solid_local_faces(model, a, 0, Some(iface.fa), None);
            let faces_b = solid_local_faces(model, b, na, Some(iface.fb), Some(&iface.remap));
            let faces = merge_coincident_fuse_faces(model, iface, faces_a, faces_b, &planes);
            assemble_fuse_cut(model, a, b, &planes, &[], &faces)
        }
    }
}

/// Fuse a boss: `small` sits on `big`'s face inside its boundary. `big`'s contact face keeps
/// its boundary and gains `small`'s footprint as a hole; `small`'s contact face is dropped;
/// the shared footprint edges stitch the hole to `small`'s walls (`assemble_fuse_cut`'s
/// `edge_for` dedups them). `small`'s footprint region becomes interior — solid on both sides,
/// so no face there — which is exactly what the dropped pair leaves (cell coplanar-contact-boss).
fn contained_boss_fuse(
    model: &mut Model,
    cc: &ContainedContact,
) -> Result<Handle<Solid>, BoolError> {
    let planes_big = collect_planes(model, cc.big_solid)?;
    let na = planes_big.len();
    let planes_small = collect_planes(model, cc.small_solid)?;
    let mut planes = planes_big;
    planes.extend(planes_small);

    // Big's faces, skipping its contact face — re-added below as an annulus.
    let mut faces = solid_local_faces(model, cc.big_solid, 0, Some(cc.big_face), None);
    let big_shell = model.solids.get(cc.big_solid).outer;
    let big_pos = model
        .shells
        .get(big_shell)
        .faces
        .iter()
        .position(|&f| f == cc.big_face)
        .expect("contact face is on big's shell");
    let ring_nodes = |l: &Loop| -> Vec<Node> {
        l.half_edges
            .iter()
            .map(|&he| Node::Orig(he_start(model, he)))
            .collect()
    };
    let big_f = model.faces.get(cc.big_face);
    let mut inner: Vec<Vec<Node>> = big_f.inner.iter().map(&ring_nodes).collect();
    // `small`'s outer loop is CCW about its own (opposite) normal, i.e. CW about big's normal —
    // exactly the hole winding, used as-is so its edges stay opposed to small's walls.
    inner.push(ring_nodes(&model.faces.get(cc.small_face).outer));
    faces.push(LocalFace {
        plane_idx: big_pos,
        loop_nodes: ring_nodes(&big_f.outer),
        inner,
        flip: false,
    });
    // Small's faces except the dropped contact face.
    faces.extend(solid_local_faces(
        model,
        cc.small_solid,
        na,
        Some(cc.small_face),
        None,
    ));

    assemble_fuse_cut(model, cc.big_solid, cc.small_solid, &planes, &[], &faces)
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use proptest::prelude::*;

    /// The single-ring successor `(i + 1) % n` — what every caller but a multi-ring `∂f`
    /// hands [`arrange::stitch_cycles`].
    fn ident_next(n: usize) -> Vec<usize> {
        (0..n).map(|i| (i + 1) % n).collect()
    }

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

    /// The same L-prism, its profile started one vertex earlier so the reflex corner
    /// `(1,1)` lands at index 1 of the cap's loop. Geometrically identical.
    fn rotated_l_prism() -> (Model, Handle<Solid>) {
        let l = Profile2d {
            points: vec![
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
                p2(0.0, 0.0),
                p2(2.0, 0.0),
            ],
        };
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// `PlaneInfo::n_out` is documented as the single source of "outward". Two
    /// independent sources say which way that is: the ring, which winds CCW about the
    /// outward normal, and the b-rep's own `Surface` plus `Orientation`. They must agree
    /// on every face of every solid.
    ///
    /// `outer_tri` used to read the turn at the first non-collinear corner, which is the
    /// ring's winding **only when that corner is convex**. Nothing enforced that. The
    /// four fixtures below were safe by accident — none starts its cap loop one vertex
    /// before a reflex corner. `rotated_l_prism` does, and it is the same solid.
    #[test]
    fn outward_normals_agree_with_their_orientation() {
        let mut cube = Model::new();
        let c = cube.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let (ml, sl) = l_prism();
        let (mu, su) = u_prism();
        let (mr, sr) = rotated_l_prism();
        for (name, m, s) in [
            ("cube", &cube, c),
            ("l_prism", &ml, sl),
            ("u_prism", &mu, su),
            ("rotated_l_prism", &mr, sr),
        ] {
            for pi in &collect_planes(m, s).unwrap() {
                let dot = pi.plane.normal().dot(pi.n_out);
                assert!(
                    dot.abs() > 0.5,
                    "{name}: n_out is not parallel to its plane"
                );
                assert_eq!(
                    dot > 0.0,
                    pi.orient == Orientation::Forward,
                    "{name}: n_out disagrees with the face's orientation"
                );
            }
        }
    }

    #[test]
    fn a_rotated_profile_is_the_same_solid_to_the_boolean() {
        // The silent half of the same fault. `orient_seam_loop` multiplies `orient_sign`
        // of two planes, and nothing pairs that with a `tri` whose sign would cancel it
        // (as `order_along` does — there both factors flip together). So an inward `n_out`
        // on the L's cap reverses the hole, and only `validate` and the signed mesh volume
        // would ever say so, in release where the `debug_assert` is gone.
        //
        // Rotating a profile cannot change a solid. Here the cap is the face that carries
        // the dimple's hole, so this is the shortest path from the fault to a wrong b-rep.
        let (mut m, l) = rotated_l_prism();
        let stub = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 0.5]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        let r = boolean(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.08)).abs() < 1e-9, "volume {vol}");
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

    /// Cell 3g opened non-convex `Common`: the same seam the corner bite arranges for
    /// `Cut`/`Fuse` (`cut_non_convex_overlap_corner_bite` says `Cut = 3.0 − 0.224`, so the
    /// overlap is `0.224`), now kept as the intersection. Only the keep/flip table differs.
    #[test]
    fn common_non_convex_overlap_is_their_intersection() {
        let (mut m, l, bx) = l_and_corner_box();
        let r = boolean(&mut m, BoolKind::Common, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.224).abs() < 1e-9, "volume {vol}");
    }

    /// The first `Common` whose seam closes into a loop — where cell 3g's `Inside` keep
    /// first reaches `orient_seam_loop`. A bar drilled through a cube: `∩ = [1,2]²×[0,3]`,
    /// the middle segment of the bar. On the cube's `z=0` and `z=3` caps the kept square
    /// `[1,2]²` is bounded entirely by seam (an island face, no `∂f`), so the loop's
    /// orientation runs through `orient_seam_loop` with `material_outside = false` — the
    /// sign a hole (Cut) never exercised. `1·1·3`.
    #[test]
    fn a_common_can_leave_a_closed_seam_loop() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let r = boolean(&mut m, BoolKind::Common, cube, bar).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 3.0).abs() < 1e-9, "volume {vol}");
    }

    /// The seam path answers a *convex* overlap `Common` — `two_boxes ∩ = [0.5,1]³` with
    /// the same 8 verts / 12 edges / 6 faces the deleted convex `common` gave. This was
    /// the measurement cell 3g's n1 took (through `overlap_fuse_cut` directly) before n2
    /// deleted `common`; `common_of_two_cubes_is_their_overlap` now exercises the same
    /// result through the dispatcher.
    #[test]
    fn the_seam_path_answers_common_overlap() {
        let (mut m, a, b) = two_boxes();
        let r = overlap_fuse_cut(&mut m, BoolKind::Common, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let reach = m.reachable();
        assert_eq!(reach.vertices.len(), 8);
        assert_eq!(reach.edges.len(), 12);
        assert_eq!(reach.faces.len(), 6);
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.125).abs() < 1e-12, "volume {vol}");
    }

    /// The L with a box straddling its reflex corner (1,1): a *single* chord with one
    /// reflex bend. The box vertex `(1.6,1.6,·)` sits in the L's notch — inside the
    /// convex hull, outside the L — exactly where a convex half-space test would
    /// misclassify it `Inside`. Overlap = `xy(1.0 − notch 0.36) · z(0.8)` = `0.512`.
    fn l_and_reflex_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.6, 0.6, 0.2]),
            Point3::from_array([1.6, 1.6, 1.4]),
        );
        (m, l, bx)
    }

    #[test]
    fn cut_across_reflex_corner_bite() {
        // Still transitions==2, so it reconstructs into a correct L-shaped face. The
        // strongest classification test: only exact `point_in_solid` gets the volume
        // right, a convex half-space test would not.
        let (mut m, l, bx) = l_and_reflex_box();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        // Overlap = xy(1.0 − notch 0.36 = 0.64) · z(0.8) = 0.512.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
    }

    /// Drilling. A rod skewers the L's bottom bar, both its ends outside, and what comes
    /// out is a solid of genus 1 — the first this kernel has ever made.
    ///
    /// The guard that stood here said a through-hole was "not yet representable". It had
    /// been representable since cell 3f-1 put an inner loop on a face; what could not
    /// follow were the seam list, which built no vertex for an edge whose endpoints agree,
    /// and the boundary model, which indexed `∂f` by edge. Both are cell 3e-3's.
    ///
    /// `3 − 0.2·0.3·1`. Each cap takes the rod's cross-section as a hole, and the rod's
    /// four walls contribute the tunnel — each of them a face whose kept region is bounded
    /// by two arcs and two runs holding no vertex at all.
    #[test]
    fn drill_through_the_l() {
        let (mut m, l, rod) = l_and_rod();
        let r = boolean(&mut m, BoolKind::Cut, l, rod).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "genus 1 is clean: {vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 2.94).abs() < 1e-9, "{}", props.volume);
        assert!((props.area - 14.88).abs() < 1e-9, "{}", props.area);
        assert_eq!(holed_faces(&m, r).len(), 2, "both caps");
    }

    /// The `Fuse`: the rod stands proud on both faces of the bar, so the caps take the
    /// same two holes and each rod wall splits into the stub above and the stub below.
    /// `3 + 0.12 − 0.06`.
    #[test]
    fn fuse_the_l_and_the_rod() {
        let (mut m, l, rod) = l_and_rod();
        let r = boolean(&mut m, BoolKind::Fuse, l, rod).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 3.06).abs() < 1e-9, "{}", props.volume);
        assert!((props.area - 15.0).abs() < 1e-9, "{}", props.area);
        assert_eq!(holed_faces(&m, r).len(), 2);
    }

    /// The same two solids the other way round: the bar severs the rod, and `Cut` must
    /// answer with two solids. It cannot — a `Solid` has one outer shell.
    ///
    /// Measured with the guard removed: `Ok`, volume `0.06`, and `validate` reporting
    /// `NegativeGenus { v: 16, e: 24, f: 12, genus: -1 }` — two disjoint boxes in one
    /// shell. `boolean` never runs `validate`, so nothing else would have said a word.
    /// `pierced_multi` had been hiding this: severing A takes an edge of A through B.
    #[test]
    fn cut_rod_by_l_disconnects() {
        let (mut m, l, rod) = l_and_rod();
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, rod, l),
            tag::DISCONNECTED_RESULT,
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
    /// misses. No edge threads the other solid, so the arcs are the whole story.
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
    ) -> (Vec<PlaneInfo>, HashMap<Handle<Face>, usize>) {
        let mut planes = collect_planes(m, a).unwrap();
        planes.extend(collect_planes(m, b).unwrap());
        let surf_ix = planes
            .iter()
            .enumerate()
            .map(|(i, pi)| (pi.face, i))
            .collect();
        (planes, surf_ix)
    }

    /// The face of `solid` whose outward normal is `n` (there is exactly one, for
    /// the axis-aligned fixtures here).
    fn face_facing(
        m: &Model,
        solid: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Face>, usize>,
        n: [f64; 3],
    ) -> Handle<Face> {
        let want = Vector3::from_array(n);
        let shell = m.solids.get(solid).outer;
        let mut hit = None;
        for &fh in &m.shells.get(shell).faces {
            let pi = &planes[surf_ix[&fh]];
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
        surf_ix: &HashMap<Handle<Face>, usize>,
    ) -> Vec<arrange::SeamSegment> {
        let inc_x = arrange::edge_planes(m, x, surf_ix).unwrap();
        let inc_y = arrange::edge_planes(m, y, surf_ix).unwrap();
        arrange::seam_segments_on(m, f, y, planes, surf_ix, &inc_x, &inc_y).unwrap()
    }

    fn paths(
        m: &Model,
        f: Handle<Face>,
        x: Handle<Solid>,
        y: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Face>, usize>,
    ) -> Vec<arrange::SeamPath> {
        let inc_x = arrange::edge_planes(m, x, surf_ix).unwrap();
        let inc_y = arrange::edge_planes(m, y, surf_ix).unwrap();
        arrange::seam_paths_on(m, f, y, planes, surf_ix, &inc_x, &inc_y).unwrap()
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

        /// The assembly's structure, checked without repeating the assembly. A node is
        /// a boundary node (on `∂f`, two A-planes in its triple) or an interior node
        /// (a B-edge piercing `f`, two B-planes); an arc has exactly the former at its
        /// two ends and the latter within; a loop has none.
        ///
        /// Like the proptest above, an honest `Err` is skipped, so this says nothing
        /// about rejected inputs. And it says nothing about closed loops: two boxes
        /// sharing a corner cannot make one. Measured over a 144-instance sweep of
        /// this generator: 859 faces accepted, 5 rejected, 381 arcs, **0 loops**. The
        /// loop and the doubly-crossed edge are covered by hand-built fixtures, not
        /// by chance.
        #[test]
        fn prop_seam_paths_have_manifold_structure(
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
            let len_a = collect_planes(&m, a).unwrap().len();
            let (planes, surf_ix) = combined(&m, a, b);
            let (Ok(inc_x), Ok(inc_y)) = (
                arrange::edge_planes(&m, a, &surf_ix),
                arrange::edge_planes(&m, b, &surf_ix),
            ) else {
                return Ok(());
            };
            for &f in &m.shells.get(m.solids.get(a).outer).faces {
                let Ok(out) = arrange::seam_paths_on(&m, f, b, &planes, &surf_ix, &inc_x, &inc_y)
                else {
                    continue;
                };
                for path in &out {
                    let n = path.nodes();
                    for e in n {
                        // A plane index < len_a belongs to A. Boundary ⇔ two of them.
                        let from_a = e.triple.iter().filter(|&&i| i < len_a).count();
                        prop_assert_eq!(e.on_edge.is_some(), from_a == 2, "{:?}", e);
                        prop_assert!(from_a == 1 || from_a == 2, "{e:?}");
                    }
                    match path {
                        arrange::SeamPath::Open(_) => {
                            prop_assert!(n.len() >= 2);
                            prop_assert!(n[0].on_edge.is_some() && n[n.len() - 1].on_edge.is_some());
                            prop_assert!(n[1..n.len() - 1].iter().all(|e| e.on_edge.is_none()));
                        }
                        arrange::SeamPath::Closed(_) => {
                            prop_assert!(n.len() >= 3);
                            prop_assert!(n.iter().all(|e| e.on_edge.is_none()));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn seam_segments_on_two_overlapping_boxes() {
        // A = [0,2]×[0,2.2]×[0,2.4] and B = [1,3.5]×[1,3.2]×[1,3.4] share a corner.
        // Extents were unequal because with two cubes the seam lands on a face's *centre*,
        // where both fan diagonals cross; `segment_crosses_face` grazed from every apex
        // and said `contact_degenerate`. Cell (5a) deleted the fan, and `two_boxes` now
        // runs this shape at the centre. The coordinates stay: one variable at a time.
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
        // A face split into two chords. The slab meets the U-prism's base cap
        // (z=0) along {z=0, y=1.5}; the U has material there only for x∈[0,1] and
        // x∈[2,3], so the seam is *two* segments — which cell 3e-2 taught the seam path
        // to stitch (the old `multichord` guard retired there).
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

    #[test]
    fn seam_paths_split_a_face_into_two_open_arcs() {
        // The same base cap as `seam_segments_split_a_face_into_two_chords`, now
        // assembled. Two disjoint arcs, each a lone segment whose both ends sit on
        // `∂cap` — the two-chord shape the retired `multichord` guard once rejected.
        let (m, u, slab) = u_and_slab();
        let (planes, surf_ix) = combined(&m, u, slab);
        let cap = face_facing(&m, u, &planes, &surf_ix, [0.0, 0.0, -1.0]);

        let out = paths(&m, cap, u, slab, &planes, &surf_ix);
        assert_eq!(out.len(), 2, "two arcs on one face: {out:?}");
        for p in &out {
            assert!(matches!(p, arrange::SeamPath::Open(_)), "{p:?}");
            let n = p.nodes();
            assert_eq!(n.len(), 2);
            // Both ends are boundary nodes, so both carry the cap edge they lie on.
            assert!(n.iter().all(|e| e.on_edge.is_some()), "{n:?}");
        }
        let mut spans: Vec<(f64, f64)> = out
            .iter()
            .map(|p| {
                let (a, b) = (p.nodes()[0].point[0], p.nodes()[1].point[0]);
                (a.min(b), a.max(b))
            })
            .collect();
        spans.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        assert!(
            (spans[0].0).abs() < 1e-9 && (spans[0].1 - 1.0).abs() < 1e-9,
            "{spans:?}"
        );
        assert!(
            (spans[1].0 - 2.0).abs() < 1e-9 && (spans[1].1 - 3.0).abs() < 1e-9,
            "{spans:?}"
        );
    }

    /// Canonicalize a node sequence: an arc up to reversal, a loop up to rotation and
    /// reversal. `seam_paths_on` does fix a start and a direction — that is a contract
    /// and `seam_paths_are_deterministic` pins it — but these goldens are about the
    /// *arc order*, and should not also freeze which end the walk entered from.
    fn canon(pts: &[[f64; 3]], closed: bool) -> Vec<[f64; 3]> {
        let rots = if closed { pts.len() } else { 1 };
        (0..rots)
            .flat_map(|i| {
                let rot: Vec<[f64; 3]> = pts[i..].iter().chain(&pts[..i]).copied().collect();
                let mut rev = rot.clone();
                rev.reverse();
                [rot, rev]
            })
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or_default()
    }

    fn assert_points(p: &arrange::SeamPath, want: &[[f64; 3]]) {
        let closed = matches!(p, arrange::SeamPath::Closed(_));
        let got: Vec<[f64; 3]> = p.nodes().iter().map(|n| n.point.as_array()).collect();
        let (got, want) = (canon(&got, closed), canon(want, closed));
        assert_eq!(got.len(), want.len(), "{got:?} vs {want:?}");
        assert!(
            got.iter()
                .zip(&want)
                .all(|(g, w)| near(Point3::from_array(*g), *w)),
            "{got:?} vs {want:?}"
        );
    }

    #[test]
    fn seam_paths_follow_the_true_arc_order_through_a_bend() {
        // The staircase of `cut_staircase_seam_arc`, assembled. On the box's bottom face
        // the arc runs (2,0.5) → (2,1) → (1,1) → (1,1.5): two bends turning opposite
        // ways. The old projection sort could not tell a reflex turn from an arc folded
        // back on itself, so `strict` rejected both. Adjacency reads this order off the
        // 1-manifold instead, and cell 3e-1 retired the guard.
        let (m, l, bx) = l_and_popup_box();
        let (planes, surf_ix) = combined(&m, l, bx);
        let floor = face_facing(&m, bx, &planes, &surf_ix, [0.0, 0.0, -1.0]);

        let out = paths(&m, floor, bx, l, &planes, &surf_ix);
        assert_eq!(out.len(), 1, "one arc: {out:?}");
        assert!(matches!(out[0], arrange::SeamPath::Open(_)), "{:?}", out[0]);
        assert_points(
            &out[0],
            &[
                [2.0, 0.5, 0.2],
                [2.0, 1.0, 0.2],
                [1.0, 1.0, 0.2],
                [1.0, 1.5, 0.2],
            ],
        );
        // Ends on `∂floor`, bends where the L's two vertical edges pierce it.
        let n = out[0].nodes();
        assert!(n[0].on_edge.is_some() && n[3].on_edge.is_some());
        assert!(n[1].on_edge.is_none() && n[2].on_edge.is_none());
    }

    #[test]
    fn seam_paths_close_a_loop_inside_a_face() {
        // The dimple of `cut_blind_dimple_is_unsupported`. The stub's footprint never
        // reaches the top face's boundary, so the seam is a closed ring in the face
        // interior — an inner loop. Zero boundary nodes, four interior ones, one per
        // vertical stub edge.
        let (m, l, stub) = l_and_dimple();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);

        let out = paths(&m, top, l, stub, &planes, &surf_ix);
        assert_eq!(out.len(), 1, "one loop: {out:?}");
        assert!(
            matches!(out[0], arrange::SeamPath::Closed(_)),
            "{:?}",
            out[0]
        );
        assert!(out[0].nodes().iter().all(|n| n.on_edge.is_none()));
        assert_points(
            &out[0],
            &[
                [0.3, 0.3, 1.0],
                [0.3, 0.7, 1.0],
                [0.7, 0.7, 1.0],
                [0.7, 0.3, 1.0],
            ],
        );
    }

    /// A big cube whose `y=0, z=0` edge is crossed **twice** by the seam: the notch
    /// spans `x∈[3,7]` and hangs below both `y=0` and `z=0`, so that one edge enters
    /// and leaves it. Extents are asymmetric so no crossing lands on a face centre or a
    /// fan diagonal — a debt to the fan, which cell (5a) deleted. Kept: one variable at
    /// a time.
    ///
    /// `edge_seam` maps an edge to *one* seam triple. This input is what that map
    /// cannot represent — and `boolean` rejects it today (see
    /// `an_edge_crossed_twice_is_rejected`).
    fn cube_and_notch() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let y = m.add_cuboid(
            Point3::from_array([3.0, -1.0, -1.0]),
            Point3::from_array([7.0, 1.4, 1.2]),
        );
        (m, a, y)
    }

    #[test]
    fn seam_path_crosses_one_face_edge_twice() {
        // Both boundary nodes lie on the *same* edge of `f`. Any splice keyed by
        // `Handle<Edge>` keeps one and drops the other — a manifold face that
        // `validate` accepts and that is wrong. Cell 3d must key by the seam triple
        // and order the two crossings along the edge with `order_along`.
        let (m, a, y) = cube_and_notch();
        let (planes, surf_ix) = combined(&m, a, y);
        let floor = face_facing(&m, a, &planes, &surf_ix, [0.0, 0.0, -1.0]);

        let out = paths(&m, floor, a, y, &planes, &surf_ix);
        assert_eq!(out.len(), 1, "one arc: {out:?}");
        assert_points(
            &out[0],
            &[
                [3.0, 0.0, 0.0],
                [3.0, 1.4, 0.0],
                [7.0, 1.4, 0.0],
                [7.0, 0.0, 0.0],
            ],
        );
        let n = out[0].nodes();
        assert_eq!(
            n[0].on_edge, n[3].on_edge,
            "both ends ride the cube's y=0,z=0 edge"
        );
        assert!(n[0].on_edge.is_some());
    }

    /// A thin rod skewering the L's bottom bar in `z`, both ends outside. Each of its four
    /// vertical edges pierces the L's two caps, so the caps take a closed seam loop and the
    /// rod's walls take two chords apiece — and each wall's vertical edges are crossed
    /// **twice**, leaving runs with no vertex at all.
    fn l_and_rod() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let rod = m.add_cuboid(
            Point3::from_array([0.3, 0.3, -0.5]),
            Point3::from_array([0.5, 0.6, 1.5]),
        );
        (m, l, rod)
    }

    /// Group a face's boundary seam nodes by the index of the edge each rides.
    fn by_edge(m: &Model, f: Handle<Face>, paths: &[arrange::SeamPath]) -> Vec<Vec<[usize; 3]>> {
        let hes = &m.faces.get(f).outer.half_edges;
        let mut out = vec![Vec::new(); hes.len()];
        for nd in paths.iter().flat_map(|p| p.nodes()) {
            if let Some(e) = nd.on_edge {
                let i = hes
                    .iter()
                    .position(|he| he.edge == e)
                    .expect("an outer edge");
                out[i].push(nd.triple);
            }
        }
        out
    }

    /// Cell 3e-3's core object. The cube's floor has all four corners outside the notch and
    /// yet the seam enters and leaves across one edge, so an edge-indexed model sees one
    /// transition slot where there are two crossings. The run model sees two crossings and
    /// two runs, one of which holds **no vertex** — the piece of that edge inside the notch.
    #[test]
    fn two_crossings_on_one_edge_make_a_run_with_no_vertex() {
        let (m, a, y) = cube_and_notch();
        let (planes, surf_ix) = combined(&m, a, y);
        let floor = face_facing(&m, a, &planes, &surf_ix, [0.0, 0.0, -1.0]);
        let p = surf_ix[&floor];
        let inc = arrange::edge_planes(&m, a, &surf_ix).unwrap();
        let ps = paths(&m, floor, a, y, &planes, &surf_ix);
        let be = by_edge(&m, floor, &ps);
        assert_eq!(be.iter().map(|v| v.len()).max(), Some(2));

        let bnd = arrange::face_vertex_triples(&m, floor, p, &inc).unwrap();
        let br = arrange::boundary_runs(&planes, p, &bnd, &be).unwrap();
        assert_eq!(br.crossings.len(), 2);
        let mut lens: Vec<usize> = br.runs.iter().map(|r| r.len()).collect();
        assert_eq!(br.runs.len(), 2);
        lens.sort_unstable();
        assert_eq!(
            lens,
            vec![0, 4],
            "one run holds nothing, the other the whole floor"
        );

        // The empty run runs from the first crossing the walk meets to the second, so the
        // pair must be ordered the way the doubly-crossed edge is directed. That edge is
        // `y = 0, z = 0`; the notch cuts it at `x = 3` and `x = 7`.
        let empty = br.runs.iter().position(|r| r.is_empty()).unwrap();
        let x_of = |t: [usize; 3]| {
            three_planes(
                &planes[t[0]].plane,
                &planes[t[1]].plane,
                &planes[t[2]].plane,
            )
            .unwrap()
            .as_array()[0]
        };
        let hes = &m.faces.get(floor).outer.half_edges;
        let on_y0 = hes
            .iter()
            .position(|he| {
                let b = m.edges.get(he.edge).bounds.unwrap();
                b.iter()
                    .all(|&v| m.vertices.get(v).point.as_array()[1] == 0.0)
            })
            .unwrap();
        let b = m.edges.get(hes[on_y0].edge).bounds.unwrap();
        let start = m.vertices.get(he_start(&m, hes[on_y0])).point.as_array()[0];
        let end = b
            .iter()
            .map(|&v| m.vertices.get(v).point.as_array()[0])
            .find(|&x| x != start)
            .unwrap();
        let (first, second) = if end > start { (3.0, 7.0) } else { (7.0, 3.0) };
        assert!((x_of(br.crossings[empty]) - first).abs() < 1e-9);
        assert!((x_of(br.crossings[(empty + 1) % 2]) - second).abs() < 1e-9);
    }

    /// The rod's wall: four crossings on two edges, two apiece, and **two** runs with no
    /// vertex. The kept region of that face (for a `Cut`) is bounded by arcs alone.
    #[test]
    fn a_rod_wall_has_two_runs_with_no_vertex() {
        let (m, l, rod) = l_and_rod();
        let (planes, surf_ix) = combined(&m, l, rod);
        let wall = face_facing(&m, rod, &planes, &surf_ix, [-1.0, 0.0, 0.0]);
        let p = surf_ix[&wall];
        let inc = arrange::edge_planes(&m, rod, &surf_ix).unwrap();
        let ps = paths(&m, wall, rod, l, &planes, &surf_ix);
        assert_eq!(ps.len(), 2, "two chords, one per cap: {ps:?}");

        let be = by_edge(&m, wall, &ps);
        let bnd = arrange::face_vertex_triples(&m, wall, p, &inc).unwrap();
        let br = arrange::boundary_runs(&planes, p, &bnd, &be).unwrap();
        assert_eq!(br.crossings.len(), 4);
        let mut lens: Vec<usize> = br.runs.iter().map(|r| r.len()).collect();
        lens.sort_unstable();
        assert_eq!(lens, vec![0, 0, 2, 2]);

        // `Cut(l, rod)` keeps the rod's inside-L piece, so the vertex-bearing runs (which
        // hold the wall's corners, all outside L) are dropped and the empty ones are kept.
        let kept_vert: Vec<bool> = m
            .faces
            .get(wall)
            .outer
            .half_edges
            .iter()
            .map(|&he| {
                point_in_solid(&m, m.vertices.get(he_start(&m, he)).point, l).unwrap()
                    == Side::Inside
            })
            .collect();
        assert_eq!(kept_vert, vec![false; 4]);
        let kept = arrange::run_classes(&br.runs, &kept_vert, &[br.runs.len()]).unwrap();
        assert_eq!(kept.iter().filter(|&&k| k).count(), 2);
        for (j, r) in br.runs.iter().enumerate() {
            assert_eq!(kept[j], r.is_empty(), "run {j}");
        }
    }

    /// `run_classes` on paper. Alternation shapes it, `classof` anchors it, and a run with
    /// no vertex has only the alternation to go by.
    #[test]
    fn run_classes_alternate_and_are_anchored_by_a_vertex() {
        // `cube_and_notch`'s floor: run 0 empty (inside the notch), run 1 the four corners.
        let runs = vec![vec![], vec![0, 1, 2, 3]];
        let kept = arrange::run_classes(&runs, &[true; 4], &[2]).unwrap();
        assert_eq!(kept, vec![false, true]);

        // A vertex disagreeing with the propagation is the two machineries in conflict.
        assert_rejects(
            || arrange::run_classes(&[vec![0], vec![1]], &[true, true], &[2]),
            tag::SEAM_COUNT_MISMATCH,
        );
        // Crossings alternate enter/exit around a closed curve, so an odd count is a lie.
        assert_rejects(
            || arrange::run_classes(&[vec![0]], &[true], &[1]),
            tag::SEAM_COUNT_MISMATCH,
        );
        // Nothing anchors a boundary made only of vertex-free runs.
        assert_rejects(
            || arrange::run_classes(&[vec![], vec![]], &[], &[2]),
            tag::SEAM_COUNT_MISMATCH,
        );

        // Two rings, seeded apart (cell 3f-6). The outer runs `[_, {0,1}]` are kept, the rim
        // runs `[_, {2,3}]` dropped — the firing lid, whose outer boundary lies in the box and
        // whose pocket rim lies out of it. A single global alternation would carry the outer
        // seed onto the rim and conflict with the rim's own vertices; per-ring seeding is what
        // lets the two loops disagree. Result `[F, T, T, F]`.
        let runs = vec![vec![], vec![0, 1], vec![], vec![2, 3]];
        let kept_vert = [true, true, false, false];
        assert_eq!(
            arrange::run_classes(&runs, &kept_vert, &[2, 2]).unwrap(),
            vec![false, true, true, false],
            "each ring seeded by its own vertex"
        );
        // A ring with an odd run count: its crossings do not alternate around a closed loop.
        assert_rejects(
            || arrange::run_classes(&[vec![0], vec![1], vec![2], vec![3]], &[true; 4], &[1, 3]),
            tag::SEAM_COUNT_MISMATCH,
        );
        // A ring of only vertex-free runs, beside a well-anchored one: no anchor, still a lie.
        assert_rejects(
            || arrange::run_classes(&[vec![], vec![], vec![0], vec![1]], &[true; 2], &[2, 2]),
            tag::SEAM_COUNT_MISMATCH,
        );
        // A vertex inside one ring disagreeing with that ring's propagation — the two machines
        // in conflict, caught by the global check even when the other ring is clean.
        assert_rejects(
            || arrange::run_classes(&[vec![], vec![0, 1], vec![2], vec![3]], &[true; 4], &[2, 2]),
            tag::SEAM_COUNT_MISMATCH,
        );
    }

    /// The integer algebra cell 3e-3 rests on, before a line of it is wired: with crossings
    /// as the index space, `kept[t] && !kept[t+1]` says "crossing `t+1` is kept→dropped",
    /// so an arc's `kd` slot is its kd **crossing minus one**. `stitch_cycles` is then the
    /// same function it always was — it never knew what its indices meant.
    #[test]
    fn crossing_indices_feed_stitch_cycles_unchanged() {
        // The cube's floor: crossings `[c0 (x=3), c1 (x=7)]`, runs `[empty, corners]`.
        // The single arc runs c1 → c0 through the notch; c1 is dropped→kept, c0 the other.
        let kept = [false, true];
        let (kd_cross, dk_cross) = (0usize, 1usize);
        let n = 2;
        let kd = vec![(kd_cross + n - 1) % n];
        let dk = vec![(dk_cross + n - 1) % n];
        assert_eq!((kd.clone(), dk.clone()), (vec![1], vec![0]));
        assert_eq!(
            arrange::stitch_cycles(&kept, &kd, &dk, &ident_next(kept.len())).unwrap(),
            vec![vec![0]]
        );

        // The rod's wall: crossings `c0,c1` on one vertical edge and `c2,c3` on the other,
        // the chords `c1–c2` (at `z = 1`) and `c3–c0` (at `z = 0`). An arc's kd end is the
        // crossing whose *incoming* run is kept — `is_kd(c_j) = kept[j-1]`.
        let n = 4;
        let kd_dk = |kept: [bool; 4], arcs: [[usize; 2]; 2]| {
            let is_kd = |c: usize| kept[(c + n - 1) % n];
            let (mut kd, mut dk) = (vec![], vec![]);
            for [x, y] in arcs {
                let (k, d) = if is_kd(x) { (x, y) } else { (y, x) };
                assert!(is_kd(k) && !is_kd(d), "an arc's ends oppose");
                kd.push((k + n - 1) % n);
                dk.push((d + n - 1) % n);
            }
            (kd, dk)
        };
        let arcs = [[1, 2], [3, 0]];

        // `Cut(l, rod)` keeps the vertex-free runs: one region, the rectangle inside L.
        let (kd, dk) = kd_dk([true, false, true, false], arcs);
        assert_eq!((kd.clone(), dk.clone()), (vec![0, 2], vec![1, 3]));
        let kept = [true, false, true, false];
        assert_eq!(
            arrange::stitch_cycles(&kept, &kd, &dk, &ident_next(kept.len())).unwrap(),
            vec![vec![0, 1]],
            "one cycle: the wall's middle"
        );

        // `Fuse` keeps the other two: the stub above L and the stub below, two faces.
        let (kd, dk) = kd_dk([false, true, false, true], arcs);
        assert_eq!((kd.clone(), dk.clone()), (vec![1, 3], vec![0, 2]));
        let kept = [false, true, false, true];
        assert_eq!(
            arrange::stitch_cycles(&kept, &kd, &dk, &ident_next(kept.len())).unwrap(),
            vec![vec![0], vec![1]],
            "two cycles: one face becomes two"
        );
    }

    /// Two rings, threaded by arcs, before the multi-ring `∂f` is wired (cell 3f-6). The
    /// firing fixture's lid: outer loop crossed twice, the pocket rim crossed twice, and two
    /// arcs each running outer↔rim. In the crossing/run index space the outer ring takes
    /// `{0, 1}` and the rim `{2, 3}`, so `next` cycles *within* each block — `[1, 0, 3, 2]`,
    /// not the single ring's `(i + 1) % 4`.
    ///
    /// The kept runs are `0` (outer) and `2` (rim). Arc `0` joins the outer kd to the rim dk,
    /// arc `1` the rim kd to the outer dk; `succ` walking each arc's `dk` ring lands on the
    /// other ring's kd, so the two arcs close into **one** cycle spanning both rings. Feed the
    /// same runs the single-ring `next` and each arc closes on *itself* — two cycles, the one
    /// connected region wrongly split. That split is the whole of what the successor map fixes.
    #[test]
    fn two_rings_thread_into_one_cycle() {
        let run_kept = [true, false, true, false]; // outer {0,1}, rim {2,3}, each alternating
        let next = [1, 0, 3, 2]; // per-ring cyclic successor
        // `prev` is `next`'s inverse — the caller derives an arc's slot as `prev[crossing]`,
        // the multi-ring form of cell 3e-3's `crossing − 1`.
        let mut prev = [0usize; 4];
        for (i, &j) in next.iter().enumerate() {
            prev[j] = i;
        }
        assert_eq!(prev, [1, 0, 3, 2]);
        // kd crossings `1` (outer) and `3` (rim); dk crossings `0` (outer) and `2` (rim). Slots
        // via `prev`: arc 0 = (kd slot `prev[1]=0`, dk slot `prev[2]=3`), arc 1 = (kd `prev[3]=2`,
        // dk `prev[0]=1`).
        let kd = [prev[1], prev[3]]; // [0, 2]
        let dk = [prev[2], prev[0]]; // [3, 1]
        assert_eq!((kd, dk), ([0, 2], [3, 1]));
        assert_eq!(
            arrange::stitch_cycles(&run_kept, &kd, &dk, &next).unwrap(),
            vec![vec![0, 1]],
            "one cycle threads the outer ring and the rim"
        );
        // The single-ring successor severs the crossing: `succ` cannot leave the ring `dk`
        // sits on, so each arc closes on itself and the one region splits into two.
        assert_eq!(
            arrange::stitch_cycles(&run_kept, &kd, &dk, &ident_next(4)).unwrap(),
            vec![vec![0], vec![1]],
            "single-ring next wrongly splits the region"
        );
    }

    /// `reconstruct_face_paths`'s kept/dropped alternation edges on `f`'s outer loop, as
    /// indices into `f.outer.half_edges`.
    ///
    /// Independent of `keep`: flipping it flips every `kept[i]`, and this only
    /// compares neighbours. So no `BoolKind` need be chosen here.
    fn transitions_on(m: &Model, f: Handle<Face>, other: Handle<Solid>) -> Vec<usize> {
        let hes = &m.faces.get(f).outer.half_edges;
        let inside: Vec<bool> = hes
            .iter()
            .map(|&he| {
                let v = m.vertices.get(he_start(m, he)).point;
                point_in_solid(m, v, other).unwrap() == Side::Inside
            })
            .collect();
        (0..inside.len())
            .filter(|&i| inside[i] != inside[(i + 1) % inside.len()])
            .collect()
    }

    #[test]
    fn boundary_nodes_match_the_transition_count() {
        // An independent oracle: `reconstruct_face_paths` finds the seam's ends by counting
        // kept/dropped alternations around `∂f`; the arrangement finds them as
        // degree-1 nodes. They must agree — through entirely different machinery
        // (ray-cast classification vs. exact segment/face crossings).
        //
        // The identity needs one hypothesis: no edge of `f` is crossed twice. An edge
        // crossed twice has equal endpoint classes, contributing no transition, while
        // the arrangement sees both nodes. `cube_and_notch` is exactly that case, and
        // it is checked below as the documented counterexample rather than hidden.
        let cases: Vec<(Model, Handle<Solid>, Handle<Solid>)> = vec![
            l_and_corner_box(),
            u_and_slab(),
            l_and_popup_box(),
            l_and_dimple(),
        ];
        for (m, x, y) in cases {
            let (planes, surf_ix) = combined(&m, x, y);
            for &f in &m.shells.get(m.solids.get(x).outer).faces {
                let boundary = paths(&m, f, x, y, &planes, &surf_ix)
                    .iter()
                    .flat_map(|p| p.nodes())
                    .filter(|n| n.on_edge.is_some())
                    .count();
                assert_eq!(boundary, transitions_on(&m, f, y).len(), "face {f:?}");
            }
        }
    }

    /// Every non-convex overlap input that reaches `reconstruct_face_paths`: the three that
    /// succeed and the three that it rejects from *inside*. The rejecting three matter
    /// most — it calls `seam_paths_on`, which can reject for reasons the old convex code
    /// never could (`vertex_on_face_plane`, `point_on_ring`, `fourplane`). Should that fire
    /// on a face visited *before* the intended one, the reject tag silently changes.
    /// Measuring only the accepted inputs would not see it.
    struct OverlapCase {
        name: &'static str,
        /// `(most open arcs on any one face, closed loops over all faces)`.
        ///
        /// The arrangement's own reading of each fixture, and it moves as the ladder
        /// climbs. Arcs stopped mattering at cell 3e-2, and loops at 3f-3 — however many,
        /// so long as they share a winding. What is left for `pokehole` is a loop *beside*
        /// an arc. `(1, 0)` says the arrangement never had anything to say about this
        /// fixture at all.
        expect: (usize, usize),
        m: Model,
        x: Handle<Solid>,
        y: Handle<Solid>,
    }

    fn overlap_fixtures() -> Vec<OverlapCase> {
        let case = |name, expect, (m, x, y)| OverlapCase {
            name,
            expect,
            m,
            x,
            y,
        };
        vec![
            // `l_and_corner_box` serves both the Cut and the Fuse test: the arrangement
            // never sees `BoolKind`.
            case("l_and_corner_box", (1, 0), l_and_corner_box()),
            case("l_and_reflex_box", (1, 0), l_and_reflex_box()),
            // Two loops, measured, not guessed: the U pierces the slab's `y=1.5` face
            // twice — its two prongs cut rectangles wholly inside that face. So this
            // fixture needs *both* multi-chord and inner loops before it can pass;
            // `multichord` merely happens to fire first, on the U's base cap.
            case("u_and_slab (two flat loops)", (2, 2), u_and_slab()),
            // `(1, 0)`: nothing about the arrangement ever blocked this one. What
            // rejected it was the shape of a single arc, until cell 3e-1.
            case("l_and_popup_box (folded arc)", (1, 0), l_and_popup_box()),
            case("l_and_dimple (inner loop)", (1, 1), l_and_dimple()),
            // `(2, 0)`: two chords on one face and not a single closed loop — the only
            // fixture that asks for multi-chord alone. Measured.
            case("l_and_notch_bar (two chords)", (2, 0), l_and_notch_bar()),
            // The first non-convex inner loop: six nodes, one of them reflex.
            case("l_and_ell_stub (non-convex loop)", (1, 1), l_and_ell_stub()),
            // An arc and a loop on one face — the only shape `pokehole` keeps after 3f-3.
            case("l_and_staple (arc beside a loop)", (2, 1), l_and_staple()),
        ]
    }

    #[test]
    fn arrangement_agrees_with_todays_seam_bookkeeping() {
        // Cell 3d swaps `reconstruct_face`'s projection sort for the arrangement's arc
        // order. This proves the swap is output-preserving *before* any production code
        // moves. Every claim below is one the wired code will rely on.
        let (mut n_faces, mut n_arcs, mut n_loops) = (0usize, 0usize, 0usize);
        for OverlapCase {
            name,
            expect,
            m,
            x,
            y,
        } in overlap_fixtures()
        {
            let (mut max_open, mut tot_closed) = (0usize, 0usize);
            for (f_solid, o_solid) in [(x, y), (y, x)] {
                let (planes, surf_ix) = combined(&m, x, y);
                for &f in &m.shells.get(m.solids.get(f_solid).outer).faces {
                    let inc_f = arrange::edge_planes(&m, f_solid, &surf_ix).unwrap();
                    let inc_o = arrange::edge_planes(&m, o_solid, &surf_ix).unwrap();

                    // (1) The precondition of the whole cell: no face rejects.
                    let out =
                        arrange::seam_paths_on(&m, f, o_solid, &planes, &surf_ix, &inc_f, &inc_o)
                            .unwrap_or_else(|e| panic!("{name}: face {f:?} rejected: {e:?}"));

                    let hes = &m.faces.get(f).outer.half_edges;
                    let opens: Vec<&arrange::SeamPath> = out
                        .iter()
                        .filter(|p| matches!(p, arrange::SeamPath::Open(_)))
                        .collect();
                    let closed = out.len() - opens.len();

                    // (2) Boundary nodes ride exactly the transition edges. This is the
                    // set equality the wired code checks in production; a count would
                    // pass even with both ends on one edge.
                    let bnd: std::collections::BTreeSet<usize> = out
                        .iter()
                        .flat_map(|p| p.nodes())
                        .filter_map(|n| n.on_edge)
                        .map(|e| {
                            hes.iter()
                                .position(|he| he.edge == e)
                                .expect("a boundary node rides an edge of its own face")
                        })
                        .collect();
                    let trans: std::collections::BTreeSet<usize> =
                        transitions_on(&m, f, o_solid).into_iter().collect();
                    assert_eq!(bnd, trans, "{name}: face {f:?}");

                    max_open = max_open.max(opens.len());
                    tot_closed += closed;

                    // (4) On a single chord the *bends* — the interior nodes, the only
                    // thing today sorts — are strictly increasing along the chord
                    // direction, so today's stable projection sort reproduces the arc
                    // order exactly. Ties would not: the sort would then fall back on
                    // `seam`'s insertion order, which the arrangement knows nothing of.
                    //
                    // The endpoints are *not* part of this. A reflex bend projects
                    // outside the chord interval — on `l_and_reflex_box` the bends land
                    // at −0.24 and 0.96 against endpoints 0.0 and 0.72. "Bends lie
                    // between the endpoints" is false, and does not need to be true.
                    n_faces += 1;
                    n_arcs += opens.len();
                    n_loops += closed;
                    if opens.len() == 1 {
                        let pts: Vec<Point3> = opens[0].nodes().iter().map(|n| n.point).collect();
                        let dir = pts[pts.len() - 1] - pts[0];
                        let proj: Vec<f64> = pts[1..pts.len() - 1]
                            .iter()
                            .map(|&p| (p - pts[0]).dot(dir))
                            .collect();
                        assert!(
                            proj.windows(2).all(|w| w[0] < w[1]),
                            "{name}: face {f:?} bends not strictly monotone: {proj:?}"
                        );
                    }
                }
            }
            // (3) What the arrangement sees, per fixture. `l_and_popup_box` claiming
            // `(1, 0)` is what says only the arc's shape blocks it — necessary, not
            // sufficient: the wired code still checks the seam list agrees.
            assert_eq!((max_open, tot_closed), expect, "{name}");
        }
        // Pin the exercise. A fixture edit that quietly stops reaching the arrangement
        // would otherwise leave every assertion above vacuously true.
        assert_eq!((n_faces, n_arcs, n_loops), (122, 64, 5));
    }

    #[test]
    fn a_doubly_crossed_edge_breaks_the_transition_oracle() {
        // The domain boundary of the oracle above, kept as a test rather than a
        // comment. The cube's floor has zero transitions — all four of its corners
        // sit outside the notch — yet the seam enters and leaves across one edge.
        let (m, a, y) = cube_and_notch();
        let (planes, surf_ix) = combined(&m, a, y);
        let floor = face_facing(&m, a, &planes, &surf_ix, [0.0, 0.0, -1.0]);

        assert_eq!(transitions_on(&m, floor, y).len(), 0);
        let boundary = paths(&m, floor, a, y, &planes, &surf_ix)
            .iter()
            .flat_map(|p| p.nodes())
            .filter(|n| n.on_edge.is_some())
            .count();
        assert_eq!(boundary, 2);
    }

    #[test]
    fn seam_paths_are_deterministic() {
        // The assembly queries its adjacency map but never iterates it: `HashMap`
        // order varies per instance, and a path's start and direction must not. Cell
        // 3d builds faces from these paths, and the operation log replays them.
        let (m, l, bx) = l_and_popup_box();
        let (planes, surf_ix) = combined(&m, l, bx);
        let floor = face_facing(&m, bx, &planes, &surf_ix, [0.0, 0.0, -1.0]);
        let once = paths(&m, floor, bx, l, &planes, &surf_ix);
        for _ in 0..8 {
            assert_eq!(paths(&m, floor, bx, l, &planes, &surf_ix), once);
        }
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
        let (mut m, l, bx) = l_and_reflex_box();
        let from_arrange = seam_endpoint_triples(&m, l, bx);
        assert_eq!(from_arrange.len(), 8);
        boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // And the folded arc (cell 3e-1). `discovered_triples` reads `reachable()`, so
        // a seam vertex the reconstruction built but no result face used would drop out
        // of the set and break this — which is what a lost face looks like. Volume
        // alone could coincide; this cannot.
        let (mut m, l, bx) = l_and_popup_box();
        let from_arrange = seam_endpoint_triples(&m, l, bx);
        boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // And the inner loop (cell 3f-1): its four nodes are the hole's rim, used once
        // by the lid and once by a pocket wall. Lose the loop and they leave `reachable`.
        let (mut m, l, stub) = l_and_dimple();
        let from_arrange = seam_endpoint_triples(&m, l, stub);
        assert_eq!(from_arrange.len(), 4);
        boolean(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // Two chords (cell 3e-2). Six seam nodes: two boundary ends and a bend for each of
        // the cap's arcs, and the same six seen again from the bar's floor. Lose one of the
        // bar's two floor faces and its four nodes leave `reachable`.
        let (mut m, l, bar) = l_and_notch_bar();
        let from_arrange = seam_endpoint_triples(&m, l, bar);
        boolean(&mut m, BoolKind::Cut, l, bar).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // The non-convex hole (cell 3h): six rim nodes, none of them on `∂f`.
        let (mut m, l, stub) = l_and_ell_stub();
        let from_arrange = seam_endpoint_triples(&m, l, stub);
        assert_eq!(from_arrange.len(), 6);
        boolean(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // An arc beside a loop (cell 3f-4). Lose the loop, or place it on the wrong region,
        // and its four rim nodes leave `reachable`.
        let (mut m, l, st) = l_and_staple();
        let from_arrange = seam_endpoint_triples(&m, l, st);
        boolean(&mut m, BoolKind::Cut, l, st).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // Two flat loops (cell 3f-3): sixteen nodes, and the two prong cross-sections must
        // both survive as faces or their eight rim nodes leave `reachable`.
        let (mut m, u, slab) = u_and_slab();
        let from_arrange = seam_endpoint_triples(&m, u, slab);
        boolean(&mut m, BoolKind::Cut, u, slab).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // Swapped, the same four nodes are the island's whole outer ring (cell 3f-2).
        // The arrangement does not know which solid is `a`, so it offers the same set;
        // the result has to still contain all of it. If the island face were dropped,
        // the box would lose its floor and every one of the four would go with it.
        let (mut m, l, stub) = l_and_dimple();
        let from_arrange = seam_endpoint_triples(&m, stub, l);
        assert_eq!(from_arrange.len(), 4);
        boolean(&mut m, BoolKind::Cut, stub, l).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));
    }

    /// The lid's hole survives a seam that crosses the lid's outer ring. `stitch_cycles`
    /// leaves one region — the lid with a corner bitten off — and the rim, being a ring of
    /// three-plane points, is placed inside it by the same containment test a seam loop
    /// gets. Measured before cell 3f-5: `Ok`, volume 0.96996 for a correct 0.916625, with
    /// `validate` reporting five violations.
    ///
    /// `0.92 − 0.15³`. A boolean composing on a shape a boolean can make.
    ///
    /// The box was asymmetric until cell (5a), because the fan was. `[0.85,1.15]³` pierces
    /// the lid at `(0.85, 0.85)`, dead on the diagonal of every fan the lid and its rim
    /// admit, and `contact_degenerate` spoke — honestly, but about the triangulation rather
    /// than about the geometry. Exact containment has no diagonals, so the coordinate comes
    /// back and brings with it the `0.916625` design.md §9 has called the right answer since
    /// before cell 3f-5.
    #[test]
    fn cut_a_pocket_at_a_corner() {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(Point3::from_array([0.85; 3]), Point3::from_array([1.15; 3]));
        let r = boolean(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 0.916625).abs() < 1e-9, "{}", props.volume);
        // A rectangular bite at a convex corner trades three quarter-faces for three more.
        assert!((props.area - 6.8).abs() < 1e-9, "{}", props.area);
        assert_eq!(holed_faces(&m, r).len(), 1, "the lid kept its hole");
    }

    /// The same bite, mirrored in `z`, so the lid never meets the seam and its hole rides
    /// out through `whole()` rather than through `place_loops`. Two paths, one number.
    #[test]
    fn cut_a_pocket_at_a_bottom_corner() {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.85, 0.85, -0.15]),
            Point3::from_array([1.15, 1.15, 0.15]),
        );
        let r = boolean(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 0.916625).abs() < 1e-9, "{}", props.volume);
        assert_eq!(holed_faces(&m, r).len(), 1, "the lid kept its hole");
    }

    /// The box crosses the pocket rim, so a boundary node rides a rim edge and `∂f` is no
    /// longer one ring. Rejected where the rule applies, in `edge_ix`, and named for what
    /// it is. The rim edge itself is fine now: it reaches this far only because
    /// `edge_incidence` gives it both incident planes and a seam vertex gets built on it.
    ///
    /// The obvious box, which threads the pocket void and leaves through a wall. Cell 3f-5
    /// had to hang a smaller one over the rim's corner instead, because `pierced_multi`
    /// spoke here; cell 3e-3 retired that guard and `contact_degenerate` spoke underneath
    /// it; cell (5a) retired that one too, and `SEAM_ACROSS_HOLE_RIM` the last — cell 3f-6.
    /// The box crosses the lid's pocket rim, so `∂(lid)` is two rings the seam threads into
    /// one corner notch; the pocket opens to the outside and the lid is left hole-free.
    ///
    /// `V = 0.92 − (cube∩box − pocket∩box) = 0.92 − (0.030375 − 0.003375) = 0.893`.
    #[test]
    fn cut_across_a_hole_rim() {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.55, 0.55, 0.85]),
            Point3::from_array([1.15; 3]),
        );
        let r = boolean(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 0.893).abs() < 1e-9, "{}", props.volume);
        // The absorbed rim leaves genus 0: no face carries an inner loop.
        assert!(holed_faces(&m, r).is_empty());
    }

    /// Same two solids, opposite order — a face of the *box* is arranged first, and the
    /// pocket's rim edge pierces it. The result must not depend on which operand is named
    /// first.
    ///
    /// Both orders also pin `seam_segments_on`'s hole-aware sweep. Revert it to walk outer
    /// rings only and the rim crossings go uncounted, the crossing count on the face turns
    /// odd, and `arrangement_degenerate` speaks in place of the honest one — measured, both
    /// ways. Detectors are not interchangeable.
    ///
    /// `V(box) − (cube∩box − pocket∩box) = 0.108 − 0.027 = 0.081`.
    #[test]
    fn cut_across_a_hole_rim_either_way() {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.55, 0.55, 0.85]),
            Point3::from_array([1.15; 3]),
        );
        let r = boolean(&mut m, BoolKind::Cut, bx, pc).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 0.081).abs() < 1e-9, "{}", props.volume);
        assert!(holed_faces(&m, r).is_empty());
    }

    /// An imprint is the one holed operand the convex path can see, and `has_coplanar_pair`
    /// catches it — the region face is coplanar with the lid it was cut from. This pins the
    /// tag that took over when the door guard came down.
    #[test]
    fn an_imprinted_convex_operand_rejects_as_coplanar_pair() {
        let (mut m, ic) = imprinted_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, ic, bx),
            tag::COPLANAR_PAIR,
        );
    }

    /// A slab over the pocketed cube, its underside at height `z0`. The rectangle is
    /// asymmetric so that the cube's four vertical edges, which pierce the underside at
    /// `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`, miss its fan diagonals; a square slab has all
    /// four apexes degenerate at once.
    ///
    /// `z0 = 0.3` runs below the pocket floor, so the slab's underside meets only the
    /// cube's outer walls. `z0 = 0.7` runs between the floor and the lid and meets the
    /// pocket walls as well — two loops, nested.
    fn pocket_and_slab(z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, pc) = pocketed_cube();
        let slab = m.add_cuboid(
            Point3::from_array([-0.2, -0.25, z0]),
            Point3::from_array([1.3, 1.2, 1.5]),
        );
        (m, slab, pc)
    }

    fn holed_faces(m: &Model, s: Handle<Solid>) -> Vec<Handle<Face>> {
        solid_faces(m, s)
            .into_iter()
            .filter(|&f| !m.faces.get(f).inner.is_empty())
            .collect()
    }

    /// `flip` + `inner`, finally exercised. The pocket's lid is strictly inside the slab
    /// and the seam misses it, so it rides out through `whole()` as a `Cut`'s inside-B
    /// piece: every ring reversed, orientation toggled, and the hole still a hole. Written
    /// since cell 3f-1, believed but never run — no operand could carry a hole in.
    ///
    /// Two holed faces come out: that lid at `z = 1`, and the slab's underside at
    /// `z = 0.3`, where the cube's cross-section is a seam loop placed as a hole. So the
    /// two roads to an inner loop — carried in, and discovered — meet on one solid.
    ///
    /// `V(slab) − V(B ∩ {z ≥ 0.3}) = 2.61 − (0.7 − 0.08)`.
    #[test]
    fn cut_slab_by_pocket() {
        let (mut m, slab, pc) = pocket_and_slab(0.3);
        let r = boolean(&mut m, BoolKind::Cut, slab, pc).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 1.99).abs() < 1e-9, "{}", props.volume);
        assert!((props.area - 15.03).abs() < 1e-9, "{}", props.area);

        let mut z: Vec<f64> = holed_faces(&m, r)
            .into_iter()
            .map(|f| {
                let l = &m.faces.get(f).outer.half_edges[0];
                m.vertices.get(m.edges.get(l.edge).bounds.unwrap()[0]).point[2]
            })
            .collect();
        z.sort_by(f64::total_cmp);
        assert_eq!(z, vec![0.3, 1.0]);
    }

    /// The union. `2.61 + 0.92 − 0.62`, closing inclusion–exclusion with the two cuts.
    /// Here the lid is inside the slab and dropped, so only the discovered hole survives.
    #[test]
    fn fuse_slab_and_pocket() {
        let (mut m, slab, pc) = pocket_and_slab(0.3);
        let r = boolean(&mut m, BoolKind::Fuse, slab, pc).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 2.91).abs() < 1e-9, "{}", props.volume);
        assert!((props.area - 12.63).abs() < 1e-9, "{}", props.area);
        assert_eq!(holed_faces(&m, r).len(), 1);
    }

    /// The other cut: `0.92 − 0.62 = 0.30`, a plain `1 × 1 × 0.3` box. The pocket is
    /// entirely above the slab's underside, so nothing of it survives — a strong check,
    /// because a hole leaking through here would show up in the volume at once.
    #[test]
    fn cut_pocket_by_slab() {
        let (mut m, slab, pc) = pocket_and_slab(0.3);
        let r = boolean(&mut m, BoolKind::Cut, pc, slab).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 0.3).abs() < 1e-9, "{}", props.volume);
        assert!((props.area - 3.2).abs() < 1e-9, "{}", props.area);
        assert!(holed_faces(&m, r).is_empty());
    }

    /// Nested loops, at last (cell 3f-7). Cell 3f-3 argued nesting was reachable and reached
    /// for a polyhedral torus; a pocket sliced between its floor and its lid is enough. The
    /// slab's underside carries the cube's cross-section as one loop and the pocket's inside
    /// it — a loop within a loop, which the pairwise guard over-rejected.
    ///
    /// `Cut` opens, both orders: `A ∩ B = B ∩ {z ≥ 0.7}` = `0.30 − 0.048 = 0.252`, slab `1.74`,
    /// pocketed cube `0.92`, so `Cut(slab,pc) = 1.488` and `Cut(pc,slab) = 0.668`. `Cut(pc,slab)`
    /// is the one that makes the cube cross-section an *island with a hole* (the slab is B, its
    /// underside dropped, so no region survives to own the loops); `Cut(slab,pc)` keeps the slab
    /// region and hangs the pocket loop in it as a hole beside an island.
    ///
    /// `Fuse` seals the pocket (`[0.3,0.7]² × [0.5,0.7]`, capped by the slab at `z = 0.7`) into
    /// an enclosed cavity — a second shell. Cell (5c) assembles it: the union material is
    /// `2.408` and the result carries one cavity of volume `0.032`, two shells, `validate`
    /// clean (the void's inward orientation). Both orders (Fuse is commutative).
    #[test]
    fn a_slab_between_the_lid_and_the_floor_nests_two_loops() {
        for (swap, expect) in [(false, 1.488), (true, 0.668)] {
            let (mut m, slab, pc) = pocket_and_slab(0.7);
            let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
            let r = boolean(&mut m, BoolKind::Cut, x, y).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "Cut swap={swap}: {vs:?}");
            let props = nacre_props::mass_props(&m, r).unwrap();
            assert!(
                (props.volume - expect).abs() < 1e-9,
                "Cut swap={swap}: {} vs {expect}",
                props.volume
            );
        }
        // Fuse seals the pocket into a cavity — a second shell, assembled by cell (5c).
        for swap in [false, true] {
            let (mut m, slab, pc) = pocket_and_slab(0.7);
            let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
            let r = boolean(&mut m, BoolKind::Fuse, x, y).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "Fuse swap={swap}: {vs:?}");
            let props = nacre_props::mass_props(&m, r).unwrap();
            assert!(
                (props.volume - 2.408).abs() < 1e-9,
                "Fuse swap={swap}: {}",
                props.volume
            );
            assert_eq!(m.solids.get(r).cavities.len(), 1, "Fuse swap={swap}");
            assert_eq!(m.reachable().shells.len(), 2, "Fuse swap={swap}");
        }
    }

    /// The sign convention `assemble_fuse_cut` classifies components by (cell 5c): an outward,
    /// material-enclosing shell has positive signed flux; an inward void shell (a cavity) has
    /// negative. `reversed_shell` flips one into the other.
    #[test]
    fn shell_signed_flux_is_positive_outward_and_negative_inward() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let outer = m.solids.get(cube).outer;
        let out_faces = m.shells.get(outer).faces.clone();
        assert!(shell_signed_flux(&m, &out_faces) > 0.0); // 3 · 8 = 24
        let void = m.reversed_shell(outer);
        let void_faces = m.shells.get(void).faces.clone();
        assert!(shell_signed_flux(&m, &void_faces) < 0.0);
    }

    /// `detect_coincident_interface` counts only cross-solid opposite-normal coplanar
    /// pairs, so imprinting a *different* face leaves the stack looking clean and routes
    /// into `coincident_merge` — the one path an imprinted operand can take. Measured
    /// before `solid_local_faces` learned to carry `face.inner`: `Ok`, volume 2.0533.
    #[test]
    fn coincident_merge_keeps_an_imprinted_hole() {
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
        let r = boolean(&mut m, BoolKind::Fuse, solid, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 2.0).abs() < 1e-9, "volume {vol}");
        assert!(
            solid_faces(&m, r)
                .into_iter()
                .any(|fh| !m.faces.get(fh).inner.is_empty()),
            "the imprinted face kept its hole"
        );
    }

    #[test]
    fn cut_u_by_slab() {
        // Two loops on one face, resolved. The slab's `y = 1.5` face has its whole boundary
        // dropped, so both loops are islands: the two prong cross-sections become the
        // result's new end caps, one input face giving two output faces, both flipped.
        //
        // `5.3 − 1.3`. Area 18: the clipped U's perimeter is 10, its caps 4 apiece.
        let (mut m, u, slab) = u_and_slab();
        let r = boolean(&mut m, BoolKind::Cut, u, slab).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 4.0).abs() < 1e-9, "volume {}", props.volume);
        assert!((props.area - 18.0).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn fuse_u_and_slab() {
        // The same face, the other way: `Fuse` keeps both outsides, so its boundary is kept
        // and both loops are *holes*. `5.3 + 8 − 1.3`.
        //
        // Area `28 + 18 − 2·2`: each island's 1.0 is buried on both sides of the interface.
        let (mut m, u, slab) = u_and_slab();
        let r = boolean(&mut m, BoolKind::Fuse, u, slab).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 12.0).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 42.0).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn cut_slab_by_u() {
        // Operands swapped: the holed face is now on **A**, and the U's caps split into two
        // cycles apiece under `flip`. Two blind pockets, `8 − 1.3`.
        //
        // Area `28 − 2 + (4·0.5 + 1) + (4·0.8 + 1)`: the face gives up its two 1.0 holes and
        // the pockets hand back walls and floors.
        let (mut m, u, slab) = u_and_slab();
        let r = boolean(&mut m, BoolKind::Cut, slab, u).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 6.7).abs() < 1e-9, "volume {}", props.volume);
        assert!((props.area - 33.2).abs() < 1e-9, "area {}", props.area);
    }

    /// The L with a box biting its reflex corner and poking out the top. The box top
    /// (z=1.2) clears the L's z=1 **deliberately**: sunk inside the L's slab, the L's
    /// vertical edges at (2,1) and (1,1) would pierce the box's bottom *and* top face,
    /// which `pierced_multi` used to reject. Cell 3e-3 supports it; the fixture keeps its
    /// clearance so that it goes on testing one thing.
    fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.2]),
            Point3::from_array([2.5, 1.5, 1.2]),
        );
        (m, l, bx)
    }

    #[test]
    fn cut_staircase_seam_arc() {
        // On the box's bottom face the seam runs (2,0.5) → (2,1) → (1,1) → (1,1.5): a
        // staircase whose two bends turn opposite ways. `strict` used to reject it,
        // unable to tell a reflex turn from an arc folded back on itself.
        //
        // The reconstructed face there is (0.5,0.5) → (0.5,1.5) → (1,1.5) → (1,1) →
        // (2,1) → (2,0.5): the overlap footprint, area 1.0, reflex at (1,1). A correct
        // simple polygon. `strict` was rejecting a right answer.
        //
        // Overlap = footprint 1.0 × z∈[0.2,1] = 0.8.
        let (mut m, l, bx) = l_and_popup_box();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.8)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn fuse_staircase_seam_arc() {
        // The `Fuse` counterpart, closing the inclusion–exclusion: V_L + V_box − 0.8.
        // `validate` cannot see a self-intersecting face (it stays manifold, Euler
        // holds), so the volume is what pins the folded arc — with OCCT alongside.
        let (mut m, l, bx) = l_and_popup_box();
        let r = boolean(&mut m, BoolKind::Fuse, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 2.0 - 0.8)).abs() < 1e-9, "volume {vol}");
    }

    /// The L-prism and an L-shaped bar lying in its notch, biting two convex corners of
    /// the L's top face. The bar spans `z ∈ [0.5, 1.5]`, so its body clears the cap.
    ///
    /// Each bite crosses **two different** edges of the cap, which is exactly why no edge
    /// is pierced twice — the bar takes corners, not edges. Two chords, no closed loop.
    fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bar = Profile2d {
            points: vec![
                p2(1.8, 0.8),
                p2(2.1, 0.8),
                p2(2.1, 2.1),
                p2(0.8, 2.1),
                p2(0.8, 1.8),
                p2(1.8, 1.8),
            ],
        };
        let OpOutput::Extrude { solid: b, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane {
                    origin: Point3::from_array([0.0, 0.0, 0.5]),
                    ..SketchPlane::world_xy()
                },
                profile: bar,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, b)
    }

    #[test]
    fn one_arc_stitches_to_the_old_splice() {
        // The single-arc case must come out `[[0]]`, because that cycle *is* today's
        // splice: run `dk+1 ..= kd`, then the arc. A hexagon with one chord.
        let kept = [true, true, true, false, false, false];
        assert_eq!(
            arrange::stitch_cycles(&kept, &[2], &[5], &ident_next(6)).unwrap(),
            vec![vec![0]]
        );
    }

    #[test]
    fn two_chords_make_one_ring_or_two() {
        // Which it is depends on nothing but where the kept vertices sit.
        //
        // The L's cap: `[T,T,F,T,F,T]`, arcs on transitions `(kd=1, dk=2)` and
        // `(kd=3, dk=4)`. From `dk=2` the run is vertex 3 alone and ends at `kd=3`, arc 1;
        // from `dk=4` the run is `5,0,1` and ends at `kd=1`, arc 0. One cycle, two arcs —
        // the kept region is the cap minus two corners.
        let cap = [true, true, false, true, false, true];
        assert_eq!(
            arrange::stitch_cycles(&cap, &[1, 3], &[2, 4], &ident_next(6)).unwrap(),
            vec![vec![0, 1]]
        );

        // The bar's floor seen from `Cut`'s B-piece: `[T,F,F,F,T,F]`, arcs `(kd=0, dk=5)`
        // and `(kd=4, dk=3)`. Each run is one vertex and closes on its own arc: two
        // cycles, two faces — the two bites.
        let floor = [true, false, false, false, true, false];
        assert_eq!(
            arrange::stitch_cycles(&floor, &[0, 4], &[5, 3], &ident_next(6)).unwrap(),
            vec![vec![0], vec![1]]
        );
    }

    #[test]
    fn cycles_start_at_the_lowest_unused_arc() {
        // Replay rests on the order of `LocalFace`s, which rests on this. Feed the same
        // two cycles with the arcs swapped and the output order swaps with them — it
        // follows `opens`, which `seam_paths_on` already orders deterministically.
        let floor = [true, false, false, false, true, false];
        assert_eq!(
            arrange::stitch_cycles(&floor, &[4, 0], &[3, 5], &ident_next(6)).unwrap(),
            vec![vec![0], vec![1]]
        );
    }

    #[test]
    fn a_broken_successor_map_is_rejected_not_looped() {
        // Where `classof`'s ray casting and the arrangement's exact crossings would
        // disagree. Today's single-arc code *assumes* all three; with several arcs the
        // assumption can be silently false, so it becomes the same reject that already
        // arbitrates those two machines.
        let cap = [true, true, false, true, false, true];
        // A `kd` that is not a kept→dropped transition.
        assert_rejects(
            || arrange::stitch_cycles(&cap, &[0, 3], &[2, 4], &ident_next(6)).map(|_| ()),
            tag::SEAM_COUNT_MISMATCH,
        );
        // Two arcs claiming the same `kd`.
        assert_rejects(
            || arrange::stitch_cycles(&cap, &[1, 1], &[2, 4], &ident_next(6)).map(|_| ()),
            tag::SEAM_COUNT_MISMATCH,
        );
        // Four transitions but only one arc: the `kd`s no longer cover them.
        assert_rejects(
            || arrange::stitch_cycles(&cap, &[1], &[2], &ident_next(6)).map(|_| ()),
            tag::SEAM_COUNT_MISMATCH,
        );
    }

    #[test]
    fn the_notch_bar_cuts_two_chords_in_one_face() {
        // What `multichord` guards on `l_and_notch_bar`, pinned before the guard moves.
        //
        // The bar bites two *corners* of the L's cap, so each of its two arcs ends on two
        // **different** edges of `∂f`. That is the whole trick: an arc with both ends on
        // one edge means that edge is pierced twice, which `pierced_multi` took before
        // `multichord` ever saw it (cell 3e-3 retired that guard, and the fixture keeps its
        // corners so it goes on testing one thing). Every multi-chord shape I tried first fell
        // into exactly that trap.
        //
        // The cap keeps its corners `(0,0) (2,0) (1,1) (0,2)` and drops `(2,1) (1,2)`.
        let (m, l, bar) = l_and_notch_bar();
        let (planes, surf_ix) = combined(&m, l, bar);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);

        let ps = paths(&m, top, l, bar, &planes, &surf_ix);
        assert_eq!(ps.len(), 2);
        assert!(ps.iter().all(|p| matches!(p, arrange::SeamPath::Open(_))));
        // Two boundary ends and one bend apiece; the bend is the bar's vertical edge
        // piercing the cap, and the ends ride the cap's own edges.
        for p in &ps {
            let nodes = p.nodes();
            assert_eq!(nodes.len(), 3);
            assert!(nodes[0].on_edge.is_some() && nodes[2].on_edge.is_some());
            assert!(nodes[1].on_edge.is_none());
        }

        let hes = &m.faces.get(top).outer.half_edges;
        let inside: Vec<bool> = hes
            .iter()
            .map(|&he| {
                point_in_solid(&m, m.vertices.get(he_start(&m, he)).point, bar).unwrap()
                    == Side::Inside
            })
            .collect();
        assert_eq!(inside, [false, false, true, false, true, false]);

        // Four transitions on four distinct edges: one crossing apiece, so the run model
        // and the old edge-indexed one agree here, which is the point of the fixture.
        let ts = transitions_on(&m, top, bar);
        assert_eq!(ts.len(), 4);
        let ends: BTreeSet<Handle<Edge>> = ps
            .iter()
            .flat_map(|p| p.nodes())
            .filter_map(|nd| nd.on_edge)
            .collect();
        assert_eq!(ends.len(), 4);
    }

    #[test]
    fn cut_notch_bar() {
        // Two chords on one face, resolved. The bar bites the cap's corners `(2,1)` and
        // `(1,2)`; the kept region is the cap minus both, one ring using both arcs.
        // `3 − 2·(0.2 · 0.2 · 0.5) = 2.96`.
        //
        // The area is unchanged at 14: a corner bite removes `0.04` of cap, `0.1` of the
        // `x=2` wall and `0.1` of the `y=1` wall, and hands back exactly those three as
        // the bar's own faces. So area alone would not have noticed the bite at all — the
        // volume and `validate` are what score it here.
        let (mut m, l, bar) = l_and_notch_bar();
        let r = boolean(&mut m, BoolKind::Cut, l, bar).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 2.96).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 14.0).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn fuse_notch_bar() {
        // The `Fuse` counterpart, closing inclusion–exclusion: `3 + 0.69 − 0.04`. The bar
        // measures `0.3·1.3 + 1.0·0.3` in section, `1.0` tall.
        //
        // Area `14 + 6.58 − 0.96`: the bar's own surface is `5.2·1.0 + 2·0.69`, and each
        // bite buries `0.24` on either side of the interface.
        let (mut m, l, bar) = l_and_notch_bar();
        let r = boolean(&mut m, BoolKind::Fuse, l, bar).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (3.0 + 0.69 - 0.04)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 19.62).abs() < 1e-9, "area {}", props.area);
    }

    /// The L-prism with an **L-shaped** stub standing wholly inside its top face,
    /// `z ∈ [0.5, 1.5]`. The seam on the cap is a closed loop with a reflex node — the
    /// suite's first non-convex inner loop, and the shape a winding must be read from.
    ///
    /// Its coordinates dodge the cap's fan diagonals from `(0,0)` (`y = x`, `y = x/2`,
    /// `y = 2x`), which `segment_crosses_face` would graze.
    fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let ell = Profile2d {
            points: vec![
                p2(0.2, 0.25),
                p2(0.85, 0.25),
                p2(0.85, 0.4),
                p2(0.35, 0.4), // reflex
                p2(0.35, 0.9),
                p2(0.2, 0.9),
            ],
        };
        let OpOutput::Extrude { solid: stub, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane {
                    origin: Point3::from_array([0.0, 0.0, 0.5]),
                    ..SketchPlane::world_xy()
                },
                profile: ell,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, stub)
    }

    /// The ring's nodes as points, in ring order — tests only, to name a node by where it is.
    fn ring_points(nodes: &[arrange::SeamEnd], ring: &[[usize; 3]]) -> Vec<[f64; 3]> {
        ring.iter()
            .map(|t| {
                nodes
                    .iter()
                    .find(|nd| nd.triple == *t)
                    .unwrap()
                    .point
                    .as_array()
            })
            .collect()
    }

    #[test]
    fn a_loop_that_winds_the_wrong_way_for_where_it_sits_is_rejected() {
        // Containment says hole or island; the winding must agree. Two of the four
        // combinations are contradictions, and neither `props`, OCCT, nor the mesh gate
        // would notice — only this, and `validate` downstream.
        assert!(check_loop_class(true, -1).is_ok(), "a hole runs clockwise");
        assert!(
            check_loop_class(false, 1).is_ok(),
            "an island runs counter-clockwise"
        );
        assert_rejects(|| check_loop_class(true, 1), tag::LOOP_CLASS_MISMATCH);
        assert_rejects(|| check_loop_class(false, -1), tag::LOOP_CLASS_MISMATCH);

        // On the real dimple: reversing the ring reverses the winding, so containment and
        // the ring stop agreeing. Nothing else in the kernel would notice.
        let (m, l, stub) = l_and_dimple();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let p = surf_ix[&top];
        let ps = paths(&m, top, l, stub, &planes, &surf_ix);
        let arrange::SeamPath::Closed(nodes) = &ps[0] else {
            panic!("closed loop")
        };
        let mut ring = arrange::orient_seam_loop(&planes, p, nodes, true).unwrap();
        ring.reverse();
        let w = arrange::loop_winding(&planes, p, &ring).unwrap();
        assert_rejects(|| check_loop_class(true, w), tag::LOOP_CLASS_MISMATCH);
    }

    /// The containment forest the placement rests on, measured directly (cell 3f-7). The
    /// staple cap's nested pair: the ring sits inside the cycle, the cycle inside nothing.
    /// Cells 3f-3/3f-4 reached for a polyhedral torus to fire nesting; the staple cap is a
    /// real pair of rings, one inside the other, and nesting is now supported, not rejected.
    #[test]
    fn nest_loops_reads_the_containment_forest() {
        let (planes, p, _, cycle, ring) = staple_cap_rings(true, &[4, 5, 0, 1, 2], false);
        let nest = nest_loops(&planes, p, &[], &[cycle.clone(), ring.clone()]).unwrap();
        assert_eq!(nest[0].containers, Vec::<usize>::new()); // cycle: contained by nothing
        assert_eq!(nest[1].containers, vec![0]); // ring: inside the cycle
        assert!(nest[0].region.is_none() && nest[1].region.is_none());
        // With the cycle as a region, the ring nests in no loop and that region holds it.
        let owned = nest_loops(&planes, p, std::slice::from_ref(&cycle), &[ring]).unwrap();
        assert_eq!(owned[0].containers, Vec::<usize>::new());
        assert_eq!(owned[0].region, Some(0));
    }

    /// `classify_nesting` on paper (cell 3f-7): the combinatorics, with no coordinate read, so
    /// every branch is a synthetic forest. Depth parity classes each seam loop, containment
    /// gives each hole its owner, and the winding is the independent third check.
    #[test]
    fn classify_nesting_places_the_forest() {
        let nst = |region: Option<usize>, containers: &[usize]| Nesting {
            region,
            containers: containers.to_vec(),
        };
        use Owner::{Island as IslandOf, Region};
        use Placement::*;

        // Depth 0, all four: a seam loop is a hole in a kept region and an island in a dropped
        // one; an `f` hole (index >= n_rings) is always a hole, owned if a region holds it and
        // discarded otherwise — never reclassified an island.
        let nest = [
            nst(Some(0), &[]),
            nst(None, &[]),
            nst(Some(0), &[]),
            nst(None, &[]),
        ];
        assert_eq!(
            classify_nesting(&nest, 2, &[-1, 1, -1, -1]).unwrap(),
            vec![Hole(Region(0)), Island, Hole(Region(0)), Dropped]
        );

        // Depth 1 in a kept region: hole C, island P inside it.
        let nest = [nst(Some(0), &[]), nst(Some(0), &[0])];
        assert_eq!(
            classify_nesting(&nest, 2, &[-1, 1]).unwrap(),
            vec![Hole(Region(0)), Island]
        );
        // Depth 1 in a dropped region (the `Cut(pc,slab)` order): island C, hole P owned by it.
        let nest = [nst(None, &[]), nst(None, &[0])];
        assert_eq!(
            classify_nesting(&nest, 2, &[1, -1]).unwrap(),
            vec![Island, Hole(IslandOf(0))]
        );

        // A deep concentric chain in a dropped region — island, hole, island, hole. The last
        // hole is owned by the *deepest* island (L2), not the outermost (L0). A kept region
        // cannot start this chain: a depth-0 loop in it is a hole, so deep islands live only in
        // a dropped one.
        let nest = [
            nst(None, &[]),
            nst(None, &[0]),
            nst(None, &[0, 1]),
            nst(None, &[0, 1, 2]),
        ];
        assert_eq!(
            classify_nesting(&nest, 4, &[1, -1, 1, -1]).unwrap(),
            vec![Island, Hole(IslandOf(0)), Island, Hole(IslandOf(2))]
        );

        // The same lone loop is an island as a seam ring but a discarded hole as an `f` hole —
        // `n_rings` is what withholds the parity from the hole.
        assert_eq!(
            classify_nesting(&[nst(None, &[])], 1, &[1]).unwrap(),
            vec![Island]
        );
        assert_eq!(
            classify_nesting(&[nst(None, &[])], 0, &[-1]).unwrap(),
            vec![Dropped]
        );

        // The winding disagreeing with the parity is the two machines in conflict.
        assert_rejects(
            || classify_nesting(&[nst(Some(0), &[])], 1, &[1]).map(|_| ()),
            tag::LOOP_CLASS_MISMATCH,
        );
    }

    #[test]
    fn a_hole_winds_clockwise_and_an_island_counter_clockwise() {
        // The dimple's square hole. Walked with the material on its left it runs clockwise
        // about the cap's `+z`, so the winding is `-1`. Flip `material_outside` — which is
        // exactly how cell 3f-2 reads the same face as an island — and it must be `+1`.
        let (m, l, stub) = l_and_dimple();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let p = surf_ix[&top];
        let ps = paths(&m, top, l, stub, &planes, &surf_ix);
        let arrange::SeamPath::Closed(nodes) = &ps[0] else {
            panic!("closed loop")
        };

        let hole = arrange::orient_seam_loop(&planes, p, nodes, true).unwrap();
        assert_eq!(arrange::loop_winding(&planes, p, &hole).unwrap(), -1);

        let island = arrange::orient_seam_loop(&planes, p, nodes, false).unwrap();
        assert_eq!(arrange::loop_winding(&planes, p, &island).unwrap(), 1);

        // Nothing but the ring's direction went into that. Reversing it by hand agrees.
        let mut reversed = hole.clone();
        reversed.reverse();
        assert_eq!(arrange::loop_winding(&planes, p, &reversed).unwrap(), 1);

        // A square turns the same way everywhere, so `nodes[0]` would have done. That is
        // precisely what the next test refutes.
        for i in 0..hole.len() {
            assert_eq!(arrange::turn_at(&planes, p, &hole, i).unwrap(), -1);
        }
    }

    #[test]
    fn a_reflex_node_turns_against_its_ring() {
        // The whole reason `loop_winding` hunts for a hull vertex. On the L-shaped hole the
        // turn is `-1` at five nodes and `+1` at the reflex one, `(0.35, 0.4)` — read
        // `turn_at` there and the ring looks counter-clockwise, which it is not.
        //
        // This is the `outer_tri` bug restated: the turn at one corner is the ring's winding
        // only when that corner is convex. There a fixture found it after the fact; here the
        // test finds it before there is any code to be wrong.
        let (m, l, stub) = l_and_ell_stub();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let p = surf_ix[&top];
        let ps = paths(&m, top, l, stub, &planes, &surf_ix);
        let arrange::SeamPath::Closed(nodes) = &ps[0] else {
            panic!("closed loop")
        };
        let hole = arrange::orient_seam_loop(&planes, p, nodes, true).unwrap();
        assert_eq!(arrange::loop_winding(&planes, p, &hole).unwrap(), -1);

        let pts = ring_points(nodes, &hole);
        let reflex = pts
            .iter()
            .position(|q| near(Point3::from_array(*q), [0.35, 0.4, 1.0]))
            .expect("the reflex node");
        assert_eq!(arrange::turn_at(&planes, p, &hole, reflex).unwrap(), 1);
        let turns: Vec<i8> = (0..hole.len())
            .map(|i| arrange::turn_at(&planes, p, &hole, i).unwrap())
            .collect();
        assert_eq!(turns.iter().filter(|&&t| t == 1).count(), 1);

        // And the node the search lands on is the lexicographically least, `(0.2, 0.25)` —
        // a hull vertex, where the turn is the winding. The test finds it by reading
        // coordinates; `loop_winding` finds it with an exact predicate.
        let lo = (0..pts.len())
            .min_by(|&i, &j| pts[i].partial_cmp(&pts[j]).unwrap())
            .unwrap();
        assert!(near(Point3::from_array(pts[lo]), [0.2, 0.25, 1.0]));
        assert_eq!(arrange::turn_at(&planes, p, &hole, lo).unwrap(), -1);
        assert_ne!(lo, reflex);

        // ★ The teeth. A ring is a cycle, so its winding cannot depend on where the walk
        // began. Start it at the reflex node and a `turn_at(ring[0])` implementation reads
        // `+1` — the exact fault `outer_tri` shipped. Measured: without this rotation, such
        // an implementation passes every assertion above.
        let mut rotated = hole.clone();
        rotated.rotate_left(reflex);
        assert_eq!(arrange::turn_at(&planes, p, &rotated, 0).unwrap(), 1);
        assert_eq!(arrange::loop_winding(&planes, p, &rotated).unwrap(), -1);
    }

    #[test]
    fn the_ell_stub_cuts_a_non_convex_loop() {
        // Every closed seam loop in the suite so far has been a rectangle, and a convex
        // ring turns the same way at every node. Cell 3h's hull-vertex search would never
        // be exercised by one. This loop has a reflex node, at `(0.35, 0.4)`.
        let (m, l, stub) = l_and_ell_stub();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let ps = paths(&m, top, l, stub, &planes, &surf_ix);
        assert_eq!(ps.len(), 1);
        let arrange::SeamPath::Closed(nodes) = &ps[0] else {
            panic!("closed loop")
        };
        assert_eq!(nodes.len(), 6);
        assert!(nodes.iter().all(|nd| nd.on_edge.is_none()));
    }

    #[test]
    fn cut_ell_dimple() {
        // The blind pocket is the stub's L-shaped section, `0.65·0.15 + 0.15·0.5 = 0.1725`,
        // half a unit deep. Area `14 − 0.1725 + 2.6·0.5 + 0.1725`: the lid gives up exactly
        // what the floor hands back, so only the walls move it.
        let (mut m, l, stub) = l_and_ell_stub();
        let r = boolean(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (3.0 - 0.1725 * 0.5)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 15.3).abs() < 1e-9, "area {}", props.area);
    }

    /// The L-prism and a П-shaped staple straddling the L's reflex corner. The profile
    /// lives in the **XZ** sketch plane and extrudes along `−y`, so the L's cap (`z = 1`)
    /// is *parallel* to the extrusion axis and the staple's section there falls into two
    /// pieces: one wholly inside the cap, one wrapping the corner `(1,1)`.
    ///
    /// That parallelism is the whole point. A prism cut by a plane **perpendicular** to
    /// its axis meets a face in the profile, which is connected — so every component of
    /// `profile ∩ f` reaches `∂f`, and a face can never carry both an arc and a loop. Every
    /// earlier attempt at such a fixture died on that.
    ///
    /// Leg bottoms sit at `z = 0.5` and `z = 0.45`: two coplanar faces of *one* operand
    /// trip `coplanar_pair` at the door, exactly as `u_prism`'s staggered prongs avoid.
    /// And the legs span `y ∈ [0.65, 1.3]`, not `[0.7, 1.3]`, because `(1.4, 0.7)` lies on
    /// the cap's fan diagonal `y = x/2` and `segment_crosses_face` would graze it.
    fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let staple = Profile2d {
            points: vec![
                p2(0.1, 0.5),
                p2(0.6, 0.5),
                p2(0.6, 1.3),
                p2(0.8, 1.3),
                p2(0.8, 0.45),
                p2(1.4, 0.45),
                p2(1.4, 1.5),
                p2(0.1, 1.5),
            ],
        };
        let OpOutput::Extrude { solid: st, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane {
                    origin: Point3::from_array([0.0, 1.3, 0.0]),
                    x_axis: Vector3::from_array([1.0, 0.0, 0.0]),
                    y_axis: Vector3::from_array([0.0, 0.0, 1.0]),
                },
                profile: staple,
                dist: 0.65,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, st)
    }

    /// The staple cap's `(∂f ring, cycle ring, oriented loop ring)`, as plane triples.
    ///
    /// The cycle is built the way `reconstruct_face_paths` builds it: `Cut(L, staple)` keeps
    /// `[T,T,T,F,T,T]`, so the run is `4,5,0,1,2` and the arc is spliced starting at its
    /// `kd` end (edge 2). `material_outside` picks the loop's direction.
    #[allow(clippy::type_complexity)]
    fn staple_cap_rings(
        material_outside: bool,
        kept_run: &[usize],
        arc_forward: bool,
    ) -> (
        Vec<PlaneInfo>,
        usize,
        Vec<[usize; 3]>,
        Vec<[usize; 3]>,
        Vec<[usize; 3]>,
    ) {
        let (m, l, st) = l_and_staple();
        let (planes, surf_ix) = combined(&m, l, st);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let p = surf_ix[&top];
        let ps = paths(&m, top, l, st, &planes, &surf_ix);
        let arrange::SeamPath::Open(arc) = &ps[0] else {
            panic!("an arc")
        };
        let arrange::SeamPath::Closed(lp) = &ps[1] else {
            panic!("a loop")
        };
        let inc_f = arrange::edge_planes(&m, l, &surf_ix).unwrap();
        let bnd = arrange::face_vertex_triples(&m, top, p, &inc_f).unwrap();

        let mut cycle: Vec<[usize; 3]> = kept_run.iter().map(|&i| bnd[i]).collect();
        if arc_forward {
            cycle.extend(arc.iter().map(|nd| nd.triple));
        } else {
            cycle.extend(arc.iter().rev().map(|nd| nd.triple));
        }
        let ring = arrange::orient_seam_loop(&planes, p, lp, material_outside).unwrap();
        (planes, p, bnd, cycle, ring)
    }

    #[test]
    fn a_loop_is_inside_the_face_it_was_found_on() {
        // `SeamPath::Closed` has claimed since cell 3f-1 that a closed seam loop never
        // touches `∂f`. Nothing checked it. Every node of the staple's loop is strictly
        // inside the cap's own ring, and `point_on_ring` is what would say otherwise.
        //
        // The `∂f` ring is a hexagon with a reflex corner, so this is not a convex test.
        let (planes, p, bnd, _, ring) = staple_cap_rings(true, &[4, 5, 0, 1, 2], false);
        for t in &ring {
            assert!(arrange::point_in_ring(&planes, p, *t, &bnd).unwrap());
        }
    }

    #[test]
    fn the_older_loops_are_inside_their_faces_too() {
        // The two fixtures cell 3f-4 must not move: a square hole and an L-shaped one. Both
        // rings sit strictly inside the same reflex hexagon, and every clear ray agrees.
        for (name, f) in [
            (
                "dimple",
                l_and_dimple as fn() -> (Model, Handle<Solid>, Handle<Solid>),
            ),
            ("ell stub", l_and_ell_stub),
        ] {
            let (m, l, stub) = f();
            let (planes, surf_ix) = combined(&m, l, stub);
            let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
            let p = surf_ix[&top];
            let ps = paths(&m, top, l, stub, &planes, &surf_ix);
            let arrange::SeamPath::Closed(lp) = &ps[0] else {
                panic!("{name}: a loop")
            };
            let inc_f = arrange::edge_planes(&m, l, &surf_ix).unwrap();
            let bnd = arrange::face_vertex_triples(&m, top, p, &inc_f).unwrap();
            let ring = arrange::orient_seam_loop(&planes, p, lp, true).unwrap();
            for t in &ring {
                let rays = arrange::every_ray(&planes, p, *t, &bnd).unwrap();
                assert!(!rays.is_empty(), "{name}: no clear ray");
                assert!(rays.iter().all(|&x| x), "{name}: {rays:?}");
            }
        }
    }

    #[test]
    fn a_loop_is_placed_by_where_it_is_not_by_how_it_winds() {
        // The shape cell 3f-3 could not decide. One face, one arc, one loop — and the loop
        // is a *hole* or an *island* depending only on which side of the arc it lies.
        //
        // `Cut(L, staple)`: the kept region is the cap minus the corner bite, and the loop
        // sits inside it. `Cut(staple, L)`: the kept region *is* the corner bite, and the
        // same loop sits outside it. Neither the winding nor `kept[0]` can tell those apart.
        let (planes, p, _, cycle, ring) = staple_cap_rings(true, &[4, 5, 0, 1, 2], false);
        for t in &ring {
            assert!(
                arrange::point_in_ring(&planes, p, *t, &cycle).unwrap(),
                "the loop is inside the kept region: a hole"
            );
        }
        // And the ring does not contain the region: containment is not symmetric.
        assert!(!arrange::point_in_ring(&planes, p, cycle[0], &ring).unwrap());
        assert!(!arrange::point_in_ring(&planes, p, cycle[7], &ring).unwrap());

        let (planes, p, _, cycle, ring) = staple_cap_rings(false, &[3], true);
        for t in &ring {
            assert!(
                !arrange::point_in_ring(&planes, p, *t, &cycle).unwrap(),
                "the loop is outside the kept region: an island"
            );
        }
    }

    #[test]
    fn every_clear_ray_agrees_and_half_of_them_are_not_clear() {
        // The ring is simple, so the parity cannot depend on which ray was cast. That is a
        // second machine, free.
        //
        // And the candidates are not interchangeable, which is why both of a node's planes
        // must be tried: the staple's legs are extruded from the same `y = 0.65` and
        // `y = 1.3` caps, so the loop and the arc *share* those planes. A ray along the
        // loop's `y` plane runs straight through the arc's nodes `(0.8, 0.65)` and
        // `(1.4, 0.65)`. Only the `x` lines survive — two of four candidates.
        let (planes, p, bnd, cycle, ring) = staple_cap_rings(true, &[4, 5, 0, 1, 2], false);

        let vs_face = arrange::every_ray(&planes, p, ring[0], &bnd).unwrap();
        assert_eq!(vs_face.len(), 4, "all four candidates are clear of `∂f`");
        assert!(vs_face.iter().all(|&x| x), "and all agree: inside");

        let vs_cycle = arrange::every_ray(&planes, p, ring[0], &cycle).unwrap();
        assert_eq!(
            vs_cycle.len(),
            2,
            "the loop's `y` plane meets the arc's nodes"
        );
        assert!(vs_cycle.iter().all(|&x| x), "the clear ones agree: inside");
    }

    #[test]
    fn a_ring_inside_a_ring_is_what_nesting_looks_like() {
        // `nested_loops` has no operand in the suite that produces it — a polyhedral torus
        // would. The detector can still be aimed at real geometry: the staple cap's kept
        // region contains its loop, so feeding that pair to the containment test as though
        // they were two loops fires exactly the condition cell 3f-4 rejects. It is the
        // detector under test, not the fixture.
        let (planes, p, _, cycle, ring) = staple_cap_rings(true, &[4, 5, 0, 1, 2], false);
        assert!(arrange::point_in_ring(&planes, p, ring[0], &cycle).unwrap());
        assert!(!arrange::point_in_ring(&planes, p, cycle[0], &ring).unwrap());
    }

    #[test]
    fn the_staple_leaves_an_arc_beside_a_loop() {
        // The one shape `pokehole` will still guard after cell 3f-3. The cap carries both:
        // the near leg cuts a rectangle wholly inside it, the far leg wraps the reflex
        // corner `(1,1)` and leaves an arc whose two ends ride *different* edges. That was
        // chosen to keep `pierced_multi` quiet; cell 3e-3 retired it, and the fixture stays
        // as it is so that it goes on testing one thing.
        let (m, l, st) = l_and_staple();
        let (planes, surf_ix) = combined(&m, l, st);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let ps = paths(&m, top, l, st, &planes, &surf_ix);
        assert_eq!(ps.len(), 2);

        let arrange::SeamPath::Open(arc) = &ps[0] else {
            panic!("an open arc")
        };
        assert_eq!(arc.len(), 5);
        let ends: BTreeSet<Handle<Edge>> = arc.iter().filter_map(|nd| nd.on_edge).collect();
        assert_eq!(ends.len(), 2, "the arc's ends ride two distinct edges");

        let arrange::SeamPath::Closed(loop_) = &ps[1] else {
            panic!("a closed loop")
        };
        assert_eq!(loop_.len(), 4);
        assert!(loop_.iter().all(|nd| nd.on_edge.is_none()));
    }

    #[test]
    fn cut_l_staple() {
        // A loop beside an arc, resolved. The cap's kept region is the hexagon minus the
        // corner bite, and the near leg's rectangle sits inside it — so the loop is a
        // **hole** of that region. The winding says clockwise, which is only a cross-check
        // now: containment decided.
        //
        // `V_∩ = 0.5·0.5·0.65 + 0.27·0.55 = 0.311`, the two legs' parts inside the L.
        let (mut m, l, st) = l_and_staple();
        let r = boolean(&mut m, BoolKind::Cut, l, st).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 2.689).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 15.755).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn fuse_l_staple() {
        // `3 + 0.7605 − 0.311`. The staple measures `1.17` in section, `0.65` deep.
        let (mut m, l, st) = l_and_staple();
        let r = boolean(&mut m, BoolKind::Fuse, l, st).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (3.0 + 0.7605 - 0.311)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 16.72).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn cut_staple_by_l() {
        // ★ The same loop, on the same face, is now an **island**. Swap the operands and the
        // cap's kept region becomes the corner bite alone; the loop lies outside it, in a
        // dropped region, so its interior is what survives.
        //
        // Cell 3f-3's counterexample made flesh: a hole and an island wind oppositely without
        // nesting, so no winding could have told these two apart. Position did.
        //
        // `0.7605 − 0.311`, and the three volumes close inclusion–exclusion exactly.
        let (mut m, l, st) = l_and_staple();
        let r = boolean(&mut m, BoolKind::Cut, st, l).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (0.7605 - 0.311)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 4.68).abs() < 1e-9, "area {}", props.area);
    }

    /// The L with a stub rising out of its top face, footprint strictly inside that
    /// face. Unlike the rod of `drill_through_the_l`, the stub enters the L from within,
    /// so each of its vertical edges crosses exactly one face and its bottom ring stays
    /// inside — one chord, no tunnel, and the seam loop is the whole story.
    fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let stub = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 0.5]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        (m, l, stub)
    }

    #[test]
    fn the_dimple_face_matches_the_shape_cell_3f_1_opens() {
        // What `pokehole` guards on `l_and_dimple` is exactly one face of the L, and
        // exactly one shape: no open arc, one closed loop, and the face's whole outer
        // boundary outside the stub. That last part is `kept[0] == true` written
        // without choosing a `keep` — for both `Cut` and `Fuse`, A keeps its outside.
        //
        // Crossing the seam flips the class, so a lone loop on an all-kept boundary
        // must bound a *dropped* interior: a hole, never an island. That is the whole
        // proof, and it needs no area.
        //
        // Read the same geometric fact with the operands swapped and it says the other
        // thing: with the L as B, `keep == Inside`, so every one of those boundary
        // vertices is *dropped* and the loop's interior is what survives — the island of
        // cell 3f-2. One measurement, two shapes; there is nothing further to measure.
        //
        // The arrangement's own preconditions are already measured — the five-fixture
        // sweep in `arrangement_agrees_with_todays_seam_bookkeeping` unwraps
        // `seam_paths_on` on every face of both solids.
        let (m, l, stub) = l_and_dimple();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);

        let out = paths(&m, top, l, stub, &planes, &surf_ix);
        assert_eq!(out.len(), 1);
        let arrange::SeamPath::Closed(nodes) = &out[0] else {
            panic!("the dimple's seam is a closed loop: {out:?}");
        };
        assert_eq!(nodes.len(), 4);

        for &he in &m.faces.get(top).outer.half_edges {
            let v = m.vertices.get(he_start(&m, he)).point;
            assert_eq!(point_in_solid(&m, v, stub).unwrap(), Side::Outside);
        }
    }

    #[test]
    fn a_seam_loop_keeps_material_on_its_left() {
        // The global sign, pinned before anything is wired. `loop_orient_mismatch`
        // cannot catch a loop that is flipped *as a whole* — every edge would be wrong
        // together — and `validate` is a caller's option, not the kernel's. So the
        // derivation is checked here, on the real fixture, against a sequence computed
        // by hand.
        //
        // The L's top face has `n_out = +z`; the hole is the stub's footprint. Walk it
        // with the material (outside the stub) on the left and you go clockwise seen
        // from +z. On the edge (0.3,0.3) → (0.3,0.7) the left is `z × y = −x`: outside
        // the stub. Only the test reads a coordinate; the rule reads three signs.
        //
        // The rule is local: it asks which side of the seam the kept material lies on,
        // never whether the loop bounds a hole or an island. So the second call below
        // is not a symmetry curiosity — `material_outside == false` is exactly the
        // island's outer ring, the same L face read with `keep == Inside` when the
        // operands are swapped. Cell 3f-2 wires it, and this golden is what pins its
        // direction: downstream, only `validate` and the signed mesh volume look.
        let (m, l, stub) = l_and_dimple();
        let (planes, surf_ix) = combined(&m, l, stub);
        let top = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, 1.0]);
        let p = surf_ix[&top];

        let out = paths(&m, top, l, stub, &planes, &surf_ix);
        let arrange::SeamPath::Closed(nodes) = &out[0] else {
            panic!("closed loop")
        };
        let point_of = |t: &[usize; 3]| nodes.iter().find(|nd| nd.triple == *t).unwrap().point;

        let outside = arrange::orient_seam_loop(&planes, p, nodes, true).unwrap();
        let got: Vec<[f64; 3]> = outside.iter().map(|t| point_of(t).as_array()).collect();
        assert_points_cycle(
            &got,
            &[
                [0.3, 0.3, 1.0],
                [0.3, 0.7, 1.0],
                [0.7, 0.7, 1.0],
                [0.7, 0.3, 1.0],
            ],
        );

        // The rule's only degree of freedom is that bit. Flip it and the loop must
        // reverse exactly — if it does not, the derivation is wrong somewhere.
        let inside = arrange::orient_seam_loop(&planes, p, nodes, false).unwrap();
        let mut rev = outside.clone();
        rev.reverse();
        assert_eq!(inside, rev);
    }

    /// Compare two closed sequences up to rotation (not reflection — the direction is
    /// the whole point).
    fn assert_points_cycle(got: &[[f64; 3]], want: &[[f64; 3]]) {
        assert_eq!(got.len(), want.len(), "{got:?} vs {want:?}");
        let hit = (0..want.len()).any(|r| {
            got.iter()
                .zip(want[r..].iter().chain(&want[..r]))
                .all(|(g, w)| near(Point3::from_array(*g), *w))
        });
        assert!(hit, "{got:?} vs {want:?} up to rotation");
    }

    #[test]
    fn cut_the_stub_by_the_l_leaves_an_island_face() {
        // Swap the operands of `cut_blind_dimple` and the same loop lands on a face whose
        // boundary is *all* dropped: the kept region is the loop's interior alone. That
        // face has no `∂f` at all — its outer loop *is* the seam ring, four `Discovered`
        // vertices and nothing else. The answer is the `0.4 × 0.4 × 0.5` box above `z = 1`.
        //
        // This is where `flip` first meets a ring that came from `orient_seam_loop`. The
        // rule wound it CCW about the L's `+z`, keeping the material (inside the stub) on
        // its left; `flip` reverses it and the face becomes the box's downward-facing
        // floor. Nothing but `validate` and the signed mesh volume can see that go wrong.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean(&mut m, BoolKind::Cut, stub, l).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.08).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn fuse_the_stub_and_the_l() {
        // Not a new branch — the same hole, on the same face of the L, reached with the
        // operands the other way round. `Fuse` keeps both outsides and flips neither, so
        // what this pins is that the answer does not depend on which solid is `a`: the
        // hole now lands on B, and 3.08 is 3.08.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean(&mut m, BoolKind::Fuse, stub, l).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 0.16 - 0.08)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn cut_blind_dimple() {
        // The stub's footprint never reaches the L's top-face boundary, so the seam is
        // a closed ring in the face interior. It is the face's inner loop, and the
        // result is a blind pocket: `3 − 0.4² × 0.5`.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.08)).abs() < 1e-9, "volume {vol}");
    }

    /// A slotted bar — a cuboid with a full-width groove cut across its top — has two
    /// coplanar top strips that share one Surface. Before cell coplanar-narrow the door
    /// guard rejected any second boolean on it as `COPLANAR_PAIR`; the strips are disjoint
    /// (they share no edge), so it now chains. Cut a blind pocket into one strip: the result
    /// is the bar minus the groove (`3 − 0.5`) minus the pocket (`0.3³`).
    #[test]
    fn a_slotted_bar_chains_through_a_cut() {
        let mut m = Model::new();
        let bar = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        // Oversized in y so the groove walls do not coincide with the bar's y-faces
        // (a flush cutter would be a cross-operand coplanar pair).
        let groove = m.add_cuboid(
            Point3::from_array([1.0, -0.5, 0.5]),
            Point3::from_array([2.0, 1.5, 1.5]),
        );
        let slotted = boolean(&mut m, BoolKind::Cut, bar, groove).unwrap();
        m.rebuild_adjacency();
        // The two z=1 strips are coplanar (share one Surface) but disjoint (no shared edge):
        // `has_coplanar_pair` sees the pair, the narrowed guard passes it.
        let planes = collect_planes(&m, slotted).unwrap();
        assert!(has_coplanar_pair(&planes));
        assert!(!solid_has_coplanar_neighbour_edge(&m, slotted, &planes));

        let pocket = m.add_cuboid(
            Point3::from_array([0.2, 0.2, 0.7]),
            Point3::from_array([0.5, 0.5, 1.5]),
        );
        let r = boolean(&mut m, BoolKind::Cut, slotted, pocket).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.5 - 0.027)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn a_flipped_hole_loop_is_caught() {
        // `loop_orient_mismatch` cannot see a loop flipped as a whole, and `f2`'s golden
        // pins the derivation — but the two downstream detectors must actually fire on
        // *this* shape, not merely exist. They are different in kind: `validate` sees a
        // rim edge used twice the same way, `tessellate` sees a hole wound like its
        // outer ring. Volume and OCCT see neither: `props` sums `|area|`.
        //
        // `Store` is append-only, so the face cannot be edited. Push a replacement with
        // the hole reversed, swap it into a fresh shell and solid, and move the live
        // handle: the old face falls out of `reachable()`.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();

        let faces = m.shells.get(m.solids.get(r).outer).faces.clone();
        let holed = *faces
            .iter()
            .find(|&&f| !m.faces.get(f).inner.is_empty())
            .expect("the L's top face carries the hole");
        let f = m.faces.get(holed).clone();
        let mut hole = f.inner[0].clone();
        hole.half_edges.reverse();
        for he in &mut hole.half_edges {
            he.forward = !he.forward;
        }
        let bad = m.faces.push(Face {
            inner: vec![hole],
            ..f
        });
        let swapped = faces
            .iter()
            .map(|&x| if x == holed { bad } else { x })
            .collect();
        let shell = m.shells.push(Shell { faces: swapped });
        let solid = m.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        m.live_solids = vec![solid];
        m.rebuild_adjacency();

        let vs = nacre_validate::validate(&m);
        assert!(
            vs.iter()
                .any(|v| matches!(v, nacre_validate::Violation::NonOpposedEdge { .. })),
            "validate stayed quiet: {vs:?}"
        );
        assert!(matches!(
            nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()),
            Err(nacre_tess::TessError::HoleWinding)
        ));
    }

    #[test]
    fn fuse_blind_dimple() {
        // The `Fuse` counterpart — a boss on the L — closing the inclusion–exclusion:
        // `V_L + V_stub − V_overlap`. Both put a hole in the same face of the L.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean(&mut m, BoolKind::Fuse, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 0.16 - 0.08)).abs() < 1e-9, "volume {vol}");
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

    /// The unit cube with a 0.4-square imprinted on its top face: no material moves,
    /// the face is merely split into a holed lid and a coplanar region face.
    fn imprinted_cube() -> (Model, Handle<Solid>) {
        let (mut m, top) = cube_with_top();
        let OpOutput::ImprintSketch { solid, .. } = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: top,
                profile: small_square(),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        (m, solid)
    }

    /// A face with an inner loop implies one of exactly two things about the shell: the
    /// rim is bounded by inward walls, making the solid non-convex (a pocket), or by a
    /// coplanar region face (an imprint). Nothing else closes. These two tests measure
    /// both halves, and together they say where a holed operand can arrive: an imprint is
    /// rejected as `COPLANAR_PAIR` (its region face is coplanar with the lid it was cut
    /// from); everything else takes the general seam path.
    #[test]
    fn a_pocketed_cube_is_not_convex() {
        let (m, pc) = pocketed_cube();
        assert!(!m.faces.get(solid_faces(&m, pc)[1]).inner.is_empty()); // the lid is holed
        let planes = collect_planes(&m, pc).unwrap();
        assert!(!is_convex(&m, &planes, &solid_vertex_handles(&m, pc)));
    }

    #[test]
    fn an_imprinted_cube_is_convex() {
        let (m, ic) = imprinted_cube();
        let planes = collect_planes(&m, ic).unwrap();
        assert!(is_convex(&m, &planes, &solid_vertex_handles(&m, ic)));
        assert!(has_coplanar_pair(&planes)); // the lid and its region face
        // The rim shares an edge, so the production guard also rejects it — unlike a slot's
        // two disjoint strips (cell coplanar-narrow).
        assert!(solid_has_coplanar_neighbour_edge(&m, ic, &planes));
    }

    /// A boolean *result* (here `Common` of two overlapping cubes) carries
    /// `Discovered` corners — the mixed-plane meets. This is the first operand
    /// that drives `is_convex`'s definition-based branch (cell 5d-3): those
    /// corners must be judged by their plane triple (`three_plane_orient3d`), not
    /// their coordinate cache, and the overlap box still reads convex.
    #[test]
    fn is_convex_judges_a_discovered_operand_by_its_triple() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 1.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let c = boolean(&mut m, BoolKind::Common, a, b).unwrap(); // the box [1,2]³
        m.rebuild_adjacency();
        let vhs = solid_vertex_handles(&m, c);
        // Six of the eight corners are mixed A/B-plane meets ⇒ Discovered; the
        // branch is genuinely exercised (the two original corners stay Constructed).
        let n_disc = vhs
            .iter()
            .filter(|&&vh| matches!(m.vertices.get(vh).origin, Origin::Discovered { .. }))
            .count();
        assert_eq!(n_disc, 6, "expected 6 Discovered corners, got {n_disc}");
        let planes = collect_planes(&m, c).unwrap();
        assert!(is_convex(&m, &planes, &vhs));
    }

    /// The first time a boolean result is fed back as an operand: the overlap box
    /// (Discovered corners) stacked on a third box merges through the coincident-
    /// interface path — which runs `is_convex` on that Discovered-cornered
    /// operand. Volume is the sum and the shell stays closed.
    #[test]
    fn a_boolean_result_stacks_as_an_operand() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 1.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let c = boolean(&mut m, BoolKind::Common, a, b).unwrap(); // [1,2]³
        m.rebuild_adjacency();
        let d = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 2.0]),
            Point3::from_array([2.0, 2.0, 3.0]),
        );
        let r = boolean(&mut m, BoolKind::Fuse, c, d).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 2.0).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0);
    }

    /// The pocketed cube with a second square imprinted on the pocket's floor: non-convex,
    /// so it takes the seam-free path, and carrying a coplanar pair (that floor and the
    /// region face cut from it), so exact containment cannot read its rings as triples.
    fn imprinted_pocketed_cube() -> (Model, Handle<Solid>) {
        let (mut m, top) = cube_with_top();
        let OpOutput::PocketOnFace { bottom_face, .. } =
            apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        let OpOutput::ImprintSketch { solid, .. } = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: bottom_face,
                profile: Profile2d {
                    points: vec![p2(-0.1, -0.1), p2(0.1, -0.1), p2(0.1, 0.1), p2(-0.1, 0.1)],
                },
            },
        )
        .unwrap() else {
            unreachable!()
        };
        (m, solid)
    }

    /// The price of exact containment, paid at the last door. `general_boolean` had no
    /// coplanar guard — `contained_result` reuses whole shells and never looks at a ring —
    /// so an imprinted non-convex operand used to sail through containment and disjointness.
    /// Now `edge_crosses_face` asks a face for its rings as three-plane triples, and an
    /// imprinted face has none: the rim's two neighbours are the same plane. Rejected, and
    /// `containment_boolean_already_keeps_a_pocket` measures that a *pocket* still rides
    /// through — the guard costs the imprint alone, not every hole.
    #[test]
    fn an_imprinted_nonconvex_operand_rejects_as_coplanar_pair() {
        let (mut m, ipc) = imprinted_pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.05, 0.05, 0.05]),
            Point3::from_array([0.15, 0.15, 0.15]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, ipc, bx),
            tag::COPLANAR_PAIR,
        );
    }

    /// A holed operand already survives the seam-free path — `general_boolean` never
    /// had a hole guard, and `contained_result` reuses whole shells, so the pocket rides
    /// through untouched. Nothing tested it. Pin it before the guard comes down.
    #[test]
    fn containment_boolean_already_keeps_a_pocket() {
        for (kind, want) in [(BoolKind::Cut, 0.919), (BoolKind::Fuse, 0.92)] {
            let (mut m, pc) = pocketed_cube();
            let bx = m.add_cuboid(
                Point3::from_array([0.05, 0.05, 0.05]),
                Point3::from_array([0.15, 0.15, 0.15]),
            );
            let r = boolean(&mut m, kind, pc, bx).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{kind:?} {vs:?}");
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            assert!((vol - want).abs() < 1e-9, "{kind:?} volume {vol}");
        }
    }

    /// How many of `x`'s edges pierce a face of `y`, by exact containment.
    fn crossings_of(
        m: &Model,
        x: Handle<Solid>,
        y: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Face>, usize>,
    ) -> Result<usize, BoolError> {
        let inc_y = arrange::edge_planes(m, y, surf_ix)?;
        let rings = solid_face_rings(m, y, surf_ix, &inc_y)?;
        let mut hits = 0;
        for (_, bounds, inc) in edge_incidence(m, x, surf_ix)? {
            let (p0, p1) = (
                m.vertices.get(bounds[0]).point,
                m.vertices.get(bounds[1]).point,
            );
            hits += pierced_faces(planes, inc, p0, p1, &rings)?.len();
        }
        Ok(hits)
    }

    /// The census: every edge against every face, on every fixture, both directions. The
    /// counts are the ones the fan gave before cell (5a) deleted it — the exact test
    /// reproduces them all, and adds `two_boxes`, which the fan could not answer at all.
    ///
    /// It also measures what cell (5a) risked. `face_rings` calls `face_vertex_triples` on
    /// every face now, not only the ones the seam touches, so a straight angle anywhere would
    /// reject; and `point_in_ring` needs a clear ray for every piercing point, inside the face
    /// or outside it. Neither fires. `no_clear_ray` and `loop_orient_mismatch` stay unfired.
    #[test]
    fn every_edge_against_every_face_on_every_fixture() {
        let extra = [
            ("cube_and_notch", (2, 4), cube_and_notch()),
            ("l_and_rod", (0, 8), l_and_rod()),
            ("pocket_and_slab(0.3)", (0, 4), pocket_and_slab(0.3)),
            ("pocket_and_slab(0.7)", (0, 8), pocket_and_slab(0.7)),
            ("two_boxes", (3, 3), two_boxes()),
        ];
        let want: HashMap<&str, (usize, usize)> = HashMap::from([
            ("l_and_corner_box", (3, 3)),
            ("l_and_reflex_box", (3, 5)),
            ("u_and_slab (two flat loops)", (8, 0)),
            ("l_and_popup_box (folded arc)", (4, 4)),
            ("l_and_dimple (inner loop)", (0, 4)),
            ("l_and_notch_bar (two chords)", (6, 6)),
            ("l_and_ell_stub (non-convex loop)", (0, 6)),
            ("l_and_staple (arc beside a loop)", (3, 9)),
        ]);
        for (name, expect, m, x, y) in overlap_fixtures()
            .into_iter()
            .map(|c| (c.name, want[c.name], c.m, c.x, c.y))
            .chain(extra.into_iter().map(|(n, e, (m, x, y))| (n, e, m, x, y)))
        {
            let (planes, surf_ix) = combined(&m, x, y);
            let got = (
                crossings_of(&m, x, y, &planes, &surf_ix).unwrap(),
                crossings_of(&m, y, x, &planes, &surf_ix).unwrap(),
            );
            assert_eq!(got, expect, "{name}");
        }
    }

    /// The one fixture the fan will not speak on, and the reason cell (5a) exists.
    ///
    /// `B`'s edge along `x` at `(y, z) = (0.5, 0.5)` meets `A`'s `x = 1` face at `(1, 0.5,
    /// 0.5)` — that square's **centre**, where both its diagonals cross. A fan from any apex
    /// lays a diagonal through it, so the fan grazes from all four and rejects. Exact
    /// containment does not care: the piercing point is `{B_y05, B_z05, A_x1}`, a three-plane
    /// point, and cell 3f-4 built its containment test for a different reason entirely.
    ///
    /// design.md §9 line 447 calls this the gate on unifying the convex path.
    #[test]
    fn the_centre_of_a_square_face_is_pierced_not_grazed() {
        let (m, a, b) = two_boxes();
        let (planes, surf_ix) = combined(&m, a, b);
        assert_eq!(crossings_of(&m, b, a, &planes, &surf_ix).unwrap(), 3);

        // Down to the one pair, so the point is named rather than counted.
        let x_face = face_facing(&m, a, &planes, &surf_ix, [1.0, 0.0, 0.0]);
        let q = surf_ix[&x_face];
        let inc_a = arrange::edge_planes(&m, a, &surf_ix).unwrap();
        let rings = arrange::face_rings(&m, x_face, q, &inc_a).unwrap();
        let (bounds, along_x) = *arrange::edge_planes(&m, b, &surf_ix)
            .unwrap()
            .values()
            .find(|(bd, _)| {
                bd.iter().all(|&v| {
                    let p = m.vertices.get(v).point.as_array();
                    p[1] == 0.5 && p[2] == 0.5
                })
            })
            .expect("B's edge through (0.5, 0.5)");
        let (p0, p1) = (
            m.vertices.get(bounds[0]).point,
            m.vertices.get(bounds[1]).point,
        );
        assert!(arrange::edge_crosses_face(&planes, along_x, p0, p1, q, &rings).unwrap());

        // And the whole boolean runs on it through the seam path — the measurement
        // design.md §9 line 447 asked for. Cell (5b) then deleted the convex path, so
        // `boolean` reaches the same code; this keeps the direct call as the §447 record.
        let (mut m, a, b) = two_boxes();
        let r = overlap_fuse_cut(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (1.0 - 0.125)).abs() < 1e-9, "volume {vol}");

        // The `Fuse` leg of the same corner overlap, which cell (5a) never measured on the
        // seam path — cell (5b) routes it here. `1 + 1 − 0.125`.
        let (mut m, a, b) = two_boxes();
        let r = overlap_fuse_cut(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.875).abs() < 1e-9, "volume {vol}");
    }

    /// A hole subtracts. An edge dropped straight through the pocket's mouth crosses the lid's
    /// plane inside its outer ring and inside its rim, so it pierces no material.
    #[test]
    fn an_edge_through_a_pocket_mouth_pierces_nothing() {
        let (mut m, pc) = pocketed_cube();
        let rod = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.6]),
            Point3::from_array([0.6, 0.6, 1.4]),
        );
        let (planes, surf_ix) = combined(&m, pc, rod);
        let inc_pc = arrange::edge_planes(&m, pc, &surf_ix).unwrap();
        let lid = solid_faces(&m, pc)
            .into_iter()
            .find(|&f| !m.faces.get(f).inner.is_empty())
            .unwrap();
        let q = surf_ix[&lid];
        let rings = arrange::face_rings(&m, lid, q, &inc_pc).unwrap();
        let mut vertical = 0;
        for (_, bounds, inc) in edge_incidence(&m, rod, &surf_ix).unwrap() {
            let (p0, p1) = (
                m.vertices.get(bounds[0]).point,
                m.vertices.get(bounds[1]).point,
            );
            if p0.as_array()[2] == p1.as_array()[2] {
                continue;
            }
            vertical += 1;
            assert!(!arrange::edge_crosses_face(&planes, inc, p0, p1, q, &rings).unwrap());
        }
        assert_eq!(vertical, 4);
    }

    /// The two inputs the convex `Fuse`/`Cut` path rejected but the seam path answers — both
    /// convex, both blocked only by `edge_seam`'s one-triple-per-edge, both already green on
    /// their non-convex twins (`cut_notch_bar`, `drill_through_the_l`, `fuse_the_l_and_the_rod`).
    /// Cell (5b) deleted the convex path, so `boolean` routes these here now.
    ///
    /// The drilled cube is genus 1, and `validate` clean does not prove that (a tunnel is
    /// clean too). Its two holed caps — the drill's entry and exit — do: `holed_faces == 2`.
    #[test]
    fn the_seam_path_answers_both_convex_pokes() {
        // A notch bitten out of one edge: the edge is crossed twice, so it threads the notch
        // rather than straddling. `10³ − 4·1.4·1.2` cut, `+2·4·1.4·1.2 − 6.72` fused.
        let (mut m, a, y) = cube_and_notch();
        let cut = boolean(&mut m, BoolKind::Cut, a, y).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let v = nacre_props::mass_props(&m, cut).unwrap().volume;
        assert!((v - 993.28).abs() < 1e-9, "notch cut {v}");

        let (mut m, a, y) = cube_and_notch();
        let fuse = boolean(&mut m, BoolKind::Fuse, a, y).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let v = nacre_props::mass_props(&m, fuse).unwrap().volume;
        assert!((v - 1014.4).abs() < 1e-9, "notch fuse {v}");

        // A bar straight through a cube: a genus-1 solid, and the canonical thing a user
        // drills. `27 − 1·1·3` cut, `27 + (5·1·1) − 3` fused.
        let drill = |kind| {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
            let bar = m.add_cuboid(
                Point3::from_array([1.0, 1.0, -1.0]),
                Point3::from_array([2.0, 2.0, 4.0]),
            );
            let r = boolean(&mut m, kind, a, bar).unwrap();
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
            (m, r)
        };
        let (m, cut) = drill(BoolKind::Cut);
        assert_eq!(holed_faces(&m, cut).len(), 2, "the drill's two caps");
        assert!((nacre_props::mass_props(&m, cut).unwrap().volume - 24.0).abs() < 1e-9);
        let (m, fuse) = drill(BoolKind::Fuse);
        assert!((nacre_props::mass_props(&m, fuse).unwrap().volume - 29.0).abs() < 1e-9);
    }

    /// An edge of one convex solid, threading the other, severs it. The bar runs through the
    /// cube and out both ends, so `Cut(bar, cube)` leaves the bar in two pieces — two solids,
    /// which one handle cannot answer. The convex path rejected this as `poke_through`; the
    /// seam path names it for what it is. First convex firing of `disconnected_result`
    /// (cell 3e-3's was the non-convex `Cut(rod, L)`).
    #[test]
    fn a_convex_cut_can_sever_its_operand() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, bar, a),
            tag::DISCONNECTED_RESULT,
        );
    }

    /// The two contacts, named. Both were `contact_degenerate` — "grazed a fan diagonal" — and
    /// neither ever was a graze.
    ///
    /// The vertex case is why `point_on_ring` is asked **before** `point_in_ring`: when the
    /// piercing point is a ring node, both of its rays carry that node, every candidate is
    /// skipped, and `point_in_ring` alone comes back `no_clear_ray`. Measured below.
    #[test]
    fn an_edge_touching_a_face_is_a_contact_not_a_graze() {
        let (mut m, l) = l_prism();
        // Its floor lies in the L's own bottom plane.
        let flat = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.5, 1.5, 0.5]),
        );
        // A corner exactly on the L cap's reflex vertex `(1,1)`, and one on an edge's interior.
        let at_vertex = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -0.5]),
            Point3::from_array([1.2, 1.2, 0.5]),
        );
        let at_edge = m.add_cuboid(
            Point3::from_array([1.5, 0.0, -0.5]),
            Point3::from_array([1.7, 0.2, 0.5]),
        );

        let pierce = |other: Handle<Solid>, want: [f64; 2]| {
            let (planes, surf_ix) = combined(&m, l, other);
            let inc_l = arrange::edge_planes(&m, l, &surf_ix).unwrap();
            let bottom = face_facing(&m, l, &planes, &surf_ix, [0.0, 0.0, -1.0]);
            let q = surf_ix[&bottom];
            let rings = arrange::face_rings(&m, bottom, q, &inc_l).unwrap();
            let (bounds, inc) = *arrange::edge_planes(&m, other, &surf_ix)
                .unwrap()
                .values()
                .find(|(bd, _)| {
                    bd.iter().all(|&v| {
                        let p = m.vertices.get(v).point.as_array();
                        [p[0], p[1]] == want
                    }) && bd
                        .iter()
                        .any(|&v| m.vertices.get(v).point.as_array()[2] < 0.0)
                })
                .expect("the vertical edge");
            let (p0, p1) = (
                m.vertices.get(bounds[0]).point,
                m.vertices.get(bounds[1]).point,
            );
            (planes, q, rings, inc, p0, p1)
        };

        // (a) The whole edge lies in the other face's plane. Both endpoint signs are zero.
        {
            let (planes, surf_ix) = combined(&m, l, flat);
            let inc_flat = arrange::edge_planes(&m, flat, &surf_ix).unwrap();
            let floor = face_facing(&m, flat, &planes, &surf_ix, [0.0, 0.0, -1.0]);
            let q = surf_ix[&floor];
            let rings = arrange::face_rings(&m, floor, q, &inc_flat).unwrap();
            let (bounds, inc) = *arrange::edge_planes(&m, l, &surf_ix)
                .unwrap()
                .values()
                .find(|(bd, _)| {
                    bd.iter()
                        .all(|&v| m.vertices.get(v).point.as_array()[2] == 0.0)
                })
                .expect("an edge of the L's bottom");
            let (p0, p1) = (
                m.vertices.get(bounds[0]).point,
                m.vertices.get(bounds[1]).point,
            );
            assert_rejects(
                || arrange::edge_crosses_face(&planes, inc, p0, p1, q, &rings),
                tag::VERTEX_ON_FACE_PLANE,
            );
        }

        // (b) The piercing point is a ring node.
        {
            let (planes, q, rings, inc, p0, p1) = pierce(at_vertex, [1.0, 1.0]);
            assert_rejects(
                || arrange::edge_crosses_face(&planes, inc, p0, p1, q, &rings),
                tag::POINT_ON_RING,
            );
            let mut x = [inc[0], inc[1], q];
            x.sort_unstable();
            assert_rejects(
                || arrange::point_in_ring(&planes, q, x, &rings[0]),
                tag::NO_CLEAR_RAY, // the tag the pre-check exists to prevent
            );
        }

        // (c) The piercing point is inside a ring edge. `point_in_ring` names this one itself.
        {
            let (planes, q, rings, inc, p0, p1) = pierce(at_edge, [1.5, 0.0]);
            assert_rejects(
                || arrange::edge_crosses_face(&planes, inc, p0, p1, q, &rings),
                tag::POINT_ON_RING,
            );
        }
    }

    /// The pocketed cube's planes, its lid, that lid's plane index, and its edge index.
    fn pocket_lid_setup() -> (
        Model,
        Vec<PlaneInfo>,
        Handle<Face>,
        usize,
        arrange::EdgePlanes,
    ) {
        let (m, pc) = pocketed_cube();
        let planes = collect_planes(&m, pc).unwrap();
        let surf_ix: HashMap<Handle<Face>, usize> = planes
            .iter()
            .enumerate()
            .map(|(i, pi)| (pi.face, i))
            .collect();
        // The lid and the pocket floor share an outward normal, so pick by the hole.
        let lid = solid_faces(&m, pc)
            .into_iter()
            .find(|&fh| !m.faces.get(fh).inner.is_empty())
            .expect("a holed face");
        let p = surf_ix[&lid];
        let inc = arrange::edge_planes(&m, pc, &surf_ix).unwrap();
        (m, planes, lid, p, inc)
    }

    /// A hole rim is a ring of three-plane points, exactly like a seam loop: the lid's
    /// plane and the two pocket walls meeting at each rim vertex. So `place_loops` can
    /// place it with no new machinery — cell 3f-5 just appends it to the loop list.
    ///
    /// The winding is the load-bearing measurement. `f.inner` is stored clockwise about
    /// the face's *outward* normal, and `loop_winding` reads `orient_sign(P)`; whether
    /// those two conventions agree was assumed, never measured. They do: `−1`, the value
    /// `check_loop_class(is_hole = true, ·)` demands, and it does not depend on `flip`.
    #[test]
    fn a_pocket_rim_is_a_clockwise_ring_of_three_plane_points() {
        let (m, planes, lid, p, inc) = pocket_lid_setup();
        let rims = arrange::hole_rings(&m, lid, p, &inc).unwrap();
        assert_eq!(rims.len(), 1);
        let rim = &rims[0];
        assert_eq!(rim.len(), 4);
        for t in rim {
            assert!(t.contains(&p), "every rim node lies on the lid");
            // Its other two planes are pocket walls: neither is the lid, and each pair
            // is distinct (a straight angle would already have rejected).
            assert_eq!(t.iter().filter(|&&x| x != p).count(), 2);
        }
        assert_eq!(arrange::loop_winding(&planes, p, rim).unwrap(), -1);
    }

    /// The rim sits inside the lid's outer ring, and the outer ring's vertices sit outside
    /// the rim. Both are ordinary `point_in_ring` questions once the rim is triples.
    #[test]
    fn a_pocket_rim_is_inside_its_lid_and_the_lid_is_not_inside_it() {
        let (m, planes, lid, p, inc) = pocket_lid_setup();
        let rim = &arrange::hole_rings(&m, lid, p, &inc).unwrap()[0];
        let bnd = arrange::face_vertex_triples(&m, lid, p, &inc).unwrap();

        for t in rim {
            assert!(arrange::point_in_ring(&planes, p, *t, &bnd).unwrap());
            // Every clear ray must agree — the ring is simple. A second machine, free.
            let rays = arrange::every_ray(&planes, p, *t, &bnd).unwrap();
            assert!(!rays.is_empty() && rays.iter().all(|&r| r), "{rays:?}");
        }
        for t in &bnd {
            assert!(!arrange::point_in_ring(&planes, p, *t, rim).unwrap());
        }
    }

    /// A lid's hole rim, placed by `nest_loops`: the lid's outer boundary is its region, and
    /// the rim nests in no loop, so it comes out owned by region `0` and depth `0`.
    #[test]
    fn nest_loops_owns_a_hole_ring() {
        let (m, planes, lid, p, inc) = pocket_lid_setup();
        let rim = arrange::hole_rings(&m, lid, p, &inc).unwrap().remove(0);
        let bnd = arrange::face_vertex_triples(&m, lid, p, &inc).unwrap();
        let nest = nest_loops(&planes, p, &[bnd], &[rim]).unwrap();
        assert_eq!(nest[0].region, Some(0));
        assert_eq!(nest[0].containers, Vec::<usize>::new());
        // A region ring the seam has bitten a corner from is exercised end to end by
        // `cut_a_pocket_at_a_corner`; it cannot be faked here, because dropping a node
        // from `bnd` leaves two nodes sharing only the lid's plane and `ring_edge`
        // rightly refuses to invent an edge between them.
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
    fn cut_with_a_hollow_operand_far_from_the_void() {
        // A cavitied operand whose seam misses the void: the corner cut is far from
        // the [1,2]³ void, so the void is carried through and preserved (cell (5c-in)).
        // The seam front-end walks all shells, so the void no longer silently vanishes
        // (it used to read as convex and return vol 26.875, cavities 0). Now: correct
        // 25.875 (27 − 1 void − 0.125 corner) with the cavity intact.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        m.rebuild_adjacency();
        let cutter = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
        let r = boolean(&mut m, BoolKind::Cut, hollow, cutter).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
    }

    #[test]
    fn blind_hole_drills_into_a_void() {
        // A cut reaching *into* a void: a stub drilled from below the hollow box up into
        // its void. The void loses its enclosure and merges with the outer shell (an
        // open pocket, cavities 0). The seam machinery (all-shells) + the (5c) component
        // split reconstruct it exactly: remove the channel [1.4,1.6]²×[0,1]=0.04 through
        // the floor, void interior removes nothing ⇒ 26 − 0.04 = 25.96.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        m.rebuild_adjacency();
        let stub = m.add_cuboid(
            Point3::from_array([1.4, 1.4, -0.5]),
            Point3::from_array([1.6, 1.6, 1.5]),
        );
        let r = boolean(&mut m, BoolKind::Cut, hollow, stub).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.96).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0); // the void opened to outside
    }

    #[test]
    fn a_tunnel_drilled_through_a_void() {
        // A cut passing all the way through the void (bottom to top). Removes the floor
        // and ceiling channels [1.4,1.6]²×([0,1]∪[2,3]) = 0.08; the void interior removes
        // nothing ⇒ 26 − 0.08 = 25.92. Result is a genus-1 solid (a straight tunnel),
        // one shell, no cavity.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let tunnel = m.add_cuboid(
            Point3::from_array([1.4, 1.4, -0.5]),
            Point3::from_array([1.6, 1.6, 3.5]),
        );
        let r = boolean(&mut m, BoolKind::Cut, hollow, tunnel).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.92).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0);
    }

    #[test]
    fn a_slab_that_splits_a_hollow_box_disconnects() {
        // A slab cut through the whole box (and its void) severs it into two solids —
        // honestly rejected as `DISCONNECTED_RESULT`, not silently mis-assembled.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.4, -0.5]),
            Point3::from_array([3.5, 1.6, 3.5]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, hollow, slab),
            tag::DISCONNECTED_RESULT,
        );
    }

    #[test]
    fn a_hollow_part_takes_a_second_far_cut() {
        // The headline: keep cutting a part after it is hollow. A bore at the corner
        // opposite the void — the void survives, and the result is fed back as an
        // operand (chaining past the first cavity-producing op, which the door used to
        // block).
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let bore = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let r = boolean(&mut m, BoolKind::Cut, hollow, bore).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
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
    fn a_profile_reaching_past_the_face_is_rejected() {
        // The frame origin is the top face centroid (0.5, 0.5). Each profile reaches past
        // the [0,1]² face region, so all three face ops reject rather than build a
        // silently-invalid inner loop (cell imprint-containment; n0 measured that apply was
        // Ok, validate passed, props integrated garbage — including a negative volume — and
        // only tessellate's NoEar caught it). Boundary contact ("touches") rejects too.
        let crosses = Profile2d {
            points: vec![p2(-0.7, -0.2), p2(0.7, -0.2), p2(0.7, 0.2), p2(-0.7, 0.2)],
        };
        let outside = Profile2d {
            points: vec![p2(1.8, 1.8), p2(2.2, 1.8), p2(2.2, 2.2), p2(1.8, 2.2)],
        };
        let overhang = Profile2d {
            points: vec![p2(-2.0, -2.0), p2(2.0, -2.0), p2(2.0, 2.0), p2(-2.0, 2.0)],
        };
        let touches = Profile2d {
            points: vec![p2(-0.5, -0.2), p2(0.3, -0.2), p2(0.3, 0.2), p2(-0.5, 0.2)],
        };
        for (name, profile) in [
            ("crosses", crosses),
            ("outside", outside),
            ("overhang", overhang),
            ("touches", touches),
        ] {
            for kind in ["imprint", "pad", "pocket"] {
                let (mut m, top) = cube_with_top();
                let op = match kind {
                    "imprint" => Operation::ImprintSketch {
                        face: top,
                        profile: profile.clone(),
                    },
                    "pad" => Operation::PadOnFace {
                        face: top,
                        profile: profile.clone(),
                        dist: 0.3,
                    },
                    _ => Operation::PocketOnFace {
                        face: top,
                        profile: profile.clone(),
                        dist: 0.3,
                    },
                };
                assert_eq!(
                    apply(&mut m, &op),
                    Err(OpError::ProfileNotContainedInFace),
                    "{name}/{kind}"
                );
            }
        }
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
    fn fuse_a_boss_onto_a_face() {
        // A small box sits on the base's top face inside its boundary (contained footprint).
        // Coplanar contact: before this cell it rejected as vertex_on_face_plane. The base
        // top face gains the boss footprint as a hole, the boss contact face is dropped, and
        // the footprint edges stitch the two. Volume 1 + 0.5²·1 = 1.25.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 1.0]),
            Point3::from_array([0.75, 0.75, 2.0]),
        );
        let r = boolean(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.25).abs() < 1e-12, "volume {vol}");
        // No coplanar-adjacent faces (the boss walls are perpendicular to the base top), so
        // the bossed solid chains into a further boolean.
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn fuse_stacked_cubes() {
        // A=[0,1]³ and B=[0,1]²×[1,2] share the z=1 face ⇒ merge into a clean 1×1×2 box. The
        // four coplanar side pairs are spliced into 4 faces, and the four interface corners
        // (each a straight angle on a split vertical edge) are dissolved, fusing the split
        // edges — a canonical 6-face / 8-vertex / 12-edge box (cell fuse-coplanar-merge).
        let (mut m, a, b) = stacked_cubes();
        let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let reach = m.reachable();
        assert_eq!(reach.faces.len(), 6);
        assert_eq!(reach.vertices.len(), 8);
        assert_eq!(reach.edges.len(), 12);
        // No coplanar-adjacent faces and no straight angles, so the fused solid chains.
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 2.0).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn a_fused_stack_chains_through_a_cut() {
        // The dissolved 1×1×2 box (cell fuse-coplanar-merge) feeds a second boolean. Before
        // the merge/dissolve this rejected — first as COPLANAR_PAIR (the flat edges), then as
        // LOOP_ORIENT_MISMATCH (the straight-angle interface corners). A clean box cuts.
        let (mut m, a, b) = stacked_cubes();
        let stack = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        // A cutter straddling z=1 (the fused interface) — the seam runs where the split
        // vertical edges used to be. Result: 2 − 0.5·0.5·1.0.
        let cutter = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let r = boolean(&mut m, BoolKind::Cut, stack, cutter).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.75).abs() < 1e-12, "volume {vol}");
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
        // Cell (5b) deleted `fuse_cut`, whose combined-plane coplanar check caught this;
        // the seam path's per-operand check passes and `edge_crosses_face` finds A's top
        // edge lying in B's z=1 plane, so the tag moves to `vertex_on_face_plane`. Still
        // rejected, one plane shy of honest (§9 line 451): the cause is two coplanar faces.
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, a, b),
            tag::VERTEX_ON_FACE_PLANE,
        );
    }

    #[test]
    fn same_ground_overlap_is_unsupported() {
        // Two boxes sharing the z=0 ground with overlapping footprints: only
        // same-normal coplanar contact + 3D overlap ⇒ needs the general 2D
        // coplanar path (next unit) ⇒ Unsupported.
        //
        // Until cell (5b) this was `coplanar_pair`, from `fuse_cut`'s combined-plane check.
        // With that path gone, the seam path's per-operand check passes and `plane_side`
        // meets A's floor edge on B's z=0 plane first — `vertex_on_face_plane`. Still an
        // honest rejection, less specific about why (§9 line 451). One of two fixtures that
        // moved this way; `coincident_merge_rejects_offset_footprint` is the other.
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.0]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, a, b),
            tag::VERTEX_ON_FACE_PLANE,
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
        // half-spaces. It never reaches `common`'s `is_convex` guard, though: the box's
        // floor is **coplanar with the L's**, so the L's bottom edges lie in that plane and
        // `boundaries_intersect` bails on a tangential contact first. Measured, not assumed.
        //
        // The fan used to call this `contact_degenerate` — "grazed a diagonal from every
        // apex". It never was a graze. Cell (5a) reads the endpoints' exact sides of the
        // plane, finds both zero, and says what is actually true.
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
            tag::VERTEX_ON_FACE_PLANE,
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
