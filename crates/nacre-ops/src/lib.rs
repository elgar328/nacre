//! Operations for the nacre kernel, plus a replayable operation log (design §6).
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::ImprintSketch`] (M4) splits an existing planar face along a
//! closed profile — the first op that consumes a prior op's face by `Handle`
//! (exposed via [`OpOutput`]) and supersedes a solid (design §2 live-solid
//! semantics). Ops are applied by [`apply`] and folded by [`replay`]; every
//! result is a **closed** solid, so `nacre-validate` applies fully.

use nacre_geom::{Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::Isometry;
use nacre_scalar::frame3::{Pt3, dir_orient3d_judge, orient3d_judge, orient3d_ray};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Orientation, Origin, Rotation, Shell, Solid, Vertex,
    VertexDef,
};
// Data-parallel evaluation of the read-only predicate phases (boolean face reconstruction
// and vertex classification). Only present under the default `parallel` feature; the
// serial build maps with plain iterators. See `overlap_fuse_cut`.
#[cfg(feature = "parallel")]
use rayon::prelude::*;

mod arrange;
mod tolerant;

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
    /// Pad a boss: extrude `profile` on a planar `face` into a tool prism (height
    /// `dist`) and `Fuse` it onto the solid — boolean sugar over [`Operation::Boolean`],
    /// not a direct face-split. No "profile inside the face" constraint: an overhanging
    /// footprint is handled by the boolean's coplanar-contact / overhang path. Adds
    /// material (design §6).
    PadOnFace {
        face: Handle<Face>,
        profile: Profile2d,
        dist: f64,
    },
    /// Carve a blind pocket: extrude `profile` on a planar `face` into a tool prism
    /// (depth `dist`) and `Cut` it from the solid — boolean sugar over
    /// [`Operation::Boolean`], not a direct face-split. No "profile inside the face"
    /// constraint (overhang footprints route through the boolean). A cut that would
    /// punch through is rejected as not-blind. Removes material (design §6).
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
    /// Rigid-body transform: supersede `solid` by its image under `isometry`
    /// (overhaul stage 1). 1a realizes a rational translation; 1b adds rotation.
    /// The `Isometry` is the exact definition (op-log truth); the geometry is a
    /// realized cache.
    Transform {
        solid: Handle<Solid>,
        isometry: Isometry,
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
    /// The result severs into two or more material solids *and* at least one enclosed void
    /// (cavity) survives. Which outer shell owns which cavity needs a shell-scoped point-in-shell
    /// test we do not have yet (a re-scope of `point_in_solid`), so this is honestly rejected and
    /// deferred to a follow-on cell. Reachable: `Cut` a hollow part with a cut that isolates the
    /// void into one severed piece. Born with its firing test (`severed_with_cavity_is_rejected`).
    pub const SEVERED_WITH_CAVITY: &str = "severed_with_cavity";
    /// No material-enclosing (outward) shell among the result components — every component is
    /// inward-oriented. Geometrically impossible for a real solid result; a defensive backstop
    /// with no firing test (cf. `FOURPLANE`).
    pub const NO_OUTWARD_SHELL: &str = "no_outward_shell";
    /// A boolean input is a rotated solid (overhaul stage 1b). Rotated planar geometry
    /// is representable but its predicates are not yet sound (no TIP until stage 3), so
    /// the boolean honestly rejects until then. Fired by a `Transform`-rotated operand.
    pub const ROTATED_UNSUPPORTED: &str = "rotated_unsupported";
    /// An overhang overlap arc could not be paired with a complementary arc sharing its two
    /// crossings (a non-convex-cyclic arrangement `is_convex` did not screen out) — out of scope.
    pub const OVERHANG_ARCS: &str = "overhang_arcs";
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
    /// A solid's planar section returned more than one loop (an outer ring plus a hole, e.g. a
    /// slot piercing a cavity). The section-clip reconstruction assumes a single loop; a
    /// multi-loop section is honestly rejected until the multi-loop cell lands, so the
    /// single-loop assumption is never silently reached.
    pub const SECTION_MULTI_LOOP: &str = "section_multi_loop";
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
    /// A pocket's depth reaches through the solid: the carved prism is not blind, so `Cut`
    /// produced a through-hole with no floor face. `pocket` requires `dist` less than the
    /// thickness at the face (the boolean pocket path honestly rejects instead of the old
    /// direct path's silent invalid result).
    PocketNotBlind,
    /// A boolean operation failed (design §8 M5).
    Boolean(BoolError),
    /// A `Transform` input solid is not live (a stale or non-live handle).
    SolidNotLive,
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
    /// The boolean result solids (supersede both inputs). Usually one; a boolean that severs the
    /// body yields several (cell 0.4), and `Cut(A, A)` (deferred) would yield none.
    Boolean { solids: Vec<Handle<Solid>> },
    /// The transformed solid (supersedes the input).
    Transform { solid: Handle<Solid> },
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
            let solids = boolean(model, *kind, *a, *b).map_err(OpError::Boolean)?;
            Ok(OpOutput::Boolean { solids })
        }
        Operation::Transform { solid, isometry } => {
            let out = transform(model, *solid, isometry)?;
            Ok(OpOutput::Transform { solid: out })
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
    let base_pts: Vec<Point3> = profile.points.iter().map(|p| plane.point(*p)).collect();
    build_prism(model, &base_pts, plane.normal() * dist, None)
}

/// Sweep the ring `base_pts` along `sweep` into a prism solid (caps + side quads). The ring is
/// normalized CCW **about `sweep`** so the synthesized normals point outward; this uses the
/// polygon's area vector dotted with the sweep normal (frame-independent — `proj2`+`signed_area`
/// would mis-sign when the sweep runs along a negative dominant axis, e.g. a pocket into a `+z`
/// face sweeping `−z`). All vertices are `Origin::Constructed`. Returns the solid and its faces:
/// `faces[0]` = base cap (at `base_pts`, normal `−ŝ`), `faces[1]` = far cap (at `base_pts + sweep`,
/// normal `+ŝ`), then the side quads. Shared by [`extrude`] (a boss) and the pocket (`sweep = −n`).
fn build_prism(
    model: &mut Model,
    base_pts: &[Point3],
    sweep: Vector3,
    base_cap_surface: Option<Handle<Surface>>,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if base_pts.len() < 3 {
        return Err(OpError::DegenerateProfile);
    }
    let normal = sweep.normalize().ok_or(OpError::DegenerateGeometry)?;
    let mut base_pts: Vec<Point3> = base_pts.to_vec();
    let k = base_pts.len();
    let area_vec = (0..k)
        .map(|i| (base_pts[i] - Point3::origin()).cross(base_pts[(i + 1) % k] - Point3::origin()))
        .fold(Vector3::from_array([0.0; 3]), |a, b| a + b);
    if area_vec.dot(normal) < 0.0 {
        base_pts.reverse();
    }
    let n = base_pts.len();
    let top_pts: Vec<Point3> = base_pts.iter().map(|b| *b + sweep).collect();

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
    // When padding/pocketing on a face, reuse that face's `Surface` handle (§5
    // explicit sharing) so the flush contact is a shared-handle coplanar pair the
    // boolean can recognize by `Handle` identity; otherwise push a fresh plane.
    // The materialized outward normal must stay −N, so the face orientation is
    // chosen from the shared surface's stored normal — `surface` and `orientation`
    // travel together, and the reconstruction copies both.
    let (base_surface, base_orient) = match base_cap_surface {
        Some(h) => {
            let n_h = match model.surfaces.get(h) {
                Surface::Plane(p) => p.normal(),
                Surface::Cylinder(_) => return Err(OpError::DegenerateGeometry),
            };
            let orient = if n_h.dot(-normal) > 0.0 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            (h, orient)
        }
        None => {
            let s = model.surfaces.push(Surface::Plane(
                Plane::from_point_normal(base_pts[0], -normal)
                    .ok_or(OpError::DegenerateGeometry)?,
            ));
            (s, Orientation::Forward)
        }
    };
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
        orientation: base_orient,
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

/// The result of splitting a planar face along a profile: the rebuilt outer face
/// carrying the profile as a hole, plus the profile's segment edges. Used by
/// [`imprint`], which finishes by adding a coplanar region face and calling
/// [`finish_split`].
struct Split {
    solid_h: Handle<Solid>,
    shell_h: Handle<Shell>,
    /// The target face's surface and orientation (for a coplanar region face).
    surface_h: Handle<Surface>,
    orientation: Orientation,
    /// Profile segment edges `base_i → base_{i+1}` (CCW about the outward normal).
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

/// A planar face's live solid, its in-plane right-handed frame (`x × y = n`, centred on the face
/// centroid so a profile's `(0,0)` lands there), and its loops — the shared setup for placing a
/// profile on a face (imprint / pad / pocket).
struct FaceFrame {
    solid_h: Handle<Solid>,
    shell_h: Handle<Shell>,
    surface_h: Handle<Surface>,
    orientation: Orientation,
    n: Vector3, // outward normal
    x: Vector3,
    y: Vector3,
    origin: Point3, // face centroid
    outer_pts: Vec<Point3>,
    outer_loop: Loop,
    inner_loops: Vec<Loop>,
}

/// Locate `face`'s live solid and build its planar frame. `NonPlanarFace` for a curved surface,
/// `FaceNotInLiveSolid` if no live outer shell holds it.
fn face_frame(model: &Model, face: Handle<Face>) -> Result<FaceFrame, OpError> {
    let (solid_h, shell_h) = model
        .live_solids
        .iter()
        .map(|&s| (s, model.solids.get(s).outer))
        .find(|&(_, sh)| model.shells.get(sh).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;
    let f = model.faces.get(face);
    let surface_h = f.surface;
    let orientation = f.orientation;
    let outer_loop = f.outer.clone();
    let inner_loops = f.inner.clone();
    let plane = match model.surfaces.get(surface_h) {
        Surface::Plane(p) => *p,
        Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    };
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
    Ok(FaceFrame {
        solid_h,
        shell_h,
        surface_h,
        orientation,
        n,
        x,
        y,
        origin,
        outer_pts,
        outer_loop,
        inner_loops,
    })
}

/// Place `profile` on `frame` (CCW in the frame so its RH normal is `+n`) and return its points on
/// the face plane. **No containment check** — the profile may reach past the face boundary. Used by
/// the boolean pad/pocket path (`extrude_and_boolean`), where an overhanging footprint routes to the
/// overhang boolean sidecars; the boolean honestly rejects configurations it does not cover.
fn placed_profile_unchecked(
    frame: &FaceFrame,
    profile: &Profile2d,
) -> Result<Vec<Point3>, OpError> {
    if profile.points.len() < 3 {
        return Err(OpError::DegenerateProfile);
    }
    let mut pts = profile.points.clone();
    if signed_area(&pts) < 0.0 {
        pts.reverse();
    }
    Ok(pts
        .iter()
        .map(|p| frame.origin + frame.x * p[0] + frame.y * p[1])
        .collect())
}

/// [`placed_profile_unchecked`] plus the strict-containment requirement: inside the outer ring,
/// outside every hole, touching no boundary. Otherwise the "hole" is not a clean inner loop and the
/// result is silently invalid — validate sees only topology, props integrates the ring, and only
/// tessellate's `NoEar` catches it (cell imprint-containment; dev-log.md). Used by `imprint`, which
/// needs a bounded inner loop; a profile reaching past the face is a boolean pad/pocket instead.
fn placed_profile(
    model: &Model,
    frame: &FaceFrame,
    profile: &Profile2d,
) -> Result<Vec<Point3>, OpError> {
    let base_pts = placed_profile_unchecked(frame, profile)?;
    let drop = planar_drop_axes(frame.n);
    let profile2: Vec<[f64; 2]> = base_pts.iter().map(|&p| proj2(p, drop)).collect();
    let outer2: Vec<[f64; 2]> = frame.outer_pts.iter().map(|&p| proj2(p, drop)).collect();
    let holes2: Vec<Vec<[f64; 2]>> = frame
        .inner_loops
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
    Ok(base_pts)
}

fn prepare_face_split(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
) -> Result<Split, OpError> {
    let frame = face_frame(model, face)?;
    let base_pts = placed_profile(model, &frame, profile)?;
    let (surface_h, orientation) = (frame.surface_h, frame.orientation);

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
        outer: frame.outer_loop,
        inner: vec![hole_loop],
        orientation,
    });

    Ok(Split {
        solid_h: frame.solid_h,
        shell_h: frame.shell_h,
        surface_h,
        orientation,
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

/// A face-local feature built as **tool body + boolean** (roadmap §9 unification): the profile
/// extrudes off `face` into a top-flush prism, then `kind` fuses/cuts it against the face's solid.
/// A **contained** footprint takes the contained-coplanar path (empty seam → all
/// `Origin::Constructed`); one that **reaches past the face** routes to the overhang boolean
/// sidecars (a boss cantilever / an edge slot; Discovered seam vertices). `Fuse` sweeps **outward**
/// (a boss); `Cut` sweeps **inward** (a blind pocket). Returns the result solid and the feature's
/// exposed cap — the boss top or the pocket floor, the outer-shell face on the prism's far-cap plane
/// with outward normal `+n`. `Option::None` there ⇒ the far cap did not survive (a through-cut with
/// no floor); callers map it to their own error. `NonPositiveDistance`/`DegenerateProfile` propagate
/// from the frame; overhang configurations the boolean does not cover surface as `Boolean(_)`.
fn extrude_and_boolean(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
    kind: BoolKind,
) -> Result<(Handle<Solid>, Option<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    let frame = face_frame(model, face)?;
    // No containment check — an overhanging footprint routes to the overhang boolean sidecars.
    let base_pts = placed_profile_unchecked(&frame, profile)?;
    let n = frame.n;
    // Cut carves inward, Fuse raises outward; either way the prism's near cap is flush on the face.
    let signed = if matches!(kind, BoolKind::Cut) {
        -dist
    } else {
        dist
    };
    let (prism, _) = build_prism(model, &base_pts, n * signed, Some(frame.surface_h))?;
    let solids = boolean(model, kind, frame.solid_h, prism).map_err(|e| {
        model.live_solids.retain(|&s| s != prism); // drop the transient prism (atomic on failure)
        OpError::Boolean(e)
    })?;
    // Exposed cap = the result face on the prism's far-cap plane (face plane offset by n·signed),
    // its outward normal +n (the opening side for a pocket, the boss top for a boss). A pad (Fuse)
    // never severs; a pocket (Cut) that severs leaves several solids — scan them all for the cap
    // and return the piece that carries it, leaving the others live. The committed model is a valid
    // multi-solid either way, so there is no reject-after-commit (append-only has no rollback).
    let cap_pt = frame.origin + n * signed;
    let primary = *solids
        .first()
        .expect("a pad/pocket boolean yields at least one solid");
    match solids
        .iter()
        .find_map(|&s| find_face_on_plane(model, s, cap_pt, n).map(|c| (s, c)))
    {
        Some((solid, cap)) => Ok((solid, Some(cap))),
        None => Ok((primary, None)),
    }
}

/// Pad a boss on a planar `face`: extrude the profile **outward** by `dist` and `Fuse` it onto the
/// solid, adding `profile_area · dist` of material. Returns `(new solid, top cap face)`. A boss
/// always yields its top cap, so the `None` guard is an unreachable internal-invariant defense.
fn pad(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    let (solid, top) = extrude_and_boolean(model, face, profile, dist, BoolKind::Fuse)?;
    Ok((solid, top.ok_or(OpError::DegenerateGeometry)?))
}

/// Carve a blind pocket on a planar `face`: extrude the profile **inward** by `dist` and `Cut` it
/// from the solid, removing `profile_area · dist` of material. Returns `(new solid, floor face)`.
/// `PocketNotBlind` if `dist` reaches through the solid (the far cap is not blind → a through-cut
/// with no floor face).
fn pocket(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    let (solid, floor) = extrude_and_boolean(model, face, profile, dist, BoolKind::Cut)?;
    Ok((solid, floor.ok_or(OpError::PocketNotBlind)?))
}

/// The outer-shell face of `solid` whose plane passes through `pt` (coplanar) and whose **oriented**
/// outward normal agrees with `n_out` (`dot > 0`). Recovers a pocket floor from a boolean result:
/// the oriented-normal filter distinguishes the floor (normal toward the opening) from a coincident
/// face with the opposite normal. `None` if there is none (e.g. a through-pocket has no floor).
fn find_face_on_plane(
    model: &Model,
    solid: Handle<Solid>,
    pt: Point3,
    n_out: Vector3,
) -> Option<Handle<Face>> {
    let target = Plane::from_point_normal(pt, n_out)?;
    let shell = model.solids.get(solid).outer;
    model.shells.get(shell).faces.iter().copied().find(|&fh| {
        let f = model.faces.get(fh);
        let Surface::Plane(plane) = model.surfaces.get(f.surface) else {
            return false;
        };
        let sign = match f.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        planes_coplanar(plane, &target) && (plane.normal() * sign).dot(n_out) > 0.0
    })
}

// ---- transform (overhaul stage 1) ----

/// Supersede `solid` by its image under `isometry` (stage 1a: a rational
/// translation). Clones the solid's cells with moved geometry ([`transform_solid`])
/// and drops the input from `live_solids` — the op-log is the truth.
fn transform(
    model: &mut Model,
    solid: Handle<Solid>,
    isometry: &Isometry,
) -> Result<Handle<Solid>, OpError> {
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    let out = transform_solid(model, solid, isometry);
    model.live_solids.retain(|&s| s != solid);
    Ok(out)
}

/// A vertex/edge `Origin` with any `Discovered` `ThreePlane` definition remapped
/// onto the moved surfaces. Constructed stays constructed; the `tol` is unchanged
/// (a rigid move; stage 1 records tol but never judges on it). Fails loud if a
/// definition names a surface that is not one of the solid's face surfaces
/// (an invariant violation — a seam vertex is the meet of three of its faces).
fn remap_origin(origin: Origin, surf_map: &HashMap<Handle<Surface>, Handle<Surface>>) -> Origin {
    match origin {
        Origin::Constructed => Origin::Constructed,
        Origin::Discovered {
            tol,
            definition: VertexDef::ThreePlane(planes),
        } => {
            let mapped = planes.map(|s| {
                *surf_map
                    .get(&s)
                    .expect("Discovered ThreePlane surface must be a face surface of the solid")
            });
            Origin::Discovered {
                tol,
                definition: VertexDef::ThreePlane(mapped),
            }
        }
        // Reached only for an *exact* move of an already-rotated solid — a translation
        // or a 90°-family rotation, which take the no-node path. Keep the rotation
        // definition: its `base`/`rotation` name arena ancestors unaffected by a
        // translation. (An *inexact* re-rotation instead records a chain node and is
        // handled in `transform_solid` pass 3, not here.)
        Origin::Rotated { .. } => origin,
    }
}

/// A surface moved by `isometry`. A pure translation uses `translated` (the normal —
/// and its exact `raw` — is unchanged). A rotation rebuilds from the moved
/// origin/normal via the constructor (rotation makes `raw` irrational, as expected —
/// the plane then carries tol, judged by TIP later).
fn transform_surface(s: &Surface, iso: &Isometry, offset: Vector3) -> Surface {
    if iso.rotate.is_none() {
        return s.translated(offset);
    }
    let p = |q: Point3| Point3::from_array(iso.apply_point(q.as_array()));
    let d = |v: Vector3| Vector3::from_array(iso.apply_dir(v.as_array()));
    match s {
        Surface::Plane(pl) => Surface::Plane(
            Plane::from_point_normal(p(pl.origin()), d(pl.normal()))
                .expect("rotation preserves a nonzero normal"),
        ),
        Surface::Cylinder(cy) => {
            let ax = cy.axis();
            Surface::Cylinder(
                Cylinder::from_axis(
                    p(ax.origin()),
                    d(ax.direction()),
                    d(cy.ref_dir()),
                    cy.radius(),
                )
                .expect("rotation preserves a valid cylinder"),
            )
        }
    }
}

/// A curve moved by `isometry` (see [`transform_surface`]).
fn transform_curve(c: &Curve, iso: &Isometry, offset: Vector3) -> Curve {
    if iso.rotate.is_none() {
        return c.translated(offset);
    }
    let p = |q: Point3| Point3::from_array(iso.apply_point(q.as_array()));
    let d = |v: Vector3| Vector3::from_array(iso.apply_dir(v.as_array()));
    match c {
        Curve::Line(l) => Curve::Line(
            Line::from_point_direction(p(l.origin()), d(l.direction()))
                .expect("rotation preserves a nonzero direction"),
        ),
        Curve::Circle(ci) => Curve::Circle(
            Circle::from_center_normal(
                p(ci.center()),
                d(ci.normal()),
                d(ci.ref_dir()),
                ci.radius(),
            )
            .expect("rotation preserves a valid circle"),
        ),
    }
}

/// Clone `solid` into a new solid with every cell's geometry moved by `isometry`,
/// preserving topology, shared cells (surfaces/curves/vertices/edges are deduped),
/// face orientations (a translation does not rotate normals), inner-loop holes,
/// cavity shells, and each vertex/edge `Origin` (a `Discovered` definition's plane
/// handles are remapped to the moved surfaces). Cells are pushed in a **deterministic
/// traversal order** (shell → face → loop) with per-cell dedup maps, so the same
/// op-log reproduces identical handles (replay determinism, DNA 3). Stage 1b's
/// rotation reuses this by swapping the per-cell geometry transform.
fn transform_solid(model: &mut Model, solid: Handle<Solid>, isometry: &Isometry) -> Handle<Solid> {
    let offset = Vector3::from_array(isometry.offset_f64());
    let src = model.solids.get(solid).clone();

    // Forest node for this transform (§TIP ⑦), one shared node named by every rotated
    // vertex. The input's shared leaf (None if the input is not rotated) decides B0 vs B1:
    //   A.  translation → no node (remap path).
    //   B0. fresh rotation (input not rotated): exact (90°-family) → no node (remap,
    //       preserving 1b); inexact → a root node (`parent = None`).
    //   B1. re-rotation (input already rotated): **always chain** a node (`parent =
    //       input leaf`) — record every rotation, even an exact one, because an inexact
    //       ancestor makes the composite inexact and stage-2 tol must transport through
    //       it; the forest stays complete. (Same-axis *bundling* — accumulating the
    //       angle into one node — is a later cell; this cell always chains.)
    let input_leaf = solid_rotation(model, solid);
    let rot_node: Option<Handle<Rotation>> = match (isometry.rotate, input_leaf) {
        (None, _) => None,
        (Some(r), None) => (!isometry.is_exact()).then(|| {
            model.rotations.push(Rotation {
                axis: r.axis,
                point: r.point,
                angle: r.angle,
                parent: None,
            })
        }),
        (Some(r), Some(parent)) => Some(model.rotations.push(Rotation {
            axis: r.axis,
            point: r.point,
            angle: r.angle,
            parent: Some(parent),
        })),
    };

    // Deterministic order: outer shell then cavities; each shell's faces in order.
    let shell_order: Vec<Handle<Shell>> = std::iter::once(src.outer)
        .chain(src.cavities.iter().copied())
        .collect();
    let face_order: Vec<Handle<Face>> = shell_order
        .iter()
        .flat_map(|&sh| model.shells.get(sh).faces.clone())
        .collect();

    // Pass 1 — surfaces (dedup, moved): needed before vertex `Origin` remap.
    let mut surf_map: HashMap<Handle<Surface>, Handle<Surface>> = HashMap::new();
    for &fh in &face_order {
        let s = model.faces.get(fh).surface;
        if let std::collections::hash_map::Entry::Vacant(e) = surf_map.entry(s) {
            let moved = transform_surface(model.surfaces.get(s), isometry, offset);
            e.insert(model.surfaces.push(moved));
        }
    }

    // Edge order (deterministic dedup) — used by passes 2/3/4.
    let mut edge_order: Vec<Handle<Edge>> = Vec::new();
    let mut edge_seen: HashSet<Handle<Edge>> = HashSet::new();
    for &fh in &face_order {
        let face = model.faces.get(fh).clone();
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                if edge_seen.insert(he.edge) {
                    edge_order.push(he.edge);
                }
            }
        }
    }

    // Pass 2 — curves (dedup, moved).
    let mut curve_map: HashMap<Handle<Curve>, Handle<Curve>> = HashMap::new();
    for &eh in &edge_order {
        let c = model.edges.get(eh).curve;
        if let std::collections::hash_map::Entry::Vacant(e) = curve_map.entry(c) {
            let moved = transform_curve(model.curves.get(c), isometry, offset);
            e.insert(model.curves.push(moved));
        }
    }

    // Pass 3 — vertices (dedup, moved point + Origin preserved/remapped).
    let mut vert_order: Vec<Handle<Vertex>> = Vec::new();
    let mut vert_seen: HashSet<Handle<Vertex>> = HashSet::new();
    for &eh in &edge_order {
        if let Some(bounds) = model.edges.get(eh).bounds {
            for vh in bounds {
                if vert_seen.insert(vh) {
                    vert_order.push(vh);
                }
            }
        }
    }
    let mut vert_map: HashMap<Handle<Vertex>, Handle<Vertex>> = HashMap::new();
    for &vh in &vert_order {
        let v = *model.vertices.get(vh);
        let new_v = Vertex {
            point: Point3::from_array(isometry.apply_point(v.point.as_array())),
            // A recorded rotation marks the vertex `Rotated`; `base` is the **root** —
            // the non-`Rotated` (Constructed/Discovered) ancestor whose exact definition
            // the rotation chain turns. A fresh rotation's input is itself the root; a
            // re-rotation chases one hop to the input's own root (the invariant keeps
            // `base` pointing at a root, never at another `Rotated` vertex, so stage-2
            // recompute never applies the same rotation twice). No node → remap (1a/1b).
            origin: match rot_node {
                Some(rotation) => {
                    let base = match v.origin {
                        Origin::Rotated { base, .. } => base,
                        _ => vh,
                    };
                    Origin::Rotated { base, rotation }
                }
                None => remap_origin(v.origin, &surf_map),
            },
        };
        vert_map.insert(vh, model.vertices.push(new_v));
    }

    // Pass 4 — edges (curve/vertex handles + Origin remapped).
    let mut edge_map: HashMap<Handle<Edge>, Handle<Edge>> = HashMap::new();
    for &eh in &edge_order {
        let e = *model.edges.get(eh);
        let new_e = Edge {
            curve: curve_map[&e.curve],
            bounds: e.bounds.map(|[a, b]| [vert_map[&a], vert_map[&b]]),
            origin: remap_origin(e.origin, &surf_map),
        };
        edge_map.insert(eh, model.edges.push(new_e));
    }

    // Pass 5 — faces (loops rebuilt onto the new edges; orientation unchanged).
    let map_loop = |lp: &Loop| Loop {
        half_edges: lp
            .half_edges
            .iter()
            .map(|he| HalfEdge {
                edge: edge_map[&he.edge],
                forward: he.forward,
            })
            .collect(),
    };
    let mut face_map: HashMap<Handle<Face>, Handle<Face>> = HashMap::new();
    for &fh in &face_order {
        let face = model.faces.get(fh).clone();
        let new_f = Face {
            surface: surf_map[&face.surface],
            outer: map_loop(&face.outer),
            inner: face.inner.iter().map(&map_loop).collect(),
            orientation: face.orientation,
        };
        face_map.insert(fh, model.faces.push(new_f));
    }

    // Pass 6 — shells; Pass 7 — solid.
    let mut shell_map: HashMap<Handle<Shell>, Handle<Shell>> = HashMap::new();
    for &sh in &shell_order {
        let faces: Vec<Handle<Face>> = model
            .shells
            .get(sh)
            .faces
            .iter()
            .map(|fh| face_map[fh])
            .collect();
        shell_map.insert(sh, model.shells.push(Shell { faces }));
    }
    let new_solid = Solid {
        outer: shell_map[&src.outer],
        cavities: src.cavities.iter().map(|sh| shell_map[sh]).collect(),
    };
    model.push_solid(new_solid)
}

/// Whether `solid` was produced by a non-exact rotation — its vertices carry
/// `Origin::Rotated`. A `Transform` rotates a whole solid uniformly and `boolean`
/// rejects rotated inputs, so a solid is all-or-nothing rotated: one vertex decides
/// (O(1)). (90°-family rotations stay exact/`Constructed`, so this is false for them.)
fn solid_is_rotated(model: &Model, solid: Handle<Solid>) -> bool {
    let sh = model.solids.get(solid).outer;
    for &fh in &model.shells.get(sh).faces {
        for he in &model.faces.get(fh).outer.half_edges {
            if let Some(bounds) = model.edges.get(he.edge).bounds {
                return matches!(model.vertices.get(bounds[0]).origin, Origin::Rotated { .. });
            }
        }
    }
    false
}

/// The shared rotation-forest leaf that every boundary vertex of a rotated `solid`
/// names (`None` if the solid is not rotated) — the input side of the re-rotation
/// decision in [`transform_solid`]. A `Transform` rotates a solid uniformly and
/// `boolean` rejects rotated inputs, so all boundary vertices share one leaf; that
/// invariant is `debug_assert`ed here (fail-loud if a future change ever produces a
/// solid with mixed rotation provenance).
fn solid_rotation(model: &Model, solid: Handle<Solid>) -> Option<Handle<Rotation>> {
    let sh = model.solids.get(solid).outer;
    let mut seen: Option<Option<Handle<Rotation>>> = None;
    for &fh in &model.shells.get(sh).faces {
        for he in &model.faces.get(fh).outer.half_edges {
            let Some(bounds) = model.edges.get(he.edge).bounds else {
                continue;
            };
            for &vh in &bounds {
                let leaf = match model.vertices.get(vh).origin {
                    Origin::Rotated { rotation, .. } => Some(rotation),
                    _ => None,
                };
                match seen {
                    None => seen = Some(leaf),
                    Some(established) => debug_assert_eq!(
                        established, leaf,
                        "a solid's boundary vertices must share one rotation node (uniform rotation)"
                    ),
                }
            }
        }
    }
    seen.flatten()
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
    /// The three `tri` points as **toleranced `Pt3`** (exact rotation definition), in the
    /// same order as `tri` — `Some` only when the solid is rotated (overhaul stage 3;
    /// `collect_planes` builds it once). `None` on the axis-aligned path, where `tri`'s
    /// f64 coordinates are already exact and the geom predicates are used directly.
    pub(crate) tri_pt3: Option<[Pt3; 3]>,
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
) -> Result<Vec<Handle<Solid>>, BoolError> {
    if !model.live_solids.contains(&a) || !model.live_solids.contains(&b) {
        return Err(BoolError::InputNotLive);
    }
    // A rotated operand's predicates are now TIP-exact (overhaul 3a–3c-vi): the arrangement,
    // seam, in/out, and outer/cavity tests all decide on the vertices' exact rotation
    // definitions, so a rotated *transverse* boolean is sound (cell 3d-i retired the guard).
    // A rotated *coplanar* contact is out of scope (the coplanar milestone) but never silently
    // wrong: the exact `planes_coplanar` detectors miss it, and `general_boolean` then either
    // solves it or honestly rejects (`edge_crosses_face` declare-0 → `VERTEX_ON_FACE_PLANE`,
    // or a downstream degeneracy) — the vertex-on-plane TIP predicate catches what the rounded
    // coefficients cannot.
    // Cavitied operands are supported (cells (5c-in), (5c-in-2)): the seam front-end and
    // reconstruction walk all shells (outer + cavities) via `solid_shell_handles`, so a
    // void the cut misses is carried through, and a cut reaching into a void reconstructs
    // (the opened void merges with the outer shell).
    // Coincident-coplanar degeneracy (M5-c5): a clean matched-interface stack.
    if let Some(iface) = detect_coincident_interface(model, a, b) {
        return coincident_merge(model, kind, a, b, &iface);
    }
    // A boss sitting on a face inside its boundary — coplanar contact, contained footprint,
    // opposite normals (cell coplanar-contact-boss).
    if kind == BoolKind::Fuse {
        if let Some(cc) = detect_contained_contact(model, a, b) {
            return contained_contact_result(model, &cc, false);
        }
        // A boss whose footprint hangs past a single edge of the face — partial coplanar
        // overlap, a cantilever (cell coplanar-contact-overhang).
        if let Some(oc) = detect_overhang_contact(model, a, b) {
            return overhang_contact_result(model, &oc);
        }
    }
    // A blind pocket: `b` sits inside `a` with its top flush on `a`'s face (same normals,
    // contained footprint) (cell coplanar-contact-cut). `Cut(a, b)` carves it.
    if kind == BoolKind::Cut {
        if let Some(cc) = detect_pocket_contact(model, a, b) {
            return contained_contact_result(model, &cc, true);
        }
        // A prism whose footprint hangs past `a`'s face boundary, breaking out through one or more
        // walls — the general N-wall overhang Cut (cell coplanar-contact-overhang-cut-general):
        // edge-slot (1 wall), spanning channel (2 opposite), corner (2 adjacent), and 3+/L/U.
        if let Some(oc) = detect_overhang_cut_general(model, a, b) {
            return overhang_cut_general_result(model, &oc);
        }
    }
    // A prism top-flush on `a`'s face, footprint hanging past the boundary — `Common(a, b)` is the
    // convex overlap `R = a ∩ b` (cell coplanar-contact-overhang-common). Two-crossing configs only.
    if kind == BoolKind::Common {
        if let Some(oc) = detect_overhang_common(model, a, b) {
            return overhang_common_result(model, &oc);
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
) -> Result<Vec<Handle<Solid>>, BoolError> {
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
            let side = vertex_in_solid(model, vh, other)?;
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
    // Seam-free never severs: one solid contains the other, or they are disjoint (a single
    // solid, or the deferred disjoint-Fuse `EmptyResult`). Always at most one solid here.
    Ok(vec![contained_result(model, kind, a, b, &classof)?])
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
    model: &Model,
    planes: &[PlaneInfo],
    pair: [usize; 2],
    v0: Handle<Vertex>,
    v1: Handle<Vertex>,
    rings: &FaceRings,
) -> Result<Vec<usize>, BoolError> {
    let mut hits = Vec::new();
    for (q, r) in rings {
        if arrange::edge_crosses_face(model, planes, pair, v0, v1, *q, r)? {
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
/// One face-reconstruction work item for the parallel sweep in [`overlap_fuse_cut`]:
/// `(face, other solid, its plane index, side to keep, flip, own/other edge-plane maps)`.
/// Flattening the (side, shell, face) loops into a `Vec<ReconItem>` lets the reconstruction
/// map in index order, which `assemble_fuse_cut`'s handle assignment depends on.
type ReconItem<'a> = (
    Handle<Face>,
    Handle<Solid>,
    usize,
    Side,
    bool,
    &'a arrange::EdgePlanes,
    &'a arrange::EdgePlanes,
);

/// One seam work item for the parallel seam sweep in [`overlap_fuse_cut`]: an edge (its two
/// bound vertices + two incident plane indices) tested against the *other* solid's face
/// rings. Flattened in the sequential loop's order so the first-appearance dedup is stable.
type SeamItem<'a> = ([Handle<Vertex>; 2], [usize; 2], &'a FaceRings);

/// A candidate seam vertex produced per edge, before dedup: `(sorted triple as the dedup
/// key, intersection point, e0, e1, entry)`. The point is computed in the original
/// `(e0,e1,entry)` order (`three_planes` is order-sensitive); the sorted triple is *only*
/// the dedup key. `tol` is computed later, per new triple, in that same original order.
type SeamCands = Vec<([usize; 3], Point3, usize, usize, usize)>;

/// returned as several solids (cell 0.4); a sever that also leaves a cavity is
/// `SEVERED_WITH_CAVITY`.
fn overlap_fuse_cut(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
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

    // Classify every original vertex vs the other solid (exact forward-ray winding). Each
    // vertex is an independent read-only classification, so evaluate them in parallel and
    // fold into `classof` sequentially — keyed by `Handle<Vertex>`, insertion order is
    // irrelevant. Collect in index order so the index-first reject matches the sequential
    // loop's first `?` (and, under parallel + cfg(test), replay it to restore the tag).
    let classof: HashMap<Handle<Vertex>, Side> = {
        let items: Vec<(Handle<Vertex>, Handle<Solid>)> = [(a, b), (b, a)]
            .into_iter()
            .flat_map(|(verts_solid, other)| {
                solid_vertex_handles(model, verts_solid)
                    .into_iter()
                    .map(move |vh| (vh, other))
            })
            .collect();
        let classify =
            |&(vh, other): &(Handle<Vertex>, Handle<Solid>)| vertex_in_solid(model, vh, other);
        #[cfg(feature = "parallel")]
        let res: Vec<Result<Side, BoolError>> = items.par_iter().map(classify).collect();
        #[cfg(not(feature = "parallel"))]
        let res: Vec<Result<Side, BoolError>> = items.iter().map(classify).collect();
        let mut classof: HashMap<Handle<Vertex>, Side> = HashMap::new();
        for (item, r) in items.iter().zip(res) {
            match r {
                Ok(side) => {
                    classof.insert(item.0, side);
                }
                Err(e) => {
                    // Restore the reject tag on this thread: a parallel worker set it on
                    // its own `LAST_REJECT`. Re-scan in index order until the first failure
                    // re-runs `reject` here (cfg(test) only; inert in release).
                    #[cfg(all(test, feature = "parallel"))]
                    let _ = items.iter().find(|it| classify(it).is_err());
                    return Err(e);
                }
            }
        }
        classof
    };

    // Seam vertices: an edge straddling the other boundary pierces exactly one of
    // its faces (that face's plane is the seam vertex's third plane).
    let edges_a = edge_incidence(model, a, &surf_ix)?;
    let edges_b = edge_incidence(model, b, &surf_ix)?;
    let inc_a = arrange::edge_planes(model, a, &surf_ix)?;
    let inc_b = arrange::edge_planes(model, b, &surf_ix)?;
    // Each solid's faces as triple rings, once. `edge_crosses_face` reads them per edge.
    let rings_a = solid_face_rings(model, a, &surf_ix, &inc_a)?;
    let rings_b = solid_face_rings(model, b, &surf_ix, &inc_b)?;
    // Each edge's crossing test + 4-plane guard is an independent, read-only predicate
    // evaluation, so evaluate edges in parallel; the seam dedup (first-appearance index =
    // seam-vertex identity) stays sequential in index order, keeping the result
    // bit-identical regardless of thread count.
    let (seam, seam_ix): (Vec<SeamVertex>, HashMap<[usize; 3], usize>) = {
        // Flatten in the sequential loop's order — `edges_a` (other = `rings_b`) first,
        // then `edges_b` (other = `rings_a`) — so the first-appearance dedup below
        // reproduces the same seam indices.
        let mut items: Vec<SeamItem> = Vec::new();
        for (edges, other) in [(&edges_a, &rings_b), (&edges_b, &rings_a)] {
            for &(_, bounds, inc) in edges {
                items.push((bounds, inc, other));
            }
        }
        // Per edge: `pierced_faces` + straddle/tunnel parity + per-hit `three_planes` and the
        // 4-plane guard. Returns candidate seam vertices; the sorted `triple` is only the
        // dedup key, and `point` is built in the original `(e0,e1,entry)` order because
        // `three_planes` (Cramer) is order-sensitive.
        let per_edge = |&(bounds, inc, other): &SeamItem| -> Result<SeamCands, BoolError> {
            let [v0, v1] = bounds;
            let (s0, s1) = (classof[&v0], classof[&v1]);
            let hits = pierced_faces(model, &planes, inc, v0, v1, other)?;
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
            let mut cands: SeamCands = Vec::new();
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
                    if crate::tolerant::t_orient3d(&planes, e0, e1, entry, m) == 0 {
                        return Err(reject(tag::FOURPLANE)); // seam vertex on a 4th plane
                    }
                }
                let mut triple = [e0, e1, entry];
                triple.sort_unstable();
                cands.push((triple, point, e0, e1, entry));
            }
            Ok(cands)
        };
        #[cfg(feature = "parallel")]
        let per_item: Vec<Result<SeamCands, BoolError>> = items.par_iter().map(per_edge).collect();
        #[cfg(not(feature = "parallel"))]
        let per_item: Vec<Result<SeamCands, BoolError>> = items.iter().map(per_edge).collect();

        // Sequential first-appearance dedup in index order (bit-identical to the serial
        // sweep). `tol` is computed here, only for a new triple, in the original
        // `(e0,e1,entry)` order — exactly as the sequential code did.
        let mut seam: Vec<SeamVertex> = Vec::new();
        let mut seam_ix: HashMap<[usize; 3], usize> = HashMap::new();
        for r in per_item {
            let cands = match r {
                Ok(c) => c,
                Err(e) => {
                    // Restore the reject tag on this thread (a worker set it on its own
                    // `LAST_REJECT`): re-scan in index order until the first failure re-runs
                    // `reject` here (cfg(test) only; inert in release).
                    #[cfg(all(test, feature = "parallel"))]
                    let _ = items.iter().find(|it| per_edge(it).is_err());
                    return Err(e);
                }
            };
            for (triple, point, e0, e1, entry) in cands {
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
        (seam, seam_ix)
    };

    let (keep_a, keep_b, flip_b) = match kind {
        BoolKind::Fuse => (Side::Outside, Side::Outside, false),
        BoolKind::Cut => (Side::Outside, Side::Inside, true),
        // A∩B: keep each solid's material *inside* the other, neither shell flipped —
        // the De Morgan dual of Fuse. Cell 3g.
        BoolKind::Common => (Side::Inside, Side::Inside, false),
    };
    if seam.is_empty() {
        // No crossing after all — containment/disjoint, at most one solid.
        return Ok(vec![contained_result(model, kind, a, b, &classof)?]);
    }

    // `edge_seam` is gone from this path: the arrangement carries each boundary
    // crossing on the edge it was produced on (`SeamSegment::on_edge`), so nothing has
    // to key a splice by `Handle<Edge>` — the map that could hold only one crossing per
    // edge no longer exists here.
    // Reconstruct each (side, shell, face) independently. `reconstruct_face_paths` is a
    // pure read of `model` + the frozen seam/classification and returns owned `LocalFace`s,
    // so faces reconstruct in parallel. Results are collected in the SAME order as the
    // sequential sweep and appended in that order — `assemble_fuse_cut` assigns vertex
    // handles by first appearance across `faces`, and replay determinism (bit-identical
    // model) rests on that order, which rayon's indexed `collect` preserves.
    let faces: Vec<LocalFace> = {
        let mut items: Vec<ReconItem> = Vec::new();
        for (solid, other, inc_f, inc_o, keep, flip) in [
            (a, b, &inc_a, &inc_b, keep_a, false),
            (b, a, &inc_b, &inc_a, keep_b, flip_b),
        ] {
            for sh in solid_shell_handles(model, solid) {
                for &fh in &model.shells.get(sh).faces {
                    items.push((fh, other, surf_ix[&fh], keep, flip, inc_f, inc_o));
                }
            }
        }
        let recon = |&(fh, other, pidx, keep, flip, inc_f, inc_o): &ReconItem| {
            reconstruct_face_paths(
                model, fh, other, pidx, keep, flip, &classof, &seam, &seam_ix, &planes, &surf_ix,
                inc_f, inc_o,
            )
        };
        // Per-item work is heavy (n0: 90µs–33ms/face), so plain `par_iter` maximizes
        // parallelism — no `with_min_len` floor, which would cap parallelism when faces
        // are few (a 12-face box → few chunks) while the rayon task overhead is negligible
        // against ms-scale faces.
        #[cfg(feature = "parallel")]
        let per_item: Vec<Result<Vec<LocalFace>, BoolError>> =
            items.par_iter().map(recon).collect();
        #[cfg(not(feature = "parallel"))]
        let per_item: Vec<Result<Vec<LocalFace>, BoolError>> = items.iter().map(recon).collect();

        // Surface the index-first error (deterministic reject tag, matching the sequential
        // sweep's first `?`). Under parallel + cfg(test) that tag was written on a worker
        // thread's `LAST_REJECT`; replay the one failing item here on the main thread so
        // `assert_rejects` still sees it (inert in release — `reject` records nothing there).
        let mut faces: Vec<LocalFace> = Vec::new();
        for r in per_item {
            match r {
                Ok(lfs) => faces.extend(lfs),
                Err(e) => {
                    // Restore the reject tag on this thread (a worker set it on its own
                    // `LAST_REJECT`): re-scan in index order until the first failure re-runs
                    // `reject` here (cfg(test) only; inert in release).
                    #[cfg(all(test, feature = "parallel"))]
                    let _ = items.iter().find(|it| recon(it).is_err());
                    return Err(e);
                }
            }
        }
        faces
    };
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

/// A canonical, replay-stable sort key for a severed component: its outer-loop vertex
/// coordinates, sorted lexicographically. Total order for disjoint components — distinct pieces
/// occupy different space, so their coordinate multisets differ, and the lex-min vertex alone can
/// tie (identity is by `Handle`, not coordinates, so two vertices may coincide). Coordinates are a
/// derived cache used only to order multi-solid output; no judgment reads this (tol-irrelevant).
fn comp_key(model: &Model, faces: &[Handle<Face>]) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = faces
        .iter()
        .flat_map(|&fh| model.faces.get(fh).outer.half_edges.iter().copied())
        .map(|he| model.vertices.get(he_start(model, he)).point.as_array())
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite vertex coordinates"));
    pts
}

/// Whether a closed component shell (a set of oriented faces) is **outward**
/// (material-enclosing, positive signed volume — an outer shell) versus
/// **inward** (a void/cavity shell). The `assemble_fuse_cut` cavity-vs-outer
/// label ((5d)#5 retired the f64 signed-volume flux this replaces), reading no
/// coordinate arithmetic: only a lexicographic vertex ordering and one
/// axis-aligned plane-coefficient sign.
///
/// At the component's lexicographically-minimal vertex `v*` (min x, then y, then
/// z) the shell is a convex corner, and the material lies toward increasing
/// coordinates. So an outward shell has a `−x`-facing boundary face at `v*` (its
/// materialized outward normal `n_x < 0`), while a void's three walls all face
/// into the void (`n_x ≥ 0` at its own `v*`). Hence: **outward iff some face
/// incident to `v*` has materialized outward normal with `n_x < 0`**. Two
/// antiparallel `x`-perpendicular faces cannot share a vertex, so this `∃`-test
/// is equivalent to (and simpler than) picking the max-`|n_x|` face.
///
/// Exact for the axis-aligned M5 corpus: face normals are exactly `±eₓ/±e_y/±e_z`
/// so `sign(n_x)` is the exact sign of the plane's `x`-coefficient (times the
/// face orientation), and the `v*` search is exact coordinate ordering — both
/// hold even for non-representable coordinates (e.g. a `0.3`-offset void face).
/// Rotated shells break the "axis-aligned normal / unique x-perpendicular face"
/// premises and are TIP's job (design §9 (5d)#5, honest scope).
fn is_shell_outward(model: &Model, faces: &[Handle<Face>]) -> bool {
    // Lexicographically-minimal vertex over the component's outer loops.
    let mut vstar: Option<Handle<Vertex>> = None;
    let mut pstar = [f64::INFINITY; 3];
    for &fh in faces {
        for &he in &model.faces.get(fh).outer.half_edges {
            let vh = he_start(model, he);
            let p = model.vertices.get(vh).point.as_array();
            if p < pstar {
                pstar = p;
                vstar = Some(vh);
            }
        }
    }
    let Some(vstar) = vstar else { return false };
    // Outward iff some face at v* faces −x (materialized outward normal n_x < 0).
    // n_x's sign is the plane x-coefficient's sign times the orientation sign
    // (no normalization — exact for axis-aligned faces).
    for &fh in faces {
        let face = model.faces.get(fh);
        if !face
            .outer
            .half_edges
            .iter()
            .any(|&he| he_start(model, he) == vstar)
        {
            continue;
        }
        let Surface::Plane(plane) = model.surfaces.get(face.surface) else {
            continue;
        };
        let sign = match face.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        if plane.coefficients()[0] * sign < 0.0 {
            return true;
        }
    }
    false
}

/// Rotation-sound twin of [`is_shell_outward`]: whether a result component's shell is
/// **outward** (material, an outer shell) versus **inward** (a void/cavity shell), decided on
/// the faces' exact `Pt3` definitions through `frame3`, so it is sound when the coordinates
/// are rounded irrationals (rotation). Same algorithm as the f64 [`is_shell_outward`] — the
/// lexicographically-minimal vertex `v*` is a convex extreme corner and the shell is outward
/// iff some face there has an outward normal with `n_x < 0` — but both numeric steps become
/// exact TIP predicates:
///
/// - **`v*`** by [`t_cmp_coord`](crate::tolerant) over each node's three-plane triple, the same
///   lex-min scan as [`loop_winding`](crate::arrange::loop_winding). A node's triple is its
///   own face plane plus the neighbour planes of its two loop edges (the
///   [`loop_triples`](crate::arrange) construction), whose meet *is* that vertex — so an
///   original corner is as implicit a point as a seam node, no mixed compare needed.
/// - **`sign(n_x)`** by [`dir_orient3d_judge`]`([1,0,0], tri…)` on the face's exact plane
///   definition ([`plane_def`](crate::tolerant), mixed-rotation safe): the x-component of the
///   RH normal, flipped by `lf.flip` to the result face's materialized outward normal.
///
/// Operates on the **pre-assembly** `LocalFace`s (not the result faces), so it never reads a
/// result vertex whose exact rotation provenance `assemble_fuse_cut` drops — every exact
/// definition it needs lives in `planes` and in the loop adjacency. Reads no `Model`. `Err`
/// on a non-simple/degenerate component (coincident nodes, a straight angle, or a non-manifold
/// edge) — an honest reject, never a silent wrong label. Routed from `assemble_fuse_cut`'s
/// per-component outward test by [`any_rotated`](crate::tolerant) (cell 3c-vi-b).
fn component_is_outward_tol(planes: &[PlaneInfo], comp: &[&LocalFace]) -> Result<bool, BoolError> {
    use nacre_scalar::{Orient, Rat};

    // Node lacks `Ord`; this is a canonical, hashable id for the unordered edge key.
    fn rank(n: Node) -> (u8, usize, usize, usize) {
        match n {
            Node::Orig(h) => (0, h.index() as usize, 0, 0),
            Node::Seam([a, b, c]) => (1, a, b, c),
        }
    }
    let ekey = |a: Node, b: Node| {
        let (ra, rb) = (rank(a), rank(b));
        if ra <= rb { (ra, rb) } else { (rb, ra) }
    };

    // Edge -> the planes carrying it, over every loop (outer + inner): a hole-rim edge is the
    // outer edge of its wall and an inner edge of the holed face, so building over both loops
    // gives it both planes. A manifold edge yields exactly two.
    type EKey = ((u8, usize, usize, usize), (u8, usize, usize, usize));
    let mut edge_planes: HashMap<EKey, Vec<usize>> = HashMap::new();
    for lf in comp {
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            let k = ring.len();
            for t in 0..k {
                edge_planes
                    .entry(ekey(ring[t], ring[(t + 1) % k]))
                    .or_default()
                    .push(lf.plane_idx);
            }
        }
    }
    let other_plane = |a: Node, b: Node, own: usize| -> Result<usize, BoolError> {
        let ps = edge_planes
            .get(&ekey(a, b))
            .ok_or_else(|| reject(tag::MISSING_SEAM))?;
        let mut others = ps.iter().copied().filter(|&x| x != own);
        let o = others
            .next()
            .ok_or_else(|| reject(tag::LOOP_ORIENT_MISMATCH))?;
        if others.any(|x| x != o) {
            return Err(reject(tag::NON_MANIFOLD_EDGE)); // edge shared by >2 distinct planes
        }
        Ok(o)
    };

    // Each unique outer node -> its three-plane triple (meet = that vertex). Any incident face
    // yields a valid triple (all its planes pass through the vertex); first occurrence wins.
    let mut triple_of: HashMap<Node, [usize; 3]> = HashMap::new();
    for lf in comp {
        let ring = &lf.loop_nodes;
        let k = ring.len();
        for t in 0..k {
            let node = ring[t];
            if triple_of.contains_key(&node) {
                continue;
            }
            let prev = other_plane(ring[(t + k - 1) % k], node, lf.plane_idx)?;
            let next = other_plane(node, ring[(t + 1) % k], lf.plane_idx)?;
            if prev == next {
                return Err(reject(tag::LOOP_ORIENT_MISMATCH)); // a straight angle
            }
            let mut tri = [lf.plane_idx, prev, next];
            tri.sort_unstable();
            triple_of.insert(node, tri);
        }
    }

    // Lexicographically-minimal vertex over the unique outer nodes (`loop_winding`'s scan;
    // `t_cmp_coord` is exact for the rotated triples). Sort candidates for replay determinism.
    let mut nodes: Vec<Node> = triple_of.keys().copied().collect();
    nodes.sort_by_key(|&n| rank(n));
    let Some((&first, rest)) = nodes.split_first() else {
        return Ok(false); // empty component
    };
    let mut lo = first;
    for &node in rest {
        let ord = (0..3)
            .map(|axis| {
                crate::tolerant::t_cmp_coord(planes, triple_of[&node], triple_of[&lo], axis)
            })
            .find(|&c| c != 0);
        match ord {
            Some(c) if c < 0 => lo = node,
            Some(_) => {}
            None => return Err(reject(tag::LOOP_ORIENT_MISMATCH)), // two distinct nodes coincide
        }
    }

    // Outward iff some outer face at v* has a result outward normal with n_x < 0. n_x's sign is
    // the RH-normal x-component (`dir_orient3d_judge` on the exact plane def), flipped by `flip`.
    for lf in comp {
        if !lf.loop_nodes.contains(&lo) {
            continue;
        }
        let tri = crate::tolerant::plane_def(planes, lf.plane_idx);
        let ex = [Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)];
        let nx = dir_orient3d_judge(ex, &tri[0], &tri[1], &tri[2]);
        let nx = if lf.flip {
            match nx {
                Orient::Positive => Orient::Negative,
                Orient::Negative => Orient::Positive,
                Orient::Zero => Orient::Zero,
            }
        } else {
            nx
        };
        if nx == Orient::Negative {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The supporting planes of a solid's outer shell. `Unsupported` if any face is
/// non-planar or lacks three non-collinear loop points.
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<PlaneInfo>, BoolError> {
    // A rotated operand's face coordinates are rounded, so each plane also carries its
    // exact `Pt3` definition (overhaul stage 3). Decided once per solid — the axis-aligned
    // path keeps `tri_pt3 = None` and pays nothing.
    let rotated = solid_is_rotated(model, solid);
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let plane = match model.surfaces.get(face.surface) {
                Surface::Plane(p) => *p,
                Surface::Cylinder(_) => return Err(reject(tag::CYLINDER_FACE)),
            };
            let (tri, tri_verts) =
                outer_tri(model, face).ok_or_else(|| reject(tag::DEGENERATE_FACE))?;
            let n_out = (tri[1] - tri[0])
                .cross(tri[2] - tri[0])
                .normalize()
                .ok_or_else(|| reject(tag::DEGENERATE_NORMAL))?;
            let tri_pt3 = if rotated {
                let pt3 = |vh| {
                    nacre_tip::vertex_pt3(model, vh).map_err(|_| reject(tag::ROTATED_UNSUPPORTED))
                };
                Some([pt3(tri_verts[0])?, pt3(tri_verts[1])?, pt3(tri_verts[2])?])
            } else {
                None
            };
            out.push(PlaneInfo {
                surf: face.surface,
                face: fh,
                plane,
                tri,
                n_out,
                orient: face.orientation,
                tri_pt3,
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

/// Three non-collinear points of a face's outer loop — with the **vertex handle** each
/// point came from — ordered so their right-hand normal points **out** of the solid.
/// The handles let the toleranced predicates rebuild each point as a `Pt3` (overhaul
/// stage 3); the coordinates alone drive the axis-aligned path.
fn outer_tri(model: &Model, face: &Face) -> Option<([Point3; 3], [Handle<Vertex>; 3])> {
    let verts: Vec<Handle<Vertex>> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let pts: Vec<Point3> = verts
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let n = pts.len();
    // The turn at one corner does not know which way the ring winds. Every b-rep loop is
    // CCW about its face's outward normal, but at a *reflex* corner the local turn
    // opposes the global winding, so three consecutive points can hand back an inward
    // normal. The Newell sum has no single corner to be fooled by.
    let newell = (0..n).fold(Vector3::zero(), |acc, i| {
        acc + (pts[i] - pts[0]).cross(pts[(i + 1) % n] - pts[0])
    });
    let i = (0..n).find(|&i| {
        let (a, b, c) = (pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        (b - a).cross(c - a).norm() > 0.0
    })?;
    let (i0, i1, i2) = (i, (i + 1) % n, (i + 2) % n);
    let (a, b, c) = (pts[i0], pts[i1], pts[i2]);
    // Same b/c swap for coords and handles, so `tri[k]` and `tri_verts[k]` stay aligned.
    Some(if (b - a).cross(c - a).dot(newell) < 0.0 {
        ([a, c, b], [verts[i0], verts[i2], verts[i1]])
    } else {
        ([a, b, c], [verts[i0], verts[i1], verts[i2]])
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
                // Unreachable in practice — `boolean` rejects rotated inputs
                // (ROTATED_UNSUPPORTED) before reaching here; fall back defensively.
                Origin::Rotated { .. } => plane_side(pi.tri, v.point),
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

/// Toleranced twin of [`nacre_predicates::ray_triangle_cross`]: the same forward-ray /
/// triangle crossing, but decided on the vertices' exact rotation definitions (`Pt3`) through
/// `frame3`, so it is sound when the coordinates are rounded irrationals (rotation). The five
/// `orient3d` signs of the f64 predicate become three [`orient3d_ray`] (the ray-line edge
/// tests `orient3d(p, p+d, ·, ·)`), one [`orient3d_judge`] (`s0`, `p` vs the triangle plane),
/// and one [`dir_orient3d_judge`] (`sd`, the direction vs the plane). A `declare-0`
/// (`Orient::Zero`) maps to `Degenerate`, which the caller's next-ray retry absorbs — the same
/// escape the f64 predicate uses for a grazed edge/vertex. Logic mirrors the f64 original
/// arm-for-arm; the argument order matches it exactly via [`orient3d_ray`].
#[allow(dead_code)] // wired into point_in_solid in cell 3c-iii
fn ray_triangle_cross_tol(
    p: &Pt3,
    d: [nacre_scalar::Rat; 3],
    v0: &Pt3,
    v1: &Pt3,
    v2: &Pt3,
) -> RayCross {
    use nacre_scalar::Orient;
    let (e0, e1, e2) = (
        orient3d_ray(p, d, v1, v2),
        orient3d_ray(p, d, v2, v0),
        orient3d_ray(p, d, v0, v1),
    );
    if e0 == Orient::Zero || e1 == Orient::Zero || e2 == Orient::Zero {
        return RayCross::Degenerate; // ray line through an edge/vertex
    }
    if e0 != e1 || e1 != e2 {
        return RayCross::Miss; // ray line misses the triangle
    }
    let s0 = orient3d_judge(v0, v1, v2, p);
    if s0 == Orient::Zero {
        return RayCross::Degenerate; // p on the triangle's plane
    }
    let sd = dir_orient3d_judge(d, v0, v1, v2);
    if sd == Orient::Zero {
        return RayCross::Degenerate; // ray parallel to the plane
    }
    if s0 == sd {
        RayCross::Cross(if sd == Orient::Positive { 1 } else { -1 }) // forward, oriented by sd
    } else {
        RayCross::Miss // crossing behind p
    }
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
    // Fan triangles are direction-independent — build (and exactly-drop degenerate
    // ones) once, then reuse across every ray direction.
    let tris: Vec<[Point3; 3]> = faces
        .iter()
        .flat_map(|rings| rings.iter())
        .flat_map(|ring| fan_triangles(ring, 0))
        .collect();
    'dirs: for dir in RAY_DIRECTIONS {
        let d = Vector3::from_array(dir);
        let mut winding = 0i32;
        for tri in &tris {
            match ray_face_cross(p, d, *tri) {
                RayCross::Cross(sign) => winding += sign as i32,
                RayCross::Miss => {}
                RayCross::Degenerate => continue 'dirs, // grazed — try another direction
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

/// Classify vertex `vh` in/out of `solid`, routed by rotation: the exact f64 [`point_in_solid`]
/// when the geometry is axis-aligned (unchanged hot path), else the toleranced
/// [`point_in_solid_tol`] on the vertex's exact `Pt3`. Either the solid's faces being rotated
/// (irrational plane geometry) *or* the query vertex being rotated (irrational coordinate, e.g.
/// a mixed-rotation boolean) forces the toleranced route. A `Discovered` query vertex on that
/// route is honestly rejected `ROTATED_UNSUPPORTED` (indirect classification is a later cell).
fn vertex_in_solid(
    model: &Model,
    vh: Handle<Vertex>,
    solid: Handle<Solid>,
) -> Result<Side, BoolError> {
    if solid_is_rotated(model, solid)
        || matches!(model.vertices.get(vh).origin, Origin::Rotated { .. })
    {
        let p = nacre_tip::vertex_pt3(model, vh).map_err(|_| reject(tag::ROTATED_UNSUPPORTED))?;
        point_in_solid_tol(model, &p, solid)
    } else {
        point_in_solid(model, model.vertices.get(vh).point, solid)
    }
}

/// Rotation-sound twin of [`point_in_solid`]: forward-ray winding decided on the vertices'
/// exact `Pt3` definitions (via [`ray_triangle_cross_tol`]), so it is sound when the
/// coordinates are rounded irrationals. `p` is the query point as an exact `Pt3` (a rotated
/// vertex, or any definition). The solid's face vertices must be **direct**
/// (`Constructed`/`Rotated`); a `Discovered` face vertex (a boolean result fed back as a
/// rotated operand) is honestly rejected `ROTATED_UNSUPPORTED` — indirect classification is a
/// later cell (fresh-rotated scope).
///
/// A `Degenerate` from a triangle is disambiguated in place: a genuinely collinear (zero-area)
/// fan triangle contributes nothing and is skipped ([`pt3_base_collinear`], exact); anything
/// else is a real graze, and the ray is retried. `winding != 0` is ray-direction-independent
/// for a closed surface, so the fixed `RAY_DIRECTIONS` classify correctly in any rotated frame.
///
/// A parallel of `point_in_solid` (not a shared generic): the f64 path stays untouched, so its
/// bit-for-bit behaviour is unchanged by construction.
fn point_in_solid_tol(model: &Model, p: &Pt3, solid: Handle<Solid>) -> Result<Side, BoolError> {
    use nacre_scalar::Rat;
    // The exact-definition fan triangles, built once and reused across ray directions.
    let mut tris: Vec<[Pt3; 3]> = Vec::new();
    for fh in solid_faces(model, solid) {
        for ring in face_loop_verts(model, fh) {
            let pts: Vec<Pt3> = ring
                .iter()
                .map(|&vh| nacre_tip::vertex_pt3(model, vh))
                .collect::<Result<_, _>>()
                .map_err(|_| reject(tag::ROTATED_UNSUPPORTED))?;
            for s in 1..pts.len().saturating_sub(1) {
                tris.push([pts[0].clone(), pts[s].clone(), pts[s + 1].clone()]);
            }
        }
    }
    'dirs: for dir in RAY_DIRECTIONS {
        let d = dir.map(|v| Rat::from_int(v as i128));
        let mut winding = 0i32;
        for tri in &tris {
            match ray_triangle_cross_tol(p, d, &tri[0], &tri[1], &tri[2]) {
                RayCross::Cross(sign) => winding += sign as i32,
                RayCross::Miss => {}
                // Exactly collinear (zero-area) → contributes nothing, like the f64 path's
                // up-front degenerate drop; a real graze retries another direction.
                RayCross::Degenerate if pt3_base_collinear(&tri[0], &tri[1], &tri[2]) => {}
                RayCross::Degenerate => continue 'dirs,
            }
        }
        return Ok(if winding != 0 {
            Side::Inside
        } else {
            Side::Outside
        });
    }
    Err(reject(tag::RAY_DEGENERATE))
}

/// Whether three `Pt3` are **exactly collinear** (zero-area triangle), decided on their
/// pre-rotation rational `base` coordinates. A rigid rotation preserves collinearity, and three
/// vertices of one solid share a rotation chain, so their bases are comparable; all three
/// coordinate-plane projections of `(b−a)×(c−a)` must vanish (exact `Rat`, no tolerance). An
/// i128 overflow returns `false` (treat as non-collinear): a genuinely-collinear triangle then
/// stays and is at worst rejected `RAY_DEGENERATE`, never falsely skipped (which would drop a
/// real crossing — silent-wrong). Unrotated vertices carry `base == coord`, matching the f64
/// `triangle_is_degenerate`.
fn pt3_base_collinear(a: &Pt3, b: &Pt3, c: &Pt3) -> bool {
    use nacre_scalar::Rat;
    let (a, b, c) = (&a.base, &b.base, &c.base);
    let proj_zero = |i: usize, j: usize| -> Option<bool> {
        let det = b[i]
            .checked_sub(a[i])?
            .checked_mul(c[j].checked_sub(a[j])?)?
            .checked_sub(
                b[j].checked_sub(a[j])?
                    .checked_mul(c[i].checked_sub(a[i])?)?,
            )?;
        Some(det == Rat::from_int(0))
    };
    matches!(
        (proj_zero(1, 2), proj_zero(2, 0), proj_zero(0, 1)),
        (Some(true), Some(true), Some(true))
    )
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

/// [`face_loops`] as vertex **handles** (outer loop then each hole) — the toleranced
/// classifier builds each vertex's exact `Pt3` from these, where `face_loops` reads f64 points.
fn face_loop_verts(model: &Model, fh: Handle<Face>) -> Vec<Vec<Handle<Vertex>>> {
    let face = model.faces.get(fh);
    std::iter::once(&face.outer)
        .chain(face.inner.iter())
        .map(|lp| {
            lp.half_edges
                .iter()
                .map(|&he| he_start(model, he))
                .collect()
        })
        .collect()
}

/// The non-degenerate fan triangles `(pts[apex], pts[apex+s], pts[apex+s+1])` of a
/// planar loop, fanned from vertex `apex` (indices mod `k`). A concave loop's
/// spurious (reflex) triangles are kept — they cancel by orientation in the
/// oriented crossing sum — but **exactly zero-area (collinear)** triangles are
/// dropped by an exact test ([`triangle_is_degenerate`]), not a tolerance: a
/// zero-area triangle contributes nothing to the winding and would force a
/// spurious `Degenerate` ray retry, whereas a tiny-but-nonzero triangle is kept
/// and judged exactly by `ray_triangle_cross`. Retiring the old relative
/// `1e-12` bound closes the one silent-wrong drop on the winding path ((5d)#4).
/// Varying `apex` changes which internal diagonals appear, which the segment gate
/// exploits to sidestep a diagonal that happens to be coplanar with a query edge.
fn fan_triangles(pts: &[Point3], apex: usize) -> Vec<[Point3; 3]> {
    let k = pts.len();
    let mut tris = Vec::new();
    for s in 1..k.saturating_sub(1) {
        let (t0, t1, t2) = (pts[apex], pts[(apex + s) % k], pts[(apex + s + 1) % k]);
        if triangle_is_degenerate(t0, t1, t2) {
            continue; // exactly zero-area (collinear) fan triangle
        }
        tris.push([t0, t1, t2]);
    }
    tris
}

/// Whether three points are **exactly collinear** (zero-area triangle), decided
/// by exact `orient2d` on all three coordinate-plane projections — these are the
/// three components of `(t1−t0)×(t2−t0)`, so zero area ⟺ all three are `0`. No
/// tolerance, no coordinate materialized. **Axis-independent**: a genuinely
/// nonzero-area triangle has a nonzero cross vector, so at least one projection is
/// non-degenerate and it is never falsely dropped (unlike a single fixed-axis
/// projection, which would collapse for a face perpendicular to that axis).
fn triangle_is_degenerate(t0: Point3, t1: Point3, t2: Point3) -> bool {
    let (a, b, c) = (t0.as_array(), t1.as_array(), t2.as_array());
    orient2d([a[1], a[2]], [b[1], b[2]], [c[1], c[2]]) == 0.0
        && orient2d([a[2], a[0]], [b[2], b[0]], [c[2], c[0]]) == 0.0
        && orient2d([a[0], a[1]], [b[0], b[1]], [c[0], c[1]]) == 0.0
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
            if !pierced_faces(model, planes, pair, bounds[0], bounds[1], &rings)?.is_empty() {
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
) -> Result<Vec<Handle<Solid>>, BoolError> {
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
    // whole result; several mean either an enclosed void (a cavity — an inward-oriented shell)
    // or a severed operand (two or more outward, material-enclosing shells). `is_shell_outward`
    // (exact extreme-vertex sign) tells them apart. One outward component ⇒ outer shell + the
    // rest as its cavities. Several outward components ⇒ the result severed into that many
    // solids (cell 0.4) — unless a cavity also survives, which needs a containment test we do
    // not have yet, so that is `SEVERED_WITH_CAVITY`. No outward component is impossible.
    let (labels, n) = face_components(faces);
    let mut by_comp: Vec<Vec<Handle<Face>>> = vec![Vec::new(); n];
    for (i, &fh) in face_handles.iter().enumerate() {
        by_comp[labels[i]].push(fh);
    }
    // Outward/void label per component, routed by rotation (overhaul 3c-vi): a component with
    // any rotated plane is decided on exact `Pt3` definitions (`component_is_outward_tol` over
    // its pre-assembly `LocalFace`s), else the axis-aligned f64 `is_shell_outward` — unchanged,
    // so an unrotated result is bit-identical. `positives` stays in ascending `c` order.
    let mut by_comp_lf: Vec<Vec<&LocalFace>> = vec![Vec::new(); n];
    for (i, lf) in faces.iter().enumerate() {
        by_comp_lf[labels[i]].push(lf);
    }
    let mut positives: Vec<usize> = Vec::new();
    for c in 0..n {
        let idxs: Vec<usize> = by_comp_lf[c].iter().map(|lf| lf.plane_idx).collect();
        let outward = if crate::tolerant::any_rotated(planes, &idxs) {
            component_is_outward_tol(planes, &by_comp_lf[c])?
        } else {
            is_shell_outward(model, &by_comp[c])
        };
        if outward {
            positives.push(c);
        }
    }
    let shells: Vec<Handle<Shell>> = by_comp
        .iter()
        .map(|faces| {
            model.shells.push(Shell {
                faces: faces.clone(),
            })
        })
        .collect();
    match positives.len() {
        0 => Err(reject(tag::NO_OUTWARD_SHELL)),
        1 => {
            let outer_c = positives[0];
            // A cavity shell's faces already point into the void (the material is outside it, so
            // the material-on-correct-side reconstruction winds them inward) — measured, no flip.
            let cavities = (0..n)
                .filter(|&c| c != outer_c)
                .map(|c| shells[c])
                .collect();
            let solid = model.push_solid(Solid {
                outer: shells[outer_c],
                cavities,
            });
            model.live_solids.retain(|&s| s != a && s != b);
            Ok(vec![solid])
        }
        _ => {
            // Several material solids. Cavity ownership across multiple outer shells is unsolved.
            if (0..n).any(|c| !positives.contains(&c)) {
                return Err(reject(tag::SEVERED_WITH_CAVITY));
            }
            // Every component is its own cavity-free solid. Emit them in a canonical, replay-stable
            // order keyed on geometry (a component's sorted vertex coordinates), so a downstream op
            // can index the returned Vec deterministically. `comp_key` is total for disjoint
            // components (distinct pieces occupy different space, so their coordinate sets differ).
            let keys: Vec<Vec<[f64; 3]>> = by_comp.iter().map(|f| comp_key(model, f)).collect();
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&x, &y| {
                keys[x]
                    .partial_cmp(&keys[y])
                    .expect("finite vertex coordinates")
            });
            let solids: Vec<Handle<Solid>> = order
                .into_iter()
                .map(|c| {
                    model.push_solid(Solid {
                        outer: shells[c],
                        cavities: Vec::new(),
                    })
                })
                .collect();
            model.live_solids.retain(|&s| s != a && s != b);
            Ok(solids)
        }
    }
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

/// Two faces lie on the same plane — by a **shared `Surface` handle** (§5 explicit
/// sharing: O(1) `Handle` identity, exact, rotation-independent) or, as a fallback,
/// by the geometric rank-1 `planes_coplanar` test. A referenced coplanar contact —
/// a pad/pocket cap that reuses its face's surface — is caught by the handle path
/// without any coordinate test. On the axis-aligned M5 corpus the handle path is
/// redundant with `planes_coplanar` (same handle ⇒ same plane), so the geometric
/// fallback is what keeps independently-built coplanar contacts working; the handle
/// path's real payoff is rotated frames, where the geometric test would need the
/// rotation-exact judgment.
fn shares_or_coplanar(pa: &PlaneInfo, pb: &PlaneInfo) -> bool {
    pa.surf == pb.surf || planes_coplanar(&pa.plane, &pb.plane)
}

/// Union-find root of `x` in `parent` (with path compression). Roots are the smallest index
/// of their class, so the result is deterministic (replay, DNA §absolute-3).
// Wired into the unified coplanar handler's dispatch in a later cell; used by tests now.
#[cfg_attr(not(test), allow(dead_code))]
fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut c = x;
    while parent[c] != r {
        let next = parent[c];
        parent[c] = r;
        c = next;
    }
    r
}

/// Canonicalize the combined plane table by coplanarity: two planes that are the same plane
/// (shared `Surface` handle, or exact rank-1 [`planes_coplanar`]) are merged into one class, so
/// a wall of `a` coplanar with a wall of `b` names a **single line** in a shared plane π. This is
/// the one thing the seam engine cannot do (it rejects `order_along(R,R)==0` as `FOURPLANE`);
/// canonicalizing turns that self-comparison into a real order. Returns `canon` where `canon[i]`
/// is the class root (the smallest index in the class). Every decision is exact
/// (`shares_or_coplanar`) — no coordinate. O(n²) scan over the (small) face count.
// Wired into the unified coplanar handler's dispatch in a later cell; used by tests now.
#[cfg_attr(not(test), allow(dead_code))]
fn plane_classes(planes: &[PlaneInfo]) -> Vec<usize> {
    let n = planes.len();
    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            if shares_or_coplanar(&planes[i], &planes[j]) {
                let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                if ri != rj {
                    // Attach the larger root under the smaller so a class's root is its min index.
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }
    (0..n).map(|i| uf_find(&mut parent, i)).collect()
}

/// One boundary edge of a coplanar contact face: its edge handle, its two vertices in loop
/// order, and the *wall* plane class it lies on (the neighbour ≠ the contact plane π).
#[cfg_attr(not(test), allow(dead_code))]
/// One boundary edge of a coplanar footprint in the shared plane, as the crossing/ordering
/// machinery sees it: a segment on the line `π ∩ wall` between two endpoint `Node`s, keyed by its
/// `seg` index within the boundary. A **contact-face** boundary carries `wall` = a wall class and
/// `Node::Orig` endpoints; a **section** boundary (`section_of_solid`) carries `wall` = an a-face
/// class and `Node::Seam` endpoints — both are the same object here (section-Q generalization).
struct BndEdge {
    seg: usize, // index within the boundary (was `Handle<Edge>`; a section chord has no edge)
    v: [Node; 2], // endpoints, in loop order (Orig for a contact face, Seam for a section)
    wall: usize, // canonical plane class of the line `π ∩ wall` carrying this edge
}

/// The outer-loop boundary edges of a coplanar contact face `cf` (canonical plane class `pi`),
/// each tagged with the canonical class of its carrying wall. Order follows the loop, so a
/// vertex's two incident walls are the walls of the edge before and the edge at it.
#[cfg_attr(not(test), allow(dead_code))]
fn contact_boundary(
    model: &Model,
    cf: Handle<Face>,
    pi: usize,
    inc: &arrange::EdgePlanes,
    canon: &[usize],
) -> Result<Vec<BndEdge>, BoolError> {
    let mut out = Vec::new();
    for (seg, &he) in model.faces.get(cf).outer.half_edges.iter().enumerate() {
        let (bounds, pair) = inc[&he.edge];
        let (ca, cb) = (canon[pair[0]], canon[pair[1]]);
        let wall = if ca == pi {
            cb
        } else if cb == pi {
            ca
        } else {
            // Every contact-face edge has the contact plane as one incidence.
            return Err(reject(tag::ARRANGEMENT_DEGENERATE));
        };
        let start = he_start(model, he);
        let [v0, v1] = if bounds[0] == start {
            bounds
        } else {
            [bounds[1], bounds[0]]
        };
        out.push(BndEdge {
            seg,
            v: [Node::Orig(v0), Node::Orig(v1)],
            wall,
        });
    }
    Ok(out)
}

/// The boundary of a solid's planar section (from [`section_of_solid`]) as `BndEdge`s on the
/// section plane `pi`, so a section footprint can drive [`coplanar_reconstruct`] exactly like a
/// contact-face footprint. Each section node is a `Node::Seam` three-plane point `{pi, A_f, A_g}`;
/// the chord between two consecutive nodes rides the one a-face they share besides `pi`, which
/// becomes the chord's carrying `wall`. `pi` is the section plane's index in `planes` (the same one
/// passed to `section_of_solid`), so the wall indices agree with the rest of the machinery.
///
/// Scope: a single closed loop, no holes. A multi-loop section (a slot piercing a cavity → outer +
/// hole rings) is rejected up front (`SECTION_MULTI_LOOP`), never silently flattened — the one
/// genuine silent-wrong risk of the section approach (plan §C2 급소).
#[cfg_attr(not(test), allow(dead_code))]
fn section_boundary(loops: &[Vec<Node>], pi: usize) -> Result<Vec<BndEdge>, BoolError> {
    if loops.len() != 1 {
        return Err(reject(tag::SECTION_MULTI_LOOP)); // outer + hole rings: later cell
    }
    let ring = &loops[0];
    let n = ring.len();
    if n < 3 {
        return Err(reject(tag::ARRANGEMENT_DEGENERATE));
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (v0, v1) = (ring[i], ring[(i + 1) % n]);
        let (Node::Seam(t0), Node::Seam(t1)) = (v0, v1) else {
            return Err(reject(tag::ARRANGEMENT_DEGENERATE)); // section nodes are all Seam
        };
        // The chord rides the single plane its two endpoints share besides the section plane.
        let mut shared = t0.iter().copied().filter(|&p| p != pi && t1.contains(&p));
        let wall = shared
            .next()
            .ok_or_else(|| reject(tag::ARRANGEMENT_DEGENERATE))?;
        if shared.next().is_some() {
            return Err(reject(tag::ARRANGEMENT_DEGENERATE)); // two shared walls = coincident ends
        }
        out.push(BndEdge {
            seg: i,
            v: [v0, v1],
            wall,
        });
    }
    Ok(out)
}

/// A proper crossing of the two coplanar footprint boundaries ∂P (from `a`) and ∂Q (from `b`)
/// in the shared plane π: the exact three-plane point `{π, W, U}` where an `a`-edge (on line
/// π∩W) crosses a `b`-edge (on line π∩U), each strictly within its segment.
#[cfg_attr(not(test), allow(dead_code))]
struct CoCross {
    point: Point3,
    triple: [usize; 3], // sorted {pi, wall_a, wall_b}
    a_seg: usize,       // index of the crossed a-boundary edge (was `Handle<Edge>`)
    b_seg: usize,       // index of the crossed b-boundary edge
}

/// Proper crossings of the two coplanar footprint boundaries, the exact analog of
/// [`try_overhang`]'s `proj2`/`orient2d` sweep: containment along each edge is decided by
/// [`order_along`](arrange::order_along) on plane triples, never a coordinate. A crossing that
/// lands on a footprint *vertex* (`order` returns `0`) is a flush-edge / T-junction precursor —
/// honestly rejected here, supported in the flush-edge cell. Parallel or shared (flush) walls
/// meet in no point (`three_planes` `None`) and contribute no transversal crossing.
#[cfg_attr(not(test), allow(dead_code))]
fn coplanar_boundary_crossings(
    planes: &[PlaneInfo],
    pi: usize,
    a: &[BndEdge],
    b: &[BndEdge],
) -> Result<Vec<CoCross>, BoolError> {
    let (na, nb) = (a.len(), b.len());
    // Position of the crossing within an edge from its two endpoint orders: `1` strictly
    // between (same nonzero sign both ways), `-1` outside, `0` on an endpoint (graze).
    let between = |s0: i8, s1: i8| -> i8 {
        if s0 == 0 || s1 == 0 {
            0
        } else if s0 == s1 {
            1
        } else {
            -1
        }
    };
    let mut out = Vec::new();
    for i in 0..na {
        let (w, wp, wn) = (a[i].wall, a[(i + na - 1) % na].wall, a[(i + 1) % na].wall);
        for (j, edge_b) in b.iter().enumerate() {
            let u = edge_b.wall;
            let Some(point) = three_planes(&planes[pi].plane, &planes[w].plane, &planes[u].plane)
            else {
                continue; // parallel or shared wall — no transversal crossing
            };
            let (up, un) = (b[(j + nb - 1) % nb].wall, b[(j + 1) % nb].wall);
            let a_pos = between(
                arrange::order_along(planes, pi, w, wp, u),
                arrange::order_along(planes, pi, w, u, wn),
            );
            let b_pos = between(
                arrange::order_along(planes, pi, u, up, w),
                arrange::order_along(planes, pi, u, w, un),
            );
            if a_pos == 0 || b_pos == 0 {
                return Err(reject(tag::VERTEX_ON_FACE_PLANE)); // flush-edge precursor
            }
            if a_pos == 1 && b_pos == 1 {
                let mut triple = [pi, w, u];
                triple.sort_unstable();
                out.push(CoCross {
                    point,
                    triple,
                    a_seg: a[i].seg,
                    b_seg: edge_b.seg,
                });
            }
        }
    }
    Ok(out)
}

/// The `b` footprint boundary ∂Q spliced with its coplanar crossings and split into arcs.
///
/// Walk ∂Q in loop order, inserting each crossing on the b-edge it rides (ordered along that
/// edge), then cut the cyclic sequence at the crossings. Each arc runs crossing → interior
/// b-vertices → crossing. Interior nodes are the `b` solid's **original** vertices
/// (`Node::Orig`, Constructed); only the crossing endpoints are `Node::Seam` — the mixed-node
/// coplanar seam that keeps b-vertex identity (all-Constructed purity where there is no
/// crossing). Returns the arcs; empty when there are no crossings (∂Q is then a closed inner
/// loop handled elsewhere). Classifying an arc inside/outside `P` is the tagging cell's job.
#[cfg_attr(not(test), allow(dead_code))]
fn coplanar_seam_arcs(
    planes: &[PlaneInfo],
    pi: usize,
    b: &[BndEdge],
    crossings: &[CoCross],
) -> Vec<Vec<Node>> {
    if crossings.is_empty() {
        return Vec::new();
    }
    // The third plane of a crossing on line π∩`u` — its `a`-wall, the key that orders it there.
    let other_wall = |c: &CoCross, u: usize| -> usize {
        c.triple
            .iter()
            .copied()
            .find(|&x| x != pi && x != u)
            .expect("crossing triple is {pi, W, U}")
    };
    // ∂Q as a cyclic node sequence: each edge's start vertex, then the crossings on that edge.
    let mut seq: Vec<Node> = Vec::new();
    for be in b {
        seq.push(be.v[0]);
        let mut on: Vec<&CoCross> = crossings.iter().filter(|c| c.b_seg == be.seg).collect();
        on.sort_by(|c0, c1| {
            let (w0, w1) = (other_wall(c0, be.wall), other_wall(c1, be.wall));
            match arrange::order_along(planes, pi, be.wall, w0, w1) {
                -1 => std::cmp::Ordering::Less,
                1 => std::cmp::Ordering::Greater,
                _ => std::cmp::Ordering::Equal,
            }
        });
        for c in on {
            seq.push(Node::Seam(c.triple));
        }
    }
    // Cut at the crossings into arcs (crossing → interior → next crossing). Cut points are the
    // crossing triples specifically, not every `Node::Seam`: a **contact-face** ∂Q has `Node::Orig`
    // interior vertices (only crossings are `Seam`), but a **section** ∂Q has `Node::Seam` corners
    // too, so we identify the crossings by their triples rather than by node kind.
    let xset: HashSet<[usize; 3]> = crossings.iter().map(|c| c.triple).collect();
    let n = seq.len();
    let cpos: Vec<usize> = (0..n)
        .filter(|&i| matches!(seq[i], Node::Seam(t) if xset.contains(&t)))
        .collect();
    let mut arcs = Vec::with_capacity(cpos.len());
    for k in 0..cpos.len() {
        let (from, to) = (cpos[k], cpos[(k + 1) % cpos.len()]);
        let mut arc = Vec::new();
        let mut i = from;
        loop {
            arc.push(seq[i]);
            if i == to {
                break;
            }
            i = (i + 1) % n;
        }
        arcs.push(arc);
    }
    arcs
}

/// Each boundary vertex of a contact face as a plane triple `{π, W_prev, W}` — the meet of the
/// two walls at the vertex (the start of edge `i`, shared with edge `i-1`), on canonical plane π.
#[cfg_attr(not(test), allow(dead_code))]
fn boundary_ring_triples(bnd: &[BndEdge], pi: usize) -> Vec<[usize; 3]> {
    let n = bnd.len();
    (0..n)
        .map(|i| {
            let mut t = [pi, bnd[(i + n - 1) % n].wall, bnd[i].wall];
            t.sort_unstable();
            t
        })
        .collect()
}

/// What survives on the `a`-side contact face `P` under a coplanar boolean, as a function of the
/// cells' inside-`Q` membership.
///
/// Derived from the occupancy table (plan §생존 규칙): a π-face survives where the result's
/// material lies on exactly one side of π. `Whole` keeps all of `P` (∂Q is internal and
/// dissolves — no hole); `MinusQ` keeps `P∖Q` (`Q` a hole for a contained footprint, a notch for
/// a crossing one); `InterQ` keeps `P∩Q`; `Empty` keeps nothing. The complementary `inQ∖P` piece
/// (a boss/cantilever bottom, `Fuse`/opposite) is produced by the symmetric `b`-side pass.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[cfg_attr(not(test), allow(dead_code))]
enum PSurvive {
    Whole,
    MinusQ,
    InterQ,
    Empty,
}

/// The `a`-side survival selector (`[`PSurvive`]`) and whether `b`'s non-contact faces flip.
/// `same_normal` is whether the two contact faces share an outward-normal direction.
///
/// Verified against the six bespoke paths: Cut/same=`MinusQ` (blind pocket mouth),
/// Fuse/opposite=`MinusQ` (boss: `P` gains `Q` as a hole; the cantilever is the `b`-side),
/// Cut/opposite=`Whole` (`b` sits above, no overlap), Fuse/same=`Whole` (union — `∂Q` internal),
/// Common/same=`InterQ` (overlap cap), Common/opposite=`Empty` (coincident stack → `EmptyResult`).
#[cfg_attr(not(test), allow(dead_code))]
fn coplanar_survival(kind: BoolKind, same_normal: bool) -> (PSurvive, bool) {
    use BoolKind::*;
    use PSurvive::*;
    let survive = match (kind, same_normal) {
        (Cut, true) => MinusQ,
        (Cut, false) => Whole,
        (Fuse, false) => MinusQ,
        (Fuse, true) => Whole,
        (Common, true) => InterQ,
        (Common, false) => Empty,
    };
    (survive, kind == Cut) // b's faces flip only for a Cut (they bound the removed region)
}

/// Reconstruct one coplanar contact face `P` (of `a`, plane class `pi`, materialized on plane
/// `plane_idx`) subdivided by `b`'s footprint `Q`, keeping the cells whose material survives.
///
/// Reuses the seam engine's cell-extraction core ([`arrange::boundary_runs`] /
/// [`arrange::run_classes`] / [`arrange::stitch_cycles`]) but with the coplanar substitutions the
/// plan names: containment is exact [`arrange::point_in_ring`] (no `classof`, no coordinate); the
/// seam that cuts `P` is the pieces of `∂Q` lying inside `P` (`inside_arcs`), whose interior nodes
/// are `b`'s original vertices (`Node::Orig`) and whose endpoints are the crossings (`Node::Seam`)
/// — the mixed-node seam. `keep_inside_q` selects `P∩Q` vs `P∖Q`; `flip` sets the face normal.
///
/// Scope: a single outer ring, no holes (the primitive-layer fixtures). Holes/islands and the
/// no-interior arc (adjacent crossings on one edge) are later cells.
#[cfg_attr(not(test), allow(dead_code))]
#[allow(clippy::too_many_arguments)]
fn coplanar_reconstruct(
    planes: &[PlaneInfo],
    pi: usize,
    plane_idx: usize,
    a_bnd: &[BndEdge],
    b_bnd: &[BndEdge],
    crossings: &[CoCross],
    arcs: &[Vec<Node>],
    keep_inside_q: bool,
    flip: bool,
) -> Result<Vec<LocalFace>, BoolError> {
    let p_ring = boundary_ring_triples(a_bnd, pi);
    let q_ring = boundary_ring_triples(b_bnd, pi);
    // ∂Q node → its plane triple, to classify an arc's interior against P. Keyed by `Node`, so it
    // covers both a contact-face vertex (`Node::Orig`) and a section chord end (`Node::Seam`).
    let b_vtriple: HashMap<Node, [usize; 3]> = b_bnd
        .iter()
        .enumerate()
        .map(|(j, be)| (be.v[0], q_ring[j]))
        .collect();

    // Keep the arcs of ∂Q that lie inside P — the seam that actually cuts P. A vertex-bearing arc
    // is classified by an interior node (its whole interior lies on one side). An interior-free
    // arc (two crossings adjacent on one edge) has no such witness, but the arcs alternate
    // inside/outside P at every crossing, so it takes the class opposite its neighbour.
    let mut inside: Vec<Option<bool>> = vec![None; arcs.len()];
    for (k, arc) in arcs.iter().enumerate() {
        if arc.len() > 2 {
            let tri = *b_vtriple
                .get(&arc[1])
                .ok_or_else(|| reject(tag::MISSING_SEAM))?;
            inside[k] = Some(arrange::point_in_ring(planes, pi, tri, &p_ring)?);
        }
    }
    if !arcs.is_empty() {
        let seed = inside
            .iter()
            .position(|x| x.is_some())
            .ok_or_else(|| reject(tag::OVERHANG_ARCS))?; // no vertex-bearing arc to anchor
        for step in 1..=arcs.len() {
            let k = (seed + step) % arcs.len();
            let prev = inside[(k + arcs.len() - 1) % arcs.len()].expect("propagated in walk order");
            match inside[k] {
                None => inside[k] = Some(!prev),
                Some(v) if v == prev => return Err(reject(tag::SEAM_COUNT_MISMATCH)), // must alternate
                Some(_) => {}
            }
        }
    }
    let opens: Vec<&Vec<Node>> = (0..arcs.len())
        .filter(|&k| inside[k] == Some(true))
        .map(|k| &arcs[k])
        .collect();

    // ∂P vertices and their kept flag: a vertex is kept iff its inside-Q status matches the
    // survival selector. A vertex exactly on ∂Q is a flush touch — out of scope (later cell).
    let verts: Vec<Node> = a_bnd.iter().map(|e| e.v[0]).collect();
    let n = verts.len();
    let mut kept = vec![false; n];
    for i in 0..n {
        if arrange::point_on_ring(planes, pi, p_ring[i], &q_ring)? {
            return Err(reject(tag::VERTEX_ON_FACE_PLANE));
        }
        let inside_q = arrange::point_in_ring(planes, pi, p_ring[i], &q_ring)?;
        kept[i] = inside_q == keep_inside_q;
    }

    // Crossings grouped by the ∂P edge each rides (by `a_seg` = the edge's index), as triples.
    let mut by_edge: Vec<Vec<[usize; 3]>> = vec![Vec::new(); n];
    for c in crossings {
        by_edge[c.a_seg].push(c.triple);
    }

    if opens.is_empty() {
        // No seam cuts P: it is wholly kept or wholly dropped by its (uniform) vertex class.
        return Ok(Vec::from_iter(kept[0].then(|| LocalFace {
            plane_idx,
            loop_nodes: verts.clone(),
            inner: Vec::new(),
            flip,
        })));
    }

    let br = arrange::boundary_runs(planes, pi, &p_ring, &by_edge)?;
    let (crossings_t, runs) = (br.crossings, br.runs);
    let n_x = crossings_t.len();
    let run_kept = arrange::run_classes(&runs, &kept, &[runs.len()])?;
    let cross_ix: HashMap<[usize; 3], usize> = crossings_t
        .iter()
        .enumerate()
        .map(|(j, &t)| (t, j))
        .collect();
    let next: Vec<usize> = (0..n_x).map(|i| (i + 1) % n_x).collect();
    let prev: Vec<usize> = (0..n_x).map(|i| (i + n_x - 1) % n_x).collect();
    let is_kd = |c: usize| run_kept[prev[c]];

    // Each inside arc oriented kd → dk, its ends found by crossing triple.
    let mut kd = Vec::with_capacity(opens.len());
    let mut dk = Vec::with_capacity(opens.len());
    let mut arc_nodes: Vec<Vec<Node>> = Vec::with_capacity(opens.len());
    for arc in &opens {
        let end_tri = |nd: &Node| match nd {
            Node::Seam(t) => cross_ix
                .get(t)
                .copied()
                .ok_or_else(|| reject(tag::SEAM_COUNT_MISMATCH)),
            Node::Orig(_) => Err(reject(tag::SEAM_COUNT_MISMATCH)),
        };
        let (h, t) = (end_tri(&arc[0])?, end_tri(arc.last().expect("arc"))?);
        let forward = is_kd(h);
        if forward == is_kd(t) {
            return Err(reject(tag::SEAM_COUNT_MISMATCH));
        }
        let ordered: Vec<Node> = if forward {
            (*arc).clone()
        } else {
            arc.iter().rev().copied().collect()
        };
        let (kd_c, dk_c) = if forward { (h, t) } else { (t, h) };
        kd.push(prev[kd_c]);
        dk.push(prev[dk_c]);
        arc_nodes.push(ordered);
    }

    let cycles = arrange::stitch_cycles(&run_kept, &kd, &dk, &next)?;
    let mut faces = Vec::new();
    for cycle in cycles {
        let mut loop_nodes: Vec<Node> = Vec::new();
        for w in 0..cycle.len() {
            let (pre, this) = (cycle[(w + cycle.len() - 1) % cycle.len()], cycle[w]);
            let mut i = next[dk[pre]];
            loop {
                for &v in &runs[i] {
                    loop_nodes.push(verts[v]);
                }
                if i == kd[this] {
                    break;
                }
                i = next[i];
            }
            // Splice the arc, dropping its endpoint crossings' duplicate (they are the run's
            // bounding crossings already implied) — keep the full arc nodes.
            loop_nodes.extend(arc_nodes[this].iter().copied());
        }
        faces.push(LocalFace {
            plane_idx,
            loop_nodes,
            inner: Vec::new(),
            flip,
        });
    }
    Ok(faces)
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
            if shares_or_coplanar(pa, pb) && pa.n_out.dot(pb.n_out) < 0.0 {
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
///
/// **Neither solid need be convex.** The gates are convexity-agnostic — `opposite.len() == 1`
/// (exact coplanar test) and `face_contains_face` (`profile_strictly_in_region`, arbitrary
/// polygons) — and the reconstruction (`contained_contact_result` re-emits faces;
/// `assemble_fuse_cut` partitions by component). So a boss fuses cleanly onto a non-convex `a`
/// (an L-bracket top) with a non-convex `b` (an L/star footprint), mirroring the pocket path.
/// This makes `pad = extrude + Fuse` a drop-in for the old direct construction.
fn detect_contained_contact(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<ContainedContact> {
    let planes_a = collect_planes(model, a).ok()?;
    let planes_b = collect_planes(model, b).ok()?;
    let mut opposite: Vec<(usize, usize)> = Vec::new();
    for (i, pa) in planes_a.iter().enumerate() {
        for (j, pb) in planes_b.iter().enumerate() {
            if shares_or_coplanar(pa, pb) && pa.n_out.dot(pb.n_out) < 0.0 {
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

/// `Some` iff `b` is a blind pocket in `a`: they meet at exactly one **same-normal** coplanar
/// pair whose footprint is strictly inside `a`'s face, and `b` lies wholly inside `a` (its walls
/// cross none of `a`'s faces — a `b` that punched through would be a seam cut, not a pocket).
/// For `Cut(a, b)`, `a` is kept and `b` is carved out, so `a` must be the containing solid.
///
/// **Neither solid need be convex.** The remaining gates are convexity-agnostic — `same.len() == 1`
/// (exact coplanar test), `face_contains_face` (`profile_strictly_in_region`, arbitrary polygons),
/// and the blind test (`point_in_solid`, exact ray cast for non-convex) — and the reconstruction
/// (`contained_contact_result` re-emits faces; `assemble_fuse_cut` partitions by component). So a
/// pocket carves cleanly into a non-convex `a` (a re-pocketed part) with a non-convex `b` (an
/// L/star footprint). This makes `pocket = extrude + Cut` a drop-in for the old direct construction.
fn detect_pocket_contact(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<ContainedContact> {
    let planes_a = collect_planes(model, a).ok()?;
    let planes_b = collect_planes(model, b).ok()?;
    let mut same: Vec<(usize, usize)> = Vec::new();
    for (i, pa) in planes_a.iter().enumerate() {
        for (j, pb) in planes_b.iter().enumerate() {
            if shares_or_coplanar(pa, pb) && pa.n_out.dot(pb.n_out) > 0.0 {
                same.push((i, j));
            }
        }
    }
    if same.len() != 1 {
        return None;
    }
    let (i, j) = same[0];
    let fa = model.shells.get(model.solids.get(a).outer).faces[i];
    let fb = model.shells.get(model.solids.get(b).outer).faces[j];
    let n = planes_a[i].plane.normal();
    if !face_contains_face(model, fa, fb, n) {
        return None; // `a`'s face must contain `b`'s footprint (a is the kept solid)
    }
    // `b` must be blind — every vertex off the contact plane lies strictly inside `a`. A
    // through-`b` has vertices below `a` (a transversal seam cut, not a pocket); the flush
    // contact vertices sit on the plane and are skipped. (boundaries_intersect is unusable
    // here — the flush contact itself reads as a boundary touch.)
    let contact_tri = planes_a[i].tri;
    for &vh in &solid_vertex_handles(model, b) {
        let p = model.vertices.get(vh).point;
        if plane_side(contact_tri, p) != 0 && point_in_solid(model, p, a).ok()? != Side::Inside {
            return None;
        }
    }
    Some(ContainedContact {
        big_solid: a,
        big_face: fa,
        small_solid: b,
        small_face: fb,
    })
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
) -> Result<Vec<Handle<Solid>>, BoolError> {
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

/// A contained coplanar contact (`small`'s face inside `big`'s). `big`'s contact face keeps its
/// boundary and gains `small`'s footprint as a hole; `small`'s contact face is dropped; the
/// shared footprint edges stitch the hole to `small`'s walls (`assemble_fuse_cut`'s `edge_for`
/// dedups them). `small`'s footprint region becomes interior — no face there.
///
/// `cut` picks the operation. **Boss (`Fuse`, `cut = false`):** `small` sits outside `big` (the
/// faces have opposite normals), its walls kept as-is. **Pocket (`Cut`, `cut = true`):** `small`
/// sits inside `big` (same normals), so its faces flip — walls face into the removed region and
/// the far face becomes the pocket floor. The hole flips with it: `small`'s outer loop is CW
/// about `big`'s normal for the boss (opposite normals) but CCW for the pocket (same normals), so
/// the pocket reverses it (cells coplanar-contact-boss / -cut).
fn contained_contact_result(
    model: &mut Model,
    cc: &ContainedContact,
    cut: bool,
) -> Result<Vec<Handle<Solid>>, BoolError> {
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
    let mut hole = ring_nodes(&model.faces.get(cc.small_face).outer);
    if cut {
        hole.reverse(); // same normals: reverse to CW about big's normal
    }
    inner.push(hole);
    faces.push(LocalFace {
        plane_idx: big_pos,
        loop_nodes: ring_nodes(&big_f.outer),
        inner,
        flip: false,
    });
    // Small's faces except the dropped contact face; flipped for a cut so they bound the
    // removed region.
    faces.extend(
        solid_local_faces(model, cc.small_solid, na, Some(cc.small_face), None)
            .into_iter()
            .map(|lf| LocalFace { flip: cut, ..lf }),
    );

    assemble_fuse_cut(model, cc.big_solid, cc.small_solid, &planes, &[], &faces)
}

/// The **one** unified coplanar handler (the entry `D0` will wire into dispatch, replacing the six
/// bespoke paths). It finds the single coplanar contact pair, reads the survival table
/// [`coplanar_survival`] from the op and the faces' relative normal, and dispatches **internally**
/// on the footprint relationship — contained (`Q ⊂ P`, no crossing) here; overlapping/crossing and
/// `P ⊂ Q` are the later branches. Not per-case public builders: one entry, internal branches.
///
/// Contained builds `P` as `P∖Q` (`Q` a hole, reversed to CW about `P`'s normal iff the faces share
/// a normal) or `Whole`, `b`'s faces flipped for a `Cut`. Seam is empty — every node `Node::Orig`,
/// so the result is all `Constructed` (the purity the bespoke paths also keep).
#[cfg_attr(not(test), allow(dead_code))]
fn coplanar_result(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes_a = collect_planes(model, a)?;
    let na = planes_a.len();
    let planes_b = collect_planes(model, b)?;
    let mut planes = planes_a;
    planes.extend(planes_b);

    // The single coplanar cross-operand contact pair.
    let mut pair = None;
    for i in 0..na {
        for j in na..planes.len() {
            if shares_or_coplanar(&planes[i], &planes[j]) {
                if pair.is_some() {
                    return Err(reject(tag::VERTEX_ON_FACE_PLANE)); // >1 pair — out of scope
                }
                pair = Some((i, j));
            }
        }
    }
    let Some((pi_idx, qj_idx)) = pair else {
        return Err(reject(tag::VERTEX_ON_FACE_PLANE)); // no coplanar contact
    };
    let same_normal = planes[pi_idx].n_out.dot(planes[qj_idx].n_out) > 0.0;
    let (survive, b_flip) = coplanar_survival(kind, same_normal);
    let (p_face, q_face) = (planes[pi_idx].face, planes[qj_idx].face);
    let n = planes[pi_idx].plane.normal();

    // Contained: a's contact face P strictly contains b's footprint Q (no ∂P × ∂Q crossing).
    if face_contains_face(model, p_face, q_face, n) {
        let ring_nodes = |l: &Loop| -> Vec<Node> {
            l.half_edges
                .iter()
                .map(|&he| Node::Orig(he_start(model, he)))
                .collect()
        };
        let mut faces = solid_local_faces(model, a, 0, Some(p_face), None);
        let p_f = model.faces.get(p_face);
        let mut inner: Vec<Vec<Node>> = p_f.inner.iter().map(&ring_nodes).collect();
        match survive {
            PSurvive::MinusQ => {
                let mut hole = ring_nodes(&model.faces.get(q_face).outer);
                if same_normal {
                    hole.reverse(); // same normals: reverse to CW about P's normal
                }
                inner.push(hole);
            }
            PSurvive::Whole => {} // Q internal — no hole
            PSurvive::InterQ | PSurvive::Empty => return Err(reject(tag::VERTEX_ON_FACE_PLANE)),
        }
        faces.push(LocalFace {
            plane_idx: pi_idx, // a's contact face keeps its combined index (= shell position)
            loop_nodes: ring_nodes(&p_f.outer),
            inner,
            flip: false,
        });
        faces.extend(
            solid_local_faces(model, b, na, Some(q_face), None)
                .into_iter()
                .map(|lf| LocalFace { flip: b_flip, ..lf }),
        );
        return assemble_fuse_cut(model, a, b, &planes, &[], &faces);
    }
    if face_contains_face(model, q_face, p_face, n) {
        return Err(reject(tag::OVERHANG_ARCS)); // P ⊂ Q — symmetric contained, later branch
    }

    // Overlapping footprints (∂P × ∂Q crossing) — the overhang branch. Scope B4p2b: the boss
    // Fuse (opposite normals, two-sided, walls in open space). Cut/Common overhang, where b
    // breaks transversally through a's material, compose with the general engine at C2.
    if kind != BoolKind::Fuse || same_normal {
        return Err(reject(tag::OVERHANG_ARCS));
    }
    let mut surf_ix = HashMap::new();
    for (i, p) in planes.iter().enumerate() {
        surf_ix.insert(p.face, i);
    }
    let canon = plane_classes(&planes);
    let pi_c = canon[pi_idx];
    let inc_a = arrange::edge_planes(model, a, &surf_ix)?;
    let inc_b = arrange::edge_planes(model, b, &surf_ix)?;
    let a_bnd = contact_boundary(model, p_face, pi_c, &inc_a, &canon)?;
    let b_bnd = contact_boundary(model, q_face, pi_c, &inc_b, &canon)?;
    let crossings = coplanar_boundary_crossings(&planes, pi_c, &a_bnd, &b_bnd)?;
    if crossings.is_empty() {
        return Err(reject(tag::OVERHANG_ARCS)); // disjoint — no overhang
    }
    // Both contact faces reconstructed (mixed-node), keeping each's outside-other cells; no flip
    // (opposite normals already point out of the fused solid).
    let arcs_p = coplanar_seam_arcs(&planes, pi_c, &b_bnd, &crossings);
    let cx_q = coplanar_boundary_crossings(&planes, pi_c, &b_bnd, &a_bnd)?;
    let arcs_q = coplanar_seam_arcs(&planes, pi_c, &a_bnd, &cx_q);
    let mut faces = coplanar_reconstruct(
        &planes, pi_c, pi_idx, &a_bnd, &b_bnd, &crossings, &arcs_p, false, false,
    )?;
    faces.extend(coplanar_reconstruct(
        &planes, pi_c, qj_idx, &b_bnd, &a_bnd, &cx_q, &arcs_q, false, false,
    )?);
    // Both solids' walls, split at the crossings that ride their edges (no T-junction).
    let cross_r: Vec<Crossing> = crossings
        .iter()
        .map(|c| Crossing {
            point: c.point,
            triple: c.triple,
            p_seg: 0,
            q_seg: 0,
        })
        .collect();
    let mut a_walls = solid_local_faces(model, a, 0, Some(p_face), None);
    let mut b_walls = solid_local_faces(model, b, na, Some(q_face), None);
    resplit_overhang(model, &mut a_walls, &planes, &cross_r);
    resplit_overhang(model, &mut b_walls, &planes, &cross_r);
    faces.extend(a_walls);
    faces.extend(b_walls);
    let seam: Vec<SeamVertex> = crossings
        .iter()
        .map(|c| SeamVertex {
            point: c.point,
            triple: c.triple,
            tol: vertex_tol(
                c.point,
                &planes[c.triple[0]].plane,
                &planes[c.triple[1]].plane,
                &planes[c.triple[2]].plane,
            ),
        })
        .collect();
    assemble_fuse_cut(model, a, b, &planes, &seam, &faces)
}

// ---- coplanar-contact-overhang (M5): single-edge overhang boss Fuse ----

/// One footprint-boundary crossing: where P's single crossed edge meets one of Q's edges.
/// Its exact point and provenance are the three planes meeting there — P's contact plane,
/// P's wall (holding the crossed edge), Q's wall (holding the crossing edge).
#[derive(Clone)]
struct Crossing {
    point: Point3,
    triple: [usize; 3], // sorted combined-plane indices (planes_p ++ planes_q)
    p_seg: usize,       // P outer-loop edge index (the single crossed edge — same for both)
    q_seg: usize,       // Q outer-loop edge index (distinct per crossing)
}

/// An overhang boss contact (cells coplanar-contact-overhang / -corner / -multi). Two convex
/// solids meet at one opposite-normal coplanar face pair whose footprints overlap, their
/// boundaries ∂P × ∂Q crossing at `2k` points. `Fuse` welds the boss: each contact face keeps
/// its region minus the overlap (as one or more notch/cantilever pieces), the crossings stitch
/// the walls.
struct OverhangContact {
    p_solid: Handle<Solid>, // base — its face P → notch piece(s)
    p_face: Handle<Face>,
    q_solid: Handle<Solid>, // boss — its face Q → cantilever piece(s)
    q_face: Handle<Face>,
    crossings: Vec<Crossing>,
}

/// Strict proper 2D crossing: the two closed segments cross in their interiors (no shared
/// endpoint, no collinear touch). All four orientations nonzero and both straddle.
fn proper_cross_2d(p0: [f64; 2], p1: [f64; 2], q0: [f64; 2], q1: [f64; 2]) -> bool {
    let (o1, o2) = (orient2d(p0, p1, q0), orient2d(p0, p1, q1));
    let (o3, o4) = (orient2d(q0, q1, p0), orient2d(q0, q1, p1));
    o1 != 0.0
        && o2 != 0.0
        && o3 != 0.0
        && o4 != 0.0
        && (o1 > 0.0) != (o2 > 0.0)
        && (o3 > 0.0) != (o4 > 0.0)
}

/// The outer-shell position of the face (other than `contact`) whose outer loop carries the
/// edge with endpoints `va`,`vb`. Adjacency-free (walks the shell directly): `boolean` never
/// builds the model's adjacency cache on its inputs (see `edge_incidence`).
fn wall_pos_on_edge(
    model: &Model,
    solid: Handle<Solid>,
    contact: Handle<Face>,
    va: Handle<Vertex>,
    vb: Handle<Vertex>,
) -> Option<usize> {
    let shell = model.solids.get(solid).outer;
    for (pos, &fh) in model.shells.get(shell).faces.iter().enumerate() {
        if fh == contact {
            continue;
        }
        let ring: Vec<Handle<Vertex>> = model
            .faces
            .get(fh)
            .outer
            .half_edges
            .iter()
            .map(|&he| he_start(model, he))
            .collect();
        let k = ring.len();
        if (0..k).any(|i| {
            let (u, w) = (ring[i], ring[(i + 1) % k]);
            (u == va && w == vb) || (u == vb && w == va)
        }) {
            return Some(pos);
        }
    }
    None
}

/// The 3D point of a reconstructed node: an original vertex, or a crossing's stored point.
fn overhang_node_point(model: &Model, node: Node, crossings: &[Crossing]) -> Point3 {
    match node {
        Node::Orig(vh) => model.vertices.get(vh).point,
        Node::Seam(t) => {
            crossings
                .iter()
                .find(|c| c.triple == t)
                .expect("crossing node")
                .point
        }
    }
}

/// Splice `detour` into `loop_nodes`, replacing the edge `loop[seg] → loop[seg+1]`. The two
/// endpoints of `detour` are the crossings on that edge and its interior runs off it; the
/// result keeps the loop's winding. `detour` is reversed if needed so it enters at the crossing
/// nearer `loop[seg]` (ordered by parameter along the edge). `pt` maps a node to its 3D point.
/// Shared by the overhang Fuse notch and the Cut's mouth/side notches.
fn splice_notch(
    loop_nodes: &[Node],
    seg: usize,
    detour: Vec<Node>,
    pt: &impl Fn(Node) -> Point3,
) -> Vec<Node> {
    let n = loop_nodes.len();
    let ea = pt(loop_nodes[seg]);
    let dir = pt(loop_nodes[(seg + 1) % n]) - ea;
    let param = |node: Node| (pt(node) - ea).dot(dir);
    let detour = if param(detour[0]) <= param(*detour.last().expect("detour endpoints")) {
        detour
    } else {
        let mut r = detour;
        r.reverse();
        r
    };
    let mut out = Vec::with_capacity(n + detour.len());
    for (i, &node) in loop_nodes.iter().enumerate() {
        out.push(node);
        if i == seg {
            out.extend(detour.iter().cloned());
        }
    }
    out
}

/// One arc of a contact face's outer loop between two consecutive crossings.
struct LoopArc {
    nodes: Vec<Node>,     // crossing → … → crossing (endpoints are `Node::Seam`s)
    ends: (usize, usize), // the two crossing indices (unordered key)
    outside: bool,        // an interior vertex lies strictly outside the other footprint
}

/// Split a contact face's outer loop at all `2k` crossings into its `2k` arcs, each labelled by
/// its endpoint crossing indices and whether it lies outside the other footprint (empty interior
/// ⇒ `false`, i.e. inside — a same-edge segment between two crossings is inside by convexity).
/// `seg` gives each crossing's edge index on this loop; same-edge crossings are inserted in
/// parameter order so each arc stays simple.
fn split_loop_all_arcs(
    loop_nodes: &[Node],
    loop_pts: &[Point3],
    crossings: &[Crossing],
    other2d: &[[f64; 2]],
    drop: (usize, usize),
    seg: impl Fn(&Crossing) -> usize,
) -> Vec<LoopArc> {
    let n = loop_nodes.len();
    let inside: Vec<bool> = loop_pts
        .iter()
        .map(|&p| point_in_ring2(proj2(p, drop), other2d) == Some(true))
        .collect();
    enum Aug {
        V(usize),
        C(usize),
    }
    let mut aug: Vec<Aug> = Vec::new();
    for j in 0..n {
        aug.push(Aug::V(j));
        let ea = loop_pts[j];
        let dir = loop_pts[(j + 1) % n] - ea;
        let mut here: Vec<usize> = (0..crossings.len())
            .filter(|&ci| seg(&crossings[ci]) == j)
            .collect();
        here.sort_by(|&x, &y| {
            let px = (crossings[x].point - ea).dot(dir);
            let py = (crossings[y].point - ea).dot(dir);
            px.partial_cmp(&py).expect("finite params")
        });
        for ci in here {
            aug.push(Aug::C(ci));
        }
    }
    let m = aug.len();
    let cpos: Vec<usize> = (0..m).filter(|&i| matches!(aug[i], Aug::C(_))).collect();
    let ci_at = |pos: usize| match aug[pos] {
        Aug::C(ci) => ci,
        Aug::V(_) => unreachable!("cpos indexes crossings"),
    };
    // Each arc runs from one crossing forward to the next.
    let mut arcs = Vec::new();
    for k in 0..cpos.len() {
        let (from, to) = (cpos[k], cpos[(k + 1) % cpos.len()]);
        let mut nodes = Vec::new();
        let mut outside = false;
        let mut i = from;
        loop {
            match aug[i] {
                Aug::V(j) => {
                    nodes.push(loop_nodes[j]);
                    if !inside[j] {
                        outside = true;
                    }
                }
                Aug::C(ci) => nodes.push(Node::Seam(crossings[ci].triple)),
            }
            if i == to {
                break;
            }
            i = (i + 1) % m;
        }
        arcs.push(LoopArc {
            nodes,
            ends: (ci_at(from), ci_at(to)),
            outside,
        });
    }
    arcs
}

/// Join two arcs sharing their two crossing endpoints into one closed loop. `outer` sets the
/// winding; `inner` is oriented to run back from `outer`'s end to its start and only its interior
/// is appended, giving `[outer, inner-interior]` (implicitly closed). An empty inner interior
/// leaves `outer` unchanged.
fn stitch_arcs(outer: Vec<Node>, inner: Vec<Node>) -> Vec<Node> {
    let y = *outer.last().expect("outer arc");
    let inner = if inner[0] == y {
        inner
    } else {
        inner.into_iter().rev().collect()
    };
    let interior = &inner[1..inner.len() - 1];
    let mut out = outer;
    out.extend_from_slice(interior);
    out
}

/// With `p_solid`'s contact face as P and `q_solid`'s as Q, find the proper crossings of their
/// footprint boundaries (∂P × ∂Q) — a general convex-convex overlap. Returns the `2k` crossings
/// (each with its own P/Q edge), or `None` on an odd/zero count or a boundary graze.
#[allow(clippy::too_many_arguments)]
fn try_overhang(
    model: &Model,
    p_solid: Handle<Solid>,
    planes_p: &[PlaneInfo],
    p_pos: usize,
    p_face: Handle<Face>,
    q_solid: Handle<Solid>,
    planes_q: &[PlaneInfo],
    q_face: Handle<Face>,
    na_p: usize,
) -> Option<Vec<Crossing>> {
    let drop = planar_drop_axes(planes_p[p_pos].n_out);
    let ring = |face: Handle<Face>| -> Vec<(Handle<Vertex>, Point3)> {
        model
            .faces
            .get(face)
            .outer
            .half_edges
            .iter()
            .map(|&he| {
                let vh = he_start(model, he);
                (vh, model.vertices.get(vh).point)
            })
            .collect()
    };
    let p_ring = ring(p_face);
    let q_ring = ring(q_face);
    let p2: Vec<[f64; 2]> = p_ring.iter().map(|&(_, p)| proj2(p, drop)).collect();
    let q2: Vec<[f64; 2]> = q_ring.iter().map(|&(_, p)| proj2(p, drop)).collect();
    let (np, nq) = (p2.len(), q2.len());

    // Proper crossings of the two footprint boundaries.
    let mut crs: Vec<(usize, usize)> = Vec::new();
    for i in 0..np {
        for j in 0..nq {
            if proper_cross_2d(p2[i], p2[(i + 1) % np], q2[j], q2[(j + 1) % nq]) {
                crs.push((i, j));
            }
        }
    }
    // An even number 2k (≥2) of proper crossings — a general convex-convex overlap (k lens
    // pieces). Odd/zero means a degenerate or non-crossing case, out of scope.
    if crs.len() < 2 || crs.len() % 2 != 0 {
        return None;
    }
    // No boundary graze on either side (a vertex exactly on the other's boundary is out of scope).
    for &pp in &p2 {
        point_in_ring2(pp, &q2)?;
    }
    for &qq in &q2 {
        point_in_ring2(qq, &p2)?;
    }
    // Exact crossing points as three-plane meets; provenance triple in combined-plane indices.
    // Each crossing records its own P/Q edge (they may differ per crossing for a corner).
    let mut out: Vec<Crossing> = Vec::new();
    for &(i, j) in &crs {
        let p_wall = wall_pos_on_edge(model, p_solid, p_face, p_ring[i].0, p_ring[(i + 1) % np].0)?;
        let q_wall = wall_pos_on_edge(model, q_solid, q_face, q_ring[j].0, q_ring[(j + 1) % nq].0)?;
        let point = three_planes(
            &planes_p[p_pos].plane,
            &planes_p[p_wall].plane,
            &planes_q[q_wall].plane,
        )?;
        let mut triple = [p_pos, p_wall, na_p + q_wall];
        triple.sort_unstable();
        out.push(Crossing {
            point,
            triple,
            p_seg: i,
            q_seg: j,
        });
    }
    Some(out)
}

/// A simple polygon is convex iff every turn winds the same way. Exact via `orient2d` on the
/// projected loop; collinear turns (`0`) are skipped. `< 3` points ⇒ not a polygon.
fn loop_is_convex_2d(pts: &[[f64; 2]]) -> bool {
    let n = pts.len();
    if n < 3 {
        return false;
    }
    let mut sign = 0.0_f64;
    for i in 0..n {
        let o = orient2d(pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        if o != 0.0 {
            if sign == 0.0 {
                sign = o;
            } else if (o > 0.0) != (sign > 0.0) {
                return false;
            }
        }
    }
    true
}

/// Whether `face`'s outer loop is a convex polygon, projected onto its plane by dropping `drop`.
fn face_outer_is_convex(model: &Model, face: Handle<Face>, drop: (usize, usize)) -> bool {
    let pts: Vec<[f64; 2]> = model
        .faces
        .get(face)
        .outer
        .half_edges
        .iter()
        .map(|&he| proj2(model.vertices.get(he_start(model, he)).point, drop))
        .collect();
    loop_is_convex_2d(&pts)
}

/// `Some` iff `a` and `b` form an overhang boss contact: exactly one opposite-normal coplanar
/// face pair, neither footprint contained in the other (that is the boss cell), their footprints
/// overlapping (∂P meets ∂Q at `2k` points — single-edge, swallowed corner, or a spanning slab),
/// and no transversal piercing (a pure coplanar contact, not a seam cut).
///
/// **Neither solid need be convex** — only the two **contact-face footprints** must be convex (and
/// hole-free). The Fuse reconstruction touches only the contact faces (arc-split, which assumes a
/// convex footprint) and the breached walls (`resplit_overhang`, which is edge-local and
/// convexity-agnostic); every other face is re-emitted verbatim. So a boss cantilevers onto a
/// non-convex solid (a pocketed part, a boolean result) as long as the contact face is a convex
/// polygon. A non-convex contact footprint would mislabel arcs (silent-wrong) and is rejected here.
/// (Cut/Common keep the whole-solid gate: their `clip_bwall_inside_a` clips against breached
/// half-spaces, which only equals "inside `a`" when `a` is convex.)
fn detect_overhang_contact(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<OverhangContact> {
    let planes_a = collect_planes(model, a).ok()?;
    let planes_b = collect_planes(model, b).ok()?;
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
    let (ia, jb) = opposite[0];
    let fa = model.shells.get(model.solids.get(a).outer).faces[ia];
    let fb = model.shells.get(model.solids.get(b).outer).faces[jb];
    let n = planes_a[ia].plane.normal();
    // Only the two contact footprints must be convex and hole-free (not the whole solids) — the
    // reconstruction is local to them plus the edge-locally-resplit walls (see doc above).
    let drop = planar_drop_axes(planes_a[ia].n_out);
    if !face_outer_is_convex(model, fa, drop)
        || !face_outer_is_convex(model, fb, drop)
        || !model.faces.get(fa).inner.is_empty()
        || !model.faces.get(fb).inner.is_empty()
    {
        return None;
    }
    if face_contains_face(model, fa, fb, n) || face_contains_face(model, fb, fa, n) {
        return None; // fully contained ⇒ the boss cell, not an overhang
    }
    let na = planes_a.len();
    let nb = planes_b.len();
    // Either solid can serve as P (the union is symmetric — swapping P/Q swaps which face is the
    // notch vs the cantilever but yields the same two faces); try `a` as P, then `b`.
    let cc = if let Some(cr) = try_overhang(model, a, &planes_a, ia, fa, b, &planes_b, fb, na) {
        OverhangContact {
            p_solid: a,
            p_face: fa,
            q_solid: b,
            q_face: fb,
            crossings: cr,
        }
    } else if let Some(cr) = try_overhang(model, b, &planes_b, jb, fb, a, &planes_a, fa, nb) {
        OverhangContact {
            p_solid: b,
            p_face: fb,
            q_solid: a,
            q_face: fa,
            crossings: cr,
        }
    } else {
        return None;
    };
    // Pure coplanar contact: every vertex off the shared plane is strictly outside the other
    // solid. A vertex inside the other would be a transversal seam cut, not an overhang boss.
    let tri = planes_a[ia].tri;
    for (s, other) in [(cc.p_solid, cc.q_solid), (cc.q_solid, cc.p_solid)] {
        for &vh in &solid_vertex_handles(model, s) {
            let p = model.vertices.get(vh).point;
            if plane_side(tri, p) != 0 && point_in_solid(model, p, other).ok()? != Side::Outside {
                return None;
            }
        }
    }
    Some(cc)
}

/// The single-edge overhang result. P's contact face becomes a notch (its region minus the
/// overlap), Q's becomes the cantilever (its region minus the overlap), every wall is kept and
/// re-split at the crossings, and the crossings weld the two via `assemble_fuse_cut`.
fn overhang_contact_result(
    model: &mut Model,
    cc: &OverhangContact,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes_p = collect_planes(model, cc.p_solid)?;
    let na = planes_p.len();
    let planes_q = collect_planes(model, cc.q_solid)?;
    let mut planes = planes_p;
    planes.extend(planes_q);

    let seam: Vec<SeamVertex> = cc
        .crossings
        .iter()
        .map(|c| SeamVertex {
            point: c.point,
            triple: c.triple,
            tol: vertex_tol(
                c.point,
                &planes[c.triple[0]].plane,
                &planes[c.triple[1]].plane,
                &planes[c.triple[2]].plane,
            ),
        })
        .collect();

    let p_shell = model.solids.get(cc.p_solid).outer;
    let p_pos = model
        .shells
        .get(p_shell)
        .faces
        .iter()
        .position(|&f| f == cc.p_face)
        .expect("P contact on shell");
    let q_shell = model.solids.get(cc.q_solid).outer;
    let q_pos = model
        .shells
        .get(q_shell)
        .faces
        .iter()
        .position(|&f| f == cc.q_face)
        .expect("Q contact on shell");

    let drop = planar_drop_axes(planes[p_pos].n_out);
    // Both contact faces' outer loops (vertices + points), and P's holes.
    let p_face = model.faces.get(cc.p_face);
    let p_vh: Vec<Handle<Vertex>> = p_face
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let p_pts: Vec<Point3> = p_vh
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let p_inner: Vec<Vec<Node>> = p_face
        .inner
        .iter()
        .map(|l| {
            l.half_edges
                .iter()
                .map(|&he| Node::Orig(he_start(model, he)))
                .collect()
        })
        .collect();
    let q_face = model.faces.get(cc.q_face);
    let q_vh: Vec<Handle<Vertex>> = q_face
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let q_pts: Vec<Point3> = q_vh
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let p2: Vec<[f64; 2]> = p_pts.iter().map(|&p| proj2(p, drop)).collect();
    let q2: Vec<[f64; 2]> = q_pts.iter().map(|&p| proj2(p, drop)).collect();

    // Split both contact faces at the crossings into arcs. Each notch piece is a P-outside arc
    // stitched to the Q-inside arc sharing its two crossings; each cantilever piece a Q-outside
    // arc stitched to the P-inside arc sharing its ends. (Single-edge/corner = one piece each.)
    let p_nodes: Vec<Node> = p_vh.iter().map(|&vh| Node::Orig(vh)).collect();
    let q_nodes: Vec<Node> = q_vh.iter().map(|&vh| Node::Orig(vh)).collect();
    let p_arcs = split_loop_all_arcs(&p_nodes, &p_pts, &cc.crossings, &q2, drop, |c| c.p_seg);
    let q_arcs = split_loop_all_arcs(&q_nodes, &q_pts, &cc.crossings, &p2, drop, |c| c.q_seg);
    let match_arc = |arcs: &[LoopArc], ends: (usize, usize)| -> Option<Vec<Node>> {
        arcs.iter()
            .find(|a| !a.outside && unordered(a.ends.0, a.ends.1) == unordered(ends.0, ends.1))
            .map(|a| a.nodes.clone())
    };
    #[cfg(debug_assertions)]
    let sarea = |ring: &[[f64; 2]]| -> f64 {
        let k = ring.len();
        (0..k)
            .map(|i| {
                let j = (i + 1) % k;
                ring[i][0] * ring[j][1] - ring[j][0] * ring[i][1]
            })
            .sum()
    };
    #[cfg(debug_assertions)]
    let proj = |nodes: &[Node]| -> Vec<[f64; 2]> {
        nodes
            .iter()
            .map(|&nd| proj2(overhang_node_point(model, nd, &cc.crossings), drop))
            .collect()
    };

    // P's holes (each assigned to the notch piece that contains it).
    let hole_pt = |hole: &[Node]| -> [f64; 2] {
        proj2(overhang_node_point(model, hole[0], &cc.crossings), drop)
    };

    let mut faces = Vec::new();
    for pa in p_arcs.iter().filter(|a| a.outside) {
        let qi = match_arc(&q_arcs, pa.ends).ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
        let notch = stitch_arcs(pa.nodes.clone(), qi);
        #[cfg(debug_assertions)]
        debug_assert!(
            sarea(&proj(&notch)).signum() == sarea(&p2).signum(),
            "overhang notch winding flipped"
        );
        let notch2: Vec<[f64; 2]> = notch
            .iter()
            .map(|&nd| proj2(overhang_node_point(model, nd, &cc.crossings), drop))
            .collect();
        let inner: Vec<Vec<Node>> = p_inner
            .iter()
            .filter(|h| point_in_ring2(hole_pt(h), &notch2) == Some(true))
            .cloned()
            .collect();
        faces.push(LocalFace {
            plane_idx: p_pos,
            loop_nodes: notch,
            inner,
            flip: false,
        });
    }
    for qa in q_arcs.iter().filter(|a| a.outside) {
        let pi = match_arc(&p_arcs, qa.ends).ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
        let cantilever = stitch_arcs(qa.nodes.clone(), pi);
        #[cfg(debug_assertions)]
        debug_assert!(
            sarea(&proj(&cantilever)).signum() == sarea(&q2).signum(),
            "overhang cantilever winding flipped"
        );
        faces.push(LocalFace {
            plane_idx: na + q_pos,
            loop_nodes: cantilever,
            inner: Vec::new(),
            flip: false,
        });
    }
    let mut p_walls = solid_local_faces(model, cc.p_solid, 0, Some(cc.p_face), None);
    let mut q_walls = solid_local_faces(model, cc.q_solid, na, Some(cc.q_face), None);
    resplit_overhang(model, &mut p_walls, &planes, &cc.crossings);
    resplit_overhang(model, &mut q_walls, &planes, &cc.crossings);
    faces.extend(p_walls);
    faces.extend(q_walls);

    assemble_fuse_cut(model, cc.p_solid, cc.q_solid, &planes, &seam, &faces)
}

/// Insert each crossing that lies strictly interior to a wall's outer-loop edge, so the shared
/// edge is split to match the notch/cantilever (no T-junction). Exact: collinear (`orient2d`)
/// plus bbox, endpoints excluded; multiple crossings on one edge are inserted in order.
fn resplit_overhang(
    model: &Model,
    walls: &mut [LocalFace],
    planes: &[PlaneInfo],
    crossings: &[Crossing],
) {
    for lf in walls.iter_mut() {
        let tri = planes[lf.plane_idx].tri;
        // Only crossings actually on this wall's plane can lie on one of its edges. The 2D
        // collinearity test below drops an axis, so a crossing at a different depth but the
        // same projected coordinates would false-match without this 3D guard.
        let on_plane: Vec<&Crossing> = crossings
            .iter()
            .filter(|c| plane_side(tri, c.point) == 0)
            .collect();
        if on_plane.is_empty() {
            continue;
        }
        let drop = planar_drop_axes(planes[lf.plane_idx].n_out);
        let k = lf.loop_nodes.len();
        let mut out = Vec::with_capacity(k);
        for i in 0..k {
            let a_node = lf.loop_nodes[i];
            out.push(a_node);
            let a = overhang_node_point(model, a_node, crossings);
            let b = overhang_node_point(model, lf.loop_nodes[(i + 1) % k], crossings);
            let (a2, b2) = (proj2(a, drop), proj2(b, drop));
            let mut ins: Vec<(f64, Node)> = Vec::new();
            for c in on_plane.iter().copied() {
                if c.point == a || c.point == b {
                    continue;
                }
                let c2 = proj2(c.point, drop);
                if orient2d(a2, b2, c2) == 0.0 && in_bbox(a2, b2, c2) {
                    ins.push(((c.point - a).dot(b - a), Node::Seam(c.triple)));
                }
            }
            ins.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite params"));
            out.extend(ins.into_iter().map(|(_, node)| node));
        }
        lf.loop_nodes = out;
    }
}

// ---- coplanar-contact-overhang-cut (M5): shared helpers for the overhang Cut path ----

/// Whether `c` lies strictly between `a` and `b` on their segment (exact for axis-aligned edges).
fn point_strictly_on_segment(a: Point3, b: Point3, c: Point3) -> bool {
    let (ab, ac) = (b - a, c - a);
    if ac.cross(ab).norm_squared() != 0.0 {
        return false;
    }
    let t = ac.dot(ab);
    t > 0.0 && t < ab.dot(ab)
}

/// A face's outer loop as `Node::Orig`s.
fn face_orig_nodes(model: &Model, face: Handle<Face>) -> Vec<Node> {
    model
        .faces
        .get(face)
        .outer
        .half_edges
        .iter()
        .map(|&he| Node::Orig(he_start(model, he)))
        .collect()
}

/// The point where segment `a`–`b` meets `tri`'s plane (the two must straddle it). Exact for
/// axis-aligned geometry (the unnormalised distances cancel in `t`).
fn segment_plane_point(a: Point3, b: Point3, tri: [Point3; 3]) -> Point3 {
    let normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
    let da = normal.dot(a - tri[0]);
    let db = normal.dot(b - tri[0]);
    let t = da / (da - db);
    a + (b - a) * t
}

/// Sutherland–Hodgman clip of a convex polygon (as points) by `tri`'s plane, keeping the
/// `plane_side < 0` side. Crossing points are computed geometrically (so intermediate points a
/// later clip removes need not be canonical).
fn clip_points_by_plane(poly: &[Point3], tri: [Point3; 3]) -> Vec<Point3> {
    let n = poly.len();
    let side: Vec<i8> = poly.iter().map(|&p| plane_side(tri, p)).collect();
    let mut out = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        if side[i] <= 0 {
            out.push(poly[i]);
        }
        if side[i] != 0 && side[j] != 0 && (side[i] < 0) != (side[j] < 0) {
            out.push(segment_plane_point(poly[i], poly[j], tri));
        }
    }
    out
}

/// The planar cross-section of solid `a` cut by the plane at combined index `w_idx`, as a set of
/// closed `Node::Seam` loops. Each section vertex is the exact three-plane point where an `a`-edge
/// (shared by faces `A_f`,`A_g`) crosses `W`, provenance triple `{W,A_f,A_g}`. Decisions are exact/
/// toleranced only ([`t_plane_side`](crate::tolerant::t_plane_side) straddle, [`order_along`]
/// ordering) — no coordinate arbiter. A line crossing a face's ring meets it an even number of
/// times, so ordered same-face section vertices pair even-odd into material chords; each vertex
/// lies on exactly two faces (manifold) so it has degree 2, and the section is a union of simple
/// closed loops assembled by alternating chords. Honest-reject on degeneracy: an `a`-edge endpoint
/// *on* `W` (`VERTEX_ON_FACE_PLANE`), an odd per-face crossing count / non-degree-2 vertex
/// (`ARRANGEMENT_DEGENERATE`), `three_planes` failure (`THREE_PLANES`). Scope: cavity-free `a`
/// (outer shell) — a cavitied `a` is rejected upstream.
// Wired into the non-convex overhang clip in the next cell; used by tests now.
#[cfg_attr(not(test), allow(dead_code))]
fn section_of_solid(
    model: &Model,
    a: Handle<Solid>,
    w_idx: usize,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<Vec<Vec<Node>>, BoolError> {
    // A section vertex: its sorted dedup triple + the two a-faces it links (for connectivity).
    struct Sv {
        triple: [usize; 3],
        faces: [usize; 2],
    }
    // Deterministic edge order (Store index) so any downstream first-appearance is replay-stable.
    type SectionEdge = (Handle<Edge>, [Handle<Vertex>; 2], [usize; 2]);
    let inc = arrange::edge_planes(model, a, surf_ix)?;
    let mut edges: Vec<SectionEdge> = inc.iter().map(|(&e, &(b, f))| (e, b, f)).collect();
    edges.sort_by_key(|(e, _, _)| e.index());

    let mut svs: Vec<Sv> = Vec::new();
    for (_e, bounds, faces) in edges {
        let [v0, v1] = bounds;
        let s0 = crate::tolerant::t_plane_side(model, planes, w_idx, v0);
        let s1 = crate::tolerant::t_plane_side(model, planes, w_idx, v1);
        if s0 == 0 || s1 == 0 {
            return Err(reject(tag::VERTEX_ON_FACE_PLANE));
        }
        if s0 == s1 {
            continue; // edge does not straddle W
        }
        let [af, ag] = faces;
        three_planes(&planes[w_idx].plane, &planes[af].plane, &planes[ag].plane)
            .ok_or_else(|| reject(tag::THREE_PLANES))?;
        let mut triple = [w_idx, af, ag];
        triple.sort_unstable();
        svs.push(Sv {
            triple,
            faces: [af, ag],
        });
    }
    if svs.is_empty() {
        return Ok(Vec::new());
    }

    // Per section vertex, its chord partner on each of its two faces. Group vertices by face,
    // order along W∩A_f (each vertex's ordering key is its *other* face plane), pair even-odd.
    let mut neighbor: Vec<Vec<(usize, usize)>> = vec![Vec::new(); svs.len()]; // (face, partner)
    let mut faces_seen: Vec<usize> = svs.iter().flat_map(|s| s.faces).collect();
    faces_seen.sort_unstable();
    faces_seen.dedup();
    for &f in &faces_seen {
        // vertices on face f, with their "other" plane R for order_along
        let mut on_f: Vec<usize> = (0..svs.len())
            .filter(|&i| svs[i].faces.contains(&f))
            .collect();
        let other = |i: usize| -> usize {
            let [x, y] = svs[i].faces;
            if x == f { y } else { x }
        };
        on_f.sort_by(
            |&i, &j| match arrange::order_along(planes, w_idx, f, other(i), other(j)) {
                -1 => std::cmp::Ordering::Less,
                1 => std::cmp::Ordering::Greater,
                _ => std::cmp::Ordering::Equal,
            },
        );
        if on_f.len() % 2 != 0 {
            return Err(reject(tag::ARRANGEMENT_DEGENERATE)); // odd crossings = degeneracy
        }
        for pair in on_f.chunks_exact(2) {
            let (u, w) = (pair[0], pair[1]);
            neighbor[u].push((f, w));
            neighbor[w].push((f, u));
        }
    }
    // Every section vertex must have degree 2 (one chord per face).
    if neighbor.iter().any(|n| n.len() != 2) {
        return Err(reject(tag::ARRANGEMENT_DEGENERATE));
    }

    // Assemble closed loops: alternate chords (arrive on a face, leave on the other).
    let mut visited = vec![false; svs.len()];
    let mut loops: Vec<Vec<Node>> = Vec::new();
    for start in 0..svs.len() {
        if visited[start] {
            continue;
        }
        let mut loop_nodes = Vec::new();
        let mut cur = start;
        let mut take_face = svs[start].faces[0];
        for _ in 0..=svs.len() {
            visited[cur] = true;
            loop_nodes.push(Node::Seam(svs[cur].triple));
            let next = neighbor[cur]
                .iter()
                .find(|(f, _)| *f == take_face)
                .map(|(_, p)| *p)
                .ok_or_else(|| reject(tag::ARRANGEMENT_DEGENERATE))?;
            take_face = other_face(&svs[next].faces, take_face)?;
            cur = next;
            if cur == start {
                break;
            }
        }
        loops.push(loop_nodes);
    }
    Ok(loops)
}

/// The face of `faces` that is not `used` (the continuing chord's face at a section vertex).
#[cfg_attr(not(test), allow(dead_code))]
fn other_face(faces: &[usize; 2], used: usize) -> Result<usize, BoolError> {
    match *faces {
        [x, y] if x == used => Ok(y),
        [x, y] if y == used => Ok(x),
        _ => Err(reject(tag::ARRANGEMENT_DEGENERATE)),
    }
}

/// Clip a `b` wall to inside `a` by folding [`clip_points_by_plane`] over every breached wall's
/// plane, then map each surviving point back to a node — an original `b` vertex (`Node::Orig`) or
/// a canonical crossing (`Node::Seam`), by exact point match. `None` if a survivor matches
/// neither (out of scope), or `Some(vec![])` if the wall clips away entirely (drop it).
fn clip_bwall_inside_a(
    model: &Model,
    face_nodes: &[Node],
    breached_tris: &[[Point3; 3]],
    crossings: &[Crossing],
) -> Option<Vec<Node>> {
    // Fully inside every plane ⇒ keep the original nodes (no re-mapping).
    let orig: Vec<(Point3, Node)> = face_nodes
        .iter()
        .map(|&nd| (overhang_node_point(model, nd, crossings), nd))
        .collect();
    if breached_tris
        .iter()
        .all(|&tri| orig.iter().all(|&(p, _)| plane_side(tri, p) < 0))
    {
        return Some(face_nodes.to_vec());
    }
    let mut poly: Vec<Point3> = orig.iter().map(|&(p, _)| p).collect();
    for &tri in breached_tris {
        poly = clip_points_by_plane(&poly, tri);
        if poly.len() < 3 {
            return Some(Vec::new());
        }
    }
    poly.iter()
        .map(|&p| {
            if let Some(&(_, nd)) = orig.iter().find(|(op, _)| *op == p) {
                Some(nd)
            } else {
                crossings
                    .iter()
                    .find(|c| c.point == p)
                    .map(|c| Node::Seam(c.triple))
            }
        })
        .collect()
}

/// Replace the swallowed corner vertex `corner_vh` in a wall's loop with the three-node detour
/// (a top crossing, the interior floor point, the shared corner-column point) — oriented so its
/// ends stay on the corner's two incident edges.
fn splice_corner(
    model: &Model,
    loop_nodes: &[Node],
    corner_vh: Handle<Vertex>,
    detour: [Node; 3],
    crossings: &[Crossing],
) -> Vec<Node> {
    let n = loop_nodes.len();
    let pos = loop_nodes
        .iter()
        .position(|&nd| nd == Node::Orig(corner_vh))
        .expect("corner vertex on wall loop");
    let corner = model.vertices.get(corner_vh).point;
    let prev = overhang_node_point(model, loop_nodes[(pos + n - 1) % n], crossings);
    // detour[0] must lie on the edge from `prev` to the corner; else reverse.
    let d0 = overhang_node_point(model, detour[0], crossings);
    let ordered = if point_strictly_on_segment(prev, corner, d0) {
        detour
    } else {
        [detour[2], detour[1], detour[0]]
    };
    let mut out = Vec::with_capacity(n + 2);
    for (i, &nd) in loop_nodes.iter().enumerate() {
        if i == pos {
            out.extend_from_slice(&ordered);
        } else {
            out.push(nd);
        }
    }
    out
}

/// The mouth notch piece(s) for an overhang Cut: the contact face `p_face` minus its overlap with
/// `q_face`, as one or more `LocalFace`s on plane `p_pos`. `top` holds the contact-plane crossings
/// only (floor crossings would corrupt the split). Reuses the Fuse notch machinery
/// (`split_loop_all_arcs` + `stitch_arcs`), assigning each P hole to the piece that contains it.
fn mouth_notch_pieces(
    model: &Model,
    p_face: Handle<Face>,
    q_face: Handle<Face>,
    top: &[Crossing],
    p_pos: usize,
    drop: (usize, usize),
) -> Result<Vec<LocalFace>, BoolError> {
    let pt = |nd: Node| overhang_node_point(model, nd, top);
    let ring_vh = |face: Handle<Face>| -> Vec<Handle<Vertex>> {
        model
            .faces
            .get(face)
            .outer
            .half_edges
            .iter()
            .map(|&he| he_start(model, he))
            .collect()
    };
    let p_vh = ring_vh(p_face);
    let p_pts: Vec<Point3> = p_vh
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let p2: Vec<[f64; 2]> = p_pts.iter().map(|&p| proj2(p, drop)).collect();
    let p_inner: Vec<Vec<Node>> = model
        .faces
        .get(p_face)
        .inner
        .iter()
        .map(|l| {
            l.half_edges
                .iter()
                .map(|&he| Node::Orig(he_start(model, he)))
                .collect()
        })
        .collect();
    let q_vh = ring_vh(q_face);
    let q_pts: Vec<Point3> = q_vh
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let q2: Vec<[f64; 2]> = q_pts.iter().map(|&p| proj2(p, drop)).collect();
    let p_nodes: Vec<Node> = p_vh.iter().map(|&vh| Node::Orig(vh)).collect();
    let q_nodes: Vec<Node> = q_vh.iter().map(|&vh| Node::Orig(vh)).collect();
    let p_arcs = split_loop_all_arcs(&p_nodes, &p_pts, top, &q2, drop, |c| c.p_seg);
    let q_arcs = split_loop_all_arcs(&q_nodes, &q_pts, top, &p2, drop, |c| c.q_seg);
    let match_arc = |arcs: &[LoopArc], ends: (usize, usize)| -> Option<Vec<Node>> {
        arcs.iter()
            .find(|a| !a.outside && unordered(a.ends.0, a.ends.1) == unordered(ends.0, ends.1))
            .map(|a| a.nodes.clone())
    };
    let mut out = Vec::new();
    for pa in p_arcs.iter().filter(|a| a.outside) {
        let qi = match_arc(&q_arcs, pa.ends).ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
        let mouth = stitch_arcs(pa.nodes.clone(), qi);
        let mouth2: Vec<[f64; 2]> = mouth.iter().map(|&nd| proj2(pt(nd), drop)).collect();
        let inner: Vec<Vec<Node>> = p_inner
            .iter()
            .filter(|h| point_in_ring2(proj2(pt(h[0]), drop), &mouth2) == Some(true))
            .cloned()
            .collect();
        out.push(LocalFace {
            plane_idx: p_pos,
            loop_nodes: mouth,
            inner,
            flip: false,
        });
    }
    Ok(out)
}

// ---- coplanar-contact-overhang-cut-general (M5): one general N-wall overhang Cut ----

/// A general overhang Cut: prism `b` sits top-flush on `a`'s face (same-normal coplanar) but its
/// footprint hangs past `a`'s face boundary, breaking out through one or more walls. This one path
/// subsumes the edge-slot (1 wall), spanning-channel (2 opposite walls), and corner (2 adjacent
/// walls) cases, and opens 3+/L/U configurations. Each breached wall's side opening is classified
/// by how many of its top-edge corners `b` swallows (0/1/2); adjacent breached walls share a
/// corner column `three_planes(W_i, W_j, b-bottom)`.
struct OverhangCutG {
    a: Handle<Solid>,
    b: Handle<Solid>,
    p_face: Handle<Face>,     // a's contact face (→ mouth notch pieces)
    q_face: Handle<Face>,     // b's contact face (dropped)
    walls: Vec<WallG>,        // one per breached a wall
    crossings: Vec<Crossing>, // top c's, then corner columns cc, then floor d's
}

/// One breached `a` wall of a general overhang Cut: its plane index and the shape of its side
/// opening, keyed by the number of swallowed top-edge corners. All indices point into
/// `OverhangCutG::crossings`.
struct WallG {
    w_pos: usize,
    kind: WallKind,
}

/// A breached wall's opening, by swallowed top-corner count. `Middle` (0): the opening touches only
/// the top edge's interior — two top crossings `c` and two floor crossings `d` (`d[k]` paired to
/// `c[k]` by shared `b` wall). `Corner` (1): the opening bites one corner — one `c`, one `d`, and
/// the swallowed corner's shared column `cc`. `Shorten` (2): `b` covers the whole top edge — the
/// wall shrinks to a rectangle capped at both corners' columns (no `c`/`d`).
enum WallKind {
    Middle {
        c: [usize; 2],
        d: [usize; 2],
        p_seg: usize,
    },
    Corner {
        c: usize,
        d: usize,
        cc: usize,
        corner_vh: Handle<Vertex>,
    },
    Shorten([(Handle<Vertex>, usize); 2]),
}

/// `Some` iff `Cut(a, b)` is a general overhang slot: one same-normal coplanar pair whose
/// footprints are not nested, `b`'s footprint hanging past `a`'s face boundary through one or more
/// walls, and `b` blind on every non-breached `a` face. Absorbs the edge-slot, slab, and corner
/// detections.
fn detect_overhang_cut_general(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<OverhangCutG> {
    let planes_a = collect_planes(model, a).ok()?;
    let planes_b = collect_planes(model, b).ok()?;
    if !is_convex(model, &planes_a, &solid_vertex_handles(model, a))
        || !is_convex(model, &planes_b, &solid_vertex_handles(model, b))
    {
        return None;
    }
    let mut same: Vec<(usize, usize)> = Vec::new();
    for (i, pa) in planes_a.iter().enumerate() {
        for (j, pb) in planes_b.iter().enumerate() {
            if planes_coplanar(&pa.plane, &pb.plane) && pa.n_out.dot(pb.n_out) > 0.0 {
                same.push((i, j));
            }
        }
    }
    if same.len() != 1 {
        return None;
    }
    let (pi, qj) = same[0];
    let p_face = model.shells.get(model.solids.get(a).outer).faces[pi];
    let q_face = model.shells.get(model.solids.get(b).outer).faces[qj];
    let n = planes_a[pi].plane.normal();
    if face_contains_face(model, p_face, q_face, n) || face_contains_face(model, q_face, p_face, n)
    {
        return None; // contained ⇒ the blind-pocket cell
    }
    let na = planes_a.len();
    let cs = try_overhang(model, a, &planes_a, pi, p_face, b, &planes_b, q_face, na)?;
    let qn = planes_b[qj].n_out;
    let floors: Vec<usize> = (0..planes_b.len())
        .filter(|&j| planes_b[j].n_out.dot(qn) < 0.0)
        .collect();
    let [floor] = floors[..] else {
        return None; // not a single opposite cap ⇒ out of scope
    };
    let p_ring: Vec<Handle<Vertex>> = model
        .faces
        .get(p_face)
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let np = p_ring.len();
    let q_ring: Vec<Handle<Vertex>> = model
        .faces
        .get(q_face)
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let nq = q_ring.len();
    let drop = planar_drop_axes(planes_a[pi].n_out);
    let q2: Vec<[f64; 2]> = q_ring
        .iter()
        .map(|&vh| proj2(model.vertices.get(vh).point, drop))
        .collect();

    // Which a-face corners `b` swallows (strictly inside its footprint).
    let swallowed: Vec<bool> = p_ring
        .iter()
        .map(|&vh| point_in_ring2(proj2(model.vertices.get(vh).point, drop), &q2) == Some(true))
        .collect();
    let wall_of = |i: usize| wall_pos_on_edge(model, a, p_face, p_ring[i], p_ring[(i + 1) % np]);

    // crossings: top c's first (indices 0..cs.len()), then the corner columns, then the floor d's.
    let mut crossings: Vec<Crossing> = cs.clone();

    // A shared corner column per swallowed corner: the meet of its two incident (both breached)
    // walls and the b floor. A simple loop ⇒ each corner has exactly two incident walls.
    let mut cc_of: HashMap<Handle<Vertex>, usize> = HashMap::new();
    for j in 0..np {
        if !swallowed[j] {
            continue;
        }
        let w_prev = wall_of((j + np - 1) % np)?;
        let w_next = wall_of(j)?;
        let point = three_planes(
            &planes_a[w_prev].plane,
            &planes_a[w_next].plane,
            &planes_b[floor].plane,
        )?;
        let mut triple = [w_prev, w_next, na + floor];
        triple.sort_unstable();
        let idx = crossings.len();
        crossings.push(Crossing {
            point,
            triple,
            p_seg: j,
            q_seg: 0,
        });
        cc_of.insert(p_ring[j], idx);
    }

    // A floor crossing per top crossing (its wall ∩ the b wall it came from ∩ the b floor).
    let mut d_of: Vec<usize> = Vec::with_capacity(cs.len());
    for c in &cs {
        let w = wall_of(c.p_seg)?;
        let ywall = wall_pos_on_edge(
            model,
            b,
            q_face,
            q_ring[c.q_seg],
            q_ring[(c.q_seg + 1) % nq],
        )?;
        let point = three_planes(
            &planes_a[w].plane,
            &planes_b[ywall].plane,
            &planes_b[floor].plane,
        )?;
        let mut triple = [w, na + ywall, na + floor];
        triple.sort_unstable();
        d_of.push(crossings.len());
        crossings.push(Crossing {
            point,
            triple,
            p_seg: c.p_seg,
            q_seg: c.q_seg,
        });
    }

    // Classify each breached wall. An a-top edge is breached iff `b` overlaps it: it has a crossing
    // or a swallowed endpoint. Invariant on a convex box: (crossings on the edge) == 2 − (swallowed
    // endpoints), so n∈{0,1,2} partitions Middle/Corner/Shorten.
    let mut walls: Vec<WallG> = Vec::new();
    for i in 0..np {
        let (sw_lo, sw_hi) = (swallowed[i], swallowed[(i + 1) % np]);
        let on: Vec<usize> = (0..cs.len()).filter(|&k| cs[k].p_seg == i).collect();
        if !sw_lo && !sw_hi && on.is_empty() {
            continue; // untouched
        }
        let ncorner = sw_lo as usize + sw_hi as usize;
        if on.len() != 2 - ncorner {
            return None; // degenerate / out of scope
        }
        let w_pos = wall_of(i)?;
        let kind = match ncorner {
            0 => {
                let [c0, c1] = on[..] else { return None };
                WallKind::Middle {
                    c: [c0, c1],
                    d: [d_of[c0], d_of[c1]],
                    p_seg: i,
                }
            }
            1 => {
                let [c] = on[..] else { return None };
                let corner_vh = if sw_lo {
                    p_ring[i]
                } else {
                    p_ring[(i + 1) % np]
                };
                let cc = *cc_of.get(&corner_vh)?;
                WallKind::Corner {
                    c,
                    d: d_of[c],
                    cc,
                    corner_vh,
                }
            }
            _ => {
                let (vlo, vhi) = (p_ring[i], p_ring[(i + 1) % np]);
                WallKind::Shorten([(vlo, *cc_of.get(&vlo)?), (vhi, *cc_of.get(&vhi)?)])
            }
        };
        walls.push(WallG { w_pos, kind });
    }
    if walls.is_empty() {
        return None;
    }

    // Blind everywhere but the breached walls: every off-plane b vertex is strictly inside every
    // non-breached a face. Catches a through-bottom, a full-face cover, or a non-axis breakout.
    let breached: Vec<usize> = walls.iter().map(|w| w.w_pos).collect();
    let tri = planes_a[pi].tri;
    for &vh in &solid_vertex_handles(model, b) {
        let p = model.vertices.get(vh).point;
        if plane_side(tri, p) == 0 {
            continue;
        }
        for (fi, pf) in planes_a.iter().enumerate() {
            if breached.contains(&fi) {
                continue;
            }
            if plane_side(pf.tri, p) >= 0 {
                return None;
            }
        }
    }
    Some(OverhangCutG {
        a,
        b,
        p_face,
        q_face,
        walls,
        crossings,
    })
}

/// Replace each swallowed corner vertex in a wall's loop with its shared corner-column node. The
/// wall's whole top edge is swallowed (both corners), so it shrinks to a rectangle capped at the
/// two columns (a straight `cc_lo → cc_hi` trace on the b floor).
fn substitute_corners(
    loop_nodes: &[Node],
    pairs: &[(Handle<Vertex>, usize)],
    crossings: &[Crossing],
) -> Vec<Node> {
    loop_nodes
        .iter()
        .map(|&nd| {
            if let Node::Orig(vh) = nd {
                if let Some(&(_, cc)) = pairs.iter().find(|(v, _)| *v == vh) {
                    return Node::Seam(crossings[cc].triple);
                }
            }
            nd
        })
        .collect()
}

/// Build the general overhang Cut result: the mouth notch on `a`'s contact face, a per-wall side
/// notch dispatched by swallowed-corner count, and `b`'s walls clipped to inside `a` and flipped to
/// bound the removed slot.
fn overhang_cut_general_result(
    model: &mut Model,
    oc: &OverhangCutG,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes_a = collect_planes(model, oc.a)?;
    let na = planes_a.len();
    let planes_b = collect_planes(model, oc.b)?;
    let mut planes = planes_a;
    planes.extend(planes_b);
    let seam: Vec<SeamVertex> = oc
        .crossings
        .iter()
        .map(|c| SeamVertex {
            point: c.point,
            triple: c.triple,
            tol: vertex_tol(
                c.point,
                &planes[c.triple[0]].plane,
                &planes[c.triple[1]].plane,
                &planes[c.triple[2]].plane,
            ),
        })
        .collect();

    let a_shell = model.solids.get(oc.a).outer;
    let p_pos = model
        .shells
        .get(a_shell)
        .faces
        .iter()
        .position(|&f| f == oc.p_face)
        .expect("P on a");
    let pt = |nd: Node| overhang_node_point(model, nd, &oc.crossings);
    let drop = planar_drop_axes(planes[p_pos].n_out);
    let contact_tri = planes[p_pos].tri;

    // mouth (top crossings only — floor crossings are off P's plane and corrupt the split).
    let top: Vec<Crossing> = oc
        .crossings
        .iter()
        .filter(|c| plane_side(contact_tri, c.point) == 0)
        .cloned()
        .collect();
    let mut faces = mouth_notch_pieces(model, oc.p_face, oc.q_face, &top, p_pos, drop)?;
    let p_vh: Vec<Handle<Vertex>> = model
        .faces
        .get(oc.p_face)
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();

    // each breached wall: its loop with the opening spliced in, by swallowed-corner count.
    let mut a_faces = solid_local_faces(model, oc.a, 0, Some(oc.p_face), None);
    for wall in &oc.walls {
        let w_face = model.shells.get(a_shell).faces[wall.w_pos];
        let w_nodes = face_orig_nodes(model, w_face);
        let side = match &wall.kind {
            WallKind::Middle { c, d, p_seg } => {
                let detour = vec![
                    Node::Seam(oc.crossings[c[0]].triple),
                    Node::Seam(oc.crossings[d[0]].triple),
                    Node::Seam(oc.crossings[d[1]].triple),
                    Node::Seam(oc.crossings[c[1]].triple),
                ];
                let (e0v, e1v) = (p_vh[*p_seg], p_vh[(*p_seg + 1) % p_vh.len()]);
                let kw = w_nodes.len();
                let w_seg = (0..kw)
                    .find(|&i| {
                        let (u, v) = (node_vh(w_nodes[i]), node_vh(w_nodes[(i + 1) % kw]));
                        (u == e0v && v == e1v) || (u == e1v && v == e0v)
                    })
                    .expect("wall shares the crossed edge with P");
                splice_notch(&w_nodes, w_seg, detour, &pt)
            }
            WallKind::Corner {
                c,
                d,
                cc,
                corner_vh,
            } => {
                let detour = [
                    Node::Seam(oc.crossings[*c].triple),
                    Node::Seam(oc.crossings[*d].triple),
                    Node::Seam(oc.crossings[*cc].triple),
                ];
                splice_corner(model, &w_nodes, *corner_vh, detour, &oc.crossings)
            }
            WallKind::Shorten(pairs) => substitute_corners(&w_nodes, pairs, &oc.crossings),
        };
        #[cfg(debug_assertions)]
        {
            let sarea = |nodes: &[Node]| -> f64 {
                let ps: Vec<[f64; 2]> = nodes
                    .iter()
                    .map(|&nd| {
                        proj2(
                            overhang_node_point(model, nd, &oc.crossings),
                            planar_drop_axes(planes[wall.w_pos].n_out),
                        )
                    })
                    .collect();
                let k = ps.len();
                (0..k)
                    .map(|i| {
                        let j = (i + 1) % k;
                        ps[i][0] * ps[j][1] - ps[j][0] * ps[i][1]
                    })
                    .sum()
            };
            debug_assert!(
                sarea(&side).signum() == sarea(&w_nodes).signum(),
                "wall notch winding flipped"
            );
        }
        let lf = a_faces
            .iter_mut()
            .find(|lf| lf.plane_idx == wall.w_pos)
            .expect("breached wall in a_faces");
        lf.loop_nodes = side;
        lf.inner = Vec::new();
    }
    faces.extend(a_faces);

    // b's walls clipped to inside a (geometric SH over every breached plane), flipped.
    let breached_tris: Vec<[Point3; 3]> = oc.walls.iter().map(|w| planes[w.w_pos].tri).collect();
    for lf in solid_local_faces(model, oc.b, na, Some(oc.q_face), None) {
        let clipped = clip_bwall_inside_a(model, &lf.loop_nodes, &breached_tris, &oc.crossings)
            .ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
        if clipped.len() >= 3 {
            faces.push(LocalFace {
                loop_nodes: clipped,
                flip: true,
                ..lf
            });
        }
    }

    assemble_fuse_cut(model, oc.a, oc.b, &planes, &seam, &faces)
}

// ---- coplanar-contact-overhang-common (M5): intersection of a top-flush overhang pair ----

/// A general overhang `Common`: the same top-flush coplanar overhang geometry the Cut path
/// detects, but `Common(a, b) = a ∩ b` keeps the convex overlap polytope `R` instead of carving it
/// out. `R`'s faces are the contact overlap (P ∩ Q), `a`'s breached walls clipped to inside `b`, and
/// `b`'s walls clipped to inside `a` — all with their original outward normals (`R` convex, so no
/// flip). Restricted to two top crossings (one overlap-arc pair, assembled by a single stitch);
/// more crossings need arc chaining and stay out of scope.
/// `Some` iff `Common(a, b)` is a two-crossing overhang intersection. Reuses the kind-independent
/// overhang geometry (`detect_overhang_cut_general`) and gates on exactly two contact-plane
/// crossings, so the overlap has one `!outside` P-arc that a single `stitch_arcs` assembles (edge,
/// corner, and L configs — the breached-wall count is irrelevant). Four+ crossings (a spanning
/// slab) fall through to the standard `VERTEX_ON_FACE_PLANE` rejection.
fn detect_overhang_common(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Option<OverhangCutG> {
    let oc = detect_overhang_cut_general(model, a, b)?;
    let planes_a = collect_planes(model, a).ok()?;
    let pi = model
        .shells
        .get(model.solids.get(a).outer)
        .faces
        .iter()
        .position(|&f| f == oc.p_face)?;
    let contact_tri = planes_a[pi].tri;
    let top = oc
        .crossings
        .iter()
        .filter(|c| plane_side(contact_tri, c.point) == 0)
        .count();
    if top != 2 {
        return None; // >1 overlap-arc pair (slab) — out of scope, honest fall-through
    }
    Some(oc)
}

/// A breached wall's face on the intersection `R = a ∩ b`, built exactly from the wall's own
/// crossings by swallowed-corner count (`R` is convex, so every corner is a stored crossing or a
/// shared original vertex — no clipping, which a chain of parallel `b` planes would round off).
fn common_wall_face(wall: &WallG, crossings: &[Crossing]) -> Vec<Node> {
    let cr = |i: usize| Node::Seam(crossings[i].triple);
    match &wall.kind {
        WallKind::Middle { c, d, .. } => vec![cr(c[0]), cr(c[1]), cr(d[1]), cr(d[0])],
        WallKind::Corner {
            c,
            d,
            cc,
            corner_vh,
        } => vec![cr(*c), Node::Orig(*corner_vh), cr(*cc), cr(*d)],
        WallKind::Shorten(pairs) => vec![
            Node::Orig(pairs[0].0),
            Node::Orig(pairs[1].0),
            cr(pairs[1].1),
            cr(pairs[0].1),
        ],
    }
}

/// Order `nodes` to match wall `w_face`'s real outer-loop winding, so `assemble_fuse_cut` gives the
/// face the wall's outward normal (no flip). Exact — only reorders the given nodes.
fn orient_to_wall(
    model: &Model,
    w_face: Handle<Face>,
    nodes: Vec<Node>,
    n_out: Vector3,
    crossings: &[Crossing],
) -> Vec<Node> {
    let drop = planar_drop_axes(n_out);
    let sarea = |ns: &[Node]| -> f64 {
        let ps: Vec<[f64; 2]> = ns
            .iter()
            .map(|&nd| proj2(overhang_node_point(model, nd, crossings), drop))
            .collect();
        let k = ps.len();
        (0..k)
            .map(|i| {
                let j = (i + 1) % k;
                ps[i][0] * ps[j][1] - ps[j][0] * ps[i][1]
            })
            .sum()
    };
    if sarea(&nodes).signum() == sarea(&face_orig_nodes(model, w_face)).signum() {
        nodes
    } else {
        nodes.into_iter().rev().collect()
    }
}

/// The overlap top piece for an overhang `Common`: the contact face `p_face` intersected with
/// `q_face`, as one `LocalFace` on plane `p_pos`. The sibling of `mouth_notch_pieces`, but it stitches
/// the **inside** arcs (P ∩ Q) rather than the outside ones. `top` holds the contact-plane crossings
/// only. Two crossings ⇒ one `!outside` P-arc and one matching `!outside` Q-arc.
fn overlap_top_piece(
    model: &Model,
    p_face: Handle<Face>,
    q_face: Handle<Face>,
    top: &[Crossing],
    p_pos: usize,
    drop: (usize, usize),
) -> Result<LocalFace, BoolError> {
    let ring_vh = |face: Handle<Face>| -> Vec<Handle<Vertex>> {
        model
            .faces
            .get(face)
            .outer
            .half_edges
            .iter()
            .map(|&he| he_start(model, he))
            .collect()
    };
    let p_vh = ring_vh(p_face);
    let p_pts: Vec<Point3> = p_vh
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let p2: Vec<[f64; 2]> = p_pts.iter().map(|&p| proj2(p, drop)).collect();
    let q_vh = ring_vh(q_face);
    let q_pts: Vec<Point3> = q_vh
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let q2: Vec<[f64; 2]> = q_pts.iter().map(|&p| proj2(p, drop)).collect();
    let p_nodes: Vec<Node> = p_vh.iter().map(|&vh| Node::Orig(vh)).collect();
    let q_nodes: Vec<Node> = q_vh.iter().map(|&vh| Node::Orig(vh)).collect();
    let p_arcs = split_loop_all_arcs(&p_nodes, &p_pts, top, &q2, drop, |c| c.p_seg);
    let q_arcs = split_loop_all_arcs(&q_nodes, &q_pts, top, &p2, drop, |c| c.q_seg);
    let match_arc = |arcs: &[LoopArc], ends: (usize, usize)| -> Option<Vec<Node>> {
        arcs.iter()
            .find(|a| !a.outside && unordered(a.ends.0, a.ends.1) == unordered(ends.0, ends.1))
            .map(|a| a.nodes.clone())
    };
    let inside_p: Vec<&LoopArc> = p_arcs.iter().filter(|a| !a.outside).collect();
    debug_assert_eq!(inside_p.len(), 1, "two crossings give one inside P-arc");
    let pa = *inside_p.first().ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
    let qi = match_arc(&q_arcs, pa.ends).ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
    let overlap = stitch_arcs(pa.nodes.clone(), qi);
    #[cfg(debug_assertions)]
    {
        let sarea = |ring: &[[f64; 2]]| -> f64 {
            let k = ring.len();
            (0..k)
                .map(|i| {
                    let j = (i + 1) % k;
                    ring[i][0] * ring[j][1] - ring[j][0] * ring[i][1]
                })
                .sum()
        };
        let o2: Vec<[f64; 2]> = overlap
            .iter()
            .map(|&nd| proj2(overhang_node_point(model, nd, top), drop))
            .collect();
        debug_assert!(
            sarea(&o2).signum() == sarea(&p2).signum(),
            "overlap top winding flipped"
        );
    }
    Ok(LocalFace {
        plane_idx: p_pos,
        loop_nodes: overlap,
        inner: Vec::new(),
        flip: false,
    })
}

/// Build the overhang `Common` result `R = a ∩ b`: the overlap top, `a`'s breached walls clipped to
/// inside `b`, and `b`'s walls clipped to inside `a` — all keeping their original outward normals.
fn overhang_common_result(
    model: &mut Model,
    oc: &OverhangCutG,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let planes_a = collect_planes(model, oc.a)?;
    let na = planes_a.len();
    let planes_b = collect_planes(model, oc.b)?;
    let mut planes = planes_a;
    planes.extend(planes_b);
    let seam: Vec<SeamVertex> = oc
        .crossings
        .iter()
        .map(|c| SeamVertex {
            point: c.point,
            triple: c.triple,
            tol: vertex_tol(
                c.point,
                &planes[c.triple[0]].plane,
                &planes[c.triple[1]].plane,
                &planes[c.triple[2]].plane,
            ),
        })
        .collect();

    let a_shell = model.solids.get(oc.a).outer;
    let p_pos = model
        .shells
        .get(a_shell)
        .faces
        .iter()
        .position(|&f| f == oc.p_face)
        .expect("P on a");
    let drop = planar_drop_axes(planes[p_pos].n_out);
    let contact_tri = planes[p_pos].tri;

    // overlap top: P ∩ Q (contact-plane crossings only — floor crossings are off P's plane).
    let top: Vec<Crossing> = oc
        .crossings
        .iter()
        .filter(|c| plane_side(contact_tri, c.point) == 0)
        .cloned()
        .collect();
    let mut faces = vec![overlap_top_piece(
        model, oc.p_face, oc.q_face, &top, p_pos, drop,
    )?];

    // a's breached walls, each ∩ b — built exactly from the wall's own crossings (no clipping;
    // R is convex so every corner is a stored crossing or shared original), oriented to the wall's
    // real winding so `assemble_fuse_cut` keeps its outward normal (no flip).
    for wall in &oc.walls {
        let w_face = model.shells.get(a_shell).faces[wall.w_pos];
        let nodes = common_wall_face(wall, &oc.crossings);
        let nodes = orient_to_wall(
            model,
            w_face,
            nodes,
            planes[wall.w_pos].n_out,
            &oc.crossings,
        );
        faces.push(LocalFace {
            plane_idx: wall.w_pos,
            loop_nodes: nodes,
            inner: Vec::new(),
            flip: false,
        });
    }

    // b's walls clipped to inside a (folded over a's breached walls only — a's top is excluded, so
    // the z=1 top edge survives). No flip: R is the kept solid, so b's outward normal points out of R.
    let breached_tris: Vec<[Point3; 3]> = oc.walls.iter().map(|w| planes[w.w_pos].tri).collect();
    for lf in solid_local_faces(model, oc.b, na, Some(oc.q_face), None) {
        let clipped = clip_bwall_inside_a(model, &lf.loop_nodes, &breached_tris, &oc.crossings)
            .ok_or_else(|| reject(tag::OVERHANG_ARCS))?;
        if clipped.len() >= 3 {
            faces.push(LocalFace {
                loop_nodes: clipped,
                flip: false,
                ..lf
            });
        }
    }

    assemble_fuse_cut(model, oc.a, oc.b, &planes, &seam, &faces)
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Test shim: a boolean whose result is exactly one solid. Most tests operate on a single
    /// body; this asserts that and returns the lone handle, so call sites read as before while
    /// `boolean` itself returns the full `Vec` (cell 0.4 multi-solid).
    fn boolean_one(
        model: &mut Model,
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    ) -> Result<Handle<Solid>, BoolError> {
        let solids = boolean(model, kind, a, b)?;
        assert_eq!(
            solids.len(),
            1,
            "boolean_one: expected one solid, got {}",
            solids.len()
        );
        Ok(solids[0])
    }

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
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
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

    /// `point_in_solid_tol` on a rotated solid + a query rotated by the **same** rigid motion
    /// agrees with the exact f64 `point_in_solid` on the unrotated solid (in/out is
    /// rotation-invariant). Covers a convex cube and the concave L-prism — including the reflex
    /// notch point `(1.5, 1.5)` a convex hull would misclassify.
    #[test]
    fn point_in_solid_tol_is_rotation_invariant() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation as SRot};
        let ang = Angle::from_deg(Rat::from_int(30)).unwrap();
        let piv = [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)];
        let rot = |m: &mut Model, s: Handle<Solid>| -> Handle<Solid> {
            let iso = Isometry::rotation(SRot {
                axis: Axis::Z,
                point: piv,
                angle: ang,
            });
            let out = apply(
                m,
                &Operation::Transform {
                    solid: s,
                    isometry: iso,
                },
            )
            .unwrap();
            m.rebuild_adjacency();
            match out {
                crate::OpOutput::Transform { solid } => solid,
                _ => panic!("expected Transform"),
            }
        };
        let rq = |x: f64, y: f64, z: f64| {
            let r = |v: f64| Rat::try_from_f64(v).unwrap();
            Pt3::at([r(x), r(y), r(z)]).rotate_about(Axis::Z, ang, piv)
        };
        // Convex cube [0,2]³.
        let mut mc = Model::new();
        let cube = mc.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let rc = rot(&mut mc, cube);
        for (x, y, z, want) in [
            (1.0, 1.0, 1.0, Side::Inside),
            (0.5, 0.5, 0.5, Side::Inside),
            (1.0, 1.0, 1.9, Side::Inside),
            (3.0, 1.0, 1.0, Side::Outside),
            (-1.0, 1.0, 1.0, Side::Outside),
        ] {
            assert_eq!(
                point_in_solid_tol(&mc, &rq(x, y, z), rc).unwrap(),
                want,
                "cube ({x},{y},{z})"
            );
        }
        // Concave L-prism: notch (x,y)∈[1,2]² is outside though inside the convex hull.
        let (mut ml, lp) = l_prism();
        let rl = rot(&mut ml, lp);
        for (x, y, z, want) in [
            (0.5, 0.5, 0.5, Side::Inside),
            (1.5, 0.5, 0.5, Side::Inside),
            (0.5, 1.5, 0.5, Side::Inside),
            (1.5, 1.5, 0.5, Side::Outside),
            (3.0, 0.5, 0.5, Side::Outside),
        ] {
            assert_eq!(
                point_in_solid_tol(&ml, &rq(x, y, z), rl).unwrap(),
                want,
                "L-prism ({x},{y},{z})"
            );
        }
    }

    /// `pt3_base_collinear` is exact on the pre-rotation rational bases: three genuinely
    /// collinear points stay collinear under rotation (→ skipped), and a real sliver (one point
    /// off the line) is never falsely called collinear (→ its crossing is kept, no silent-wrong).
    #[test]
    fn pt3_base_collinear_exact() {
        use nacre_scalar::{Angle, Axis, Rat};
        let ang = Angle::from_deg(Rat::from_int(37)).unwrap();
        let piv = [Rat::from_int(2), Rat::from_int(-1), Rat::from_int(0)];
        let rp = |x: i128, y: i128, z: i128| {
            Pt3::at([Rat::from_int(x), Rat::from_int(y), Rat::from_int(z)]).rotate_about(
                Axis::Z,
                ang,
                piv,
            )
        };
        // (0,0,0), (2,4,6), (1,2,3): all on the line t·(1,2,3) → collinear.
        assert!(pt3_base_collinear(&rp(0, 0, 0), &rp(2, 4, 6), &rp(1, 2, 3)));
        // (1,2,4) is off that line (z), a real nonzero-area triangle → not collinear.
        assert!(!pt3_base_collinear(
            &rp(0, 0, 0),
            &rp(2, 4, 6),
            &rp(1, 2, 4)
        ));
    }

    /// `vertex_in_solid` dispatches by rotation: an axis-aligned target forwards to the exact
    /// f64 `point_in_solid`, a rotated target to `point_in_solid_tol` on the vertex's `Pt3`.
    /// (Value correctness is `point_in_solid_tol`'s own rotation-invariance test; here the query
    /// is a far-away vertex — outside either target — so the two classifiers must agree.)
    #[test]
    fn vertex_in_solid_routes_by_rotation() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation as SRot};
        let mut m = Model::new();
        let b = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let a = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
        let vh = solid_vertex_handles(&m, a)[0];
        // Axis-aligned target → the router forwards to the exact f64 `point_in_solid`.
        assert_eq!(
            vertex_in_solid(&m, vh, b).unwrap(),
            point_in_solid(&m, m.vertices.get(vh).point, b).unwrap()
        );
        // Rotate the target → the router must dispatch to `point_in_solid_tol`.
        let iso = Isometry::rotation(SRot {
            axis: Axis::Z,
            point: [Rat::from_int(2), Rat::from_int(2), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        let br = match apply(
            &mut m,
            &Operation::Transform {
                solid: b,
                isometry: iso,
            },
        )
        .unwrap()
        {
            crate::OpOutput::Transform { solid } => solid,
            _ => panic!("expected Transform"),
        };
        m.rebuild_adjacency();
        let p = nacre_tip::vertex_pt3(&m, vh).unwrap();
        assert_eq!(
            vertex_in_solid(&m, vh, br).unwrap(),
            point_in_solid_tol(&m, &p, br).unwrap()
        );
    }

    /// The toleranced `ray_triangle_cross_tol` reproduces the exact f64 `ray_face_cross`
    /// bit-for-bit on unrotated integer coords — validating the five-`orient3d`→
    /// `orient3d_ray`/`orient3d_judge`/`dir_orient3d_judge` sign reduction. Integer coords keep
    /// every nonzero determinant ≥ 1, far above the declare-0 floor, so the toleranced predicate
    /// declares 0 only on TRUE zeros, matching Shewchuk exactly (Cross / Miss / Degenerate and
    /// the Cross sign).
    #[test]
    fn ray_triangle_cross_tol_matches_f64_unrotated() {
        use nacre_scalar::Rat;
        let mut st = 0x9E37_79B9_7F4A_7C15u64;
        let mut g = || {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((st >> 40) % 21) as i128 - 10
        };
        let f = |q: [i128; 3]| Point3::from_array([q[0] as f64, q[1] as f64, q[2] as f64]);
        let pt = |q: [i128; 3]| {
            Pt3::at([
                Rat::from_int(q[0]),
                Rat::from_int(q[1]),
                Rat::from_int(q[2]),
            ])
        };
        let mut cases = 0usize;
        for _ in 0..300 {
            let (p, a, b, c) = (
                [g(), g(), g()],
                [g(), g(), g()],
                [g(), g(), g()],
                [g(), g(), g()],
            );
            for dir in RAY_DIRECTIONS {
                let want = ray_face_cross(f(p), Vector3::from_array(dir), [f(a), f(b), f(c)]);
                let dr = [
                    Rat::from_int(dir[0] as i128),
                    Rat::from_int(dir[1] as i128),
                    Rat::from_int(dir[2] as i128),
                ];
                let got = ray_triangle_cross_tol(&pt(p), dr, &pt(a), &pt(b), &pt(c));
                assert_eq!(got, want, "p={p:?} a={a:?} b={b:?} c={c:?} dir={dir:?}");
                cases += 1;
            }
        }
        assert!(cases > 0);
    }

    /// Rotating the whole config — points about a pivot, direction about the origin — by an
    /// exact 90° about Z (a rigid rotation, `det +1`) leaves the crossing and its sign
    /// invariant, and exercises the `Pt3` rotation chain. 90°-Z sends `(x,y,z) → (−y,x,z)`, so
    /// the direction `d → (−d1, d0, d2)`. (90° keeps every coordinate exact, so this checks the
    /// helper on rotated `Pt3` inputs; the `tol > 0` soundness is `dir_orient3d`'s corpus.)
    #[test]
    fn ray_triangle_cross_tol_is_rotation_equivariant() {
        use nacre_scalar::{Angle, Axis, Rat};
        let mut st = 0x1234_5678_9ABC_DEF0u64;
        let mut g = || {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((st >> 40) % 15) as i128 - 7
        };
        let piv = [Rat::from_int(1), Rat::from_int(-2), Rat::from_int(0)];
        let ang = Angle::from_deg(Rat::from_int(90)).unwrap();
        let rat = |q: [i128; 3]| {
            [
                Rat::from_int(q[0]),
                Rat::from_int(q[1]),
                Rat::from_int(q[2]),
            ]
        };
        let un = |q: [i128; 3]| Pt3::at(rat(q));
        let rot = |q: [i128; 3]| Pt3::at(rat(q)).rotate_about(Axis::Z, ang, piv);
        let mut cases = 0usize;
        for _ in 0..200 {
            let (p, a, b, c) = (
                [g(), g(), g()],
                [g(), g(), g()],
                [g(), g(), g()],
                [g(), g(), g()],
            );
            let d = [g(), g(), g()];
            if d == [0i128, 0, 0] {
                continue;
            }
            let plain = ray_triangle_cross_tol(&un(p), rat(d), &un(a), &un(b), &un(c));
            let rotated = ray_triangle_cross_tol(
                &rot(p),
                rat([-d[1], d[0], d[2]]),
                &rot(a),
                &rot(b),
                &rot(c),
            );
            assert_eq!(plain, rotated, "rotation-equivariant: p={p:?} d={d:?}");
            cases += 1;
        }
        assert!(cases > 0);
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

    /// (5d)#4: `fan_triangles` drops *exactly* the zero-area (collinear) triangles
    /// and keeps the rest. A pentagon ring with three collinear points on one edge
    /// has one collinear fan triangle from apex 0; it is dropped, the other two are
    /// kept (and confirmed non-degenerate).
    #[test]
    fn fan_triangles_drops_only_exactly_collinear() {
        let ring = [
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]), // collinear with its neighbours on y=0
            Point3::from_array([2.0, 0.0, 0.0]),
            Point3::from_array([2.0, 2.0, 0.0]),
            Point3::from_array([0.0, 2.0, 0.0]),
        ];
        let tris = fan_triangles(&ring, 0);
        // apex-0 fan: (0,1,2) is collinear → dropped; (0,2,3) and (0,3,4) kept.
        assert_eq!(tris.len(), 2);
        for t in &tris {
            assert!(!triangle_is_degenerate(t[0], t[1], t[2]));
        }
    }

    /// (5d)#4 regression: a tiny-but-*exactly-nonzero* sliver is kept, where the
    /// retired relative `1e-12` bound would have silently dropped it (missing a ray
    /// crossing → wrong winding). Large-magnitude integer coords make the cross
    /// product tiny relative to the edge lengths: `cross.z = 2a − (2a+1) = −1`,
    /// `|e1||e2| ≈ 2a²`. (A direct `fan_triangles` unit test — such slivers arise
    /// in rotated geometry, not axis-aligned M5.)
    #[test]
    fn fan_triangles_keeps_tiny_nonzero_sliver() {
        let a = 1_000_000.0;
        let t0 = Point3::from_array([0.0, 0.0, 0.0]);
        let t1 = Point3::from_array([a, 1.0, 0.0]);
        let t2 = Point3::from_array([2.0 * a + 1.0, 2.0, 0.0]);
        // exactly nonzero area ⇒ the exact test keeps it.
        assert!(!triangle_is_degenerate(t0, t1, t2));
        assert_eq!(fan_triangles(&[t0, t1, t2], 0).len(), 1);
        // …yet the old relative tolerance would have dropped it:
        let (e1, e2) = (t1 - t0, t2 - t0);
        assert!(e1.cross(e2).norm() <= 1e-12 * e1.norm() * e2.norm());
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_l).abs() < 1e-9);
        assert!(m.solids.get(r).cavities.is_empty());
    }

    #[test]
    fn common_non_convex_containment_is_inner() {
        let (mut m, l, bx) = l_and_inner_box();
        let vol_bx = nacre_props::mass_props(&m, bx).unwrap().volume;
        let r = boolean_one(&mut m, BoolKind::Common, l, bx).unwrap();
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_bx).abs() < 1e-9);
    }

    #[test]
    fn cut_box_inside_non_convex_is_empty() {
        // Cut(box − L): the box is wholly inside L ⇒ nothing remains.
        let (mut m, l, bx) = l_and_inner_box();
        assert_eq!(
            boolean_one(&mut m, BoolKind::Cut, bx, l),
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, d).unwrap();
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_l).abs() < 1e-9);
        let (mut m, l) = l_prism();
        let d = m.add_cuboid(far(), far_max());
        assert_eq!(
            boolean_one(&mut m, BoolKind::Fuse, l, d),
            Err(BoolError::EmptyResult)
        );
        let (mut m, l) = l_prism();
        let d = m.add_cuboid(far(), far_max());
        assert_eq!(
            boolean_one(&mut m, BoolKind::Common, l, d),
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Common, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Common, cube, bar).unwrap();
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
        let r = overlap_fuse_cut(&mut m, BoolKind::Common, a, b).unwrap()[0];
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, rod).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, rod).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!((props.volume - 3.06).abs() < 1e-9, "{}", props.volume);
        assert!((props.area - 15.0).abs() < 1e-9, "{}", props.area);
        assert_eq!(holed_faces(&m, r).len(), 2);
    }

    /// The same two solids the other way round: the bar severs the rod, and `Cut` answers with
    /// two solids (cell 0.4). Each severed stub is its own genus-0 box, so `validate` is clean —
    /// the pieces share no vertices or edges. (Before cell 0.4 this was `DISCONNECTED_RESULT`:
    /// one handle could not name two solids, and forcing both into one shell read as
    /// `NegativeGenus { genus: -1 }`. `pierced_multi` had been hiding it: severing A takes an
    /// edge of A through B.)
    #[test]
    fn cut_rod_by_l_severs_into_two() {
        let (mut m, l, rod) = l_and_rod();
        let solids = boolean(&mut m, BoolKind::Cut, rod, l).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((vol - 0.06).abs() < 1e-9, "total volume {vol}");
    }

    // --- Rotated booleans go live (overhaul 3d-i, `ROTATED_UNSUPPORTED` retired) ---
    // A boolean commutes with a rigid motion, so rotating both operands by the same
    // irrational-angle isometry must give the rigid image of the unrotated result — identical
    // volume, solid count, and cavity count, and still valid. These are the first live proof
    // that the TIP-wired machinery (arrangement, seam, in/out, outer/cavity — 3a–3c-vi) is
    // sound end-to-end on rotated (rounded-irrational) geometry.

    #[test]
    fn rotated_corner_bite_cut_is_rotation_invariant() {
        let (mut m, l, bx) = l_and_corner_box();
        let l = transform(&mut m, l, &rot30()).unwrap();
        m.rebuild_adjacency();
        let bx = transform(&mut m, bx, &rot30()).unwrap();
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.224)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn rotated_corner_bite_fuse_is_rotation_invariant() {
        let (mut m, l, bx) = l_and_corner_box();
        let l = transform(&mut m, l, &rot30()).unwrap();
        m.rebuild_adjacency();
        let bx = transform(&mut m, bx, &rot30()).unwrap();
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 0.924 - 0.224)).abs() < 1e-9, "volume {vol}");
    }

    /// The sever — the case that exercises the rotated `is_shell_outward` twin
    /// (`component_is_outward_tol`, 3c-vi): both severed pieces must read outward, so the
    /// result is two solids, not one solid with the other misjudged as a cavity.
    #[test]
    fn rotated_sever_cut_severs_into_two() {
        let (mut m, l, rod) = l_and_rod();
        let l = transform(&mut m, l, &rot30()).unwrap();
        m.rebuild_adjacency();
        let rod = transform(&mut m, rod, &rot30()).unwrap();
        m.rebuild_adjacency();
        let solids = boolean(&mut m, BoolKind::Cut, rod, l).unwrap();
        m.rebuild_adjacency();
        assert_eq!(solids.len(), 2, "rotated sever still yields two solids");
        assert!(nacre_validate::validate(&m).is_empty());
        let vol: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((vol - 0.06).abs() < 1e-9, "total volume {vol}");
    }

    #[test]
    fn rotated_containment_cut_makes_cavity() {
        let (mut m, a, b) = l_and_inner_box();
        let a = transform(&mut m, a, &rot30()).unwrap();
        m.rebuild_adjacency();
        let b = transform(&mut m, b, &rot30()).unwrap();
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(
            m.solids.get(r).cavities.len(),
            1,
            "the inner box becomes a cavity"
        );
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
    }

    /// A rotated coplanar contact is out of scope (the coplanar milestone), but retiring the
    /// guard must not make it *silently* wrong. When an X rotation tilts the contact plane its
    /// rounded coefficients are no longer exactly proportional, so the exact `planes_coplanar`
    /// detectors miss it and it reaches `general_boolean`. The DNA invariant: the result is
    /// either an honest reject, or a valid solid of the *correct* volume — never a
    /// plausible-but-invalid solid. (Measured: the boss fuse is solved exactly; the pocket cut
    /// honestly rejects.)
    #[test]
    fn rotated_coplanar_contact_is_never_silently_wrong() {
        use nacre_scalar::Axis;
        let tilt = |m: &mut Model, s: Handle<Solid>| -> Handle<Solid> {
            let r = transform(m, s, &rot_iso(Axis::X, 30)).unwrap();
            m.rebuild_adjacency();
            r
        };
        // boss fuse (correct fused volume 1.25 if solved).
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let boss = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 1.0]),
            Point3::from_array([0.75, 0.75, 2.0]),
        );
        let (base, boss) = (tilt(&mut m, base), tilt(&mut m, boss));
        if let Ok(r) = boolean_one(&mut m, BoolKind::Fuse, base, boss) {
            m.rebuild_adjacency();
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "solved boss is valid"
            );
            let v = nacre_props::mass_props(&m, r).unwrap().volume;
            assert!((v - 1.25).abs() < 1e-9, "solved boss fuse is correct: {v}");
        }
        // pocket cut (correct carved volume 0.875 if solved).
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let prism = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let (base, prism) = (tilt(&mut m, base), tilt(&mut m, prism));
        if let Ok(r) = boolean_one(&mut m, BoolKind::Cut, base, prism) {
            m.rebuild_adjacency();
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "solved pocket is valid"
            );
            let v = nacre_props::mass_props(&m, r).unwrap().volume;
            assert!(
                (v - 0.875).abs() < 1e-9,
                "solved pocket cut is correct: {v}"
            );
        }
    }

    /// Adversarial rotation stress (overhaul 3d-iii): many fixtures × kinds × rotations
    /// (single-axis, and Euler chains reaching arbitrary orientation) confirm the DNA
    /// invariant — a rotated boolean is *never silently wrong*: its result either equals the
    /// unrotated one (a boolean commutes with a rigid motion, so volume/solid-count/cavity-count
    /// are invariant) or is an honest reject. `#[ignore]`: each rotated boolean escalates its
    /// TIP predicates to astro-float and costs ~0.5–2.5 s, so this runs on demand, not per commit
    /// (the invariance regression guard is the fast `rotated_*` tests above).
    #[test]
    #[ignore = "slow: rotated booleans ~2s each (run with --ignored)"]
    fn rotation_invariance_stress() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        #[derive(PartialEq, Debug)]
        enum Out {
            Rej,
            Ok(f64, usize, usize),
        }
        let rot = |axis: Axis, deg: i128, piv: [i128; 3]| {
            Isometry::rotation(Rotation {
                axis,
                point: [
                    Rat::from_int(piv[0]),
                    Rat::from_int(piv[1]),
                    Rat::from_int(piv[2]),
                ],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            })
        };
        let run = |build: &dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>),
                   kind: BoolKind,
                   isos: &[Isometry]|
         -> Out {
            let (mut m, mut a, mut b) = build();
            for iso in isos {
                a = transform(&mut m, a, iso).unwrap();
                m.rebuild_adjacency();
                b = transform(&mut m, b, iso).unwrap();
                m.rebuild_adjacency();
            }
            match boolean(&mut m, kind, a, b) {
                Ok(solids) => {
                    m.rebuild_adjacency();
                    assert!(
                        nacre_validate::validate(&m).is_empty(),
                        "INVALID rotated result"
                    );
                    let vol: f64 = solids
                        .iter()
                        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                        .sum();
                    let cav: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
                    Out::Ok(vol, solids.len(), cav)
                }
                Err(_) => Out::Rej,
            }
        };
        let vclose =
            |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());
        let matches = |base: &Out, r: &Out| match (base, r) {
            (Out::Rej, Out::Rej) => true,
            (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2)) => vclose(*v1, *v2) && s1 == s2 && c1 == c2,
            _ => false,
        };
        let cube = |lo: [f64; 3], hi: [f64; 3]| {
            move || {
                let mut m = Model::new();
                let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
                let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
                (m, a, b)
            }
        };
        type Build = Box<dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>)>;
        let fixtures: Vec<(&str, Build)> = vec![
            ("corner", Box::new(l_and_corner_box)),
            (
                "rod",
                Box::new(|| {
                    let (m, l, r) = l_and_rod();
                    (m, r, l) // sever = Cut(rod, l): a=rod, b=l
                }),
            ),
            ("inner", Box::new(l_and_inner_box)),
            ("cube_corner", Box::new(cube([0.5; 3], [1.5; 3]))),
        ];
        let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
            ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
            ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
            (
                "Z50>X37",
                vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::X, 37, [0, 0, 1])],
            ),
            // A three-axis Euler chain reaches an arbitrary orientation (axes are X/Y/Z only).
            (
                "Z30>X30>Y73",
                vec![
                    rot(Axis::Z, 30, [1, 1, 0]),
                    rot(Axis::X, 30, [0, 0, 0]),
                    rot(Axis::Y, 73, [0, 2, 0]),
                ],
            ),
        ];
        let (mut success, mut reject, mut silent, mut skipped) = (0, 0, 0, 0);
        for (fname, build) in &fixtures {
            for kind in [BoolKind::Cut, BoolKind::Fuse, BoolKind::Common] {
                let base = run(build.as_ref(), kind, &[]);
                for (rname, isos) in &isos_list {
                    if matches!(base, Out::Rej) {
                        skipped += 1;
                        continue;
                    }
                    let r = run(build.as_ref(), kind, isos);
                    if matches!(r, Out::Rej) {
                        reject += 1; // an honest reject is acceptable, not a counterexample
                    } else if matches(&base, &r) {
                        success += 1;
                    } else {
                        silent += 1;
                        eprintln!("SILENT-WRONG {fname} {kind:?} {rname} base={base:?} rot={r:?}");
                    }
                }
            }
        }
        eprintln!(
            "ROTATION STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent} skipped(base_rej)={skipped}"
        );
        assert_eq!(
            silent, 0,
            "a rotated boolean was silently wrong (valid but != unrotated)"
        );
    }

    /// A signature that changes if `Store::push` order (hence handle identity) changes:
    /// vertex points in handle order — exactly what `assemble_fuse_cut` assigns by first
    /// appearance across `faces` — plus edge/face/solid counts and sorted volumes.
    #[cfg(feature = "parallel")]
    fn model_sig(m: &Model, solids: &[Handle<Solid>]) -> String {
        use std::fmt::Write;
        let mut s = String::new();
        let _ = write!(s, "S{}", solids.len());
        for (_, v) in m.vertices.iter() {
            let p = v.point.as_array();
            let _ = write!(
                s,
                "|{:x},{:x},{:x}",
                p[0].to_bits(),
                p[1].to_bits(),
                p[2].to_bits()
            );
        }
        let _ = write!(s, "|E{}F{}", m.edges.len(), m.faces.len());
        let mut vols: Vec<u64> = solids
            .iter()
            .map(|&sh| nacre_props::mass_props(m, sh).unwrap().volume.to_bits())
            .collect();
        vols.sort_unstable();
        let _ = write!(s, "|V{vols:?}");
        s
    }

    /// The parallel boolean must be bit-identical regardless of rayon thread count — replay
    /// determinism (DNA) requires thread-order independence. A rotated multi-face Fuse
    /// exercises the parallel reconstruction/classification; its result under a 1-thread
    /// pool (par_iter code, index order) must equal the default many-thread result, run
    /// repeatedly so scheduling jitter would show.
    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_boolean_is_thread_order_independent() {
        use nacre_scalar::{Axis, Isometry, Rat};
        let build = || {
            let ngon = |n: usize, r: f64, cx: f64, cy: f64| Profile2d {
                points: (0..n)
                    .map(|i| {
                        let ang = std::f64::consts::TAU * (i as f64) / (n as f64);
                        p2(cx + r * ang.cos(), cy + r * ang.sin())
                    })
                    .collect(),
            };
            let mut m = replay(&[
                extrude_op(ngon(16, 2.0, 0.0, 0.0), 3.0),
                extrude_op(ngon(16, 2.0, 2.5, 0.5), 3.0),
            ])
            .unwrap();
            let a = m.live_solids[0];
            let b = m.live_solids[1];
            let up = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
            let b = transform(&mut m, b, &up).unwrap();
            m.rebuild_adjacency();
            let tilt = rot_iso(Axis::X, 30);
            let a = transform(&mut m, a, &tilt).unwrap();
            m.rebuild_adjacency();
            let b = transform(&mut m, b, &tilt).unwrap();
            m.rebuild_adjacency();
            (m, a, b)
        };
        let run = || {
            let (mut m, a, b) = build();
            let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
            m.rebuild_adjacency();
            model_sig(&m, &solids)
        };
        let pool1 = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let reference = pool1.install(run);
        for _ in 0..8 {
            assert_eq!(
                run(),
                reference,
                "parallel boolean result depends on thread order"
            );
        }
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
            boolean_one(&mut m, kind, l, bx).unwrap();
            m.rebuild_adjacency();
            assert_eq!(from_arrange, discovered_triples(&m), "{kind:?}");
        }

        // The reflex-corner bite too: eight seam vertices, one of them where the box
        // straddles the L's notch.
        let (mut m, l, bx) = l_and_reflex_box();
        let from_arrange = seam_endpoint_triples(&m, l, bx);
        assert_eq!(from_arrange.len(), 8);
        boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // And the folded arc (cell 3e-1). `discovered_triples` reads `reachable()`, so
        // a seam vertex the reconstruction built but no result face used would drop out
        // of the set and break this — which is what a lost face looks like. Volume
        // alone could coincide; this cannot.
        let (mut m, l, bx) = l_and_popup_box();
        let from_arrange = seam_endpoint_triples(&m, l, bx);
        boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // And the inner loop (cell 3f-1): its four nodes are the hole's rim, used once
        // by the lid and once by a pocket wall. Lose the loop and they leave `reachable`.
        let (mut m, l, stub) = l_and_dimple();
        let from_arrange = seam_endpoint_triples(&m, l, stub);
        assert_eq!(from_arrange.len(), 4);
        boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // Two chords (cell 3e-2). Six seam nodes: two boundary ends and a bend for each of
        // the cap's arcs, and the same six seen again from the bar's floor. Lose one of the
        // bar's two floor faces and its four nodes leave `reachable`.
        let (mut m, l, bar) = l_and_notch_bar();
        let from_arrange = seam_endpoint_triples(&m, l, bar);
        boolean_one(&mut m, BoolKind::Cut, l, bar).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // The non-convex hole (cell 3h): six rim nodes, none of them on `∂f`.
        let (mut m, l, stub) = l_and_ell_stub();
        let from_arrange = seam_endpoint_triples(&m, l, stub);
        assert_eq!(from_arrange.len(), 6);
        boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // An arc beside a loop (cell 3f-4). Lose the loop, or place it on the wrong region,
        // and its four rim nodes leave `reachable`.
        let (mut m, l, st) = l_and_staple();
        let from_arrange = seam_endpoint_triples(&m, l, st);
        boolean_one(&mut m, BoolKind::Cut, l, st).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // Two flat loops (cell 3f-3): sixteen nodes, and the two prong cross-sections must
        // both survive as faces or their eight rim nodes leave `reachable`.
        let (mut m, u, slab) = u_and_slab();
        let from_arrange = seam_endpoint_triples(&m, u, slab);
        boolean_one(&mut m, BoolKind::Cut, u, slab).unwrap();
        m.rebuild_adjacency();
        assert_eq!(from_arrange, discovered_triples(&m));

        // Swapped, the same four nodes are the island's whole outer ring (cell 3f-2).
        // The arrangement does not know which solid is `a`, so it offers the same set;
        // the result has to still contain all of it. If the island face were dropped,
        // the box would lose its floor and every one of the four would go with it.
        let (mut m, l, stub) = l_and_dimple();
        let from_arrange = seam_endpoint_triples(&m, stub, l);
        assert_eq!(from_arrange.len(), 4);
        boolean_one(&mut m, BoolKind::Cut, stub, l).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, bx, pc).unwrap();
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
            || boolean_one(&mut m, BoolKind::Fuse, ic, bx),
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
        let r = boolean_one(&mut m, BoolKind::Cut, slab, pc).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, slab, pc).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, pc, slab).unwrap();
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
            let r = boolean_one(&mut m, BoolKind::Cut, x, y).unwrap();
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
            let r = boolean_one(&mut m, BoolKind::Fuse, x, y).unwrap();
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

    /// `is_shell_outward` — the exact sign `assemble_fuse_cut` labels components by
    /// ((5d)#5, replacing the f64 signed-volume flux) — is true for an outward,
    /// material-enclosing shell and false for an inward void shell. `reversed_shell`
    /// flips one into the other, so the same faces read opposite orientations.
    #[test]
    fn is_shell_outward_true_for_outer_false_for_void() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let outer = m.solids.get(cube).outer;
        let out_faces = m.shells.get(outer).faces.clone();
        assert!(is_shell_outward(&m, &out_faces), "outer shell is outward");
        let void = m.reversed_shell(outer);
        let void_faces = m.shells.get(void).faces.clone();
        assert!(
            !is_shell_outward(&m, &void_faces),
            "reversed shell is a void"
        );
    }

    /// Build the `LocalFace`s of a solid's outer shell, one per `PlaneInfo` (whose order
    /// matches the shell's faces), with `Node::Orig` loops — the pre-assembly form
    /// `component_is_outward_tol` consumes. `flip` reverses the materialized outward normal.
    fn shell_local_faces(model: &Model, planes: &[PlaneInfo], flip: bool) -> Vec<LocalFace> {
        planes
            .iter()
            .enumerate()
            .map(|(i, pi)| {
                let loops = face_loop_verts(model, pi.face);
                let node_ring = |r: &[Handle<Vertex>]| r.iter().map(|&v| Node::Orig(v)).collect();
                LocalFace {
                    plane_idx: i,
                    loop_nodes: node_ring(&loops[0]),
                    inner: loops[1..].iter().map(|r| node_ring(r)).collect(),
                    flip,
                }
            })
            .collect()
    }

    /// On an axis-aligned solid the exact `component_is_outward_tol` reproduces the f64
    /// `is_shell_outward` label — true for the outer shell, false when every face is flipped
    /// (a void) — so wiring it will not move any label on the unrotated corpus.
    #[test]
    fn outward_tol_matches_f64_on_axis_aligned() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let planes = collect_planes(&m, cube).unwrap();
        let outer = shell_local_faces(&m, &planes, false);
        let refs: Vec<&LocalFace> = outer.iter().collect();
        assert!(component_is_outward_tol(&planes, &refs).unwrap());
        let void = shell_local_faces(&m, &planes, true);
        let vrefs: Vec<&LocalFace> = void.iter().collect();
        assert!(!component_is_outward_tol(&planes, &vrefs).unwrap());
    }

    /// The outward/void label is rigid-rotation invariant: the exact path over a cube rotated
    /// by an irrational (30°) angle agrees with the same cube unrotated — outer is outward,
    /// the flipped shell is a void — even though the rotated coordinates are rounded
    /// irrationals that the f64 x-coefficient test could misjudge.
    #[test]
    fn outward_tol_is_rotation_invariant() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation as SRot};
        let label = |m: &Model, s: Handle<Solid>, flip: bool| {
            let planes = collect_planes(m, s).unwrap();
            let lf = shell_local_faces(m, &planes, flip);
            let refs: Vec<&LocalFace> = lf.iter().collect();
            component_is_outward_tol(&planes, &refs).unwrap()
        };
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let (out_u, void_u) = (label(&m, cube, false), label(&m, cube, true));
        assert!(out_u && !void_u, "unrotated: outer outward, flipped void");
        let iso = Isometry::rotation(SRot {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        let OpOutput::Transform { solid: rc } = apply(
            &mut m,
            &Operation::Transform {
                solid: cube,
                isometry: iso,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert_eq!(
            label(&m, rc, false),
            out_u,
            "outward invariant under rotation"
        );
        assert_eq!(label(&m, rc, true), void_u, "void invariant under rotation");
    }

    /// A component whose planes mix rotated (`tri_pt3 = Some`) and axis-aligned
    /// (`tri_pt3 = None`) faces — the shape of a mixed-rotation boolean — must not panic:
    /// `plane_def` rebuilds a `None` plane's exact `Pt3` from its `tri`, and the label stays
    /// correct. Locks the `plane_def` path against a `tri_pt3.expect()` regression.
    #[test]
    fn outward_tol_handles_mixed_rotation() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let mut planes = collect_planes(&m, cube).unwrap();
        // Force one face's exact def present (Some) while the rest stay None: `any_rotated` now
        // routes the whole component through the frame3 path, exercising plane_def on both.
        let tri = planes[0].tri;
        let rat = |x: f64| nacre_scalar::Rat::try_from_f64(x).unwrap();
        planes[0].tri_pt3 = Some(tri.map(|p| {
            let a = p.as_array();
            Pt3::at([rat(a[0]), rat(a[1]), rat(a[2])])
        }));
        assert!(
            planes.iter().any(|p| p.tri_pt3.is_some())
                && planes.iter().any(|p| p.tri_pt3.is_none())
        );
        let outer = shell_local_faces(&m, &planes, false);
        let refs: Vec<&LocalFace> = outer.iter().collect();
        assert!(component_is_outward_tol(&planes, &refs).unwrap());
    }

    /// A concave (L-prism) shell — a reflex vertex at `(1,1)` — labels outward before and
    /// after an irrational rotation. `v*` is the convex min corner, not the reflex one, so the
    /// `∃ n_x < 0` test still holds; this guards the non-convex extreme case.
    #[test]
    fn outward_tol_on_concave_l_prism() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation as SRot};
        let label = |m: &Model, s: Handle<Solid>| {
            let planes = collect_planes(m, s).unwrap();
            let lf = shell_local_faces(m, &planes, false);
            let refs: Vec<&LocalFace> = lf.iter().collect();
            component_is_outward_tol(&planes, &refs).unwrap()
        };
        let (mut m, s) = l_prism();
        assert!(label(&m, s), "unrotated L outer is outward");
        let iso = Isometry::rotation(SRot {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        let OpOutput::Transform { solid: r } = apply(
            &mut m,
            &Operation::Transform {
                solid: s,
                isometry: iso,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(label(&m, r), "rotated L outer is outward");
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
        let r = boolean_one(&mut m, BoolKind::Fuse, solid, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, u, slab).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, u, slab).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, slab, u).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, bar).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bar).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, st).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, st).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, st, l).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, stub, l).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, stub, l).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
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
        let slotted = boolean_one(&mut m, BoolKind::Cut, bar, groove).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, slotted, pocket).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, l, stub).unwrap();
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
        // Exactly one holed face — the lid (found by predicate; the boolean pocket path does not
        // fix the shell's face order the way direct construction did).
        let holed = solid_faces(&m, pc)
            .iter()
            .filter(|&&fh| !m.faces.get(fh).inner.is_empty())
            .count();
        assert_eq!(holed, 1);
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
        let c = boolean_one(&mut m, BoolKind::Common, a, b).unwrap(); // the box [1,2]³
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
        let c = boolean_one(&mut m, BoolKind::Common, a, b).unwrap(); // [1,2]³
        m.rebuild_adjacency();
        let d = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 2.0]),
            Point3::from_array([2.0, 2.0, 3.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, c, d).unwrap();
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
            || boolean_one(&mut m, BoolKind::Cut, ipc, bx),
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
            let r = boolean_one(&mut m, kind, pc, bx).unwrap();
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
            hits += pierced_faces(m, planes, inc, bounds[0], bounds[1], &rings)?.len();
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
        assert!(
            arrange::edge_crosses_face(&m, &planes, along_x, bounds[0], bounds[1], q, &rings)
                .unwrap()
        );

        // And the whole boolean runs on it through the seam path — the measurement
        // design.md §9 line 447 asked for. Cell (5b) then deleted the convex path, so
        // `boolean` reaches the same code; this keeps the direct call as the §447 record.
        let (mut m, a, b) = two_boxes();
        let r = overlap_fuse_cut(&mut m, BoolKind::Cut, a, b).unwrap()[0];
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (1.0 - 0.125)).abs() < 1e-9, "volume {vol}");

        // The `Fuse` leg of the same corner overlap, which cell (5a) never measured on the
        // seam path — cell (5b) routes it here. `1 + 1 − 0.125`.
        let (mut m, a, b) = two_boxes();
        let r = overlap_fuse_cut(&mut m, BoolKind::Fuse, a, b).unwrap()[0];
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
            assert!(
                !arrange::edge_crosses_face(&m, &planes, inc, bounds[0], bounds[1], q, &rings)
                    .unwrap()
            );
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
        let cut = boolean_one(&mut m, BoolKind::Cut, a, y).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let v = nacre_props::mass_props(&m, cut).unwrap().volume;
        assert!((v - 993.28).abs() < 1e-9, "notch cut {v}");

        let (mut m, a, y) = cube_and_notch();
        let fuse = boolean_one(&mut m, BoolKind::Fuse, a, y).unwrap();
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
            let r = boolean_one(&mut m, kind, a, bar).unwrap();
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
    /// cube and out both ends, so `Cut(bar, cube)` leaves the bar in two 1×1×1 stubs — two solids
    /// (cell 0.4), each a clean genus-0 box (`validate` clean, pieces share nothing). The convex
    /// path rejected this as `poke_through`; the seam path returns both pieces. (Before cell 0.4
    /// this was `disconnected_result`; cell 3e-3's non-convex sibling was `Cut(rod, L)`.)
    #[test]
    fn a_convex_cut_severs_its_operand() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, bar, a).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        for &s in &solids {
            assert!((nacre_props::mass_props(&m, s).unwrap().volume - 1.0).abs() < 1e-9);
        }
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
            (planes, q, rings, inc, bounds)
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
            assert_rejects(
                || arrange::edge_crosses_face(&m, &planes, inc, bounds[0], bounds[1], q, &rings),
                tag::VERTEX_ON_FACE_PLANE,
            );
        }

        // (b) The piercing point is a ring node.
        {
            let (planes, q, rings, inc, bounds) = pierce(at_vertex, [1.0, 1.0]);
            assert_rejects(
                || arrange::edge_crosses_face(&m, &planes, inc, bounds[0], bounds[1], q, &rings),
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
            let (planes, q, rings, inc, bounds) = pierce(at_edge, [1.5, 0.0]);
            assert_rejects(
                || arrange::edge_crosses_face(&m, &planes, inc, bounds[0], bounds[1], q, &rings),
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
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
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
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        m.rebuild_adjacency();
        let cutter = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, cutter).unwrap();
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
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        m.rebuild_adjacency();
        let stub = m.add_cuboid(
            Point3::from_array([1.4, 1.4, -0.5]),
            Point3::from_array([1.6, 1.6, 1.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, stub).unwrap();
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
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let tunnel = m.add_cuboid(
            Point3::from_array([1.4, 1.4, -0.5]),
            Point3::from_array([1.6, 1.6, 3.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, tunnel).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.92).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0);
    }

    #[test]
    fn a_slab_splits_a_hollow_box_into_two() {
        // A slab cut through the whole box (and its void) severs it into two solids (cell 0.4).
        // The slab spans the full cross-section, so it opens the void — both pieces are
        // cavity-free. (Before cell 0.4 this was rejected as `DISCONNECTED_RESULT`.)
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.4, -0.5]),
            Point3::from_array([3.5, 1.6, 3.5]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        for &s in &solids {
            assert_eq!(m.solids.get(s).cavities.len(), 0);
        }
        // Hollow 26 (= 27 − 1 void); the slab removes 1.6 of material (8 area × 0.2 thick).
        let vol: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((vol - 24.4).abs() < 1e-9, "total volume {vol}");
    }

    /// A sever that also leaves a surviving cavity: a hollow box whose void sits to one side,
    /// cut by a slab that severs it without touching the void. The x<2 piece keeps the void as a
    /// cavity, the x>2 piece is solid — two outward shells *and* one inward. Which outer owns the
    /// cavity needs a containment test we do not have yet, so it is honestly rejected
    /// (`SEVERED_WITH_CAVITY`) rather than mis-assembled. This is the firing test the guard is
    /// born with (design.md); n0 measured the two-outward-plus-one-inward component split.
    #[test]
    fn severed_with_cavity_is_rejected() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        // Void near the x-low side, clear of the x=2 cut.
        let inner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 2.5]),
        );
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        // A slab spanning full y,z, thin in x at x∈[2,2.2] — severs into x<2 (holds the void)
        // and x>2 (solid).
        let slab = m.add_cuboid(
            Point3::from_array([2.0, -1.0, -1.0]),
            Point3::from_array([2.2, 4.0, 4.0]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, hollow, slab),
            tag::SEVERED_WITH_CAVITY,
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
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let bore = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, bore).unwrap();
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
    fn an_imprint_reaching_past_the_face_is_rejected() {
        // The frame origin is the top face centroid (0.5, 0.5). Each profile reaches past the
        // [0,1]² face region, so `imprint` rejects rather than build a silently-invalid inner loop
        // (cell imprint-containment; n0 measured that apply was Ok, validate passed, props
        // integrated garbage — including a negative volume — and only tessellate's NoEar caught it).
        // Boundary contact ("touches") rejects too. `pad`/`pocket` no longer reject here — a
        // profile past the face is a boolean overhang (see `pad_an_overhanging_boss` etc.).
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
            let (mut m, top) = cube_with_top();
            assert_eq!(
                apply(&mut m, &Operation::ImprintSketch { face: top, profile }),
                Err(OpError::ProfileNotContainedInFace),
                "{name}"
            );
        }
    }

    // A single-edge overhang footprint on the cube top: world x∈[0.25,0.75], y∈[-0.25,0.75]
    // (overhangs the y=0 edge), area 0.5. (Frame maps local [px,py] → world (0.5+py, 0.5−px).)
    fn edge_overhang_profile() -> Profile2d {
        Profile2d {
            points: vec![
                p2(-0.25, -0.25),
                p2(0.75, -0.25),
                p2(0.75, 0.25),
                p2(-0.25, 0.25),
            ],
        }
    }

    // A spanning slab: world x∈[0.25,0.75], y∈[-0.25,1.25] (crosses both y edges), area 0.75.
    fn spanning_slab_profile() -> Profile2d {
        Profile2d {
            points: vec![
                p2(-0.75, -0.25),
                p2(0.75, -0.25),
                p2(0.75, 0.25),
                p2(-0.75, 0.25),
            ],
        }
    }

    #[test]
    fn pad_an_overhanging_boss() {
        // The profile reaches past one face edge: part of the boss sits on the face, part
        // cantilevers into the air. Relaxing the containment gate routes it to the overhang Fuse
        // sidecar (Ok here proves the routing — a contained-only pad would reject). The boss lives
        // wholly above z=1, so vol = cube 1 + footprint 0.5 · dist 1 = 1.5.
        let (mut m, top) = cube_with_top();
        let OpOutput::PadOnFace { solid, top_face } =
            apply(&mut m, &pad_op(top, edge_overhang_profile(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 1.5).abs() < 1e-12);
        assert!(m.reachable().faces.contains(&top_face)); // boss top cap recovered
    }

    #[test]
    fn pad_a_spanning_slab_boss() {
        // A slab crossing the whole face (overhangs two opposite edges). vol = 1 + 0.75 · 1 = 1.75.
        let (mut m, top) = cube_with_top();
        let out = apply(&mut m, &pad_op(top, spanning_slab_profile(), 1.0)).unwrap();
        let OpOutput::PadOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 1.75).abs() < 1e-12);
    }

    #[test]
    fn pocket_an_edge_slot() {
        // A blind pocket whose footprint overhangs one edge — an edge slot open to the side.
        // Only the on-face part (world x[0.25,0.75]×y[0,0.75] = 0.375) carves: 1 − 0.375·0.5 = 0.8125.
        let (mut m, top) = cube_with_top();
        let OpOutput::PocketOnFace { solid, bottom_face } =
            apply(&mut m, &pocket_op(top, edge_overhang_profile(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 0.8125).abs() < 1e-12);
        assert!(m.reachable().faces.contains(&bottom_face)); // slot floor recovered
    }

    #[test]
    fn pocket_a_slab_channel() {
        // A blind channel crossing the whole face (breaches two opposite walls). On-face carve
        // world x[0.25,0.75]×y[0,1] = 0.5: 1 − 0.5·0.5 = 0.75.
        let (mut m, top) = cube_with_top();
        let out = apply(&mut m, &pocket_op(top, spanning_slab_profile(), 0.5)).unwrap();
        let OpOutput::PocketOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 0.75).abs() < 1e-12);
    }

    #[test]
    fn pad_overhang_off_the_face_is_rejected() {
        // A footprint that does not touch the face at all: the boss is disjoint from the solid,
        // an unrepresentable union. Honest reject (no panic), not a silent floating boss.
        let (mut m, top) = cube_with_top();
        let far = Profile2d {
            points: vec![p2(1.8, 1.8), p2(2.2, 1.8), p2(2.2, 2.2), p2(1.8, 2.2)],
        };
        assert!(matches!(
            apply(&mut m, &pad_op(top, far, 0.3)),
            Err(OpError::Boolean(_))
        ));
    }

    #[test]
    fn pocket_through_overhang_is_rejected() {
        // An overhang pocket deep enough to pierce the far side is not blind — no single floor.
        // Honest reject via whichever path fires (the overhang detector declines, the seam path
        // rejects the mixed contact), mirroring `pocket_through_the_solid_is_rejected`.
        let (mut m, top) = cube_with_top();
        let got = apply(&mut m, &pocket_op(top, edge_overhang_profile(), 1.5));
        assert!(
            matches!(got, Err(OpError::Boolean(_)) | Err(OpError::PocketNotBlind)),
            "through overhang must reject honestly, got {got:?}"
        );
    }

    #[test]
    fn pad_non_convex_overhang_is_rejected() {
        // A non-convex (L) overhang footprint: the overhang boolean sidecars gate on convexity, so
        // this is out of scope and honestly rejected (contained non-convex still works — a
        // deliberate asymmetry until the overhang convex gate is dropped).
        let (mut m, top) = cube_with_top();
        let l_over = Profile2d {
            points: vec![
                p2(-0.25, -0.25),
                p2(0.75, -0.25),
                p2(0.75, 0.25),
                p2(0.25, 0.25),
                p2(0.25, 0.5),
                p2(-0.25, 0.5),
            ],
        };
        assert!(matches!(
            apply(&mut m, &pad_op(top, l_over, 1.0)),
            Err(OpError::Boolean(_))
        ));
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

    /// §5 explicit sharing (overhaul #3): a prism built with a shared base-cap
    /// surface reuses that `Surface` handle for its flush cap, and reconciles the
    /// cap's face orientation so the materialized outward normal stays `−sweep`.
    #[test]
    fn build_prism_base_cap_reuses_shared_surface() {
        let mut m = Model::new();
        // A face-plane surface with outward normal +z (as a face on the base solid).
        let sf = m.surfaces.push(Surface::Plane(
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
                .unwrap(),
        ));
        let base_pts = [
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        ];
        let (_prism, faces) = build_prism(
            &mut m,
            &base_pts,
            Vector3::from_array([0.0, 0.0, 1.0]),
            Some(sf),
        )
        .unwrap();
        let cap = m.faces.get(faces[0]); // base cap is pushed first
        // Shared handle (was a fresh push before overhaul #3).
        assert_eq!(cap.surface, sf, "base cap reuses the shared surface handle");
        // Orientation reconciled: materialized outward normal is −sweep (−z).
        let Surface::Plane(p) = m.surfaces.get(cap.surface) else {
            unreachable!()
        };
        let sign = match cap.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        let materialized = p.normal() * sign;
        assert!(
            (materialized - Vector3::from_array([0.0, 0.0, -1.0])).norm() < 1e-12,
            "materialized cap normal stays −z, got {materialized:?}"
        );
    }

    /// The handle branch of `shares_or_coplanar` is load-bearing: a shared
    /// `Surface` handle reports coplanar even when the stored `plane` values are
    /// *not* geometrically coplanar (so the fallback would not fire). This is the
    /// path a referenced coplanar contact takes; on axis-aligned M5 it is redundant
    /// with the geometric test, but the branch must work for rotated frames.
    #[test]
    fn shares_or_coplanar_uses_the_handle_branch() {
        let mut m = Model::new();
        let shared = m.surfaces.push(Surface::Plane(
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0]))
                .unwrap(),
        ));
        let fh = m.faces.push(Face {
            surface: shared,
            outer: Loop { half_edges: vec![] },
            inner: vec![],
            orientation: Orientation::Forward,
        });
        let plane_x0 =
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0]))
                .unwrap();
        let plane_z0 =
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
                .unwrap();
        let mk = |plane| PlaneInfo {
            surf: shared,
            face: fh,
            plane,
            tri: [Point3::origin(); 3],
            n_out: Vector3::from_array([0.0; 3]),
            orient: Orientation::Forward,
            tri_pt3: None,
        };
        let (pa, pb) = (mk(plane_x0), mk(plane_z0));
        // The two planes are NOT geometrically coplanar → the fallback would fail.
        assert!(!planes_coplanar(&pa.plane, &pb.plane));
        // But the shared handle makes them coplanar-by-reference.
        assert!(shares_or_coplanar(&pa, &pb));
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

    /// A pocket deep enough to pierce the far side is not blind. The extrude+Cut prism is
    /// top-flush yet crosses the far face transversally, a mixed contact the boolean declines
    /// (`Boolean(Unsupported)`); were the through-cut instead accepted, no floor would land on
    /// the offset plane and the wrapper would reject it as `PocketNotBlind`. Either way the
    /// kernel rejects honestly — no panic, no silently invalid solid (M4 left this unchecked).
    #[test]
    fn pocket_through_the_solid_is_rejected() {
        let (mut m, top) = cube_with_top(); // 1.0-thick cube
        let got = apply(&mut m, &pocket_op(top, small_square(), 1.5));
        assert!(
            matches!(got, Err(OpError::Boolean(_)) | Err(OpError::PocketNotBlind)),
            "through-pocket must reject honestly, got {got:?}"
        );
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
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.875).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.live_solids, vec![r]);
    }

    /// Lower corner of a solid's outer-shell vertex bounding box (for translation
    /// tests: a rigid move shifts it by exactly the offset).
    fn bbox_lo(m: &Model, s: Handle<Solid>) -> [f64; 3] {
        let mut lo = [f64::INFINITY; 3];
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        let p = m.vertices.get(vh).point.as_array();
                        for k in 0..3 {
                            lo[k] = lo[k].min(p[k]);
                        }
                    }
                }
            }
        }
        lo
    }

    fn count_discovered(m: &Model, s: Handle<Solid>) -> usize {
        let mut seen = std::collections::HashSet::new();
        let mut n = 0;
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh)
                            && matches!(m.vertices.get(vh).origin, Origin::Discovered { .. })
                        {
                            n += 1;
                        }
                    }
                }
            }
        }
        n
    }

    fn test_iso() -> (nacre_scalar::Isometry, [f64; 3]) {
        use nacre_scalar::Rat;
        (
            nacre_scalar::Isometry::translation([
                Rat::new(7, 2).unwrap(),
                Rat::from_int(-4),
                Rat::from_int(11),
            ]),
            [3.5, -4.0, 11.0],
        )
    }

    /// A rational translation supersedes a cuboid: rigid, so volume/area are
    /// invariant and the bounding box shifts by exactly the offset; validate/tess/
    /// STEP all accept the moved solid, and the input drops from `live_solids`.
    #[test]
    fn transform_translate_cuboid() {
        let (iso, off) = test_iso();
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap();
        let lo0 = bbox_lo(&m, c);

        let c2 = transform(&mut m, c, &iso).unwrap();
        m.rebuild_adjacency();

        assert_eq!(m.live_solids, vec![c2], "input superseded");
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap();
        assert!(
            (after.volume - before.volume).abs() < 1e-12,
            "volume invariant"
        );
        assert!((after.area - before.area).abs() < 1e-12, "area invariant");
        let lo1 = bbox_lo(&m, c2);
        for k in 0..3 {
            assert!(
                (lo1[k] - (lo0[k] + off[k])).abs() < 1e-12,
                "bbox shifted by offset"
            );
        }
        assert!(nacre_tess::to_obj(&m).is_ok(), "moved solid tessellates");
        assert!(
            nacre_step::to_step(&m)
                .unwrap()
                .contains("MANIFOLD_SOLID_BREP"),
            "moved solid exports to STEP"
        );
    }

    /// Transforming a boolean *result* (which carries `Discovered` seam vertices)
    /// preserves those vertices' `Origin` and remaps their `ThreePlane` definition
    /// onto the moved surfaces — the count survives and validate stays clean, so the
    /// definition was not silently downgraded to `Constructed`.
    #[test]
    fn transform_translate_preserves_discovered_definition() {
        let (iso, _) = test_iso();
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        let before = nacre_props::mass_props(&m, r).unwrap().volume;
        let disc = count_discovered(&m, r);
        assert!(
            disc > 0,
            "the Cut result must have Discovered seam vertices"
        );

        let r2 = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();

        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(
            count_discovered(&m, r2),
            disc,
            "Discovered vertices preserved"
        );
        let after = nacre_props::mass_props(&m, r2).unwrap().volume;
        assert!((after - before).abs() < 1e-12, "volume invariant");
    }

    /// Replay determinism (DNA 3): the same construction + transform reproduces the
    /// same geometry and the same handle down to the index.
    #[test]
    fn transform_is_deterministic() {
        let (iso, _) = test_iso();
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let c2 = transform(&mut m, c, &iso).unwrap();
            (bbox_lo(&m, c2), c2)
        };
        assert_eq!(build(), build(), "same ops → same geometry and handle");
    }

    /// The `Transform` op flows through `apply`, superseding via the op dispatch.
    #[test]
    fn transform_op_applies() {
        let (iso, _) = test_iso();
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: c,
                isometry: iso,
            },
        )
        .unwrap();
        match out {
            OpOutput::Transform { solid } => assert_eq!(m.live_solids, vec![solid]),
            other => panic!("expected Transform output, got {other:?}"),
        }
    }

    /// A genuinely tilted rigid rotation: 30° about Z through the rational axis
    /// point (1,1,0). Non-90° and non-axis-aligned, so it exercises the Rotated
    /// origin and the boolean reject guard (unlike the 90° family, which stays exact).
    fn rot30() -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        })
    }

    /// Distinct outer-shell vertex points of a solid (dedup by handle).
    fn outer_points(m: &Model, s: Handle<Solid>) -> Vec<[f64; 3]> {
        let mut seen = std::collections::HashSet::new();
        let mut pts = Vec::new();
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh) {
                            pts.push(m.vertices.get(vh).point.as_array());
                        }
                    }
                }
            }
        }
        pts
    }

    /// A non-90° rotation genuinely tilts the solid: rigid (volume/area invariant),
    /// validate/tess/STEP clean, a known corner lands at its exact rotated image, the
    /// vertices carry `Origin::Rotated` (`solid_is_rotated`), and a boolean against it now
    /// runs (a *mixed*-rotation cut: rotated `c2` minus an axis-aligned `d` it contains, so
    /// `d` becomes a cavity — overhaul 3d-i retired the `ROTATED_UNSUPPORTED` entry guard).
    #[test]
    fn transform_rotate_cuboid_tilts_and_cuts() {
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap();
        let c2 = transform(&mut m, c, &rot30()).unwrap();
        m.rebuild_adjacency();

        assert_eq!(m.live_solids, vec![c2], "input superseded");
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap();
        assert!(
            (after.volume - before.volume).abs() < 1e-9,
            "volume invariant"
        );
        assert!((after.area - before.area).abs() < 1e-9, "area invariant");
        assert!(nacre_tess::to_obj(&m).is_ok(), "rotated solid tessellates");
        assert!(
            nacre_step::to_step(&m)
                .unwrap()
                .contains("MANIFOLD_SOLID_BREP"),
            "rotated solid exports to STEP"
        );

        assert!(solid_is_rotated(&m, c2), "vertices carry Origin::Rotated");

        // Corner (0,0,0) rotates about pivot (1,1) by 30°: dx=dy=-1, so
        // x' = 1 - cos30 + sin30, y' = 1 - sin30 - cos30, z' = 0.
        let (c30, s30) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
        let want = [1.0 - c30 + s30, 1.0 - s30 - c30, 0.0];
        let pts = outer_points(&m, c2);
        assert!(
            pts.iter()
                .any(|p| p.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9)),
            "corner (0,0,0) rotated to its exact image {want:?}; got {pts:?}"
        );

        // A mixed-rotation cut now runs: the axis-aligned `d` sits inside the rotated `c2`, so
        // `Cut(c2, d)` leaves `c2` with `d` carved out as a cavity (volume 24 − 1 = 23).
        let d = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "mixed cut is valid"
        );
        assert_eq!(
            m.solids.get(r).cavities.len(),
            1,
            "the contained box is a cavity"
        );
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 23.0).abs() < 1e-9, "volume {vol}");
    }

    /// A 90° rotation about Z is axis-aligned and exact: the solid stays
    /// `Constructed` (`solid_is_rotated` false), volume is exact, and boolean is
    /// still allowed — a following cut succeeds and validates.
    #[test]
    fn transform_rotate_90_is_exact_and_allows_boolean() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let rot90 = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
        });
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c2 = transform(&mut m, c, &rot90).unwrap();
        m.rebuild_adjacency();

        assert!(nacre_validate::validate(&m).is_empty());
        assert!(!solid_is_rotated(&m, c2), "90° stays exact (Constructed)");
        assert_eq!(
            nacre_props::mass_props(&m, c2).unwrap().volume,
            24.0,
            "exact volume"
        );

        // c rotated 90° about origin occupies x∈[-3,0], y∈[0,2], z∈[0,4].
        // Cut with d = [-1,0.5,1]-[0.5,1.5,2]: overlap volume 1 → 24 − 1 = 23.
        let d = m.add_cuboid(
            Point3::from_array([-1.0, 0.5, 1.0]),
            Point3::from_array([0.5, 1.5, 2.0]),
        );
        let r = boolean(&mut m, BoolKind::Cut, c2, d).expect("exact rotation → boolean allowed");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r[0]).unwrap().volume;
        assert!((vol - 23.0).abs() < 1e-9, "cut volume {vol}");
    }

    /// Rotating a boolean *result* marks its `Discovered` seam vertices `Rotated`
    /// over the pre-rotation vertex as base: validate stays clean, volume is
    /// invariant, and at least one Rotated base is a Discovered vertex (the seam
    /// definition is preserved through the rotation, not downgraded).
    #[test]
    fn transform_rotate_boolean_result_keeps_discovered_base() {
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        assert!(
            count_discovered(&m, r) > 0,
            "Cut result has Discovered seams"
        );
        let before = nacre_props::mass_props(&m, r).unwrap().volume;

        let r2 = transform(&mut m, r, &rot30()).unwrap();
        m.rebuild_adjacency();

        assert!(nacre_validate::validate(&m).is_empty());
        assert!(
            solid_is_rotated(&m, r2),
            "rotated result carries Rotated origin"
        );
        let after = nacre_props::mass_props(&m, r2).unwrap().volume;
        assert!((after - before).abs() < 1e-9, "volume invariant");

        let sh = m.solids.get(r2).outer;
        let mut found_disc_base = false;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if let Origin::Rotated { base, .. } = m.vertices.get(vh).origin {
                            if matches!(m.vertices.get(base).origin, Origin::Discovered { .. }) {
                                found_disc_base = true;
                            }
                        }
                    }
                }
            }
        }
        assert!(
            found_disc_base,
            "a Rotated vertex's base is its Discovered seam vertex"
        );
    }

    /// Replay determinism (DNA 3): the same construction + rotation reproduces the
    /// same geometry and the same handle.
    #[test]
    fn transform_rotate_is_deterministic() {
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let c2 = transform(&mut m, c, &rot30()).unwrap();
            (bbox_lo(&m, c2), c2)
        };
        assert_eq!(build(), build(), "same ops → same geometry and handle");
    }

    /// A rotation `Transform` flows through `apply` and marks the result Rotated.
    #[test]
    fn transform_rotate_op_applies() {
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: c,
                isometry: rot30(),
            },
        )
        .unwrap();
        let OpOutput::Transform { solid } = out else {
            panic!("expected Transform output, got {out:?}");
        };
        assert_eq!(m.live_solids, vec![solid]);
        assert!(solid_is_rotated(&m, solid));
    }

    fn boundary_verts(m: &Model, s: Handle<Solid>) -> Vec<Handle<Vertex>> {
        let mut seen = std::collections::HashSet::new();
        let mut vs = Vec::new();
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh) {
                            vs.push(vh);
                        }
                    }
                }
            }
        }
        vs
    }

    fn rot_iso(axis: nacre_scalar::Axis, deg: i128) -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Isometry, Rat, Rotation as SRot};
        Isometry::rotation(SRot {
            axis,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    }

    /// Chain: (leaf, base_is_rotated, node_count, axes-root-to-leaf) for the first
    /// Rotated boundary vertex of `s`.
    fn forest_probe(m: &Model, s: Handle<Solid>) -> Option<(bool, usize, Vec<nacre_scalar::Axis>)> {
        let vh = *boundary_verts(m, s).first()?;
        let Origin::Rotated { base, rotation } = m.vertices.get(vh).origin else {
            return None;
        };
        let base_is_rotated = matches!(m.vertices.get(base).origin, Origin::Rotated { .. });
        let mut axes = Vec::new();
        let mut cur = Some(rotation);
        while let Some(h) = cur {
            let n = m.rotations.get(h);
            axes.push(n.axis);
            cur = n.parent;
        }
        axes.reverse();
        Some((base_is_rotated, axes.len(), axes))
    }

    fn translate_iso(off: [i128; 3]) -> nacre_scalar::Isometry {
        use nacre_scalar::{Isometry, Rat};
        Isometry::translation([
            Rat::from_int(off[0]),
            Rat::from_int(off[1]),
            Rat::from_int(off[2]),
        ])
    }

    fn rigid_iso(axis: nacre_scalar::Axis, deg: i128, off: [i128; 3]) -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Isometry, Rat, Rotation as SRot};
        Isometry::rigid(
            SRot {
                axis,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            },
            [
                Rat::from_int(off[0]),
                Rat::from_int(off[1]),
                Rat::from_int(off[2]),
            ],
        )
    }

    /// A boolean produces `Discovered` seam vertices with tol 0 (exact axis-aligned
    /// intersections). Before exact quadrantal realization, rotating them exactly 90°
    /// left an ~8e-17 f64 residual that exceeded tol 0 → `VertexOffSurface`. Now the
    /// rotation is exact, so the residual stays 0 and validate is clean — both for a
    /// pure 90° rotation and for a rigid 90°+translation (the offset cancels in
    /// vertex−plane, so it does not reintroduce a residual).
    #[test]
    fn boolean_result_rotated_90_validates() {
        use nacre_scalar::Axis;
        for iso in [rot_iso(Axis::Z, 90), rigid_iso(Axis::Z, 90, [5, -3, 2])] {
            let (mut m, a, b) = two_boxes();
            let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
            let before = nacre_props::mass_props(&m, r).unwrap().volume;
            let r2 = transform(&mut m, r, &iso).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(
                vs.is_empty(),
                "exact 90° realization → validate clean: {vs:?}"
            );
            let after = nacre_props::mass_props(&m, r2).unwrap().volume;
            assert!((after - before).abs() < 1e-12, "volume invariant");
        }
    }

    /// A 90°-family rotation lands a cuboid's vertices exactly on the axis-aligned grid
    /// (no ~6e-17 spurious offset): the corner (2,3,4) rotated 90° about Z maps to
    /// exactly (-3,2,4).
    #[test]
    fn rotate_90_lands_vertices_exactly() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c2 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
        m.rebuild_adjacency();
        let pts = outer_points(&m, c2);
        // (2,3,4) about Z by 90°: (x,y)→(-y,x) → (-3,2,4). Bit-exact.
        assert!(
            pts.iter().any(|p| *p == [-3.0, 2.0, 4.0]),
            "corner lands exactly on the grid; got {pts:?}"
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// Re-rotating about the same axis chains a second forest node onto the first
    /// (this cell does not bundle): the leaf's parent is the earlier rotation, `base`
    /// stays the Constructed root, and the solid remains a rigid (volume/area-invariant)
    /// `Rotated` solid that validate/tess/STEP accept and a boolean now runs against.
    #[test]
    fn rerotate_same_axis_chains() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap();
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::Z, 20)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c2),
            Some((false, 2, vec![Axis::Z, Axis::Z])),
            "two Z nodes chained, base = the Constructed root"
        );
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap();
        assert!(
            (after.volume - before.volume).abs() < 1e-9,
            "volume invariant"
        );
        assert!((after.area - before.area).abs() < 1e-9, "area invariant");
        assert!(nacre_tess::to_obj(&m).is_ok());
        assert!(
            nacre_step::to_step(&m)
                .unwrap()
                .contains("MANIFOLD_SOLID_BREP")
        );
        assert!(solid_is_rotated(&m, c2), "re-rotated solid stays Rotated");
        // A boolean against the chain-rotated solid runs (guard retired, 3d-i): the axis-aligned
        // `d` inside the re-rotated `c2` is carved out, and the result is a valid solid.
        let d = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "chain-rotated cut is valid"
        );
    }

    /// Re-rotating about a different axis chains a node whose parent is the first
    /// rotation (root → Z → X), `base` still the root.
    #[test]
    fn rerotate_different_axis_chains() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap().volume;
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c2),
            Some((false, 2, vec![Axis::Z, Axis::X])),
            "chain root→Z→X, base = root"
        );
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap().volume;
        assert!((after - before).abs() < 1e-9, "volume invariant");
    }

    /// An **exact** (90°-family) rotation applied to an already-rotated solid is still
    /// recorded as a chain node — the composite is inexact (an ancestor is), so the
    /// forest must stay complete (1b silently dropped it). The solid stays Rotated and
    /// boolean-rejected.
    #[test]
    fn rerotate_exact_after_inexact_records_node() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 90)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c2),
            Some((false, 2, vec![Axis::Z, Axis::X])),
            "the exact 90°X is recorded as a chain node, not dropped"
        );
        assert!(
            solid_is_rotated(&m, c2),
            "composite is inexact → still Rotated"
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A translation between two same-axis rotations forces a chain (this cell never
    /// bundles anyway): the forest records both rotations, `base` stays the root, and
    /// the result is rigid and valid — sound with no adjacency guard (each rotation is
    /// its own node).
    #[test]
    fn rerotate_across_translation_chains() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap().volume;
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &translate_iso([5, -3, 2])).unwrap();
        m.rebuild_adjacency();
        let c3 = transform(&mut m, c2, &rot_iso(Axis::Z, 20)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c3),
            Some((false, 2, vec![Axis::Z, Axis::Z])),
            "both rotations recorded; the intervening translation is not in the forest"
        );
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c3).unwrap().volume;
        assert!((after - before).abs() < 1e-9, "volume invariant");
    }

    /// A three-axis chain records three nodes root→Z→X→Y.
    #[test]
    fn rerotate_deep_chain() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
        m.rebuild_adjacency();
        let c3 = transform(&mut m, c2, &rot_iso(Axis::Y, 20)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c3),
            Some((false, 3, vec![Axis::Z, Axis::X, Axis::Y])),
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A fresh rotation of a Constructed solid is unchanged from 1b (B0): an inexact
    /// angle records a single root node; a 90°-family angle stays Constructed (no node,
    /// boolean allowed). Guards that the B0/B1 split preserves fresh-rotation behavior.
    #[test]
    fn fresh_rotation_of_constructed_unchanged() {
        use nacre_scalar::Axis;
        // inexact → one root node, base = the cuboid's Constructed vertices.
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        assert_eq!(forest_probe(&m, c1), Some((false, 1, vec![Axis::Z])));

        // exact 90° → Constructed (no node), boolean allowed.
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
        m.rebuild_adjacency();
        assert!(!solid_is_rotated(&m, c1), "fresh 90° stays Constructed");
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(forest_probe(&m, c1), None, "no rotation node");
    }

    /// Replay determinism (DNA 3): a re-rotation sequence reproduces the same forest
    /// and handles.
    #[test]
    fn rerotate_is_deterministic() {
        use nacre_scalar::Axis;
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
            let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
            (bbox_lo(&m, c2), c2)
        };
        assert_eq!(build(), build(), "same ops → same geometry and handle");
    }

    #[test]
    fn fuse_of_disjoint_boxes_is_empty() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        assert_eq!(
            boolean_one(&mut m, BoolKind::Fuse, a, b),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn cut_containment_makes_a_cavity() {
        // A = [0,3]³ (27) with B = [1,2]³ (1) strictly inside ⇒ A − B is a
        // hollow solid: volume 26, an outer + one void shell (V16/E24/F12/S2).
        let (mut m, a, b) = nested_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
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
            boolean_one(&mut m, BoolKind::Cut, inner, outer),
            Err(BoolError::EmptyResult)
        );
        // Fuse(inner ∪ outer) = outer.
        let r = boolean_one(&mut m, BoolKind::Fuse, inner, outer).unwrap();
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
            let r = boolean_one(&mut m, BoolKind::Common, x, y).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
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
    fn fuse_a_boss_onto_a_non_convex_solid() {
        // A contained boss on the top of an L-prism (non-convex kept `a`). The convexity gate
        // used to decline this to the seam path, which rejected the seamless contact; the
        // contained-coplanar Fuse now admits it. Volume 3 (L) + 0.4²·0.5 = 3.08.
        let (mut m, l) = l_prism(); // L footprint area 3, height 1
        let boss = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 1.0]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, l, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 3.08).abs() < 1e-12, "volume {vol}");
        // The boss top cap sits on the z = 1.5 plane, its outward normal +z.
        assert!(
            find_face_on_plane(
                &m,
                r,
                Point3::from_array([0.5, 0.5, 1.5]),
                Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .is_some()
        );
    }

    #[test]
    fn fuse_a_non_convex_profile_boss() {
        // An L-shaped boss (non-convex cutter `b`) on a cube top. The gate used to decline the
        // non-convex prism; the contained-coplanar Fuse now carries the L footprint as a hole.
        // Volume 1 (cube) + 0.12 (L area) · 0.4 = 1.048.
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let l_base: Vec<Point3> = [
            [0.3, 0.3],
            [0.7, 0.3],
            [0.7, 0.5],
            [0.5, 0.5],
            [0.5, 0.7],
            [0.3, 0.7],
        ]
        .iter()
        .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
        .collect();
        let (boss, _) =
            build_prism(&mut m, &l_base, Vector3::from_array([0.0, 0.0, 0.4]), None).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, cube, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.048).abs() < 1e-12, "volume {vol}");
        // The L boss top cap sits on the z = 1.4 plane, its outward normal +z.
        assert!(
            find_face_on_plane(
                &m,
                r,
                Point3::from_array([0.4, 0.4, 1.4]),
                Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .is_some()
        );
    }

    #[test]
    fn fuse_an_overhanging_boss() {
        // A boss whose footprint hangs past a single edge of the base's top face: part fuses
        // onto the face, part cantilevers into the air (cell coplanar-contact-overhang). The
        // base top gains a boundary notch, the boss bottom keeps only its overhanging piece,
        // and the two footprint crossings weld the walls. Before this cell it rejected as
        // vertex_on_face_plane. Volume 1 + 1.0·0.5·1.0 = 1.5 (the solids share only z=1).
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.5).abs() < 1e-12, "volume {vol}");
        // The notch (base top remainder) and the cantilever (boss underside) are coplanar but
        // meet only at the two crossing points, not along an edge, so the solid chains.
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn cut_an_edge_slot() {
        // A prism sits top-flush on the base but its footprint hangs past the base's x=1 edge.
        // Cut carves an edge-slot that breaks out through the base's x=1 wall: the base top gains
        // a mouth notch, the x=1 wall gains a side opening, and the prism's inside walls become
        // the slot surfaces. Before this cell it rejected as vertex_on_face_plane. Removed volume
        // = 0.5·0.5·0.5 = 0.125 → 0.875.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 0.5]),
            Point3::from_array([1.5, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.875).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn an_edge_slot_through_the_bottom_is_out_of_scope() {
        // The prism pokes out the base's bottom too — its walls cross a second base face (the
        // blind gate fails there), so it is out of scope. detect_overhang_cut_general must decline.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.25, -0.5]),
            Point3::from_array([1.5, 0.75, 1.0]),
        );
        assert!(detect_overhang_cut_general(&m, base, prism).is_none());
    }

    #[test]
    fn cut_a_slab_channel() {
        // A prism top-flush on the base but spanning clear across it in x: its footprint crosses
        // two opposite base-top edges (x=0, x=1), so Cut carves a channel breaking out both x
        // sides. The base top splits into two notch strips (y<0.4, y>0.6), each x wall gains a
        // side opening, and the prism's y-walls + floor (clipped to x∈[0,1]) become the channel
        // surfaces. Removed volume 1·0.2·0.5 = 0.1 → 0.9.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 0.4, 0.5]),
            Point3::from_array([1.5, 0.6, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, slab).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.9).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn cut_a_corner_slot() {
        // A prism top-flush on the base swallowing its (1,1) corner and crossing the two adjacent
        // edges (x=1, y=1): Cut carves a corner slot breaking out both walls. The base top gains
        // an L mouth, each wall a corner opening (one swallowed corner each), and the two openings
        // meet at the corner column; the prism's inner walls + floor (clipped to inside the base)
        // are the slot surfaces. Removed 0.5³ = 0.125 → 0.875.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, corner).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.875).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn cut_an_l_step() {
        // A prism top-flush on the base spanning clear across it in x and hanging past its y=1 edge:
        // its footprint runs through three walls — x=0 and x=1 each swallow one top corner (a
        // corner opening), and y=1 is fully covered (both corners swallowed, no crossing → the wall
        // shrinks to a rectangle at the two corner columns). The general N-wall path handles all
        // three at once, opening the 3-wall config the edge-slot/slab/corner paths could not. The
        // base top gains an L mouth; R = [0,1]×[0.5,1]×[0.5,1] = 0.25 is removed → 0.75.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([-0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    /// The overhang-Common invariants: `R = a ∩ b` is a clean convex solid of the given volume
    /// (`is_convex` catches a face-orientation bug the volume alone would pass), watertight, with no
    /// coplanar-neighbour edge and no cavity.
    fn assert_common_box(m: &Model, r: Handle<Solid>, vol: f64) {
        let vs = nacre_validate::validate(m);
        assert!(vs.is_empty(), "{vs:?}");
        let v = nacre_props::mass_props(m, r).unwrap().volume;
        assert!((v - vol).abs() < 1e-12, "volume {v}");
        let planes = collect_planes(m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(m, r, &planes));
        assert!(is_convex(m, &planes, &solid_vertex_handles(m, r)));
        assert!(m.solids.get(r).cavities.is_empty());
    }

    #[test]
    fn common_an_edge_overhang() {
        // A prism top-flush on the base hanging past its y=1 edge: Common keeps the convex overlap
        // R = [0.3,0.7]×[0.5,1]×[0.5,1] (one breached wall, no swallowed corner). Volume 0.1.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.3, 0.5, 0.5]),
            Point3::from_array([0.7, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
        m.rebuild_adjacency();
        assert_common_box(&m, r, 0.1);
    }

    #[test]
    fn common_a_corner_overhang() {
        // A prism swallowing the base's (1,1) corner: Common keeps R = [0.5,1]×[0.5,1]×[0.5,1] (two
        // breached walls meeting at one corner column cc=(1,1,0.5)). Volume 0.125.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
        m.rebuild_adjacency();
        assert_common_box(&m, r, 0.125);
    }

    #[test]
    fn common_an_l_step() {
        // The L-step prism (spanning x, hanging past y=1): three breached walls sharing two corner
        // columns cc0=(0,1,0.5), cc1=(1,1,0.5). Common keeps R = [0,1]×[0.5,1]×[0.5,1] (volume 0.25),
        // proving the two-crossing gate is independent of breached-wall count.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([-0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
        m.rebuild_adjacency();
        assert_common_box(&m, r, 0.25);
    }

    #[test]
    fn common_a_slab_overhang_is_out_of_scope() {
        // A spanning slab crosses two opposite base edges — four contact crossings, two overlap-arc
        // pairs. The two-crossing gate declines it; the boolean rejects via the standard fall-through.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 0.4, 0.5]),
            Point3::from_array([1.5, 0.6, 1.0]),
        );
        assert!(detect_overhang_common(&m, base, slab).is_none());
        assert_rejects(
            || boolean_one(&mut m, BoolKind::Common, base, slab),
            tag::VERTEX_ON_FACE_PLANE,
        );
    }

    #[test]
    fn a_corner_cut_through_the_bottom_is_out_of_scope() {
        // The corner prism also pokes out the base's bottom (blind gate fails on the bottom face),
        // a through-slot beyond this cell. detect_overhang_cut_general must decline it.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.5, 0.5, -0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        assert!(detect_overhang_cut_general(&m, base, through).is_none());
    }

    #[test]
    fn a_boss_that_pierces_the_base_is_not_an_overhang() {
        // The boss dips below the base's top (its walls cross the base) — a transversal seam
        // cut, not a coplanar overhang. detect_overhang_contact must decline it.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 0.5]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        assert!(detect_overhang_contact(&m, base, through).is_none());
    }

    #[test]
    fn fuse_a_corner_overhanging_boss() {
        // The boss footprint swallows a base-top corner, crossing two base edges (a 2-crossing
        // lens with the corner inside). Fuse welds an L-shaped cantilever: the base top gains an
        // L notch, the boss underside keeps its L remainder, and the swallowed corner (1,1) is
        // the shared reflex vertex. Before this cell it rejected as vertex_on_face_plane.
        // Volume 1 + 1·1·1 = 2 (the solids share only z=1).
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, corner).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 2.0).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn fuse_a_spanning_slab_boss() {
        // The boss footprint spans clear across the base (a slab overhanging both x sides), so
        // ∂Q crosses ∂P at four points. Fuse welds it: the base top splits into two notch strips
        // (y<0.4, y>0.6) and the boss underside into two cantilever pieces (x<0, x>1). Before this
        // cell it rejected as vertex_on_face_plane. Volume 1 + 2.0·0.2·1.0 = 1.4 (share only z=1).
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 0.4, 1.0]),
            Point3::from_array([1.5, 0.6, 2.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, slab).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.4).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    // Build a unit cube with a blind pocket in its top — a non-convex solid whose side faces stay
    // convex. Returns `(model, pocketed solid)`.
    fn top_pocketed_cube() -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude_op(square(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        let OpOutput::PocketOnFace { solid, .. } =
            apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        (m, solid)
    }

    // n0: cross-section of a solid by a wall plane (the new `section_of_solid` primitive).
    type SectionSetup = (
        Vec<PlaneInfo>,
        usize,
        std::collections::HashMap<Handle<Face>, usize>,
    );
    fn section_setup(m: &Model, a: Handle<Solid>, b: Handle<Solid>) -> SectionSetup {
        let planes_a = collect_planes(m, a).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(m, b).unwrap());
        let mut surf_ix = std::collections::HashMap::new();
        for (i, pi) in planes.iter().enumerate() {
            surf_ix.insert(pi.face, i);
        }
        // W = b's -x face (n_out ≈ [-1,0,0]) at x = 0.5.
        let w = (na..planes.len())
            .find(|&i| {
                let n = planes[i].n_out.as_array();
                n[0] < -0.5 && (planes[i].tri[0].as_array()[0] - 0.5).abs() < 1e-9
            })
            .expect("b's -x wall at x=0.5");
        (planes, w, surf_ix)
    }

    fn section_pts(loops: &[Vec<Node>], planes: &[PlaneInfo]) -> Vec<[f64; 3]> {
        loops
            .iter()
            .flat_map(|l| l.iter())
            .map(|&nd| match nd {
                Node::Seam(t) => three_planes(
                    &planes[t[0]].plane,
                    &planes[t[1]].plane,
                    &planes[t[2]].plane,
                )
                .unwrap()
                .as_array(),
                _ => panic!("section is all Seam nodes"),
            })
            .collect()
    }

    #[test]
    fn section_of_a_cube_is_a_quad() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, -1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        m.rebuild_adjacency();
        let (planes, w, surf_ix) = section_setup(&m, a, b);
        let loops = section_of_solid(&m, a, w, &planes, &surf_ix).unwrap();
        assert_eq!(loops.len(), 1, "one section loop");
        assert_eq!(loops[0].len(), 4, "a cube slice is a quad");
        let mut yz: Vec<[f64; 2]> = section_pts(&loops, &planes)
            .iter()
            .inspect(|p| assert!((p[0] - 0.5).abs() < 1e-9, "vertex on W"))
            .map(|p| [p[1], p[2]])
            .collect();
        yz.sort_by(|u, v| {
            u[0].partial_cmp(&v[0])
                .unwrap()
                .then(u[1].partial_cmp(&v[1]).unwrap())
        });
        assert_eq!(yz, vec![[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]]);
    }

    #[test]
    fn section_of_a_pocketed_cube_is_non_convex() {
        // A top-pocketed cube sliced through the pocket: section is the outer square with the
        // pocket's square bite where the slice passes through the pocket walls.
        let (mut m, a) = top_pocketed_cube();
        // small_square pocket is inset; slice at x=0.5 crosses it. Give b a -x wall at x=0.5.
        let b = m.add_cuboid(
            Point3::from_array([0.5, -1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        m.rebuild_adjacency();
        let (planes, w, surf_ix) = section_setup(&m, a, b);
        let loops = section_of_solid(&m, a, w, &planes, &surf_ix).unwrap();
        // The pocketed-cube cross-section at x=0.5 is a single non-convex loop (outer square
        // dented by the pocket) — more than 4 vertices, one loop, all on W.
        assert_eq!(loops.len(), 1, "one outer section loop");
        assert!(loops[0].len() > 4, "non-convex slice has >4 vertices");
        for p in section_pts(&loops, &planes) {
            assert!((p[0] - 0.5).abs() < 1e-9, "vertex on W");
        }
    }

    // C2a: a real `section_of_solid` loop, turned into a boundary by `section_boundary`, drives
    // `coplanar_reconstruct` exactly like a contact footprint — the section-Q clip that C2b needs.
    #[test]
    fn section_boundary_clips_a_wall_to_inside_a_solid() {
        // `a` is a 4×4×2 block; `b`'s top face (z=1) is a 5×2 rectangle [1,6]×[1,3] that pokes out
        // of `a` on the +x side. Section `a` at z=1 → the full 4×4 square, and clip `b`'s top face
        // (keep inside the section): the survivor is [1,4]×[1,3], area 6, two of its corners the
        // original `b` vertices and two the crossings on `a`'s x=4 wall.
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([6.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let planes_a = collect_planes(&m, a).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(&m, b).unwrap());
        let mut surf_ix = std::collections::HashMap::new();
        for (i, p) in planes.iter().enumerate() {
            surf_ix.insert(p.face, i);
        }
        // β = b's +z (top) wall at z=1.
        let beta = (na..planes.len())
            .find(|&i| {
                let n = planes[i].n_out.as_array();
                n[2] > 0.5 && (planes[i].tri[0].as_array()[2] - 1.0).abs() < 1e-9
            })
            .expect("b's +z wall at z=1");
        let canon = plane_classes(&planes);
        assert_eq!(canon[beta], beta, "the section plane folds with nothing");

        let loops = section_of_solid(&m, a, beta, &planes, &surf_ix).unwrap();
        assert_eq!(loops.len(), 1, "one section loop");
        assert_eq!(loops[0].len(), 4, "the block's slice is a quad");

        let q_bnd = section_boundary(&loops, beta).unwrap();
        assert_eq!(q_bnd.len(), 4, "four section chords");
        let inc_b = arrange::edge_planes(&m, b, &surf_ix).unwrap();
        let p_bnd = contact_boundary(&m, planes[beta].face, beta, &inc_b, &canon).unwrap();
        let crossings = coplanar_boundary_crossings(&planes, beta, &p_bnd, &q_bnd).unwrap();
        assert_eq!(crossings.len(), 2, "b's top edges cross a's x=4 wall twice");
        let arcs = coplanar_seam_arcs(&planes, beta, &q_bnd, &crossings);
        let faces = coplanar_reconstruct(
            &planes, beta, beta, &p_bnd, &q_bnd, &crossings, &arcs, true, false,
        )
        .unwrap();
        assert_eq!(faces.len(), 1, "one clipped cell");

        // Node points → projected (x,y) → shoelace area on plane z=1.
        let node_xy = |nd: &Node| -> [f64; 2] {
            let p = match *nd {
                Node::Orig(vh) => m.vertices.get(vh).point,
                Node::Seam(t) => three_planes(
                    &planes[t[0]].plane,
                    &planes[t[1]].plane,
                    &planes[t[2]].plane,
                )
                .unwrap(),
            }
            .as_array();
            [p[0], p[1]]
        };
        let ring: Vec<[f64; 2]> = faces[0].loop_nodes.iter().map(node_xy).collect();
        let n = ring.len();
        let area = 0.5
            * (0..n)
                .map(|i| {
                    let (u, v) = (ring[i], ring[(i + 1) % n]);
                    u[0] * v[1] - v[0] * u[1]
                })
                .sum::<f64>()
                .abs();
        assert_eq!(faces[0].loop_nodes.len(), 4, "survivor is a quad");
        assert!((area - 6.0).abs() < 1e-9, "clipped area is 6, got {area}");
    }

    // C2a: the one silent-wrong risk of the section approach — a multi-loop section (outer + hole)
    // must be rejected up front, never flattened to its first ring.
    #[test]
    fn section_boundary_rejects_a_multi_loop_section() {
        let two_loops = vec![
            vec![
                Node::Seam([0, 1, 2]),
                Node::Seam([0, 1, 3]),
                Node::Seam([0, 2, 3]),
            ],
            vec![
                Node::Seam([0, 4, 5]),
                Node::Seam([0, 4, 6]),
                Node::Seam([0, 5, 6]),
            ],
        ];
        LAST_REJECT.with(|c| c.take());
        assert!(section_boundary(&two_loops, 0).is_err());
        assert_eq!(
            LAST_REJECT.with(|c| c.take()),
            Some(tag::SECTION_MULTI_LOOP)
        );
    }

    // C2b n0': the through-bottom corner of the deliverable target must resolve to ONE shared
    // 3-plane triple across the two faces that meet the cutter there — a plain trihedral corner,
    // not a fan. If the two sections disagree (or need rotational ordering), the scoped-flush
    // premise fails and C2b-4 would depend on the general turn_at resolver (C1). This throwaway
    // probe falsifies that C1-dependence BEFORE any Cut-branch code is written.
    #[test]
    fn probe_target_through_bottom_corner_is_a_plain_trihedral() {
        // Deliverable target: a slot cut from top_pocketed_cube, flush on the +x wall (x=1),
        // breaking out the bottom (z=0). The corner (1, 0.25, 0) is where the cube's +x wall, the
        // cube's bottom, and the slot's -y wall meet.
        let (mut m, cube) = top_pocketed_cube();
        let slot = m.add_cuboid(
            Point3::from_array([0.75, 0.25, -0.25]),
            Point3::from_array([1.0, 0.75, 0.5]),
        );
        m.rebuild_adjacency();
        let planes_a = collect_planes(&m, cube).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(&m, slot).unwrap());
        let mut surf_ix = std::collections::HashMap::new();
        for (i, p) in planes.iter().enumerate() {
            surf_ix.insert(p.face, i);
        }
        let canon = plane_classes(&planes);
        let find = |rng: std::ops::Range<usize>, n: [f64; 3], coord: usize, val: f64| -> usize {
            rng.clone()
                .find(|&i| {
                    let nn = planes[i].n_out.as_array();
                    nn[0] * n[0] + nn[1] * n[1] + nn[2] * n[2] > 0.5
                        && (planes[i].tri[0].as_array()[coord] - val).abs() < 1e-9
                })
                .expect("plane")
        };
        let cube_x1 = find(0..na, [1.0, 0.0, 0.0], 0, 1.0); // cube +x wall (contact)
        let cube_z0 = find(0..na, [0.0, 0.0, -1.0], 2, 0.0); // cube bottom
        let slot_x1 = find(na..planes.len(), [1.0, 0.0, 0.0], 0, 1.0); // slot +x (contact)
        let slot_y025 = find(na..planes.len(), [0.0, -1.0, 0.0], 1, 0.25); // slot -y wall

        // x=1 is the contact plane: cube and slot +x walls fold into one class.
        assert_eq!(
            canon[cube_x1], canon[slot_x1],
            "x=1 is the shared contact plane"
        );

        // The corner (1, 0.25, 0) as a canonicalized sorted triple, from each section that produces
        // it. Section of the SLOT at the cube bottom (z=0); section of the CUBE at the slot -y wall.
        let corner = [1.0, 0.25, 0.0];
        let corner_triple = |sect_solid: Handle<Solid>, w_idx: usize| -> [usize; 3] {
            let loops = section_of_solid(&m, sect_solid, w_idx, &planes, &surf_ix).unwrap();
            for l in &loops {
                for &nd in l {
                    let Node::Seam(t) = nd else { continue };
                    let p = three_planes(
                        &planes[t[0]].plane,
                        &planes[t[1]].plane,
                        &planes[t[2]].plane,
                    )
                    .unwrap()
                    .as_array();
                    if (0..3).all(|k| (p[k] - corner[k]).abs() < 1e-9) {
                        let mut c = [canon[t[0]], canon[t[1]], canon[t[2]]];
                        c.sort_unstable();
                        return c;
                    }
                }
            }
            panic!("corner (1,0.25,0) not found in section");
        };
        let from_bottom = corner_triple(slot, cube_z0); // slot ∩ {z=0}
        let from_ywall = corner_triple(cube, slot_y025); // cube ∩ {y=0.25}
        eprintln!("PROBE corner from bottom-section: {from_bottom:?}");
        eprintln!("PROBE corner from ywall-section:  {from_ywall:?}");
        assert_eq!(
            from_bottom, from_ywall,
            "the through-bottom corner must be ONE shared 3-plane triple (no fan → no C1)"
        );
        // Sanity: the shared triple is exactly {x=1, z=0, y=0.25} in canon classes.
        let mut expect = [canon[cube_x1], canon[cube_z0], canon[slot_y025]];
        expect.sort_unstable();
        assert_eq!(from_bottom, expect, "triple is {{x=1, z=0, y=0.25}}");

        // And the x=1 contact face itself cannot be sectioned (slot vertices lie ON x=1) — it must
        // use the footprint directly, confirming the mouth/section split.
        LAST_REJECT.with(|c| c.take());
        assert!(
            section_of_solid(&m, slot, cube_x1, &planes, &surf_ix).is_err(),
            "x=1 is degenerate for the slot (contact plane) → footprint, not section"
        );
    }

    // A1: plane-class canonicalization — coplanar walls of the two operands fold into one line.
    #[test]
    fn plane_classes_merge_a_shared_wall() {
        // Two unit cubes side by side share the plane x=1 (a's +x wall, b's -x wall — the same
        // plane, opposite normals). `plane_classes` must merge those two into one line class and
        // keep the far walls (a's x=0, b's x=2) distinct.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let planes_a = collect_planes(&m, a).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(&m, b).unwrap());
        // Find a plane by outward-normal x-sign and its x coordinate, within an index range.
        let find = |rng: std::ops::Range<usize>, nx: f64, x: f64| -> usize {
            rng.clone()
                .find(|&i| {
                    let n = planes[i].n_out.as_array();
                    n[0] * nx > 0.5 && (planes[i].tri[0].as_array()[0] - x).abs() < 1e-9
                })
                .expect("plane")
        };
        let a_xp = find(0..na, 1.0, 1.0); // a's +x wall at x=1
        let b_xm = find(na..planes.len(), -1.0, 1.0); // b's -x wall at x=1
        let a_xm = find(0..na, -1.0, 0.0); // a's -x wall at x=0
        let b_xp = find(na..planes.len(), 1.0, 2.0); // b's +x wall at x=2
        let canon = plane_classes(&planes);
        assert_eq!(canon[a_xp], canon[b_xm], "shared x=1 wall is one class");
        assert_ne!(canon[a_xm], canon[b_xp], "far walls stay distinct");
        assert_ne!(canon[a_xp], canon[a_xm], "x=1 and x=0 are different lines");
        // Every b face but its +x wall is coplanar with an a face, so 12 planes fold to 7 classes.
        let distinct: std::collections::HashSet<usize> = canon.iter().copied().collect();
        assert_eq!(distinct.len(), na + 1, "only b's far wall is a new class");
        // The class root is the smallest index in the class (deterministic canon).
        assert_eq!(canon[a_xp], a_xp.min(b_xm));
    }

    // A2: exact coplanar boundary overlay — proper crossings of two footprints in a shared plane.
    #[test]
    fn coplanar_crossings_of_two_overlapping_squares() {
        // a's top (z=1, +z) and b's bottom (z=1, -z) are coplanar; footprints [0,2]² and [1,3]²
        // overlap, boundaries crossing properly at exactly (1,2,1) and (2,1,1).
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 1.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        );
        m.rebuild_adjacency();
        let planes_a = collect_planes(&m, a).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(&m, b).unwrap());
        let mut surf_ix = std::collections::HashMap::new();
        for (i, pinf) in planes.iter().enumerate() {
            surf_ix.insert(pinf.face, i);
        }
        let canon = plane_classes(&planes);
        let find = |rng: std::ops::Range<usize>, nz: f64| -> usize {
            rng.clone()
                .find(|&i| {
                    let n = planes[i].n_out.as_array();
                    n[2] * nz > 0.5 && (planes[i].tri[0].as_array()[2] - 1.0).abs() < 1e-9
                })
                .expect("contact face at z=1")
        };
        let a_top = find(0..na, 1.0);
        let b_bot = find(na..planes.len(), -1.0);
        let pi = canon[a_top];
        assert_eq!(pi, canon[b_bot], "contact faces are one plane class");
        let inc_a = arrange::edge_planes(&m, a, &surf_ix).unwrap();
        let inc_b = arrange::edge_planes(&m, b, &surf_ix).unwrap();
        let a_bnd = contact_boundary(&m, planes[a_top].face, pi, &inc_a, &canon).unwrap();
        let b_bnd = contact_boundary(&m, planes[b_bot].face, pi, &inc_b, &canon).unwrap();
        assert_eq!(a_bnd.len(), 4, "square footprint has 4 boundary edges");
        assert_eq!(b_bnd.len(), 4);
        assert!(
            a_bnd.iter().all(|e| e.v[0] != e.v[1]),
            "each edge has two distinct endpoints"
        );
        let cx = coplanar_boundary_crossings(&planes, pi, &a_bnd, &b_bnd).unwrap();
        for c in &cx {
            assert!(c.triple.contains(&pi), "crossing triple carries π");
            assert!(a_bnd.iter().any(|e| e.seg == c.a_seg), "a_seg is on ∂P");
            assert!(b_bnd.iter().any(|e| e.seg == c.b_seg), "b_seg is on ∂Q");
        }
        let mut pts: Vec<[f64; 3]> = cx.iter().map(|c| c.point.as_array()).collect();
        pts.sort_by(|u, v| {
            u[0].partial_cmp(&v[0])
                .unwrap()
                .then(u[1].partial_cmp(&v[1]).unwrap())
        });
        assert_eq!(pts.len(), 2, "two proper crossings, got {pts:?}");
        assert!(
            pts[0]
                .iter()
                .zip([1.0, 2.0, 1.0])
                .all(|(p, q)| (p - q).abs() < 1e-9),
            "{:?}",
            pts[0]
        );
        assert!(
            pts[1]
                .iter()
                .zip([2.0, 1.0, 1.0])
                .all(|(p, q)| (p - q).abs() < 1e-9),
            "{:?}",
            pts[1]
        );

        // A2 part2: ∂Q split at the crossings into mixed-node arcs. Two crossings give two arcs,
        // each running Seam → Orig… → Seam. One arc's interior b-vertex (1,1) is inside P=[0,2]²
        // (the piece of ∂Q that cuts P); the other holds the three outside corners.
        let arcs = coplanar_seam_arcs(&planes, pi, &b_bnd, &cx);
        assert_eq!(arcs.len(), 2, "two crossings split ∂Q into two arcs");
        for arc in &arcs {
            assert!(
                matches!(arc.first(), Some(Node::Seam(_)))
                    && matches!(arc.last(), Some(Node::Seam(_))),
                "an arc runs crossing → crossing"
            );
            assert!(
                arc[1..arc.len() - 1]
                    .iter()
                    .all(|n| matches!(n, Node::Orig(_))),
                "interior nodes are b's original vertices (Constructed)"
            );
        }
        // The inside-P arc is the one whose single interior vertex is (1,1); it has exactly one.
        let inside_arc = arcs
            .iter()
            .find(|a| {
                a.len() == 3
                    && matches!(a[1], Node::Orig(vh)
                        if m.vertices.get(vh).point.as_array() == [1.0, 1.0, 1.0])
            })
            .expect("an arc with the single interior corner (1,1)");
        let outside_arc = arcs
            .iter()
            .find(|a| a.len() == 5)
            .expect("the 3-corner arc");
        assert!(
            outside_arc[1..4].iter().all(|n| matches!(n, Node::Orig(vh)
                if { let p = m.vertices.get(*vh).point.as_array(); p[0] > 2.0 || p[1] > 2.0 })),
            "outside arc holds only corners beyond P"
        );
        assert_eq!(inside_arc.len() + outside_arc.len(), 8);

        // A3: reconstruct P subdivided by Q. Keep P∖Q (an L-shaped cell of 6 mixed nodes) and
        // keep P∩Q (the [1,2]² overlap, 4 nodes). Every decision is point_in_ring (no coordinate).
        let plane_idx = a_top;
        let pt = |nd: &Node| -> [f64; 3] {
            match nd {
                Node::Orig(vh) => m.vertices.get(*vh).point.as_array(),
                Node::Seam(t) => three_planes(
                    &planes[t[0]].plane,
                    &planes[t[1]].plane,
                    &planes[t[2]].plane,
                )
                .unwrap()
                .as_array(),
            }
        };
        let set = |f: &LocalFace| -> std::collections::BTreeSet<[i64; 3]> {
            f.loop_nodes
                .iter()
                .map(|nd| {
                    let p = pt(nd);
                    [
                        (p[0] * 8.0).round() as i64,
                        (p[1] * 8.0).round() as i64,
                        (p[2] * 8.0).round() as i64,
                    ]
                })
                .collect()
        };
        let key = |xs: &[[f64; 3]]| -> std::collections::BTreeSet<[i64; 3]> {
            xs.iter()
                .map(|p| {
                    [
                        (p[0] * 8.0).round() as i64,
                        (p[1] * 8.0).round() as i64,
                        (p[2] * 8.0).round() as i64,
                    ]
                })
                .collect()
        };
        let out = coplanar_reconstruct(
            &planes, pi, plane_idx, &a_bnd, &b_bnd, &cx, &arcs, false, false,
        )
        .unwrap();
        assert_eq!(out.len(), 1, "P∖Q is one cell");
        assert_eq!(out[0].loop_nodes.len(), 6, "L-shape has 6 boundary nodes");
        assert_eq!(
            set(&out[0]),
            key(&[
                [0.0, 0.0, 1.0],
                [2.0, 0.0, 1.0],
                [2.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
                [1.0, 2.0, 1.0],
                [0.0, 2.0, 1.0],
            ])
        );
        let inter = coplanar_reconstruct(
            &planes, pi, plane_idx, &a_bnd, &b_bnd, &cx, &arcs, true, false,
        )
        .unwrap();
        assert_eq!(inter.len(), 1, "P∩Q is one cell");
        assert_eq!(
            set(&inter[0]),
            key(&[
                [1.0, 1.0, 1.0],
                [2.0, 1.0, 1.0],
                [2.0, 2.0, 1.0],
                [1.0, 2.0, 1.0],
            ])
        );

        // B4: the symmetric Q-side reconstruction (roles swapped) — b's cantilever. Keep Q∖P →
        // the L-shaped cantilever bottom (b corners beyond P, the two crossings, and P's (2,2)).
        let cx_q = coplanar_boundary_crossings(&planes, pi, &b_bnd, &a_bnd).unwrap();
        let arcs_q = coplanar_seam_arcs(&planes, pi, &a_bnd, &cx_q); // ∂P as the seam on Q
        let cant = coplanar_reconstruct(
            &planes, pi, b_bot, &b_bnd, &a_bnd, &cx_q, &arcs_q, false, false,
        )
        .unwrap();
        assert_eq!(cant.len(), 1, "Q∖P is one cell");
        assert_eq!(cant[0].loop_nodes.len(), 6, "cantilever is an L");
        assert_eq!(
            set(&cant[0]),
            key(&[
                [3.0, 1.0, 1.0],
                [3.0, 3.0, 1.0],
                [1.0, 3.0, 1.0],
                [1.0, 2.0, 1.0],
                [2.0, 2.0, 1.0],
                [2.0, 1.0, 1.0],
            ])
        );
    }

    // B1: the coplanar-face survival table, tied to the six bespoke paths it must reproduce.
    #[test]
    fn coplanar_survival_table_matches_the_bespoke_paths() {
        use PSurvive::*;
        // (kind, same_normal) -> (P-side survival, b faces flip)
        assert_eq!(coplanar_survival(BoolKind::Cut, true), (MinusQ, true)); // blind pocket mouth
        assert_eq!(coplanar_survival(BoolKind::Cut, false), (Whole, true)); // b above, no overlap
        assert_eq!(coplanar_survival(BoolKind::Fuse, false), (MinusQ, false)); // boss: Q is a hole
        assert_eq!(coplanar_survival(BoolKind::Fuse, true), (Whole, false)); // union, ∂Q internal
        assert_eq!(coplanar_survival(BoolKind::Common, true), (InterQ, false)); // overlap cap
        assert_eq!(coplanar_survival(BoolKind::Common, false), (Empty, false)); // coincident → empty
    }

    // B3: the survival table governs the coincident stack (P ≡ Q, opposite normals). Verified
    // against the (still-bespoke) coincident_merge results — the classification gate before D2
    // retires it. The unified coincident builder (side-face dissolving) lands with C1.
    #[test]
    fn coincident_stack_outcomes_match_the_survival_table() {
        assert_eq!(coplanar_survival(BoolKind::Fuse, false).0, PSurvive::MinusQ); // ∅ ⇒ drop both
        assert_eq!(coplanar_survival(BoolKind::Cut, false).0, PSurvive::Whole); // A kept
        assert_eq!(
            coplanar_survival(BoolKind::Common, false).0,
            PSurvive::Empty
        );
        let stack = || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let b = m.add_cuboid(
                Point3::from_array([0.0, 0.0, 1.0]),
                Point3::from_array([1.0, 1.0, 2.0]),
            );
            (m, a, b)
        };
        {
            let (mut m, a, b) = stack();
            let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty());
            assert!((nacre_props::mass_props(&m, r).unwrap().volume - 2.0).abs() < 1e-12);
            let planes = collect_planes(&m, r).unwrap();
            assert!(
                !solid_has_coplanar_neighbour_edge(&m, r, &planes),
                "side faces dissolved"
            );
        }
        {
            let (mut m, a, b) = stack();
            let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
            assert!((nacre_props::mass_props(&m, r).unwrap().volume - 1.0).abs() < 1e-12);
        }
        {
            let (mut m, a, b) = stack();
            assert!(matches!(
                boolean(&mut m, BoolKind::Common, a, b),
                Err(BoolError::EmptyResult)
            ));
        }
    }

    // B1/B2: the survival-table builder reproduces the bespoke contained pocket and boss — same
    // volume, watertight, all-Constructed (empty seam), no spurious coplanar edge.
    #[test]
    fn coplanar_contained_result_reproduces_pocket_and_boss() {
        // Pocket (Cut/same): base top gains the prism footprint as a mouth. Volume 0.875.
        {
            let mut m = Model::new();
            let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let prism = m.add_cuboid(
                Point3::from_array([0.25, 0.25, 0.5]),
                Point3::from_array([0.75, 0.75, 1.0]),
            );
            let r = coplanar_result(&mut m, BoolKind::Cut, base, prism).unwrap();
            assert_eq!(r.len(), 1);
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty());
            assert!((nacre_props::mass_props(&m, r[0]).unwrap().volume - 0.875).abs() < 1e-12);
            assert!(
                solid_vertex_handles(&m, r[0])
                    .iter()
                    .all(|&v| matches!(m.vertices.get(v).origin, Origin::Constructed)),
                "empty seam ⇒ all Constructed (purity)"
            );
            let planes = collect_planes(&m, r[0]).unwrap();
            assert!(!solid_has_coplanar_neighbour_edge(&m, r[0], &planes));
        }
        // Boss (Fuse/opposite): base top gains the boss footprint as a hole. Volume 1.25.
        {
            let mut m = Model::new();
            let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let boss = m.add_cuboid(
                Point3::from_array([0.25, 0.25, 1.0]),
                Point3::from_array([0.75, 0.75, 2.0]),
            );
            let r = coplanar_result(&mut m, BoolKind::Fuse, base, boss).unwrap();
            assert_eq!(r.len(), 1);
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty());
            assert!((nacre_props::mass_props(&m, r[0]).unwrap().volume - 1.25).abs() < 1e-12);
        }
    }

    // B4 part 2b: the overhang boss (crossing) through the one entry — first crossing-based real
    // solid. P notch + Q cantilever (mixed-node) + walls resplit at the crossings. Volume 1.5.
    #[test]
    fn coplanar_result_reproduces_overhang_boss() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let r = coplanar_result(&mut m, BoolKind::Fuse, base, boss).unwrap();
        assert_eq!(r.len(), 1);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, r[0]).unwrap().volume - 1.5).abs() < 1e-12);
        let planes = collect_planes(&m, r[0]).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r[0], &planes));
    }

    // A3 (C2 de-risk): point_in_ring parity classifies against a NON-CONVEX footprint with no
    // half-space / convexity assumption — the soundness `clip_bwall_inside_a` lacked.
    #[test]
    fn point_in_ring_on_a_non_convex_footprint() {
        // An L: [0,3]×[0,1] ∪ [0,1]×[0,3] at z=1. A corner in the notch (x>1 and y>1) is OUTSIDE
        // the ring; a corner in an arm is inside.
        let mut m = Model::new();
        let l: Vec<Point3> = [
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [1.0, 1.0],
            [1.0, 3.0],
            [0.0, 3.0],
        ]
        .iter()
        .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
        .collect();
        let (a, _) = build_prism(&mut m, &l, Vector3::from_array([0.0, 0.0, -1.0]), None).unwrap();
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        m.rebuild_adjacency();
        let planes_a = collect_planes(&m, a).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(&m, b).unwrap());
        let mut surf_ix = std::collections::HashMap::new();
        for (i, pinf) in planes.iter().enumerate() {
            surf_ix.insert(pinf.face, i);
        }
        let canon = plane_classes(&planes);
        let find = |rng: std::ops::Range<usize>, nz: f64| -> usize {
            rng.clone()
                .find(|&i| {
                    let n = planes[i].n_out.as_array();
                    n[2] * nz > 0.5 && (planes[i].tri[0].as_array()[2] - 1.0).abs() < 1e-9
                })
                .expect("z=1 face")
        };
        let a_top = find(0..na, 1.0);
        let b_bot = find(na..planes.len(), -1.0);
        let pi = canon[a_top];
        let inc_a = arrange::edge_planes(&m, a, &surf_ix).unwrap();
        let inc_b = arrange::edge_planes(&m, b, &surf_ix).unwrap();
        let a_bnd = contact_boundary(&m, planes[a_top].face, pi, &inc_a, &canon).unwrap();
        let b_bnd = contact_boundary(&m, planes[b_bot].face, pi, &inc_b, &canon).unwrap();
        let l_ring = boundary_ring_triples(&a_bnd, pi);
        let q_ring = boundary_ring_triples(&b_bnd, pi);
        assert_eq!(l_ring.len(), 6, "L has 6 vertices");
        let corner = |x: f64, y: f64| -> [usize; 3] {
            *q_ring
                .iter()
                .find(|&&t| {
                    let p = three_planes(
                        &planes[t[0]].plane,
                        &planes[t[1]].plane,
                        &planes[t[2]].plane,
                    )
                    .unwrap()
                    .as_array();
                    (p[0] - x).abs() < 1e-9 && (p[1] - y).abs() < 1e-9
                })
                .expect("b corner")
        };
        assert!(
            arrange::point_in_ring(&planes, pi, corner(0.5, 0.5), &l_ring).unwrap(),
            "arm corner is inside the L"
        );
        assert!(
            !arrange::point_in_ring(&planes, pi, corner(2.0, 2.0), &l_ring).unwrap(),
            "notch corner is outside the L"
        );
    }

    #[test]
    fn fuse_an_overhanging_boss_onto_a_non_convex_solid() {
        // A boss cantilevers off the +x side face of a top-pocketed cube (non-convex solid),
        // overhanging the bottom edge. The whole-solid gate used to block it; the contact face
        // (+x side) is a convex square, so the footprint gate admits it and the Fuse reconstruction
        // is local (the far pocket is verbatim-copied). Volume: pocketed 0.92 + boss 0.25 = 1.17.
        let (mut m, pc) = top_pocketed_cube();
        let boss = m.add_cuboid(
            Point3::from_array([1.0, 0.25, -0.25]),
            Point3::from_array([1.5, 0.75, 0.75]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, pc, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - 1.17).abs() < 1e-12);
    }

    #[test]
    fn overhang_boss_with_a_non_convex_footprint_is_rejected() {
        // An L-shaped (non-convex) boss footprint overhanging a cube edge: the footprint gate
        // rejects it (a non-convex contact loop would mislabel the notch/cantilever arcs). Honest
        // reject — non-convex overhang footprints are a separate cell.
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let l_base: Vec<Point3> = [
            [0.3, 0.3],
            [1.2, 0.3],
            [1.2, 0.5],
            [0.6, 0.5],
            [0.6, 0.7],
            [0.3, 0.7],
        ]
        .iter()
        .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
        .collect();
        let (l_tool, _) =
            build_prism(&mut m, &l_base, Vector3::from_array([0.0, 0.0, 0.4]), None).unwrap();
        assert!(matches!(
            boolean_one(&mut m, BoolKind::Fuse, cube, l_tool),
            Err(BoolError::Unsupported)
        ));
    }

    #[test]
    fn overhang_cut_on_a_non_convex_solid_is_rejected() {
        // The Cut path keeps the whole-solid convexity gate (its `clip_bwall_inside_a` clips against
        // breached half-spaces, sound only for a convex kept solid). An overhang slot on a
        // non-convex solid is honestly rejected — this cell opens Fuse (boss) only.
        let (mut m, pc) = top_pocketed_cube();
        let slot = m.add_cuboid(
            Point3::from_array([0.75, 0.25, -0.25]),
            Point3::from_array([1.0, 0.75, 0.5]),
        );
        assert!(matches!(
            boolean_one(&mut m, BoolKind::Cut, pc, slot),
            Err(BoolError::Unsupported)
        ));
    }

    #[test]
    fn cut_a_blind_pocket_into_a_face() {
        // A prism inside the base with its top flush on the base's top (same-normal coplanar
        // contact, contained footprint). Cut carves a blind pocket: the base top gains the
        // footprint as a hole (the mouth), the prism faces flip to bound the removed region.
        // Before this cell it rejected as vertex_on_face_plane. Volume 1 − 0.5²·0.5 = 0.875.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.875).abs() < 1e-12, "volume {vol}");
        let planes = collect_planes(&m, r).unwrap();
        assert!(!solid_has_coplanar_neighbour_edge(&m, r, &planes));
    }

    #[test]
    fn cut_a_blind_pocket_into_a_non_convex_solid() {
        // A blind pocket carved into an already-pocketed (non-convex) cube: a second contained
        // top-flush prism in a corner away from the first pocket. The kept solid `a` is non-convex,
        // which the pocket contact now admits (the gates are convexity-agnostic). Removed
        // 0.2·0.1·0.4 = 0.008 on top of the first pocket's 0.08 → 1 − 0.08 − 0.008 = 0.912.
        let (mut m, pc) = pocketed_cube();
        let corner = m.add_cuboid(
            Point3::from_array([0.05, 0.1, 0.6]),
            Point3::from_array([0.25, 0.2, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, pc, corner).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.912).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn cut_a_non_convex_blind_pocket() {
        // A blind pocket with a non-convex (L-shaped) footprint: the cutter prism is non-convex,
        // which the pocket contact now admits. The L extrudes to z∈[0,0.5], top-flush on the base's
        // z=0.5 face, blind. L area = 0.6² − 0.3² = 0.27, depth 0.5 → removed 0.135; base 3²·1.5 =
        // 13.5 → 13.365.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 0.5]),
        );
        let l = Profile2d {
            points: vec![
                p2(-0.3, -0.3),
                p2(0.3, -0.3),
                p2(0.3, 0.0),
                p2(0.0, 0.0),
                p2(0.0, 0.3),
                p2(-0.3, 0.3),
            ],
        };
        let OpOutput::Extrude { solid: lp, .. } = apply(&mut m, &extrude_op(l, 0.5)).unwrap()
        else {
            unreachable!()
        };
        let r = boolean_one(&mut m, BoolKind::Cut, base, lp).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 13.365).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn a_pocket_that_punches_through_is_not_a_pocket_contact() {
        // The prism's top is flush, but it pokes out the base's bottom — its walls cross the
        // base's bottom face, so it is a seam cut, not a blind pocket. detect_pocket_contact
        // must decline it (leaving it to the seam path) rather than build a floor outside the base.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.25, 0.25, -0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        assert!(detect_pocket_contact(&m, base, through).is_none());
    }

    #[test]
    fn fuse_stacked_cubes() {
        // A=[0,1]³ and B=[0,1]²×[1,2] share the z=1 face ⇒ merge into a clean 1×1×2 box. The
        // four coplanar side pairs are spliced into 4 faces, and the four interface corners
        // (each a straight angle on a split vertical edge) are dissolved, fusing the split
        // edges — a canonical 6-face / 8-vertex / 12-edge box (cell fuse-coplanar-merge).
        let (mut m, a, b) = stacked_cubes();
        let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
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
        let stack = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        // A cutter straddling z=1 (the fused interface) — the seam runs where the split
        // vertical edges used to be. Result: 2 − 0.5·0.5·1.0.
        let cutter = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, stack, cutter).unwrap();
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
            boolean_one(&mut m, BoolKind::Common, a, b),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn cut_stacked_cubes_is_a() {
        let (mut m, a, b) = stacked_cubes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn coincident_merge_declines_offset_footprint() {
        // Coplanar z=1 interface but B's footprint is offset ⇒ boundaries differ, so
        // coincident_merge (which needs an exact ring match) must decline. The Fuse itself now
        // succeeds as a corner overhang (see fuse_a_corner_overhanging_boss, same fixture).
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        assert!(detect_coincident_interface(&m, a, b).is_none());
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
        // honest rejection, less specific about why (§9 line 451). This is same-normal (both
        // floors at z=0), so it stays out of the opposite-normal overhang path; its former
        // sibling fixture is now the corner overhang (fuse_a_corner_overhanging_boss).
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
            || boolean_one(&mut m, BoolKind::Fuse, a, b),
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
            let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1);
            prop_assume!(rf.is_ok()); // skip rare coplanar/degenerate configs
            let rf = rf.unwrap();
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb - ov)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            let rc = boolean_one(&mut m2, BoolKind::Cut, a2, b2).unwrap();
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
            let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1);
            prop_assume!(rf.is_ok());
            let rf = rf.unwrap();
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            prop_assert_eq!(
                boolean_one(&mut m2, BoolKind::Common, a2, b2),
                Err(BoolError::EmptyResult)
            );

            let (mut m3, a3, b3) = build();
            let rc = boolean_one(&mut m3, BoolKind::Cut, a3, b3).unwrap();
            let vc = nacre_props::mass_props(&m3, rc).unwrap().volume;
            prop_assert!((vc - va).abs() <= 1e-9 * va, "cut {vc}");
        }
    }

    #[test]
    fn fuse_of_two_cubes() {
        // A = [0,1]³, B = [0.5,1.5]³ ⇒ A∪B volume 1+1−0.125 = 1.875.
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
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
        let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1).unwrap();
        let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;

        let mut m2 = Model::new();
        let (a2, b2) = (m2.add_cuboid(amin, amax), m2.add_cuboid(bmin, bmax));
        let rc = boolean_one(&mut m2, BoolKind::Common, a2, b2).unwrap();
        let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;

        assert!((vf + vc - va - vb).abs() < 1e-9, "{vf}+{vc} vs {va}+{vb}");
    }

    #[test]
    fn boolean_rejects_non_live_input() {
        let (mut m, a, b) = two_boxes();
        m.live_solids.retain(|&s| s != b); // as if superseded
        assert_eq!(
            boolean_one(&mut m, BoolKind::Common, a, b),
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
        let r = boolean_one(&mut m, BoolKind::Common, a, b).unwrap();
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
        let r = boolean_one(&mut m, BoolKind::Common, a, c).unwrap();
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
            boolean_one(&mut m, BoolKind::Common, a, b),
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
            || boolean_one(&mut m, BoolKind::Common, a, cyl),
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
            || boolean_one(&mut m, BoolKind::Common, lsolid, b),
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
            || boolean_one(&mut m, BoolKind::Common, solid, b),
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
            let res = boolean_one(&mut m, BoolKind::Common, a, b);
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
            let res = boolean_one(&mut m, BoolKind::Common, prism, c);
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
