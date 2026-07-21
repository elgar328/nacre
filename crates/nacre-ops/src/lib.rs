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
mod trace;

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
    /// The trace arrangement on one plane class nested a hole whose containment depth exceeds one,
    /// or produced more than one unbounded contour (several disjoint bodies on the plane).
    /// `nest_cells` resolves any number of holes at depth one inside one outer loop; deeper nesting
    /// and multiple bodies are honestly rejected until the general nesting cell lands. Two distinct
    /// tags so a refactor cannot silently merge the conditions.
    pub const HOLE_DEPTH: &str = "hole_depth";
    pub const HOLE_ROOTS: &str = "hole_roots";
    /// A face whose boundary never crosses the seam, yet the seam lies on its plane — the
    /// convex path only.
    ///
    pub const MISSING_SEAM: &str = "missing_seam";
    /// A solid's planar section returned more than one loop (an outer ring plus a hole, e.g. a
    /// slot piercing a cavity). The section-clip reconstruction assumes a single loop; a
    /// multi-loop section is honestly rejected until the multi-loop cell lands, so the
    /// single-loop assumption is never silently reached.
    pub const SECTION_MULTI_LOOP: &str = "section_multi_loop";
    /// Canonicalizing a section's vertex triples through the plane classes folded two *distinct*
    /// section vertices onto one triple (e.g. `{α,β,cutter_top}` and `{α,β,base_top}` when both
    /// caps fold to the contact class π). Merging them would silently drop a point, so it is
    /// rejected instead of collapsed (DNA: never silently wrong).
    pub const SECTION_TRIPLE_COLLISION: &str = "section_triple_collision";
    /// A ∂P vertex coincides exactly with a section corner on the contact plane (a four-plane
    /// point). The scoped contact-plane flush (R2) covers a vertex strictly *on* a contact chord,
    /// but a vertex *at* a chord endpoint is an ambiguous fan — honestly rejected until the general
    /// T-junction resolver (C1) lands.
    pub const FLUSH_VERTEX_COINCIDENT: &str = "flush_vertex_coincident";
    /// A face is coplanar-overlapping with *two or more* faces of the other solid on one plane
    /// class (a non-convex `other` seating twice on one plane). The per-face classifier (B4-R0)
    /// splits against a single coincident face; multiple overlaps are out of scope until a later
    /// cell, honestly rejected rather than picking one arbitrarily.
    pub const COPLANAR_OVERLAP_MULTI: &str = "coplanar_overlap_multi";
    /// The unified driver found no shared plane class between the two solids — it was invoked on a
    /// pair with no coplanar contact at all (should be gated out upstream).
    pub const NO_COPLANAR_CONTACT: &str = "no_coplanar_contact";
    /// A `Whole`-survival contact face whose footprint OVERLAPS the other's (∂P × ∂Q cross) rather
    /// than nesting, in the one such case still unbuilt. `Whole` has two entries: `Fuse`/same-normal,
    /// which the E1 union cell now builds, and `Cut`/opposite-normal, which is exact whenever the
    /// contact plane separates the two solids (nothing to remove). What is left is a `Cut` whose tool
    /// reaches back across that plane — a pin below its own contact face — where the cut owes a notch
    /// this path cannot yet cut. Honest reject rather than a whole cap that ignores the pin.
    pub const COPLANAR_MERGE: &str = "coplanar_merge";
    /// The winding driver (M-C) was handed a pair with no seam (containment/disjoint) — the
    /// seam-free branch (`contained_result`) is not yet wired into the unified driver.
    pub const WINDING_NO_SEAM: &str = "winding_no_seam";
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
    RayCross, orient2d, plane_plane, plane_side, planes_coplanar, three_plane_orient3d,
    three_planes,
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
/// **Coverage:** planar solids — all three kinds go through one exact seam path
/// (`general_boolean`) for transversal contact, and one coplanar path
/// (`coplanar_result_unified`) when the operands share a contact plane. Anything else is
/// rejected with [`BoolError`]. Transactional: each variant computes its result in local
/// structures and pushes only after every degeneracy check passes, so a rejected boolean
/// leaves the model untouched.
///
/// **Dispatch (F2 collapse):** two routes, chosen by one exact question — *is there a genuine
/// coplanar contact?* ([`coplanar_contact_count`]). The seven bespoke detectors that used to
/// gate this (contained / pocket / overhang boss·cut·common / single-shared-plane /
/// multi-contact) all funnelled into the same builder, so they were pure routing; collapsing
/// them to the count both retired ~775 lines and opened cases their narrow gates declined
/// (a non-convex overhang footprint, a slot or corner cut punching through the far face).
pub fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    if !model.live_solids.contains(&a) || !model.live_solids.contains(&b) {
        return Err(BoolError::InputNotLive);
    }
    // The arrangement engine (`trace.rs`) is the sole boolean path: one per-plane-class 2D
    // arrangement handles transverse, coplanar-contact, coincident, contained and disjoint cases,
    // and cleans its own output (coplanar-face merge) so results are chainable.
    crate::trace::boolean_via_trace(model, kind, a, b)
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
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
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

/// The minimal per-op plane table two solids share for [`arrange::point_in_solid_idx`]: the
/// concatenated plane list (`a`'s then `b`'s), the face→index map, and each solid's
/// [`arrange::EdgePlanes`]. Detectors that classify vertices exactly but lack
/// `overlap_fuse_cut`'s setup build it once, then classify each vertex without rebuilding.
/// Indices into the returned `planes`/`surf_ix` are shared, so a vertex of `a` and a face of
/// `b` compose in one space.
type PlaneSetup = (
    Vec<PlaneInfo>,
    HashMap<Handle<Face>, usize>,
    arrange::EdgePlanes,
    arrange::EdgePlanes,
    Vec<usize>,
);

fn plane_index_setup(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<PlaneSetup, BoolError> {
    let mut planes = collect_planes(model, a)?;
    planes.extend(collect_planes(model, b)?);
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        surf_ix.insert(pi.face, i);
    }
    let inc_a = arrange::edge_planes(model, a, &surf_ix)?;
    let inc_b = arrange::edge_planes(model, b, &surf_ix)?;
    // Coplanar planes → one class, so `point_in_solid_idx` never forms a `det=0` triple from two
    // coplanar faces meeting a vertex (the cantilever-step degeneracy).
    let canon = plane_classes(&planes);
    Ok((planes, surf_ix, inc_a, inc_b, canon))
}

/// Whether three `Pt3` are **exactly collinear** (zero-area triangle), decided on their
/// pre-rotation rational `base` coordinates. A rigid rotation preserves collinearity, and three
/// vertices of one solid share a rotation chain, so their bases are comparable; all three
/// coordinate-plane projections of `(b−a)×(c−a)` must vanish (exact `Rat`, no tolerance). An
/// i128 overflow returns `false` (treat as non-collinear): a genuinely-collinear triangle then
/// stays and is at worst rejected `RAY_DEGENERATE`, never falsely skipped (which would drop a
/// real crossing — silent-wrong). Unrotated vertices carry `base == coord`, so this is the exact
/// zero-area (collinear) test on the vertices' rotation definitions.
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
    let mut node_handle = |model: &mut Model, node: Node| -> Result<Handle<Vertex>, BoolError> {
        if let Some(&h) = vh.get(&node) {
            return Ok(h);
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
                // A face references a seam node whose triple was not welded into `seam` — a
                // reconstruction dropped a crossing. Reject (never panic): an unmodeled flush
                // topology must decline honestly, not abort the kernel (DNA).
                let sv = seam
                    .iter()
                    .find(|s| s.triple == triple)
                    .ok_or_else(|| reject(tag::MISSING_SEAM))?;
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
        Ok(handle)
    };
    // Materialize all vertex handles first. Outer then inner, rings in order: `vh`'s
    // first-appearance order fixes the vertex handles, and replay depends on it.
    for lf in faces {
        for &node in lf.loop_nodes.iter().chain(lf.inner.iter().flatten()) {
            node_handle(model, node)?;
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
    // Closed-shell guard: a 2-manifold b-rep uses every edge exactly twice (once from each of the
    // two faces that share it). A reconstruction that emits a face set with a dangling edge (use
    // count 1) or a pinched one (>2) is not a solid — `validate` would call it `NonManifoldEdge`,
    // but `boolean` never runs `validate` on its own output, so without this the caller receives a
    // silently invalid solid. Honest-reject instead ("honest-reject > silent-wrong", overview §1).
    // Counted on the welded `Handle<Edge>`s, so it is exact and coordinate-free.
    {
        let mut uses: HashMap<Handle<Edge>, usize> = HashMap::new();
        for &fh in &face_handles {
            let f = model.faces.get(fh);
            for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &l.half_edges {
                    *uses.entry(he.edge).or_default() += 1;
                }
            }
        }
        if uses.values().any(|&n| n != 2) {
            return Err(reject(tag::NON_MANIFOLD_EDGE));
        }
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
/// Union-find root with path compression. Drives component grouping in [`unify_coplanar_faces`].
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

/// An unordered edge key: the two nodes in a fixed order, so `{a,b}` and `{b,a}` collide.
fn norm_edge(a: Node, b: Node) -> (Node, Node) {
    fn rank(n: Node) -> (u8, usize, usize, usize) {
        match n {
            Node::Orig(v) => (0, v.index() as usize, 0, 0),
            Node::Seam(t) => (1, t[0], t[1], t[2]),
        }
    }
    if rank(a) <= rank(b) { (a, b) } else { (b, a) }
}

/// Splice loop `b` into loop `a` across their shared edge `{u, v}` (present in `a` as `u→v`
/// or `v→u`, and in `b` with the opposite orientation), dropping that edge — the loop-join
/// core of the coplanar merge. Node identity throughout, agnostic to `Orig`/`Seam` away from
/// the shared edge. The shared edge's two endpoints remain as (collinear, for a flat merge)
/// boundary points; [`unify_coplanar_faces`] dissolves them globally afterwards.
fn splice_along(a: &[Node], b: &[Node], u: Node, v: Node) -> Vec<Node> {
    let m = a.len();
    let ia = (0..m)
        .find(|&i| {
            let (p, q) = (a[i], a[(i + 1) % m]);
            (p == u && q == v) || (p == v && q == u)
        })
        .expect("shared edge lies on A's loop");
    let (x, y) = (a[ia], a[(ia + 1) % m]);
    let n = b.len();
    // B traverses the shared edge the opposite way: `b[ib] = Y`, `b[ib+1] = X`.
    let ib = (0..n)
        .find(|&j| b[j] == y && b[(j + 1) % n] == x)
        .expect("B's face shares the edge, opposite orientation");
    let mut nodes = Vec::with_capacity(m + n - 2);
    // A from Y (ia+1) all the way to X (ia): the whole A loop, less the dropped edge.
    for t in 0..m {
        nodes.push(a[(ia + 1 + t) % m]);
    }
    // B's interior, strictly between X and Y (skip the shared edge's two endpoints).
    for t in 2..n {
        nodes.push(b[(ib + t) % n]);
    }
    nodes
}

/// Merge coplanar, same-normal, hole-free faces that share an original edge into one face,
/// then dissolve straight-angle (globally degree-2, exactly collinear) original vertices — a
/// T-junction-free defeature. Generalizes the coincident Fuse merge off the interface ring to
/// any plane class via an edge→faces map: two faces meeting flat across an original edge fuse,
/// and the now-redundant edge and its straight-angle vertices disappear.
///
/// Reused, not reinvented: `plane_classes` (`canon`) decides coplanarity, `n_out.dot > 0` the
/// shared-normal guard — a real dihedral, or an opposite-normal cantilever step, is kept — and
/// `pt3_base_collinear` the exact straight-angle test. `assemble_fuse_cut` downstream is
/// unchanged. Detection is deterministic (faces in index order, edges in loop order) for replay.
///
/// Deferred, each a safe no-op: holed faces (`inner`), seam-shared edges (detection needs `Orig`
/// endpoints), opposite-normal coplanar pairs, and non-disk components (a `debug_assert` guards
/// the disk assumption). The one caller today is the coincident Fuse; a future coplanar-contact
/// path is the second.
fn unify_coplanar_faces(
    model: &Model,
    faces: Vec<LocalFace>,
    planes: &[PlaneInfo],
    canon: &[usize],
) -> Vec<LocalFace> {
    let n = faces.len();
    // 1. Edge (unordered node pair) → the faces carrying it, over outer loops. Built in face
    //    order, so each value is ascending face indices.
    let mut edge_faces: HashMap<(Node, Node), Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        let k = lf.loop_nodes.len();
        for i in 0..k {
            let e = norm_edge(lf.loop_nodes[i], lf.loop_nodes[(i + 1) % k]);
            edge_faces.entry(e).or_default().push(fi);
        }
    }
    // 2. Mergeable edges (deterministic order): shared by exactly two faces that are coplanar,
    //    same-normal, hole-free, with both endpoints `Orig`. Union the incident faces.
    let mut mergeable: Vec<(Node, Node)> = Vec::new();
    let mut comp: Vec<usize> = (0..n).collect();
    for (fi, lf) in faces.iter().enumerate() {
        let k = lf.loop_nodes.len();
        for i in 0..k {
            let (u, v) = (lf.loop_nodes[i], lf.loop_nodes[(i + 1) % k]);
            let fs = &edge_faces[&norm_edge(u, v)];
            if fs.len() != 2 || fs[0] != fi {
                continue; // process each edge once, from its lower-index face
            }
            if !matches!(u, Node::Orig(_)) || !matches!(v, Node::Orig(_)) {
                continue; // seam-shared edge — deferred
            }
            let (li, lj) = (&faces[fs[0]], &faces[fs[1]]);
            if canon[li.plane_idx] != canon[lj.plane_idx] {
                continue; // not coplanar
            }
            if planes[li.plane_idx].n_out.dot(planes[lj.plane_idx].n_out) <= 0.0 {
                continue; // opposite normal — a genuine fold/step, keep
            }
            if !li.inner.is_empty() || !lj.inner.is_empty() {
                continue; // holed — deferred
            }
            mergeable.push((u, v));
            let (ri, rj) = (uf_find(&mut comp, fs[0]), uf_find(&mut comp, fs[1]));
            if ri != rj {
                comp[ri] = rj;
            }
        }
    }
    if mergeable.is_empty() {
        return faces; // nothing coplanar-adjacent
    }
    // 3. Fold each component into one loop by splicing across its mergeable edges.
    let mut active: Vec<Option<LocalFace>> = faces.into_iter().map(Some).collect();
    let mut slot: Vec<usize> = (0..n).collect(); // orig face → active slot currently holding it
    for (u, v) in mergeable {
        let fs = &edge_faces[&norm_edge(u, v)];
        let (si, sj) = (uf_find(&mut slot, fs[0]), uf_find(&mut slot, fs[1]));
        if si == sj {
            continue; // already one face (joined via another edge of this component)
        }
        let lf_b = active[sj].take().expect("active slot");
        let lf_a = active[si].as_ref().expect("active slot");
        let loop_nodes = splice_along(&lf_a.loop_nodes, &lf_b.loop_nodes, u, v);
        debug_assert!(
            loop_nodes.len() >= 3,
            "coplanar merge across a non-disk component"
        );
        active[si] = Some(LocalFace {
            plane_idx: lf_a.plane_idx,
            loop_nodes,
            inner: Vec::new(),
            flip: false,
        });
        slot[sj] = si;
    }
    let mut out: Vec<LocalFace> = active.into_iter().flatten().collect();
    // 4. Dissolve straight-angle Orig vertices: globally degree-2 (one edge line through them)
    //    and exactly collinear. Drop from every incident loop at once, so a vertex that is a
    //    real corner on any face survives (no T-junction).
    let mut nbrs: HashMap<Node, HashSet<Node>> = HashMap::new();
    for lf in &out {
        let k = lf.loop_nodes.len();
        for i in 0..k {
            let (a, b) = (lf.loop_nodes[i], lf.loop_nodes[(i + 1) % k]);
            nbrs.entry(a).or_default().insert(b);
            nbrs.entry(b).or_default().insert(a);
        }
    }
    let mut drop: HashSet<Node> = HashSet::new();
    for (&node, ns) in &nbrs {
        let Node::Orig(v) = node else { continue };
        if ns.len() != 2 {
            continue;
        }
        let mut it = ns.iter();
        let (Node::Orig(a), Node::Orig(b)) = (*it.next().unwrap(), *it.next().unwrap()) else {
            continue; // a seam neighbour — leave the vertex (deferred)
        };
        let (Ok(pa), Ok(pv), Ok(pb)) = (
            nacre_tip::vertex_pt3(model, a),
            nacre_tip::vertex_pt3(model, v),
            nacre_tip::vertex_pt3(model, b),
        ) else {
            continue;
        };
        if pt3_base_collinear(&pa, &pv, &pb) {
            drop.insert(node);
        }
    }
    if !drop.is_empty() {
        for lf in &mut out {
            lf.loop_nodes.retain(|nd| !drop.contains(nd));
        }
    }
    out
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

    /// The exact index-plane classifier [`arrange::point_in_solid_idx`] on hand-known answers.
    /// The query is always a *vertex* (it must carry plane identity), classified against the
    /// other solid — exactly the live use (`a`'s vertices vs `b`). Not an f64 cross-check: the
    /// answers are computed by hand, so this is an independent oracle.
    #[test]
    fn point_in_solid_idx_matches_known_answers() {
        let check = |amin: [f64; 3],
                     amax: [f64; 3],
                     bmin: [f64; 3],
                     bmax: [f64; 3],
                     inside: &dyn Fn([f64; 3]) -> bool| {
            let mut m = Model::new();
            let qa = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
            let ob = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
            m.rebuild_adjacency();
            let mut planes = collect_planes(&m, qa).unwrap();
            planes.extend(collect_planes(&m, ob).unwrap());
            let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
            for (i, pi) in planes.iter().enumerate() {
                surf_ix.insert(pi.face, i);
            }
            let inc_a = arrange::edge_planes(&m, qa, &surf_ix).unwrap();
            let inc_b = arrange::edge_planes(&m, ob, &surf_ix).unwrap();
            let canon = plane_classes(&planes);
            for vh in solid_vertex_handles(&m, qa) {
                let p = m.vertices.get(vh).point.as_array();
                let want = if inside(p) {
                    Side::Inside
                } else {
                    Side::Outside
                };
                let got = arrange::point_in_solid_idx(
                    &m, vh, &inc_a, ob, &inc_b, &planes, &surf_ix, &canon,
                )
                .unwrap_or_else(|e| panic!("vertex {p:?}: {e:?}"));
                assert_eq!(got, want, "vertex {p:?} vs B[{bmin:?}..{bmax:?}]");
            }
        };
        // Overlap: A[0,2]³ vs B[1,3]³ — a corner is inside iff every coord (∈{0,2}) is in (1,3),
        // i.e. only (2,2,2).
        check([0.; 3], [2.; 3], [1.; 3], [3.; 3], &|p| {
            p.iter().all(|&c| c > 1.0 && c < 3.0)
        });
        // Containment: B[-1,3]³ ⊃ A[0,2]³ — every corner Inside.
        check([0.; 3], [2.; 3], [-1.; 3], [3.; 3], &|_| true);
        // Disjoint: B[5,6]³ — every corner Outside.
        check([0.; 3], [2.; 3], [5.; 3], [6.; 3], &|_| false);
        // Coplanar-disjoint (shares the z=0/z=1 planes, footprints apart): A[0,1]³ vs
        // B[1,2]×[3,4]×[0,1]. A's corners lie *on* B's z-planes yet off its footprint — all
        // Outside, and crucially decided (not NO_CLEAR_RAY). This is the 4078 query species.
        check([0.; 3], [1.; 3], [1., 3., 0.], [2., 4., 1.], &|_| false);
    }

    /// The classifier on a genuinely **non-convex** `other` (an L-prism, whose notch a
    /// convex all-half-spaces test gets wrong) and a **holed** `other` (a hollow box, whose
    /// void must read Outside via the outer+cavity shell sum). These are the cases boxes cannot
    /// exercise — the whole reason classification is winding-parity, not per-face plane-side.
    #[test]
    fn point_in_solid_idx_on_non_convex_and_cavity() {
        let classify_all =
            |m: &Model, qs: Handle<Solid>, os: Handle<Solid>, want: &dyn Fn([f64; 3]) -> bool| {
                let mut planes = collect_planes(m, qs).unwrap();
                planes.extend(collect_planes(m, os).unwrap());
                let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
                for (i, pi) in planes.iter().enumerate() {
                    surf_ix.insert(pi.face, i);
                }
                let inc_q = arrange::edge_planes(m, qs, &surf_ix).unwrap();
                let inc_o = arrange::edge_planes(m, os, &surf_ix).unwrap();
                let canon = plane_classes(&planes);
                for vh in solid_vertex_handles(m, qs) {
                    let p = m.vertices.get(vh).point.as_array();
                    let got = arrange::point_in_solid_idx(
                        m, vh, &inc_q, os, &inc_o, &planes, &surf_ix, &canon,
                    )
                    .unwrap_or_else(|e| panic!("vertex {p:?}: {e:?}"));
                    assert_eq!(got == Side::Inside, want(p), "at {p:?}");
                }
            };

        // Non-convex: L-prism (bottom bar x∈[0,2]×y∈[0,1] + left column x∈[0,1]×y∈[0,2],
        // z∈[0,1]; notch x>1∧y>1 is empty). Query cube corners at {0.5,1.5}²×{0.25,0.75} —
        // the (1.5,1.5) corners sit in the notch (Outside), the rest inside the L.
        {
            let (mut m, lp) = l_prism();
            let q = m.add_cuboid(
                Point3::from_array([0.5, 0.5, 0.25]),
                Point3::from_array([1.5, 1.5, 0.75]),
            );
            m.rebuild_adjacency();
            classify_all(&m, q, lp, &|p| !(p[0] > 1.0 && p[1] > 1.0));
        }

        // Cavity: hollow box, outer [0,4]³ with a concentric void [1,3]³. Query cube corners at
        // {0.5,1.5}³ — (1.5,1.5,1.5) is in the void (Outside), the rest in the material wall.
        {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.; 3]), Point3::from_array([4.; 3]));
            let b = m.add_cuboid(Point3::from_array([1.; 3]), Point3::from_array([3.; 3]));
            let void = m.reversed_shell(m.solids.get(b).outer);
            let a_outer = m.solids.get(a).outer;
            let hollow = m.push_solid(Solid {
                outer: a_outer,
                cavities: vec![void],
            });
            let q = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
            m.rebuild_adjacency();
            // Inside the wall iff not strictly inside the void (all coords in (1,3)).
            classify_all(&m, q, hollow, &|p| !p.iter().all(|&c| c > 1.0 && c < 3.0));
        }
    }

    /// canon-awareness (Cell 2.9): `Fuse(base, boss)` leaves two coplanar faces on z=1 (base-top
    /// `+z`, boss-bottom `−z` — a cantilever step, not mergeable). Their shared incident vertices
    /// give `det=0` triples that panic the predicate without canon. With canon those planes
    /// collapse to one class, so the classifier runs **panic-free**, agrees with the f64 ray where
    /// both decide, and honestly rejects (never panics) a degenerate step vertex.
    #[test]
    fn point_in_solid_idx_canon_handles_coplanar_step() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        let probe = m.add_cuboid(
            Point3::from_array([1.1, 0.35, 0.5]),
            Point3::from_array([1.4, 0.65, 2.5]),
        );
        m.rebuild_adjacency();
        let mut planes = collect_planes(&m, probe).unwrap();
        planes.extend(collect_planes(&m, overhung).unwrap());
        let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
        for (i, pi) in planes.iter().enumerate() {
            surf_ix.insert(pi.face, i);
        }
        let inc_p = arrange::edge_planes(&m, probe, &surf_ix).unwrap();
        let inc_o = arrange::edge_planes(&m, overhung, &surf_ix).unwrap();
        let canon = plane_classes(&planes);
        // Hand-known membership (no f64 needed): overhung = base[0,1]³ ∪ boss[.5,.25,1]-[1.5,.75,2];
        // probe is a box. idx must classify every vertex (both directions) panic-free and match.
        let in_overhung = |p: [f64; 3]| {
            let inb = |lo: [f64; 3], hi: [f64; 3]| (0..3).all(|k| p[k] >= lo[k] && p[k] <= hi[k]);
            inb([0.0; 3], [1.0; 3]) || inb([0.5, 0.25, 1.0], [1.5, 0.75, 2.0])
        };
        let in_probe = |p: [f64; 3]| {
            p[0] >= 1.1 && p[0] <= 1.4 && p[1] >= 0.35 && p[1] <= 0.65 && p[2] >= 0.5 && p[2] <= 2.5
        };
        for vh in solid_vertex_handles(&m, probe) {
            let p = m.vertices.get(vh).point.as_array();
            let got = arrange::point_in_solid_idx(
                &m, vh, &inc_p, overhung, &inc_o, &planes, &surf_ix, &canon,
            )
            .unwrap_or_else(|e| panic!("probe vertex {p:?} vs overhung: {e:?}"));
            assert_eq!(
                got == Side::Inside,
                in_overhung(p),
                "probe vertex {p:?} vs overhung"
            );
        }
        for vh in solid_vertex_handles(&m, overhung) {
            let p = m.vertices.get(vh).point.as_array();
            let got = arrange::point_in_solid_idx(
                &m, vh, &inc_o, probe, &inc_p, &planes, &surf_ix, &canon,
            )
            .unwrap_or_else(|e| panic!("overhung vertex {p:?} vs probe: {e:?}"));
            assert_eq!(
                got == Side::Inside,
                in_probe(p),
                "overhung vertex {p:?} vs probe"
            );
        }
    }

    /// Go/no-go for the classifier-unification track (plan R7): on grid-aligned, shared-
    /// coordinate two-box configs — where the exact ray, constrained to the query's own axis
    /// planes, is most prone to grazing — does `point_in_solid_idx` ever reject
    /// (`NO_CLEAR_RAY`) a strictly-in/out vertex? A strict vertex's axis rays cross the other
    /// box's faces at footprint-interior points, never its edges, so the answer should be *zero*
    /// regressions. Each answer is also checked against the hand-computed box membership.
    #[test]
    fn point_in_solid_idx_no_reject_regression_on_aligned() {
        // (A box, B box) sharing coordinates in the adversarial grid.
        type BoxPair = ([f64; 3], [f64; 3], [f64; 3], [f64; 3]);
        let configs: &[BoxPair] = &[
            ([0.; 3], [2.; 3], [1.; 3], [3.; 3]),           // corner overlap
            ([0.; 3], [2.; 3], [1., 1., 1.], [2., 2., 2.]), // B is A's octant (shares corner)
            ([0.; 3], [2.; 3], [-1.; 3], [3.; 3]),          // B ⊃ A
            ([0.; 3], [2.; 3], [0.5, 0.5, 0.5], [1.5, 1.5, 1.5]), // B ⊂ A
            ([0.; 3], [2.; 3], [1., 0., 0.], [3., 2., 2.]), // face-flush slab (shares x=... none; y,z flush)
            ([0.; 3], [2.; 3], [5.; 3], [6.; 3]),           // disjoint
            ([0.; 3], [1.; 3], [0., 0., 1.], [1., 1., 2.]), // stacked on shared z=1 face
            ([0.; 3], [3.; 3], [1., 1., -1.], [2., 2., 4.]), // thin bar piercing through
        ];
        // Strict membership of `p` in axis box [lo,hi]: Some(true/false), None on its boundary.
        let strict = |p: [f64; 3], lo: [f64; 3], hi: [f64; 3]| -> Option<bool> {
            let mut on = false;
            for k in 0..3 {
                if p[k] < lo[k] || p[k] > hi[k] {
                    return Some(false); // a coord beyond the slab ⇒ strictly outside
                }
                if p[k] == lo[k] || p[k] == hi[k] {
                    on = true;
                }
            }
            if on { None } else { Some(true) }
        };
        let (mut compared, mut regress) = (0u32, 0u32);
        for &(amin, amax, bmin, bmax) in configs {
            let mut m = Model::new();
            let sa = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
            let sb = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
            m.rebuild_adjacency();
            let mut planes = collect_planes(&m, sa).unwrap();
            planes.extend(collect_planes(&m, sb).unwrap());
            let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
            for (i, pi) in planes.iter().enumerate() {
                surf_ix.insert(pi.face, i);
            }
            let inc_a = arrange::edge_planes(&m, sa, &surf_ix).unwrap();
            let inc_b = arrange::edge_planes(&m, sb, &surf_ix).unwrap();
            let canon = plane_classes(&planes);
            for (qs, os, inc_q, inc_o, olo, ohi) in [
                (sa, sb, &inc_a, &inc_b, bmin, bmax),
                (sb, sa, &inc_b, &inc_a, amin, amax),
            ] {
                for vh in solid_vertex_handles(&m, qs) {
                    let p = m.vertices.get(vh).point.as_array();
                    let Some(want) = strict(p, olo, ohi) else {
                        continue; // on the other box's boundary — honest-reject territory, skip
                    };
                    let idx = arrange::point_in_solid_idx(
                        &m, vh, inc_q, os, inc_o, &planes, &surf_ix, &canon,
                    );
                    match idx {
                        Ok(s) => {
                            compared += 1;
                            let inside = s == Side::Inside;
                            assert_eq!(
                                inside, want,
                                "idx wrong at {p:?} vs box [{olo:?}..{ohi:?}]"
                            );
                        }
                        Err(_) => regress += 1, // NO_CLEAR_RAY on a strict vertex = R7 regression
                    }
                }
            }
        }
        assert!(compared > 0, "measured nothing");
        assert_eq!(
            regress, 0,
            "R7: idx rejected {regress} strict vertices f64 would classify"
        );
        eprintln!("go/no-go: compared={compared} regress=0 — R7 benign on aligned grid");
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

    /// Classify vertex `vh` (of `vh_solid`) in/out of `target` on the exact index-plane
    /// substrate — the test replacement for the retired f64 `point_in_solid`.
    fn classify_in(
        m: &Model,
        vh: Handle<Vertex>,
        vh_solid: Handle<Solid>,
        target: Handle<Solid>,
    ) -> Side {
        let (planes, surf_ix, inc_v, inc_o, canon) =
            plane_index_setup(m, vh_solid, target).unwrap();
        arrange::point_in_solid_idx(m, vh, &inc_v, target, &inc_o, &planes, &surf_ix, &canon)
            .unwrap()
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

    // A hollow operand in a COPLANAR contact — the combination nothing covered until now. The
    // cavity goldens above all take the seam path (transversal cuts) and every coplanar golden uses
    // solid operands, so the intersection of the two was a blind spot, and the coplanar driver had
    // never received the all-shell patch the seam front-end got in cell (5c-in). It emitted only
    // outer-shell faces, so the void vanished: the fuse read 27.0625 — the *un-hollowed* cube plus
    // the boss — with `cavities: 0` and a clean `validate`, because what remained was still a
    // closed shell. Silent-wrong, invisible to every guard. Now the driver walks all shells.
    #[test]
    fn a_hollow_part_takes_a_coplanar_boss() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        // Top-flush boss on z=3: a genuine coplanar contact, clear of the void's planes.
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 3.0]),
            Point3::from_array([0.75, 0.75, 4.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Fuse, hollow, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        // 27 − 1 void + 0.25·0.25·1 boss.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 26.0625).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
    }

    #[test]
    fn a_hollow_part_takes_a_coplanar_pocket() {
        // The Cut twin of the boss case: a top-flush pocket sunk into a hollow part. Same blind
        // spot, same silent-wrong before the fix (26.96875 with the void gone).
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let tool = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 2.5]),
            Point3::from_array([0.75, 0.75, 3.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, tool).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        // 27 − 1 void − 0.25·0.25·0.5 pocket.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.96875).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
    }

    #[test]
    fn a_coplanar_boss_over_a_void_plane_is_rejected() {
        // The boss straddles the void's x=1 and y=1 planes (it spans [0.9,1.1]²), so classifying
        // the void's walls against it is no longer the clean whole-face case. Measured: honestly
        // rejected, before and after the all-shell fix — the fix does not widen this one. Pinned so
        // that a future change either keeps the reject or turns it into a correct 26.04, never into
        // a silent answer.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let boss = m.add_cuboid(
            Point3::from_array([0.9, 0.9, 3.0]),
            Point3::from_array([1.1, 1.1, 4.0]),
        );
        m.rebuild_adjacency();
        assert_eq!(
            boolean_one(&mut m, BoolKind::Fuse, hollow, boss),
            Err(BoolError::Unsupported)
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
    fn an_edge_slot_through_the_bottom() {
        // The prism pokes out the base's bottom too, so the old convex/blind overhang-cut gate
        // declined it and the seam path could not build it either (honest reject). The F2 dispatch
        // collapse hands it to the unified coplanar driver, which carves the slot exactly:
        // base 1.0 − (x∈[0.5,1] · y∈[0.25,0.75] · z∈[0,1]) = 1 − 0.25 = 0.75.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.25, -0.5]),
            Point3::from_array([1.5, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn a_corner_cut_through_the_bottom() {
        // The corner prism pokes out the base's bottom, so the old blind gate declined it. The F2
        // collapse routes it to the unified driver: base 1.0 − corner column (x,y ∈ [0.5,1], full
        // height) = 1 − 0.25 = 0.75.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.5, 0.5, -0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn a_boss_that_pierces_the_base_is_not_an_overhang() {
        // The boss dips below the base's top (its walls cross the base) — a transversal seam cut,
        // not a coplanar overhang. The old overhang detector declined it and the seam path built
        // it; after the F2 collapse the coplanar arm declines it and the same seam path still
        // does. Union = 1.0 + boss 0.75 − overlap 0.125 = 1.625.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 0.5]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.625).abs() < 1e-12, "volume {vol}");
    }

    // The Cut and Common twins of `fuse_a_corner_overhanging_boss` below — the same two solids,
    // the same shared z=1 plane. The boss sits entirely above it, so it removes nothing and shares
    // nothing: Cut is the base untouched and Common is empty. Both used to be rejected
    // (`coplanar_merge` for Cut; Common's every face dropped, which assembly reported as
    // `no_outward_shell`). The `Whole` survival cell now checks whether the contact plane actually
    // separates the solids, which is what makes the whole cap correct here.
    #[test]
    fn cut_by_a_corner_overhanging_boss_removes_nothing() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, base, corner).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
        // Structure, not just volume: the base comes through as itself. Six faces means the cap was
        // not split along ∂Q and the boss contributed nothing.
        let s = m.solids.get(r);
        assert!(s.cavities.is_empty());
        assert_eq!(m.shells.get(s.outer).faces.len(), 6, "a clean cube");
    }

    #[test]
    fn cut_a_seated_block_by_the_part_below_it() {
        // The operands swapped: now the canonical contact face is the upper block's *lower* cap, so
        // the separation test runs with the plane's normal the other way round. Same answer — the
        // block keeps its volume.
        let mut m = Model::new();
        let block = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, block, base).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn common_with_a_corner_overhanging_boss_is_empty() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        m.rebuild_adjacency();
        assert_eq!(
            boolean(&mut m, BoolKind::Common, base, corner),
            Err(BoolError::EmptyResult)
        );
    }

    #[test]
    fn cut_by_an_overhanging_boss_carrying_a_pin_is_rejected() {
        // Same seating, but the tool carries a pin reaching below the contact plane, so the plane
        // no longer separates the solids and the cut owes a real notch (1 − 0.2·0.2·0.5 = 0.98).
        // We cannot build that yet, and emitting the whole cap would be silently wrong — so this
        // pins the guard that keeps the widening above honest.
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let block = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        let pin = m.add_cuboid(
            Point3::from_array([0.55, 0.55, 0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        m.rebuild_adjacency();
        let tool = boolean_one(&mut m, BoolKind::Fuse, block, pin).unwrap();
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, tool).unwrap().volume - 1.02).abs() < 1e-12,
            "the pinned tool itself"
        );
        assert_eq!(
            boolean_one(&mut m, BoolKind::Cut, base, tool),
            Err(BoolError::Unsupported)
        );
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

    /// Every three-plane point plane class `p` puts on the line `p ∩ q`, gathered the way a
    /// per-plane-class arrangement would: both solids' faces seated on `p`, plus both solids'
    /// sections cut by `p`. Triples are canonicalized, then kept only if they name `q` as well —
    /// those are exactly the points on the shared line. `section_of_solid` rejecting (it declines a
    /// plane it has vertices on) contributes nothing rather than failing: a plane class with no
    /// section still has its faces.
    /// One boundary segment of a plane class, named entirely in plane classes. It rides `wall`, so
    /// it lies on the line `p ∩ wall`; on that line a point is named by its *third* plane alone
    /// (`order_along`'s convention, arrange.rs:204-214), and the segment's two ends are named by
    /// `end[0]` / `end[1]`.
    #[derive(Clone, Copy, Debug)]
    struct Chord {
        wall: usize,
        end: [usize; 2],
    }

    /// What the chord builder had to drop. Both are ways a class ends up with *less* boundary than
    /// it really has, and downstream neither is distinguishable from "the section was never
    /// computed" unless it is counted here.
    #[derive(Default, Debug, Clone, Copy)]
    struct ChordStats {
        /// A ring edge whose two endpoint triples do not share exactly one plane besides `p` —
        /// what a canonized four-plane vertex looks like. Produces no chord.
        skipped_naming: usize,
        /// A canonized triple with a repeated index (`[c, w, w]`), which is not a point.
        degenerate_triples: usize,
    }

    /// The arrangement vertices class `p` **mints**: for every pair of chords riding different
    /// walls, the three-plane point `{p, w1, w2}` when it falls strictly inside both chords.
    ///
    /// `order_along(planes, p, w, x, e)` orders two points of the line `p ∩ w`, each named by its
    /// third plane — so `x` is strictly between the chord's ends exactly when the two comparisons
    /// come back non-zero and opposite. **Strict is an invariant, not a choice:** admitting an
    /// endpoint mints a point that already exists, and a zero-length result edge panics at
    /// `Line::through_points(..).expect("distinct")` (lib.rs:3576) — an abort, not a reject.
    ///
    /// Walls whose normals are coplanar with `p`'s meet it in no point; `t_plane_pair_dir_sign`
    /// returns 0 there and the pair contributes nothing (that is the collinear-overlap case, E5).
    fn minted_crossings(planes: &[PlaneInfo], p: usize, chords: &[Chord]) -> HashSet<[usize; 3]> {
        let strictly_inside = |c: &Chord, x: usize| -> bool {
            let (s0, s1) = (
                arrange::order_along(planes, p, c.wall, x, c.end[0]),
                arrange::order_along(planes, p, c.wall, x, c.end[1]),
            );
            s0 != 0 && s1 != 0 && s0 != s1
        };
        let mut out = HashSet::new();
        for (i, c1) in chords.iter().enumerate() {
            for c2 in &chords[i + 1..] {
                if c1.wall == c2.wall
                    || tolerant::t_plane_pair_dir_sign(planes, p, c1.wall, c2.wall) == 0
                {
                    continue;
                }
                if strictly_inside(c1, c2.wall) && strictly_inside(c2, c1.wall) {
                    let mut t = [p, c1.wall, c2.wall];
                    t.sort_unstable();
                    out.insert(t);
                }
            }
        }
        out
    }

    /// One fixture's verdict: how the two sides of every intersecting class pair compare on their
    /// shared line, once crossings are minted.
    #[derive(Default, Debug)]
    struct SweepTally {
        pairs: usize,
        agreed: usize,
        outside: usize, // disagreement, but the missing point is outside every host face
        /// Inside a host face, but the lacking class's section producer declined — it never had
        /// the input, so this is a `section_of_solid` limitation, not a minting disagreement.
        inside_missing_section: usize,
        /// Inside a host face with both sections computed. **This is the one that would refute
        /// per-plane minting**, and the assertion below pins it at zero.
        inside_unexplained: usize,
        /// **The dominant category, and the increment's real finding.** `point_in_ring` declines
        /// because the point lies *on* the host face's boundary — neither strictly in nor out. Far
        /// from undecidable noise, this is the case that matters most: one class puts a vertex on
        /// the shared line where the other class's boundary passes through with no vertex at all.
        /// That is a T-junction, the thing `resplit_overhang` exists to repair in today's engine.
        /// On the boundary, and the lacking class had **both sections computed** — the only
        /// sub-population that can speak about minting rather than about absent input.
        on_boundary_complete: usize,
        /// On the boundary, but the lacking class's section producer declined. Same confound that
        /// accounted for 42/42 of the inside-disagreements; measured here rather than modelled.
        on_boundary_missing_section: usize,
        /// Boundary a class could not build: four-plane naming skips and canon-degenerate triples.
        /// Both look like absent input downstream, so they are surfaced, not folded in.
        chord_skips: usize,
        chord_degenerate: usize,
        /// Containment genuinely undecidable and not on any boundary. Counted, never read as
        /// "outside" — reading a reject as a negative is a mistake this work already made once.
        undecided: usize,
    }

    #[allow(clippy::too_many_arguments)]
    fn on_a_host_ring(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Face>, usize>,
        inc_a: &arrange::EdgePlanes,
        inc_b: &arrange::EdgePlanes,
        canon: &[usize],
        c: usize,
        t: [usize; 3],
    ) -> bool {
        let canonize = |x: [usize; 3]| {
            let mut v = [canon[x[0]], canon[x[1]], canon[x[2]]];
            v.sort_unstable();
            v
        };
        for (solid, inc) in [(a, inc_a), (b, inc_b)] {
            for sh in solid_shell_handles(m, solid) {
                for &fh in &m.shells.get(sh).faces {
                    let fi = surf_ix[&fh];
                    if canon[fi] != c {
                        continue;
                    }
                    let Ok(ts) = arrange::face_vertex_triples(m, fh, fi, inc) else {
                        continue;
                    };
                    let mut rings: Vec<Vec<[usize; 3]>> =
                        vec![ts.into_iter().map(canonize).collect()];
                    // Hole rims bound the face too — a point on a hole rim is on the boundary.
                    if let Ok(hs) = arrange::hole_rings(m, fh, fi, inc) {
                        rings.extend(
                            hs.into_iter()
                                .map(|r| r.into_iter().map(canonize).collect()),
                        );
                    }
                    for ring in &rings {
                        if arrange::point_on_ring(planes, c, t, ring) == Ok(true) {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// Is the three-plane point `t` inside any face either solid seats on class `c`? Rejections are
    /// reported separately rather than folded into "outside" — reading a reject as a negative is a
    /// mistake this line of work has already made once.
    #[allow(clippy::too_many_arguments)]
    fn inside_a_host_face(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        planes: &[PlaneInfo],
        surf_ix: &HashMap<Handle<Face>, usize>,
        inc_a: &arrange::EdgePlanes,
        inc_b: &arrange::EdgePlanes,
        canon: &[usize],
        c: usize,
        t: [usize; 3],
    ) -> Result<bool, BoolError> {
        let canonize = |x: [usize; 3]| {
            let mut v = [canon[x[0]], canon[x[1]], canon[x[2]]];
            v.sort_unstable();
            v
        };
        let mut rejected = None;
        for (solid, inc) in [(a, inc_a), (b, inc_b)] {
            for sh in solid_shell_handles(m, solid) {
                for &fh in &m.shells.get(sh).faces {
                    let fi = surf_ix[&fh];
                    if canon[fi] != c {
                        continue;
                    }
                    let Ok(ts) = arrange::face_vertex_triples(m, fh, fi, inc) else {
                        continue;
                    };
                    let ring: Vec<[usize; 3]> = ts.into_iter().map(canonize).collect();
                    match arrange::point_in_ring(planes, c, t, &ring) {
                        Ok(true) => {
                            // Inside the outer loop, but a point inside a hole is outside the face.
                            let in_hole = arrange::hole_rings(m, fh, fi, inc)
                                .map(|hs| {
                                    hs.into_iter().any(|r| {
                                        let h: Vec<[usize; 3]> =
                                            r.into_iter().map(canonize).collect();
                                        arrange::point_in_ring(planes, c, t, &h) == Ok(true)
                                    })
                                })
                                .unwrap_or(false);
                            if !in_hole {
                                return Ok(true);
                            }
                        }
                        Ok(false) => {}
                        Err(e) => rejected = Some(e),
                    }
                }
            }
        }
        match rejected {
            Some(e) => Err(e),
            None => Ok(false),
        }
    }

    /// **Differential coverage sweep — the cutover roadmap + silent-wrong guard.**
    ///
    /// Runs the isolated trace engine (`crate::trace::boolean_via_trace`, unwired) against the
    /// OCCT-validated production `boolean` over the whole two-solid fixture corpus. The load-bearing
    /// invariant: wherever BOTH engines succeed, trace's result is manifold and its volume+area
    /// match production — else it is a trace silent-wrong. Everything trace declines is the recorded
    /// cutover roadmap (`trace_gap`), not a failure; a cell only prod declines is UNVERIFIED (no
    /// oracle). `Model` has no `Clone`, so each op rebuilds the fixture fresh (mirrors
    /// `rotation_invariance_stress`). Compared on `vclose` (prod's own harness bound) so f64
    /// accumulation order never fakes a disagreement.
    #[test]
    fn trace_vs_production_coverage_sweep() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        #[derive(Clone, Copy, Debug)]
        enum Out {
            Ok { vol: f64, area: f64, valid: bool },
            Err,
        }
        // `Err(())` = the engine panicked (caught). `Ok(Out::Err)` = an honest `Err(BoolError)`.
        let run = |build: &dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>),
                   kind: BoolKind,
                   trace: bool|
         -> Result<Out, ()> {
            catch_unwind(AssertUnwindSafe(|| {
                let (mut m, a, b) = build();
                let r = if trace {
                    crate::trace::boolean_via_trace(&mut m, kind, a, b)
                } else {
                    boolean(&mut m, kind, a, b)
                };
                match r {
                    Ok(solids) => {
                        m.rebuild_adjacency();
                        let vol = solids
                            .iter()
                            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                            .sum();
                        let area = solids
                            .iter()
                            .map(|&s| nacre_props::mass_props(&m, s).unwrap().area)
                            .sum();
                        let valid = nacre_validate::validate(&m).is_empty();
                        Out::Ok { vol, area, valid }
                    }
                    Err(_) => Out::Err,
                }
            }))
            .map_err(|_| ())
        };
        let vclose =
            |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());

        type Build = Box<dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>)>;
        let fixtures: Vec<(&str, Build)> = vec![
            (
                "same_ground",
                Box::new(|| {
                    let mut m = Model::new();
                    let a = m.add_cuboid(
                        Point3::from_array([0.0; 3]),
                        Point3::from_array([1.0, 1.0, 1.0]),
                    );
                    let b = m.add_cuboid(
                        Point3::from_array([0.5, 0.5, 0.0]),
                        Point3::from_array([1.5, 1.5, 1.0]),
                    );
                    (m, a, b)
                }),
            ),
            ("l_and_inner_box", Box::new(l_and_inner_box)),
            ("l_and_corner_box", Box::new(l_and_corner_box)),
            ("l_and_reflex_box", Box::new(l_and_reflex_box)),
            ("u_and_slab", Box::new(u_and_slab)),
            ("cube_and_notch", Box::new(cube_and_notch)),
            ("l_and_rod", Box::new(l_and_rod)),
            ("l_and_popup_box", Box::new(l_and_popup_box)),
            ("l_and_notch_bar", Box::new(l_and_notch_bar)),
            ("l_and_ell_stub", Box::new(l_and_ell_stub)),
            ("l_and_staple", Box::new(l_and_staple)),
            ("l_and_dimple", Box::new(l_and_dimple)),
            ("two_boxes", Box::new(two_boxes)),
            ("nested_boxes", Box::new(nested_boxes)),
        ];

        let (mut agree, mut trace_wrong, mut trace_gap, mut both_err, mut prod_only_err, mut panic) =
            (0, 0, 0, 0, 0, 0);
        let mut notes: Vec<String> = Vec::new();
        for (name, build) in &fixtures {
            for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                let prod = run(build.as_ref(), kind, false).expect("production boolean panicked");
                let tr = run(build.as_ref(), kind, true);
                let class = match (prod, tr) {
                    (_, Err(())) => {
                        panic += 1;
                        notes.push(format!("{name} {kind:?}: TRACE PANIC"));
                        "TRACE_PANIC"
                    }
                    (
                        Out::Ok {
                            vol: pv, area: pa, ..
                        },
                        Ok(Out::Ok {
                            vol: tv,
                            area: ta,
                            valid,
                        }),
                    ) => {
                        if valid && vclose(pv, tv) && vclose(pa, ta) {
                            agree += 1;
                            "agree"
                        } else {
                            trace_wrong += 1;
                            notes.push(format!(
                                "{name} {kind:?}: WRONG prod(v={pv},a={pa}) trace(v={tv},a={ta},valid={valid})"
                            ));
                            "TRACE_WRONG"
                        }
                    }
                    (Out::Ok { .. }, Ok(Out::Err)) => {
                        trace_gap += 1;
                        "trace_gap"
                    }
                    (Out::Err, Ok(Out::Ok { valid, .. })) => {
                        prod_only_err += 1;
                        notes.push(format!(
                            "{name} {kind:?}: prod-only-err, trace Ok (UNVERIFIED, valid={valid})"
                        ));
                        "prod_only_err"
                    }
                    (Out::Err, Ok(Out::Err)) => {
                        both_err += 1;
                        "both_err"
                    }
                };
                eprintln!("{name:18} {kind:?}: {class}");
            }
        }
        eprintln!(
            "=== agree {agree} · gap {trace_gap} · wrong {trace_wrong} · both_err {both_err} · prod_only_err {prod_only_err} · panic {panic} ==="
        );
        for n in &notes {
            eprintln!("  {n}");
        }

        assert_eq!(panic, 0, "trace panicked (an un-honest reject): {notes:?}");
        assert_eq!(trace_wrong, 0, "trace is silently wrong on: {notes:?}");
        // Measured 2026-07-21 after nest_cells learned multi-hole faces (dropping the HOLE_MULTI
        // reject — the U's two prongs are two holes in the slab's one face): agree 42, trace_gap 0,
        // wrong/panic/prod_only_err/both_err 0. u_and_slab graduated. **The whole corpus now agrees
        // with production** — the cutover threshold: every two-solid fixture × 3 ops matches on
        // volume+area+manifold+solid-count. Floor at the measured 42 (the maximum): any drop is a
        // regression.
        assert!(
            agree >= 42,
            "only {agree} agreements (was 42, the whole corpus); trace regressed a case: {notes:?}"
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
    fn overhang_boss_with_a_non_convex_footprint() {
        // An L-shaped (non-convex) boss footprint overhanging a cube edge. The old convexity-gated
        // overhang detector declined this, so it used to be an honest reject; the F2 dispatch
        // collapse hands it to the unified coplanar driver, which builds it exactly. Volume =
        // cube 1.0 + L-prism (area 0.9·0.2 + 0.3·0.2 = 0.24) · height 0.4 = 1.096.
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
        let r = boolean_one(&mut m, BoolKind::Fuse, cube, l_tool).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.096).abs() < 1e-12, "volume {vol}");
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
    fn a_pocket_that_punches_through_drills_a_bore() {
        // A top-flush tool that pokes out the base's bottom: the pocket becomes a through hole.
        // The exit face has to come out annular, and until the coplanar reconstruct learned to
        // emit a hole it came out whole instead, leaving the bore's walls nothing to close
        // against — an open shell the assembly guard rejected. (Honest reject, never a wrong
        // answer; the previous cell pinned it as such.)
        //
        // Area is the assertion that matters here: volume alone cannot tell a bore from a shape
        // that merely displaces the same material. 0.75 (top) + 0.75 (bottom) + 4 (sides) +
        // 2.0 (the bore's four inner walls) = 7.5, against 6.0 for the cube.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.25, 0.25, -0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, r).unwrap();
        assert!((p.volume - 0.75).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 7.5).abs() < 1e-12, "area {}", p.area);
        // A bore, not a void: no cavity shell, and both caps carry the hole (the top from the
        // coincident contact, the bottom from the section the tool cuts through it).
        let s = m.solids.get(r);
        assert!(s.cavities.is_empty(), "a through hole is not a cavity");
        let faces = &m.shells.get(s.outer).faces;
        assert_eq!(faces.len(), 10, "6 base faces + the bore's 4 walls");
        assert_eq!(
            faces
                .iter()
                .filter(|&&fh| !m.faces.get(fh).inner.is_empty())
                .count(),
            2,
            "both caps are annular"
        );
    }

    #[test]
    fn a_boss_that_punches_through_keeps_the_stub() {
        // The Fuse twin, and the same emission path: the tool's a-side face keeps `P∖Q`, so the
        // base's bottom needs the same hole for the stub below it to join on. Volume
        // 1 + 0.5·0.5·0.5 = 1.125; area 1.0 (top, the flush tool cap dissolves into it) + 0.75
        // (bottom) + 4 (sides) + 1.0 (stub walls) + 0.25 (stub floor) = 7.0.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.25, 0.25, -0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, r).unwrap();
        assert!((p.volume - 1.125).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 7.0).abs() < 1e-12, "area {}", p.area);
        let s = m.solids.get(r);
        assert_eq!(
            m.shells
                .get(s.outer)
                .faces
                .iter()
                .filter(|&&fh| !m.faces.get(fh).inner.is_empty())
                .count(),
            1,
            "only the bottom is annular — the flush top merges away"
        );
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

    // ---- `unify_coplanar_faces`: the general (interface-free) coplanar merge, on hand-built
    // `LocalFace` lists the coincident goldens above never reach — chains, opposite normals,
    // holes, seam edges, and the asymmetric T-junction the global dissolve exists to prevent. ----

    fn mk_vert(m: &mut Model, x: f64, y: f64, z: f64) -> Handle<Vertex> {
        m.vertices.push(Vertex {
            point: Point3::from_array([x, y, z]),
            origin: Origin::Constructed,
        })
    }

    /// A `PlaneInfo` whose only field `unify_coplanar_faces` reads is `n_out`; the rest is a
    /// valid-but-unreferenced dummy (`surf`/`face`/`plane` are never dereferenced there).
    fn mk_plane(m: &mut Model, n: [f64; 3]) -> PlaneInfo {
        let normal = Vector3::from_array(n);
        let plane = Plane::from_point_normal(Point3::origin(), normal).unwrap();
        let surf = m.surfaces.push(Surface::Plane(plane));
        let face = m.faces.push(Face {
            surface: surf,
            outer: Loop { half_edges: vec![] },
            inner: vec![],
            orientation: Orientation::Forward,
        });
        PlaneInfo {
            surf,
            face,
            plane,
            tri: [Point3::origin(); 3],
            n_out: normal,
            orient: Orientation::Forward,
            tri_pt3: None,
        }
    }

    fn oloop(vs: &[Handle<Vertex>]) -> Vec<Node> {
        vs.iter().map(|&v| Node::Orig(v)).collect()
    }

    fn face(plane_idx: usize, nodes: Vec<Node>, inner: Vec<Vec<Node>>) -> LocalFace {
        LocalFace {
            plane_idx,
            loop_nodes: nodes,
            inner,
            flip: false,
        }
    }

    #[test]
    fn unify_merges_a_coplanar_chain() {
        // Three unit squares on z=0 (+z), tiled in x, each sharing a vertical edge with the
        // next. One plane class, one normal ⇒ all fuse into a single face; the four
        // straight-angle mid-edge vertices dissolve, leaving one 4-corner rectangle.
        let mut m = Model::new();
        let p = vec![mk_plane(&mut m, [0.0, 0.0, 1.0])];
        let canon = vec![0];
        let g = |m: &mut Model, x, y| mk_vert(m, x, y, 0.0);
        let (c00, c10, c20, c30) = (
            g(&mut m, 0., 0.),
            g(&mut m, 1., 0.),
            g(&mut m, 2., 0.),
            g(&mut m, 3., 0.),
        );
        let (c31, c21, c11, c01) = (
            g(&mut m, 3., 1.),
            g(&mut m, 2., 1.),
            g(&mut m, 1., 1.),
            g(&mut m, 0., 1.),
        );
        let faces = vec![
            face(0, oloop(&[c00, c10, c11, c01]), vec![]),
            face(0, oloop(&[c10, c20, c21, c11]), vec![]),
            face(0, oloop(&[c20, c30, c31, c21]), vec![]),
        ];
        let out = unify_coplanar_faces(&m, faces, &p, &canon);
        assert_eq!(out.len(), 1, "three coplanar faces fuse into one");
        let l = &out[0].loop_nodes;
        assert_eq!(l.len(), 4, "straight-angle mid vertices dissolved: {l:?}");
        for c in [c00, c30, c31, c01] {
            assert!(l.contains(&Node::Orig(c)), "corner kept");
        }
        for c in [c10, c20, c11, c21] {
            assert!(!l.contains(&Node::Orig(c)), "mid vertex dropped");
        }
    }

    #[test]
    fn unify_keeps_opposite_normal_coplanar() {
        // Two coplanar faces sharing an edge but with opposite outward normals — a genuine
        // fold / cantilever step (cell canon-step), not redundant. The shared-normal guard
        // (`n_out.dot > 0`) keeps them separate.
        let mut m = Model::new();
        let p = vec![
            mk_plane(&mut m, [0., 0., 1.]),
            mk_plane(&mut m, [0., 0., -1.]),
        ];
        let canon = vec![0, 0]; // same plane class, opposite normal
        let a = mk_vert(&mut m, 0., 0., 0.);
        let b = mk_vert(&mut m, 1., 0., 0.);
        let c = mk_vert(&mut m, 1., 1., 0.);
        let d = mk_vert(&mut m, 0., 1., 0.);
        let h = mk_vert(&mut m, 0., -1., 0.);
        let e = mk_vert(&mut m, 1., -1., 0.);
        let faces = vec![
            face(0, oloop(&[a, b, c, d]), vec![]),
            face(1, oloop(&[b, a, h, e]), vec![]),
        ];
        let out = unify_coplanar_faces(&m, faces, &p, &canon);
        assert_eq!(out.len(), 2, "opposite-normal pair stays separate");
    }

    #[test]
    fn unify_keeps_holed_faces() {
        // A holed face is left separate (merging holes is deferred).
        let mut m = Model::new();
        let p = vec![mk_plane(&mut m, [0., 0., 1.])];
        let canon = vec![0];
        let g = |m: &mut Model, x, y| mk_vert(m, x, y, 0.0);
        let (a, b, c, d) = (
            g(&mut m, 0., 0.),
            g(&mut m, 1., 0.),
            g(&mut m, 1., 1.),
            g(&mut m, 0., 1.),
        );
        let (e, f) = (g(&mut m, 2., 0.), g(&mut m, 2., 1.));
        let (h0, h1, h2) = (
            g(&mut m, 0.2, 0.2),
            g(&mut m, 0.4, 0.2),
            g(&mut m, 0.3, 0.4),
        );
        let faces = vec![
            face(0, oloop(&[a, b, c, d]), vec![oloop(&[h0, h1, h2])]),
            face(0, oloop(&[b, e, f, c]), vec![]),
        ];
        let out = unify_coplanar_faces(&m, faces, &p, &canon);
        assert_eq!(out.len(), 2, "a holed face is not merged");
    }

    #[test]
    fn unify_skips_seam_shared_edges() {
        // Two coplanar same-normal faces sharing an edge whose endpoints are seam nodes.
        // Detection requires `Orig` endpoints (splice reads original vertices), so the pair is
        // left separate — and, crucially, `splice_along`'s `Orig`-only path is never entered
        // (no panic).
        let mut m = Model::new();
        let p = vec![mk_plane(&mut m, [0., 0., 1.])];
        let canon = vec![0];
        let a = Node::Orig(mk_vert(&mut m, 0., 0., 0.));
        let b = Node::Orig(mk_vert(&mut m, 1., 1., 0.));
        let (s1, s2) = (Node::Seam([1, 2, 3]), Node::Seam([2, 3, 4]));
        let faces = vec![
            face(0, vec![s1, s2, a], vec![]),
            face(0, vec![s2, s1, b], vec![]),
        ];
        let out = unify_coplanar_faces(&m, faces, &p, &canon);
        assert_eq!(out.len(), 2, "seam-shared edge is not merged");
    }

    #[test]
    fn unify_keeps_a_vertex_that_is_a_corner_elsewhere() {
        // F0,F1 on z=0 merge; their shared-edge endpoint (1,0,0) is a straight angle on the
        // merged face but a real corner on a perpendicular face G (plane y=0). Global degree 3
        // ⇒ it is NOT dissolved — a per-face local rule would have, opening a T-junction. Its
        // twin (1,1,0), on the merged face only, IS dissolved.
        let mut m = Model::new();
        let p = vec![
            mk_plane(&mut m, [0., 0., 1.]),
            mk_plane(&mut m, [0., -1., 0.]),
        ];
        let canon = vec![0, 1];
        let v000 = mk_vert(&mut m, 0., 0., 0.);
        let v100 = mk_vert(&mut m, 1., 0., 0.);
        let v200 = mk_vert(&mut m, 2., 0., 0.);
        let v210 = mk_vert(&mut m, 2., 1., 0.);
        let v110 = mk_vert(&mut m, 1., 1., 0.);
        let v010 = mk_vert(&mut m, 0., 1., 0.);
        let v201 = mk_vert(&mut m, 2., 0., 1.);
        let v101 = mk_vert(&mut m, 1., 0., 1.);
        let faces = vec![
            face(0, oloop(&[v000, v100, v110, v010]), vec![]),
            face(0, oloop(&[v100, v200, v210, v110]), vec![]),
            face(1, oloop(&[v200, v100, v101, v201]), vec![]), // perpendicular, not coplanar
        ];
        let out = unify_coplanar_faces(&m, faces, &p, &canon);
        assert_eq!(out.len(), 2, "z=0 pair merges; G stays");
        let merged = out.iter().find(|lf| lf.plane_idx == 0).unwrap();
        assert!(
            merged.loop_nodes.contains(&Node::Orig(v100)),
            "corner-elsewhere vertex kept (no T-junction)"
        );
        assert!(
            !merged.loop_nodes.contains(&Node::Orig(v110)),
            "pure straight-angle vertex dropped"
        );
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
