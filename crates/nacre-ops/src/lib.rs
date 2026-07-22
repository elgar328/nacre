//! Operations for the nacre kernel, plus a replayable operation log (design §6).
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::PadOnFace`]/[`Operation::PocketOnFace`] (M4) consume a prior op's face by
//! `Handle` (exposed via [`OpOutput`]) and supersede a solid (design §2 live-solid
//! semantics) — each is a tool prism plus a boolean, not a direct face-split. Ops are
//! applied by [`apply`] and folded by [`replay`]; every result is a **closed** solid, so
//! `nacre-validate` applies fully.

use nacre_geom::{Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::Isometry;
use nacre_scalar::frame3::{Pt3, dir_orient3d_judge};
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
}

/// Names of the `Unsupported` reject sites, shared by the guard that raises one
/// and the test that asserts it — a renamed tag then cannot silently drift out
/// of a test's expectation. Not `#[cfg(test)]`: the guards name these in release
/// builds too. They are `const`, so they inline away where the tag is unused.
pub(crate) mod tag {
    /// An outer-shell edge used by other than two face loops. A backstop with no
    /// firing test: `validate` calls this `NonOpposedEdge` and every shell the
    /// operations build is manifold — but `boolean` never runs `validate` on its
    /// inputs, so a direct caller could still hand one in.
    pub const NON_MANIFOLD_EDGE: &str = "non_manifold_edge";
    /// A closed seam loop's edges disagree about which side its material lies on, or
    /// two of its nodes coincide, or two neighbours share no plane pair.
    ///
    /// A cross-check rather than a defence: the exact
    /// predicate (`order_along`) and the orientation bookkeeping (`PlaneInfo::n_out`
    /// vs its plane's normal) must agree, edge by edge, and neither is assumed right.
    /// It survives release on purpose. Downstream, `validate` catches a wrong loop and
    /// `tessellate` catches a wrong *hole* — but a caller may run neither.
    ///
    /// It cannot catch a *globally* flipped loop: every edge would be wrong together.
    /// A globally flipped loop would need a golden on the loop producer itself.
    ///
    /// Unreachable today: "material on the left" is a global property, so a consistent
    /// loop makes every edge agree. Not dead code — relax `FOURPLANE` or cell 3c's
    /// node-identity argument and this is what speaks first.
    pub const LOOP_ORIENT_MISMATCH: &str = "loop_orient_mismatch";
    pub const COPLANAR_PAIR: &str = "coplanar_pair";
    /// The result severs into two or more material solids *and* at least one enclosed void
    /// (cavity) survives. Which outer shell owns which cavity needs a shell-scoped point-in-shell
    /// test we do not have yet, so this is honestly rejected and deferred to a follow-on cell.
    /// Reachable: `Cut` a hollow part with a cut that isolates the void into one severed piece.
    /// Born with its firing test (`severed_with_cavity_is_rejected`).
    ///
    /// ★ The starting point for that cell is `arrange::point_in_solid_idx`, removed by the
    /// seam-engine excise (2026-07-22 dev-log cell) once the arrangement stopped needing it (it seeds the unbounded cell and propagates inward
    /// instead). It answers *point in solid* by ray casting; what this needs is *point in shell*,
    /// so recover it from that commit and re-scope it rather than deriving one from scratch.
    pub const SEVERED_WITH_CAVITY: &str = "severed_with_cavity";
    /// No material-enclosing (outward) shell among the result components — every component is
    /// inward-oriented. Geometrically impossible for a real solid result; a defensive backstop
    /// with no firing test (cf. `FOURPLANE`).
    pub const NO_OUTWARD_SHELL: &str = "no_outward_shell";
    /// A boolean input is a rotated solid (overhaul stage 1b). Rotated planar geometry
    /// is representable but its predicates are not yet sound (no TIP until stage 3), so
    /// the boolean honestly rejects until then. Fired by a `Transform`-rotated operand.
    pub const ROTATED_UNSUPPORTED: &str = "rotated_unsupported";
    pub const THREE_PLANES: &str = "three_planes";
    /// Two **different** arrangement vertices (distinct plane triples) materialized to the same
    /// coordinate. The triple is the truth and the coordinate only its cache (overview §5), so this
    /// says the exact substrate and the f64 cache disagree about how many vertices exist — always a
    /// defect upstream, never a property of the input. Raised where the seam table is built, while
    /// both triples are still in hand; without it the disagreement surfaces much later as a
    /// zero-length edge. The known cause is a **split plane table** (one geometric plane carried by
    /// two classes); a genuine 4-plane concurrency would do the same.
    pub const SEAM_ALIAS: &str = "seam_alias";
    /// A result loop asked for an edge between two vertices at the same coordinate. Every ring node
    /// is a distinct arrangement vertex, so this cannot happen for well-named input — it is the
    /// backstop that keeps a degenerate one from aborting the kernel (`Line::through_points` used to
    /// `expect`). `SEAM_ALIAS` catches the known cause earlier, so this has no firing test.
    pub const ZERO_LENGTH_EDGE: &str = "zero_length_edge";
    pub const FOURPLANE: &str = "fourplane";
    pub const CYLINDER_FACE: &str = "cylinder_face";
    pub const DEGENERATE_FACE: &str = "degenerate_face";
    pub const DEGENERATE_NORMAL: &str = "degenerate_normal";
    /// Every candidate ray from a loop's nodes has a ring node on its line.
    ///
    /// `point_in_ring` casts along `P ∩ Q_a` for a node's own plane `Q_a`; a ring node on
    /// that line makes the crossing parity ambiguous. Candidates are `2 · |loop|` lines and
    /// two directions, and half of them can be spoiled at once — `l_and_staple`'s loop and
    /// arc share both `y` planes, so only the `x` lines are clear there. Unfired today.
    pub const NO_CLEAR_RAY: &str = "no_clear_ray";
    /// A loop's node lies *on* the ring it is being tested against.
    ///
    /// A hole ring never touches the outer ring it sits in, and `point_in_ring` checks that
    /// exactly: the ray's line meets an edge at `X`, and `X == v` strictly inside that edge means
    /// `v` is on the ring. Unfired.
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
    /// A `Whole`-survival contact face whose footprint OVERLAPS the other's (∂P × ∂Q cross) rather
    /// than nesting, in the one such case still unbuilt. `Whole` has two entries: `Fuse`/same-normal,
    /// which the E1 union cell now builds, and `Cut`/opposite-normal, which is exact whenever the
    /// contact plane separates the two solids (nothing to remove). What is left is a `Cut` whose tool
    /// reaches back across that plane — a pin below its own contact face — where the cut owes a notch
    /// this path cannot yet cut. Honest reject rather than a whole cap that ignores the pin.
    pub const COPLANAR_MERGE: &str = "coplanar_merge";
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
    /// A pad/pocket target face is not planar (only planar faces carry a sketch frame;
    /// curved-face features arrive with the quadric milestones).
    NonPlanarFace,
    /// A pad/pocket target face belongs to no live solid's outer shell (a stale
    /// or non-live handle).
    FaceNotInLiveSolid,
    /// A pocket's depth reaches through the solid: the carved prism is not blind, so `Cut`
    /// produced a through-hole with no floor face. `pocket` requires `dist` less than the
    /// thickness at the face (the boolean pocket path honestly rejects instead of the old
    /// direct path's silent invalid result).
    PocketNotBlind,
    /// A pad's footprint does not meet the face at all: the `Fuse` came back severed, which two
    /// one-shell solids can only do if they never touched. Like [`OpError::PocketNotBlind`] this is
    /// the *operation's* premise breaking, not a boolean failure — the boolean answered correctly
    /// (a base and a detached boss). Use `Operation::Boolean` directly if two disjoint solids are
    /// what you want. An overhanging footprint still touches and is not this error.
    PadMissesFace,
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
/// Every half-edge walked here belongs to a valid solid, so its edge is bounded.
pub(crate) fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
    let [a, b] = model
        .edges
        .get(he.edge)
        .bounds
        .expect("a solid's loop edge is bounded");
    if he.forward { a } else { b }
}

/// A planar face's live solid, its in-plane right-handed frame (`x × y = n`, centred on the face
/// centroid so a profile's `(0,0)` lands there), and its loops — the shared setup for placing a
/// profile on a face (pad / pocket).
struct FaceFrame {
    solid_h: Handle<Solid>,
    surface_h: Handle<Surface>,
    n: Vector3, // outward normal
    x: Vector3,
    y: Vector3,
    origin: Point3, // face centroid
}

/// Locate `face`'s live solid and build its planar frame. `NonPlanarFace` for a curved surface,
/// `FaceNotInLiveSolid` if no live outer shell holds it.
fn face_frame(model: &Model, face: Handle<Face>) -> Result<FaceFrame, OpError> {
    let (solid_h, _) = model
        .live_solids
        .iter()
        .map(|&s| (s, model.solids.get(s).outer))
        .find(|&(_, sh)| model.shells.get(sh).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;
    let f = model.faces.get(face);
    let surface_h = f.surface;
    let orientation = f.orientation;
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
    let outer_pts: Vec<Point3> = f
        .outer
        .half_edges
        .iter()
        .map(|he| model.vertices.get(he_start(model, *he)).point)
        .collect();
    let origin = Point3::centroid(&outer_pts).ok_or(OpError::DegenerateGeometry)?;
    Ok(FaceFrame {
        solid_h,
        surface_h,
        n,
        x,
        y,
        origin,
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
    let (prism, prism_faces) = build_prism(model, &base_pts, n * signed, Some(frame.surface_h))?;
    let solids = boolean(model, kind, frame.solid_h, prism).map_err(|e| {
        model.live_solids.retain(|&s| s != prism); // drop the transient prism (atomic on failure)
        OpError::Boolean(e)
    })?;
    // A pad's prism must actually meet the face. Two live solids are each one outer shell, so a
    // `Fuse` of them can only come back severed if they never touched — the footprint missed the
    // face entirely. Returning the piece that carries the cap would hand back a floating boss and
    // silently drop the base, so this is `PadMissesFace`: the boolean succeeded and answered
    // correctly (two solids); it is the *pad's* premise that broke. A footprint that merely
    // overhangs still touches, fuses into one solid, and takes the normal path.
    //
    // `assemble_fuse_cut` already retired the inputs, so restoring `live_solids` is what keeps the
    // "no reject-after-commit" contract true from the outside: the model the caller sees is the one
    // it had before. Only `live_solids` is touched — the store stays append-only.
    if matches!(kind, BoolKind::Fuse) && solids.len() > 1 {
        model.live_solids.retain(|s| !solids.contains(s));
        model.live_solids.push(frame.solid_h);
        return Err(OpError::PadMissesFace);
    }
    // Exposed cap = the result face on the prism's far-cap plane (face plane offset by n·signed),
    // its outward normal +n (the opening side for a pocket, the boss top for a boss). A pocket (Cut)
    // that severs leaves several solids — scan them all for the cap and return the piece that
    // carries it, leaving the others live; that is a valid multi-solid model, not a failure.
    // `build_prism` returns the far cap as `faces[1]`; the store is append-only, so it is still
    // readable after the boolean retired the prism, and it names the cap plane exactly.
    let far_cap = prism_faces[1];
    match solids
        .iter()
        .find_map(|&s| find_face_coplanar_with(model, s, far_cap, n).map(|c| (s, c)))
    {
        Some((solid, cap)) => Ok((solid, Some(cap))),
        // Nothing carries the cap. If anything survived at all, hand it back capless and let
        // `pad`/`pocket` decide; if the boolean came back empty the prism removed the whole solid,
        // which is `PocketNotBlind` taken to its limit — not merely floorless, but nothing left.
        // Only `Cut` can empty a result: `Fuse` of two non-empty solids is never empty.
        None => match solids.first() {
            Some(&primary) => Ok((primary, None)),
            None => {
                debug_assert!(
                    matches!(kind, BoolKind::Cut),
                    "a Fuse cannot produce an empty result"
                );
                Err(OpError::PocketNotBlind)
            }
        },
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

/// The outer-shell face of `solid` that lies on `reference`'s plane with its outward normal on
/// `want`'s side — how a pad/pocket recovers its own exposed cap (the boss top, the pocket floor)
/// from the boolean result. `None` if there is none (a through-pocket has no floor).
///
/// **`reference` is a real face, not a `(point, normal)` pair, and that is the point.** Naming the
/// plane by coefficients meant comparing `d = −n·origin` computed at *different* points of the same
/// plane: exact only when the dot happens to reproduce bit for bit, which for an axis-aligned frame
/// it does (`n·p` is one coordinate) and for a slanted one it does not. With a face in hand the
/// question is answered the way the kernel answers identity everywhere else:
///
/// 1. **the same `Surface` handle** — integers, not coordinates (overview §2). `assemble_fuse_cut`
///    gives a result face the surface of the operand plane it came from, so the surviving cap
///    normally lands here.
/// 2. **the faces' own coordinates, exactly** — every `outer_tri` point of the candidate lies on
///    `reference`'s tri plane (`plane_side`, an exact `orient3d` on the points the user gave).
///    Needed because `plane_idx` names a *class representative*: if the cap plane merged with a
///    coplanar face of the other operand, the survivor can carry that operand's surface instead.
///    A `None` from `outer_tri` (no non-collinear triple) means no evidence — that face is skipped,
///    and a degenerate `reference` leaves only branch 1.
///
/// The direction filter reads the **candidate's** outward normal against `want`, never
/// `reference`'s: a pocket's tool cap faces along the sweep (`−n`) while the floor it becomes faces
/// back into the void (`+n`). Coplanarity is settled by then, so the two are parallel and the dot
/// is a full magnitude away from zero — an f64 read whose sign cannot round the wrong way.
///
/// If the cap survives as several faces they all satisfy this, and the first is returned; the
/// coefficient test had the same ambiguity.
///
/// **Measured (2026-07-22): the corpus does not separate the two branches** — disabling either one
/// leaves the whole suite at 207 passed / 23 failed. So branch 2 has no firing test today and is a
/// documented backstop (cf. `NON_MANIFOLD_EDGE`); branch 1 is kept because handle identity is the
/// strongest answer available and is the path a surviving cap normally takes.
fn find_face_coplanar_with(
    model: &Model,
    solid: Handle<Solid>,
    reference: Handle<Face>,
    want: Vector3,
) -> Option<Handle<Face>> {
    let ref_surf = model.faces.get(reference).surface;
    let ref_tri = outer_tri(model, model.faces.get(reference)).map(|(tri, _)| tri);
    let shell = model.solids.get(solid).outer;
    model.shells.get(shell).faces.iter().copied().find(|&fh| {
        let Some((tri, _)) = outer_tri(model, model.faces.get(fh)) else {
            return false;
        };
        let coplanar = model.faces.get(fh).surface == ref_surf
            || ref_tri.is_some_and(|r| tri.iter().all(|&q| plane_side(r, q) == 0));
        coplanar && (tri[1] - tri[0]).cross(tri[2] - tri[0]).dot(want) > 0.0
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

use nacre_geom::intersect::{plane_plane, plane_side, planes_coplanar, three_planes};
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
    /// The **plane class** this face's plane belongs to — `plane_classes`' root index, filled by
    /// `plane_index_setup` once the classes are known; `usize::MAX` until then (and in a hand-built
    /// `PlaneInfo`, where no class table exists).
    ///
    /// `planes` is a **per-face** table while the arrangement reasons **per plane**, so the same
    /// `usize` means "face" in one place and "plane" in another. Four bugs came from comparing the
    /// two: an index that names a plane must be a class root, and this field is what lets a
    /// predicate *check* that (`debug_assert!(planes[i].class == i)`) instead of trusting it.
    pub(crate) class: usize,
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
                class: usize::MAX, // filled by `plane_index_setup` once the classes exist
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
    /// about the face's outward normal. Only the non-convex path ever fills this.
    inner: Vec<Vec<Node>>,
    flip: bool,
}

/// The minimal per-op plane table two solids share: the
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
    let canon = fill_classes(&mut planes);
    Ok((planes, surf_ix, inc_a, inc_b, canon))
}

/// Union-find the coplanar faces into plane classes, and tell each `PlaneInfo` which class it is on.
///
/// Coplanar planes → one class, so no exact predicate ever forms a `det=0` triple from two
/// coplanar faces meeting a vertex (the cantilever-step degeneracy).
///
/// The classes exist only once the whole table is assembled, so this is where a `PlaneInfo` learns
/// which plane it is on. Every index that *names a plane* (a triple's element, `wc`, a wall, a
/// predicate argument) must be one of these roots; `class` is what lets a consumer assert that
/// rather than assume it. **Anything that builds a plane table calls this** — a table whose `class`
/// is left unfilled is invisible to that check.
pub(crate) fn fill_classes(planes: &mut [PlaneInfo]) -> Vec<usize> {
    let canon = plane_classes(planes);
    for (i, pi) in planes.iter_mut().enumerate() {
        pi.class = canon[i];
    }
    canon
}

/// Whether three `Pt3` are **exactly collinear** (zero-area triangle), decided on their
/// pre-rotation rational `base` coordinates. A rigid rotation preserves collinearity, and three
/// vertices of one solid share a rotation chain, so their bases are comparable; all three
/// coordinate-plane projections of `(b−a)×(c−a)` must vanish (exact `Rat`, no tolerance). An
/// i128 overflow returns `false` (treat as non-collinear): a genuinely-collinear triangle then
/// stays and is at worst rejected `RAY_DEGENERATE`, never falsely skipped (which would drop a
/// real crossing — silent-wrong). Unrotated vertices carry `base == coord`, so this is the exact
/// zero-area (collinear) test on the vertices' rotation definitions.
#[cfg(test)]
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
    // No faces means no result — `Common` of two solids that miss each other, `Cut` of a box that
    // is wholly inside what cuts it. That is an answer, not a failure: a solid is bounded by faces,
    // so a non-empty result cannot have none. The inputs are still consumed, exactly as they are on
    // any other successful boolean — the retire below sits inside the `positives` match, which this
    // early return skips, so it has to happen here too or the operands stay live.
    if faces.is_empty() {
        model.live_solids.retain(|&s| s != a && s != b);
        return Ok(Vec::new());
    }
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
    // A ring's two consecutive nodes are distinct arrangement vertices, so their points differ and
    // the line through them exists. Reject rather than panic if it does not: an aborting kernel is
    // below the floor (`overview.md`: out-of-coverage input declines honestly). `SEAM_ALIAS`
    // catches the known way this happens — two triples on one point — at the seam table, where the
    // names are still in hand, so this is a backstop with no firing test (cf. `NON_MANIFOLD_EDGE`).
    let mut edge_for = |model: &mut Model,
                        va: Handle<Vertex>,
                        vb: Handle<Vertex>|
     -> Result<Handle<Edge>, BoolError> {
        let key = unordered(va.index() as usize, vb.index() as usize);
        if let Some(&e) = edge_of.get(&key) {
            return Ok(e);
        }
        let pa = model.vertices.get(va).point;
        let pb = model.vertices.get(vb).point;
        let line = Line::through_points(pa, pb).ok_or_else(|| reject(tag::ZERO_LENGTH_EDGE))?;
        let curve = model.curves.push(Curve::Line(line));
        let e = model.edges.push(Edge {
            curve,
            bounds: Some([va, vb]),
            origin: Origin::Constructed,
        });
        edge_of.insert(key, e);
        Ok(e)
    };

    let mut face_handles = Vec::new();
    for lf in faces {
        let mut ring = |model: &mut Model, nodes: &[Node]| -> Result<Loop, BoolError> {
            let handles: Vec<Handle<Vertex>> = nodes.iter().map(|nd| vh[nd]).collect();
            let k = handles.len();
            let mut half_edges: Vec<HalfEdge> = (0..k)
                .map(|t| {
                    let (va, vb) = (handles[t], handles[(t + 1) % k]);
                    let e = edge_for(model, va, vb)?;
                    let forward = model.edges.get(e).bounds.expect("bounded")[0] == va;
                    Ok(HalfEdge { edge: e, forward })
                })
                .collect::<Result<Vec<_>, BoolError>>()?;
            if lf.flip {
                // Cut's inside-A B-pieces: reverse every loop and toggle the
                // orientation below, so the outward normal points into the removed
                // region and each loop still keeps material on its left.
                half_edges.reverse();
                for he in &mut half_edges {
                    he.forward = !he.forward;
                }
            }
            Ok(Loop { half_edges })
        };
        let outer = ring(model, &lf.loop_nodes)?;
        let inner: Vec<Loop> = lf
            .inner
            .iter()
            .map(|h| ring(model, h))
            .collect::<Result<Vec<_>, BoolError>>()?;
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
fn shares_or_coplanar(planes: &[PlaneInfo], i: usize, j: usize) -> bool {
    let (pa, pb) = (&planes[i], &planes[j]);
    // Three independent witnesses, OR-ed, so this can only ever merge *more* than before:
    //  1. the same `Surface` handle — coplanar by reference (what an ops-built tool's base cap and
    //     its target face share, and what a chained operand's split coplanar faces share);
    //  2. exactly proportional coefficients — the original test, kept;
    //  3. the faces' own coordinates, exactly (`t_planes_coplanar`) — the only one of the three
    //     that does not read a *derived* value, and the one that catches two independently built
    //     solids whose walls coincide (`add_cuboid` stacked on `add_cuboid`), where the rounded
    //     coefficients of differently-sized faces are not exactly proportional.
    pa.surf == pb.surf
        || planes_coplanar(&pa.plane, &pb.plane)
        || tolerant::t_planes_coplanar(planes, i, j)
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
            if shares_or_coplanar(planes, i, j) {
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
    if node_rank(a) <= node_rank(b) {
        (a, b)
    } else {
        (b, a)
    }
}

/// A total order on `Node`. There is no `Ord` derive — `Orig` and `Seam` carry different shapes —
/// but every iteration over nodes has to be deterministic for replay, so they are ranked by hand.
fn node_rank(n: Node) -> (u8, usize, usize, usize) {
    match n {
        Node::Orig(v) => (0, v.index() as usize, 0, 0),
        Node::Seam(t) => (1, t[0], t[1], t[2]),
    }
}

/// Merge every group of coplanar, same-facing result faces into one face per connected piece, then
/// dissolve straight-angle vertices — the mandatory post-boolean defeature.
///
/// **Erase the interior, do not stitch the exterior.** When two faces become one region the boundary
/// between them stops being a boundary, so the merge is: take every directed ring edge in the group,
/// drop the ones that appear as an opposed pair (`a→b` together with `b→a`), and re-thread what is
/// left. Nothing has to be spliced, which is what lets one rule cover every way the pieces can meet
/// — sharing an edge, a hole filled exactly by a neighbour, a hole filled by *several* neighbours,
/// and any chain of those (one erase settles them all at once). The earlier version stitched loops
/// with `splice_along` and so had to special-case "exactly two hole-free faces across one edge",
/// leaving `// holed — deferred` for the rest; a flush tool cap then stayed two faces forever.
///
/// A group is one plane class, one outward direction, one `flip` — mixing any of those would fold
/// material the wrong way — split further into **edge-connected components**, because faces that
/// merely lie on the same plane without touching must each survive on their own.
///
/// Reused, not reinvented: `canon` for coplanarity, `n_out.dot > 0` for facing,
/// [`arrange::loop_winding`] to tell an outer ring from a hole, [`arrange::point_in_ring`] to give
/// each hole its owner — the same two exact predicates `trace::nest_cells` uses for the same
/// question, neither of which reads a coordinate.
///
/// Rejects rather than guesses: a directed edge appearing twice the same way (two faces claiming the
/// same side), an undirected edge on three or more rings (non-manifold in the plane), or a node with
/// two outgoing edges after erasure (pieces meeting at a single point, where the cycle is not
/// unique). A group containing an `Orig` node is passed through untouched — `loop_winding` needs the
/// plane triples, and the arrangement emits `Seam` for everything, so this is unreachable in
/// production.
pub(crate) fn unify_coplanar_faces(
    faces: Vec<LocalFace>,
    planes: &[PlaneInfo],
    canon: &[usize],
) -> Result<Vec<LocalFace>, BoolError> {
    let n = faces.len();
    let eligible = all_seam;
    // One plane class, one outward direction, one flip.
    let group_key = |lf: &LocalFace| -> (usize, bool, bool) {
        let c = canon[lf.plane_idx];
        let facing = planes[lf.plane_idx].n_out.dot(planes[c].n_out) > 0.0;
        (c, facing, lf.flip)
    };

    // Edge-connected components within a group.
    let mut comp: Vec<usize> = (0..n).collect();
    let mut carriers: HashMap<(Node, Node), Vec<usize>> = HashMap::new();
    for (fi, lf) in faces.iter().enumerate() {
        if !eligible(lf) {
            continue;
        }
        // Holes count here too: a tool cap sitting flush inside another face touches it only
        // along that hole, so leaving `inner` out would put the two in different components and
        // nothing would merge at all.
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            for (a, b) in ring_edges(ring) {
                carriers.entry(norm_edge(a, b)).or_default().push(fi);
            }
        }
    }
    for fs in carriers.values() {
        for w in fs.windows(2) {
            if group_key(&faces[w[0]]) == group_key(&faces[w[1]]) {
                let (ri, rj) = (uf_find(&mut comp, w[0]), uf_find(&mut comp, w[1]));
                if ri != rj {
                    comp[ri] = rj;
                }
            }
        }
    }

    // Group members by component, in face order so the result is replay-stable.
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (fi, lf) in faces.iter().enumerate() {
        if eligible(lf) {
            let r = uf_find(&mut comp, fi);
            members[r].push(fi);
        }
    }

    let mut merged: Vec<LocalFace> = Vec::new();
    let mut kept: Vec<Option<LocalFace>> = faces.into_iter().map(Some).collect();
    for mem in &members {
        if mem.len() < 2 {
            continue; // nothing to merge; the face (if any) is emitted as-is below
        }
        let group: Vec<&LocalFace> = mem
            .iter()
            .map(|&fi| kept[fi].as_ref().expect("member present"))
            .collect();
        let rings = merge_component(&group, planes, canon)?;
        let (plane_idx, flip) = (group[0].plane_idx, group[0].flip);
        merged.extend(rings.into_iter().map(|(outer, inner)| LocalFace {
            plane_idx,
            loop_nodes: outer,
            inner,
            flip,
        }));
        for &fi in mem {
            kept[fi] = None;
        }
    }
    let mut out: Vec<LocalFace> = kept.into_iter().flatten().collect();
    out.extend(merged);
    dissolve_straight_angles(&mut out);
    Ok(out)
}

/// Every node of every ring is a `Seam` triple — `loop_winding` and `point_in_ring` name their
/// arguments by plane, so an `Orig` node has nothing to give them.
fn all_seam(lf: &LocalFace) -> bool {
    lf.loop_nodes
        .iter()
        .chain(lf.inner.iter().flatten())
        .all(|nd| matches!(nd, Node::Seam(_)))
}

/// A ring's directed edges, `i → i+1` around.
fn ring_edges(ring: &[Node]) -> impl Iterator<Item = (Node, Node)> + '_ {
    (0..ring.len()).map(move |i| (ring[i], ring[(i + 1) % ring.len()]))
}

/// The `[usize; 3]` form a ring's nodes carry, for the exact predicates.
fn seam_ring(ring: &[Node]) -> Vec<[usize; 3]> {
    ring.iter()
        .map(|nd| match nd {
            Node::Seam(t) => *t,
            Node::Orig(_) => unreachable!("callers filter on `all_seam`"),
        })
        .collect()
}

/// An outer ring with the holes that belong to it — what one merged region looks like before it
/// becomes a `LocalFace`.
type RegionRings = (Vec<Node>, Vec<Vec<Node>>);

/// One edge-connected group → its faces after erasing the interior boundary: each outer ring with
/// the holes that belong to it.
fn merge_component(
    group: &[&LocalFace],
    planes: &[PlaneInfo],
    canon: &[usize],
) -> Result<Vec<RegionRings>, BoolError> {
    // 1. Collect directed edges. A repeat in the same direction means two faces claim the same side.
    let mut dirs: HashMap<(Node, Node), usize> = HashMap::new();
    for lf in group {
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            for e in ring_edges(ring) {
                *dirs.entry(e).or_insert(0) += 1;
            }
        }
    }
    if dirs.values().any(|&c| c > 1) {
        return Err(reject(tag::COPLANAR_MERGE));
    }
    // 2. An edge carried in both directions is interior — it separates nothing. Anything carried
    //    three or more times (either direction) is non-manifold in the plane.
    let mut undirected: HashMap<(Node, Node), usize> = HashMap::new();
    for &(a, b) in dirs.keys() {
        *undirected.entry(norm_edge(a, b)).or_insert(0) += 1;
    }
    if undirected.values().any(|&c| c > 2) {
        return Err(reject(tag::COPLANAR_MERGE));
    }
    // 3. Re-thread what survives. Two outgoing edges at one node means the pieces meet at a point
    //    and the cycles are not determined.
    let mut next: HashMap<Node, Node> = HashMap::new();
    for &(a, b) in dirs.keys() {
        if dirs.contains_key(&(b, a)) {
            continue; // interior
        }
        if next.insert(a, b).is_some() {
            return Err(reject(tag::COPLANAR_MERGE));
        }
    }
    let mut starts: Vec<Node> = next.keys().copied().collect();
    starts.sort_by_key(|&nd| node_rank(nd));
    let mut seen: HashSet<Node> = HashSet::new();
    let mut cycles: Vec<Vec<Node>> = Vec::new();
    for start in starts {
        if seen.contains(&start) {
            continue;
        }
        let mut cyc = vec![start];
        seen.insert(start);
        let mut cur = next[&start];
        while cur != start {
            if !seen.insert(cur) {
                return Err(reject(tag::COPLANAR_MERGE)); // walk re-entered another cycle
            }
            cyc.push(cur);
            cur = *next.get(&cur).ok_or_else(|| reject(tag::COPLANAR_MERGE))?;
        }
        if cyc.len() < 3 {
            return Err(reject(tag::COPLANAR_MERGE));
        }
        cycles.push(cyc);
    }
    if cycles.is_empty() {
        return Err(reject(tag::COPLANAR_MERGE)); // everything erased: not a region
    }
    // 4. Winding tells an outer ring from a hole; the class root is the frame both are read in.
    let wc = canon[group[0].plane_idx];
    let mut outers: Vec<Vec<Node>> = Vec::new();
    let mut holes: Vec<Vec<Node>> = Vec::new();
    for cyc in cycles {
        match arrange::loop_winding(planes, wc, &seam_ring(&cyc))? {
            1 => outers.push(cyc),
            -1 => holes.push(cyc),
            _ => return Err(reject(tag::COPLANAR_MERGE)),
        }
    }
    // 5. Each hole belongs to the outer ring that contains it — the same question `nest_cells` asks
    //    of the arrangement's cells, answered by the same predicate.
    let mut faces: Vec<RegionRings> = outers.into_iter().map(|o| (o, Vec::new())).collect();
    for hole in holes {
        let probe = seam_ring(&hole)[0];
        let mut owner = None;
        for (i, (outer, _)) in faces.iter().enumerate() {
            if arrange::point_in_ring(planes, wc, probe, &seam_ring(outer))? {
                if owner.is_some() {
                    return Err(reject(tag::COPLANAR_MERGE)); // nested deeper than this brick names
                }
                owner = Some(i);
            }
        }
        faces[owner.ok_or_else(|| reject(tag::COPLANAR_MERGE))?]
            .1
            .push(hole);
    }
    Ok(faces)
}

/// Drop straight-angle vertices: a node whose only two neighbours across **all** rings lie with it
/// on one line. Every node is a `Seam` triple `{p, q, r}`, so the test is combinatorial and exact —
/// the node is on a line iff some pair of its planes is shared by both neighbours. Dropping from
/// every incident ring at once keeps a vertex that is a real corner somewhere (degree > 2), which is
/// what stops a T-junction from opening.
fn dissolve_straight_angles(out: &mut [LocalFace]) {
    let mut nbrs: HashMap<Node, HashSet<Node>> = HashMap::new();
    for lf in out.iter() {
        for ring in std::iter::once(&lf.loop_nodes).chain(lf.inner.iter()) {
            for (a, b) in ring_edges(ring) {
                nbrs.entry(a).or_default().insert(b);
                nbrs.entry(b).or_default().insert(a);
            }
        }
    }
    let mut drop: HashSet<Node> = HashSet::new();
    for (&node, ns) in &nbrs {
        let Node::Seam(t) = node else { continue };
        if ns.len() != 2 {
            continue;
        }
        let mut it = ns.iter();
        let (Node::Seam(a), Node::Seam(b)) = (*it.next().unwrap(), *it.next().unwrap()) else {
            continue;
        };
        let both_have = |p: usize| a.contains(&p) && b.contains(&p);
        if (both_have(t[0]) && both_have(t[1]))
            || (both_have(t[0]) && both_have(t[2]))
            || (both_have(t[1]) && both_have(t[2]))
        {
            drop.insert(node);
        }
    }
    if drop.is_empty() {
        return;
    }
    for lf in out.iter_mut() {
        lf.loop_nodes.retain(|nd| !drop.contains(nd));
        for ring in &mut lf.inner {
            ring.retain(|nd| !drop.contains(nd));
        }
    }
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

    fn p2(x: f64, y: f64) -> Point2 {
        Point2::from_array([x, y])
    }

    /// Is there an outer-shell face on the plane through `pt` with normal `n`, oriented that way?
    /// The production path names a cap by the *face* that made it (`find_face_coplanar_with`); a
    /// test that wants to say "a face sits on z = 1.5 facing +z" has no such face in hand, and
    /// asserting geometry from coordinates is exactly what a test may do.
    fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
        let Some(target) = Plane::from_point_normal(pt, n) else {
            return false;
        };
        let shell = m.solids.get(solid).outer;
        m.shells.get(shell).faces.iter().any(|&fh| {
            let f = m.faces.get(fh);
            let Surface::Plane(plane) = m.surfaces.get(f.surface) else {
                return false;
            };
            let sign = match f.orientation {
                Orientation::Forward => 1.0,
                Orientation::Reversed => -1.0,
            };
            planes_coplanar(plane, &target) && (plane.normal() * sign).dot(n) > 0.0
        })
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
        // The silent half of the same fault. A loop's orientation multiplies `orient_sign`
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
        // Cut(box − L): the box is wholly inside L ⇒ nothing remains. Nothing is an answer, so the
        // boolean succeeds with no solids — and consumes both operands like any other success.
        let (mut m, l, bx) = l_and_inner_box();
        assert!(boolean(&mut m, BoolKind::Cut, bx, l).unwrap().is_empty());
        assert!(m.live_solids.is_empty(), "both operands are consumed");
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
        // Fusing things that never touch does not merge them — it keeps both, whole.
        let (mut m, l) = l_prism();
        let vol_d = 1.0; // far()..far_max() is the unit box
        let d = m.add_cuboid(far(), far_max());
        let both = boolean(&mut m, BoolKind::Fuse, l, d).unwrap();
        assert_eq!(both.len(), 2, "disjoint operands stay two solids");
        m.rebuild_adjacency();
        let mut vols: Vec<f64> = both
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!(
            (vols[0] - vol_d).abs() < 1e-9 && (vols[1] - vol_l).abs() < 1e-9,
            "{vols:?}"
        );
        // Their intersection, on the other hand, really is empty.
        let (mut m, l) = l_prism();
        let d = m.add_cuboid(far(), far_max());
        assert!(boolean(&mut m, BoolKind::Common, l, d).unwrap().is_empty());
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

    /// The first `Common` whose kept region closes into a loop. A bar drilled through a cube:
    /// `∩ = [1,2]²×[0,3]`,
    /// the middle segment of the bar. On the cube's `z=0` and `z=3` caps the kept square
    /// `[1,2]²` is bounded entirely by the cut (an island face, no `∂f`), so the loop is
    /// oriented with material *inside* it — the sign a hole (Cut) never exercised. `1·1·3`.
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

    fn near(a: Point3, b: [f64; 3]) -> bool {
        (a - Point3::from_array(b)).norm() < 1e-9
    }

    proptest! {}

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

    /// A ring's nodes as coordinates — each triple is three planes, so its point is their meet.
    fn ring_points(planes: &[PlaneInfo], ring: &[[usize; 3]]) -> Vec<[f64; 3]> {
        ring.iter()
            .map(|t| {
                three_planes(
                    &planes[t[0]].plane,
                    &planes[t[1]].plane,
                    &planes[t[2]].plane,
                )
                .unwrap()
                .as_array()
            })
            .collect()
    }

    #[test]
    fn a_hole_winds_clockwise_and_an_island_counter_clockwise() {
        // The two rings of a holed face are stored with opposite windings — that is what makes one
        // a hole and the other its outer boundary — and `loop_winding` must read exactly that.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        assert_eq!(arrange::loop_winding(&planes, p, &outer).unwrap(), 1);
        assert_eq!(arrange::loop_winding(&planes, p, &hole).unwrap(), -1);

        // Nothing but the ring's direction went into that. Reversing it by hand agrees.
        for (name, ring, want) in [("outer", &outer, -1i8), ("hole", &hole, 1)] {
            let mut reversed = ring.clone();
            reversed.reverse();
            assert_eq!(
                arrange::loop_winding(&planes, p, &reversed).unwrap(),
                want,
                "{name} reversed"
            );
        }
    }

    #[test]
    fn a_reflex_node_turns_against_its_ring() {
        // The L cap's outer ring is a hexagon with exactly one reflex corner, at `(1, 1)`. A
        // convex node turns with the ring and the reflex one turns against it, so the turn signs
        // are *not* all equal — which is why the winding cannot be read off an arbitrary node.
        let (planes, p, outer, _) = holed_face_rings("dimple");
        let pts = ring_points(&planes, &outer);
        let reflex = pts
            .iter()
            .position(|q| near(Point3::from_array(*q), [1.0, 1.0, 1.0]))
            .expect("the reflex node");
        assert_eq!(arrange::turn_at(&planes, p, &outer, reflex).unwrap(), -1);
        let turns: Vec<i8> = (0..outer.len())
            .map(|i| arrange::turn_at(&planes, p, &outer, i).unwrap())
            .collect();
        assert_eq!(
            turns.iter().filter(|&&t| t == -1).count(),
            1,
            "one reflex corner: {turns:?}"
        );

        // The node `loop_winding` lands on is the lexicographically least — a hull vertex, where
        // the turn *is* the winding. The test finds it by reading coordinates; `loop_winding`
        // finds it with an exact predicate.
        let lo = (0..pts.len())
            .min_by(|&i, &j| pts[i].partial_cmp(&pts[j]).unwrap())
            .unwrap();
        assert_ne!(lo, reflex);
        assert_eq!(arrange::turn_at(&planes, p, &outer, lo).unwrap(), 1);

        // ★ The teeth. A ring is a cycle, so its winding cannot depend on where the walk began.
        // Start it at the reflex node and a `turn_at(ring[0])` implementation reads the reflex
        // sign — the exact fault `outer_tri` shipped. Measured: without this rotation, such an
        // implementation passes every assertion above.
        let mut rotated = outer.clone();
        rotated.rotate_left(reflex);
        assert_eq!(arrange::turn_at(&planes, p, &rotated, 0).unwrap(), -1);
        assert_eq!(arrange::loop_winding(&planes, p, &rotated).unwrap(), 1);
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

    /// One face ring, as plane triples.
    type Ring = Vec<[usize; 3]>;

    /// A **holed reflex face** from the live engine, as plane triples: `(planes, p, outer, hole)`.
    ///
    /// `Cut(L-prism, stub)` leaves the L's top cap carrying a hole — `"dimple"` a square one,
    /// `"ell"` an L-shaped one. Both rings come from [`arrange::face_vertex_triples`] and
    /// [`arrange::hole_rings`], which the boolean itself uses, so the fixture exercises only code
    /// the kernel runs.
    ///
    /// This replaces a helper that built its rings from `seam_paths_on`/`orient_seam_loop` — the
    /// retired seam engine. The *properties* below are about `point_in_ring`/`every_ray`, which are
    /// live and load-bearing (`nest_cells` picks a hole's host with them, `unify_coplanar_faces`
    /// groups by them), so they had to be re-homed rather than deleted with their old fixture.
    fn holed_face_rings(which: &str) -> (Vec<PlaneInfo>, usize, Ring, Ring) {
        let (mut m, l, stub) = if which == "dimple" {
            l_and_dimple()
        } else {
            l_and_ell_stub()
        };
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).expect("the cut");
        m.rebuild_adjacency();
        let planes = collect_planes(&m, r).unwrap();
        let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
        for (i, pi) in planes.iter().enumerate() {
            surf_ix.insert(pi.face, i);
        }
        let canon = plane_classes(&planes);
        let inc = arrange::edge_planes(&m, r, &surf_ix).unwrap();
        for &fh in &m.shells.get(m.solids.get(r).outer).faces {
            let p = surf_ix[&fh];
            let holes = arrange::hole_rings(&m, fh, p, &inc, &planes, &canon).unwrap();
            if let Some(hole) = holes.into_iter().next() {
                let outer = arrange::face_vertex_triples(&m, fh, p, &inc, &planes, &canon).unwrap();
                assert_eq!(outer.len(), 6, "{which}: the L's cap is a reflex hexagon");
                return (planes, p, outer, hole);
            }
        }
        panic!("{which}: no holed face");
    }

    #[test]
    fn a_loop_is_inside_the_face_it_was_found_on() {
        // A hole ring never touches its face's outer ring: every node is *strictly* inside, so
        // `point_in_ring` must say so for all of them. The outer ring is the L's cap, a hexagon
        // with a reflex corner, so this is not a convex test.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        for t in &hole {
            assert!(arrange::point_in_ring(&planes, p, *t, &outer).unwrap());
        }
    }

    #[test]
    fn the_older_loops_are_inside_their_faces_too() {
        // Two hole shapes that must not move: a square one and an L-shaped one. Both sit
        // strictly inside the same reflex hexagon, and **every** clear ray agrees — the parity
        // cannot depend on which ray was cast, which is a second machine for free.
        for which in ["dimple", "ell"] {
            let (planes, p, outer, hole) = holed_face_rings(which);
            for t in &hole {
                let rays = arrange::every_ray(&planes, p, *t, &outer).unwrap();
                assert!(!rays.is_empty(), "{which}: no clear ray");
                assert!(rays.iter().all(|&x| x), "{which}: {rays:?}");
            }
        }
    }

    #[test]
    fn a_loop_is_placed_by_where_it_is_not_by_how_it_winds() {
        // Containment is about **where** a ring is, never which way it runs. A hole ring is stored
        // clockwise about the face normal and its outer ring counter-clockwise, and neither
        // direction may enter the answer: reversing either must change nothing.
        //
        // And containment is **not symmetric** — the classic way to get this wrong is a test that
        // only ever asks it one way round.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        let (mut rev_outer, mut rev_hole) = (outer.clone(), hole.clone());
        rev_outer.reverse();
        rev_hole.reverse();
        for t in &hole {
            assert!(arrange::point_in_ring(&planes, p, *t, &outer).unwrap());
            assert!(
                arrange::point_in_ring(&planes, p, *t, &rev_outer).unwrap(),
                "reversing the outer ring must not move the hole"
            );
        }
        for t in &outer {
            assert!(!arrange::point_in_ring(&planes, p, *t, &hole).unwrap());
            assert!(
                !arrange::point_in_ring(&planes, p, *t, &rev_hole).unwrap(),
                "nor may reversing the hole swallow the outer ring"
            );
        }
    }

    #[test]
    fn every_clear_ray_agrees() {
        // The ring is simple, so the parity cannot depend on the ray. `every_ray` returns one
        // answer per usable candidate and they must be unanimous; a disagreement means the ray
        // choice leaked into the result.
        //
        // (Its ancestor also pinned that *half* the candidates were unusable — that count came
        // from the retired staple fixture, whose loop and arc shared a plane. The holed L cap has
        // no such sharing, so only the unanimity survives the move.)
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        let rays = arrange::every_ray(&planes, p, hole[0], &outer).unwrap();
        assert!(!rays.is_empty(), "at least one candidate is clear");
        assert!(rays.iter().all(|&x| x), "and they agree: inside — {rays:?}");
    }

    #[test]
    fn a_ring_inside_a_ring_is_what_nesting_looks_like() {
        // `nested_loops` has no operand in the suite that produces it — a polyhedral torus would.
        // The detector can still be aimed at real geometry: a holed face *is* a ring inside a ring,
        // which fires exactly the condition the nesting brick asks about. It is the detector under
        // test, not the fixture.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        assert!(arrange::point_in_ring(&planes, p, hole[0], &outer).unwrap());
        assert!(!arrange::point_in_ring(&planes, p, outer[0], &hole).unwrap());
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
    fn cut_the_stub_by_the_l_leaves_an_island_face() {
        // Swap the operands of `cut_blind_dimple` and the same loop lands on a face whose
        // boundary is *all* dropped: the kept region is the loop's interior alone. That
        // face has no `∂f` at all — its outer loop *is* the seam ring, four `Discovered`
        // vertices and nothing else. The answer is the `0.4 × 0.4 × 0.5` box above `z = 1`.
        //
        // This is where `flip` first meets a discovered hole ring. It is wound CCW about the
        // L's `+z`, keeping the material (inside the stub) on
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
    fn a_coplanar_boss_over_a_void_plane_is_solved() {
        // The boss straddles the void's x=1 and y=1 planes (it spans [0.9,1.1]²), so classifying
        // the void's walls against it is not the clean whole-face case. This was pinned as an
        // honest reject with the standing instruction that a future change may "turn it into a
        // correct 26.04, never a silent answer" — and family #3 did: once one geometric plane is
        // one class whatever the two faces' sizes, the arrangement solves it. 26 (hollow) + 0.04
        // (boss); area 60 + 0.8 (boss sides) + 0.04 (its top) − 0.04 (its footprint); void intact.
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
        let r = boolean_one(&mut m, BoolKind::Fuse, hollow, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 26.04).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 60.8).abs() < 1e-9, "area {}", props.area);
        assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
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
        // A footprint that does not touch the face at all. The boolean is not what fails here — it
        // fuses the two into a base plus a detached boss, which is the right answer (see
        // `a_touchless_boss_fuses_into_two_solids`). What breaks is the *pad's* premise, so the
        // error names that, and the model the caller is left holding is the one it started with.
        let (mut m, top) = cube_with_top();
        let far = Profile2d {
            points: vec![p2(1.8, 1.8), p2(2.2, 1.8), p2(2.2, 2.2), p2(1.8, 2.2)],
        };
        let before = m.live_solids.clone();
        assert_eq!(
            apply(&mut m, &pad_op(top, far, 0.3)),
            Err(OpError::PadMissesFace)
        );
        let (mut a, mut b) = (before, m.live_solids.clone());
        a.sort_by_key(|h| h.index());
        b.sort_by_key(|h| h.index());
        assert_eq!(a, b, "a rejected pad must leave the live model untouched");
    }

    /// **Chaining onto a fused boss.** The fuse leaves the base's `z=1` face a *ring* — a face with
    /// a hole where the boss sits — and the second boolean cuts through both. Every plane class the
    /// cut opens then meets that ring along the **hole's own edge**, which is the case that used to
    /// label inconsistently: the ring's neighbouring vertices there point *into* the hole, so
    /// reading the occupied side off a flank put the material on the wrong side of `W`. The side
    /// now comes from the ring's travel ([`trace::run_body_above`]), and the run leaves as its own
    /// homogeneous segment rather than being swallowed by the straddling stretch beside it.
    ///
    /// Hand volume: `1 + 0.5·0.5·1` fused, less the cutter's `0.2·0.2` column over `z ∈ [0.5, 2]`
    /// — `1.25 − 0.06 = 1.19`. The same shape is scored against OCCT by
    /// `boss_fuse_then_cut_matches_occt`, but that oracle is `#[ignore]`d, so this is the copy that
    /// runs on every `cargo test`.
    #[test]
    fn a_boss_fused_then_cut_through() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 1.0]),
            Point3::from_array([0.75, 0.75, 2.0]),
        );
        let bossed = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, bossed).unwrap().volume - 1.25).abs() < 1e-12,
            "the fused boss itself"
        );
        let cutter = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([0.6, 0.6, 2.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, bossed, cutter).expect("the chained cut");
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, r).unwrap().volume - 1.19).abs() < 1e-12,
            "base + boss less the drilled column: {}",
            nacre_props::mass_props(&m, r).unwrap().volume
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "a chained result is still a clean model"
        );
    }

    /// The kernel's answer for a boss that misses the face, stated on its own so nobody "fixes" the
    /// boolean to reject it: fusing two solids that do not touch **is** two solids, and both are
    /// whole. Only `pad` refuses that outcome, because a pad is defined as material joined to a face
    /// (`pad_overhang_off_the_face_is_rejected`).
    #[test]
    fn a_touchless_boss_fuses_into_two_solids() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        // Shares the z = 1 plane class with the base's top, but sits far away in x/y — so the plane
        // carries two separate bodies, which is exactly what `hole_roots` used to refuse.
        let boss = m.add_cuboid(
            Point3::from_array([1.8, 1.8, 1.0]),
            Point3::from_array([2.2, 2.2, 1.3]),
        );
        let solids = boolean(&mut m, BoolKind::Fuse, base, boss).unwrap();
        assert_eq!(solids.len(), 2, "disjoint operands stay two solids");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let mut vols: Vec<f64> = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!(
            (vols[0] - 0.048).abs() < 1e-12 && (vols[1] - 1.0).abs() < 1e-12,
            "both pieces whole: {vols:?}"
        );
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
    fn a_non_convex_pad_cantilevers_and_runs_flush() {
        // **Three hard properties at once**, which is what makes this footprint worth keeping. The
        // frame maps local `[px, py]` to world `(0.5 + py, 0.5 − px)`, so the L below lands on
        // `(0.25,0.75) (0.25,−0.25) (0.75,−0.25) (0.75,0.25) (1.0,0.25) (1.0,0.75)`:
        //   1. **non-convex** — the L has a reflex corner at `(0.75, 0.25)`;
        //   2. **overhanging** — `y < 0` cantilevers past the cube's `y = 0` edge;
        //   3. **flush** — the edge `x = 1.0, y∈[0.25,0.75]` lies *exactly* on the face's `x = 1`
        //      boundary, the "profile rim shares the face rim" case.
        // It used to reject because the overhang sidecars gated on convexity; that gate is gone.
        //
        // Hand-checked shape: footprint `0.25 + 0.375 = 0.625`, prism wholly above `z = 1`, so
        //   volume 1 + 0.625 = 1.625
        //   area   5 (cube minus its top) + 0.5 (top left uncovered) + 3.5 (prism sides)
        //          + 0.625 (prism cap) + 0.125 (the cantilever's underside) = 9.75
        // The underside term is the cantilever: a contained pad would not have one.
        //
        // OCCT cannot score this directly — `pad` builds its tool prism internally, and rebuilding
        // it here would lean on the same frame mapping the assertion is testing.
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
        let OpOutput::PadOnFace { solid, top_face } =
            apply(&mut m, &pad_op(top, l_over, 1.0)).expect("the cantilevered L pad")
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, solid).unwrap();
        assert!((p.volume - 1.625).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 9.75).abs() < 1e-12, "area {}", p.area);
        assert!(m.reachable().faces.contains(&top_face)); // boss top cap recovered
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
        // Each `tri` is three NON-collinear points of its own plane. A degenerate `tri` (three
        // equal points) would make every `orient3d` vanish, so the coordinate branch would report
        // coplanar and this test would pass without the handle branch ever mattering.
        let mk = |plane, tri| PlaneInfo {
            surf: shared,
            face: fh,
            plane,
            tri,
            n_out: Vector3::from_array([0.0; 3]),
            orient: Orientation::Forward,
            tri_pt3: None,
            class: usize::MAX,
        };
        let p = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
        let planes = vec![
            mk(plane_x0, [p(0., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]), // in x = 0
            mk(plane_z0, [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.)]), // in z = 0
        ];
        // Neither fallback fires: the coefficients are not proportional, and the coordinates say
        // these really are two different planes.
        assert!(!planes_coplanar(&planes[0].plane, &planes[1].plane));
        assert!(!tolerant::t_planes_coplanar(&planes, 0, 1));
        // The shared handle alone makes them coplanar-by-reference.
        assert!(shares_or_coplanar(&planes, 0, 1));
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

    /// A prism raised on a `(1,1,1)`-slanted sketch plane, its far cap holding a blind pocket. The
    /// cap and its two anti-parallel side walls meet in a triple whose `raw` coefficients are
    /// exactly dependent (`det = 0`); before family #3's dir-sign fix the guard read `sqrt`-rounded
    /// unit normals, called that triple non-degenerate, and the consumer aborted on `D = 0`. Now
    /// the guard reads the same coefficients the consumer does, so the arrangement runs. Volume:
    /// a `2×2` base × `2` deep block is `8`, less the `0.4²×0.5` pocket.
    #[test]
    fn a_pocket_on_a_slanted_face() {
        let plane =
            SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
                .unwrap();
        let mut m = Model::new();
        let big = Profile2d {
            points: vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)],
        };
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane,
                profile: big,
                dist: 2.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let out = apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)).unwrap();
        let OpOutput::PocketOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
        assert!((vol - (8.0 - 0.16 * 0.5)).abs() < 1e-9, "volume {vol}");
    }

    /// The same slanted cap, but a boss (pad, Fuse) instead of a pocket — the sweep runs the other
    /// way, a different code path. Volume: the `8` block plus a `0.4²×0.5` stub.
    #[test]
    fn a_pad_on_a_slanted_face() {
        let plane =
            SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
                .unwrap();
        let mut m = Model::new();
        let big = Profile2d {
            points: vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)],
        };
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane,
                profile: big,
                dist: 2.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let out = apply(
            &mut m,
            &Operation::PadOnFace {
                face: faces[1],
                profile: small_square(),
                dist: 0.5,
            },
        )
        .unwrap();
        let OpOutput::PadOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
        assert!((vol - (8.0 + 0.16 * 0.5)).abs() < 1e-9, "volume {vol}");
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

        /// A blind pocket on a randomly-slanted face: the arrangement must give a valid solid of the
        /// right volume or reject honestly — **never panic**. Drives general (non-axis) plane normals
        /// through the dir-sign guard and the `angular_order`/`turn_at` consumers that read its zeros.
        ///
        /// Was `#[ignore]`d for a residual `D = 0` panic on general normals: a triple naming one
        /// geometric plane through two coincident faces. Symmetric normals like `(1,1,1)` cleared
        /// it, `(0.446, 0.737, 0.990)` did not. **Un-ignored 2026-07-22** — naming every plane by
        /// its class made those two faces one index, so the degenerate triple can no longer form.
        #[test]
        fn pocket_on_a_random_slanted_face_is_valid_or_rejects(
            nx in -1.0f64..1.0,
            ny in -1.0f64..1.0,
            nz in 0.2f64..1.0, // keep the normal clear of the sketch's degenerate zero
        ) {
            let normal = Vector3::from_array([nx, ny, nz]);
            prop_assume!(normal.norm() > 0.3);
            let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
            let mut m = Model::new();
            let big = Profile2d {
                points: vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)],
            };
            let OpOutput::Extrude { faces, .. } =
                apply(&mut m, &Operation::Extrude { plane, profile: big, dist: 2.0 }).unwrap()
            else { unreachable!() };
            match apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)) {
                Ok(OpOutput::PocketOnFace { solid, .. }) => {
                    m.rebuild_adjacency();
                    prop_assert!(nacre_validate::validate(&m).is_empty());
                    let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
                    prop_assert!((vol - (8.0 - 0.16 * 0.5)).abs() <= 1e-9 * 8.0, "volume {}", vol);
                }
                Ok(_) => prop_assert!(false, "unexpected op output"),
                Err(_) => {} // an honest reject is acceptable; a panic is not (and would fail the test)
            }
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

    /// **Every producer states its side in the label frame** — the invariant family #2 restored,
    /// swept over the whole two-solid corpus (prints the interesting classes with `--nocapture`).
    ///
    /// The arrangement states its cell labels as `[*_above, *_below]` about one direction per plane
    /// class: the class root's **stored surface normal** (`Seated{body_above}` and `emit_faces`'
    /// `flip` are written against it). `arrange::side_of` answers in the root's **outward** frame
    /// instead, and the two are opposite exactly when the root face is `Reversed`
    /// (`orient_sign == -1`) — which no `add_cuboid` face ever is, but a face an earlier boolean
    /// re-emitted flipped is. `graze_above` read `side_of` raw, so on a pocket wall it flipped the
    /// wrong label bit. It no longer reads a point's side at all — [`trace::run_body_above`] derives
    /// the occupied side from the ring's travel, and the frame term cancels there because
    /// `order_along`'s direction and the label frame are defined by the same stored normal — but
    /// this class stays the corpus's only crossed-frame witness, so it is what would catch a
    /// producer that regresses to a raw `side_of`.
    ///
    /// What this pins, measured before the fix (2026-07-22):
    /// - the pocket fixture has 5 `orient_sign == -1` classes carrying seated *and* graze segments
    ///   (4 walls + the floor); every other fixture has **no** `Reversed` root at all, which is why
    ///   the whole corpus passed with the frames crossed and why converting cannot regress it;
    /// - the four wall classes stopped at `loop_orient_mismatch`; the floor class did **not** — its
    ///   four rim edges all carry a graze, so seated and graze were wrong *together*, consistently,
    ///   and the label survived verification while being inverted (a silent wrong, not a reject);
    /// - the hand-derived sides on those classes, so a re-crossed frame fails here first.
    #[test]
    fn every_producer_states_its_side_in_the_label_frame() {
        let mut boxed: Vec<(&str, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
        macro_rules! fixture {
            ($name:ident) => {{
                let (m, a, b) = $name();
                boxed.push((stringify!($name), m, a, b));
            }};
        }
        fixture!(two_boxes);
        fixture!(nested_boxes);
        fixture!(cube_and_notch);
        fixture!(stacked_cubes);
        fixture!(l_and_corner_box);
        fixture!(l_and_reflex_box);
        fixture!(l_and_inner_box);
        fixture!(l_and_popup_box);
        fixture!(l_and_notch_bar);
        fixture!(l_and_ell_stub);
        fixture!(l_and_staple);
        fixture!(l_and_dimple);
        fixture!(l_and_rod);
        fixture!(u_and_slab);
        {
            // The pocket family: `pocketed_cube` is itself a boolean result, so its pocket walls
            // are `Reversed` faces. Box coordinates are `pocket_corner_cut`'s.
            let (mut m, pc) = pocketed_cube();
            let bx = m.add_cuboid(
                Point3::from_array([0.85, 0.85, 0.85]),
                Point3::from_array([1.15, 1.15, 1.15]),
            );
            boxed.push(("pocket_corner_cut", m, pc, bx));
        }

        let mut reversed_with_graze: Vec<String> = Vec::new();
        let mut reversed_seated_only: Vec<String> = Vec::new();
        for (name, m, a, b) in &boxed {
            let audits = trace::frame_audit(m, BoolKind::Cut, *a, *b).unwrap();
            for au in &audits {
                let interesting = au.orient_sign < 0 || au.failed_at.is_some();
                if !interesting {
                    continue;
                }
                let where_ = format!(
                    "{name}: wc={} pt={:?} n={:?} orient_sign={} seated={:?} graze={:?} trans={} \
                     declined={:?} failed_at={:?}",
                    au.wc,
                    au.root_point,
                    au.root_normal,
                    au.orient_sign,
                    au.seated,
                    au.grazes,
                    au.transversals,
                    au.declined,
                    au.failed_at,
                );
                println!("{where_}");
                if au.orient_sign < 0 {
                    if au.grazes.is_empty() {
                        if !au.seated.is_empty() {
                            reversed_seated_only.push(where_);
                        }
                    } else {
                        reversed_with_graze.push(where_);
                    }
                }
            }
        }
        println!("--- reversed-root classes carrying a graze (the set the fix moves) ---");
        for r in &reversed_with_graze {
            println!("  {r}");
        }
        println!("--- reversed-root classes with seated but no graze (the alternative's risk) ---");
        for r in &reversed_seated_only {
            println!("  {r}");
        }
        // The pocket fixture must keep supplying such classes, or this test has stopped exercising
        // the crossed-frame configuration and would pass vacuously.
        assert_eq!(
            reversed_with_graze.len(),
            5,
            "the pocket's 4 walls + floor are the corpus's only reversed-root classes with a graze"
        );
        assert!(
            reversed_with_graze.iter().all(|r| r.starts_with("pocket")),
            "no other fixture may have one: {reversed_with_graze:?}"
        );

        // Hand-derived sides on the pocket wall class `x = 0.7` (root = the pocket's +x wall, whose
        // outward normal points into the void, so the stored normal `+x` makes "above" the material
        // side `x > 0.7`). The wall is seated with its body above; the two side walls and the floor
        // graze it from `x < 0.7`, i.e. below. Crossed frames invert the grazes.
        //
        // The **box's top face** grazes it too, from `x > 0.7`: the pocket's opening makes that face
        // a notched region whose edge rides this plane with the material outside the pocket. That is
        // a run whose flanks *differ*, which the engine used to read as a straddling transversal —
        // this class is the corpus's only `Reversed` root, so it is also the only place the frame
        // handling of `trace::run_body_above` is exercised against a crossed frame: the three
        // `false` entries below are the pre-existing answers, unchanged by the new rule.
        let (m, pc, bx) = boxed
            .iter()
            .find_map(|(n, m, a, b)| (*n == "pocket_corner_cut").then_some((m, *a, *b)))
            .unwrap();
        let wall = trace::frame_audit(m, BoolKind::Cut, pc, bx)
            .unwrap()
            .into_iter()
            .find(|au| au.root_point == [0.7, 0.7, 1.0] && au.root_normal == [1.0, 0.0, 0.0])
            .expect("the x=0.7 pocket wall class");
        assert_eq!(wall.orient_sign, -1, "the pocket wall is a Reversed face");
        assert_eq!(
            wall.seated,
            vec![true; 4],
            "body above = material at x > 0.7"
        );
        assert_eq!(
            wall.grazes,
            vec![true, false, false, false],
            "the top face grazes from x > 0.7 (its pocket-opening edge, material outside); \
             the two side walls and the floor graze from x < 0.7"
        );
        assert_eq!(wall.failed_at, None, "the class labels consistently");
    }

    /// `pocket_corner_cut` by hand, so the pocket family keeps a regression net that runs without
    /// OCCT: the unit cube less a `0.4²×0.5` pocket is `0.92`, and the corner box `[0.85,1.15]³`
    /// bites `0.15³` of solid (it clears the pocket, whose footprint stops at `x = 0.7`).
    #[test]
    fn a_corner_cut_off_a_pocketed_cube() {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.85, 0.85, 0.85]),
            Point3::from_array([1.15, 1.15, 1.15]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (0.92 - 0.15 * 0.15 * 0.15)).abs() < 1e-12,
            "volume {}",
            props.volume
        );
        // A corner bite replaces three 0.15² squares with three more: the area is unchanged at
        // 6 − 0.16 (the lid's hole) + 0.8 (four pocket walls) + 0.16 (its floor).
        assert!((props.area - 6.8).abs() < 1e-12, "area {}", props.area);
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
    /// Renamed from `..._is_empty`: that name recorded the old engine's limit, not the answer.
    /// A union of things that never touch is both of them, whole.
    fn fuse_of_disjoint_boxes_is_two_solids() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        for &s in &solids {
            let v = nacre_props::mass_props(&m, s).unwrap().volume;
            assert!((v - 1.0).abs() < 1e-12, "each unit box survives whole: {v}");
        }
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
        // Cut(inner − outer): inner is wholly removed ⇒ empty, which is an answer, not an error —
        // and a successful boolean consumes its operands, so the Fuse below needs a fresh model
        // (it used to reuse this one only because the empty Cut was an error that consumed nothing).
        let (mut m, outer, inner) = nested_boxes();
        assert!(
            boolean(&mut m, BoolKind::Cut, inner, outer)
                .unwrap()
                .is_empty()
        );
        assert!(m.live_solids.is_empty(), "both operands are consumed");

        // Fuse(inner ∪ outer) = outer.
        let (mut m, outer, inner) = nested_boxes();
        let vol_outer = nacre_props::mass_props(&m, outer).unwrap().volume;
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
        assert!(has_face_on_plane(
            &m,
            r,
            Point3::from_array([0.5, 0.5, 1.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        ));
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
        assert!(has_face_on_plane(
            &m,
            r,
            Point3::from_array([0.4, 0.4, 1.4]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        ));
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
        // They meet only along the base's top face — a contact of zero volume.
        assert!(
            boolean(&mut m, BoolKind::Common, base, corner)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn cut_by_an_overhanging_boss_carrying_a_pin_owes_a_notch() {
        // Same seating, but the tool carries a pin reaching below the contact plane, so the plane
        // no longer separates the solids and the cut owes a real notch (1 − 0.2·0.2·0.5 = 0.98).
        //
        // This used to be an honest reject: the tool's z=1 cap is an annulus-like face whose
        // *inner* edge rides the pin's walls, and reading its occupancy off the ring's flank put
        // the material on the wrong side, so the class would not label. With the side read from
        // the ring's travel instead (`trace::run_body_above`), the notch comes out at the
        // hand-computed volume with a clean model.
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
        let notched =
            boolean_one(&mut m, BoolKind::Cut, base, tool).expect("the notch is buildable");
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, notched).unwrap().volume - 0.98).abs() < 1e-12,
            "the notch the tool owes: {}",
            nacre_props::mass_props(&m, notched).unwrap().volume
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "a notched result is still a clean model"
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

    /// An axis-aligned `PlaneInfo` at `d` along its normal, with a **non-degenerate `tri`** whose
    /// right-hand normal is `n_out`. The merge reads more than `n_out` now — `loop_winding` and
    /// `point_in_ring` name their arguments by plane and evaluate exact predicates on `tri` — so a
    /// dummy triangle would make those answers meaningless.
    fn mk_axis_plane(m: &mut Model, axis: usize, d: f64, positive: bool) -> PlaneInfo {
        let mut n = [0.0; 3];
        n[axis] = if positive { 1.0 } else { -1.0 };
        let normal = Vector3::from_array(n);
        let mut at = [0.0; 3];
        at[axis] = d;
        let origin = Point3::from_array(at);
        let plane = Plane::from_point_normal(origin, normal).unwrap();
        let surf = m.surfaces.push(Surface::Plane(plane));
        let face = m.faces.push(Face {
            surface: surf,
            outer: Loop { half_edges: vec![] },
            inner: vec![],
            orientation: Orientation::Forward,
        });
        // Two in-plane directions whose cross product is `+normal`, so `tri` winds outward.
        let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
        let (i, j) = if positive { (i, j) } else { (j, i) };
        let step = |k: usize| {
            let mut q = at;
            q[k] += 1.0;
            Point3::from_array(q)
        };
        PlaneInfo {
            surf,
            face,
            plane,
            tri: [origin, step(i), step(j)],
            n_out: normal,
            orient: Orientation::Forward,
            tri_pt3: None,
            class: usize::MAX,
        }
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
            class: usize::MAX,
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
        //
        // Named the way the arrangement names things — every vertex is the meeting of three
        // planes — because the straight-angle test reads those triples. (The old `Orig` fixture
        // exercised a path the engine stopped producing when it went all-`Seam`.)
        let mut m = Model::new();
        let p = vec![
            mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0, the shared class
            mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
            mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
            mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
            mk_axis_plane(&mut m, 0, 3.0, true),  // 4: x=3
            mk_axis_plane(&mut m, 1, 0.0, false), // 5: y=0
            mk_axis_plane(&mut m, 1, 1.0, true),  // 6: y=1
        ];
        let canon: Vec<usize> = (0..p.len()).collect();
        let v = |x: usize, y: usize| Node::Seam([0, x, y]); // sorted: class, x-plane, y-plane
        let (c00, c10, c20, c30) = (v(1, 5), v(2, 5), v(3, 5), v(4, 5));
        let (c01, c11, c21, c31) = (v(1, 6), v(2, 6), v(3, 6), v(4, 6));
        let faces = vec![
            face(0, vec![c00, c10, c11, c01], vec![]),
            face(0, vec![c10, c20, c21, c11], vec![]),
            face(0, vec![c20, c30, c31, c21], vec![]),
        ];
        let out = unify_coplanar_faces(faces, &p, &canon).unwrap();
        assert_eq!(out.len(), 1, "three coplanar faces fuse into one");
        let l = &out[0].loop_nodes;
        assert_eq!(l.len(), 4, "straight-angle mid vertices dissolved: {l:?}");
        for c in [c00, c30, c31, c01] {
            assert!(l.contains(&c), "corner kept");
        }
        for c in [c10, c20, c11, c21] {
            assert!(!l.contains(&c), "mid vertex dropped");
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
        let out = unify_coplanar_faces(faces, &p, &canon).unwrap();
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
        let out = unify_coplanar_faces(faces, &p, &canon).unwrap();
        assert_eq!(out.len(), 2, "a holed face is not merged");
    }

    #[test]
    fn a_hole_filled_by_two_faces_still_merges() {
        // Replaces `unify_skips_seam_shared_edges`, whose premise is gone twice over: the `Orig`
        // gate it pinned was removed when the engine went all-`Seam`, and `splice_along`, whose
        // panic it guarded against, no longer exists.
        //
        // What matters now is that the merge is not special-cased to "a hole filled by exactly one
        // neighbour". A [0,3]² face with a [1,2]² hole, and that hole filled by **two** pieces split
        // at x=1.5: every ring edge between them is carried in both directions, so erasing interior
        // boundary leaves only the outer square — one face, no hole, whatever the filling is cut
        // into. This is the case that separates a general rule from a bespoke one.
        let mut m = Model::new();
        let p = vec![
            mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
            mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
            mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
            mk_axis_plane(&mut m, 0, 1.5, true),  // 3: x=1.5, where the filling is split
            mk_axis_plane(&mut m, 0, 2.0, true),  // 4: x=2
            mk_axis_plane(&mut m, 0, 3.0, true),  // 5: x=3
            mk_axis_plane(&mut m, 1, 0.0, false), // 6: y=0
            mk_axis_plane(&mut m, 1, 1.0, true),  // 7: y=1
            mk_axis_plane(&mut m, 1, 2.0, true),  // 8: y=2
            mk_axis_plane(&mut m, 1, 3.0, true),  // 9: y=3
        ];
        let canon: Vec<usize> = (0..p.len()).collect();
        let v = |x: usize, y: usize| Node::Seam([0, x, y]);
        let (o00, o30, o33, o03) = (v(1, 6), v(5, 6), v(5, 9), v(1, 9));
        let (h11, h12, h22, h21) = (v(2, 7), v(2, 8), v(4, 8), v(4, 7));
        let (m12, m11) = (v(3, 8), v(3, 7)); // the split points on the hole's top and bottom
        let faces = vec![
            // Outer square with the hole, wound the way `emit_faces` states it: outer CCW, hole CW.
            face(
                0,
                vec![o00, o30, o33, o03],
                vec![vec![h11, h12, m12, h22, h21, m11]],
            ),
            face(0, vec![h11, m11, m12, h12], vec![]), // left filler
            face(0, vec![m11, h21, h22, m12], vec![]), // right filler
        ];
        let out = unify_coplanar_faces(faces, &p, &canon).unwrap();
        assert_eq!(out.len(), 1, "the hole is filled, so one face remains");
        assert!(out[0].inner.is_empty(), "and it has no hole left");
        assert_eq!(out[0].loop_nodes.len(), 4, "just the outer square");
        for c in [o00, o30, o33, o03] {
            assert!(out[0].loop_nodes.contains(&c), "outer corner kept");
        }
    }

    #[test]
    fn unify_keeps_a_vertex_that_is_a_corner_elsewhere() {
        // F0,F1 on z=0 merge; their shared-edge endpoint (1,0,0) is a straight angle on the
        // merged face but a real corner on a perpendicular face G (plane y=0). Global degree 3
        // ⇒ it is NOT dissolved — a per-face local rule would have, opening a T-junction. Its
        // twin (1,1,0), on the merged face only, IS dissolved.
        let mut m = Model::new();
        let p = vec![
            mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
            mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
            mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
            mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
            mk_axis_plane(&mut m, 1, 0.0, false), // 4: y=0, the perpendicular face's plane
            mk_axis_plane(&mut m, 1, 1.0, true),  // 5: y=1
            mk_axis_plane(&mut m, 2, 1.0, true),  // 6: z=1
        ];
        let canon: Vec<usize> = (0..p.len()).collect();
        let (v000, v100, v200) = (
            Node::Seam([0, 1, 4]),
            Node::Seam([0, 2, 4]),
            Node::Seam([0, 3, 4]),
        );
        let (v010, v110, v210) = (
            Node::Seam([0, 1, 5]),
            Node::Seam([0, 2, 5]),
            Node::Seam([0, 3, 5]),
        );
        let (v101, v201) = (Node::Seam([2, 4, 6]), Node::Seam([3, 4, 6]));
        let faces = vec![
            face(0, vec![v000, v100, v110, v010], vec![]),
            face(0, vec![v100, v200, v210, v110], vec![]),
            face(4, vec![v200, v100, v101, v201], vec![]), // perpendicular, not coplanar
        ];
        let out = unify_coplanar_faces(faces, &p, &canon).unwrap();
        assert_eq!(out.len(), 2, "z=0 pair merges; G stays");
        let merged = out.iter().find(|lf| lf.plane_idx == 0).unwrap();
        assert!(
            merged.loop_nodes.contains(&v100),
            "corner-elsewhere vertex kept (no T-junction)"
        );
        assert!(
            !merged.loop_nodes.contains(&v110),
            "pure straight-angle vertex dropped"
        );
    }

    #[test]
    fn common_stacked_cubes_is_empty() {
        // The stack shares only its interface plane, so the intersection has no volume.
        let (mut m, a, b) = stacked_cubes();
        assert!(boolean(&mut m, BoolKind::Common, a, b).unwrap().is_empty());
        assert!(m.live_solids.is_empty(), "both operands are consumed");
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
            // Not `prop_assume!`: a matched-footprint stack is squarely in coverage whatever the
            // dimensions are, so a reject here is a defect, not an uninteresting sample. Assuming
            // it away is how this property went on passing while the kernel aborted on 2% of the
            // space and rejected 95% of it (family #3).
            let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1)
                .expect("stacked boxes fuse at any dimensions");
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            // The stack shares only its interface plane ⇒ no volume in common, at any dimensions.
            prop_assert!(boolean(&mut m2, BoolKind::Common, a2, b2).unwrap().is_empty());

            let (mut m3, a3, b3) = build();
            let rc = boolean_one(&mut m3, BoolKind::Cut, a3, b3).unwrap();
            let vc = nacre_props::mass_props(&m3, rc).unwrap().volume;
            prop_assert!((vc - va).abs() <= 1e-9 * va, "cut {vc}");
        }
    }

    /// **The topology of a boolean does not depend on whether the coordinates are
    /// f64-representable.** The same shape is built twice — once on tidy integers, once on
    /// dimensions that are not exact binary fractions — and both must give the same b-rep counts,
    /// with each volume matching its own formula.
    ///
    /// This is the invariant family #3 restored. Plane identity used to be read from the faces'
    /// *derived* coefficients, which are not exactly proportional for two differently-sized faces
    /// on one plane, so one plane became two classes and the arrangement named one point twice —
    /// but only when the arithmetic did not happen to cancel, which tidy coordinates hid
    /// (measured before the fix: 200/200 random stacked pairs under-merged, 12/600 ops aborted).
    #[test]
    fn boolean_topology_is_the_same_on_untidy_coordinates() {
        let counts = |dx: f64, dy: f64, z0: f64, h1: f64, h2: f64, kind: BoolKind| {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0, 0.0, z0]),
                Point3::from_array([dx, dy, z0 + h1]),
            );
            let b = m.add_cuboid(
                Point3::from_array([0.0, 0.0, z0 + h1]),
                Point3::from_array([dx, dy, z0 + h1 + h2]),
            );
            let r = boolean_one(&mut m, kind, a, b).expect("stacked boxes fuse/cut");
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty());
            let s = m.solids.get(r);
            let sh = m.shells.get(s.outer);
            let faces = sh.faces.len();
            let mut edges = std::collections::HashSet::new();
            let mut verts = std::collections::HashSet::new();
            for &fh in &sh.faces {
                let f = m.faces.get(fh);
                for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for he in &l.half_edges {
                        edges.insert(he.edge);
                        if let Some(bd) = m.edges.get(he.edge).bounds {
                            verts.extend(bd);
                        }
                    }
                }
            }
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            ((faces, edges.len(), verts.len(), s.cavities.len()), vol)
        };
        // The untidy dimensions are the minimal case the proptest shrank to when the kernel aborted.
        let (tidy_dx, tidy_dy, tidy_z0, tidy_h1, tidy_h2) = (2.0, 0.5, 0.0, 0.5, 2.0);
        let (dx, dy, z0, h1, h2) = (
            1.628165457453874,
            0.5,
            0.11200046228159026,
            0.5,
            2.07926124157585,
        );
        for kind in [BoolKind::Fuse, BoolKind::Cut] {
            let (tidy_shape, tidy_vol) = counts(tidy_dx, tidy_dy, tidy_z0, tidy_h1, tidy_h2, kind);
            let (shape, vol) = counts(dx, dy, z0, h1, h2, kind);
            assert_eq!(tidy_shape, shape, "{kind:?}: same topology either way");
            let want = |a: f64, b: f64| match kind {
                BoolKind::Fuse => a + b,
                _ => a,
            };
            let (tw, w) = (
                want(tidy_dx * tidy_dy * tidy_h1, tidy_dx * tidy_dy * tidy_h2),
                want(dx * dy * h1, dx * dy * h2),
            );
            assert!((tidy_vol - tw).abs() < 1e-9, "{kind:?} tidy {tidy_vol}");
            assert!((vol - w).abs() < 1e-9 * w, "{kind:?} untidy {vol}");
        }
    }

    /// One geometric plane is one class **whatever the two faces' sizes**, and the coefficient test
    /// alone still cannot say so — the two walls' un-normalized coefficient 4-vectors are not
    /// exactly proportional. Pins that the coordinate branch is what earns the merge.
    #[test]
    fn one_plane_is_one_class_whatever_the_face_size() {
        let (dx, dy) = (1.628165457453874f64, 0.5f64);
        let (z0, h1, h2) = (0.11200046228159026f64, 0.5f64, 2.07926124157585f64);
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, z0]),
            Point3::from_array([dx, dy, z0 + h1]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, z0 + h1]),
            Point3::from_array([dx, dy, z0 + h1 + h2]),
        );
        let (planes, _, _, _, canon) = plane_index_setup(&m, a, b).unwrap();
        // The two `+X` walls: same plane x = dx, different face sizes (heights h1 vs h2).
        let x_walls: Vec<usize> = (0..planes.len())
            .filter(|&i| {
                planes[i].n_out.as_array() == [1.0, 0.0, 0.0]
                    && (planes[i].tri[0].as_array()[0] - dx).abs() < 1e-12
            })
            .collect();
        assert_eq!(x_walls.len(), 2, "one wall from each box: {x_walls:?}");
        let (i, j) = (x_walls[0], x_walls[1]);
        assert!(
            !planes_coplanar(&planes[i].plane, &planes[j].plane),
            "the coefficient test still cannot prove these coplanar — that is the whole point"
        );
        assert!(
            tolerant::t_planes_coplanar(&planes, i, j),
            "coordinates can"
        );
        assert_eq!(canon[i], canon[j], "so they are one class");
    }

    /// A vertex where one plane is split between two faces is named by **the planes that touch it**,
    /// not by the loop's two neighbours.
    ///
    /// The overhang chain puts the base's exposed top and the cantilever's underside on one plane
    /// (`z = 1`, opposite normals — `unify` rightly keeps them apart, `canon` rightly calls them one
    /// class). At `(1, 0.25, 1)` the side wall's loop runs straight through their shared line, so the
    /// old rule named that vertex with the same plane twice: a triple defining no point, which the
    /// exact predicates — whose precondition is `D ≠ 0` — aborted on. Measured 2026-07-22 as the only
    /// path a degenerate triple reached them (16 arrivals in the OCCT suite, now 0).
    #[test]
    fn a_vertex_is_named_by_the_planes_that_touch_it() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let cutter = m.add_cuboid(
            Point3::from_array([1.1, 0.35, 0.5]),
            Point3::from_array([1.4, 0.65, 2.5]),
        );
        let (planes, surf_ix, inc_a, _, canon) = plane_index_setup(&m, overhung, cutter).unwrap();
        let mut checked = 0usize;
        for sh in solid_shell_handles(&m, overhung) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let tris =
                    arrange::face_vertex_triples(&m, fh, p, &inc_a, &planes, &canon).unwrap();
                for t in &tris {
                    let c: Vec<usize> = t.iter().map(|&k| canon[k]).collect();
                    assert!(
                        c[0] != c[1] && c[1] != c[2] && c[0] != c[2],
                        "vertex triple {t:?} names one class twice ({c:?})"
                    );
                    // A name that denotes three distinct classes must denote a real point.
                    assert!(
                        three_planes(
                            &planes[t[0]].plane,
                            &planes[t[1]].plane,
                            &planes[t[2]].plane
                        )
                        .is_some(),
                        "triple {t:?} defines no point"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "the chained operand has vertices to name");
    }

    /// **A plane triple is always in class form.** `planes` is a per-face table, so the same
    /// `usize` could mean "face" or "plane"; producers settle it by emitting class roots, and a
    /// consumer's raw `==` then means "same plane". Four silent-wrong bugs on this branch came from
    /// the two meanings meeting in one comparison, so the invariant is asserted, not assumed.
    #[test]
    fn plane_triples_are_always_canon() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        // An *overhanging* boss splits `z = 1` between two faces with opposite normals (the base's
        // exposed top and the boss underside) — the shape that makes "face index" and "plane index"
        // differ at all. A boss sitting wholly inside the top merges into one holed face instead.
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let probe = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([0.6, 0.6, 2.5]),
        );
        let (planes, surf_ix, inc_a, _, canon) = plane_index_setup(&m, chained, probe).unwrap();
        // The fixture must actually merge two faces into one class, or this proves nothing.
        assert!(
            canon.iter().enumerate().any(|(i, &c)| c != i),
            "fixture has no split plane — the invariant would be vacuous"
        );
        let mut checked = 0usize;
        for sh in solid_shell_handles(&m, chained) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let mut rings =
                    vec![arrange::face_vertex_triples(&m, fh, p, &inc_a, &planes, &canon).unwrap()];
                rings.extend(arrange::hole_rings(&m, fh, p, &inc_a, &planes, &canon).unwrap());
                for t in rings.iter().flatten() {
                    for &k in t {
                        assert_eq!(canon[k], k, "triple {t:?} names face {k}, not its class");
                    }
                    assert!(
                        t[0] < t[1] && t[1] < t[2],
                        "triple {t:?} is not three distinct classes in sorted order"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "the chained operand has vertices to name");
    }

    /// **A point on the cut plane reads zero no matter which face names it.** `t_orient3d`'s
    /// on-plane shortcut is a raw `==` against the triple, so before the triples were canon a
    /// vertex named by face 6 of the `z = 1` class was invisible to a query about face 1 of that
    /// same class, and the numeric branch answered ±1 for a point lying exactly on the plane.
    #[test]
    fn a_vertex_on_the_cut_plane_reads_zero_whichever_face_names_it() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let probe = m.add_cuboid(
            Point3::from_array([1.1, 0.35, 0.5]),
            Point3::from_array([1.4, 0.65, 2.5]),
        );
        let (planes, surf_ix, inc_a, _, canon) = plane_index_setup(&m, chained, probe).unwrap();
        assert!(
            canon.iter().enumerate().any(|(i, &c)| c != i),
            "fixture has no split plane — a sibling face is what this test is about"
        );
        let mut on_plane = 0usize;
        for sh in solid_shell_handles(&m, chained) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let tris =
                    arrange::face_vertex_triples(&m, fh, p, &inc_a, &planes, &canon).unwrap();
                for t in &tris {
                    // Ask about every face of every class the vertex names — including the sibling
                    // faces that are not the class root, which is where the old bug lived.
                    for (q, _) in planes.iter().enumerate() {
                        if !t.contains(&canon[q]) {
                            continue;
                        }
                        assert_eq!(
                            arrange::side_of(&planes, *t, canon[q]),
                            0,
                            "vertex {t:?} lies on plane class {} (face {q}) but does not read 0",
                            canon[q]
                        );
                        on_plane += 1;
                    }
                }
            }
        }
        assert!(on_plane > 0, "some vertex lies on some queried plane");
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
        // A failing boolean's error is surfaced as `OpError::Boolean`. This used to be driven by a
        // disjoint `Common`, but that is no longer an error (it is an empty result, see
        // `boolean_op_passes_an_empty_result_through`), so the wrapping is exercised with a boolean
        // that genuinely fails: a handle that is not live.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
        m.live_solids.retain(|&s| s != b); // retire `b` behind the op's back
        assert_eq!(
            apply(
                &mut m,
                &Operation::Boolean {
                    kind: BoolKind::Common,
                    a,
                    b
                }
            ),
            Err(OpError::Boolean(BoolError::InputNotLive))
        );
    }

    /// **An empty result is an answer.** Two solids that miss each other have no intersection, and
    /// that is what `Common` reports: `Ok` with no solids, both operands consumed like any other
    /// successful boolean. Stated on its own because the name is the contract — if someone makes
    /// this an error again, the failure points straight at what was decided (2026-07-22), and the
    /// `live_solids` assertion pins the retire that an early return would otherwise skip.
    #[test]
    fn a_disjoint_common_is_empty_not_an_error() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        let solids = boolean(&mut m, BoolKind::Common, a, b).expect("empty is not a failure");
        assert!(solids.is_empty());
        assert!(
            m.live_solids.is_empty(),
            "a successful boolean consumes its operands"
        );
    }

    /// An empty boolean reaches the caller as an empty solid list, not an error — the op layer
    /// passes the kernel's answer through rather than reinterpreting it.
    #[test]
    fn boolean_op_passes_an_empty_result_through() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        // Offset in all axes so no faces are coplanar with A.
        let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
        let out = apply(
            &mut m,
            &Operation::Boolean {
                kind: BoolKind::Common,
                a,
                b,
            },
        )
        .unwrap();
        assert_eq!(out, OpOutput::Boolean { solids: vec![] });
        assert!(m.live_solids.is_empty(), "both operands are consumed");
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
        assert!(boolean(&mut m, BoolKind::Common, a, b).unwrap().is_empty());
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
    fn a_corner_flush_common_keeps_the_non_convex_overlap() {
        // A **corner-flush** `Common`: the L-prism and the box both start at the origin, so
        // **three** of their face planes coincide — `z = 0` (both floors), `x = 0`, `y = 0`. Every
        // vertex of the shared corner lies exactly on the other solid's face planes, which is what
        // the old reject tag said: `VERTEX_ON_FACE_PLANE`. The arrangement engine names such a
        // point by its plane triple like any other, so the configuration is no longer special.
        //
        // Not covered by the other two non-convex `Common` locks:
        // `common_non_convex_overlap_is_their_intersection` (l_and_corner_box) and
        // `common_non_convex_containment_is_inner` both meet transversally, with no coplanar pair.
        //
        // Hand-checked shape, not just volume. The overlap is the L
        // `x∈[0,1.5]×y∈[0,1]` (1.5) plus `x∈[0,1]×y∈[1,1.5]` (0.5) = 2.0, over `z∈[0,0.5]`:
        //   volume 2.0 · 0.5 = 1.0
        //   area   2 · 2.0 (caps) + 6.0 (the L's perimeter) · 0.5 = 7.0
        // and the L has six sides, so eight faces.
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
        let r = boolean_one(&mut m, BoolKind::Common, lsolid, b).expect("corner-flush Common");
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, r).unwrap();
        assert!((p.volume - 1.0).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 7.0).abs() < 1e-12, "area {}", p.area);
        assert_eq!(m.solids.get(r).cavities.len(), 0);
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
