//! Feature operations (design §7): the public sketch/extrude/pad/pocket API and the `apply`/
//! `replay` driver. The top layer — it composes the boolean engine ([`crate::boolean`]) and rigid
//! transform ([`crate::transform`]) over the plane substrate below.

use crate::BoolError;
use crate::boolean::boolean;
use crate::exact::Swept;
use crate::planes::outer_tri;
use crate::transform::transform;
use nacre_geom::intersect::{
    RingSide, drop_collinear_midpoints, plane_side, point_in_ring_2d_rat,
    ring_self_intersection_rat, rings_cross_rat,
};
use nacre_geom::{Curve, Line, Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::{Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, MotionNode, Orientation, Origin, Shell, Solid, SurfaceDef,
    Vertex, VertexDef,
};

/// A sketch-plane frame: a 2-D point `(u, v)` maps to `origin + u·x + v·y`.
/// The axes are unit and orthogonal (the constructors ensure it); the normal is `x × y`.
///
/// ★★★★★ **The fields are private, and that is the whole point.** They used to be `pub`, so a
/// caller handed the kernel three *normalized* f64 vectors — and normalizing is where the
/// exactness dies: a plane with normal `(1, 1, 1)` has coefficients `[1, 1, 1, 0]`, three
/// integers, but its unit axes square to `0.9999999999999999…` and no exact form survives. The
/// kernel then had nothing to build on and dropped the whole prism to f64.
///
/// So a plane is built through a constructor that **keeps what the caller stated**
/// ([`PlaneDef`]), and the axes below are the *realization* of that. The two cannot describe
/// different planes because only one of them is written down.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SketchPlane {
    // `pub(crate)`, not `pub`: the constructors live in this crate's `lib.rs`, and what has to be
    // closed is the **public** surface — a caller outside must not be able to hand in three
    // normalized axes and call that a plane.
    pub(crate) origin: Point3,
    pub(crate) x_axis: Vector3,
    pub(crate) y_axis: Vector3,
    /// What the caller stated, exactly — `None` when they stated only axes (a frame the kernel
    /// cannot reconstruct, which then takes the f64 path it always took).
    pub(crate) def: Option<PlaneDef>,
}

/// **A sketch plane as its author stated it** — the exact truth behind [`SketchPlane`]'s f64 axes.
///
/// ★★★ **One field, and the invariants are structural.** This used to carry coefficients, an
/// origin, and a reference direction as three halves that every constructor had to keep agreeing
/// ("an origin that is not on `coeffs` is a definition describing two different planes", which
/// cost 0.04 of volume the one time it happened). Now the definition is the three points alone:
///
/// - the sketch's `(0, 0)` **is** `points[0]`,
/// - `+u` **is** `points[1] − points[0]` — a difference of two points of the plane, so it lies in
///   the plane by definition,
/// - the normal's direction is `(p1 − p0) × (p2 − p0)` — the point order carries the polarity.
///
/// Nothing is left to check, and nothing can disagree. The canonical coefficients are *derived*
/// (`nacre_scalar::plane_name_exact` — total since S2, `Narrow | Wide`), which also removes the
/// old constructors' failure class "the coefficients do not fit `i128`": three in-window points
/// always name their plane, however wide its canonical form.
///
/// ★ `ref_dir()` is **not** a unit vector and is not projected; the normalization a frame needs
/// is exactly one `1/√(rational)` at realization time — never something stored.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneDef {
    /// Three points of the plane, non-collinear (`plane_name_exact` is what verified it — every
    /// constructor rejects a collinear triple as "no plane"). `points[0]` is the sketch origin,
    /// `points[1] − points[0]` the `+u` direction, and the order fixes the normal's sign.
    pub(crate) points: [[nacre_scalar::Rat; 3]; 3],
}

impl PlaneDef {
    /// The three defining points — the sketch origin first, then the point `+u` runs toward,
    /// then the point fixing the normal's side.
    pub fn points(&self) -> [[nacre_scalar::Rat; 3]; 3] {
        self.points
    }

    /// Where the sketch's `(0, 0)` sits — the first defining point.
    pub fn origin(&self) -> [nacre_scalar::Rat; 3] {
        self.points[0]
    }

    /// The `+u` direction, in the plane, not unit length — `points[1] − points[0]`.
    ///
    /// The subtraction cannot overflow: both points passed through a constructor, and every
    /// constructor either lifted decimals (narrow) or added one lifted decimal to another —
    /// widths nowhere near `i128`'s ceiling.
    pub fn ref_dir(&self) -> [nacre_scalar::Rat; 3] {
        core::array::from_fn(|i| {
            self.points[1][i]
                .checked_sub(self.points[0][i])
                .expect("constructor-bounded widths")
        })
    }
}

/// One closed ring of a [`Profile2d`], stored as its rational truth.
///
/// The coordinates are what the author's decimals *spelled* (`Rat::from_decimal`), not the f64s
/// that carried them — the same truth/cache split every dimension in the kernel gets
/// (`docs/truth-and-cache.md`). [`Ring2d::realized`] is the f64 cache, and the round-trip is
/// bit-preserving (`from_decimal(x).to_f64() == x`), so the realization is exactly the f64 the
/// caller handed in.
#[derive(Clone, Debug, PartialEq)]
pub struct Ring2d {
    points: Vec<[Rat; 2]>,
}

impl Ring2d {
    /// The ring's rational truth — normalized (flat corners dissolved), in author order.
    pub fn points(&self) -> &[[Rat; 2]] {
        &self.points
    }

    /// The f64 realization of [`points`](Ring2d::points) — the cache the fallback placement and
    /// diagnostics consume. Bit-identical to what the caller passed, point for point that
    /// survived normalization.
    pub fn realized(&self) -> Vec<Point2> {
        self.points
            .iter()
            .map(|p| Point2::from_array([p[0].to_f64(), p[1].to_f64()]))
            .collect()
    }
}

/// A closed planar region: one outer ring and any number of hole rings (straight segments only,
/// at least 3 points each).
///
/// **Winding is not the caller's business, and not this type's either.** The rings are stored in
/// author order; the prism builder is the single place that decides orientation, because it is
/// the only one that knows the sweep direction (a pocket sweeps *into* a face, which flips what
/// "counter-clockwise" means). Fixing the winding here as well would put that decision in two
/// places, which is exactly how the holes and the outer ring come to disagree.
///
/// **The constructors take the boundary f64s and keep the truth** (`docs/truth-and-cache.md`):
/// each coordinate becomes the rational its shortest decimal spells, a dimension outside the
/// decimal window (`~1e38` above, `~1e-22` below for a 17-digit value) is a named error
/// **at construction** rather than a silent f64 fallback downstream, and every **flat corner**
/// (a vertex strictly mid-run on a straight edge) is dissolved — its two walls would be one
/// plane, so the corner it names has no three-plane definition, and deleting it changes no
/// geometry. What survives construction is the profile's normal form.
///
/// The constructors do **not** run [`Profile2d::check`]; that contract (simplicity, disjoint
/// rings, hole containment) is `O(n²)` and **every operation that consumes a profile runs it
/// first**, so the kernel never works from an unverified one. `sketch::from_rings` is the
/// checking constructor for loose rings.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile2d {
    outer: Ring2d,
    holes: Vec<Ring2d>,
}

/// Which ring of a [`Profile2d`] an error is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProfileRing {
    /// The outer ring.
    Outer,
    /// The hole at this index in [`Profile2d::holes`].
    Hole(usize),
}

impl Profile2d {
    /// A simple polygon — no holes. `Err` when a coordinate has no rational truth (outside the
    /// decimal window); see the type doc.
    pub fn polygon(points: Vec<Point2>) -> Result<Profile2d, OpError> {
        Profile2d::with_holes(points, Vec::new())
    }

    /// An outer ring with holes. Lifts every coordinate to its rational truth and dissolves flat
    /// corners; the *contract* ([`Profile2d::check`]) is still enforced by every operation that
    /// consumes the profile.
    pub fn with_holes(outer: Vec<Point2>, inners: Vec<Vec<Point2>>) -> Result<Profile2d, OpError> {
        let lift = |ring: &[Point2], id: ProfileRing| -> Result<Vec<[Rat; 2]>, OpError> {
            ring.iter()
                .enumerate()
                .map(|(i, p)| {
                    let point = |x: f64| Rat::from_decimal(x);
                    match (point(p[0]), point(p[1])) {
                        (Some(x), Some(y)) => Ok([x, y]),
                        _ => Err(OpError::ProfileOutsideDecimalWindow { ring: id, point: i }),
                    }
                })
                .collect()
        };
        let outer = lift(&outer, ProfileRing::Outer)?;
        let holes = inners
            .iter()
            .enumerate()
            .map(|(h, ring)| lift(ring, ProfileRing::Hole(h)))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Profile2d::from_rat_rings(outer, holes))
    }

    /// The constructor behind the f64 boundary — rings already in their rational truth
    /// (`sketch::from_rings` arrives here after lifting and classifying on that truth).
    /// Normalization (the flat-corner dissolve) happens here, once, for every entry path.
    pub(crate) fn from_rat_rings(outer: Vec<[Rat; 2]>, holes: Vec<Vec<[Rat; 2]>>) -> Profile2d {
        let ring = |points: Vec<[Rat; 2]>| Ring2d {
            points: drop_collinear_midpoints(points),
        };
        Profile2d {
            outer: ring(outer),
            holes: holes.into_iter().map(ring).collect(),
        }
    }

    /// Verify the contract: every ring is a **simple polygon** of at least three points, the rings
    /// are pairwise disjoint, and each hole lies inside the outer ring and inside no other hole.
    ///
    /// Every operation that consumes a profile calls this first, so an unverified profile never
    /// reaches the topology. It is public because an app can ask before it builds.
    ///
    /// **Simplicity is a contract on authored input, not an invariant of kernel data.** A drawn
    /// ring's inside is *defined* by even-odd, which needs simplicity to mean anything; the
    /// contours a boolean produces get their inside from the arrangement instead, and those may
    /// legitimately be non-simple (a figure-8 pinch — see `loop_winding`). Do not "unify" the two.
    ///
    /// Exact — on the truth: every decision is an `orient2d` sign on the rationals the author's
    /// decimals spelled, not on their f64 carriers. The two can disagree (three points collinear
    /// in decimal sit a hair off the line in binary), and it is the decimal that is the author's
    /// meaning. Edge indices in the errors refer to the normalized rings [`outer`](Profile2d::outer)
    /// returns — construction dissolved flat corners first.
    ///
    /// `O(n²)` in the ring size, and it runs on every extrude and every replay of one. Measured
    /// on a convex ring of 17-digit coordinates (the worst case — nothing short-circuits;
    /// `profile_check_wall_clock` in `tests/perf.rs`): **34 ms at 100 points, 3.2 s at 1 000** —
    /// each rational sign costs ~1.7 µs (the narrow route's gcd reductions; such coordinates stay
    /// inside `i128`) against the old f64 predicate's nanoseconds. Hand-drawn sketches are tens
    /// of points and pay well under a millisecond; a generator emitting thousands of points per
    /// ring is now firmly outside this function's comfort, and the named follow-ups — an exact
    /// bounding-box prefilter for the segment pairs, or an f64 filter with a sound error bound
    /// escalating to `Rat` — live in `docs/truth-and-cache.md`'s open items, deliberately
    /// unbuilt until that population exists.
    pub fn check(&self) -> Result<(), OpError> {
        let rings = || {
            std::iter::once((ProfileRing::Outer, &self.outer)).chain(
                self.holes
                    .iter()
                    .enumerate()
                    .map(|(i, r)| (ProfileRing::Hole(i), r)),
            )
        };
        for (id, r) in rings() {
            if r.points().len() < 3 {
                return Err(OpError::DegenerateProfile);
            }
            // Simplicity first: `point_in_ring_2d_rat` below is only meaningful on a simple
            // ring. The predicate reports a zero-length edge as a pair with itself.
            if let Some((a, b)) = ring_self_intersection_rat(r.points()) {
                return Err(if a == b {
                    OpError::ZeroLengthProfileEdge { ring: id, edge: a }
                } else {
                    OpError::SelfIntersectingProfile {
                        ring: id,
                        edges: (a, b),
                    }
                });
            }
        }
        let all: Vec<(ProfileRing, &Ring2d)> = rings().collect();
        for (i, (ida, a)) in all.iter().enumerate() {
            for (idb, b) in &all[i + 1..] {
                if rings_cross_rat(a.points(), b.points()) {
                    return Err(OpError::ProfileRingsMeet { a: *ida, b: *idb });
                }
            }
        }
        // The rings are disjoint, so any one vertex answers for a whole ring.
        for (h, hole) in self.holes.iter().enumerate() {
            if point_in_ring_2d_rat(hole.points()[0], self.outer.points()) != RingSide::Inside {
                return Err(OpError::HoleNotInsideOuter { hole: h });
            }
            for (k, other) in self.holes.iter().enumerate() {
                if k != h
                    && point_in_ring_2d_rat(hole.points()[0], other.points()) == RingSide::Inside
                {
                    // A ring inside a hole is an island — material again, so it belongs to a
                    // profile of its own. `sketch::from_rings` is what splits those out.
                    return Err(OpError::NestedHole { outer: k, inner: h });
                }
            }
        }
        Ok(())
    }

    /// The outer ring — the normalized rational truth.
    pub fn outer(&self) -> &Ring2d {
        &self.outer
    }

    /// The hole rings — the normalized rational truth.
    pub fn holes(&self) -> &[Ring2d] {
        &self.holes
    }
}

/// A modelling operation.
///
/// ★ **`Extrude` is the large variant and it is not boxed.** It carries a [`SketchPlane`], which
/// now holds the plane's exact rational definition (~400 bytes). An op log is tens of entries
/// long and is walked once per replay, so the wasted space is measured in kilobytes; boxing would
/// buy that back at the cost of an indirection on the one type a caller constructs by hand.
#[allow(clippy::large_enum_variant)]
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
    /// Reflect `solid` in the coordinate plane `axis = offset`, superseding it (pair with
    /// [`Operation::Copy`] to keep the original — the usual move, since mirroring exists to build
    /// the other half of a symmetric part). Lengths are preserved and handedness is reversed:
    /// this is a reflection, not a negative scale.
    Mirror {
        solid: Handle<Solid>,
        axis: Axis,
        offset: Rat,
    },
    /// Duplicate `solid` in place, keeping the original live — **the only operation that adds to
    /// `live_solids` without removing anything** (design §2). Every other edit supersedes its
    /// input, so this is what makes "cut with the same tool twice", "keep the original and a moved
    /// copy", and pattern/mirror sugar expressible at all.
    Copy { solid: Handle<Solid> },
}

/// Which boolean to compute.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolKind {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

/// A failure while applying an operation.
#[derive(Debug, PartialEq)]
pub enum OpError {
    /// A profile ring with fewer than 3 points.
    DegenerateProfile,
    /// A profile coordinate outside the decimal window (`Rat::from_decimal` — `~1e38` above,
    /// `~1e-22` below for a full-width value), so it has no rational truth for the kernel to
    /// keep. Named **at construction**: the alternative was a silent fall to the f64 path, whose
    /// prism records no exact points and whose surfaces cannot survive a motion undemoted.
    /// CAD dimensions live nowhere near the window's edges; reaching this is an input mistake.
    ProfileOutsideDecimalWindow {
        /// Which ring of the profile.
        ring: ProfileRing,
        /// The offending point's index within that ring, as given.
        point: usize,
    },
    /// A profile ring meets itself — it crosses, touches, doubles back over one of its own edges,
    /// or doubles back over one of its own edges. Such a ring has no unambiguous inside, so the
    /// region it is supposed to bound is undefined. `edges` are the two offending edge indices
    /// within that ring. The input is wrong; this is not a missing capability.
    SelfIntersectingProfile {
        /// Which ring of the profile.
        ring: ProfileRing,
        /// The offending pair of edge indices within that ring.
        edges: (usize, usize),
    },
    /// A profile ring repeats a point, so one of its edges has zero length. It would leave a
    /// degenerate edge in the topology. (Reported apart from
    /// [`OpError::SelfIntersectingProfile`] because "you typed the same point twice" and "your
    /// outline crosses itself" are different mistakes to the author, though the same predicate
    /// finds them.)
    ZeroLengthProfileEdge {
        /// Which ring of the profile.
        ring: ProfileRing,
        /// The zero-length edge's index within that ring.
        edge: usize,
    },
    /// Two rings of a profile touch or cross. A hole must lie strictly inside the outer ring and
    /// strictly outside its siblings; anything else has no unambiguous inside.
    ProfileRingsMeet {
        /// The first ring.
        a: ProfileRing,
        /// The second ring.
        b: ProfileRing,
    },
    /// A hole ring lies outside the outer ring, so it cuts nothing. Left unchecked this is a
    /// *silent* wrong: the prism builds, `validate` is clean, and the volume comes out reduced by
    /// a hole that is not there.
    HoleNotInsideOuter {
        /// The index in [`Profile2d::holes`].
        hole: usize,
    },
    /// A hole ring lies inside another hole. That region is material again — an island — and
    /// belongs to a profile of its own; `sketch::from_rings` splits those out by nesting depth.
    NestedHole {
        /// The containing hole's index.
        outer: usize,
        /// The contained hole's index.
        inner: usize,
    },
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
    /// A `Discovered` vertex or edge of a `Transform`/`Copy` input names a surface that is not one
    /// of that solid's face surfaces, so its exact definition cannot be carried onto the duplicate.
    /// The solid's provenance is inconsistent — declined rather than aborting the kernel, and
    /// unreached in the suite (the definitions and the result faces both name the plane class's
    /// representative surface).
    OriginNotOnSolid,
    /// A `Mirror` input carries curved geometry (a cylindrical face, a circular edge). A
    /// reflection reverses a circle's parametrisation, and which convention a mirrored quadric
    /// should take is a curved-geometry decision, so it is declined rather than guessed.
    MirrorNotPlanar,
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
    /// The duplicate. Unlike every other output, the input stays live alongside it.
    Copy { solid: Handle<Solid> },
    /// The reflected solid (supersedes the input).
    Mirror { solid: Handle<Solid> },
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
        Operation::Copy { solid } => {
            let out = crate::transform::copy(model, *solid)?;
            Ok(OpOutput::Copy { solid: out })
        }
        Operation::Mirror {
            solid,
            axis,
            offset,
        } => {
            let out = crate::transform::mirror(model, *solid, *axis, *offset)?;
            Ok(OpOutput::Mirror { solid: out })
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
pub(crate) fn extrude(
    model: &mut Model,
    plane: &SketchPlane,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    profile.check()?;
    // ★★★★★ **A plane the caller stated is drawn in its own frame.**
    //
    // In world coordinates a tilted plane's axes are irrational and `exact()` declines, so the
    // whole prism used to drop to f64. Inside the plane's frame those axes are `x̂`/`ŷ` and the
    // rational path applies unchanged — the walls come straight from the profile's decimals and
    // the far cap is `w = dist`.
    //
    // ★★ **The gate is `plane.exact()`, the same expression the face path uses.** A plane that
    // lifts in the world takes the world, definition or not; otherwise every axis-aligned model
    // would gain motions it does not need and leave the exact predicate path for nothing.
    //
    // ★ **The frame is the caller's own**: `def.origin` and `def.ref_dir` are what they wrote, so
    // the sketch's `(0, 0)` and `+u` land where they asked. `flip` is measured, not derived — the
    // plane's canonical coefficients carry no direction, and `ŵ` has to face the way the sweep does.
    let frame = plane.def.filter(|_| plane.exact().is_none()).and_then(|d| {
        // ★★ **Built with `−normal`, the sense `build_prism` gives a base cap.** The key is
        // direction-free, so this surface and the base cap intern together — and matching the
        // sense means the face does not have to be spelled `Reversed` to compensate. Measured:
        // pushing `+normal` instead flipped the stored normal on 781 base caps, which
        // `Face::orientation` absorbed correctly but for no reason.
        //
        // ★ Through the sketch origin rather than a ring point, which is also *more* accurate
        // here: the caller's origin is exact, so `d` comes out exact where the ring point's dot
        // product rounds (measured `5.55e-17` against `0`).
        let pl = Plane::from_point_normal(plane.origin(), -plane.normal())?;
        let (h, _) = model.push_surface_with_points(
            Surface::Plane(pl),
            SurfaceDef::Constructed,
            Some(d.points()),
        );
        // The caller stated the pair, so the placement is `Named` (S4) — `Canonical` is for
        // frames nobody named, like a face's.
        let placement = nacre_topo::FramePlacement::Named {
            origin: d.origin(),
            ref_dir: d.ref_dir(),
        };
        let (_, _, _, w) = crate::rotated_vertex::frame_world_basis(model, h, &placement, false)?;
        let n = plane.normal().as_array();
        let flip = (0..3).map(|k| w[k] * n[k]).sum::<f64>() < 0.0;
        Some(model.push_motion(
            nacre_topo::Motion::Frame {
                plane: h,
                placement,
                flip,
            },
            None,
        ))
    });
    let (outer, holes) = swept_profile(model, plane, profile, dist, frame);
    // ★★★★ **The base cap *is* the plane the caller named**, so where they stated it exactly
    // (`PlaneDef`) it can record that plane's own points even when nothing else about the prism
    // can — the walls and far cap are irrational in the world unless `plane.exact()` allows.
    //
    // ★★★ **The points travel, not a `Surface`.** Pre-creating one here would have to guess the
    // f64 plane `build_prism` builds, and it builds it through the *oriented* ring's first point —
    // which reversing the ring can change. Two spellings of one plane whose `d` differs by an ulp
    // is precisely the defect this whole line of work removes, so the geometry stays exactly where
    // it was and only the exact record is added. (Measured: building it from the sketch origin
    // instead moved a volume.)
    build_prism(
        model,
        outer,
        holes,
        plane.normal(),
        None,
        plane.def.as_ref().map(|d| d.points()),
    )
}

/// Place a profile on its plane and sweep it — **exactly where the plane admits it**.
///
/// The rational path is not an optimization: it is what makes `extrude(7.7)` and
/// `extrude(1.1)` then `extrude(6.6)` put their caps on the same plane rather than an
/// ulp apart. Where it does not apply — a frame without an exact form, an i128 overflow
/// in the ring arithmetic — this falls back to the f64 arithmetic that was here before,
/// which is no worse than it was. (A dimension outside the decimal window no longer
/// arrives here at all: the profile constructor names it.)
///
/// **Mapping only — no winding decision, and no containment check.** Forcing the outer
/// ring CCW here would be a second opinion on a question `build_prism` already answers
/// from the sweep, and two opinions is how an outer ring and its holes end up wound the
/// same way (the pocket case, where the sweep runs `−n` and flips the outer ring). A
/// profile may also reach past the face boundary; an overhanging footprint routes to
/// the overhang boolean sidecars, which reject honestly what they do not cover.
fn swept_profile(
    model: &Model,
    plane: &SketchPlane,
    profile: &Profile2d,
    dist: f64,
    frame: Option<Handle<MotionNode>>,
) -> (Swept, Vec<Swept>) {
    if let Some(rings) = crate::exact::prism_rings(model, plane, profile, dist, frame) {
        return rings;
    }
    let sweep = plane.normal() * dist;
    let place = |ring: &Ring2d| {
        Swept::along(
            ring.realized().iter().map(|p| plane.point(*p)).collect(),
            sweep,
        )
    };
    (
        place(profile.outer()),
        profile.holes().iter().map(place).collect(),
    )
}

/// Sweep a profile's rings along `sweep` into a prism solid: caps, side walls, and — for each
/// hole ring — a wall of its own plus an inner loop on each cap.
///
/// **This is the only place winding is decided**, because it is the only place that knows the
/// sweep. The outer ring is normalized CCW **about `sweep`** (area vector dotted with the sweep
/// normal — frame-independent, unlike `proj2`+`signed_area`, which mis-signs when the sweep runs
/// along a negative dominant axis, e.g. a pocket into a `+z` face sweeping `−z`). Each hole is
/// then normalized to the **opposite** sense *of that normalized outer ring*, never of the input:
/// get that backwards and a pocket — where the sweep flips the outer ring — silently produces
/// holes wound the same way as the outer, which is not a hole at all.
///
/// Opposite winding is all a hole needs. The wall quads are built from the ring's own traversal,
/// so a reversed ring yields walls facing into the hole (out of the material), and the cap loops
/// come out opposed to the cap's outer loop, which is what makes them holes.
///
/// All vertices are `Origin::Constructed`. Returns the solid and its faces: `faces[0]` = base cap
/// (at the ring, normal `−ŝ`), `faces[1]` = far cap, then the outer walls, then each hole's walls.
/// Shared by [`extrude`] (a boss) and the pocket (`sweep = −n`).
pub(crate) fn build_prism(
    model: &mut Model,
    outer_ring: Swept,
    inner_rings: Vec<Swept>,
    normal: Vector3,
    base_cap_surface: Option<Handle<Surface>>,
    // ★ The world points of the plane the caller named, when they named one. Its canonical name is
    // derived from these, so there is no second half that could travel separately.
    base_cap_points: Option<[[nacre_scalar::Rat; 3]; 3]>,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if outer_ring.base.len() < 3 || inner_rings.iter().any(|h| h.base.len() < 3) {
        return Err(OpError::DegenerateProfile);
    }

    let outer_pts = oriented_ring(outer_ring, normal, true);
    // `false` = opposite to the normalized outer ring, whichever way that ended up.
    let hole_pts: Vec<Swept> = inner_rings
        .into_iter()
        .map(|h| oriented_ring(h, normal, false))
        .collect();

    // ★★★ **Surfaces before topology.** A vertex is defined by the three faces that meet at it,
    // and `Store` is append-only, so the handles have to exist before the vertex does. The two
    // arenas are separate, so interleaving them differently does not shift either one's numbering
    // — but the surfaces' order *among themselves* is what the numbering depends on, and it is
    // preserved exactly: base cap, then top cap, then walls, outer ring before the holes.

    // Base cap: outward normal −N, loops reversed.
    // When padding/pocketing on a face, reuse that face's `Surface` handle (explicit sharing) so
    // the flush contact is a shared-handle coplanar pair the boolean can recognize by `Handle`
    // identity; otherwise push a fresh plane. The materialized outward normal must stay −N, so the
    // face orientation is chosen from the shared surface's stored normal — `surface` and
    // `orientation` travel together, and the reconstruction copies both.
    let (base_surface, base_orient) = match base_cap_surface {
        Some(h) => {
            let n_h = match model.surface(h) {
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
            // ★★★★ **The caller's statement wins here, and the frame's is the fallback.**
            //
            // The base cap *is* the plane they named, and they named it in the world — so stating
            // it that way keeps it `Constructed`, keeps its judgment exact, and lets two extrudes
            // on one plane share it **whatever frames they chose**. Writing it as `[0,0,1,0]` in
            // this prism's frame instead would be a second exact description of one plane, under a
            // different `SurfaceKey` — the very duplication this work removes.
            //
            // The frame's answer is what a plane with no caller statement gets (an axis-aligned
            // sketch, where the two agree anyway).
            //
            // ★★★★★ **One frame or the other, and the `def` goes with the points.**
            // A plane's points and its `SurfaceDef` are two halves of one statement: `Constructed`
            // means *"these speak about the world"* and `Moved` means *"about the pre-motion
            // frame"*. There used to be a third half — coefficients, chosen by their own `or_else`,
            // so a caller whose points overflowed while their coefficients did not got world
            // coefficients beside the prism's **frame** ring. That half no longer exists: the name
            // is derived from whichever points are recorded, so the two can no longer come apart.
            //
            // ★ The frame's `surface_def()` is the same one the top cap takes, and it agrees with
            // `Constructed` wherever there is no motion — which is every case a missing caller
            // statement can produce.
            let (base_def, cap_pts) = match base_cap_points {
                Some(p) => (SurfaceDef::Constructed, Some(p)),
                None => {
                    let e = outer_pts.exact.as_ref();
                    (
                        e.map_or(SurfaceDef::Constructed, |e| e.surface_def()),
                        e.and_then(|e| e.cap_points(false)),
                    )
                }
            };
            let (s, flipped) = model.push_surface_with_points(
                Surface::Plane(
                    Plane::from_point_normal(outer_pts.base[0], -normal)
                        .ok_or(OpError::DegenerateGeometry)?,
                ),
                base_def,
                cap_pts,
            );
            // The plane was built with `−N` as its normal, so `Forward` is what states an outward
            // `−N` — unless a shared surface points the other way, which `flipped` reports.
            let orient = if flipped {
                Orientation::Forward.flipped()
            } else {
                Orientation::Forward
            };
            (s, orient)
        }
    };
    // Top cap: outward normal +N.
    let top_def = match outer_pts.exact.as_ref() {
        Some(e) if e.top.len() >= 3 => e.surface_def(),
        _ => SurfaceDef::Constructed,
    };
    let top_points = outer_pts.exact.as_ref().and_then(|e| e.cap_points(true));
    let (top_surface, top_flipped) = model.push_surface_with_points(
        Surface::Plane(
            Plane::from_point_normal(outer_pts.top[0], normal)
                .ok_or(OpError::DegenerateGeometry)?,
        ),
        top_def,
        top_points,
    );
    let top_orient = if top_flipped {
        Orientation::Forward.flipped()
    } else {
        Orientation::Forward
    };

    // Wall surfaces, in the same order the faces will be emitted: the outer ring's, then each
    // hole's (facing into the hole).
    let outer_walls = wall_surfaces(model, &outer_pts)?;
    let hole_walls: Vec<Vec<(Handle<Surface>, bool)>> = hole_pts
        .iter()
        .map(|h| wall_surfaces(model, h))
        .collect::<Result<_, _>>()?;

    // ── Topology. Every surface it needs already exists.
    let caps = (base_surface, top_surface);
    let outer = sweep_ring(model, &outer_pts, &outer_walls, caps)?;
    let holes: Vec<RingCells> = hole_pts
        .iter()
        .zip(hole_walls.iter())
        .map(|(h, w)| sweep_ring(model, h, w, caps))
        .collect::<Result<_, _>>()?;

    let mut faces =
        Vec::with_capacity(2 + outer.len() + holes.iter().map(|h| h.len()).sum::<usize>());
    faces.push(model.faces.push(Face {
        surface: base_surface,
        outer: outer.cap_loop(Cap::Base),
        inner: holes.iter().map(|h| h.cap_loop(Cap::Base)).collect(),
        orientation: base_orient,
    }));
    faces.push(model.faces.push(Face {
        surface: top_surface,
        outer: outer.cap_loop(Cap::Top),
        inner: holes.iter().map(|h| h.cap_loop(Cap::Top)).collect(),
        orientation: top_orient,
    }));
    for (ring, walls) in
        std::iter::once((&outer, &outer_walls)).chain(holes.iter().zip(hole_walls.iter()))
    {
        ring.push_walls(model, walls, &mut faces);
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

/// Which cap a ring's loop is being built for.
#[derive(Clone, Copy)]
enum Cap {
    Base,
    Top,
}

/// One ring swept into cells: the two rings of vertices and the three edge families that join
/// them. Built the same way for the outer ring and for a hole — the difference is only which way
/// the ring runs, which the caller has already decided.
struct RingCells {
    /// Kept for [`RingCells::len`]; the geometry itself is read off the `Swept` these came from,
    /// which is also where the wall planes were built (`wall_surfaces`).
    base_pts: Vec<Point3>,
    be: Vec<Handle<Edge>>, // base  B_i -> B_{i+1}
    te: Vec<Handle<Edge>>, // top   T_i -> T_{i+1}
    ve: Vec<Handle<Edge>>, // riser B_i -> T_i
}

impl RingCells {
    fn len(&self) -> usize {
        self.base_pts.len()
    }

    /// The cap loop for this ring. The base cap faces `−ŝ`, so its loops run backwards.
    fn cap_loop(&self, cap: Cap) -> Loop {
        let n = self.len();
        match cap {
            Cap::Base => Loop {
                half_edges: (0..n)
                    .rev()
                    .map(|i| HalfEdge {
                        edge: self.be[i],
                        forward: false,
                    })
                    .collect(),
            },
            Cap::Top => Loop {
                half_edges: (0..n)
                    .map(|i| HalfEdge {
                        edge: self.te[i],
                        forward: true,
                    })
                    .collect(),
            },
        }
    }

    /// One quad per ring segment. The quad's winding follows the ring's, so a ring wound against
    /// the outer one yields walls whose normals point into the hole.
    fn push_walls(
        &self,
        model: &mut Model,
        walls: &[(Handle<Surface>, bool)],
        faces: &mut Vec<Handle<Face>>,
    ) {
        let n = self.len();
        for (i, &(surface, flipped)) in walls.iter().enumerate().take(n) {
            let j = (i + 1) % n;
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: self.be[i],
                        forward: true,
                    },
                    HalfEdge {
                        edge: self.ve[j],
                        forward: true,
                    },
                    HalfEdge {
                        edge: self.te[i],
                        forward: false,
                    },
                    HalfEdge {
                        edge: self.ve[i],
                        forward: false,
                    },
                ],
            };
            faces.push(model.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                // The quad's winding is the ring's, which is what the plane above was built from;
                // a shared surface pointing the other way spells the same outward as `Reversed`.
                orientation: if flipped {
                    Orientation::Forward.flipped()
                } else {
                    Orientation::Forward
                },
            }));
        }
    }
}

/// One plane per ring segment, pushed **before** any of the ring's topology exists — see
/// `build_prism`. The same three points `push_walls` used to build them from, and in the same
/// order, so the surface arena's numbering is untouched.
fn wall_surfaces(model: &mut Model, ring: &Swept) -> Result<Vec<(Handle<Surface>, bool)>, OpError> {
    let n = ring.base.len();
    (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            // The witness is the same three points **in the frame the coefficients are written
            // in** — the world's own points when there is no frame, so this is unchanged there.
            let def = match ring.exact.as_ref() {
                Some(e) => e.surface_def(),
                None => SurfaceDef::Constructed,
            };
            Ok(model.push_surface_with_points(
                Surface::Plane(
                    Plane::through_points(ring.base[i], ring.base[j], ring.top[i])
                        .ok_or(OpError::DegenerateGeometry)?,
                ),
                def,
                // ★ The same three points the f64 plane above is built through, in rationals — so
                // this wall and any other face of the same plane record one array and derive one
                // name. `None` here is the f64 path, where there are no rationals at all.
                ring.exact.as_ref().map(|e| e.wall_points(i)),
            ))
        })
        .collect()
}

/// `pts` wound counter-clockwise about `normal` when `ccw`, clockwise when not. The test is the
/// polygon's area vector against `normal`, so it does not care which axis dominates.
/// **Requires a nonzero area.** The winding is read from the sign of the area vector, and a ring
/// that encloses nothing gives zero — the comparison below would then pick a side by accident.
/// [`Profile2d::check`] is what guarantees it: a simple polygon cannot have zero area, and a ring
/// that folds back on itself (a symmetric bowtie cancels to exactly zero) is not simple.
fn oriented_ring(ring: Swept, normal: Vector3, ccw: bool) -> Swept {
    let v = &ring.base;
    let k = v.len();
    let area_vec = (0..k)
        .map(|i| (v[i] - Point3::origin()).cross(v[(i + 1) % k] - Point3::origin()))
        .fold(Vector3::from_array([0.0; 3]), |a, b| a + b);
    if (area_vec.dot(normal) < 0.0) == ccw {
        ring.reversed()
    } else {
        ring
    }
}

/// Push one ring's vertices and edges (base ring, top ring, risers).
///
/// The top ring arrives already computed rather than being derived here as
/// `base + sweep`: where the frame allows it that arithmetic is done in exact rationals
/// (see [`crate::exact`]), and a dimension split into two then lands on the same points
/// as the undivided one instead of an ulp away.
fn sweep_ring(
    model: &mut Model,
    ring: &Swept,
    walls: &[(Handle<Surface>, bool)],
    caps: (Handle<Surface>, Handle<Surface>),
) -> Result<RingCells, OpError> {
    let n = ring.base.len();
    let base_pts: Vec<Point3> = ring.base.clone();
    let top_pts: Vec<Point3> = ring.top.clone();
    // Corner `i` is where the wall before it, the wall after it, and the cap meet.
    //
    // ★ **Two walls that are one plane would name a line, not a point** — and since surfaces are
    // interned, "one plane" *is* "one handle", so the check is a comparison. Since S3 this is a
    // guard on an invariant, not a live case: the only producer of adjacent same-plane walls was
    // a profile with a collinear midpoint, and `Profile2d`'s constructor now dissolves those, so
    // every corner gets its three-plane definition (locked by
    // `a_collinear_midpoint_profile_builds_its_clean_twin_bit_for_bit`). Non-adjacent walls may
    // still legitimately share a plane (a notch), which never lands `prev == here`.
    let define = |i: usize, cap: Handle<Surface>| -> Option<VertexDef> {
        let prev = walls[(i + n - 1) % n].0;
        let here = walls[i].0;
        (prev != here && prev != cap && here != cap)
            .then_some(VertexDef::ThreePlane([prev, here, cap]))
    };
    // ★★★★ **A vertex drawn in a frame is `Moved`, not `Constructed`.**
    //
    // `Constructed` means *"this f64 coordinate is the truth"*, which is exactly right in the
    // world and exactly wrong here: the truth is a rational `(u, v, w)` in the plane's frame, and
    // the world coordinate is what realizing it produced. Recording it as `Moved` against the
    // frame is what lets a judge replay the definition and check the coordinate against it —
    // and what lets `shared_base` cancel the whole sketch's frame and judge it *exactly*.
    //
    // The base vertex is the sketch point itself, at its frame coordinate. It belongs to no
    // shell — the store is append-only and a vertex nothing references is simply a definition
    // that outlives its use, which is what `Origin::Moved` needs one of.
    let motion = ring.exact.as_ref().and_then(|e| e.motion);
    let push_verts = |model: &mut Model,
                      ps: &[Point3],
                      frame_pt: &dyn Fn(usize) -> Option<Point3>,
                      cap: Handle<Surface>|
     -> Vec<Handle<Vertex>> {
        ps.iter()
            .enumerate()
            .map(|(i, p)| {
                let origin = match (motion, frame_pt(i)) {
                    (Some(m), Some(fp)) => {
                        let base = model.vertices.push(Vertex {
                            point: fp,
                            origin: Origin::Constructed,
                            definition: None,
                        });
                        Origin::Moved { base, motion: m }
                    }
                    _ => Origin::Constructed,
                };
                model.vertices.push(Vertex {
                    point: *p,
                    origin,
                    definition: define(i, cap),
                })
            })
            .collect()
    };
    let ex = ring.exact.as_ref();
    let bv = push_verts(model, &base_pts, &|i| ex.map(|e| e.base_f64(i)), caps.0);
    let tv = push_verts(model, &top_pts, &|i| ex.map(|e| e.top_f64(i)), caps.1);

    let (mut be, mut te, mut ve) = (Vec::new(), Vec::new(), Vec::new());
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
    Ok(RingCells {
        base_pts,
        be,
        te,
        ve,
    })
}

/// A planar face's live solid, its in-plane right-handed frame (`x × y = n`, centred on the face
/// centroid so a profile's `(0,0)` lands there), and its loops — the shared setup for placing a
/// profile on a face (pad / pocket).
/// **Which frame a face's sketch lives in** — the plane, its [`nacre_topo::FramePlacement`]
/// (a face has no caller to name one, so it is always `Canonical` today), and whether the
/// plane's canonical coefficients need negating to face the way the face does. Everything a
/// [`nacre_topo::Motion::Frame`] node needs, before the model has one.
#[derive(Clone, Copy, Debug)]
struct SketchFrame {
    plane: Handle<Surface>,
    placement: nacre_topo::FramePlacement,
    flip: bool,
}

struct FaceFrame {
    solid_h: Handle<Solid>,
    surface_h: Handle<Surface>,
    n: Vector3, // outward normal
    x: Vector3,
    y: Vector3,
    origin: Point3,
    /// ★ Set when this face's sketch lives in its plane's **own frame** rather than the world:
    /// the plane to take the frame from, and which way round. `x`/`y`/`origin` above are then that
    /// frame's, realized — so what a caller is told and what the operation builds are one thing.
    sketch_frame: Option<SketchFrame>,
}

/// **Which way is "right" and "up" on a face pointing `n`** — the `(u, v)` a sketch frame takes.
///
/// This is the **arbitrary-axis convention** (DXF/AutoCAD, and what most CAD puts on a face):
/// cross the world `ẑ` into the normal, unless the normal *is* vertical, in which case cross `ŷ`.
///
/// ★ **Not [`nacre_math::Vector::any_perpendicular`], and the difference is the point.** That one
/// answers a question of fact — *"give me a unit vector perpendicular to this"* — by crossing in
/// whichever world axis the normal is **least** aligned with, which keeps the cross far from zero
/// and is exactly right for a cylinder's seam or a STEP `ref_dir`. It is the wrong answer for a
/// frame a person draws in, because the axis it picks changes with the normal: a box lid comes out
/// `u = −ŷ, v = +x̂`, ninety degrees from world XY and from [`SketchPlane::world_xy`] itself.
///
/// What this convention buys, on top of the lid agreeing with `world_xy`:
///
/// * **On every face that is not horizontal, `v` points up.** `u = ẑ × n` is horizontal, so
///   `v·ẑ = (n × (ẑ × n))·ẑ = 1 − n_z²`, which is positive unless `n` is vertical. Sketching on a
///   wall, "up" is up.
/// * **Axis-aligned faces keep axes in `{0, ±1}`**, so `SketchPlane::exact` still fires and the
///   rational construction path is not lost.
///
/// ★★ **The branch is exact, not toleranced.** DXF switches on `|n_x| < 1/64`, a threshold only a
/// float-only kernel needs; `n_x == 0 && n_y == 0` is the real question and this kernel can ask it.
/// Nor is the non-vertical branch fragile near vertical: `ẑ × n = (−n_y, n_x, 0)` is a plain
/// rotation of `(n_x, n_y)` into the plane, with no cancellation to lose digits to.
///
/// ★ **A discontinuity is unavoidable and this convention chooses where to put it** — no continuous
/// tangent frame exists on the sphere. `any_perpendicular` breaks along whole arcs (wherever two
/// components tie for smallest, e.g. `(0.5, 0.5, 0.707)`, far from any pole); this breaks at the
/// two poles only. It is not smooth *at* the poles either — approaching `+ẑ` from different sides
/// gives different limits — but the set where that happens is two points instead of three arcs.
///
/// Returns both axes rather than just `u`, so the two call sites cannot disagree about which way
/// `v` runs. `None` only for the zero vector.
pub(crate) fn frame_axes(n: Vector3) -> Option<(Vector3, Vector3)> {
    let [nx, ny, nz] = n.as_array();
    let raw = if nx == 0.0 && ny == 0.0 {
        // ŷ × n for a vertical normal: (n_z, 0, 0).
        Vector3::from_array([nz, 0.0, 0.0])
    } else {
        // ẑ × n.
        Vector3::from_array([-ny, nx, 0.0])
    };
    // `-0.0` is worth nothing to keep: the negations above produce it whenever a component is
    // zero, and a frame reported to a caller as `[-0.0, 1.0, 0.0]` invites a double-take for no
    // reason. Adding zero is the identity on every other value.
    let tidy = |v: Vector3| Vector3::from_array(v.as_array().map(|c| c + 0.0));
    let u = tidy(raw.normalize()?);
    Some((u, tidy(n.cross(u))))
}

/// The sketch plane of a planar face — **the very frame [`Operation::PadOnFace`] and
/// [`Operation::PocketOnFace`] place their profile in**, so a caller can work out where its
/// `(0, 0)` will land before it builds anything.
///
/// That equality is the contract, not a coincidence: this is a projection of the frame those
/// operations use, never a second derivation. A test pins a hand-placed profile against a pad to
/// keep it that way.
///
/// `NonPlanarFace` for a curved surface (only a plane carries a frame); `FaceNotInLiveSolid` if no
/// live solid's outer shell holds the face.
pub fn face_plane(model: &Model, face: Handle<Face>) -> Result<SketchPlane, OpError> {
    let f = face_frame(model, face)?;
    Ok(realized_plane(f.origin, f.x, f.y))
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
    let plane = match model.surface(surface_h) {
        Surface::Plane(p) => *p,
        Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    };
    let sign = match orientation {
        Orientation::Forward => 1.0,
        Orientation::Reversed => -1.0,
    };
    let n = plane.normal() * sign;
    let (x, y) = frame_axes(n).ok_or(OpError::DegenerateGeometry)?;
    // ★ The origin is the **world origin projected onto the face's plane** — a property of the
    // plane, not of the face.
    //
    // It used to be the face region's area centroid, chosen over the mean of the outer loop's
    // corners because a centroid does not move when a vertex is added along a straight edge. Both
    // are computed in `f64` from the face's own vertices, though, and `exact.rs` lifts the frame
    // origin with `Rat::from_decimal` — so a rounded cache became the truth. Padding one footprint
    // twice then placed the second profile an ulp from the first and left faces of area `2.2e-16`
    // that `validate` did not report.
    //
    // The projection is `(−d / n·n)·n` from the plane's rational coefficients: one division, no
    // f64 in the derivation, and invariant under negating or scaling those coefficients — so two
    // faces of one plane cannot disagree about where `(0, 0)` is. Measured over the suite: for
    // every `Constructed` surface the realized point lies on the f64 plane at distance exactly `0`
    // (1121/1121), which the centroid did not always manage.
    //
    // Only `SurfaceDef::Constructed` coefficients are world truth. A `Moved` surface records its
    // **pre-motion** frame, so projecting those gives a pre-motion point — measured `0.29` away
    // from the world plane, not a rounding but a different place. Those keep the f64 projection,
    // which is the same rule computed from the description that is available.
    // `narrow()` gates the wide vessel out: a `Wide` name (S2) carries identity only, so it
    // keeps the f64 projection exactly as a missing name did.
    let origin = match (
        model.surface_defs.get(&surface_h),
        model.surface_name.get(&surface_h).and_then(|n| n.narrow()),
    ) {
        (Some(SurfaceDef::Constructed), Some(&c)) => nacre_scalar::plane_origin_projection(c)
            .map(|p| Point3::from_array([p[0].to_f64(), p[1].to_f64(), p[2].to_f64()]))
            .unwrap_or_else(|| plane.project(Point3::origin())),
        _ => plane.project(Point3::origin()),
    };
    // ★★★★★ **When the operation will sketch in the plane's own frame, report *that* frame.**
    //
    // `face_plane`'s contract is that it names the frame `PadOnFace` places a profile in, and a
    // tilted face is about to be sketched in its plane's frame rather than in world coordinates.
    // Reporting the world axes here and using the frame's there would put a caller's profile a
    // quarter turn from where it asked for it — measured, as a boss that missed its own face.
    //
    // ★★ **The gate is the same one the operation uses**, and it has to be the same expression,
    // not the same intent: take the frame only when the world axes do not lift to exact
    // orthonormal rationals. Axis-aligned faces therefore never go near it and are untouched.
    //
    // ★ **`flip` is measured, not derived** — see `frame_world_basis`. `flip = true` negates `ŵ`
    // and `û` together and leaves `v̂`, so the second reading is a sign change rather than a
    // second realization.
    // ★★ A face has no caller to name a frame, so its placement is `Canonical` (S4) — derived
    // when the chain is flattened, stored nowhere. That is also what opens this branch for a
    // plane whose name is `Wide` or whose canonical values overflow `i128`: `frame_world_basis`
    // succeeds through the arbitrary-precision road where the old narrow derivation declined.
    let world = realized_plane(origin, x, y);
    let sketch = (world.exact().is_none())
        .then(|| {
            crate::rotated_vertex::frame_world_basis(
                model,
                surface_h,
                &nacre_topo::FramePlacement::Canonical,
                false,
            )
        })
        .flatten()
        .map(|(o, u, v, w)| {
            // ★ `flip = true` negates `ŵ` and `û` together and leaves `v̂` — a half-turn about
            // `v` — so the second reading is a sign change rather than a second realization.
            let flip = (0..3).map(|k| w[k] * n.as_array()[k]).sum::<f64>() < 0.0;
            let sgn = if flip { -1.0 } else { 1.0 };
            (
                flip,
                Point3::from_array(o),
                Vector3::from_array(u.map(|c| c * sgn)),
                // ★★ **`v̂` as realized, not as `ŵ × û` recomputed here.** It has its own exact
                // rational form (`plane_frame`), so realizing it costs one rounding where a cross
                // product costs two that do not cancel — measured, a wall whose `v` is exactly
                // `ẑ` came back three ulps short of `1.0` through the cross product.
                Vector3::from_array(v),
            )
        });
    let (x, y, origin, sketch_frame) = match sketch {
        Some((flip, o, u, v)) => (
            u,
            v,
            o,
            Some(SketchFrame {
                plane: surface_h,
                placement: nacre_topo::FramePlacement::Canonical,
                flip,
            }),
        ),
        None => (x, y, origin, None),
    };
    Ok(FaceFrame {
        solid_h,
        surface_h,
        n,
        x,
        y,
        origin,
        sketch_frame,
    })
}

/// A [`SketchPlane`] from a **realized** (rounded) frame basis — kernel-internal, and
/// deliberately without a definition.
///
/// ★★ **Do not lift these axes.** The public [`SketchPlane::from_axes`] lifts what a *caller*
/// wrote, because written decimals are a statement. A realized basis is a cache of a plane that
/// already has an exact definition (points + motion); lifting it would mint a second,
/// ulp-different "truth" for the same wall — two exact descriptions of one plane, the defect
/// class the frame work exists to remove — and a sketch built on that lift would land on a
/// non-interned plane a hair off the face it means.
fn realized_plane(origin: Point3, x_axis: Vector3, y_axis: Vector3) -> SketchPlane {
    SketchPlane {
        origin,
        x_axis,
        y_axis,
        def: None,
    }
}

/// A face-local feature built as **tool body + boolean**: the profile
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
    profile.check()?;
    let frame = face_frame(model, face)?;
    // No containment check — an overhanging footprint routes to the overhang boolean sidecars.
    let n = frame.n;
    // Cut carves inward, Fuse raises outward; either way the prism's near cap is flush on the face.
    let signed = if matches!(kind, BoolKind::Cut) {
        -dist
    } else {
        dist
    };
    // The face's own frame, so a pad or pocket takes the same exact-rational path an
    // extrude does: its axes are `{0, ±1}` exactly whenever the face is axis-aligned.
    let plane = realized_plane(frame.origin, frame.x, frame.y);
    // ★★★ **The frame is not decided here — `face_frame` already decided it**, and that is the
    // point: `face_plane` promises a caller the frame this operation will use, so there must be
    // exactly one place that picks it. All that is left is to name it as a motion node.
    //
    // ★★ **`push_motion` interns**, so two sketches on one face name the *same* node — which is
    // what makes their surfaces intern too (`SurfaceKey` is `(name, motion)`) and is the
    // whole point of the exercise: two routes to one height become one `Handle<Surface>` at
    // construction, with no f64 comparison anywhere. With `Canonical` placement the node is
    // `(plane, Canonical, flip)` — nothing per-sketch in the key at all.
    let sketch_frame = frame.sketch_frame.map(|f| {
        model.push_motion(
            nacre_topo::Motion::Frame {
                plane: f.plane,
                placement: f.placement,
                flip: f.flip,
            },
            None,
        )
    });
    let (outer, holes) = swept_profile(model, &plane, profile, signed, sketch_frame);
    let (prism, prism_faces) = build_prism(
        model,
        outer,
        holes,
        n * signed.signum(),
        Some(frame.surface_h),
        // The pad reuses the face's own surface, so it pushes no base cap and states nothing.
        None,
    )?;
    let (solids, class_of) =
        crate::boolean::boolean_with_classes(model, kind, frame.solid_h, prism).map_err(|e| {
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
    // ★★★ **Which surface the cap's plane became** — the boolean's own answer, not a guess made
    // afterwards. Its plane classes are decided with evidence (an exact `orient3d`, a composed
    // rotation proof, or a coincidence within the limit), and a later comparison of handles or
    // coordinates can see none of that: on a tilted face the cap merges with a face of the other
    // operand and the survivor carries *that* surface, which is exactly what `find_face_coplanar_with`
    // was left to guess at. `pad` used to throw away a correct solid when the guess missed.
    let want_surf = class_of
        .get(&far_cap)
        .copied()
        .unwrap_or_else(|| model.faces.get(far_cap).surface);
    match solids
        .iter()
        .find_map(|&s| find_face_coplanar_with(model, s, far_cap, want_surf, n).map(|c| (s, c)))
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
pub(crate) fn pad(
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
pub(crate) fn pocket(
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
pub(crate) fn find_face_coplanar_with(
    model: &Model,
    solid: Handle<Solid>,
    reference: Handle<Face>,
    ref_surf: Handle<Surface>,
    want: Vector3,
) -> Option<Handle<Face>> {
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
