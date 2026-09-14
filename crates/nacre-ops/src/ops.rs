//! Feature operations (design §7): the public sketch/extrude/pad/pocket API and the `apply`/
//! `replay` driver. The top layer — it composes the boolean engine ([`crate::boolean`]) and rigid
//! transform ([`crate::transform`]) over the plane substrate below.

use crate::BoolError;
use crate::boolean::boolean;
use crate::exact::{Seg3, Swept};
use crate::planes::outer_tri;
use crate::transform::transform;
use nacre_geom::Plane;
use nacre_geom::intersect::{RingSide, orient2d_rat, plane_side};
use nacre_geom::mixed::{
    Seg2d, mixed_ring_self_intersection, mixed_rings_cross, point_in_mixed_ring,
};
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::{Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::CylinderDef;
use nacre_topo::PointCache;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, MotionNode, Orientation, Shell, Solid, Surface, Vertex,
    VertexDef,
};
use std::borrow::Cow;

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

/// One closed ring of a [`Profile2d`], stored as its rational truth: its vertices and the step
/// leaving each — straight to the next vertex, or an arc around a stated centre ([`Seg2d`],
/// `nacre-geom`'s vocabulary, which is also what its predicates read). `segs[i]` runs
/// `vertices[i] → vertices[(i + 1) % n]`.
///
/// The coordinates are what the author's decimals *spelled* (`Rat::from_decimal`), not the f64s
/// that carried them — the same truth/cache split every dimension in the kernel gets
/// (`docs/truth-and-cache.md`). A computed coordinate (an arc's far end) is computed in `Rat` and
/// never rounds through f64.
///
/// **Normal form** ([`Ring2d::normalized`]). A vertex strictly mid-run on a straight edge is
/// dissolved — its two walls would be one plane, so the corner it names has no three-plane
/// definition, and deleting it changes no geometry. Two adjacent arcs of one circle in one
/// direction are one arc — their meeting vertex names no corner either: a pair of fillets that ate
/// their whole edge is a half circle, two half circles are the circle. One vertex with one arc is a
/// whole circle, and that vertex is its seam.
#[derive(Clone, Debug, PartialEq)]
pub struct Ring2d {
    vertices: Vec<[Rat; 2]>,
    segs: Vec<Seg2d>,
}

impl Ring2d {
    /// A polygon ring, in normal form.
    pub(crate) fn polygon(points: Vec<[Rat; 2]>) -> Ring2d {
        let n = points.len();
        Ring2d::normalized(points, vec![Seg2d::Line; n])
    }

    /// The normal form of `segs[i]: vertices[i] → vertices[i + 1]` — see the type doc. Runs to a
    /// fixpoint; each pass removes at most one vertex, and a ring of one or two vertices is left
    /// for [`Profile2d::check`] to judge (one vertex with an arc is a whole circle).
    pub(crate) fn normalized(mut vertices: Vec<[Rat; 2]>, mut segs: Vec<Seg2d>) -> Ring2d {
        debug_assert_eq!(vertices.len(), segs.len(), "one step leaves each vertex");
        loop {
            let n = vertices.len();
            if n < 2 {
                break;
            }
            let between = |a: [Rat; 2], b: [Rat; 2], p: [Rat; 2]| {
                p[0] >= a[0].min(b[0])
                    && p[0] <= a[0].max(b[0])
                    && p[1] >= a[1].min(b[1])
                    && p[1] <= a[1].max(b[1])
            };
            let mut dropped = false;
            // Later vertices first, the wrap-around pair last, so the ring keeps its first vertex
            // whenever it can (two half circles merge into the circle seamed at the first).
            for k in (1..n).chain(std::iter::once(0)) {
                let prev = (k + n - 1) % n;
                let next = (k + 1) % n;
                let flat = match (segs[prev], segs[k]) {
                    // Strictly mid-run: collinear and between its neighbours, and neither of
                    // them — a repeated point is a zero-length edge for `check` to name, not a
                    // corner to dissolve.
                    (Seg2d::Line, Seg2d::Line) => {
                        n >= 3
                            && vertices[k] != vertices[prev]
                            && vertices[k] != vertices[next]
                            && orient2d_rat(vertices[prev], vertices[k], vertices[next]) == 0
                            && between(vertices[prev], vertices[next], vertices[k])
                    }
                    (
                        Seg2d::Arc {
                            center: c1,
                            radius: r1,
                            ccw: w1,
                        },
                        Seg2d::Arc {
                            center: c2,
                            radius: r2,
                            ccw: w2,
                        },
                    ) => c1 == c2 && r1 == r2 && w1 == w2,
                    _ => false,
                };
                if flat {
                    // The step arriving at `k` runs on to `next`; the vertex and the step that
                    // left it go.
                    vertices.remove(k);
                    segs.remove(k);
                    dropped = true;
                    break;
                }
            }
            if !dropped {
                break;
            }
        }
        Ring2d { vertices, segs }
    }

    pub fn vertices(&self) -> &[[Rat; 2]] {
        &self.vertices
    }

    pub fn segs(&self) -> &[Seg2d] {
        &self.segs
    }

    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Straight steps only.
    pub fn is_polygon(&self) -> bool {
        self.segs.iter().all(|s| matches!(s, Seg2d::Line))
    }

    /// **Which way the ring runs**, exactly: `Positive` is counter-clockwise (`+x → +y`) — the
    /// signed-area question over straight steps and quarter-turn arcs, answered in integers by
    /// [`nacre_scalar::winding_sign_quarter_arcs`] on the profile's own coordinates (the author's
    /// decimals, not a lifted frame's widths). `None` for an arc that is not a quarter-turn
    /// multiple. Reachable by type, unreached by any producer today: `Edge2d::arc_rat` accepts any
    /// angle, but its one production caller (the kit's fillet) only rounds axis-aligned corners,
    /// so every arc a profile brings here is a quarter-turn multiple. A `None` reaches
    /// `prism_rings_in` and is refused as `PlaneWithoutExactForm` — a wall with no one at it.
    pub(crate) fn winding_sign(&self) -> Option<nacre_scalar::Orient> {
        let n = self.vertices.len();
        let (mut lines, mut arcs) = (Vec::new(), Vec::new());
        for i in 0..n {
            let (s0, e0) = (self.vertices[i], self.vertices[(i + 1) % n]);
            match self.segs[i] {
                Seg2d::Line => lines.push([s0, e0]),
                Seg2d::Arc {
                    center,
                    radius,
                    ccw,
                } => arcs.push(nacre_scalar::QuarterArc {
                    center,
                    radius,
                    start: s0,
                    end: e0,
                    ccw,
                }),
            }
        }
        nacre_scalar::winding_sign_quarter_arcs(&lines, &arcs)
    }
}

/// A closed planar region: one outer ring and any number of hole rings — straight steps and
/// circular arcs ([`Ring2d`]).
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
        Profile2d {
            outer: Ring2d::polygon(outer),
            holes: holes.into_iter().map(Ring2d::polygon).collect(),
        }
    }

    /// Rings already in normal form — `sketch`'s door once it has chained and sorted them.
    pub(crate) fn from_normalized_rings(outer: Ring2d, holes: Vec<Ring2d>) -> Profile2d {
        Profile2d { outer, holes }
    }

    /// Whether any ring carries an arc.
    pub fn has_arcs(&self) -> bool {
        !self.outer.is_polygon() || self.holes.iter().any(|h| !h.is_polygon())
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
        let undecidable = |_| OpError::ProfileUndecidable;
        for (id, r) in rings() {
            if r.is_empty() || (r.is_polygon() && r.len() < 3) {
                return Err(OpError::DegenerateProfile);
            }
            // Simplicity first: the parity below is only meaningful on a simple ring. The
            // predicate reports a zero-length edge as a pair with itself.
            if let Some((a, b)) = mixed_ring_self_intersection(r.mixed()).map_err(undecidable)? {
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
                if mixed_rings_cross(a.mixed(), b.mixed()).map_err(undecidable)? {
                    return Err(OpError::ProfileRingsMeet { a: *ida, b: *idb });
                }
            }
        }
        for (h, hole) in self.holes.iter().enumerate() {
            let probe = hole.vertices()[0];
            if point_in_mixed_ring(probe, self.outer.mixed()).map_err(undecidable)?
                != RingSide::Inside
            {
                return Err(OpError::HoleNotInsideOuter { hole: h });
            }
            for (k, other) in self.holes.iter().enumerate() {
                if k != h
                    && point_in_mixed_ring(probe, other.mixed()).map_err(undecidable)?
                        == RingSide::Inside
                {
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
/// ★ **`DatumPlane` is the large variant and it is not boxed.** It carries a [`SketchPlane`] in
/// its `Stated` arm, which holds a plane's exact rational definition (~400 bytes). An op log is
/// tens of entries long and is walked once per replay, so the wasted space is measured in
/// kilobytes; boxing would buy that back at the cost of an indirection on the one type a caller
/// constructs by hand. (`Extrude` used to be that variant; naming its plane instead of carrying
/// it made it small.)
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
pub enum Operation {
    /// **Put a plane in the model because the caller named it** — not because a face lies on it.
    ///
    /// Until this existed a log could name only two kinds of plane: the three seeded world planes
    /// and the plane of a face it had already built ([`face_sketch_frame`]). Everything else was
    /// stated *by value* inside [`Operation::Extrude`], which is why that variant is the last one
    /// carrying a `SketchPlane` at all.
    ///
    /// ★ **It is an operation, not a plain function, because `replay` must reproduce the handle.**
    /// A plane minted outside the log leaves a model that is not self-contained, and the
    /// "handles in a log are index vocabulary" contract (`docs/design.md` §2) is false for it.
    ///
    /// The plane may already exist — planes are the one thing interned at construction — in which
    /// case this returns the handle that exists and the arena does not grow. That is the intended
    /// answer, not a special case: *same plane, same handle*.
    DatumPlane { def: DatumDef },
    /// Extrude `profile`, drawn in `frame`, by `dist` along that frame's `ŵ`.
    ///
    /// ★ **The plane is named, not carried.** `frame` holds a `Handle<Surface>`, so the plane it
    /// sketches on is one the model already has — a seeded world plane
    /// ([`SketchFrame::world`]), the plane of a face ([`face_sketch_frame`]), or one a
    /// [`Operation::DatumPlane`] put there. That is what makes the base cap a *shared* handle
    /// rather than a second statement of the same plane.
    ///
    /// `dist` is a **thickness** and must be positive; which way it goes is the frame's, measured
    /// by whoever built the frame. (Contrast [`DatumDef::Offset`], whose `dist` is a signed
    /// displacement because there the sign is the only thing that says a side.)
    Extrude {
        frame: SketchFrame,
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
    /// Checked arithmetic overflowed while classifying the profile — `SketchError::Undecidable`'s
    /// twin. Refused by name; nothing guesses.
    ProfileUndecidable,
    /// An arc that is not a whole number of quarter turns. The builder's exact winding reads a
    /// ring's area as `a + b·π`, which needs every arc's angle to be a multiple of `π/2`; the
    /// vocabulary states such an arc exactly (a `3-4-5` lens is a valid region), the builder does
    /// not stand it yet — that family opens with named points.
    ArcSweepNotQuarterTurn,
    /// Two arcs of different circles meet at a profile vertex. That corner is a point on two
    /// cylinders and a plane, which no [`nacre_topo::VertexDef`] states yet, and the ruling between
    /// the two cylinder walls has no curve the kernel derives (`derive_edge_curve` declines two
    /// distinct cylinders) — the same frontier as the cylinder–cylinder boolean (M6b). A lens, a
    /// cam lobe; a straight step between the arcs is what builds today.
    ArcsMeetAtVertex,
    /// A non-positive extrusion distance.
    NonPositiveDistance,
    /// An extrusion distance outside the decimal window (`Rat::from_decimal` — `~1e38` above,
    /// `~1e-22` below for a full-width value): it has no rational truth for the sweep to be
    /// computed in. The sibling of [`OpError::ProfileOutsideDecimalWindow`], named at the
    /// operation's door — this used to fall silently to a point-less f64 prism.
    DistOutsideDecimalWindow,
    /// The sketch plane (or the frame chain carrying it) has no exact form to build in: its
    /// axes fell outside the decimal window or were degenerate, or the placement arithmetic
    /// overflowed `i128`. The prism this used to build silently in f64 recorded no exact
    /// points, could not survive a motion, and is the population `Inexact` grew from — a named
    /// reject is the honest answer (S6b).
    PlaneWithoutExactForm,
    /// The face's sketch frame **exists** — [`face_plane`] reports it and the pad sketches in
    /// it — but no [`SketchFrame`] realizes to it, so [`face_sketch_frame`] has nothing true
    /// to return. Not [`Self::PlaneWithoutExactForm`]: that says the exact form is missing,
    /// and here it is present — what is missing is a *spelling*.
    ///
    /// The population: world-branch faces whose surface carries a motion. The pad elides the
    /// frame node and sketches in world axes there, and every `SketchFrame` on that surface
    /// means "the pre-motion frame, then the motion" — stating the world axes in it would need
    /// the motion's inverse image, which is irrational. Defined by **verification failure**
    /// (no candidate's realization matches), not by that condition — which is how the
    /// invariant-plane restatement (2026-08-17) shrank it without a code change here: a plane
    /// its motion *fixes* keeps the world statement and verifies. What remains recorded, and
    /// so still lands here: exactly-statable-but-shifted images (a normal-wise translation
    /// after a turn), mirror chains, and moved sources re-moved (stage-1 boundaries).
    FrameNotRepresentable,
    /// A [`SketchFrame::named`] coordinate (origin or `ref_dir`) outside the decimal window
    /// (`Rat::from_decimal`), so the frame claim has no exact statement to check. The frame
    /// sibling of [`OpError::ProfileOutsideDecimalWindow`] and [`OpError::DistOutsideDecimalWindow`].
    FrameOutsideDecimalWindow,
    /// A [`SketchFrame::named`] origin whose exact residual against the plane's name is nonzero —
    /// the stated sketch origin does not lie on the stated plane. Rejected at construction:
    /// silently projecting it (or falling back to `Canonical`) would move the caller's sketch.
    OriginNotOnPlane,
    /// A [`SketchFrame::named`] `ref_dir` whose projection into the plane vanishes (parallel to
    /// the normal, or zero), so it names no `+u` direction.
    RefDirParallelToNormal,
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
    /// The log names a cell this model does not have. A log's handles are an **index
    /// vocabulary** (`docs/design.md` §2): [`replay`] re-anchors each one onto the model it is
    /// building, and an index past the end of the store means the log is not the one that built
    /// this arena — a hand-written handle, a truncated log, a log spliced from another session,
    /// or a session that kept recording after a *late* reject left arena cells the log does not
    /// account for (see [`replay`]'s note on self-containment).
    ///
    /// ★ **Existence, not legality.** Whether the cell is a *legal target* is still the
    /// operation's own question — [`OpError::SolidNotLive`] / [`OpError::FaceNotInLiveSolid`].
    /// And an index that is merely *wrong* rather than out of range cannot be caught here at
    /// all: it names a real cell, just not the intended one. That is why re-anchoring requires
    /// the self-containment premise rather than replacing it.
    /// A datum offset of zero. **It is a reject, not a no-op**: the plane it names is one the
    /// model already holds, but stated in a different coordinate system — and `push_plane` keys on
    /// `(name, motion)`, so accepting it would hand that plane a *second* handle. Planes are the
    /// one thing interned at construction precisely so that cannot happen, and flush contact is
    /// decided by comparing handles. Every nonzero offset of a tilted plane is irrational in the
    /// world and therefore has no rival statement; zero is the only case that collides.
    ZeroOffset,
    /// A datum named the same vertex twice — sorting reveals it, and two points do not fix a
    /// plane. Distinct from [`OpError::CollinearVertices`]: three different points on a line is a
    /// different mistake from two points and a typo.
    DuplicateVertex,
    /// A datum named three vertices that lie on one line — no unique plane through them.
    CollinearVertices,
    /// A datum named a vertex that is not a three-plane point (today: a cylinder seam vertex,
    /// whose pair pins a curve rather than a point). Its exact coordinate is not derivable from
    /// its definition, so the plane through it would not be exact either.
    VertexNotThreePlane,
    /// **The named vertices' carriers do not all share one motion**, so no single frame holds a
    /// rational coordinate for them: reaching the world means realizing a rotation, and cos/sin
    /// are irrational.
    ///
    /// ★ This is not a caller mistake. The dominant shape is not two solids combined by hand but a
    /// **single boolean result**: a cut between a turned operand and a still one leaves corners
    /// where an unmoved wall meets two turned ones, and measured, 12 of that solid's 20 vertices
    /// are in that state (`tests/point_width.rs`, `a_datum_on_straddling_carriers_has_no_name`).
    ///
    /// ★★★★ **Corrected 2026-08-08 — what opens this is not the judging layer.** This used to say
    /// the next stage's homogeneous lift opened the population. It does not, and the lift's
    /// machinery was built long before anything could reach it. Such
    /// a plane has **no exact name**, and from there:
    ///
    /// ```text
    /// no name → frame_chain declines → no SketchFrame → never a base cap → never in a plane table
    /// ```
    ///
    /// so it never gets as far as being judged. What would open it is a **frame for a nameless
    /// plane** — a realization that takes the plane's coefficients as intervals at a precision
    /// rather than as exact rationals or `BigInt`s, which is what `MoveNode::Frame` and
    /// `FrameWide` both require today. `docs/truth-and-cache.md`'s open item says so, and
    /// `a_plane_with_no_name_cannot_host_a_sketch` runs the chain.
    ///
    /// It is named separately so that "how much does this cost us today" stays countable.
    VerticesInMixedFrames,
    /// A mixed-frame datum's **judged frame could not be decided** at the fixed rung: no
    /// arbitrary-axis branch's squared length — or the normal's, or the origin's denominator —
    /// could be bounded away from zero (`nacre_cip::FrameThrough::of`).
    ///
    /// ★ Deliberately **not** [`OpError::CollinearVertices`] and not
    /// [`OpError::DegenerateGeometry`]: both claim the construction *is* degenerate, and an
    /// interval that fails to clear zero proves nothing of the kind — the points may be exactly
    /// collinear or merely too close to call. Failing to prove health is its own cause, so it
    /// gets its own name (C7).
    ThroughFrameUndecided,
    LogHandleOutOfRange {
        /// Which store the index was meant for.
        cell: LogCell,
        /// The index the log named.
        index: u32,
    },
}

/// Which store an operation-log handle indexes. (`Surface` joins when datum ops arrive — S5.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogCell {
    Face,
    Solid,
    Surface,
    Vertex,
}

/// **How a datum plane is stated** — the variant is the kind of statement, the way
/// `PlanePoints` splits truth into what is written and what is pointed at.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
pub enum DatumDef {
    /// **Stated in world coordinates.** Every [`SketchPlane`] constructor produces this, so the
    /// vocabulary a caller already has — `world_xy`, `from_origin_normal`, `through_points`,
    /// `from_axes`, `with_origin` — is the vocabulary of a datum, unchanged. It is also the whole
    /// of today's population: every plane any operation states is a rational point triple.
    ///
    /// [`OpError::PlaneWithoutExactForm`] when the plane carries no exact statement (a coordinate
    /// outside the decimal window, a zero normal) — the same proposition that error already names.
    Stated(SketchPlane),
    /// **The plane through three vertices the model already holds** — named, not measured.
    ///
    /// This is the one thing the stated vocabulary cannot say. A caller who wants "the plane
    /// through those three corners" can only read their coordinates today, and a discovered
    /// vertex's coordinate is rounded: measured on tilted geometry, **every one of 220 triples**
    /// produced a plane with a *different name* than the plane actually through them
    /// (`tests/point_width.rs`). Naming the vertices keeps the statement exact.
    ///
    /// ★ **The order is the direction.** The three are sorted before they are stored — the same
    /// vertices are the same plane in any order — but the caller's order fixes a normal by the
    /// right-hand rule, and `flip` is measured against it exactly as it is for a stated plane.
    /// Reversing two of them returns *the same handle with the opposite frame*, which is the only
    /// way a caller can choose a side (`dist` is positive-only).
    ///
    /// Rejects by cause rather than by one blanket failure, because the causes have different
    /// futures: [`OpError::VerticesInMixedFrames`] waits on a frame for a nameless plane (see
    /// there — this used to say "what the next stage opens", and that stage turned out not to be
    /// the one that opens it), and the rest are the caller's. (`VertexPointTooWide` retired with
    /// open item 17: a meet wider than `Rat` names its plane through `plane_name_from_meets`,
    /// so width stopped being a cause.)
    ThroughVertices([Handle<Vertex>; 3]),
    /// **`dist` away from a plane the model already holds**, stated inside that plane's own frame
    /// as the rational triple `(0,0,d), (1,0,d), (0,1,d)`.
    ///
    /// ★ **This is the exact form of an offset, and the frame is why.** In the world, "d along the
    /// normal" of a tilted plane is `p + d·n̂` — irrational, because `n̂` carries a square root.
    /// Inside the frame the same plane is `w = d` and every coordinate is a written decimal; the
    /// irrationality lives in the frame's realization, which is machinery the kernel already has.
    /// It is also what `docs/truth-and-cache.md` prescribes instead of floating a sketch origin off
    /// its plane.
    ///
    /// `frame` is taken rather than a bare handle so a caller can say **which side**: `+dist` runs
    /// along the frame's `ŵ`, which for a face's frame is its outward normal. The plane's
    /// *identity*, though, depends only on `(plane, signed distance)` — an origin and a `+u` do not
    /// move a parallel plane — so the operation normalizes to the plane's canonical frame and
    /// folds `flip` into the sign. Two callers holding different frames of one plane and meaning
    /// the same side get **one handle**.
    ///
    /// Rejects: [`OpError::ZeroOffset`] (`dist == 0` names the plane you already hold),
    /// [`OpError::DistOutsideDecimalWindow`], [`OpError::PlaneWithoutExactForm`] (the frame cannot
    /// be realized exactly, or its world form overflows).
    Offset { frame: SketchFrame, dist: f64 },
}

/// The handles an operation produced. Not `Copy`: `Extrude` carries a `Vec`.
///
/// ★ **`DatumPlane` is the large variant and it is not boxed**, for the same reason
/// [`Operation::Extrude`] is not: it carries a [`SketchFrame`], whose `Named` placement is six
/// `Rat`s. An `OpOutput` is a return value — one per `apply`, read and dropped, never accumulated
/// — so the width costs a stack copy, and boxing would buy that back by putting an indirection on
/// the value every caller destructures.
#[allow(clippy::large_enum_variant)]
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
    /// The datum plane, and **the frame the caller's statement implies**.
    ///
    /// ★ The frame is not a convenience. A caller who had to rebuild it with
    /// [`SketchFrame::named`] would hand back `f64`, taking `Rat::from_decimal` a second time —
    /// and a computed value can lift to a different rational than the one the statement already
    /// holds. Returning it keeps the caller's sketch origin and `+u` exact, which is what
    /// `extrude` does internally today.
    ///
    /// `flip` is `false`: which way `ŵ` must face is measured by the operation that consumes the
    /// frame, never stated here ([`SketchFrame`]'s constructor contract).
    DatumPlane {
        plane: Handle<Surface>,
        frame: SketchFrame,
    },
}

/// Apply one operation to `model`, returning the handles it created. Does not
/// rebuild the adjacency cache (do that once after a batch — see [`replay`]).
pub fn apply(model: &mut Model, op: &Operation) -> Result<OpOutput, OpError> {
    match op {
        Operation::DatumPlane { def } => {
            let (plane, frame) = datum_plane(model, def)?;
            Ok(OpOutput::DatumPlane { plane, frame })
        }
        Operation::Extrude {
            frame,
            profile,
            dist,
        } => {
            let (solid, faces) = extrude_on_frame(model, frame, profile, *dist)?;
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
///
/// ★★ **A log's handles are an index vocabulary, and this is where they are re-anchored.**
/// Six of [`Operation`]'s variants name a cell by `Handle`, and those handles belong to the
/// model the log was *recorded* against — a different arena from the one being built here. A
/// `Handle`'s identity is its index (`Store`'s manual `Eq`/`Hash` use nothing else), so the
/// index is the part that carries meaning across models, and [`nacre_store::Store::handle_at`]
/// turns it back
/// into a handle of *this* model, one operation at a time.
///
/// **Why one operation at a time, and not a pre-pass**: operation *N*'s handle names a cell
/// operation *N−1* created, so nothing can be checked before the walk. That also makes it
/// atomic for free — a re-anchoring failure happens before the operation pushes anything, and
/// the half-built `model` is local, so the caller never sees a partly-mutated arena.
///
/// ★ **The premise this rests on: the log is the whole history.** Cells put into a model
/// outside the log (`Model::add_cuboid`, a direct `push_*`) shift every later index, and so
/// does a *late* reject — one that pushed cells before declining (`PadMissesFace`,
/// `PocketNotBlind`, a boolean's reject) leaves them in the append-only arena while the log has
/// no entry for them. [`OpError::LogHandleOutOfRange`] catches only the case where the index
/// runs off the end; an index that lands on a real-but-wrong cell cannot be detected here. So a
/// session that keeps recording after a late reject must rebuild from its log first.
///
/// **`apply` does not re-anchor**, deliberately: its model belongs to the caller, so its
/// handles do too, and quietly re-anchoring there would launder a genuinely foreign handle and
/// destroy the cross-model guard that catches it.
pub fn replay(ops: &[Operation]) -> Result<Model, OpError> {
    let mut model = Model::new();
    for op in ops {
        let op = rebind(&model, op)?;
        apply(&mut model, &op)?;
    }
    model.rebuild_adjacency();
    Ok(model)
}

/// The operation with its handles re-anchored onto `model` — borrowed when there is nothing to
/// re-anchor, which since the plane-handle vocabulary landed is only a `Stated` datum.
fn rebind<'a>(model: &Model, op: &'a Operation) -> Result<Cow<'a, Operation>, OpError> {
    let surface = |h: Handle<Surface>| {
        model
            .surface_handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Surface,
                index: h.index(),
            })
    };
    let vertex = |h: Handle<Vertex>| {
        model
            .vertex_handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Vertex,
                index: h.index(),
            })
    };
    let face = |h: Handle<Face>| {
        model
            .faces
            .handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Face,
                index: h.index(),
            })
    };
    let solid = |h: Handle<Solid>| {
        model
            .solids
            .handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Solid,
                index: h.index(),
            })
    };
    Ok(match op {
        // ★ `Extrude` carries a handle now — it is no longer the one variant that does not.
        Operation::Extrude {
            frame,
            profile,
            dist,
        } => Cow::Owned(Operation::Extrude {
            frame: frame.rebound(surface(frame.plane())?),
            profile: profile.clone(),
            dist: *dist,
        }),
        Operation::PadOnFace {
            face: f,
            profile,
            dist,
        } => Cow::Owned(Operation::PadOnFace {
            face: face(*f)?,
            profile: profile.clone(),
            dist: *dist,
        }),
        Operation::PocketOnFace {
            face: f,
            profile,
            dist,
        } => Cow::Owned(Operation::PocketOnFace {
            face: face(*f)?,
            profile: profile.clone(),
            dist: *dist,
        }),
        Operation::Boolean { kind, a, b } => Cow::Owned(Operation::Boolean {
            kind: *kind,
            a: solid(*a)?,
            b: solid(*b)?,
        }),
        Operation::Transform { solid: s, isometry } => Cow::Owned(Operation::Transform {
            solid: solid(*s)?,
            isometry: *isometry,
        }),
        Operation::Mirror {
            solid: s,
            axis,
            offset,
        } => Cow::Owned(Operation::Mirror {
            solid: solid(*s)?,
            axis: *axis,
            offset: *offset,
        }),
        Operation::Copy { solid: s } => Cow::Owned(Operation::Copy { solid: solid(*s)? }),
        // A `Stated` datum names no cell — it is three rational points and nothing else.
        Operation::DatumPlane {
            def: DatumDef::Stated(_),
        } => Cow::Borrowed(op),
        Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        } => {
            let mut out = *vs;
            for v in out.iter_mut() {
                *v = vertex(*v)?;
            }
            Cow::Owned(Operation::DatumPlane {
                def: DatumDef::ThroughVertices(out),
            })
        }
        Operation::DatumPlane {
            def: DatumDef::Offset { frame, dist },
        } => Cow::Owned(Operation::DatumPlane {
            def: DatumDef::Offset {
                frame: frame.rebound(model.surface_handle_at(frame.plane.index()).ok_or(
                    OpError::LogHandleOutOfRange {
                        cell: LogCell::Surface,
                        index: frame.plane.index(),
                    },
                )?),
                dist: *dist,
            },
        }),
    })
}

/// **The definition of a corner where a straight wall meets an arc wall on a cap**: the wall's
/// plane and the cap's plane meet in a line, and that line crosses the arc's cylinder at this
/// point — [`VertexDef::Pierce`], the same definition the boolean mints for its pierce corners.
///
/// The root is read the way `QuadRoot` is defined: the two planes in **ascending handle order**,
/// each by its **stored canonical name** (`surface_name`, the sign convention that fixes the meet
/// line's direction), fed to [`nacre_scalar::quad::plane_plane_cylinder`] with the cylinder's own
/// statement; the root whose parameter equals this point's is the name. A wall tangent to the
/// cylinder — a fillet's, a slot's straight side — is the double root, one point. Every statement
/// has to live in one frame for the meet to mean anything: a cap borrowed from another body in
/// another frame (a pad on a turned face) declines by name rather than reading a frame line
/// against a world cylinder.
fn pierce_def(
    model: &Model,
    plane: Handle<Surface>,
    cap: Handle<Surface>,
    cylinder: Handle<Surface>,
    at: &[nacre_scalar::Rat; 3],
) -> Result<VertexDef, OpError> {
    use nacre_scalar::quad::CylinderMeet;
    use nacre_topo::QuadRoot;
    let (a, b) = if plane.index() < cap.index() {
        (plane, cap)
    } else {
        (cap, plane)
    };
    let name = |h: Handle<Surface>| -> Result<[nacre_scalar::Rat; 4], OpError> {
        model
            .surface_name
            .get(&h)
            .and_then(|n| n.narrow().copied())
            .ok_or(OpError::PlaneWithoutExactForm)
    };
    let (pa, pb) = (name(a)?, name(b)?);
    let nacre_topo::Surface::Cylinder { def, motion } = model.surface_truth(cylinder) else {
        return Err(OpError::DegenerateGeometry);
    };
    if model.plane_motion(a) != *motion || model.plane_motion(b) != *motion {
        return Err(OpError::PlaneWithoutExactForm);
    }
    let meet =
        nacre_scalar::quad::plane_plane_cylinder(&pa, &pb, &def.origin(), &def.dir(), def.radius())
            .ok_or(OpError::PlaneWithoutExactForm)?;
    // This point's parameter along the meet line: `(at − base)·dir / dir·dir`.
    let param = |line: &nacre_scalar::quad::MeetLine| -> Option<nacre_scalar::Rat> {
        let (bse, d) = (line.base(), line.dir());
        let mut num = nacre_scalar::Rat::from_int(0);
        let mut den = nacre_scalar::Rat::from_int(0);
        for k in 0..3 {
            num = num.checked_add(at[k].checked_sub(bse[k])?.checked_mul(d[k])?)?;
            den = den.checked_add(d[k].checked_mul(d[k])?)?;
        }
        num.checked_mul(nacre_scalar::Rat::new(den.denom(), den.numer())?)
    };
    let root = match meet {
        CylinderMeet::Tangent { line, s } => {
            debug_assert_eq!(param(&line), Some(s), "the tangent point is this corner");
            QuadRoot::Double
        }
        CylinderMeet::Pair { line, s } => {
            let t = param(&line).ok_or(OpError::PlaneWithoutExactForm)?;
            // Compared as values: a rational root still arrives as `mid ± k·√disc` when the
            // discriminant is a perfect square, so `b == 0` is not the test — the difference's
            // sign is.
            let is = |q: &nacre_scalar::quad::QuadVal| {
                q.checked_sub(&nacre_scalar::quad::QuadVal::from_rat(t))
                    .is_some_and(|d| d.sign() == nacre_scalar::Orient::Zero)
            };
            if is(&s[0]) {
                QuadRoot::Lo
            } else if is(&s[1]) {
                QuadRoot::Hi
            } else {
                return Err(OpError::PlaneWithoutExactForm);
            }
        }
        _ => return Err(OpError::DegenerateGeometry),
    };
    Ok(VertexDef::Pierce {
        planes: [a, b],
        cylinder,
        root,
    })
}

fn push_line_edge(
    model: &mut Model,
    a: Handle<Vertex>,
    b: Handle<Vertex>,
    carriers: [Handle<Surface>; 2],
) -> Result<Handle<Edge>, OpError> {
    model
        .push_edge(carriers, [a, b])
        .ok_or(OpError::DegenerateGeometry)
}

/// Put a stated plane in the arena and hand back its handle and the frame the statement implies.
///
/// ★★ **The cache faces `−normal`, and that is a convention with two independent reasons.**
/// (i) It is the sense a base cap gets: `extrude` pushes `−plane.normal()` and `Model::new` seeds
/// the world planes along `−axis`; S9 measured that seeding `+axis` instead flipped 781 stored
/// cap normals for nothing. (ii) `WorkingPlane::frame_sign` records whether a plane's *stored* normal
/// agrees with its root face's outward normal, and a base cap's outward is `−N` — so `−normal`
/// leaves that sign exactly where it is today. A datum that later becomes a base cap therefore
/// interns with `flipped == false` and nothing downstream has to compensate.
///
/// Which point anchors the cache is a choice the arena keeps (interning discards the newcomer's
/// cache), and `tests/plane_anchor.rs` measures what that choice costs: on the tilted population
/// where anchors can disagree at all, the judgment path does not read the disagreeing part and the
/// result holds to 1.1e-15 at the worst anchor.
fn datum_plane(
    model: &mut Model,
    def: &DatumDef,
) -> Result<(Handle<Surface>, SketchFrame), OpError> {
    /// The same solve `Model::through_points_rat` performs, with its single `None` split into the
    /// causes that have different owners.
    ///
    /// ★ The split is the point. `three_planes_rat`'s one `None` hid two different facts for as
    /// long as it existed, and the measurement that found it is the reason this variant is being
    /// added at all — repeating the shape here would make the next stage unable to read how much
    /// of the population it is opening.
    ///
    /// ★★ A carrier triple that does not meet is **not** on the list: those three planes met, or
    /// the vertex would not exist. That is asserted, not rejected, so a broken invariant cannot
    /// arrive disguised as a user error.
    ///
    /// ★★★★★ **The motion comes back with the points, because it *is* which frame they are
    /// written in.** Returning the points alone is what produced the defect this function was
    /// rewritten for: the caller then chose a motion, chose `None`, and a plane stated in a frame
    /// was filed as a world plane. `DatumDef::Offset` has always returned `(points, motion)` from
    /// one decision for exactly this reason.
    /// What the three vertices state, once the causes are told apart.
    ///
    /// (`Named` is 288 B beside a unit variant — the same shape and the same verdict as
    /// `PlanePoints`: it lives for one call on one stack frame, and boxing would buy nothing.)
    #[allow(clippy::large_enum_variant)]
    enum ThroughStatement {
        /// One shared frame — the named road (S5(ii)-1): the vertices' meets in that frame
        /// (**any width** since open item 17 — a meet wider than `Rat` still names its plane
        /// through `plane_name_from_meets`), and the frame.
        Named(
            [nacre_scalar::MeetPoint; 3],
            Option<Handle<nacre_topo::MotionNode>>,
        ),
        /// ★ No one frame holds all three: either a vertex the door cannot place at all (16-2's
        /// straddle) or three placeable vertices whose frames differ (16-1's first wall). No
        /// rational triple, no name, and the plane takes the judged road.
        Nameless,
    }

    #[allow(clippy::type_complexity)]
    fn through_points_by_cause(
        model: &Model,
        vs: [Handle<Vertex>; 3],
    ) -> Result<ThroughStatement, OpError> {
        // ★★★★ **Which frame a vertex is solvable in is [`nacre_topo::Model::vertex_meet`]'s
        // answer, not a second copy of its rule.** This function used to compare
        // `plane_motion(tri[0])` against the other two itself. That reading calls a turned
        // solid's corner a straddle — the invariant-plane restatement leaves its cap
        // world-stated while the walls carry a node — and, worse, it can now *disagree* with
        // the door: `push_plane_through` derives the interning name through `vertex_meet`, so a
        // producer that says "no frame" while the door says "this one" files a frame-local name
        // as a world plane. That is the defect `a_datum_through_frame_local_vertices_is_not_a
        // _world_plane` exists to catch. One decision, one place.
        let mut pts: [Option<nacre_scalar::MeetPoint>; 3] = [None, None, None];
        let mut frames = [None; 3];
        for (i, vh) in vs.iter().enumerate() {
            let tri = match model.vertices.get(*vh).def {
                nacre_topo::VertexDef::ThreePlane(tri) => tri,
                // A through-vertices datum needs three-plane meets; a seam vertex has no
                // point-meet at all and a pierce point has no rational one — the same honest
                // reject, spelled per variant.
                nacre_topo::VertexDef::OnSeam(_) | nacre_topo::VertexDef::Pierce { .. } => {
                    return Err(OpError::VertexNotThreePlane);
                }
            };
            if tri.iter().any(|h| !model.surface_name.contains_key(h)) {
                // ★★ The one population still refused here: a carrier that is itself a nameless
                // datum (a datum on a nameless datum — depth). Its triangle needs the judged
                // machinery recursively, which is open item 16-3's depth question. This is now
                // the **only** thing `VerticesInMixedFrames` names; measure before splitting
                // the label further.
                return Err(OpError::VerticesInMixedFrames);
            }
            let Some((meet, frame)) = model.vertex_meet(*vh) else {
                // ★ A straddling vertex — no rational coordinate anywhere, but a complete
                // definition (the meet of its carriers). The judged road takes it (16-2), as
                // long as every carrier can hand out a witness triangle of its own.
                //
                // ★★★ **"The carriers meet" is not assertable here, and finding that out cost a
                // red run.** The old pass two `expect`ed it, and rightly: it ran only after the
                // frames agreed. Solving three names *stated in different frames* is not the
                // same question — they are not three planes in one coordinate system, so
                // `three_planes_big` declines for a straddling vertex as a matter of course. An
                // assertion here would measure the proposition next door. The invariant now
                // lives where it is meaningful: past this `else`, every `pts[i]` is `Some`, so
                // the named road below has the meets by construction rather than by `expect`.
                return Ok(ThroughStatement::Nameless);
            };
            pts[i] = Some(meet);
            frames[i] = Some(frame);
        }
        // ★★ Pure vertices in differing frames — the judged road. The collinearity question is
        // *not* asked there: with no shared frame there is no exact solve to ask it in, and the
        // judged constructor's failure is reported as undecided, never as proven collinear.
        if !(frames[1] == frames[0] && frames[2] == frames[0]) {
            return Ok(ThroughStatement::Nameless);
        }
        // ★ The meets are kept at whatever width they need (open item 17) — the name is a
        // function of the *plane*, and `plane_name_from_meets` derives it without ever asking a
        // coordinate to fit `Rat`.
        let meets = pts.map(|p| p.expect("filled above"));
        if nacre_scalar::plane_name_from_meets([&meets[0], &meets[1], &meets[2]]).is_none() {
            return Err(OpError::CollinearVertices);
        }
        Ok(ThroughStatement::Named(meets, frames[0].flatten()))
    }

    match def {
        DatumDef::Stated(sp) => {
            let d = sp.def.as_ref().ok_or(OpError::PlaneWithoutExactForm)?;
            let cache = Plane::from_point_normal(sp.origin(), -sp.normal())
                .ok_or(OpError::DegenerateGeometry)?;
            let (plane, _flipped) = model.push_plane(cache, d.points(), None);
            // ★ `Named`, unconditionally — never derived. The canonical frame of the ZX plane has
            // `+u = −x̂` while the script convention (and `SketchPlane::world_zx`) says `+ẑ`, so a
            // placement inferred from the plane would silently turn some sketches. The values are
            // the `PlaneDef`'s own rationals, and `SketchFrame::named`'s checks hold structurally:
            // `origin = points[0]` is on the plane, and `ref_dir = points[1] − points[0]` lies in
            // it and is nonzero because the triple is not collinear.
            let placement = nacre_topo::FramePlacement::Named {
                origin: d.origin(),
                ref_dir: d.ref_dir(),
            };
            // ★★★ **`flip` is measured here, against the normal the caller stated** — and that is
            // the only place the caller's *direction* can survive.
            //
            // A plane's canonical name has no direction, and planes intern: state `z = 0` facing
            // `+ẑ` and state it facing `−ẑ`, and both come back as **one handle whose realized `ŵ`
            // is `+ẑ`** (measured). So a frame built with `flip: false` would silently answer "up"
            // to a caller who said "down". Measuring against `sp.normal()` — the direction their
            // own point order fixes — makes the returned frame mean what they said, which is what
            // lets an operation take a frame where it used to take a plane and sweep the same way.
            //
            // This is still S9's rule, not an exception to it: `flip` is *measured*, never stated,
            // and `measured_frame` is the one place that measures.
            let frame = measured_frame(model, plane, placement, sp.normal())
                .ok_or(OpError::PlaneWithoutExactForm)?;
            Ok((plane, frame))
        }
        DatumDef::ThroughVertices(vs) => {
            // ★★★ **Every reject here happens before anything is pushed.** The causes are told
            // apart first, then the plane is built — so a refusal leaves the model exactly as it
            // found it, and the name of the refusal says which stage owns it.
            let mut sorted = *vs;
            sorted.sort_by_key(|v| v.index());
            if sorted[0] == sorted[1] || sorted[1] == sorted[2] {
                return Err(OpError::DuplicateVertex);
            }
            let statement = through_points_by_cause(model, sorted)?;

            // ★ **The caller's order is the stated normal**, by the right-hand rule — the one
            // place their choice of side can survive, since the stored triple is sorted and a
            // canonical name carries no direction. `measured_frame` then measures `flip` against
            // it, exactly as the `Stated` arm does against `sp.normal()`.
            //
            // The vertices' world caches are the right input here even though they are rounded:
            // this asks only for a *direction*, and `measured_frame` compares it against the
            // plane's own world realization.
            let world = vs.map(|v| model.vertex_point(v));
            let stated = (world[1] - world[0])
                .cross(world[2] - world[0])
                .normalize()
                .ok_or(OpError::CollinearVertices)?;

            let ThroughStatement::Named(meets, motion) = statement else {
                // ★★★ **The judged road** (open item 16, first wall): every vertex pure, frames
                // differing — no name exists, and the frame is derived from the defining points
                // as intervals. **Validation comes before the push**: the judged constructor can
                // refuse (an undecidable basis), and rejecting after `push_plane_through` would
                // leave a frameless nameless plane in the arena — the reject-after-commit shape
                // this arm's own header forbids. The shared derivation
                // (`through_judged_points`) is the one `frame_chain` will re-run, so what was
                // validated is what gets framed.
                let jpts = crate::rotated_vertex::through_judged_points(model, sorted)
                    .ok_or(OpError::PlaneWithoutExactForm)?; // unreachable: causes told apart above
                let ft = nacre_cip::FrameThrough::of(jpts, false)
                    .ok_or(OpError::ThroughFrameUndecided)?;
                // The cache anchors at the validated realization of the first stored vertex —
                // the definition's own replay, same rule as the named road below.
                let anchor =
                    Point3::from_array(ft.anchor_coord().ok_or(OpError::ThroughFrameUndecided)?);
                let cache =
                    Plane::from_point_normal(anchor, -stated).ok_or(OpError::DegenerateGeometry)?;
                let (plane, _flipped) = model.push_plane_through(cache, sorted, None);
                let frame =
                    measured_frame(model, plane, nacre_topo::FramePlacement::Canonical, stated)
                        .ok_or(OpError::PlaneWithoutExactForm)?;
                return Ok((plane, frame));
            };

            // ★★★ **The cache is a world description, and the meets are not world coordinates
            // unless the carriers share no motion.** Realizing them means walking the very chain
            // the definition names — `exact.rs` states the rule ("the realization must be the
            // definition's own replay, not a second route to the same real number"), and taking
            // any other road here is how a plane's cache and its truth end up describing
            // different planes.
            //
            // ★ **Split by width so the narrow bits stay put**: a `Narrow` meet takes the road
            // this arm always took, letter for letter. A `Wide` meet (open item 17) has no
            // `Rat` triple to replay — and no f64 realization of it survives an exact lift
            // either (a value like 5⁻⁴⁰ rounds to a dyadic whose denominator leaves `i128`) —
            // so its cache anchors at the **first stored vertex's own world cache**: the
            // definition's replay already performed by whoever made the vertex (the 8/8 lock is
            // what says re-solving and replaying lands on it), and a rounded cache like every
            // anchor (`plane_anchor.rs` measured what anchor wobble costs — nothing the judging
            // reads).
            let anchor = match (&meets[0], motion) {
                (nacre_scalar::MeetPoint::Narrow(p), None) => {
                    Point3::from_array(p.map(|r| r.to_f64()))
                }
                (nacre_scalar::MeetPoint::Narrow(p), Some(leaf)) => {
                    let chain = crate::rotated_vertex::motion_chain(model, leaf)
                        .ok_or(OpError::PlaneWithoutExactForm)?;
                    let w = crate::rotated_vertex::replay(nacre_cip::WitnessPoint::at(*p), &chain)
                        .ok_or(OpError::PlaneWithoutExactForm)?;
                    Point3::from_array(w.coord())
                }
                (nacre_scalar::MeetPoint::Wide(_), _) => model.vertex_point(sorted[0]),
            };
            let cache =
                Plane::from_point_normal(anchor, -stated).ok_or(OpError::DegenerateGeometry)?;
            let (plane, _flipped) = model.push_plane_through(cache, sorted, motion);
            let frame = measured_frame(model, plane, nacre_topo::FramePlacement::Canonical, stated)
                .ok_or(OpError::PlaneWithoutExactForm)?;
            Ok((plane, frame))
        }
        DatumDef::Offset { frame, dist } => {
            // ★ Zero is the one offset that would duplicate a plane the model already holds; see
            // `OpError::ZeroOffset`. Checked before the lift so the reject names the real fault.
            if *dist == 0.0 {
                return Err(OpError::ZeroOffset);
            }
            let d =
                nacre_scalar::Rat::from_decimal(*dist).ok_or(OpError::DistOutsideDecimalWindow)?;
            let base = frame.plane;

            // ★★ **Normalize to the plane's canonical frame, folding `flip` into the sign.** A
            // parallel plane is fixed by `(plane, signed distance)` — a placement's origin and
            // `+u` do not move it — so building in the caller's frame would give one geometric
            // plane as many handles as there are frames naming it. The caller's `ŵ` is compared
            // against the canonical one to recover which side they meant; the dot is between two
            // realizations of the same unit normal, so it is a full magnitude from zero.
            let canonical = nacre_topo::FramePlacement::Canonical;
            let cb = crate::rotated_vertex::frame_world_basis(model, base, &canonical, false)
                .ok_or(OpError::PlaneWithoutExactForm)?;
            let sb = crate::rotated_vertex::frame_world_basis(
                model,
                base,
                frame.placement(),
                frame.flip(),
            )
            .ok_or(OpError::PlaneWithoutExactForm)?;
            let same_side: f64 = (0..3).map(|k| cb.3[k] * sb.3[k]).sum();
            let d = if same_side < 0.0 {
                nacre_scalar::Rat::from_int(0)
                    .checked_sub(d)
                    .ok_or(OpError::DistOutsideDecimalWindow)?
            } else {
                d
            };

            // ★★ **Say it in the world when the world can hold it** — the same node-omission
            // normalization S9 froze for frames. A plane stated under a frame node lives at the
            // key `(name, Some(node))`, so an offset of the world XY plane would *not* intern with
            // a box's cap on the same plane. Where the canonical basis lifts to exact rational
            // orthonormal axes, the offset plane is rational in the world and is stated there.
            // (The gate expression is `exact()`, the same one the extrude road uses — asked here
            // of the *realized* basis. The existing call site is untouched, so no node population
            // moves.)
            let realized = realized_plane(
                Point3::from_array(cb.0),
                Vector3::from_array(cb.1),
                Vector3::from_array(cb.2),
            );
            let (points, motion) = match realized.exact() {
                // ★ An overflowing pullback is a **named reject**, not a quiet switch to the frame
                // road — that switch is exactly where the duplicate handle would appear.
                Some(rf) => (
                    rf.offset_plane_points(d)
                        .ok_or(OpError::PlaneWithoutExactForm)?,
                    None,
                ),
                None => {
                    let zero = nacre_scalar::Rat::from_int(0);
                    let one = nacre_scalar::Rat::from_int(1);
                    let node = push_frame_node(
                        model,
                        SketchFrame {
                            plane: base,
                            placement: canonical,
                            flip: false,
                        },
                    );
                    (
                        [[zero, zero, d], [one, zero, d], [zero, one, d]],
                        Some(node),
                    )
                }
            };

            // The cache rides the canonical realization: anchor `o + d·ŵ`, facing `−ŵ` — the same
            // `−normal` convention `Stated` and every base cap keep.
            //
            // ★ **`d` here is the *signed* distance, not the caller's `dist`.** Using the raw one
            // puts the cache on the far side of the plane its own truth names whenever the
            // caller's frame is flipped — an incoherent surface, and one that a `flip` positive
            // control caught rather than any amount of reading.
            let signed = d.to_f64();
            let anchor = Point3::from_array(core::array::from_fn(|k| cb.0[k] + signed * cb.3[k]));
            let w = Vector3::from_array(cb.3);
            let cache = Plane::from_point_normal(anchor, -w).ok_or(OpError::DegenerateGeometry)?;
            let (plane, _flipped) = model.push_plane(cache, points, motion);
            Ok((plane, SketchFrame::canonical(plane)))
        }
    }
}

/// **The frame's basis in exact rationals, when it has one** — the gate that decides whether a
/// sketch is written in world coordinates or inside a motion node.
///
/// Only a frame on a plane with **no motion of its own** can be rational in the world: a moved
/// plane's axes carry the motion's irrational part, which is precisely why that road exists. So a
/// chain of one narrow frame node is the whole population, and anything longer (or wide) declines.
///
/// ★ Declining is not a failure — it is the frame-node road, the same one a tilted face takes.
fn exact_frame(model: &Model, frame: &SketchFrame) -> Option<crate::exact::RatFrame> {
    let chain =
        crate::rotated_vertex::frame_chain(model, frame.plane(), frame.placement(), frame.flip())?;
    match chain.as_slice() {
        [nacre_cip::MoveNode::Frame { frame: pf }] => crate::exact::RatFrame::of_plane_frame(pf),
        _ => None,
    }
}

/// **Extrude a profile on a frame the model already holds** — the handle vocabulary of
/// [`Operation::Extrude`], and the same three steps `extrude_and_boolean` takes for a pad.
///
/// The base cap **is** the frame's plane, so it is handed to [`build_prism`] as a surface rather
/// than as points: the flush contact then reads as one shared handle, which is what the boolean
/// recognizes. Nothing new is pushed for it.
///
/// `dist > 0` is a **thickness**; which way it goes is the frame's `ŵ`, measured by whoever built
/// the frame (a datum against the caller's stated normal, a face against its outward). That is why
/// this can take a frame where the operation used to take a plane and sweep the same way.
/// The builder's exact winding needs quarter-turn arcs (`Ring2d::winding_sign`); an arc of any
/// other angle is refused by name here, before the builder reads the ring.
fn refuse_non_quarter_arcs(profile: &Profile2d) -> Result<(), OpError> {
    let zero = Rat::from_int(0);
    let quarter = |r: &Ring2d| -> Result<(), OpError> {
        let n = r.len();
        for i in 0..n {
            let Seg2d::Arc { center, .. } = r.segs()[i] else {
                continue;
            };
            let (s0, e0) = (r.vertices()[i], r.vertices()[(i + 1) % n]);
            if s0 == e0 {
                continue; // a whole circle
            }
            let v = |p: [Rat; 2]| -> Option<[Rat; 2]> {
                Some([p[0].checked_sub(center[0])?, p[1].checked_sub(center[1])?])
            };
            let (a, b) = (
                v(s0).ok_or(OpError::ProfileUndecidable)?,
                v(e0).ok_or(OpError::ProfileUndecidable)?,
            );
            let dot = a[0]
                .checked_mul(b[0])
                .and_then(|x| x.checked_add(a[1].checked_mul(b[1])?));
            let cross = a[0]
                .checked_mul(b[1])
                .and_then(|x| x.checked_sub(a[1].checked_mul(b[0])?));
            let (Some(dot), Some(cross)) = (dot, cross) else {
                return Err(OpError::ProfileUndecidable);
            };
            // A quarter, a half, three quarters: perpendicular radii, or opposite ones.
            if !((dot == zero && cross != zero) || (dot < zero && cross == zero)) {
                return Err(OpError::ArcSweepNotQuarterTurn);
            }
        }
        Ok(())
    };
    quarter(profile.outer())?;
    profile.holes().iter().try_for_each(quarter)
}

pub(crate) fn extrude_on_frame(
    model: &mut Model,
    frame: &SketchFrame,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    if nacre_scalar::Rat::from_decimal(dist).is_none() {
        return Err(OpError::DistOutsideDecimalWindow);
    }
    profile.check()?;
    refuse_non_quarter_arcs(profile)?;
    // ★ A `SketchFrame` may name any surface — `SketchFrame::canonical` makes no claim and checks
    // nothing — so a cylinder can reach here, which a `SketchPlane` never could. Reject it by name
    // rather than letting the frame derivation fail later for a reason that reads as something
    // else ("no exact form" when the truth is "not a plane").
    match model.surface(frame.plane()) {
        nacre_geom::Surface::Plane(_) => {}
        nacre_geom::Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    }
    let (_, _, _, w) = crate::rotated_vertex::frame_world_basis(
        model,
        frame.plane(),
        frame.placement(),
        frame.flip(),
    )
    .ok_or(OpError::PlaneWithoutExactForm)?;

    // World road when the frame's basis is rational, frame-node road otherwise — the same
    // question `SketchPlane::exact` asks a caller's axes, asked of the frame in `Rat`.
    let (rat, node) = match exact_frame(model, frame) {
        Some(f) => (f, None),
        None => (
            crate::exact::RatFrame::identity(),
            Some(push_frame_node(model, *frame)),
        ),
    };
    let (outer, holes) = crate::exact::prism_rings_in(model, rat, profile, dist, node)
        .ok_or(OpError::PlaneWithoutExactForm)?;
    build_prism(
        model,
        outer,
        holes,
        Vector3::from_array(w),
        Some(frame.plane()),
        None,
    )
}

/// Place a profile on its plane and sweep it — **exactly, or not at all** (S6b).
///
/// The rational path is not an optimization: it is what makes `extrude(7.7)` and
/// `extrude(1.1)` then `extrude(6.6)` put their caps on the same plane rather than an
/// ulp apart. Where it does not apply — a plane with no exact statement (axes outside the
/// decimal window, a degenerate pair), a frame the chain cannot realize, an i128 overflow in
/// the placement arithmetic — the answer is a **named reject**, not the silent f64 prism this
/// used to build: a point-less solid cannot state itself, cannot survive a motion, and is the
/// population `Inexact` grew from. (Profile coordinates and `dist` outside the decimal window
/// are named before this runs.)
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
) -> Result<(Swept, Vec<Swept>), OpError> {
    crate::exact::prism_rings(model, plane, profile, dist, frame)
        .ok_or(OpError::PlaneWithoutExactForm)
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
/// No vertex carries a measured tolerance. Returns the solid and its faces: `faces[0]` = base cap
/// (at the ring, normal `−ŝ`), `faces[1]` = far cap, then the outer walls, then each hole's walls.
/// Shared by [`extrude_on_frame`] (a boss) and the pocket (`sweep = −n`).
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
    // A polygon needs three corners; a ring with an arc bounds area with one (a circle) or two
    // (a half disk, a slot's end drawn alone).
    let thin = |r: &Swept| r.base.len() < 3 && !r.exact.segs.iter().any(Seg3::is_arc);
    if thin(&outer_ring) || inner_rings.iter().any(thin) {
        return Err(OpError::DegenerateProfile);
    }

    // Whether the sweep runs along the frame normal (a boss) or against it (a pocket): the rings
    // are normalized about the *sweep*, the arcs' `ccw` is stated about the *normal*, and the
    // lateral orientation rule below compares the two.
    let sweep_up = outer_ring.normal.dot(normal) > 0.0;
    let outer_pts = oriented_ring(outer_ring, sweep_up, true)?;
    let hole_pts: Vec<Swept> = inner_rings
        .into_iter()
        .map(|h| oriented_ring(h, sweep_up, false))
        .collect::<Result<_, _>>()?;

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
                nacre_geom::Surface::Plane(p) => p.normal(),
                nacre_geom::Surface::Cylinder(_) => return Err(OpError::DegenerateGeometry),
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
            let (base_motion, cap_pts) = match base_cap_points {
                Some(p) => (None, p),
                None => {
                    let e = &outer_pts.exact;
                    (
                        e.motion,
                        e.cap_points(false).ok_or(OpError::DegenerateGeometry)?,
                    )
                }
            };
            let (s, flipped) = model.push_plane(
                Plane::from_point_normal(outer_pts.base[0], -normal)
                    .ok_or(OpError::DegenerateGeometry)?,
                cap_pts,
                base_motion,
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
    let top_motion = outer_pts.exact.motion;
    let top_points = outer_pts
        .exact
        .cap_points(true)
        .ok_or(OpError::DegenerateGeometry)?;
    let (top_surface, top_flipped) = model.push_plane(
        Plane::from_point_normal(outer_pts.top[0], normal).ok_or(OpError::DegenerateGeometry)?,
        top_points,
        top_motion,
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
    for (ring, walls, swept_pts) in std::iter::once((&outer, &outer_walls, &outer_pts)).chain(
        holes
            .iter()
            .zip(hole_walls.iter())
            .zip(hole_pts.iter())
            .map(|((r, w), p)| (r, w, p)),
    ) {
        ring.push_walls(model, walls, &swept_pts.exact.segs, sweep_up, &mut faces);
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
    /// Whether step `i`'s edge, walked in its own direction, follows the ring. A straight edge is
    /// pushed in ring order; a circle edge's own direction is fixed by the kernel's convention
    /// (`[A, B]` counter-clockwise about the axis — for a whole circle `A == B`, so the vertex
    /// order cannot carry it), so a clockwise arc's edge runs *against* the ring and every loop
    /// that walks it flips `forward`.
    along: Vec<bool>,
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
                        forward: !self.along[i],
                    })
                    .collect(),
            },
            Cap::Top => Loop {
                half_edges: (0..n)
                    .map(|i| HalfEdge {
                        edge: self.te[i],
                        forward: self.along[i],
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
        segs: &[Seg3],
        sweep_up: bool,
        faces: &mut Vec<Handle<Face>>,
    ) {
        let n = self.len();
        for (i, &(surface, flipped)) in walls.iter().enumerate().take(n) {
            // A plane wall's sense is the plane's (`flipped` from `push_plane`). A cylinder wall
            // is `Forward` iff the material lies **inside** the cylinder. Seen from `+ŵ`, the
            // material is on the ring's left iff the sweep runs along `ŵ` (the rings are
            // normalized about the sweep — outer and holes alike, since a hole runs the other way
            // and the material is outside it), and the arc's centre is on its left iff the arc is
            // counter-clockwise about `ŵ`. Material inside ⟺ the two sides agree: a boss's convex
            // circle is `Forward`, a bore's circle and a notch's concave arc `Reversed`.
            let flipped = match &segs[i] {
                Seg3::Line => flipped,
                Seg3::Arc { ccw, .. } => *ccw != sweep_up,
            };
            let j = (i + 1) % n;
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: self.be[i],
                        forward: self.along[i],
                    },
                    HalfEdge {
                        edge: self.ve[j],
                        forward: true,
                    },
                    HalfEdge {
                        edge: self.te[i],
                        forward: !self.along[i],
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
            match &ring.exact.segs[i] {
                Seg3::Line => Ok(model.push_plane(
                    Plane::through_points(ring.base[i], ring.base[j], ring.top[i])
                        .ok_or(OpError::DegenerateGeometry)?,
                    ring.exact.wall_points(i),
                    ring.exact.motion,
                )),
                // The wall is the cylinder about the arc's centre along the frame normal, stated
                // exactly and interned by that statement: every arc of one circle in one sketch
                // lands on one surface. The `bool` is a plane's flip; a cylinder's sense is
                // decided by `push_walls` from the arc's turn.
                Seg3::Arc {
                    center,
                    radius,
                    ref_dir,
                    cache,
                    ..
                } => {
                    let def = CylinderDef::new(*center, ring.exact.normal, *ref_dir, *radius)
                        .ok_or(OpError::DegenerateGeometry)?;
                    Ok((model.push_cylinder(*cache, def, ring.exact.motion), false))
                }
            }
        })
        .collect()
}

/// `pts` wound counter-clockwise about `normal` when `ccw`, clockwise when not. The test is the
/// polygon's area vector against `normal`, so it does not care which axis dominates.
/// **Requires a nonzero area.** The winding is read from the sign of the area vector, and a ring
/// that encloses nothing gives zero — the comparison below would then pick a side by accident.
/// [`Profile2d::check`] is what guarantees it: a simple polygon cannot have zero area, and a ring
/// that folds back on itself (a symmetric bowtie cancels to exactly zero) is not simple.
/// Normalize a ring's direction about the **sweep**: `ccw` for the outer ring, its opposite for a
/// hole. The winding is read exactly ([`crate::exact::SweptRat::winding`], about the ring's world
/// normal — the frame's motion, a reflection included, already folded in) and turned to the
/// sweep's sense by `sweep_up`; a ring with arcs has no f64 polygon area to read, and a polygon's
/// reads the same as before (`debug_assert`ed — the check that caught the mirrored pad).
fn oriented_ring(ring: Swept, sweep_up: bool, ccw: bool) -> Result<Swept, OpError> {
    let about_normal = ring.exact.winding;
    if about_normal == nacre_scalar::Orient::Zero {
        return Err(OpError::DegenerateProfile);
    }
    let ccw_about_normal = about_normal == nacre_scalar::Orient::Positive;
    debug_assert!(
        ring.exact.segs.iter().any(Seg3::is_arc) || {
            let v = &ring.base;
            let k = v.len();
            let area_vec = (0..k)
                .map(|i| (v[i] - Point3::origin()).cross(v[(i + 1) % k] - Point3::origin()))
                .fold(Vector3::from_array([0.0; 3]), |a, b| a + b);
            (area_vec.dot(ring.normal) > 0.0) == ccw_about_normal
        },
        "the exact winding agrees with the f64 polygon area"
    );
    let ccw_about_sweep = ccw_about_normal == sweep_up;
    Ok(if ccw_about_sweep != ccw {
        ring.reversed()
    } else {
        ring
    })
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
    // ★ Total since S7: the guard's `None` (adjacent same-plane walls) is unreachable — S3's
    // profile constructor dissolves collinear midpoints — and a wall can never equal a cap (a
    // wall contains the sweep direction, a cap has it as normal). The honest reject stands in
    // for the unreachable arm; the suite is what would refute "unreachable".
    //
    // (The frame base vertex died here (S7, Q2): a frame-drawn corner is the intersection of
    // three planes sharing the frame's motion, and solving them in that frame and replaying
    // the chain reproduces the stored coordinate bit for bit — measured 8/8, re-measured
    // suite-wide by the reuse differential. The shell-less base-vertex scaffolding went
    // with it.)
    let define = |model: &Model,
                  i: usize,
                  cap: Handle<Surface>,
                  at: &[nacre_scalar::Rat; 3]|
     -> Result<VertexDef, OpError> {
        let prev = walls[(i + n - 1) % n].0;
        let here = walls[i].0;
        if here == cap || prev == cap {
            return Err(OpError::DegenerateGeometry);
        }
        let (prev_arc, here_arc) = (
            ring.exact.segs[(i + n - 1) % n].is_arc(),
            ring.exact.segs[i].is_arc(),
        );
        match (prev_arc, here_arc) {
            // Two straight walls and the cap: the corner every polygon prism has.
            (false, false) => {
                if prev == here {
                    return Err(OpError::DegenerateGeometry);
                }
                Ok(VertexDef::ThreePlane([prev, here, cap]))
            }
            // A whole circle: one wall, one vertex — the rim's point at `+ref_dir`, the seam.
            (true, true) if n == 1 => Ok(VertexDef::OnSeam([here, cap])),
            // Two arcs of one circle merged in the profile's normal form; two of different
            // circles are a corner the kernel does not define.
            (true, true) => Err(OpError::ArcsMeetAtVertex),
            // A straight wall meets an arc wall on the cap: the wall's plane and the cap's plane
            // meet in a line that crosses the arc's cylinder there.
            (true, false) => pierce_def(model, here, cap, prev, at),
            (false, true) => pierce_def(model, prev, cap, here, at),
        }
    };
    let push_verts = |model: &mut Model,
                      ps: &[Point3],
                      exact: &[[nacre_scalar::Rat; 3]],
                      cap: Handle<Surface>|
     -> Result<Vec<Handle<Vertex>>, OpError> {
        ps.iter()
            .zip(exact.iter())
            .enumerate()
            .map(|(i, (p, at))| {
                let def = define(model, i, cap, at)?;
                Ok(model.push_vertex(def, PointCache::Unmeasured(*p)))
            })
            .collect()
    };
    let bv = push_verts(model, &base_pts, &ring.exact.base, caps.0)?;
    let tv = push_verts(model, &top_pts, &ring.exact.top, caps.1)?;

    let along: Vec<bool> = ring
        .exact
        .segs
        .iter()
        .map(|s| !matches!(s, Seg3::Arc { ccw: false, .. }))
        .collect();
    let (mut be, mut te, mut ve) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        let j = (i + 1) % n;
        // Carriers: the same expression `define` uses for the corner triples — a base/top edge
        // runs between wall `i` and its cap, a riser between the walls either side of corner `i`.
        // A clockwise arc's edge is stored the other way round so its own direction — the
        // kernel's `[A, B]` counter-clockwise — names the same points; the loops then walk it
        // backwards (`along`).
        let (p, q) = if along[i] { (i, j) } else { (j, i) };
        be.push(push_line_edge(model, bv[p], bv[q], [walls[i].0, caps.0])?);
        te.push(push_line_edge(model, tv[p], tv[q], [walls[i].0, caps.1])?);
        ve.push(push_line_edge(
            model,
            bv[i],
            tv[i],
            [walls[(i + n - 1) % n].0, walls[i].0],
        )?);
    }
    Ok(RingCells {
        base_pts,
        be,
        te,
        ve,
        along,
    })
}

/// A planar face's live solid, its in-plane right-handed frame (`x × y = n`, centred on the face
/// centroid so a profile's `(0,0)` lands there), and its loops — the shared setup for placing a
/// profile on a face (pad / pocket).
/// **Which frame a sketch lives in** — the plane (a handle: one statement of the plane, shared
/// with every face on it), its [`nacre_topo::FramePlacement`], and whether the plane's canonical
/// coefficients need negating to face the way the sketch does. Everything a
/// [`nacre_topo::Motion::Frame`] node needs, before the model has one (S9).
///
/// ★★ **The fields are private and the constructors validate** — the reason this type is not a
/// plain record. A `Named` placement is a *claim*: "this origin lies on that plane, this
/// direction crosses its normal". [`SketchFrame::named`] checks the claim exactly, at
/// construction, and rejects by name — silently substituting `Canonical` would move a caller's
/// sketch and answer a question they did not ask (`docs/truth-and-cache.md`'s rule). A public
/// field would let a literal walk around the check.
///
/// ★ **`flip` is not the caller's to state** — the canonical coefficients carry no direction, so
/// which way `ŵ` must face is a fact about the *use* (a sweep's sense, a face's outward normal),
/// measured by the consuming operation against the realized basis. Constructors set `false`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SketchFrame {
    plane: Handle<Surface>,
    placement: nacre_topo::FramePlacement,
    flip: bool,
}

impl SketchFrame {
    /// The frame nobody named: origin at the world origin's projection, axes by the
    /// arbitrary-axis convention — derived from the plane when the chain is flattened, stored
    /// nowhere. No claim is made, so there is nothing to validate.
    pub fn canonical(plane: Handle<Surface>) -> SketchFrame {
        SketchFrame {
            plane,
            placement: nacre_topo::FramePlacement::Canonical,
            flip: false,
        }
    }

    /// **One of the three world planes, framed the way the convention names them** — the frame
    /// twin of [`SketchPlane::world_xy`] / [`world_yz`](SketchPlane::world_yz) /
    /// [`world_zx`](SketchPlane::world_zx), and the door most sketches come through.
    ///
    /// ★★ **XY and YZ take `Canonical`; ZX must take `Named`.** The arbitrary-axis rule gives the
    /// ZX plane `+u = −x̂`, while the convention — and `SketchPlane::world_zx` — says `+u = +ẑ`.
    /// Derived and stated differ there and only there, so a caller who reached for
    /// [`SketchFrame::canonical`] on the ZX seed would find their sketch a quarter turn from where
    /// they asked. This function is where that one exception lives, spelled once.
    ///
    /// ★ **`flip` is `false` and that is a measurement, not a derivation**: each seed's canonical
    /// `ŵ` realizes to `+axis` because of how `Model::new` writes the seeds' points and how the
    /// canonical name normalizes — pinned by `a_world_frame_faces_its_axis`, and by S9's
    /// `a_seeded_planes_canonical_frame_is_the_world_basis_exactly` underneath it. Re-seed the
    /// world planes differently and this turns silently; the tests are what stop that.
    ///
    /// To sketch facing the *other* way, state a plane facing that way
    /// ([`Operation::DatumPlane`] measures `flip` from the normal you state) — the same move as
    /// writing `SketchPlane::from_origin_normal(o, -ẑ)` today.
    pub fn world(model: &Model, axis: nacre_scalar::Axis) -> SketchFrame {
        let plane = model.world_plane(axis);
        match axis {
            // The two whose derived frame already is the convention.
            nacre_scalar::Axis::Z | nacre_scalar::Axis::X => SketchFrame::canonical(plane),
            // ZX: `+u = +ẑ`, stated because it cannot be derived.
            nacre_scalar::Axis::Y => SketchFrame {
                plane,
                placement: nacre_topo::FramePlacement::Named {
                    origin: [nacre_scalar::Rat::from_int(0); 3],
                    ref_dir: [
                        nacre_scalar::Rat::from_int(0),
                        nacre_scalar::Rat::from_int(0),
                        nacre_scalar::Rat::from_int(1),
                    ],
                },
                flip: false,
            },
        }
    }

    /// A frame whose origin and `+u` direction the caller states — validated **now**, against
    /// the plane's exact name, so a bad claim is a named reject at the door rather than a sketch
    /// somewhere else:
    ///
    /// * [`OpError::FrameOutsideDecimalWindow`] — a coordinate of `origin`/`ref_dir` has no
    ///   decimal truth (`Rat::from_decimal` — `~1e38` above, `~1e-22` below for a full-width
    ///   value); the claim cannot even be stated exactly.
    /// * [`OpError::PlaneWithoutExactForm`] — the plane carries no name: the same proposition as
    ///   everywhere this error fires, "there is no exact statement to check against". ★ Not a
    ///   test-only population, as this used to say — [`Model::push_plane`] is public and leaves a
    ///   nameless plane behind for a collinear triple, so a caller outside this crate can reach
    ///   it. Every production push passes a non-collinearity gate first, which is why
    ///   `every_plane_still_has_a_name` holds; that is *we do not*, not *it cannot*.
    ///   [`SketchFrame::canonical`] makes no claim and so cannot answer here — the same plane
    ///   reaches this error one step later, inside the operation, and
    ///   `a_plane_with_no_name_cannot_host_a_sketch` pins that the two doors agree.
    /// * [`OpError::OriginNotOnPlane`] — the stated origin's residual against the plane's name
    ///   is nonzero. Exact, total ([`nacre_scalar::plane_residual_sign`]): a `Wide` name checks
    ///   through arbitrary precision, never a shrug.
    /// * [`OpError::RefDirParallelToNormal`] — `ref_dir`'s projection into the plane vanishes
    ///   (parallel to the normal, or zero), so it picks no `+u`.
    ///
    /// `ref_dir` need not lie in the plane — its in-plane part is what names `+u` (projected
    /// exactly where the frame is built). It need not be unit either; only its direction speaks.
    pub fn named(
        model: &Model,
        plane: Handle<Surface>,
        origin: Point3,
        ref_dir: Vector3,
    ) -> Result<SketchFrame, OpError> {
        let lift = |v: [f64; 3]| -> Result<[nacre_scalar::Rat; 3], OpError> {
            let mut out = [nacre_scalar::Rat::from_int(0); 3];
            for (o, c) in out.iter_mut().zip(v) {
                *o =
                    nacre_scalar::Rat::from_decimal(c).ok_or(OpError::FrameOutsideDecimalWindow)?;
            }
            Ok(out)
        };
        let (origin, ref_dir) = (lift(origin.as_array())?, lift(ref_dir.as_array())?);
        let name = model
            .surface_name
            .get(&plane)
            .ok_or(OpError::PlaneWithoutExactForm)?;
        if nacre_scalar::plane_residual_sign(name, origin) != 0 {
            return Err(OpError::OriginNotOnPlane);
        }
        // The existing frame machinery is the judge — `WideFrame::named_of` is total in width
        // (arbitrary precision), so its only `None` is the projected `u_raw` vanishing: exactly
        // the parallel-or-zero claim this constructor rejects. (The narrow `plane_frame_named`
        // is not consulted here: its `None` can mean `i128` overflow, which is a width fact,
        // not a defect in the claim.)
        if nacre_cip::WideFrame::named_of(name, &origin, &ref_dir, false).is_none() {
            return Err(OpError::RefDirParallelToNormal);
        }
        Ok(SketchFrame {
            plane,
            placement: nacre_topo::FramePlacement::Named { origin, ref_dir },
            flip: false,
        })
    }

    /// The same frame, with its plane re-anchored onto another model's arena — `replay`'s index
    /// vocabulary, and nothing else. Crate-private because it is not a claim a caller can make:
    /// only the code that owns the model being built knows the index means anything there.
    pub(crate) fn rebound(self, plane: Handle<Surface>) -> SketchFrame {
        SketchFrame { plane, ..self }
    }

    /// The plane this frame sketches on — one handle, shared with every face on that plane.
    pub fn plane(&self) -> Handle<Surface> {
        self.plane
    }

    /// Where the frame's origin and `+u` come from: `Canonical` (derived) or `Named` (stated).
    pub fn placement(&self) -> &nacre_topo::FramePlacement {
        &self.placement
    }

    /// Whether the plane's canonical coefficients are negated to face the way the sketch does —
    /// measured by the consuming operation, `false` as constructed.
    pub fn flip(&self) -> bool {
        self.flip
    }
}

/// A [`SketchFrame`] with its `flip` measured — the placement's realized `ŵ` dotted against the
/// direction the sketch must face (a sweep's sense, a face's outward normal). This is the one
/// place `flip` is decided (S9): every road calls it, so no two can measure differently. `None`
/// when the chain cannot realize a basis (a plane with no name).
///
/// ★★★ **Only the frame comes back — deliberately.** The basis realized here to measure `flip`
/// is the *unflipped* one, and it used to ride along "so deciding and looking cost one
/// realization, not two". That saving is what broke `face_plane`'s contract: the one caller who
/// wanted the axes combined the measured `flip` with the unflipped basis **by hand**, as a
/// half-turn about `v̂` — while the realization (`frame_chain`) half-turns about `û` — and every
/// flip=true face was reported a frame point-symmetric to the one the pad actually built in
/// (measured: a footprint centred on the face through `face_plane`'s own coordinates landed
/// outside it, `PadMissesFace` on 2 of 6 faces of a turned block). A caller that needs the axes
/// asks [`crate::rotated_vertex::frame_world_basis`] *with the measured flip*, so the geometry of
/// `flip` is written in exactly one place; the second 4-point replay is one plain f64 chain per
/// user operation, which is what the hand-combination was saving.
fn measured_frame(
    model: &Model,
    plane: Handle<Surface>,
    placement: nacre_topo::FramePlacement,
    toward: Vector3,
) -> Option<SketchFrame> {
    let basis = crate::rotated_vertex::frame_world_basis(model, plane, &placement, false)?;
    let n = toward.as_array();
    let flip = (0..3).map(|k| basis.3[k] * n[k]).sum::<f64>() < 0.0;
    Some(SketchFrame {
        plane,
        placement,
        flip,
    })
}

/// Name a [`SketchFrame`] as the [`nacre_topo::Motion::Frame`] node the sweep writes coordinates
/// against — the one road from the frame value to a node, shared by the extrude and face paths.
///
/// ★★ **`push_motion` interns**, so two sketches in one frame name the *same* node — which is
/// what makes their surfaces intern too (`SurfaceKey` is `(name, motion)`): two routes to one
/// height become one `Handle<Surface>` at construction, with no f64 comparison anywhere. With
/// `Canonical` placement the node is `(plane, Canonical, flip)` — nothing per-sketch in the key.
fn push_frame_node(model: &mut Model, frame: SketchFrame) -> Handle<nacre_topo::MotionNode> {
    model.push_motion(
        nacre_topo::Motion::Frame {
            plane: frame.plane,
            placement: frame.placement,
            flip: frame.flip,
        },
        None,
    )
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

/// **Where a frame is, in space** — its origin and its two axes, realized as f64.
///
/// [`face_plane`] answers this for a face; this answers it for any frame a caller holds, which
/// is what a viewer needs to draw a sketch where it was drawn: the sketch's coordinates are
/// `(u, v)` in this frame, and `origin + u·x + v·y` is the point.
///
/// ★ It is a **report**, not a truth. The frame's statement is the truth — this is that
/// statement realized, with all the rounding a realization carries, and nothing exact should be
/// decided from it.
///
/// `None` when the frame's chain cannot be realized at all (a plane with no name). A caller
/// that cannot place a thing should decline to draw it rather than draw it somewhere wrong.
pub fn frame_plane(model: &Model, frame: &SketchFrame) -> Option<SketchPlane> {
    let basis = crate::rotated_vertex::frame_world_basis(
        model,
        frame.plane(),
        frame.placement(),
        frame.flip(),
    )?;
    Some(realized_plane(
        Point3::from_array(basis.0),
        Vector3::from_array(basis.1),
        Vector3::from_array(basis.2),
    ))
}

/// A planar face's sketch frame **as a [`SketchFrame`]** — the plane handle, placement, and
/// measured flip that [`Operation::PadOnFace`] / [`Operation::PocketOnFace`] sketch in. Where
/// [`face_plane`] projects that frame to realized f64 axes for a caller to *look at*, this is
/// the exact vocabulary itself (S9): the same value `face_frame` builds internally, no longer
/// thrown away at the boundary.
///
/// ★★★ **What comes back is verified against the pad's frame, by realization, to the bit.** On a
/// world-branch face (axes lifting exactly) the operation elides the frame node and sketches in
/// `face_frame`'s world axes — so this transcribes *those* axes into the vocabulary and returns a
/// candidate only if realizing it lands bit-identically on them. The old fallback assumed the
/// canonical frame realizes to the same axes ("the node-omission normalization"); that holds only
/// for flip=false faces of motion-free planes, and everywhere else the returned frame put a
/// sketch somewhere the pad does not (measured: point-symmetric on every flip=true axis-aligned
/// face). Verification is the contract now — no candidate can be returned wrong, whatever
/// population shows up next.
///
/// A face has no caller to name a placement, so the canonical frame is tried first (the stronger
/// normal form); where it realizes elsewhere, the pad's axes are transcribed as a `Named`
/// placement (origin + `ref_dir`, S9's vocabulary) with `flip` measured as everywhere else.
///
/// Errors as [`face_plane`]: `NonPlanarFace`, `FaceNotInLiveSolid`; `PlaneWithoutExactForm` when
/// the plane carries no name to derive a frame from (a test-only unregistered surface); and
/// [`OpError::FrameNotRepresentable`] when the frame exists but no spelling realizes to it —
/// world-branch faces whose surface carries a motion; since the invariant-plane restatement
/// a plane its motion fixes carries none, so the residual population is the recorded one
/// (see the variant's doc).
pub fn face_sketch_frame(model: &Model, face: Handle<Face>) -> Result<SketchFrame, OpError> {
    let f = face_frame(model, face)?;
    if let Some(sf) = f.sketch_frame {
        return Ok(sf);
    }
    // The world-branch population: the operation will build no node and sketch in `f`'s axes.
    // Every candidate below must prove itself by realizing to exactly those axes — bits, not a
    // tolerance: both sides come from exact roads ({0,±1} axes, rational projections), so
    // agreement is exact when it holds and a threshold would only paper over a third derivation.
    let pad_frame =
        [f.origin.as_array(), f.x.as_array(), f.y.as_array()].map(|c| c.map(f64::to_bits));
    let verified = |sf: SketchFrame| -> Option<SketchFrame> {
        let (o, u, v, _) =
            crate::rotated_vertex::frame_world_basis(model, f.surface_h, sf.placement(), sf.flip)?;
        ([o, u, v].map(|c| c.map(f64::to_bits)) == pad_frame).then_some(sf)
    };
    let canonical = || {
        measured_frame(
            model,
            f.surface_h,
            nacre_topo::FramePlacement::Canonical,
            f.n,
        )
    };
    // The transcription: the pad's own origin and +u, said as a `Named` placement. `named`'s
    // exact checks (on-plane origin, non-degenerate ref_dir) ride along; any failure just drops
    // the candidate — the refusal below is the answer, never a silent wrong frame.
    let transcribed = || {
        let sf = SketchFrame::named(model, f.surface_h, f.origin, f.x).ok()?;
        measured_frame(model, f.surface_h, sf.placement, f.n)
    };
    if !model.surface_name.contains_key(&f.surface_h) {
        // No name at all: nothing can realize. The distinct, older proposition.
        return Err(OpError::PlaneWithoutExactForm);
    }
    canonical()
        .and_then(&verified)
        .or_else(|| transcribed().and_then(&verified))
        .ok_or(OpError::FrameNotRepresentable)
}

/// Locate `face`'s live solid and build its planar frame. `NonPlanarFace` for a curved surface,
/// `FaceNotInLiveSolid` if no live outer shell holds it.
fn face_frame(model: &Model, face: Handle<Face>) -> Result<FaceFrame, OpError> {
    let (solid_h, _) = model
        .live_solids
        .iter()
        .map(|&s| (s, model.solids.get(s).outer))
        // Index-only equality again: a face handle from another model can match here. The
        // shell lookup that follows is where the cross-store guard fires.
        .find(|&(_, sh)| model.shells.get(sh).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;
    let f = model.faces.get(face);
    let surface_h = f.surface;
    let orientation = f.orientation;
    let plane = match model.surface(surface_h) {
        nacre_geom::Surface::Plane(p) => *p,
        nacre_geom::Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    };
    let sign = f64::from(orientation.sign());
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
    // Only a world-stated plane's coefficients are world truth (`motion: None`). A moved
    // surface records its **pre-motion** frame, so projecting those gives a pre-motion point —
    // measured `0.29` away from the world plane, not a rounding but a different place. Those
    // keep the f64 projection, which is the same rule computed from the description available.
    // `narrow()` gates the wide vessel out: a `Wide` name (S2) carries identity only, so it
    // keeps the f64 projection exactly as a missing name did.
    let world_stated = matches!(
        model.surface_truth(surface_h),
        nacre_topo::Surface::Plane { motion: None, .. }
    );
    let origin = match (
        world_stated,
        model.surface_name.get(&surface_h).and_then(|n| n.narrow()),
    ) {
        (true, Some(&c)) => nacre_scalar::plane_origin_projection(c)
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
    // ★ **`flip` is measured, not derived** — by `measured_frame`, the one measuring place (S9).
    // ★★★ **The axes are then realized *with* that flip, by the same function the operation's
    // prism replays through.** They used to be read off the unflipped basis with the sign applied
    // by hand here, as a half-turn about `v̂` — but the realization (`frame_chain`) half-turns
    // about `û` (its `ref_dir` is derived from the unflipped coefficients and survives the sign),
    // so every flip=true face was reported a frame point-symmetric to the one the pad built in.
    // Asking `frame_world_basis` with the measured flip leaves the geometry of `flip` written in
    // exactly one place; for flip=false the call is bit-identical to the measuring one.
    // ★★ A face has no caller to name a frame, so its placement is `Canonical` (S4) — derived
    // when the chain is flattened, stored nowhere. That is also what opens this branch for a
    // plane whose name is `Wide` or whose canonical values overflow `i128`: `frame_world_basis`
    // succeeds through the arbitrary-precision road where the old narrow derivation declined.
    let world = realized_plane(origin, x, y);
    let sketch = (world.exact().is_none())
        .then(|| measured_frame(model, surface_h, nacre_topo::FramePlacement::Canonical, n))
        .flatten()
        .and_then(|sf| {
            let (o, u, v, _) = crate::rotated_vertex::frame_world_basis(
                model,
                surface_h,
                &nacre_topo::FramePlacement::Canonical,
                sf.flip,
            )?;
            Some((
                sf,
                Point3::from_array(o),
                Vector3::from_array(u),
                // ★★ **`v̂` as realized, not as `ŵ × û` recomputed here.** It has its own exact
                // rational form (`plane_frame`), so realizing it costs one rounding where a cross
                // product costs two that do not cancel — measured, a wall whose `v` is exactly
                // `ẑ` came back three ulps short of `1.0` through the cross product.
                Vector3::from_array(v),
            ))
        });
    let (x, y, origin, sketch_frame) = match sketch {
        Some((sf, o, u, v)) => (u, v, o, Some(sf)),
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
/// no measured tolerance); one that **reaches past the face** routes to the overhang boolean
/// sidecars (a boss cantilever / an edge slot; Discovered seam vertices). `Fuse` sweeps **outward**
/// (a boss); `Cut` sweeps **inward** (a blind pocket). Returns the result solid and the feature's
/// exposed cap — the boss top or the pocket floor, the outer-shell face on the prism's far-cap plane
/// with outward normal `+n`. A cap that did not survive (a through-cut has no floor) is `no_cap`,
/// **the caller's error raised here** rather than an `Option` the caller turns into one: only this
/// scope holds the boolean's result solids, and restoring `live_solids` from them is what keeps a
/// reject from committing. `NonPositiveDistance`/`DegenerateProfile` propagate from the frame;
/// overhang configurations the boolean does not cover surface as `Boolean(_)`.
fn extrude_and_boolean(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
    kind: BoolKind,
    no_cap: OpError,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    if nacre_scalar::Rat::from_decimal(dist).is_none() {
        return Err(OpError::DistOutsideDecimalWindow);
    }
    profile.check()?;
    refuse_non_quarter_arcs(profile)?;
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
    let sketch_frame = frame.sketch_frame.map(|f| push_frame_node(model, f));
    let (outer, holes) = swept_profile(model, &plane, profile, signed, sketch_frame)?;
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
        Some((solid, cap)) => Ok((solid, cap)),
        // Nothing carries the cap — either the prism reached through (a pocket with no floor) or
        // it removed the solid outright, which is the same verdict taken to its limit. Only `Cut`
        // can empty a result: a `Fuse` of two non-empty solids is never empty.
        //
        // ★ **The restore is the whole reason this arm lives here.** `assemble_fuse_cut` already
        // retired the operand and installed its own results, so returning an error now would hand
        // the caller a failure *and* a model it never asked for: its solid gone, a through-cut in
        // its place. Putting `live_solids` back is what makes the reject true from the outside —
        // the same move `PadMissesFace` makes above, for the same reason. (The arena keeps the
        // prism's cells; the store is append-only. That residue is why a session must rebuild
        // from its log before recording again — see `tests/replay.rs`.)
        None => {
            debug_assert!(
                !solids.is_empty() || matches!(kind, BoolKind::Cut),
                "a Fuse cannot produce an empty result"
            );
            model.live_solids.retain(|s| !solids.contains(s));
            model.live_solids.push(frame.solid_h);
            Err(no_cap)
        }
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
    extrude_and_boolean(
        model,
        face,
        profile,
        dist,
        BoolKind::Fuse,
        OpError::DegenerateGeometry,
    )
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
    extrude_and_boolean(
        model,
        face,
        profile,
        dist,
        BoolKind::Cut,
        OpError::PocketNotBlind,
    )
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
///    A `None` from `outer_tri` (no non-collinear triple) means no evidence *for this branch* —
///    such a candidate can still match by handle, and a degenerate `reference` leaves only
///    branch 1.
///
/// The direction filter reads the **candidate's** outward normal against `want`, never
/// `reference`'s: a pocket's tool cap faces along the sweep (`−n`) while the floor it becomes faces
/// back into the void (`+n`). Outward is the face's *stated* one — `plane.normal()` ×
/// `orientation`, the same cutover `collect_planes` made — not a re-derivation from its loop.
/// Coplanarity is settled by then, so the two are parallel and the dot is a full magnitude away
/// from zero — an f64 read whose sign cannot round the wrong way.
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
        let face = model.faces.get(fh);
        let nacre_geom::Surface::Plane(pl) = model.surface(face.surface) else {
            return false;
        };
        let coplanar = face.surface == ref_surf
            || ref_tri.is_some_and(|r| {
                outer_tri(model, face)
                    .is_some_and(|(tri, _)| tri.iter().all(|&q| plane_side(r, q) == 0))
            });
        let sign = f64::from(face.orientation.sign());
        coplanar && pl.normal().dot(want) * sign > 0.0
    })
}

/// **Do the two roads to a sketch frame agree?** — the measurement the vocabulary swap rests on.
///
/// Today an extrude reads its frame from the `SketchPlane` a caller handed in: `exact()` lifts the
/// caller's f64 axes and, when they lift to exact orthonormal rationals, the whole prism is built
/// in world coordinates with no motion node. When `Operation::Extrude` starts naming a
/// [`SketchFrame`] instead, the axes will come from **realizing the frame's exact form** — a
/// different route to the same real vectors.
///
/// ★★★ **A one-ulp disagreement there is not an ulp of error.** `exact()` would flip to `None`,
/// the plane would silently take the frame-node road, and the arena would gain motion nodes and
/// write its points in frame coordinates: a *different but still valid* model. So the assertion
/// order below matters — **same road first**, values second. This crate has been bitten by exactly
/// this shape before: `nacre_scalar::plane_frame_named` records `v̂` realized as `ŵ × û` coming out
/// `0.999999999999999_7`, "an exact path quietly lost".
///
/// The prediction is agreement, and it is structural rather than lucky: a datum's `ref_dir` is
/// `points[1] − points[0]`, which *is* the caller's `+u`; for a unit rational axis `|u_raw|² = 1`
/// so `inv_sqrt_exact` returns exactly one; and `v_raw` has its own exact form. But an argument is
/// not a gate.
#[cfg(test)]
mod frame_road {
    use super::*;
    use crate::exact::RatFrame;

    /// The frame road's realized basis, wrapped so `exact()` can be asked of it — the same two
    /// lines `extrude_on_frame` will run.
    fn realized(model: &Model, frame: &SketchFrame) -> Option<SketchPlane> {
        let (o, u, v, _) = crate::rotated_vertex::frame_world_basis(
            model,
            frame.plane(),
            frame.placement(),
            frame.flip(),
        )?;
        Some(realized_plane(
            Point3::from_array(o),
            Vector3::from_array(u),
            Vector3::from_array(v),
        ))
    }

    /// State `sp` as a datum and hand back both roads' frames for it.
    fn both(model: &mut Model, sp: SketchPlane) -> (Option<RatFrame>, Option<RatFrame>) {
        let OpOutput::DatumPlane { frame, .. } = apply(
            model,
            &Operation::DatumPlane {
                def: DatumDef::Stated(sp),
            },
        )
        .expect("every fixture here states a plane the kernel can hold") else {
            unreachable!()
        };
        let by_value = sp.exact();
        // ★ The repaired road: ask in `Rat`, never realize.
        let by_frame = exact_frame(model, &frame);
        (by_value, by_frame)
    }

    /// The road as it was before the repair — realize the axes, then lift them back. Kept so the
    /// repair is visibly a different *question* rather than a refactor of the same one.
    fn both_realized(model: &mut Model, sp: SketchPlane) -> (Option<RatFrame>, Option<RatFrame>) {
        let OpOutput::DatumPlane { frame, .. } = apply(
            model,
            &Operation::DatumPlane {
                def: DatumDef::Stated(sp),
            },
        )
        .expect("every fixture here states a plane the kernel can hold") else {
            unreachable!()
        };
        (sp.exact(), realized(model, &frame).and_then(|p| p.exact()))
    }

    fn population() -> Vec<(&'static str, SketchPlane)> {
        let p3 = Point3::from_array;
        let v3 = Vector3::from_array;
        vec![
            ("world_xy", SketchPlane::world_xy()),
            ("world_yz", SketchPlane::world_yz()),
            // ★ The axis whose derived frame is not the convention.
            ("world_zx", SketchPlane::world_zx()),
            (
                "offset_xy",
                SketchPlane::world_xy().with_origin(p3([0.0, 0.0, 0.5])),
            ),
            (
                "far_offset_xy",
                SketchPlane::world_xy().with_origin(p3([50.0, -37.25, 0.5])),
            ),
            // ★ Rational tilt: axes lift exactly, so this takes the world road today — the only
            // family where the two routes could disagree and it would matter.
            (
                "rational_tilt_wf",
                SketchPlane::from_axes(
                    p3([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
                    v3([0.6, 0.8, 0.0]),
                    v3([-0.48, 0.36, 0.8]),
                ),
            ),
            (
                "rational_tilt_345",
                SketchPlane::from_axes(
                    p3([1.0, 2.0, 3.0]),
                    v3([0.6, 0.8, 0.0]),
                    v3([0.0, 0.0, 1.0]),
                ),
            ),
            // Irrational tilt: `exact()` declines on both roads — the agreement that matters here
            // is that they decline *together*.
            (
                "irrational_tilt",
                SketchPlane::from_origin_normal(p3([0.0; 3]), v3([1.0, 1.0, 1.0]))
                    .expect("a plane"),
            ),
            // ★ A real population: three call sites extrude on a plane facing −ŷ.
            (
                "negative_normal",
                SketchPlane::from_origin_normal(p3([0.0, 1.3, 0.0]), v3([0.0, -1.0, 0.0]))
                    .expect("a plane"),
            ),
            (
                "through_points",
                SketchPlane::through_points(
                    p3([1.0, 0.0, 0.0]),
                    p3([1.0, 2.0, 0.0]),
                    p3([1.0, 0.0, 3.0]),
                )
                .expect("a plane"),
            ),
        ]
    }

    /// ★★ **They agree — once the question is asked in rationals.**
    ///
    /// Before the repair two families diverged (`rational_tilt_wf`, `rational_tilt_345`): an axis
    /// survived the frame road only when it needed no normalizing, because `reduce_direction`
    /// turns a `(0.6, 0.8, 0)` axis into the primitive `(3, 4, 0)` with `uu = 25` and a
    /// *realization* multiplies by a numerically computed `1/5`, landing on `0.6000000000000001`.
    /// `Rat::from_decimal` lifted that, orthonormality failed, and a perfectly rational plane
    /// took the frame-node road — a different arena, not an ulp.
    ///
    /// ★ Two hypotheses died on the way, and both were mine. **Origin cancellation** (`1.6 − 1.0`)
    /// is not the cause: a unit axis at `(1, 2, 3)` survives, which
    /// [`a_realized_axis_survives_only_when_it_needs_no_normalizing`] pins. And the planning
    /// argument "`|u_raw|² = 1`, so `inv_sqrt_exact` returns exactly one" held only for axes that
    /// are *already* unit; I generalized it to every rational axis and predicted agreement.
    ///
    /// [`RatFrame::of_plane_frame`] asks in `Rat` and never realizes — the same rule
    /// `plane_frame_named` already states for `v̂`, one level up.
    #[test]
    fn the_two_roads_take_the_same_road() {
        for (name, sp) in population() {
            let mut m = Model::new();
            let (by_value, by_frame) = both(&mut m, sp);
            println!(
                "stat frame_road {name} exact_by_value={} exact_by_frame={}",
                by_value.is_some(),
                by_frame.is_some()
            );
            assert_eq!(
                by_value.is_some(),
                by_frame.is_some(),
                "{name}: the two roads disagree about whether the axes are rational. That is not \
                 an ulp — it moves the plane between the world road and the frame-node road, and \
                 the arena differs by whole motion nodes."
            );
        }
    }

    /// ★★ **The road the repair replaced, kept as a positive control.**
    ///
    /// Realizing the axes and lifting them back still loses exactly the two families whose axes
    /// need normalizing. Without this the repair would read as a refactor of one question; with
    /// it, the two routes are visibly different questions and the fix is visibly load-bearing.
    #[test]
    fn realizing_the_axes_first_still_loses_them() {
        let mut lost = Vec::new();
        for (name, sp) in population() {
            let mut m = Model::new();
            let (by_value, by_realized) = both_realized(&mut m, sp);
            if by_value.is_some() != by_realized.is_some() {
                lost.push(name);
            }
        }
        assert_eq!(
            lost,
            vec!["rational_tilt_wf", "rational_tilt_345"],
            "the realized route is what the repair stopped using; if this list changed, the \
             measurement the repair was built on moved with it"
        );
    }

    /// **The mechanism, isolated** — so the repair is checked against the cause, not the symptom.
    ///
    /// On the realized route an axis survives exactly when its reduced form is already unit;
    /// moving the frame's origin changes nothing, which is what ruled out the "differencing two
    /// realized points cancels the origin" explanation. Asked in `Rat`, all four are rational.
    #[test]
    fn a_realized_axis_survives_only_when_it_needs_no_normalizing() {
        let p3 = Point3::from_array;
        let v3 = Vector3::from_array;
        for (name, origin, u) in [
            ("unit_axis_at_origin", p3([0.0; 3]), v3([1.0, 0.0, 0.0])),
            ("unit_axis_moved", p3([1.0, 2.0, 3.0]), v3([1.0, 0.0, 0.0])),
            ("scaled_axis_at_origin", p3([0.0; 3]), v3([0.6, 0.8, 0.0])),
            (
                "scaled_axis_moved",
                p3([1.0, 2.0, 3.0]),
                v3([0.6, 0.8, 0.0]),
            ),
        ] {
            let mut m = Model::new();
            let sp = SketchPlane::from_axes(origin, u, v3([0.0, 0.0, 1.0]));
            assert!(sp.exact().is_some(), "{name}: caller axes lift");
            let (_, by_realized) = both_realized(&mut m, sp);
            let (_, by_rat) = both(&mut m, sp);
            let unit_axis = u.as_array().iter().filter(|c| **c != 0.0).count() == 1;
            println!(
                "stat frame_axis_loss {name} realized={} rational={}",
                by_realized.is_some(),
                by_rat.is_some()
            );
            assert!(
                by_rat.is_some(),
                "{name}: asked in Rat, every one of these frames is rational"
            );
            assert_eq!(
                by_realized.is_some(),
                unit_axis,
                "{name}: the realized route keeps an axis it does not have to normalize and \
                 loses one it does — whatever the origin is"
            );
        }
    }

    /// Where both roads *do* reach the exact road, the rational frames are identical — so the
    /// loss above is the only thing standing between the two vocabularies.
    #[test]
    fn the_two_roads_realize_the_same_axes_where_they_agree_on_the_road() {
        let mut agreed = 0;
        for (name, sp) in population() {
            let mut m = Model::new();
            let (by_value, by_frame) = both(&mut m, sp);
            let (Some(a), Some(b)) = (by_value, by_frame) else {
                println!("stat frame_axes {name} both_declined");
                continue;
            };
            assert_eq!(
                a, b,
                "{name}: same road, different axes — the prism would be built on a different \
                 rational frame and every vertex would move"
            );
            agreed += 1;
            println!("stat frame_axes {name} identical");
        }
        assert!(
            agreed >= 5,
            "only {agreed} fixtures reached the exact road; the population stopped measuring what \
             it was chosen to measure"
        );
    }

    /// ★★ **The invariant the agreement above quietly rests on, stated.**
    ///
    /// `exact()` lifts the frame's *realized* origin with `Rat::from_decimal`, so a `Named`
    /// placement's origin has to survive a round trip through `f64` to come back as the same
    /// rational. Every `Named` origin today came from `Rat::from_decimal` in the first place
    /// (`SketchFrame::named` lifts an `f64`; a datum takes `PlaneDef`'s lifted points), so it
    /// does — but that is a property of the current producers, **not of the type**. A future
    /// producer that *computes* an origin would break the round trip, and the plane would move to
    /// the frame-node road without anything failing.
    #[test]
    fn every_named_origin_survives_the_round_trip_it_is_relied_on_for() {
        for (name, sp) in population() {
            let mut m = Model::new();
            let OpOutput::DatumPlane { frame, .. } = apply(
                &mut m,
                &Operation::DatumPlane {
                    def: DatumDef::Stated(sp),
                },
            )
            .expect("stated") else {
                unreachable!()
            };
            let nacre_topo::FramePlacement::Named { origin, .. } = frame.placement() else {
                panic!("{name}: a stated datum always names its placement");
            };
            for (k, r) in origin.iter().enumerate() {
                let round = nacre_scalar::Rat::from_decimal(r.to_f64());
                assert_eq!(
                    round,
                    Some(*r),
                    "{name}: origin[{k}] does not survive f64 — `exact()` lifts the realized \
                     origin, so this is what keeps the frame road on the world road"
                );
            }
        }
    }
}

/// **The handle road and the value road build the same prism.**
///
/// [`Operation::Extrude`] still carries a [`SketchPlane`]; [`extrude_on_frame`] is the road it is
/// about to take, and nothing calls it yet. Before the call sites move, this asks the only
/// question that matters: given the *same* plane said two ways, do the two roads leave the same
/// arena behind?
///
/// ★ The first assertion is **which road was taken**, not whether the models match. A frame whose
/// basis is not rational is written inside a motion node with its points in frame coordinates — a
/// different but still valid model, and an equality check would report that without naming it.
/// [`frame_road`] measured that this used to happen for rational-tilt planes; here it must not
/// happen at all.
#[cfg(test)]
mod frame_differential {
    use super::*;

    /// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
    /// when the plane is not one the model already holds (a seed, or a face's).
    fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
        match apply(
            m,
            &Operation::DatumPlane {
                def: DatumDef::Stated(plane),
            },
        ) {
            Ok(OpOutput::DatumPlane { frame, .. }) => frame,
            other => panic!("stating a plane: {other:?}"),
        }
    }
    use nacre_topo::Surface;

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Profile2d {
        Profile2d::polygon(vec![
            Point2::from_array([x0, y0]),
            Point2::from_array([x1, y0]),
            Point2::from_array([x1, y1]),
            Point2::from_array([x0, y1]),
        ])
        .expect("a rectangle")
    }

    /// A ring with a hole, so the differential covers the inner-loop plumbing too.
    fn washer() -> Profile2d {
        let ring = |a: f64, b: f64| {
            vec![
                Point2::from_array([a, a]),
                Point2::from_array([b, a]),
                Point2::from_array([b, b]),
                Point2::from_array([a, b]),
            ]
        };
        Profile2d::with_holes(ring(0.0, 4.0), vec![ring(1.0, 3.0)]).expect("a washer")
    }

    /// Everything about a model that a road could change, named and indexed.
    fn arena(m: &Model) -> Vec<(String, String)> {
        let mut out = vec![(
            "len".into(),
            format!(
                "{} {} {} {} {} {}",
                m.vertices.len(),
                m.edges.len(),
                m.faces.len(),
                m.shells.len(),
                m.solids.len(),
                m.surface_count()
            ),
        )];
        for (h, _) in m.vertices.iter() {
            let p = m.vertex_point(h).as_array();
            out.push((
                format!("v{}", h.index()),
                format!(
                    "{:x},{:x},{:x}",
                    p[0].to_bits(),
                    p[1].to_bits(),
                    p[2].to_bits()
                ),
            ));
        }
        for (h, e) in m.edges.iter() {
            out.push((
                format!("e{}", h.index()),
                format!(
                    "{},{} {},{}",
                    e.surfaces[0].index(),
                    e.surfaces[1].index(),
                    e.vertices[0].index(),
                    e.vertices[1].index()
                ),
            ));
        }
        for (h, f) in m.faces.iter() {
            let loops: Vec<String> = std::iter::once(&f.outer)
                .chain(f.inner.iter())
                .map(|lp| {
                    lp.half_edges
                        .iter()
                        .map(|he| {
                            format!("{}{}", he.edge.index(), if he.forward { "+" } else { "-" })
                        })
                        .collect::<Vec<_>>()
                        .join(" ")
                })
                .collect();
            out.push((
                format!("f{}", h.index()),
                format!(
                    "s{} {:?} {}",
                    f.surface.index(),
                    f.orientation,
                    loops.join(" | ")
                ),
            ));
        }
        out.push((
            "live".into(),
            m.live_solids
                .iter()
                .map(|h| h.index().to_string())
                .collect::<Vec<_>>()
                .join(","),
        ));
        out
    }

    /// How many motion nodes a model holds — the "which road" signal. A world-road prism makes none.
    fn nodes(m: &Model) -> usize {
        // Motion handles are only reachable through surfaces' truth; count the distinct leaves.
        let mut seen = std::collections::BTreeSet::new();
        for i in 0..m.surface_count() as u32 {
            let h = m.surface_handle_at(i).expect("in range");
            let motion = match m.surface_truth(h) {
                Surface::Plane { motion, .. } | Surface::Cylinder { motion, .. } => *motion,
            };
            if let Some(node) = motion {
                seen.insert(node.index());
            }
        }
        seen.len()
    }

    fn population() -> Vec<(&'static str, SketchPlane)> {
        let p3 = Point3::from_array;
        let v3 = Vector3::from_array;
        vec![
            ("world_xy", SketchPlane::world_xy()),
            ("world_yz", SketchPlane::world_yz()),
            ("world_zx", SketchPlane::world_zx()),
            (
                "offset_xy",
                SketchPlane::world_xy().with_origin(p3([0.0, 0.0, 0.5])),
            ),
            (
                "far_offset_xy",
                SketchPlane::world_xy().with_origin(p3([50.0, -37.25, 0.5])),
            ),
            (
                "rational_tilt_wf",
                SketchPlane::from_axes(
                    p3([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
                    v3([0.6, 0.8, 0.0]),
                    v3([-0.48, 0.36, 0.8]),
                ),
            ),
            (
                "rational_tilt_345",
                SketchPlane::from_axes(
                    p3([1.0, 2.0, 3.0]),
                    v3([0.6, 0.8, 0.0]),
                    v3([0.0, 0.0, 1.0]),
                ),
            ),
            (
                "irrational_tilt",
                SketchPlane::from_origin_normal(p3([0.0; 3]), v3([1.0, 1.0, 1.0]))
                    .expect("a plane"),
            ),
            (
                "negative_normal",
                SketchPlane::from_origin_normal(p3([0.0, 1.3, 0.0]), v3([0.0, -1.0, 0.0]))
                    .expect("a plane"),
            ),
        ]
    }

    /// ★★ **The gate.** Same plane, said two ways, one arena.
    #[test]
    fn the_two_roads_build_the_same_prism() {
        for profile_name in ["rect", "washer"] {
            for (name, sp) in population() {
                let profile = || {
                    if profile_name == "rect" {
                        rect(0.1, 0.2, 2.3, 1.7)
                    } else {
                        washer()
                    }
                };
                let what = format!("{name}/{profile_name}");

                // The value road: today's operation.
                let mut by_value = Model::new();
                let __frame0 = datum_frame(&mut by_value, sp);
                apply(
                    &mut by_value,
                    &Operation::Extrude {
                        frame: __frame0,
                        profile: profile(),
                        dist: 1.5,
                    },
                )
                .unwrap_or_else(|e| panic!("{what}: value road failed: {e:?}"));

                // The handle road: state the plane, then extrude on the frame it hands back.
                let mut by_frame = Model::new();
                let OpOutput::DatumPlane { frame, .. } = apply(
                    &mut by_frame,
                    &Operation::DatumPlane {
                        def: DatumDef::Stated(sp),
                    },
                )
                .unwrap_or_else(|e| panic!("{what}: datum failed: {e:?}")) else {
                    unreachable!()
                };
                extrude_on_frame(&mut by_frame, &frame, &profile(), 1.5)
                    .unwrap_or_else(|e| panic!("{what}: frame road failed: {e:?}"));

                // ★ Which road, first. A frame-node road is a valid model too, so an equality check
                // alone would report the difference without naming its cause.
                assert_eq!(
                    nodes(&by_frame),
                    nodes(&by_value),
                    "{what}: the two roads disagree about whether a motion node is needed — the \
                     arena differs by whole nodes, not by an ulp"
                );

                let (a, b) = (arena(&by_value), arena(&by_frame));
                for (x, y) in a.iter().zip(&b) {
                    assert_eq!(x, y, "{what}: arenas first differ at {}", x.0);
                }
                assert_eq!(a.len(), b.len(), "{what}: different cell counts");
                println!(
                    "stat frame_differential {what} identical nodes={}",
                    nodes(&by_value)
                );
            }
        }
    }

    /// ★★ **The negative control: the ZX trap is real, and the sugar is what avoids it.**
    ///
    /// The arbitrary-axis rule gives the ZX plane `+u = −x̂` while the convention — and
    /// `SketchPlane::world_zx` — says `+u = +ẑ`. So reaching for `SketchFrame::canonical` on the
    /// ZX seed puts a caller's profile a quarter turn from where they asked, and the test above,
    /// which uses the datum's own `Named` frame, would never notice.
    ///
    /// Without this, "the sugar is just `canonical`" is a simplification that passes everything.
    #[test]
    fn the_zx_seed_without_the_sugar_turns_the_sketch() {
        let profile = rect(0.0, 0.0, 2.0, 1.0);

        let mut stated = Model::new();
        let __w0 = SketchFrame::world(&stated, Axis::Y);
        apply(
            &mut stated,
            &Operation::Extrude {
                frame: __w0,
                profile: profile.clone(),
                dist: 1.0,
            },
        )
        .expect("the convention");

        let mut derived = Model::new();
        let zx = derived.world_plane(nacre_scalar::Axis::Y);
        extrude_on_frame(&mut derived, &SketchFrame::canonical(zx), &profile, 1.0)
            .expect("the derivation is a perfectly good frame — it is just a different one");

        let (a, b) = (arena(&stated), arena(&derived));
        assert_ne!(
            a, b,
            "the derived ZX frame must differ from the convention — if these ever agree, either \
             the seeding or the arbitrary-axis rule moved, and SketchFrame::world is dead weight"
        );
        // And the sugar is what closes it: same convention, same arena.
        let mut sugared = Model::new();
        let f = SketchFrame::world(&sugared, nacre_scalar::Axis::Y);
        extrude_on_frame(&mut sugared, &f, &profile, 1.0).expect("the sugar");
        let c = arena(&sugared);
        for (x, y) in a.iter().zip(&c) {
            assert_eq!(x, y, "the sugar must reproduce the convention: {}", x.0);
        }
        println!("stat frame_differential world_zx_without_sugar differs=true sugar_matches=true");
    }
}
