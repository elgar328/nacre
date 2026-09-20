//! Operations for the nacre kernel, plus a replayable operation log.
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::PadOnFace`]/[`Operation::PocketOnFace`] (M4) consume a prior op's face by
//! `Handle` (exposed via [`OpOutput`]) and supersede a solid (live-solid
//! semantics) — each is a tool prism plus a boolean, not a direct face-split. Ops are
//! applied by [`apply`] and folded by [`replay`]; every result is a **closed** solid, so
//! `nacre-validate` applies fully.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::Rat;
use nacre_store::Handle;
use nacre_topo::{Face, HalfEdge, Model, Solid, Vertex};

mod arrangement;
/// The cylinder-band pass — see the module docs.
mod bands;
mod boolean;
mod combinatorics;
/// The cylinder chart: the lateral faces' arrangement, emitted from cells (capability D).
mod cyl_chart;
mod exact;
mod nesting;
mod ops;
mod par;
mod planes;
mod realize;
/// Public because the measurements that read it live in other crates — see the module doc.
pub mod reject_census;
mod reuse;
mod rotated_vertex;
mod sketch;
mod tolerant;
mod transform;

pub use boolean::{BoolReport, boolean, boolean_with_report};
// The report's vocabulary: what a judgement was asked about, and what it established. Re-exported
// so a consumer reads one crate, not two.
pub use nacre_cip::Decision;
pub use nacre_cip::predicate::{Evidence, Site};
pub use nacre_geom::mixed::Edge2d;
pub use ops::{
    BoolKind, DatumDef, LogCell, OpError, OpOutput, Operation, PlaneDef, Profile2d, ProfileRing,
    Ring2d, SketchFrame, SketchPlane, apply, face_plane, face_sketch_frame, frame_plane, replay,
};
pub use realize::{
    CacheDecline, Precision, RealizeError, Realized, RefineReport, realize_cache, realize_vertex,
    realize_vertex_decimal, refine_vertex_cache,
};
pub use sketch::{SketchError, arc_to_rat, arc_turns, arc_turns_rat, from_paths, from_rings};

impl SketchPlane {
    /// The world XY plane: `+u = x̂`, `+v = ŷ`, normal `+ẑ`.
    pub fn world_xy() -> Self {
        Self::axis_plane([0, 0, 1], [1, 0, 0], [0, 1, 0])
    }

    /// The world YZ plane: `+u = ŷ`, `+v = ẑ`, normal `+x̂`.
    pub fn world_yz() -> Self {
        Self::axis_plane([1, 0, 0], [0, 1, 0], [0, 0, 1])
    }

    /// The world ZX plane: `+u = ẑ`, `+v = x̂`, normal `+ŷ`.
    ///
    /// ★★ **The axes are named, not derived.** `ẑ × n` would give `−x̂` here; the convention a
    /// person expects (and the one the script layer documents) is `+u = ẑ`. A named plane gets to
    /// say, which is exactly what [`PlaneDef::ref_dir`] is for.
    pub fn world_zx() -> Self {
        Self::axis_plane([0, 1, 0], [0, 0, 1], [1, 0, 0])
    }

    /// One of the three world planes, stated exactly: normal, `+u`, `+v` as integer triples.
    /// The defining points are `[0, u, v]` — `u × v = n` for all three world planes, so the
    /// point order carries the same normal the coefficients used to state, and `points[1]`
    /// carries the *named* `+u` (the `world_zx` convention `+u = ẑ` included).
    fn axis_plane(n: [i128; 3], u: [i128; 3], v: [i128; 3]) -> Self {
        let r = |a: [i128; 3]| a.map(Rat::from_int);
        let f = |a: [i128; 3]| Vector3::from_array(a.map(|c| c as f64));
        debug_assert_eq!(
            {
                let (u, v) = (f(u), f(v));
                u.cross(v).as_array()
            },
            f(n).as_array(),
            "axis_plane point order must reproduce the stated normal"
        );
        Self {
            origin: Point3::origin(),
            x_axis: f(u),
            y_axis: f(v),
            def: Some(PlaneDef {
                points: [[Rat::from_int(0); 3], r(u), r(v)],
            }),
        }
    }

    /// A plane through `origin` with the given `normal`, its axes synthesized by the same
    /// convention a face's frame uses (`ops::frame_axes` — cross the world `ẑ` into the normal,
    /// or `ŷ` when the normal is vertical). `None` if `normal` is zero.
    ///
    /// ★ It has to be the same convention: this and [`face_plane`] answer the same question, and a
    /// caller that builds a frame here and compares it with one read off a face would otherwise
    /// find them ninety degrees apart.
    ///
    /// ★★★★ **`normal` is kept, not just consumed.** Normalizing it is what destroys the exact
    /// form — the caller's `(1, 1, 1)` is coefficients `[1, 1, 1, 0]`, while `normalize` of it
    /// squares to `0.9999999999999999…`. The unit axes below stay the f64 cache; the definition
    /// records what was handed in. A normal outside the decimal window simply leaves `def` empty,
    /// so this constructor is **never stricter than it was**.
    pub fn from_origin_normal(origin: Point3, normal: Vector3) -> Option<Self> {
        let n = normal.normalize()?;
        let (x_axis, y_axis) = ops::frame_axes(n)?;
        Some(Self {
            origin,
            x_axis,
            y_axis,
            def: Self::normal_def(origin, normal),
        })
    }

    /// `from_origin_normal`'s exact half: the plane through `origin` with normal `normal`, and the
    /// same `ẑ × n` convention spelled in rationals (un-normalized — a cross product is already in
    /// the plane, so nothing needs projecting).
    ///
    /// The definition is the three points `[o, o + u, o + w]`, and **both in-plane directions are
    /// basis crosses** — `u = ẑ × n` (or `ŷ × n` for a vertical normal), `w = x̂ × n` (or a
    /// sibling), which are component *shuffles* of the written decimals: no products, no
    /// divisions, nothing to overflow. The only declines left are a coordinate outside the
    /// decimal window and a zero normal.
    ///
    /// ★★ **Polarity is algebraic**: `(a × n) × (b × n) = det[a, b, n] · n`, so `u × w` is a
    /// *known scalar* times `n` — `n₁` for the `(ẑ, x̂)` pair — and flipping `w`'s sign when that
    /// scalar is negative makes the point order face the caller's normal, both signs, exactly.
    ///
    /// ★★★ **Why not `w = n × u`, the "obvious" second direction: its components are products.**
    /// S6a shipped that and the checked arithmetic looked like a formality; the day the silent
    /// f64 fallback stopped absorbing failures (S6b), a proptest found the window's
    /// small-exponent corner — a `10²¹` denominator squares to `10⁴²`, and even the *primitive*
    /// direction of that cross needs 137 bits. The retired `named_plane_points` solved an axis
    /// for the same reason. Basis crosses stay inside the inputs' own widths.
    fn normal_def(origin: Point3, normal: Vector3) -> Option<PlaneDef> {
        let o = origin.as_array().map(Rat::from_decimal);
        let n = normal.as_array().map(Rat::from_decimal);
        let (o, n) = ([o[0]?, o[1]?, o[2]?], [n[0]?, n[1]?, n[2]?]);
        let zero = Rat::from_int(0);
        let neg = |x: Rat| {
            zero.checked_sub(x)
                .expect("negation cannot overflow a lifted decimal")
        };
        // (u, w0, s): two independent in-plane basis crosses and the scalar with
        // `u × w0 = s · n`, per the identity above (arms chosen so `s ≠ 0`).
        let (u, w0, s) = if n[0] == zero && n[1] == zero {
            if n[2] == zero {
                return None; // a zero normal names no plane
            }
            // Vertical: u = ŷ × n, w0 = x̂ × n; det[ŷ, x̂, n] = −n₂.
            ([n[2], zero, zero], [zero, neg(n[2]), zero], neg(n[2]))
        } else if n[1] != zero {
            // The stated convention: u = ẑ × n; w0 = x̂ × n; det[ẑ, x̂, n] = n₁.
            ([neg(n[1]), n[0], zero], [zero, neg(n[2]), n[1]], n[1])
        } else {
            // n₁ = 0, n₀ ≠ 0: u = ẑ × n; w0 = ŷ × n; det[ẑ, ŷ, n] = −n₀.
            ([neg(n[1]), n[0], zero], [n[2], zero, neg(n[0])], neg(n[0]))
        };
        let w = if s > zero { w0 } else { w0.map(neg) };
        let add = |a: [Rat; 3], b: [Rat; 3]| -> Option<[Rat; 3]> {
            Some([
                a[0].checked_add(b[0])?,
                a[1].checked_add(b[1])?,
                a[2].checked_add(b[2])?,
            ])
        };
        Some(PlaneDef {
            points: [o, add(o, u)?, add(o, w)?],
        })
    }

    /// **A plane through three written points**: `origin` is the sketch's `(0, 0)`, `+u` runs
    /// toward `x_point`, and `+v` leans toward `y_hint`.
    ///
    /// ★★★ **Everything here is exact by construction.** The written points, lifted to their
    /// decimal truth, **are** the definition — origin first, so the sketch `(0, 0)` and the `+u`
    /// direction (`x_point − origin`) fall out of the structure. `None` if the three are
    /// collinear or fall outside the decimal window.
    pub fn through_points(origin: Point3, x_point: Point3, y_hint: Point3) -> Option<Self> {
        let x = (x_point - origin).normalize()?;
        let v = y_hint - origin;
        let y = (v - x * v.dot(x)).normalize()?;
        let lift = |p: Point3| {
            let a = p.as_array().map(Rat::from_decimal);
            Some([a[0]?, a[1]?, a[2]?])
        };
        let def = (|| {
            let (o, xp, yh) = (lift(origin)?, lift(x_point)?, lift(y_hint)?);
            // The written points ARE the definition; the only thing to verify is that they name
            // a plane at all. `plane_name_exact` is total (Narrow | Wide), so `None` means
            // exactly one thing: collinear. There is no
            // "answer does not fit i128" failure class.
            nacre_scalar::plane_name_exact(o, xp, yh)?;
            Some(PlaneDef {
                points: [o, xp, yh],
            })
        })();
        Some(Self {
            origin,
            x_axis: x,
            y_axis: y,
            def,
        })
    }

    /// The same plane **moved to pass through `p`**, with `p` as the sketch's `(0, 0)`.
    ///
    /// ★★★★ **The plane travels with the origin** — `plane(ZX, { origin: … })` sets the position
    /// as well as the 2-D origin, and `SketchPlane { origin, ..world_xy() }` always meant that.
    /// Keeping the plane in place while moving the origin would leave the definition describing
    /// one plane and its origin sitting on another: measured, `world_xy().with_origin([0, 0, 0.5])`
    /// recorded `z = 0` for a cap at `z = 0.5`, and a boolean built on that lost 0.04 of volume.
    /// Since the definition is three points and the origin is the first of them, the move is a
    /// translation of the whole triple: differences (`ref_dir`) and the normal are untouched,
    /// exactly the "translation does not turn `+u`" the old form promised.
    pub fn with_origin(mut self, p: Point3) -> Self {
        self.origin = p;
        self.def = self.def.and_then(|d| {
            let a = p.as_array().map(Rat::from_decimal);
            let origin = [a[0]?, a[1]?, a[2]?];
            let shift = |q: [Rat; 3]| -> Option<[Rat; 3]> {
                Some([
                    q[0].checked_sub(d.points[0][0])?.checked_add(origin[0])?,
                    q[1].checked_sub(d.points[0][1])?.checked_add(origin[1])?,
                    q[2].checked_sub(d.points[0][2])?.checked_add(origin[2])?,
                ])
            };
            Some(PlaneDef {
                points: [origin, shift(d.points[1])?, shift(d.points[2])?],
            })
        });
        self
    }

    /// A frame from axes the caller already holds — **and the axes' decimal truth is its
    /// definition**. The boundary rule that `Profile2d` applies to coordinates applies to
    /// axes too: what the caller wrote *is* the statement, so the plane through
    /// `[o, o + x, o + y]` and the `+u` direction `x` are recorded exactly. A 45°-rotated frame
    /// — whose axes never lift to exact orthonormal rationals — extrudes through the frame
    /// road (wide names included) instead of falling silently to f64.
    ///
    /// ★ The **world lift** still comes first: axes whose decimals square and cross to exact
    /// `1`/`0` — a Pythagorean frame like `(0.6, 0.8, 0)`/`(−0.48, 0.36, 0.8)` — pass `exact()`
    /// and take the world-rational path, definition or not. That population is how the
    /// `n·n`-overflow walls (the census `wf` family) are built.
    ///
    /// ★★ **What the definition states — and what it does not.** Three points, `+u`, and the
    /// polarity (point order); the realized frame is *orthonormal*, exactly as for every other
    /// definition (`ref_dir` is any length, `+v` is derived on `y`'s side). Exact decimal
    /// orthogonality is deliberately **not** required — a rotated pair `[c, s, 0]/[−s, c, 0]`
    /// cancels to exactly zero, but two independently rounded axes need not, and requiring it
    /// would strand exactly the callers this lift exists for. A skewed (non-orthogonal) pair is
    /// outside [`SketchPlane`]'s contract; the frame realization drops the skew.
    ///
    /// `def` stays `None` only for axes outside the decimal window or a degenerate pair
    /// (`x × y = 0` in the decimals) — those keep the f64 path they had.
    pub fn from_axes(origin: Point3, x_axis: Vector3, y_axis: Vector3) -> Self {
        let def = (|| {
            let lift = |p: [f64; 3]| {
                let a = p.map(Rat::from_decimal);
                Some([a[0]?, a[1]?, a[2]?])
            };
            let o = lift(origin.as_array())?;
            let add = |a: [Rat; 3], b: [Rat; 3]| -> Option<[Rat; 3]> {
                Some([
                    a[0].checked_add(b[0])?,
                    a[1].checked_add(b[1])?,
                    a[2].checked_add(b[2])?,
                ])
            };
            let px = add(o, lift(x_axis.as_array())?)?;
            let py = add(o, lift(y_axis.as_array())?)?;
            // Total: `None` means exactly one thing — the axes are parallel (or zero) in their
            // decimal truth, and name no plane.
            nacre_scalar::plane_name_exact(o, px, py)?;
            Some(PlaneDef {
                points: [o, px, py],
            })
        })();
        Self {
            origin,
            x_axis,
            y_axis,
            def,
        }
    }

    /// The sketch's `(0, 0)` in space.
    #[inline]
    pub fn origin(&self) -> Point3 {
        self.origin
    }

    /// The `+u` direction.
    #[inline]
    pub fn x_axis(&self) -> Vector3 {
        self.x_axis
    }

    /// The `+v` direction.
    #[inline]
    pub fn y_axis(&self) -> Vector3 {
        self.y_axis
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

/// Why a boolean could not be computed. The engine rejects out-of-coverage
/// input honestly rather than returning a plausibly-wrong solid (overview
/// 불리언 전략).
/// Exhaustive on purpose: new failure modes become [`RejectReason`] variants, not new variants
/// here, so a consumer can handle this enum completely and still not be broken by growth.
///
/// `Eq` is deliberately absent: [`RejectWhere`] carries `f64` coordinates, whose equality is
/// not an equivalence ([`Point3`]'s own rule). Compare rejects by projecting `reason` out —
/// never by comparing whole errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoolError {
    /// The engine declined to answer, and [`RejectReason`] names *which* guard spoke.
    /// The reason travels in the value so a consumer can say why, and so the site that
    /// returned is the site that is reported (guards that are raised and then swallowed
    /// by an alternative path cannot be mistaken for the surfaced one).
    ///
    /// `at` is where the guard was looking when it spoke — a witness for diagnostics
    /// ([`RejectWhere`]), `None` when the reason has no meaningful single location.
    ///
    /// ★ Not `Unsupported`: that name would assert every reject is a
    /// coverage limit, which is false for two of the three [`RejectClass`]es —
    /// `Impossible` (no milestone will build this input) and `SuspectedDefect` (ours, not
    /// the caller's). The variant states what happened; *what kind* of answer it is stays
    /// where it always was, in [`RejectReason::class`].
    Rejected {
        reason: RejectReason,
        at: Option<RejectWhere>,
    },
    /// An input solid handle is not in [`Model::live_solids`].
    InputNotLive,
}

/// Where a reject was looking — the location payload that rides **beside** the reason in
/// [`BoolError::Rejected`], for a consumer that wants to point at the failure (a viewer
/// drawing a marker), never inside [`RejectReason`] itself: the reason is categorical
/// vocabulary (the census key, the thing tests name), and a coordinate is a measurement.
///
/// Two rules about what these values are:
/// - **A witness, not a census.** A reject with several offending entities carries one,
///   chosen deterministically (the minimum by the entity's own order), so identical inputs
///   yield identical errors in every build.
/// - **A diagnostic realization, cache-grade.** The coordinates are rounded `f64` world
///   positions realized at the raise site. Draw with them, report them — never judge with
///   them, and never compare them with `==` (use a distance against a tolerance; this type
///   has no `Eq` for the same reason [`Point3`] has none).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RejectWhere {
    /// A single point — e.g. the pinch vertex of a non-manifold result.
    Point(Point3),
    /// A segment between two points — e.g. the edge along which a result touches itself.
    Segment([Point3; 2]),
}

/// What kind of answer a rejection is — **the API a consumer should branch on**.
///
/// [`RejectReason`]'s variant names are engine vocabulary (plane triples, seam runs, ring
/// naming); they are stable identifiers for logs and bug reports, not something an
/// application author should have to understand. This classification is knowledge only the
/// kernel has, so it is exposed rather than left for every consumer to guess at.
///
/// Exhaustive on purpose — a consumer switching on it should be forced to decide about every
/// class, and the taxonomy is meant to stay this small (the reasons grow, not the classes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RejectClass {
    /// **The kernel's current coverage ends here** — a fact about the kernel, and nothing is
    /// judged about the input itself. Epistemically the same input *may* succeed at a later
    /// milestone (the kernel cannot rule it out from where this class is assigned), but that is
    /// a possibility, not a promise: some members will only ever earn a more precise rejection.
    /// (The 45° fold is not an example: the merge abstains and the fold reaches its
    /// `Impossible` truth, `SelfTouchingResult`.)
    ///
    /// ★ Not `NotSupportedYet`: a "Yet" reads as a prediction that
    /// support is coming, and a
    /// deliberate refusal (the DNA: reject honestly rather than guess) is not an unfinished
    /// feature. The name states the present fact.
    NotSupported,
    /// No valid solid exists for this input, at any milestone: the operands are not valid
    /// 2-manifolds, or the requested combination pinches. A statement of fact, not advice —
    /// which coincidence was unintended, and what the author does about it, is theirs
    /// (the same rule [`RejectReason::SelfTouchingResult`]'s doc states: a diagnosis, not a
    /// prompt).
    Impossible,
    /// An engine invariant broke: the arrangement built something malformed, or a backstop
    /// that should be unreachable spoke. Reported rather than returned (DNA: never silently
    /// wrong), and worth a bug report. Say "could not produce a valid result", not "your fault".
    ///
    /// The classification of variants that have never been observed to fire is provisional —
    /// tighten it once the reason census says which are reachable.
    SuspectedDefect,
}

/// Which guard raised an [`BoolError::Rejected`].
///
/// Named so a guard and the test that asserts it share one identifier — a renamed reason then
/// cannot silently drift out of a test's expectation. `#[non_exhaustive]`: reasons are added
/// and refined as coverage grows, so match with a wildcard arm and branch on [`Self::class`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum RejectReason {
    /// An **operand** has an edge used by other than two face loops, so it is not a valid
    /// 2-manifold. `validate` calls this `NonOpposedEdge` and every shell the operations build
    /// is manifold — but `boolean` never runs `validate` on its inputs, so a direct caller could
    /// still hand one in. Raised only where an operand is read (its plane table and its seam
    /// neighbours); the *result*-side closure check is [`Self::OpenResultShell`], which is a
    /// different situation. Unfired across the suite.
    NonManifoldEdge,
    /// One result solid uses an edge **more than twice**: its own surface meets itself along that
    /// line, so it would pinch there and no 2-manifold solid contains it. The edge twin of
    /// [`Self::NonManifoldVertex`].
    ///
    /// ★★ **It is about one solid, not two.** Two bodies that meet along a line are two bodies and
    /// come back as such — handles are minted per output solid, so their contact is two edges of two
    /// uses each and there is nothing here to count. What is left is the case where the material
    /// **runs around** the contact: cut the pinch and one piece remains, so no pair of solids
    /// exists to hand back and the reject is `Impossible` at any milestone. A square block with two
    /// square voids meeting along a line is the smallest example (`tests/contact_separates.rs`).
    ///
    /// Two sampled boxes sharing exactly an edge are **not** this: those are answers,
    /// and the grid proptest scores them against inclusion-exclusion instead of skipping them.
    NonManifoldResultEdge,
    /// **The result's own surface touches itself**, leaving the material no thickness where it
    /// does. An embedded boundary cannot do that — two pieces of it would occupy the same points —
    /// so what came back is not a solid, however plausible its volume.
    ///
    /// ★★★★★ **Two witnesses, one proposition.** The first is planar and combinatorial: an *edge*
    /// of the solid lies in the interior of one of that same solid's faces. The second is curved
    /// and has no edge at all: a **lateral tangent to one of the solid's own plane
    /// faces along a line**. Nothing splits at such a contact, so there is no edge to count and no
    /// vertex to link — the topology sees a perfectly ordinary solid — and it is caught instead by
    /// asking what the material near the line is: three regions (the lens inside the cylinder, the
    /// **two** wedges beside it, the far half-space), adjacent only `L–W2` and `F–W2`, so the kept
    /// ones fall apart exactly when `(W2 ∧ ¬L ∧ ¬F)` or `(L ∧ F ∧ ¬W2)` — and the pieces land in
    /// one body. See `boolean::tangency_reject`. ★ The second disjunct is *not* a defect on its
    /// own: two bodies touching along a line are two valid solids, and the kernel returns them.
    ///
    /// Distinct from [`Self::NonManifoldResultEdge`], whose proposition ("uses an edge more than
    /// twice") is **false** here: the face is not split at the contact, so every edge is still used
    /// exactly twice and the topology count sees nothing. (That sentence was written for the planar
    /// witness and is the reason the curved one needed a rule of its own — ☑ measured: a tangent
    /// `Cut` returns `Ok` with `validate` clean and every net silent.) `Impossible` for the stronger reason —
    /// that one is now about a single solid pinching *itself* along an edge — the day it named an
    /// output shape the kernel might allow arrived, and several solids is what two touching bodies
    /// get. Both survive because each states a different proposition about one solid's surface: an
    /// edge inside a face here, an edge four faces share there.
    ///
    /// What produces it is two surfaces of the model meeting exactly, leaving the material between
    /// them no thickness at all: a pocket whose wedge tip lands on the far wall, a boss flush with
    /// a neighbouring face. Parasolid refuses the same bodies
    /// (`PK_FACE_state_bad_face_face_c`); other kernels name the condition "zero thickness
    /// geometry".
    ///
    /// ★ **This says what is wrong, not what to do about it.** Which coincidence was unintended —
    /// and whether the answer is a different dimension, a different operation, or two bodies
    /// instead of one — is the author's design intent, which the kernel cannot see. Suggesting a
    /// nudge would be guessing at it, and [`crate::BoolReport`] already writes the rule this
    /// follows: *a diagnosis, not a prompt*.
    SelfTouchingResult,
    /// The assembled boundary leaves an edge used **once** — a dangling edge, so the face set is
    /// not closed. Unlike [`Self::NonManifoldResultEdge`] this says nothing bad about the input:
    /// the assembly dropped a face, which is ours to fix. Unfired in the suite.
    OpenResultShell,
    // ---- the arrangement's own consistency checks ----
    //
    // These ten questions used to share one name, `LoopOrientMismatch`, raised from **22
    // places** across three files — and its doc claimed it was unreachable, which measurement
    // refuted. One label over that many questions makes a reject name a dead end: it says the
    // engine is unhappy without saying about what, and a census cannot tell two causes apart.
    // Each name below is one question, asked at one kind of place.
    /// **Two distinct plane triples name the same point.** A vertex is named by the three planes
    /// meeting there, so two names for one point means a **fourth** plane passes through it — the
    /// arrangement then holds an edge with no direction, or a loop whose extreme vertex recurs.
    ///
    /// The same substrate limit [`Self::FourPlane`] names; this is where it surfaces on the
    /// paths that do not run that check.
    CoincidentNodes,
    /// A ring's node names do not chain: two consecutive nodes share other than exactly the
    /// class plane and one wall, or a node has no third plane to be named by.
    RingNaming,
    /// A ring is shorter than a triangle, so it bounds nothing and has no winding to read.
    DegenerateRing,
    /// **A corner with no turn**: the two edges meeting at a ring node run along one line
    /// (their wall planes' determinant vanishes, or both edges name the same wall), so the loop
    /// has no left or right there.
    ///
    /// ★ **No fixture in the suite, and the reason is worth keeping.** Its one subject was a pair of
    /// overlapping unit boxes with one turned 60° about Z, whose `Fuse` stopped in the assembly's
    /// vertex naming — a node that is a corner of one body and a straight run of the other, named
    /// by borrowing the other body's plane. Scoping the naming per solid removed the borrowing, and
    /// the pair (which only *touches* at 60°) now comes back as two bodies. The sweep was re-run
    /// over that family — three axes, seventeen angles, four overlaps, all three kinds, **612
    /// booleans, none declining**. Like [`Self::DegenerateWitness`] the check stays wired: it is the
    /// honest answer if a ring ever does run straight through a node, and a corpus is not a proof.
    StraightAngle,
    /// **Two edges leave one arrangement vertex at the same angle**, so the cyclic order around
    /// that vertex has no answer — and the face walk is built from exactly that order.
    ///
    /// ★ A line and an arc **tangent** at the vertex are not this when the arc
    /// leaves the other way (a fillet's smooth corner is a half turn, read by geometry). What is
    /// still this: the same way — two edges tangent *and* co-directed, whose order is a matter of
    /// **curvature** (two tangent circles, the shape M6b's cylinder pairs will bring).
    ///
    /// Their lines are then identical (parallel plus a shared point), which upstream is supposed
    /// to have already resolved: `merge_coincident` folds an edge traced twice, and `Aliases`
    /// folds two walls that carry one line. Reaching here means one of those did not, and the
    /// honest answer is that this arrangement cannot be ordered rather than an order picked
    /// arbitrarily — the ordering feeds `next`, so a guess there is a wrong face, silently.
    UnorderedEdges,
    /// **Two facts about one edge disagree**: a graze coincident with a same-solid crossing, two
    /// seated faces claiming opposite sides, or one solid crossing the same edge twice. Whichever
    /// is right, nothing on that edge says which.
    EdgeOccupancyConflict,
    /// **A circle carrying a contribution over only part of it was never cut.**
    ///
    /// A lateral face with a hole marks a class over an *arc*, and `split_circles` feeds both ends
    /// of every such extent into the split — so a partial contribution forces its own cut and this
    /// cannot arise. It is checked where a whole-circle mask would otherwise be built from a
    /// partial trace, because that mask is wrong in silence: the arc inside the hole would flip the
    /// bits of a face that is not there.
    PartialCircleUncut,
    /// **Neither traversal of the class's rings produced a subdivision.** Two things must hold and
    /// both are checked: exactly one outer boundary per component, and `V − E + F = 2C` — the
    /// Euler relation the walk's own cell count must satisfy.
    ///
    /// ★ The second is not a stricter reading of the first. A walk can attach the wrong edges to
    /// the right number of contours: measured, one wrong sign in the arc turn merged four cells
    /// into one eight-half-edge orbit and the contour count accepted it. Two numbers, two failure
    /// modes, one proposition — this reason's.
    RingOrientation,
    /// **A closed ring met a line an odd number of times.** Crossings of a closed curve with a
    /// line come in pairs, so the alternation `segment_meets_face` reads along that line —
    /// outside, inside, outside — has no consistent end.
    ///
    /// Raised rather than assumed away: the leftover would otherwise be read as "outside", and
    /// this judge's whole job is to notice a contact, so a wrong "outside" is a body whose
    /// surface meets itself shipped in silence. A `SuspectedDefect` because the rings it walks
    /// are the arrangement's own output — nothing a user writes can make a closed ring odd.
    RingParity,
    /// **A cell the inside/outside labels never reached.** Labels spread across shared edges, and
    /// a component sharing none is bridged by its nesting host instead — so this says the host
    /// was not found, and the cell has no side.
    UnreachedCell,
    /// Labels reached a cell two ways and disagreed, so the flip relation does not hold across
    /// the whole complex.
    LabelConflict,
    /// One face's trace on a plane class came back incomplete, so the arrangement cannot
    /// conclude — a consumer must not read "no segments" as "the plane misses the solid".
    /// `kind` says what the tracer could not do and `face` is the operand face it gave up on
    /// (an *input* face: a rejected boolean restores the live set, so the handle stays valid).
    /// `None` names a **synthetic** face — the cap a half-space clip puts on an operand — which
    /// has no handle to give.
    ///
    /// A class can decline several faces; this names the first. The full list is the audit's
    /// business, not the error's.
    TraceDeclined {
        kind: DeclineKind,
        face: Option<Handle<Face>>,
    },
    /// The result severs into two or more material solids *and* at least one enclosed void
    /// (cavity) survives. Which outer shell owns which cavity needs a shell-scoped point-in-shell
    /// test we do not have yet, so this is honestly rejected and deferred to a follow-on cell.
    /// Reachable: `Cut` a hollow part with a cut that isolates the void into one severed piece.
    /// Born with its firing test (`severed_with_cavity_is_rejected`).
    /// A surviving cavity that no material component contains — geometrically impossible for a
    /// valid boolean result (a void lies inside exactly one piece). A defensive backstop; cavity
    /// ownership is otherwise decided exactly by [`combinatorics::point_in_component`] containment.
    CavityNoOwner,
    /// No material-enclosing (outward) shell among the result components — every component is
    /// inward-oriented. Geometrically impossible for a real solid result; a defensive backstop
    /// with no firing test.
    NoOutwardShell,
    /// The assembled result has an **odd Euler characteristic** (`V − E + F − L_i`), which no
    /// closed 2-manifold can have (it must equal the even `2(S − G)`) — so the arrangement produced
    /// a malformed solid and the boolean rejects rather than return it (DNA: never silently wrong).
    /// This is the Euler-parity backstop for malformity that is *not* a pinch (see
    /// `NonManifoldVertex`); e.g. a dropped face. Rotation-independent. Checked post-assembly in
    /// `boolean`, per solid.
    EulerParity,
    /// One result solid has a **non-manifold vertex** — a "pinch" where two or more of its face-fans
    /// meet at one point (a cutter's convex corner exactly on the target's concave corner), even
    /// though every edge of it is manifold. No valid 2-manifold solid has one, so the boolean
    /// rejects with this clear reason rather than the incidental `EulerParity` (which also misses an
    /// *even* number of pinches). Rotation-independent; checked per solid post-assembly via
    /// `nacre_topo::nonmanifold_vertices`.
    ///
    /// ★ **Two solids touching only at a corner are not this.** They are two solids, and each gets
    /// its own copy of the point; what remains is a body whose material runs around the contact, so
    /// there is no pair to part into.
    NonManifoldVertex,
    /// The assembled result has an even Euler characteristic but a **negative genus** (`S − χ/2 < 0`)
    /// — more handles than a solid can have, so it is not a valid closed 2-manifold. A count-based
    /// backstop below the pinch and parity checks; checked per solid post-assembly.
    NegativeGenus,
    /// Three planes that should meet in a point do not (a parallel pair), so an arrangement
    /// vertex has no name.
    ThreePlanes,
    /// Two **different** arrangement vertices (distinct plane triples) materialized to the same
    /// coordinate. The triple is the truth and the coordinate only its cache, so this
    /// says the exact substrate and the f64 cache disagree about how many vertices exist — always a
    /// defect upstream, never a property of the input. Raised where the seam table is built, while
    /// both triples are still in hand; without it the disagreement surfaces much later as a
    /// zero-length edge. The known cause is a **split plane table** (one geometric plane carried by
    /// two classes).
    ///
    /// ★ A genuine 4-plane concurrency used to reach here too, and that was **this reject naming
    /// the wrong thing** — a property of the input reported at the class that says "report a bug".
    /// The arrangement folds those names now, including the pierce that was missing: a plane which
    /// *carries* an arrangement line rather than crossing it (`Aliases::wall_family`, and
    /// `tests/concurrent_line.rs` for the shape). What remains here is meant to be a real defect.
    SeamAlias,
    /// A result loop asked for an edge between two vertices at the same coordinate. Every ring node
    /// is a distinct arrangement vertex, so this cannot happen for well-named input — it is the
    /// backstop that keeps a degenerate one from aborting the kernel (`Line::through_points` used to
    /// `expect`). `SeamAlias` catches the known cause earlier, so this has no firing test.
    ZeroLengthEdge,
    /// Four planes concurrent at one point: two distinct plane triples name the same arrangement
    /// vertex, which the substrate cannot express.
    ///
    /// **What builds one.** A tool *edge* lying inside one of the target's planes — tangential
    /// contact rather than a crossing. Three planes then share that edge's line (the tool's two
    /// faces and the target's one), so every point of the line already lies on three planes and a
    /// fourth turns it into a vertex: each plane crossing the line yields a four-plane point.
    ///
    /// **★ The model this was written for now builds, and this has no firing test.** That model —
    /// a bar spun **45°** about an axis in the plane, whose bottom corner edge lands back in it
    /// because half-width equals the pivot-to-bottom offset — is `nacre-oracle`'s
    /// `the_four_plane_cut_matches_occt`, where OCCT scores the result and agrees on volume and
    /// centroid. Vertex identity was normalised (names identify, structures carry geometry), and
    /// `coverage/rotation.rs` keeps the near misses around it building.
    ///
    /// What survives is the guard, at one site: a run whose candidate handles are all parallel to
    /// the line they must cut (`arrangement`). It has no reproduction. It stays because a
    /// substrate that cannot name a point must say so rather than pick one of the names — and
    /// because an unfired reject costs nothing, while a missing one costs a wrong solid.
    ///
    /// ★ There is no **operand**-side guard — a vertex found on more than three plane classes
    /// (a gusset's apex landing on a wall's top edge) is not refused. A point's name is
    /// one function (`combinatorics::canonical_triple`), and an operand vertex names itself
    /// from the classes its topology knows.
    ///
    /// **Exact, not toleranced.** Measured on that model: perturbing an operand coordinate by
    /// **one ULP** in either direction removed the concurrency. The judgement returned zero
    /// because the determinant *is* zero at 200 bits, not because it fell below the coincidence
    /// limit (~55 orders of magnitude lower). Exact rational input would produce the same
    /// concurrency — this was never an artefact of `f64` construction, which is why exact
    /// rational construction (`crate::exact`) left it exactly where it was.
    FourPlane,
    /// A ring holds a vertex the arrangement names as a `plane ∩ plane ∩ cylinder` **pierce
    /// point**, on a path that speaks only three-plane names — the ring walks, the wall-and-handle
    /// derivation, the seam table. The point is exactly named; what is missing is that these paths
    /// have no other name to carry it by, and dropping it from a ring would silently answer about a
    /// different polygon.
    ///
    /// ★ **Not [`Self::RingNaming`].** That one's sentence is "these names do not chain" — a fact
    /// about a *three-plane* naming that came out degenerate. A pierce point has no such name to
    /// begin with, so reporting it there would put two causes under one label.
    ///
    /// ★ Its per-face sibling is `DeclineKind::PierceNode`: the tracer's decline structure carries
    /// a face handle, so the paths inside it say *where* instead of raising this.
    ///
    /// ★★ **A ring may now hold a pierce node, and this no longer refuses one.**
    /// `loop_winding` compares pierce coordinates through the quad tower
    /// (`combinatorics::CoordKey`); what is left here is the narrower sentence — a pierce node
    /// whose *cylinder* is not on any of the ring's own carriers, so the point cannot be re-solved
    /// from its name. That arm has no fixture, and neither does the ray-cast one next to it; both
    /// are recorded rather than hidden. (`DegenerateWitness` is the precedent: no fixture, written
    /// down, and not taken as licence to delete the guard.)
    PierceVertexUnnamed,
    /// A **result** vertex's definition names a surface the finished solid keeps no face on, so
    /// the point could not be re-solved from the solid's own geometry — the property transform
    /// and replay stand on (`defs_are_remappable`). Measured population: the
    /// **contact-cut** — a plate cut by a boss that only touches its top face. The cut
    /// removes no material, but the arrangement minted pierce vertices where the boss's rim
    /// crossed the plate's edge, and the assembly kept them with definitions still saying
    /// `wall ∩ boss cylinder` after every boss face was gone. The volume was already right; the
    /// names were not. A shipped check rather than a debug_assert, because this population
    /// reaches it — until the assembly learns to shed the stale corners, refusing is the floor.
    ///
    /// ★ **Two bodies meeting on a full wall are not this**: the wall is interior, so the
    /// result keeps no face on it, but the coplanar merge carries a component's circle holes,
    /// the corners that sat on the wall dissolve, and that
    /// population builds. What is left under this name is the contact-cut: a definition naming a
    /// surface the result has **no face on at all**, which no amount of merging repairs.
    VertexNamesAbsentSurface,
    /// A plane class neither perpendicular nor parallel to a cylinder's axis whose lateral faces
    /// could not be shown to miss the plane — where they meet, the intersection is an ellipse
    /// (vocabulary the kernel does not have yet). ★ The gate asks the faces first: a slanted
    /// plane that
    /// runs past every lateral face of the cylinder (a gusset beside a plate's holes) passes.
    ObliqueCylinderCut,
    /// **A cell called a circle sits on a class that is not perpendicular to its axis** — so its
    /// boundary is an ellipse and every rule that reads a *radius* as that boundary's width is
    /// wrong there.
    ///
    /// ★ **A backstop, deliberately, and it is classified as one.** [`Self::ObliqueCylinderCut`]
    /// says *the kernel's coverage ends here* and judges nothing about the input; this says
    /// something stronger and different — the gate promised no oblique class would carry a circle
    /// (see [`crate::combinatorics::class_carries_circle`]) and the promise did not hold. Sharing
    /// the older name would have put two sentences under one label, which is the shape this
    /// kernel keeps re-learning. The vessel-guard precedent is [`Self::PartialCircleUncut`], and
    /// it takes the same class.
    ///
    /// ☑ Nothing raises it today, and that is the point: the alternative at its two sites is not a
    /// refusal but a **silently wrong containment answer** — a rim witness off the cell's own
    /// plane, or two radii compared as if an ellipse's width were its radius.
    ObliqueCircleClass,
    /// Two cylinder classes **of different operands may share a face**: one surface stated by
    /// both (one handle with rows of both solids, or one surface under two handles — a
    /// `translate`d twin), or axes within the radius sum whose lateral faces could not be shown
    /// to miss each other along **any** direction the pair can state (`planes::separating_dirs`:
    /// each axis, and for skew axes the common perpendicular). Where the faces do meet, their
    /// intersection is a quartic curve, which is M6b's.
    ///
    /// ★ What the check must deliver for the name to be true, four ways. The distance is
    /// measured for skew axes as well as *parallel* ones ([`nacre_scalar::cylinders_clear`]), so
    /// a drill crossing a bore with room to spare passes. The distance is not the whole rule —
    /// it is a fact about two infinite
    /// surfaces: a stud fused through a cube and a second stud across it have crossing
    /// axes, though no *face* of one reaches the other, so non-parallel pairs ask the faces.
    /// Pairs **within one solid** (a plate's two fillets, a slot's two half cylinders) are not
    /// asked, parallel pairs read their
    /// spans, and the coincident surface under two handles is the one parallel refusal left.
    /// ★★ And the face test asks more than the two **axes**: a filleted plate and a drill
    /// laid across it are apart along the rulings' cross product and nothing else
    /// (the rational candidate set is closed).
    CylinderPairContact,
    /// The population gate could not decide a (plane, cylinder) pair **exactly** — a plane
    /// with no narrow rational description, a rotated class, a moved cylinder (its def is
    /// pre-motion), or checked-`Rat` overflow. Conservative honest refusal, never a guess.
    ///
    /// ★ The chart's emitter raises it too, for the
    /// `chamber` sentence: a lateral cell whose two ends **disagree** about its chamber (☑ the
    /// (0,0)-corner boss, where the plate classes' disk cells carry no B material while the
    /// caps do — an arrangement label defect refused here instead of assembling an open shell),
    /// and a cylinder class with no row or rows of both solids. The cell reader
    /// raises it for a cell with a face whose **two ends both say nothing** — answering
    /// «present» there would be a guess; measured 0 with stations placed by name.
    CylinderGateUndecided,
    /// **The trace does not determine whether a lateral face is present over a sector.**
    ///
    /// A cell label says where *material* is; whether the cylinder's own face bounds it there is a
    /// different proposition, answered by the contributions covering the rim arc (`bands`'
    /// `face_spans`). This name is what that answer being *contradictory* is called, and it has
    /// two spellings — both of them "the evidence points both ways", which is why they share one
    /// name rather than reporting the shape that produced it:
    ///
    /// - The interval's **two rims disagree.** Every cut rim is a band boundary, so a sector
    ///   exists over an interval's whole height or over none of it; ends that disagree mean that
    ///   is false here, and choosing an end to believe is a guess.
    /// - **Two faces of one solid** cover the same arc and disagree about reaching in. That takes
    ///   one cylinder class holding several faces which share a rim; the eventual answer is
    ///   likely "any of them is enough", but it has never been measured.
    ///
    /// ★ Not [`Self::CylinderGateUndecided`]: nothing here failed to be *computed*. The
    /// arrangement decided every piece exactly and the pieces contradict each other.
    CylinderFaceUndecided,
    /// An operand face has no three non-collinear outer-loop points, so it spans no plane.
    /// (Its sibling `DegenerateNormal` — a zero-length triangle normal — died when `n_out`
    /// moved to the stored orientation: no triangle cross is taken, so there is nothing
    /// left to degenerate.)
    DegenerateFace,
    /// A plane's own frame cannot be stated exactly, so a sketch built in it has no exact
    /// definition to judge from. Either the plane recorded no rational coefficients, or the
    /// squared lengths the frame's realization divides by do not fit `i128`.
    ///
    /// ★ **Distinct from the retired `InexactSurface` on purpose** (both ended with "no exact
    /// definition here", but they say different things about *why*: that one names a surface with
    /// no provenance, this one names a frame the rationals cannot hold. Sharing a label would
    /// leave a firing pointing at the wrong producer.
    FrameOutOfRange,
    /// **The model's rotation history is longer than the judging budget.**
    ///
    /// A rotated point's realization carries an error `C · 2⁻ᵖʳᵉᶜ`, and `C` grows about one bit
    /// per turn (measured). The kernel therefore sizes `prec` from the model, and this says that
    /// size exceeded [`crate::planes::JUDGE_PREC_CAP`]: `needed` bits to separate what the
    /// operation must separate, against a cap of `cap`.
    ///
    /// **Nothing is wrong with the model** — this is a cost limit, not a resolution one. The same
    /// solid turned fewer times builds, and raising the cap would build this one too, slowly.
    /// It exists so the cause has a name: at a fixed 256 bits a solid turned 245 times used to
    /// fail as `LoopOrientMismatch`, a symptom three layers away from the reason.
    PrecisionBudget { needed: usize, cap: usize },
    /// **A judgement ran out of bits.** A determinant could not be separated from zero even at
    /// the judging cap, and the separation it *did* bound is wider than the coincidence limit —
    /// so calling the two things one would be a guess, and the kernel says so instead.
    ///
    /// Sibling of [`Self::PrecisionBudget`], which is the same shortage seen before any work
    /// starts: that one is the whole model being too deep, this one is a single judgement being
    /// harder than the model's own depth suggested (a very thin witness). More bits would answer
    /// it. Nothing is wrong with the model.
    JudgeExhausted,
    /// **A judgement had no distance to measure.** Turning a determinant into a length means
    /// dividing by a cofactor, and this one came out exactly zero: a witness triangle that has
    /// collapsed to a line, or three planes with no meeting point to speak of.
    ///
    /// **A different cause from [`Self::JudgeExhausted`], and the difference is what to do about
    /// it**: no amount of precision creates a distance that is not there (measured — the
    /// determinant was still bit-exactly zero 8192 bits deeper). The arrangement asked for a
    /// point that does not exist.
    ///
    /// ★ **Unfired across the suite since the outwardness label became nesting parity.** The one
    /// judgement the corpus ever came back degenerate on was the component material/void test,
    /// which asked for the sign of a rotated plane's normal; the parity label asks containment
    /// instead, which the substrate always answers. A sweep over three rotation axes, ten angles,
    /// two operand shapes, four overlaps and all three kinds (720 booleans) found none, and none
    /// of `JudgeExhausted` either. `undecided_reject` is still wired and still the right shape —
    /// an undecided judgement must never reach the geometry as a zero — it simply has no fixture,
    /// which is recorded rather than taken as licence to delete the guard.
    ///
    /// ★ What it stands between: a triple of three planes sharing a line (a four-plane operand
    /// vertex named per face loop would produce one) makes the judge, asked that «point»
    /// for a side, come back degenerate, and everything downstream would read the `0` as «on every
    /// class» — false alias sets, parallel planes taken for a shared line, one merged wall family,
    /// every name folded onto a corner elsewhere. This reject is the only thing between that
    /// and a solid with a vertex in the wrong place. No producer mints such a dependent name.
    DegenerateWitness,
    /// Every candidate ray from a loop's nodes has a ring node on its line — or, in 3D, every node
    /// of a component grazes the boundary it is being classified against.
    ///
    /// `point_in_ring` casts along `P ∩ Q_a` for a node's own plane `Q_a`; a ring node on
    /// that line makes the crossing parity ambiguous. Candidates are `2 · |loop|` lines and
    /// two directions, and half of them can be spoiled at once — `l_and_staple`'s loop and
    /// arc share both `y` planes, so only the `x` lines are clear there.
    ///
    /// ★ **The population that can reach it.** A component's
    /// nesting depth is decided by `point_in_component` from one of that component's own nodes.
    /// Contacts separate, so a touching void and its host are two components, and a void whose
    /// *every* corner sits on a wall — a
    /// diamond inscribed in a square block — leaves no candidate that does not graze. The retry
    /// over the other nodes is what keeps ordinary contacts (a cube corner: one node of eight)
    /// clear. Locked in `tests/contact_separates.rs`.
    ///
    /// ★★ **Raised only where the retries actually run out.** A single node
    /// failing to decide is an *abstention*, not an error — `point_in_component` says so in its
    /// type (`Ok(None)`); saying it with this reason instead measured 122 of the whole suite's
    /// 151 raises as that guard being caught and swallowed by
    /// its own retry loop. What raises this reason today is the caller whose node supply is
    /// exhausted (the 3D depth/cavity classification in `boolean.rs`) — and, same shape one
    /// dimension down, `point_in_ring`'s rayless case, the nesting engine's exhausted offer, and the
    /// coplanar cleaning pass's own mirror of that road (`boolean.rs`, which names the empty
    /// list [`Self::RingHasNoWitness`] and the exhausted one this, as the arrangement's
    /// `cell_in_cell` does; neither has a population there).
    ///
    /// ☑ **What keeps its 3D population empty.** A wall boss's Common parted by a slab is two
    /// components of half-disc caps, a panel and a lateral, with no vertex among them: reading a
    /// cut cap as the
    /// disk it is (`combinatorics::face_circle`) gives every one of them a witness. The
    /// **corner** families, whose axis stands on the plate's corner edge, have probes that meet a
    /// **corner on the ray**: the half-open rule is spelled once and read by the mixed ring
    /// parity as well as the planar roads, so their Fuse and Cut decide. In the corner
    /// Commons parted by a slab,
    /// measured per attempt, two of each probe's three rays meet the other half's **lateral**;
    /// a lateral's boundary loops answer whether the crossing is on it on
    /// the cylinder's own chart, and the crossing census raises this
    /// reason nowhere. What can still exhaust a probe list is the probe on a ring, a tangent
    /// ray, or a seam-incident root with an arc above it.
    ///
    /// ★★★★★ **Exhaustion is not emptiness.** A probe list that starts **empty** — a
    /// ring cut out of a cylinder names its corners with the quadric and `three_plane_probes`
    /// keeps only plane triples — is not a list of rays that were cast and blocked.
    /// That fact has its own name ([`Self::RingHasNoWitness`]).
    ///
    /// ★ For the nesting question the two are told apart by one rule and not by which of
    /// four spellings a caller happened to hold — see [`Self::RingHasNoWitness`].
    ///
    /// ☑ **A component is offered more than its corners.** A rounded outline under a quarter
    /// turn (six of the
    /// motion group's rotations) and a void whose every corner rides a wall are the supply
    /// being thin, not the geometry being hard: given the points its **edges** name
    /// ([`crate::combinatorics::edge_interior_points`]) every one of those questions decides. The
    /// six rotations **build**, with one body, a clean `validate` and the exact volume; the void
    /// gets the name its shape has one grazing corner down
    /// ([`Self::SelfTouchingResult`]). What is left under this name raises **0 or 1 times per
    /// whole-suite sweep** — a backstop with no population, and no fixture surfaces it.
    NoClearRay,
    /// **A ring — or a component — offered no point to ask about** — not a ray that was blocked,
    /// and not a value that could not be formed: the containment roads draw their witnesses from
    /// a ring's *corners*, and a ring cut out of a cylinder has only pierce-named ones; the 3D
    /// depth classification draws them from a component's vertices and, failing those, from its
    /// cut caps' interiors, and a thin segment of a disk can offer neither.
    ///
    /// The three names beside each other, once: [`Self::NoClearRay`] is «every witness we had was
    /// blocked», [`Self::WitnessNotRational`] is «a value we needed could not be formed
    /// exactly», and this is «there was no witness to begin with». Reading the first for the
    /// third sends a diagnosis to the wrong layer.
    ///
    /// ☑ **Nothing in today's corpus raises it.** ★ The three refusals of the nesting
    /// question follow one rule ([`crate::nesting::cell_inside`]): a value that could not be
    /// **formed** is
    /// [`Self::WitnessNotRational`], an offer that was **empty** is this one, and an offer every
    /// member of which **abstained** is [`Self::NoClearRay`]. So this name means exactly what
    /// it says — the cell had nothing to offer — and there is one supply. On the **component**
    /// road, an offset boss's Common is a 0.2-deep segment prism whose halves
    /// have pierce-named corners only and no cap candidate from the centre inside — so the
    /// cut cap offers two points per **chord** (`ring_interior_candidates`). What would still
    /// reach it is a segment cut again along its chord's normal line, in a multi-body result.
    /// The remedy for every such shape at once is to widen the probe's **type** so a pierce
    /// corner is itself a witness.
    ///
    /// ★★ **A disk answers the converse through its rim.**
    /// [`crate::nesting`]'s reverse pass skips an interior witness without setting a flag, so a
    /// `Cell::Disk` whose whole supply is its centre would offer nothing and fall straight to
    /// this name, every time. Nothing in the corpus walks that road (the `== 0` above is an
    /// emptiness of *this* corpus); the disk's rim witnesses close
    /// it, and `a_disk_inside_a_disk_is_decided_by_its_rim` is that road's only coverage.
    RingHasNoWitness,
    /// **An exact *value* could not be formed** — a class with no narrow rational description (a
    /// rotated one, say) or a coordinate past `Rat`'s ceiling.
    ///
    /// Raised where the cylinder work needs a number rather than a sign: a plane's axis parameter
    /// against a cylinder (`planes::axis_param_of_plane`, read by the band pass and the
    /// transversal-circle test), the rational chart the circle nesting projects a ring into
    /// (`nesting::cell_inside`'s coordinate road), the order of two points on a meet line when one of
    /// them has no exact description (`combinatorics::order_pinned` and the interval overlay that
    /// calls it), and the arrangement's split passes, where a point's own description — its
    /// `(line, s)` or its `dir_sign` — could not be formed.
    ///
    /// ★ **It used to name the arc split's demand for both endpoints' coordinates, and that is
    /// gone.** The arc split asked for a coordinate where it wanted an order; it asks for the order
    /// now, so what is left under this name there is the honest cause — a description past `Rat` —
    /// and not a shape the road cannot spell. A chained cylinder operand no longer stops here.
    ///
    /// ★ **Distinct from the gate's own name on purpose.** The population gate asks *signs*, and
    /// those were made total (they clear denominators and answer in `BigInt`), so
    /// [`Self::CylinderGateUndecided`] means the geometry — a rotated class, a moved cylinder, a
    /// surface contact. This one means the arithmetic: the wall a wide model meets after the gate
    /// has already said yes. Sharing one name would put a width limit inside a geometric verdict.
    ///
    /// ☑ **The chart half of that wall is closed**: the chart reads a **primitive** normal, so
    /// a plain sub-millimetre model does not meet the width of its plane offset's denominator
    /// squared by the chart's second axis. What is left under this name on that road is a
    /// **primitive**
    /// normal past ~2⁶³, and the corpus has none. ⚠ It has raises but no frozen
    /// *surfacing* fixture; `reject-trace` is what keeps it visible.
    WitnessNotRational,
    /// **An arc-bounded boundary this assembly cannot spell yet** — a backstop, no longer a
    /// stopper.
    ///
    /// The arc population builds (the straddling boss — `bands`' fences), so no stopper
    /// carries this name. What keeps the name alive are its honest backstops:
    ///
    /// * a **band rim spelled `Rim::Circle` on a cut circle** — the emitter spells a cut
    ///   rim once, as the chain of its arcs, so `band_loop` reaching a whole-circle rim that the
    ///   split cut is a producer inconsistency, spelled there rather than assumed away;
    /// * the region walk's own two: a run along an **uncut** rim that is not the whole circle,
    ///   and cycles whose winding does not classify into one lower and one upper rim with holes
    ///   the assembly can bridge (`classify_cycles`' abstentions — a chain with more than two
    ///   seam contacts, a hole meeting the seam at other than zero or two, two bridging holes;
    ///   ☑ all measured 0 across the suite);
    /// * a **whole-disk bound on a cut circle** (`circle_loop`) — a producer inconsistency (the
    ///   trace subdivides a cut disk into cells), named honestly rather than as a dropped
    ///   crossing;
    /// * a **hole meeting the seam at other than two contacts** — the merge that produces such a
    ///   hole abstains on the same count, so this is a producer inconsistency too;
    /// * a **contact whose cut circle has no world description** — `band_loop` orders the two
    ///   contacts by their circles' exact axial parameters, and the population gate demanded a
    ///   world description of every plane class long before a band could be emitted, so this too
    ///   is spelled rather than assumed away;
    /// * a **chain rim with no seam contact** — a chain winds once about the axis, so it
    ///   passes the seam somewhere and the split named that point; a chain the assembly cannot
    ///   attach its slit to is a producer inconsistency, named.
    ///
    /// "Not yet" is still the literal truth for all five.
    ArcBoundNotYet,
    /// **An *operand* face is bounded by a cylinder, and the tracer names rings by planes.**
    ///
    /// The tracer reads each operand face's loops as three-plane triples with a carried wall
    /// class beside each edge (`combinatorics::loop_triples`). Both of those are plane-only,
    /// and both are total until a boolean's *result* is fed back in as an operand: a boss standing
    /// on a wall leaves the plate's caps bitten by an arc and the wall split by two rulings, so
    /// those faces' rings run along the boss's lateral surface. There is exactly one such loop the
    /// road already speaks — a **full circle**, one rim edge whose far face is the cylinder, which
    /// is named by that cylinder's class — and everything else lands here.
    ///
    /// ★ **Distinct from [`Self::ArcBoundNotYet`] and [`Self::RulingBoundNotYet`], which are the
    /// *assembly's*.** Those two name configurations the result side cannot **spell**; this one
    /// names an operand the trace cannot **read**. Sharing either would put an input-side coverage
    /// limit under an output-side name, and the census keys on these.
    ///
    /// ★★★★★ **Two sites, and neither of them fires today** — measured over the whole suite
    /// (`--features reject-trace`: 86 raises across 15 reasons, this one **absent**). What a caller
    /// meets on a chained operand has kept moving outward as the rungs went in, and the sentence
    /// here has had to move with it four times: [`DeclineKind::PierceNode`] at the tracer's ring
    /// naming, then [`Self::PierceVertexUnnamed`] at the arrangement's then-plane-only overlay,
    /// then [`Self::WitnessNotRational`] at the arc split, which asked every segment for its two
    /// ends' exact coordinates — and now [`Self::CurvedStraightRun`], **one stage further into the
    /// arrangement**: its split passes all take a cylinder-pinned end, so a chained operand reaches
    /// the *cell* stage, where `arrangement::walk_cells` asks `combinatorics::loop_winding` for a
    /// ring's winding and its extreme-node walk steps past two tangent arcs of one circle.
    /// ☑ Measured by backtrace, not inferred: the raise arrives through `walk_cells`/`per_class`.
    ///
    /// * `planes::face_clears_footprint` — the population **gate**. It used to refuse the whole
    ///   chained population here, because its clearance scan stopped at the first corner a
    ///   cylinder made. That corner is now **read exactly** (a pierce point carries one radicand
    ///   and everything it is measured against is rational), so what is left under this name is
    ///   only what the scan genuinely cannot spell: an **arc** edge — whose bulge breaks the
    ///   convex-hull argument the scan rests on — and an `OnSeam` vertex, which pins a curve
    ///   rather than a point. The gate's *arithmetic* failures are not this name: they wear
    ///   [`Self::CylinderGateUndecided`], whose sentence is true of them and false of these.
    /// * `combinatorics::loop_triples` — the road itself, **swallowed**: `trace_input` maps a
    ///   loop it cannot name to `None` and the tracer reports which *loop* failed
    ///   ([`DeclineKind::OuterRing`] / [`DeclineKind::HoleRing`]), which is the finer fact when a
    ///   face has several. The reject census still records the raise, so the cause is in the ledger
    ///   even where the surfaced label is the loop's.
    ///
    ///   ★★ **That second site is now a backstop with no firings, and deliberately so.** Once the
    ///   ring naming learned to restate an operand's pierce corner
    ///   (`combinatorics::pierce_name_from_def`) and to write a curved carrier, the curved rings of
    ///   today's population are **described** rather than declined — measured, zero raises from
    ///   this site across the whole suite. What can still reach it: a corner whose def is not
    ///   `Pierce` (a ruling ending at a seam vertex — a seam *joint* between two legs of one
    ///   arc is read as one step, not as a corner), a plane with no *narrow* world
    ///   description (a wide or rotated chain), a handle that answers to both candidate classes
    ///   or to neither, two laterals meeting at one corner (M6b's), and a loop whose every joint
    ///   is a seam. ★ Measured with the restatement switched off: the ring
    ///   **declines by this name** and nothing panics, which is the whole reason the backstop is
    ///   under the road rather than trusted away.
    ///
    /// ★ The two ask the same question of a face in two vocabularies — the gate reads the edge's
    /// **carriers** from the model, the road reads the far face's **class** — and they were
    /// measured to pick the same faces (a face-by-face probe over wall, corner, bore and top-boss
    /// operands). The road's side of that is now
    /// `an_operand_bounded_by_a_cylinder_is_named_in_class_space`, which locks the naming rather
    /// than the decline.
    CurvedOperandBoundary,
    /// **A circle and a ruling of one plane class cross** — two cylinders' curves crossing on one
    /// plane, at a point this kernel cannot name.
    ///
    /// ★ **Crossing, not merely meeting.** A tangency divides nothing — as for a
    /// tangency at a vertex — so a curve that only touches another mints no
    /// node and needs none. The name says so because reading it the other way refuses
    /// valid parts.
    ///
    /// A class may carry a circle (⊥ one cylinder's axis) and rulings (∥ another's) at once, and
    /// that alone is fine: neither split cuts
    /// the other's kind, so the two populations never interfere — **while they do not meet**. A
    /// crossing they both walk past is a node no road mints, and the point itself is
    /// `plane ∩ cylinder ∩ cylinder`, a nested radical with no name.
    ///
    /// ★★ **The argument that this population cannot arrive**: a crossing means two lateral
    /// *faces* share a point, which `planes::lateral_faces_clear` denies. It holds because the
    /// question is asked of the **arc**: a circle does not become an edge **whole** while the
    /// face on that cylinder uses
    /// only a quarter of it — what a contribution states is the extent its face covers, and
    /// *that* arc is the edge.
    ///
    /// ★★★ **So this is a backstop, and deliberately so.** A ruling sits nearer its cylinder's
    /// axis than that cylinder's own radius (`√(r² − h²) < r`), so an arc reaching one is inside
    /// the drill's reach and the pair rule refuses first — measured on a plate whose corner
    /// fillets grow until their arcs do reach (`f = 6, 7, 8`: all `CylinderPairContact`). It
    /// ships rather than living under `debug_assertions` because a node no road mints is a
    /// **silently wrong solid**, not a crash, and the argument above is easy to read wrong.
    ///
    /// Not [`RejectReason::CurvedOperandBoundary`] (the road behind can read these rings) and not
    /// [`RejectReason::CylinderPairContact`] (the two lateral faces do **not** touch).
    CircleCrossesRuling,
    /// The **rulings ladder's** own refusal — a configuration its machinery does not arrange
    /// yet. The assembly's edge road is open (ruling edges mint with their own key and
    /// carriers) and so is the gate's record-and-pass arm, so this name is
    /// **reachable from production**. ★ A boss whose **cap sits
    /// inside the other body's material** (the half-height variant) builds,
    /// and so does the crossing census's whole column of 15 — a lateral region
    /// that crosses a z-line transversally in one sector while ending on it in another
    /// (the emitter walks regions of the chart; `cyl_chart::regions`).
    ///
    /// What raises it today is a **producer inconsistency**, never a shape: in the region walk,
    /// a boundary run along a rim whose end is not one of the rim's nodes, a run along a ruling
    /// the wall class's pieces do not tile, a cycle whose pieces do not chain end to end, or a
    /// component with no outer cycle (☑ all measured 0 across the suite);
    /// and the assembly's standing guards — a segment lying *on* the lateral, and ruling end
    /// names that share no single plane. [`Self::ArcBoundNotYet`]'s straight sibling.
    ///
    /// ★ **«A class carrying both circles and rulings» is not one of those guards**:
    /// it would refuse a population it has no reason to, since the two only conflict when they
    /// *meet*
    /// and a meeting means two lateral faces share a point — which the cylinder-pair gate
    /// refuses first, by `CylinderPairContact`.
    RulingBoundNotYet,
    /// **A loop's winding had to be read across a *curved* straight stretch.**
    ///
    /// `loop_winding` reads the turn at the ring's extreme node, walking back past nodes the loop
    /// runs straight through. The direction it carries out of that walk is the far edge's direction
    /// **at its own start**, which is the extreme node's arriving direction only because a line's
    /// tangent does not change along it. Two arcs of one circle are tangent-continuous too — so the
    /// walk steps past them just the same — but their tangents differ, and the winding read there
    /// would be the turn at some *other* point of the ring, confidently wrong.
    ///
    /// ★★★★★ **Not an unfired guard.** Where every extreme node is a
    /// real corner — a disk cut by chords is convex, a polygon minus disks turns at every arc end —
    /// the walk never takes a step. The arc split also cuts circles on a population that
    /// produces the other shape (a **chained** cylinder operand), and all four chained fixtures
    /// stop here. ☑ Measured.
    ///
    /// ★ The answer it named still stands, and it is not a wider version of this walk: read the
    /// winding at the extremum of the **region**, which may lie in an arc's interior.
    CurvedStraightRun,
    /// A loop's node lies *on* the ring it is being tested against.
    ///
    /// A hole ring never touches the outer ring it sits in, and `point_in_ring` checks that
    /// exactly: the ray's line meets an edge at `X`, and `X == v` strictly inside that edge means
    /// `v` is on the ring. Unfired.
    PointOnRing,
    /// The trace arrangement on one plane class nested a hole whose containment depth exceeds one.
    /// `nest_cells` resolves any number of holes at depth one inside one outer loop; deeper nesting
    /// is honestly rejected until the general nesting cell lands. Distinct from [`Self::HoleRoots`]
    /// so a refactor cannot silently merge the conditions.
    HoleDepth,
    /// A plane class's arrangement produced no unbounded contour, or more than one. A closed figure
    /// always has an outside, so "none" is impossible; "several" means several disjoint bodies on
    /// the plane, which the nesting resolver does not cover yet.
    HoleRoots,
    /// A face whose boundary never crosses the seam, yet the seam lies on its plane — the
    /// convex path only.
    MissingSeam,
    /// A `Whole`-survival contact face whose footprint OVERLAPS the other's (∂P × ∂Q cross) rather
    /// than nesting, in the one such case still unbuilt. `Whole` has two entries: `Fuse`/same-normal,
    /// which the same-ground union builds, and `Cut`/opposite-normal, which is exact whenever the
    /// contact plane separates the two solids (nothing to remove). What is left is a `Cut` whose tool
    /// reaches back across that plane — a pin below its own contact face — where the cut owes a notch
    /// this path cannot yet cut. Honest reject rather than a whole cap that ignores the pin.
    ///
    /// ★ **Unfired across the suite (census).** Coplanar pieces meeting at a point during the
    /// merge are an **abstention**, not this (the merge emits a pinching
    /// group unmerged and the whole-result judgement, hoisted before minting, names the shape —
    /// `SelfTouchingResult` on every input that reaches it). What remains under this name are
    /// the nine defensive guards of the coplanar merge, none with a known input.
    CoplanarMerge,
}

/// What a face's trace on one plane class could not do — the detail behind
/// [`RejectReason::TraceDeclined`].
///
/// These name arrangement steps, not user-facing situations; branch on
/// [`RejectReason::class`] and keep these for logs and bug reports. They all mean the same
/// thing to a caller ("this configuration is beyond the tracer"), and they are kept apart so a
/// refactor cannot silently merge two different degeneracies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum DeclineKind {
    /// A ring vertex's plane triple collapses (two of its planes coincide), so it names no point.
    CollapsedTriple,
    /// A ring vertex is a `plane ∩ plane ∩ cylinder` **pierce point**, and something on the
    /// tracer's road could not take it.
    ///
    /// ★★★★★ **It has narrowed to nothing, and that is the shape of two rungs.** Last cell this
    /// read "four places could raise it and exactly one does, 20 raises — the seated road's
    /// `emit_ring`". That road stopped deriving an edge's carrier from its corners' plane sets and
    /// the raise went with it, leaving **three** sites: a walk with no exact side, a run whose body
    /// it cannot place, an order it cannot form. ☑ Measured **0** across the workspace suite and
    /// the ignored sweep — the name is now entirely unexercised, kept because each of those three
    /// can still say it.
    ///
    /// ★ Distinct from [`Self::CollapsedTriple`] on purpose: that one's sentence is "two of its
    /// planes coincide", which is simply not what happened here. The census keys its `detail` on
    /// this name, so reusing the other would put a false cause in the ledger.
    PierceNode,
    /// The scan's crossing arm met an edge riding a cylinder that it cannot name a point on: an
    /// **arc** (its crossing waits on the lateral's ruling sweep — see
    /// `arrangement::crossing_on_ruling`), or a **ruling** whose crossing has no exact statement
    /// (a class that is not ⊥ to the axis, a plane not through it, a tangency, both roots on one
    /// side). A crossing on a ruling with a statement is named as a pierce node and passes.
    ///
    /// ★ Distinct from [`Self::PierceNode`], which is about a *corner*. The two travel together on
    /// the population that produced them (a ruling's ends lie on the cylinder, so they are pierce
    /// points), but they are different sentences, and a ring whose names are perfectly good while a
    /// **carrier** is curved is a fact worth seeing on its own.
    CurvedRingWall,
    /// The face's outer ring could not be named as plane triples.
    OuterRing,
    /// One of the face's hole rings could not be named. Not "no hole": swallowing it would trace
    /// the face as if it were solid.
    HoleRing,
    /// Every vertex of a ring lies on the class plane — a ring lying in the cut plane.
    AllOnPlane,
    /// An on-plane run's bounding node has no nameable wall plane.
    ///
    /// (There is no sibling for a crossing edge's wall being unreadable from the endpoint
    /// names: the tracer reads the wall the producer carries
    /// (`NamedRing`), so there is no name derivation to fail.)
    RunName,
    /// An on-plane run's node lies on the cut plane yet is not named by it — four planes meet
    /// there. Raised as [`RejectReason::FourPlane`] rather than as a `TraceDeclined`, since the
    /// substrate limit is the cause and the naming failure only the symptom.
    FourPlane,
    /// Two arrangement features on the class line order as equal — they coincide.
    CoincidentFeatures,
    /// A run's two nodes did not end up adjacent after ordering, so the run is not one interval.
    RunSplit,
    /// The sweep along the class line entered and left unequally (an unbalanced parity), which a
    /// closed boundary cannot do.
    OddParity,
    /// **A seated edge's carrier is the face's own plane class**, so `fc ∩ wall` is not a line and
    /// the segment cannot be named on it.
    ///
    /// ★★ **Its sentence changed when the road stopped guessing.** It used to mean "the two
    /// endpoint names share no single plane besides `fc`" — an artefact of *deriving* the carrier,
    /// measured **0** and gone with the derivation. What replaces it is a fact about the carrier
    /// the producer hands over: the face across this edge is in the same class as the face itself,
    /// which the coplanar merge should have folded. One name, one proposition.
    ///
    /// ☑ Measured unexercised (`unreachable!()`, whole suite and ignored sweep green). The guard
    /// exists because the derivation it replaced could not produce this case — it filtered
    /// `c != fc` — and taking the carrier is what makes the case spellable at all.
    SeatedEdgeNaming,
    /// **A seated face's ring rode a cylinder this class never received an element for.**
    ///
    /// The seated walk emits no segment for an arc or a ruling — a [`Seg`](crate::arrangement) is
    /// a straight edge on the class, and that element is already there, contributed by the
    /// cylinder's own lateral face of the same solid. This is the net under that sentence: a
    /// missing element does not decline anywhere, it silently relabels every cell on the class.
    ///
    /// ☑ Measured unexercised — the population gate leaves only the two shapes that do produce
    /// one (a ⊥ class cuts a circle, a class through the axis leaves two rulings).
    SeatedCurveUnbacked,
    /// **A point's name carries no plane that cuts the line it sits on**, so nothing in that name
    /// pins it there — see [`combinatorics::pin_on_line`].
    ///
    /// ★ Distinct from [`Self::FourPlane`], which `third_on_l`'s alias arm raises for the same
    /// missing pin: there the arm has *already* established a fourth plane through the point, so
    /// the substrate limit is the cause and the missing pin only its symptom. Where no such fact is
    /// in hand, saying the symptom is the honest answer.
    ///
    /// ☑ Measured unexercised at **both** its sites: `third_on_l`'s ordinary arm, 0 of 1,622,692
    /// calls where the lone off-line class fails the cut test; and the seated road's, by
    /// `unreachable!()` with the whole suite and the ignored sweep green.
    NoPinOnLine,
    /// A cylinder face could not state its shape exactly — its outer loop could not be cut at
    /// its slits or named, a class has no exact station, or its whole rims are more than two or
    /// not at the ends of its range (a panel and a chain rim are stated, not
    /// declined).
    CylSpan,
    /// **A lateral face has a hole here and this road could not *read* it.**
    ///
    /// A fuse can bury part of a lateral in the other body — a boss straddling a plate's edge keeps
    /// its band but loses the angles inside the plate — and what is left is a band with a hole in
    /// the `(t, θ)` chart. The trace **speaks in arcs** now (`circle_on_class` walks the hole and
    /// answers per angular extent), so what this name is left saying is narrower than it was: the
    /// hole's loop could not be named at all, or a corner of it lies on a second cylinder, or the
    /// ring walk could not decide a node's side.
    ///
    /// ★★★★★ **The name exists because the alternative was being silently wrong.** The face table
    /// described a lateral by its **outer** span, so the trace used to answer "the class cuts a
    /// full
    /// circle" — false for the angles inside the hole. That falsehood did not fail here: it flowed
    /// on, and two stages later `label_cells` found the flip relation broken and said
    /// [`RejectReason::LabelConflict`] — a symptom, not the cause. The refusal belongs where the
    /// false sentence is made, and once the road could state the truth the refusal shrank to what
    /// it still cannot describe.
    ///
    /// ★ A feature the walk *found* but the arc road has no extent for is
    /// [`Self::CylHoleFeature`], not this — the two say different things about the same face.
    CylFaceHole,
    /// **A cycle's boundary — a hole's, a panel's, a chain rim's — met this class in a shape the
    /// arc road has no extent for.**
    ///
    /// Distinct from [`Self::CylFaceHole`], whose proposition is "the hole could not be *read*":
    /// here the ring walked and its features came back, and it is turning one of them into a
    /// counter-clockwise extent that fails — a rim run the hole continues straight through
    /// (an extent and a parity toggle at once), a crossing on a carrier that is not a ruling, an
    /// odd number of crossings, two holes claiming one extent, or boundaries that do not alternate
    /// around the circle. Every one of them is a fact about the *shape*, and none is produced by
    /// today's population; the name is what will say so when one is.
    CylHoleFeature,
    /// A class runs **through a lateral's axis** and the ruling trace could not be stated
    /// exactly — a rim without a ⊥ class to name its ends, or checked arithmetic past `Rat`
    /// (the rulings ladder). Declined whole rather than contributed partially: a
    /// half-contributed rectangle leaves the class's 1-skeleton dangling.
    ///
    /// ★ Not its population: a class plane that holds a **tangent ruling** together with
    /// the wall tangent there — a gusset's side plane exactly through a fillet's axis (`arcwalls
    /// p1p3`), or a box face coplanar with a slot's tangent wall and running past the tangent
    /// point (`rrect-box`). One line on two planes and the cylinder: the corner's names are
    /// joined in the alias table before any trace
    /// (`seed_from_operands` → `Aliases::record_on_cylinder`), the lateral states its own side of
    /// the shared line as a graze segment in the plane vocabulary, and a crossing with a class
    /// circle is a cut point only on an arc. Both populations build; what remains here is the
    /// exact-arithmetic decline the variant is for.
    Ruling,
}

impl DeclineKind {
    /// The stable kebab-case identifier used in logs and the class audit.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CollapsedTriple => "collapsed-triple",
            Self::PierceNode => "pierce-node",
            Self::CurvedRingWall => "curved-ring-wall",
            Self::OuterRing => "outer-ring",
            Self::HoleRing => "hole-ring",
            Self::AllOnPlane => "all-on-plane",
            Self::RunName => "run-name",
            Self::FourPlane => "four-plane",
            Self::CoincidentFeatures => "coincident-features",
            Self::RunSplit => "run-split",
            Self::OddParity => "odd-parity",
            Self::SeatedEdgeNaming => "seated-edge-naming",
            Self::NoPinOnLine => "no-pin-on-line",
            Self::SeatedCurveUnbacked => "seated-curve-unbacked",
            Self::CylSpan => "cyl-span",
            Self::CylFaceHole => "cyl-face-hole",
            Self::CylHoleFeature => "cyl-hole-feature",
            Self::Ruling => "ruling",
        }
    }
}

impl std::fmt::Display for DeclineKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl RejectReason {
    /// The stable snake_case identifier — the same string the reject tags used, so logs and
    /// issue reports do not change meaning across this refactor.
    ///
    /// [`Self::TraceDeclined`] answers `"trace_declined"`; its [`DeclineKind`] carries the
    /// detail (and `Display` prints both).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NonManifoldEdge => "non_manifold_edge",
            Self::NonManifoldResultEdge => "non_manifold_result_edge",
            Self::SelfTouchingResult => "self_touching_result",
            Self::OpenResultShell => "open_result_shell",
            Self::CoincidentNodes => "coincident_nodes",
            Self::RingNaming => "ring_naming",
            Self::DegenerateRing => "degenerate_ring",
            Self::StraightAngle => "straight_angle",
            Self::UnorderedEdges => "unordered_edges",
            Self::EdgeOccupancyConflict => "edge_occupancy_conflict",
            Self::RingOrientation => "ring_orientation",
            Self::RingParity => "ring_parity",
            Self::UnreachedCell => "unreached_cell",
            Self::LabelConflict => "label_conflict",
            Self::TraceDeclined { .. } => "trace_declined",
            Self::CavityNoOwner => "cavity_no_owner",
            Self::NoOutwardShell => "no_outward_shell",
            Self::EulerParity => "euler_parity",
            Self::NonManifoldVertex => "non_manifold_vertex",
            Self::NegativeGenus => "negative_genus",
            Self::ThreePlanes => "three_planes",
            Self::SeamAlias => "seam_alias",
            Self::ZeroLengthEdge => "zero_length_edge",
            Self::FourPlane => "fourplane",
            Self::PierceVertexUnnamed => "pierce_vertex_unnamed",
            Self::VertexNamesAbsentSurface => "vertex_names_absent_surface",
            Self::ObliqueCylinderCut => "oblique_cylinder_cut",
            Self::ObliqueCircleClass => "oblique_circle_class",
            Self::CylinderPairContact => "cylinder_pair_contact",
            Self::CylinderGateUndecided => "cylinder_gate_undecided",
            Self::CylinderFaceUndecided => "cylinder_face_undecided",
            Self::DegenerateFace => "degenerate_face",
            Self::FrameOutOfRange => "frame_out_of_range",
            Self::PrecisionBudget { .. } => "precision_budget",
            Self::JudgeExhausted => "judge_exhausted",
            Self::DegenerateWitness => "degenerate_witness",
            Self::NoClearRay => "no_clear_ray",
            Self::RingHasNoWitness => "ring_has_no_witness",
            Self::WitnessNotRational => "witness_not_rational",
            Self::ArcBoundNotYet => "arc_bound_not_yet",
            Self::CurvedOperandBoundary => "curved_operand_boundary",
            Self::CircleCrossesRuling => "circle_crosses_ruling",
            Self::RulingBoundNotYet => "ruling_bound_not_yet",
            Self::CurvedStraightRun => "curved_straight_run",
            Self::PointOnRing => "point_on_ring",
            Self::HoleDepth => "hole_depth",
            Self::HoleRoots => "hole_roots",
            Self::MissingSeam => "missing_seam",
            Self::CoplanarMerge => "coplanar_merge",
            Self::PartialCircleUncut => "partial_circle_uncut",
        }
    }

    /// What kind of answer this is — see [`RejectClass`]. **Pierce on this, not on the variant.**
    ///
    /// The split follows what each guard's own documentation says it detects: invalid operands
    /// (no valid solid exists) are `Impossible`, coverage limits are `NotSupported`, and
    /// "the arrangement built something malformed" backstops are `SuspectedDefect`.
    pub fn class(self) -> RejectClass {
        match self {
            // The operands are not valid 2-manifolds, or the combination genuinely pinches.
            Self::NonManifoldEdge
            | Self::NonManifoldVertex
            | Self::NonManifoldResultEdge
            | Self::SelfTouchingResult
            | Self::DegenerateFace => RejectClass::Impossible,
            // Built later: quadrics, deeper nesting, rotated-chain witnesses, degenerate
            // arrangements the substrate cannot name yet.
            Self::TraceDeclined { .. }
            | Self::ThreePlanes
            | Self::FourPlane
            | Self::PierceVertexUnnamed
            | Self::VertexNamesAbsentSurface
            | Self::ObliqueCylinderCut
            | Self::CylinderPairContact
            | Self::CylinderGateUndecided
            // A coordinate outside `Rat`'s range: the *kernel* cannot represent it exactly, not
            // that no answer exists — a wider rational would lift this.
            // Likewise a frame past `i128`: a wider rational would lift it.
            | Self::FrameOutOfRange
            // A cost limit, not a resolution one: more bits would answer it.
            | Self::PrecisionBudget { .. }
            | Self::JudgeExhausted
            // The arrangement named a point that is not a point — a degenerate configuration the
            // substrate cannot describe, not a defect in the assembly.
            | Self::DegenerateWitness
            // The substrate cannot name what this configuration asks it to name: a point with
            // four planes through it, a ring it cannot chain, a corner with no turn, or an edge
            // two facts disagree about.
            | Self::CoincidentNodes
            | Self::RingNaming
            | Self::StraightAngle
            | Self::UnorderedEdges
            | Self::EdgeOccupancyConflict
            // Two exact facts about one lateral face's presence, pointing opposite ways: the
            // configuration is outside what this road covers, not a broken arrangement.
            | Self::CylinderFaceUndecided
            | Self::NoClearRay
            | Self::RingHasNoWitness
            | Self::WitnessNotRational
            | Self::ArcBoundNotYet
            | Self::RulingBoundNotYet
            | Self::CurvedOperandBoundary
            | Self::CircleCrossesRuling
            | Self::CurvedStraightRun
            | Self::PointOnRing
            | Self::HoleDepth
            | Self::CoplanarMerge => RejectClass::NotSupported,
            // An invariant broke: malformed assembly, or a backstop that should be unreachable.
            Self::OpenResultShell
            | Self::EulerParity
            | Self::NegativeGenus
            | Self::CavityNoOwner
            | Self::NoOutwardShell
            | Self::SeamAlias
            | Self::ZeroLengthEdge
            // The arrangement built something that is not a subdivision: a ring that bounds
            // nothing, an edge used once, faces that will not close, a cell with no side.
            | Self::DegenerateRing
            | Self::RingOrientation
            | Self::RingParity
            | Self::UnreachedCell
            | Self::LabelConflict
            | Self::HoleRoots
            | Self::PartialCircleUncut
            | Self::ObliqueCircleClass
            | Self::MissingSeam => RejectClass::SuspectedDefect,
        }
    }
}

impl std::fmt::Display for RejectReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TraceDeclined { kind, .. } => write!(f, "{}({kind})", self.as_str()),
            _ => f.write_str(self.as_str()),
        }
    }
}

/// Build an `Rejected` carrying *which* guard raised it.
///
/// Every `Rejected` in this crate is built here or in [`reject_at`] (this, with a location).
/// The reason rides in the returned value, so a guard that is raised and then swallowed by an
/// alternative path (several sites try another route on `Err`) can never be mistaken for the
/// one that actually surfaced.
///
/// The pair being the only funnel is also what makes [`reject_census`] possible:
/// `#[track_caller]` here records *which* guard rang, at its own line, for every call site
/// at once.
#[inline]
#[track_caller]
pub(crate) fn reject(reason: RejectReason) -> BoolError {
    reject_census::raised(reason, std::panic::Location::caller());
    BoolError::Rejected { reason, at: None }
}

/// [`reject`], carrying where the guard was looking — same census instrumentation
/// (`#[track_caller]` sees through to the raise site).
#[inline]
#[track_caller]
pub(crate) fn reject_at(reason: RejectReason, at: RejectWhere) -> BoolError {
    reject_census::raised(reason, std::panic::Location::caller());
    BoolError::Rejected {
        reason,
        at: Some(at),
    }
}

/// Assert that `f` rejects *through the intended guard* — a reject test whose fixture drifts
/// onto a different guard then fails instead of silently passing. Compares the projected
/// reason only: the location payload is a measurement, asserted (approximately) by the
/// per-reason payload tests, not here.
#[cfg(test)]
fn assert_rejects<T: std::fmt::Debug + PartialEq>(
    f: impl FnOnce() -> Result<T, BoolError>,
    expect: RejectReason,
) {
    match f() {
        Err(BoolError::Rejected { reason, .. }) => assert_eq!(reason, expect),
        other => panic!("expected a reject with {expect:?}, got {other:?}"),
    }
}

/// The start vertex of a half-edge (`vertices[0]` if forward, else `vertices[1]`).
pub(crate) fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
    model.he_start(he)
}

use std::collections::HashMap;

fn unordered(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

#[cfg(test)]
#[path = "tests/suite/mod.rs"]
pub mod tests;

/// **Ledger totals, printed** — a measurement, not an assertion. Run last, single-threaded, so
/// the sums cover the whole lib suite:
/// `cargo test -p nacre-ops --lib -- --include-ignored --test-threads=1 --nocapture zzz_ledger`
/// is *not* it (a filter runs only this); the whole-suite run is
/// `cargo test -p nacre-ops --lib -- --include-ignored --test-threads=1 --nocapture --skip stress --skip spike --skip direction_families`.
#[cfg(test)]
#[path = "tests/zzz_ledger.rs"]
mod zzz_ledger;
