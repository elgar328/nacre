//! Operations for the nacre kernel, plus a replayable operation log (design §6).
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::PadOnFace`]/[`Operation::PocketOnFace`] (M4) consume a prior op's face by
//! `Handle` (exposed via [`OpOutput`]) and supersede a solid (design §2 live-solid
//! semantics) — each is a tool prism plus a boolean, not a direct face-split. Ops are
//! applied by [`apply`] and folded by [`replay`]; every result is a **closed** solid, so
//! `nacre-validate` applies fully.

use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::Rat;
use nacre_store::Handle;
use nacre_topo::{Face, HalfEdge, Model, Solid, Vertex};

mod arrangement;
mod boolean;
mod combinatorics;
mod exact;
mod ops;
mod par;
mod planes;
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
pub use ops::{
    BoolKind, OpError, OpOutput, Operation, PlaneDef, Profile2d, ProfileRing, SketchPlane, apply,
    face_plane, replay,
};
pub use sketch::{Curve2d, Edge2d, SketchError, from_edges, from_rings};

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
    fn axis_plane(n: [i128; 3], u: [i128; 3], v: [i128; 3]) -> Self {
        let r = |a: [i128; 3]| a.map(Rat::from_int);
        let f = |a: [i128; 3]| Vector3::from_array(a.map(|c| c as f64));
        let zero = Rat::from_int(0);
        Self {
            origin: Point3::origin(),
            x_axis: f(u),
            y_axis: f(v),
            def: Some(PlaneDef {
                coeffs: [r(n)[0], r(n)[1], r(n)[2], zero],
                origin: [zero; 3],
                ref_dir: r(u),
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
    fn normal_def(origin: Point3, normal: Vector3) -> Option<PlaneDef> {
        let o = origin.as_array().map(Rat::from_decimal);
        let n = normal.as_array().map(Rat::from_decimal);
        let (o, n) = ([o[0]?, o[1]?, o[2]?], [n[0]?, n[1]?, n[2]?]);
        let coeffs = nacre_scalar::plane_from_point_normal(n, o)?;
        let zero = Rat::from_int(0);
        let ref_dir = if n[0] == zero && n[1] == zero {
            [n[2], zero, zero] // ŷ × n for a vertical normal
        } else {
            [zero.checked_sub(n[1])?, n[0], zero] // ẑ × n
        };
        Some(PlaneDef {
            coeffs,
            origin: o,
            ref_dir,
        })
    }

    /// **A plane through three written points**: `origin` is the sketch's `(0, 0)`, `+u` runs
    /// toward `x_point`, and `+v` leans toward `y_hint`.
    ///
    /// ★★★ **Everything here is exact by construction.** The plane is
    /// [`nacre_scalar::plane_through_points`] of the three; `ref_dir` is `x_point − origin`, a
    /// difference of written points that **already lies in the plane**. `None` if the three are
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
            let coeffs = nacre_scalar::plane_through_points(o, xp, yh)?;
            let d = |i: usize| xp[i].checked_sub(o[i]);
            Some(PlaneDef {
                coeffs,
                origin: o,
                ref_dir: [d(0)?, d(1)?, d(2)?],
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
    /// Keeping the coefficients while moving the origin would leave the definition describing one
    /// plane and its origin sitting on another: measured, `world_xy().with_origin([0, 0, 0.5])`
    /// recorded `z = 0` for a cap at `z = 0.5`, and a boolean built on that lost 0.04 of volume.
    /// The invariant that stops it — *the origin satisfies the coefficients* — is asserted in
    /// `a_named_plane_records_what_its_caller_stated`.
    ///
    /// The direction is untouched: a translation does not turn `+u`.
    pub fn with_origin(mut self, p: Point3) -> Self {
        self.origin = p;
        self.def = self.def.and_then(|d| {
            let a = p.as_array().map(Rat::from_decimal);
            let origin = [a[0]?, a[1]?, a[2]?];
            Some(PlaneDef {
                coeffs: nacre_scalar::plane_from_point_normal(
                    [d.coeffs[0], d.coeffs[1], d.coeffs[2]],
                    origin,
                )?,
                origin,
                ref_dir: d.ref_dir,
            })
        });
        self
    }

    /// A frame from axes the caller already holds — **with no exact definition**, so anything
    /// built on it takes the f64 path.
    ///
    /// ★ Not public. Normalized axes are exactly what this type exists to stop being handed, and
    /// the callers that legitimately have only axes are internal: a face's own frame, and
    /// `exact()`'s own tests.
    pub(crate) fn from_axes(origin: Point3, x_axis: Vector3, y_axis: Vector3) -> Self {
        Self {
            origin,
            x_axis,
            y_axis,
            def: None,
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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolError {
    /// The engine declined to answer, and [`RejectReason`] names *which* guard spoke.
    /// The reason travels in the value so a consumer can say why, and so the site that
    /// returned is the site that is reported (guards that are raised and then swallowed
    /// by an alternative path cannot be mistaken for the surfaced one).
    Unsupported { reason: RejectReason },
    /// An input solid handle is not in `model.live_solids`.
    InputNotLive,
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
    /// Inside the kernel's remit but outside what is built yet — the same input may succeed
    /// at a later milestone. "Not supported yet."
    NotSupportedYet,
    /// No valid solid exists for this input, at any milestone: the operands are not valid
    /// 2-manifolds, or the requested combination pinches. The design has to change.
    Impossible,
    /// An engine invariant broke: the arrangement built something malformed, or a backstop
    /// that should be unreachable spoke. Reported rather than returned (DNA: never silently
    /// wrong), and worth a bug report. Say "could not produce a valid result", not "your fault".
    ///
    /// The classification of variants that have never been observed to fire is provisional —
    /// tighten it once the reason census (dev-log) says which are reachable.
    SuspectedDefect,
}

/// Which guard raised an [`BoolError::Unsupported`].
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
    /// different situation and used to share this name. Unfired across the suite (2026-07-26).
    NonManifoldEdge,
    /// The assembled boundary uses an edge **more than twice**: the two bodies meet exactly along
    /// that edge, so the result would pinch there and no 2-manifold solid contains it. The edge
    /// twin of [`Self::NonManifoldVertex`], and `Impossible` for the same reason — a design that
    /// leans on an exact edge-to-edge touch has no valid answer at any milestone.
    ///
    /// **The most common reject in the suite** (2026-07-26 census): the grid proptest lands on it
    /// whenever two sampled boxes share exactly an edge, e.g. `[2,4]×[2,4]×[1,4]` fused with
    /// `[2,4]×[4,8]×[4,6]`, which touch only along the line `y = 4, z = 4`.
    NonManifoldResultEdge,
    /// The assembled boundary leaves an edge used **once** — a dangling edge, so the face set is
    /// not closed. Unlike [`Self::NonManifoldResultEdge`] this says nothing bad about the input:
    /// the assembly dropped a face, which is ours to fix. Unfired in the suite (2026-07-26).
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
    StraightAngle,
    /// Two ring edges meet at more than one vertex, or at none, so which vertex the corner *is*
    /// cannot be decided (a two-gon, or a self-bounded rim).
    AmbiguousCorner,
    /// **Two edges leave one arrangement vertex at the same angle**, so the cyclic order around
    /// that vertex has no answer — and the face walk is built from exactly that order.
    ///
    /// Their lines are then identical (parallel plus a shared point), which upstream is supposed
    /// to have already resolved: `merge_coincident` folds an edge traced twice, and `Aliases`
    /// folds two walls that carry one line. Reaching here means one of those did not, and the
    /// honest answer is that this arrangement cannot be ordered rather than an order picked
    /// arbitrarily — the ordering feeds `next`, so a guess there is a wrong face, silently.
    UnorderedEdges,
    /// A result face's edge is used by only one plane, so the boundary is open there. The
    /// under-used twin of [`Self::NonManifoldEdge`], which is the over-used case.
    UnpairedSeamEdge,
    /// **Two facts about one edge disagree**: a graze coincident with a same-solid crossing, two
    /// seated faces claiming opposite sides, or one solid crossing the same edge twice. Whichever
    /// is right, nothing on that edge says which.
    EdgeOccupancyConflict,
    /// No assignment of the class's rings to cells leaves exactly one outer boundary per
    /// component — the planar subdivision does not close into faces the way a subdivision must.
    RingOrientation,
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
    /// A face lies on a surface whose exact definition the kernel cannot state
    /// (`SurfaceDef::Inexact`). Judging on rounded coefficients is what produces two plane classes
    /// for one wall, so the operation declines instead of pretending.
    ///
    /// **It no longer means "rotated, then translated"** — the motion history names that now. What
    /// is left is a surface with no recorded provenance at all (`nacre-validate` reports such a
    /// model) and one defensive branch in the mirror path, so a firing of this today points at a
    /// producer that skipped `Model::push_surface`, not at a motion the kernel cannot describe.
    InexactSurface,
    /// The assembled result has an **odd Euler characteristic** (`V − E + F − L_i`), which no
    /// closed 2-manifold can have (it must equal the even `2(S − G)`) — so the arrangement produced
    /// a malformed solid and the boolean rejects rather than return it (DNA: never silently wrong).
    /// This is the Euler-parity backstop for malformity that is *not* a pinch (see
    /// `NonManifoldVertex`); e.g. a dropped face. Rotation-independent. Checked post-assembly in
    /// `boolean`, per solid.
    EulerParity,
    /// The assembled result has a **non-manifold vertex** — a "pinch" where two or more face-fans
    /// meet at one point (a cutter's convex corner exactly on the target's concave corner; two
    /// solids touching only at a corner), even though every edge is manifold. No valid 2-manifold
    /// solid has one, so the boolean rejects with this clear reason rather than the incidental
    /// `EulerParity` (which also misses an *even* number of pinches). Rotation-independent; checked
    /// per solid post-assembly via `nacre_topo::nonmanifold_vertices`.
    NonManifoldVertex,
    /// The assembled result has an even Euler characteristic but a **negative genus** (`S − χ/2 < 0`)
    /// — more handles than a solid can have, so it is not a valid closed 2-manifold. A count-based
    /// backstop below the pinch and parity checks; checked per solid post-assembly.
    NegativeGenus,
    /// Three planes that should meet in a point do not (a parallel pair), so an arrangement
    /// vertex has no name.
    ThreePlanes,
    /// Two **different** arrangement vertices (distinct plane triples) materialized to the same
    /// coordinate. The triple is the truth and the coordinate only its cache (overview §5), so this
    /// says the exact substrate and the f64 cache disagree about how many vertices exist — always a
    /// defect upstream, never a property of the input. Raised where the seam table is built, while
    /// both triples are still in hand; without it the disagreement surfaces much later as a
    /// zero-length edge. The known cause is a **split plane table** (one geometric plane carried by
    /// two classes); a genuine 4-plane concurrency would do the same.
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
    /// What survives is the guard, at two sites: an **operand** vertex found on more than three
    /// plane classes (`combinatorics`), and a run whose candidate handles are all parallel to the
    /// line they must cut (`arrangement`). Neither has a reproduction. It stays because a
    /// substrate that cannot name a point must say so rather than pick one of the names — and
    /// because an unfired reject costs nothing, while a missing one costs a wrong solid.
    ///
    /// **Exact, not toleranced.** Measured on that model: perturbing an operand coordinate by
    /// **one ULP** in either direction removed the concurrency. The judgement returned zero
    /// because the determinant *is* zero at 200 bits, not because it fell below the coincidence
    /// limit (~55 orders of magnitude lower). Exact rational input would produce the same
    /// concurrency — this was never an artefact of `f64` construction, which is why exact
    /// rational construction (`crate::exact`) left it exactly where it was.
    FourPlane,
    /// An operand carries a cylindrical face. The planar engine covers planes only (M6 adds
    /// quadrics).
    CylinderFace,
    /// An operand face has no three non-collinear outer-loop points, so it spans no plane.
    DegenerateFace,
    /// An operand face's outer triangle has a zero-length normal, so it has no outward direction.
    DegenerateNormal,
    /// An operand coordinate lies outside the exact rational scalar's range, so its plane has no
    /// exact definition to reason with. `Rat` is `Ratio<i128>`: the numerator `mantissa · 2^exp`
    /// must fit (`|x| ≲ 1.7e38`) and so must the denominator `2^k` (`|x| ≳ 2^-74 ≈ 5.3e-23`);
    /// exact zero is always fine. **No CAD model lives at either extreme** — this exists so the
    /// kernel says so by name. (It used to panic on the mixed-rotation path and silently succeed
    /// on the axis-aligned one, because the definition was built lazily, per judgment.)
    CoordinateOutOfRange,
    /// A plane's own frame cannot be stated exactly, so a sketch built in it has no exact
    /// definition to judge from. Either the plane recorded no rational coefficients, or the
    /// squared lengths the frame's realization divides by do not fit `i128`.
    ///
    /// ★ **Distinct from [`Self::InexactSurface`] on purpose.** Both end with "no exact
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
    DegenerateWitness,
    /// Every candidate ray from a loop's nodes has a ring node on its line.
    ///
    /// `point_in_ring` casts along `P ∩ Q_a` for a node's own plane `Q_a`; a ring node on
    /// that line makes the crossing parity ambiguous. Candidates are `2 · |loop|` lines and
    /// two directions, and half of them can be spoiled at once — `l_and_staple`'s loop and
    /// arc share both `y` planes, so only the `x` lines are clear there. Unfired today.
    NoClearRay,
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
    /// which the E1 union cell now builds, and `Cut`/opposite-normal, which is exact whenever the
    /// contact plane separates the two solids (nothing to remove). What is left is a `Cut` whose tool
    /// reaches back across that plane — a pin below its own contact face — where the cut owes a notch
    /// this path cannot yet cut. Honest reject rather than a whole cap that ignores the pin.
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
    /// The face's outer ring could not be named as plane triples.
    OuterRing,
    /// One of the face's hole rings could not be named. Not "no hole": swallowing it would trace
    /// the face as if it were solid.
    HoleRing,
    /// Every vertex of a ring lies on the class plane — a ring lying in the cut plane.
    AllOnPlane,
    /// A strictly crossing ring edge has no nameable wall plane beside the face's own.
    CrossingName,
    /// An on-plane run's bounding node has no nameable wall plane.
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
    /// A seated (on-plane) edge has no unique wall plane, so its segment cannot be named.
    SeatedEdgeNaming,
}

impl DeclineKind {
    /// The stable kebab-case identifier used in logs and the class audit.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CollapsedTriple => "collapsed-triple",
            Self::OuterRing => "outer-ring",
            Self::HoleRing => "hole-ring",
            Self::AllOnPlane => "all-on-plane",
            Self::CrossingName => "crossing-name",
            Self::RunName => "run-name",
            Self::FourPlane => "four-plane",
            Self::CoincidentFeatures => "coincident-features",
            Self::RunSplit => "run-split",
            Self::OddParity => "odd-parity",
            Self::SeatedEdgeNaming => "seated-edge-naming",
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
            Self::OpenResultShell => "open_result_shell",
            Self::CoincidentNodes => "coincident_nodes",
            Self::RingNaming => "ring_naming",
            Self::DegenerateRing => "degenerate_ring",
            Self::StraightAngle => "straight_angle",
            Self::AmbiguousCorner => "ambiguous_corner",
            Self::UnorderedEdges => "unordered_edges",
            Self::UnpairedSeamEdge => "unpaired_seam_edge",
            Self::EdgeOccupancyConflict => "edge_occupancy_conflict",
            Self::RingOrientation => "ring_orientation",
            Self::UnreachedCell => "unreached_cell",
            Self::LabelConflict => "label_conflict",
            Self::TraceDeclined { .. } => "trace_declined",
            Self::CavityNoOwner => "cavity_no_owner",
            Self::NoOutwardShell => "no_outward_shell",
            Self::InexactSurface => "inexact_surface",
            Self::EulerParity => "euler_parity",
            Self::NonManifoldVertex => "non_manifold_vertex",
            Self::NegativeGenus => "negative_genus",
            Self::ThreePlanes => "three_planes",
            Self::SeamAlias => "seam_alias",
            Self::ZeroLengthEdge => "zero_length_edge",
            Self::FourPlane => "fourplane",
            Self::CylinderFace => "cylinder_face",
            Self::DegenerateFace => "degenerate_face",
            Self::DegenerateNormal => "degenerate_normal",
            Self::CoordinateOutOfRange => "coordinate_out_of_range",
            Self::FrameOutOfRange => "frame_out_of_range",
            Self::PrecisionBudget { .. } => "precision_budget",
            Self::JudgeExhausted => "judge_exhausted",
            Self::DegenerateWitness => "degenerate_witness",
            Self::NoClearRay => "no_clear_ray",
            Self::PointOnRing => "point_on_ring",
            Self::HoleDepth => "hole_depth",
            Self::HoleRoots => "hole_roots",
            Self::MissingSeam => "missing_seam",
            Self::CoplanarMerge => "coplanar_merge",
        }
    }

    /// What kind of answer this is — see [`RejectClass`]. **Branch on this, not on the variant.**
    ///
    /// The split follows what each guard's own documentation says it detects: invalid operands
    /// (no valid solid exists) are `Impossible`, coverage limits are `NotSupportedYet`, and
    /// "the arrangement built something malformed" backstops are `SuspectedDefect`.
    pub fn class(self) -> RejectClass {
        match self {
            // The operands are not valid 2-manifolds, or the combination genuinely pinches.
            Self::NonManifoldEdge
            | Self::NonManifoldVertex
            | Self::NonManifoldResultEdge
            | Self::DegenerateFace
            | Self::DegenerateNormal => RejectClass::Impossible,
            // Built later: quadrics, deeper nesting, rotated-chain witnesses, degenerate
            // arrangements the substrate cannot name yet.
            Self::TraceDeclined { .. }
            | Self::InexactSurface
            | Self::ThreePlanes
            | Self::FourPlane
            | Self::CylinderFace
            // A coordinate outside `Rat`'s range: the *kernel* cannot represent it exactly, not
            // that no answer exists — a wider rational would lift this.
            | Self::CoordinateOutOfRange
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
            | Self::AmbiguousCorner
            | Self::UnorderedEdges
            | Self::EdgeOccupancyConflict
            | Self::NoClearRay
            | Self::PointOnRing
            | Self::HoleDepth
            | Self::CoplanarMerge => RejectClass::NotSupportedYet,
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
            | Self::UnpairedSeamEdge
            | Self::RingOrientation
            | Self::UnreachedCell
            | Self::LabelConflict
            | Self::HoleRoots
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

/// Build an `Unsupported` carrying *which* guard raised it.
///
/// Every `Unsupported` in this crate is built here. The reason rides in the returned value, so
/// a guard that is raised and then swallowed by an alternative path (several sites try another
/// route on `Err`) can never be mistaken for the one that actually surfaced.
#[inline]
pub(crate) fn reject(reason: RejectReason) -> BoolError {
    BoolError::Unsupported { reason }
}

/// Assert that `f` rejects *through the intended guard* — a reject test whose fixture drifts
/// onto a different guard then fails instead of silently passing.
#[cfg(test)]
fn assert_rejects<T: std::fmt::Debug + PartialEq>(
    f: impl FnOnce() -> Result<T, BoolError>,
    expect: RejectReason,
) {
    assert_eq!(f(), Err(BoolError::Unsupported { reason: expect }));
}

/// The start vertex of a half-edge (`bounds[0]` if forward, else `bounds[1]`).
/// Every half-edge walked here belongs to a valid solid, so its edge is bounded.
pub(crate) fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
    // A solid's loop edge is always bounded; only the standalone full circle of
    // design §4 is not, and that is never part of a face's loop.
    model.he_start(he).expect("a solid's loop edge is bounded")
}

use std::collections::HashMap;

fn unordered(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

#[cfg(test)]
pub mod tests {

    use super::*;
    use crate::tolerant::Judge;
    use crate::transform::transform;
    use crate::{boolean::*, ops::*, planes::*};
    use nacre_cip::Pt3;
    use nacre_geom::intersect::{planes_coplanar, three_planes};
    use nacre_geom::{Plane, Surface};
    use nacre_topo::{Loop, Orientation, Origin, VertexDef};
    use proptest::prelude::*;
    use std::collections::HashMap;

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
            let Surface::Plane(plane) = m.surface(f.surface) else {
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
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)])
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
        Profile2d::polygon(points)
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
        let tri = Profile2d::polygon(vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(1.0, 1.5)]);
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
        let l = Profile2d::polygon(vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ]);
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 12);
        assert_eq!(m.faces.len(), 8);
    }

    /// The L-prism: profile `[(0,0),(2,0),(2,1),(1,1),(1,2),(0,2)]` extruded to
    /// z ∈ [0,1]. Material = bottom bar (x∈[0,2],y∈[0,1]) ∪ left bar (x∈[0,1],
    /// y∈[1,2]); the notch (x∈[1,2],y∈[1,2]) is empty.
    fn l_prism() -> (Model, Handle<Solid>) {
        let l = Profile2d::polygon(vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ]);
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// The same L-prism, its profile started one vertex earlier so the reflex corner
    /// `(1,1)` lands at index 1 of the cap's loop. Geometrically identical.
    fn rotated_l_prism() -> (Model, Handle<Solid>) {
        let l = Profile2d::polygon(vec![
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
            p2(0.0, 0.0),
            p2(2.0, 0.0),
        ]);
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// `FaceInfo::n_out` is documented as the single source of "outward". Two
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
        // `collect_planes` debug_asserts, per face, that `sign(normal·n_out)` equals the topo
        // `face.orientation` — the invariant this test used to check by reading a `.orient` field
        // it kept alongside. That field is gone (its production role is `orient_sign`), so the
        // check lives at construction now; running collect_planes on shapes with reversed faces
        // (the L/U notch, a rotated solid) exercises it. Here we add the parallel invariant, which
        // is not debug_asserted the same way.
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
            }
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

    /// The L-prism with a `[0.1,0.9]³` box strictly inside its bottom bar
    /// (non-coplanar coordinates ⇒ no shared face planes). `V_L = 3`, `V_box =
    /// 0.512`.
    fn l_and_inner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
        (m, l, bx)
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

    // --- Rotated booleans go live (overhaul 3d-i, `ROTATED_UNSUPPORTED` retired) ---
    // A boolean commutes with a rigid motion, so rotating both operands by the same
    // irrational-angle isometry must give the rigid image of the unrotated result — identical
    // volume, solid count, and cavity count, and still valid. These are the first live proof
    // that the CIP-wired machinery (arrangement, seam, in/out, outer/cavity — 3a–3c-vi) is
    // sound end-to-end on rotated (rounded-irrational) geometry.

    /// A cutter's convex corner landing exactly on the target's concave corner — three planes
    /// (x=1,y=1,z=1) meeting at one point (1,1,1), six faces there — makes a **non-manifold pinch**
    /// (two face-fans touching at the point). No valid 2-manifold solid exists, so `boolean` rejects
    /// with the clear `NON_MANIFOLD_VERTEX` reason (not the incidental `EULER_PARITY`), leaving the
    /// live model untouched. **Rotation-independent**: both the axis-aligned and the rotated framings
    /// (exact and CIP-kernel paths) hit the same pinch and reject. R = a cube minus a far-corner
    /// octant; C = a cube whose +corner is that removed octant's inner corner. (Coplanar contact away
    /// from a corner works — [`a_rotated_boolean_result_can_be_cut_again`].)
    #[test]
    fn a_corner_coincident_cut_is_rejected_not_silently_wrong() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        // `rotate`: None = axis-aligned (exact path); Some = the result and cutter tilted (CIP path).
        let run = |rotate: bool| {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
            let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
            // C's +corner is (1,1,1) = R's concave corner → three shared planes meet there.
            let c = m.add_cuboid(
                Point3::from_array([-1.0, -1.0, -1.0]),
                Point3::from_array([1.0; 3]),
            );
            m.rebuild_adjacency();
            let (r, c) = if rotate {
                let r = transform(&mut m, r, &iso).unwrap();
                m.rebuild_adjacency();
                let c = transform(&mut m, c, &iso).unwrap();
                m.rebuild_adjacency();
                (r, c)
            } else {
                (r, c)
            };
            let live = m.live_solids.clone();
            assert_rejects(
                || boolean(&mut m, BoolKind::Cut, r, c),
                RejectReason::NonManifoldVertex,
            );
            assert_eq!(m.live_solids, live, "reject must not mutate the live set");
        };
        run(false); // axis-aligned
        run(true); // rotated
    }

    /// Two cubes touching only at the corner (1,1,1): their Fuse pinches two solids at a single
    /// vertex (non-manifold), so it is rejected with the clear `NON_MANIFOLD_VERTEX` — the direct,
    /// minimal pinch (every edge is manifold; only the vertex is the defect).
    #[test]
    fn two_cubes_touching_at_a_corner_fuse_is_non_manifold() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        m.rebuild_adjacency();
        let live = m.live_solids.clone();
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, a, b),
            RejectReason::NonManifoldVertex,
        );
        assert_eq!(m.live_solids, live, "reject must not mutate the live set");
    }

    /// The kernel's validator now catches the pinch too (it previously only exposed it indirectly as
    /// `EulerParity`). Bypass `boolean`'s reject via `arrangement::boolean` to obtain the malformed
    /// corner-touch Fuse solid, then confirm `validate` reports a `NonManifoldVertex`.
    #[test]
    fn validate_reports_the_non_manifold_pinch() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        m.rebuild_adjacency();
        crate::arrangement::boolean(&mut m, BoolKind::Fuse, a, b).unwrap();

        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(
            issues
                .iter()
                .any(|v| matches!(v, nacre_validate::Violation::NonManifoldVertex { .. })),
            "validate must flag the pinch: {issues:?}"
        );
    }

    /// Predicates over a rotated result's witness planes are rotation-invariant against the same
    /// result unrotated — the provenance witness (with outward winding) defines the exact plane,
    /// so `Judge::orient3d` agrees on every definite triple. A regression guard on the witness itself.
    #[test]
    fn rotated_result_witness_predicates_are_invariant() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let table = |m: &Model, s: Handle<Solid>| {
            let f = collect_planes(m, s).unwrap();
            let c = plane_classes(&crate::planes::test_judge(&f));
            dense_planes(&f, &c).0
        };
        let pu = table(&m, r);
        let r2 = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();
        let pr = table(&m, r2);
        assert_eq!(pu.len(), pr.len(), "rotation preserves the plane count");
        let n = pu.len();
        let indep = |p: &[PlaneGeom], a: usize, b: usize, c: usize| {
            let nrm = |k: usize| p[k].plane.normal();
            nrm(a).dot(nrm(b).cross(nrm(c))).abs() > 0.3
        };
        let mut disagree = 0;
        for p in 0..n {
            for q in (p + 1)..n {
                for rr in (q + 1)..n {
                    if !indep(&pu, p, q, rr) {
                        continue;
                    }
                    for j in 0..n {
                        if j == p || j == q || j == rr {
                            continue;
                        }
                        let su = crate::planes::test_judge(&pu).orient3d(p, q, rr, j);
                        let sr = crate::planes::test_judge(&pr).orient3d(p, q, rr, j);
                        if su != 0 && sr != 0 && su != sr {
                            eprintln!("DISAGREE orient3d ({p},{q},{rr},{j}): u={su} r={sr}");
                            disagree += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(
            disagree, 0,
            "{disagree} predicate disagreements (witness wrong)"
        );
    }

    /// Rotating a boolean *result* and feeding it back into a boolean (was `ROTATED_UNSUPPORTED`):
    /// `collect_planes` now witnesses each rotated seam face's plane through provenance — its
    /// plane is `R(π)` for an operand plane `π`, recovered from the operand face still on `π` and
    /// rotated by the face's own chain. A boolean commutes with a rigid motion, so the rotated
    /// chain's result matches the unrotated chain's (volume, solid count) and stays valid.
    #[test]
    fn a_rotated_boolean_result_can_be_cut_again() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        // Chain: R = Cut(A, B) removes a far-corner octant; then Cut(R, C) removes a near one.
        let build = |m: &mut Model| -> (Handle<Solid>, Handle<Solid>) {
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
            let r = boolean_one(m, BoolKind::Cut, a, b).unwrap();
            // A clean slab that severs R at x = 0.5 — no plane of C coincides with any of R's
            // (avoids the separate rotated-coplanar-contact gap; isolates the witness).
            let c = m.add_cuboid(
                Point3::from_array([-1.0, -1.0, -1.0]),
                Point3::from_array([0.5, 4.0, 4.0]),
            );
            (r, c)
        };
        // Unrotated reference (reuse already works when nothing is rotated).
        let mut m0 = Model::new();
        let (r0, c0) = build(&mut m0);
        m0.rebuild_adjacency();
        let ref_out = boolean(&mut m0, BoolKind::Cut, r0, c0).unwrap();
        m0.rebuild_adjacency();
        let ref_vol: f64 = ref_out
            .iter()
            .map(|&s| nacre_props::mass_props(&m0, s).unwrap().volume)
            .sum();
        // Rotated: turn the *result* R (and C) by the same isometry, then reuse R.
        let mut m = Model::new();
        let (r, c) = build(&mut m);
        m.rebuild_adjacency();
        let r = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();
        let c = transform(&mut m, c, &iso).unwrap();
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Cut, r, c).unwrap();
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(
            issues.is_empty(),
            "rotated-result reuse must be valid: {issues:?}"
        );
        let vol: f64 = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert_eq!(
            out.len(),
            ref_out.len(),
            "solid count invariant under rotation"
        );
        assert!(
            (vol - ref_vol).abs() < 1e-6,
            "rotated reuse volume {vol} vs unrotated {ref_vol}"
        );
    }

    /// Result-reuse rotation stress: build R with a first boolean, then feed R into a second
    /// boolean with a fresh cutter C — once unrotated, once with R and C rotated by the same
    /// isometry. A boolean commutes with a rigid motion, so the rotated reuse must equal the
    /// unrotated one (volume, solid count, cavity count) or be an honest reject — never silently
    /// wrong. This is the invariant on the newly-enabled rotated-*Discovered* geometry (every
    /// vertex of a boolean result is `Discovered`, so its rotation exercises the provenance
    /// witness on every face). `#[ignore]`: rotated booleans escalate to astro-float (~1–3 s).
    #[test]
    #[ignore = "slow: rotated result-reuse booleans (run with --ignored)"]
    fn rotated_result_reuse_stress() {
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
        // (first kind, second kind, |m| -> (a, b, c)). R = kind1(a, b); out = kind2(R, c).
        type Build = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>, Handle<Solid>)>;
        let cuboid = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
            m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
        };
        let fixtures: Vec<(&str, BoolKind, BoolKind, Build)> = vec![
            (
                "corner_then_slab",
                BoolKind::Cut,
                BoolKind::Cut,
                Box::new(move |m: &mut Model| {
                    let a = cuboid(m, [0.0; 3], [2.0; 3]);
                    let b = cuboid(m, [1.0; 3], [3.0; 3]);
                    let c = cuboid(m, [-1.0, -1.0, -1.0], [0.5, 4.0, 4.0]);
                    (a, b, c)
                }),
            ),
            (
                "fuse_then_bite",
                BoolKind::Fuse,
                BoolKind::Cut,
                Box::new(move |m: &mut Model| {
                    let a = cuboid(m, [0.0; 3], [2.0, 1.0, 1.0]);
                    let b = cuboid(m, [0.0, 0.0, 0.0], [1.0, 2.0, 1.0]);
                    let c = cuboid(m, [1.3, 1.3, -1.0], [3.0, 3.0, 2.0]);
                    (a, b, c)
                }),
            ),
            (
                "cut_then_fuse",
                BoolKind::Cut,
                BoolKind::Fuse,
                Box::new(move |m: &mut Model| {
                    let a = cuboid(m, [0.0; 3], [2.0; 3]);
                    let b = cuboid(m, [1.3, 1.3, 1.3], [3.0, 3.0, 3.0]);
                    let c = cuboid(m, [-0.7, 0.4, 0.4], [0.3, 1.4, 1.4]);
                    (a, b, c)
                }),
            ),
        ];
        let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
            ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
            ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
            (
                "Z50>Y37",
                vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::Y, 37, [0, 0, 1])],
            ),
        ];
        let run = |k1: BoolKind, k2: BoolKind, build: &Build, isos: &[Isometry]| -> Out {
            let mut m = Model::new();
            let (a, b, c) = build(&mut m);
            m.rebuild_adjacency();
            let Ok(r) = boolean(&mut m, k1, a, b) else {
                return Out::Rej;
            };
            assert_eq!(r.len(), 1, "first boolean is a single solid");
            let mut r = r[0];
            let mut c = c;
            m.rebuild_adjacency();
            for iso in isos {
                r = transform(&mut m, r, iso).unwrap();
                m.rebuild_adjacency();
                c = transform(&mut m, c, iso).unwrap();
                m.rebuild_adjacency();
            }
            match boolean(&mut m, k2, r, c) {
                Ok(solids) => {
                    m.rebuild_adjacency();
                    assert!(
                        nacre_validate::validate(&m).is_empty(),
                        "INVALID rotated reuse result"
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
        let (mut success, mut reject, mut silent) = (0, 0, 0);
        for (fname, k1, k2, build) in &fixtures {
            let base = run(*k1, *k2, build, &[]);
            for (rname, isos) in &isos_list {
                let r = run(*k1, *k2, build, isos);
                match (&base, &r) {
                    (_, Out::Rej) => reject += 1,
                    (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2))
                        if vclose(*v1, *v2) && s1 == s2 && c1 == c2 =>
                    {
                        success += 1
                    }
                    _ => {
                        silent += 1;
                        eprintln!("SILENT-WRONG {fname} {rname} base={base:?} rot={r:?}");
                    }
                }
            }
        }
        eprintln!("REUSE STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent}");
        assert_eq!(
            silent, 0,
            "a rotated result-reuse boolean was silently wrong"
        );
        assert!(
            success >= 1,
            "at least one rotated reuse must actually succeed"
        );
    }

    /// Adversarial rotation stress (overhaul 3d-iii): many fixtures × kinds × rotations
    /// (single-axis, and Euler chains reaching arbitrary orientation) confirm the DNA
    /// invariant — a rotated boolean is *never silently wrong*: its result either equals the
    /// unrotated one (a boolean commutes with a rigid motion, so volume/solid-count/cavity-count
    /// are invariant) or is an honest reject. `#[ignore]`: each rotated boolean escalates its
    /// CIP predicates to astro-float and costs ~0.5–2.5 s, so this runs on demand, not per commit
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
    /// determinism (DNA) requires thread-order independence. Each fixture's result under a
    /// 1-thread pool must equal the default many-thread result, run repeatedly so scheduling
    /// jitter would show.
    ///
    /// **★ Two things this test has to keep honest about itself.**
    ///
    /// First, it passed for the whole stretch when there was *no* parallelism — the cutover
    /// took the old engine's `par_iter` calls with it, and a one-thread pool is trivially
    /// equal to a many-thread one when nothing forks. So the fixtures must have real
    /// parallel width: the fin fold below reaches into the dozens of plane classes, where
    /// the two-ngon fuse has barely a dozen.
    ///
    /// Second, `model_sig` compares coordinates, counts and volumes — **not the report**.
    /// `Notes::sorted` is a stable sort by site, so any two entries sharing a site keep
    /// their *arrival* order, and under `parallel` that is the schedule's to decide. So the
    /// signature includes `BoolReport`, which is the only thing that would catch it.
    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_boolean_is_thread_order_independent() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

        // A hub with fins arrayed around it — every fin is turned by an angle with no exact
        // f64, so every judgement is on the toleranced path and the report is non-empty.
        let fin_fold = |n: i128| -> String {
            let mut m = replay(&[extrude_op(
                Profile2d::polygon(vec![
                    p2(-3.0, -3.0),
                    p2(3.0, -3.0),
                    p2(3.0, 3.0),
                    p2(-3.0, 3.0),
                ]),
                2.0,
            )])
            .unwrap();
            let mut acc = m.live_solids[0];
            let mut sig = String::new();
            for i in 0..n {
                let fin = {
                    let out = ops::apply(
                        &mut m,
                        &Operation::Extrude {
                            plane: SketchPlane::world_xy(),
                            profile: Profile2d::polygon(vec![
                                p2(2.0, -0.4),
                                p2(8.0, -0.4),
                                p2(8.0, 0.4),
                                p2(2.0, 0.4),
                            ]),
                            dist: 1.0,
                        },
                    )
                    .unwrap();
                    m.rebuild_adjacency();
                    match out {
                        OpOutput::Extrude { solid, .. } => solid,
                        o => panic!("{o:?}"),
                    }
                };
                let fin = transform(
                    &mut m,
                    fin,
                    &Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::new(360 * i, n).unwrap()).unwrap(),
                    }),
                )
                .unwrap();
                m.rebuild_adjacency();
                let (solids, report) =
                    boolean_with_report(&mut m, BoolKind::Fuse, acc, fin).unwrap();
                m.rebuild_adjacency();
                acc = solids[0];
                // **The report travels in the signature**, not just the geometry. Measured on
                // this fixture: up to 841 coincidences per boolean, so `loosest`'s `max_by_key`
                // is choosing among hundreds of candidates — which is the tie-break that a
                // schedule could otherwise decide.
                assert!(
                    i == 0 || report.coincidences > 0,
                    "the report is empty, so comparing it proves nothing"
                );
                sig.push_str(&format!("{report:?}"));
            }
            sig.push_str(&model_sig(&m, &[acc]));
            sig
        };

        let build = || {
            let ngon = |n: usize, r: f64, cx: f64, cy: f64| {
                Profile2d::polygon(
                    (0..n)
                        .map(|i| {
                            let ang = std::f64::consts::TAU * (i as f64) / (n as f64);
                            p2(cx + r * ang.cos(), cy + r * ang.sin())
                        })
                        .collect(),
                )
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
        let two_ngons = || {
            let (mut m, a, b) = build();
            let (solids, report) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).unwrap();
            m.rebuild_adjacency();
            format!("{report:?}{}", model_sig(&m, &solids))
        };
        let pool1 = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let seven_fins = || fin_fold(7);
        // **The only fixture here whose alias table is not empty.** `Aliases::record` returns
        // immediately below four planes, so every other model leaves the table at zero and
        // never exercises the snapshot-and-absorb a parallel round is built on. Measured on
        // this one: 36 aliases, settled over two rounds.
        let four_plane = || {
            let (mut m, target, bar) = four_plane_model(0.2, 45);
            let (solids, report) = boolean_with_report(&mut m, BoolKind::Cut, target, bar).unwrap();
            m.rebuild_adjacency();
            format!("{report:?}{}", model_sig(&m, &solids))
        };
        for (name, run) in [
            (
                "two rotated ngons",
                &two_ngons as &(dyn Fn() -> String + Sync),
            ),
            (
                "a seven-fin fold",
                &seven_fins as &(dyn Fn() -> String + Sync),
            ),
            (
                "a four-plane concurrency",
                &four_plane as &(dyn Fn() -> String + Sync),
            ),
        ] {
            let reference = pool1.install(run);
            // A signature that came back empty would make this vacuous whatever the schedule.
            assert!(!reference.is_empty(), "{name}: nothing to compare");
            for _ in 0..4 {
                assert_eq!(
                    run(),
                    reference,
                    "{name}: the result depends on thread order"
                );
            }
        }
    }

    /// A U-prism: a bottom bar `y∈[0,1]` with two prongs rising from it. The prong
    /// tops sit at *different* heights (y=2.3 and y=2.0) on purpose — level tops
    /// would be coplanar faces, which the pre-cutover `has_coplanar_pair` door guard
    /// rejected before the seam machinery ran. That guard is gone; the staggering stays
    /// as this fixture's pinned shape. Area 3 + 1 + 1.3, extruded 1.0 ⇒ volume 5.3.
    fn u_prism() -> (Model, Handle<Solid>) {
        let u = Profile2d::polygon(vec![
            p2(0.0, 0.0),
            p2(3.0, 0.0),
            p2(3.0, 2.3),
            p2(2.0, 2.3),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ]);
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

    /// The L-prism and an L-shaped bar lying in its notch, biting two convex corners of
    /// the L's top face. The bar spans `z ∈ [0.5, 1.5]`, so its body clears the cap.
    ///
    /// Each bite crosses **two different** edges of the cap, which is exactly why no edge
    /// is pierced twice — the bar takes corners, not edges. Two chords, no closed loop.
    fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bar = Profile2d::polygon(vec![
            p2(1.8, 0.8),
            p2(2.1, 0.8),
            p2(2.1, 2.1),
            p2(0.8, 2.1),
            p2(0.8, 1.8),
            p2(1.8, 1.8),
        ]);
        let OpOutput::Extrude { solid: b, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
                profile: bar,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, b)
    }

    /// The L-prism with an **L-shaped** stub standing wholly inside its top face,
    /// `z ∈ [0.5, 1.5]`. The seam on the cap is a closed loop with a reflex node — the
    /// suite's first non-convex inner loop, and the shape a winding must be read from.
    ///
    /// Its coordinates dodge the cap's fan diagonals from `(0,0)` (`y = x`, `y = x/2`,
    /// `y = 2x`), which `segment_crosses_face` would graze.
    fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let ell = Profile2d::polygon(vec![
            p2(0.2, 0.25),
            p2(0.85, 0.25),
            p2(0.85, 0.4),
            p2(0.35, 0.4), // reflex
            p2(0.35, 0.9),
            p2(0.2, 0.9),
        ]);
        let OpOutput::Extrude { solid: stub, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
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
    fn ring_points(planes: &[PlaneGeom], ring: &[[usize; 3]]) -> Vec<[f64; 3]> {
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
        assert_eq!(
            combinatorics::loop_winding(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap(),
            1
        );
        assert_eq!(
            combinatorics::loop_winding(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &hole).unwrap()
            )
            .unwrap(),
            -1
        );

        // Nothing but the ring's direction went into that. Reversing it by hand agrees.
        for (name, ring, want) in [("outer", &outer, -1i8), ("hole", &hole, 1)] {
            let mut reversed = ring.clone();
            reversed.reverse();
            assert_eq!(
                combinatorics::loop_winding(
                    &crate::planes::test_judge(&planes),
                    p,
                    &combinatorics::ring_from_names(p, &reversed).unwrap()
                )
                .unwrap(),
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
        assert_eq!(
            combinatorics::turn_at(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &outer).unwrap(),
                reflex
            )
            .unwrap(),
            -1
        );
        let turns: Vec<i8> = (0..outer.len())
            .map(|i| {
                combinatorics::turn_at(
                    &crate::planes::test_judge(&planes),
                    p,
                    &combinatorics::ring_from_names(p, &outer).unwrap(),
                    i,
                )
                .unwrap()
            })
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
        assert_eq!(
            combinatorics::turn_at(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &outer).unwrap(),
                lo
            )
            .unwrap(),
            1
        );

        // ★ The teeth. A ring is a cycle, so its winding cannot depend on where the walk began.
        // Start it at the reflex node and a `turn_at(ring[0])` implementation reads the reflex
        // sign — the exact fault `outer_tri` shipped. Measured: without this rotation, such an
        // implementation passes every assertion above.
        let mut rotated = outer.clone();
        rotated.rotate_left(reflex);
        assert_eq!(
            combinatorics::turn_at(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &rotated).unwrap(),
                0
            )
            .unwrap(),
            -1
        );
        assert_eq!(
            combinatorics::loop_winding(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &rotated).unwrap()
            )
            .unwrap(),
            1
        );
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
    /// tripped the pre-cutover door guard, exactly as `u_prism`'s staggered prongs avoid.
    /// And the legs span `y ∈ [0.65, 1.3]`, not `[0.7, 1.3]`, because `(1.4, 0.7)` lies on
    /// the cap's fan diagonal `y = x/2` and `segment_crosses_face` would graze it.
    fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let staple = Profile2d::polygon(vec![
            p2(0.1, 0.5),
            p2(0.6, 0.5),
            p2(0.6, 1.3),
            p2(0.8, 1.3),
            p2(0.8, 0.45),
            p2(1.4, 0.45),
            p2(1.4, 1.5),
            p2(0.1, 1.5),
        ]);
        let OpOutput::Extrude { solid: st, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::from_origin_normal(
                    Point3::from_array([0.0, 1.3, 0.0]),
                    Vector3::from_array([0.0, -1.0, 0.0]),
                )
                .expect("a unit normal"),
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
    /// `"ell"` an L-shaped one. Both rings come from [`combinatorics::face_vertex_triples`] and
    /// [`combinatorics::hole_rings`], which the boolean itself uses, so the fixture exercises only code
    /// the kernel runs.
    ///
    /// This replaces a helper that built its rings from `seam_paths_on`/`orient_seam_loop` — the
    /// retired seam engine. The *properties* below are about `point_in_ring`/`every_ray`, which are
    /// live and load-bearing (`nest_cells` picks a hole's host with them, `unify_coplanar_faces`
    /// groups by them), so they had to be re-homed rather than deleted with their old fixture.
    fn holed_face_rings(which: &str) -> (Vec<PlaneGeom>, usize, Ring, Ring) {
        let (mut m, l, stub) = if which == "dimple" {
            l_and_dimple()
        } else {
            l_and_ell_stub()
        };
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).expect("the cut");
        m.rebuild_adjacency();
        let faces_tab = collect_planes(&m, r).unwrap();
        let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
        for (i, pi) in faces_tab.iter().enumerate() {
            surf_ix.insert(pi.face.expect("a real face table"), i);
        }
        let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
        let (planes, plane_ix) = dense_planes(&faces_tab, &canon);
        let inc = combinatorics::edge_faces(&m, r, &surf_ix).unwrap();
        for &fh in &m.shells.get(m.solids.get(r).outer).faces {
            let fp = surf_ix[&fh];
            let holes = combinatorics::hole_rings(
                &m,
                fh,
                fp,
                &inc,
                &crate::planes::test_judge(&planes),
                &plane_ix,
            )
            .unwrap();
            if let Some(hole) = holes.into_iter().next() {
                let outer = combinatorics::face_vertex_triples(
                    &m,
                    fh,
                    fp,
                    &inc,
                    &crate::planes::test_judge(&planes),
                    &plane_ix,
                )
                .unwrap();
                assert_eq!(outer.len(), 6, "{which}: the L's cap is a reflex hexagon");
                return (planes, plane_ix[fp], outer, hole);
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
            assert!(
                combinatorics::point_in_ring(
                    &crate::planes::test_judge(&planes),
                    p,
                    *t,
                    &combinatorics::ring_from_names(p, &outer).unwrap()
                )
                .unwrap()
            );
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
                let rays = combinatorics::every_ray(
                    &crate::planes::test_judge(&planes),
                    p,
                    *t,
                    &combinatorics::ring_from_names(p, &outer).unwrap(),
                )
                .unwrap();
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
            assert!(
                combinatorics::point_in_ring(
                    &crate::planes::test_judge(&planes),
                    p,
                    *t,
                    &combinatorics::ring_from_names(p, &outer).unwrap()
                )
                .unwrap()
            );
            assert!(
                combinatorics::point_in_ring(
                    &crate::planes::test_judge(&planes),
                    p,
                    *t,
                    &combinatorics::ring_from_names(p, &rev_outer).unwrap()
                )
                .unwrap(),
                "reversing the outer ring must not move the hole"
            );
        }
        for t in &outer {
            assert!(
                !combinatorics::point_in_ring(
                    &crate::planes::test_judge(&planes),
                    p,
                    *t,
                    &combinatorics::ring_from_names(p, &hole).unwrap()
                )
                .unwrap()
            );
            assert!(
                !combinatorics::point_in_ring(
                    &crate::planes::test_judge(&planes),
                    p,
                    *t,
                    &combinatorics::ring_from_names(p, &rev_hole).unwrap()
                )
                .unwrap(),
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
        let rays = combinatorics::every_ray(
            &crate::planes::test_judge(&planes),
            p,
            hole[0],
            &combinatorics::ring_from_names(p, &outer).unwrap(),
        )
        .unwrap();
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
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                hole[0],
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
        assert!(
            !combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                outer[0],
                &combinatorics::ring_from_names(p, &hole).unwrap()
            )
            .unwrap()
        );
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

    /// **`Origin` no longer tells result faces apart.** The arrangement names every vertex
    /// it emits by the three planes meeting there, so an operand corner the cut never touched
    /// comes back as `Discovered`, exactly like a seam vertex. Nothing carries over as
    /// `Constructed` (the arrangement builds no vertex from an original handle).
    ///
    /// This is a contract, not a curiosity: `pipeline.rs`'s island test selected a face by
    /// "all its vertices are `Discovered`", which was unique under the old engine and is
    /// now true of *every* face. It flipped the wrong face and only the last assertion
    /// noticed. Selecting a face by provenance is what this locks out.
    ///
    /// The subject is `Cut(l, stub)` — the **holed** result, so `face_half_edges` walks
    /// `inner` rings too (`count_discovered` walks only `outer` and would miss them).
    ///
    /// **Unrotated only.** Rotating a boolean result re-marks these vertices `Rotated` over
    /// a `Discovered` base — that is `transform_rotate_boolean_result_keeps_discovered_base`,
    /// and this lock must not be read as contradicting it.
    #[test]
    fn an_unrotated_boolean_names_every_vertex_by_its_plane_triple() {
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();

        let mut seen = std::collections::HashSet::new();
        let mut holed = 0;
        for sh in solid_shell_handles(&m, r) {
            for &fh in &m.shells.get(sh).faces {
                let face = m.faces.get(fh);
                holed += usize::from(!face.inner.is_empty());
                for he in face_half_edges(face) {
                    for vh in m.edges.get(he.edge).bounds.into_iter().flatten() {
                        if !seen.insert(vh) {
                            continue;
                        }
                        assert!(
                            matches!(
                                m.vertices.get(vh).origin,
                                Origin::Discovered {
                                    definition: VertexDef::ThreePlane(_),
                                    ..
                                }
                            ),
                            "vertex {:?} is {:?}, not a plane triple",
                            m.vertices.get(vh).point.as_array(),
                            m.vertices.get(vh).origin
                        );
                    }
                }
            }
        }
        assert_eq!(holed, 1, "the blind dimple leaves exactly one holed face");
        assert_eq!(seen.len(), 20, "the L's 12 corners + the dimple's 8");
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

    /// A sever that also leaves a surviving cavity: a hollow box whose void sits to one side,
    /// cut by a slab that severs it without touching the void. The x<2 piece keeps the void as a
    /// cavity, the x>2 piece is solid — two outward shells *and* one inward. `point_in_component`
    /// assigns the void to the x<2 piece that nests it (a containment test), rather than rejecting.
    #[test]
    fn severed_with_cavity_assigns_the_void_to_its_piece() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        // Void near the x-low side (1×2×2 = 4), clear of the x=2 cut.
        let inner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 2.5]),
        );
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        // A slab spanning full y,z, thin in x at x∈[2,2.2] — severs into x<2 (holds the void, vol
        // 2·3·3 − 4 = 14) and x>2 (solid, vol 0.8·3·3 = 7.2).
        let slab = m.add_cuboid(
            Point3::from_array([2.0, -1.0, -1.0]),
            Point3::from_array([2.2, 4.0, 4.0]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
        assert_eq!(
            solids.len(),
            2,
            "the slab severs the hollow box into two pieces"
        );
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // Exactly one piece owns the void; volumes match the hand calculation.
        let with_cav: Vec<_> = solids
            .iter()
            .filter(|&&s| !m.solids.get(s).cavities.is_empty())
            .collect();
        assert_eq!(
            with_cav.len(),
            1,
            "the void is assigned to exactly one piece"
        );
        let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
        let hollow_piece = *with_cav[0];
        assert!(
            (vol(hollow_piece) - 14.0).abs() < 1e-9,
            "hollow piece {}",
            vol(hollow_piece)
        );
        let total: f64 = solids.iter().map(|&s| vol(s)).sum();
        assert!((total - 21.2).abs() < 1e-9, "total {total}");
    }

    /// The adjacent case: a cut that passes *through* the void opens it — the void wall becomes
    /// exterior boundary, so no cavity survives. Handled by the plain sever path (no cavity to
    /// assign), not the containment code, but pinned so a regression there is caught.
    #[test]
    fn a_cut_through_the_void_leaves_no_cavity() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3])); // void 1³
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        // Slab x∈[1.4,1.6] passes through the void (x∈[1,2]) → severs AND opens the void.
        let slab = m.add_cuboid(
            Point3::from_array([1.4, -1.0, -1.0]),
            Point3::from_array([1.6, 4.0, 4.0]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // The void is opened, so neither piece keeps a cavity; material = 26 − (1.8 − 0.2) = 24.4.
        let total_cavities: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
        assert_eq!(
            total_cavities, 0,
            "the cut opened the void — no surviving cavity"
        );
        let total: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((total - 24.4).abs() < 1e-9, "total {total}");
    }

    /// Nested cavities: a hollow box A ([0,6]³ − [1,5]³ void) with a smaller hollow box B
    /// ([2,4]³ − [2.5,3.5]³ void) floating inside A's void. `Fuse(A,B)` is one arrangement with
    /// four components (two materials, two voids); B's void is contained by **both** A's outer
    /// shell and B's own, so the containment assignment must pick the **innermost** (B), not A.
    /// The result is two solids, each keeping its own void (A: 216−64 = 152, B: 8−1 = 7).
    #[test]
    fn a_void_nested_in_a_floating_island_goes_to_the_inner_solid() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([6.0; 3]));
        let void = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([5.0; 3]));
        let a = boolean_one(&mut m, BoolKind::Cut, big, void).unwrap();
        m.rebuild_adjacency();
        let bbig = m.add_cuboid(Point3::from_array([2.0; 3]), Point3::from_array([4.0; 3]));
        let bvoid = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
        let b = boolean_one(&mut m, BoolKind::Cut, bbig, bvoid).unwrap();
        m.rebuild_adjacency();
        let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // Each solid keeps exactly one void — B's void was assigned to B (innermost), not A.
        let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
        for &s in &solids {
            assert_eq!(
                m.solids.get(s).cavities.len(),
                1,
                "each piece keeps its own void"
            );
        }
        let mut vols: Vec<f64> = solids.iter().map(|&s| vol(s)).collect();
        vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((vols[0] - 7.0).abs() < 1e-9, "inner {}", vols[0]);
        assert!((vols[1] - 152.0).abs() < 1e-9, "outer {}", vols[1]);
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
        let far = SketchPlane::world_xy().with_origin(Point3::from_array([5.0, 0.0, 0.0]));
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

    /// The `0.4` square on `[0.3, 0.7]²` of the unit cube's lid.
    ///
    /// ★ **On a lid these are world coordinates.** The sketch origin is the world origin projected
    /// onto the plane, and the arbitrary-axis convention gives `n = ẑ` the axes `u = +x̂, v = +ŷ`,
    /// so a frame point `(a, b)` is world `(a, b, 1)`.
    fn small_square() -> Profile2d {
        Profile2d::polygon(vec![p2(0.3, 0.7), p2(0.3, 0.3), p2(0.7, 0.3), p2(0.7, 0.7)])
    }

    /// **Chaining onto a fused boss.** The fuse leaves the base's `z=1` face a *ring* — a face with
    /// a hole where the boss sits — and the second boolean cuts through both. Every plane class the
    /// cut opens then meets that ring along the **hole's own edge**, which is the case that used to
    /// label inconsistently: the ring's neighbouring vertices there point *into* the hole, so
    /// reading the occupied side off a flank put the material on the wrong side of `W`. The side
    /// now comes from the ring's travel ([`arrangement::run_body_above`]), and the run leaves as its own
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

    /// explicit sharing (overhaul #3): a prism built with a shared base-cap
    /// surface reuses that `Surface` handle for its flush cap, and reconciles the
    /// cap's face orientation so the materialized outward normal stays `−sweep`.
    #[test]
    fn build_prism_base_cap_reuses_shared_surface() {
        let mut m = Model::new();
        // A face-plane surface with outward normal +z (as a face on the base solid).
        let sf = m.push_surface(
            Surface::Plane(
                Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
                    .unwrap(),
            ),
            nacre_topo::SurfaceDef::Constructed,
        );
        let base_pts = [
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        ];
        let (_prism, faces) = build_prism(
            &mut m,
            crate::exact::Swept::along(base_pts.to_vec(), Vector3::from_array([0.0, 0.0, 1.0])),
            vec![],
            Vector3::from_array([0.0, 0.0, 1.0]),
            Some(sf),
            None,
        )
        .unwrap();
        let cap = m.faces.get(faces[0]); // base cap is pushed first
        // Shared handle (was a fresh push before overhaul #3).
        assert_eq!(cap.surface, sf, "base cap reuses the shared surface handle");
        // Orientation reconciled: materialized outward normal is −sweep (−z).
        let Surface::Plane(p) = m.surface(cap.surface) else {
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

    /// ★★★★★ **A plane's record and its `SurfaceDef` are one statement.**
    ///
    /// `Constructed` means *"these points speak about the world"*, `Moved` means *"about the
    /// pre-motion frame"* — so a base cap that takes the caller's world triple must be
    /// `Constructed`, and one that falls back to the prism's own ring must take the frame's `def`
    /// with it. There is nothing else to keep in step: the plane's canonical name is derived from
    /// whichever triple is recorded.
    ///
    /// ★★ **It used to be three halves and they could come apart.** Coefficients were supplied
    /// beside the points and chosen by a **separate** `or_else`, so a caller whose points
    /// overflowed while their coefficients did not got world coefficients recorded beside the
    /// prism's own frame ring — one plane stated two ways. The agreement filter could not catch it
    /// either, since `c · p` overflows at exactly those widths. Removing the coefficient parameter
    /// is what made the pairing structural; this pins what is left of the choice.
    #[test]
    fn a_prisms_base_cap_records_the_frame_its_def_names() {
        let r = nacre_scalar::Rat::from_int;
        let base_pts = [
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        ];
        let prism = |pts: Option<[[nacre_scalar::Rat; 3]; 3]>| {
            let mut m = Model::new();
            let (_prism, faces) = build_prism(
                &mut m,
                crate::exact::Swept::along(base_pts.to_vec(), Vector3::from_array([0.0, 0.0, 1.0])),
                vec![],
                Vector3::from_array([0.0, 0.0, 1.0]),
                None,
                pts,
            )
            .unwrap();
            let surf = m.faces.get(faces[0]).surface; // base cap is pushed first
            (m, surf)
        };

        // ★ The caller's triple, when there is one — here a plane nowhere near the prism, so a
        // record that ignored it would be visibly different rather than coincidentally equal.
        let far_pts = [[r(0), r(0), r(3)], [r(1), r(0), r(3)], [r(0), r(1), r(3)]];
        let (m, surf) = prism(Some(far_pts));
        assert_eq!(
            m.surface_points.get(&surf),
            Some(&far_pts),
            "the caller's triple was not the one recorded"
        );
        assert_eq!(
            m.surface_name.get(&surf),
            Some(&nacre_scalar::PlaneName::Narrow([r(0), r(0), r(1), r(-3)])),
            "the name was not derived from the triple that was recorded"
        );
        assert_eq!(
            m.surface_defs.get(&surf),
            Some(&nacre_topo::SurfaceDef::Constructed)
        );

        // ★ And with no caller statement, the ring answers — `Swept::along` carries no rationals,
        // so there is nothing to record and no name to derive.
        let (m, surf) = prism(None);
        assert!(!m.surface_points.contains_key(&surf));
        assert!(!m.surface_name.contains_key(&surf));
    }

    /// The handle branch of `shares_or_coplanar` is load-bearing: a shared
    /// `Surface` handle reports coplanar even when the stored `plane` values are
    /// *not* geometrically coplanar (so the fallback would not fire). This is the
    /// path a referenced coplanar contact takes; on axis-aligned M5 it is redundant
    /// with the geometric test, but the branch must work for rotated frames.
    #[test]
    fn shares_or_coplanar_uses_the_handle_branch() {
        let mut m = Model::new();
        let shared = m.push_surface(
            Surface::Plane(
                Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0]))
                    .unwrap(),
            ),
            nacre_topo::SurfaceDef::Constructed,
        );
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
        let mk = |plane, tri: [Point3; 3]| FaceInfo {
            // Unmoved and hand-built: nothing to record, and the base frame is unused anyway.
            base_rat: None,
            motion: None,
            surf: shared,
            face: Some(fh),
            plane,
            tri,
            n_out: Vector3::from_array([0.0; 3]),
            // Unread: this table only ever reaches `Judge::planes_coplanar`, which decides on `tri`.
            orient_sign: 1,
            tri_pt3: tri.map(|p| nacre_cip::Pt3::exact(p.as_array()).expect("exact")),
            rotated: false,
        };
        let p = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
        let planes = vec![
            mk(plane_x0, [p(0., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]), // in x = 0
            mk(plane_z0, [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.)]), // in z = 0
        ];
        // Neither fallback fires: the coefficients are not proportional, and the coordinates say
        // these really are two different planes.
        assert!(!planes_coplanar(&planes[0].plane, &planes[1].plane));
        assert!(!crate::planes::test_judge(&planes).planes_coplanar(0, 1));
        // The shared handle alone makes them coplanar-by-reference.
        assert!(shares_or_coplanar(
            &crate::planes::test_judge(&planes),
            0,
            1
        ));
    }

    fn pocket_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
        Operation::PocketOnFace {
            face,
            profile,
            dist,
        }
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
            let big = Profile2d::polygon(vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)]);
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
            let rect = Profile2d::polygon(vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)]);
            let mut m = Model::new();
            let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: rect,
                dist: sz,
            }).unwrap() else { unreachable!() };
            let hole = Profile2d::polygon(vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)]);
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
            let rect = Profile2d::polygon(vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)]);
            let mut m = Model::new();
            let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: rect,
                dist: sz,
            }).unwrap() else { unreachable!() };
            let hole = Profile2d::polygon(vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)]);
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

    /// The reported four-plane model: a unit cube with a block fused on its top, and a bar spun
    /// `deg`° about an axis in the block's `x = 0.5` plane. At 45° with `half_z == 0.2` the bar's
    /// half-width and its pivot-to-bottom offset are equal, so its bottom corner edge lands in that
    /// plane and the y-planes cutting the edge become four-plane vertices.
    fn four_plane_model(half_z: f64, deg: i128) -> (Model, Handle<Solid>, Handle<Solid>) {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let block = m.add_cuboid(
            Point3::from_array([0.5, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let target = boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses on")[0];
        m.rebuild_adjacency();
        let bar = m.add_cuboid(
            Point3::from_array([0.3, -0.5, 1.0 - half_z]),
            Point3::from_array([0.7, 1.5, 1.0 + half_z]),
        );
        m.rebuild_adjacency();
        let bar = transform(
            &mut m,
            bar,
            &Isometry::rotation(Rotation {
                axis: Axis::Y,
                point: [Rat::new(1, 2).unwrap(), Rat::from_int(0), Rat::from_int(1)],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            }),
        )
        .unwrap();
        m.rebuild_adjacency();
        (m, target, bar)
    }

    /// **What the corpus actually contains in the way of concurrent vertices — and whether the
    /// engine's rule for noticing them is complete.**
    ///
    /// The four-plane work rests on one premise: a point's identity can be made a function of the
    /// *set* of planes through it, because **every producer that meets the point derives the same
    /// set**. Until now that was checked on a single model. This measures it against ground truth
    /// (every plane asked, not the rule asking itself) over the whole fixture corpus.
    ///
    /// Two things are asserted, and the second is the load-bearing one:
    ///
    /// 1. **Exactly four.** No corpus point has five or more planes through it. That matters
    ///    because it is what makes both discovery rules complete: the trace learns `{wc} ∪ t`, and
    ///    with `|S| = 4` that *is* `S`. A five-plane point would leave it one short — so if this
    ///    ever fires, the identity rule needs the full set from somewhere else, and the message
    ///    says which model found it.
    /// 2. **The trace's rule reproduces ground truth.** Wherever the trace would notice a
    ///    concurrency (a name on class `wc` that does not mention `wc`), the set it would record
    ///    equals the set every plane agrees on.
    ///
    /// The count is printed rather than pinned: this is a *measurement*, and a number here would
    /// only pin today's fixture list.
    #[test]
    fn concurrent_vertices_are_four_planes_and_the_trace_sees_all_of_them() {
        let mut boxed: Vec<(String, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
        macro_rules! fixture {
            ($name:ident) => {{
                let (m, a, b) = $name();
                boxed.push((stringify!($name).to_string(), m, a, b));
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
        // ...and the same shapes tilted, since a rotated operand is where concurrencies actually
        // turn up: an axis-aligned corpus would measure the easy half and call it the whole.
        macro_rules! tilted {
            ($name:ident) => {{
                let (mut m, a, b) = $name();
                let iso = rot_iso(nacre_scalar::Axis::Z, 30);
                let a = transform(&mut m, a, &iso).unwrap();
                m.rebuild_adjacency();
                let b = transform(&mut m, b, &iso).unwrap();
                m.rebuild_adjacency();
                boxed.push((format!("{} (tilted)", stringify!($name)), m, a, b));
            }};
        }
        tilted!(two_boxes);
        tilted!(nested_boxes);
        tilted!(cube_and_notch);
        tilted!(stacked_cubes);
        tilted!(l_and_corner_box);
        tilted!(l_and_reflex_box);
        tilted!(l_and_inner_box);
        tilted!(l_and_popup_box);
        tilted!(l_and_notch_bar);
        tilted!(l_and_ell_stub);
        tilted!(l_and_staple);
        tilted!(l_and_dimple);
        tilted!(l_and_rod);
        tilted!(u_and_slab);

        // ★ And the models that actually have one. Without these the sweep asserts nothing: the
        // corpus above turns out to carry **no** concurrent vertex at all, so it can measure the
        // blast radius of the coming stages but not the discovery rule. These are the reported
        // model (a bar spun 45° whose bottom corner edge lands in the block's x = 0.5 plane) and
        // its variants; `rejects.rs` pins the same shape from outside.
        for (label, half_z, deg) in [
            ("four_plane (the reported model)", 0.2, 45),
            ("four_plane at 44 deg (a near miss)", 0.2, 44),
            ("four_plane, taller bar (a near miss)", 0.3, 45),
        ] {
            let (m, t, bar) = four_plane_model(half_z, deg);
            boxed.push((label.to_string(), m, t, bar));
        }

        let (mut with_any, mut total, mut trace_rule_checked) = (0usize, 0usize, 0usize);
        let mut lines_found = 0usize;
        for (name, m, a, b) in &boxed {
            let found = arrangement::concurrency_audit(m, *a, *b).unwrap();
            if !found.is_empty() {
                with_any += 1;
            }
            total += found.len();
            for c in &found {
                lines_found += c.lines.len();
                assert_eq!(
                    c.planes.len(),
                    4,
                    "{name}: a {}-plane point at {:?} — the trace learns only `{{wc}} ∪ t`, \
                     which is four names at most, so this one would be discovered incomplete",
                    c.planes.len(),
                    c.planes
                );
                if !c.triple.contains(&c.wc) {
                    trace_rule_checked += 1;
                    let mut derived = c.triple.to_vec();
                    derived.push(c.wc);
                    derived.sort_unstable();
                    assert_eq!(
                        derived, c.planes,
                        "{name}: on class {} the trace would record {derived:?} for the point \
                         named {:?}, but every plane says {:?}",
                        c.wc, c.triple, c.planes
                    );
                }
            }
        }
        println!(
            "concurrency audit: {with_any}/{} fixtures carry a concurrent vertex, {total} in total",
            boxed.len()
        );
        // A sweep that quietly stops finding anything reads as agreement, so say what was actually
        // exercised — and fail if the load-bearing assertion never ran.
        assert!(
            trace_rule_checked > 0,
            "no observation reached the trace's discovery condition, so the rule was not measured"
        );
        // The same discovery must also yield the *line* aliases — three planes sharing a line show
        // up as a sub-triple of `S` that names no point. In this corpus every concurrency comes
        // from exactly that (a tool edge lying in a target plane), so each one carries one.
        assert!(
            lines_found > 0,
            "no line-sharing triple was derived, yet every concurrency here comes from one"
        );
        println!("  and {lines_found} carried a line-sharing triple");
        println!("  of which {trace_rule_checked} exercised the trace's discovery rule");
    }

    /// **Every producer states its side in the label frame** — the invariant family #2 restored,
    /// swept over the whole two-solid corpus (prints the interesting classes with `--nocapture`).
    ///
    /// The arrangement states its cell labels as `[*_above, *_below]` about one direction per plane
    /// class: the class root's **stored surface normal** (`Seated{body_above}` and `emit_faces`'
    /// `flip` are written against it). `combinatorics::side_of` answers in the root's **outward** frame
    /// instead, and the two are opposite exactly when the root face is `Reversed`
    /// (`orient_sign == -1`) — which no `add_cuboid` face ever is, but a face an earlier boolean
    /// re-emitted flipped is. `graze_above` read `side_of` raw, so on a pocket wall it flipped the
    /// wrong label bit. It no longer reads a point's side at all — [`arrangement::run_body_above`] derives
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
            let audits = arrangement::frame_audit(m, BoolKind::Cut, *a, *b).unwrap();
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
        // handling of `arrangement::run_body_above` is exercised against a crossed frame: the three
        // `false` entries below are the pre-existing answers, unchanged by the new rule.
        let (m, pc, bx) = boxed
            .iter()
            .find_map(|(n, m, a, b)| (*n == "pocket_corner_cut").then_some((m, *a, *b)))
            .unwrap();
        let wall = arrangement::frame_audit(m, BoolKind::Cut, pc, bx)
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

    /// Transforming a boolean *result* (which carries `Discovered` seam vertices) does not
    /// **downgrade** those vertices to `Constructed` — the failure this guards.
    ///
    /// A motion that records a forest node supersedes the origin with `Moved { base, motion }`,
    /// and the definition is preserved *through the base*: the base is the seam vertex, still in
    /// the arena with its `ThreePlane` definition, and the node says how it moved. A motion that
    /// records nothing instead remaps the definition's planes in place. Either way the truth
    /// survives; only its spelling depends on whether the motion was worth recording.
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
        // Every seam vertex still names its three-plane definition — directly, or through the
        // base of the motion that moved it. None fell back to `Constructed`.
        let named = boundary_verts(&m, r2)
            .into_iter()
            .filter(|&vh| {
                let root = match m.vertices.get(vh).origin {
                    Origin::Moved { base, .. } => m.vertices.get(base).origin,
                    other => other,
                };
                matches!(root, Origin::Discovered { .. })
            })
            .count();
        assert_eq!(named, disc, "seam definitions preserved");
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

    /// The bit-identity guard, on the one producer that can break it.
    ///
    /// A definition and its cached coordinate must agree **exactly**, and they do only because the
    /// replay performs the very same float operations the producer did. A reflection is the step
    /// where that is easiest to lose — `Pt3::mirror` must walk the same `2c − x` that
    /// `AxisMirror::point` just walked — and a chain that ends in one is what this checks.
    #[test]
    fn a_mirrored_rotated_vertex_reconstructs_from_its_definition() {
        use nacre_scalar::{Axis, Rat};
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let r = transform(&mut m, c, &rot30()).unwrap();
        m.rebuild_adjacency();
        let mirrored = crate::transform::mirror(&mut m, r, Axis::X, Rat::from_int(0)).unwrap();
        m.rebuild_adjacency();

        let shell = m.solids.get(mirrored).outer;
        let mut checked = 0;
        for &fh in &m.shells.get(shell).faces.clone() {
            for he in &m.faces.get(fh).outer.half_edges.clone() {
                for vh in m.edges.get(he.edge).bounds.iter().flatten() {
                    assert!(
                        matches!(m.vertices.get(*vh).origin, Origin::Moved { .. }),
                        "the image keeps its rotation definition"
                    );
                    let Origin::Moved {
                        base,
                        motion: rotation,
                    } = m.vertices.get(*vh).origin
                    else {
                        unreachable!("just asserted Rotated")
                    };
                    let replayed = crate::rotated_vertex::replay_chain_coord(
                        &m,
                        m.vertices.get(base).point.as_array(),
                        rotation,
                    )
                    .expect("the root coordinate lifts to an exact rational");
                    assert_eq!(
                        replayed,
                        m.vertices.get(*vh).point.as_array(),
                        "definition reproduces the stored coordinate bit for bit"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "walked no vertices");
    }

    /// Replay determinism (DNA 3) for the one additive operation: a copy reproduces the same
    /// geometry *and* the same handle index, and leaves the same live set behind it.
    #[test]
    fn copy_is_deterministic() {
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let twin = crate::transform::copy(&mut m, c).unwrap();
            (bbox_lo(&m, twin), twin, m.live_solids.clone())
        };
        assert_eq!(build(), build(), "same ops → same geometry and handles");
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
    /// vertices carry `Origin::Moved` (`solid_is_rotated`), and a boolean against it now
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

        assert!(solid_is_rotated(&m, c2), "vertices carry Origin::Moved");

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
                        if let Origin::Moved { base, .. } = m.vertices.get(vh).origin {
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

    /// Whether any wall of `s` is a **rotated image** — what `planes::solid_is_rotated` used to
    /// ask of the vertices, now asked of the surfaces that actually record it.
    fn solid_is_rotated(m: &Model, s: Handle<Solid>) -> bool {
        m.shells.get(m.solids.get(s).outer).faces.iter().any(|&fh| {
            matches!(
                m.surface_defs.get(&m.faces.get(fh).surface),
                Some(nacre_topo::SurfaceDef::Moved { .. })
            )
        })
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
        let Origin::Moved {
            base,
            motion: rotation,
        } = m.vertices.get(vh).origin
        else {
            return None;
        };
        let base_is_rotated = matches!(m.vertices.get(base).origin, Origin::Moved { .. });
        let mut axes = Vec::new();
        let mut cur = Some(rotation);
        while let Some(h) = cur {
            let n = m.motion(h);
            if let nacre_topo::Motion::Rotate { axis, .. } = n.motion {
                axes.push(axis);
            }
            cur = n.parent;
        }
        axes.reverse();
        Some((base_is_rotated, axes.len(), axes))
    }

    /// **A chained boolean must not lose exactness.**
    ///
    /// A rotated solid's face coordinates are rounded, so its planes are truthful only through an
    /// exact *definition* (`FaceInfo::tri_pt3` built from a rotation history). A boolean's *result*
    /// is just as rotated as its operands — but the result carries no rotation provenance, so
    /// `collect_planes` describes every one of its faces by `Pt3::exact` of the rounded triangle
    /// and the kernel starts treating a rounded copy as the truth. That is what makes one wall
    /// become two plane classes on the next operation.
    ///
    /// The invariant: **every face of a boolean between rotated operands is described by a
    /// rotation definition, not by its rounded coordinates.**
    #[test]
    fn a_chained_boolean_keeps_its_faces_exact() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, 0.0]),
            Point3::from_array([1.0, 1.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, -0.2, 1.0]),
            Point3::from_array([4.0, 0.2, 3.0]),
        );
        m.rebuild_adjacency();
        let a = transform(&mut m, a, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let b = transform(&mut m, b, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        // The operands hold up: a rotated solid's faces do carry definitions.
        for (name, s) in [("operand a", a), ("operand b", b)] {
            let planes = crate::planes::collect_planes(&m, s).expect("planes");
            assert!(
                planes.iter().all(|f| f.rotated),
                "{name}: a rotated operand's faces must be described by their rotation"
            );
        }
        let r = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse")[0];
        m.rebuild_adjacency();
        let planes = crate::planes::collect_planes(&m, r).expect("planes");
        let described = planes.iter().filter(|f| f.rotated).count();
        assert_eq!(
            described,
            planes.len(),
            "the result of a rotated boolean is rotated too: {described}/{} faces carry a \
             definition, the rest are rounded coordinates declared exact",
            planes.len()
        );
    }

    /// **The same motion, applied twice, is the same node.**
    ///
    /// The motion handle is the canonical name of "which motion", and judgments use it to decide
    /// whether a whole judgement can be answered exactly in the pre-motion frame. Without
    /// interning, turning two solids by the same 30° makes two nodes, their shared motion stops
    /// cancelling, and rotating a model turns its exact questions into assumed ones — which is
    /// what `a_shared_rotation_still_assumes_nothing` measures. (The identity it replaced was a
    /// 64-bit hash of the chain's contents, where a collision would have answered a *different*
    /// question with full confidence.)
    #[test]
    fn the_same_motion_applied_twice_is_one_node() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([2.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let a = transform(&mut m, a, &rot_iso(Axis::X, 30)).unwrap();
        m.rebuild_adjacency();
        let b = transform(&mut m, b, &rot_iso(Axis::X, 30)).unwrap();
        m.rebuild_adjacency();
        fn leaf(m: &Model, s: Handle<Solid>) -> Handle<nacre_topo::MotionNode> {
            let sh = m.solids.get(s).outer;
            let fh = m.shells.get(sh).faces[0];
            let vh = m
                .edges
                .get(m.faces.get(fh).outer.half_edges[0].edge)
                .bounds
                .unwrap()[0];
            match m.vertices.get(vh).origin {
                Origin::Moved { motion, .. } => motion,
                other => panic!("a rotated solid's vertices are moved, got {other:?}"),
            }
        }
        assert_eq!(leaf(&m, a), leaf(&m, b), "one motion, one node");
        // …and a *different* motion is a different node, or the identity would be worthless.
        let c = m.add_cuboid(
            Point3::from_array([5.0, 0.0, 0.0]),
            Point3::from_array([6.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let c = transform(&mut m, c, &rot_iso(Axis::X, 31)).unwrap();
        m.rebuild_adjacency();
        assert_ne!(
            leaf(&m, a),
            leaf(&m, c),
            "different motions must not share a node"
        );
    }

    /// **A boolean result rotated again continues its history — per wall.**
    ///
    /// One solid does not have one rotation history. A result's vertices are all `Discovered`, so
    /// asking them "what rotation is this solid at" answers `None` and the next rotation would
    /// start a fresh root — replaying a pre-first-rotation witness through only the *second*
    /// rotation, which is a plane that does not exist. And its walls can come from operands
    /// rotated by different angles, so there is no single answer to give. Each surface therefore
    /// chains from its own leaf.
    #[test]
    fn a_rerotated_boolean_result_continues_each_walls_history() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        );
        m.rebuild_adjacency();
        // Different angles, so the two operands' walls carry genuinely different histories.
        let a = transform(&mut m, a, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let b = transform(&mut m, b, &rot_iso(Axis::Z, 50)).unwrap();
        m.rebuild_adjacency();
        let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap()[0];
        m.rebuild_adjacency();
        let r = transform(&mut m, r, &rot_iso(Axis::Z, 20)).unwrap();
        m.rebuild_adjacency();

        let mut leaves = std::collections::HashSet::new();
        let sh = m.solids.get(r).outer;
        for &fh in &m.shells.get(sh).faces {
            let s = m.faces.get(fh).surface;
            let nacre_topo::SurfaceDef::Moved {
                motion: rotation, ..
            } = m
                .surface_defs
                .get(&s)
                .copied()
                .expect("every face's surface is defined")
            else {
                panic!("a rotated result's walls must carry a rotation");
            };
            assert_eq!(
                crate::rotated_vertex::motion_chain(&m, rotation)
                    .expect("an axis-aligned history holds no frame")
                    .len(),
                2,
                "both rotations, once each"
            );
            leaves.insert(rotation);
        }
        assert_eq!(leaves.len(), 2, "the two operands' histories stay apart");
    }

    /// **A copy of a rotated solid is still exactly defined.**
    ///
    /// `copy` is `transform` under a *zero* translation, so a surface rule that reads "any
    /// translation makes a rotated plane inexpressible" swallows it — and then a copy of a rotated
    /// solid cannot take part in a boolean at all, though its geometry is bit-identical to the
    /// original's. The identity is not a translation.
    #[test]
    fn a_copy_of_a_rotated_solid_answers_like_the_original() {
        use nacre_scalar::Axis;
        let build = |use_copy: bool| {
            let mut m = Model::new();
            let hub = m.add_cuboid(
                Point3::from_array([-1.0, -1.0, 0.0]),
                Point3::from_array([1.0, 1.0, 3.0]),
            );
            let fin = m.add_cuboid(
                Point3::from_array([0.5, -0.2, 1.0]),
                Point3::from_array([4.0, 0.2, 3.0]),
            );
            m.rebuild_adjacency();
            let mut fin = transform(&mut m, fin, &rot_iso(Axis::Z, 30)).unwrap();
            m.rebuild_adjacency();
            if use_copy {
                fin = crate::transform::copy(&mut m, fin).unwrap();
                m.rebuild_adjacency();
            }
            let r = boolean(&mut m, BoolKind::Fuse, hub, fin).expect("fuse");
            nacre_props::mass_props(&m, r[0]).expect("props").volume
        };
        // Bit-identical, not merely close: the copy walks the same definitions.
        assert_eq!(build(true), build(false));
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
        let (boss, _) = build_prism(
            &mut m,
            crate::exact::Swept::along(l_base, Vector3::from_array([0.0, 0.0, 0.4])),
            vec![],
            Vector3::from_array([0.0, 0.0, 1.0]),
            None,
            None,
        )
        .unwrap();
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
    fn cut_by_an_overhanging_boss_carrying_a_pin_owes_a_notch() {
        // Same seating, but the tool carries a pin reaching below the contact plane, so the plane
        // no longer separates the solids and the cut owes a real notch (1 − 0.2·0.2·0.5 = 0.98).
        //
        // This used to be an honest reject: the tool's z=1 cap is an annulus-like face whose
        // *inner* edge rides the pin's walls, and reading its occupancy off the ring's flank put
        // the material on the wrong side, so the class would not label. With the side read from
        // the ring's travel instead (`arrangement::run_body_above`), the notch comes out at the
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
        let canon = plane_classes(&crate::planes::test_judge(&planes));
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
        let (l_tool, _) = build_prism(
            &mut m,
            crate::exact::Swept::along(l_base, Vector3::from_array([0.0, 0.0, 0.4])),
            vec![],
            Vector3::from_array([0.0, 0.0, 1.0]),
            None,
            None,
        )
        .unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, cube, l_tool).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.096).abs() < 1e-12, "volume {vol}");
    }

    // ---- `unify_coplanar_faces`: the general (interface-free) coplanar merge, on hand-built
    // `LocalFace` lists the coincident goldens above never reach — chains, opposite normals,
    // holes, seam edges, and the asymmetric T-junction the global dissolve exists to prevent. ----

    /// An axis-aligned `FaceInfo` at `d` along its normal, with a **non-degenerate `tri`** whose
    /// right-hand normal is `n_out`. The merge reads more than `n_out` now — `loop_winding` and
    /// `point_in_ring` name their arguments by plane and evaluate exact predicates on `tri` — so a
    /// dummy triangle would make those answers meaningless.
    fn mk_axis_plane(m: &mut Model, axis: usize, d: f64, positive: bool) -> PlaneGeom {
        let mut n = [0.0; 3];
        n[axis] = if positive { 1.0 } else { -1.0 };
        let normal = Vector3::from_array(n);
        let mut at = [0.0; 3];
        at[axis] = d;
        let origin = Point3::from_array(at);
        let plane = Plane::from_point_normal(origin, normal).unwrap();
        let surf = m.push_surface(Surface::Plane(plane), nacre_topo::SurfaceDef::Constructed);
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
        let _ = face;
        let tri = [origin, step(i), step(j)];
        PlaneGeom {
            // A hand-built table has no recorded coefficients; the composed-rotation route
            // declines and the fixture takes the same escalating path it always did.
            base_rat: None,
            base: crate::planes::BaseFrame::none(),
            surf,
            plane,
            tri,
            tri_pt3: tri.map(|p| nacre_cip::Pt3::exact(p.as_array()).expect("exact")),
            rotated: false,
            frame_sign: 1, // `plane` is built from `normal`, so the two agree
            exact_coeffs: PlaneGeom::reconcile(&plane, tri, false).0,
            exact_normal: PlaneGeom::reconcile(&plane, tri, false).1,
        }
    }

    /// A `FaceInfo` for the `unify_coplanar_faces` tests, which read none of its geometry; the rest is a
    /// valid-but-unreferenced dummy (`surf`/`face`/`plane` are never dereferenced there).
    fn face(plane_idx: usize, nodes: Vec<Node>, inner: Vec<Vec<Node>>) -> LocalFace {
        LocalFace {
            plane_idx,
            loop_nodes: crate::boolean::Ring::from_clean_names(plane_idx, nodes),
            inner: inner
                .into_iter()
                .map(|r| crate::boolean::Ring::from_clean_names(plane_idx, r))
                .collect(),
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
        let _canon: Vec<usize> = (0..p.len()).collect();
        let v = |x: usize, y: usize| Node::Seam([0, x, y]); // sorted: class, x-plane, y-plane
        let (c00, c10, c20, c30) = (v(1, 5), v(2, 5), v(3, 5), v(4, 5));
        let (c01, c11, c21, c31) = (v(1, 6), v(2, 6), v(3, 6), v(4, 6));
        let faces = vec![
            face(0, vec![c00, c10, c11, c01], vec![]),
            face(0, vec![c10, c20, c21, c11], vec![]),
            face(0, vec![c20, c30, c31, c21], vec![]),
        ];
        let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p)).unwrap();
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
    fn an_overhang_fuse_keeps_the_two_z1_caps_separate() {
        // An overhanging boss splits `z = 1` between two coplanar faces with **opposite** outward
        // normals — the base's exposed top (`+z`) and the boss underside (`-z`). They must not be
        // fused into one face: their `flip` differs, so `unify`'s `(plane_idx, flip)` group key
        // keeps them apart.
        //
        // This replaces two retired tests (`unify_keeps_opposite_normal_coplanar`,
        // `unify_keeps_holed_faces`) that built `Node::Orig` faces to exercise a passthrough the
        // arrangement never triggers — it emits all-`Seam`. The real invariant is exercised here on
        // the production path, in the default `cargo test` run: `overhang_fuse_then_cut_matches_occt`
        // proves it against OCCT but is `#[ignore]`, so this hand-computed volume is the non-ignored
        // guard. A wrong merge collapses the topology — the volume shifts or `validate` speaks.
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
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!(
            (vol - 1.5).abs() < 1e-12,
            "base 1 + boss 0.5, no overlap: {vol}"
        );
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
        let _canon: Vec<usize> = (0..p.len()).collect();
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
        let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p)).unwrap();
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
        let _canon: Vec<usize> = (0..p.len()).collect();
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
        let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p)).unwrap();
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

    /// One geometric plane is one class **whatever the two faces' sizes**.
    ///
    /// ★ **The reason changed, and that is the news.** This used to assert that the coefficient
    /// test *could not* prove these coplanar — two walls of one plane at different face sizes have
    /// un-normalized 4-vectors that are not exactly proportional — and that the coordinate branch
    /// was what earned the merge. Surfaces are interned on their canonical rational coefficients
    /// now, so the two walls are handed **one handle**, and the merge is a handle comparison
    /// before any geometry is asked. The f64 non-proportionality is still real and still pinned,
    /// on the planes themselves, in `nacre_topo`'s
    /// `two_faces_of_one_plane_disagree_in_f64_and_agree_in_the_rationals`.
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
        let PlaneSetup {
            planes: faces_tab,
            geom: _planes,
            plane_ix,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        // The two `+X` walls: same plane x = dx, different face sizes (heights h1 vs h2). Two
        // *faces*, so this searches the face table — the plane table holds one entry for both,
        // which is the property under test.
        let x_walls: Vec<usize> = (0..faces_tab.len())
            .filter(|&i| {
                faces_tab[i].n_out.as_array() == [1.0, 0.0, 0.0]
                    && (faces_tab[i].tri[0].as_array()[0] - dx).abs() < 1e-12
            })
            .collect();
        assert_eq!(x_walls.len(), 2, "one wall from each box: {x_walls:?}");
        let (i, j) = (x_walls[0], x_walls[1]);
        assert_eq!(
            faces_tab[i].surf, faces_tab[j].surf,
            "two faces, one surface — interning collapsed them at construction"
        );
        assert!(
            crate::planes::test_judge(&faces_tab).planes_coplanar(i, j),
            "and the geometry agrees, so nothing rests on the handle alone"
        );
        assert_eq!(plane_ix[i], plane_ix[j], "so they are one plane-table row");
    }

    /// ★★★ **Solving a vertex's three planes lands on its coordinate.**
    ///
    /// This is the premise the whole of `docs/truth-and-cache.md` rests on: a point *is* the
    /// meeting of three surfaces, and the stored coordinate is a rounded answer to that question.
    /// Stage 0 measured it by reading the census; this asserts it on data the kernel itself built,
    /// which is a different claim — the definitions have to be *right*, not merely present.
    ///
    /// ★ **Split by origin, because one row proves nothing.** A `Discovered` vertex's coordinate
    /// was produced by solving exactly this triple, so agreement there is an identity. The rows
    /// that carry weight are `Constructed` and `Moved`, where the coordinate came from somewhere
    /// else entirely — construction arithmetic, or a motion replayed on a base point.
    #[test]
    fn a_vertex_definition_solves_to_its_own_coordinate() {
        // ★ **Three solids, because one would not exercise three kinds of vertex.** The first
        // draft of this test used a fused-then-turned-then-cut solid and reported
        // `Constructed (0,0) Discovered (24,24) Moved (0,0)` with a worst error of exactly zero —
        // a boolean recomputes every vertex it emits, so the only row present was the tautological
        // one. A plain box keeps its constructed corners; a turned box keeps them as `Moved`.
        let mut m = Model::new();
        let plain = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let to_turn = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let OpOutput::Transform { solid: turned } = apply(
            &mut m,
            &Operation::Transform {
                solid: to_turn,
                isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                    axis: nacre_scalar::Axis::Z,
                    point: [nacre_scalar::Rat::from_int(0); 3],
                    angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
                }),
            },
        )
        .expect("turn") else {
            unreachable!("transform yields Transform output")
        };
        m.rebuild_adjacency();
        let a = m.add_cuboid(
            Point3::from_array([10.0, 0.0, 0.0]),
            Point3::from_array([12.0, 3.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([11.0, 1.0, 0.5]),
            Point3::from_array([14.0, 2.0, 2.5]),
        );
        m.rebuild_adjacency();
        let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse")[0];
        m.rebuild_adjacency();

        let mut counts = [(0usize, 0usize); 3]; // (with definition, total) by origin
        let mut worst = [0.0f64; 3];
        let mut diam = 0.0f64;
        for solid in [plain, turned, fused] {
            let shell = m.solids.get(solid).outer;
            let mut seen: Vec<Handle<Vertex>> = Vec::new();
            for &fh in &m.shells.get(shell).faces {
                let face = m.faces.get(fh);
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        for &vh in m.edges.get(he.edge).bounds.iter().flatten() {
                            if !seen.contains(&vh) {
                                seen.push(vh);
                            }
                        }
                    }
                }
            }
            for &vh in &seen {
                let v = *m.vertices.get(vh);
                for c in v.point.as_array() {
                    diam = diam.max(c.abs());
                }
                let kind = match v.origin {
                    Origin::Constructed => 0,
                    Origin::Discovered { .. } => 1,
                    Origin::Moved { .. } => 2,
                };
                counts[kind].1 += 1;
                let Some(VertexDef::ThreePlane(planes)) = v.definition else {
                    continue;
                };
                counts[kind].0 += 1;
                // Free cross-check: where `Origin` carries the same answer, the two must agree.
                if let Origin::Discovered {
                    definition: VertexDef::ThreePlane(od),
                    ..
                } = v.origin
                {
                    assert_eq!(
                        od, planes,
                        "the two records of one vertex's planes disagree"
                    );
                }
                let coeffs = planes.map(|s| match m.surface(s) {
                    nacre_geom::Surface::Plane(p) => p.coefficients(),
                    nacre_geom::Surface::Cylinder(_) => panic!("ThreePlane named a cylinder"),
                });
                let solved = solve_three_planes(coeffs).expect("three planes meeting at a point");
                let d = solved
                    .iter()
                    .zip(v.point.as_array())
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0f64, f64::max);
                worst[kind] = worst[kind].max(d);
            }
        }

        let limit = diam * 2f64.powi(-40);
        eprintln!(
            "[definition coverage] Constructed {:?} worst {:e} | Discovered {:?} worst {:e} | \
             Moved {:?} worst {:e} | limit {:e}",
            counts[0], worst[0], counts[1], worst[1], counts[2], worst[2], limit
        );
        for (i, name) in ["Constructed", "Discovered", "Moved"].iter().enumerate() {
            assert!(
                counts[i].1 > 0,
                "{name} is not exercised — the row proves nothing"
            );
            assert_eq!(
                counts[i].0, counts[i].1,
                "{name} vertices without a definition: {:?}",
                counts[i]
            );
            assert!(
                worst[i] <= limit,
                "{name}: a definition solved {:e} away from its own coordinate (limit {limit:e})",
                worst[i]
            );
        }
        // ★ `Discovered` is the tautological row — its coordinate *came from* solving this very
        // triple, so a zero there says nothing. `Constructed` and `Moved` are the claims.
        assert_eq!(
            worst[1], 0.0,
            "a Discovered vertex is its own solve, exactly"
        );
    }

    /// The negative control for the assertion above: point a definition at the wrong plane and the
    /// solve must land somewhere else. Without this, an agreement test passes on any model whose
    /// planes happen to be near each other.
    #[test]
    fn a_wrong_plane_in_a_definition_is_caught() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        let shell = m.solids.get(s).outer;
        let surfaces: Vec<_> = m
            .shells
            .get(shell)
            .faces
            .iter()
            .map(|&fh| m.faces.get(fh).surface)
            .collect();
        // Take a real corner's triple and swap one plane for another face of the same box.
        let corner = m
            .shells
            .get(shell)
            .faces
            .iter()
            .find_map(|&fh| {
                let face = m.faces.get(fh);
                m.edges.get(face.outer.half_edges[0].edge).bounds
            })
            .expect("a bounded edge")[0];
        let Some(VertexDef::ThreePlane(mut planes)) = m.vertices.get(corner).definition else {
            panic!("a constructed corner has a definition")
        };
        let good = solve_three_planes(planes.map(|h| match m.surface(h) {
            nacre_geom::Surface::Plane(p) => p.coefficients(),
            nacre_geom::Surface::Cylinder(_) => unreachable!(),
        }))
        .expect("meets at a point");
        // Any face not already in the triple. Every one of them moves the point (or leaves the
        // three not meeting at all, which is just as good a refutation).
        let other = *surfaces
            .iter()
            .find(|h| !planes.contains(h))
            .expect("a fourth face");
        planes[0] = other;
        let bad = solve_three_planes(planes.map(|h| match m.surface(h) {
            nacre_geom::Surface::Plane(p) => p.coefficients(),
            nacre_geom::Surface::Cylinder(_) => unreachable!(),
        }));
        let moved = bad.is_none_or(|bad| {
            good.iter()
                .zip(bad)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max)
                > 0.1
        });
        assert!(
            moved,
            "swapping a plane must move the solved point, or the assertion above is vacuous"
        );
    }

    /// Three planes by Cramer, or `None` when they do not meet in a point.
    fn solve_three_planes(p: [[f64; 4]; 3]) -> Option<[f64; 3]> {
        let n = |i: usize| [p[i][0], p[i][1], p[i][2]];
        let (a, b, c) = (n(0), n(1), n(2));
        let det3 = |m: [[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d = det3([a, b, c]);
        if d == 0.0 {
            return None;
        }
        let rhs = [-p[0][3], -p[1][3], -p[2][3]];
        Some(core::array::from_fn(|j| {
            let mut m = [a, b, c];
            for (row, r) in m.iter_mut().zip(rhs) {
                row[j] = r;
            }
            det3(m) / d
        }))
    }

    /// ★ **A profile with a collinear vertex makes its two walls one surface.**
    ///
    /// The two segments either side of a straight-through vertex lie on the *same* plane, so the
    /// vertex they would define is named by two planes, not three — `[S, S, cap]` determines a
    /// line, not a point. Interning is what makes that detectable at all: the two walls share one
    /// handle, so the check is `s[i] == s[j]`, one comparison. Without it they are two handles
    /// naming one plane and the test leaks silently.
    ///
    /// This is the secondary reason surfaces are interned (`docs/truth-and-cache.md`, 「남은 것」 3),
    /// and it is worth pinning because nothing else in the suite builds such a profile.
    #[test]
    fn a_collinear_profile_vertex_gives_its_two_walls_one_surface() {
        let mut m = Model::new();
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        // A unit square whose bottom edge carries a redundant midpoint.
        let profile = Profile2d::polygon(vec![
            p(0.0, 0.0),
            p(0.5, 0.0), // collinear with its neighbours
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
        ]);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile,
                dist: 1.0,
            },
        )
        .expect("extrude") else {
            unreachable!("extrude yields Extrude output")
        };
        m.rebuild_adjacency();
        let shell = m.solids.get(solid).outer;
        let surfaces: Vec<_> = m
            .shells
            .get(shell)
            .faces
            .iter()
            .map(|&fh| m.faces.get(fh).surface)
            .collect();
        // Five profile points ⇒ five wall quads, plus two caps.
        assert_eq!(surfaces.len(), 7, "five wall quads and two caps");
        let mut distinct = surfaces.clone();
        distinct.sort_unstable_by_key(|h| h.index());
        distinct.dedup();
        assert_eq!(
            distinct.len(),
            6,
            "the split bottom edge's two walls are one plane, so one handle: {surfaces:?}"
        );
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
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, overhung, cutter).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let _ = &faces_tab;
        let mut checked = 0usize;
        for sh in solid_shell_handles(&m, overhung) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let tris =
                    combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix).unwrap();
                for t in &tris {
                    // The triple is already dense plane ids: distinct means three real planes.
                    assert!(
                        t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                        "vertex triple {t:?} names one plane twice"
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

    /// **Dense plane ids are order-isomorphic to the sparse roots.** `dense_planes` ranks the class
    /// roots, so any comparison, sort or lex-min over plane indices reads the same either way.
    ///
    /// This is a **migration gate, not a permanent invariant**: it exists so the claim is measured
    /// before the split rides on it, and it retires with `canon` — its subject, not its coverage,
    /// is what goes away.
    #[test]
    fn dense_plane_ids_are_monotone_in_canon() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        // An overhanging boss splits `z = 1` between two faces, so classes really do merge and the
        // ranking really does compress — without that the map is the identity and proves nothing.
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
        // Rebuild the pieces `dense_planes` consumes, so this locks its contract without needing
        // `canon` to escape `plane_index_setup`. `plane_classes` is the same union-find the setup
        // runs; `dense_planes` the same ranking.
        let mut faces = collect_planes(&m, chained).unwrap();
        faces.extend(collect_planes(&m, probe).unwrap());
        let canon = plane_classes(&crate::planes::test_judge(&faces));
        let (geom, plane_ix) = dense_planes(&faces, &canon);
        assert!(
            canon.iter().enumerate().any(|(i, &c)| c != i),
            "fixture has no split plane — the invariant would be vacuous"
        );
        assert!(
            geom.len() < canon.len(),
            "the ranking must actually compress"
        );
        for i in 0..canon.len() {
            for j in 0..canon.len() {
                assert_eq!(
                    canon[i].cmp(&canon[j]),
                    plane_ix[i].cmp(&plane_ix[j]),
                    "faces {i}/{j}: canon {}/{} vs dense {}/{}",
                    canon[i],
                    canon[j],
                    plane_ix[i],
                    plane_ix[j]
                );
            }
        }
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
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, chained, probe).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        // The fixture must actually merge two faces into one plane, or this proves nothing.
        assert!(
            planes.len() < faces_tab.len(),
            "fixture has no split plane — the invariant would be vacuous"
        );
        // A producer hands out dense plane ids (`loop_triples` maps face indices through
        // `plane_ix`), so "every element is a plane, not a face" is now the type, not a runtime
        // check. What remains testable is that the ids are in range and sorted-distinct.
        let mut checked = 0usize;
        for sh in solid_shell_handles(&m, chained) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let mut rings = vec![
                    combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix).unwrap(),
                ];
                rings.extend(combinatorics::hole_rings(&m, fh, p, &inc_a, &jd, &plane_ix).unwrap());
                for t in rings.iter().flatten() {
                    for &k in t {
                        assert!(
                            k < planes.len(),
                            "triple {t:?} names {k}, out of the plane table"
                        );
                    }
                    assert!(
                        t[0] < t[1] && t[1] < t[2],
                        "triple {t:?} is not three distinct planes in sorted order"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "the chained operand has vertices to name");
    }

    /// **A point reads zero on each of its three defining planes.** `Judge::orient3d`'s on-plane
    /// shortcut is a raw `==` against the triple, so a vertex on the query plane must name it by the
    /// same id the query uses. The face/plane split makes that automatic — a plane has exactly one
    /// id now, so the old failure (a vertex named by face 6 of the `z = 1` class invisible to a
    /// query about face 1 of it) cannot be expressed. What is left to check is the identity itself.
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
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, chained, probe).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        assert!(
            planes.len() < faces_tab.len(),
            "fixture has no split plane — the sibling faces this used to distinguish"
        );
        let mut on_plane = 0usize;
        for sh in solid_shell_handles(&m, chained) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let tris =
                    combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix).unwrap();
                for t in &tris {
                    // The vertex lies on exactly its three defining planes; each must read 0.
                    for &q in t {
                        assert_eq!(
                            combinatorics::side_of(&jd, *t, q),
                            0,
                            "vertex {t:?} lies on plane {q} but does not read 0"
                        );
                        on_plane += 1;
                    }
                }
            }
        }
        assert!(on_plane > 0, "some vertex lies on some queried plane");
    }

    // ---- boolean Common algorithm (M5-c3 commit 2) ----

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
            RejectReason::CylinderFace,
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
    /// ★★★★ **Padding the same footprint twice must not leave zero-area faces.**
    ///
    /// It did, and `validate` reported nothing. Padding `1.1` and then `6.6` on the top of a
    /// cuboid left two faces of area `2.2e-16` on the plane `z = 2.1`, whose long edges sat one
    /// ulp apart (`-0.6` against `-0.6000000000000001`).
    ///
    /// The ulp came from the sketch frame's **origin**. `face_frame` takes it from the face's area
    /// centroid, and `ring_area_centroid` was rounding twice more than it needed to, so the first
    /// pad's top face reported its centre as `-1.11e-16` instead of `0`. The second pad then placed
    /// the same profile one ulp away from where the first had placed it, and the kernel — correctly
    /// — built the one-ulp-wide faces that answer describes.
    ///
    /// ★ `SketchPlane::exact()` does not catch this: it checks only that the axes are orthonormal,
    /// never the origin.
    #[test]
    fn padding_one_footprint_twice_leaves_no_zero_area_face() {
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let rect =
            || Profile2d::polygon(vec![p(-1.0, -0.6), p(1.0, -0.6), p(1.0, 0.6), p(-1.0, 0.6)]);
        let mut m = Model::new();
        let mut solid = m.add_cuboid(
            Point3::from_array([-2.0, -2.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        m.rebuild_adjacency();
        // The topmost face pointing up. Everything here is axis-aligned, so this is unambiguous.
        let top = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
            let shell = m.solids.get(s).outer;
            *m.shells
                .get(shell)
                .faces
                .iter()
                .max_by(|&&a, &&b| {
                    let h = |f: Handle<Face>| {
                        crate::ops::face_plane(m, f)
                            .ok()
                            .filter(|sp| sp.normal().as_array()[2] > 0.5)
                            .map(|sp| sp.origin.as_array()[2])
                            .unwrap_or(f64::NEG_INFINITY)
                    };
                    h(a).partial_cmp(&h(b)).unwrap()
                })
                .expect("a face")
        };
        for dist in [1.1, 6.6] {
            let face = top(&m, solid);
            let OpOutput::PadOnFace { solid: out, .. } = apply(
                &mut m,
                &Operation::PadOnFace {
                    face,
                    profile: rect(),
                    dist,
                },
            )
            .expect("pad") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            solid = out;
        }
        let shell = m.solids.get(solid).outer;
        let degenerate: Vec<_> = m
            .shells
            .get(shell)
            .faces
            .iter()
            .filter_map(|&f| {
                let area = nacre_props::face_props(&m, f).ok()?.area;
                (area < 1e-9).then_some((f, area))
            })
            .collect();
        assert!(
            degenerate.is_empty(),
            "zero-area faces survived: {degenerate:?}"
        );
        // A plain box with one rib on top: 6 + 5 walls/cap, no leftovers from the seam.
        assert_eq!(m.shells.get(shell).faces.len(), 11);
    }
    /// ★★★★★ **The target: two ways of reaching one height land on one plane, far from the
    /// origin.**
    ///
    /// The frame's origin is where the drift used to enter — `face_frame` took it from the face's
    /// area centroid, computed in f64 from the face's own vertices, and `exact.rs` lifted that as
    /// truth. Padding the same footprint `1.1` then `6.6` put the second profile an ulp from the
    /// first and left faces of area `2.2e-16`.
    ///
    /// Placed **far from the world origin**, because that is where the projected origin is least
    /// like the old centroid — if anything about the new rule were fragile with distance, a
    /// hundred units of it would show here.
    #[test]
    fn two_routes_to_one_height_share_a_plane_far_from_the_origin() {
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        // On a lid, frame coordinates are world x and y — the origin is the world origin projected
        // onto the plane and the axes are `u = +x̂`, `v = +ŷ`. So this ring is world
        // x ∈ [99, 101], y ∈ [99.4, 100.6].
        let rect = || {
            Profile2d::polygon(vec![
                p(99.0, 100.6),
                p(99.0, 99.4),
                p(101.0, 99.4),
                p(101.0, 100.6),
            ])
        };
        let top = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
            let shell = m.solids.get(s).outer;
            *m.shells
                .get(shell)
                .faces
                .iter()
                .max_by(|&&a, &&b| {
                    let h = |f: Handle<Face>| {
                        crate::ops::face_plane(m, f)
                            .ok()
                            .filter(|sp| sp.normal().as_array()[2] > 0.5)
                            .map(|sp| sp.origin.as_array()[2])
                            .unwrap_or(f64::NEG_INFINITY)
                    };
                    h(a).partial_cmp(&h(b)).unwrap()
                })
                .expect("a face")
        };
        let build = |dists: &[f64]| -> (Model, Handle<Solid>) {
            let mut m = Model::new();
            let mut solid = m.add_cuboid(
                Point3::from_array([98.0, 98.0, 0.0]),
                Point3::from_array([102.0, 102.0, 1.0]),
            );
            m.rebuild_adjacency();
            for &dist in dists {
                let face = top(&m, solid);
                let OpOutput::PadOnFace { solid: out, .. } = apply(
                    &mut m,
                    &Operation::PadOnFace {
                        face,
                        profile: rect(),
                        dist,
                    },
                )
                .expect("pad") else {
                    unreachable!()
                };
                m.rebuild_adjacency();
                solid = out;
            }
            (m, solid)
        };
        let (m1, one) = build(&[7.7]);
        let (m2, two) = build(&[1.1, 6.6]);
        // The two boss tops are the same plane, to the bit.
        let z = |m: &Model, s| {
            crate::ops::face_plane(m, top(m, s))
                .unwrap()
                .origin
                .as_array()[2]
        };
        assert_eq!(
            z(&m1, one),
            z(&m2, two),
            "7.7 and 1.1+6.6 disagree on the cap plane"
        );
        // And the two-step route left nothing degenerate behind.
        let shell = m2.solids.get(two).outer;
        let degenerate: Vec<_> = m2
            .shells
            .get(shell)
            .faces
            .iter()
            .filter_map(|&f| {
                let area = nacre_props::face_props(&m2, f).ok()?.area;
                (area < 1e-9).then_some((f, area))
            })
            .collect();
        assert!(
            degenerate.is_empty(),
            "zero-area faces survived: {degenerate:?}"
        );
        assert_eq!(
            m1.shells.get(m1.solids.get(one).outer).faces.len(),
            m2.shells.get(shell).faces.len(),
            "the two routes did not build the same solid"
        );
    }

    /// ★★★ **One plane, one origin** — even when two faces of it were made by different operations.
    ///
    /// This is what the area centroid could not promise: it was a property of the *face*, so a
    /// boolean that reshaped one face moved its sketch origin away from its coplanar neighbour's.
    /// The projection is a property of the plane, and surfaces are interned, so the two cannot
    /// disagree. Only the origin is asserted — the axes follow the face's `Orientation`, which is
    /// today's behaviour and a separate question.
    #[test]
    fn two_faces_of_one_plane_share_a_sketch_origin() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        // A notch out of one end: the lid `z = 1` becomes two faces of the same plane.
        let cutter = m.add_cuboid(
            Point3::from_array([0.8, -1.0, 0.5]),
            Point3::from_array([1.2, 3.0, 2.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, a, cutter).expect("cut");
        m.rebuild_adjacency();
        let lids: Vec<_> = m
            .shells
            .get(m.solids.get(r).outer)
            .faces
            .iter()
            .filter(|&&f| {
                crate::ops::face_plane(&m, f).is_ok_and(|sp| {
                    sp.normal().as_array()[2] > 0.5 && sp.origin.as_array()[2] == 1.0
                })
            })
            .copied()
            .collect();
        assert_eq!(lids.len(), 2, "the notch should leave two lid faces");
        let o = |f| crate::ops::face_plane(&m, f).unwrap().origin.as_array();
        assert_eq!(
            o(lids[0]),
            o(lids[1]),
            "coplanar faces disagree on the origin"
        );
        assert_eq!(
            o(lids[0]),
            [0.0, 0.0, 1.0],
            "the lid's origin is the world origin projected"
        );
    }
    /// ★★★★★ **The convention, face by face.** This table *is* the rule — every property below is a
    /// consequence of it, and pinning the consequences without pinning the table would let a
    /// different rule that happens to satisfy them slip in.
    #[test]
    fn the_six_axis_directions_get_the_frames_the_convention_names() {
        let v = |a: [f64; 3]| Vector3::from_array(a);
        for (n, u, w) in [
            // The lid: this is the row that must equal `SketchPlane::world_xy`.
            ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ([0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            // Every wall: `v` is +ẑ.
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
            ([-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
            ([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            ([0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ] {
            let (gu, gv) = crate::ops::frame_axes(v(n)).expect("a unit normal has frame axes");
            assert_eq!(gu.as_array(), u, "u for normal {n:?}");
            assert_eq!(gv.as_array(), w, "v for normal {n:?}");
            // ★ And the axes stay exactly representable, so the rational construction path still
            // fires — losing that would drop every axis-aligned model to f64 silently.
            let plane = SketchPlane::from_axes(Point3::origin(), gu, gv);
            assert!(plane.exact().is_some(), "exact path lost for normal {n:?}");
        }
    }

    /// The lid's frame and [`SketchPlane::world_xy`] name the same plane, so they must name it the
    /// same way. They did not: `any_perpendicular` gave the lid `u = −ŷ, v = +x̂`, ninety degrees
    /// round, and the kernel carried both spellings at once.
    #[test]
    fn a_lid_gets_the_same_frame_as_the_world_xy_plane() {
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        let (u, v) = crate::ops::frame_axes(up).unwrap();
        let w = SketchPlane::world_xy();
        assert_eq!(u.as_array(), w.x_axis.as_array());
        assert_eq!(v.as_array(), w.y_axis.as_array());
    }

    /// **On anything but a horizontal face, `v` points up.** `u = ẑ × n` is horizontal, so
    /// `v·ẑ = 1 − n_z² > 0` whenever `n` is not vertical — the reason a sketch on a wall has "up"
    /// where a person expects it. Checked on tilts the axis-aligned table cannot reach.
    #[test]
    fn every_non_horizontal_face_has_its_v_pointing_up() {
        for n in [
            [1.0, 1.0, 0.0],
            [0.6, 0.0, 0.8],
            [-0.3, 0.5, -0.81],
            [0.0, 1.0, 0.001],
            [7.0, -13.0, 5.0],
        ] {
            let n = Vector3::from_array(n).normalize().unwrap();
            let (u, v) = crate::ops::frame_axes(n).unwrap();
            assert_eq!(u.as_array()[2], 0.0, "u must be horizontal for {n:?}");
            assert!(v.as_array()[2] > 0.0, "v points down for {n:?}");
        }
    }

    /// `u ⊥ n`, `|u| = 1`, and `(u, v, n)` right-handed — for the tilted normals too, where the
    /// table above says nothing.
    #[test]
    fn the_reference_axis_is_a_unit_normal_perpendicular() {
        for n in [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 1.0, 1.0],
            [-2.0, 0.5, 3.25],
            [1e-9, 0.0, 1.0],
        ] {
            let n = Vector3::from_array(n).normalize().unwrap();
            let (u, v) = crate::ops::frame_axes(n).unwrap();
            assert!((u.norm() - 1.0).abs() < 1e-15, "|u| for {n:?}");
            assert!(u.dot(n).abs() < 1e-15, "u·n for {n:?}");
            assert!((u.cross(v) - n).norm() < 1e-15, "handedness for {n:?}");
        }
    }

    /// ★★ **The jump at the poles is intended, not a bug to be fixed later.**
    ///
    /// No continuous tangent frame exists on the sphere, so some set of normals must jump; this
    /// convention spends that budget on the two poles and nowhere else. A normal a billionth off
    /// vertical takes the other branch and lands ninety degrees away — pinned here so the next
    /// reader can see it was chosen.
    #[test]
    fn the_frame_jumps_at_the_poles_and_that_is_the_deal() {
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        let tilted = Vector3::from_array([1e-9, 0.0, 1.0]).normalize().unwrap();
        assert_eq!(
            crate::ops::frame_axes(up).unwrap().0.as_array(),
            [1.0, 0.0, 0.0]
        );
        assert_eq!(
            crate::ops::frame_axes(tilted).unwrap().0.as_array(),
            [0.0, 1.0, 0.0]
        );
    }

    /// No `-0.0` reaches a caller. It compares equal to `0.0` and lifts to the same rational, so
    /// this is presentation only — but a frame printed as `[-0.0, 1.0, 0.0]` reads like a defect.
    #[test]
    fn no_axis_component_is_negative_zero() {
        for n in [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ] {
            let n = Vector3::from_array(n);
            let (u, v) = crate::ops::frame_axes(n).unwrap();
            for c in u.as_array().iter().chain(v.as_array().iter()) {
                assert!(
                    !(*c == 0.0 && c.is_sign_negative()),
                    "negative zero in {n:?}'s frame"
                );
            }
        }
    }

    /// The zero vector has no frame — the only `None`.
    #[test]
    fn the_zero_vector_has_no_frame_axes() {
        assert!(crate::ops::frame_axes(Vector3::zero()).is_none());
    }
    /// ★★★★★ **A boss on a tilted face, and another beside it — `pad` used to lose the second one.**
    ///
    /// The prism's far cap and the first boss's cap are one plane, and the boolean says so: its
    /// classes are decided with evidence, and on a tilted face that evidence is a composed-rotation
    /// proof or a coincidence within the limit, never a handle match — the cap's surface has no
    /// rational coefficients to intern by, so it is minted fresh.
    ///
    /// `find_face_coplanar_with` then had to guess which surface the class had collapsed to, from
    /// handles and an exact `plane_side` on f64 points. Both miss: the survivor carries the *other*
    /// operand's surface, and the two f64 planes sit `1.8e-15` apart. `pad` turned "cap not found"
    /// into a hard error and **threw away a correct solid** — the volume was already right.
    ///
    /// Now it asks the boolean instead.
    #[test]
    fn a_second_boss_on_a_tilted_face_keeps_its_cap() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let mut m = Model::new();
        let mut s = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        m.rebuild_adjacency();
        // Two turns, so the face's normal is off every world axis and its frame is not exact.
        for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid: s,
                    isometry: Isometry::rotation(Rotation {
                        axis,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                    }),
                },
            )
            .expect("turn") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            s = solid;
        }
        // Where the original +z went, so the face can be picked without a heuristic tie.
        let (sy, cy) = (53f64).to_radians().sin_cos();
        let (sz, cz) = (17f64).to_radians().sin_cos();
        let up = Vector3::from_array([cz * sy, sz * sy, cy]);
        // The **original** tilted face, not a boss raised on it: lowest along `up` among the faces
        // that point that way. Taking the highest would stack the second boss on the first.
        let facing = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
            *m.shells
                .get(m.solids.get(s).outer)
                .faces
                .iter()
                .filter(|&&f| {
                    crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99)
                })
                .min_by(|&&a, &&b| {
                    let h = |f: Handle<Face>| {
                        (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                    };
                    h(a).partial_cmp(&h(b)).unwrap()
                })
                .expect("a face along up")
        };
        // The frame is a function of the plane, so it survives the face being reshaped by the first
        // pad — both columns are placed from one reading.
        let (cu, cv) = {
            let f = facing(&m, s);
            let sp = crate::ops::face_plane(&m, f).expect("planar");
            let d = nacre_props::face_props(&m, f).unwrap().centroid - sp.origin;
            (d.dot(sp.x_axis), d.dot(sp.y_axis))
        };
        let before = nacre_props::mass_props(&m, s).unwrap().volume;
        let mut caps = Vec::new();
        for (lo, hi) in [(-1.0f64, -0.4f64), (0.4, 1.0)] {
            let f = facing(&m, s);
            let profile = Profile2d::polygon(vec![
                p(cu + lo, cv - 0.5),
                p(cu + hi, cv - 0.5),
                p(cu + hi, cv + 0.5),
                p(cu + lo, cv + 0.5),
            ]);
            let OpOutput::PadOnFace { solid, top_face } = apply(
                &mut m,
                &Operation::PadOnFace {
                    face: f,
                    profile,
                    dist: 7.7,
                },
            )
            .expect("a boss on a tilted face keeps its cap") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            s = solid;
            caps.push(top_face);
        }
        // Each column is 0.6 × 1.0 × 7.7 = 4.62 of material.
        let after = nacre_props::mass_props(&m, s).unwrap().volume;
        assert!(
            (after - before - 2.0 * 4.62).abs() < 1e-9,
            "volume {before} -> {after}"
        );
        // ★ And the handle it returned is the boss top, not some other face that happened to pass:
        // area 0.6, outward along the face's normal.
        for cap in caps {
            let props = nacre_props::face_props(&m, cap).expect("a planar cap");
            assert!((props.area - 0.6).abs() < 1e-9, "cap area {}", props.area);
            assert!(
                props.normal.is_some_and(|n| n.dot(up) > 0.99),
                "cap faces the wrong way"
            );
        }
        assert!(nacre_validate::validate(&m).is_empty());
    }

    // ---- tilted-face sketches: what the kernel does today (2b's corpus) ----
    //
    // ★★★ These pin **today's** behaviour, not a wish. The suite had almost none of them, so the
    // stage that makes tilted frames exact would otherwise be built with nothing to measure against.
    // Each says which part 2b is expected to change and which part must not move.

    /// A profile edge running `(1, 2)` sweeps a wall whose normal is `(2, 1, 0)` — and **no
    /// rational-degree rotation reaches it** (the angle is `atan(1/2)`). That is the case
    /// `docs/truth-and-cache.md` names for `Motion::Frame`.
    ///
    /// What holds today and must keep holding:
    ///
    /// * the wall's **world coefficients are rational** — `[2, 1, 0, −10]`, so a frame built on it
    ///   has exact data to derive from;
    /// * its sketch frame follows the arbitrary-axis convention — origin at the world origin's
    ///   projection `(4, 2, 0) = (10/5)·(2,1,0)`, `u = ẑ × n` normalized, `v = +ẑ`;
    /// * a pad on it **works**, through the f64 path.
    ///
    /// What 2b changes: `exact()` is `false` here, so the prism it raises records no rational
    /// coefficients. ★ The axes must **not** move — this plane's own frame is the world, so
    /// `ẑ × n` is the same vector before and after.
    #[test]
    fn a_sketch_on_a_prism_side_wall_takes_the_f64_path_today() {
        let (m, wall) = prism_with_a_slanted_wall();
        let sp = crate::ops::face_plane(&m, wall).expect("planar");
        let c = m
            .surface_name
            .get(&m.faces.get(wall).surface)
            .expect("a world-frame wall has a name")
            .narrow()
            .expect("a world-frame wall's name is narrow");
        assert_eq!(c.map(|r| r.to_f64()), [2.0, 1.0, 0.0, -10.0]);
        assert_eq!(
            sp.origin.as_array(),
            [4.0, 2.0, 0.0],
            "the projected origin"
        );
        assert_eq!(
            sp.y_axis.as_array(),
            [0.0, 0.0, 1.0],
            "v points up on a wall"
        );
        assert!(
            (sp.x_axis - Vector3::from_array([-1.0, 2.0, 0.0]).normalize().unwrap()).norm() < 1e-15,
            "u is ẑ × n normalized, got {:?}",
            sp.x_axis.as_array()
        );
        // ★ The gap 2b closes: the frame is not exact, so nothing built here records coefficients.
        assert!(
            sp.exact().is_none(),
            "a tilted frame has no rational form today"
        );
    }

    /// The same wall, actually used: a boss on it comes out right through the f64 path. Pinned so
    /// that making the frame exact cannot change the **answer**, only how it is recorded.
    #[test]
    fn a_boss_on_a_slanted_wall_is_correct_today() {
        let (mut m, wall) = prism_with_a_slanted_wall();
        let profile = centred_on(&m, wall, 0.5);
        let OpOutput::PadOnFace { solid, top_face } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: wall,
                profile,
                dist: 1.0,
            },
        )
        .expect("pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        // Base prism 15 × 3 = 45, plus a 1 × 1 × 1 boss.
        let props = nacre_props::mass_props(&m, solid).unwrap();
        assert!(
            (props.volume - 46.0).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!(
            (nacre_props::face_props(&m, top_face).unwrap().area - 1.0).abs() < 1e-9,
            "the boss top is 1 × 1"
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// ★★★ **A sketch on a wall raised from a sketch on a wall** — the nesting `Motion::Frame` was
    /// redesigned for. The second wall's *world* normal is irrational, so its frame cannot be
    /// written down by naming a normal; only by naming the plane.
    ///
    /// It works today, through f64. Pinned because nesting is where a frame that names its plane
    /// by handle must terminate its recursion.
    #[test]
    fn a_sketch_on_a_wall_raised_from_a_slanted_wall_works_today() {
        let (mut m, wall) = prism_with_a_slanted_wall();
        let profile = centred_on(&m, wall, 0.5);
        let OpOutput::PadOnFace { solid, top_face } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: wall,
                profile,
                dist: 1.0,
            },
        )
        .expect("pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        // A side face of that boss: not the cap, not on the original wall's plane.
        let cap_n = nacre_props::face_props(&m, top_face)
            .unwrap()
            .normal
            .unwrap();
        let side = *m
            .shells
            .get(m.solids.get(solid).outer)
            .faces
            .iter()
            .find(|&&f| {
                f != top_face
                    && nacre_props::face_props(&m, f).is_ok_and(|p| {
                        (p.area - 1.0).abs() < 1e-9
                            && p.normal.is_some_and(|n| n.dot(cap_n).abs() < 0.5)
                    })
            })
            .expect("a boss side face");
        let profile = centred_on(&m, side, 0.2);
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: side,
                profile,
                dist: 0.5,
            },
        )
        .expect("nested pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let props = nacre_props::mass_props(&m, solid).unwrap();
        // The nested boss is 0.4 × 0.4 × 0.5 = 0.08.
        assert!(
            (props.volume - 46.08).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A pocket on the slanted wall — the other face-based operation, so the sweep runs inward.
    #[test]
    fn a_pocket_in_a_slanted_wall_is_correct_today() {
        let (mut m, wall) = prism_with_a_slanted_wall();
        let profile = centred_on(&m, wall, 0.5);
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: wall,
                profile,
                dist: 0.5,
            },
        )
        .expect("pocket") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let props = nacre_props::mass_props(&m, solid).unwrap();
        assert!(
            (props.volume - 44.5).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A pentagonal prism whose fourth wall is slanted, and that wall's handle. Footprint area 15
    /// (a 4×4 square less the 1×2 triangle the slant cuts off), swept 3.
    fn prism_with_a_slanted_wall() -> (Model, Handle<Face>) {
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let mut m = Model::new();
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d::polygon(vec![
                    p(0.0, 0.0),
                    p(4.0, 0.0),
                    p(4.0, 2.0),
                    p(3.0, 4.0),
                    p(0.0, 4.0),
                ]),
                dist: 3.0,
            },
        )
        .expect("extrude") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let wall = *m
            .shells
            .get(m.solids.get(solid).outer)
            .faces
            .iter()
            .find(|&&f| {
                crate::ops::face_plane(&m, f).is_ok_and(|sp| {
                    let n = sp.normal().as_array();
                    n[0].abs() > 0.1 && n[1].abs() > 0.1 && n[2].abs() < 1e-12
                })
            })
            .expect("a slanted wall");
        (m, wall)
    }

    /// A `2·half` square centred on `face`, in that face's own sketch frame. The frame is a
    /// function of the plane, so reading it once is enough even if the face is reshaped later.
    /// ★★★★★ **The payoff, and its limit — both measured.**
    ///
    /// Two bosses of the same height on one tilted face used to be two plane records that agreed
    /// only if their f64 coefficients happened to. Written in the plane's own frame they are both
    /// `w = 7.7`, and `SurfaceKey` is `(coefficients, motion)` — so they are **one
    /// `Handle<Surface>` at construction**, before anything is compared. And every face of the
    /// result states itself exactly, where before a tilted sketch recorded nothing at all.
    ///
    /// ★★★ **What this does *not* buy, stated plainly**: `7.7` against `1.1 + 6.6` — the target
    /// the plan named. Stacking sketches the second boss on the **first boss's cap**, which is a
    /// different plane and therefore a different frame, so the two caps come out `w = 7.7` and
    /// `w = 6.6` — two exact descriptions of one plane that `SurfaceKey` cannot equate. That is
    /// not a regression (before this they had no descriptions at all, and the merge still happens
    /// through the judge), but the plan's headline claim only holds for sketches sharing a frame,
    /// which is what this pins instead.
    #[test]
    fn two_bosses_on_one_tilted_face_share_a_cap_plane_by_name() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let mut m = Model::new();
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d::polygon(vec![
                    p(0.0, 0.0),
                    p(4.0, 0.0),
                    p(4.0, 4.0),
                    p(0.0, 4.0),
                ]),
                dist: 3.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let mut s = solid;
        for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid: s,
                    isometry: Isometry::rotation(Rotation {
                        axis,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                    }),
                },
            )
            .unwrap() else {
                unreachable!()
            };
            m.rebuild_adjacency();
            s = solid;
        }
        let (sy, cy) = (53f64).to_radians().sin_cos();
        let (sz, cz) = (17f64).to_radians().sin_cos();
        let up = Vector3::from_array([cz * sy, sz * sy, cy]);
        let facing = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
            *m.shells
                .get(m.solids.get(s).outer)
                .faces
                .iter()
                .filter(|&&f| {
                    crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99)
                })
                .min_by(|&&a, &&b| {
                    let h = |f: Handle<Face>| {
                        (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                    };
                    h(a).partial_cmp(&h(b)).unwrap()
                })
                .unwrap()
        };
        let mut caps = Vec::new();
        for (lo, hi) in [(-1.0f64, -0.4f64), (0.4, 1.0)] {
            let f = facing(&m, s);
            let sp = crate::ops::face_plane(&m, f).unwrap();
            let d = nacre_props::face_props(&m, f).unwrap().centroid - sp.origin;
            let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
            let profile = Profile2d::polygon(vec![
                p(cu + lo, cv - 0.5),
                p(cu + hi, cv - 0.5),
                p(cu + hi, cv + 0.5),
                p(cu + lo, cv + 0.5),
            ]);
            let OpOutput::PadOnFace { solid, top_face } = apply(
                &mut m,
                &Operation::PadOnFace {
                    face: f,
                    profile,
                    dist: 7.7,
                },
            )
            .expect("boss") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            s = solid;
            let su = m.faces.get(top_face).surface;
            assert_eq!(
                m.surface_name
                    .get(&su)
                    .and_then(|n| n.narrow())
                    .map(|c| c.map(|r| r.to_f64())),
                Some([0.0, 0.0, 10.0, -77.0]),
                "★ a cap raised in a frame records `w = 7.7` there, exactly"
            );
            caps.push(su);
        }
        assert_eq!(caps[0], caps[1], "one plane, one handle, by name");
        // Rule audit: how many of the result's face surfaces state themselves exactly.
        let sh = m.solids.get(s).outer;
        let (mut with, mut without) = (0, 0);
        for &f in &m.shells.get(sh).faces {
            let su = m.faces.get(f).surface;
            if m.surface_name.contains_key(&su) {
                with += 1
            } else {
                without += 1
            }
        }
        assert_eq!(
            (with, without),
            (16, 0),
            "★ every face of a twice-turned, twice-bossed result states itself exactly"
        );
    }

    /// ★★★★★ **A plane states itself exactly, and its f64 axes are the realization of that.**
    ///
    /// This is the whole point of closing the struct: `(1, 1, 1)` is coefficients `[1, 1, 1, 0]`,
    /// three integers, while the unit axes derived from it square to `0.9999999999999999…`. The
    /// old API stored only the axes and threw the normal away, so nothing exact survived the door.
    #[test]
    fn a_named_plane_records_what_its_caller_stated() {
        let f = |d: &PlaneDef| d.coeffs.map(|r| r.to_f64());
        // The three world planes, with the axes the script layer documents.
        for (p, want, u) in [
            (
                SketchPlane::world_xy(),
                [0.0, 0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ),
            (
                SketchPlane::world_yz(),
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ),
            // ★ `ẑ × n` would give `−x̂` here; a named plane says `+u = ẑ` instead.
            (
                SketchPlane::world_zx(),
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ),
        ] {
            let d = p.def.expect("a world plane states itself");
            assert_eq!(f(&d), want);
            assert_eq!(d.ref_dir.map(|r| r.to_f64()), u);
            assert_eq!(p.x_axis().as_array(), u, "the axis follows the definition");
        }
        // A tilted normal the caller wrote: exact coefficients, though its axes never can be.
        let tilt =
            SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
                .unwrap();
        assert_eq!(f(&tilt.def.unwrap()), [1.0, 1.0, 1.0, 0.0]);
        assert!(
            tilt.exact().is_none(),
            "★ the axes still have no exact form — that is what the definition exists to replace"
        );
        // Three written points: the plane is exact and `+u` runs toward `x_point`.
        let tp = SketchPlane::through_points(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([1.0, 2.0, 0.0]),
            Point3::from_array([1.0, 0.0, 3.0]),
        )
        .unwrap();
        let d = tp.def.unwrap();
        assert_eq!(f(&d), [1.0, 0.0, 0.0, -1.0], "the plane x = 1");
        assert_eq!(
            d.ref_dir.map(|r| r.to_f64()),
            [0.0, 2.0, 0.0],
            "x_point − origin"
        );
        assert_eq!(d.origin.map(|r| r.to_f64()), [1.0, 0.0, 0.0]);
        // Moving the sketch origin keeps the plane and moves only `(0, 0)`.
        let moved = tp.with_origin(Point3::from_array([1.0, 5.0, 5.0]));
        let m = moved.def.unwrap();
        assert_eq!(f(&m), f(&d), "same plane");
        assert_eq!(m.origin.map(|r| r.to_f64()), [1.0, 5.0, 5.0]);
        // ★★★★★ **The invariant that ties the two halves together: the origin is *on* the plane.**
        // Without it a definition can describe one plane while its sketch sits on another —
        // exactly what `with_origin` did until a boolean lost 0.04 of volume over it.
        for p in [
            SketchPlane::world_xy(),
            SketchPlane::world_yz(),
            SketchPlane::world_zx(),
            SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
            SketchPlane::world_zx().with_origin(Point3::from_array([10.0, 20.0, 5.0])),
            SketchPlane::from_origin_normal(
                Point3::from_array([0.0, 1.3, 0.0]),
                Vector3::from_array([0.0, -1.0, 0.0]),
            )
            .unwrap(),
            tilt,
            tp,
            moved,
        ] {
            let d = p.def.expect("stated");
            let mut s = d.coeffs[3];
            for k in 0..3 {
                s = s
                    .checked_add(d.coeffs[k].checked_mul(d.origin[k]).unwrap())
                    .unwrap();
            }
            assert_eq!(
                s,
                nacre_scalar::Rat::from_int(0),
                "the origin must lie on the plane it names: {:?} vs {:?}",
                d.coeffs.map(|r| r.to_f64()),
                d.origin.map(|r| r.to_f64())
            );
        }
        // ★ And the axes-only route is honest about having no definition.
        assert!(
            SketchPlane::from_axes(
                Point3::origin(),
                Vector3::from_array([1.0, 0.0, 0.0]),
                Vector3::from_array([0.0, 1.0, 0.0]),
            )
            .def
            .is_none()
        );
    }

    /// ★★★★★ **A prism raised on a named tilted plane states every one of its faces.**
    ///
    /// In world coordinates that plane's axes are irrational, so the whole prism used to drop to
    /// f64 and record nothing. Two things fixed it: the base cap **is** the plane the caller
    /// named, so it states itself in the world; and the walls and far cap are built **inside that
    /// plane's frame**, where the axes are `x̂`/`ŷ` and the profile's own decimals are the truth.
    ///
    /// ★★★ **The base cap stays in the world on purpose.** Writing it as `[0,0,1,0]` in this
    /// prism's frame would be a second exact description of one plane under a different
    /// `SurfaceKey` — the duplication this work exists to remove. Stated in the world it is
    /// `Constructed`, its judgment stays exact, and two extrudes share it whatever frames they chose.
    #[test]
    fn a_prism_on_a_named_tilted_plane_states_all_of_its_faces() {
        let plane =
            SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
                .unwrap();
        assert!(plane.exact().is_none(), "the axes have no exact form");
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane,
                profile: square(),
                dist: 1.0,
            },
        )
        .expect("extrude on a tilted plane") else {
            unreachable!()
        };
        let coeffs = |f: Handle<Face>, m: &Model| {
            m.surface_name
                .get(&m.faces.get(f).surface)
                .and_then(|n| n.narrow())
                .map(|c| c.map(|r| r.to_f64()))
        };
        assert_eq!(
            coeffs(faces[0], &m),
            Some([1.0, 1.0, 1.0, 0.0]),
            "★ the base cap is the caller's plane, in the world"
        );
        assert!(
            matches!(
                m.surface_defs.get(&m.faces.get(faces[0]).surface),
                Some(nacre_topo::SurfaceDef::Constructed)
            ),
            "★ and it carries no motion, so its judgment stays exact"
        );
        assert_eq!(
            coeffs(faces[1], &m),
            Some([0.0, 0.0, 1.0, -1.0]),
            "★ the far cap is `w = dist` in the frame"
        );
        // ★ Every face now states itself — that is the whole measurement.
        for &f in &faces {
            assert!(coeffs(f, &m).is_some(), "a face with no exact plane");
        }

        // ★★★★★ **Two extrudes on one named plane put their far caps on one handle** — by name,
        // at construction, with no f64 comparison. That is what the frame buys over the f64 path,
        // where the two would agree only if their rounded coefficients happened to.
        let cap_of = |m: &mut Model, d: f64| -> Handle<nacre_geom::Surface> {
            let OpOutput::Extrude { faces, .. } = apply(
                m,
                &Operation::Extrude {
                    plane,
                    profile: square(),
                    dist: d,
                },
            )
            .expect("extrude") else {
                unreachable!()
            };
            m.faces.get(faces[1]).surface
        };
        let a = cap_of(&mut m, 2.5);
        let b = cap_of(&mut m, 2.5);
        assert_eq!(a, b, "one height on one plane is one plane");
        assert_ne!(a, cap_of(&mut m, 2.6), "and a different height is not");
    }

    fn centred_on(m: &Model, face: Handle<Face>, half: f64) -> Profile2d {
        let sp = crate::ops::face_plane(m, face).expect("planar");
        let d = nacre_props::face_props(m, face).unwrap().centroid - sp.origin;
        let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        Profile2d::polygon(vec![
            p(cu - half, cv - half),
            p(cu + half, cv - half),
            p(cu + half, cv + half),
            p(cu - half, cv + half),
        ])
    }
}
