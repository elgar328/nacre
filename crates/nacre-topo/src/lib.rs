//! b-rep topology and the truth-only `Model` aggregate (design.md §2, §4).
//!
//! The topology (`Vertex`/`Edge`/`Face`/`Loop`/`HalfEdge`/`Shell`/`Solid`)
//! references exact geometry only by `Handle` — geometry never knows about
//! topology, topology never inspects coordinates (overview 절대원칙 2).
//!
//! [`Model`] is the **truth**: exact geometry stores + topology stores + the
//! derived [`Adjacency`] cache. It holds no tessellation and no operation log —
//! a mesh cache and the op log are companions owned at higher layers (§0: "the
//! model is the replay result"; putting them here would make topo depend on
//! tess/ops and break the truth/cache split).

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
mod adjacency;
mod topology;

pub use adjacency::{Adjacency, nonmanifold_vertices};
pub use topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

use nacre_geom::{Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Angle, Axis, Rat};
use nacre_store::{Handle, Store};
use std::collections::{HashMap, HashSet};

/// How a discovered vertex is *defined* — the primitives whose intersection it
/// is (design §4). This definition is the **truth**; the vertex's `f64` point is
/// a cache derived from it. Sign decisions (in/out, orientation) feed this
/// definition to the indirect predicates rather than the cached coordinate
/// (design §3, §8 M5), so a discovered vertex must carry it.
///
/// M5 polyhedral vertices are three-plane intersections; later milestones add
/// variants (a line∩plane point, quadric intersections). `Handle<Surface>` is
/// `Copy` regardless of `Surface`, so this stays `Copy`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VertexDef {
    /// The intersection of three planes (their surface handles).
    ThreePlane([Handle<Surface>; 3]),
    /// A point on the intersection **curve** of two surfaces — the M3 cylinder seam vertex:
    /// the rim circle (lateral cylinder ∩ cap plane) at parameter `θ = 0` (S7).
    ///
    /// ★ The pair pins a curve; what picks the point on it is the cylinder's `ref_dir`, which
    /// sits in the cylinder's truth since M6-0 ([`CylinderDef`]): `OnSeam([cylinder, cap])`
    /// **is** "the rim ∩ the `+ref_dir` ray" — a unique point, exactly designated. The
    /// definition is complete; what remains deferred (with 3b, recorded honestly) is the
    /// machinery that *regenerates* the cached coordinate from it.
    /// M6 grows the vocabulary by variants, each stating its own truth — the invariants are
    /// per-variant (Q5's doctrine); [`VertexDef::Branch`] (M6-1) is the first.
    OnSeam([Handle<Surface>; 2]),
    /// One of the (at most two) points where two planes' meet line crosses a cylinder's
    /// lateral surface (M6-1) — the structure says the carrier kinds, deliberately not an
    /// array of three lookalike handles (validate's carrier check becomes structural).
    ///
    /// ★ **`root` picks the point; its meaning is a convention.** The meet line's direction
    /// is `ℓ = n₁ × n₂` where `n₁`/`n₂` are the **canonical-name normals** of `planes[0]`/
    /// `planes[1]` **in stored (ascending-handle) order** — canonicalization fixes each
    /// normal's sign ("first nonzero component positive"), so ℓ is deterministic; `Lo`/`Hi`
    /// is ascending parameter along ℓ (`nacre_scalar::quad`'s pair order). A tangency
    /// (double root) is one point and spells it `Lo`. ★★ Anything that re-sorts the two
    /// plane handles (a transform remapping them) must **toggle `root` when they swap** —
    /// swapping flips ℓ and with it the meaning of `Lo`/`Hi`.
    ///
    /// The producer arrives with M6-2's boolean; until then hand-built fixtures and validate
    /// are the consumers (the `FaceMisoriented`-control precedent).
    Branch {
        /// The two cutting planes, ascending handle order (the `ThreePlane` precedent).
        planes: [Handle<Surface>; 2],
        /// The cylinder whose lateral surface the meet line crosses.
        cylinder: Handle<Surface>,
        /// Which of the two crossings, along the canonical line direction.
        root: QuadRoot,
    },
}

/// Which root of the two-point plane·plane·cylinder crossing a [`VertexDef::Branch`] means —
/// ascending parameter along the canonical meet-line direction (see `Branch`'s doc for the
/// full convention). Definition vocabulary, so it lives here beside [`VertexDef`], not in
/// scalar (whose pair is positional).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum QuadRoot {
    /// The smaller parameter — and the spelling a tangency's single point takes.
    Lo,
    /// The larger parameter.
    Hi,
}

impl QuadRoot {
    /// The other root — what a re-sort that swaps the two plane handles must apply.
    #[inline]
    pub fn flipped(self) -> QuadRoot {
        match self {
            QuadRoot::Lo => QuadRoot::Hi,
            QuadRoot::Hi => QuadRoot::Lo,
        }
    }
}

impl VertexDef {
    /// Every carrier handle the definition references, in stored order — the one spelling of
    /// "the surfaces this vertex is defined by" (reference integrity, the off-definition
    /// check and the remappability gate all ask exactly this; carrier *kinds* stay
    /// per-variant checks).
    pub fn carriers(&self) -> impl Iterator<Item = Handle<Surface>> {
        let (arr, n): ([Handle<Surface>; 3], usize) = match *self {
            VertexDef::ThreePlane(s) => (s, 3),
            VertexDef::OnSeam([a, b]) => ([a, b, b], 2),
            VertexDef::Branch {
                planes: [a, b],
                cylinder,
                ..
            } => ([a, b, cylinder], 3),
        };
        arr.into_iter().take(n)
    }
}

/// One motion in a history — what a [`MotionNode`] carries.
///
/// Motions do **not** commute, so a history is one ordered chain, never a store per kind: "turn
/// then place" and "place then turn" are different motions and must stay tellable apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Motion {
    /// One axis-aligned rotation about the line through the rational pivot `point`.
    Rotate {
        axis: Axis,
        point: [Rat; 3],
        angle: Angle,
    },
    /// One exact rational translation.
    Translate { offset: [Rat; 3] },
    /// One reflection in the coordinate plane `axis = offset` (`x ↦ 2·offset − x` on that axis).
    ///
    /// **Improper** — the only motion here with `det = −1`. It preserves lengths, distances and
    /// incidence like the others, but it *negates* every determinant of its images rather than
    /// leaving them alone. Judgments that answer a determinant question in a shared pre-motion
    /// frame must account for that; see `nacre-ops`' `BaseFrame` and `nacre-cip`'s `shared_base`.
    Mirror { axis: Axis, offset: Rat },
    /// A change of basis **into a plane's own frame** — what makes a sketch on a tilted face
    /// exact. Coordinates written against this node are read as `(u, v, w)` in that plane's
    /// frame; the motion carries them out into the frame the plane itself lives in.
    ///
    /// ★★★★ **The normal is named, not spelled.** A wall raised on a tilted face has an
    /// irrational world normal — there is no `[Rat; 3]` for it — so the node points at the
    /// *plane*, whose data is rational **in its own frame**. The recursion terminates
    /// because it walks down to a plane with no frame, which is where the world is.
    ///
    /// ★★★ **Where the frame sits on the plane is [`FramePlacement`]** (S4): the derived
    /// `Canonical` convention by default, or the caller's `Named` values. See its own doc.
    ///
    /// ★★★ **`flip` is what makes the node a whole frame and not half of one.** Canonical plane
    /// coefficients carry **no direction** — the first nonzero component is forced positive,
    /// because their question is *"are these the same plane"*. A frame's `ŵ` is a direction, and
    /// two faces of one plane can face opposite ways; `flip` says to negate the coefficients, so
    /// the node names the sense as well as the plane. Without it a sketch on a reversed face comes
    /// out mirrored in `u` with its sweep running inward (measured: a tilted second boss came back
    /// `PadMissesFace`). It is measured (realized `ŵ` against the face's outward normal), never
    /// chosen by a caller.
    ///
    /// ★ **Proper** (`det = +1`) — `(u, v, w)` is right-handed by construction (`v = w × u`), so
    /// unlike [`Mirror`](Self::Mirror) it contributes nothing to a chain's parity. Negating the
    /// coefficients flips `ŵ` and `û` together and leaves `v̂`, which is a half-turn: still proper.
    Frame {
        plane: Handle<Surface>,
        placement: FramePlacement,
        flip: bool,
    },
}

/// Where a sketch frame sits on its plane — the derived convention, or the caller's values.
/// The dichotomy mirrors `PlanePoints` (`docs/truth-and-cache.md`): state a value when it can
/// be written, derive when it cannot.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FramePlacement {
    /// **The default.** The frame is a pure function of the plane — origin at the world
    /// origin's projection (`(−d/n·n)·n`), axes by the arbitrary-axis convention (`u = ẑ × n`,
    /// `ŷ × n` when the normal is exactly vertical) — derived when the chain is flattened and
    /// stored nowhere. That is what lets a plane whose canonical values overflow `i128` (or
    /// whose name is `Wide`, S2) take the **same convention through arbitrary-precision
    /// realization** instead of falling to f64: the S4 opening. One plane, one node — sketches
    /// on one face share it automatically. The derivation convention is frozen spec (changing
    /// it would silently turn every stored sketch).
    Canonical,
    /// The caller's own origin and `+u` direction, stated in the frame the plane's data lives
    /// in — what a caller-named plane (`PlaneDef`) supplies, validated at its construction
    /// (origin on the plane, `ref_dir` not parallel to the normal). `ref_dir` need not be unit
    /// length nor exactly in-plane; the realization projects it exactly.
    Named { origin: [Rat; 3], ref_dir: [Rat; 3] },
}

/// A node in the motion-history forest (design §CIP ⑦): one [`Motion`] applied to a solid, with a
/// parent link so several points can share a history's tail.
///
/// Stored in [`Model::motions`]; a moved surface's
/// moved surfaces name their leaf node ([`SurfaceTruth`]'s motion slot). The tol a motion contributes is
/// application-point-dependent, so it is **not** stored here — judgment computes it by traversing
/// to the root.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MotionNode {
    pub motion: Motion,
    pub parent: Option<Handle<MotionNode>>,
}

/// What makes two surfaces the same plane, for `Model::surface_ids`: the canonical name
/// ([`nacre_scalar::PlaneName`] — `Narrow | Wide`, S2) and **the motion it is stated in**.
/// Identical names under different motions are different planes, because `Constructed` names
/// speak about the world and `Moved` ones about the pre-motion frame.
pub type SurfaceKey = (nacre_scalar::PlaneName, Option<Handle<MotionNode>>);

/// The statement key for a plane the name key cannot hold (open item 16): the sorted defining
/// triple and the motion it is stated under. See `Model::surface_through_ids`.
type ThroughKey = ([Handle<Vertex>; 3], Option<Handle<MotionNode>>);

/// The interning key for a cylinder — **deliberately conservative** (M6-0): the whole exact
/// statement, `ref_dir` included, plus the motion it is stated under. Two statements of one
/// geometric cylinder with different `ref_dir`s stay two handles, because merging them would
/// split the seam (seam vertices and the seam edge cite the surface as their carrier). A key
/// this literal cannot merge wrongly; geometric identity across different statements is the
/// predicates' to answer per question (rule 6), starting M6-1. No `flipped` report either — a
/// literal-identical statement realizes to a literal-identical cache.
type CylinderKey = (CylinderDef, Option<Handle<MotionNode>>);

/// How many planes were named **`Wide`** — the canonical answer exceeded `i128` and took the
/// arbitrary-precision vessel (S2). Before S2 these were the *unnamed* planes; now they intern
/// and carry identity like any other.
///
/// ★ **What `Wide` still cannot do**: host a sketch frame (`ops::frame_chain` reads
/// [`nacre_scalar::PlaneName::narrow`]) or ride the exact shortcuts — those decline exactly as
/// they declined on a missing name, until S4's `Canonical` placement opens frames without names.
/// The non-rotated coplanarity test reads the faces' own triangles and never looks at a name.
///
/// ★★ **Bounded by the type, not by the corpus.** A canonical name is a product of two point
/// differences, so `Rat = Ratio<i128>` inputs admit answers to roughly `2^2291`; today's models
/// stay far inside `i128` because they are written in short decimals, and the count reflects that
/// rather than any guarantee.
pub static WIDE_PLANES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// How many pushes interned onto a **seeded world plane** (handles 0–2) — the census stat that
/// explains a plane-digest diff the wide counter cannot (S9): a seeding-shaped change moves
/// survivors by *name collision*, not by width, and a falsifiability bridge needs a number for
/// that population too.
pub static SEEDED_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Whether a face uses its surface normal as-is (`Forward`) or flipped
/// (`Reversed`). A pure tag — full derives.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reversed,
}

impl Orientation {
    /// The opposite orientation.
    #[inline]
    pub fn flipped(self) -> Orientation {
        match self {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
        }
    }

    /// The flag as a sign: stored surface normal × `sign()` = the face's **stated outward**
    /// — the one reading props, STEP (`same_sense`) and the boolean engine all share. One
    /// spelling, because the ±1 map used to live inline at six production sites and a copy
    /// drifting is exactly how a face comes to lie about which way it faces.
    #[inline]
    pub fn sign(self) -> i8 {
        match self {
            Orientation::Forward => 1,
            Orientation::Reversed => -1,
        }
    }
}

/// **A surface's exact truth** — what the surface *is*, as opposed to the f64
/// [`Surface`](nacre_geom::Surface) beside it, which is its realization (S6b,
/// `docs/truth-and-cache.md`).
///
/// The `motion` field is what the retired `SurfaceDef` used to say from a side table: `None` is the
/// world (`Constructed`), `Some` names the motion history the data is stated *before*
/// (`Moved`). Holding it inside the variant is the point — a surface whose provenance is
/// unrecorded, or whose exact form does not exist (`Inexact`), is **unrepresentable** here,
/// which is what retires both.
/// ★ **`Plane` is the large variant and it is not boxed** (the `Operation::Extrude`
/// precedent): planes dominate the arena — a prism is all planes, a cylinder contributes one
/// curved surface — so boxing the points would put an allocation and a pointer chase on the
/// common case to shrink the rare one.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
pub enum SurfaceTruth {
    Plane {
        /// The plane's three exact points, stated in the frame `motion` names (the world when
        /// `None`) — the value `Model::surface_points` used to carry.
        points: PlanePoints,
        /// The motion history carrying the points out to the world; `None` = the world itself.
        motion: Option<Handle<MotionNode>>,
    },
    /// A cylinder's exact truth (M6-0): the rational statement of its lateral surface, beside
    /// the same motion slot a plane carries — a *moved* cylinder records its history instead of
    /// silently degrading (the old side-table path demoted it to `Inexact`).
    Cylinder {
        /// The exact statement, in the frame `motion` names (the world when `None`).
        def: CylinderDef,
        /// See [`SurfaceTruth::Plane::motion`].
        motion: Option<Handle<MotionNode>>,
    },
}

/// A plane's three points — the kind of statement is the variant. `Known` carries values
/// (construction planes — walls, caps, caller-stated planes); `Through` points at model vertices
/// (datum planes — S5, `docs/truth-and-cache.md`).
///
/// ★ **`Known` is nine `Rat` (288 B) beside three handles (24 B), and it stays unboxed.** The
/// size gap is not new — every plane already pays it — and `Known` is very nearly the whole
/// population, so boxing it would put an allocation and an indirection on the common read in
/// order to shrink the rare one. Same reasoning `PlaneName` records for keeping `Wide` inline.
/// (Not measured; the shape of the trade is what decides it, as it did there.)
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
pub enum PlanePoints {
    /// Three exact rational points, non-collinear, in the pre-motion frame.
    Known([[nacre_scalar::Rat; 3]; 3]),
    /// **Three model vertices the plane passes through.** What a caller means by "the plane
    /// through those corners" — a coordinate read off a discovered vertex is rounded, and the
    /// plane built from rounded coordinates is a *different* plane (measured: on tilted geometry
    /// every one of 220 triples, `nacre-ops/tests/point_width.rs`).
    ///
    /// ★ **Sorted**, so that the same three vertices are the same statement whichever order they
    /// arrive in. Direction is *not* lost by sorting: a plane's canonical name carries none
    /// either, and the caller's order is measured into the frame's `flip` where every other
    /// producer already puts it.
    ///
    /// ★★ **The motion composes rather than replacing.** The vertices pin the plane in the frame
    /// they are stated in and `motion` carries it out — the two add, never multiply. A mover must
    /// therefore keep these handles verbatim and record a node; transporting *and* keeping them
    /// would move the plane twice, and doing neither would leave the truth behind while the cache
    /// moved.
    Through([nacre_store::Handle<Vertex>; 3]),
}

/// **A cylinder's exact truth** (M6-0): the lateral surface as its producer stated it — origin,
/// axis direction, seam reference direction and radius, all rational, in the pre-motion frame.
///
/// ★ `dir` and `ref_dir` are **raw, unnormalized** — the `normal_def` precedent: normalizing
/// divides by an irrational length and would destroy the exact form. The f64 cache
/// ([`nacre_geom::Cylinder`]) holds the realized unit frame; this holds what the cylinder *is*.
/// The component form is the cylinder's analogue of a plane's *points* — a direct lift of the
/// caller's vocabulary plus component shuffles, never a derived product (the shape the
/// falsification table forbids for coefficients).
///
/// ★★ **`ref_dir` is model geometry, not a chart choice.** It fixes the seam permanently
/// (`θ = 0` on the `+ref_dir` side of the axis); the rational half-angle chart (M6-1) places its
/// own excluded point *on* this seam — the chart adapts to the seam, never the reverse.
///
/// Construction is checked ([`CylinderDef::new`]); fields are private so a literal cannot bypass
/// the check (the `SketchFrame` precedent).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct CylinderDef {
    origin: [Rat; 3],
    dir: [Rat; 3],
    ref_dir: [Rat; 3],
    radius: Rat,
}

impl CylinderDef {
    /// The checked constructor — `None` when the statement means no cylinder, and **only then**:
    /// a zero `dir`, a non-positive `radius`, or a `ref_dir` with no component perpendicular to
    /// the axis (`ref_dir × dir = 0`, which a zero `ref_dir` satisfies too).
    ///
    /// ★ **Width is not a cause.** The parallelism test runs in
    /// [`nacre_scalar::parallel_rat`], which clears denominators and answers in integers, so it
    /// cannot decline. It used to run in checked `Rat` and answer `None` on overflow — a
    /// "conservative refusal" that conflated *no cylinder* with *the arithmetic ran out*, and
    /// the callers below read it as the first: a statement whose axis carries a small component
    /// with a long decimal (denominator ~10²⁰, whose square leaves `i128`) crashed the
    /// constructor's `expect`. Measured population: 80% of computed near-axis-aligned
    /// directions, 0% of hand-written short decimals.
    pub fn new(origin: [Rat; 3], dir: [Rat; 3], ref_dir: [Rat; 3], radius: Rat) -> Option<Self> {
        let zero = Rat::from_int(0);
        if dir.iter().all(|c| *c == zero) || radius <= zero {
            return None;
        }
        if nacre_scalar::parallel_rat(&ref_dir, &dir) {
            return None;
        }
        Some(CylinderDef {
            origin,
            dir,
            ref_dir,
            radius,
        })
    }

    /// A point on the axis, exact.
    pub fn origin(&self) -> [Rat; 3] {
        self.origin
    }

    /// The axis direction, raw (unnormalized), exact. Nonzero by construction.
    pub fn dir(&self) -> [Rat; 3] {
        self.dir
    }

    /// The seam reference direction, raw, exact — not parallel to `dir` by construction. The
    /// seam is the lateral line on the `+ref_dir` side of the axis.
    pub fn ref_dir(&self) -> [Rat; 3] {
        self.ref_dir
    }

    /// The radius, exact. Positive by construction.
    pub fn radius(&self) -> Rat {
        self.radius
    }
}

/// Why [`Model::add_cylinder_exact`] refused a statement — **every way in, named**.
///
/// The entry promises no panics: an application's numbers are input, not a caller bug, and a
/// panic in wasm is a dead session rather than a sentence. The first three are the caller's
/// statement; the last two cannot happen for a statement that passed them (see the entry's doc)
/// and exist so that "cannot happen" never has to be spelled `expect`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CylinderError {
    /// `axis` or `ref_dir` is not a unit vector, or the two are not perpendicular. This is the
    /// precondition that makes every point the entry derives exact.
    FrameNotOrthonormal,
    NonPositiveRadius,
    NonPositiveHeight,
    /// An exact product left `i128` — the statement is representable, its derived points are not.
    Overflow,
    /// The truth, a cache, or an edge's curve refused what the statement said.
    Degenerate,
}

/// Everything a cylinder solid is assembled from, already derived: the **exact truth** (`def`,
/// the caps' points) and its **realization** (the caches, the seam coordinates).
///
/// ★ The two cylinder constructors differ only in *how they reach here* — a caller's f64
/// statement (`Model::add_cylinder`) or an exact orthonormal frame
/// ([`Model::add_cylinder_exact`]) — so the b-rep below is written once. The f64 road cannot
/// simply delegate to the exact one: its axis comes out of `normalize()`, and a normalized f64
/// direction has no rational form to hand over.
struct CylinderParts {
    def: CylinderDef,
    lateral: Cylinder,
    /// Bottom then top: each cap's realized plane beside its three exact points.
    caps: [(Plane, [[Rat; 3]; 3]); 2],
    /// The rims' `θ = 0` points (bottom, top) — where the seam vertices sit.
    seam_pts: [Point3; 2],
    motion: Option<Handle<MotionNode>>,
}

/// The truth-only aggregate: exact geometry + topology stores + the derived
/// adjacency cache. No `tess`, no `ops` (see the crate docs).
#[derive(Debug)]
pub struct Model {
    // exact geometry (truth)
    /// ★ Private (stage S1, `docs/truth-and-cache.md`): a surface can only enter through
    /// [`Model::push_plane`]/[`Model::push_cylinder`], which state its truth —
    /// a surface **without** a record is unrepresentable from outside this crate. Read through
    /// [`Model::surface`]/[`Model::surface_count`]; there is deliberately no whole-store
    /// iterator (the arena keeps superseded surfaces — consumers walk the live faces).
    surfaces: Store<Surface>,
    /// ★★ **The truth beside the cache** (S6b): index-parallel to [`Model::surfaces`], so a
    /// `Handle<Surface>` names both — the store above holds the f64 *realization*, this holds
    /// what the surface *is*. Total: a surface cannot enter without its truth, which is what
    /// retired `SurfaceDef`/`Inexact` and the point-less population.
    ///
    /// (Why the truth is the `Vec` and the cache the `Store`, when the doc draws it the other
    /// way: `Handle<T>`'s type parameter. Retyping `Face::surface` would ripple through every
    /// crate for a distinction the shared index already erases — the flip happens with the
    /// final rename, when `SurfaceCache` gets its real shape.)
    surface_truths: Vec<SurfaceTruth>,
    /// The three seeded world planes, in normal-axis order Z(XY)·X(YZ)·Y(ZX) — captured at
    /// [`Model::new`] so [`Model::world_plane`] needs no handle minting. Always length 3.
    world_planes: Vec<Handle<Surface>>,
    /// Per-edge curve caches, index-parallel to `edges` (S8) — **cache, not truth**: the
    /// carriers and endpoints decide the curve ([`Model::derive_edge_curve`]), and
    /// [`Model::rebuild_edge_cache`] discards and regenerates the lot. Filled eagerly by
    /// [`Model::push_edge`]; read through [`Model::edge_curve`]. A raw `edges.push` without a
    /// cache entry desyncs the two — the accessor's debug_assert and validate's parallelism
    /// check watch for that (the store stays `pub` per the doc's final shape).
    edge_cache: Vec<EdgeCache>,
    /// Per-vertex coordinate caches, index-parallel to `vertices` (S7) — the realized `coord`
    /// and, for discovered vertices, the measured `tol`. Filled by [`Model::push_vertex`];
    /// read through [`Model::vertex_point`] / [`Model::vertex_tol`]. No rebuild exists (3b ⏸):
    /// for discovered and seam vertices the coordinate carries information the definition
    /// cannot yet reproduce.
    vertex_cache: Vec<PointCache>,
    /// The motion-history forest (§CIP ⑦): motion definitions named by moved surfaces. Not geometry — a definition store.
    ///
    /// **Interned** — private (S1), so writing through [`Model::push_motion`] is enforced by the
    /// type, not by discipline. Read through [`Model::motion`].
    motions: Store<MotionNode>,
    /// Interning table for [`Model::push_motion`]: the handle already issued for a given
    /// `(motion, parent)`. Not iterated (a `HashMap`'s order must never reach a result).
    motion_ids: HashMap<MotionNode, Handle<MotionNode>>,
    /// Each surface's **canonical name** ([`nacre_scalar::PlaneName`]), derived from its points —
    /// present for every surface whose producer had a rational description to record.
    ///
    /// ★ **The point is that two statements of one plane get the same value.** `nacre_geom::Plane`
    /// keeps an un-normalized normal whose length follows the *face's size*, so the same plane
    /// reaches `coefficients()` as `[2.2, 0, 0, −6.6000000000000005]` from one face and
    /// `[13.2, 0, 0, −39.599999999999994]` from another — not exactly proportional, because `d`
    /// is a rounded product. Canonical names have no scale to disagree about, and since S2 the
    /// vessel is arbitrary-precision (`Narrow | Wide`) so **width cannot lose a name either**.
    ///
    /// ★★ **Built from the dimensions the user wrote, never lifted from the f64 coefficients
    /// above.** Lifting is lossless and useless here: it preserves the rounding, so the two
    /// vectors stay different (`nacre_scalar::canonical_plane_coeffs`).
    ///
    /// ★★★ **The name is stated in the frame this surface's truth names** — the world for
    /// `motion: None`, and the **pre-motion** frame for a moved surface, whose world
    /// coefficients are irrational and so cannot be written
    /// down at all. A moved surface therefore inherits its source's name unchanged: the motion
    /// is recorded beside it, not folded into it. (That pairing — an exact name plus a
    /// motion — is the shape `docs/truth-and-cache.md` builds toward.)
    ///
    /// Absent is ordinary — an f64 construction path (no points to derive from).
    /// Iterate through the faces, never over the map. Arithmetic consumers (frames, `base_rat`,
    /// exact transports) read [`nacre_scalar::PlaneName::narrow`]; `Wide` carries identity only.
    pub surface_name: HashMap<Handle<Surface>, nacre_scalar::PlaneName>,
    /// Interning table for [`Model::push_plane`]: the handle already issued for a
    /// plane, keyed by its canonical coefficients **and the motion they are stated in**. The twin
    /// of [`Model::motion_ids`]; not iterated (a `HashMap`'s order must never reach a result).
    ///
    /// ★ **The motion belongs in the key.** `Constructed` coefficients speak about the world and
    /// `Moved` ones about the pre-motion frame, so two identical arrays under different motions
    /// are different planes. `Inexact` has no coefficients and never interns.
    ///
    /// ★★ **The witness does not.** Two faces of one plane sharing one moved
    /// witness is already how this works — `collect_planes` says the witness "was captured from
    /// whichever face first reached this surface" and winds it to *each* face's own outward
    /// normal — so it is not part of what makes two planes the same.
    surface_ids: HashMap<SurfaceKey, Handle<Surface>>,
    /// Interning for the planes the name key **cannot** hold: a `Through` statement whose exact
    /// world coefficients are irrational (mixed-frame datum — open item 16) has no canonical
    /// name, so it interns by the **statement itself**: the sorted vertex triple and the motion
    /// it is stated under.
    ///
    /// ★ This is statement identity, not geometric identity. Two *different* triples on one
    /// geometric plane get two handles here, and rule 6's qualification is exactly that:
    /// a nameless plane's geometric identity is the predicates' to answer, per question.
    /// What this table guarantees is the same thing construction-time sorting guarantees one
    /// level down — **the same statement never becomes two handles.**
    surface_through_ids: HashMap<ThroughKey, Handle<Surface>>,
    /// Interning for cylinders (M6-0) — by the whole exact statement plus motion; see
    /// [`CylinderKey`] for why the key is deliberately this literal.
    cylinder_ids: HashMap<CylinderKey, Handle<Surface>>,
    // topology (references geometry by Handle only)
    pub vertices: Store<Vertex>,
    pub edges: Store<Edge>,
    pub faces: Store<Face>,
    pub shells: Store<Shell>,
    pub solids: Store<Solid>,
    /// The live solids — the "current model" (design §2 supersede semantics).
    /// Editing ops supersede topology by pushing new cells and updating this
    /// list; the old cells stay in the append-only arena but, unreferenced by
    /// any live solid, drop out of the reachable closure. Producers register
    /// through [`Model::push_solid`]; `validate`/`Adjacency`/`nacre-step`
    /// traverse [`Model::reachable`], not the whole store.
    pub live_solids: Vec<Handle<Solid>>,
    // derived cache (rebuilt on demand)
    pub adj: Adjacency,
}

/// One vertex's realized coordinate and measured tolerance — a **cache** beside the vertex
/// store (index-parallel), filled by [`Model::push_vertex`]. `tol: Some` is a discovered
/// vertex's measured residual (kept exact — `0.0` means exactly zero); `None` is a constructed
/// vertex (checkers apply their own epsilon). Unlike [`EdgeCache`] there is **no rebuild**:
/// the coordinate is truth-bearing for seam vertices and hard-won for discovered ones (3b ⏸).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointCache {
    coord: Point3,
    tol: Option<f64>,
}

/// One edge's realized curve — a **cache** beside the edge store (index-parallel), derived
/// from the edge's carriers and endpoints by [`Model::derive_edge_curve`]. Discard and
/// regenerate with [`Model::rebuild_edge_cache`]; the field is private so the only writers
/// are the derivation itself.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeCache {
    curve: Curve,
}

/// The handles reachable from a model's live solids — the live model (design §2).
///
/// Only the sets consumers need today: `validate`'s Euler counts vertices/edges/
/// faces/shells, and `Adjacency`/loop/incidence walk faces/edges. Surfaces and
/// curves are not tracked: the M4 face ops never orphan a *used* one, and the three
/// seeded world planes (S9) are deliberately face-less — traversal through
/// definitions (world planes, datum references) is 열린 항목 4.
#[derive(Debug, Default)]
pub struct Reachable {
    pub vertices: HashSet<Handle<Vertex>>,
    pub edges: HashSet<Handle<Edge>>,
    pub faces: HashSet<Handle<Face>>,
    pub shells: HashSet<Handle<Shell>>,
}

/// Whether a handle indexes inside its store (guards the reachable traversal
/// against dangling/out-of-range handles).
#[inline]
fn in_bounds<T>(h: Handle<T>, store: &Store<T>) -> bool {
    (h.index() as usize) < store.len()
}

impl Loop {
    /// This loop with its winding reversed: half-edges in reverse order, each traversed the
    /// opposite way. Both parts matter — reversing the *order* flips the normal the winding
    /// implies, and flipping each `forward` keeps every edge used once in each direction, so two
    /// faces that share an edge stay opposed when both are reversed (still a valid 2-manifold).
    ///
    /// Two callers want it for opposite reasons. [`Model::reversed_shell`] pairs it with an
    /// [`Orientation`] toggle to turn a boundary inward (a cavity). A **reflection** pairs it with
    /// nothing: mirroring negates the normal a winding implies, so rewinding restores it and the
    /// orientation flag stays as it was.
    pub fn reversed(&self) -> Loop {
        Loop {
            half_edges: self
                .half_edges
                .iter()
                .rev()
                .map(|he| HalfEdge {
                    edge: he.edge,
                    forward: !he.forward,
                })
                .collect(),
        }
    }
}

/// `Default` is [`Model::new`] — **seeded**. The derive used to hand out an *empty* model,
/// which after S9 would be one without the world planes: a public back door to the state the
/// seeding exists to remove. There is deliberately no unseeded constructor.
impl Default for Model {
    fn default() -> Self {
        Self::new()
    }
}

impl Model {
    /// A model with the three **world axis planes pre-seeded** — surface handles 0 (XY, z = 0),
    /// 1 (YZ, x = 0), 2 (ZX, y = 0), deterministic so a replayed log and a live session name the
    /// same planes (`docs/truth-and-cache.md` rule: 세계 축 평면 셋은 `Model::new()` 가 심는다).
    ///
    /// Each seed's truth is the canonical triple `[0, u, v]` — the very points
    /// `SketchPlane::world_*` states — so any later producer of the same plane (a cuboid face on
    /// an axis, a `z = 0` sketch's base cap) **interns onto the seed**: one plane, one handle,
    /// stated once.
    ///
    /// ★ **The seed's f64 cache points down the −axis**, not up. This is the sense the dominant
    /// producers push — `extrude` builds base caps with `−normal` (its comment records that
    /// pushing `+normal` "flipped the stored normal on 781 base caps … for no reason"), and an
    /// origin cuboid's bottom/left/front faces cross to `−axis` raws — so a −axis seed keeps
    /// their `flipped` bits false and their `Orientation` spellings unchanged. The *sketch*
    /// convention ("`world_xy`'s normal is `+ẑ`") is about the frame, not the stored cache;
    /// frame derivation goes through the canonical name plus a measured `flip`, so the cache's
    /// direction never leaks to a caller.
    pub fn new() -> Self {
        let mut m = Model {
            surfaces: Store::default(),
            surface_truths: Vec::new(),
            edge_cache: Vec::new(),
            vertex_cache: Vec::new(),
            motions: Store::default(),
            motion_ids: HashMap::new(),
            surface_name: HashMap::new(),
            surface_ids: HashMap::new(),
            surface_through_ids: HashMap::new(),
            cylinder_ids: HashMap::new(),
            vertices: Store::default(),
            edges: Store::default(),
            faces: Store::default(),
            shells: Store::default(),
            solids: Store::default(),
            live_solids: Vec::new(),
            adj: Adjacency::default(),
            world_planes: Vec::new(),
        };
        let r = nacre_scalar::Rat::from_int;
        // (normal axis, +u, +v) for XY / YZ / ZX — the `SketchPlane::axis_plane` triples.
        let seeds: [([f64; 3], [i128; 3], [i128; 3]); 3] = [
            ([0.0, 0.0, 1.0], [1, 0, 0], [0, 1, 0]),
            ([1.0, 0.0, 0.0], [0, 1, 0], [0, 0, 1]),
            ([0.0, 1.0, 0.0], [0, 0, 1], [1, 0, 0]),
        ];
        let world_planes: Vec<Handle<Surface>> = seeds
            .into_iter()
            .map(|(n, u, v)| {
                let cache = nacre_geom::Plane::from_point_normal(
                    nacre_math::Point3::origin(),
                    nacre_math::Vector3::from_array(n.map(|c| -c)),
                )
                .expect("a unit axis");
                let points = [[r(0); 3], u.map(r), v.map(r)];
                let (h, flipped) = m.push_plane(cache, points, None);
                debug_assert!(!flipped, "an empty model cannot intern a seed");
                h
            })
            .collect();
        m.world_planes = world_planes;
        debug_assert_eq!(
            m.world_planes.iter().map(|h| h.index()).collect::<Vec<_>>(),
            [0, 1, 2],
            "seed handles are deterministic"
        );
        m
    }

    /// Recompute the [`Adjacency`] cache from the current topology stores.
    /// Call once after a batch of additions (the cache is otherwise stale).
    pub fn rebuild_adjacency(&mut self) {
        // Build against an immutable borrow, then move into place — avoids
        // borrowing `self.adj` mutably while iterating the other stores.
        let adj = Adjacency::rebuild(&*self);
        self.adj = adj;
    }

    /// Push a solid into the store **and mark it live** (design §2). This is the
    /// blessed way for a producer to add a solid; the reachable closure
    /// ([`Model::reachable`]) grows to include it. Editing ops instead mutate
    /// [`Model::live_solids`] directly (drop the superseded solid, add the new).
    pub fn push_solid(&mut self, solid: Solid) -> Handle<Solid> {
        let h = self.solids.push(solid);
        self.live_solids.push(h);
        h
    }

    /// Push a motion node, **interned**: the same `(motion, parent)` always yields the same
    /// handle.
    ///
    /// The handle is the canonical name of "which motion", and judgments use it to decide whether
    /// a set of points shares one rigid motion — which lets the whole judgement be answered
    /// exactly in the pre-motion frame. That identity has to be *both* collision-free and free of
    /// false misses:
    ///
    /// - it used to be a 64-bit hash of the chain's contents, and a collision would hand the exact
    ///   predicate two incompatible frames and answer a different question with full confidence;
    /// - a raw `Store::push` per transform is collision-free but *misses*: turning two solids by
    ///   the same 30° would make two nodes, and their shared motion would stop cancelling — so
    ///   rotating a model would turn its exact questions into assumed ones (measured: it breaks
    ///   `a_shared_rotation_still_assumes_nothing`).
    ///
    /// Interning gives both. It also keeps the forest small, since a chain shared by many solids
    /// is stored once.
    pub fn push_motion(
        &mut self,
        motion: Motion,
        parent: Option<Handle<MotionNode>>,
    ) -> Handle<MotionNode> {
        let node = MotionNode { motion, parent };
        if let Some(&h) = self.motion_ids.get(&node) {
            return h;
        }
        let h = self.motions.push(node);
        self.motion_ids.insert(node, h);
        h
    }

    /// The one push everything funnels through (private): the f64 cache and the exact truth,
    /// index-parallel, in one motion — so the two stores cannot come apart.
    fn push_raw(&mut self, surface: Surface, truth: SurfaceTruth) -> Handle<Surface> {
        let h = self.surfaces.push(surface);
        self.surface_truths.push(truth);
        debug_assert_eq!(self.surface_truths.len(), self.surfaces.len());
        h
    }

    /// The surface's exact truth — what it *is*, beside the f64 realization
    /// [`Model::surface`] returns. Total: a surface without a truth entry is unrepresentable
    /// (S6b), which is what retired `SurfaceDef::Inexact` and the `UndefinedSurface` violation.
    #[inline]
    pub fn surface_truth(&self, h: Handle<Surface>) -> &SurfaceTruth {
        &self.surface_truths[h.index() as usize]
    }

    /// The seeded world plane whose **normal** runs along `axis` — `Z` names the XY plane
    /// (z = 0), `X` the YZ plane, `Y` the ZX plane. Deterministic (handles 0–2, pushed by
    /// [`Model::new`]), so a handle in an op log and one from a live session agree.
    pub fn world_plane(&self, axis: nacre_scalar::Axis) -> Handle<Surface> {
        let ix = match axis {
            nacre_scalar::Axis::Z => 0usize,
            nacre_scalar::Axis::X => 1,
            nacre_scalar::Axis::Y => 2,
        };
        let h = self.world_planes[ix];
        debug_assert!(
            matches!(
                self.surface_truth(h),
                SurfaceTruth::Plane { motion: None, .. }
            ),
            "seed handles must stay the world planes"
        );
        h
    }

    /// The surface a handle names — the **f64 cache** of [`Model::surface_truth`]'s answer.
    ///
    /// Reading is open; **writing is not** — the store is private (S1), so a surface can only
    /// enter through [`Model::push_plane`]/[`Model::push_cylinder`], which state its truth. The
    /// lock:
    ///
    /// ```compile_fail,E0616
    /// let m = nacre_topo::Model::new();
    /// let _ = m.surfaces.len(); // private field — read through `surface`/`surface_count`
    /// ```
    #[inline]
    pub fn surface(&self, h: Handle<Surface>) -> &Surface {
        self.surfaces.get(h)
    }

    /// How many surfaces the arena holds — live and superseded alike. Handle-validity checks
    /// and the `a_boolean_mints_no_surface` lock read this; nothing iterates the store
    /// (superseded surfaces are still in it — consumers walk the live faces).
    #[inline]
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// The surface an **index** names — how a log's index vocabulary is re-anchored onto the model
    /// being built (`docs/design.md` §2).
    ///
    /// A handle in an operation log carries only its index across models, so `replay` turns that
    /// index back into a handle of its own arena before applying the operation. `None` past the
    /// end: existence, not legality.
    ///
    /// **This does not open the store.** Reading was already open ([`Model::surface`],
    /// [`Model::surface_count`]); writing still goes only through [`Model::push_plane`] /
    /// [`Model::push_cylinder`], which state the truth. The `compile_fail` lock on
    /// [`Model::surface`] is untouched.
    ///
    /// The one legitimate shape is "re-anchor a log's index onto the model I am building" — using
    /// it to quiet a cross-model panic hides the bug instead of fixing it. Its consumer today is
    /// `nacre-ops`' `rebind`, for the operations that name a plane.
    #[inline]
    #[must_use]
    pub fn surface_handle_at(&self, index: u32) -> Option<Handle<Surface>> {
        self.surfaces.handle_at(index)
    }

    /// The vertex a log's index names — [`Model::surface_handle_at`]'s twin, for the same reason
    /// (`replay` re-anchors an operation's handles onto the model it is rebuilding).
    ///
    /// ★ A bounds check is all this can be, and for vertices that is a weaker guarantee than it
    /// looks: a rejected operation still leaves cells behind (measured 63–84), so a log recorded
    /// across a reject can re-anchor **in range and onto the wrong vertex**. The discipline that
    /// answers it — rebuild from the log before recording again — lives with `replay`, and the
    /// generated-session proptest is what holds it.
    #[inline]
    pub fn vertex_handle_at(&self, index: u32) -> Option<Handle<Vertex>> {
        self.vertices.handle_at(index)
    }

    /// The motion node a handle names. Same seal as [`Model::surface`]: writing goes through
    /// [`Model::push_motion`] (interned), reading through here.
    #[inline]
    pub fn motion(&self, h: Handle<MotionNode>) -> &MotionNode {
        self.motions.get(h)
    }

    /// Push a **plane**, stating its truth outright: the f64 cache, three exact points, and the
    /// motion they are written before (`None` = the world). The truth is not optional — that is
    /// the S6b point: a point-less plane, and with it `SurfaceDef::Inexact`, stopped existing.
    ///
    /// ★★★★★ **The points are the only thing a producer states.** The canonical name
    /// ([`Model::surface_name`]) is *derived* here, from those points, by
    /// [`nacre_scalar::plane_name_exact`] — so a plane cannot be described two ways, because there
    /// is only one place to describe it. (Producers used to hand in coefficients beside the
    /// points; before that parameter went away, every producer was compared against this
    /// derivation across the suite: **83,883 agreements, 0 disagreements**.)
    ///
    /// ★ **Stated in the frame `motion` names** — the world for `None`, the pre-motion frame
    /// otherwise. An interned plane keeps the first pusher's triple.
    ///
    /// ★ **The `bool` says the returned surface's cache normal points the *other* way** from the
    /// one handed in, and a caller that meets it must record its face `Orientation::flipped()`.
    /// It exists because a plane's canonical form has no direction — `[0,0,1,−3]` and
    /// `[0,0,−1,3]` are one plane — so once identical planes share a handle the direction has to
    /// be reconciled somewhere, and the honest place is where the caller still knows what it
    /// asked for.
    ///
    /// ★★ **It is not a dormant path.** Measured across the suite, interning hits 7,095 times
    /// and **1,916 of those report `flipped`** — a boss meeting the plate it sits on is one
    /// plane approached from both sides, which is as ordinary as it sounds.
    pub fn push_plane(
        &mut self,
        cache: nacre_geom::Plane,
        points: [[nacre_scalar::Rat; 3]; 3],
        motion: Option<Handle<MotionNode>>,
    ) -> (Handle<Surface>, bool) {
        // ★★★★★ **The name is derived, so it cannot disagree with the thing it names.**
        // `plane_name_exact` computes the canonical form at unbounded precision — `None` only
        // for collinear points (which no production `PlaneDef` can supply); since S2 the vessel
        // (`PlaneName::Narrow | Wide`) always holds the answer, so every plane interns, wide
        // ones included. [`WIDE_PLANES`] counts the names that took the wide vessel.
        // ★★★★★ **The name is derived, so it cannot disagree with the thing it names.**
        // `plane_name_exact` computes the canonical form at unbounded precision — `None` only
        // for collinear points (which no production `PlaneDef` can supply); since S2 the vessel
        // (`PlaneName::Narrow | Wide`) always holds the answer, so every plane interns, wide
        // ones included.
        let name = nacre_scalar::plane_name_exact(points[0], points[1], points[2]);
        self.intern_plane(cache, name, PlanePoints::Known(points), motion)
    }

    /// **Interning, once — the half every plane producer shares.**
    ///
    /// A producer differs only in *how it derives the name* and *which `PlanePoints` it stores*.
    /// Everything after that — the key, the already-issued reply and its `flipped`, the arena
    /// push, the two side tables, the two counters — is the same, and was duplicated once, which
    /// promptly cost both counters on the new road ([`WIDE_PLANES`] and [`SEEDED_HITS`] were
    /// simply absent from it). Sharing the tail makes losing them structurally impossible rather
    /// than a thing to remember.
    fn intern_plane(
        &mut self,
        cache: nacre_geom::Plane,
        name: Option<nacre_scalar::PlaneName>,
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,
    ) -> (Handle<Surface>, bool) {
        if name.as_ref().is_some_and(|n| n.narrow().is_none()) {
            WIDE_PLANES.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
        let key = name.map(|n| (n, motion));
        if let Some(k) = &key {
            if let Some(&h) = self.surface_ids.get(k) {
                if (h.index() as usize) < 3 {
                    SEEDED_HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                // Same plane, already issued. The canonical form says nothing about direction, so
                // report whether the survivor points the other way and let the caller spell its
                // outward the other way round.
                return (h, self.flipped_against(h, &cache));
            }
        }
        let h = self.push_raw(
            Surface::Plane(cache),
            SurfaceTruth::Plane { points, motion },
        );
        if let Some((n, _)) = &key {
            // One clone per push — the name is derived once here, never on a judging loop.
            self.surface_name.insert(h, n.clone());
        }
        if let Some(k) = key {
            self.surface_ids.insert(k, h);
        }
        (h, false)
    }

    /// Whether the already-issued surface's cache normal points the other way from the one the
    /// caller just built — the `flipped` report both interning roads share. The f64 cache dot is
    /// exact here: two caches of one plane have parallel normals, so the sign cannot be lost to
    /// rounding.
    fn flipped_against(&self, h: Handle<Surface>, cache: &nacre_geom::Plane) -> bool {
        match self.surfaces.get(h) {
            Surface::Plane(p) => p.normal().dot(cache.normal()) < 0.0,
            Surface::Cylinder(_) => false,
        }
    }

    /// **Push a plane stated as the three vertices it passes through** — [`push_plane`]'s twin
    /// for the datum vocabulary, with the same interning contract and the same `flipped` report.
    ///
    /// The name is derived the same way, one step further back: solve each vertex from its three
    /// carriers, then take the canonical form of the plane through those points. So a `Through`
    /// plane and a `Known` plane that *are* the same plane share one handle, which is the whole
    /// point of interning — the variant records how this plane's existence is grounded, not who
    /// asked for it first.
    ///
    /// ★ `vertices` is **sorted** by the caller before it gets here (the same three vertices are
    /// the same statement in any order). Direction is not lost: it lives in the frame's measured
    /// `flip`, exactly as it does for a stated plane.
    ///
    /// ★★ **A statement the name key cannot hold still interns — by the statement itself**
    /// (open item 16). A mixed-frame datum's exact world coefficients are irrational, so
    /// [`Model::plane_name_through`] answers `None`; such a plane takes the second key
    /// (`surface_through_ids`) — the sorted triple and the motion. That is *statement*
    /// identity: the same three vertices under the same motion are one handle, and geometric
    /// identity across different statements is the predicates' to answer per question (rule 6's
    /// own qualification — *"interning 없이 술어가 매번 답한다"*). This is **not** the
    /// record-less population S2 drained: the truth (handles + motion) is complete; what does
    /// not exist is a rational description of it.
    ///
    /// ★ The producer remains responsible for rejecting **before** pushing whatever it cannot
    /// frame — this door stores; it does not validate framability.
    ///
    /// [`push_plane`]: Model::push_plane
    pub fn push_plane_through(
        &mut self,
        cache: nacre_geom::Plane,
        vertices: [Handle<Vertex>; 3],
        motion: Option<Handle<MotionNode>>,
    ) -> (Handle<Surface>, bool) {
        debug_assert!(
            vertices[0].index() < vertices[1].index() && vertices[1].index() < vertices[2].index(),
            "a Through statement must arrive sorted and duplicate-free"
        );
        let name = self.plane_name_through(vertices);
        if name.is_some() {
            return self.intern_plane(cache, name, PlanePoints::Through(vertices), motion);
        }
        if let Some(&h) = self.surface_through_ids.get(&(vertices, motion)) {
            return (h, self.flipped_against(h, &cache));
        }
        let h = self.push_raw(
            Surface::Plane(cache),
            SurfaceTruth::Plane {
                points: PlanePoints::Through(vertices),
                motion,
            },
        );
        self.surface_through_ids.insert((vertices, motion), h);
        (h, false)
    }

    /// **The name a `Through` statement derives**, and the one place that derivation lives — the
    /// producer's check and [`Model::push_plane_through`] read the same answer, so "we rejected
    /// what we could not name" is structural rather than two functions agreeing by habit.
    ///
    /// `None` when any vertex is not a three-plane point, when the carriers do not share one
    /// motion (no frame holds a rational coordinate then), or when the three points are
    /// collinear. ★ A meet too wide for `Rat` is **not** on the list since open item 17: the
    /// name is derived from the meets at whatever width they need
    /// ([`nacre_scalar::plane_name_from_meets`]) — width was the arithmetic's problem, never
    /// the statement's.
    pub fn plane_name_through(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<nacre_scalar::PlaneName> {
        let m = self.through_meets(vertices)?;
        nacre_scalar::plane_name_from_meets([&m[0], &m[1], &m[2]])
    }

    /// **The three vertices' exact meeting points, in the one frame they share** — the single
    /// solve behind [`Model::plane_name_through`] (at push, width-free) and, through the
    /// all-narrow projection [`Model::through_points_rat`], the judging table's witness
    /// triangle. One spelling, so the name a plane interns under and the points a predicate
    /// reasons about cannot describe different planes.
    ///
    /// `None` on any of: a vertex that is not a three-plane point, carriers that do not share
    /// one motion (then no frame holds a rational coordinate at all), or a carrier with no
    /// recorded name.
    pub fn through_meets(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<[nacre_scalar::MeetPoint; 3]> {
        let mut pts: [Option<nacre_scalar::MeetPoint>; 3] = [None, None, None];
        let mut frame = None;
        for (i, vh) in vertices.iter().enumerate() {
            let tri = match self.vertices.get(*vh).def {
                VertexDef::ThreePlane(tri) => tri,
                // OnSeam pins a curve, not a point; a Branch *is* a point but its coordinates
                // are quadratic-irrational — neither has the rational meet a datum statement
                // needs, so both decline here (honest, and spelled per variant so the next
                // variant is a compile error, not a silent fall-through).
                VertexDef::OnSeam(_) | VertexDef::Branch { .. } => return None,
            };
            let mine = self.plane_motion(tri[0]);
            if tri.iter().any(|h| self.plane_motion(*h) != mine) {
                return None; // this vertex's carriers straddle frames
            }
            match frame {
                None if i == 0 => frame = Some(mine),
                f if f == Some(mine) => {}
                _ => return None, // the three vertices do not share one frame
            }
            let names = tri.map(|h| self.surface_name.get(&h));
            let [Some(a), Some(b), Some(c)] = names else {
                return None;
            };
            pts[i] = Some(nacre_scalar::three_planes_big([a, b, c])?);
        }
        Some([pts[0].take()?, pts[1].take()?, pts[2].take()?])
    }

    /// [`Model::through_meets`]' all-narrow projection — the form a **witness triangle** takes,
    /// since a witness base is a `[Rat; 3]` by type. `None` additionally when any meet is
    /// [`nacre_scalar::MeetPoint::Wide`]; the judging table then builds its witness another way
    /// (the plane's own frame probes), so this is a road fork, not a refusal.
    pub fn through_points_rat(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<[[nacre_scalar::Rat; 3]; 3]> {
        let meets = self.through_meets(vertices)?;
        let mut pts = [[nacre_scalar::Rat::from_int(0); 3]; 3];
        for (o, m) in pts.iter_mut().zip(&meets) {
            *o = *m.narrow()?;
        }
        Some(pts)
    }

    /// The motion a surface's truth records, whichever variant it is.
    #[inline]
    pub fn plane_motion(&self, h: Handle<Surface>) -> Option<Handle<MotionNode>> {
        match self.surface_truth(h) {
            SurfaceTruth::Plane { motion, .. } | SurfaceTruth::Cylinder { motion, .. } => *motion,
        }
    }

    /// Push a **cylinder** — the lateral surface, stating its exact truth (M6-0), with
    /// interning by the whole statement (see [`CylinderKey`] for why the key is deliberately
    /// that literal: merging two `ref_dir`s would split the seam). No `flipped` report — a
    /// literal-identical statement realizes to a literal-identical cache, so there is no other
    /// way round to report.
    pub fn push_cylinder(
        &mut self,
        cache: nacre_geom::Cylinder,
        def: CylinderDef,
        motion: Option<Handle<MotionNode>>,
    ) -> Handle<Surface> {
        let key = (def.clone(), motion);
        if let Some(&h) = self.cylinder_ids.get(&key) {
            return h;
        }
        let h = self.push_raw(
            Surface::Cylinder(cache),
            SurfaceTruth::Cylinder { def, motion },
        );
        self.cylinder_ids.insert(key, h);
        h
    }

    /// Push a plane with truth but **no name and no interning** — test-only.
    ///
    /// Two fixture populations need this door: hand-built merge fixtures that deliberately hold
    /// *one geometric plane as two handles* (interning would collapse them), and dummy planes
    /// whose handles are never dereferenced. The successor of the retired
    /// `push_surface_unrecorded`, minus the unrecordedness — the truth is still stated, so
    /// nothing point-less enters the arena even from tests.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_plane_unregistered(
        &mut self,
        cache: nacre_geom::Plane,
        points: [[nacre_scalar::Rat; 3]; 3],
    ) -> Handle<Surface> {
        self.push_raw(
            Surface::Plane(cache),
            SurfaceTruth::Plane {
                points: PlanePoints::Known(points),
                motion: None,
            },
        )
    }

    /// Overwrite a plane's truth points — **test-only**, and deliberately incoherence-capable:
    /// the surface's derived name and interning key stay whatever the original points said, so
    /// this door exists for fixtures that need adversarial point widths on an existing surface
    /// (the overflowing-move probe) and must never grow a production caller.
    #[cfg(any(test, feature = "test-util"))]
    pub fn set_plane_points_for_test(
        &mut self,
        h: Handle<Surface>,
        pts: [[nacre_scalar::Rat; 3]; 3],
    ) {
        match &mut self.surface_truths[h.index() as usize] {
            SurfaceTruth::Plane { points, .. } => *points = PlanePoints::Known(pts),
            SurfaceTruth::Cylinder { .. } => panic!("a cylinder has no plane points"),
        }
    }

    /// A new shell whose faces are copies of `src`'s with their outward normals
    /// flipped inward: every loop's winding is reversed and every face
    /// `orientation` is toggled. Pushes fresh [`Face`] cells and a fresh
    /// [`Shell`], but **reuses** `src`'s surfaces, edges, curves, and vertices —
    /// which stay valid handles after the source solid is superseded
    /// (append-only). This is the building block for a cavity (void) shell: the
    /// boundary of a solid whose interior becomes empty space (design §8 M5
    /// containment). The reversed winding keeps each shared edge used with
    /// opposed half-edges (a valid 2-manifold), and the toggled orientation makes
    /// the outward normal point into the void.
    pub fn reversed_shell(&mut self, src: Handle<Shell>) -> Handle<Shell> {
        let src_faces = self.shells.get(src).faces.clone();
        let faces: Vec<Handle<Face>> = src_faces
            .iter()
            .map(|&fh| {
                let face = self.faces.get(fh).clone();
                self.faces.push(Face {
                    surface: face.surface,
                    outer: face.outer.reversed(),
                    inner: face.inner.iter().map(Loop::reversed).collect(),
                    orientation: face.orientation.flipped(),
                })
            })
            .collect();
        self.shells.push(Shell { faces })
    }

    /// The vertex a half-edge starts at: its edge's `vertices[0]` when the use runs
    /// forward, `bounds[1]` when it runs back.
    ///
    /// A traversal accessor, not an analysis — the same kind of thing as
    /// [`Model::reachable`], and the reason it lives here: walking a loop's
    /// corners is the first thing every consumer above does, and it was written
    /// twice (with two different failure policies) before this existed.
    ///
    /// `None` for an edge with no endpoints — the standalone full circle of §4,
    /// which is a legitimate form but has no start. Callers pick their own
    /// policy: a solid's loop edge is always bounded, so `nacre-ops` unwraps with
    /// A vertex's realized coordinate — **the one road to a coordinate from a vertex** (S7),
    /// read from the index-parallel cache [`Model::push_vertex`] fills.
    #[inline]
    pub fn vertex_point(&self, vh: Handle<Vertex>) -> Point3 {
        debug_assert_eq!(
            self.vertex_cache.len(),
            self.vertices.len(),
            "vertex cache out of step with the vertex store — push vertices through Model::push_vertex"
        );
        self.vertex_cache[vh.index() as usize].coord
    }

    /// A vertex's measured coordinate tolerance: `Some` for a discovered vertex (the residual
    /// the arrangement measured when it made the coordinate — `0.0` means exactly zero, kept
    /// exact), `None` for a constructed one (no measurement — a checker applies its own
    /// construction epsilon).
    #[inline]
    pub fn vertex_tol(&self, vh: Handle<Vertex>) -> Option<f64> {
        debug_assert_eq!(self.vertex_cache.len(), self.vertices.len());
        self.vertex_cache[vh.index() as usize].tol
    }

    /// Push a vertex: its definition (the truth) plus the realized coordinate and measured
    /// tolerance (the cache, moved verbatim — S7 moves the field, it does not re-derive; the
    /// coordinate re-derivation question is 3b, deliberately on hold). The one write road.
    ///
    /// ★ There is **no `rebuild_vertex_cache`**: a discovered coordinate is the arrangement's
    /// carefully-made value (measured: 238 of 1,992 differ from a naive re-solve), and a seam
    /// vertex's coordinate is load-bearing (M6). The «discard and regenerate» warrant S8 gave
    /// edges is honestly absent here until then.
    pub fn push_vertex(
        &mut self,
        def: VertexDef,
        coord: Point3,
        tol: Option<f64>,
    ) -> Handle<Vertex> {
        match def {
            VertexDef::ThreePlane([a, b, c]) => debug_assert!(
                a != b && b != c && a != c,
                "a three-plane definition needs three distinct planes"
            ),
            VertexDef::OnSeam([a, b]) => {
                debug_assert!(a != b, "a seam vertex needs two distinct carriers")
            }
            VertexDef::Branch {
                planes: [a, b],
                cylinder,
                ..
            } => debug_assert!(
                a.index() < b.index() && cylinder != a && cylinder != b,
                "a branch definition needs two sorted distinct planes and a distinct cylinder"
            ),
        }
        let h = self.vertices.push(Vertex { def });
        self.vertex_cache.push(PointCache { coord, tol });
        h
    }

    /// An edge's curve — **the one road to a curve from an edge** (S8), read from the
    /// index-parallel cache [`Model::push_edge`] fills.
    #[inline]
    pub fn edge_curve(&self, e: Handle<Edge>) -> &Curve {
        debug_assert_eq!(
            self.edge_cache.len(),
            self.edges.len(),
            "edge cache out of step with the edge store — push edges through Model::push_edge"
        );
        &self.edge_cache[e.index() as usize].curve
    }

    /// Push an edge: canonicalize the carrier pair, derive its curve, fill the cache — the one
    /// write road (S8). `None` when the curve does not derive, which for the line arms means
    /// coincident endpoints (a zero-length edge); the caller maps that to its own reject
    /// (`DegenerateGeometry` / `ZeroLengthEdge`). ★ A rim is `[v, v]` and NOT degenerate — the
    /// circle arm never reads the endpoints (see [`Model::derive_edge_curve`]).
    pub fn push_edge(
        &mut self,
        surfaces: [Handle<Surface>; 2],
        vertices: [Handle<Vertex>; 2],
    ) -> Option<Handle<Edge>> {
        let surfaces = Edge::carrier_pair(surfaces[0], surfaces[1]);
        let curve = self.derive_edge_curve(surfaces, vertices)?;
        let h = self.edges.push(Edge { surfaces, vertices });
        self.edge_cache.push(EdgeCache { curve });
        Some(h)
    }

    /// Discard every edge-curve cache and derive it afresh — the «cache, not truth» warrant
    /// (`docs/truth-and-cache.md`): nothing is lost, because nothing there was truth.
    pub fn rebuild_edge_cache(&mut self) {
        self.edge_cache = self
            .edges
            .iter()
            .map(|(_, e)| EdgeCache {
                curve: self
                    .derive_edge_curve(e.surfaces, e.vertices)
                    .expect("every stored edge derives its curve"),
            })
            .collect();
    }

    /// that invariant while `nacre-props` reports it as unsupported input.
    #[inline]
    pub fn he_start(&self, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = self.edges.get(he.edge).vertices;
        if he.forward { a } else { b }
    }

    /// **An edge's curve, derived from what the model already holds** — the S8 shape of the
    /// truth/cache split: the carriers and the endpoints decide the curve, so the stored one
    /// is a cache that can be discarded and regenerated.
    ///
    /// Dispatch by carrier type:
    /// * **Plane × Plane** (and the self-adjacent cylinder **seam**): the line through the two
    ///   endpoint coordinates — the very expression every producer used to build the stored
    ///   curve, so the derivation is bit-identical, and an endpoint pair that coincides is the
    ///   `None` (a degenerate line — the check lives in this arm only).
    /// * **Plane × Cylinder** (a rim): the circle centred where the cylinder's axis crosses the
    ///   cap plane, with the **cylinder's** frame (`axis direction`, `ref_dir`, `radius`) — the
    ///   same parameters `add_cylinder` builds the stored rims from, so tessellation's `θ`
    ///   parameterization is preserved. The endpoints are not read: a rim is a closed edge
    ///   (`[v, v]`), which is not a degeneracy.
    /// * **Cylinder × Cylinder**: no producer builds one before M6 — `None`, honestly.
    ///
    /// ★ The M3 rim population is axis-perpendicular by construction; a *tilted* plane over a
    /// cylinder would cross in an ellipse, which this arm cannot express (M6). The debug
    /// assertion keeps that boundary visible.
    pub fn derive_edge_curve(
        &self,
        surfaces: [Handle<Surface>; 2],
        vertices: [Handle<Vertex>; 2],
    ) -> Option<Curve> {
        let endpoints_line = || -> Option<Curve> {
            let p0 = self.vertex_point(vertices[0]);
            let p1 = self.vertex_point(vertices[1]);
            Some(Curve::Line(Line::through_points(p0, p1)?))
        };
        match (self.surface(surfaces[0]), self.surface(surfaces[1])) {
            (Surface::Plane(_), Surface::Plane(_)) => endpoints_line(),
            (Surface::Cylinder(_), Surface::Cylinder(_)) if surfaces[0] == surfaces[1] => {
                endpoints_line() // the seam — a parameterization joint, straight along the axis
            }
            (Surface::Plane(p), Surface::Cylinder(c))
            | (Surface::Cylinder(c), Surface::Plane(p)) => {
                let axis = c.axis();
                debug_assert!(
                    {
                        let n = p.normal();
                        let d = axis.direction();
                        n.cross(d).norm_squared() <= 1e-18 * n.norm_squared()
                    },
                    "a tilted plane over a cylinder crosses in an ellipse — M6, no producer yet"
                );
                let center = nacre_geom::intersect::line_plane(&axis, p)?;
                Some(Curve::Circle(Circle::from_center_normal(
                    center,
                    axis.direction(),
                    c.ref_dir(),
                    c.radius(),
                )?))
            }
            (Surface::Cylinder(_), Surface::Cylinder(_)) => None, // two distinct cylinders: M6
        }
    }

    /// The handles reachable from the live solids — the live model (design §2).
    ///
    /// Superseded cells left in the append-only arena are excluded (nothing live
    /// references them). Every step is bounds-guarded, so this is safe even on a
    /// corrupt or partially-built model (a dangling handle simply prunes that
    /// branch; `validate`'s reference-integrity check reports it separately).
    pub fn reachable(&self) -> Reachable {
        let mut r = Reachable::default();
        for &solid_h in &self.live_solids {
            if !in_bounds(solid_h, &self.solids) {
                continue;
            }
            let solid = self.solids.get(solid_h);
            for &shell_h in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
                if !in_bounds(shell_h, &self.shells) || !r.shells.insert(shell_h) {
                    continue;
                }
                for &face_h in &self.shells.get(shell_h).faces {
                    if !in_bounds(face_h, &self.faces) || !r.faces.insert(face_h) {
                        continue;
                    }
                    let face = self.faces.get(face_h);
                    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                        for he in &lp.half_edges {
                            if !in_bounds(he.edge, &self.edges) || !r.edges.insert(he.edge) {
                                continue;
                            }
                            for v in self.edges.get(he.edge).vertices {
                                if in_bounds(v, &self.vertices) {
                                    r.vertices.insert(v);
                                }
                            }
                        }
                    }
                }
            }
        }
        r
    }

    /// Add an axis-aligned box `min`..`max` to this model and return its solid.
    ///
    /// Requires `max[i] > min[i]` on every axis (a degenerate box is a caller
    /// bug → panic); each
    /// face is wound so its plane normal points outward, so all faces are
    /// [`Orientation::Forward`]. Does **not** rebuild adjacency — call
    /// [`Model::rebuild_adjacency`] once after all additions.
    #[cfg(any(test, feature = "test-util"))]
    pub fn add_cuboid(&mut self, min: Point3, max: Point3) -> Handle<Solid> {
        let [x0, y0, z0] = min.as_array();
        let [x1, y1, z1] = max.as_array();
        debug_assert!(
            x1 > x0 && y1 > y0 && z1 > z0,
            "cuboid needs max > min on every axis"
        );

        // 8 corners, numbered by bits (i·x, j·y, k·z).
        let corners = [
            Point3::from_array([x0, y0, z0]), // V0
            Point3::from_array([x1, y0, z0]), // V1
            Point3::from_array([x1, y1, z0]), // V2
            Point3::from_array([x0, y1, z0]), // V3
            Point3::from_array([x0, y0, z1]), // V4
            Point3::from_array([x1, y0, z1]), // V5
            Point3::from_array([x1, y1, z1]), // V6
            Point3::from_array([x0, y1, z1]), // V7
        ];
        // 6 faces: (first 3 loop vertices for the outward plane, half-edges as
        // (edge index, forward)). Loops wound CCW seen from outside → outward
        // normal. See the design plan's winding table (hand-verified).
        type FaceDef = ([usize; 3], [(usize, bool); 4]);
        let faces_def: [FaceDef; 6] = [
            ([0, 3, 2], [(3, false), (2, false), (1, false), (0, false)]), // Bottom −Z
            ([4, 5, 6], [(4, true), (5, true), (6, true), (7, true)]),     // Top +Z
            ([0, 1, 5], [(0, true), (9, true), (4, false), (8, false)]),   // Front −Y
            ([2, 3, 7], [(2, true), (11, true), (6, false), (10, false)]), // Back +Y
            ([0, 4, 7], [(8, true), (7, false), (11, false), (3, true)]),  // Left −X
            ([1, 2, 6], [(1, true), (10, true), (5, false), (9, false)]),  // Right +X
        ];
        // ★ **Surfaces before vertices**, so each corner can name the three it lies on. Separate
        // arenas, so the interleaving does not move any handle; the surfaces' order among
        // themselves is what matters and it is unchanged.
        let surf: [(Handle<Surface>, bool); 6] = core::array::from_fn(|i| {
            let (tri, _) = &faces_def[i];
            // The same three corners in rationals. `from_decimal` because a corner is a value the
            // caller *wrote* — lifting the f64 bit pattern instead would carry its drift in and
            // defeat the whole point (see `Model::surface_name`).
            let rat_corner = |k: usize| -> Option<[nacre_scalar::Rat; 3]> {
                let c = corners[k].as_array();
                Some([
                    nacre_scalar::Rat::from_decimal(c[0])?,
                    nacre_scalar::Rat::from_decimal(c[1])?,
                    nacre_scalar::Rat::from_decimal(c[2])?,
                ])
            };
            self.push_plane(
                Plane::through_points(corners[tri[0]], corners[tri[1]], corners[tri[2]])
                    .expect("non-degenerate box"),
                // ★ The same three corners, exactly. The name is derived from them, so a corner
                // whose decimals are wide enough to overflow the narrow derivation still gets
                // one. A corner outside the decimal window is a caller bug (the radius/height
                // precedent in `add_cylinder`): the truth is not optional any more.
                core::array::from_fn(|k| {
                    rat_corner(tri[k]).expect("cuboid corners inside the decimal window")
                }),
                None,
            )
        });

        let vh: [Handle<Vertex>; 8] = core::array::from_fn(|i| {
            // Which three of the six faces meet at corner `i`. The corner order is
            // `V0..V3` round the bottom then `V4..V7` round the top, so the high bit picks the cap
            // and the position round the ring picks the two walls.
            let m = i % 4;
            let cap = if i < 4 { 0 } else { 1 }; // Bottom −Z / Top +Z
            let along_x = if m == 1 || m == 2 { 5 } else { 4 }; // Right +X / Left −X
            let along_y = if m >= 2 { 3 } else { 2 }; // Back +Y / Front −Y
            self.push_vertex(
                VertexDef::ThreePlane([surf[cap].0, surf[along_y].0, surf[along_x].0]),
                corners[i],
                None,
            )
        });

        // 12 edges as (start, end) vertex indices plus the two faces each edge runs between
        // (indices into `faces_def` — its carriers): bottom ring, top ring, verticals.
        const EDGES: [(usize, usize, [usize; 2]); 12] = [
            (0, 1, [0, 2]), // Bottom·Front
            (1, 2, [0, 5]), // Bottom·Right
            (2, 3, [0, 3]), // Bottom·Back
            (3, 0, [0, 4]), // Bottom·Left
            (4, 5, [1, 2]), // Top·Front
            (5, 6, [1, 5]), // Top·Right
            (6, 7, [1, 3]), // Top·Back
            (7, 4, [1, 4]), // Top·Left
            (0, 4, [2, 4]), // Front·Left
            (1, 5, [2, 5]), // Front·Right
            (2, 6, [3, 5]), // Back·Right
            (3, 7, [3, 4]), // Back·Left
        ];
        // The carrier columns restate what `faces_def` already says (which loops use which
        // edge) — keep the two tables from drifting apart.
        debug_assert!(EDGES.iter().enumerate().all(|(e, &(_, _, fs))| {
            faces_def
                .iter()
                .enumerate()
                .all(|(f, (_, hes))| hes.iter().any(|&(he, _)| he == e) == fs.contains(&f))
        }));
        let eh: [Handle<Edge>; 12] = core::array::from_fn(|i| {
            let (a, b, [fa, fb]) = EDGES[i];
            self.push_edge([surf[fa].0, surf[fb].0], [vh[a], vh[b]])
                .expect("non-degenerate box")
        });

        let fh: [Handle<Face>; 6] = core::array::from_fn(|i| {
            let (_, hes) = &faces_def[i];
            let (surface, flipped) = surf[i];
            let outer = Loop {
                half_edges: hes
                    .iter()
                    .map(|&(e, forward)| HalfEdge {
                        edge: eh[e],
                        forward,
                    })
                    .collect(),
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                // The winding table below is written for a surface whose normal this face uses as-is;
                // a shared surface may point the other way, and then the same outward direction is
                // spelled `Reversed`.
                orientation: if flipped {
                    Orientation::Forward.flipped()
                } else {
                    Orientation::Forward
                },
            })
        });

        let shell = self.shells.push(Shell { faces: fh.to_vec() });
        self.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        })
    }

    /// Add a closed cylinder solid and return it: the axis runs from `base` along
    /// `axis` for `height`, with the given `radius`. The b-rep is the shared seam form
    /// (`Model::cylinder_solid`); this entry's own work is **deriving the exact truth from a
    /// caller's f64 statement**.
    ///
    /// `radius`/`height` must be positive, `axis` nonzero, and every stated coordinate inside
    /// the decimal window (caller bug → panic). Does **not** rebuild adjacency — call
    /// [`Model::rebuild_adjacency`] once after all additions.
    ///
    /// ★ **The panics are why this is a test convenience.** A statement whose *computed*
    /// coordinates leave the decimal window is not a caller bug when the caller is an
    /// application — it is an input. The production road states a cylinder through
    /// [`Model::add_cylinder_exact`], whose frame makes those coordinates rational by
    /// construction and which names every refusal instead of panicking.
    #[cfg(any(test, feature = "test-util"))]
    pub fn add_cylinder(
        &mut self,
        base: Point3,
        axis: Vector3,
        radius: f64,
        height: f64,
    ) -> Handle<Solid> {
        debug_assert!(
            radius > 0.0 && height > 0.0,
            "cylinder needs positive radius and height"
        );
        let d = axis.normalize().expect("cylinder axis must be nonzero");
        let u = d
            .any_perpendicular()
            .expect("a unit axis has a perpendicular");
        let c0 = base;
        let c1 = base + d * height;
        let p_bot = c0 + u * radius; // seam point on the bottom rim (angle 0)
        let p_top = c1 + u * radius; // seam point on the top rim

        // The exact truth (M6-0): a direct lift of the caller's statement — origin and the raw,
        // unnormalized axis (normalizing would destroy the exact form; the `normal_def`
        // precedent). `ref_dir` replicates `any_perpendicular`'s own rule in rationals: cross
        // the axis with the basis axis of its smallest |component| (ties X→Y→Z).
        //
        // ★ **The basis choice reads `d` — the very components `any_perpendicular` reads — so
        // agreement is structural, not order-theoretic.** The first spelling compared the *raw*
        // components and argued "one positive scale preserves |·| order"; that is a real-number
        // argument, and f64 division rounds: a strict `|x| > |y|` can collapse to equality in
        // `d`, flipping which side of the `<=` tie-break each rule lands on (measured — axis
        // `[0.34, 0.33999999999999997, 1.0]`: raw picks Y, `d` picks X, seam ~90° apart, the
        // validate net fires on a healthy model; pinned in `tests/cylinder_truth.rs`).
        //
        // The *cross* still uses the raw exact components — `ê_k × raw` is positively parallel
        // to `ê_k × d` whichever values chose `k` — so the seam direction stays exact.
        // A statement outside the decimal window is a caller bug → panic (the radius/height
        // precedent above).
        let lift = |x: f64| -> nacre_scalar::Rat {
            nacre_scalar::Rat::from_decimal(x)
                .expect("cylinder statement inside the decimal window")
        };
        let def = {
            let zero = nacre_scalar::Rat::from_int(0);
            // Lift then negate (not lift the negated f64): `-0.0` has no decimal of its own.
            let neg = |x: f64| {
                zero.checked_sub(lift(x))
                    .expect("negating a lifted decimal cannot overflow")
            };
            let a = axis.as_array();
            let ax = d.as_array().map(f64::abs);
            // ê_k × axis, k = the smallest-|component| basis axis of `d` — the identical
            // comparison chain `any_perpendicular` runs on the identical inputs.
            let ref_dir = if ax[0] <= ax[1] && ax[0] <= ax[2] {
                [zero, neg(a[2]), lift(a[1])]
            } else if ax[1] <= ax[2] {
                [lift(a[2]), zero, neg(a[0])]
            } else {
                [neg(a[1]), lift(a[0]), zero]
            };
            // ★ **Unreachable, and now provably so.** `CylinderDef::new` refuses a zero axis,
            // a non-positive radius, and a `ref_dir` parallel to the axis. The first two are
            // the caller's debug_asserts above; the third cannot happen here: `ref_dir` is
            // `ê_k × a` for the basis axis `k` of *smallest* |component|, and `a ∥ ê_k` would
            // need `|a_k|` to be both the largest and the smallest component — true only for
            // the zero axis. (Positive scaling preserves parallelism, so picking `k` from the
            // normalized `d` while crossing the raw `a` does not disturb the argument.)
            // Before the parallelism test became total, this `expect` also fired on statements
            // it had no business rejecting — a long-decimal axis component whose square left
            // `i128`.
            CylinderDef::new(
                base.as_array().map(lift),
                a.map(lift),
                ref_dir,
                lift(radius),
            )
            .expect("non-degenerate cylinder")
        };

        // A cap plane's three exact points: the decimal truth of the realized center and two
        // rim-direction offsets the construction already computed (the `add_cuboid` precedent —
        // the producer's own f64 is its statement). For an axis whose normalization is exact
        // (`ẑ`, a Pythagorean triple) these are exact by construction; for an irrational axis
        // they are the truth of what was *built*, which is all a `Constructed` surface ever
        // claims. A cap outside the decimal window is a caller bug (the radius/height
        // precedent above): the truth is not optional any more.
        let w = d.cross(u);
        let cap_points = |c: Point3| -> [[nacre_scalar::Rat; 3]; 3] {
            let lift = |p: Point3| -> [nacre_scalar::Rat; 3] {
                p.as_array().map(|x| {
                    nacre_scalar::Rat::from_decimal(x)
                        .expect("cylinder caps inside the decimal window")
                })
            };
            [lift(c), lift(c + u * radius), lift(c + w * radius)]
        };
        self.cylinder_solid(CylinderParts {
            def,
            lateral: Cylinder::from_axis(c0, d, u, radius).expect("non-degenerate cylinder"),
            caps: [
                (
                    Plane::from_point_normal(c0, -d).expect("nonzero axis"),
                    cap_points(c0),
                ),
                (
                    Plane::from_point_normal(c1, d).expect("nonzero axis"),
                    cap_points(c1),
                ),
            ],
            seam_pts: [p_bot, p_top],
            motion: None,
        })
        .expect("a rim derives its circle from the cap and the cylinder")
        .0
    }

    /// **State a cylinder exactly** — the production road, and the one with no panics in it.
    ///
    /// `axis` and `ref_dir` must be **unit and perpendicular** (checked, exactly). That single
    /// precondition is what makes everything below exact rather than realized-then-lifted: the
    /// far centre `base + axis·height`, both seam points `c + ref_dir·radius`, and each cap
    /// plane's third point `c + (axis × ref_dir)·radius` are rational products of rational
    /// inputs. `Model::add_cylinder` (the test entry) cannot do this — it normalizes an f64 axis, so its caps
    /// and seams have to be lifted back out of computed floats, and a statement whose computed
    /// coordinates leave the decimal window panics there.
    ///
    /// ★ The precondition is not a restriction on *what can be built*: an application states a
    /// cylinder on a sketch frame, and a frame either has an exact orthonormal basis or is
    /// carried by `motion` — in which case the cylinder is stated in the frame's own coordinates,
    /// where the basis is `{0, ±1}`. So the tilted case is not the irrational case.
    ///
    /// Returns the solid beside its three faces (lateral, bottom cap, top cap). Does **not**
    /// rebuild adjacency.
    pub fn add_cylinder_exact(
        &mut self,
        base: [Rat; 3],
        axis: [Rat; 3],
        ref_dir: [Rat; 3],
        radius: Rat,
        height: Rat,
        motion: Option<Handle<MotionNode>>,
    ) -> Result<(Handle<Solid>, [Handle<Face>; 3]), CylinderError> {
        let zero = Rat::from_int(0);
        let one = Rat::from_int(1);
        if radius <= zero {
            return Err(CylinderError::NonPositiveRadius);
        }
        if height <= zero {
            return Err(CylinderError::NonPositiveHeight);
        }
        let dot = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<Rat> {
            let mut acc = zero;
            for k in 0..3 {
                acc = acc.checked_add(a[k].checked_mul(b[k])?)?;
            }
            Some(acc)
        };
        let (aa, rr, ar) = (
            dot(&axis, &axis).ok_or(CylinderError::Overflow)?,
            dot(&ref_dir, &ref_dir).ok_or(CylinderError::Overflow)?,
            dot(&axis, &ref_dir).ok_or(CylinderError::Overflow)?,
        );
        if aa != one || rr != one || ar != zero {
            return Err(CylinderError::FrameNotOrthonormal);
        }
        // The third direction of the frame, so a cap plane gets a second rim point rather than a
        // second statement of the same one — three points on a circle name its plane.
        let cross = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<[Rat; 3]> {
            let term =
                |i: usize, j: usize| a[i].checked_mul(b[j])?.checked_sub(a[j].checked_mul(b[i])?);
            Some([term(1, 2)?, term(2, 0)?, term(0, 1)?])
        };
        let w = cross(&axis, &ref_dir).ok_or(CylinderError::Overflow)?;
        // `p + dir·s`, exactly — the only arithmetic this entry does, and the reason the
        // orthonormal precondition is worth checking.
        let step = |p: &[Rat; 3], dir: &[Rat; 3], s: Rat| -> Option<[Rat; 3]> {
            let mut out = *p;
            for k in 0..3 {
                out[k] = p[k].checked_add(dir[k].checked_mul(s)?)?;
            }
            Some(out)
        };
        let over = CylinderError::Overflow;
        let c0 = base;
        let c1 = step(&c0, &axis, height).ok_or(over)?;
        let cap_points = |c: &[Rat; 3]| -> Option<[[Rat; 3]; 3]> {
            Some([*c, step(c, &ref_dir, radius)?, step(c, &w, radius)?])
        };
        let (bottom_points, top_points) =
            (cap_points(&c0).ok_or(over)?, cap_points(&c1).ok_or(over)?);
        let (p_bot, p_top) = (bottom_points[1], top_points[1]);
        let def = CylinderDef::new(base, axis, ref_dir, radius).ok_or(CylinderError::Degenerate)?;

        // The caches are the realization of exactly these statements — nothing here is measured
        // or re-derived, so a cache cannot disagree with the truth beside it.
        let point = |p: [Rat; 3]| Point3::from_array(p.map(Rat::to_f64));
        let vector = |p: [Rat; 3]| Vector3::from_array(p.map(Rat::to_f64));
        let (d, u) = (vector(axis), vector(ref_dir));
        let deg = CylinderError::Degenerate;
        self.cylinder_solid(CylinderParts {
            def,
            lateral: Cylinder::from_axis(point(c0), d, u, radius.to_f64()).ok_or(deg)?,
            caps: [
                (
                    Plane::from_point_normal(point(c0), -d).ok_or(deg)?,
                    bottom_points,
                ),
                (
                    Plane::from_point_normal(point(c1), d).ok_or(deg)?,
                    top_points,
                ),
            ],
            seam_pts: [point(p_bot), point(p_top)],
            motion,
        })
        .ok_or(deg)
    }

    /// Assemble a cylinder's b-rep from parts already derived — **the one spelling** of the
    /// V2/E3/F3 seam form, shared by both constructors.
    ///
    /// Two seam vertices, two full-circle rim edges (`bounds: Some([seam, seam])`,
    /// start == end), one straight seam edge, a cylindrical lateral face whose loop uses the
    /// seam edge twice (opposite orientation), and two planar caps (each a single-half-edge rim
    /// loop). This forms a valid CW-complex (Euler χ = 2) that
    /// [`validate`](../nacre_validate/fn.validate.html) accepts — a closed periodic surface
    /// needs a seam vertex, so the rims are `Some([v, v])`, not `bounds: None` (that form is for
    /// a standalone full circle; §4).
    ///
    /// Returns the solid beside its three faces in push order (lateral, bottom cap, top cap).
    /// `None` if an edge cannot derive its curve from the carriers it states — each caller says
    /// what that means for it. Does **not** rebuild adjacency.
    fn cylinder_solid(
        &mut self,
        parts: CylinderParts,
    ) -> Option<(Handle<Solid>, [Handle<Face>; 3])> {
        let CylinderParts {
            def,
            lateral: lateral_cache,
            caps: [(bottom_plane, bottom_points), (top_plane, top_points)],
            seam_pts: [p_bot, p_top],
            motion,
        } = parts;

        // ★ **Surfaces before edges** (S8): an edge states its two carriers, so the lateral
        // cylinder and both cap planes must exist first. Separate arenas — the interleaving
        // moves no handle; the surfaces' order among themselves (lateral → bottom cap → top
        // cap) is what matters and it is unchanged (the `add_cuboid` precedent).
        let lateral_surface = self.push_cylinder(lateral_cache, def, motion);
        let (bottom_cap_surface, _) = self.push_plane(bottom_plane, bottom_points, motion);
        let (top_cap_surface, _) = self.push_plane(top_plane, top_points, motion);

        // A seam vertex lies on two surfaces only — the rim circle's `θ = 0` point. `OnSeam`
        // states exactly that (S7), and since M6-0 the designation is complete: the cylinder's
        // truth carries `ref_dir`, so the pair means "rim ∩ the `+ref_dir` ray" — one point
        // (see `VertexDef::OnSeam`; regenerating the cached coordinate stays deferred with 3b).
        let v_bot = self.push_vertex(
            VertexDef::OnSeam([lateral_surface, bottom_cap_surface]),
            p_bot,
            None,
        );
        let v_top = self.push_vertex(
            VertexDef::OnSeam([lateral_surface, top_cap_surface]),
            p_top,
            None,
        );

        // Rims are full circles seamed at their vertex (start == end); the seam is
        // a straight edge joining the two rim seam points.
        let bottom = self.push_edge([lateral_surface, bottom_cap_surface], [v_bot, v_bot])?;
        let top = self.push_edge([lateral_surface, top_cap_surface], [v_top, v_top])?;
        // Self-adjacent: a seam is a parameterization joint of ONE surface, not an
        // intersection of two — the confirmed spelling (M6-0; see `Edge::surfaces`), guarded
        // by validate's "self-adjacent ⇔ cylinder" carrier rule.
        let seam = self.push_edge([lateral_surface, lateral_surface], [v_bot, v_top])?;

        // Lateral cylindrical face: one loop wrapping the seam twice (opposite).
        let lateral = {
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: bottom,
                        forward: true,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: true,
                    },
                    HalfEdge {
                        edge: top,
                        forward: false,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: false,
                    },
                ],
            };
            self.faces.push(Face {
                surface: lateral_surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Bottom cap: outward normal −d, the bottom rim reversed.
        let bottom_cap = {
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: bottom,
                    forward: false,
                }],
            };
            self.faces.push(Face {
                surface: bottom_cap_surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Top cap: outward normal +d, the top rim forward.
        let top_cap = {
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: top,
                    forward: true,
                }],
            };
            self.faces.push(Face {
                surface: top_cap_surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };

        let shell = self.shells.push(Shell {
            faces: vec![lateral, bottom_cap, top_cap],
        });
        let solid = self.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        Some((solid, [lateral, bottom_cap, top_cap]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;
    use proptest::prelude::*;

    /// ★★★ **The defect the rational coefficients exist to remove, pinned from both sides.**
    ///
    /// Two boxes meet on the plane `x = 3` with faces of different size. `Plane` stores an
    /// un-normalized normal whose length follows that size, so `d = −raw·origin` is a differently
    /// rounded product on each side and the two coefficient vectors are **not exactly
    /// proportional** — the f64 test says "different planes" about one plane. Measured across the
    /// census, 18 pairs are merged only because a second test looks at the faces' coordinates
    /// instead (`docs/dev-log.md`).
    ///
    /// The rational coefficients are built from the corners the caller wrote and canonicalized, so
    /// they have no scale to disagree about and come out **equal**.
    ///
    /// Both halves are load-bearing. If the first assertion ever fails the f64 defect was fixed
    /// somewhere else and this test should be re-read, not deleted; if the second fails the
    /// rational path stopped reaching these surfaces.
    #[test]
    fn two_faces_of_one_plane_disagree_in_f64_and_agree_in_the_rationals() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 2.2, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([3.0, 0.0, 0.0]),
            Point3::from_array([5.0, 13.2, 1.0]),
        );
        // Each solid's face on x = 3: a's outward +X, b's outward −X.
        let face_on_x3 = |s: Handle<Solid>, want_x: f64| -> Handle<Surface> {
            let shell = m.solids.get(s).outer;
            *m.shells
                .get(shell)
                .faces
                .iter()
                .map(|&fh| &m.faces.get(fh).surface)
                .find(|&&sh| match m.surfaces.get(sh) {
                    Surface::Plane(p) => {
                        let [a, b, c, d] = p.coefficients();
                        b == 0.0 && c == 0.0 && a != 0.0 && (-d / a - want_x).abs() < 1e-12
                    }
                    Surface::Cylinder(_) => false,
                })
                .expect("a face on x = 3")
        };
        let (sa, sb) = (face_on_x3(a, 3.0), face_on_x3(b, 3.0));

        // ★ **They are one handle now** — that is what the rational coefficients bought.
        assert_eq!(sa, sb, "one plane, one surface");
        assert!(m.surface_name.contains_key(&sa), "and it is recorded");

        // The f64 defect that made this necessary, shown on the planes themselves rather than
        // through the model, since the model no longer holds two of them. `Plane` keeps an
        // un-normalized normal whose length follows the face's size, so `d = −raw·origin` is a
        // differently rounded product on each side and the two vectors are not exactly
        // proportional — the f64 test says "different planes" about one plane.
        let wall = |dy: f64| {
            Plane::through_points(
                Point3::from_array([3.0, 0.0, 0.0]),
                Point3::from_array([3.0, dy, 0.0]),
                Point3::from_array([3.0, 0.0, 1.0]),
            )
            .expect("non-degenerate")
        };
        let (pa, pb) = (wall(2.2), wall(13.2));
        assert!(
            !nacre_geom::intersect::planes_coplanar(&pa, &pb),
            "f64 coefficients of one plane at two face sizes: {:?} vs {:?}",
            pa.coefficients(),
            pb.coefficients()
        );
        let rat = |dy: f64| {
            let r = |x: f64| nacre_scalar::Rat::from_decimal(x).expect("decimal");
            nacre_scalar::plane_through_points(
                [r(3.0), r(0.0), r(0.0)],
                [r(3.0), r(dy), r(0.0)],
                [r(3.0), r(0.0), r(1.0)],
            )
        };
        assert_eq!(
            rat(2.2),
            rat(13.2),
            "the rationals have no scale to disagree about"
        );
        assert!(rat(2.2).is_some());
    }

    fn build(min: [f64; 3], max: [f64; 3]) -> Model {
        let mut m = Model::new();
        m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
        m.rebuild_adjacency();
        m
    }

    fn he_start(m: &Model, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = m.edges.get(he.edge).vertices;
        if he.forward { a } else { b }
    }
    fn he_end(m: &Model, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = m.edges.get(he.edge).vertices;
        if he.forward { b } else { a }
    }
    fn face_plane_normal(m: &Model, f: &Face) -> Vector3 {
        match m.surfaces.get(f.surface) {
            Surface::Plane(p) => p.normal(),
            // Planar-only helper: callers filter to plane faces (caps), never cylinders.
            Surface::Cylinder(_) => unreachable!("face_plane_normal called on a curved face"),
        }
    }
    fn face_centroid(m: &Model, f: &Face) -> Point3 {
        let pts: Vec<Point3> = f
            .outer
            .half_edges
            .iter()
            .map(|he| m.vertex_point(he_start(m, *he)))
            .collect();
        Point3::centroid(&pts).unwrap()
    }

    // --- golden ---

    #[test]
    fn cuboid_counts() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.edges.len(), 12);
        assert_eq!(m.faces.len(), 6);
        assert_eq!(m.shells.len(), 1);
        assert_eq!(m.solids.len(), 1);
        // Still 6 — but differently composed since S9: the origin box's bottom/left/front
        // faces intern onto the three seeded world planes (same name, same handle), so the
        // arena holds 3 seeds + 3 fresh (top/back/right). Seeding adds nothing here precisely
        // because the seeds are these planes.
        assert_eq!(m.surfaces.len(), 6);
        assert_eq!(
            m.edge_cache.len(),
            m.edges.len(),
            "the curve cache stays index-parallel"
        );
    }

    #[test]
    fn cuboid_corner_points() {
        let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        let expected = vec![
            [-2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 4.0, 0.0],
            [-2.0, 4.0, 0.0],
            [-2.0, 1.0, 10.0],
            [3.0, 1.0, 10.0],
            [3.0, 4.0, 10.0],
            [-2.0, 4.0, 10.0],
        ];
        let got: Vec<[f64; 3]> = m
            .vertices
            .iter()
            .map(|(vh, _)| m.vertex_point(vh).as_array())
            .collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn face_normals_point_outward() {
        let m = build([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let center = Point3::origin().lerp(Point3::from_array([2.0, 3.0, 4.0]), 0.5);
        for (_, f) in m.faces.iter() {
            let outward = face_plane_normal(&m, f).dot(face_centroid(&m, f) - center);
            assert!(outward > 0.0);
        }
    }

    #[test]
    fn every_edge_used_twice_opposite() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.adj.edge_uses.len(), 12);
        for uses in m.adj.edge_uses.values() {
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1);
        }
    }

    #[test]
    fn every_vertex_incident_to_three_edges() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.adj.vertex_edges.len(), 8);
        for edges in m.adj.vertex_edges.values() {
            assert_eq!(edges.len(), 3);
        }
    }

    #[test]
    fn outer_loops_are_closed() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        for (_, f) in m.faces.iter() {
            let hes = &f.outer.half_edges;
            assert_eq!(hes.len(), 4);
            for i in 0..hes.len() {
                assert_eq!(he_end(&m, hes[i]), he_start(&m, hes[(i + 1) % hes.len()]));
            }
        }
    }

    /// ★★ S8, the «discard and regenerate» warrant: the edge-curve cache rebuilt from the
    /// carriers and endpoints is bit-identical to the one `push_edge` filled eagerly — proof
    /// that nothing in it was truth.
    #[test]
    fn edge_cache_discard_and_regenerate_bit_identical() {
        let mut m = Model::new();
        m.add_cuboid(
            Point3::from_array([-2.0, 1.0, 0.0]),
            Point3::from_array([3.0, 4.0, 10.0]),
        );
        m.add_cylinder(
            Point3::from_array([8.0, 0.0, 0.0]),
            Vector3::from_array([0.3, -0.4, 1.0]),
            1.25,
            2.5,
        );
        let snapshot = |m: &Model| -> Vec<Curve> {
            m.edges
                .iter()
                .map(|(eh, _)| m.edge_curve(eh).clone())
                .collect()
        };
        let before = snapshot(&m);
        m.rebuild_edge_cache();
        assert_eq!(
            before,
            snapshot(&m),
            "regeneration must reproduce the cache bit for bit"
        );
    }

    #[test]
    fn edge_endpoints_lie_on_their_curve() {
        let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        for (eh, e) in m.edges.iter() {
            let curve = m.edge_curve(eh);
            let [a, b] = e.vertices;
            assert!(curve.contains(m.vertex_point(a), 1e-9));
            assert!(curve.contains(m.vertex_point(b), 1e-9));
        }
    }

    #[test]
    fn euler_poincare_holds() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (v, e, f) = (
            m.vertices.len() as i64,
            m.edges.len() as i64,
            m.faces.len() as i64,
        );
        // V − E + F = 2(S − G) + L_i, with S=1, G=0, L_i=0. Formal validate:
        // nacre-validate (next unit).
        assert_eq!(v - e + f, 2);
    }

    // --- cylinder (seam b-rep) --- (validate-clean lives in nacre-validate)

    fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
        let mut m = Model::new();
        m.add_cylinder(Point3::from_array(base), Vector3::from_array(axis), r, h);
        m.rebuild_adjacency();
        m
    }

    #[test]
    fn cylinder_counts_and_euler() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        assert_eq!(m.vertices.len(), 2);
        assert_eq!(m.edges.len(), 3);
        assert_eq!(m.faces.len(), 3);
        assert_eq!(m.shells.len(), 1);
        assert_eq!(m.solids.len(), 1);
        // Euler χ = V − E + F = 2 (one shell, genus 0, no inner loops).
        assert_eq!(
            m.vertices.len() as i64 - m.edges.len() as i64 + m.faces.len() as i64,
            2
        );
    }

    #[test]
    fn cylinder_geometry() {
        // +Z axis, r=2, h=5. Seam direction is X (least-aligned axis of +Z).
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let pts: Vec<[f64; 3]> = m
            .vertices
            .iter()
            .map(|(vh, _)| m.vertex_point(vh).as_array())
            .collect();
        // Seam direction for +Z is any_perpendicular([0,0,1]) = X×Z = [0,-1,0], so
        // the seam vertices sit at radius 2 along −Y, at z=0 and z=5.
        assert_eq!(pts, vec![[0.0, -2.0, 0.0], [0.0, -2.0, 5.0]]);
        // Two rim circles carry a Circle; the straight seam a Line.
        let mut circles = 0;
        for (eh, _) in m.edges.iter() {
            if let Curve::Circle(c) = m.edge_curve(eh) {
                assert_eq!(c.radius(), 2.0);
                circles += 1;
            }
        }
        assert_eq!(circles, 2);
    }

    #[test]
    fn cylinder_caps_point_outward() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        // Planar caps only (the lateral cylindrical face has no single normal).
        let mut caps = 0;
        for (_, f) in m.faces.iter() {
            if matches!(m.surfaces.get(f.surface), Surface::Plane(_)) {
                let n = face_plane_normal(&m, f).as_array();
                // Bottom cap → −Z, top cap → +Z (outward along the axis).
                assert!(n == [0.0, 0.0, -1.0] || n == [0.0, 0.0, 1.0]);
                caps += 1;
            }
        }
        assert_eq!(caps, 2);
    }

    /// ★★ S7: **a seam vertex states its two carriers** — `OnSeam([lateral, its own cap])`.
    /// The pair pins the rim circle; the coordinate pins the point until M6's `ref_dir` truth
    /// (see `VertexDef::OnSeam`). The lateral surface must be the cylinder, the other its cap.
    #[test]
    fn a_seam_vertex_states_its_rim_carriers() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let defs: Vec<_> = m.vertices.iter().map(|(_, v)| v.def).collect();
        assert_eq!(
            defs.len(),
            2,
            "a cylinder has exactly its two seam vertices"
        );
        for (i, d) in defs.iter().enumerate() {
            let VertexDef::OnSeam([a, b]) = d else {
                panic!("a seam vertex carries OnSeam, got {d:?}")
            };
            assert!(
                matches!(m.surface(*a), Surface::Cylinder(_)),
                "first carrier is the lateral cylinder"
            );
            let Surface::Plane(p) = m.surface(*b) else {
                panic!("second carrier is the cap plane")
            };
            // Bottom vertex names the bottom cap (through z = 0), top the top cap (z = 5).
            let z = if i == 0 { 0.0 } else { 5.0 };
            assert_eq!(p.distance(Point3::from_array([0.0, 0.0, z])), 0.0);
        }
    }

    #[test]
    fn cylinder_seam_edge_is_self_adjacent() {
        // The novel topology: the seam edge is used twice by the SAME lateral
        // face with opposite orientation (a valid non-manifold-looking but
        // manifold seam). Every edge is still used exactly twice, opposite.
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let mut seam_uses = None;
        for (eh, _) in m.edges.iter() {
            let uses = &m.adj.edge_uses[&eh];
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1); // opposite orientation
            // The seam's discriminator IS the new invariant: self-adjacent carriers (S8).
            if m.edges.get(eh).surfaces[0] == m.edges.get(eh).surfaces[1] {
                seam_uses = Some(uses.clone());
            }
        }
        let uses = seam_uses.expect("a seam line edge exists");
        assert_eq!(uses[0].0, uses[1].0); // both uses are the same (lateral) face
    }

    // --- proptest ---

    fn box_strategy() -> impl Strategy<Value = ([f64; 3], [f64; 3])> {
        (
            prop::array::uniform3(-1e3f64..1e3),
            prop::array::uniform3(1e-2f64..1e3),
        )
            .prop_map(|(min, ext)| {
                let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
                (min, max)
            })
    }

    proptest! {
        #[test]
        fn prop_cuboid_structural((min, max) in box_strategy()) {
            let m = build(min, max);
            prop_assert_eq!(m.vertices.len(), 8);
            prop_assert_eq!(m.edges.len(), 12);
            prop_assert_eq!(m.faces.len(), 6);
            prop_assert_eq!(m.adj.edge_uses.len(), 12);
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
            prop_assert_eq!(m.adj.vertex_edges.len(), 8);
            for edges in m.adj.vertex_edges.values() {
                prop_assert_eq!(edges.len(), 3);
            }
        }

        #[test]
        fn prop_cuboid_outward_normals((min, max) in box_strategy()) {
            let m = build(min, max);
            let center = Point3::from_array(min).lerp(Point3::from_array(max), 0.5);
            for (_, f) in m.faces.iter() {
                prop_assert!(face_plane_normal(&m, f).dot(face_centroid(&m, f) - center) > 0.0);
            }
        }

        #[test]
        fn prop_cuboid_endpoints_on_curves((min, max) in box_strategy()) {
            let m = build(min, max);
            let scale = 1e-6 * (max.iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
            for (eh, e) in m.edges.iter() {
                let curve = m.edge_curve(eh);
                let [a, b] = e.vertices;
                prop_assert!(curve.contains(m.vertex_point(a), scale));
                prop_assert!(curve.contains(m.vertex_point(b), scale));
            }
        }

        #[test]
        fn prop_cylinder_structural(
            base in prop::array::uniform3(-1e3f64..1e3),
            axis in prop::array::uniform3(-1.0f64..1.0),
            r in 0.5f64..10.0,
            h in 0.1f64..10.0,
        ) {
            let axis = Vector3::from_array(axis);
            prop_assume!(axis.norm() > 0.1); // skip near-zero axes
            let mut m = Model::new();
            m.add_cylinder(Point3::from_array(base), axis, r, h);
            m.rebuild_adjacency();
            prop_assert_eq!(m.vertices.len(), 2);
            prop_assert_eq!(m.edges.len(), 3);
            prop_assert_eq!(m.faces.len(), 3);
            // Every edge used exactly twice, opposite orientation (seam included).
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
        }

        /// The same structure over the population that **actually stressed the exact
        /// arithmetic**: an axis a hair off `ẑ`, whose tiny components carry a full f64's worth
        /// of decimal digits and so lift to rationals with ~10²⁰ denominators. Squaring one of
        /// those leaves `i128`, which is what the parallelism test used to refuse — and
        /// `add_cylinder` read that refusal as "degenerate cylinder" and panicked.
        ///
        /// ★ The sibling above cannot stand in for this: its uniform axis reaches this family
        /// about **0.01%** of the time (measured), which is why the defect sat green for a
        /// milestone and then surfaced from one unlucky seed. Here it is ~80%.
        #[test]
        fn prop_cylinder_near_axis_aligned_is_built_not_refused(
            u in -1.0f64..1.0,
            v in -1.0f64..1.0,
            k in 1i32..9,
            j in 1i32..9,
            r in 0.5f64..10.0,
            h in 0.1f64..10.0,
        ) {
            let axis = Vector3::from_array([u * 10f64.powi(-k), v * 10f64.powi(-j), 1.0]);
            let mut m = Model::new();
            m.add_cylinder(Point3::origin(), axis, r, h);
            m.rebuild_adjacency();
            prop_assert_eq!(m.vertices.len(), 2);
            prop_assert_eq!(m.edges.len(), 3);
            prop_assert_eq!(m.faces.len(), 3);
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
        }
    }

    // --- reversed_shell (M5 containment building block) ---

    #[test]
    fn orientation_flip_is_involution() {
        assert_eq!(Orientation::Forward.flipped(), Orientation::Reversed);
        assert_eq!(Orientation::Reversed.flipped(), Orientation::Forward);
        assert_eq!(
            Orientation::Forward.flipped().flipped(),
            Orientation::Forward
        );
    }

    #[test]
    fn reversed_shell_toggles_orientation_and_reverses_loops() {
        let mut m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let outer = m.solids.get(m.live_solids[0]).outer;
        let f0 = m.shells.get(outer).faces[0];
        let orig = m.faces.get(f0).clone();

        let rev_shell = m.reversed_shell(outer);
        // Fresh cells (not reused faces), fresh shell.
        assert_ne!(rev_shell, outer);
        let rf0 = m.shells.get(rev_shell).faces[0];
        assert_ne!(rf0, f0);
        let rev = m.faces.get(rf0);

        assert_eq!(rev.surface, orig.surface); // surface reused
        assert_eq!(rev.orientation, orig.orientation.flipped());
        let n = orig.outer.half_edges.len();
        assert_eq!(rev.outer.half_edges.len(), n);
        // Reversed winding: he[i] mirrors orig[n-1-i] with the edge reused and
        // the traversal direction flipped.
        for i in 0..n {
            let o = orig.outer.half_edges[n - 1 - i];
            let r = rev.outer.half_edges[i];
            assert_eq!(r.edge, o.edge);
            assert_ne!(r.forward, o.forward);
        }
    }

    #[test]
    fn reversed_shell_is_a_valid_manifold() {
        // Reversing every face's winding preserves the b-rep manifold: each edge
        // is still used by exactly two faces with opposed half-edges. The
        // reversed shell reuses the cube's edges, but the source solid is
        // superseded, so `Adjacency` (reachable-scoped) counts only the reversed
        // faces — no 4-use false positive.
        let mut m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let outer = m.solids.get(m.live_solids[0]).outer;
        let rev = m.reversed_shell(outer);
        let s = m.push_solid(Solid {
            outer: rev,
            cavities: vec![],
        });
        m.live_solids.retain(|&h| h == s); // supersede the original cube
        m.rebuild_adjacency();
        assert_eq!(m.adj.edge_uses.len(), 12);
        for uses in m.adj.edge_uses.values() {
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1); // opposite forward
        }
    }

    /// ★★★★★ **A plane is named by its points, and by nothing else.**
    ///
    /// The three assertions are the three ways this can go wrong. A name has to come out for an
    /// ordinary plane; it has to be the plane the points are **on**, not some other; and it has to
    /// come out even where the derivation's narrow route cannot reach — that last one is what a
    /// producer used to lose a name to, and losing a name closes the exact road for everything
    /// built on that plane.
    ///
    /// ★ There is no "wrong name" case left to test. Coefficients are no longer something a
    /// producer can hand in, so a plane cannot be stated twice — the state the old agreement filter
    /// watched for is now unspellable.
    #[test]
    fn a_plane_is_named_by_its_points() {
        use nacre_scalar::Rat;
        let r = Rat::from_int;
        let pl = |z: f64| {
            nacre_geom::Plane::from_point_normal(
                Point3::from_array([0.0, 0.0, z]),
                Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap()
        };

        let mut m = Model::new();
        let pts = [[r(0), r(0), r(0)], [r(1), r(0), r(0)], [r(0), r(1), r(0)]];
        let (ok, _) = m.push_plane(pl(0.0), pts, None);
        assert_eq!(
            m.surface_name.get(&ok),
            Some(&nacre_scalar::PlaneName::Narrow([r(0), r(0), r(1), r(0)])),
            "the plane z = 0 was not named, or was named as something else"
        );

        // ★★★★ **Wide enough that the narrow derivation gives up.** Decimal arithmetic reduces to
        // denominators that are powers of two on one coordinate and powers of five on another, and
        // a triple of *coprime* denominators needs their product to state its plane — which
        // `(b − a) × (c − a)` then needs squared.
        //
        // ★ An earlier spelling here used points chosen to overflow the old **agreement check**
        // (`c · p`) and asserted the derivation gave up on them too. It does not: those are
        // different products, and the derivation went through. Two propositions, one fixture.
        let q = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
        let wide_pts = [
            [q(1, 1 << 53), r(0), r(0)],
            [r(0), q(1, 5i128.pow(23)), r(0)],
            [r(0), r(0), q(1, (1 << 40) * 5i128.pow(11))],
        ];
        assert_eq!(
            nacre_scalar::plane_through_points(wide_pts[0], wide_pts[1], wide_pts[2]),
            None,
            "the narrow route was expected to overflow here — the case has stopped being the case"
        );
        let (wide, _) = m.push_plane(pl(2.0), wide_pts, None);
        // ★ The stored value carries the proposition — `Narrow` says "fits i128" directly.
        // (This used to compare the global counter before/after, which races against other
        // tests pushing wide planes in parallel; the value cannot.)
        let name = *m
            .surface_name
            .get(&wide)
            .expect("a plane the narrow route cannot reach went unnamed")
            .narrow()
            .expect("this fixture's canonical answer is small — it must be stored Narrow");
        // ★★ And it names *these* points' plane — exact rationals, no tolerance.
        for p in &wide_pts {
            let mut acc = name[3];
            for k in 0..3 {
                acc = acc
                    .checked_add(name[k].checked_mul(p[k]).expect("no overflow"))
                    .expect("no overflow");
            }
            assert_eq!(
                acc,
                r(0),
                "a recorded point is off the name derived from it"
            );
        }

        // ★ (A plane with no points cannot be pushed at all any more — the "no points, no
        // name" arm retired with the old API; the type is the assertion now.)
    }

    /// ★★★★★ **A wide name interns — and opens no shortcut** (S2's whole behavioral change).
    ///
    /// The fixture is a triple whose **canonical answer** exceeds `i128` (cross-product terms
    /// multiply two ~2^90 coprime numerators — the same triple `nacre-scalar` locks as `Wide`).
    /// Before S2 such a plane got no name at all, so two statements of it were two handles;
    /// now it interns like any other.
    ///
    /// ★ Note the population: an axis-aligned cuboid — even on seventeen-digit corners — is NOT
    /// wide, because its canonical answers are tiny (`[10^13, 0, 0, −c]`); only the *narrow
    /// route's intermediates* overflow there, and the big route has named those since it exists.
    /// ★★ Measured while building S4's end-to-end lock: today's sketch→extrude walls cap out
    /// around ~115 bits (the decimal window bounds the products), so a wide name currently
    /// arises only from hand-stated triples like this one — datum planes (S5) are the coming
    /// production source.
    ///
    /// What a wide name still does not do: `narrow()` is `None`, so the **narrow shortcuts**
    /// (`base_rat`, integer Shewchuk, `Isometry` transport) decline exactly as they did on a
    /// missing name. Since S4 a wide name **does** host a sketch frame — through the
    /// arbitrary-precision road (`FramePlacement::Canonical`, locked in nacre-ops) — so what
    /// this pins is the shortcut boundary, not a frame one.
    #[test]
    fn a_wide_plane_interns_but_opens_no_narrow_shortcut() {
        let q = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
        let big1 = (1i128 << 90) + 1;
        let big2 = (1i128 << 90) + 3;
        let a = [q(big1, 3), q(big2, 7), q(0, 1)];
        let b = [q(-big2, 5), q(big1, 11), q(0, 1)];
        let c = [q(1, 13), q(1, 17), q(1, 19)];
        // Fixture qualification: the narrow route gives up on these points.
        assert_eq!(
            nacre_scalar::plane_through_points(a, b, c),
            None,
            "the narrow route was expected to overflow here — the case has stopped being the case"
        );

        let pl = nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap();
        let mut m = Model::new();
        let before = WIDE_PLANES.load(std::sync::atomic::Ordering::Relaxed);
        let (first, _) = m.push_plane(pl, [a, b, c], None);
        assert!(
            WIDE_PLANES.load(std::sync::atomic::Ordering::Relaxed) > before,
            "an answer past i128 must be counted as wide"
        );
        let name = m
            .surface_name
            .get(&first)
            .expect("a plane too wide for i128 must still be named");
        assert!(
            name.narrow().is_none(),
            "a wide name must not open the narrow shortcuts (base_rat, Shewchuk)"
        );

        // ★ The same plane stated again — permuted, even — is the same handle now.
        let (second, _) = m.push_plane(pl, [a, b, c], None);
        assert_eq!(first, second, "one wide plane, one handle");
        let (permuted, _) = m.push_plane(pl, [b, c, a], None);
        assert_eq!(
            first, permuted,
            "two spellings of one wide plane must intern"
        );
    }

    /// ★★ **A statement the name key cannot hold interns by the statement** (open item 16).
    ///
    /// The fixture reaches namelessness through `OnSeam` carriers — the cheapest population
    /// `plane_name_through` declines inside this crate. Semantically an ops producer would
    /// refuse this particular datum by cause; the door's contract is narrower ("store the
    /// statement, once") and holds for every nameless reason identically, which is what is
    /// pinned here. The real mixed-frame population is exercised end-to-end in `nacre-ops`.
    #[test]
    fn a_nameless_through_statement_interns_by_its_statement() {
        let mut m = Model::new();
        m.add_cylinder(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.0,
            2.0,
        );
        m.add_cuboid(
            Point3::from_array([4.0, 0.0, 0.0]),
            Point3::from_array([5.0, 1.0, 1.0]),
        );
        let mut vs: Vec<Handle<Vertex>> = Vec::new();
        let mut i = 0u32;
        while let Some(h) = m.vertex_handle_at(i) {
            i += 1;
            if matches!(m.vertices.get(h).def, VertexDef::OnSeam(_)) && vs.len() < 2 {
                vs.push(h);
            } else if matches!(m.vertices.get(h).def, VertexDef::ThreePlane(_)) && vs.len() == 2 {
                vs.push(h);
                break;
            }
        }
        let mut triple: [Handle<Vertex>; 3] = [vs[0], vs[1], vs[2]];
        triple.sort_by_key(|v| v.index());
        assert!(
            m.plane_name_through(triple).is_none(),
            "the fixture must be nameless, or this test measures the name road"
        );

        let cache = nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.5]),
            Vector3::from_array([0.0, 0.0, -1.0]),
        )
        .unwrap();
        let before = m.surface_count();
        let (h, flipped) = m.push_plane_through(cache, triple, None);
        assert!(!flipped);
        assert!(
            !m.surface_name.contains_key(&h),
            "a nameless statement must not invent a name"
        );
        assert_eq!(m.surface_count(), before + 1, "stored once");

        // The same statement again — one handle; and a cache built facing the other way is the
        // same plane with `flipped` reported, exactly as the name road reports it.
        let (again, flipped_same) = m.push_plane_through(cache, triple, None);
        assert_eq!(h, again, "one statement, one handle");
        assert!(!flipped_same);
        let reversed = nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap();
        let (still, flipped_now) = m.push_plane_through(reversed, triple, None);
        assert_eq!(h, still, "direction is not part of the statement");
        assert!(
            flipped_now,
            "but the survivor's other-way cache is reported"
        );
        assert_eq!(m.surface_count(), before + 1, "and nothing new was stored");

        // A different motion is a different statement — the same rule the name key keeps.
        let node = m.push_motion(
            Motion::Translate {
                offset: [
                    nacre_scalar::Rat::from_int(1),
                    nacre_scalar::Rat::from_int(0),
                    nacre_scalar::Rat::from_int(0),
                ],
            },
            None,
        );
        let (moved, _) = m.push_plane_through(cache, triple, Some(node));
        assert_ne!(h, moved, "the motion belongs in the statement key");
    }

    /// ★★★ S9: **the three world planes are born with the model** — deterministic handles,
    /// canonical truth, −axis caches — and `Default` is the same seeded model (the unseeded
    /// back door is closed).
    #[test]
    fn a_new_model_carries_the_three_world_planes() {
        let r = nacre_scalar::Rat::from_int;
        for m in [Model::new(), Model::default()] {
            assert_eq!(m.surface_count(), 3);
            let want = [
                (
                    nacre_scalar::Axis::Z,
                    [r(0), r(0), r(1), r(0)],
                    [0.0, 0.0, -1.0],
                ),
                (
                    nacre_scalar::Axis::X,
                    [r(1), r(0), r(0), r(0)],
                    [-1.0, 0.0, 0.0],
                ),
                (
                    nacre_scalar::Axis::Y,
                    [r(0), r(1), r(0), r(0)],
                    [0.0, -1.0, 0.0],
                ),
            ];
            for (i, (axis, name, cache_n)) in want.into_iter().enumerate() {
                let h = m.world_plane(axis);
                assert_eq!(h.index() as usize, i, "deterministic seed handles");
                assert_eq!(
                    m.surface_name.get(&h),
                    Some(&nacre_scalar::PlaneName::Narrow(name))
                );
                assert!(matches!(
                    m.surface_truth(h),
                    SurfaceTruth::Plane { motion: None, .. }
                ));
                let Surface::Plane(pl) = m.surface(h) else {
                    panic!("a seed is a plane")
                };
                assert_eq!(
                    pl.normal().as_array(),
                    cache_n,
                    "the cache points down the −axis (the base-cap sense)"
                );
            }
        }
    }

    /// ★★ S9: the seeds are the interning survivors — an origin cuboid's bottom/left/front
    /// faces carry the seed handles, so "the world plane" and "that face's plane" are one
    /// surface, stated once.
    #[test]
    fn an_origin_cuboids_axis_faces_intern_onto_the_seeds() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let face_surfaces: Vec<_> = m
            .shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .map(|&fh| m.faces.get(fh).surface)
            .collect();
        for axis in [
            nacre_scalar::Axis::Z,
            nacre_scalar::Axis::X,
            nacre_scalar::Axis::Y,
        ] {
            assert!(
                face_surfaces.contains(&m.world_plane(axis)),
                "the {axis:?}-normal face at 0 must be the seed itself"
            );
        }
        // And nothing pointless was minted: 3 seeds + the 3 off-origin faces.
        assert_eq!(m.surface_count(), 6);
    }

    /// ★★★ S6a: **a cylinder's caps record their three exact points** — the `add_cuboid`
    /// precedent applied to the last un-gated production path that minted point-less planes.
    /// The direct evidence that the record is a real name and not a dead entry: a cap that
    /// shares a plane with a box face **interns to the same handle**, which no point-less
    /// surface could ever do.
    /// `surface_handle_at` gives back the handle the arena already issued, and nothing more.
    ///
    /// The pair to this is [`Model::surface`]'s `compile_fail` lock, which still refuses to open
    /// the store: this adds a *name* for an index round-trip that `iter`-style access could
    /// already express, not a new power. Past the end is `None`, because the question it answers
    /// is existence.
    #[test]
    fn surface_handle_at_is_the_handle_the_arena_issued() {
        let mut m = Model::new();
        for axis in [
            nacre_scalar::Axis::Z,
            nacre_scalar::Axis::X,
            nacre_scalar::Axis::Y,
        ] {
            let h = m.world_plane(axis);
            assert_eq!(
                m.surface_handle_at(h.index()),
                Some(h),
                "a seeded plane's index names it back"
            );
        }
        let n = m.surface_count() as u32;
        assert_eq!(m.surface_handle_at(n), None, "past the end is None");
        assert_eq!(m.surface_handle_at(u32::MAX), None);

        // A surface pushed after the seeds is reachable by its own index too.
        let cuboid = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let s = m
            .faces
            .get(m.shells.get(m.solids.get(cuboid).outer).faces[0])
            .surface;
        assert_eq!(m.surface_handle_at(s.index()), Some(s));
    }

    #[test]
    fn a_cylinders_caps_record_points_and_intern_with_a_coplanar_face() {
        let mut m = Model::new();
        let cuboid = m.add_cuboid(
            Point3::from_array([-3.0, -3.0, 0.0]),
            Point3::from_array([-1.0, -1.0, 2.0]),
        );
        let cyl = m.add_cylinder(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.0,
            2.0,
        );
        m.rebuild_adjacency();
        // Every planar face of the cylinder carries points and a derived name.
        let planar: Vec<_> = m
            .shells
            .get(m.solids.get(cyl).outer)
            .faces
            .iter()
            .map(|&fh| m.faces.get(fh).surface)
            .filter(|&s| matches!(m.surface(s), Surface::Plane(_)))
            .collect();
        assert_eq!(planar.len(), 2, "two caps");
        for s in &planar {
            assert!(
                matches!(m.surface_truth(*s), SurfaceTruth::Plane { .. }),
                "a cap without truth"
            );
            assert!(m.surface_name.contains_key(s), "a cap without a name");
        }
        // The top cap lies on `z = 2`, the same plane as the box's top face — one plane, one
        // handle, across two producers.
        let box_top = m
            .shells
            .get(m.solids.get(cuboid).outer)
            .faces
            .iter()
            .map(|&fh| m.faces.get(fh).surface)
            .find(|s| {
                m.surface_name
                    .get(s)
                    .and_then(|n| n.narrow())
                    .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, -2.0])
            })
            .expect("the box top names z = 2");
        assert!(
            planar.contains(&box_top),
            "the coplanar cap and box face must intern to one handle"
        );
    }
    /// ★★★★★ **A plane can point at three vertices, and it is the same plane as the values.**
    ///
    /// The corners of a unit box name three coordinate planes each, so a datum through them is
    /// derivable — and the plane through `(1,0,0)`, `(0,1,0)`, `(0,0,1)` has a name like any
    /// other. That the name comes out at all is what makes every road below (frames, far caps,
    /// predicates) work unchanged: they read the name, and the name arrives at push.
    #[test]
    fn a_plane_stated_through_vertices_is_named_like_any_other() {
        let (mut m, vs) = box_corner_vertices();
        let cache = nacre_geom::Plane::from_point_normal(
            Point3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
        )
        .unwrap();
        let (h, _) = m.push_plane_through(cache, vs, None);
        let name = m
            .surface_name
            .get(&h)
            .expect("a Through plane must be named");
        assert_eq!(
            name.narrow().map(|c| c.map(|r| r.to_f64())),
            Some([1.0, 1.0, 1.0, -1.0]),
            "x + y + z = 1 is the plane through those three corners"
        );
        assert!(matches!(
            m.surface_truth(h),
            SurfaceTruth::Plane {
                points: PlanePoints::Through(_),
                ..
            }
        ));
    }

    /// ★★ **Interning does not care which way the plane was stated.** A `Known` push of the same
    /// plane returns the handle the `Through` push made — one plane, one handle, which is the
    /// property `PlanePoints` must not break by adding a variant.
    #[test]
    fn a_through_plane_and_a_known_plane_that_are_one_plane_share_a_handle() {
        let (mut m, vs) = box_corner_vertices();
        let cache = nacre_geom::Plane::from_point_normal(
            Point3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
        )
        .unwrap();
        let (through, _) = m.push_plane_through(cache, vs, None);
        let r = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
        let (known, _) = m.push_plane(
            cache,
            [
                [r(1, 1), r(0, 1), r(0, 1)],
                [r(0, 1), r(1, 1), r(0, 1)],
                [r(0, 1), r(0, 1), r(1, 1)],
            ],
            None,
        );
        assert_eq!(through, known, "one plane, two statements, one handle");
        // ★ And the first statement wins, as interning is documented to: the truth still points
        // at vertices. Which variant a plane ends up with is arena history, not a promise.
        assert!(matches!(
            m.surface_truth(known),
            SurfaceTruth::Plane {
                points: PlanePoints::Through(_),
                ..
            }
        ));
    }

    /// ★★★ **The reference goes backwards in time only (C5).** A datum is pushed after the
    /// vertices it names, and those vertices name surfaces pushed before *them* — so the walk
    /// plane → vertex → surface strictly descends in index and cannot cycle.
    #[test]
    fn a_through_plane_points_only_at_older_cells() {
        let (mut m, vs) = box_corner_vertices();
        let cache = nacre_geom::Plane::from_point_normal(
            Point3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
        )
        .unwrap();
        let (h, _) = m.push_plane_through(cache, vs, None);
        let SurfaceTruth::Plane {
            points: PlanePoints::Through(named),
            ..
        } = m.surface_truth(h)
        else {
            unreachable!()
        };
        for v in named {
            assert!(v.index() < m.vertices.len() as u32);
            let VertexDef::ThreePlane(tri) = m.vertices.get(*v).def else {
                unreachable!()
            };
            for s in tri {
                assert!(
                    s.index() < h.index(),
                    "a vertex's carrier must be older than the datum that names the vertex"
                );
            }
        }
    }

    /// Three corners of the unit box, sorted — each the meeting of three coordinate planes.
    fn box_corner_vertices() -> (Model, [Handle<Vertex>; 3]) {
        let mut m = Model::new();
        m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let mut want = Vec::new();
        for i in 0..m.vertices.len() as u32 {
            let vh = m.vertices.handle_at(i).unwrap();
            let c = m.vertex_point(vh).as_array();
            if c.iter().filter(|x| **x == 1.0).count() == 1
                && c.iter().all(|x| *x == 0.0 || *x == 1.0)
            {
                want.push(vh);
            }
        }
        want.sort_by_key(|v| v.index());
        let vs = [want[0], want[1], want[2]];
        (m, vs)
    }
    /// ★★★★ **Both counters see the second producer too.**
    ///
    /// `push_plane_through` began as a copy of `push_plane`'s tail and the copy silently dropped
    /// [`WIDE_PLANES`] and [`SEEDED_HITS`] — the census bridges `stat wide_planes` and
    /// `stat seeded_hits` went blind on the road the design predicts will *feed* them (open item
    /// 2b: *"the first producer of a `Wide` name is the datum"*). Sharing one interning tail is
    /// the fix; this is the observation that it worked, and it is two assertions because a
    /// counter that cannot move is indistinguishable from a population that never arrives.
    #[test]
    fn a_through_plane_is_counted_by_both_bridges() {
        use std::sync::atomic::Ordering::Relaxed;

        // ① A `Through` datum landing on a **seed** — three corners of a unit box on `z = 0`.
        let (mut m, vs) = seed_plane_corners();
        let cache = nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap();
        let seeded_before = SEEDED_HITS.load(Relaxed);
        let (h, _) = m.push_plane_through(cache, vs, None);
        assert!(
            (h.index() as usize) < 3,
            "three corners of the box's z = 0 face are the world XY seed"
        );
        assert!(
            SEEDED_HITS.load(Relaxed) > seeded_before,
            "a Through statement that interns onto a seed must be counted like any other"
        );

        // ② A `Through` datum whose name needs the wide vessel. The carriers are stated with
        // coprime ~2^90 numerators, so the canonical coefficients leave `i128` with no content to
        // divide out — the same construction `nacre-scalar` uses to reach that arm.
        let wide_before = WIDE_PLANES.load(Relaxed);
        let (mut m2, vs2) = wide_named_corner();
        let cache2 = nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
        )
        .unwrap();
        let (h2, _) = m2.push_plane_through(cache2, vs2, None);
        assert!(
            m2.surface_name
                .get(&h2)
                .is_some_and(|n| n.narrow().is_none()),
            "this fixture must actually produce a wide name, or ② measures nothing"
        );
        assert!(
            WIDE_PLANES.load(Relaxed) > wide_before,
            "a Through statement with a wide name must reach the wide-name bridge"
        );
    }

    /// Three corners of the unit box that lie on `z = 0` — the world XY seed.
    fn seed_plane_corners() -> (Model, [Handle<Vertex>; 3]) {
        let mut m = Model::new();
        m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let mut want = Vec::new();
        for i in 0..m.vertices.len() as u32 {
            let vh = m.vertices.handle_at(i).unwrap();
            if m.vertex_point(vh).as_array()[2] == 0.0 {
                want.push(vh);
            }
        }
        want.sort_by_key(|v| v.index());
        let vs = [want[0], want[1], want[2]];
        (m, vs)
    }

    /// A model whose three planes meet at a corner and whose *plane through* those corners needs
    /// the wide vessel. Built through the test door so the fixture is exactly three planes.
    fn wide_named_corner() -> (Model, [Handle<Vertex>; 3]) {
        let r = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
        let (b1, b2) = ((1i128 << 90) + 1, (1i128 << 90) + 3);
        let mut m = Model::new();
        // Nine planes: three per vertex, each triple meeting at a point whose coordinates carry
        // the coprime numerators, so the plane through the three points is wide.
        let mut vs = Vec::new();
        for (k, off) in [(0i128, 0i128), (1, 1), (2, 3)].iter().enumerate() {
            let _ = k;
            let (a, c) = (b1 + off.0, b2 + off.1);
            let tri = [
                axis_plane_at(&mut m, 0, r(1, a)),
                axis_plane_at(&mut m, 1, r(1, c)),
                axis_plane_at(&mut m, 2, r(off.0 + 1, b1)),
            ];
            let coord = Point3::from_array([
                1.0 / a as f64,
                1.0 / c as f64,
                (off.0 + 1) as f64 / b1 as f64,
            ]);
            vs.push(m.push_vertex(VertexDef::ThreePlane(tri), coord, None));
        }
        vs.sort_by_key(|v| v.index());
        let out = [vs[0], vs[1], vs[2]];
        (m, out)
    }

    /// The plane `x_axis = value`, pushed with its exact triple.
    fn axis_plane_at(m: &mut Model, axis: usize, value: nacre_scalar::Rat) -> Handle<Surface> {
        let z = nacre_scalar::Rat::from_int(0);
        let one = nacre_scalar::Rat::from_int(1);
        let mut pts = [[z; 3]; 3];
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        for (i, p) in pts.iter_mut().enumerate() {
            p[axis] = value;
            if i == 1 {
                p[u] = one;
            }
            if i == 2 {
                p[v] = one;
            }
        }
        let mut n = [0.0; 3];
        n[axis] = 1.0;
        let cache = nacre_geom::Plane::from_point_normal(
            Point3::from_array(core::array::from_fn(|k| {
                if k == axis { value.to_f64() } else { 0.0 }
            })),
            nacre_math::Vector3::from_array(n),
        )
        .unwrap();
        m.push_plane(cache, pts, None).0
    }
}
