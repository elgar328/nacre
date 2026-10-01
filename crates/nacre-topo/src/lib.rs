//! b-rep topology and the truth-only `Model` aggregate.
//!
//! The topology (`Vertex`/`Edge`/`Face`/`Loop`/`HalfEdge`/`Shell`/`Solid`)
//! references exact geometry only by `Handle` — geometry never knows about
//! topology, topology never inspects coordinates.
//!
//! [`Model`] is the **truth**: exact geometry stores + topology stores + the
//! derived [`Adjacency`] cache. It holds no tessellation and no operation log —
//! a mesh cache and the op log are companions owned at higher layers ("the
//! model is the replay result"; putting them here would make topo depend on
//! tess/ops and break the truth/cache split).

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
mod adjacency;
mod model;
mod topology;

pub use adjacency::{Adjacency, nonmanifold_vertices};
pub use topology::{Edge, EdgeDecline, Face, HalfEdge, Loop, Shell, Solid};

use nacre_exact::{Angle, Axis, HpBounded, Mag, Rat};
use nacre_geom::{Circle, Curve, Cylinder, Line, Plane};
use nacre_math::{Point3, Vector3};
use nacre_store::{Handle, Store};
use std::collections::{HashMap, HashSet};

/// How a discovered vertex is *defined* — the primitives whose intersection it
/// is. This definition is the **truth**; the vertex's `f64` point is
/// a cache derived from it. Sign decisions (in/out, orientation) feed this
/// definition to the indirect predicates rather than the cached coordinate
/// (M5), so a discovered vertex must carry it.
///
/// M5 polyhedral vertices are three-plane intersections; later milestones add
/// variants (a line∩plane point, quadric intersections). `Handle<Surface>` is
/// `Copy` regardless of `Surface`, so this stays `Copy`.
///
/// ★ **A 0-cell is its definition and nothing else.** The
/// coordinate, and what has been proven about it, live in the index-parallel point cache
/// ([`Model::vertex_point`] / [`Model::vertex_cache`]), filled by [`Model::push_vertex`] —
/// definition first, coordinate second. There is no origin tag (Constructed / Discovered /
/// Moved): such a tag never means exactness, and a moved vertex's motion is always the *faces'*
/// motion, which they record themselves. Nor is there a one-field `Vertex { def }` wrapper —
/// it would make every reader write `.def` to reach the only thing there.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Vertex {
    /// The intersection of three planes (their surface handles).
    ThreePlane([Handle<Surface>; 3]),
    /// A point on the intersection **curve** of two surfaces — the M3 cylinder seam vertex:
    /// the rim circle (lateral cylinder ∩ cap plane) at parameter `θ = 0`.
    ///
    /// ★ The pair pins a curve; what picks the point on it is the cylinder's `ref_dir`, which
    /// sits in the cylinder's truth ([`CylinderDef`]): `OnSeam([cylinder, cap])`
    /// **is** "the rim ∩ the `+ref_dir` ray" — a unique point, exactly designated. The
    /// definition is complete and `nacre_ops::realize_vertex` regenerates the coordinate from
    /// it; what remains is writing that back into the cache.
    /// M6 grows the vocabulary by variants, each stating its own truth — the invariants are
    /// per-variant; [`Vertex::Pierce`] is the first.
    OnSeam([Handle<Surface>; 2]),
    /// One of the (at most two) points where two planes' meet line crosses a cylinder's
    /// lateral surface — the structure says the carrier kinds, deliberately not an
    /// array of three lookalike handles (validate's carrier check becomes structural).
    ///
    /// ★ **`root` picks the point; its meaning is a convention.** The meet line's direction
    /// is `ℓ = n₁ × n₂` where `n₁`/`n₂` are the **canonical-name normals** of `planes[0]`/
    /// `planes[1]` **in stored (ascending-handle) order** — canonicalization fixes each
    /// normal's sign ("first nonzero component positive"), so ℓ is deterministic; `Lo`/`Hi`
    /// is ascending parameter along ℓ (`nacre_exact::quad`'s pair order). A tangency
    /// (double root) is one point and spells it [`QuadRoot::Double`] — **not `Lo`**, which
    /// would leave `Lo` unable to say which of the two it means. ★★ Anything
    /// that re-sorts the two plane handles (a transform remapping them, a mint site working in
    /// class indices) must restate `root` with it, and **[`QuadRoot::canonical`] is the one
    /// place that rule lives** — do not spell it again at the site.
    ///
    /// Minted in production by `nacre-ops`' arrangement (15 sites):
    /// the engine names the point as `NodeId::Pierce` in class-index space and the
    /// assembler restates it here in handle space — two canonical orders, one correspondence,
    /// established at the mint site through [`QuadRoot::canonical`].
    ///
    /// The name says what the point *is* — a line piercing a cylinder — not which root was
    /// picked (that is `root`); a name like "branch" would say only the latter.
    Pierce {
        /// The two cutting planes, ascending handle order (the `ThreePlane` precedent).
        planes: [Handle<Surface>; 2],
        /// The cylinder whose lateral surface the meet line crosses.
        cylinder: Handle<Surface>,
        /// Which of the two crossings, along the canonical line direction.
        root: QuadRoot,
    },
}

/// Which root of the two-point plane·plane·cylinder crossing a [`Vertex::Pierce`] means —
/// ascending parameter along the canonical meet-line direction (see `Pierce`'s doc for the
/// full convention). Definition vocabulary, so it lives here beside [`Vertex`], not in
/// scalar (whose pair is positional).
///
/// ★ **The declaration order is load-bearing.** The derived `Ord` is `Lo < Hi`, which *is* the
/// ascending-parameter convention, and comparison keys downstream read it (the arrangement's
/// reject witness is chosen as the smallest name). Reordering these two variants would silently
/// reverse those choices, so the order is a decision, not a listing.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum QuadRoot {
    /// The smaller parameter along the canonical meet-line direction.
    Lo,
    /// The larger parameter.
    Hi,
    /// **A tangency: the two roots coincide, so there is one point.**
    ///
    /// ★ It is a stored variant, not just a solver's report, because otherwise the stored value
    /// cannot say which it is: spelling a tangency `Lo` leaves
    /// `Lo` meaning either "the smaller of two" or "the only one", and then nothing holding a
    /// definition can tell whether a re-sort should toggle it. That ambiguity would make the
    /// remap rule below wrong for a tangency.
    Double,
}

impl QuadRoot {
    /// **The name the same point takes once the two planes are swapped.**
    ///
    /// Not "the other root" — the *other name for this point*, which is what a swap needs. For a
    /// genuine pair the swap sends `lo ↦ −hi` and `hi ↦ −lo` (see [`Self::canonical`]), so the
    /// two names trade places. For a tangency `s = mid` and the swap sends `mid ↦ −mid` along
    /// `ℓ ↦ −ℓ`: **the point goes to itself**, so its name does too.
    #[inline]
    pub fn flipped(self) -> QuadRoot {
        match self {
            QuadRoot::Lo => QuadRoot::Hi,
            QuadRoot::Hi => QuadRoot::Lo,
            QuadRoot::Double => QuadRoot::Double,
        }
    }

    /// **The one place "which root is it, once the pair is put in canonical order" is answered.**
    ///
    /// A [`Vertex::Pierce`] stores its two plane carriers ascending, and its root is defined
    /// against the meet line `ℓ = n₁ × n₂` of *that* order. So anything holding a solver's
    /// `(pair, root)` in some other order has to restate the root — and the restatement is a
    /// **derivation, not a convention**. Swapping the two planes in
    /// `nacre_exact::quad::plane_plane_cylinder` sends `ℓ ↦ −ℓ` while `base` is invariant (the
    /// Cramer system's row swap and its row-2 sign change cancel in numerator and denominator
    /// alike), leaves `a`, `c` and the discriminant alone and sends `b ↦ −b`. So the roots become
    /// `lo′ = −hi` and `hi′ = −lo`, and `base + lo′·(−ℓ) = base + hi·ℓ`: the swapped `Lo` and the
    /// original `Hi` are **the same point**.
    ///
    /// ★★ **A tangency is the exception, and it is not written here.** With `disc = 0` there is
    /// one point and the swap fixes it, so its name must not change — which is exactly what
    /// [`Self::flipped`] already says about [`QuadRoot::Double`]. The exception lives in the
    /// primitive, so this function has no special case to get wrong.
    ///
    /// ★ **Generic over the index space on purpose.** `Handle<Surface>` orders by its `index` and
    /// the arrangement names planes by bare `usize` **class** indices; both must answer the same
    /// way, and one function answering both is what keeps a second copy from being written.
    ///
    /// ★★★★★ **The rule this states is not about swapping — it is about `ℓ`.** `Lo`/`Hi` are the
    /// order along `ℓ = n₁ × n₂`, so *anything* that reverses `ℓ` trades the two names, and the
    /// swap is only the case this function can see. **Negating either normal is another**: the
    /// plane is the same set, and `plane_plane_cylinder` fixes the base by
    /// `{n₁·x = −d₁, n₂·x = −d₂, ℓ·x = 0}` — a condition `−ℓ` satisfies identically — so the base
    /// does not move and `ℓ` alone reverses. ⇒ **flip once per reversal; an even number is no flip
    /// at all.** A caller restating a name across two *spellings* of the same planes (rather than
    /// two orders of them) owes the sign half — `nacre-ops`' `pierce_name_from_def` is the one that
    /// does, coming back from handle space into class space.
    ///
    /// ☑ **The swap half is measured.** Where `assemble_fuse_cut` calls it, the population that
    /// exercises it is a boss whose classes arrive in the other order — a boolean's **result used
    /// as the next operand** (a second boolean builds its classes afresh) — and dropping the
    /// call there turns `an_operand_bounded_by_a_cylinder_is_named_in_class_space` red.
    #[inline]
    pub fn canonical<T: Ord>(pair: [T; 2], root: QuadRoot) -> ([T; 2], QuadRoot) {
        let [first, second] = pair;
        if second < first {
            ([second, first], root.flipped())
        } else {
            ([first, second], root)
        }
    }
}

impl Vertex {
    /// Every carrier handle the definition references, in stored order — the one spelling of
    /// "the surfaces this vertex is defined by" (reference integrity, the off-definition
    /// check and the remappability gate all ask exactly this; carrier *kinds* stay
    /// per-variant checks).
    pub fn carriers(&self) -> impl Iterator<Item = Handle<Surface>> {
        let (arr, n): ([Handle<Surface>; 3], usize) = match *self {
            Vertex::ThreePlane(s) => (s, 3),
            Vertex::OnSeam([a, b]) => ([a, b, b], 2),
            Vertex::Pierce {
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
    /// One axis-aligned rotation about the line through the rational `pivot`.
    ///
    /// ★ The field was called `point` while this very sentence called it the pivot — `Axis` gives
    /// only a direction, so a rotation needs *which line*, and that is what this names.
    Rotate {
        axis: Axis,
        pivot: [Rat; 3],
        angle: Angle,
    },
    /// One exact rational translation.
    Translate { offset: [Rat; 3] },
    /// One reflection in the coordinate plane `axis = offset` (`x ↦ 2·offset − x` on that axis).
    ///
    /// **Improper** — the only motion here with `det = −1`. It preserves lengths, distances and
    /// incidence like the others, but it *negates* every determinant of its images rather than
    /// leaving them alone. Judgments that answer a determinant question in a shared pre-motion
    /// frame must account for that; see `nacre-ops`' `BaseFrame` and `nacre-judge`'s `shared_base`.
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
    /// ★★★ **Where the frame sits on the plane is [`FramePlacement`]**: the derived
    /// `Canonical` convention by default, or the caller's `Named` values. See its own doc.
    ///
    /// ★★★ **`flip` is what makes the node a whole frame and not half of one.** Canonical plane
    /// coefficients carry **no direction** — the first nonzero component is forced positive,
    /// because their question is *"are these the same plane"*. A frame's `ŵ` is a direction, and
    /// two faces of one plane can face opposite ways; `flip` says to negate the coefficients, so
    /// the node names the sense as well as the plane. Without it a sketch on a reversed face comes
    /// out mirrored in `u` with its sweep running inward (measured: a tilted second boss came back
    /// `PadMissesFace`). It is decided from the truth by the consuming operation (the plane's
    /// facing, the motion's handedness, the side the use faces), never chosen by a caller.
    ///
    /// ★ **Proper** (`det = +1`) — `(u, v, w)` is right-handed by construction (`v = w × u`), so
    /// unlike [`Mirror`](Self::Mirror) it contributes nothing to a chain's parity. The coefficients
    /// are negated before the placement is read, so under `Canonical` (whose `+u` is derived from
    /// the normal) `ŵ` and `û` turn together and `v̂` stays, and under `Named` (whose `+u` is
    /// stated) `ŵ` and `v̂` turn and `û` stays — a half-turn either way: still proper.
    Frame {
        plane: Handle<Surface>,
        placement: FramePlacement,
        flip: bool,
    },
}

/// Where a sketch frame sits on its plane — the derived convention, or the caller's values.
/// The dichotomy mirrors `PlanePoints`: state a value when it can
/// be written, derive when it cannot.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FramePlacement {
    /// **The default.** The frame is a pure function of the plane — origin at the world
    /// origin's projection (`(−d/n·n)·n`), axes by the arbitrary-axis convention (`u = ẑ × n`,
    /// `ŷ × n` when the normal is exactly vertical, where `n` is the normal the frame faces — the
    /// coefficients with the node's `flip` spent) — derived when the chain is flattened and
    /// stored nowhere. That is what lets a plane whose canonical values overflow `i128` (or
    /// whose name is `Wide`) take the **same convention through arbitrary-precision
    /// realization** instead of falling to f64. One plane, one node — sketches
    /// on one face share it automatically. The derivation convention is frozen spec (changing
    /// it would silently turn every stored sketch).
    Canonical,
    /// The caller's own origin and `+u` direction, stated in the frame the plane's data lives
    /// in — what a caller-named plane (`PlaneDef`) supplies, validated at its construction
    /// (origin on the plane, `ref_dir` not parallel to the normal). `ref_dir` need not be unit
    /// length nor exactly in-plane; the realization projects it exactly.
    Named { origin: [Rat; 3], ref_dir: [Rat; 3] },
}

/// A node in the motion-history forest: one [`Motion`] applied to a solid, with a
/// parent link so several points can share a history's tail.
///
/// Stored in [`Model::motions`]; a moved surface's
/// moved surfaces name their leaf node ([`Surface`]'s motion slot). The tol a motion contributes is
/// application-point-dependent, so it is **not** stored here — judgment computes it by traversing
/// to the root.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MotionNode {
    pub motion: Motion,
    pub parent: Option<Handle<MotionNode>>,
}

/// What makes two surfaces the same plane, for `Model::surface_ids`: the canonical name
/// ([`nacre_exact::PlaneName`] — `Narrow | Wide`) and **the motion it is stated in**.
/// Identical names under different motions are different planes, because a name with no
/// motion speaks about the world and one with a motion about the pre-motion frame.
pub type SurfaceKey = (nacre_exact::PlaneName, Option<Handle<MotionNode>>);

/// The statement key for a plane the name key cannot hold: the sorted defining
/// triple and the motion it is stated under. See `Model::surface_through_ids`.
type ThroughKey = ([Handle<Vertex>; 3], Option<Handle<MotionNode>>);

/// The interning key for a cylinder — **deliberately conservative**: the whole exact
/// statement, `ref_dir` included, plus the motion it is stated under. Two statements of one
/// geometric cylinder with different `ref_dir`s stay two handles, because merging them would
/// split the seam (seam vertices and the seam edge cite the surface as their carrier). A key
/// this literal cannot merge wrongly; geometric identity across different statements is the
/// predicates' to answer per question. No `flipped` report either — a
/// literal-identical statement realizes to a literal-identical cache.
type CylinderKey = (CylinderDef, Option<Handle<MotionNode>>);

/// How many planes were named **`Wide`** — the canonical answer exceeded `i128` and took the
/// arbitrary-precision vessel. They intern
/// and carry identity like any other.
///
/// ★ **What `Wide` cannot do**: ride the narrow shortcuts —
/// [`nacre_exact::PlaneName::narrow`] is `None`, so they decline exactly as they decline on a
/// missing name. It does host a sketch frame, through the arbitrary-precision frame road.
/// The non-rotated coplanarity test reads the faces' own triangles and never looks at a name.
///
/// ★★ **Bounded by the type, not by the corpus.** A canonical name is a product of two point
/// differences, so `Rat = Ratio<i128>` inputs admit answers to roughly `2^2291`; today's models
/// stay far inside `i128` because they are written in short decimals, and the count reflects that
/// rather than any guarantee.
#[cfg(any(test, feature = "test-util"))]
pub static WIDE_PLANES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// How many pushes interned onto a **seeded world plane** (handles 0–2) — the census stat that
/// explains a plane-digest diff the wide counter cannot: a seeding-shaped change moves
/// survivors by *name collision*, not by width, and a falsifiability bridge needs a number for
/// that population too.
#[cfg(any(test, feature = "test-util"))]
pub static SEEDED_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// **How far the realization road reaches for surfaces** — the one road to a realization,
/// counted push by push. A vertex reaches it whole; a plane reaches it here wherever it has a world
/// name (anchor and unit normal from the truth), and elsewhere keeps the cache its pusher brought —
/// `nacre-ops` realizes what the truth answers and hands over the producer's figure for the rest,
/// which these counters do not see.
///
/// ⚠ `Model::apply_derivation` (private, so not a link) counts before it overwrites: `differs`
/// reads "how far the
/// producer's value was from the truth's", not "how far the cache moved after the fact".
///
/// Read through [`surface_derive_counts`]; the census prints them as `stat` rows, the same
/// falsifiability bridge [`WIDE_PLANES`] and [`SEEDED_HITS`] are. They are process-global and
/// never reset, so a number means "over everything this process built".
///
/// ★ **A test build's instrument** (`test-util`), like the two counters beside it: a product
/// build pushes surfaces without counting, and computes no derivation it does not apply. It
/// retires with the migration it measures — when `declined` is zero or justified, the producer's
/// cache leaves the push doors and nothing is left to count.
#[cfg(any(test, feature = "test-util"))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SurfaceDeriveCounts {
    /// Pushes whose cache the truth could derive (`Model::derive_surface_cache`, private — so
    /// this is deliberately not a link: a public field's doc cannot point inside the crate).
    pub derived: u64,
    /// Pushes where it declined — no world name (no name, a moved `Wide` name, a motion chain that
    /// does not fold — a frame, a turn off the quarters — or an overflow in the carry) or a moved
    /// cylinder. **This is the population that must reach zero (or be justified) before `fallback`
    /// can leave the push doors' signatures.**
    pub declined: u64,
    /// `declined`, split by cause. The split is what says *which* work removing the parameter
    /// needs, and the causes are not interchangeable: an unnamed plane is a fixture door, a moved
    /// `Wide` name wants a carry wider than `Rat`, a motion that still declines holds a frame or
    /// a turn off the quarters (no rational map exists — the rest fold, `Model::chain_plane_coeffs`),
    /// and arithmetic is an `i128` ceiling.
    pub declined_unnamed: u64,
    /// See [`SurfaceDeriveCounts::declined_unnamed`].
    pub declined_wide: u64,
    /// See [`SurfaceDeriveCounts::declined_unnamed`].
    pub declined_motion: u64,
    /// See [`SurfaceDeriveCounts::declined_unnamed`].
    pub declined_arith: u64,
    /// A cylinder whose statement is written before a motion — see
    /// [`SurfaceDeriveCounts::declined_unnamed`]. (A statement `Cylinder::from_axis` itself
    /// refuses would land here too; none has been seen.)
    pub declined_cylinder: u64,
    /// Of `derived`, how many differ bit-for-bit from what the producer stated — the census
    /// `in:` diff this cell predicts before it causes it.
    pub differs: u64,
}

#[cfg(any(test, feature = "test-util"))]
static SURFACE_DERIVED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DECLINED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DIFFERS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DECLINED_UNNAMED: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DECLINED_WIDE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DECLINED_MOTION: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DECLINED_ARITH: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(any(test, feature = "test-util"))]
static SURFACE_DECLINED_CYLINDER: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// A snapshot of [`SurfaceDeriveCounts`] — see its doc for what each number means.
#[cfg(any(test, feature = "test-util"))]
pub fn surface_derive_counts() -> SurfaceDeriveCounts {
    use std::sync::atomic::Ordering::Relaxed;
    SurfaceDeriveCounts {
        derived: SURFACE_DERIVED.load(Relaxed),
        declined: SURFACE_DECLINED.load(Relaxed),
        differs: SURFACE_DIFFERS.load(Relaxed),
        declined_unnamed: SURFACE_DECLINED_UNNAMED.load(Relaxed),
        declined_wide: SURFACE_DECLINED_WIDE.load(Relaxed),
        declined_motion: SURFACE_DECLINED_MOTION.load(Relaxed),
        declined_arith: SURFACE_DECLINED_ARITH.load(Relaxed),
        declined_cylinder: SURFACE_DECLINED_CYLINDER.load(Relaxed),
    }
}

/// A realized surface's bits, in one array, so "identical" means *identical* and not
/// `PartialEq`'s f64 equality (`-0.0 == 0.0`, and a `NaN` that never compares equal to itself).
/// A plane is its origin and unit normal — what the cache realizes, and all of it, so the last four
/// slots are zero; a cylinder is its axis origin and direction, its `ref_dir` and its radius.
#[cfg(any(test, feature = "test-util"))]
fn surface_bits(s: &nacre_geom::Surface) -> [u64; 10] {
    let pack = |a: [f64; 3], b: [f64; 3], c: [f64; 3], d: f64| {
        [
            a[0].to_bits(),
            a[1].to_bits(),
            a[2].to_bits(),
            b[0].to_bits(),
            b[1].to_bits(),
            b[2].to_bits(),
            c[0].to_bits(),
            c[1].to_bits(),
            c[2].to_bits(),
            d.to_bits(),
        ]
    };
    match s {
        nacre_geom::Surface::Plane(p) => {
            pack(p.origin().as_array(), p.normal().as_array(), [0.0; 3], 0.0)
        }
        nacre_geom::Surface::Cylinder(cy) => pack(
            cy.axis().origin().as_array(),
            cy.axis().direction().as_array(),
            cy.ref_dir().as_array(),
            cy.radius(),
        ),
    }
}

/// A sign against a reference direction: a face's against its surface's normal, and a plane's
/// ([`Surface::Plane::sense`]) against its points' direction. A pure tag — full derives.
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
    /// spelling: a ±1 map inline at each production site is a copy that can drift, which is
    /// exactly how a face comes to lie about which way it faces.
    #[inline]
    pub fn sign(self) -> i8 {
        match self {
            Orientation::Forward => 1,
            Orientation::Reversed => -1,
        }
    }
}

/// **A surface's exact truth** — what the surface *is*, as opposed to the f64
/// [`Surface`](nacre_geom::Surface) beside it, which is its realization.
///
/// The `motion` field says where the data is stated: `None` is the world, `Some` names the motion
/// history the data is stated *before*. Holding it inside the variant
/// is the point — a surface whose provenance is unrecorded, or whose exact form does not exist, is
/// **unrepresentable** here.
/// ★ **`Plane` is the large variant and it is not boxed** (the `Operation::Extrude`
/// precedent): planes dominate the arena — a prism is all planes, a cylinder contributes one
/// curved surface — so boxing the points would put an allocation and a pointer chase on the
/// common case to shrink the rare one.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    Plane {
        /// The plane's three exact points, stated in the frame `motion` names (the world when
        /// `None`).
        points: PlanePoints,
        /// The motion history carrying the points out to the world; `None` = the world itself.
        motion: Option<Handle<MotionNode>>,
        /// **Which way the plane faces, against its points.** The points' own direction is the
        /// world image of `n = (p₁ − p₀) × (p₂ − p₀)` — `chain_parity · L(n)`, where `L` is the
        /// linear part of the unfolded motion chain and the parity is `−1` per reflection
        /// (a reflection carries points, so it reverses the cross product it does not reverse
        /// as a direction). `Forward` says the plane's normal is that direction, `Reversed` its
        /// negation; a face's outward is the plane's normal times the face's own
        /// [`Orientation`].
        ///
        /// ★ **The sense is truth, not cache.** It used to live only in the f64 cache's normal,
        /// so the faces' flags leaned on a cache. The producer states it exactly — it cannot be
        /// read back off the cache without an `f64 → truth` arrow.
        ///
        /// ★ **Not part of any interning key.** One plane is one handle whichever way it was
        /// asked for; the first statement's sense stands and later askers learn how theirs
        /// relates through the `flipped` a push door returns.
        sense: Orientation,
    },
    /// A cylinder's exact truth: the rational statement of its lateral surface, beside
    /// the same motion slot a plane carries — a *moved* cylinder records its history.
    Cylinder {
        /// The exact statement, in the frame `motion` names (the world when `None`).
        def: CylinderDef,
        /// See [`Surface::Plane::motion`].
        motion: Option<Handle<MotionNode>>,
    },
}

/// A plane's three points — the kind of statement is the variant. `Known` carries values
/// (construction planes — walls, caps, caller-stated planes); `Through` points at model vertices
/// (datum planes).
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
    Known([[nacre_exact::Rat; 3]; 3]),
    /// **Three model vertices the plane passes through.** What a caller means by "the plane
    /// through those corners" — a coordinate read off a discovered vertex is rounded, and the
    /// plane built from rounded coordinates is a *different* plane (measured: on tilted geometry
    /// every one of 220 triples, `nacre-ops/tests/instruments/point_width.rs`).
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

/// **A cylinder's exact truth**: the lateral surface as its producer stated it — origin,
/// axis direction, seam reference direction and **squared** radius, all rational, in the
/// pre-motion frame.
///
/// ★ `dir` and `ref_dir` are **raw, unnormalized** — the `normal_def` precedent: normalizing
/// divides by an irrational length and would destroy the exact form. The f64 cache
/// ([`nacre_geom::Cylinder`]) holds the realized unit frame; this holds what the cylinder *is*.
/// The component form is the cylinder's analogue of a plane's *points* — a direct lift of the
/// caller's vocabulary plus component shuffles, never a derived product (the shape
/// forbidden for coefficients).
///
/// ★★ **`ref_dir` is model geometry, not a chart choice.** It fixes the seam permanently
/// (`θ = 0` on the `+ref_dir` side of the axis); the rational half-angle chart places its
/// own excluded point *on* this seam — the chart adapts to the seam, never the reverse.
///
/// Construction is checked ([`CylinderDef::new`]); fields are private so a literal cannot bypass
/// the check (the `SketchFrame` precedent).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct CylinderDef {
    origin: [Rat; 3],
    dir: [Rat; 3],
    ref_dir: [Rat; 3],
    /// `r²`, not `r`: the square of the distance from the axis to any rational point on the
    /// surface is rational, the radius itself need not be, and no exact predicate ever asks for
    /// the radius unsquared. A radius wanted as a number is a realization ([`Self::radius_f64`]),
    /// or the exact rational when the square has one ([`Self::radius_exact`]). **Wide**, because
    /// the square of a stated `Rat` need not fit `i128` (a 16-digit decimal below `1e-4`) and a
    /// statement is never refused for the width of its square.
    r2: nacre_exact::BigRat,
}

impl CylinderDef {
    /// The checked constructor — `None` when the statement means no cylinder, and **only then**:
    /// a zero `dir`, a non-positive squared radius `r2`, or a `ref_dir` with no component
    /// perpendicular to the axis (`ref_dir × dir = 0`, which a zero `ref_dir` satisfies too).
    ///
    /// ★ **Width is not a cause.** The parallelism test runs in
    /// [`nacre_exact::parallel_rat`], which clears denominators and answers in integers, so it
    /// cannot decline. Checked `Rat` would answer `None` on overflow — conflating *no cylinder*
    /// with *the arithmetic ran out*, and the callers below read `None` as the first: a statement
    /// whose axis carries a small component with a long decimal (denominator ~10²⁰, whose square
    /// leaves `i128`) would crash the constructor's `expect`. Measured population: 80% of computed
    /// near-axis-aligned
    /// directions, 0% of hand-written short decimals.
    pub fn new(
        origin: [Rat; 3],
        dir: [Rat; 3],
        ref_dir: [Rat; 3],
        r2: nacre_exact::BigRat,
    ) -> Option<Self> {
        let zero = Rat::from_int(0);
        if dir.iter().all(|c| *c == zero) || !r2.is_positive() {
            return None;
        }
        if nacre_exact::parallel_rat(&ref_dir, &dir) {
            return None;
        }
        Some(CylinderDef {
            origin,
            dir,
            ref_dir,
            r2,
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

    /// The squared radius, exact. Positive by construction.
    pub fn r2(&self) -> &nacre_exact::BigRat {
        &self.r2
    }

    /// The radius as an exact rational, when the squared radius has one — true of every radius a
    /// caller stated as a number, and of any arc through rational points whose `|start − centre|²`
    /// is a rational's square. `None` is not a refusal: the cylinder is as well-defined as any
    /// other, its radius merely has no rational spelling.
    pub fn radius_exact(&self) -> Option<Rat> {
        nacre_exact::rat_sqrt_exact_big(&self.r2)
    }

    /// The radius realized as an `f64`, correctly rounded — a rational radius' own `to_f64`,
    /// otherwise `√r²` at 128 → 256 bits (`nacre_exact::sqrt_f64`). The cache's number.
    pub fn radius_f64(&self) -> f64 {
        nacre_exact::sqrt_f64(&self.r2).expect("r² > 0 by construction")
    }
}

/// The truth-only aggregate: exact geometry + topology stores + the derived
/// adjacency cache. No `tess`, no `ops` (see the crate docs).
#[derive(Debug)]
pub struct Model {
    // exact geometry (truth)
    /// ★★★ **The truth, in the arena** — so a `Handle<Surface>` names what the surface *is*,
    /// not a realization of it. Total: a surface cannot enter without its truth, so there is no
    /// point-less surface.
    ///
    /// ★ Private: a surface can only enter through
    /// [`Model::push_plane`]/[`Model::push_cylinder`], which state its truth —
    /// a surface **without** a record is unrepresentable from outside this crate. Read through
    /// [`Model::surface`]/[`Model::surface_cache`]/[`Model::surface_count`]; there is
    /// deliberately no whole-store iterator (the arena keeps superseded surfaces — consumers
    /// walk the live faces).
    surfaces: Store<Surface>,
    /// Per-surface f64 caches, index-parallel to `surfaces` — **cache, not truth**: the truth
    /// above decides the realization, and a refinement pass may discard and regenerate the lot.
    /// Filled eagerly by [`Model::push_plane_raw`]/[`Model::push_cylinder_raw`]; read through
    /// [`Model::surface_cache`].
    ///
    /// ★★ **Why this is a private `Vec` and not a second `Store`**: a `Store` is append-only and
    /// sealed, so nothing could ever rewrite a cache entry at a higher precision. The cache has
    /// to be writable to be a cache at all — `edge_cache` and its
    /// [`Model::rebuild_edge_cache`] set that precedent.
    surface_cache: Vec<SurfaceCache>,
    /// The three seeded world planes, in normal-axis order Z(XY)·X(YZ)·Y(ZX) — captured at
    /// [`Model::new`] so [`Model::world_plane`] needs no handle minting. Always length 3.
    world_planes: Vec<Handle<Surface>>,
    /// Per-edge curve caches, index-parallel to `edges` — **cache, not truth**: the
    /// carriers and endpoints decide the curve ([`Model::derive_edge_curve`]), and
    /// [`Model::rebuild_edge_cache`] discards and regenerates the lot. Filled eagerly by
    /// [`Model::push_edge`]; read through [`Model::edge_curve`]. A raw `edges.push` without a
    /// cache entry desyncs the two — the accessor's debug_assert and validate's parallelism
    /// check watch for that (the store stays `pub`).
    edge_cache: Vec<EdgeCache>,
    /// Per-vertex coordinate caches, index-parallel to `vertices` — the coordinate and what
    /// has been proven about it. Filled by [`Model::push_vertex`]; read through
    /// [`Model::vertex_cache`] or its coordinate piece [`Model::vertex_point`].
    vertex_cache: Vec<PointCache>,
    /// The motion-history forest: motion definitions named by moved surfaces. Not geometry — a
    /// definition store.
    ///
    /// **Interned** — private, so writing through [`Model::push_motion`] is enforced by the
    /// type, not by discipline. Read through [`Model::motion`].
    motions: Store<MotionNode>,
    /// Interning table for [`Model::push_motion`]: the handle already issued for a given
    /// `(motion, parent)`. Not iterated (a `HashMap`'s order must never reach a result).
    motion_ids: HashMap<MotionNode, Handle<MotionNode>>,
    /// **Each motion node's whole chain, folded** — index-parallel to `motions`, filled by
    /// [`Model::push_motion`] as the node is born (its parent's fold, then its own), so asking a
    /// chain what it amounts to costs one read at any depth. `None` where the chain holds a node
    /// the rationals cannot state (a frame, a turn off the quarters), and for every child after.
    ///
    /// A cache: derived from the nodes, never an interning key, and no node is made or merged for
    /// it. Read through [`Model::chain_plane_coeffs`], [`Model::chain_point_rat`] and
    /// [`Model::chain_dir_rat`].
    ///
    /// ★ **Folded once, not walked per question** (measured): walking the chain at every question
    /// made 2,000 recorded quarter turns cost 6.7 s against 0.20 s, and a 4,200-turn irrational
    /// history 2.14 s against 1.44 s — the questions are asked per plane per transform, so a walk
    /// is quadratic in the history.
    motion_folds: Vec<Option<nacre_exact::AxisAffine>>,
    /// Each surface's **canonical name** ([`nacre_exact::PlaneName`]), derived from its points —
    /// present for every surface whose producer had a rational description to record.
    ///
    /// ★ **The point is that two statements of one plane get the same value.** `nacre_geom::Plane`
    /// keeps an un-normalized normal whose length follows the *face's size*, so the same plane
    /// reaches `coefficients()` as `[2.2, 0, 0, −6.6000000000000005]` from one face and
    /// `[13.2, 0, 0, −39.599999999999994]` from another — not exactly proportional, because `d`
    /// is a rounded product. Canonical names have no scale to disagree about, and the
    /// vessel is arbitrary-precision (`Narrow | Wide`) so **width cannot lose a name either**.
    ///
    /// ★★ **Built from the dimensions the user wrote, never lifted from the f64 coefficients
    /// above.** Lifting is lossless and useless here: it preserves the rounding, so the two
    /// vectors stay different (`nacre_exact::canonical_plane_coeffs`).
    ///
    /// ★★★ **The name is stated in the frame this surface's truth names** — the world for
    /// `motion: None`, and the **pre-motion** frame for a moved surface, whose world
    /// coefficients are irrational and so cannot be written
    /// down at all. A moved surface therefore inherits its source's name unchanged: the motion
    /// is recorded beside it, not folded into it.
    ///
    /// Absent is ordinary — an f64 construction path (no points to derive from).
    /// Iterate through the faces, never over the map. Arithmetic consumers (frames, `base_rat`,
    /// exact transports) read [`nacre_exact::PlaneName::narrow`]; `Wide` carries identity only.
    pub surface_name: HashMap<Handle<Surface>, nacre_exact::PlaneName>,
    /// Interning table for [`Model::push_plane`]: the handle already issued for a
    /// plane, keyed by its canonical coefficients **and the motion they are stated in**. The twin
    /// of [`Model::motion_ids`]; not iterated (a `HashMap`'s order must never reach a result).
    ///
    /// ★ **The motion belongs in the key.** Coefficients with no motion speak about the world and
    /// those with one about the pre-motion frame, so two identical arrays under different motions
    /// are different planes.
    ///
    /// ★★ **The witness does not.** Two faces of one plane sharing one moved
    /// witness is already how this works — `collect_planes` says the witness "was captured from
    /// whichever face first reached this surface" and winds it to *each* face's own outward
    /// normal — so it is not part of what makes two planes the same.
    surface_ids: HashMap<SurfaceKey, Handle<Surface>>,
    /// Interning for the planes the name key **cannot** hold: a `Through` statement whose three
    /// vertices no one frame solves rationally (mixed-frame datum) has no canonical
    /// name, so it interns by the **statement itself**: the sorted vertex triple and the motion
    /// it is stated under. The same plane may be named elsewhere — two handles, one plane — and
    /// the class discovery merges them by predicate.
    ///
    /// ★ This is statement identity, not geometric identity. Two *different* triples on one
    /// geometric plane get two handles here, by design: a nameless plane's geometric identity is the predicates' to answer, per question.
    /// What this table guarantees is the same thing construction-time sorting guarantees one
    /// level down — **the same statement never becomes two handles.**
    surface_through_ids: HashMap<ThroughKey, Handle<Surface>>,
    /// Interning for cylinders — by the whole exact statement plus motion; see
    /// [`CylinderKey`] for why the key is deliberately this literal.
    cylinder_ids: HashMap<CylinderKey, Handle<Surface>>,
    // topology (references geometry by Handle only)
    vertices: Store<Vertex>,
    edges: Store<Edge>,
    faces: Store<Face>,
    shells: Store<Shell>,
    solids: Store<Solid>,
    /// The live solids — the "current model" (supersede semantics).
    /// Editing ops supersede topology by pushing new cells and updating this
    /// list; the old cells stay in the append-only arena but, unreferenced by
    /// any live solid, drop out of the reachable closure. Producers register
    /// through [`Model::push_solid`]; `validate`/`Adjacency`/`nacre-step`
    /// traverse [`Model::reachable`], not the whole store.
    live_solids: Vec<Handle<Solid>>,
    // derived cache (rebuilt on demand)
    adj: Adjacency,
    /// **A pure accelerator — the one thing here that is neither truth nor cache.**
    ///
    /// A vertex's coordinate is realized by folding its motion chain from its rational base, so a
    /// solid whose history is `n` deep pays `1 + 2 + … + n` to build: every generation re-walks the
    /// prefix the generation before it already walked. This table holds the high-precision value a
    /// chain reaches at one node, so the next motion can carry on from it instead of starting over.
    ///
    /// ★ **The key is the definition** — `(base, leaf, prec)` — and everything else follows from
    /// that. The value is a pure function of the key, and motion nodes are interned append-only, so
    /// an entry **cannot go stale**: there is no invalidation rule to get wrong. A motion the
    /// statements absorbed moves the base, so it lands on a *different* key and simply misses; a
    /// motion that fixes its plane records no node; a world-stated triple has no leaf. All three
    /// "do not extend the prefix" cases fall out of the key rather than out of a test.
    ///
    /// ☑ **Empty is always correct.** Nothing here is needed for an answer — drop it, and the next
    /// realization folds from the base as it always did, to the same bits. That is what makes it
    /// safe to evict on a hit, and why [`Model::clear_prefix_hp`] owes no one an explanation.
    prefix_hp: HashMap<PrefixKey, PrefixValue>,
}

/// What a remembered prefix is filed under — the **definition**, never a vertex: the rational base,
/// the chain node it was folded to, and the precision it was folded at.
///
/// ★ Spelled here rather than in `nacre-ops` because this is the table's own shape, and the reader
/// and the writer live in different crates. Two spellings of one key are two keys that drift.
pub type PrefixKey = ([Rat; 3], Handle<MotionNode>, usize);

/// A remembered prefix: **how many chain nodes were folded to reach it**, and the point.
///
/// ⚠ The count is not "how many motion nodes up". One `Motion::Frame` expands into several
/// `MoveNode`s (`nacre_ops`'s `motion_chain` appends a whole sub-chain for it), so counting parent
/// steps would slice the wrong suffix the moment a frame is in the history — and a solid built on a
/// sketch frame and then moved has exactly that shape. The length is what the reader needs and what
/// the writer already knows, so it travels with the value.
pub type PrefixValue = (usize, [HpBounded; 3]);

/// One vertex's realized coordinate — a **cache** beside the vertex store (index-parallel),
/// filled by [`Model::push_vertex`]. The variant says what the cache **knows** about the
/// coordinate, never where the vertex came from, and there are only three things it can know:
///
/// - [`Self::Bounded`] — the coordinate was **realized from the definition** and rounded once
///   (`nacre_ops::realize_vertex`), so `bound` is a *proven* per-axis containment: the truth lies
///   within `coord ± bound`. Every vertex an operation makes carries this unless the cache road
///   declined.
/// - [`Self::Ceiling`] — the cache road **stopped**, and a road willing to pay more can still
///   answer. Two walls stop it and they say the same thing to a reader: the ladder's first two
///   rungs did not decide the coordinate, or a chain the road would have to replay ran past what the cache
///   road pays for (`nacre_ops`'s cost cap — a chain whose fold answers is read at any depth).
///   `nacre_ops::refine_vertex_cache` is what lifts these.
/// - [`Self::Unrealized`] — **no realization stands behind the coordinate**: the kernel has no road
///   to it (the realization declines by name), or nobody asked (a hand-built fixture). The figure
///   is the construction's own, and a checker applies its construction epsilon.
///
/// ⚠★★★ **There is no stored tolerance.** A residual is **not** a bound: it is one distance,
/// from the point to its carriers, and says nothing about how
/// far each coordinate sits from the truth — a near-degenerate crossing can be close to every
/// carrier and far from the exact corner. Nothing carries it: what a cache is for is saying what
/// has been **proven** about the coordinate, and the arrangement's self-touch sieve reads this
/// bound too.
///
/// Unlike [`EdgeCache`] there is no discard-and-rebuild: a vertex is realized when it is pushed, so
/// the cache is the realization's memo from the start. What can happen later is a **refinement** —
/// `Ceiling` raised to `Bounded` — which only ever makes a coordinate more accurate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PointCache {
    /// Realized from the definition: the truth lies within `coord ± bound` on every axis.
    Bounded { coord: Point3, bound: [Mag; 3] },
    /// The cache road stopped — too few bits on its one rung, or a history past its cost cap.
    /// A paid realization can still answer; until one does, `coord` is the construction's figure.
    Ceiling { coord: Point3 },
    /// No realization stands behind `coord`: no road to it, or nobody asked.
    Unrealized { coord: Point3 },
}

impl PointCache {
    /// The cached coordinate, whichever of the three the cache knows about it.
    #[inline]
    pub fn coord(&self) -> Point3 {
        match *self {
            PointCache::Bounded { coord, .. }
            | PointCache::Ceiling { coord }
            | PointCache::Unrealized { coord } => coord,
        }
    }

    /// The proven per-axis bound, `None` unless the coordinate was realized from the definition.
    #[inline]
    pub fn bound(&self) -> Option<&[Mag; 3]> {
        match self {
            PointCache::Bounded { bound, .. } => Some(bound),
            PointCache::Ceiling { .. } | PointCache::Unrealized { .. } => None,
        }
    }
}

/// One edge's realized curve — a **cache** beside the edge store (index-parallel), derived
/// from the edge's carriers and endpoints by [`Model::derive_edge_curve`]. Discard and
/// regenerate with [`Model::rebuild_edge_cache`]; the field is private so the only writers
/// are the derivation itself.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeCache {
    curve: Curve,
}

/// One surface's realized geometry — a **cache** beside the surface store (index-parallel),
/// the f64 answer to what [`Model::surface`] states exactly.
///
/// ★ **[`EdgeCache`]'s mirror**: topo wraps, and the geometry's own methods stay in
/// `nacre-geom` — `distance`, `normal_at`, `translated` and the rest dispatch on
/// [`nacre_geom::Surface`], which is geom's vocabulary and stays there. The wrapper exists so
/// the cache is a named thing with room to grow: a measured `tol` (the surface analogue of
/// [`PointCache`]'s) arrives with the refinement pass that can produce it.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceCache {
    realized: nacre_geom::Surface,
}

/// The handles reachable from a model's live solids — the live model.
///
/// Only the sets consumers need today: `validate`'s Euler counts vertices/edges/
/// faces/shells, and `Adjacency`/loop/incidence walk faces/edges. Surfaces and
/// curves are not tracked: the M4 face ops never orphan a *used* one, and the three
/// seeded world planes are deliberately face-less — traversal through
/// definitions (world planes, datum references) is not built.
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
    /// Two callers want it for opposite reasons. `Model::reversed_shell` (a `test-util` door)
    /// pairs it with an
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

/// `Default` is [`Model::new`] — **seeded**. A derived `Default` would hand out an *empty* model —
/// one without the world planes: a public back door to the state the
/// seeding exists to remove. There is deliberately no unseeded constructor.
impl Default for Model {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[path = "tests/lib.rs"]
mod tests;
