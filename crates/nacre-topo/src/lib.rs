//! b-rep topology and the truth-only `Model` aggregate.
//!
//! The topology (`Vertex`/`Edge`/`Face`/`Loop`/`HalfEdge`/`Shell`/`Solid`)
//! references exact geometry only by `Handle` — geometry never knows about
//! topology, topology never inspects coordinates (overview, principle 2).
//!
//! [`Model`] is the **truth**: exact geometry stores + topology stores + the
//! derived [`Adjacency`] cache. It holds no tessellation and no operation log —
//! a mesh cache and the op log are companions owned at higher layers ("the
//! model is the replay result"; putting them here would make topo depend on
//! tess/ops and break the truth/cache split).

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
mod adjacency;
mod topology;

pub use adjacency::{Adjacency, nonmanifold_vertices};
pub use topology::{Edge, Face, HalfEdge, Loop, Shell, Solid};

use nacre_geom::{Circle, Curve, Cylinder, Line, Plane};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Angle, Axis, HpBounded, Mag, Rat};
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
    /// is ascending parameter along ℓ (`nacre_scalar::quad`'s pair order). A tangency
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
    /// `nacre_scalar::quad::plane_plane_cylinder` sends `ℓ ↦ −ℓ` while `base` is invariant (the
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
    /// ★★★ **Where the frame sits on the plane is [`FramePlacement`]**: the derived
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
/// The dichotomy mirrors `PlanePoints`: state a value when it can
/// be written, derive when it cannot.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FramePlacement {
    /// **The default.** The frame is a pure function of the plane — origin at the world
    /// origin's projection (`(−d/n·n)·n`), axes by the arbitrary-axis convention (`u = ẑ × n`,
    /// `ŷ × n` when the normal is exactly vertical) — derived when the chain is flattened and
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

/// A node in the motion-history forest (design §CIP ⑦): one [`Motion`] applied to a solid, with a
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
/// ([`nacre_scalar::PlaneName`] — `Narrow | Wide`) and **the motion it is stated in**.
/// Identical names under different motions are different planes, because `Constructed` names
/// speak about the world and `Moved` ones about the pre-motion frame.
pub type SurfaceKey = (nacre_scalar::PlaneName, Option<Handle<MotionNode>>);

/// The statement key for a plane the name key cannot hold: the sorted defining
/// triple and the motion it is stated under. See `Model::surface_through_ids`.
type ThroughKey = ([Handle<Vertex>; 3], Option<Handle<MotionNode>>);

/// The interning key for a cylinder — **deliberately conservative**: the whole exact
/// statement, `ref_dir` included, plus the motion it is stated under. Two statements of one
/// geometric cylinder with different `ref_dir`s stay two handles, because merging them would
/// split the seam (seam vertices and the seam edge cite the surface as their carrier). A key
/// this literal cannot merge wrongly; geometric identity across different statements is the
/// predicates' to answer per question (rule 6). No `flipped` report either — a
/// literal-identical statement realizes to a literal-identical cache.
type CylinderKey = (CylinderDef, Option<Handle<MotionNode>>);

/// How many planes were named **`Wide`** — the canonical answer exceeded `i128` and took the
/// arbitrary-precision vessel. They intern
/// and carry identity like any other.
///
/// ★ **What `Wide` cannot do**: ride the narrow shortcuts —
/// [`nacre_scalar::PlaneName::narrow`] is `None`, so they decline exactly as they decline on a
/// missing name. It does host a sketch frame, through the arbitrary-precision frame road.
/// The non-rotated coplanarity test reads the faces' own triangles and never looks at a name.
///
/// ★★ **Bounded by the type, not by the corpus.** A canonical name is a product of two point
/// differences, so `Rat = Ratio<i128>` inputs admit answers to roughly `2^2291`; today's models
/// stay far inside `i128` because they are written in short decimals, and the count reflects that
/// rather than any guarantee.
pub static WIDE_PLANES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// How many pushes interned onto a **seeded world plane** (handles 0–2) — the census stat that
/// explains a plane-digest diff the wide counter cannot: a seeding-shaped change moves
/// survivors by *name collision*, not by width, and a falsifiability bridge needs a number for
/// that population too.
pub static SEEDED_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// **How far the realization road reaches for surfaces** — principle 4's "one road to a
/// realization", counted push by push. A vertex reaches it whole; a plane reaches it
/// **halfway**: the anchor is a function of the truth and the row and the
/// sense stay the producer's.
///
/// ⚠ `Model::apply_derivation` (private, so not a link) counts before it overwrites: `differs`
/// reads "how far the
/// producer's value was from the truth's", not "how far the cache moved after the fact".
///
/// Read through [`surface_derive_counts`]; the census prints them as `stat` rows, the same
/// falsifiability bridge [`WIDE_PLANES`] and [`SEEDED_HITS`] are. They are process-global and
/// never reset, so a number means "over everything this process built".
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SurfaceDeriveCounts {
    /// Pushes whose cache the truth could derive (`Model::derive_surface_cache`, private — so
    /// this is deliberately not a link: a public field's doc cannot point inside the crate).
    pub derived: u64,
    /// Pushes where it declined — no name, a `Wide` name, a motion that is not a rational
    /// translation chain, a moved cylinder, or an overflow. **This is the population that must
    /// reach zero (or be justified) before `cache` can leave the push doors' signatures.**
    pub declined: u64,
    /// `declined`, split by cause. The split is what says *which* work removing the parameter
    /// needs, and the causes are not interchangeable: an unnamed plane is a fixture door, a
    /// `Wide` name wants an arbitrary-precision arm, a motion wants a chain folded to an
    /// `Isometry` (rotations included — `nacre_scalar::Isometry::plane_coeffs` already carries
    /// a plane through the 90° family), and arithmetic is an `i128` ceiling.
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
    /// Interning hits whose **discarded** cache differs bit-for-bit from the survivor's — the
    /// only measurement of the determinism claim ("the same geometry writes the same file") on
    /// the *product* population rather than on an experiment built to show it.
    pub discarded_differing: u64,
}

static SURFACE_DERIVED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static SURFACE_DECLINED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static SURFACE_DIFFERS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static SURFACE_DISCARDED_DIFFERING: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);
static SURFACE_DECLINED_UNNAMED: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);
static SURFACE_DECLINED_WIDE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static SURFACE_DECLINED_MOTION: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static SURFACE_DECLINED_ARITH: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static SURFACE_DECLINED_CYLINDER: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// A snapshot of [`SurfaceDeriveCounts`] — see its doc for what each number means.
pub fn surface_derive_counts() -> SurfaceDeriveCounts {
    use std::sync::atomic::Ordering::Relaxed;
    SurfaceDeriveCounts {
        derived: SURFACE_DERIVED.load(Relaxed),
        declined: SURFACE_DECLINED.load(Relaxed),
        differs: SURFACE_DIFFERS.load(Relaxed),
        discarded_differing: SURFACE_DISCARDED_DIFFERING.load(Relaxed),
        declined_unnamed: SURFACE_DECLINED_UNNAMED.load(Relaxed),
        declined_wide: SURFACE_DECLINED_WIDE.load(Relaxed),
        declined_motion: SURFACE_DECLINED_MOTION.load(Relaxed),
        declined_arith: SURFACE_DECLINED_ARITH.load(Relaxed),
        declined_cylinder: SURFACE_DECLINED_CYLINDER.load(Relaxed),
    }
}

/// A realized surface's bits, in one array, so "identical" means *identical* and not
/// `PartialEq`'s f64 equality (`-0.0 == 0.0`, and a `NaN` that never compares equal to itself).
/// Ten numbers either way: a plane's origin, unit normal and four coefficients; a cylinder's
/// axis origin and direction, its `ref_dir` and its radius.
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
            let c = p.coefficients();
            pack(
                p.origin().as_array(),
                p.normal().as_array(),
                [c[0], c[1], c[2]],
                c[3],
            )
        }
        nacre_geom::Surface::Cylinder(cy) => pack(
            cy.axis().origin().as_array(),
            cy.axis().direction().as_array(),
            cy.ref_dir().as_array(),
            cy.radius(),
        ),
    }
}

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
/// [`Surface`](nacre_geom::Surface) beside it, which is its realization.
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
pub enum Surface {
    Plane {
        /// The plane's three exact points, stated in the frame `motion` names (the world when
        /// `None`) — the value `Model::surface_points` used to carry.
        points: PlanePoints,
        /// The motion history carrying the points out to the world; `None` = the world itself.
        motion: Option<Handle<MotionNode>>,
    },
    /// A cylinder's exact truth: the rational statement of its lateral surface, beside
    /// the same motion slot a plane carries — a *moved* cylinder records its history instead of
    /// silently degrading (the old side-table path demoted it to `Inexact`).
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
    r2: nacre_scalar::BigRat,
}

impl CylinderDef {
    /// The checked constructor — `None` when the statement means no cylinder, and **only then**:
    /// a zero `dir`, a non-positive squared radius `r2`, or a `ref_dir` with no component
    /// perpendicular to the axis (`ref_dir × dir = 0`, which a zero `ref_dir` satisfies too).
    ///
    /// ★ **Width is not a cause.** The parallelism test runs in
    /// [`nacre_scalar::parallel_rat`], which clears denominators and answers in integers, so it
    /// cannot decline. It used to run in checked `Rat` and answer `None` on overflow — a
    /// "conservative refusal" that conflated *no cylinder* with *the arithmetic ran out*, and
    /// the callers below read it as the first: a statement whose axis carries a small component
    /// with a long decimal (denominator ~10²⁰, whose square leaves `i128`) crashed the
    /// constructor's `expect`. Measured population: 80% of computed near-axis-aligned
    /// directions, 0% of hand-written short decimals.
    pub fn new(
        origin: [Rat; 3],
        dir: [Rat; 3],
        ref_dir: [Rat; 3],
        r2: nacre_scalar::BigRat,
    ) -> Option<Self> {
        let zero = Rat::from_int(0);
        if dir.iter().all(|c| *c == zero) || !r2.is_positive() {
            return None;
        }
        if nacre_scalar::parallel_rat(&ref_dir, &dir) {
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
    pub fn r2(&self) -> &nacre_scalar::BigRat {
        &self.r2
    }

    /// The radius as an exact rational, when the squared radius has one — true of every radius a
    /// caller stated as a number, and of any arc through rational points whose `|start − centre|²`
    /// is a rational's square. `None` is not a refusal: the cylinder is as well-defined as any
    /// other, its radius merely has no rational spelling.
    pub fn radius_exact(&self) -> Option<Rat> {
        nacre_scalar::rat_sqrt_exact_big(&self.r2)
    }

    /// The radius realized as an `f64`, correctly rounded — a rational radius' own `to_f64`,
    /// otherwise `√r²` at 128 → 256 bits (`nacre_scalar::sqrt_f64`). The cache's number.
    pub fn radius_f64(&self) -> f64 {
        nacre_scalar::sqrt_f64(&self.r2).expect("r² > 0 by construction")
    }
}

/// Why `Model::add_cylinder_exact` (a `test-util` door) refused a statement — **every way in,
/// named**.
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
#[cfg(any(test, feature = "test-util"))]
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
    /// ★★★ **The truth, in the arena** — so a `Handle<Surface>` names what the surface *is*,
    /// not a realization of it. Total: a surface cannot enter without its truth, which is what
    /// retired `SurfaceDef`/`Inexact` and the point-less population.
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
    /// Each surface's **canonical name** ([`nacre_scalar::PlaneName`]), derived from its points —
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
    /// vectors stay different (`nacre_scalar::canonical_plane_coeffs`).
    ///
    /// ★★★ **The name is stated in the frame this surface's truth names** — the world for
    /// `motion: None`, and the **pre-motion** frame for a moved surface, whose world
    /// coefficients are irrational and so cannot be written
    /// down at all. A moved surface therefore inherits its source's name unchanged: the motion
    /// is recorded beside it, not folded into it.
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
    /// world coefficients are irrational (mixed-frame datum) has no canonical
    /// name, so it interns by the **statement itself**: the sorted vertex triple and the motion
    /// it is stated under.
    ///
    /// ★ This is statement identity, not geometric identity. Two *different* triples on one
    /// geometric plane get two handles here, and rule 6's qualification is exactly that:
    /// a nameless plane's geometric identity is the predicates' to answer, per question.
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
///   answer. Two walls stop it and they say the same thing to a reader: the ladder's first rung did
///   not decide the coordinate, or the motion history ran past what the cache road pays for
///   (`nacre_ops`'s cost cap). `nacre_ops::refine_vertex_cache` is what lifts these.
/// - [`Self::Unrealized`] — **no realization stands behind the coordinate**: the kernel has no road
///   to it (the realization declines by name), or nobody asked (a hand-built fixture). The figure
///   is the construction's own, and a checker applies its construction epsilon.
///
/// ⚠★★★ **There is no stored tolerance, and the variant that held one is gone.** An earlier cut
/// kept the residual the arrangement measured (`Measured { residual }`). That number is **not** a
/// bound: a residual is one distance, from the point to its carriers, and says nothing about how
/// far each coordinate sits from the truth — a near-degenerate crossing can be close to every
/// carrier and far from the exact corner. The arrangement still measures it (`SeamVertex.tol`,
/// which the self-touch sieve reads); the *cache* no longer stores it, because what a cache is for
/// is saying what has been **proven** about the coordinate.
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

impl Model {
    /// A model with the three **world axis planes pre-seeded** — surface handles 0 (XY, z = 0),
    /// 1 (YZ, x = 0), 2 (ZX, y = 0), deterministic so a replayed log and a live session name the
    /// same planes.
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
            surface_cache: Vec::new(),
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
            prefix_hp: HashMap::new(),
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

    /// The prefix accelerator's value at one chain node — `(nodes folded, the point)` — if it is
    /// remembered. `None` is never an error: the caller folds from the base instead.
    #[inline]
    pub fn prefix_hp(
        &self,
        base: [Rat; 3],
        leaf: Handle<MotionNode>,
        prec: usize,
    ) -> Option<&PrefixValue> {
        self.prefix_hp.get(&(base, leaf, prec))
    }

    /// **Take one prefix value and leave another — "consumed" is how this table evicts.**
    ///
    /// A remembered prefix is read by exactly one successor (a transform maps one vertex to one
    /// vertex), so removing what was used and inserting what was produced keeps the table at a
    /// single live generation without tracking generations at all. `used` is the key a reader hit,
    /// `None` when it folded from the base.
    ///
    /// ⚠ **Only a caller that is extending a chain may insert**, which is why this is one door and
    /// not two: a vertex minted fresh by an arrangement is nobody's prefix, and remembering it
    /// would grow the table by one entry per boolean, forever.
    pub fn hand_over_prefix_hp(
        &mut self,
        used: Option<PrefixKey>,
        key: PrefixKey,
        value: PrefixValue,
    ) {
        if let Some(used) = used {
            self.prefix_hp.remove(&used);
        }
        self.prefix_hp.insert(key, value);
    }

    /// Drop every remembered prefix. Costs nothing but time: the next realization folds from the
    /// base and reaches the same bits (`Model::prefix_hp`'s contract).
    ///
    /// ☑ Production-facing on purpose even though nothing in this workspace calls it yet: it is
    /// the door that makes "empty is always correct" usable rather than merely true, and a
    /// consumer holding a long-lived model is the caller it is for.
    #[cfg(any(test, feature = "test-util"))]
    pub fn clear_prefix_hp(&mut self) {
        self.prefix_hp.clear();
    }

    /// How many prefixes are remembered — the instrument behind "the table stays bounded by the
    /// live generation, not by history length".
    ///
    /// ⚠ **Test-gated because it has no production consumer**, the same reason
    /// [`Model::push_plane_unregistered`] is: this counts an accelerator's internals, which is a
    /// thing to assert about, not a thing to build on. An ungated `pub` here would ship a
    /// permanent public API through the `nacre` facade for a caller that does not exist.
    #[cfg(any(test, feature = "test-util"))]
    #[inline]
    pub fn prefix_hp_len(&self) -> usize {
        self.prefix_hp.len()
    }

    /// Recompute the [`Adjacency`] cache from the current topology stores.
    /// Call once after a batch of additions (the cache is otherwise stale).
    pub fn rebuild_adjacency(&mut self) {
        // Build against an immutable borrow, then move into place — avoids
        // borrowing `self.adj` mutably while iterating the other stores.
        let adj = Adjacency::rebuild(&*self);
        self.adj = adj;
    }

    /// Push a solid into the store **and mark it live**. This is the
    /// blessed way for a producer to add a solid; the reachable closure
    /// ([`Model::reachable`]) grows to include it. Editing ops instead mutate
    /// [`Model::live_solids`] directly (drop the superseded solid, add the new).
    pub fn push_solid(&mut self, solid: Solid) -> Handle<Solid> {
        let h = self.solids.push(solid);
        self.live_solids.push(h);
        h
    }

    /// The live solids — the "current model" (supersede semantics).
    #[inline]
    pub fn live_solids(&self) -> &[Handle<Solid>] {
        &self.live_solids
    }

    /// **Drop these solids from the live set** — supersede, the editing ops' half.
    ///
    /// ⚠★★★ **Order-preserving on the survivors, and that is load-bearing, not incidental.**
    /// `nacre_step::to_step` exports the live set **in order**, so permuting it here would
    /// silently permute the exported STEP entities — and nothing would catch that: the census
    /// reads the arena (it never sees live order) and there is no golden STEP text anywhere.
    /// `Vec::retain` keeps relative order, which is why this is spelled with it.
    ///
    /// ★ Takes what to **drop**, not a keep-predicate: every caller reads "supersede these", and a
    /// predicate would make the call site say the opposite of the name.
    pub fn supersede_live(&mut self, drop: &[Handle<Solid>]) {
        self.live_solids.retain(|h| !drop.contains(h));
    }

    /// **Make a solid that is already in the arena live again.**
    ///
    /// ★ [`Model::push_solid`]'s doc has always described this move — *"editing ops instead mutate
    /// live_solids directly (drop the superseded solid, **add the new**)"* — but there was no door
    /// for the second half, so callers reached for the field. The population is the reject paths
    /// in `ops`, which retire the operands through a boolean and then have to put the original
    /// back when the op itself refuses.
    pub fn make_live(&mut self, h: Handle<Solid>) {
        debug_assert!(
            (h.index() as usize) < self.solids.len(),
            "a solid the arena does not hold cannot be live"
        );
        self.live_solids.push(h);
    }

    /// **Put the live set back** — the rollback half of a rejected operation.
    ///
    /// ★ This exists so the transaction has a name. The idiom it replaces (`let snapshot =
    /// …clone()` … `model.live_solids = snapshot`) is syntactically just an assignment and says
    /// nothing about what it is for; it is there because a rejected op once left the model changed.
    pub fn restore_live(&mut self, snapshot: Vec<Handle<Solid>>) {
        self.live_solids = snapshot;
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

    /// A plane's push (private): the exact truth and its f64 cache, index-parallel, in one
    /// motion — so the two cannot come apart.
    ///
    /// ★★★★ **The truth comes first because the arena holds it, and the cache is derived from it
    /// here**: this door calls [`Model::apply_derivation`], so what the producer hands
    /// in survives only in the parts the truth does not decide — the row (`raw`) and with it the
    /// sense — and wholesale where [`Model::derive_surface_cache`] declines.
    /// ⚠ Only the anchor is derived; «the door takes only the truth» is **not** reached
    /// while `cache` is still a parameter.
    ///
    /// ★★★ **Two doors split by kind, rather than one taking both enums.** A single
    /// `push_raw(truth: Surface, cache: nacre_geom::Surface)` could be handed a plane truth
    /// beside a cylinder cache, and **four** sites downstream would need an `unreachable!` to
    /// say the pairing holds. A typed door makes the mismatch unspellable, and those four cite
    /// the door instead of asserting the fact.
    ///
    /// ⚠ **No lock here, and this is why.** The proposition is the signature itself, and a
    /// `compile_fail` doc-test cannot reach a private function to demonstrate it. The public
    /// doors ([`Model::push_plane`], [`Model::push_cylinder`]) already took narrow types; what
    /// was wide was this crate-internal one. Recorded rather than locked.
    ///
    /// ★★★ **The name enters here too**, for the same reason the cache does: this is
    /// the one place a plane reaches the arena, so it is the one place that can derive a cache
    /// from the truth — and the derivation reads the name. Inserting it one line later
    /// (at interning) would leave a window in which a surface exists without the name that
    /// describes it.
    fn push_plane_raw(
        &mut self,
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,
        name: Option<nacre_scalar::PlaneName>,
        cache: nacre_geom::Plane,
    ) -> Handle<Surface> {
        let h = self.surfaces.push(Surface::Plane { points, motion });
        self.surface_cache.push(SurfaceCache {
            realized: nacre_geom::Surface::Plane(cache),
        });
        if let Some(n) = name {
            self.surface_name.insert(h, n);
        }
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "the truth and its cache enter together or not at all"
        );
        self.apply_derivation(h);
        h
    }

    /// **Realize this surface's cache from its truth**, keeping the producer's value where the
    /// truth cannot say ([`Model::derive_surface_cache`] declines).
    ///
    /// ★ Counting happens **first**, against the value the producer stated, so the census `stat`
    /// rows keep meaning "how far the producer's value was from the truth's".
    ///
    /// ⚠ Planes only. [`Model::push_cylinder_raw`] still calls [`Model::measure_derivation`]:
    /// the cylinder arm is derived and **thrown away**, deliberately, because a moved cylinder's
    /// world statement lives in `nacre-ops` and this door cannot reach it.
    fn apply_derivation(&mut self, h: Handle<Surface>) {
        self.measure_derivation(h);
        if let Some(realized) = self.derive_surface_cache(h) {
            self.surface_cache[h.index() as usize] = SurfaceCache { realized };
        }
    }

    /// A cylinder's push (private) — [`Model::push_plane_raw`]'s twin, and the other half of the
    /// reason neither takes a `nacre_geom::Surface`.
    fn push_cylinder_raw(
        &mut self,
        def: CylinderDef,
        motion: Option<Handle<MotionNode>>,
        cache: nacre_geom::Cylinder,
    ) -> Handle<Surface> {
        let h = self.surfaces.push(Surface::Cylinder { def, motion });
        self.surface_cache.push(SurfaceCache {
            realized: nacre_geom::Surface::Cylinder(cache),
        });
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "the truth and its cache enter together or not at all"
        );
        self.measure_derivation(h);
        h
    }

    /// **The f64 cache this surface's truth realizes to** — principle 4's road, for surfaces.
    ///
    /// A vertex has a whole one (`push_vertex_realized` realizes the definition at
    /// birth). A surface has **half** of one: the anchor is derived here, the row and
    /// the sense are whatever the producer handed in.
    ///
    /// ★ **What that half buys** (measured): a producer-stated anchor on a tilted plane
    /// varies by up to **22 ulps** with which face asked for the plane first, and for **20 of 29**
    /// moved planes it does not satisfy the plane's own coefficients. The anchor is
    /// a function of the truth, so neither varies with who asked.
    ///
    /// What it derives:
    /// * **Plane** — the **anchor**, and nothing else: the truth's first point, carried to the
    ///   world and realized. The row (`raw`) and with it the sense are copied from the value the
    ///   producer stated.
    /// * **Cylinder** — `origin` and `radius` descend exactly from [`CylinderDef`]; the axis
    ///   direction and `ref_dir` do not (`Cylinder::from_axis` normalizes both). ⚠ Nothing
    ///   **applies** this arm today — [`Model::push_cylinder_raw`] only measures it.
    ///
    /// ★★★★★ **Why the anchor and not the row** (measured). A canonical row is the
    /// tidier answer, but it moves 32 census result rows, costs an exact-coefficient
    /// road (`WorkingPlane::reconcile` gates on `Plane::spans_exactly`, which a rescaled row
    /// fails), and makes `plane_origin_projection` square the coefficients — which overflows for
    /// 8 of this corpus's planes and spends half the `i128` width budget. The anchor costs none
    /// of that: **0** census rows move, the exact road is untouched, and the arithmetic is one
    /// addition. What it buys is the same thing: the cache becomes **recomputable from the
    /// arena**, so two statements of one plane no longer disagree about where it is anchored.
    ///
    /// ⚠ **This is not «the cache is a function of the geometry».** The anchor is the *first*
    /// point of the *first* pusher's triple; two files whose first pushers state the plane
    /// differently still differ. Inside one model interning makes that unreachable — one name,
    /// one handle, one truth.
    ///
    /// ☑ **The sense cannot come from the truth** (measured): 342 non-seed `Known` planes carry a
    /// cache normal opposing their own point order, so the point order does not name a direction.
    /// Copying the row sidesteps the question — `flipped` and every face's outward spelling are
    /// bit-unchanged by this derivation.
    ///
    /// `None` — the caller keeps the cache it has — for an unnamed plane, a `Wide` name, a
    /// motion that is not a rational translation chain, a `Through` truth, a moved cylinder, or
    /// an overflow. ⚠ The name is still required even though the anchor does not read it: every
    /// number above was measured with that gate on, and widening it is its own measurement.
    ///
    /// ⚠★★★ **One door still bypasses this entirely** — [`Model::push_plane_unregistered`], which
    /// skips the name and therefore the derivation. It is the only remaining way for a surface's
    /// truth and its cache to disagree about where the plane is anchored, and it is `cfg(test)`:
    /// every call site is a fixture that wants one geometric plane held as two handles.
    fn derive_surface_cache(&self, h: Handle<Surface>) -> Option<nacre_geom::Surface> {
        let rat3 = |v: [Rat; 3]| [v[0].to_f64(), v[1].to_f64(), v[2].to_f64()];
        match self.surface(h) {
            Surface::Plane { points, motion } => {
                // The gate, kept verbatim: a plane the model cannot name in the world is one this
                // derivation declines, and every measured number above assumes that population.
                self.world_plane_name(h)?.narrow()?;
                let stated = match self.surface_cache(h) {
                    nacre_geom::Surface::Plane(p) => *p,
                    nacre_geom::Surface::Cylinder(_) => return None,
                };
                let PlanePoints::Known(pts) = points else {
                    // A `Through` truth names vertices, whose meet may not fit `Rat` at all.
                    return None;
                };
                let t = match motion {
                    None => [Rat::from_int(0); 3],
                    Some(leaf) => self.chain_translation(*leaf)?,
                };
                let mut anchor = [0.0f64; 3];
                for (k, a) in anchor.iter_mut().enumerate() {
                    *a = pts[0][k].checked_add(t[k])?.to_f64();
                }
                // The row verbatim — `coefficients()` is `[raw, −raw·origin]`, so its first three
                // are the `raw` the producer built, and copying them keeps the sense with it.
                let c = stated.coefficients();
                let raw = Vector3::from_array([c[0], c[1], c[2]]);
                Plane::from_point_normal(Point3::from_array(anchor), raw)
                    .map(nacre_geom::Surface::Plane)
            }
            Surface::Cylinder { def, motion } => {
                motion.is_none().then_some(())?;
                Cylinder::from_axis(
                    Point3::from_array(rat3(def.origin())),
                    Vector3::from_array(rat3(def.dir())),
                    Vector3::from_array(rat3(def.ref_dir())),
                    def.radius_f64(),
                )
                .map(nacre_geom::Surface::Cylinder)
            }
        }
    }

    /// Count what [`Model::derive_surface_cache`] would do at this push, and change nothing —
    /// the measuring half of [`Model::apply_derivation`], for the doors that do not apply it.
    fn measure_derivation(&self, h: Handle<Surface>) {
        use std::sync::atomic::Ordering::Relaxed;
        match self.derive_surface_cache(h) {
            None => {
                SURFACE_DECLINED.fetch_add(1, Relaxed);
                self.decline_reason(h).fetch_add(1, Relaxed);
            }
            Some(d) => {
                if surface_bits(&d) != surface_bits(self.surface_cache(h)) {
                    SURFACE_DIFFERS.fetch_add(1, Relaxed);
                }
                SURFACE_DERIVED.fetch_add(1, Relaxed);
            }
        };
    }

    /// Which counter a decline belongs to — a **diagnosis of the road already taken**, not a
    /// second derivation. [`Model::derive_surface_cache`] asks [`Model::world_plane_name`],
    /// which folds every cause into one `None`; the populations have to be told apart because
    /// they need different work, and only one of them (`Wide`) is arithmetic at all.
    fn decline_reason(&self, h: Handle<Surface>) -> &'static std::sync::atomic::AtomicU64 {
        match self.surface(h) {
            Surface::Cylinder { .. } => &SURFACE_DECLINED_CYLINDER,
            Surface::Plane { .. } => match self.surface_name.get(&h) {
                None => &SURFACE_DECLINED_UNNAMED,
                Some(n) => match (self.plane_motion(h), n.narrow()) {
                    (_, None) => &SURFACE_DECLINED_WIDE,
                    (Some(leaf), Some(_)) if self.chain_translation(leaf).is_none() => {
                        &SURFACE_DECLINED_MOTION
                    }
                    _ => &SURFACE_DECLINED_ARITH,
                },
            },
        }
    }

    /// Count an interning hit whose incoming cache is **discarded** in favour of the survivor's,
    /// when the two are not the same bits — the product-population measurement of "the same
    /// geometry, described twice, writes the same file".
    fn count_discarded_cache(&self, h: Handle<Surface>, discarded: &nacre_geom::Plane) {
        if surface_bits(self.surface_cache(h))
            != surface_bits(&nacre_geom::Surface::Plane(*discarded))
        {
            SURFACE_DISCARDED_DIFFERING.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }

    /// The surface's exact truth — what it *is*, beside the f64 realization
    /// [`Model::surface_cache`] returns. Total: a surface without a truth is unrepresentable — the
    /// truth **is** the arena entry a handle names, which is what retired `SurfaceDef::Inexact`
    /// and the `UndefinedSurface` violation.
    ///
    /// ★★★★ **Which door answers which question** (measured).
    /// **Classification, comparison and branching ask the truth**; display, tessellation,
    /// bounding and measurement ask [`Model::surface_cache`]. A *kind* — "is this face planar"
    /// — is a fact about what the surface **is**, so it is asked here even though the cache
    /// would answer it correctly (the two can never disagree: `push_plane_raw` and
    /// `push_cylinder_raw` pair them by type). The difference is when a wrong answer becomes
    /// possible: a new surface kind lands in *this* enum first, so a `match` here goes red a
    /// step before one on the cache does.
    #[inline]
    pub fn surface(&self, h: Handle<Surface>) -> &Surface {
        self.surfaces.get(h)
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
            matches!(self.surface(h), Surface::Plane { motion: None, .. }),
            "seed handles must stay the world planes"
        );
        h
    }

    /// **Debug-only: the cross-store guard for an index-parallel cache read.**
    ///
    /// A cache is a plain `Vec` indexed by `h.index()`, so reading it alone accepts a handle
    /// minted by *another* model and silently answers with the wrong cell. [`Store::get`] is
    /// where that guard lives (*"Handle was minted by a different Store"*), so a cache read asks
    /// its store first and throws the answer away.
    ///
    /// ★ Measured: a `Store::get` carries this guard for free; indexing a cache alone drops
    /// it, and a foreign handle reaches `surface_cache`/`vertex_point`/`edge_curve` without a
    /// sound. Release builds pay nothing — `Store::get`'s assertion is `cfg(debug_assertions)`
    /// and so is this call.
    #[inline]
    #[cfg(debug_assertions)]
    fn debug_guard<T>(store: &Store<T>, h: Handle<T>) {
        let _ = store.get(h);
    }

    /// The surface a handle names, realized — the **f64 cache** of [`Model::surface`]'s
    /// answer.
    ///
    /// ★★★★ **What belongs here**: display, tessellation, bounding and
    /// measurement — every question whose answer is a number a rounded copy can carry, plus the
    /// one that compares a cached coordinate against its cached carrier. A *kind* question does
    /// not, even though it would answer correctly; [`Model::surface`] says why.
    ///
    /// Reading is open; **writing is not** — the store is private, so a surface can only
    /// enter through [`Model::push_plane`]/[`Model::push_cylinder`], which state its truth. The
    /// lock:
    ///
    /// ```compile_fail,E0616
    /// let m = nacre_topo::Model::new();
    /// let _ = m.surfaces.len(); // private field — read through `surface_cache`/`surface_count`
    /// ```
    #[inline]
    pub fn surface_cache(&self, h: Handle<Surface>) -> &nacre_geom::Surface {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.surfaces, h);
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "surface cache out of step with the surface store — push surfaces through \
             Model::push_plane / push_cylinder"
        );
        &self.surface_cache[h.index() as usize].realized
    }

    /// How many surfaces the arena holds — live and superseded alike. Handle-validity checks
    /// and the `a_boolean_mints_no_surface` lock read this; nothing iterates the store
    /// (superseded surfaces are still in it — consumers walk the live faces).
    #[inline]
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// The surface an **index** names — how a log's index vocabulary is re-anchored onto the model
    /// being built.
    ///
    /// A handle in an operation log carries only its index across models, so `replay` turns that
    /// index back into a handle of its own arena before applying the operation. `None` past the
    /// end: existence, not legality.
    ///
    /// **This does not open the store.** Reading was already open ([`Model::surface_cache`],
    /// [`Model::surface_count`]); writing still goes only through [`Model::push_plane`] /
    /// [`Model::push_cylinder`], which state the truth. The `compile_fail` lock on
    /// [`Model::surface_cache`] is untouched.
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

    /// The vertex a handle names — the **truth**, which for a vertex is its definition.
    ///
    /// ★ One door per entity per side:
    /// `x(h)` is the arena entry itself, `x_cache(h)` is what was realized from it, and the pieces
    /// underneath are reached by chaining on the returned type rather than by more doors here.
    #[inline]
    pub fn vertex(&self, h: Handle<Vertex>) -> &Vertex {
        self.vertices.get(h)
    }

    /// The edge a handle names — the truth beside [`Model::edge_curve`]'s cache.
    #[inline]
    pub fn edge(&self, h: Handle<Edge>) -> &Edge {
        self.edges.get(h)
    }

    /// The face a handle names. ★ There is no `face_cache`: a face has no realized twin, which is
    /// why [`Model::push_face`] is worth an invariant rather than a cache pairing.
    #[inline]
    pub fn face(&self, h: Handle<Face>) -> &Face {
        self.faces.get(h)
    }

    /// The shell a handle names — cache-less for the same reason as [`Model::face`].
    #[inline]
    pub fn shell(&self, h: Handle<Shell>) -> &Shell {
        self.shells.get(h)
    }

    /// The solid a handle names. Which solids are *live* is a separate question —
    /// [`Model::live_solids`] answers it; the arena keeps superseded ones forever.
    #[inline]
    pub fn solid(&self, h: Handle<Solid>) -> &Solid {
        self.solids.get(h)
    }

    /// How many vertices the arena holds — live and superseded alike, like
    /// [`Model::surface_count`].
    #[inline]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// How many edges the arena holds — live and superseded alike.
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// How many faces the arena holds — live and superseded alike.
    #[inline]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// How many shells the arena holds — live and superseded alike.
    #[inline]
    pub fn shell_count(&self) -> usize {
        self.shells.len()
    }

    /// How many solids the arena holds — live and superseded alike. Most of them are dead:
    /// measured, a box moved 120 times leaves 120 superseded solids behind one live one.
    #[inline]
    pub fn solid_count(&self) -> usize {
        self.solids.len()
    }

    /// The edge an **index** names — [`Model::surface_handle_at`]'s twin for edges.
    ///
    /// ★★★★ **This family is also how the whole arena is walked, and that is deliberate.** There
    /// is no `iter()` door and there will not be one: exposing one invites treating the arena as
    /// the model, when most of what it holds is superseded (measured: 99% of the vertices after a
    /// box is moved 120 times). A consumer that wants the live model walks `live_solids` and the
    /// faces under it.
    ///
    /// But a consumer that wants the **arena** — `validate`'s reference-integrity and vertex-def
    /// checks — genuinely needs every cell, dead ones included: a torn page is torn whether or not
    /// anything still points at it. Those checks walk `0..x_count()` through these doors. Filtering
    /// them to the reachable set would make four planted-corruption tests pass while measuring
    /// nothing, because each plants its corruption on an unreachable cell.
    ///
    /// ☑ Walking by index is not a weaker `iter()`: `nacre-store` locks `handle_at` and `iter`
    /// to the same order, by unit test and by proptest, under the note that the door "grants no
    /// new power".
    #[inline]
    #[must_use]
    pub fn edge_handle_at(&self, index: u32) -> Option<Handle<Edge>> {
        self.edges.handle_at(index)
    }

    /// The face an index names — see [`Model::edge_handle_at`] for why this family exists.
    #[inline]
    #[must_use]
    pub fn face_handle_at(&self, index: u32) -> Option<Handle<Face>> {
        self.faces.handle_at(index)
    }

    /// The shell an index names — see [`Model::edge_handle_at`].
    #[inline]
    #[must_use]
    pub fn shell_handle_at(&self, index: u32) -> Option<Handle<Shell>> {
        self.shells.handle_at(index)
    }

    /// The solid an index names — see [`Model::edge_handle_at`].
    #[inline]
    #[must_use]
    pub fn solid_handle_at(&self, index: u32) -> Option<Handle<Solid>> {
        self.solids.handle_at(index)
    }

    /// The motion node a handle names. Same seal as [`Model::surface_cache`]: writing goes through
    /// [`Model::push_motion`] (interned), reading through here.
    #[inline]
    pub fn motion(&self, h: Handle<MotionNode>) -> &MotionNode {
        self.motions.get(h)
    }

    /// Whether more than `n` recorded motions stand between `leaf` and the world — the length of
    /// the chain a replay would walk (a `Frame` node counts once; its own short expansion is not a
    /// history). Walks at most `n + 1` nodes, so asking about a 4,000-deep history costs `n`.
    pub fn motion_deeper_than(&self, leaf: Handle<MotionNode>, n: usize) -> bool {
        let mut depth = 0;
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            depth += 1;
            if depth > n {
                return true;
            }
            cur = self.motion(h).parent;
        }
        false
    }

    /// Whether every node of `leaf`'s recorded chain fixes the plane `coeffs` **as a set** —
    /// the consumer-side twin of the producer's `Isometry::fixes_plane`, on the same scalar
    /// atoms (one rule per motion kind, two thin composers).
    ///
    /// What it licenses: a world-stated plane among motion-carrying carriers is usable in
    /// the carriers' *pre-motion* frame **iff** the chain fixes it — then its world equation
    /// is the same equation there, and the corner solves as if all three shared the chain.
    /// The invariant-plane restatement mints exactly this shape (a turned block's corner =
    /// restated cap × two chained walls), and this check is what keeps the licence honest:
    /// a frame-hosted datum's cap hits the `Frame` arm and the caller stays declined.
    ///
    /// ★ **It lives here, below `nacre-ops`, because [`Model::vertex_meet`] needs it.** Without
    /// the rescue that door reads a turned solid's corner as straddling frames — every carrier
    /// but the fixed cap moved — and the datum road it gates loses its named form. Keeping the
    /// rule in `nacre-ops` would mean a second walk down here, and this rule already cost the
    /// kernel a regression by existing in more than one spelling.
    ///
    /// Conservative by construction — `Frame` nodes are never fixed (their basis is
    /// irrational), and the scalar atoms answer `false` on overflow.
    pub fn chain_fixes_plane(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_scalar::Rat; 4],
    ) -> bool {
        self.chain_fixes(leaf, coeffs, false)
    }

    /// **The single rational translation `leaf`'s whole chain amounts to**, or `None` when the
    /// chain is anything else — the third question of this family, and the one that lets a
    /// *moved* statement be restated in the world exactly.
    ///
    /// A rational translation maps a rational statement to a rational statement, so a chain made
    /// only of [`Motion::Translate`] nodes loses nothing: a plane's `d` shifts by `−n·t`
    /// ([`nacre_scalar::Isometry::plane_coeffs`]), a cylinder's origin by `+t`. What such a move
    /// *does* lose is the exactness of the `f64` **cache** — which is why the producer still
    /// records the node (`transform`'s `carry_of`, and `nacre-ops`' reuse road reads a
    /// world-stated carrier's coordinate as the statement itself). So this answers a question
    /// about *descriptions*, for the per-operation mirrors that carry them; it does not license
    /// dropping the history.
    ///
    /// ★ **The parent chain is walked raw, so a [`Motion::Frame`] node refuses outright.** A
    /// frame's expansion can come back as translations, and a statement written *under* a frame
    /// is in that frame's coordinates — folding those as world translations is a different
    /// question. Translations commute and compose by addition, so no order is implied here.
    /// `None` on overflow too (checked throughout).
    pub fn chain_translation(&self, leaf: Handle<MotionNode>) -> Option<[nacre_scalar::Rat; 3]> {
        let mut total = [nacre_scalar::Rat::from_int(0); 3];
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            let node = self.motion(h);
            let Motion::Translate { offset } = node.motion else {
                return None;
            };
            for (o, t) in total.iter_mut().zip(offset) {
                *o = o.checked_add(t)?;
            }
            cur = node.parent;
        }
        Some(total)
    }

    /// **A plane's canonical name in the world**, whatever frame its truth is written in — the
    /// door between "how this surface got here" and "where it is".
    ///
    /// Unmoved: the name itself, which is already world. Moved by a chain that folds to a
    /// rational translation ([`Model::chain_translation`]): that name carried out exactly
    /// (`d' = d − n·t`, [`nacre_scalar::Isometry::plane_coeffs`], canonicalized). `None` for
    /// anything else — a rotation, a frame node, a **moved** [`nacre_scalar::PlaneName::Wide`]
    /// (no narrow vessel for the transport to take), or an overflow — and the caller declines
    /// rather than guessing. An *unmoved* `Wide` name comes back verbatim: it is already world,
    /// and refusing it would switch off a capability that is locked elsewhere.
    ///
    /// **Planes only.** A cylinder surface has no entry in `surface_name`, so it answers `None`;
    /// its own world statement is `nacre-ops`' `world_cylinder_def`.
    pub fn world_plane_name(&self, surf: Handle<Surface>) -> Option<nacre_scalar::PlaneName> {
        let name = self.surface_name.get(&surf)?;
        match self.plane_motion(surf) {
            None => Some(name.clone()),
            Some(leaf) => {
                let t = self.chain_translation(leaf)?;
                let c = name.narrow()?;
                Some(nacre_scalar::PlaneName::Narrow(
                    nacre_scalar::Isometry::translation(t).plane_coeffs(*c)?,
                ))
            }
        }
    }

    /// The strict twin: every node carries the plane's **coefficient row verbatim**, not merely
    /// the set. The one place they part is a mirror whose plane is the carrier itself (`n ∥ axis`,
    /// on-plane): the set maps to itself but the row comes back negated — harmless to a Cramer
    /// solve (both determinants negate, the ratio stands), fatal to a determinant *sign* read.
    /// So [`Model::chain_fixes_plane`] licenses solving, and this licenses judging-table
    /// descriptions (`nacre_ops`' mirror re-chaining).
    pub fn chain_preserves_plane_row(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_scalar::Rat; 4],
    ) -> bool {
        self.chain_fixes(leaf, coeffs, true)
    }

    fn chain_fixes(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_scalar::Rat; 4],
        rows_verbatim: bool,
    ) -> bool {
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            let n: &MotionNode = self.motion(h);
            let fixed = match &n.motion {
                Motion::Rotate { axis, .. } => {
                    nacre_scalar::axis_rotation_fixes_plane(*axis, coeffs)
                }
                Motion::Translate { offset } => {
                    nacre_scalar::translation_fixes_plane(offset, coeffs)
                }
                Motion::Mirror { axis, offset } => {
                    if rows_verbatim {
                        coeffs[match axis {
                            nacre_scalar::Axis::X => 0,
                            nacre_scalar::Axis::Y => 1,
                            nacre_scalar::Axis::Z => 2,
                        }] == nacre_scalar::Rat::from_int(0)
                    } else {
                        nacre_scalar::mirror_fixes_plane(*axis, *offset, coeffs)
                    }
                }
                Motion::Frame { .. } => false,
            };
            if !fixed {
                return false;
            }
            cur = n.parent;
        }
        true
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
        // for collinear points (which no production `PlaneDef` can supply); the vessel
        // (`PlaneName::Narrow | Wide`) always holds the answer, so every plane interns, wide
        // ones included. [`WIDE_PLANES`] counts the names that took the wide vessel.
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
        let key = name.clone().map(|n| (n, motion));
        if let Some(k) = &key {
            if let Some(&h) = self.surface_ids.get(k) {
                if (h.index() as usize) < 3 {
                    SEEDED_HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                // ★ The cache the caller built is **dropped here** — the survivor's stands. That
                // is what makes a plane's realization depend on which face asked first, and
                // [`Model::count_discarded_cache`] is the only measurement of how often the two
                // actually differ.
                self.count_discarded_cache(h, &cache);
                // Same plane, already issued. The canonical form says nothing about direction, so
                // report whether the survivor points the other way and let the caller spell its
                // outward the other way round.
                return (h, self.flipped_against(h, &cache));
            }
        }
        // One clone per push — the name is derived once here, never on a judging loop.
        let h = self.push_plane_raw(points, motion, name, cache);
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
        match self.surface_cache(h) {
            nacre_geom::Surface::Plane(p) => p.normal().dot(cache.normal()) < 0.0,
            nacre_geom::Surface::Cylinder(_) => false,
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
    /// ★★ **A statement the name key cannot hold still interns — by the statement itself.**
    /// A mixed-frame datum's exact world coefficients are irrational, so
    /// [`Model::plane_name_through`] answers `None`; such a plane takes the second key
    /// (`surface_through_ids`) — the sorted triple and the motion. That is *statement*
    /// identity: the same three vertices under the same motion are one handle, and geometric
    /// identity across different statements is the predicates' to answer per question (rule 6's
    /// own qualification — *"without interning, the predicate answers every time"*). This is
    /// **not** the
    /// record-less population: the truth (handles + motion) is complete; what does
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
        let h = self.push_plane_raw(PlanePoints::Through(vertices), motion, None, cache);
        self.surface_through_ids.insert((vertices, motion), h);
        (h, false)
    }

    /// **The name a `Through` statement derives**, and the one place that derivation lives — the
    /// producer's check and [`Model::push_plane_through`] read the same answer, so "we rejected
    /// what we could not name" is structural rather than two functions agreeing by habit.
    ///
    /// `None` when any vertex is not a three-plane point, when the carriers do not share one
    /// motion (no frame holds a rational coordinate then), or when the three points are
    /// collinear. ★ A meet too wide for `Rat` is **not** on the list: the
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
    /// `None` on any of: a vertex that is not a three-plane point, a vertex [`Model::vertex_meet`]
    /// cannot place in one frame, a carrier with no recorded name, or three vertices whose frames
    /// disagree.
    pub fn through_meets(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<[nacre_scalar::MeetPoint; 3]> {
        let mut pts: [Option<nacre_scalar::MeetPoint>; 3] = [None, None, None];
        let mut frame = None;
        for (i, vh) in vertices.iter().enumerate() {
            let (p, mine) = self.vertex_meet(*vh)?;
            match frame {
                None if i == 0 => frame = Some(mine),
                f if f == Some(mine) => {}
                _ => return None, // the three vertices do not share one frame
            }
            pts[i] = Some(p);
        }
        Some([pts[0].take()?, pts[1].take()?, pts[2].take()?])
    }

    /// **One vertex's exact meeting point, and the frame it is stated in** — the per-vertex half
    /// of [`Model::through_meets`], which is its caller for the three-vertex case.
    ///
    /// The frame comes back beside the point because a coordinate means nothing without it: `None`
    /// is the world, `Some(node)` a pre-motion frame. A consumer comparing this against anything
    /// stated in world coordinates must **require `None`** — the three-vertex caller above only
    /// needs the three to *agree*, which is a weaker demand and would silently mix frames if it
    /// were copied.
    ///
    /// ★★★ **"Carriers in one frame" is not "carriers with one `motion` field".** The sentence
    /// that used to stand here — *`None` when the carriers do not share one motion, since then no
    /// frame holds a rational coordinate at all* — became false with the invariant-plane
    /// restatement, and reading it as still true cost every turned solid's corners their named
    /// datum road for 169 commits. A motion that **fixes** a plane restates nothing, so that plane
    /// stays world-stated beside carriers that moved; its world equation *is* its equation in
    /// their pre-motion frame, and [`Model::chain_fixes_plane`] is what proves it.
    ///
    /// `None` on any of: a vertex that is not a three-plane point; carriers carrying **two**
    /// motion histories; a world-stated carrier the shared chain does not fix (or whose name is
    /// `Wide`, since the licence reads narrow coefficients — conservative, and recorded); a
    /// carrier with no recorded name; or three carriers that meet in no point.
    pub fn vertex_meet(
        &self,
        v: Handle<Vertex>,
    ) -> Option<(nacre_scalar::MeetPoint, Option<Handle<MotionNode>>)> {
        self.vertex_meet_of(self.vertices.get(v))
    }

    /// [`Model::vertex_meet`] on a definition that has not been pushed yet — what an operation
    /// asks before it states a vertex, so the cache it pushes is already the realization.
    pub fn vertex_meet_of(
        &self,
        def: &Vertex,
    ) -> Option<(nacre_scalar::MeetPoint, Option<Handle<MotionNode>>)> {
        let tri = match *def {
            Vertex::ThreePlane(tri) => tri,
            // OnSeam pins a curve, not a point; a Pierce *is* a point but its coordinates
            // are quadratic-irrational — neither has the rational meet a datum statement
            // needs, so both decline here (honest, and spelled per variant so the next
            // variant is a compile error, not a silent fall-through).
            Vertex::OnSeam(_) | Vertex::Pierce { .. } => return None,
        };
        // ★★★ **The third door — solve in the world**. The two below want *one*
        // frame: the world (nothing moved) or one shared chain. A second-generation array breaks
        // both — an array fused in x and then moved in y puts carriers with chains `T2` and
        // `T1·T2` on one corner — and yet every one of those planes states the world exactly,
        // because a rational translation carries a rational name ([`Model::world_plane_name`]).
        // So when the shared-frame roads decline, the triple is solved from the **world names**
        // and the answer carries no leaf: the caller must not replay anything.
        //
        // ★ **It is a fallback, deliberately.** Running it first would re-spell points the two
        // doors already answer, and those spellings are what the corpus is pinned on. Reached
        // only where today's answer is `None`, it can open a population and cannot move one.
        //
        // ★ It transports the **name**, never the point: `world_plane_name` moves each carrier's
        // equation into the world, and `three_planes_big` then meets three world planes. Nothing
        // here realizes a coordinate and shifts it.
        let world_road = || -> Option<nacre_scalar::MeetPoint> {
            let (a, b, c) = (
                self.world_plane_name(tri[0])?,
                self.world_plane_name(tri[1])?,
                self.world_plane_name(tri[2])?,
            );
            nacre_scalar::three_planes_big([&a, &b, &c])
        };
        let motions = tri.map(|h| self.plane_motion(h));
        // One shared leaf among the **moved** carriers; no moved carrier means the world.
        let mut leaf = None;
        for m in motions.iter().flatten() {
            match leaf {
                None => leaf = Some(*m),
                Some(l) if l == *m => {}
                // Two histories — no shared frame, so ask whether both state the world.
                Some(_) => return world_road().map(|p| (p, None)),
            }
        }
        let names = tri.map(|h| self.surface_name.get(&h));
        let [Some(a), Some(b), Some(c)] = names else {
            return None;
        };
        // ★★ **A world-stated carrier among chained ones is not a straddle if the chain fixes
        // it.** The invariant-plane restatement mints exactly that shape — a turned block's
        // corner is a restated cap × two chained walls — and a fixed plane's world equation *is*
        // its equation in the pre-motion frame, so the corner solves as if all three shared the
        // chain. Reading the mismatch as a straddle is what cost every turned solid's corners
        // their named datum road; `chain_fixes_plane` is the licence that tells them apart.
        //
        // ★ `narrow()` is asked **only of the carriers being licensed**, never of the other two:
        // this function's answer is a `MeetPoint` of any width, and refusing a `Wide` carrier
        // here would quietly switch off the capability `a_datum_through_wide_meets_keeps_its_name`
        // locks. A `Wide` *fixed* carrier declines the licence — the same conservatism the atoms
        // have, recorded rather than papered over.
        if let Some(leaf) = leaf {
            for (n, m) in [a, b, c].iter().zip(&motions) {
                if m.is_none() && !n.narrow().is_some_and(|c| self.chain_fixes_plane(leaf, c)) {
                    // The carriers straddle frames — the world road is the remaining question.
                    return world_road().map(|p| (p, None));
                }
            }
        }
        Some((nacre_scalar::three_planes_big([a, b, c])?, leaf))
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
        match self.surface(h) {
            Surface::Plane { motion, .. } | Surface::Cylinder { motion, .. } => *motion,
        }
    }

    /// Push a **cylinder** — the lateral surface, stating its exact truth, with
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
        let h = self.push_cylinder_raw(def, motion, cache);
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
        self.push_plane_raw(PlanePoints::Known(points), None, None, cache)
    }

    /// A new shell whose faces are copies of `src`'s with their outward normals
    /// flipped inward: every loop's winding is reversed and every face
    /// `orientation` is toggled. Pushes fresh [`Face`] cells and a fresh
    /// [`Shell`], but **reuses** `src`'s surfaces, edges, curves, and vertices —
    /// which stay valid handles after the source solid is superseded
    /// (append-only). This is the building block for a cavity (void) shell: the
    /// boundary of a solid whose interior becomes empty space (M5
    /// containment). The reversed winding keeps each shared edge used with
    /// opposed half-edges (a valid 2-manifold), and the toggled orientation makes
    /// the outward normal point into the void.
    #[cfg(any(test, feature = "test-util"))]
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

    /// A vertex's cache — **the one road to a coordinate from a vertex**, read from the
    /// index-parallel store [`Model::push_vertex`] fills. The coordinate piece is
    /// [`Model::vertex_point`] and the proven bound is [`PointCache::bound`]; read the whole when
    /// the variant itself — what the cache *knows* — is the question.
    #[inline]
    pub fn vertex_cache(&self, vh: Handle<Vertex>) -> &PointCache {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.vertices, vh);
        debug_assert_eq!(
            self.vertex_cache.len(),
            self.vertices.len(),
            "vertex cache out of step with the vertex store — push vertices through Model::push_vertex"
        );
        &self.vertex_cache[vh.index() as usize]
    }

    /// A vertex's realized coordinate — [`Model::vertex_cache`]'s coordinate piece.
    #[inline]
    pub fn vertex_point(&self, vh: Handle<Vertex>) -> Point3 {
        self.vertex_cache(vh).coord()
    }

    /// Push a vertex: its definition (the truth) plus its cache. **The road that builds a vertex,
    /// and the only one** — nothing else may add to the store.
    ///
    /// ★ There is no `rebuild_vertex_cache` that *replaces* what is here: an operation realizes the
    /// definition *before* it pushes (`nacre_ops`'s push funnel), so what arrives is already the
    /// realization's memo ([`PointCache::Bounded`]) or a named reason it is not
    /// ([`PointCache::Ceiling`], [`PointCache::Unrealized`]). The old warrant for the absence
    /// (*"238 of 1,992 differ from a naive re-solve"*) compared against a naive re-solve; the
    /// exact-rounding realization is not one, and the census says so vertex by vertex.
    ///
    /// ⚠ **A second writer does exist, and it only ever raises**: `nacre_ops::refine_vertex_cache`
    /// lifts a `Ceiling` to `Bounded` by paying for the realization this road would not. It cannot
    /// reach the other variants and cannot move a coordinate anywhere but closer to the truth, so
    /// what stays true of the cache behind this door is *append-only in accuracy*, not in bytes.
    pub fn push_vertex(&mut self, def: Vertex, cache: PointCache) -> Handle<Vertex> {
        match def {
            Vertex::ThreePlane([a, b, c]) => debug_assert!(
                a != b && b != c && a != c,
                "a three-plane definition needs three distinct planes"
            ),
            Vertex::OnSeam([a, b]) => {
                debug_assert!(a != b, "a seam vertex needs two distinct carriers")
            }
            Vertex::Pierce {
                planes: [a, b],
                cylinder,
                ..
            } => debug_assert!(
                a.index() < b.index() && cylinder != a && cylinder != b,
                "a pierce definition needs two sorted distinct planes and a distinct cylinder"
            ),
        }
        let h = self.vertices.push(def);
        self.vertex_cache.push(cache);
        h
    }

    /// **Raise a [`PointCache::Ceiling`] to [`PointCache::Bounded`]** — the second writer of the
    /// vertex cache, and the only thing it can do is make a coordinate more accurate.
    ///
    /// ★ Deliberately not `set_vertex_cache`. It cannot reach the other two variants: an
    /// `Unrealized` coordinate has no realization behind it (raising it would be inventing one) and
    /// a `Bounded` one is already the realization. The debug assert is the type saying so out loud
    /// — the caller that pays for the realization is `nacre_ops::refine_vertex_cache`, and its own
    /// contract is that it looks at nothing else.
    ///
    /// ⚠ **Anything derived from this coordinate is now stale**, starting with the edge curves the
    /// endpoints decide ([`Model::rebuild_edge_cache`]). The paying caller is responsible for that,
    /// because it is the one that knows whether it moved anything at all.
    pub fn refine_vertex_cache(&mut self, vh: Handle<Vertex>, coord: Point3, bound: [Mag; 3]) {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.vertices, vh);
        let slot = &mut self.vertex_cache[vh.index() as usize];
        debug_assert!(
            matches!(slot, PointCache::Ceiling { .. }),
            "only a Ceiling is raised, not {slot:?}"
        );
        *slot = PointCache::Bounded { coord, bound };
    }

    /// An edge's curve — **the one road to a curve from an edge**, read from the
    /// index-parallel cache [`Model::push_edge`] fills.
    #[inline]
    pub fn edge_curve(&self, e: Handle<Edge>) -> &Curve {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.edges, e);
        debug_assert_eq!(
            self.edge_cache.len(),
            self.edges.len(),
            "edge cache out of step with the edge store — push edges through Model::push_edge"
        );
        &self.edge_cache[e.index() as usize].curve
    }

    /// Push an edge: canonicalize the carrier pair, derive its curve, fill the cache — the one
    /// write road. `None` when the curve does not derive, which for the line arms means
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

    /// Push a face — **the door that carries the face invariant**.
    ///
    /// ★ Unlike [`Model::push_vertex`] this pairs no cache with the truth, because a face has
    /// none (`surface_cache`/`edge_cache`/`vertex_cache` exist; a face cache does not). What this
    /// door is worth is therefore the **invariant**, not cache-sync: a face whose loops do not
    /// close is a face no consumer can walk, and until now nothing said so at the moment it was
    /// built.
    ///
    /// ☑ **Measured before it shipped** — every face the production road builds satisfies this,
    /// live or superseded, across the boolean / twice-cut / cylinder / tilted-frame fixtures. The
    /// one fixture that did not was a deliberately malformed «franken» face in a `transform` test,
    /// and it was reshaped rather than exempted, so the invariant holds with no hole behind it.
    pub fn push_face(&mut self, face: Face) -> Handle<Face> {
        debug_assert!(
            (face.surface.index() as usize) < self.surfaces.len(),
            "a face names a surface the arena does not hold"
        );
        debug_assert!(
            !face.outer.half_edges.is_empty(),
            "a face's outer loop has no half-edges"
        );
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            debug_assert!(
                !lp.half_edges.is_empty(),
                "a face's inner loop has no half-edges"
            );
            let n = lp.half_edges.len();
            for i in 0..n {
                let (he, nx) = (&lp.half_edges[i], &lp.half_edges[(i + 1) % n]);
                let [a, b] = self.edges.get(he.edge).vertices;
                let end = if he.forward { b } else { a };
                let [c, d] = self.edges.get(nx.edge).vertices;
                let start = if nx.forward { c } else { d };
                debug_assert_eq!(end, start, "a face loop does not close at half-edge {i}");
            }
        }
        self.faces.push(face)
    }

    /// Push a face **without its invariant** — test-only.
    ///
    /// `validate`'s own tests have to plant models that are wrong: a loop that does not close, a
    /// face duplicated onto another's loop, a winding deliberately reversed. Those cannot go
    /// through [`Model::push_face`], whose assert dereferences the very edge being dangled. This
    /// is the same exception [`Model::push_plane_unregistered`] is, and the only one that is ever
    /// justified: something the product cannot express, kept for the tests that must express it.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_face_unchecked(&mut self, face: Face) -> Handle<Face> {
        self.faces.push(face)
    }

    /// Push a shell **without its invariant** — test-only, see [`Model::push_face_unchecked`].
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_shell_unchecked(&mut self, shell: Shell) -> Handle<Shell> {
        self.shells.push(shell)
    }

    /// Push a solid **without making it live** — test-only.
    ///
    /// ⚠ Not [`Model::push_solid`], and the difference is the whole point: that door marks the
    /// solid live, while a planted dangling reference has to stay **unreachable**, because what it
    /// proves is that the arena-wide checks see cells nothing points at.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_solid_unlisted(&mut self, solid: Solid) -> Handle<Solid> {
        self.solids.push(solid)
    }

    /// Push a shell — [`Model::push_face`]'s twin, and a cache-less door for the same reason.
    ///
    /// ☑ **Measured before it shipped**: across the same fixtures no shell the production road
    /// builds is empty or names a face outside the arena, live or superseded (36 shells, 0 and 0).
    pub fn push_shell(&mut self, shell: Shell) -> Handle<Shell> {
        debug_assert!(
            !shell.faces.is_empty(),
            "a shell with no faces bounds nothing"
        );
        debug_assert!(
            shell
                .faces
                .iter()
                .all(|f| (f.index() as usize) < self.faces.len()),
            "a shell names a face the arena does not hold"
        );
        self.shells.push(shell)
    }

    /// Discard every edge-curve cache and derive it afresh — the «cache, not truth» warrant:
    /// nothing is lost, because nothing there was truth.
    /// ⚠★★★ **Only the reachable edges are re-derived, and a superseded one keeps what it has.**
    /// The arena is append-only, so most of what is in it is dead: measured, a boolean corner has
    /// 24 dead edges of 48, a twice-cut one 60 of 108, a thrice-moved box 36 of 48. Re-deriving
    /// those costs the work twice over and — once a caller can *move* a coordinate
    /// ([`Model::refine_vertex_cache`]) — risks a dead cell's endpoints becoming coincident, where
    /// the derivation answers `None` and this would die on the `expect`.
    ///
    /// ☑ That `None` does not happen today: measured over the same fixtures, **zero** stored edges
    /// fail to derive, dead or live. The filter is not a workaround for a live failure — it is what
    /// makes the failure structurally unreachable, because a superseded edge is never derived again.
    pub fn rebuild_edge_cache(&mut self) {
        let reach = self.reachable();
        self.edge_cache = self
            .edges
            .iter()
            .map(|(eh, e)| match reach.edges.contains(&eh) {
                true => EdgeCache {
                    curve: self
                        .derive_edge_curve(e.surfaces, e.vertices)
                        .expect("every live edge derives its curve"),
                },
                false => self.edge_cache[eh.index() as usize].clone(),
            })
            .collect();
    }

    /// The vertex a half-edge starts at: its edge's `vertices[0]` when the use runs
    /// forward, `vertices[1]` when it runs back.
    ///
    /// A traversal accessor, not an analysis — the same kind of thing as
    /// [`Model::reachable`], and the reason it lives here: walking a loop's
    /// corners is the first thing every consumer above does, and it was written
    /// twice (with two different failure policies) before this existed.
    ///
    /// Total: every edge has both endpoints by type, and a closed rim
    /// states that by repeating its seam vertex — `[v, v]`, so the start is `v` either way.
    /// A standalone full circle with no seam would have no start, but the type does not
    /// express that form: it is a wireframe/open-shell element and a v1 non-goal. Nothing
    /// guards it because nothing can build it.
    #[inline]
    pub fn he_start(&self, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = self.edges.get(he.edge).vertices;
        if he.forward { a } else { b }
    }

    /// **An edge's curve, derived from what the model already holds** — the edge's
    /// truth/cache split: the carriers and the endpoints decide the curve, so the stored one
    /// is a cache that can be discarded and regenerated.
    ///
    /// Dispatch by carrier type:
    /// * **Plane × Plane** (and the self-adjacent cylinder **seam**): the line through the two
    ///   endpoint coordinates — the very expression a producer would build the stored
    ///   curve with, so the derivation is bit-identical, and an endpoint pair that coincides is the
    ///   `None` (a degenerate line — the check lives in this arm only).
    /// * **Plane × Cylinder** (a rim, or an arc of one): the circle centred where the cylinder's
    ///   axis crosses the cap plane, with the **cylinder's** frame (`axis direction`, `ref_dir`,
    ///   `radius`) — the same parameters `add_cylinder` builds the stored rims from, so
    ///   tessellation's `θ` parameterization is preserved. The endpoints are not read: a full rim
    ///   is a closed edge (`[v, v]`), which is not a degeneracy.
    ///
    ///   ★★ **On a circle carrier, the vertex *order* says which arc**: two distinct
    ///   endpoints cut a circle into two pieces the endpoints alone cannot tell apart, so
    ///   `[A, B]` means the piece from A to B **counter-clockwise about the axis direction**,
    ///   and the two complementary arcs between one vertex pair are the two orders. Producers
    ///   uphold this (`boolean`'s edge welding keys arcs in CCW order); the curve stored here is
    ///   the whole circle either way. Both consumers read it the same way: tessellation's
    ///   `sample_edge` walks `θ(v0) → θ(v0) + Δθ` with `Circle::angle_of` as the one spelling of
    ///   θ, and `validate`'s `loop_winding` adds each arc's circular segment with Δθ from the
    ///   same order.
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
        // ★★ **Both kinds asked of the cache here, deliberately**. Every arm
        // reads cache *values* out of the very binding it matched — `p.normal()`, `c.axis()`,
        // `c.radius()` — so dispatching on the truth would double the lookups and leave two
        // matches whose agreement no reader could check. The rule sends *kind questions* to the
        // truth; this match's answer is a curve, and the kinds only choose how to derive it.
        match (
            self.surface_cache(surfaces[0]),
            self.surface_cache(surfaces[1]),
        ) {
            (nacre_geom::Surface::Plane(_), nacre_geom::Surface::Plane(_)) => endpoints_line(),
            (nacre_geom::Surface::Cylinder(_), nacre_geom::Surface::Cylinder(_))
                if surfaces[0] == surfaces[1] =>
            {
                endpoints_line() // the seam — a parameterization joint, straight along the axis
            }
            (nacre_geom::Surface::Plane(p), nacre_geom::Surface::Cylinder(c))
            | (nacre_geom::Surface::Cylinder(c), nacre_geom::Surface::Plane(p)) => {
                let axis = c.axis();
                // A plane **parallel** to the axis meets the lateral along rulings — straight,
                // so the endpoints decide, exactly like the seam arm above (the rulings
                // ladder). Same scale convention as the ⊥ assertion below, so the band between
                // the two tests is symmetric and only a genuinely tilted plane (an ellipse)
                // falls through to it.
                {
                    let n = p.normal();
                    let d = axis.direction();
                    if n.dot(d).powi(2) <= 1e-18 * n.norm_squared() * d.norm_squared() {
                        return endpoints_line();
                    }
                }
                debug_assert!(
                    {
                        let n = p.normal();
                        let d = axis.direction();
                        n.cross(d).norm_squared() <= 1e-18 * n.norm_squared()
                    },
                    "a tilted plane over a cylinder crosses in an ellipse — M6-3, no producer yet"
                );
                let center = nacre_geom::intersect::line_plane(&axis, p)?;
                Some(Curve::Circle(Circle::from_center_normal(
                    center,
                    axis.direction(),
                    c.ref_dir(),
                    c.radius(),
                )?))
            }
            (nacre_geom::Surface::Cylinder(_), nacre_geom::Surface::Cylinder(_)) => None, // two distinct cylinders: M6
        }
    }

    /// The handles reachable from the live solids — the live model.
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
                Vertex::ThreePlane([surf[cap].0, surf[along_y].0, surf[along_x].0]),
                PointCache::Unrealized { coord: corners[i] },
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

        // The exact truth: a direct lift of the caller's statement — origin and the raw,
        // unnormalized axis (normalizing would destroy the exact form; the `normal_def`
        // precedent). `ref_dir` replicates `any_perpendicular`'s own rule in rationals: cross
        // the axis with the basis axis of its smallest |component| (ties X→Y→Z).
        //
        // ★ **The basis choice reads `d` — the very components `any_perpendicular` reads — so
        // agreement is structural, not order-theoretic.** Comparing the *raw*
        // components on the argument "one positive scale preserves |·| order" is a real-number
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
                nacre_scalar::BigRat::square_of(lift(radius)),
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
    #[cfg(any(test, feature = "test-util"))]
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
        // The truth carries the squared radius; a stated `radius` is squared once, here — wide,
        // so no radius is refused for the width of its square.
        let def = CylinderDef::new(base, axis, ref_dir, nacre_scalar::BigRat::square_of(radius))
            .ok_or(CylinderError::Degenerate)?;

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
    /// needs a seam vertex, so the rims repeat it as `[v, v]` — an endpointless edge (the
    /// standalone full circle) is a form this type does not express.
    ///
    /// Returns the solid beside its three faces in push order (lateral, bottom cap, top cap).
    /// `None` if an edge cannot derive its curve from the carriers it states — each caller says
    /// what that means for it. Does **not** rebuild adjacency.
    #[cfg(any(test, feature = "test-util"))]
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

        // ★ **Surfaces before edges**: an edge states its two carriers, so the lateral
        // cylinder and both cap planes must exist first. Separate arenas — the interleaving
        // moves no handle; the surfaces' order among themselves (lateral → bottom cap → top
        // cap) is what matters (the `add_cuboid` precedent).
        let lateral_surface = self.push_cylinder(lateral_cache, def, motion);
        // ★★★ **A cap's plane may already be in the model, facing the other way** — and it is the
        // *rule*, not the exception, once cylinders come from operations: the bottom cap lies on
        // the very plane the sketch frame names. `push_plane` interns by the plane's canonical
        // name, which has no direction, so it hands back "the surface exists, its cache points
        // the other way" and the face must record `Orientation::flipped()` — the `add_cuboid`
        // spelling. Dropping the bit gave a cylinder standing on a face **two upward caps**
        // (measured: bottom cap outward `+Z` on a top-face frame), which is not a solid at all.
        let (bottom_cap_surface, bottom_flipped) =
            self.push_plane(bottom_plane, bottom_points, motion);
        let (top_cap_surface, top_flipped) = self.push_plane(top_plane, top_points, motion);
        let facing = |flipped: bool| {
            if flipped {
                Orientation::Forward.flipped()
            } else {
                Orientation::Forward
            }
        };

        // A seam vertex lies on two surfaces only — the rim circle's `θ = 0` point. `OnSeam`
        // states exactly that, and the designation is complete: the cylinder's
        // truth carries `ref_dir`, so the pair means "rim ∩ the `+ref_dir` ray" — one point
        // (see `Vertex::OnSeam`; regenerating the cached coordinate is not done here).
        let v_bot = self.push_vertex(
            Vertex::OnSeam([lateral_surface, bottom_cap_surface]),
            PointCache::Unrealized { coord: p_bot },
        );
        let v_top = self.push_vertex(
            Vertex::OnSeam([lateral_surface, top_cap_surface]),
            PointCache::Unrealized { coord: p_top },
        );

        // Rims are full circles seamed at their vertex (start == end); the seam is
        // a straight edge joining the two rim seam points.
        let bottom = self.push_edge([lateral_surface, bottom_cap_surface], [v_bot, v_bot])?;
        let top = self.push_edge([lateral_surface, top_cap_surface], [v_top, v_top])?;
        // Self-adjacent: a seam is a parameterization joint of ONE surface, not an
        // intersection of two — the confirmed spelling (see `Edge::surfaces`), guarded
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
                orientation: facing(bottom_flipped),
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
                orientation: facing(top_flipped),
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
#[path = "tests/lib.rs"]
mod tests;
