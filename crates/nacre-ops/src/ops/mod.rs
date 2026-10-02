//! Feature operations: the public sketch/extrude API and the `apply`/
//! `replay` driver. The top layer — it composes the boolean engine ([`crate::boolean`]) and rigid
//! transform ([`crate::transform`]) over the plane substrate below.

use crate::BoolError;
use crate::BoolKind;
use crate::boolean::boolean;
use crate::construct::{Seg3, Swept};
use crate::transform::transform;
use nacre_exact::{Axis, Isometry, Rat};
use nacre_geom::Plane;
use nacre_geom::intersect::{RingSide, orient2d_rat};
use nacre_geom::mixed::{
    Edge2d, mixed_ring_self_intersection, mixed_rings_cross, point_in_mixed_ring,
};
use nacre_math::{Point2, Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::CylinderDef;
use nacre_topo::PointCache;
use nacre_topo::{Edge, Face, HalfEdge, Loop, Model, Orientation, Shell, Solid, Surface, Vertex};
use std::borrow::Cow;

mod apply;
mod datum;
/// Solids built through the product's own operations — see the module doc.
#[cfg(any(test, feature = "test-util"))]
pub mod fixtures;
mod frame;
mod plane;
mod prism;
mod profile;

pub use apply::*;
use datum::*;
pub use frame::*;
pub use plane::*;
pub(crate) use prism::*;
pub use profile::*;

/// A modelling operation.
///
/// ★ **`DatumPlane` is the large variant and it is not boxed.** It carries a [`SketchPlane`] in
/// its `Stated` arm, which holds a plane's exact rational definition (~400 bytes). An op log is
/// tens of entries long and is walked once per replay, so the wasted space is measured in
/// kilobytes; boxing would buy that back at the cost of an indirection on the one type a caller
/// constructs by hand. (`Extrude` names its plane by handle, inside its [`SketchFrame`], which
/// keeps it small.)
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq)]
pub enum Operation {
    /// **Put a plane in the model because the caller named it** — not because a face lies on it.
    ///
    /// Without it a log could name only two kinds of plane: the three seeded world planes and the
    /// plane of a face it had already built ([`face_sketch_frame`]); every other plane would have
    /// to travel by value. This is the one variant that carries a plane's statement (`Stated`).
    ///
    /// ★ **It is an operation, not a plain function, because `replay` must reproduce the handle.**
    /// A plane minted outside the log leaves a model that is not self-contained, and the
    /// "handles in a log are index vocabulary" contract (see `docs/design.md`) is false for it.
    ///
    /// The plane may already exist — planes are the one thing interned at construction — in which
    /// case this returns the handle that exists and the arena does not grow. That is the intended
    /// answer, not a special case: *same plane, same handle*.
    DatumPlane { def: DatumDef },
    /// Extrude `profile`, drawn in `frame`, by `dist` along that frame's `ŵ` — or against it when
    /// `dist` is negative.
    ///
    /// ★ **The plane is named, not carried.** `frame` holds a `Handle<Surface>`, so the plane it
    /// sketches on is one the model already has — a seeded world plane
    /// ([`SketchFrame::world`]), the plane of a face ([`face_sketch_frame`]), or one a
    /// [`Operation::DatumPlane`] put there. That is what makes the base cap a *shared* handle
    /// rather than a second statement of the same plane, whichever way the prism runs.
    ///
    /// ★ **The frame says where `ŵ` points and where the sketch's `+u` lies; the sign says which
    /// side of the plane the body takes.** `flip` is the frame maker's, decided from the truth, so
    /// a face's frame always faces outward — and turning the frame over would re-derive its axes
    /// and land the same sketch mirrored. Sweeping against `ŵ` keeps the axes and is the only
    /// spelling of «drawn on this face, cut into it». `0` names no body
    /// ([`OpError::ZeroDistance`]).
    Extrude {
        frame: SketchFrame,
        profile: Profile2d,
        dist: f64,
    },
    /// Boolean of two live solids (M5). Inputs the engine cannot answer are
    /// rejected with [`BoolError`].
    Boolean {
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    },
    /// Rigid-body transform: supersede `solid` by its image under `isometry` — a rotation, then
    /// a translation. The `Isometry` is the exact definition (op-log truth); the geometry is a
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
    /// `live_solids` without removing anything**. Every other edit supersedes its
    /// input, so this is what makes "cut with the same tool twice", "keep the original and a moved
    /// copy", and pattern/mirror sugar expressible at all.
    Copy { solid: Handle<Solid> },
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
    /// cylinders and a plane, which no [`nacre_topo::Vertex`] states yet, and the ruling between
    /// the two cylinder walls has no curve the kernel derives (`derive_edge_curve` declines two
    /// distinct cylinders) — the same frontier as the cylinder–cylinder boolean (M6b). A lens, a
    /// cam lobe; a straight step between the arcs is what builds today.
    ArcsMeetAtVertex,
    /// An extrusion distance of zero: a prism of no thickness. Either sign is a side
    /// ([`Operation::Extrude`]); zero is neither.
    ZeroDistance,
    /// An extrusion distance outside the decimal window (`Rat::from_decimal` — `~1e38` above,
    /// `~1e-22` below for a full-width value): it has no rational truth for the sweep to be
    /// computed in. The sibling of [`OpError::ProfileOutsideDecimalWindow`], named at the
    /// operation's door — this used to fall silently to a point-less f64 prism.
    DistOutsideDecimalWindow,
    /// The sketch plane (or the frame chain carrying it) has no exact form to build in: its
    /// axes fell outside the decimal window or were degenerate, or the placement arithmetic
    /// overflowed `i128`. A prism built silently in f64 instead would record no exact points and
    /// could not survive a motion — a named reject is the honest answer.
    PlaneWithoutExactForm,
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
    /// A sketch asked to stand on a face or surface that is not planar (only planar faces carry a
    /// sketch frame; curved-face features arrive with the quadric milestones).
    NonPlanarFace,
    /// A face asked for its sketch frame belongs to no live solid's outer shell (a stale or
    /// non-live handle).
    FaceNotInLiveSolid,
    /// A boolean operation failed (M5).
    Boolean(BoolError),
    /// A `Transform` input solid is not live (a stale or non-live handle).
    SolidNotLive,
    /// A vertex or edge of a `Transform`/`Copy` input names a surface that is not one
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
    /// vocabulary** (see `docs/design.md`): [`replay`] re-anchors each one onto the model it is
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
    /// are in that state (`tests/instruments/point_width.rs`,
    /// `a_datum_on_straddling_carriers_has_no_name`).
    ///
    /// ★★★★ **What opens this is not the judging layer.** The homogeneous lift
    /// does not open the population. Such
    /// a plane has **no exact name**, and from there:
    ///
    /// ```text
    /// no name → frame_chain declines → no SketchFrame → never a base cap → never in a plane table
    /// ```
    ///
    /// so it never gets as far as being judged. What would open it is a **frame for a nameless
    /// plane** — a realization that takes the plane's coefficients as intervals at a precision
    /// rather than as exact rationals or `BigInt`s, which is what `MoveNode::Frame` and
    /// `FrameWide` both require today;
    /// `a_plane_with_no_name_cannot_host_a_sketch` runs the chain.
    ///
    /// It is named separately so that "how much does this cost us today" stays countable.
    VerticesInMixedFrames,
    /// A mixed-frame datum's **judged frame could not be decided** at the fixed rung: no
    /// arbitrary-axis branch's squared length — or the normal's, or the origin's denominator —
    /// could be bounded away from zero (`nacre_judge::FrameThrough::of`).
    ///
    /// ★ Deliberately **not** [`OpError::CollinearVertices`] and not
    /// [`OpError::DegenerateGeometry`]: both claim the construction *is* degenerate, and an
    /// interval that fails to clear zero proves nothing of the kind — the points may be exactly
    /// collinear or merely too close to call. Failing to prove health is its own cause, and a
    /// refusal is named by its cause.
    ThroughFrameUndecided,
    LogHandleOutOfRange {
        /// Which store the index was meant for.
        cell: LogCell,
        /// The index the log named.
        index: u32,
    },
}

/// Which store an operation-log handle indexes. (`Surface` is what datum ops name.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogCell {
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
    /// (`tests/instruments/point_width.rs`). Naming the vertices keeps the statement exact.
    ///
    /// ★ **The order is the direction.** The three are sorted before they are stored — the same
    /// vertices are the same plane in any order — but the caller's order fixes a normal by the
    /// right-hand rule, and `flip` is measured against it exactly as it is for a stated plane.
    /// Reversing two of them returns *the same handle with the opposite frame*, which is the only
    /// way a caller can choose a side (`dist` is positive-only).
    ///
    /// Rejects by cause rather than by one blanket failure, because the causes have different
    /// futures: [`OpError::VerticesInMixedFrames`] waits on a frame for a nameless plane (see
    /// there), and the rest are the caller's. (Width is not a cause:
    /// a meet wider than `Rat` names its plane through `plane_name_from_meets`.)
    ThroughVertices([Handle<Vertex>; 3]),
    /// **`dist` away from a plane the model already holds**, stated inside that plane's own frame
    /// as the rational triple `(0,0,d), (1,0,d), (0,1,d)`.
    ///
    /// ★ **This is the exact form of an offset, and the frame is why.** In the world, "d along the
    /// normal" of a tilted plane is `p + d·n̂` — irrational, because `n̂` carries a square root.
    /// Inside the frame the same plane is `w = d` and every coordinate is a written decimal; the
    /// irrationality lives in the frame's realization, which is machinery the kernel already has.
    /// It is also the alternative to floating a sketch origin off
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
    /// The boolean result solids (supersede both inputs). Usually one; a boolean that severs the
    /// body yields several, and `Cut(A, A)` (deferred) would yield none.
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
    /// `flip` makes `ŵ` face the normal the statement names — a stated plane's normal, the
    /// right-hand turn of the caller's vertex order — so a frame is returned facing the way its
    /// author said even when the plane interned onto a handle facing the other way. An offset
    /// returns its plane's canonical frame (`flip` false).
    DatumPlane {
        plane: Handle<Surface>,
        frame: SketchFrame,
    },
}

/// **Which frame a sketch lives in** — the plane (a handle: one statement of the plane, shared
/// with every face on it), its [`nacre_topo::FramePlacement`], and whether the plane's canonical
/// coefficients need negating to face the way the sketch does. Everything a
/// [`nacre_topo::Motion::Frame`] node needs, before the model has one.
///
/// ★★ **The fields are private and the constructors validate** — the reason this type is not a
/// plain record. A `Named` placement is a *claim*: "this origin lies on that plane, this
/// direction crosses its normal". [`SketchFrame::named`] checks the claim exactly, at
/// construction, and rejects by name — silently substituting `Canonical` would move a caller's
/// sketch and answer a question they did not ask. A public
/// field would let a literal walk around the check.
///
/// ★ **`flip` is not the caller's to state** — the canonical coefficients carry no direction, so
/// which way `ŵ` must face is a fact about the *use* (a sweep's sense, a face's outward normal),
/// decided by the consuming operation from the truth: which way that use faces against the
/// plane's own facing. Constructors set `false`.
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
    /// canonical name normalizes — pinned by `a_world_frame_faces_its_axis`, and by
    /// `a_seeded_planes_canonical_frame_is_the_world_basis_exactly` underneath it. Re-seed the
    /// world planes differently and this turns silently; the tests are what stop that.
    ///
    /// To sketch facing the *other* way, state a plane facing that way
    /// ([`Operation::DatumPlane`] turns `flip` to face the normal you state) — the same move as
    /// writing `SketchPlane::from_origin_normal(o, -ẑ)` today.
    pub fn world(model: &Model, axis: nacre_exact::Axis) -> SketchFrame {
        let plane = model.world_plane(axis);
        match axis {
            // The two whose derived frame already is the convention.
            nacre_exact::Axis::Z | nacre_exact::Axis::X => SketchFrame::canonical(plane),
            // ZX: `+u = +ẑ`, stated because it cannot be derived.
            nacre_exact::Axis::Y => SketchFrame {
                plane,
                placement: nacre_topo::FramePlacement::Named {
                    origin: [nacre_exact::Rat::from_int(0); 3],
                    ref_dir: [
                        nacre_exact::Rat::from_int(0),
                        nacre_exact::Rat::from_int(0),
                        nacre_exact::Rat::from_int(1),
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
    ///   is nonzero. Exact, total ([`nacre_exact::plane_residual_sign`]): a `Wide` name checks
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
        let lift = |v: [f64; 3]| -> Result<[nacre_exact::Rat; 3], OpError> {
            let mut out = [nacre_exact::Rat::from_int(0); 3];
            for (o, c) in out.iter_mut().zip(v) {
                *o = nacre_exact::Rat::from_decimal(c).ok_or(OpError::FrameOutsideDecimalWindow)?;
            }
            Ok(out)
        };
        let (origin, ref_dir) = (lift(origin.as_array())?, lift(ref_dir.as_array())?);
        let name = model
            .surface_name
            .get(&plane)
            .ok_or(OpError::PlaneWithoutExactForm)?;
        if nacre_exact::plane_residual_sign(name, origin) != 0 {
            return Err(OpError::OriginNotOnPlane);
        }
        // The existing frame machinery is the judge — `WideFrame::named_of` is total in width
        // (arbitrary precision), so its only `None` is the projected `u_raw` vanishing: exactly
        // the parallel-or-zero claim this constructor rejects. (The narrow `plane_frame_named`
        // is not consulted here: its `None` can mean `i128` overflow, which is a width fact,
        // not a defect in the claim.)
        if nacre_judge::WideFrame::named_of(name, &origin, &ref_dir, false).is_none() {
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
    /// decided by the consuming operation, `false` as constructed.
    pub fn flip(&self) -> bool {
        self.flip
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
/// [`prism::frame_road`] measured that this used to happen for rational-tilt planes; here it must not
/// happen at all.
#[cfg(test)]
#[path = "../tests/ops_frame_differential.rs"]
mod frame_differential;
