//! Rigid-body transform of a solid (design overhaul stage 1): copy a solid under an isometry,
//! remapping every surface/curve and carrying each vertex's exact `Origin` forward so a rotated
//! operand stays exactly defined. [`copy`] is the same walk with no motion at all.

use crate::OpError;
use nacre_geom::{AxisMirror, Cylinder, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Motion, MotionNode, Shell, Solid, Vertex, VertexDef,
};
use std::collections::{HashMap, HashSet};

/// Supersede `solid` by its image under `isometry` (stage 1a: a rational
/// translation). Clones the solid's cells with moved geometry ([`transform_solid`])
/// and drops the input from `live_solids` — the op-log is the truth.
pub(crate) fn transform(
    model: &mut Model,
    solid: Handle<Solid>,
    isometry: &Isometry,
) -> Result<Handle<Solid>, OpError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !defs_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let out = transform_solid(model, solid, &Xform::Rigid(isometry))?;
    model.live_solids.retain(|&s| s != solid);
    Ok(out)
}

/// An independent twin of `solid` at the same place — **the one operation that only adds to
/// `live_solids`** (design §2 supersede semantics). Every other edit supersedes: `transform` and
/// `boolean` drop their inputs, so without this the kernel can *move* a solid but never *copy* one,
/// and "cut with the same tool twice" or "keep the original and a moved copy" cannot be expressed.
///
/// It is [`transform_solid`] under a zero translation: the cells are duplicated (two live solids
/// must not share cells — a shared edge would read as four face uses and break the manifold check)
/// while the geometry is rebuilt bit-for-bit (a pure translation keeps a plane's exact `raw`).
///
/// **The carried measured tolerance is exact here, not provisional.** A copy keeps `tol`
/// unchanged, and its point and three planes are bit-identical, so the measured residual
/// cannot have changed. (A *recorded* move re-realizes coordinates and drops the measurement
/// instead — see the tol rule in pass 3.)
///
/// A non-live input is rejected rather than resurrected: reusing a superseded handle is a caller
/// bug, and letting it succeed would hide it.
pub(crate) fn copy(model: &mut Model, solid: Handle<Solid>) -> Result<Handle<Solid>, OpError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !defs_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let zero = Isometry::translation([Rat::from_int(0); 3]);
    // `push_solid` registers the twin as live; the input is *not* retained away — that missing
    // line is the whole difference from `transform`.
    transform_solid(model, solid, &Xform::Rigid(&zero))
}

/// Whether every vertex **definition** in `solid` names surfaces the walk will remap — that is,
/// face surfaces of this solid. The remap cannot proceed otherwise, and a kernel must decline
/// rather than abort ("honest-reject > silent-wrong"; `assemble_fuse_cut` takes the same line at
/// its own naming failure), so the three entry points check first and reject.
///
/// ★ **Widened from `Discovered`-only to every `Some(definition)` (S7)** — the old gate
/// inspected discovered vertices because only that road panicked; with the definition becoming
/// the vertex's identity, every def must remap. The widening is expected to fire **zero** times
/// (every producer writes definitions from its own face surfaces — the coverage gate measures
/// 100%), and the suite staying green is that evidence: a firing here is a new rejection, which
/// is a failing test. The positive control (`a_foreign_definition_is_rejected`) shows the gate
/// actually bites.
///
/// ★★ **That expectation was once false, and is now enforced where it is made.** A boolean's
/// assembly used to name a four-plane vertex by the arrangement's canonical triple, which could
/// include a plane the result kept no face on; the solid built fine and refused to move, here,
/// two operations after the mistake. `assemble_fuse_cut` now derives the name from the faces that
/// meet the vertex *and* asserts this predicate on what it returns (debug builds), so the whole
/// suite is the corpus for "every producer writes definitions from its own face surfaces" rather
/// than this gate being the first to find out.
pub(crate) fn defs_are_remappable(model: &Model, solid: Handle<Solid>) -> bool {
    let src = model.solids.get(solid);
    let shells: Vec<Handle<Shell>> = std::iter::once(src.outer)
        .chain(src.cavities.iter().copied())
        .collect();
    let mut surfs: HashSet<Handle<Surface>> = HashSet::new();
    for &sh in &shells {
        for &fh in &model.shells.get(sh).faces {
            surfs.insert(model.faces.get(fh).surface);
        }
    }
    let named = |def: &VertexDef| def.carriers().all(|s| surfs.contains(&s));
    for &sh in &shells {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    let edge = model.edges.get(he.edge);
                    for vh in edge.vertices.iter() {
                        if !named(&model.vertices.get(*vh).def) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

/// Supersede `solid` by its reflection in the coordinate plane `axis = offset`.
///
/// Lengths are preserved and **handedness is reversed** — this is a reflection, not a negative
/// scale (the kernel has no scale, and `Isometry` cannot hold an improper motion). The plane is
/// axis-aligned, the same restriction rotation already has: a general plane's unit normal is
/// irrational, so the image could not be reproduced from an exact definition.
///
/// Like `transform` this consumes its input; pair it with [`copy`] to keep the original — which
/// is the usual move, since mirroring exists to build the other half of a symmetric part.
pub(crate) fn mirror(
    model: &mut Model,
    solid: Handle<Solid>,
    axis: Axis,
    offset: Rat,
) -> Result<Handle<Solid>, OpError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !defs_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let out = transform_solid(model, solid, &Xform::mirror(axis, offset))?;
    model.live_solids.retain(|&s| s != solid);
    Ok(out)
}

/// The motion nodes this `motion` appends to `parent`, or `None` when it records nothing.
///
/// **One rule for all three motions.** An `Isometry` is "rotate about a pivot, then translate",
/// so it contributes up to two nodes in that order — **order is the definition**, since the two
/// do not commute. A reflection contributes one. A reflection used to be carried a second way
/// entirely (conjugate the input's chain, `M ∘ R = (M R M⁻¹) ∘ M`, and mirror its root vertex),
/// which meant two mechanisms for one question and a chain that could not hold the reflection it
/// had just performed. It is a motion; it goes in the chain.
///
/// A node is *omitted* only when the motion changes nothing the definition needs to say: a zero
/// translation is the identity (`copy` is `transform_solid` under one), a 90°-family rotation of
/// an exact datum keeps it exact, and a translation that lands every coordinate back on an exact
/// `f64` does too (`exact_translate`, from [`translation_is_exact`]).
///
/// **Both exceptions lapse once the datum already has a history.** "This motion kept the
/// coordinates exact" is a statement about *this step*; it says nothing about the chain, and a
/// chain missing a link describes the datum as it was *before* that link — silently. So once
/// `leaf` is `Some`, every motion is recorded.
///
/// The translation half of that used not to be reachable — an inexact rotation leaves coordinates
/// using the full 53-bit mantissa, so no nonzero translation of such a solid is exact. A
/// reflection in a dyadic plane *is* exactness-preserving and can carry a history, which is what
/// makes the rule load-bearing rather than merely correct.
fn chain_motion(
    model: &mut Model,
    parent: Option<Handle<MotionNode>>,
    motion: &Xform<'_>,
    exact: bool,
) -> Option<Handle<MotionNode>> {
    let mut leaf = parent;
    match motion {
        Xform::Rigid(iso) => {
            if let Some(r) = iso.rotate.filter(|_| !iso.is_exact() || leaf.is_some()) {
                leaf = Some(model.push_motion(
                    Motion::Rotate {
                        axis: r.axis,
                        point: r.point,
                        angle: r.angle,
                    },
                    leaf,
                ));
            }
            if iso.translate.iter().any(|r| r.numer() != 0) && (!exact || leaf.is_some()) {
                leaf = Some(model.push_motion(
                    Motion::Translate {
                        offset: iso.translate,
                    },
                    leaf,
                ));
            }
        }
        Xform::Mirror { axis, offset, .. } => {
            if !exact || leaf.is_some() {
                leaf = Some(model.push_motion(
                    Motion::Mirror {
                        axis: *axis,
                        offset: *offset,
                    },
                    leaf,
                ));
            }
        }
    }
    leaf
}

/// Whether `motion` lands **every** coordinate of `solid` back on an exact `f64`.
///
/// The twin of `Isometry::is_exact` for the data-dependent motions: a 90°-family rotation keeps any
/// exact datum exact by its cos/sin alone, but a translation's or a reflection's exactness depends
/// on the *data* — `p + t` and `2c − p` are representable only when the parameter is dyadic **and**
/// the result fits 53 bits. So it is measured, not derived, and when it holds the coefficients stay
/// the truth: no node, no tolerance, and the judgment keeps the exact `f64` predicate path that
/// `Witness::is_rotated` gates.
///
/// **Per solid, not per face.** `transform_solid` requires a solid's boundary vertices to share
/// one node ("uniform motion"), and a *dyadic* parameter does split by magnitude — measured, `0.5`
/// is exact at `p = 1.0` and rounds at `p = 1e17`. Deciding per face would then leave one solid
/// with some vertices carrying a node and some not, and break that invariant. One pass over every
/// vertex and every face's plane origin — the two things a judgment reads — and one rounding
/// anywhere puts the whole solid on the recorded path.
///
/// (The per-plane **surface-statement** exemption — pass 1's `invariant` flag, a plane the
/// motion fixes as a set — is a different question and does not touch this one: the vertices
/// of such a plane still move and still share the solid's node story; only the plane's own
/// statement needed no new spelling.)
///
/// This answers for the **translation/reflection** part only; `chain_motion` reads
/// `Isometry::is_exact` for the turn itself, and a rotation that is not exact makes the whole
/// chain recorded anyway (the exception lapses the moment a history exists).
fn motion_is_exact(model: &Model, solid: Handle<Solid>, motion: &Xform<'_>) -> bool {
    // Representable *and* reached by the `f64` arithmetic the producer actually performs — the
    // second half matters: an exactly-representable answer the producer does not land on would
    // make the definition and the cached coordinate disagree.
    let translated = |x: f64, t: Rat| match Rat::try_from_f64(x).and_then(|r| r.checked_add(t)) {
        Some(sum) => Rat::try_from_f64(sum.to_f64()) == Some(sum) && x + t.to_f64() == sum.to_f64(),
        None => false,
    };
    let reflected = |x: f64, c: Rat| match Rat::try_from_f64(x)
        .and_then(|r| c.checked_mul(Rat::from_int(2))?.checked_sub(r))
    {
        Some(v) => Rat::try_from_f64(v.to_f64()) == Some(v) && 2.0 * c.to_f64() - x == v.to_f64(),
        None => false,
    };
    let all = |p: Point3| match motion {
        Xform::Rigid(iso) => p
            .as_array()
            .iter()
            .zip(iso.translate)
            .all(|(&x, t)| translated(x, t)),
        Xform::Mirror { axis, offset, .. } => reflected(p.as_array()[axis.index()], *offset),
    };
    let src = model.solids.get(solid);
    // ★ S6a: the surfaces' exact points must survive the no-node path too. An exact motion
    // carries a `Constructed` surface's rational triple through `point_rat`/`mirror_point_rat`,
    // and that arithmetic can overflow `i128` even when every f64 above is exact (the two
    // conditions are independent). Dropping the points — the old behaviour — is what minted the
    // point-less population `Inexact` grows from; recording a node instead keeps the original
    // triple as the pre-motion truth. The probe is per solid like everything here, so one
    // overflowing surface puts the whole solid on the recorded path rather than splitting it.
    let points_move = |s: Handle<Surface>| -> bool {
        // ★★★★★ **An exhaustive `match`, deliberately — this used to be the most dangerous
        // `let`-`else` in the file, and the compiler could not see it.** Adding
        // `PlanePoints::Through` produced exactly two non-exhaustive-match errors and **not**
        // one here — the fallback arm would have swallowed the new variant silently, answering
        // for geometry it had never seen. Spelled as a match, the next variant (M6's) is a
        // compile error at exactly this decision.
        match model.surface_truth(s) {
            // The very function pass 1 will transport with — sharing it is what makes the
            // probe's promise ("this will not overflow") structural rather than a parallel
            // re-derivation.
            nacre_topo::SurfaceTruth::Plane {
                points: nacre_topo::PlanePoints::Known(p),
                ..
            } => transport_points(motion, *p).is_some(),
            // A `Through` plane's truth carries geometry **by reference** — the no-node path
            // would transport nothing while the f64 cache moved, leaving truth and cache
            // describing different planes. So it refuses that path.
            nacre_topo::SurfaceTruth::Plane {
                points: nacre_topo::PlanePoints::Through(_),
                ..
            } => false,
            // A cylinder's truth carries geometry since M6-0, so the probe asks the same
            // question it asks a plane: can the very function pass 1 will transport with
            // carry this statement? (`transport_cylinder` — shared, like `transport_points`
            // above.) A mirror answers `false` here and is then rejected by pass 1's
            // `MirrorNotPlanar` before any transport runs.
            nacre_topo::SurfaceTruth::Cylinder { def, .. } => {
                transport_cylinder(motion, def).is_some()
            }
        }
    };
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            // Exhaustive for the same reason as `points_move` above: a new `Surface` variant
            // (M6's sphere/cone) must be a compile error here, not a silently skipped probe.
            match model.surface(face.surface) {
                Surface::Plane(pl) => {
                    if !all(pl.origin()) {
                        return false;
                    }
                }
                // Licensed by this function's own charter: it probes "the two things a
                // judgment reads" — vertex coordinates and *plane* origins — and no M5
                // judgment reads a cylinder's origin (cylinders never reach a boolean).
                // When M6 admits cylinders to judgments, this arm owes an origin/axis probe.
                Surface::Cylinder(_) => {}
            }
            if !points_move(face.surface) {
                return false;
            }
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in model.edges.get(he.edge).vertices.iter() {
                        if !all(model.vertex_point(vh)) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

/// A plane's exact triple carried through an **exact** (recorded-nothing) motion — the one
/// transport both `motion_is_exact`'s probe and pass 1 use, so the probe's feasibility answer
/// and the actual transport cannot disagree. `None` on `i128` overflow, which the probe turns
/// into "record a node instead".
fn transport_points(motion: &Xform<'_>, p: [[Rat; 3]; 3]) -> Option<[[Rat; 3]; 3]> {
    let each = |q: [Rat; 3]| -> Option<[Rat; 3]> {
        match motion {
            Xform::Rigid(iso) => iso.point_rat(q),
            Xform::Mirror { axis, offset, .. } => nacre_scalar::mirror_point_rat(q, *axis, *offset),
        }
    };
    Some([each(p[0])?, each(p[1])?, each(p[2])?])
}

/// The cylinder twin of [`transport_points`] — pass 1 transports with the very function the
/// probe checked, so "this will not overflow" is structural, not a parallel re-derivation.
///
/// Rigid: the origin rides `point_rat`; the two directions ride the rotation alone
/// (`dir_rat` — a direction is a difference of points, so pivot and translation cancel); the
/// radius is invariant under a rigid motion. The image is rebuilt through the checked
/// constructor, which a rigid motion cannot fail except on overflow — then `None`, the
/// conservative answer the probe turns into "take the recorded path".
///
/// Mirror: `None`, deliberately — pass 1's `MirrorNotPlanar` rejection preempts (a mirrored
/// cylinder has no cache image, `Xform::surface`), so this arm only ever steers the probe and
/// must refuse quietly rather than panic.
fn transport_cylinder(
    motion: &Xform<'_>,
    def: &nacre_topo::CylinderDef,
) -> Option<nacre_topo::CylinderDef> {
    match motion {
        Xform::Rigid(iso) => nacre_topo::CylinderDef::new(
            iso.point_rat(def.origin())?,
            iso.dir_rat(def.dir())?,
            iso.dir_rat(def.ref_dir())?,
            def.radius(),
        ),
        Xform::Mirror { .. } => None,
    }
}

/// The motion history a moved surface's image carries — the surface twin of the vertex
/// `Origin` rules, in S6b form (the old `SurfaceDef` table collapsed to one field):
///
/// | source motion | motion records a node | motion records nothing |
/// |---|---|---|
/// | `None` (world) | `Some(leaf)` — the points stay pre-motion | `None` (points transported exactly) |
/// | `Some(m)` | `Some(leaf)`, hanging off `m` — replaying from the root applies every motion once | unreachable (the exceptions lapse once a history exists) |
///
/// The `None (world)` row has a third outcome decided **before** this function is asked: a
/// rigid motion that *fixes the plane* (`Isometry::fixes_plane` — pass 1's `invariant` flag)
/// records nothing and carries the points **verbatim**, so the image interns back onto the
/// source handle. Per plane, not per solid — the caps of a turned block take it while the
/// walls land in this table.
///
/// **One table for all three motions.** A reflection used to have its own column here, carrying a
/// moved surface by conjugating its chain over a mirrored witness; it now appends a node like
/// everything else, so the reflection is *in* the definition rather than folded into it.
///
/// A cylinder rides the same rows through its own motion slot — the old `SurfaceDef` path
/// demoted a moved cylinder to `Inexact` (a `Constructed` source with no points to move); its
/// history is simply recorded now.
fn moved_surface_motion(
    model: &mut Model,
    src: Handle<Surface>,
    motion: &Xform<'_>,
    exact: bool,
    surf_rot: &mut HashMap<Option<Handle<MotionNode>>, Option<Handle<MotionNode>>>,
) -> Option<Handle<MotionNode>> {
    // **Each surface chains from its own leaf, not the solid's.** One solid does not have one
    // history: a boolean between differently-moved operands hands back walls that came from
    // different ones. Memoized per distinct parent so surfaces that did share a history still do.
    let parent = model.plane_motion(src);
    let leaf = match surf_rot.get(&parent) {
        Some(&h) => h,
        None => {
            // `entry` cannot hold a `&mut Model` across the closure, so look up then insert.
            let h = chain_motion(model, parent, motion, exact);
            surf_rot.insert(parent, h);
            h
        }
    };
    // Nothing recorded: the motion kept the data exact and the image keeps the source's own
    // history. (`leaf` is always `Some` when `parent` is — the exceptions lapse once a history
    // exists — so `or` never resurrects a stale parent past a recorded node.)
    leaf.or(parent)
}

/// How [`transform_solid`] maps a solid's cells. One walker serves both kinds so the seven
/// passes are not duplicated; the kinds differ in exactly three places — how a point/direction
/// maps, whether a curved surface can be carried at all, and whether loops must be rewound.
pub(crate) enum Xform<'a> {
    /// A proper motion: rotation then translation. Preserves handedness.
    Rigid(&'a Isometry),
    /// A reflection in a coordinate plane. Reverses handedness, so `det = −1`.
    ///
    /// The plane is carried **twice**: `m` is the `f64` map the coordinates actually go through,
    /// `(axis, offset)` the exact statement of the same plane that the motion history records.
    /// They used to travel as separate arguments, and an `Option` that could in principle arrive
    /// empty; one variant cannot lose half of itself. [`Xform::mirror`] is the only constructor,
    /// so the two halves cannot be made to disagree either.
    Mirror {
        m: AxisMirror,
        axis: Axis,
        offset: Rat,
    },
}

impl Xform<'_> {
    /// The reflection in `axis = offset`, with its `f64` map **derived** from the exact plane
    /// rather than handed in beside it.
    fn mirror(axis: Axis, offset: Rat) -> Xform<'static> {
        Xform::Mirror {
            m: AxisMirror::new(axis.index(), offset.to_f64()).expect("axis index is 0..3"),
            axis,
            offset,
        }
    }

    fn point(&self, p: Point3) -> Point3 {
        match self {
            Xform::Rigid(iso) => Point3::from_array(iso.apply_point(p.as_array())),
            Xform::Mirror { m, .. } => m.point(p),
        }
    }

    /// `None` when the variant has no image under this motion — a mirrored cylinder, whose
    /// parametrisation handedness is a curved-geometry decision (see `Surface::mirrored`).
    fn surface(&self, s: &Surface, offset: Vector3) -> Option<Surface> {
        match self {
            Xform::Rigid(iso) => Some(transform_surface(s, iso, offset)),
            Xform::Mirror { m, .. } => s.mirrored(*m),
        }
    }

    /// A reflection negates the normal a loop's winding implies (`R(a) × R(b) = −R(a × b)`), so
    /// every loop is rewound to put it back — and then the `Orientation` flag needs no change,
    /// because a reflection preserves dot products.
    fn reverses_orientation(&self) -> bool {
        matches!(self, Xform::Mirror { .. })
    }

    /// The isometry, for the rotation-forest bookkeeping that only proper motion does.
    fn rigid(&self) -> Option<&Isometry> {
        match self {
            Xform::Rigid(iso) => Some(iso),
            Xform::Mirror { .. } => None,
        }
    }
}

/// A surface moved by `isometry`. A pure translation uses `translated` (the normal —
/// and its exact `raw` — is unchanged). A rotation rebuilds from the moved
/// origin/normal via the constructor (rotation makes `raw` irrational, as expected —
/// the plane then carries tol, judged by CIP later).
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

/// Clone `solid` into a new solid with every cell's geometry mapped by `motion`,
/// preserving topology, shared cells (surfaces/curves/vertices/edges are deduped),
/// inner-loop holes, cavity shells, and each vertex/edge `Origin` (a `Discovered`
/// definition's plane handles are remapped to the moved surfaces). Cells are pushed in a
/// **deterministic traversal order** (shell → face → loop) with per-cell dedup maps, so the same
/// op-log reproduces identical handles (replay determinism, DNA 3).
///
/// Face orientation flags are carried unchanged for **both** kinds of motion. A rigid motion
/// turns the normal and the winding together; a reflection negates the winding's implied normal,
/// which [`Xform::reverses_orientation`] undoes by rewinding every loop — and since a reflection
/// preserves dot products, `sign(plane.normal · n_out)` is then unchanged too, which is exactly
/// what the `Orientation` flag records.
///
/// `Err` only when the motion has no image for some cell's geometry (a mirrored cylinder).
fn transform_solid(
    model: &mut Model,
    solid: Handle<Solid>,
    motion: &Xform<'_>,
) -> Result<Handle<Solid>, OpError> {
    let offset = motion
        .rigid()
        .map(|iso| Vector3::from_array(iso.offset_f64()))
        .unwrap_or_else(Vector3::zero);
    let src = model.solids.get(solid).clone();

    // Decided once, for the whole solid — see [`motion_is_exact`].
    let exact = motion_is_exact(model, solid, motion);
    // Whether the source already carries a motion history (any face surface's truth records
    // one). With the vertex-side motion gone (S7 — the faces record it themselves, and a
    // vertex's motion was always its faces'), this is what remains of the old
    // "exceptions lapse once the solid has a history" test: an exact move of a fresh solid
    // records nothing and keeps the measured tolerances verbatim; anything else re-realizes
    // coordinates, so a measured tolerance no longer describes them.
    let prior_history = {
        let src = model.solids.get(solid);
        std::iter::once(src.outer)
            .chain(src.cavities.iter().copied())
            .flat_map(|sh| model.shells.get(sh).faces.clone())
            .any(|fh| {
                !matches!(
                    model.surface_truth(model.faces.get(fh).surface),
                    nacre_topo::SurfaceTruth::Plane { motion: None, .. }
                        | nacre_topo::SurfaceTruth::Cylinder { motion: None, .. }
                )
            })
    };
    let keeps_tol = exact && !prior_history;

    // Deterministic order: outer shell then cavities; each shell's faces in order.
    let shell_order: Vec<Handle<Shell>> = std::iter::once(src.outer)
        .chain(src.cavities.iter().copied())
        .collect();
    let face_order: Vec<Handle<Face>> = shell_order
        .iter()
        .flat_map(|&sh| model.shells.get(sh).faces.clone())
        .collect();

    // Pass 1 — surfaces (dedup, moved): needed before vertex `Origin` remap.
    //
    // A moved surface's coefficients are only the truth while the motion kept them exact, so each
    // one states its provenance here (`SurfaceDef`). The witness for a rotation is the face's own
    // **pre-rotation** triangle, read on the spot — the definition must not depend on anything
    // else in the model surviving, and those points are exactly the ones a rotated operand's
    // `tri_pt3` is built from today, which is why this changes no answer.
    //
    // **A surface chains from its own leaf, not from the solid's.** One solid does not have one
    // rotation history: a boolean between differently-rotated operands hands back a result whose
    // walls came from different rotations, and its vertices are all `Discovered`, so the
    // vertex-side `solid_rotation` cannot answer for it at all. `surf_rot` memoizes one new forest
    // node per distinct parent leaf, so surfaces that did share a history still share it.
    let mut surf_rot: HashMap<Option<Handle<MotionNode>>, Option<Handle<MotionNode>>> =
        HashMap::new();
    let mut surf_map: HashMap<Handle<Surface>, Handle<Surface>> = HashMap::new();
    // ★ Set when the surface the model handed back points the other way from the one built here.
    // A copied face keeps its `Orientation` because its surface moved with it — but a *shared*
    // surface did not, so the same outward direction has to be spelled the other way. This is a
    // different field from the mirror's rewind below, which turns loop *winding*; the two do not
    // interact (`Loop::reversed`'s doc spells out why a reflection touches only the winding).
    let mut surf_flip: HashMap<Handle<Surface>, bool> = HashMap::new();
    for &fh in &face_order {
        let s = model.faces.get(fh).surface;
        if surf_map.contains_key(&s) {
            continue;
        }
        let moved = motion
            .surface(model.surface(s), offset)
            .ok_or(OpError::MirrorNotPlanar)?;
        let src_truth = model.surface_truth(s).clone();
        // ★ **A motion that fixes this plane restates nothing — the source statement already
        // states the image.** The per-plane sibling of `motion_is_exact`'s whole-solid
        // exemption: a rigid motion mapping this plane onto itself *as a set* (axis ∥
        // normal, in-plane translation — `Isometry::fixes_plane`, exact) leaves the same
        // points and no node to record, and re-pushing that statement interns back onto
        // the **source handle** — the road `Copy` already takes for the identity. The caps
        // of a solid turned about their own normal are the population; their frames stay
        // world-spoken (`face_sketch_frame` verifies instead of declining) and their
        // predicates keep the exact roads.
        //
        // `Known` + narrow name only: a `Through` plane's truth is vertex handles in a
        // separate intern table, and a source that already carries a motion keeps today's
        // recorded path (stage 1 — under it a fixed plane never gains a history, so a
        // second turn about the same normal restates again). A mirror never restates here:
        // an improper motion's parity interactions deserve their own measured stage.
        let invariant = matches!(
            &src_truth,
            nacre_topo::SurfaceTruth::Plane {
                points: nacre_topo::PlanePoints::Known(_),
                motion: None,
            }
        ) && motion.rigid().is_some_and(|iso| {
            model
                .surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| iso.fixes_plane(c))
        });
        let new_motion = if invariant {
            None
        } else {
            moved_surface_motion(model, s, motion, exact, &mut surf_rot)
        };
        // ★★★★★ **Only the points move.** The image's canonical name is derived from them by
        // `Model::push_surface_with_points`, so there is no second description to keep in step —
        // the agreement between points and coefficients was checked across the suite before the
        // coefficient parameter went away (83,883 pushes, 0 disagreements) and is now structural.
        //
        // A recorded node states its plane **before** the motion, and the image's base is the
        // source's base — `moved_surface_motion` chains from the source's own leaf for the same
        // reason — so the triple is inherited verbatim. With nothing recorded the motion kept
        // everything exact, and the points are carried in the world with it; the transport is
        // the very function `motion_is_exact` probed, so it cannot fail here.
        // ★★ A `Through` plane's statement is *handles*, and handles do not move. Its image keeps
        // the same three vertices and gains the node — the definition composes as "the plane
        // through those, then this motion". Transporting them is not an option (they may not even
        // belong to this solid), and leaving them without a node would move the cache while the
        // truth stayed put. `points_move` refuses the no-node path for exactly this reason, so
        // `new_motion` is always `Some` here.
        let (new_s, flipped) = match (moved, &src_truth) {
            (
                Surface::Plane(pl),
                nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Known(p),
                    motion: src_m,
                },
            ) => {
                // ★ `invariant` must gate first: there `new_motion == *src_m == None`, and
                // the else arm's expect names a probe (`motion_is_exact`) that never ran on
                // this road — an irrational turn would panic in `transport_points`. The
                // statement is the image's verbatim, so there is nothing to transport.
                let carried = if invariant || new_motion != *src_m {
                    *p
                } else {
                    transport_points(motion, *p).expect("probed by motion_is_exact")
                };
                model.push_plane(pl, carried, new_motion)
            }
            (
                Surface::Plane(pl),
                nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Through(vs),
                    ..
                },
            ) => {
                // ★ `new_motion` may be `None`, and only for a motion that moves nothing —
                // `Copy` is `transform_solid` with the identity. Then the plane is unchanged, the
                // same statement re-pushed interns back onto the source handle, and that is the
                // right answer. `points_move` refusing the transporting path is what rules out
                // the other case, where a node is missing because the walk thought it could carry
                // points that do not exist.
                let out = model.push_plane_through(pl, *vs, new_motion);
                debug_assert!(
                    new_motion.is_some() || out.0 == s,
                    "a Through plane gained no node yet changed handle — its truth would be \
                     describing the plane it used to be"
                );
                out
            }
            (Surface::Cylinder(cy), nacre_topo::SurfaceTruth::Cylinder { def, motion: src_m }) => {
                // The same fork as the `Known` plane above, minus the invariant road (an
                // invariant-cylinder restatement — a turn about its own axis — is deliberately
                // deferred; the condition is narrower than a plane's because `ref_dir` turns).
                // A recorded node states the def **before** the motion → carried verbatim;
                // nothing recorded means the motion kept everything exact → the def rides the
                // very transport the probe checked.
                let carried = if new_motion != *src_m {
                    def.clone()
                } else {
                    transport_cylinder(motion, def).expect("probed by motion_is_exact")
                };
                (model.push_cylinder(cy, carried, new_motion), false)
            }
            (Surface::Plane(_), nacre_topo::SurfaceTruth::Cylinder { .. })
            | (Surface::Cylinder(_), nacre_topo::SurfaceTruth::Plane { .. }) => {
                unreachable!("a surface's cache and truth cannot disagree about its kind")
            }
        };
        surf_map.insert(s, new_s);
        surf_flip.insert(s, flipped);
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

    // (The old pass 2 — moving curves — died with `Store<Curve>` (S8): the moved edge's
    // curve now derives from its moved carriers and endpoints in pass 4's `push_edge`.)

    // Pass 3 — vertices (dedup): the definition's handles re-pointed onto the moved surfaces,
    // the coordinate moved, the measured tolerance kept only when nothing was re-realized.
    let mut vert_order: Vec<Handle<Vertex>> = Vec::new();
    let mut vert_seen: HashSet<Handle<Vertex>> = HashSet::new();
    for &eh in &edge_order {
        for vh in model.edges.get(eh).vertices {
            if vert_seen.insert(vh) {
                vert_order.push(vh);
            }
        }
    }
    let mut vert_map: HashMap<Handle<Vertex>, Handle<Vertex>> = HashMap::new();
    for &vh in &vert_order {
        // `defs_are_remappable` cleared every definition handle before the walk began.
        let remap = |s: Handle<Surface>| {
            *surf_map
                .get(&s)
                .expect("a definition's surface must be a face surface of the solid")
        };
        let def = match model.vertices.get(vh).def {
            VertexDef::ThreePlane(planes) => VertexDef::ThreePlane(planes.map(remap)),
            VertexDef::OnSeam(pair) => VertexDef::OnSeam(pair.map(remap)),
            // ★ `.map(remap)` alone would be wrong here: pass 1 issues new surface handles in
            // face-traversal order, so the two planes' handle order can invert. The pair is
            // re-sorted (the stored-ascending invariant) — and a swap flips the canonical
            // line direction ℓ = n₁×n₂, so `root` must toggle with it or `Lo` silently
            // starts naming the other point.
            VertexDef::Branch {
                planes: [p0, p1],
                cylinder,
                root,
            } => {
                let (a, b) = (remap(p0), remap(p1));
                let (planes, root) = if b.index() < a.index() {
                    ([b, a], root.flipped())
                } else {
                    ([a, b], root)
                };
                VertexDef::Branch {
                    planes,
                    cylinder: remap(cylinder),
                    root,
                }
            }
        };
        let coord = motion.point(model.vertex_point(vh));
        // Tolerance rule (S7, letter-preserving): an exact move of a history-less solid used to
        // keep `Discovered { tol }` verbatim (`remap_origin`); a recorded move used to demote to
        // `Moved` (checker epsilon — now `None`).
        let tol = if keeps_tol {
            model.vertex_tol(vh)
        } else {
            None
        };
        vert_map.insert(vh, model.push_vertex(def, coord, tol));
    }

    // Pass 4 — edges (carrier/vertex handles remapped; the curve cache derives from them).
    let mut edge_map: HashMap<Handle<Edge>, Handle<Edge>> = HashMap::new();
    for &eh in &edge_order {
        let e = *model.edges.get(eh);
        // The carriers move with the surfaces (pass 1 mapped every reachable one, so the
        // lookups cannot miss); `push_edge` re-canonicalizes the pair. A rigid image of a
        // non-degenerate edge cannot degenerate, so the `None` is unreachable in practice —
        // mapped to the honest reject rather than a panic all the same.
        let new_e = model
            .push_edge(
                [surf_map[&e.surfaces[0]], surf_map[&e.surfaces[1]]],
                e.vertices.map(|v| vert_map[&v]),
            )
            .ok_or(OpError::DegenerateGeometry)?;
        edge_map.insert(eh, new_e);
    }

    // Pass 5 — faces (loops rebuilt onto the new edges; orientation unchanged).
    // A reflection rewinds every loop (see the fn doc); a rigid motion keeps the winding.
    let rewind = motion.reverses_orientation();
    let map_loop = |lp: &Loop| {
        let mapped = Loop {
            half_edges: lp
                .half_edges
                .iter()
                .map(|he| HalfEdge {
                    edge: edge_map[&he.edge],
                    forward: he.forward,
                })
                .collect(),
        };
        if rewind { mapped.reversed() } else { mapped }
    };
    let mut face_map: HashMap<Handle<Face>, Handle<Face>> = HashMap::new();
    for &fh in &face_order {
        let face = model.faces.get(fh).clone();
        let new_f = Face {
            surface: surf_map[&face.surface],
            outer: map_loop(&face.outer),
            inner: face.inner.iter().map(&map_loop).collect(),
            orientation: if surf_flip[&face.surface] {
                face.orientation.flipped()
            } else {
                face.orientation
            },
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
    Ok(model.push_solid(new_solid))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_scalar::Angle;

    /// **"This step kept it exact" never excuses a chain from recording it.**
    ///
    /// A motion that leaves every coordinate on an exact `f64` needs no node *of its own* — but if
    /// the datum already has a history, the chain has to keep reproducing it, and a chain missing a
    /// link describes the datum as it was before that link. Silently.
    ///
    /// Asserted on the rule rather than on a symptom. The translation half is what the rule was
    /// found to be wrong on; the reflection half is the case that made it load-bearing, since a
    /// reflection in a dyadic plane preserves exactness *and* can carry a history.
    #[test]
    fn an_exact_motion_is_still_recorded_once_there_is_a_history() {
        let mut m = Model::new();
        let root = m.push_motion(
            Motion::Rotate {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
            },
            None,
        );
        let place = Isometry::translation([Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)]);
        let flip = Xform::mirror(Axis::X, Rat::from_int(0));

        for (what, motion, kind) in [
            ("translation", &Xform::Rigid(&place), "Translate"),
            ("reflection", &flip, "Mirror"),
        ] {
            // No history: an exact motion records nothing, and the coordinates stay the truth.
            assert_eq!(
                chain_motion(&mut m, None, motion, true),
                None,
                "an exact {what} over no history records nothing"
            );
            // With a history: recorded anyway, and the new leaf hangs off the old one.
            let leaf = chain_motion(&mut m, Some(root), motion, true)
                .unwrap_or_else(|| panic!("a {what} over a history is always recorded"));
            assert_ne!(leaf, root);
            assert_eq!(m.motion(leaf).parent, Some(root));
            // And it records **this** motion — a chain that keeps the wrong kind of node
            // reproduces the wrong datum just as silently as one that keeps none.
            let got = match m.motion(leaf).motion {
                Motion::Rotate { .. } => "Rotate",
                Motion::Translate { .. } => "Translate",
                Motion::Mirror { .. } => "Mirror",
                Motion::Frame { .. } => "Frame",
            };
            assert_eq!(got, kind, "a {what} records a {kind} node");
        }
    }
    /// ★★★ S6a: **an exact move whose rational point transport would overflow records a node
    /// instead of dropping the points.** The two exactness conditions are independent — every
    /// f64 here lands exactly (`t = 2⁻³⁰` on unit-scale corners), while one surface's stored
    /// triple has a `5⁴²` denominator, so `q + t` needs `lcm(5⁴², 2³⁰) ≈ 2.4e38 > i128`. The
    /// old behaviour kept the no-node path and silently pushed the moved surface point-less —
    /// the very population `Inexact` grows from; the `motion_is_exact` probe now puts the whole
    /// solid on the recorded path, and the original triple survives verbatim as the pre-motion
    /// truth.
    #[test]
    fn an_overflowing_exact_move_records_a_node_and_keeps_the_points() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        // Make one surface's triple adversarially deep: a 5⁴² denominator no dyadic translation
        // can share. (The map is data here — the mechanism under test reads it, nothing else.)
        let fh = m.shells.get(m.solids.get(s).outer).faces[0];
        let surf = m.faces.get(fh).surface;
        let deep = Rat::new(1, 5i128.pow(42)).unwrap();
        let pts = [
            [deep, Rat::from_int(0), Rat::from_int(0)],
            [Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)],
            [Rat::from_int(0), Rat::from_int(1), Rat::from_int(0)],
        ];
        m.set_plane_points_for_test(surf, pts);
        // Fixture qualification: the f64 side is exact, the rational side overflows.
        let t = Rat::new(1, 1 << 30).unwrap();
        assert_eq!(
            1.0f64 + t.to_f64(),
            (Rat::from_int(1).checked_add(t)).unwrap().to_f64()
        );
        assert!(deep.checked_add(t).is_none(), "the transport must overflow");

        // The z component keeps the doctored z-plane off the invariant branch (a purely-x
        // translation *fixes* it, and a fixed plane is restated with no transport at all —
        // this test is about the transport overflowing). `0 + 1` and `1 + 1` are exact, so
        // the f64-side qualification above still holds for every corner.
        let iso = Isometry::translation([t, Rat::from_int(0), Rat::from_int(1)]);
        let moved = transform_solid(&mut m, s, &Xform::Rigid(&iso)).unwrap();
        m.rebuild_adjacency();
        let moved_surf = m
            .shells
            .get(m.solids.get(moved).outer)
            .faces
            .iter()
            .map(|&f| m.faces.get(f).surface)
            .find(|&s2| {
                matches!(
                    m.surface_truth(s2),
                    nacre_topo::SurfaceTruth::Plane {
                        points: nacre_topo::PlanePoints::Known(p),
                        ..
                    } if *p == pts
                )
            })
            .expect("the deep triple must survive the move verbatim (pre-motion truth)");
        assert!(
            matches!(
                m.surface_truth(moved_surf),
                nacre_topo::SurfaceTruth::Plane {
                    motion: Some(_),
                    ..
                }
            ),
            "the overflow must force the recorded path, not drop the points"
        );
    }

    /// ★★ S6b: **a moved cylinder records its history instead of silently degrading.** The old
    /// `SurfaceDef` path read the lateral surface as `Constructed`-without-points, so a rotated
    /// cylinder's move demoted it to `Inexact` — reachable in production, pinned by nothing.
    /// The truth variant has a motion slot of its own, and this is it working.
    /// ★ S7 C2 positive control: the widened gate **bites** — a walkable solid whose vertex
    /// definition names a foreign surface is rejected with `OriginNotOnSolid`. The old
    /// (`Discovered`-only) gate passed this fixture (the vertex is `Constructed`), so this is
    /// specifically the widening's teeth; without it, "the gate fired zero times across the
    /// suite" would be indistinguishable from "the gate checks nothing".
    #[test]
    fn a_foreign_definition_is_rejected() {
        let mut m = Model::new();
        // A plane of this solid's own...
        let (own, _) = m.push_plane(
            Plane::from_point_normal(
                Point3::from_array([0.0, 0.0, 9.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            [
                [Rat::from_int(0), Rat::from_int(0), Rat::from_int(9)],
                [Rat::from_int(1), Rat::from_int(0), Rat::from_int(9)],
                [Rat::from_int(0), Rat::from_int(1), Rat::from_int(9)],
            ],
            None,
        );
        // ...and a vertex whose definition names the world seeds — surfaces this solid's face
        // set does not contain.
        let foreign = VertexDef::ThreePlane([
            m.world_plane(Axis::Z),
            m.world_plane(Axis::X),
            m.world_plane(Axis::Y),
        ]);
        let mk_v =
            |m: &mut Model, x: f64| m.push_vertex(foreign, Point3::from_array([x, 0.0, 9.0]), None);
        let v0 = mk_v(&mut m, 0.0);
        let v1 = mk_v(&mut m, 1.0);
        let e = m
            .push_edge([own, m.world_plane(Axis::Z)], [v0, v1])
            .unwrap();
        let f = m.faces.push(Face {
            surface: own,
            outer: Loop {
                half_edges: vec![HalfEdge {
                    edge: e,
                    forward: true,
                }],
            },
            inner: vec![],
            orientation: nacre_topo::Orientation::Forward,
        });
        let sh = m.shells.push(Shell { faces: vec![f] });
        let s = m.push_solid(Solid {
            outer: sh,
            cavities: vec![],
        });
        let zero = Isometry::translation([Rat::from_int(0); 3]);
        assert_eq!(
            transform(&mut m, s, &zero),
            Err(OpError::OriginNotOnSolid),
            "a foreign definition must be rejected before the walk"
        );
    }

    /// ★★ S7: a moved cylinder's seam vertices carry `OnSeam` re-pointed at the **twin's own**
    /// surfaces — the carrier pair moves with the solid, like every other definition.
    #[test]
    fn a_moved_cylinders_seam_defs_repoint_to_the_twin() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::from_array([1.0, 2.0, 0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.5,
            3.0,
        );
        m.rebuild_adjacency();
        let iso = Isometry::rotation(nacre_scalar::Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: nacre_scalar::Angle::from_deg(Rat::from_int(31)).expect("angle"),
        });
        let out = transform(&mut m, s, &iso).expect("rotate the cylinder");
        let mut twin_surfs = std::collections::HashSet::new();
        for &fh in &m.shells.get(m.solids.get(out).outer).faces {
            twin_surfs.insert(m.faces.get(fh).surface);
        }
        let mut seams = 0;
        for &fh in &m.shells.get(m.solids.get(out).outer).faces {
            let face = m.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in &m.edges.get(he.edge).vertices {
                        if let VertexDef::OnSeam(pair) = m.vertices.get(vh).def {
                            seams += 1;
                            assert!(
                                pair.iter().all(|c| twin_surfs.contains(c)),
                                "a moved seam def must name the twin's own surfaces"
                            );
                        }
                    }
                }
            }
        }
        assert!(seams > 0, "the sweep saw no seam vertices at all");
    }

    #[test]
    fn a_rotated_cylinder_records_its_motion() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.0,
            2.0,
        );
        m.rebuild_adjacency();
        let iso = Isometry::rotation(nacre_scalar::Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: nacre_scalar::Angle::from_deg(Rat::from_int(31)).expect("angle"),
        });
        let turned = transform_solid(&mut m, s, &Xform::Rigid(&iso)).unwrap();
        m.rebuild_adjacency();
        let lateral = m
            .shells
            .get(m.solids.get(turned).outer)
            .faces
            .iter()
            .map(|&f| m.faces.get(f).surface)
            .find(|&su| matches!(m.surface(su), Surface::Cylinder(_)))
            .expect("a cylinder keeps its lateral face");
        // An inexact turn records a node, and the def is carried **verbatim** — the recorded
        // node states its cylinder before the motion (the plane rule, unchanged by M6-0).
        match m.surface_truth(lateral) {
            nacre_topo::SurfaceTruth::Cylinder {
                def,
                motion: Some(_),
            } => {
                let z = Rat::from_int(0);
                assert_eq!(def.origin(), [z; 3], "pre-motion statement, verbatim");
                assert_eq!(def.dir(), [z, z, Rat::from_int(1)]);
            }
            other => panic!("a moved cylinder's truth must carry the motion, got {other:?}"),
        }
        // ★ **The measurement this file was missing.** Nothing here used to call `validate`, so
        // the *moved* cylinder population never met the net — and the net's newest rule reads a
        // cap's rim circle against its plane, two caches that a motion carries together. If they
        // ever stopped travelling together this is where it would show.
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
    }

    /// ★ M6-1: the remap of a `Branch` definition is not `.map(remap)` — pass 1 issues new
    /// surface handles in face-traversal order, so the two planes' handle order can invert,
    /// and re-sorting flips the canonical line direction ℓ = n₁×n₂, so the root must toggle
    /// with the swap or `Lo` silently names the other point.
    ///
    /// The fixture makes the swap *actually happen*: the branch planes are the cylinder's
    /// bottom cap (a fresh handle after a z-translation) and the world x = 0 seed (invariant
    /// under that translation — it keeps handle 1), so `[cap(0), x0(1)]` remaps to
    /// `[N, 1] → sorted [1, N]` — swapped. Geometry agrees with the toggle: with planes
    /// `[cap, x0]` the canonical ℓ is +y and the y = −2 point is `Lo`; with `[x0, cap′]` ℓ
    /// is −y and that same point is `Hi`.
    #[test]
    fn a_branch_definition_swap_toggles_its_root() {
        use nacre_topo::QuadRoot;
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        m.rebuild_adjacency();
        let shell = m.solids.get(s).outer;
        let faces = m.shells.get(shell).faces.clone();
        let lateral = faces
            .iter()
            .map(|&f| m.faces.get(f).surface)
            .find(|&su| matches!(m.surface(su), Surface::Cylinder(_)))
            .expect("lateral");
        let bottom = m.world_plane(Axis::Z); // the z = 0 cap interned onto the world seed
        let x0 = m.world_plane(Axis::X);
        assert!(bottom.index() < x0.index(), "the fixture's premise");
        // The two branch points of {z = 0} ∧ {x = 0} against the cylinder: (0, ∓2, 0).
        // Canonical normals (0,0,1) × (1,0,0) = +y, so y = −2 is the smaller parameter: Lo.
        let v_lo = m.push_vertex(
            VertexDef::Branch {
                planes: [bottom, x0],
                cylinder: lateral,
                root: QuadRoot::Lo,
            },
            Point3::from_array([0.0, -2.0, 0.0]),
            None,
        );
        let v_hi = m.push_vertex(
            VertexDef::Branch {
                planes: [bottom, x0],
                cylinder: lateral,
                root: QuadRoot::Hi,
            },
            Point3::from_array([0.0, 2.0, 0.0]),
            None,
        );
        // Wire the vertices into the solid (a franken-face on the x = 0 seed): transform
        // remaps only what its face walk reaches, and `defs_are_remappable` requires every
        // carrier among the face surfaces.
        let edge = m
            .push_edge([bottom, x0], [v_lo, v_hi])
            .expect("a line through distinct endpoints");
        let franken = m.faces.push(Face {
            surface: x0,
            outer: Loop {
                half_edges: vec![HalfEdge {
                    edge,
                    forward: true,
                }],
            },
            inner: vec![],
            orientation: nacre_topo::Orientation::Forward,
        });
        let mut new_faces = faces.clone();
        new_faces.push(franken);
        let sh = m.shells.push(Shell { faces: new_faces });
        let franken_solid = m.push_solid(Solid {
            outer: sh,
            cavities: vec![],
        });
        m.live_solids.retain(|&x| x == franken_solid);

        let iso = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
        let moved = transform_solid(&mut m, franken_solid, &Xform::Rigid(&iso)).expect("moves");

        // Find the two branch vertices of the moved solid and read their defs.
        let mut seen = Vec::new();
        let solid = m.solids.get(moved).clone();
        for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                let face = m.faces.get(fh).clone();
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        for &vh in m.edges.get(he.edge).vertices.iter() {
                            if let VertexDef::Branch { planes, root, .. } = m.vertices.get(vh).def {
                                seen.push((m.vertex_point(vh).as_array(), planes, root));
                            }
                        }
                    }
                }
            }
        }
        seen.sort_by(|a, b| a.0[1].partial_cmp(&b.0[1]).unwrap());
        seen.dedup_by_key(|e| e.0[1] as i64);
        assert_eq!(seen.len(), 2, "both branch vertices survive the move");
        for (p, planes, root) in &seen {
            assert!(
                planes[0].index() < planes[1].index(),
                "stored order stays ascending"
            );
            // x = 0 kept its seed handle (invariant restatement); the cap moved to a fresh
            // one — so the pair swapped, and the root must have toggled with it.
            assert_eq!(planes[0], x0, "the surviving seed now sorts first");
            let want = if p[1] < 0.0 {
                QuadRoot::Hi
            } else {
                QuadRoot::Lo
            };
            assert_eq!(
                *root, want,
                "at y = {}: ℓ flipped to −y, so the root names the same point only if it \
                 toggled",
                p[1]
            );
        }
    }

    #[test]
    fn an_exact_turn_transports_a_cylinders_truth_instead_of_recording() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::from_array([0.5, -1.25, 2.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.5,
            2.5,
        );
        m.rebuild_adjacency();
        let iso = Isometry::rotation(nacre_scalar::Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: nacre_scalar::Angle::from_deg(Rat::from_int(90)).expect("angle"),
        });
        let turned = transform_solid(&mut m, s, &Xform::Rigid(&iso)).unwrap();
        m.rebuild_adjacency();
        let lateral = m
            .shells
            .get(m.solids.get(turned).outer)
            .faces
            .iter()
            .map(|&f| m.faces.get(f).surface)
            .find(|&su| matches!(m.surface(su), Surface::Cylinder(_)))
            .expect("a cylinder keeps its lateral face");
        // A 90°-family turn is exact: nothing is recorded, and the def rides the very
        // transport the probe checked — origin through `point_rat`, directions through
        // `dir_rat` (the pivot cancels), radius invariant. (x, y) ↦ (−y, x).
        let d = |x: f64| Rat::from_decimal(x).expect("decimal");
        match m.surface_truth(lateral) {
            nacre_topo::SurfaceTruth::Cylinder { def, motion: None } => {
                assert_eq!(def.origin(), [d(1.25), d(0.5), d(2.0)]);
                assert_eq!(def.dir(), [d(0.0), d(0.0), d(1.0)]);
                assert_eq!(
                    def.ref_dir(),
                    [d(1.0), d(0.0), d(0.0)],
                    "seam turned with it"
                );
                assert_eq!(def.radius(), d(1.5), "radius is rigid-invariant");
            }
            other => panic!("an exact turn must transport, not record — got {other:?}"),
        }
    }
}
