//! Rigid-body transform of a solid (design overhaul stage 1): copy a solid under an isometry,
//! remapping every surface/curve and carrying each vertex's exact `Origin` forward so a rotated
//! operand stays exactly defined. [`copy`] is the same walk with no motion at all.

use crate::OpError;
use nacre_geom::{AxisMirror, Cylinder, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Motion, MotionNode, Origin, Shell, Solid, Vertex, VertexDef,
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
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !origins_are_remappable(model, solid) {
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
/// **The carried `Origin::Discovered { tol }` is exact here, not provisional.** `remap_origin`
/// keeps `tol` unchanged, which for a *rotated* move is a stage-1 debt that tol transport must
/// later repay; for a copy the point and its three planes are bit-identical, so the measured
/// residual cannot have changed. Tol transport will have nothing to fix here.
///
/// A non-live input is rejected rather than resurrected: reusing a superseded handle is a caller
/// bug, and letting it succeed would hide it.
pub(crate) fn copy(model: &mut Model, solid: Handle<Solid>) -> Result<Handle<Solid>, OpError> {
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !origins_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let zero = Isometry::translation([Rat::from_int(0); 3]);
    // `push_solid` registers the twin as live; the input is *not* retained away — that missing
    // line is the whole difference from `transform`.
    transform_solid(model, solid, &Xform::Rigid(&zero))
}

/// Whether every `Discovered` definition in `solid` names a surface the walk will remap — that is,
/// a face surface of this solid. [`remap_origin`] cannot proceed otherwise, and a kernel must
/// decline rather than abort ("honest-reject > silent-wrong"; `assemble_fuse_cut` takes the same
/// line at its own naming failure), so the two entry points check first and reject.
///
/// The invariant is expected to hold — a boolean's result faces and its `VertexDef::ThreePlane`
/// both name the *plane class representative* surface — but "expected" is not "proved" for severed
/// results and post-cleaning face sets, and [`copy`] makes this path routine rather than rare.
fn origins_are_remappable(model: &Model, solid: Handle<Solid>) -> bool {
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
    let named = |origin: &Origin| match origin {
        Origin::Discovered {
            definition: VertexDef::ThreePlane(planes),
            ..
        } => planes.iter().all(|s| surfs.contains(s)),
        _ => true,
    };
    for &sh in &shells {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    let edge = model.edges.get(he.edge);
                    for vh in edge.vertices.iter() {
                        if !named(&model.vertices.get(*vh).origin) {
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
    if !model.live_solids.contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !origins_are_remappable(model, solid) {
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
        Xform::Mirror { axis, offset, .. } => reflected(p.as_array()[axis_index(*axis)], *offset),
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
        let nacre_topo::SurfaceTruth::Plane {
            points: nacre_topo::PlanePoints::Known(p),
            ..
        } = model.surface_truth(s)
        else {
            return true; // a cylinder carries no points
        };
        // The very function pass 1 will transport with — sharing it is what makes the probe's
        // promise ("this will not overflow") structural rather than a parallel re-derivation.
        transport_points(motion, *p).is_some()
    };
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            if let Surface::Plane(pl) = model.surface(face.surface) {
                if !all(pl.origin()) {
                    return false;
                }
            }
            if !points_move(face.surface) {
                return false;
            }
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in model.edges.get(he.edge).vertices.iter() {
                        if !all(model.vertices.get(vh).point) {
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

/// The motion history a moved surface's image carries — the surface twin of the vertex
/// `Origin` rules, in S6b form (the old `SurfaceDef` table collapsed to one field):
///
/// | source motion | motion records a node | motion records nothing |
/// |---|---|---|
/// | `None` (world) | `Some(leaf)` — the points stay pre-motion | `None` (points transported exactly) |
/// | `Some(m)` | `Some(leaf)`, hanging off `m` — replaying from the root applies every motion once | unreachable (the exceptions lapse once a history exists) |
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
    let parent = match model.surface_truth(src) {
        nacre_topo::SurfaceTruth::Plane { motion, .. }
        | nacre_topo::SurfaceTruth::Cylinder { motion } => *motion,
    };
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

fn axis_index(axis: Axis) -> usize {
    match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
    }
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
            m: AxisMirror::new(axis_index(axis), offset.to_f64()).expect("axis index is 0..3"),
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

/// A vertex/edge `Origin` with any `Discovered` `ThreePlane` definition remapped
/// onto the moved surfaces. Constructed stays constructed; the `tol` is unchanged
/// (a rigid move; stage 1 records tol but never judges on it — for [`copy`] it is exact,
/// since the point and its planes are bit-identical). The `expect` below is unreachable:
/// both callers run [`origins_are_remappable`] first and reject when it fails.
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
        Origin::Moved { .. } => origin,
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

    // The forest nodes this transform appends (§CIP ⑦), one shared leaf named by every moved
    // vertex. `chain_motion` decides what is worth recording: a zero translation, a 90°-family
    // rotation of a still-exact datum, and an exactness-preserving reflection record nothing, and
    // every exception lapses once the solid already has a history — then every motion is recorded
    // or the chain would not reproduce the result.
    let input_leaf = solid_motion(model, solid);
    // Decided once, for the whole solid — see [`motion_is_exact`].
    let exact = motion_is_exact(model, solid, motion);
    let move_node: Option<Handle<MotionNode>> = chain_motion(model, input_leaf, motion, exact);

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
        let new_motion = moved_surface_motion(model, s, motion, exact, &mut surf_rot);
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
        let points = match &src_truth {
            nacre_topo::SurfaceTruth::Plane {
                points: nacre_topo::PlanePoints::Known(p),
                motion: src_m,
            } => {
                if new_motion != *src_m {
                    Some(*p)
                } else {
                    Some(transport_points(motion, *p).expect("probed by motion_is_exact"))
                }
            }
            nacre_topo::SurfaceTruth::Cylinder { .. } => None,
        };
        let (new_s, flipped) = match (moved, points) {
            (Surface::Plane(pl), Some(p)) => model.push_plane(pl, p, new_motion),
            (Surface::Cylinder(cy), _) => (model.push_cylinder(cy, new_motion), false),
            (Surface::Plane(_), None) => unreachable!("a plane's truth always carries points"),
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

    // Pass 3 — vertices (dedup, moved point + Origin preserved/remapped).
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
        let v = *model.vertices.get(vh);
        let new_v = Vertex {
            point: motion.point(v.point),
            // A recorded motion marks the vertex `Moved`; `base` is the **root** — the
            // non-`Moved` (Constructed/Discovered) ancestor whose exact definition the chain
            // moves. A fresh move's input is itself the root; a re-move chases one hop to the
            // input's own root (the invariant keeps `base` pointing at a root, never at another
            // `Moved` vertex, so a replay never applies the same motion twice). No node → remap.
            //
            // **A reflection needs nothing more.** It used to conjugate the chain onto a mirrored
            // copy of the root and replay the coordinate from there, because the chain had no way
            // to say "and then reflect"; now it does, `Pt3::mirror` walks the same `2c − x` the
            // producer just walked, and the root stays the root.
            origin: match move_node {
                Some(node) => {
                    let base = match v.origin {
                        Origin::Moved { base, .. } => base,
                        _ => vh,
                    };
                    Origin::Moved { base, motion: node }
                }
                None => remap_origin(v.origin, &surf_map),
            },
            // ★ **One rule, whatever the `Origin` is: map every plane, or record nothing.**
            // `remap_origin` above `expect`s its lookups because `origins_are_remappable` cleared
            // them first — but that gate only inspects `Discovered` vertices. Extending its reach
            // to cover this field would either panic here or start rejecting transforms that work
            // today, and both would make this commit change answers. A definition that cannot be
            // re-pointed is simply absent, and the coverage count says how often.
            definition: v.definition.and_then(|VertexDef::ThreePlane(planes)| {
                let mut out = [planes[0]; 3];
                for (o, s) in out.iter_mut().zip(planes) {
                    *o = *surf_map.get(&s)?;
                }
                Some(VertexDef::ThreePlane(out))
            }),
        };
        vert_map.insert(vh, model.vertices.push(new_v));
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

/// The shared rotation-forest leaf that every boundary vertex of a rotated `solid`
/// names (`None` if the solid is not rotated) — the input side of the re-rotation
/// decision in [`transform_solid`]. All boundary vertices share one leaf, and that
/// invariant is `debug_assert`ed here (fail-loud if a future change ever produces a
/// solid with mixed rotation provenance).
///
/// Two things keep it true. A `Transform` rotates a solid uniformly, so its own output is
/// uniform. And a *boolean* output cannot mix provenance either — not because rotated inputs
/// are refused (that guard, `ROTATED_UNSUPPORTED`, was retired when rotated booleans went live)
/// but because `assemble_fuse_cut` names every result vertex through `Node::Seam`, the enum's
/// only variant: even a corner that survived untouched is rebuilt as a seam vertex, so a result
/// is uniformly `Discovered` and carries no `Rotated` at all.
fn solid_motion(model: &Model, solid: Handle<Solid>) -> Option<Handle<MotionNode>> {
    let sh = model.solids.get(solid).outer;
    let mut seen: Option<Option<Handle<MotionNode>>> = None;
    for &fh in &model.shells.get(sh).faces {
        for he in &model.faces.get(fh).outer.half_edges {
            for &vh in &model.edges.get(he.edge).vertices {
                let leaf = match model.vertices.get(vh).origin {
                    Origin::Moved { motion, .. } => Some(motion),
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

        let iso = Isometry::translation([t, Rat::from_int(0), Rat::from_int(0)]);
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
        assert!(
            matches!(
                m.surface_truth(lateral),
                nacre_topo::SurfaceTruth::Cylinder { motion: Some(_) }
            ),
            "a moved cylinder's truth must carry the motion, not degrade"
        );
    }
}
