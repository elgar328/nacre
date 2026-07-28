//! Rigid-body transform of a solid (design overhaul stage 1): copy a solid under an isometry,
//! remapping every surface/curve and carrying each vertex's exact `Origin` forward so a rotated
//! operand stays exactly defined. [`copy`] is the same walk with no motion at all.

use crate::OpError;
use nacre_geom::{AxisMirror, Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Angle, Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Motion, MotionNode, Origin, Shell, Solid, SurfaceDef,
    Vertex, VertexDef,
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
    let out = transform_solid(model, solid, &Xform::Rigid(isometry), None)?;
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
    transform_solid(model, solid, &Xform::Rigid(&zero), None)
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
                    if !named(&edge.origin) {
                        return false;
                    }
                    for vh in edge.bounds.iter().flatten() {
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
    let m = AxisMirror::new(axis_index(axis), offset.to_f64()).expect("axis index is 0..3");
    let out = transform_solid(model, solid, &Xform::Mirror(m), Some((axis, offset)))?;
    model.live_solids.retain(|&s| s != solid);
    Ok(out)
}

/// What a reflection needs in order to carry a rotated solid's exact definitions: the mirror
/// itself plus memo tables, so a chain shared by every boundary vertex is conjugated once.
struct Conjugation {
    mirror: AxisMirror,
    axis: Axis,
    offset: Rat,
    nodes: HashMap<Handle<MotionNode>, Handle<MotionNode>>,
    bases: HashMap<Handle<Vertex>, Handle<Vertex>>,
    surfaces: HashMap<Handle<Surface>, Handle<Surface>>,
}

/// The reflected image of a rotation chain, as a chain.
///
/// `M ∘ R = (M R M⁻¹) ∘ M`, and for an axis-aligned mirror and the kernel's X/Y/Z rotation axes
/// the conjugate `M R M⁻¹` is again a rotation about the **same axis**, through the **mirrored
/// pivot**, by the **negated angle** — except when the mirror plane's normal *is* the rotation
/// axis, where the two commute and the angle is unchanged (`diag(1,1,−1)` commutes with `Rot_z`).
///
/// Every step is exact: the pivot is `[Rat; 3]` and the angle rational degrees, so a reflection —
/// unlike a rotation — introduces no irrational value at all. `None` only if the rational
/// arithmetic overflows.
fn conjugate_chain(
    model: &mut Model,
    leaf: Handle<MotionNode>,
    c: &mut Conjugation,
) -> Option<Handle<MotionNode>> {
    if let Some(&h) = c.nodes.get(&leaf) {
        return Some(h);
    }
    let node = *model.motions.get(leaf);
    let parent = match node.parent {
        Some(p) => Some(conjugate_chain(model, p, c)?),
        None => None,
    };
    let i = axis_index(c.axis);
    let motion = match node.motion {
        Motion::Rotate {
            axis,
            mut point,
            angle,
        } => {
            // The **pivot** is a point, so it reflects as `2c − p`.
            point[i] = c
                .offset
                .checked_mul(Rat::from_int(2))?
                .checked_sub(point[i])?;
            // `Angle` keeps `0 ≤ θ < 360`, so the negation of 0 is 0, not 360.
            let angle = if axis == c.axis || angle.deg() == Rat::from_int(0) {
                angle
            } else {
                Angle::from_deg(Rat::from_int(360).checked_sub(angle.deg())?)?
            };
            Motion::Rotate { axis, point, angle }
        }
        // An **offset is a vector, not a point** — `M ∘ T(t) ∘ M⁻¹ = T(t')` negates the mirrored
        // component and leaves the rest, with no `2c` term. The pivot formula above looks like it
        // would fit (both are `[Rat; 3]`) and would be wrong.
        Motion::Translate { mut offset } => {
            offset[i] = Rat::from_int(0).checked_sub(offset[i])?;
            Motion::Translate { offset }
        }
    };
    let h = model.push_motion(motion, parent);
    c.nodes.insert(leaf, h);
    Some(h)
}

/// The provenance a moved surface inherits — the surface twin of the vertex `Origin` rules.
///
/// `rot_node` is the forest node this transform recorded (`None` for a motion that keeps
/// coefficients exact: a translation, a 90°-family rotation, a reflection). The three source
/// definitions map as:
///
/// | source | inexact rotation (`rot_node`) | exact motion |
/// |---|---|---|
/// | `Constructed` | `Rotated { witness: this face's pre-rotation triangle, rot_node }` | `Constructed`, or the conjugated chain under a mirror |
/// | `Rotated { witness, .. }` | `Rotated { witness, rot_node }` — the node's parent *is* the source's leaf, so the chain already carries the earlier rotation | identity (`copy`): unchanged; mirror: conjugate the chain over the mirrored witness; **translation: `Inexact`** |
/// | `Inexact` | `Inexact` | `Inexact` |
///
/// The one `Inexact` producer is a translation of an already-rotated surface: the forest holds
/// rotations only, so `R` then `T` has no node to name. That is the same limit
/// `Pt3Error::TranslateInterleaved` already reports for vertices, and it is a reject, not a
/// silent rounding.
fn moved_surface_def(
    model: &mut Model,
    src: Handle<Surface>,
    face: Handle<Face>,
    motion: &Xform<'_>,
    surf_rot: &mut HashMap<Option<Handle<MotionNode>>, Handle<MotionNode>>,
    conj: &mut Option<Conjugation>,
) -> Result<SurfaceDef, OpError> {
    let rotate = motion.rigid().and_then(|i| i.rotate);
    let exact_rotation = motion.rigid().map(|i| i.is_exact()).unwrap_or(true);
    // Whether the motion moves anything at all. `copy` is [`transform_solid`] under a **zero**
    // translation, and the identity is not a translation: the geometry comes out bit-identical, so
    // a rotated surface's definition still describes it. Reading "no rotation" as "translated"
    // refused a boolean on any copy of a rotated solid.
    let displaced = motion
        .rigid()
        .map(|i| i.translate.iter().any(|r| r.numer() != 0))
        .unwrap_or(true);
    let source = model
        .surface_defs
        .get(&src)
        .copied()
        .unwrap_or(SurfaceDef::Inexact);
    // A node is recorded unless the whole history stays exact: a 90°-family rotation of a
    // `Constructed` surface keeps exact coefficients and needs none, but the same rotation of an
    // already-rotated one must still be recorded, or the chain would not reproduce the result.
    let parent = match source {
        SurfaceDef::Moved { motion, .. } => Some(motion),
        _ => None,
    };
    if let Some(r) = rotate.filter(|_| !exact_rotation || parent.is_some()) {
        // `entry` cannot hold a `&mut Model` across the closure, so look up then insert.
        let rotation = match surf_rot.get(&parent) {
            Some(&h) => h,
            None => {
                let h = model.push_motion(
                    Motion::Rotate {
                        axis: r.axis,
                        point: r.point,
                        angle: r.angle,
                    },
                    parent,
                );
                surf_rot.insert(parent, h);
                h
            }
        };
        return Ok(match source {
            SurfaceDef::Inexact => SurfaceDef::Inexact,
            // Keep the original witness: the new node's parent is this surface's own leaf, so
            // replaying the chain from the *root* witness applies every rotation, once each.
            SurfaceDef::Moved { witness, .. } => SurfaceDef::Moved {
                witness,
                motion: rotation,
            },
            SurfaceDef::Constructed => {
                let f = model.faces.get(face);
                let (witness, _) =
                    crate::planes::outer_tri(model, f).ok_or(OpError::DegenerateGeometry)?;
                SurfaceDef::Moved {
                    witness,
                    motion: rotation,
                }
            }
        });
    }
    Ok(match (source, conj) {
        (SurfaceDef::Constructed, _) => SurfaceDef::Constructed,
        (SurfaceDef::Inexact, _) => SurfaceDef::Inexact,
        // A reflection is exact (rational pivot, negated angle — `conjugate_chain`), so a mirrored
        // rotated surface stays describable: the mirrored witness under the conjugated chain.
        (SurfaceDef::Moved { witness, motion }, Some(c)) => {
            let m = c.mirror;
            let leaf = conjugate_chain(model, motion, c).ok_or(OpError::MirrorChainOverflow)?;
            SurfaceDef::Moved {
                witness: witness.map(|p| m.point(p)),
                motion: leaf,
            }
        }
        // A translation of a rotated surface: no node can name `R` then `T`. A *zero*
        // translation is not a translation — it is `copy`, and the definition still holds.
        (def @ SurfaceDef::Moved { .. }, None) if !displaced => def,
        (SurfaceDef::Moved { .. }, None) => SurfaceDef::Inexact,
    })
}

/// The reflected image of a rotation chain's root vertex, pushed as its own cell.
///
/// The chain now turns *this* point, so it has to exist. A `Constructed` root mirrors to a
/// `Constructed` root. A `Discovered` root keeps its definition with the three planes mirrored —
/// nothing consumes that today (`vertex_pt3` stops at a `Discovered` root either way), but
/// leaving the original planes would record the false claim that the mirrored point lies on them.
fn mirrored_base(model: &mut Model, base: Handle<Vertex>, c: &mut Conjugation) -> Handle<Vertex> {
    if let Some(&h) = c.bases.get(&base) {
        return h;
    }
    let v = *model.vertices.get(base);
    let origin = match v.origin {
        Origin::Discovered {
            tol,
            definition: VertexDef::ThreePlane(planes),
        } => {
            let mapped = planes.map(|s| {
                if let Some(&h) = c.surfaces.get(&s) {
                    return h;
                }
                // A root's planes belong to the pre-rotation solid, not to the one being walked,
                // so they are mirrored here rather than through the walk's surface map.
                let moved = model
                    .surfaces
                    .get(s)
                    .mirrored(c.mirror)
                    .unwrap_or_else(|| model.surfaces.get(s).clone());
                // A root is pre-rotation, so its planes are `Constructed` and a reflection keeps
                // them so. Anything else would be a rotated plane reached through a root, which
                // the chain invariant forbids — record it as inexact rather than assert.
                let def = match model.surface_defs.get(&s) {
                    Some(SurfaceDef::Constructed) => SurfaceDef::Constructed,
                    _ => SurfaceDef::Inexact,
                };
                let h = model.push_surface(moved, def);
                c.surfaces.insert(s, h);
                h
            });
            Origin::Discovered {
                tol,
                definition: VertexDef::ThreePlane(mapped),
            }
        }
        // A root is never `Rotated` (the chain always names the non-rotated ancestor).
        other => other,
    };
    let h = model.vertices.push(Vertex {
        point: c.mirror.point(v.point),
        origin,
    });
    c.bases.insert(base, h);
    h
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
    Mirror(AxisMirror),
}

impl Xform<'_> {
    fn point(&self, p: Point3) -> Point3 {
        match self {
            Xform::Rigid(iso) => Point3::from_array(iso.apply_point(p.as_array())),
            Xform::Mirror(m) => m.point(p),
        }
    }

    /// `None` when the variant has no image under this motion — a mirrored cylinder, whose
    /// parametrisation handedness is a curved-geometry decision (see `Surface::mirrored`).
    fn surface(&self, s: &Surface, offset: Vector3) -> Option<Surface> {
        match self {
            Xform::Rigid(iso) => Some(transform_surface(s, iso, offset)),
            Xform::Mirror(m) => s.mirrored(*m),
        }
    }

    fn curve(&self, c: &Curve, offset: Vector3) -> Option<Curve> {
        match self {
            Xform::Rigid(iso) => Some(transform_curve(c, iso, offset)),
            Xform::Mirror(m) => c.mirrored(*m),
        }
    }

    /// A reflection negates the normal a loop's winding implies (`R(a) × R(b) = −R(a × b)`), so
    /// every loop is rewound to put it back — and then the `Orientation` flag needs no change,
    /// because a reflection preserves dot products.
    fn reverses_orientation(&self) -> bool {
        matches!(self, Xform::Mirror(_))
    }

    /// The isometry, for the rotation-forest bookkeeping that only proper motion does.
    fn rigid(&self) -> Option<&Isometry> {
        match self {
            Xform::Rigid(iso) => Some(iso),
            Xform::Mirror(_) => None,
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

/// A curve moved by `isometry` (see [`transform_surface`]).
fn transform_curve(c: &Curve, iso: &Isometry, offset: Vector3) -> Curve {
    if iso.rotate.is_none() {
        return c.translated(offset);
    }
    let p = |q: Point3| Point3::from_array(iso.apply_point(q.as_array()));
    let d = |v: Vector3| Vector3::from_array(iso.apply_dir(v.as_array()));
    match c {
        Curve::Line(l) => Curve::Line(
            Line::from_point_direction(p(l.origin()), d(l.direction()))
                .expect("rotation preserves a nonzero direction"),
        ),
        Curve::Circle(ci) => Curve::Circle(
            Circle::from_center_normal(
                p(ci.center()),
                d(ci.normal()),
                d(ci.ref_dir()),
                ci.radius(),
            )
            .expect("rotation preserves a valid circle"),
        ),
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
    mirror_plane: Option<(Axis, Rat)>,
) -> Result<Handle<Solid>, OpError> {
    // A reflection carries a rotated input by conjugating its chain (see `conjugate_chain`).
    let mut conj = match (motion, mirror_plane) {
        (Xform::Mirror(m), Some((axis, offset))) => Some(Conjugation {
            mirror: *m,
            axis,
            offset,
            nodes: HashMap::new(),
            bases: HashMap::new(),
            surfaces: HashMap::new(),
        }),
        _ => None,
    };
    let offset = motion
        .rigid()
        .map(|iso| Vector3::from_array(iso.offset_f64()))
        .unwrap_or_else(Vector3::zero);
    let src = model.solids.get(solid).clone();

    // Forest node for this transform (§CIP ⑦), one shared node named by every rotated
    // vertex. The input's shared leaf (None if the input is not rotated) decides B0 vs B1:
    //   A.  translation → no node (remap path).
    //   B0. fresh rotation (input not rotated): exact (90°-family) → no node (remap,
    //       preserving 1b); inexact → a root node (`parent = None`).
    //   B1. re-rotation (input already rotated): **always chain** a node (`parent =
    //       input leaf`) — record every rotation, even an exact one, because an inexact
    //       ancestor makes the composite inexact and stage-2 tol must transport through
    //       it; the forest stays complete. (Same-axis *bundling* — accumulating the
    //       angle into one node — is a later cell; this cell always chains.)
    let input_leaf = solid_rotation(model, solid);
    let rot_node: Option<Handle<MotionNode>> =
        match (motion.rigid().and_then(|i| i.rotate), input_leaf) {
            (None, _) => None,
            (Some(r), None) => {
                (!motion.rigid().expect("rotate implies rigid").is_exact()).then(|| {
                    model.push_motion(
                        Motion::Rotate {
                            axis: r.axis,
                            point: r.point,
                            angle: r.angle,
                        },
                        None,
                    )
                })
            }
            (Some(r), Some(parent)) => Some(model.push_motion(
                Motion::Rotate {
                    axis: r.axis,
                    point: r.point,
                    angle: r.angle,
                },
                Some(parent),
            )),
        };

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
    let mut surf_rot: HashMap<Option<Handle<MotionNode>>, Handle<MotionNode>> = HashMap::new();
    let mut surf_map: HashMap<Handle<Surface>, Handle<Surface>> = HashMap::new();
    for &fh in &face_order {
        let s = model.faces.get(fh).surface;
        if surf_map.contains_key(&s) {
            continue;
        }
        let moved = motion
            .surface(model.surfaces.get(s), offset)
            .ok_or(OpError::MirrorNotPlanar)?;
        let def = moved_surface_def(model, s, fh, motion, &mut surf_rot, &mut conj)?;
        surf_map.insert(s, model.push_surface(moved, def));
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

    // Pass 2 — curves (dedup, moved).
    let mut curve_map: HashMap<Handle<Curve>, Handle<Curve>> = HashMap::new();
    for &eh in &edge_order {
        let c = model.edges.get(eh).curve;
        if let std::collections::hash_map::Entry::Vacant(e) = curve_map.entry(c) {
            let moved = motion
                .curve(model.curves.get(c), offset)
                .ok_or(OpError::MirrorNotPlanar)?;
            e.insert(model.curves.push(moved));
        }
    }

    // Pass 3 — vertices (dedup, moved point + Origin preserved/remapped).
    let mut vert_order: Vec<Handle<Vertex>> = Vec::new();
    let mut vert_seen: HashSet<Handle<Vertex>> = HashSet::new();
    for &eh in &edge_order {
        if let Some(bounds) = model.edges.get(eh).bounds {
            for vh in bounds {
                if vert_seen.insert(vh) {
                    vert_order.push(vh);
                }
            }
        }
    }
    let mut vert_map: HashMap<Handle<Vertex>, Handle<Vertex>> = HashMap::new();
    for &vh in &vert_order {
        let v = *model.vertices.get(vh);
        let new_v = Vertex {
            point: motion.point(v.point),
            // A recorded rotation marks the vertex `Rotated`; `base` is the **root** —
            // the non-`Rotated` (Constructed/Discovered) ancestor whose exact definition
            // the rotation chain turns. A fresh rotation's input is itself the root; a
            // re-rotation chases one hop to the input's own root (the invariant keeps
            // `base` pointing at a root, never at another `Rotated` vertex, so stage-2
            // recompute never applies the same rotation twice). No node → remap (1a/1b).
            origin: match rot_node {
                Some(rotation) => {
                    let base = match v.origin {
                        Origin::Moved { base, .. } => base,
                        _ => vh,
                    };
                    Origin::Moved {
                        base,
                        motion: rotation,
                    }
                }
                None => remap_origin(v.origin, &surf_map),
            },
        };
        // A reflection of an already-rotated vertex: the definition is the conjugated chain over
        // the mirrored root, and the *coordinate must be replayed from it* — `vertex_pt3` requires
        // the replay to match the stored point bit for bit, which reflecting the point separately
        // would not (two float routes to one real number differ in the last places).
        let new_v = match (&mut conj, v.origin) {
            (
                Some(c),
                Origin::Moved {
                    base,
                    motion: rotation,
                },
            ) => {
                let leaf =
                    conjugate_chain(model, rotation, c).ok_or(OpError::MirrorChainOverflow)?;
                let mbase = mirrored_base(model, base, c);
                let bp = model.vertices.get(mbase).point.as_array();
                let point = match crate::rotated_vertex::replay_chain_coord(model, bp, leaf) {
                    Ok(coord) => Point3::from_array(coord),
                    // A `Discovered` root is not replayable — before or after mirroring — so the
                    // reflected coordinate is the honest value and `vertex_pt3` defers as it did.
                    Err(_) => new_v.point,
                };
                Vertex {
                    point,
                    origin: Origin::Moved {
                        base: mbase,
                        motion: leaf,
                    },
                }
            }
            _ => new_v,
        };
        vert_map.insert(vh, model.vertices.push(new_v));
    }

    // Pass 4 — edges (curve/vertex handles + Origin remapped).
    let mut edge_map: HashMap<Handle<Edge>, Handle<Edge>> = HashMap::new();
    for &eh in &edge_order {
        let e = *model.edges.get(eh);
        let new_e = Edge {
            curve: curve_map[&e.curve],
            bounds: e.bounds.map(|[a, b]| [vert_map[&a], vert_map[&b]]),
            origin: remap_origin(e.origin, &surf_map),
        };
        edge_map.insert(eh, model.edges.push(new_e));
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
            orientation: face.orientation,
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
fn solid_rotation(model: &Model, solid: Handle<Solid>) -> Option<Handle<MotionNode>> {
    let sh = model.solids.get(solid).outer;
    let mut seen: Option<Option<Handle<MotionNode>>> = None;
    for &fh in &model.shells.get(sh).faces {
        for he in &model.faces.get(fh).outer.half_edges {
            let Some(bounds) = model.edges.get(he.edge).bounds else {
                continue;
            };
            for &vh in &bounds {
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
