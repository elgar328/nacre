//! Rigid-body transform of a solid (design overhaul stage 1): copy a solid under an isometry,
//! remapping every surface/curve and carrying each vertex's exact `Origin` forward so a rotated
//! operand stays exactly defined. [`copy`] is the same walk with no motion at all.

use crate::OpError;
use nacre_geom::{Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Origin, Rotation, Shell, Solid, Vertex, VertexDef,
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
    let out = transform_solid(model, solid, isometry);
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
    Ok(transform_solid(model, solid, &zero))
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
        Origin::Rotated { .. } => origin,
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

/// Clone `solid` into a new solid with every cell's geometry moved by `isometry`,
/// preserving topology, shared cells (surfaces/curves/vertices/edges are deduped),
/// face orientations (a translation does not rotate normals), inner-loop holes,
/// cavity shells, and each vertex/edge `Origin` (a `Discovered` definition's plane
/// handles are remapped to the moved surfaces). Cells are pushed in a **deterministic
/// traversal order** (shell → face → loop) with per-cell dedup maps, so the same
/// op-log reproduces identical handles (replay determinism, DNA 3). Stage 1b's
/// rotation reuses this by swapping the per-cell geometry transform.
fn transform_solid(model: &mut Model, solid: Handle<Solid>, isometry: &Isometry) -> Handle<Solid> {
    let offset = Vector3::from_array(isometry.offset_f64());
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
    let rot_node: Option<Handle<Rotation>> = match (isometry.rotate, input_leaf) {
        (None, _) => None,
        (Some(r), None) => (!isometry.is_exact()).then(|| {
            model.rotations.push(Rotation {
                axis: r.axis,
                point: r.point,
                angle: r.angle,
                parent: None,
            })
        }),
        (Some(r), Some(parent)) => Some(model.rotations.push(Rotation {
            axis: r.axis,
            point: r.point,
            angle: r.angle,
            parent: Some(parent),
        })),
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
    let mut surf_map: HashMap<Handle<Surface>, Handle<Surface>> = HashMap::new();
    for &fh in &face_order {
        let s = model.faces.get(fh).surface;
        if let std::collections::hash_map::Entry::Vacant(e) = surf_map.entry(s) {
            let moved = transform_surface(model.surfaces.get(s), isometry, offset);
            e.insert(model.surfaces.push(moved));
        }
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
            let moved = transform_curve(model.curves.get(c), isometry, offset);
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
            point: Point3::from_array(isometry.apply_point(v.point.as_array())),
            // A recorded rotation marks the vertex `Rotated`; `base` is the **root** —
            // the non-`Rotated` (Constructed/Discovered) ancestor whose exact definition
            // the rotation chain turns. A fresh rotation's input is itself the root; a
            // re-rotation chases one hop to the input's own root (the invariant keeps
            // `base` pointing at a root, never at another `Rotated` vertex, so stage-2
            // recompute never applies the same rotation twice). No node → remap (1a/1b).
            origin: match rot_node {
                Some(rotation) => {
                    let base = match v.origin {
                        Origin::Rotated { base, .. } => base,
                        _ => vh,
                    };
                    Origin::Rotated { base, rotation }
                }
                None => remap_origin(v.origin, &surf_map),
            },
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
    let map_loop = |lp: &Loop| Loop {
        half_edges: lp
            .half_edges
            .iter()
            .map(|he| HalfEdge {
                edge: edge_map[&he.edge],
                forward: he.forward,
            })
            .collect(),
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
    model.push_solid(new_solid)
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
fn solid_rotation(model: &Model, solid: Handle<Solid>) -> Option<Handle<Rotation>> {
    let sh = model.solids.get(solid).outer;
    let mut seen: Option<Option<Handle<Rotation>>> = None;
    for &fh in &model.shells.get(sh).faces {
        for he in &model.faces.get(fh).outer.half_edges {
            let Some(bounds) = model.edges.get(he.edge).bounds else {
                continue;
            };
            for &vh in &bounds {
                let leaf = match model.vertices.get(vh).origin {
                    Origin::Rotated { rotation, .. } => Some(rotation),
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
