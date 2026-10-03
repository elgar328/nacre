//! STEP (AP242 Edition 2) export for the nacre kernel.
//!
//! [`to_step`] is a thin **adapter**: it translates a kernel [`Model`] into
//! AP242 entities and delegates serialization to a swappable backend
//! (`step-io` during development). The rest of the kernel is backend-agnostic —
//! `step-io` is confined to this crate, the same isolation used for the OCCT
//! helper protocol.
//!
//! Coverage: planar + cylindrical b-rep — `Surface::{Plane, Cylinder}` bounded by
//! `Curve::{Line, Circle}` (a cylinder is the seam model: two seam vertices, two full-circle
//! rims, one seam line the lateral uses twice). The
//! surface/curve `match`es stay exhaustive, so future variants (sphere, NURBS)
//! force a compile error here. AP242 Ed2 is stamped by the backend.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use nacre_geom::{Curve, Surface};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Edge, Model, Orientation, Shell, Solid, Vertex};
use std::collections::HashMap;
use step_io::StepBuilder;
use step_io::build::Vertex as StepVertex;
use step_io::build::{
    CurveInput, FaceBoundInput, Frame, HeaderInput, SurfaceInput, VoidShellNormals,
};
use step_io::generated::model::{AdvancedFaceId, EdgeCurveId};

/// A failure while translating a [`Model`] to STEP.
#[derive(Debug)]
pub enum StepError {
    /// The `step-io` backend rejected the entity graph (its `AuthorError`,
    /// stringified so the backend type does not leak into the public API).
    Backend(String),
}

impl From<step_io::AuthorError> for StepError {
    fn from(e: step_io::AuthorError) -> Self {
        StepError::Backend(format!("{e:?}"))
    }
}

/// Export every live solid in `model` to AP242 (Ed2) STEP text.
///
/// Planar and cylindrical faces (`Surface::{Plane, Cylinder}`) bounded by lines
/// and full circles (`Curve::{Line, Circle}`); a full-circle rim is a seam edge
/// (`vertices: [v, v]`, start == end) that emits a closed STEP circle. A
/// solid with cavity shells is exported as a `BREP_WITH_VOIDS`. Coordinates are
/// emitted in millimetres (nacre is unitless; STEP needs a unit).
///
/// Every value is written from the model's caches as they stand. An operation realizes them as it
/// builds — each the nearest `f64` of the exact geometry — except behind a history longer than the
/// caches pay for at build time (a chain of more than 192 recorded motions, or a coordinate two
/// rungs do not decide), where the construction's own figure stands. `nacre_ops::refine_caches`
/// (`nacre::ops::refine_caches`) pays for those and reports what it could not settle; call it
/// before exporting such a model. It takes `&mut Model`, which an export does not.
pub fn to_step(model: &Model) -> Result<String, StepError> {
    build_step(model, model.live_solids())
}

/// Export a single solid to AP242 (Ed2) STEP text — the solid's live closure
/// only, as one STEP part. Same coverage and errors as [`to_step`].
///
/// Needed by the M5 boolean oracle: OCCT's binary `fuse`/`cut`/`common` take two
/// separate single-solid STEP files, whereas [`to_step`] emits the whole live set
/// as one file.
pub fn to_step_solid(model: &Model, solid: Handle<Solid>) -> Result<String, StepError> {
    build_step(model, &[solid])
}

/// Export the given solids to AP242 (Ed2) STEP text, one STEP part each. The
/// shared body of [`to_step`] and [`to_step_solid`].
fn build_step(model: &Model, solids: &[Handle<Solid>]) -> Result<String, StepError> {
    let mut b = StepBuilder::new()?;
    b.header(&HeaderInput {
        originating_system: Some("nacre".to_owned()),
        ..Default::default()
    });

    // Export the given (live) solids, not the whole append-only store:
    // superseded solids linger in the arena but must not reach the file.
    for &solid_h in solids {
        let solid = model.solid(solid_h);
        let part = b.part("nacre_solid")?;

        // Deduped per export: a shared vertex/edge becomes one STEP entity.
        let mut vmap: HashMap<Handle<Vertex>, StepVertex> = HashMap::new();
        let mut emap: HashMap<Handle<Edge>, EdgeCurveId> = HashMap::new();

        let outer = build_shell_faces(&mut b, model, solid.outer, &mut vmap, &mut emap)?;
        if solid.cavities.is_empty() {
            b.solid(part, "body", outer)?;
        } else {
            // A hollow solid → BREP_WITH_VOIDS. Each cavity shell is built the
            // same way as the outer shell: nacre stores a cavity with its face
            // normals pointing into the void (away from the material, like the
            // outer shell points outward), which is exactly `AwayFromMaterial`
            // — step-io keeps the authored orientation (no reversal).
            let mut voids = Vec::with_capacity(solid.cavities.len());
            for &cavity in &solid.cavities {
                voids.push(build_shell_faces(
                    &mut b, model, cavity, &mut vmap, &mut emap,
                )?);
            }
            b.solid_with_voids(
                part,
                "body",
                outer,
                voids,
                VoidShellNormals::AwayFromMaterial,
            )?;
        }
    }

    b.finish().map_err(StepError::from)
}

/// Build every face of one shell as step-io `AdvancedFaceId`s, deduping shared
/// vertices/edges through `vmap`/`emap`. Used for both the outer shell and each
/// cavity shell of a solid (a cavity is emitted identically — its faces already
/// carry the correct inward orientation).
fn build_shell_faces(
    b: &mut StepBuilder,
    model: &Model,
    shell: Handle<Shell>,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
    emap: &mut HashMap<Handle<Edge>, EdgeCurveId>,
) -> Result<Vec<AdvancedFaceId>, StepError> {
    let mut face_ids = Vec::new();
    for &fh in &model.shell(shell).faces {
        let face = model.face(fh);

        // Surface → step-io SurfaceInput. Exhaustive match: future variants
        // (sphere, NURBS) must be handled here (else compile error).
        let surface = match model.surface_cache(face.surface) {
            // ref_dir is cosmetic for a bounded planar face — any perpendicular.
            Surface::Plane(p) => SurfaceInput::Plane(frame(
                p.origin(),
                p.normal(),
                p.normal()
                    .any_perpendicular()
                    .expect("unit normal has a perpendicular"),
            )),
            Surface::Cylinder(c) => SurfaceInput::Cylinder(
                frame(c.axis().origin(), c.axis().direction(), c.ref_dir()),
                c.radius(),
            ),
        };
        let same_sense = matches!(face.orientation, Orientation::Forward);

        let mut bounds = Vec::with_capacity(1 + face.inner.len());
        bounds.push(FaceBoundInput::outer(build_loop(
            b,
            model,
            &face.outer,
            vmap,
            emap,
        )?));
        for inner in &face.inner {
            bounds.push(FaceBoundInput::inner(build_loop(
                b, model, inner, vmap, emap,
            )?));
        }

        face_ids.push(b.face(surface, same_sense, bounds)?);
    }
    Ok(face_ids)
}

/// Build a loop's edges as `(edge, forward)` for a `FaceBoundInput`.
fn build_loop(
    b: &mut StepBuilder,
    model: &Model,
    lp: &nacre_topo::Loop,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
    emap: &mut HashMap<Handle<Edge>, EdgeCurveId>,
) -> Result<Vec<(EdgeCurveId, bool)>, StepError> {
    let mut edges = Vec::with_capacity(lp.half_edges.len());
    for he in &lp.half_edges {
        edges.push((build_edge(b, model, he.edge, vmap, emap)?, he.forward));
    }
    Ok(edges)
}

fn build_edge(
    b: &mut StepBuilder,
    model: &Model,
    eh: Handle<Edge>,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
    emap: &mut HashMap<Handle<Edge>, EdgeCurveId>,
) -> Result<EdgeCurveId, StepError> {
    if let Some(&id) = emap.get(&eh) {
        return Ok(id);
    }
    let edge = model.edge(eh);
    let [v0, v1] = edge.vertices;
    // A line carries its cached direction — the truth's where the truth names one
    // (`Model::derive_edge_curve`) — written bit for bit; a circle carries its own frame.
    // A seam rim has v0 == v1, giving a closed STEP circle.
    let curve = match model.edge_curve(eh) {
        Curve::Line(l) => CurveInput::LineAlong(l.direction().as_array()),
        Curve::Circle(c) => {
            CurveInput::Circle(frame(c.center(), c.normal(), c.ref_dir()), c.radius())
        }
    };
    let sv0 = build_vertex(b, model, v0, vmap)?;
    let sv1 = build_vertex(b, model, v1, vmap)?;
    let id = b.edge(sv0, sv1, curve)?;
    emap.insert(eh, id);
    Ok(id)
}

fn build_vertex(
    b: &mut StepBuilder,
    model: &Model,
    vh: Handle<Vertex>,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
) -> Result<StepVertex, StepError> {
    if let Some(&v) = vmap.get(&vh) {
        return Ok(v);
    }
    let v = b.vertex(model.vertex_point(vh).as_array())?;
    vmap.insert(vh, v);
    Ok(v)
}

/// A step-io `Frame` (an `AXIS2_PLACEMENT_3D`) from kernel vectors.
fn frame(origin: Point3, axis: Vector3, ref_dir: Vector3) -> Frame {
    Frame {
        origin: origin.as_array(),
        axis: axis.as_array(),
        ref_dir: ref_dir.as_array(),
    }
}

#[cfg(test)]
#[path = "tests/lib.rs"]
mod tests;
