//! STEP (AP242 Edition 2) export for the nacre kernel.
//!
//! [`to_step`] is a thin **adapter**: it translates a kernel [`Model`] into
//! AP242 entities and delegates serialization to `brep-to-step`, a shape-only
//! STEP writer. The rest of the kernel is backend-agnostic — the writer is
//! confined to this crate, the same isolation used for the OCCT helper protocol.
//!
//! Coverage: planar + cylindrical b-rep — `Surface::{Plane, Cylinder}` bounded by
//! `Curve::{Line, Circle}` (a whole cylinder: two seam vertices, two full-circle rims, and a
//! lateral bounded by those two rims — one its outer bound, the other an inner one). The
//! surface/curve `match`es stay exhaustive, so future variants (sphere, NURBS)
//! force a compile error here. AP242 Ed2 is stamped by the writer.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use brep_to_step::{
    Bound, Curve as StepCurve, Edge as StepEdge, Face, Frame, Header, StepWriter,
    Surface as StepSurface, Units, Vertex as StepVertex, VoidShellNormals,
};
use nacre_geom::{Curve, Surface};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Edge, Model, Orientation, Shell, Solid, Vertex};
use std::collections::HashMap;

/// A failure while translating a [`Model`] to STEP.
#[derive(Debug)]
pub enum StepError {
    /// The writer refused an input (its `brep_to_step::Error`, as its message,
    /// so the writer's type does not leak into the public API).
    Backend(String),
}

impl From<brep_to_step::Error> for StepError {
    fn from(e: brep_to_step::Error) -> Self {
        StepError::Backend(e.to_string())
    }
}

/// Export every live solid in `model` to AP242 (Ed2) STEP text.
///
/// Planar and cylindrical faces (`Surface::{Plane, Cylinder}`) bounded by lines,
/// arcs and full circles (`Curve::{Line, Circle}`); a full-circle rim is a closed edge
/// (`vertices: [v, v]`, start == end, `v` its seam vertex) that emits a closed STEP circle. A
/// solid with cavity shells is exported as a `BREP_WITH_VOIDS`. Coordinates are
/// emitted in millimetres (nacre is unitless; STEP needs a unit).
///
/// Every value is written from the model's caches as they stand. An operation realizes them as it
/// builds — each the nearest `f64` of the exact geometry — except behind a history longer than the
/// caches pay for at build time (a chain of more than 192 recorded motions, or a coordinate two
/// rungs do not decide), where the construction's own figure stands. `nacre_ops::refine_caches`
/// (`nacre::ops::refine_caches`) pays for those and reports what it could not settle; call it
/// before exporting such a model. It takes `&mut Model`, which an export does not.
///
/// ★ **`timestamp` is the header's `FILE_NAME` time stamp, written verbatim** (an ISO 8601 string
/// such as `2026-10-05T12:00:00Z`; `""` leaves the field blank). The caller gives it because the
/// kernel reads no clock: the same model and the same time stamp give the same bytes, and
/// `wasm32-unknown-unknown` has no clock to read — asking the system for the time there panics.
/// Part 21 caps a header string at 256 characters — a longer stamp is [`StepError::Backend`] —
/// and writes a character outside ASCII as an `\X2\` escape, so verbatim means verbatim ASCII.
pub fn to_step(model: &Model, timestamp: &str) -> Result<String, StepError> {
    build_step(model, model.live_solids(), timestamp)
}

/// Export the given solids to AP242 (Ed2) STEP text, one STEP part each, in the order given.
/// Same coverage, time stamp and errors as [`to_step`].
///
/// For a caller that holds more live solids than it means to write — an application exporting
/// what it shows, while copies and unconsumed intermediates stay live in the model — and for the
/// boolean oracle, whose OCCT `fuse`/`cut`/`common` take single-solid files (`&[solid]`). Each
/// solid must be live; a superseded one would be written as if it were still there.
pub fn to_step_solids(
    model: &Model,
    solids: &[Handle<Solid>],
    timestamp: &str,
) -> Result<String, StepError> {
    build_step(model, solids, timestamp)
}

/// The shared body of [`to_step`] and [`to_step_solids`].
fn build_step(
    model: &Model,
    solids: &[Handle<Solid>],
    timestamp: &str,
) -> Result<String, StepError> {
    let mut b = StepWriter::new(Units::default())?;

    // Export the given (live) solids, not the whole append-only store:
    // superseded solids linger in the arena but must not reach the file.
    for &solid_h in solids {
        let solid = model.solid(solid_h);
        let part = b.part("nacre_solid");

        // Deduped per export: a shared vertex/edge becomes one STEP entity.
        let mut vmap: HashMap<Handle<Vertex>, StepVertex> = HashMap::new();
        let mut emap: HashMap<Handle<Edge>, StepEdge> = HashMap::new();

        let outer = build_shell_faces(&mut b, model, solid.outer, &mut vmap, &mut emap)?;
        if solid.cavities.is_empty() {
            b.solid(part, &outer)?;
        } else {
            // A hollow solid → BREP_WITH_VOIDS. Each cavity shell is built the
            // same way as the outer shell: nacre stores a cavity with its face
            // normals pointing into the void (away from the material, like the
            // outer shell points outward), which is exactly `AwayFromMaterial`
            // — the writer keeps the authored orientation (no reversal).
            let mut voids = Vec::with_capacity(solid.cavities.len());
            for &cavity in &solid.cavities {
                voids.push(build_shell_faces(
                    &mut b, model, cavity, &mut vmap, &mut emap,
                )?);
            }
            b.solid_with_voids(part, &outer, &voids, VoidShellNormals::AwayFromMaterial)?;
        }
    }

    b.finish(&Header {
        timestamp: timestamp.to_owned(),
        originating_system: "nacre".to_owned(),
        ..Header::default()
    })
    .map_err(StepError::from)
}

/// Build every face of one shell as STEP faces, deduping shared
/// vertices/edges through `vmap`/`emap`. Used for both the outer shell and each
/// cavity shell of a solid (a cavity is emitted identically — its faces already
/// carry the correct inward orientation).
fn build_shell_faces(
    b: &mut StepWriter,
    model: &Model,
    shell: Handle<Shell>,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
    emap: &mut HashMap<Handle<Edge>, StepEdge>,
) -> Result<Vec<Face>, StepError> {
    let mut face_ids = Vec::new();
    for &fh in &model.shell(shell).faces {
        let face = model.face(fh);

        // Surface → the writer's surface. Exhaustive match: future variants
        // (sphere, NURBS) must be handled here (else compile error).
        let surface = match model.surface_cache(face.surface) {
            // ref_dir is cosmetic for a bounded planar face — any perpendicular.
            Surface::Plane(p) => StepSurface::Plane(frame(
                p.origin(),
                p.normal(),
                p.normal()
                    .any_perpendicular()
                    .expect("unit normal has a perpendicular"),
            )),
            Surface::Cylinder(c) => StepSurface::Cylinder {
                frame: frame(c.axis().origin(), c.axis().direction(), c.ref_dir()),
                radius: c.radius(),
            },
        };
        let same_sense = matches!(face.orientation, Orientation::Forward);

        let mut bounds = Vec::with_capacity(1 + face.inner.len());
        bounds.push(Bound::outer(build_loop(b, model, &face.outer, vmap, emap)?));
        for inner in &face.inner {
            bounds.push(Bound::inner(build_loop(b, model, inner, vmap, emap)?));
        }

        face_ids.push(b.face(surface, same_sense, &bounds)?);
    }
    Ok(face_ids)
}

/// Build a loop's edges as `(edge, forward)` for a [`Bound`].
fn build_loop(
    b: &mut StepWriter,
    model: &Model,
    lp: &nacre_topo::Loop,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
    emap: &mut HashMap<Handle<Edge>, StepEdge>,
) -> Result<Vec<(StepEdge, bool)>, StepError> {
    let mut edges = Vec::with_capacity(lp.half_edges.len());
    for he in &lp.half_edges {
        edges.push((build_edge(b, model, he.edge, vmap, emap)?, he.forward));
    }
    Ok(edges)
}

fn build_edge(
    b: &mut StepWriter,
    model: &Model,
    eh: Handle<Edge>,
    vmap: &mut HashMap<Handle<Vertex>, StepVertex>,
    emap: &mut HashMap<Handle<Edge>, StepEdge>,
) -> Result<StepEdge, StepError> {
    if let Some(&id) = emap.get(&eh) {
        return Ok(id);
    }
    let edge = model.edge(eh);
    let [v0, v1] = edge.vertices;
    // A line carries its cached direction — the truth's where the truth names one
    // (`Model::derive_edge_curve`) — written bit for bit; a circle carries its own frame.
    // A whole rim has v0 == v1, giving a closed STEP circle.
    let curve = match model.edge_curve(eh) {
        Curve::Line(l) => StepCurve::LineAlong(l.direction().as_array()),
        Curve::Circle(c) => StepCurve::Circle {
            frame: frame(c.center(), c.normal(), c.ref_dir()),
            radius: c.radius(),
        },
    };
    let sv0 = build_vertex(b, model, v0, vmap)?;
    let sv1 = build_vertex(b, model, v1, vmap)?;
    let id = b.edge(sv0, sv1, curve)?;
    emap.insert(eh, id);
    Ok(id)
}

fn build_vertex(
    b: &mut StepWriter,
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

/// The writer's [`Frame`] (an `AXIS2_PLACEMENT_3D`) from kernel vectors.
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
