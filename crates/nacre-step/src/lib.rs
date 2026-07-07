//! STEP (AP242 Edition 2) export for the nacre kernel.
//!
//! [`to_step`] is a thin **adapter**: it translates a kernel [`Model`] into
//! AP242 entities and delegates serialization to a swappable backend
//! (`step-io` during development). The rest of the kernel is backend-agnostic —
//! `step-io` is confined to this crate, the same isolation used for
//! `BooleanEngine` and the OCCT helper protocol.
//!
//! M2 coverage is **planar b-rep only** (`Surface::Plane` + `Curve::Line`); the
//! surface/curve `match`es are exhaustive today, so when M3 adds variants the
//! compiler forces them to be handled here. AP242 Ed2 is stamped by the backend.

use nacre_geom::{Curve, Surface};
use nacre_store::Handle;
use nacre_topo::{Edge, Model, Orientation, Vertex};
use std::collections::HashMap;
use step_io::StepBuilder;
use step_io::build::Vertex as StepVertex;
use step_io::build::{CurveInput, FaceBoundInput, Frame, HeaderInput, SurfaceInput};
use step_io::generated::model::EdgeCurveId;

/// A failure while translating a [`Model`] to STEP.
#[derive(Debug)]
pub enum StepError {
    /// An edge has no endpoint vertices (a closed edge — M3). M2 expects all
    /// edges bounded.
    UnboundedEdge,
    /// A solid has inner cavity shells (`brep_with_voids` — deferred). M2
    /// exports the outer shell only.
    Cavities,
    /// The `step-io` backend rejected the entity graph (its `AuthorError`,
    /// stringified so the backend type does not leak into the public API).
    Backend(String),
    /// Curved geometry (`Curve::Circle` / `Surface::Cylinder`) — STEP emission
    /// is deferred to the next unit. Transitional: this variant and its two
    /// match arms are replaced by real emission (circle/cylinder) then.
    UnsupportedCurvedGeometry,
}

impl From<step_io::AuthorError> for StepError {
    fn from(e: step_io::AuthorError) -> Self {
        StepError::Backend(format!("{e:?}"))
    }
}

/// Export every solid in `model` to AP242 (Ed2) STEP text.
///
/// M2: planar faces (`Surface::Plane`) bounded by straight edges
/// (`Curve::Line`). Rejects closed edges ([`StepError::UnboundedEdge`]) and
/// cavities ([`StepError::Cavities`]). Coordinates are emitted in millimetres
/// (nacre is unitless; STEP needs a unit).
pub fn to_step(model: &Model) -> Result<String, StepError> {
    let mut b = StepBuilder::new()?;
    b.header(&HeaderInput {
        originating_system: Some("nacre".to_owned()),
        ..Default::default()
    });

    for (_, solid) in model.solids.iter() {
        if !solid.cavities.is_empty() {
            return Err(StepError::Cavities);
        }
        let part = b.part("nacre_solid")?;
        let shell = model.shells.get(solid.outer);

        // Deduped per export: a shared vertex/edge becomes one STEP entity.
        let mut vmap: HashMap<Handle<Vertex>, StepVertex> = HashMap::new();
        let mut emap: HashMap<Handle<Edge>, EdgeCurveId> = HashMap::new();
        let mut face_ids = Vec::new();

        for &fh in &shell.faces {
            let face = model.faces.get(fh);

            // Surface → plane frame. Exhaustive match: M3's new variants must be
            // handled here (else compile error), which is the honest rejection point.
            let frame = match model.surfaces.get(face.surface) {
                Surface::Plane(p) => Frame {
                    origin: p.origin().as_array(),
                    axis: p.normal().as_array(),
                    // Any perpendicular works — ref_dir is cosmetic for a bounded planar face.
                    ref_dir: p
                        .normal()
                        .any_perpendicular()
                        .expect("unit normal has a perpendicular")
                        .as_array(),
                },
                // Cylinder STEP emission (CYLINDRICAL_SURFACE) lands in the next unit.
                Surface::Cylinder(_) => return Err(StepError::UnsupportedCurvedGeometry),
            };
            let same_sense = matches!(face.orientation, Orientation::Forward);

            let mut bounds = Vec::with_capacity(1 + face.inner.len());
            bounds.push(FaceBoundInput::outer(build_loop(
                &mut b,
                model,
                &face.outer,
                &mut vmap,
                &mut emap,
            )?));
            for inner in &face.inner {
                bounds.push(FaceBoundInput::inner(build_loop(
                    &mut b, model, inner, &mut vmap, &mut emap,
                )?));
            }

            face_ids.push(b.face(SurfaceInput::Plane(frame), same_sense, bounds)?);
        }

        b.solid(part, "body", face_ids)?;
    }

    b.finish().map_err(StepError::from)
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
    let edge = model.edges.get(eh);
    let [v0, v1] = edge.bounds.ok_or(StepError::UnboundedEdge)?;
    // Exhaustive: CurveInput::Line derives geometry from the two vertices.
    let curve = match model.curves.get(edge.curve) {
        Curve::Line(_) => CurveInput::Line,
        // Circle STEP emission (CurveInput::Circle) lands in the next unit.
        Curve::Circle(_) => return Err(StepError::UnsupportedCurvedGeometry),
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
    let v = b.vertex(model.vertices.get(vh).point.as_array())?;
    vmap.insert(vh, v);
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::{Point3, Vector3};
    use nacre_topo::Solid;
    use step_io::read;
    use step_io::scene::geometry::{CurveKind, SurfaceKind};

    fn cuboid(min: [f64; 3], max: [f64; 3]) -> Model {
        let mut m = Model::new();
        m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
        m
    }

    #[test]
    fn cube_round_trips_through_step_io_reader() {
        let text = to_step(&cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).expect("export");
        let (model, report) = read(text.as_bytes()).expect("re-read");

        // No dropped/orphan entities — the structural lint (step-loupe's report).
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
        // A solid-bearing part promotes to ABSR, not a plain SR.
        assert_eq!(
            model.advanced_brep_shape_representation_arena.items.len(),
            1
        );
        assert_eq!(model.shape_representation_arena.items.len(), 0);

        let scene = model.scene();
        let solids: Vec<_> = scene.all_solids().collect();
        assert_eq!(solids.len(), 1);
        let faces: Vec<_> = solids[0].faces().collect();
        assert_eq!(faces.len(), 6);
        for face in &faces {
            assert!(matches!(face.surface().kind(), SurfaceKind::Plane(_)));
            let bounds: Vec<_> = face.bounds().collect();
            assert_eq!(bounds.len(), 1);
            let edges: Vec<_> = bounds[0].oriented_edges().collect();
            assert_eq!(edges.len(), 4);
            for (edge, _forward) in edges {
                assert!(matches!(edge.curve().kind(), CurveKind::Line(_)));
            }
        }
    }

    #[test]
    fn output_has_expected_entities_and_ap242e2_schema() {
        let text = to_step(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])).expect("export");
        for needle in [
            "MANIFOLD_SOLID_BREP",
            "CLOSED_SHELL",
            "ADVANCED_FACE",
            "PLANE",
            "LINE",
            "CARTESIAN_POINT",
        ] {
            assert!(text.contains(needle), "missing {needle}");
        }
        // AP242 Edition 2 schema token.
        assert!(text.contains("442 3 1 4"), "not AP242 Ed2");
    }

    #[test]
    fn asymmetric_box_round_trips() {
        let text = to_step(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])).expect("export");
        let (_, report) = read(text.as_bytes()).expect("re-read");
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
    }

    #[test]
    fn solid_with_cavities_is_rejected() {
        let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let shell = m.shells.iter().next().unwrap().0; // an existing shell handle
        m.solids.push(Solid {
            outer: shell,
            cavities: vec![shell],
        });
        assert!(matches!(to_step(&m), Err(StepError::Cavities)));
    }

    #[test]
    fn cylinder_curved_geometry_is_rejected_for_now() {
        // Curved STEP emission (CYLINDRICAL_SURFACE / CIRCLE) lands next unit; the
        // adapter rejects it honestly rather than emitting a plane/line.
        let mut m = Model::new();
        m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.0,
            2.0,
        );
        assert!(matches!(
            to_step(&m),
            Err(StepError::UnsupportedCurvedGeometry)
        ));
    }
}
