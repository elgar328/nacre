//! STEP (AP242 Edition 2) export for the nacre kernel.
//!
//! [`to_step`] is a thin **adapter**: it translates a kernel [`Model`] into
//! AP242 entities and delegates serialization to a swappable backend
//! (`step-io` during development). The rest of the kernel is backend-agnostic —
//! `step-io` is confined to this crate, the same isolation used for
//! `BooleanEngine` and the OCCT helper protocol.
//!
//! Coverage: planar + cylindrical b-rep — `Surface::{Plane, Cylinder}` bounded by
//! `Curve::{Line, Circle}` (a cylinder is the seam model of `add_cylinder`). The
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
pub fn to_step(model: &Model) -> Result<String, StepError> {
    build_step(model, &model.live_solids)
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

    // Export the given (live) solids, not the whole append-only store (design
    // §2): superseded solids linger in the arena but must not reach the file.
    for &solid_h in solids {
        let solid = model.solids.get(solid_h);
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
    for &fh in &model.shells.get(shell).faces {
        let face = model.faces.get(fh);

        // Surface → step-io SurfaceInput. Exhaustive match: future variants
        // (sphere, NURBS) must be handled here (else compile error).
        let surface = match model.surface(face.surface) {
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
    let edge = model.edges.get(eh);
    let [v0, v1] = edge.vertices;
    // CurveInput::Line derives geometry from the two vertices; a circle carries
    // its own frame. A seam rim has v0 == v1, giving a closed STEP circle.
    let curve = match model.edge_curve(eh) {
        Curve::Line(_) => CurveInput::Line,
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
    fn single_solid_export_isolates_one_solid() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([2.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );

        // The whole live model exports both solids.
        let both = to_step(&m).expect("export both");
        let (model_both, _) = read(both.as_bytes()).expect("re-read");
        assert_eq!(model_both.scene().all_solids().count(), 2);

        // A single-solid export isolates exactly one solid with its six faces.
        for h in [a, b] {
            let text = to_step_solid(&m, h).expect("export one");
            let (model, report) = read(text.as_bytes()).expect("re-read");
            assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
            let scene = model.scene();
            let solids: Vec<_> = scene.all_solids().collect();
            assert_eq!(solids.len(), 1);
            assert_eq!(solids[0].faces().count(), 6);
        }
    }

    #[test]
    fn hollow_solid_round_trips_as_brep_with_voids() {
        // A 4-cube with a concentric 2-cube void, built with the cavity
        // producer's primitive (`reversed_shell`): the inner shell reversed
        // inward. Exports as a BREP_WITH_VOIDS and reads back with one cavity.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let b_outer = m.solids.get(b).outer;
        let void = m.reversed_shell(b_outer);
        let a_outer = m.solids.get(a).outer;
        let hollow = m.push_solid(Solid {
            outer: a_outer,
            cavities: vec![void],
        });
        m.live_solids.retain(|&s| s == hollow); // supersede the two source cubes

        let text = to_step(&m).expect("export hollow");
        assert!(text.contains("BREP_WITH_VOIDS"), "no void entity in output");

        let (model, report) = read(text.as_bytes()).expect("re-read");
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
        let scene = model.scene();
        let solids: Vec<_> = scene.all_solids().collect();
        assert_eq!(solids.len(), 1);
        assert_eq!(solids[0].faces().count(), 6); // outer shell
        let voids = solids[0].voids();
        assert_eq!(voids.len(), 1);
        assert_eq!(voids[0].len(), 6); // one cavity shell, 6 faces
    }

    #[test]
    fn cylinder_round_trips_through_step_io_reader() {
        let mut m = Model::new();
        m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        let text = to_step(&m).expect("export cylinder");

        for needle in ["CYLINDRICAL_SURFACE", "CIRCLE", "442 3 1 4"] {
            assert!(text.contains(needle), "missing {needle}");
        }

        let (model, report) = read(text.as_bytes()).expect("re-read");
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);

        let scene = model.scene();
        let solids: Vec<_> = scene.all_solids().collect();
        assert_eq!(solids.len(), 1);
        let faces: Vec<_> = solids[0].faces().collect();
        assert_eq!(faces.len(), 3);

        let (mut cylindrical, mut planes) = (0, 0);
        for face in &faces {
            match face.surface().kind() {
                SurfaceKind::Cylindrical(_) => cylindrical += 1,
                SurfaceKind::Plane(_) => planes += 1,
                other => panic!("unexpected surface kind: {other:?}"),
            }
            for bound in face.bounds() {
                for (edge, _forward) in bound.oriented_edges() {
                    assert!(matches!(
                        edge.curve().kind(),
                        CurveKind::Line(_) | CurveKind::Circle(_)
                    ));
                }
            }
        }
        assert_eq!((cylindrical, planes), (1, 2));

        // The lateral (cylindrical) face's loop reuses the seam edge → 4 oriented edges.
        let lateral = faces
            .iter()
            .find(|f| matches!(f.surface().kind(), SurfaceKind::Cylindrical(_)))
            .unwrap();
        let edges: Vec<_> = lateral.bounds().next().unwrap().oriented_edges().collect();
        assert_eq!(edges.len(), 4);
    }
}
