//! Cross-crate composition: ops (replay) → validate → step/tess export.
//!
//! The per-crate unit tests each cover one layer; this guards that the layers
//! compose. STEP format details (schema token, entity list, step_io round-trip)
//! stay owned by `nacre-step`'s own tests — here we only prove the pipeline hands
//! a valid model through to real STEP/OBJ output.

use nacre_math::{Point2, Point3};
use nacre_ops::{Operation, Profile2d, SketchPlane, replay};

/// A regular hexagon of the given radius, centred on the sketch origin.
fn hexagon(r: f64) -> Profile2d {
    let points = (0..6)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / 6.0;
            Point2::from_array([r * a.cos(), r * a.sin()])
        })
        .collect();
    Profile2d { points }
}

fn hex_extrude(plane: SketchPlane) -> Operation {
    Operation::Extrude {
        plane,
        profile: hexagon(10.0),
        dist: 5.0,
    }
}

#[test]
fn hexagon_extrude_exports_to_step_and_obj() {
    let model = replay(&[hex_extrude(SketchPlane::world_xy())]).unwrap();

    // The pipeline yields a valid closed b-rep.
    assert!(nacre_validate::validate(&model).is_empty());

    // ops model → a real STEP solid.
    let step = nacre_step::to_step(&model).expect("export to STEP");
    assert!(step.contains("MANIFOLD_SOLID_BREP"));

    // ops model → an OBJ mesh: 12 vertices (2n for a hexagon prism) plus faces.
    let obj = nacre_tess::to_obj(&model);
    let v_lines = obj.lines().filter(|l| l.starts_with("v ")).count();
    let f_lines = obj.lines().filter(|l| l.starts_with("f ")).count();
    assert_eq!(v_lines, 12);
    assert!(f_lines > 0);
}

#[test]
fn multi_op_log_composes_through_export() {
    // Two extrudes on parallel planes produce two independent solids that both
    // survive validation and export.
    let far = SketchPlane {
        origin: Point3::from_array([50.0, 0.0, 0.0]),
        ..SketchPlane::world_xy()
    };
    let model = replay(&[hex_extrude(SketchPlane::world_xy()), hex_extrude(far)]).unwrap();

    assert_eq!(model.solids.len(), 2);
    assert!(nacre_validate::validate(&model).is_empty());
    assert!(nacre_step::to_step(&model).is_ok());
}
