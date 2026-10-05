//! Build a cylinder and write it to a STEP (AP242 Ed2) file for inspection.
//!
//! `cargo run -p nacre-step --example dump_cylinder_step [path]`
//! then open the file in step-loupe (its report flags dropped/orphan entities)
//! or FreeCAD (a reader outside this repo, though OCCT underneath — cross-check the
//! two-rim lateral and outward face orientation, which the step-io round-trip test
//! cannot verify).

use nacre_math::Point2;
use nacre_ops::{DatumDef, OpOutput, Operation, Ring2d, SketchPlane, apply, from_paths};
use nacre_step::to_step;
use nacre_topo::Model;

fn main() {
    let mut model = Model::new();
    // Radius 10, height 20, along +Z from the origin — stated the way an application states it:
    // a sketch plane, a whole circle on it, an extrude.
    let Ok(OpOutput::DatumPlane { frame, .. }) = apply(
        &mut model,
        &Operation::DatumPlane {
            def: DatumDef::Stated(SketchPlane::world_xy()),
        },
    ) else {
        panic!("the XY plane is stated")
    };
    let circle = Ring2d::circle(Point2::from_array([0.0, 0.0]), 10.0).expect("a circle");
    let profile = from_paths(vec![circle]).expect("one profile").remove(0);
    apply(
        &mut model,
        &Operation::Extrude {
            frame,
            profile,
            dist: 20.0,
        },
    )
    .expect("the circle extrudes");

    // The header time stamp is the caller's to give (the kernel reads no clock).
    let step = to_step(&model, "2026-10-05T00:00:00Z").expect("export cylinder to STEP");
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/cylinder.step".to_string());
    std::fs::write(&path, step).expect("write STEP file");
    println!("wrote {path} — open in step-loupe or FreeCAD");
}
