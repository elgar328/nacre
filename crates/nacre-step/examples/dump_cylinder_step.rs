//! Build a cylinder and write it to a STEP (AP242 Ed2) file for inspection.
//!
//! `cargo run -p nacre-step --example dump_cylinder_step [path]`
//! then open the file in step-loupe (its report flags dropped/orphan entities)
//! or FreeCAD (an independent OCCT reader — cross-check the seam and outward
//! face orientation, which the step-io round-trip test cannot verify).

use nacre_math::{Point3, Vector3};
use nacre_step::to_step;
use nacre_topo::Model;

fn main() {
    let mut model = Model::new();
    // Radius 10, height 20, along +Z from the origin.
    model.add_cylinder(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        10.0,
        20.0,
    );

    let step = to_step(&model).expect("export cylinder to STEP");
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/cylinder.step".to_string());
    std::fs::write(&path, step).expect("write STEP file");
    println!("wrote {path} — open in step-loupe or FreeCAD");
}
