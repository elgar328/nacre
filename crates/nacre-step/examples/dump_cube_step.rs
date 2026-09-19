//! Build a unit cube and write it to a STEP (AP242 Ed2) file for inspection.
//!
//! `cargo run -p nacre-step --example dump_cube_step [path]`
//! then open the file in step-loupe (its report flags dropped/orphan entities)
//! or FreeCAD (an independent OCCT reader — cross-check face orientation).

use nacre_math::Point3;
use nacre_step::to_step;
use nacre_topo::Model;

fn main() {
    let mut model = Model::new();
    model.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );

    let step = to_step(&model).expect("export cube to STEP");
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/cube.step".to_string());
    std::fs::write(&path, step).expect("write STEP file");
    println!("wrote {path} — open in step-loupe or FreeCAD");
}
