//! Build a unit cube and write it to an OBJ file for visual inspection.
//!
//! `cargo run -p nacre-tess --example dump_cube [path]`
//! then open the file in MeshLab / f3d / any OBJ viewer.

use nacre_math::Point3;
use nacre_tess::to_obj;
use nacre_topo::Model;

fn main() {
    let mut model = Model::new();
    model.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );

    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "cube.obj".to_string());
    let obj = to_obj(&model).expect("planar model meshes");
    std::fs::write(&path, obj).expect("write OBJ file");
    println!("wrote {path} — open in MeshLab / f3d / any OBJ viewer");
}
