//! Tessellate a cylinder and write it to OBJ files for inspection.
//!
//! `cargo run -p nacre-tess --example dump_cylinder_obj [prefix]`
//! writes `<prefix>_coarse.obj` and `<prefix>_fine.obj` (default prefix
//! `cylinder`) — two tolerances, to show that tolerance drives mesh density.
//! Open in Quick Look: the surface should be watertight and outward-facing, and
//! the fine mesh visibly rounder than the coarse one.

use nacre_math::{Point3, Vector3};
use nacre_tess::{TessConfig, tessellate};
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

    let prefix = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "cylinder".to_string());

    // ★ Both budgets, because either one alone can be the binding constraint: the
    // angular one decides for ordinary radii and the sagitta one for very large
    // circles. Varying only `tol` would leave this example silently showing the same
    // mesh twice.
    for (label, tol, angle) in [("coarse", 1.0, 20.0), ("fine", 0.02, 2.0)] {
        let obj = tessellate(
            &model,
            &TessConfig {
                tol,
                max_angle_deg: angle,
            },
        )
        .expect("cylinder meshes")
        .to_obj();
        let path = format!("{prefix}_{label}.obj");
        std::fs::write(&path, obj).expect("write OBJ file");
        println!("wrote {path} (tol {tol}, angle {angle}°) — open in Quick Look");
    }
}
