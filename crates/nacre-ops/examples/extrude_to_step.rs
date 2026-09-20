//! End-to-end demo: extrude a hexagon into a prism, then export it two ways.
//!
//! Drives the full M1+M2 stack — `nacre-ops` (replay a one-op log) → `nacre-validate`
//! (assert a clean closed b-rep) → `nacre-step`/`nacre-tess` (STEP + OBJ) — and
//! writes both files for inspection.
//!
//! `cargo run -p nacre-ops --example extrude_to_step [prefix]`
//! writes `<prefix>.step` and `<prefix>.obj` (default prefix `target/hex_prism`). Open the
//! STEP in step-loupe (its report flags dropped/orphan entities) or FreeCAD (an
//! independent OCCT reader — cross-check face orientation), and the OBJ in Quick Look.
//!
//! A regular hexagon, for no reason but familiarity: the tessellator sweeps every face, holes
//! and all, so a concave profile would mesh correctly too.
//!
//! ★ This used to read *"`to_obj` ear-clips planar faces and bridges their holes"* — both halves
//! died when the monotone sweep replaced ear clipping (it **never merges rings**, which is the
//! whole reason bridging went), and the comment outlived them by two cells.

use nacre_exact::Axis;
use nacre_math::Point2;
use nacre_ops::SketchFrame;
use nacre_ops::{Operation, Profile2d, replay};
use nacre_topo::Model;
use std::f64::consts::TAU;

/// A regular hexagon of the given radius, centred on the sketch origin.
fn hexagon(r: f64) -> Profile2d {
    let points = (0..6)
        .map(|i| {
            let a = TAU * (i as f64) / 6.0;
            Point2::from_array([r * a.cos(), r * a.sin()])
        })
        .collect();
    Profile2d::polygon(points).unwrap()
}

fn main() {
    let op = Operation::Extrude {
        frame: SketchFrame::world(&Model::new(), Axis::Z),
        profile: hexagon(10.0),
        dist: 5.0,
    };
    let model = replay(&[op]).expect("replay hexagon extrude");

    // The demo must always produce a valid closed solid — surface any regression loudly.
    let violations = nacre_validate::validate(&model);
    assert!(violations.is_empty(), "invalid model: {violations:?}");

    let prefix = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/hex_prism".to_string());

    let step = nacre_step::to_step(&model).expect("export to STEP");
    let step_path = format!("{prefix}.step");
    std::fs::write(&step_path, step).expect("write STEP file");

    let obj = nacre_tess::tessellate(&model, &nacre_tess::TessConfig::default())
        .expect("the model meshes")
        .to_obj();
    let obj_path = format!("{prefix}.obj");
    std::fs::write(&obj_path, obj).expect("write OBJ file");

    println!("wrote {step_path} — open in step-loupe or FreeCAD");
    println!("wrote {obj_path} — open in Quick Look");
}
