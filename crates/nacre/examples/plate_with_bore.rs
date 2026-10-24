//! A plate with a hole bored through it, checked and written out two ways.
//!
//! Run from the workspace root: `cargo run -p nacre --example plate_with_bore [prefix]` writes
//! `<prefix>.step` and `<prefix>.obj` (default prefix `target/plate_with_bore`, relative to the
//! current directory). Open the STEP in a STEP viewer and the OBJ in a mesh viewer (Quick Look
//! on macOS).

use nacre::prelude::*;

fn main() {
    let mut model = Model::new();
    // Every solid here is a profile on the world XY plane, extruded along +z.
    let extrude = |model: &mut Model, profile: Profile2d, dist: f64| {
        let frame = SketchFrame::world(model, Axis::Z);
        match apply(
            model,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        ) {
            Ok(OpOutput::Extrude { solid, .. }) => solid,
            other => panic!("the extrude builds: {other:?}"),
        }
    };

    // A 40 x 30 plate, 5 thick.
    let rectangle = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([40.0, 0.0]),
        Point2::from_array([40.0, 30.0]),
        Point2::from_array([0.0, 30.0]),
    ])
    .expect("a rectangle");
    let plate = extrude(&mut model, rectangle, 5.0);

    // A cylinder of radius 6 through the plate's middle, as tall as the plate.
    let circle = Ring2d::circle(Point2::from_array([20.0, 15.0]), 6.0).expect("a circle");
    let disk = from_paths(vec![circle]).expect("one profile").remove(0);
    let bore = extrude(&mut model, disk, 5.0);

    // Cut the cylinder out of the plate. Both operands are consumed; a boolean returns every
    // solid its result falls into, and here that is the one bored plate.
    let result = boolean(&mut model, BoolKind::Cut, plate, bore).expect("the bore cuts");
    let [bored] = result[..] else {
        panic!("one solid: {result:?}")
    };
    let violations = validate(&model);
    assert!(
        violations.is_empty(),
        "a valid closed b-rep: {violations:?}"
    );
    let volume = mass_props(&model, bored).expect("mass properties").volume;
    let want = 40.0 * 30.0 * 5.0 - std::f64::consts::PI * 6.0 * 6.0 * 5.0;
    assert!((volume - want).abs() < 1e-9 * want, "{volume} vs {want}");

    let prefix = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/plate_with_bore".to_string());
    // The header's time stamp is the caller's to give: the kernel reads no clock.
    let step = nacre::step::to_step(&model, "2026-10-24T00:00:00Z").expect("the STEP export");
    let step_path = format!("{prefix}.step");
    std::fs::write(&step_path, step).expect("write the STEP file");
    let obj = tessellate(&model, &TessConfig::default())
        .expect("the mesh")
        .to_obj(&model);
    let obj_path = format!("{prefix}.obj");
    std::fs::write(&obj_path, obj).expect("write the OBJ file");
    println!("volume {volume:.3}; wrote {step_path} and {obj_path}");
}
