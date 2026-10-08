//! The oracle sees orientation: a STEP written the wrong way round is refused, not repaired.

use super::*;

/// `step` with the `.T.`/`.F.` that ends the first line `pick` chooses turned over. Panics if no
/// line is chosen or the line carries no flag, so a plant that does not bite cannot pass.
fn turn_over(step: &str, what: &str, pick: impl Fn(&str) -> bool) -> String {
    let line = step
        .lines()
        .find(|l| pick(l))
        .unwrap_or_else(|| panic!("{what}: no such line"));
    let at = line
        .rfind(".T.")
        .max(line.rfind(".F."))
        .unwrap_or_else(|| panic!("{what}: no flag in {line}"));
    let flag = if &line[at..at + 3] == ".T." {
        ".F."
    } else {
        ".T."
    };
    let turned = format!("{}{flag}{}", &line[..at], &line[at + 3..]);
    let out = step.replacen(line, &turned, 1);
    assert_ne!(out, step, "{what}: the plant changed nothing");
    out
}

/// A line that states the entity `name` (`#12 = NAME(...)`).
fn states(name: &'static str) -> impl Fn(&str) -> bool {
    move |l: &str| l.contains(&format!("= {name}("))
}

#[track_caller]
fn refused_on(read: Result<OcctProps, OracleError>, shape: &str, what: &str) {
    match read {
        Err(OracleError::Faulty(faults)) => assert!(
            faults.iter().any(|f| f.shape == shape),
            "{what}: refused, but the faults {faults:?} do not name `{shape}`"
        ),
        other => panic!("{what}: expected Faulty on `{shape}`, got {other:?}"),
    }
}

/// One live solid's STEP.
fn step_of(m: &Model, s: Handle<Solid>) -> String {
    nacre_step::to_step_solids(m, &[s], crate::STAMP).expect("export")
}

/// ★★ **A face, an edge, a bound or a void written the wrong way round reads as faulty**, on a
/// curved shape too. OCCT's STEP read would turn each of them the right way round and score it
/// valid with the right volume (a cube with a face turned over still measures 1); the helper
/// turns that repair off, and the oracle refuses what `checkshape` then finds. The control is the
/// same three shapes untouched. A boolean checks its inputs, each on its own: a wrong `a` beside a
/// sound `b` is named `a`, and the other way round.
///
/// A line's `EDGE_CURVE` sense is not here: no OCCT read sees it turned over, and brep-to-step
/// writes a line's direction from its edge's end points, so nacre cannot write it wrong.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_wrong_orientation_reads_as_faulty() {
    let mut m = Model::new();
    let cube = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let tool = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5; 3]),
        Point3::from_array([1.5; 3]),
    );
    let cube_step = step_of(&m, cube);
    let tool_step = step_of(&m, tool);

    let mut c = Model::new();
    let z = Vector3::from_array([0.0, 0.0, 1.0]);
    let cylinder = nacre_ops::fixtures::cylinder(&mut c, Point3::origin(), z, 2.0, 5.0).solid;
    let cylinder_step = step_of(&c, cylinder);

    let (h, hollow) = hollow_box();
    let hollow_step = step_of(&h, hollow);

    for (what, step) in [
        ("the cube", &cube_step),
        ("the cylinder", &cylinder_step),
        ("the hollow box", &hollow_step),
    ] {
        occt_props(step).unwrap_or_else(|e| panic!("{what} untouched: {e:?}"));
    }

    // The face whose surface is the cylinder: `#n = CYLINDRICAL_SURFACE(` names it `#n`.
    let lateral = cylinder_step
        .lines()
        .find(|l| l.contains("= CYLINDRICAL_SURFACE("))
        .and_then(|l| l.split_whitespace().next())
        .expect("the cylinder's surface")
        .to_owned();
    let lateral_face =
        move |l: &str| l.contains("= ADVANCED_FACE(") && l.contains(&format!(",{lateral},"));

    for (what, turned) in [
        (
            "a face's sense",
            turn_over(&cube_step, "a face", states("ADVANCED_FACE")),
        ),
        (
            "an oriented edge",
            turn_over(&cube_step, "an edge", states("ORIENTED_EDGE")),
        ),
        (
            "a face's outer bound",
            turn_over(&cube_step, "a bound", states("FACE_OUTER_BOUND")),
        ),
        (
            "the cylinder's lateral",
            turn_over(&cylinder_step, "the lateral", lateral_face),
        ),
        (
            "the void's shell",
            turn_over(&hollow_step, "the void", states("ORIENTED_CLOSED_SHELL")),
        ),
    ] {
        refused_on(occt_props(&turned), "a", what);
    }

    let wrong_cube = turn_over(&cube_step, "a face", states("ADVANCED_FACE"));
    let wrong_tool = turn_over(&tool_step, "a face", states("ADVANCED_FACE"));
    occt_boolean(OcctBool::Cut, &cube_step, &tool_step).expect("the sound cut");
    refused_on(
        occt_boolean(OcctBool::Cut, &wrong_cube, &tool_step),
        "a",
        "a boolean's first input",
    );
    refused_on(
        occt_boolean(OcctBool::Cut, &cube_step, &wrong_tool),
        "b",
        "a boolean's second input",
    );
}
