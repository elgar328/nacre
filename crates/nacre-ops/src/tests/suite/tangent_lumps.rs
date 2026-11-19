//! **Lumps that meet only along tangent lines** — the scoreboard of booleans whose kept material
//! near a tangent line is the two wedges between the cylinder and the wall it touches. Each family
//! is run as Fuse, both Cuts and Common:
//!
//! - a plate less the pin tangent to three of its walls: the material is three lumps, the two
//!   corners each meeting the rest only along the line where the pin touches a wall;
//! - the same with a lid fused over the plate, wider than it: the lid joins the three lumps, one
//!   body whose lines of contact are each an edge four of its faces use;
//! - the same plate stepped down at its far end, the step's top crossing the lines midway: three
//!   lumps again;
//! - a plate less a pin tangent to one wall, the control: the two wedges join round the pin, one
//!   body touching itself along the line;
//! - a slab less a rod tangent to two opposite walls and longer than the slab: two lumps, the lines
//!   ending where the slab's faces end while the rod runs on. In its union the slab's top less the
//!   rod's disk is two faces meeting at the two points of contact, and so is its bottom.

use super::scoreboard::{OPS, Want, check, counted, refused, run};
use super::*;

type Build = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);

fn p(x: f64, y: f64) -> nacre_math::Point2 {
    nacre_math::Point2::from_array([x, y])
}

/// The prism over `profile`, from `z = z0` up `dist`.
fn sweep_from(m: &mut Model, profile: Profile2d, z0: f64, dist: f64) -> Handle<Solid> {
    let frame = if z0 == 0.0 {
        SketchFrame::world(m, Axis::Z)
    } else {
        datum_frame(
            m,
            SketchPlane::from_axes(
                Point3::from_array([0.0, 0.0, z0]),
                nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
                nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
            ),
        )
    };
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("the profile extrudes") else {
        unreachable!()
    };
    solid
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Profile2d {
    Profile2d::polygon(vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]).unwrap()
}

/// The plate `[0, 4] × [0, 2]`, two tall, and the unit pin about `(3, 1)` as tall: the pin touches
/// the walls `y = 0`, `y = 2` and `x = 4`.
fn plate_and_pin_on_three_walls(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep_from(m, rect(0.0, 0.0, 4.0, 2.0), 0.0, 2.0),
        sweep_from(m, circle_profile([3.0, 1.0], 1.0), 0.0, 2.0),
    )
}

/// The slab `[−1, 1] × [−2, 2]` on `z ∈ [0, 2]`, and the unit rod about the `z` axis on
/// `z ∈ [−1, 3]`: the rod touches the walls `x = ±1` and runs past the slab both ways.
fn slab_and_rod_on_two_walls(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep_from(m, rect(-1.0, -2.0, 1.0, 2.0), 0.0, 2.0),
        sweep_from(m, circle_profile([0.0, 0.0], 1.0), -1.0, 4.0),
    )
}

/// The plate `[0, 4] × [0, 4]`, two tall, and the unit pin about `(2, 3)` as tall: the pin touches
/// the wall `y = 4` alone.
fn plate_and_pin_on_one_wall(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep_from(m, rect(0.0, 0.0, 4.0, 4.0), 0.0, 2.0),
        sweep_from(m, circle_profile([2.0, 3.0], 1.0), 0.0, 2.0),
    )
}

/// [`plate_and_pin_on_three_walls`] with a lid `[−1, 5] × [−1, 3]` on `z ∈ [2, 3]` fused on: the
/// lid joins the three lumps over the lines' upper ends.
fn lidded_plate_and_pin_on_three_walls(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let (plate, pin) = plate_and_pin_on_three_walls(m);
    let lid = sweep_from(m, rect(-1.0, -1.0, 5.0, 3.0), 2.0, 1.0);
    let lidded = crate::boolean(m, BoolKind::Fuse, plate, lid).expect("the lid")[0];
    (lidded, pin)
}

/// [`plate_and_pin_on_three_walls`] stepped down to one tall on `x ∈ [0, 1]`: the step's top
/// `z = 1` crosses the three lines midway.
fn stepped_plate_and_pin_on_three_walls(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let low = sweep_from(m, rect(0.0, 0.0, 1.0, 2.0), 0.0, 1.0);
    let high = sweep_from(m, rect(1.0, 0.0, 4.0, 2.0), 0.0, 2.0);
    let stepped = crate::boolean(m, BoolKind::Fuse, low, high).expect("the step")[0];
    (
        stepped,
        sweep_from(m, circle_profile([3.0, 1.0], 1.0), 0.0, 2.0),
    )
}

const FAMILIES: [(&str, Build); 4] = [
    (
        "a lidded plate and the pin on three of its walls",
        lidded_plate_and_pin_on_three_walls,
    ),
    (
        "a stepped plate and the pin on three of its walls",
        stepped_plate_and_pin_on_three_walls,
    ),
    (
        "a plate and the pin on three of its walls",
        plate_and_pin_on_three_walls,
    ),
    (
        "a plate and the pin on one of its walls",
        plate_and_pin_on_one_wall,
    ),
];

const PI: f64 = std::f64::consts::PI;

/// Each family's answers as Fuse, the first less the second, the second less the first, Common.
fn wants(name: &str) -> [Want; 4] {
    let none = || counted(&[], &[]);
    match name {
        // The pin is inside the plate: the union is the plate, the common the pin. The plate less
        // the pin is three lumps meeting along the three lines of contact — the left one of six
        // faces and the two corners of five, `2(6 − π/2)` and `2(1 − π/4)` twice, as OCCT builds
        // them (sixteen faces, `2(8 − π)`): the line of contact is an edge of both lumps.
        "a plate and the pin on three of its walls" => [
            counted(&[16.0], &[6]),
            counted(&[2.0 - PI / 2.0, 2.0 - PI / 2.0, 12.0 - PI], &[5, 5, 6]),
            none(),
            counted(&[2.0 * PI], &[3]),
        ],
        // The two wedges beside the line of contact join round the pin: one body touching itself.
        "a plate and the pin on one of its walls" => [
            counted(&[32.0], &[6]),
            refused(RejectReason::SelfTouchingResult),
            none(),
            counted(&[2.0 * PI], &[3]),
        ],
        // The lid joins the three lumps over the lines' upper ends: one body whose three lines of
        // contact are each an edge four of its faces use (OCCT: one solid of twenty faces).
        "a lidded plate and the pin on three of its walls" => [
            counted(&[40.0], &[11]),
            refused(RejectReason::NonManifoldResultEdge),
            none(),
            counted(&[2.0 * PI], &[3]),
        ],
        // The step's top crosses the three lines midway and joins nothing: three lumps as on the
        // flat plate, the left one `10 − π` of eight faces (OCCT: the same volumes).
        "a stepped plate and the pin on three of its walls" => [
            counted(&[14.0], &[8]),
            counted(&[2.0 - PI / 2.0, 2.0 - PI / 2.0, 10.0 - PI], &[5, 5, 8]),
            none(),
            counted(&[2.0 * PI], &[3]),
        ],
        other => panic!("no answers stated for {other}"),
    }
}

/// Every family, as each of the four booleans, comes out as [`check`] reads its answer.
#[test]
fn lumps_meeting_only_on_tangent_lines() {
    for (name, build) in FAMILIES {
        for (op, want) in OPS.into_iter().zip(wants(name)) {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}

/// **A slab less a rod tangent to two of its walls and longer than it.** The union is one body of
/// twelve faces: the slab's top and bottom are two faces each, meeting at the points where the
/// rod's circle touches the walls, and the walls stay whole. OCCT builds fourteen, a count that
/// fits it also parting each wall the rod touches along the line of contact. The slab less the
/// rod is two lumps meeting along the two lines, each `2(4 − π/2)` of six faces, as OCCT builds
/// them, the lines ending where the slab's faces end. The rod less the slab is its two ends, and
/// the common the rod's middle.
#[test]
fn a_slab_and_the_rod_on_two_of_its_walls() {
    use super::scoreboard::Op;
    let every_boolean = [
        (Op::Fuse, counted(&[16.0 + 2.0 * PI], &[12])),
        (Op::CutAb, counted(&[8.0 - PI, 8.0 - PI], &[6, 6])),
        (Op::CutBa, counted(&[PI, PI], &[3, 3])),
        (Op::Common, counted(&[2.0 * PI], &[3])),
    ];
    for (op, want) in every_boolean {
        let mut m = Model::new();
        let (a, b) = slab_and_rod_on_two_walls(&mut m);
        let got = run(&mut m, op, a, b);
        check(
            &format!("a slab and the rod on two of its walls {op:?}"),
            &m,
            got,
            &want,
        );
    }
}
