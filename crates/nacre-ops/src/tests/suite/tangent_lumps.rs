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
//!   rod's disk is two faces meeting at the two points of contact, and so is its bottom;
//! - a block bored by a pin tangent to the plane `y = 0`, less a slab beyond that plane under the
//!   block's lower half — the slab's face there ending on the line — with or without a wedge on the
//!   slab whose walls run from the line: above the slab the plane has no face on the line, and the
//!   wedge's walls are a station's chords there. The block's top joins the line's two sides over
//!   the bore, so it touches itself; notched to the axis's height, with the wedge it is two bodies;
//! - that notched block with a shorter rod of its pin fused back into the notch — below the block's
//!   floor, past the slab's top, inside — and the slab, with or without the wedge: the block's top
//!   runs along the ruling where the lateral across it passes from the notch's to the rod's. With
//!   it, the filled block and a far box, and a blind notch half filled, whose top's hole is bounded
//!   by the same ruling.

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

/// [`slab_and_rod_on_two_walls`] turned off the axes: the slab between the walls `3x + 4y = ±5`
/// (the parallelogram `(−3, 1), (1, −2), (3, −1), (−1, 2)`, area ten) on `z ∈ [0, 2]`, and the same
/// rod, which touches the walls at `±(3/5, 4/5)` — points no axis-aligned sample of its circle
/// lands on.
fn slab_and_rod_on_two_slanted_walls(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let quad = Profile2d::polygon(vec![p(-3.0, 1.0), p(1.0, -2.0), p(3.0, -1.0), p(-1.0, 2.0)])
        .expect("a parallelogram");
    (
        sweep_from(m, quad, 0.0, 2.0),
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

/// The block `[−1, 3] × [−2, top]`, two tall, less the unit pin about `(1, 1)` over `pin` —
/// tangent to the plane `y = 0` along `(1, 0)` — with the pin's rod over `rod` fused back in.
fn bored_block(m: &mut Model, top: f64, pin: [f64; 2], rod: Option<[f64; 2]>) -> Handle<Solid> {
    let block = sweep_from(m, rect(-1.0, -2.0, 3.0, top), 0.0, 2.0);
    let cutter = sweep_from(m, circle_profile([1.0, 1.0], 1.0), pin[0], pin[1] - pin[0]);
    let bored = crate::boolean(m, BoolKind::Cut, block, cutter).expect("the bore")[0];
    let Some([z0, z1]) = rod else {
        return bored;
    };
    let filler = sweep_from(m, circle_profile([1.0, 1.0], 1.0), z0, z1 - z0);
    crate::boolean(m, BoolKind::Fuse, bored, filler).expect("the rod")[0]
}

/// Beyond the plane `y = 0` the slab `[−2, 4] × [−3, 0]` on `z ∈ [−1, 1]`, with the wedge
/// `(1, 0), (0, −3), (2, −3)` on `z ∈ [1, 3]` fused on it where `wedge`: its walls `3x ± y = 3`
/// run from the line down through a [`bored_block`]'s floor `y = −2`.
fn slab_beyond(m: &mut Model, wedge: bool) -> Handle<Solid> {
    let slab = sweep_from(m, rect(-2.0, -3.0, 4.0, 0.0), -1.0, 2.0);
    if !wedge {
        return slab;
    }
    let walls = Profile2d::polygon(vec![p(1.0, 0.0), p(0.0, -3.0), p(2.0, -3.0)]).expect("a wedge");
    let w = sweep_from(m, walls, 1.0, 2.0);
    crate::boolean(m, BoolKind::Fuse, slab, w).expect("the wedge")[0]
}

/// A [`bored_block`] through its height and the [`slab_beyond`] it. The block less the slab
/// keeps the two cusps beside the line below `z = 1`, and the slab's face `y = 0` ends there;
/// above it the wall's plane has no face on the line, and the wedge empties a room between its
/// walls below it.
fn bored_block_and_slab(m: &mut Model, top: f64, wedge: bool) -> (Handle<Solid>, Handle<Solid>) {
    (
        bored_block(m, top, [-1.0, 3.0], None),
        slab_beyond(m, wedge),
    )
}

/// [`bored_block_and_slab`] three high with the wedge: the block's top joins the two sides of the
/// line over the bore.
fn bored_block_and_slab_with_a_wedge(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    bored_block_and_slab(m, 3.0, true)
}

/// [`bored_block_and_slab`] three high without the wedge.
fn bored_block_and_slab_alone(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    bored_block_and_slab(m, 3.0, false)
}

/// [`bored_block_and_slab`] up to the axis's height with the wedge: the pin notches the block's
/// top, and nothing joins the two sides of the line above it.
fn notched_block_and_slab_with_a_wedge(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    bored_block_and_slab(m, 1.0, true)
}

/// [`bored_block_and_slab`] up to the axis's height without the wedge.
fn notched_block_and_slab_alone(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    bored_block_and_slab(m, 1.0, false)
}

/// The notched block with the pin's rod over `rod` fused back into the notch, and the slab: the
/// block's top `y = 1` runs along the ruling `x = 2` where the lateral across it passes from the
/// notch's to the rod's and back, and the rod's caps' rims pass there from one to the other.
fn filled_notch_and_slab(
    m: &mut Model,
    rod: [f64; 2],
    wedge: bool,
) -> (Handle<Solid>, Handle<Solid>) {
    (
        bored_block(m, 1.0, [-1.0, 3.0], Some(rod)),
        slab_beyond(m, wedge),
    )
}

fn filled_low_and_slab(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    filled_notch_and_slab(m, [-0.5, 0.7], false)
}

fn filled_low_and_slab_with_a_wedge(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    filled_notch_and_slab(m, [-0.5, 0.7], true)
}

fn filled_high_and_slab(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    filled_notch_and_slab(m, [0.3, 1.5], false)
}

fn filled_high_and_slab_with_a_wedge(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    filled_notch_and_slab(m, [0.3, 1.5], true)
}

fn filled_inside_and_slab(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    filled_notch_and_slab(m, [0.5, 1.0], false)
}

fn filled_inside_and_slab_with_a_wedge(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    filled_notch_and_slab(m, [0.5, 1.0], true)
}

/// The notch filled over `[1/2, 1]` and a box far from it: what this family asks is only whether
/// the filled block is read as an operand.
fn filled_inside_and_a_far_box(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let filled = bored_block(m, 1.0, [-1.0, 3.0], Some([0.5, 1.0]));
    let far = crate::fixtures::cuboid(m, [10.0, 10.0, 10.0], [11.0, 11.0, 11.0]);
    (filled, far)
}

/// The block notched by the pin over `[1/2, 3/2]` only — the ruling bounds a hole in the top —
/// with the rod over `[1/2, 1]` filling the notch's lower half, and the slab.
fn blind_notch_half_filled_and_slab(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        bored_block(m, 1.0, [0.5, 1.5], Some([0.5, 1.0])),
        slab_beyond(m, false),
    )
}

const FAMILIES: [(&str, Build); 16] = [
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
    (
        "a bored block and a slab with a wedge beyond its tangent wall",
        bored_block_and_slab_with_a_wedge,
    ),
    (
        "a bored block and a slab beyond its tangent wall",
        bored_block_and_slab_alone,
    ),
    (
        "a notched block and a slab with a wedge beyond its tangent wall",
        notched_block_and_slab_with_a_wedge,
    ),
    (
        "a notched block and a slab beyond its tangent wall",
        notched_block_and_slab_alone,
    ),
    (
        "a notch filled below the block's floor and a slab",
        filled_low_and_slab,
    ),
    (
        "a notch filled below the block's floor and a slab with a wedge",
        filled_low_and_slab_with_a_wedge,
    ),
    (
        "a notch filled past the slab's top and a slab",
        filled_high_and_slab,
    ),
    (
        "a notch filled past the slab's top and a slab with a wedge",
        filled_high_and_slab_with_a_wedge,
    ),
    ("a notch filled inside and a slab", filled_inside_and_slab),
    (
        "a notch filled inside and a slab with a wedge",
        filled_inside_and_slab_with_a_wedge,
    ),
    (
        "a notch filled inside and a far box",
        filled_inside_and_a_far_box,
    ),
    (
        "a blind notch half filled and a slab",
        blind_notch_half_filled_and_slab,
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
        // The block less the slab and the wedge parts along the line, below the slab's top by the
        // knife and above it by the wedge's walls, and joins over the bore: one body touching
        // itself along the line, which four of its faces use. OCCT 7.9 gives the same volumes for
        // the other three booleans and builds the cut.
        "a bored block and a slab with a wedge beyond its tangent wall" => [
            counted(&[82.0 - 2.0 * PI - 28.0 / 3.0], &[18]),
            refused(RejectReason::NonManifoldResultEdge),
            counted(&[98.0 / 3.0], &[15]),
            counted(&[28.0 / 3.0], &[10]),
        ],
        // Without the wedge the cusps the slab leaves beside the line join round the bore and above
        // the slab: one body touching itself along the line. OCCT 7.9 gives the same volumes for
        // the other three booleans and builds the cut.
        "a bored block and a slab beyond its tangent wall" => [
            counted(&[68.0 - 2.0 * PI], &[13]),
            refused(RejectReason::SelfTouchingResult),
            counted(&[28.0], &[10]),
            counted(&[8.0], &[6]),
        ],
        // Notched to the axis's height, nothing joins the line's two sides above it, and the knife
        // below the slab's top and the wedge's station above it cut the lateral along the whole
        // line: two bodies. OCCT 7.9: the same volume, eighteen faces.
        "a notched block and a slab with a wedge beyond its tangent wall" => [
            counted(&[66.0 - PI - 28.0 / 3.0], &[21]),
            counted(
                &[
                    (24.0 - PI - 28.0 / 3.0) / 2.0,
                    (24.0 - PI - 28.0 / 3.0) / 2.0,
                ],
                &[9, 9],
            ),
            counted(&[98.0 / 3.0], &[15]),
            counted(&[28.0 / 3.0], &[10]),
        ],
        // Without the wedge the block runs on round the line above the slab, joining the two sides
        // the slab's cut leaves beside it: one body touching itself along the line. OCCT 7.9 gives
        // the same volumes for the other three booleans and builds the cut.
        "a notched block and a slab beyond its tangent wall" => [
            counted(&[52.0 - PI], &[15]),
            refused(RejectReason::NonManifoldResultEdge),
            counted(&[28.0], &[10]),
            counted(&[8.0], &[6]),
        ],
        // The rod fused into the notch over `[z0, z1]` adds its whole disk there: the filled block
        // is `24 − π + π(z1 − z0)`, and it meets the slab in the block's lower half, `8`, or with
        // the wedge also its tip in the block, `28/3`; the slab alone is `36`, with the wedge `42`.
        // Where the rod leaves part of the notch below the slab's top empty, the cut keeps the two
        // cusps there touching along the line and joined through the rod: one body touching
        // itself. Below the block the rod touches the slab along the line, so the union touches
        // itself too. The slab families' volumes match OCCT 7.9, and an independent reading of
        // the sections finds the bodies counted and the touches refused.
        "a notch filled below the block's floor and a slab" => [
            refused(RejectReason::SelfTouchingResult),
            refused(RejectReason::SelfTouchingResult),
            counted(&[28.0], &[10]),
            counted(&[8.0], &[6]),
        ],
        "a notch filled below the block's floor and a slab with a wedge" => [
            refused(RejectReason::SelfTouchingResult),
            refused(RejectReason::SelfTouchingResult),
            counted(&[98.0 / 3.0], &[15]),
            counted(&[28.0 / 3.0], &[10]),
        ],
        "a notch filled past the slab's top and a slab" => [
            counted(&[52.0 + 0.2 * PI], &[19]),
            refused(RejectReason::SelfTouchingResult),
            counted(&[28.0], &[10]),
            counted(&[8.0], &[6]),
        ],
        "a notch filled past the slab's top and a slab with a wedge" => [
            counted(&[66.0 + 0.2 * PI - 28.0 / 3.0], &[25]),
            refused(RejectReason::SelfTouchingResult),
            counted(&[98.0 / 3.0], &[15]),
            counted(&[28.0 / 3.0], &[10]),
        ],
        "a notch filled inside and a slab" => [
            counted(&[52.0 - PI / 2.0], &[19]),
            refused(RejectReason::SelfTouchingResult),
            counted(&[28.0], &[10]),
            counted(&[8.0], &[6]),
        ],
        "a notch filled inside and a slab with a wedge" => [
            counted(&[66.0 - PI / 2.0 - 28.0 / 3.0], &[25]),
            refused(RejectReason::SelfTouchingResult),
            counted(&[98.0 / 3.0], &[15]),
            counted(&[28.0 / 3.0], &[10]),
        ],
        "a notch filled inside and a far box" => [
            counted(&[1.0, 24.0 - PI / 2.0], &[6, 12]),
            counted(&[24.0 - PI / 2.0], &[12]),
            counted(&[1.0], &[6]),
            none(),
        ],
        // The notch over `[1/2, 3/2]` takes `π/2` and the rod gives it back: the block is `24`.
        // Its cut keeps the rod, which fills the notch up to the slab's top, and the notch empty
        // above it: the two rooms beside the line — the slab's below `z = 1`, the notch's above —
        // meet only at `(1, 0, 1)`, and the vertex is not a manifold's.
        "a blind notch half filled and a slab" => [
            counted(&[52.0], &[17]),
            refused(RejectReason::NonManifoldVertex),
            counted(&[28.0], &[10]),
            counted(&[8.0], &[6]),
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

/// **[`a_slab_and_the_rod_on_two_of_its_walls`] between slanted walls** — the same faces and lumps,
/// the slab's area ten. The touch points sit where the rod's circle meets the walls at no axis
/// angle, so a cut the arrangement missed there is not one a mesh sample would land on: what holds
/// the top and bottom in two faces each is the touch being a vertex, and the result contract's
/// face-bounds arm (`FaceBoundsMeet`) is what refuses it if not.
#[test]
fn a_slab_and_the_rod_on_two_of_its_slanted_walls() {
    use super::scoreboard::Op;
    let every_boolean = [
        (Op::Fuse, counted(&[20.0 + 2.0 * PI], &[12])),
        (Op::CutAb, counted(&[10.0 - PI, 10.0 - PI], &[6, 6])),
        (Op::CutBa, counted(&[PI, PI], &[3, 3])),
        (Op::Common, counted(&[2.0 * PI], &[3])),
    ];
    for (op, want) in every_boolean {
        let mut m = Model::new();
        let (a, b) = slab_and_rod_on_two_slanted_walls(&mut m);
        let got = run(&mut m, op, a, b);
        check(
            &format!("a slab and the rod on two of its slanted walls {op:?}"),
            &m,
            got,
            &want,
        );
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
