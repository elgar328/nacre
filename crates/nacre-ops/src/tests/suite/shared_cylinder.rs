//! **One cylinder surface under two solids** — the scoreboard of booleans whose two operands both
//! carry lateral faces on one cylinder handle.
//!
//! Each family is built in one frame from one circle, so the two operands' laterals intern to one
//! `Handle<Surface>` (asserted per family before anything else): two identical prisms, a circle
//! swept both ways, a plate bored by the circle and the circle's pin, a short and a tall prism, a
//! half disk beside the whole one (a ruling station on the shared class), a notched plate over the
//! half pin below the notch (the material side flips across the station), a filleted corner swept
//! both ways and the fillet's pin (tangent walls on the shared class), two plates rounded on
//! opposite quarters of one circle, and two identical cylinders on a stated datum. Every family is
//! run as Fuse, both Cuts and Common.
//!
//! The gate refuses all of them as `CylinderPairContact` — one class carrying rows of both solids
//! (`planes::cylinder_gate`), because the chart reads a cylinder class as one solid's surface
//! (`cyl_chart::chart_of`).

use super::*;

type Build = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);

/// The four booleans each family is run as.
#[derive(Clone, Copy, Debug)]
enum Op {
    Fuse,
    /// The first operand minus the second.
    CutAb,
    /// The second operand minus the first.
    CutBa,
    Common,
}

const OPS: [Op; 4] = [Op::Fuse, Op::CutAb, Op::CutBa, Op::Common];

fn run(m: &mut Model, op: Op, a: Handle<Solid>, b: Handle<Solid>) -> BoolResult {
    match op {
        Op::Fuse => crate::boolean(m, BoolKind::Fuse, a, b),
        Op::CutAb => crate::boolean(m, BoolKind::Cut, a, b),
        Op::CutBa => crate::boolean(m, BoolKind::Cut, b, a),
        Op::Common => crate::boolean(m, BoolKind::Common, a, b),
    }
}

type BoolResult = Result<Vec<Handle<Solid>>, BoolError>;

fn p(x: f64, y: f64) -> nacre_math::Point2 {
    nacre_math::Point2::from_array([x, y])
}

fn profile(edges: Vec<Stated>) -> Profile2d {
    stated(edges).unwrap().remove(0)
}

fn sweep(m: &mut Model, profile: Profile2d, dist: f64) -> Handle<Solid> {
    let frame = SketchFrame::world(m, Axis::Z);
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

/// The unit circle about the origin.
fn disk() -> Profile2d {
    circle_profile([0.0, 0.0], 1.0)
}

/// The unit circle's half on `y ≥ 0`, closed by its diameter.
fn upper_half() -> Profile2d {
    profile(vec![
        line(p(-1.0, 0.0), p(1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(1.0, 0.0), 2),
    ])
}

/// The unit circle's half on `y ≤ 0`, closed by its diameter.
fn lower_half() -> Profile2d {
    profile(vec![
        line(p(1.0, 0.0), p(-1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(-1.0, 0.0), 2),
    ])
}

/// `[0, 4] × [0, 2]` with its `(4, 0)` corner rounded by the unit circle about `(3, 1)`.
fn filleted_right() -> Profile2d {
    profile(vec![
        line(p(0.0, 0.0), p(3.0, 0.0)),
        arc_turns(p(3.0, 1.0), p(3.0, 0.0), 1),
        line(p(4.0, 1.0), p(4.0, 2.0)),
        line(p(4.0, 2.0), p(0.0, 2.0)),
        line(p(0.0, 2.0), p(0.0, 0.0)),
    ])
}

/// `[2, 6] × [0, 2]` with its `(2, 0)` corner rounded by the same circle — the opposite quarter.
fn filleted_left() -> Profile2d {
    profile(vec![
        line(p(3.0, 0.0), p(6.0, 0.0)),
        line(p(6.0, 0.0), p(6.0, 2.0)),
        line(p(6.0, 2.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(2.0, 1.0)),
        arc_turns(p(3.0, 1.0), p(2.0, 1.0), 1),
    ])
}

fn identical(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, disk(), 2.0), sweep(m, disk(), 2.0))
}

fn swept_both_ways(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, disk(), 2.0), sweep(m, disk(), -2.0))
}

/// A `4 × 4 × 2` plate bored by the unit circle, and the circle's pin of the plate's height.
fn bore_and_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let square =
        Profile2d::polygon(vec![p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0)]).unwrap();
    let plate = sweep(m, square, 2.0);
    let tool = sweep(m, disk(), 2.0);
    let bored = crate::boolean(m, BoolKind::Cut, plate, tool).expect("the bore")[0];
    (bored, sweep(m, disk(), 2.0))
}

fn short_and_tall(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, disk(), 2.0), sweep(m, disk(), 3.0))
}

fn half_and_whole(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, upper_half(), 2.0), sweep(m, disk(), 2.0))
}

fn short_half_and_tall_whole(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, upper_half(), 2.0), sweep(m, disk(), 3.0))
}

/// `[−2, 2] × [0, 2]` notched by the upper half disk, and the lower half disk's pin: the shared
/// lateral is the notch's on `y > 0` (material outside) and the pin's on `y < 0` (material inside).
fn notch_and_half_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let notched = profile(vec![
        line(p(-2.0, 0.0), p(-1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(-1.0, 0.0), -2),
        line(p(1.0, 0.0), p(2.0, 0.0)),
        line(p(2.0, 0.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(-2.0, 2.0)),
        line(p(-2.0, 2.0), p(-2.0, 0.0)),
    ]);
    (sweep(m, notched, 2.0), sweep(m, lower_half(), 2.0))
}

fn fillet_swept_both_ways(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, filleted_right(), 2.0),
        sweep(m, filleted_right(), -2.0),
    )
}

/// The filleted plate and a pin on its fillet's circle, taller than the plate: the plate's two
/// walls are tangent to the pin's lateral.
fn fillet_and_its_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, filleted_right(), 2.0),
        sweep(m, circle_profile([3.0, 1.0], 1.0), 3.0),
    )
}

fn opposite_quarters(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, filleted_right(), 2.0),
        sweep(m, filleted_left(), 2.0),
    )
}

/// Two identical cylinders on a stated datum plane — the other frame road.
fn identical_on_a_datum(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let on_datum = |m: &mut Model| {
        crate::fixtures::cylinder(
            m,
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.0,
            2.0,
        )
        .solid
    };
    (on_datum(m), on_datum(m))
}

const FAMILIES: [(&str, Build); 11] = [
    ("two identical prisms", identical),
    ("one circle swept both ways", swept_both_ways),
    ("a bored plate and the bore's pin", bore_and_pin),
    ("a short and a tall prism", short_and_tall),
    ("a half disk and the whole", half_and_whole),
    (
        "a short half disk and a tall whole",
        short_half_and_tall_whole,
    ),
    ("a notched plate and the half pin", notch_and_half_pin),
    ("a fillet swept both ways", fillet_swept_both_ways),
    ("a filleted plate and its fillet's pin", fillet_and_its_pin),
    ("two plates rounded on opposite quarters", opposite_quarters),
    ("two identical cylinders on a datum", identical_on_a_datum),
];

/// The cylinder surfaces of a solid's lateral faces.
fn lateral_surfaces(m: &Model, s: Handle<Solid>) -> Vec<Handle<Surface>> {
    let mut out: Vec<Handle<Surface>> = m
        .shell(m.solid(s).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .filter(|&h| matches!(m.surface(h), Surface::Cylinder { .. }))
        .collect();
    out.sort();
    out.dedup();
    out
}

/// **Every family shares one cylinder handle between its operands** — the population this
/// scoreboard is about, checked before its answers are read.
#[test]
fn each_familys_operands_share_one_cylinder_handle() {
    for (name, build) in FAMILIES {
        let mut m = Model::new();
        let (a, b) = build(&mut m);
        let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
        assert!(
            la.iter().any(|h| lb.contains(h)),
            "{name}: the operands' laterals are {la:?} and {lb:?}"
        );
    }
}

/// **The gate refuses every family by the pair's name.**
#[test]
fn a_cylinder_shared_by_two_solids_is_refused() {
    for (name, build) in FAMILIES {
        for op in OPS {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            let got = run(&mut m, op, a, b);
            assert!(
                matches!(
                    got,
                    Err(BoolError::Rejected {
                        reason: RejectReason::CylinderPairContact,
                        ..
                    })
                ),
                "{name} {op:?}: {got:?}"
            );
        }
    }
}
