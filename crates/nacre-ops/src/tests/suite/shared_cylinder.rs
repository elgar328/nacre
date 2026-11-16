//! **One cylinder surface under two solids** — the scoreboard of booleans whose two operands both
//! carry lateral faces on one cylinder handle.
//!
//! Each family is built in one frame from one circle, so the two operands' laterals intern to one
//! `Handle<Surface>` (asserted per family before anything else): two identical prisms, a circle
//! swept both ways, a plate bored by the circle and the circle's pin, a short and a tall prism, a
//! half disk beside the whole one (a ruling station on the shared class), a notched plate over the
//! half pin below the notch (the material side flips across the station), a bored block above the
//! bore's pin (it flips across a rim), the two halves of a disk, a filleted corner swept
//! both ways and the fillet's pin (tangent walls on the shared class), two plates rounded on
//! opposite quarters of one circle, two identical cylinders on a stated datum, two chord segments
//! of the disk meeting on one ruling (each lateral ending there from its own side), two half disks
//! whose statements of the diameter's rulings overlap in part and one inside the other with both
//! ends inside (the rulings cut at each other's ends), and a bored block below the bore's pin (the
//! rim joined the other way round). Every family is run as Fuse, both Cuts and Common.
//!
//! The class carries both operands' lateral faces, and the chart reads its cells for both — the
//! inside chamber off the arrangement's labels, and the outside as the inside with every solid
//! whose lateral face is there flipped (`cyl_chart::read_cell`).
//!
//! The same cylinder stated under **two** handles is refused by the pair rule (`same_surface`):
//! the axis's sign, the seam direction, the origin on the line or a motion turning it about its
//! axis told apart. A motion that folds back to no motion at all is one statement, and one class.

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

/// A `4 × 4 × 2` block on `z ∈ [2, 4]` bored through by the unit circle swept from `z = 0`, and
/// the circle's pin on `z ∈ [0, 2]` below it: the shared lateral keeps the material outside above
/// the rim at `z = 2` and inside below it, and four faces meet on that circle.
fn bore_above_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, 2.0], [2.0, 2.0, 4.0]);
    let tool = sweep(m, disk(), 4.0);
    let bored = crate::boolean(m, BoolKind::Cut, block, tool).expect("the bore")[0];
    (bored, sweep(m, disk(), 2.0))
}

/// The unit disk's two halves: one lateral across the two stations of the diameter, and the two
/// flats on it facing each other.
fn two_halves(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, upper_half(), 2.0), sweep(m, lower_half(), 2.0))
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

/// The unit disk's segment beyond the chord `x + y = 1` and its segment beyond `x − y = 1`: the
/// two chords meet on the ruling `(1, 0)`, where one lateral face ends from `+θ` and the other's
/// from `−θ`, each against its own chord.
fn two_chord_segments(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let upper = profile(vec![
        arc_turns(p(0.0, 0.0), p(1.0, 0.0), 1),
        line(p(0.0, 1.0), p(1.0, 0.0)),
    ]);
    let lower = profile(vec![
        arc_turns(p(0.0, 0.0), p(0.0, -1.0), 1),
        line(p(1.0, 0.0), p(0.0, -1.0)),
    ]);
    (sweep(m, upper, 2.0), sweep(m, lower, 2.0))
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

/// `solid` with everything outside `z ∈ [lo, hi]` cut away by boxes wider than the unit disk — the
/// cut keeps the lateral's cylinder handle.
fn trimmed(m: &mut Model, mut solid: Handle<Solid>, lo: f64, hi: f64) -> Handle<Solid> {
    let below = crate::fixtures::cuboid(m, [-2.0, -2.0, lo - 2.0], [2.0, 2.0, lo]);
    solid = crate::boolean(m, BoolKind::Cut, solid, below).expect("the trim below")[0];
    let above = crate::fixtures::cuboid(m, [-2.0, -2.0, hi], [2.0, 2.0, hi + 2.0]);
    crate::boolean(m, BoolKind::Cut, solid, above).expect("the trim above")[0]
}

/// The upper half disk on `z ∈ [0, 2]` and on `z ∈ [1, 3]`: the two statements of each ruling of
/// the diameter overlap on `[1, 2]`.
fn halves_overlapping(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let low = sweep(m, upper_half(), 2.0);
    let high = sweep(m, upper_half(), 3.0);
    (low, trimmed(m, high, 1.0, 3.0))
}

/// The upper half disk on `z ∈ [0, 3]` and on `z ∈ [1, 2]`: one statement of each ruling lies
/// inside the other, both its ends inside.
fn half_around_half(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let tall = sweep(m, upper_half(), 3.0);
    let inner = sweep(m, upper_half(), 3.0);
    (tall, trimmed(m, inner, 1.0, 2.0))
}

/// A `4 × 4 × 2` block on `z ∈ [−4, −2]` bored through by the unit circle swept down from `z = 0`,
/// and the circle's pin on `z ∈ [−2, 0]` above it — [`bore_above_pin`] the other way up, so each
/// face meets the rim at `z = −2` from the other side.
fn bore_below_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, -4.0], [2.0, 2.0, -2.0]);
    let tool = sweep(m, disk(), -4.0);
    let bored = crate::boolean(m, BoolKind::Cut, block, tool).expect("the bore")[0];
    (bored, sweep(m, disk(), -2.0))
}

const FAMILIES: [(&str, Build); 17] = [
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
    ("a bored block above the bore's pin", bore_above_pin),
    ("the two halves of a disk", two_halves),
    ("a fillet swept both ways", fillet_swept_both_ways),
    ("a filleted plate and its fillet's pin", fillet_and_its_pin),
    ("two plates rounded on opposite quarters", opposite_quarters),
    ("two identical cylinders on a datum", identical_on_a_datum),
    ("two chord segments meeting on a ruling", two_chord_segments),
    ("two half disks on [0, 2] and [1, 3]", halves_overlapping),
    (
        "a half disk on [0, 3] around one on [1, 2]",
        half_around_half,
    ),
    ("a bored block below the bore's pin", bore_below_pin),
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

/// What a family's boolean comes out as: the volume of each output solid, ascending (none for an
/// empty result), and the solid's face count — or the refusal of a result that touches itself.
struct Want {
    volumes: Vec<f64>,
    /// Each body's face count, the bodies in the order of `volumes`.
    faces: Vec<usize>,
    refused: Option<RejectReason>,
}

/// The bodies' volumes and, body by body in the same order, their face counts.
fn counted(volumes: &[f64], faces: &[usize]) -> Want {
    assert_eq!(volumes.len(), faces.len(), "one face count per body");
    Want {
        volumes: volumes.to_vec(),
        faces: faces.to_vec(),
        refused: None,
    }
}

fn refused(reason: RejectReason) -> Want {
    Want {
        volumes: Vec::new(),
        faces: Vec::new(),
        refused: Some(reason),
    }
}

const PI: f64 = std::f64::consts::PI;

/// A unit disk's segment cut off by a chord a quarter turn long, swept 2: `(π/4 − 1/2)·2`.
const SEGMENT: f64 = PI / 2.0 - 1.0;

/// Each family's answers as Fuse, the first less the second, the second less the first, Common —
/// the volumes set out from the shapes (a fillet's quarter is `1 − π/4` short of its square's
/// corner).
fn wants(name: &str) -> [Want; 4] {
    let plate = 14.0 + PI / 2.0; // `filleted_right`, swept 2
    let none = || counted(&[], &[]);
    match name {
        "two identical prisms" | "two identical cylinders on a datum" => [
            counted(&[2.0 * PI], &[3]),
            none(),
            none(),
            counted(&[2.0 * PI], &[3]),
        ],
        "one circle swept both ways" => [
            counted(&[4.0 * PI], &[3]),
            counted(&[2.0 * PI], &[3]),
            counted(&[2.0 * PI], &[3]),
            none(),
        ],
        "a bored plate and the bore's pin" => [
            counted(&[32.0], &[6]),
            counted(&[32.0 - 2.0 * PI], &[7]),
            counted(&[2.0 * PI], &[3]),
            none(),
        ],
        "a short and a tall prism" => [
            counted(&[3.0 * PI], &[3]),
            none(),
            counted(&[PI], &[3]),
            counted(&[2.0 * PI], &[3]),
        ],
        "a half disk and the whole" => [
            counted(&[2.0 * PI], &[3]),
            none(),
            counted(&[PI], &[4]),
            counted(&[PI], &[4]),
        ],
        // The tall less the short is the lower half under a whole slab: one lateral face round
        // both, the flat, the overhang's half disk, the two caps.
        "a short half disk and a tall whole" => [
            counted(&[3.0 * PI], &[3]),
            none(),
            counted(&[2.0 * PI], &[5]),
            counted(&[PI], &[4]),
        ],
        // The notch and the pin meet on the two rulings alone: two bodies — the pin's four faces,
        // and the notched plate's eight (its bottom two pieces either side of the notch).
        "a notched plate and the half pin" => [
            counted(&[PI, 16.0 - PI], &[4, 8]),
            counted(&[16.0 - PI], &[8]),
            counted(&[PI], &[4]),
            none(),
        ],
        // The pin meets the block on the rim circle alone: two bodies — the pin's three faces, the
        // bored block's seven.
        "a bored block above the bore's pin" | "a bored block below the bore's pin" => [
            counted(&[2.0 * PI, 32.0 - 2.0 * PI], &[3, 7]),
            counted(&[32.0 - 2.0 * PI], &[7]),
            counted(&[2.0 * PI], &[3]),
            none(),
        ],
        // The two segments meet on the ruling `(1, 0)` alone: two bodies of four faces each.
        "two chord segments meeting on a ruling" => [
            counted(&[SEGMENT, SEGMENT], &[4, 4]),
            counted(&[SEGMENT], &[4]),
            counted(&[SEGMENT], &[4]),
            none(),
        ],
        "the two halves of a disk" => [
            counted(&[2.0 * PI], &[3]),
            counted(&[PI], &[4]),
            counted(&[PI], &[4]),
            none(),
        ],
        "a fillet swept both ways" => [
            counted(&[2.0 * plate], &[7]),
            counted(&[plate], &[7]),
            counted(&[plate], &[7]),
            none(),
        ],
        // The pin's circle also touches the plate's top wall, at `(3, 2)`: the plate less the pin
        // thins to nothing along that line.
        "a filleted plate and its fillet's pin" => [
            counted(&[plate + PI], &[8]),
            refused(RejectReason::SelfTouchingResult),
            counted(&[PI], &[3]),
            counted(&[2.0 * PI], &[3]),
        ],
        "two plates rounded on opposite quarters" => [
            counted(&[24.0], &[6]),
            counted(&[10.0 - PI / 2.0], &[7]),
            counted(&[10.0 - PI / 2.0], &[7]),
            counted(&[4.0 + PI], &[6]),
        ],
        // Each half disk is four faces — the lateral, the flat, two caps — and so is their union:
        // the two statements of a ruling fold where they overlap.
        "two half disks on [0, 2] and [1, 3]" => [
            counted(&[1.5 * PI], &[4]),
            counted(&[PI / 2.0], &[4]),
            counted(&[PI / 2.0], &[4]),
            counted(&[PI / 2.0], &[4]),
        ],
        // The tall less the inner is two half disks, below and above it.
        "a half disk on [0, 3] around one on [1, 2]" => [
            counted(&[1.5 * PI], &[4]),
            counted(&[PI / 2.0, PI / 2.0], &[4, 4]),
            none(),
            counted(&[PI / 2.0], &[4]),
        ],
        other => panic!("no answers stated for {other}"),
    }
}

fn face_count(m: &Model, s: Handle<Solid>) -> usize {
    m.shell(m.solid(s).outer).faces.len()
}

/// **`got` is `want`** — the refusal it names, or a valid model whose output solids have the
/// volumes and, body by body, the face counts stated, and whose mesh covers their faces.
fn check(label: &str, m: &Model, got: BoolResult, want: &Want) {
    if let Some(reason) = want.refused {
        assert!(
            matches!(&got, Err(BoolError::Rejected { reason: r, .. }) if *r == reason),
            "{label}: {got:?}"
        );
        return;
    }
    let got = got.unwrap_or_else(|e| panic!("{label}: {e:?}"));
    let bad = nacre_validate::validate(m);
    assert!(bad.is_empty(), "{label}: {bad:?}");
    let mut bodies: Vec<(f64, usize)> = got
        .iter()
        .map(|&s| {
            let v = nacre_props::mass_props(m, s).expect("props").volume;
            (v, face_count(m, s))
        })
        .collect();
    bodies.sort_by(|x, y| x.0.total_cmp(&y.0));
    let volumes: Vec<f64> = bodies.iter().map(|b| b.0).collect();
    assert_eq!(volumes.len(), want.volumes.len(), "{label}: {volumes:?}");
    for (v, w) in volumes.iter().zip(&want.volumes) {
        assert!((v - w).abs() < 1e-9, "{label}: {volumes:?}");
    }
    let faces: Vec<usize> = bodies.iter().map(|b| b.1).collect();
    assert_eq!(faces, want.faces, "{label}");
    mesh_covers_faces(label, m, &got);
}

/// **A cylinder two solids share is one class read for both.** Every family, as each of the four
/// booleans, comes out as [`check`] reads its answer.
#[test]
fn a_cylinder_shared_by_two_solids_reads_both_solids() {
    for (name, build) in FAMILIES {
        for (op, want) in OPS.into_iter().zip(wants(name)) {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}

fn extrude_in(m: &mut Model, frame: SketchFrame, profile: Profile2d, dist: f64) -> Handle<Solid> {
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

/// The frame of a datum plane through `origin` with normal `normal`.
fn datum(m: &mut Model, origin: [f64; 3], normal: [f64; 3]) -> SketchFrame {
    let plane = SketchPlane::from_origin_normal(
        nacre_math::Point3::from_array(origin),
        nacre_math::Vector3::from_array(normal),
    )
    .expect("a plane");
    let OpOutput::DatumPlane { frame, .. } = apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    )
    .expect("the datum") else {
        unreachable!()
    };
    frame
}

fn moved(m: &mut Model, solid: Handle<Solid>, isometry: nacre_exact::Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } =
        apply(m, &Operation::Transform { solid, isometry }).expect("the motion")
    else {
        unreachable!()
    };
    solid
}

fn lift(dz: i128) -> nacre_exact::Isometry {
    let z = nacre_exact::Rat::from_int;
    nacre_exact::Isometry::translation([z(0), z(0), z(dz)])
}

/// The unit circle swept `2` up from `z = 0`, and its twin swept `2` from the same plane stated
/// facing down — the same cylinder, its axis written `−ẑ`.
fn twin_with_the_axis_reversed(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 2.0);
    let down = datum(m, [0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
    (a, extrude_in(m, down, disk(), 2.0))
}

/// The twin sketched with its `+u` along `ŷ` — another seam direction.
fn twin_with_another_seam(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 2.0);
    let plane = m.seeded_plane(Axis::Z);
    let frame = SketchFrame::named(
        m,
        plane,
        nacre_math::Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .expect("the frame");
    (a, extrude_in(m, frame, disk(), 2.0))
}

/// A prism on `z ∈ [0, 4]` and a twin sketched on the datum `z = 2` — another origin on the line.
fn twin_on_another_origin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 4.0);
    let up = datum(m, [0.0, 0.0, 2.0], [0.0, 0.0, 1.0]);
    (a, extrude_in(m, up, disk(), 2.0))
}

/// A prism on `z ∈ [0, 4]` and a twin moved up `2` — another origin, carried by a motion.
fn twin_moved_up(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 4.0);
    let b = sweep(m, disk(), 2.0);
    (a, moved(m, b, lift(2)))
}

/// The twin turned a quarter about the axis — a motion that turns the seam.
fn twin_turned_a_quarter(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 2.0);
    let b = sweep(m, disk(), 2.0);
    (a, moved(m, b, super::rot_iso(Axis::Z, 90)))
}

/// The twin moved up `1` and back down.
fn twin_moved_there_and_back(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 2.0);
    let b = sweep(m, disk(), 2.0);
    let b = moved(m, b, lift(1));
    (a, moved(m, b, lift(-1)))
}

/// The twin turned a quarter and then three quarters about the axis — a whole turn.
fn twin_turned_a_whole_turn(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, disk(), 2.0);
    let b = sweep(m, disk(), 2.0);
    let b = moved(m, b, super::rot_iso(Axis::Z, 90));
    (a, moved(m, b, super::rot_iso(Axis::Z, 270)))
}

/// **One cylinder under two handles is refused by name** — the pair rule's `same_surface`
/// (`CylinderPairContact`; measured with `reject-trace`: every row here raises at the pair rule's
/// refusal). Each twin states the operand's cylinder once more, differing in one part of the
/// statement only, and each lands on a handle of its own.
#[test]
fn one_cylinder_under_two_handles_is_refused() {
    let twins: [(&str, Build); 5] = [
        ("the axis reversed", twin_with_the_axis_reversed),
        ("another seam", twin_with_another_seam),
        ("another origin", twin_on_another_origin),
        ("moved up", twin_moved_up),
        ("turned a quarter", twin_turned_a_quarter),
    ];
    for (name, build) in twins {
        for op in OPS {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
            assert_eq!((la.len(), lb.len()), (1, 1), "{name}");
            assert_ne!(la, lb, "{name}: two handles");
            let got = run(&mut m, op, a, b);
            assert!(
                matches!(
                    &got,
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

/// **A motion that folds back to none is one statement** — the twin moved there and back, or
/// turned a whole turn, interns to the operand's own handle, and the pair builds as two identical
/// prisms do.
#[test]
fn a_motion_folded_back_states_the_cylinder_once() {
    let twins: [(&str, Build); 2] = [
        ("moved there and back", twin_moved_there_and_back),
        ("turned a whole turn", twin_turned_a_whole_turn),
    ];
    for (name, build) in twins {
        for (op, want) in OPS.into_iter().zip(wants("two identical prisms")) {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            assert_eq!(
                lateral_surfaces(&m, a),
                lateral_surfaces(&m, b),
                "{name}: one handle"
            );
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}
