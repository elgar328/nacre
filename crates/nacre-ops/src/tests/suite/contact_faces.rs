//! **Faces that meet their own body only on contact lines** — the scoreboard of booleans whose
//! result has a face, or a group of faces, every edge of whose outline is a line or an arc where
//! two bodies touch: four faces use it, two of each body. Each family is run as Fuse, both Cuts and
//! Common:
//!
//! - a block bored square through and the bore's pin below it, touching on the bore's outline
//!   (two plane classes on each contact line);
//! - the same with a half-disk bore and pin (the contact is the half disk's diameter and its arc);
//! - the square pin with a blind pocket in its top (what meets its body only on contact lines is
//!   the top with the pocket's walls and floor, not one face);
//! - a half disk prism whose flat fills the window of a slab (the contact runs along two of the
//!   cylinder's rulings, where the window's side walls are tangent to it, and along two chords);
//! - a cube and a triangular prism touching along one edge line that three plane classes hold.

use super::scoreboard::{OPS, Want, check, counted, refused, run};
use super::*;

type Build = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);

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

/// The unit circle's half on `y ≥ 0`, closed by its diameter.
fn upper_half() -> Profile2d {
    profile(vec![
        line(p(-1.0, 0.0), p(1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(1.0, 0.0), 2),
    ])
}

fn cut(m: &mut Model, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
    crate::boolean(m, BoolKind::Cut, a, b).expect("the cut")[0]
}

/// A `4 × 4 × 2` block on `z ∈ [2, 4]` bored through by the square `[−1, 1]²`, and the square's
/// pin on `z ∈ [0, 2]` below it.
fn square_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, 2.0], [2.0, 2.0, 4.0]);
    let tool = crate::fixtures::cuboid(m, [-1.0, -1.0, 0.0], [1.0, 1.0, 4.0]);
    let bored = cut(m, block, tool);
    (
        bored,
        crate::fixtures::cuboid(m, [-1.0, -1.0, 0.0], [1.0, 1.0, 2.0]),
    )
}

/// The block bored through by the unit circle's upper half, and the half's pin below it — one
/// cylinder handle.
fn d_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, 2.0], [2.0, 2.0, 4.0]);
    let tool = sweep(m, upper_half(), 4.0);
    let bored = cut(m, block, tool);
    (bored, sweep(m, upper_half(), 2.0))
}

/// [`square_pin`] with a `1 × 1 × 1` pocket in the pin's top.
fn pocketed_pin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let (bored, pin) = square_pin(m);
    let pocket = crate::fixtures::cuboid(m, [-0.5, -0.5, 1.0], [0.5, 0.5, 2.0]);
    (bored, cut(m, pin, pocket))
}

/// A slab on `y ∈ [−1, 0]` with the window `[−1, 1] × [0, 2]` (in `x`, `z`) through it, and the
/// unit circle's upper half swept `z ∈ [0, 2]`, whose flat fills the window.
fn half_disk_in_a_window(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let slab = crate::fixtures::cuboid(m, [-3.0, -1.0, -1.0], [3.0, 0.0, 3.0]);
    let window = crate::fixtures::cuboid(m, [-1.0, -1.0, 0.0], [1.0, 0.0, 2.0]);
    let frame = cut(m, slab, window);
    (frame, sweep(m, upper_half(), 2.0))
}

/// The unit cube and the triangle `(1, 1)`, `(3, 3)`, `(1, 3)` swept `z ∈ [0, 1]`: they touch along
/// the line `x = y = 1`, which the planes `x = 1`, `y = 1` and `x = y` hold.
fn three_planes_on_a_line(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let cube = crate::fixtures::cuboid(m, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let triangle = profile(vec![
        line(p(1.0, 1.0), p(3.0, 3.0)),
        line(p(3.0, 3.0), p(1.0, 3.0)),
        line(p(1.0, 3.0), p(1.0, 1.0)),
    ]);
    (cube, sweep(m, triangle, 1.0))
}

const FAMILIES: [(&str, Build); 5] = [
    ("a square pin under its bore", square_pin),
    ("a half-disk pin under its bore", d_pin),
    ("a pocketed square pin under its bore", pocketed_pin),
    ("a half disk in a window", half_disk_in_a_window),
    (
        "two prisms on a line of three planes",
        three_planes_on_a_line,
    ),
];

const PI: f64 = std::f64::consts::PI;

/// Each family's answers as Fuse, the first less the second, the second less the first, Common.
fn wants(name: &str) -> [Want; 4] {
    let none = || counted(&[], &[]);
    match name {
        // The Fuse is the two bodies touching on the bore's outline — refused, todo
        // 「곧은 모서리·호에서만 만나는 면」.
        "a square pin under its bore" => [
            refused(RejectReason::NonManifoldResultEdge),
            counted(&[24.0], &[10]),
            counted(&[8.0], &[6]),
            none(),
        ],
        "a half-disk pin under its bore" => [
            refused(RejectReason::OpenResultShell),
            counted(&[32.0 - PI], &[8]),
            counted(&[PI], &[4]),
            none(),
        ],
        "a pocketed square pin under its bore" => [
            refused(RejectReason::StraightAngle),
            counted(&[24.0], &[10]),
            counted(&[7.0], &[11]),
            none(),
        ],
        "a half disk in a window" => [
            refused(RejectReason::NonManifoldResultEdge),
            counted(&[20.0], &[10]),
            counted(&[PI], &[4]),
            none(),
        ],
        "two prisms on a line of three planes" => [
            counted(&[1.0, 2.0], &[6, 5]),
            counted(&[1.0], &[6]),
            counted(&[2.0], &[5]),
            none(),
        ],
        other => panic!("no answers stated for {other}"),
    }
}

/// Every family, as each of the four booleans, comes out as [`check`] reads its answer.
#[test]
fn faces_meeting_their_body_on_contact_lines() {
    for (name, build) in FAMILIES {
        for (op, want) in OPS.into_iter().zip(wants(name)) {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}
