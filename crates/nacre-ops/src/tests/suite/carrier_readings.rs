//! **Bodies whose pieces are apart while their carriers meet** — each beside a control that
//! differs in that one element. A site that asks the infinite carrier (a plane, a line, a whole
//! circle) where it holds the bounded piece refuses the first and builds the second; a site that
//! asks the piece builds both.
//!
//! The plate is `rounded_plate(1, 0)`: `90 × 50 × 12` with one `r = 5` fillet about `(40, −20)`,
//! the quarter facing `(+x, −y)`. A line through `(43, −16)` along `z` lies on that fillet's
//! cylinder (`(3, 4)` from the axis) but off its quarter.
use super::*;
use nacre_exact::{Angle, Isometry, Rat, Rotation};
use std::f64::consts::PI;

const PLATE: f64 = 53_700.0 + 75.0 * PI;

fn q(x: f64) -> Rat {
    Rat::from_decimal(x).expect("a short decimal")
}

fn turn(m: &mut Model, s: Handle<Solid>, deg: i128, pivot: [f64; 3]) -> Handle<Solid> {
    let out = transform(
        m,
        s,
        &Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: pivot.map(q),
            angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
        }),
    )
    .expect("the solid turns");
    m.rebuild_adjacency();
    out
}

fn shift(m: &mut Model, s: Handle<Solid>, v: [f64; 3]) -> Handle<Solid> {
    let out = transform(m, s, &Isometry::translation(v.map(q))).expect("the solid moves");
    m.rebuild_adjacency();
    out
}

fn bx(m: &mut Model, min: [f64; 3], max: [f64; 3]) -> Handle<Solid> {
    let s = crate::fixtures::cuboid(m, Point3::from_array(min), Point3::from_array(max));
    m.rebuild_adjacency();
    s
}

fn zcyl(m: &mut Model, at: [f64; 3], r: f64, h: f64, seam: [f64; 3]) -> Handle<Solid> {
    let s = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array(at),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array(seam),
        r,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    s
}

fn op1(m: &mut Model, kind: BoolKind, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
    let out = boolean(m, kind, a, b).expect("the fixture's boolean builds");
    assert_eq!(out.len(), 1, "fixture boolean gives one solid");
    m.rebuild_adjacency();
    out[0]
}

/// A `5 × 5 × h` box turned `deg` about its own centre, then set at `c`.
fn turned_box(m: &mut Model, deg: i128, c: [f64; 2], z0: f64, h: f64) -> Handle<Solid> {
    let b = bx(m, [-2.5, -2.5, z0], [2.5, 2.5, z0 + h]);
    let b = turn(m, b, deg, [0.0; 3]);
    shift(m, b, [c[0], c[1], 0.0])
}

fn plate(m: &mut Model) -> Handle<Solid> {
    let p = rounded_plate(m, 1, 0);
    m.rebuild_adjacency();
    p
}

/// `40 × 40 × 10` with a slot `[15, 25] × [17, 23]` rounded at both ends (`r = 3`).
fn slotted_plate(m: &mut Model) -> Handle<Solid> {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let edges = vec![
        line(p2(0.0, 0.0), p2(40.0, 0.0)),
        line(p2(40.0, 0.0), p2(40.0, 40.0)),
        line(p2(40.0, 40.0), p2(0.0, 40.0)),
        line(p2(0.0, 40.0), p2(0.0, 0.0)),
        line(p2(15.0, 17.0), p2(25.0, 17.0)),
        arc_turns(p2(25.0, 20.0), p2(25.0, 17.0), 2),
        line(p2(25.0, 23.0), p2(15.0, 23.0)),
        arc_turns(p2(15.0, 20.0), p2(15.0, 23.0), 2),
    ];
    let profile = stated(edges).unwrap().remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 10.0,
        },
    )
    .expect("the slotted plate extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

const SLOTTED: f64 = 16_000.0 - 600.0 - 90.0 * PI;

/// The cells of [`every_kind`] for operands of volume `va`, `vb` meeting in volume `i`, every
/// result one body.
fn one_body(va: f64, vb: f64, i: f64) -> [[Want; 3]; 2] {
    [
        [
            Want::Built(1, va + vb - i),
            Want::Built(1, va - i),
            Want::Built(1, i),
        ],
        [
            Want::Built(1, va + vb - i),
            Want::Built(1, vb - i),
            Want::Built(1, i),
        ],
    ]
}

fn refused(name: RejectReason) -> [[Want; 3]; 2] {
    [[Want::Refused(name); 3]; 2]
}

// ---- a segment on a cylinder, off its faces ----

/// An L of two boxes on the plate's filleted corner: `[41, 50] × [y0, −12]` and
/// `[43, 50] × [−13, −4]`, both over `z`.
fn ell(m: &mut Model, y0: f64, z: [f64; 2]) -> (Handle<Solid>, Handle<Solid>) {
    let p = plate(m);
    let b1 = bx(m, [41.0, y0, z[0]], [50.0, -12.0, z[1]]);
    let b2 = bx(m, [43.0, -13.0, z[0]], [50.0, -4.0, z[1]]);
    let b = op1(m, BoolKind::Fuse, b1, b2);
    (p, b)
}

/// A box on the corner with a notch cut from its inner side.
fn notched_corner_box(m: &mut Model, y0: f64) -> (Handle<Solid>, Handle<Solid>) {
    let p = plate(m);
    let b = bx(m, [41.0, y0, 2.0], [50.0, -4.0, 10.0]);
    let k = bx(m, [40.0, -12.0, 1.0], [43.0, -3.0, 11.0]);
    let b = op1(m, BoolKind::Cut, b, k);
    (p, b)
}

/// **A box edge on the fillet's cylinder, beyond the quarter** — the L's edge at `(43, −16)`
/// runs along `z` on the fillet cylinder's surface, and the class `x = 43` (the arm's wall, or the
/// notch's) carries the segment. The fillet face is the quarter facing `(+x, −y)`; the edge is at `(3, 4)` from the
/// axis, off it, so the segment crosses no edge of the face. The ruling split asks the segment's
/// angle of the cylinder's lateral faces (`pieces_clear`), not its line of the cylinder, and builds.
#[test]
fn an_edge_on_the_fillet_cylinder_beyond_its_quarter_builds() {
    every_kind(
        |m| ell(m, -16.0, [2.0, 10.0]),
        one_body(PLATE, 736.0, 256.0),
    );
    every_kind(
        |m| ell(m, -16.0, [-2.0, 14.0]),
        one_body(PLATE, 1472.0, 384.0),
    );
    every_kind(
        |m| notched_corner_box(m, -16.0),
        one_body(PLATE, 736.0, 256.0),
    );
}

/// The controls: the same bodies with the edge at `y = −15.5`, off the cylinder.
#[test]
fn an_edge_off_the_fillet_cylinder_builds() {
    every_kind(
        |m| ell(m, -15.5, [2.0, 10.0]),
        one_body(PLATE, 700.0, 240.0),
    );
    every_kind(
        |m| ell(m, -15.5, [-2.0, 14.0]),
        one_body(PLATE, 1400.0, 360.0),
    );
    every_kind(
        |m| notched_corner_box(m, -15.5),
        one_body(PLATE, 700.0, 240.0),
    );
}

// ---- crossings of a circle off its arcs ----

/// **A box corner on the fillet's circle, beyond the quarter** — the box `[43, 50] × [−16, −10]`
/// puts its corner `(43, −16)` on the fillet's cap circles, off the quarter arc, where three names
/// meet on each cap. The circle split asks each crossing against the arcs' ends before it orders
/// any, so the three never meet in an order, and builds.
#[test]
fn a_corner_on_the_fillet_circle_beyond_its_quarter_builds() {
    every_kind(
        |m| {
            let p = plate(m);
            (p, bx(m, [43.0, -16.0, 2.0], [50.0, -10.0, 10.0]))
        },
        one_body(PLATE, 336.0, 96.0),
    );
}

/// The control: the corner at `y = −15.5`, off the circle.
#[test]
fn a_corner_off_the_fillet_circle_builds() {
    every_kind(
        |m| {
            let p = plate(m);
            (p, bx(m, [43.0, -15.5, 2.0], [50.0, -10.0, 10.0]))
        },
        one_body(PLATE, 308.0, 88.0),
    );
}

// ---- a ray beside an arc, inside its circle's box ----

/// `40 × 40 × 10` with a half-round notch (`r = 3`) cut from the edge `x = 40` at `y = 20`.
fn notched(m: &mut Model) -> Handle<Solid> {
    let p = bx(m, [0.0, 0.0, 0.0], [40.0, 40.0, 10.0]);
    let c = zcyl(m, [40.0, 20.0, -1.0], 3.0, 12.0, [1.0, 0.0, 0.0]);
    op1(m, BoolKind::Cut, p, c)
}

const NOTCHED: f64 = 16_000.0 - 45.0 * PI;
const PIN: f64 = 2.0 * PI;

/// A pin (`r = 1`, `z ∈ [−4, −2]`) below the notched plate at `(x, 20)`; the plate carries a
/// `37°` box, so its top face's class has no world statement.
fn pin_below(m: &mut Model, x: f64) -> (Handle<Solid>, Handle<Solid>) {
    let p = notched(m);
    let t = turned_box(m, 37, [10.0, 10.0], 5.0, 10.0);
    let a = op1(m, BoolKind::Fuse, p, t);
    (a, zcyl(m, [x, 20.0, -4.0], 1.0, 2.0, [0.0, -1.0, 0.0]))
}

/// What the pin and the plate give, apart: `plate` the plate's volume.
fn apart(plate: f64) -> [[Want; 3]; 2] {
    [
        [
            Want::Built(2, plate + PIN),
            Want::Built(1, plate),
            Want::Built(0, 0.0),
        ],
        [
            Want::Built(2, plate + PIN),
            Want::Built(1, PIN),
            Want::Built(0, 0.0),
        ],
    ]
}

/// **A pin under the notch, inside the notch circle's box** — the pin's centre `(41.5, 20)` lies
/// inside the notch circle's box `[37, 43] × [17, 23]` and outside the half arc that bounds the
/// plate, which lies in `x ≤ 40`. The Fuse's ray from the pin asks whether it misses the plate's
/// top face, whose class is nameless; the face's box takes the arc's own box (`arc_box_rat`), which
/// the ray misses, and builds. Four constructions: the pin below, the pin along `−z` above, the box
/// fused before the notch, the box cut instead.
#[test]
fn a_pin_under_the_notch_inside_its_circle_box_builds() {
    every_kind(|m| pin_below(m, 41.5), apart(NOTCHED + 125.0));
    every_kind(
        |m| {
            let p = notched(m);
            let t = turned_box(m, 37, [10.0, 10.0], 5.0, 10.0);
            let a = op1(m, BoolKind::Fuse, p, t);
            let c = crate::fixtures::cylinder(
                m,
                Point3::from_array([41.5, 20.0, 14.0]),
                Vector3::from_array([0.0, 0.0, -1.0]),
                1.0,
                2.0,
            )
            .solid;
            m.rebuild_adjacency();
            (a, c)
        },
        apart(NOTCHED + 125.0),
    );
    every_kind(
        |m| {
            let p = bx(m, [0.0, 0.0, 0.0], [40.0, 40.0, 10.0]);
            let t = turned_box(m, 37, [10.0, 10.0], 5.0, 10.0);
            let p = op1(m, BoolKind::Fuse, p, t);
            let c = zcyl(m, [40.0, 20.0, -1.0], 3.0, 12.0, [1.0, 0.0, 0.0]);
            let a = op1(m, BoolKind::Cut, p, c);
            (a, zcyl(m, [41.5, 20.0, -4.0], 1.0, 2.0, [0.0, -1.0, 0.0]))
        },
        apart(NOTCHED + 125.0),
    );
    every_kind(
        |m| {
            let p = notched(m);
            let t = turned_box(m, 37, [10.0, 10.0], 5.0, 10.0);
            let a = op1(m, BoolKind::Cut, p, t);
            (a, zcyl(m, [41.5, 20.0, -4.0], 1.0, 2.0, [0.0, -1.0, 0.0]))
        },
        apart(NOTCHED - 125.0),
    );
}

/// `40 × 40 × 10` whose edge `x = 40` bulges out in a half round (`r = 3` about `(40, 20)`),
/// carrying a `37°` box centred on `c`, `z ∈ [5, 15]`.
fn bulged_with_box(m: &mut Model, c: [f64; 2]) -> Handle<Solid> {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let edges = vec![
        line(p2(0.0, 0.0), p2(40.0, 0.0)),
        line(p2(40.0, 0.0), p2(40.0, 17.0)),
        arc_turns(p2(40.0, 20.0), p2(40.0, 17.0), 2),
        line(p2(40.0, 23.0), p2(40.0, 40.0)),
        line(p2(40.0, 40.0), p2(0.0, 40.0)),
        line(p2(0.0, 40.0), p2(0.0, 0.0)),
    ];
    let profile = stated(edges).unwrap().remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 10.0,
        },
    )
    .expect("the bulged plate extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let t = turned_box(m, 37, c, 5.0, 10.0);
    op1(m, BoolKind::Fuse, solid, t)
}

/// **A pin under the bulge, the box inside the plate's outline.** The box leaves a hole in the top
/// face whose corners its nameless walls name; the face's plane and outer ring are stated. The
/// Fuse's ray up the pin's axis crosses the top face at a rational point, read exactly against the
/// outer ring (the bulge, inside) and by its box against the hole (outside): the ray up counts the
/// top and the bottom, and the pin is apart. Read as «the hole holds the crossing», the top would go
/// uncounted and the pin would be inside the plate — the lock on that box reading.
#[test]
fn a_pin_under_a_bulging_arc_beside_a_nameless_hole_builds() {
    every_kind(
        |m| {
            let a = bulged_with_box(m, [10.0, 10.0]);
            (a, zcyl(m, [41.5, 20.0, -4.0], 1.0, 2.0, [0.0, -1.0, 0.0]))
        },
        apart(16_000.0 + 45.0 * PI + 125.0),
    );
}

/// **A pin under the bulge — the arc's own box holds it.** The box straddles the plate's edge
/// `y = 40`, so the top face's outer ring itself carries corners its nameless walls name as well
/// as the half round, and is read by its box alone. The half round bulges past the plate's straight
/// corners, so the crossing of the Fuse's ray up the pin's axis is inside the arc's box: it cannot
/// be shown outside — rightly, the ray does cross the face — and the ray down answers, the pin apart
/// from the plate. This is the lock on the arc's direction: read the other way round, the box is the
/// complementary half's, which lies inside the corners' hull, the crossing would be shown outside a
/// face it is in, the ray up would count the bottom face alone and put the pin inside the plate. The
/// lock rides on the ray up being asked before the ray down (`combinatorics::probe_dirs`'s order).
///
/// The box's centre is on the edge, and a centrally symmetric square is halved by any line through
/// its centre: half its `z ∈ [5, 10]` part is in the plate, so it adds `125 + 62.5`.
#[test]
fn a_pin_under_a_bulging_arc_is_inside_its_box_and_the_ray_down_answers() {
    every_kind(
        |m| {
            let a = bulged_with_box(m, [20.0, 40.0]);
            (a, zcyl(m, [41.5, 20.0, -4.0], 1.0, 2.0, [0.0, -1.0, 0.0]))
        },
        apart(16_000.0 + 45.0 * PI + 187.5),
    );
}

/// The controls: the pin's centre outside the circle's box, just (`x = 43.5`) and well (`x = 46`).
#[test]
fn a_pin_outside_the_notch_circle_box_builds() {
    every_kind(|m| pin_below(m, 46.0), apart(NOTCHED + 125.0));
    every_kind(|m| pin_below(m, 43.5), apart(NOTCHED + 125.0));
}

// ---- a nameless class's line off the lateral, read on the faces ----

/// **A turned box beside the fillet** — the box's walls have no world statement and their planes
/// section the fillet cylinder while their faces stay apart from it. The gate asks whether the line
/// such a class meets another recorded class in lies on a lateral face (`line_off_lateral`), and
/// asks it of the faces. Four constructions: two parallel walls of one box, which share no line;
/// three parallel walls of two boxes on a fused plate; two parallel walls beside the slotted
/// plate's slot end; and a pocketed plate with an `8°` box, whose wall meets the pocket's wall
/// `x = 43` in a line on the fillet's cylinder at `(3, 4)` from the axis, off its quarter.
#[test]
fn a_turned_box_beside_the_fillet_builds() {
    every_kind(
        |m| {
            let p = plate(m);
            (p, turned_box(m, 37, [31.0, -8.0], 6.0, 10.0))
        },
        one_body(PLATE, 250.0, 150.0),
    );
    every_kind(
        |m| {
            let p = plate(m);
            let t = turned_box(m, 37, [10.0, 16.0], 5.0, 10.0);
            let a = op1(m, BoolKind::Fuse, p, t);
            (a, turned_box(m, 37, [5.0, 19.0], 5.0, 10.0))
        },
        one_body(PLATE + 75.0, 250.0, 175.0),
    );
    every_kind(
        |m| {
            let p = slotted_plate(m);
            let b = bx(m, [-1.5, -2.5, 5.0], [1.5, 2.5, 15.0]);
            let b = turn(m, b, 37, [0.0; 3]);
            (p, shift(m, b, [18.0, 30.0, 0.0]))
        },
        one_body(SLOTTED, 150.0, 75.0),
    );
    every_kind(
        |m| {
            let p = plate(m);
            let k = bx(m, [43.0, -14.0, 6.0], [46.0, -4.0, 13.0]);
            let a = op1(m, BoolKind::Cut, p, k);
            let t = bx(m, [38.0, -6.0, 5.0], [43.0, -1.0, 15.0]);
            (a, turn(m, t, 8, [43.0, -16.0, 0.0]))
        },
        one_body(PLATE - 120.0, 250.0, 175.0),
    );
}

/// The plate and two boxes fused and turned `8°` about `(43, y0)`: `[23, 43] × [y0, y0 + 20]` over
/// `z ∈ [13, 20]`, above the plate, its corner at the pivot; and `[28, 33] × [y0 + 10, y0 + 15]`
/// over `z ∈ [5, 14]`, reaching into the plate far from the fillet.
fn corner_above_the_fillet(m: &mut Model, y0: f64) -> (Handle<Solid>, Handle<Solid>) {
    let p = plate(m);
    let pivot = [43.0, y0, 0.0];
    let t = bx(m, [23.0, y0, 13.0], [43.0, y0 + 20.0, 20.0]);
    let t = turn(m, t, 8, pivot);
    let k = bx(m, [28.0, y0 + 10.0, 5.0], [33.0, y0 + 15.0, 14.0]);
    let k = turn(m, k, 8, pivot);
    let b = op1(m, BoolKind::Fuse, t, k);
    (p, b)
}

/// The plate and a wedge above it, turned `−15°` about `p`: the triangle `p`, `p + (0, 20)`,
/// `p + (−10, 20)` over `z ∈ [13, 20]`, fused with a leg `[px − 8, px − 5] × [py + 18, py + 20]` over
/// `z ∈ [5, 14]` that reaches into the plate far from the fillet — `700 + 54 − 6 = 748`, of which the
/// leg's `42` lies in the plate. The two walls through `p` meet at about `26.6°`. Seen from
/// `(44, −17)` the quarter's box spans about `72°`, less than a right angle, so two perpendicular
/// walls through that point cannot both cross it; the wedge's two do, and both pass within the
/// radius (about `3.1` and `4.5` from the axis), so each sections the fillet cylinder. The third
/// wall and the leg's walls miss the quarter's box.
fn wedge_above_the_fillet(m: &mut Model, p: [f64; 2]) -> (Handle<Solid>, Handle<Solid>) {
    let plate = plate(m);
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let tip = p2(p[0], p[1]);
    let near = p2(p[0], p[1] + 20.0);
    let far = p2(p[0] - 10.0, p[1] + 20.0);
    let profile = stated(vec![line(tip, near), line(near, far), line(far, tip)])
        .unwrap()
        .remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 7.0,
        },
    )
    .expect("the wedge extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let w = shift(m, solid, [0.0, 0.0, 13.0]);
    let leg = bx(
        m,
        [p[0] - 8.0, p[1] + 18.0, 5.0],
        [p[0] - 5.0, p[1] + 20.0, 14.0],
    );
    let b = op1(m, BoolKind::Fuse, w, leg);
    (plate, turn(m, b, -15, [p[0], p[1], 0.0]))
}

/// What the class gate answers for the operands `build` makes: the name it refuses by, or what it
/// returned instead.
fn gate_refusal(
    build: impl Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>),
) -> Result<RejectReason, String> {
    let mut m = Model::new();
    let (a, b) = build(&mut m);
    m.rebuild_adjacency();
    match crate::arrangement::plane_index_setup(&m, a, b) {
        Err(BoolError::Rejected { reason, .. }) => Ok(reason),
        other => Err(format!("{:?}", other.map(|_| ()))),
    }
}

/// **Two nameless walls meeting on the fillet's quarter** — two walls through a corner on the
/// quarter, both sectioning the fillet cylinder, meet in a line on the lateral face, though their
/// faces stand above it. No record states that line, and the gate refuses the pair. It is held at
/// the gate: past it, the split refuses the same shape by the same name. Two constructions: the
/// upper box's corner at `(43, −24)` (`(3, −4)` from the axis) and the wedge's tip at `(44, −23)`
/// (`(4, −3)`).
#[test]
fn two_nameless_walls_meeting_on_the_quarter_are_refused_at_the_gate() {
    let bx_on = |m: &mut Model| corner_above_the_fillet(m, -24.0);
    let wedge_on = |m: &mut Model| wedge_above_the_fillet(m, [44.0, -23.0]);
    every_kind(bx_on, refused(RejectReason::CylinderGateUndecided));
    every_kind(wedge_on, refused(RejectReason::CylinderGateUndecided));
    assert_eq!(
        [gate_refusal(bx_on), gate_refusal(wedge_on)],
        [
            Ok(RejectReason::CylinderGateUndecided),
            Ok(RejectReason::CylinderGateUndecided)
        ],
        "the gate refuses both"
    );
}

/// The controls: each body with its corner on the fillet's cylinder, off its quarter. The boxes `8`
/// along `y`, the corner at `(43, −16)` (`(3, 4)` from the axis): one wall through that corner
/// stands beside the quarter's box, so no line of two sectioning walls is asked here — what this
/// holds is that the refusal above is the corner's place on the quarter. The wedge's tip at
/// `(44, −17)` (`(4, 3)`): both walls through it section the cylinder, and the gate shows their
/// line off the lateral face by its direction against the face's arc (`line_off_lateral`, the
/// interval body of `line_off_faces`).
#[test]
fn two_nameless_walls_meeting_off_the_quarter_build() {
    every_kind(
        |m| corner_above_the_fillet(m, -16.0),
        one_body(PLATE, 3000.0, 175.0),
    );
    every_kind(
        |m| wedge_above_the_fillet(m, [44.0, -17.0]),
        one_body(PLATE, 748.0, 42.0),
    );
}

/// The controls: the same boxes away from the fillet.
#[test]
fn a_turned_box_away_from_the_fillet_builds() {
    let one = |m: &mut Model| {
        let p = plate(m);
        (p, turned_box(m, 37, [10.0, 16.0], 5.0, 10.0))
    };
    every_kind(one, one_body(PLATE, 250.0, 175.0));
    every_kind(
        |m| {
            let (p, t) = one(m);
            let a = op1(m, BoolKind::Fuse, p, t);
            (a, turned_box(m, 37, [3.0, 18.0], 5.0, 10.0))
        },
        one_body(PLATE + 75.0, 250.0, 175.0),
    );
    every_kind(
        |m| {
            let p = slotted_plate(m);
            let b = bx(m, [-1.5, -2.5, 5.0], [1.5, 2.5, 15.0]);
            let b = turn(m, b, 37, [0.0; 3]);
            (p, shift(m, b, [17.0, 30.0, 0.0]))
        },
        one_body(SLOTTED, 150.0, 75.0),
    );
    every_kind(
        |m| {
            let p = plate(m);
            let k = bx(m, [43.0, -14.0, 6.0], [46.0, -4.0, 13.0]);
            let a = op1(m, BoolKind::Cut, p, k);
            let t = bx(m, [38.0, -5.0, 5.0], [43.0, 0.0, 15.0]);
            (a, turn(m, t, 8, [43.0, -15.0, 0.0]))
        },
        one_body(PLATE - 120.0, 250.0, 175.0),
    );
    every_kind(
        |m| {
            let (p, t) = one(m);
            let a = op1(m, BoolKind::Fuse, p, t);
            (a, bx(m, [-20.0, -5.0, -2.0], [-10.0, 5.0, 14.0]))
        },
        {
            let mut w = one_body(PLATE + 75.0, 1600.0, 1200.0);
            w[1][1] = Want::Built(2, 400.0);
            w
        },
    );
}

// ---- a tangent line beside a third class's plane, off its faces ----

/// `r = 1` along `z` over `[0, 2]`, its seam at `−y`.
fn unit_cylinder(m: &mut Model) -> Handle<Solid> {
    zcyl(m, [0.0, 0.0, 0.0], 1.0, 2.0, [0.0, -1.0, 0.0])
}

/// A slab `[−2, 2] × [1, 6]` whose face `y = 1` is tangent to the unit cylinder along `x = 0`, with
/// a notch `[x, 3] × [4, 7]` cut from its far side.
fn notched_slab(m: &mut Model, x: f64) -> (Handle<Solid>, Handle<Solid>) {
    let a = bx(m, [-2.0, 1.0, 0.5], [2.0, 6.0, 1.5]);
    let n = bx(m, [x, 4.0, 0.0], [3.0, 7.0, 2.0]);
    let a = op1(m, BoolKind::Cut, a, n);
    (a, unit_cylinder(m))
}

/// A slab tangent along `x = 0` and an arm `[−2, x] × [2, 6]` fused to it.
fn ell_slab(m: &mut Model, x: f64) -> (Handle<Solid>, Handle<Solid>) {
    let a = bx(m, [-2.0, 1.0, 0.5], [2.0, 3.0, 1.5]);
    let b = bx(m, [-2.0, 2.0, 0.5], [x, 6.0, 1.5]);
    let a = op1(m, BoolKind::Fuse, a, b);
    (a, unit_cylinder(m))
}

/// The tangent slab against the cylinder fused with a block `[x, 2] × [−3, −0.5]` on its far side.
fn slab_and_blocked_cylinder(m: &mut Model, x: f64) -> (Handle<Solid>, Handle<Solid>) {
    let a = bx(m, [-2.0, 1.0, 0.5], [2.0, 3.0, 1.5]);
    let c = unit_cylinder(m);
    let t = bx(m, [x, -3.0, 0.0], [2.0, -0.5, 2.0]);
    let b = op1(m, BoolKind::Fuse, c, t);
    (a, b)
}

/// What a body and the unit cylinder (or what it carries) give, touching along a line: `va`, `vb`.
fn touching(va: f64, vb: f64) -> [[Want; 3]; 2] {
    [
        [
            Want::Built(2, va + vb),
            Want::Built(1, va),
            Want::Built(0, 0.0),
        ],
        [
            Want::Built(2, va + vb),
            Want::Built(1, vb),
            Want::Built(0, 0.0),
        ],
    ]
}

/// **A tangent line beside a third class's plane, off its faces, builds** — the slab's face `y = 1`
/// is tangent to the unit cylinder along the line `x = 0, y = 1`; the notch's wall (or the arm's
/// end, or the block's side) `x = 0` holds that line as a plane, but its face stands at `y ≥ 3` (or
/// `y ≤ −1`), off the line. The gate reads the faces on the classes holding the line
/// (`line_runs_through_face` on that plane), finds none touching it, and hands the row to the
/// verdict, which reads three regions there; the rulings split carries the slab face's segment on
/// `x = 0` onto that class's ruling. Four constructions — the notch, the arm, the block on the
/// cylinder's side, the notch turned and moved — each the same answers as its control below with
/// the third plane at `x = 0.5`.
#[test]
fn a_tangent_line_beside_a_third_plane_off_its_faces_builds() {
    every_kind(|m| notched_slab(m, 0.0), touching(16.0, 2.0 * PI));
    every_kind(|m| ell_slab(m, 0.0), touching(14.0, 2.0 * PI));
    // The block meets the cylinder over the part of the disk with `x ≥ 0, y ≤ −0.5`:
    // `∫ (√(1 − x²) − 1/2) dx` over `[0, √3/2]` = `π/6 − √3/8`, over the height 2.
    let lens = PI / 6.0 - 3f64.sqrt() / 8.0;
    every_kind(
        |m| slab_and_blocked_cylinder(m, 0.0),
        touching(8.0, 2.0 * PI + 10.0 - 2.0 * lens),
    );
    every_kind(
        |m| {
            let (a, b) = notched_slab(m, 0.0);
            let a = turn(m, a, 90, [0.0; 3]);
            let b = turn(m, b, 90, [0.0; 3]);
            let a = shift(m, a, [0.5, 0.25, 0.0]);
            let b = shift(m, b, [0.5, 0.25, 0.0]);
            (a, b)
        },
        touching(16.0, 2.0 * PI),
    );
}

/// The controls: the third plane at `x = 0.5`, off the tangent line, and the slab alone.
#[test]
fn a_tangent_line_with_no_third_plane_builds() {
    every_kind(
        |m| {
            let a = bx(m, [-2.0, 1.0, 0.5], [2.0, 3.0, 1.5]);
            (a, unit_cylinder(m))
        },
        touching(8.0, 2.0 * PI),
    );
    every_kind(|m| notched_slab(m, 0.5), touching(17.0, 2.0 * PI));
    every_kind(|m| ell_slab(m, 0.5), touching(15.5, 2.0 * PI));
    // The block meets the cylinder over the part of the disk with `x ≥ 0.5, y ≤ −0.5`:
    // `∫ (√(1 − x²) − 1/2) dx` over `[1/2, √3/2]` = `π/12 − √3/4 + 1/4`, over the height 2.
    let lens = PI / 12.0 - 3f64.sqrt() / 4.0 + 0.25;
    every_kind(
        |m| slab_and_blocked_cylinder(m, 0.5),
        touching(8.0, 2.0 * PI + 7.5 - 2.0 * lens),
    );
}

/// The tangent slab with a wedge on top whose face on `x = 0` — a plane holding the tangent line —
/// is a triangle with one corner on the line, at `(0, 1, 1.75)` on the lateral, the face running off
/// it above the slab: `(y, z)` = `(1, 1.75)`, `(2.5, 1.2)`, `(2.5, 2.5)`, swept over `x ∈ [0, 2]`.
fn slab_with_a_wedge_on_the_line(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let slab = bx(m, [-2.0, 1.0, 0.5], [2.0, 3.0, 1.5]);
    let pts =
        [(1.0, 1.75), (2.5, 1.2), (2.5, 2.5)].map(|(u, v)| nacre_math::Point2::from_array([u, v]));
    let ring = crate::Ring2d::polygon_decimal(pts.to_vec()).expect("a triangle");
    let profile = crate::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = SketchFrame::world(m, Axis::X);
    let OpOutput::Extrude { solid: wedge, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    )
    .expect("the wedge extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let a = op1(m, BoolKind::Fuse, slab, wedge);
    (a, unit_cylinder(m))
}

/// **A third plane's face meeting the tangent line at a corner is refused as oblique.** The
/// wedge's face on `x = 0` meets the line `x = 0, y = 1` at one corner on the lateral; the corner's
/// other edges leave it neither along the axis nor level, so the faces across them lie on planes
/// slanted to the axis through that point of the lateral, and the gate refuses those before any
/// tangency row is read — every operation, either order. (How the face itself reads against the
/// line is the reader's own test, `a_face_meeting_the_line_at_a_corner_does_not_touch_it`.)
#[test]
fn a_third_plane_face_meeting_the_tangent_line_at_a_corner_is_refused_as_oblique() {
    every_kind(
        slab_with_a_wedge_on_the_line,
        refused(RejectReason::ObliqueCylinderCut),
    );
}

// ---- a face clear of a fillet's quarter, on a plane through its circle ----

/// The plate with a pocket `[41, 50] × [y0, −10] × [6, 13]` at its filleted corner, and an `8°`
/// box beside it.
fn pocket_beside_the_fillet(m: &mut Model, y0: f64) -> (Handle<Solid>, Handle<Solid>) {
    let p = plate(m);
    let k = bx(m, [41.0, y0, 6.0], [50.0, -10.0, 13.0]);
    let a = op1(m, BoolKind::Cut, p, k);
    let t = bx(m, [38.0, -7.0, 5.0], [44.0, -2.0, 15.0]);
    (a, turn(m, t, 8, [44.0, -17.0, 0.0]))
}

/// **A pocket wall beside the fillet's quarter** — the pocket's wall `y = −17` crosses the fillet
/// cylinder off its quarter: the plane misses the lateral face though it cuts the infinite cylinder
/// (`laterals_clear_plane`, asked before the faces), so the pair writes nothing. A record there
/// would meet the turned box's walls, which have no world statement, in a line on the cylinder's
/// surface, which the gate would ask of the faces (`line_off_lateral`). The arrangement's split
/// reads the same placement under two rules, and this test locks both: an absence proved outranks
/// an undecided «on the segment» (`crossings::collect`), and a line along the axis is off the
/// lateral when its angle is off every face's arc (`line_off_faces`). Undoing either refuses the
/// Fuse with `CylinderGateUndecided`.
#[test]
fn a_pocket_wall_beside_the_fillet_quarter_builds() {
    every_kind(
        |m| pocket_beside_the_fillet(m, -17.0),
        one_body(PLATE - 168.0, 300.0, 210.0),
    );
}

/// The controls: the pocket wall at `y = −16.5` and at `y = −14`, off the fillet cylinder.
#[test]
fn a_pocket_wall_off_the_fillet_cylinder_builds() {
    every_kind(
        |m| pocket_beside_the_fillet(m, -16.5),
        one_body(PLATE - 156.0, 300.0, 210.0),
    );
    every_kind(
        |m| pocket_beside_the_fillet(m, -14.0),
        one_body(PLATE - 96.0, 300.0, 210.0),
    );
}
