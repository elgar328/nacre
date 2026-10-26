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

/// The Fuse refused, the rest built: the pin's ray to the notched plate's top face passes inside
/// the notch circle's box but outside the half arc's.
fn fuse_refused(plate: f64) -> [[Want; 3]; 2] {
    let mut w = apart(plate);
    w[0][0] = Want::Refused(RejectReason::CylinderGateUndecided);
    w[1][0] = Want::Refused(RejectReason::CylinderGateUndecided);
    w
}

/// **A pin under the notch, inside the notch circle's box** — the pin's centre `(41.5, 20)` lies
/// inside the notch circle's box `[37, 43] × [17, 23]` and outside the half arc that bounds the
/// plate, which lies in `x ≤ 40`. The Fuse's ray from the pin asks whether it misses the plate's
/// top face, whose class is nameless; the face's box takes the arc's whole circle and cannot show
/// the miss. Four constructions: the pin below, the pin along `−z` above, the box fused before the
/// notch, the box cut instead.
#[test]
fn a_pin_under_the_notch_inside_its_circle_box_is_refused_in_the_fuse() {
    every_kind(|m| pin_below(m, 41.5), fuse_refused(NOTCHED + 125.0));
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
        fuse_refused(NOTCHED + 125.0),
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
        fuse_refused(NOTCHED + 125.0),
    );
    every_kind(
        |m| {
            let p = notched(m);
            let t = turned_box(m, 37, [10.0, 10.0], 5.0, 10.0);
            let a = op1(m, BoolKind::Cut, p, t);
            (a, zcyl(m, [41.5, 20.0, -4.0], 1.0, 2.0, [0.0, -1.0, 0.0]))
        },
        fuse_refused(NOTCHED - 125.0),
    );
}

/// The controls: the pin's centre outside the circle's box, just (`x = 43.5`) and well (`x = 46`).
#[test]
fn a_pin_outside_the_notch_circle_box_builds() {
    every_kind(|m| pin_below(m, 46.0), apart(NOTCHED + 125.0));
    every_kind(|m| pin_below(m, 43.5), apart(NOTCHED + 125.0));
}

// ---- a nameless class's line off the lateral, read on the carriers ----

/// **A turned box beside the fillet** — the box's walls have no world statement and their planes
/// section the fillet cylinder while their faces stay apart from it. The gate asks whether the line
/// such a class meets another recorded class in lies off the lateral (`line_off_lateral`), and it
/// asks the two carriers: the line is an interval box, and the angular test has no interval
/// spelling. Four constructions: two parallel walls, a second box on a fused plate, the slotted
/// plate, a pocketed plate with an `8°` box.
#[test]
fn a_turned_box_beside_the_fillet_is_refused_by_the_gate() {
    let undecided = || refused(RejectReason::CylinderGateUndecided);
    every_kind(
        |m| {
            let p = plate(m);
            (p, turned_box(m, 37, [31.0, -8.0], 6.0, 10.0))
        },
        undecided(),
    );
    every_kind(
        |m| {
            let p = plate(m);
            let t = turned_box(m, 37, [10.0, 16.0], 5.0, 10.0);
            let a = op1(m, BoolKind::Fuse, p, t);
            (a, turned_box(m, 37, [5.0, 19.0], 5.0, 10.0))
        },
        undecided(),
    );
    every_kind(
        |m| {
            let p = slotted_plate(m);
            let b = bx(m, [-1.5, -2.5, 5.0], [1.5, 2.5, 15.0]);
            let b = turn(m, b, 37, [0.0; 3]);
            (p, shift(m, b, [18.0, 30.0, 0.0]))
        },
        undecided(),
    );
    every_kind(
        |m| {
            let p = plate(m);
            let k = bx(m, [43.0, -14.0, 6.0], [46.0, -4.0, 13.0]);
            let a = op1(m, BoolKind::Cut, p, k);
            let t = bx(m, [38.0, -6.0, 5.0], [43.0, -1.0, 15.0]);
            (a, turn(m, t, 8, [43.0, -16.0, 0.0]))
        },
        undecided(),
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

/// **A tangent line beside a third class's plane, off its faces** — the slab's face `y = 1` is
/// tangent to the unit cylinder along the line `x = 0, y = 1`; the notch's wall (or the arm's end)
/// `x = 0` holds that line as a plane, but its face stands at `y ≥ 2`, off the line. The gate asks
/// which classes hold the tangent line by their infinite planes (`classes_holding_the_line`) and
/// refuses the line in another plane. Four constructions: the notch, the arm, the block on the
/// cylinder's side, the notch turned and moved.
#[test]
fn a_tangent_line_beside_a_third_plane_off_its_faces_is_refused() {
    let named = || refused(RejectReason::TangentLineInAnotherPlane);
    every_kind(|m| notched_slab(m, 0.0), named());
    every_kind(|m| ell_slab(m, 0.0), named());
    every_kind(|m| slab_and_blocked_cylinder(m, 0.0), named());
    every_kind(
        |m| {
            let (a, b) = notched_slab(m, 0.0);
            let a = turn(m, a, 90, [0.0; 3]);
            let b = turn(m, b, 90, [0.0; 3]);
            let a = shift(m, a, [0.5, 0.25, 0.0]);
            let b = shift(m, b, [0.5, 0.25, 0.0]);
            (a, b)
        },
        named(),
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

/// **A pocket floor beside the fillet's quarter** — the pocket's wall `y = −17` crosses the fillet
/// cylinder, and its face clears the fillet's quarter but not the whole circle the wall test reads
/// (`face_clears_footprint` asks the footprint's span only). The record it writes meets the turned
/// box's walls, which have no world statement, and the gate's line test asks the carriers
/// (`line_off_lateral`): refused.
#[test]
fn a_pocket_wall_beside_the_fillet_quarter_is_refused_by_the_gate() {
    every_kind(
        |m| pocket_beside_the_fillet(m, -17.0),
        refused(RejectReason::CylinderGateUndecided),
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
