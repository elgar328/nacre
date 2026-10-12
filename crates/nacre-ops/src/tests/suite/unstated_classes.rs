//! **Plane classes with no world statement beside world-stated cylinders** — the scoreboard.
//!
//! A box turned 37° about `z` has walls whose normals are irrational: no world name, so the
//! cylinder gate cannot state them in `Rat`. Each family here pairs such a class with a world
//! cylinder (a bore, a pin, a fillet) and pins two answers: what the gate says, read from
//! [`plane_index_setup`] directly — the reject census is process-global and owns no row in a
//! parallel suite — and what the boolean answers, by kind. Volumes are derived, not copied.
//!
//! The geometry, shared by the families: a plate `[0,40]²×[0,10]` (or the fillet plate of
//! [`rounded_plate`], `90×50×12`), and a `5×5` box turned about `z` around its own centre and
//! moved to `(cx, cy)`. The 37° turn is what makes a class nameless; 90° keeps every name (F0).

use super::*;
use nacre_exact::{Angle, Isometry, Rat, Rotation};

fn rat(x: i128) -> Rat {
    Rat::from_int(x)
}

fn turn(m: &mut Model, s: Handle<Solid>, axis: Axis, deg: i128, pivot: [i128; 3]) -> Handle<Solid> {
    let out = transform(
        m,
        s,
        &Isometry::rotation(Rotation {
            axis,
            pivot: pivot.map(rat),
            angle: Angle::from_deg(rat(deg)).expect("an angle"),
        }),
    )
    .expect("the solid turns");
    m.rebuild_adjacency();
    out
}

fn shift(m: &mut Model, s: Handle<Solid>, v: [i128; 3]) -> Handle<Solid> {
    let out = transform(m, s, &Isometry::translation(v.map(rat))).expect("the solid moves");
    m.rebuild_adjacency();
    out
}

fn cuboid(m: &mut Model, min: [f64; 3], max: [f64; 3]) -> Handle<Solid> {
    let s = crate::fixtures::cuboid(m, Point3::from_array(min), Point3::from_array(max));
    m.rebuild_adjacency();
    s
}

/// A `z` cylinder of radius `r` standing on `(x, y, z0)` with height `h`.
fn z_cylinder(m: &mut Model, at: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
    let s = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array(at),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        r,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    s
}

fn op(m: &mut Model, kind: BoolKind, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
    let out = boolean(m, kind, a, b).expect("the fixture's boolean builds")[0];
    m.rebuild_adjacency();
    out
}

/// The plate `[0,40]²×[0,10]` with a through bore of radius 2 at `(10, 10)` — volume `16000 − 40π`.
fn plate_with_bore(m: &mut Model) -> Handle<Solid> {
    let plate = cuboid(m, [0.0; 3], [40.0, 40.0, 10.0]);
    let bore = z_cylinder(m, [10.0, 10.0, -1.0], 2.0, 12.0);
    op(m, BoolKind::Cut, plate, bore)
}

/// The bored plate with the wall `[9,11]×[0,40]×[5,15]` fused on — its plane `x = 9` runs within
/// the bore's radius, and its face across the bore.
fn walled_bore(m: &mut Model) -> Handle<Solid> {
    let p = plate_with_bore(m);
    let wall = cuboid(m, [9.0, 0.0, 5.0], [11.0, 40.0, 15.0]);
    op(m, BoolKind::Fuse, p, wall)
}

/// The plate with a pin of radius `r`, height 3, standing on its top at `(x, y)` — volume
/// `16000 + 3πr²`.
fn plate_with_pin(m: &mut Model, x: f64, y: f64, r: f64) -> Handle<Solid> {
    let plate = cuboid(m, [0.0; 3], [40.0, 40.0, 10.0]);
    let pin = z_cylinder(m, [x, y, 10.0], r, 3.0);
    op(m, BoolKind::Fuse, plate, pin)
}

/// A `5×5` box over `z0..z0+h`, turned `deg` about `z` around its own centre, centred on
/// `(cx, cy)`.
fn turned_box(m: &mut Model, deg: i128, c: [i128; 2], z0: f64, h: f64) -> Handle<Solid> {
    let b = cuboid(m, [-2.5, -2.5, z0], [2.5, 2.5, z0 + h]);
    let b = turn(m, b, Axis::Z, deg, [0, 0, 0]);
    shift(m, b, [c[0], c[1], 0])
}

/// The plate `[0,40]²×[0,10]` with a slot through it: the stadium of half-width 3 about the segment
/// `y = 20`, `x ∈ [15, 25]` — two half-cylinder ends whose lateral faces have no whole rim, joined
/// by walls tangent to them.
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

/// [`slotted_plate`]'s volume: `(1600 − (10·6 + 9π))·10`.
fn slotted_plate_volume() -> f64 {
    (1600.0 - (60.0 + 9.0 * std::f64::consts::PI)) * 10.0
}

/// The fillet plate's top area: `90×50` less the corners its fillets of radius 5 round off, less
/// its bores of radius 3.5. The plate is 12 thick; a `5×5` box standing on `z = 6` overlaps it by
/// `25·6 = 150`.
fn rounded_area(fillets: usize, bores: usize) -> f64 {
    let pi = std::f64::consts::PI;
    90.0 * 50.0 - fillets as f64 * (25.0 - 25.0 * pi / 4.0) - bores as f64 * pi * 3.5 * 3.5
}

/// What one kind answers.
#[derive(Clone, Copy, Debug)]
enum Want {
    /// This many solids of this total volume, `validate` clean.
    Built(usize, f64),
    /// Refused by this name.
    Refused(RejectReason),
}
use Want::{Built, Refused};

struct Family {
    name: &'static str,
    build: fn(&mut Model) -> (Handle<Solid>, Handle<Solid>),
    /// The gate passes, with this many (cylinder, nameless class) pairs silent by each arm —
    /// `(beside, apart, sectioned)`: the plane clears the cylinder; only the class's faces stand
    /// apart from a cylinder with whole rims; or the plane cuts the lateral in rulings with every
    /// face apart, and the class draws them on its own arrangement — or refuses. The counts are derived from the geometry, family by family, in the comments;
    /// a family both arms could open is pinned to the arm that does, so either arm breaking shows.
    gate: Result<(usize, usize, usize), RejectReason>,
    /// Fuse, Cut, Common.
    want: [Want; 3],
}

#[test]
fn nameless_classes_beside_world_cylinders() {
    use RejectReason::CylinderGateUndecided as Undecided;
    let pi = std::f64::consts::PI;
    let bored = 16000.0 - 40.0 * pi;
    // The bored plate fused with the wall `[9,11]×[0,40]×[5,15]`: the wall's upper half stands
    // clear, and its lower half is inside the plate except where it fills the bore's strip
    // `|x − 10| ≤ 1` — `∫ 2√(4 − u²) du` over `[−1, 1]` = `2√3 + 4π/3`.
    let walled = bored + 800.0 - (400.0 - 5.0 * (2.0 * 3f64.sqrt() + 4.0 * pi / 3.0));
    // The turned box standing on `z0`: what it shares with the walled plate — `25` per unit of
    // depth inside the plate (`z` up to 10), and above it the box's section inside the wall's strip
    // `9 ≤ x ≤ 11`. The box's width across a line `x = 10 + u` is `5/cos θ` while the line passes
    // between its two near corners (`|u| ≤ a = 2.5(cos θ − sin θ)`) and falls linearly to `0` at
    // the far ones (`b = 2.5(cos θ + sin θ)`).
    let sunk = |z0: f64| {
        let (c, s) = (37f64.to_radians().cos(), 37f64.to_radians().sin());
        let (a, b) = (2.5 * (c - s), 2.5 * (c + s));
        let strip = 2.0 * (5.0 / c) * (a + ((b - a).powi(2) - (b - 1.0).powi(2)) / (2.0 * (b - a)));
        25.0 * (10.0 - z0) + strip * (15.0 - 10.0)
    };
    let refused = [Refused(Undecided); 3];
    let fams = [
        // The control: a quarter turn keeps every name, and the box sits on the plate far from
        // the bore — it overlaps the plate over `z ∈ [5, 10]`, 125.
        Family {
            name: "F0 box turned 90",
            build: |m| (plate_with_bore(m), turned_box(m, 90, [30, 30], 5.0, 10.0)),
            gate: Ok((0, 0, 0)),
            want: [
                Built(1, bored + 125.0),
                Built(1, bored - 125.0),
                Built(1, 125.0),
            ],
        },
        // Far from the bore, and every wall plane clears it as well.
        Family {
            name: "F1d box turned 37, every wall plane clear of the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 10], 5.0, 10.0)),
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, bored + 125.0),
                Built(1, bored - 125.0),
                Built(1, 125.0),
            ],
        },
        // Far from the bore, but one wall's plane, extended, runs through it (1.43 from the axis).
        Family {
            name: "F1 box turned 37, a wall plane through the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 30], 5.0, 10.0)),
            gate: Ok((3, 1, 0)),
            want: [
                Built(1, bored + 125.0),
                Built(1, bored - 125.0),
                Built(1, 125.0),
            ],
        },
        // The same box lifted off the plate: two bodies under Fuse, the plate under Cut.
        Family {
            name: "F1b box turned 37, above the plate",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 30], 20.0, 5.0)),
            gate: Ok((3, 1, 0)),
            want: [Built(2, bored + 125.0), Built(1, bored), Built(0, 0.0)],
        },
        // A cube turned 37° about `x`, centred on the plate's top: by its point symmetry half of
        // it, 62.5, lies in the plate. Its walls are oblique to the bore and clear it.
        Family {
            name: "F10 cube turned 37 about x",
            build: |m| {
                let c = cuboid(m, [-2.5; 3], [2.5; 3]);
                let c = turn(m, c, Axis::X, 37, [0, 0, 0]);
                (plate_with_bore(m), shift(m, c, [30, 30, 10]))
            },
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, bored + 62.5),
                Built(1, bored - 62.5),
                Built(1, 62.5),
            ],
        },
        // A pin of radius 2 on the plate; one wall plane of the far box cuts it.
        Family {
            name: "F7 box turned 37 beside a pin",
            build: |m| {
                (
                    plate_with_pin(m, 10.0, 10.0, 2.0),
                    turned_box(m, 37, [30, 30], 5.0, 10.0),
                )
            },
            gate: Ok((3, 1, 0)),
            want: [
                Built(1, 16000.0 + 12.0 * pi + 125.0),
                Built(1, 16000.0 + 12.0 * pi - 125.0),
                Built(1, 125.0),
            ],
        },
        // The box encloses a pin of radius 1: Common is `25·2 + 3π`, which a road that guessed
        // "outside" for the pin's disk answered as 50.
        Family {
            name: "F1e box turned 37 around a pin",
            build: |m| {
                (
                    plate_with_pin(m, 30.0, 30.0, 1.0),
                    turned_box(m, 37, [30, 30], 8.0, 10.0),
                )
            },
            gate: Ok((4, 0, 0)),
            want: refused,
        },
        // The box encloses the bore: each wall stands 2.5 from the axis, beyond the radius 2, which
        // the gate shows from the axis's distance to the plane once the box test cannot (the box of
        // the bore's lateral straddles every wall) — four beside. What refuses is the nesting, as
        // around the pin (F1e).
        Family {
            name: "F1c box turned 37 across the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [10, 10], 5.0, 10.0)),
            gate: Ok((4, 0, 0)),
            want: refused,
        },
        // A lone cylinder and a box above its top, apart.
        Family {
            name: "F11a lone cylinder, box turned 37 above, planes clear",
            build: |m| {
                (
                    z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0),
                    turned_box(m, 37, [30, 10], 20.0, 5.0),
                )
            },
            gate: Ok((4, 0, 0)),
            want: [
                Built(2, 40.0 * pi + 125.0),
                Built(1, 40.0 * pi),
                Built(0, 0.0),
            ],
        },
        Family {
            name: "F11b lone cylinder, box turned 37 above, a wall plane through it",
            build: |m| {
                (
                    z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0),
                    turned_box(m, 37, [30, 30], 20.0, 5.0),
                )
            },
            gate: Ok((3, 1, 0)),
            want: [
                Built(2, 40.0 * pi + 125.0),
                Built(1, 40.0 * pi),
                Built(0, 0.0),
            ],
        },
        // The fillet plate with one fillet, the box in the open: its rings mix pierce corners
        // and three-plane corners, which is what sees a side recorded with the wrong sign.
        Family {
            name: "F8b one-fillet plate, box turned 37 clear of the fillet",
            build: |m| {
                (
                    rounded_plate(m, 1, 0),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, rounded_area(1, 0) * 12.0 + 100.0),
                Built(1, rounded_area(1, 0) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // The playground's example 003 shape: four fillets, four bores, the box in the middle.
        Family {
            name: "F16 four fillets and four bores, box turned 37 in the middle",
            build: |m| (rounded_plate(m, 4, 4), turned_box(m, 37, [0, 0], 6.0, 10.0)),
            gate: Ok((32, 0, 0)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // Four bores and no fillet.
        Family {
            name: "V4 four bores, box turned 37",
            build: |m| {
                (
                    rounded_plate(m, 0, 4),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((14, 2, 0)),
            want: [
                Built(1, rounded_area(0, 4) * 12.0 + 100.0),
                Built(1, rounded_area(0, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // A wall plane of the box cuts a fillet, whose lateral has ruling edges: the wall's class
        // sections the fillet (its plane runs along the axis within the radius, every face of it
        // apart) and closes the plate's section on its own arrangement; the other three walls clear
        // it. The box stands on `z = 5`, so it overlaps the 12-thick plate by `25·7 = 175`.
        Family {
            name: "F8 one-fillet plate, a wall plane through the fillet",
            build: |m| {
                (
                    rounded_plate(m, 1, 0),
                    turned_box(m, 37, [10, 16], 5.0, 10.0),
                )
            },
            gate: Ok((3, 0, 1)),
            want: [
                Built(1, rounded_area(1, 0) * 12.0 + 75.0),
                Built(1, rounded_area(1, 0) * 12.0 - 175.0),
                Built(1, 175.0),
            ],
        },
        // A slot's half-cylinder end is concave from the material, so a plane can pass through its
        // lateral box beyond the arc and still cut the plate. One wall does (the axis ≈ 3.47 from its
        // plane, beyond the radius 3): the box test cannot say, the axis's distance can — beside,
        // its side asked of an axis point where every side is minted. The opposite wall cuts the
        // same end (≈ 1.53): sectioned. The other six pairs clear by their boxes. The box stands on
        // `z = 5`, overlapping the plate by `25·5 = 125`.
        Family {
            name: "S1 slotted plate, one wall beside a slot end by its distance, one sectioning it",
            build: |m| (slotted_plate(m), turned_box(m, 37, [33, 11], 5.0, 10.0)),
            gate: Ok((7, 0, 1)),
            want: [
                Built(1, slotted_plate_volume() + 125.0),
                Built(1, slotted_plate_volume() - 125.0),
                Built(1, 125.0),
            ],
        },
        // Two wall planes, each through one fillet: two sectioned pairs of the sixteen.
        Family {
            name: "V1 four fillets, a wall plane through a fillet",
            build: |m| {
                (
                    rounded_plate(m, 4, 0),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((14, 0, 2)),
            want: [
                Built(1, rounded_area(4, 0) * 12.0 + 100.0),
                Built(1, rounded_area(4, 0) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // The playground's example 003 shape (four fillets, four bores) with the box off the middle,
        // where wall planes run through fillets. Eight cylinders, four walls: thirty-two pairs,
        // counted the gate's way — a pair whose lateral box lies on one side of the wall, or whose
        // axis lies beyond the radius from the wall's plane, is beside; else a bore (whole rims)
        // with the class's faces apart is apart; else a fillet whose axis lies within the radius
        // is sectioned. The bores whose box a wall plane straddles, with the axis's distance from
        // it against the radius 3.5: P1 `(−38, 18)` 3.82 and `(−38, −18)` 3.56, both beyond; P2
        // 3.98 and 3.93, both beyond; P3 `(−38, 18)` 3.22 and 1.78 (two walls), within, and
        // `(−38, −18)` 4.35, beyond; P4 1.04 and 1.04, within.
        Family {
            name: "P1 example-003 plate, box off the middle, a wall through one fillet",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-12, -6], 6.0, 10.0),
                )
            },
            gate: Ok((31, 0, 1)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "P2 example-003 plate, walls through two fillets",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-22, -14], 6.0, 10.0),
                )
            },
            gate: Ok((30, 0, 2)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "P3 example-003 plate, walls through two fillets and three bores' boxes",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-16, -10], 6.0, 10.0),
                )
            },
            gate: Ok((28, 2, 2)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "P4 example-003 plate, walls through two fillets and two bores' boxes",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((28, 2, 2)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // A pin standing over the plate's edge keeps its top rim whole (a disk cap) while its bottom
        // rim is cut at the edge. A wall plane through it, every face apart: the class would cut the
        // disk in a chord, which is drawn only for a recorded pair, so the gate does not section a
        // cylinder with a whole rim — refused as before.
        Family {
            name: "D1 a pin over the plate's edge, a wall plane through it",
            build: |m| {
                (
                    plate_with_pin(m, 39.0, 20.0, 2.0),
                    turned_box(m, 37, [30, 28], 5.0, 10.0),
                )
            },
            gate: Err(Undecided),
            want: refused,
        },
        // A wall plane that holds the fillet's tangent ruling `(40, −25, z)`: the box's wall
        // `x = 40` turned 37° about that ruling still holds it, runs along the axis within the
        // radius (the axis is `5·sin 37° ≈ 3.0` from it), and every face of the class stands far
        // from the plate — sectioned; the other three walls clear the fillet's box. The plate's own
        // wall `y = −25` ends on that ruling, but a cylinder's own solid writes no record of its
        // tangent walls, so the gate has no other record to meet. Downstream the split pass must
        // order the class's pierce point on that ruling against the plate's corner there — one
        // point, two names, which intervals cannot order — and refuses by the gate's name.
        Family {
            name: "F8t one-fillet plate, a wall plane holding the fillet's tangent ruling",
            build: |m| {
                let plate = rounded_plate(m, 1, 0);
                let b = cuboid(m, [35.0, -80.0, 0.0], [40.0, -75.0, 10.0]);
                (plate, turn(m, b, Axis::Z, 37, [40, -25, 0]))
            },
            gate: Ok((3, 0, 1)),
            want: refused,
        },
        // A sectioned class meeting another record of the fillet on one of its rulings. The other
        // operand is one solid: the box `[41,60]×[−40,−24]×[2,8]`, whose walls `x = 41` and
        // `y = −24` cut the fillet (crossings; `y = −24` meets it at `(43, −24)`), fused with a box
        // whose wall `x = 43` is turned 37° about that ruling, overlapping the first far from the
        // plate. The turned wall's class runs along the axis within the radius (`≈ 0.011` from it)
        // with every face apart — a sectioned candidate — and it meets `x = 41` off the lateral
        // (at `y ≈ −21.35`), but meets `y = −24` on the ruling `(43, −24)` itself (`C = 0`): a line
        // on the lateral no record states, which the gate refuses. The three other turned walls
        // clear the fillet's box. With that check planted away the gate passes `(3, 0, 1)` and the
        // three operations refuse downstream as `WitnessNotRational`: the check names the refusal,
        // it does not stand between the scoreboard and a wrong volume.
        Family {
            name: "F8u a sectioned class meets a recorded wall on a fillet ruling",
            build: |m| {
                let plate = rounded_plate(m, 1, 0);
                let a0 = cuboid(m, [41.0, -40.0, 2.0], [60.0, -24.0, 8.0]);
                let t = cuboid(m, [38.0, -46.0, 0.0], [43.0, -41.0, 10.0]);
                let t = turn(m, t, Axis::Z, 37, [43, -24, 0]);
                let b = boolean(m, BoolKind::Fuse, a0, t).expect("the two boxes fuse");
                assert_eq!(b.len(), 1, "the two boxes overlap into one solid");
                (plate, b[0])
            },
            gate: Err(Undecided),
            want: refused,
        },
        // The operands the other way round: the box first. `B − A` is now the box less the plate,
        // `250 − 150`.
        Family {
            name: "P2r box first, walls through two fillets",
            build: |m| {
                let plate = rounded_plate(m, 4, 4);
                (turned_box(m, 37, [-22, -14], 6.0, 10.0), plate)
            },
            gate: Ok((30, 0, 2)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, 100.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "F8r box first, a wall plane through the fillet",
            build: |m| {
                let plate = rounded_plate(m, 1, 0);
                (turned_box(m, 37, [10, 16], 5.0, 10.0), plate)
            },
            gate: Ok((3, 0, 1)),
            want: [
                Built(1, rounded_area(1, 0) * 12.0 + 75.0),
                Built(1, 75.0),
                Built(1, 175.0),
            ],
        },
        // A rib on the plate whose wall `x = 9` runs within the bore's radius (recorded), and a
        // turned box straddling that wall's plane far from the bore. On the wall's class the
        // bore's pierce points and the box's corners share lines and rings. The plate with the
        // rib is `bored + 200 − 100`; the box, `y` 21.5..28.5, misses the rib and the bore and
        // sinks `25·5` into the plate. ★ Its rib face covers none of the cells the order between
        // those points shapes, so this family cannot see that order's sign — F14c and F14d can.
        Family {
            name: "F14 recorded rib wall, box turned 37 across its plane",
            build: |m| {
                let p = plate_with_bore(m);
                let rib = cuboid(m, [9.0, 30.0, 5.0], [11.0, 40.0, 15.0]);
                let a = op(m, BoolKind::Fuse, p, rib);
                (a, turned_box(m, 37, [10, 25], 5.0, 10.0))
            },
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, bored + 225.0),
                Built(1, bored - 25.0),
                Built(1, 125.0),
            ],
        },
        // The wall runs the plate's length, across the bore, so its face holds both the bore's
        // pierce points and the box's crossings — the order between them is what shapes it.
        Family {
            name: "F14c wall across the bore, box turned 37 across its plane",
            build: |m| (walled_bore(m), turned_box(m, 37, [10, 25], 5.0, 10.0)),
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, walled + 250.0 - sunk(5.0)),
                Built(1, walled - sunk(5.0)),
                Built(1, sunk(5.0)),
            ],
        },
        // The same, operands swapped: the class roots and their frames are the box's now. The box
        // less the walled plate is two solids — the wall's strip runs across the box's top half.
        Family {
            name: "F14c' the same, operands swapped",
            build: |m| {
                let b = turned_box(m, 37, [10, 25], 5.0, 10.0);
                (b, walled_bore(m))
            },
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, walled + 250.0 - sunk(5.0)),
                Built(2, 250.0 - sunk(5.0)),
                Built(1, sunk(5.0)),
            ],
        },
        // The box raised to stand on `z = 8`, so less of it sinks into the plate.
        Family {
            name: "F14d wall across the bore, box turned 37 standing higher",
            build: |m| (walled_bore(m), turned_box(m, 37, [10, 25], 8.0, 10.0)),
            gate: Ok((4, 0, 0)),
            want: [
                Built(1, walled + 250.0 - sunk(8.0)),
                Built(1, walled - sunk(8.0)),
                Built(1, sunk(8.0)),
            ],
        },
        // A world wall crossing the bore, fused with a turned box: the wall's face carries the
        // box's irrational corners.
        Family {
            name: "F12 named wall across the bore with irrational corners",
            build: |m| {
                let w = cuboid(m, [9.0, 0.0, 5.0], [11.0, 40.0, 15.0]);
                let r = turned_box(m, 37, [9, 32], 5.0, 10.0);
                let tool = op(m, BoolKind::Fuse, w, r);
                (plate_with_bore(m), tool)
            },
            gate: Err(Undecided),
            want: refused,
        },
        // Out of scope: the bored plate itself turned, so the cylinder has no world statement.
        Family {
            name: "F2 bored plate turned 37",
            build: |m| {
                let p = plate_with_bore(m);
                let p = turn(m, p, Axis::Z, 37, [20, 20, 0]);
                (p, turned_box(m, 0, [27, 27], 5.0, 10.0))
            },
            gate: Err(Undecided),
            want: refused,
        },
    ];
    for f in &fams {
        {
            let mut m = Model::new();
            let (a, b) = (f.build)(&mut m);
            let gate = plane_index_setup(&m, a, b)
                .map(|setup| {
                    let mut by_arm = (0, 0, 0);
                    for cy in &setup.cyls {
                        for &(_, s) in &cy.silent {
                            match s {
                                crate::planes::Silence::Beside(_) => by_arm.0 += 1,
                                crate::planes::Silence::Apart => by_arm.1 += 1,
                                crate::planes::Silence::Sectioned => by_arm.2 += 1,
                            }
                        }
                    }
                    by_arm
                })
                .map_err(|e| match e {
                    BoolError::Rejected { reason, .. } => reason,
                    other => panic!(
                        "{}: the gate failed outside its vocabulary: {other:?}",
                        f.name
                    ),
                });
            assert_eq!(gate, f.gate, "{}: the gate", f.name);
        }
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(&f.want)
        {
            let mut m = Model::new();
            let (a, b) = (f.build)(&mut m);
            match (boolean(&mut m, kind, a, b), want) {
                (Ok(out), Built(n, v)) => {
                    assert_eq!(out.len(), *n, "{} {kind:?}: solids", f.name);
                    m.rebuild_adjacency();
                    let issues = nacre_validate::validate(&m);
                    assert!(issues.is_empty(), "{} {kind:?}: {issues:?}", f.name);
                    let got: f64 = out
                        .iter()
                        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                        .sum();
                    assert!(
                        (got - v).abs() <= 1e-9 * v.max(1.0),
                        "{} {kind:?}: volume {got} against {v}",
                        f.name
                    );
                }
                (Err(BoolError::Rejected { reason, .. }), Refused(r)) => {
                    assert_eq!(&reason, r, "{} {kind:?}: the refusal's name", f.name);
                }
                (got, _) => panic!("{} {kind:?}: {got:?}, wanted {want:?}", f.name),
            }
        }
    }
}
