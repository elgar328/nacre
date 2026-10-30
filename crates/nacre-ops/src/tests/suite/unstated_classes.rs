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
//! moved to `(cx, cy)`. The 37° turn is what makes a class nameless; 90° keeps every name (the first family).

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

use Want::{Built, Refused};

#[derive(Clone, Copy)]
struct Family {
    name: &'static str,
    build: fn(&mut Model) -> (Handle<Solid>, Handle<Solid>),
    /// The gate passes, with this many (cylinder, nameless class) pairs silent by each arm —
    /// `(beside, apart, sectioned, grazes)`: the plane clears the cylinder; only the class's faces
    /// stand apart from a cylinder with whole rims; the plane cuts the lateral in rulings with every
    /// face apart, and the class draws them on its own arrangement; or the plane touches the
    /// cylinder along one ruling with every face apart — or refuses. The counts are derived from the geometry, family by family, in the comments;
    /// a family both arms could open is pinned to the arm that does, so either arm breaking shows.
    gate: Result<(usize, usize, usize, usize), RejectReason>,
    /// Fuse, Cut, Common.
    want: [Want; 3],
}

#[test]
fn nameless_classes_beside_world_cylinders() {
    use RejectReason::CylinderGateUndecided as Undecided;
    let pi = std::f64::consts::PI;
    let bored = 16000.0 - 40.0 * pi;
    let walled = walled_bore_volume();
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
            name: "box turned 90",
            build: |m| (plate_with_bore(m), turned_box(m, 90, [30, 30], 5.0, 10.0)),
            gate: Ok((0, 0, 0, 0)),
            want: [
                Built(1, bored + 125.0),
                Built(1, bored - 125.0),
                Built(1, 125.0),
            ],
        },
        // Far from the bore, and every wall plane clears it as well.
        Family {
            name: "box turned 37, every wall plane clear of the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 10], 5.0, 10.0)),
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, bored + 125.0),
                Built(1, bored - 125.0),
                Built(1, 125.0),
            ],
        },
        // Far from the bore, but one wall's plane, extended, runs through it (1.43 from the axis).
        Family {
            name: "box turned 37, a wall plane through the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 30], 5.0, 10.0)),
            gate: Ok((3, 1, 0, 0)),
            want: [
                Built(1, bored + 125.0),
                Built(1, bored - 125.0),
                Built(1, 125.0),
            ],
        },
        // The same box lifted off the plate: two bodies under Fuse, the plate under Cut.
        Family {
            name: "box turned 37, above the plate",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 30], 20.0, 5.0)),
            gate: Ok((3, 1, 0, 0)),
            want: [Built(2, bored + 125.0), Built(1, bored), Built(0, 0.0)],
        },
        // A cube turned 37° about `x`, centred on the plate's top: by its point symmetry half of
        // it, 62.5, lies in the plate. Its walls are oblique to the bore and clear it.
        Family {
            name: "cube turned 37 about x",
            build: |m| {
                let c = cuboid(m, [-2.5; 3], [2.5; 3]);
                let c = turn(m, c, Axis::X, 37, [0, 0, 0]);
                (plate_with_bore(m), shift(m, c, [30, 30, 10]))
            },
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, bored + 62.5),
                Built(1, bored - 62.5),
                Built(1, 62.5),
            ],
        },
        // A pin of radius 2 on the plate; one wall plane of the far box cuts it.
        Family {
            name: "box turned 37 beside a pin",
            build: |m| {
                (
                    plate_with_pin(m, [10.0, 10.0, 10.0], 2.0, 3.0),
                    turned_box(m, 37, [30, 30], 5.0, 10.0),
                )
            },
            gate: Ok((3, 1, 0, 0)),
            want: [
                Built(1, 16000.0 + 12.0 * pi + 125.0),
                Built(1, 16000.0 + 12.0 * pi - 125.0),
                Built(1, 125.0),
            ],
        },
        // The box encloses a pin of radius 1: Common is `25·2 + 3π`, which a road that guessed
        // "outside" for the pin's disk answered as 50.
        Family {
            name: "box turned 37 around a pin",
            build: |m| {
                (
                    plate_with_pin(m, [30.0, 30.0, 10.0], 1.0, 3.0),
                    turned_box(m, 37, [30, 30], 8.0, 10.0),
                )
            },
            gate: Ok((4, 0, 0, 0)),
            want: refused,
        },
        // The box encloses the bore: each wall stands 2.5 from the axis, beyond the radius 2, which
        // the gate shows from the axis's distance to the plane once the box test cannot (the box of
        // the bore's lateral straddles every wall) — four beside. What refuses is the nesting, as
        // with the box turned around the pin.
        Family {
            name: "box turned 37 across the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [10, 10], 5.0, 10.0)),
            gate: Ok((4, 0, 0, 0)),
            want: refused,
        },
        // A lone cylinder and a box above its top, apart.
        Family {
            name: "lone cylinder, box turned 37 above, planes clear",
            build: |m| {
                (
                    z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0),
                    turned_box(m, 37, [30, 10], 20.0, 5.0),
                )
            },
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(2, 40.0 * pi + 125.0),
                Built(1, 40.0 * pi),
                Built(0, 0.0),
            ],
        },
        Family {
            name: "lone cylinder, box turned 37 above, a wall plane through it",
            build: |m| {
                (
                    z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0),
                    turned_box(m, 37, [30, 30], 20.0, 5.0),
                )
            },
            gate: Ok((3, 1, 0, 0)),
            want: [
                Built(2, 40.0 * pi + 125.0),
                Built(1, 40.0 * pi),
                Built(0, 0.0),
            ],
        },
        // The fillet plate with one fillet, the box in the open: its rings mix pierce corners
        // and three-plane corners, which is what sees a side recorded with the wrong sign.
        Family {
            name: "one-fillet plate, box turned 37 clear of the fillet",
            build: |m| {
                (
                    rounded_plate(m, 1, 0),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, rounded_area(1, 0) * 12.0 + 100.0),
                Built(1, rounded_area(1, 0) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // The playground's example 003 shape: four fillets, four bores, the box in the middle.
        Family {
            name: "four fillets and four bores, box turned 37 in the middle",
            build: |m| (rounded_plate(m, 4, 4), turned_box(m, 37, [0, 0], 6.0, 10.0)),
            gate: Ok((32, 0, 0, 0)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        // Four bores and no fillet.
        Family {
            name: "four bores, box turned 37",
            build: |m| {
                (
                    rounded_plate(m, 0, 4),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((14, 2, 0, 0)),
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
            name: "one-fillet plate, a wall plane through the fillet",
            build: |m| {
                (
                    rounded_plate(m, 1, 0),
                    turned_box(m, 37, [10, 16], 5.0, 10.0),
                )
            },
            gate: Ok((3, 0, 1, 0)),
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
            name: "slotted plate, one wall beside a slot end by its distance, one sectioning it",
            build: |m| (slotted_plate(m), turned_box(m, 37, [33, 11], 5.0, 10.0)),
            gate: Ok((7, 0, 1, 0)),
            want: [
                Built(1, slotted_plate_volume() + 125.0),
                Built(1, slotted_plate_volume() - 125.0),
                Built(1, 125.0),
            ],
        },
        // Two wall planes, each through one fillet: two sectioned pairs of the sixteen.
        Family {
            name: "four fillets, a wall plane through a fillet",
            build: |m| {
                (
                    rounded_plate(m, 4, 0),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((14, 0, 2, 0)),
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
        // it against the radius 3.5, by the box's centre: `(−12, −6)` — `(−38, 18)` 3.82 and
        // `(−38, −18)` 3.56, both beyond; `(−22, −14)` — 3.98 and 3.93, both beyond; `(−16, −10)` —
        // `(−38, 18)` 3.22 and 1.78 (two walls), within, and `(−38, −18)` 4.35, beyond; `(−20, 0)` —
        // 1.04 and 1.04, within.
        Family {
            name: "example-003 plate, box off the middle, a wall through one fillet",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-12, -6], 6.0, 10.0),
                )
            },
            gate: Ok((31, 0, 1, 0)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "example-003 plate, walls through two fillets",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-22, -14], 6.0, 10.0),
                )
            },
            gate: Ok((30, 0, 2, 0)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "example-003 plate, walls through two fillets and three bores' boxes",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-16, -10], 6.0, 10.0),
                )
            },
            gate: Ok((28, 2, 2, 0)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, rounded_area(4, 4) * 12.0 - 150.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "example-003 plate, walls through two fillets and two bores' boxes",
            build: |m| {
                (
                    rounded_plate(m, 4, 4),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Ok((28, 2, 2, 0)),
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
            name: "a pin over the plate's edge, a wall plane through it",
            build: |m| {
                (
                    plate_with_pin(m, [39.0, 20.0, 10.0], 2.0, 3.0),
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
        // tangent walls, so the gate has no other record to meet. Downstream, on the plate's top
        // and bottom, the split pass asks whether the class's line reaches the wall's segment
        // there: the crossing `{cap, class, y = −25}` against the segment's end, the fillet's
        // tangent point (a pierce name, root `Double`) — one point, two names, which the interval
        // containment cannot order — and refuses by the gate's name.
        Family {
            name: "one-fillet plate, a wall plane holding the fillet's tangent ruling",
            build: |m| {
                let plate = rounded_plate(m, 1, 0);
                let b = cuboid(m, [35.0, -80.0, 0.0], [40.0, -75.0, 10.0]);
                (plate, turn(m, b, Axis::Z, 37, [40, -25, 0]))
            },
            gate: Ok((3, 0, 1, 0)),
            want: refused,
        },
        // The operands the other way round: the box first. `B − A` is now the box less the plate,
        // `250 − 150`.
        Family {
            name: "example-003 plate, walls through two fillets, box first",
            build: |m| {
                let plate = rounded_plate(m, 4, 4);
                (turned_box(m, 37, [-22, -14], 6.0, 10.0), plate)
            },
            gate: Ok((30, 0, 2, 0)),
            want: [
                Built(1, rounded_area(4, 4) * 12.0 + 100.0),
                Built(1, 100.0),
                Built(1, 150.0),
            ],
        },
        Family {
            name: "one-fillet plate, a wall plane through the fillet, box first",
            build: |m| {
                let plate = rounded_plate(m, 1, 0);
                (turned_box(m, 37, [10, 16], 5.0, 10.0), plate)
            },
            gate: Ok((3, 0, 1, 0)),
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
        // sinks `25·5` into the plate. Its rib face covers none of the cells the order between
        // those points shapes, so this family cannot see that order's sign — the two with the wall across the bore can.
        Family {
            name: "recorded rib wall, box turned 37 across its plane",
            build: |m| {
                let p = plate_with_bore(m);
                let rib = cuboid(m, [9.0, 30.0, 5.0], [11.0, 40.0, 15.0]);
                let a = op(m, BoolKind::Fuse, p, rib);
                (a, turned_box(m, 37, [10, 25], 5.0, 10.0))
            },
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, bored + 225.0),
                Built(1, bored - 25.0),
                Built(1, 125.0),
            ],
        },
        // The wall runs the plate's length, across the bore, so its face holds both the bore's
        // pierce points and the box's crossings — the order between them is what shapes it.
        Family {
            name: "wall across the bore, box turned 37 across its plane",
            build: |m| (walled_bore(m), turned_box(m, 37, [10, 25], 5.0, 10.0)),
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, walled + 250.0 - sunk(5.0)),
                Built(1, walled - sunk(5.0)),
                Built(1, sunk(5.0)),
            ],
        },
        // The same, operands swapped: the class roots and their frames are the box's now. The box
        // less the walled plate is two solids — the wall's strip runs across the box's top half.
        Family {
            name: "wall across the bore, box turned 37 across its plane, operands swapped",
            build: |m| {
                let b = turned_box(m, 37, [10, 25], 5.0, 10.0);
                (b, walled_bore(m))
            },
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, walled + 250.0 - sunk(5.0)),
                Built(2, 250.0 - sunk(5.0)),
                Built(1, sunk(5.0)),
            ],
        },
        // The box raised to stand on `z = 8`, so less of it sinks into the plate.
        Family {
            name: "wall across the bore, box turned 37 standing higher",
            build: |m| (walled_bore(m), turned_box(m, 37, [10, 25], 8.0, 10.0)),
            gate: Ok((4, 0, 0, 0)),
            want: [
                Built(1, walled + 250.0 - sunk(8.0)),
                Built(1, walled - sunk(8.0)),
                Built(1, sunk(8.0)),
            ],
        },
        // A world wall crossing the bore, fused with a turned box: the wall's face carries the
        // box's irrational corners.
        Family {
            name: "named wall across the bore with irrational corners",
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
            name: "bored plate turned 37",
            build: |m| {
                let p = plate_with_bore(m);
                let p = turn(m, p, Axis::Z, 37, [20, 20, 0]);
                (p, turned_box(m, 0, [27, 27], 5.0, 10.0))
            },
            gate: Err(Undecided),
            want: refused,
        },
    ];
    score(&fams, false);
    // A sectioned class meeting another record of the fillet on one of its rulings.
    score(&hinge_on_a_fillet_ruling(false), false);
    score(&hinge_on_a_fillet_ruling(true), true);
}

/// **A sectioned class meeting another record of the fillet on one of its rulings.** The other
/// operand is one solid: the box `a0 = [41,60]×[−40,−24]×[2,8]`, whose walls `x = 41` and
/// `y = −24` cut the fillet (crossings; `y = −24` meets it at `(43, −24)`), fused with a box
/// whose wall `x = 43` is turned 37° about that ruling, overlapping the first far from the
/// plate. The turned wall's class runs along the axis within the radius (`≈ 0.011` from it)
/// with every face apart — sectioned — and it meets `x = 41` off the lateral (at
/// `y ≈ −21.35`), but meets `y = −24` on the ruling `(43, −24)` itself: the line the turn is
/// about, which the gate proves the class holds and records as a line the two share. The
/// three other turned walls clear the fillet's box.
///
/// The volumes: `a0` overlaps the plate in the strip `y ∈ [−25, −24]` past `x = 41`, cut by
/// the fillet's circle about `(40, −20)`: `A = ∫ (√(25 − u²) − 1) du` over
/// `u = y + 20 ∈ [−√24, −4]`, over the height 6. The turned box misses the plate (its `y` stays
/// below −37) and overlaps `a0` in its square's band `y ∈ [−40, −24]` over the height 6 — its
/// `x` range, `52.74 ± 3.49`, inside `a0`'s.
///
/// The two boxes fuse either way round, and that decides which of the two classes the table
/// numbers first ([`the_hinge_is_read_whichever_class_is_numbered_first`]). With `swapped` the
/// tool comes first, and `B − A` is the tool less the plate.
fn hinge_on_a_fillet_ruling(swapped: bool) -> [Family; 2] {
    let want = {
        let plate = rounded_area(1, 0) * 12.0;
        let f = |u: f64| u / 2.0 * (25.0 - u * u).sqrt() + 12.5 * (u / 5.0).asin();
        let r24 = 24f64.sqrt();
        let common = 6.0 * (f(r24) - f(4.0) - (r24 - 4.0));
        let (s, c) = 37f64.to_radians().sin_cos();
        let centre = [43.0 - 2.5 * c + 19.5 * s, -24.0 - 2.5 * s - 19.5 * c];
        let tool = 1824.0 + 250.0 - 6.0 * turned_square_within(centre, 37.0, 1, -40.0, -24.0);
        let cut = if swapped {
            tool - common
        } else {
            plate - common
        };
        [
            Built(1, plate + tool - common),
            Built(1, cut),
            Built(1, common),
        ]
    };
    let family = |name, build| Family {
        name,
        build,
        gate: Ok((3, 0, 1, 0)),
        want,
    };
    [
        family(
            "a sectioned class meets a recorded wall on a fillet ruling",
            |m| (rounded_plate(m, 1, 0), hinge_tool(m, false)),
        ),
        family(
            "a sectioned class meets a recorded wall on a fillet ruling, the turned box fused first",
            |m| (rounded_plate(m, 1, 0), hinge_tool(m, true)),
        ),
    ]
}

/// [`hinge_on_a_fillet_ruling`]'s other operand: `a0` and the turned box, fused in the order
/// `turned_first` says.
fn hinge_tool(m: &mut Model, turned_first: bool) -> Handle<Solid> {
    let a0 = cuboid(m, [41.0, -40.0, 2.0], [60.0, -24.0, 8.0]);
    let t = cuboid(m, [38.0, -46.0, 0.0], [43.0, -41.0, 10.0]);
    let t = turn(m, t, Axis::Z, 37, [43, -24, 0]);
    let (x, y) = if turned_first { (t, a0) } else { (a0, t) };
    let b = boolean(m, BoolKind::Fuse, x, y).expect("the two boxes fuse");
    assert_eq!(b.len(), 1, "the two boxes overlap into one solid");
    b[0]
}

/// **The hinge's record reads the same whichever of its two classes the table numbers first** —
/// the aliases' representative is the cutting class's exact pierce name, never the one the
/// sectioned class's interval places. The two fusing orders of [`hinge_on_a_fillet_ruling`]
/// number the sectioned class after the cutting one and before it.
#[test]
fn the_hinge_is_read_whichever_class_is_numbered_first() {
    let mut first = Vec::new();
    for f in hinge_on_a_fillet_ruling(false) {
        let mut m = Model::new();
        let (a, b) = (f.build)(&mut m);
        let setup = plane_index_setup(&m, a, b).expect("the gate passes");
        let nameless = |c: usize| setup.geom[c].world.is_none();
        let rows: Vec<bool> = setup
            .cyls
            .iter()
            .flat_map(|cy| &cy.shared)
            .filter(|sr| nameless(sr.line.0) != nameless(sr.line.1))
            .map(|sr| nameless(sr.line.0))
            .collect();
        assert_eq!(rows.len(), 1, "{}: one hinge record", f.name);
        first.push(rows[0]);
    }
    first.sort_unstable();
    assert_eq!(first, [false, true], "both numberings are asked");
}

/// **Each family's gate verdict and its three answers** — what [`Family`] pins, checked. With
/// `swapped`, every question takes the operands the other way round, `(b, a)`.
fn score(fams: &[Family], swapped: bool) {
    let built = |f: &Family, m: &mut Model| {
        let (a, b) = (f.build)(m);
        if swapped { (b, a) } else { (a, b) }
    };
    for f in fams {
        {
            let mut m = Model::new();
            let (a, b) = built(f, &mut m);
            let gate = plane_index_setup(&m, a, b)
                .map(|setup| {
                    let mut by_arm = (0, 0, 0, 0);
                    for cy in &setup.cyls {
                        for &(_, r) in &cy.records {
                            let crate::planes::CylRecord::Silent(s) = r else {
                                continue;
                            };
                            match s {
                                crate::planes::Silence::Beside(_) => by_arm.0 += 1,
                                crate::planes::Silence::Apart => by_arm.1 += 1,
                                crate::planes::Silence::Sectioned => by_arm.2 += 1,
                                crate::planes::Silence::Grazes => by_arm.3 += 1,
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
            assert_eq!(gate, f.gate, "{} (swapped: {swapped}): the gate", f.name);
        }
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(&f.want)
        {
            let mut m = Model::new();
            let (a, b) = built(f, &mut m);
            match (boolean(&mut m, kind, a, b), want) {
                (Ok(out), Built(n, v)) => {
                    assert_eq!(
                        out.len(),
                        *n,
                        "{} {kind:?} (swapped: {swapped}): solids",
                        f.name
                    );
                    m.rebuild_adjacency();
                    let issues = nacre_validate::validate(&m);
                    assert!(
                        issues.is_empty(),
                        "{} {kind:?} (swapped: {swapped}): {issues:?}",
                        f.name
                    );
                    let got: f64 = out
                        .iter()
                        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                        .sum();
                    assert!(
                        (got - v).abs() <= 1e-9 * v.max(1.0),
                        "{} {kind:?} (swapped: {swapped}): volume {got} against {v}",
                        f.name
                    );
                }
                (Err(BoolError::Rejected { reason, .. }), Refused(r)) => {
                    assert_eq!(
                        &reason, r,
                        "{} {kind:?} (swapped: {swapped}): the refusal's name",
                        f.name
                    );
                }
                (got, _) => panic!(
                    "{} {kind:?} (swapped: {swapped}): {got:?}, wanted {want:?}",
                    f.name
                ),
            }
        }
    }
}

/// The area of a `5×5` square turned `deg` about its centre `c`, within `lo ≤ p[k] ≤ hi` — the
/// section a turned box or cube shares with an axis-aligned slab. The square is convex, so
/// clipping it by the two half-planes in turn is exact.
fn turned_square_within(c: [f64; 2], deg: f64, k: usize, lo: f64, hi: f64) -> f64 {
    let (s, co) = deg.to_radians().sin_cos();
    let mut poly: Vec<[f64; 2]> = [[-2.5, -2.5], [2.5, -2.5], [2.5, 2.5], [-2.5, 2.5]]
        .iter()
        .map(|&[x, y]| [c[0] + co * x - s * y, c[1] + s * x + co * y])
        .collect();
    for (bound, below) in [(lo, false), (hi, true)] {
        let inside = |p: [f64; 2]| if below { p[k] <= bound } else { p[k] >= bound };
        let mut out = Vec::new();
        for i in 0..poly.len() {
            let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
            if inside(p) {
                out.push(p);
            }
            if inside(p) != inside(q) {
                let t = (bound - p[k]) / (q[k] - p[k]);
                out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
            }
        }
        poly = out;
    }
    let twice: f64 = (0..poly.len())
        .map(|i| {
            let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
            p[0] * q[1] - q[0] * p[1]
        })
        .sum();
    twice.abs() / 2.0
}

/// A bar `[x0, x1] × [0, 40] × [z0, z1]` fused with `with` — one tool carrying a named wall and the
/// nameless walls of `with`.
fn bar_and(m: &mut Model, x: [f64; 2], z: [f64; 2], with: Handle<Solid>) -> Handle<Solid> {
    let bar = cuboid(m, [x[0], 0.0, z[0]], [x[1], 40.0, z[1]]);
    op(m, BoolKind::Fuse, bar, with)
}

/// The bored plate with a box turned 37 at `c` fused into it.
fn bored_plate_and_box(m: &mut Model, c: [i128; 2]) -> Handle<Solid> {
    let p = plate_with_bore(m);
    let t = turned_box(m, 37, c, 5.0, 10.0);
    op(m, BoolKind::Fuse, p, t)
}

/// The box `[10,15] × [20,25] × [5,15]` turned −37 about `(10, 12)`, a point of the bore's circle:
/// one wall plane turns about it, so it meets the plane `x = 10` on that ruling.
fn pivoted_box(m: &mut Model) -> Handle<Solid> {
    let b = cuboid(m, [10.0, 20.0, 5.0], [15.0, 25.0, 15.0]);
    turn(m, b, Axis::Z, -37, [10, 12, 0])
}

/// **A nameless class apart from a cylinder passes beside the records other classes write on it.**
/// The class's plane may cut the bore or pin, its faces stand apart from the lateral's box, and
/// the cylinder has only whole rims (`Silence::Apart`); a named bar cuts the same cylinder or
/// touches it, so that cylinder carries a record too. On the apart class the solid reads as the
/// bore filled or the pin removed, and that differs from the truth only inside the lateral's box,
/// where the class has no face. Every family pins, in both operand orders, the gate's tally, a
/// record on each cylinder an apart class is silent about (else the family would be one of the
/// record-free apart families of the scoreboard above), and the three answers by derived volume.
///
/// The derivations: the bar's section of the plate is its width × 40 × 5, less the bore's strip
/// `|x − 10| ∈ [1, 2]` (`5·(4π/3 − √3)`; over the pin, plus `3·(4π/3 − √3)`) or the half disk
/// (`10π`) where the bar holds the axis. A box standing on `z = 5` sinks 125 into the plate; the
/// cube sinks its section below `z = 10` (`cube_low`) over its width 5. Where a box and the bar
/// overlap is the turned square's section by the bar's slab ([`turned_square_within`]) over the
/// height 10, half of it inside the plate; the cube overlaps the bar over the width 0.5 only, by
/// `cube_band`, and inside the plate by `cube_low`.
#[test]
fn an_apart_class_passes_beside_a_recorded_cylinder() {
    use RejectReason::SelfTouchingResult;
    let pi = std::f64::consts::PI;
    let bored = 16000.0 - 40.0 * pi;
    let pinned = 16000.0 + 12.0 * pi;
    let strip = 4.0 * pi / 3.0 - 3f64.sqrt();
    let box_at = |c: [f64; 2], x: [f64; 2]| turned_square_within(c, 37.0, 0, x[0], x[1]);
    // `(|A|, |B|, |A ∩ B|)` → Fuse, Cut, Common; the other order's Cut is `|B| − |A ∩ B|`.
    let three = |a: f64, b: f64, common: f64| {
        [
            Built(1, a + b - common),
            Built(1, a - common),
            Built(1, common),
        ]
    };
    let cut11 = box_at([16.0, 18.0], [11.0, 13.0]);
    let tan12 = box_at([16.0, 18.0], [12.0, 14.0]);
    let crossing = box_at([11.0, 16.0], [11.0, 13.0]);
    let pivot_centre = {
        let (s, c) = (-37f64).to_radians().sin_cos();
        [10.0 + 2.5 * c - 10.5 * s, 12.0 + 2.5 * s + 10.5 * c]
    };
    let ruling = turned_square_within(pivot_centre, -37.0, 0, 10.0, 15.0);
    let cube_low = turned_square_within([10.0, 12.0], 37.0, 1, f64::NEG_INFINITY, 10.0);
    let cube_band = turned_square_within([10.0, 12.0], 37.0, 1, 5.0, 15.0);
    // `(|B|, |A ∩ B|)` per family.
    let on_bar_bore = (
        1050.0 - 10.0 * cut11,
        (400.0 - 5.0 * strip) + 125.0 - 5.0 * cut11,
    );
    let on_bar_tangent = (1050.0 - 10.0 * tan12, 400.0 + 125.0 - 5.0 * tan12);
    let on_bar_pin = (
        1050.0 - 10.0 * cut11,
        (400.0 + 3.0 * strip) + 125.0 - 5.0 * cut11,
    );
    let on_ruling = (
        2250.0 - 10.0 * ruling,
        (1000.0 - 10.0 * pi) + 125.0 - 5.0 * ruling,
    );
    let on_plate = (800.0, (400.0 - 5.0 * strip) + 5.0 * cut11);
    let across = (800.0, (400.0 - 5.0 * strip) + 5.0 * crossing);
    let oblique = (
        6800.0 + 125.0 - 0.5 * cube_band,
        (3400.0 - 5.0 * strip) + 5.0 * cube_low - 0.5 * cube_low,
    );
    let fams = [
        // A bar `x ∈ [11, 13]` cuts the bore. The box at `(16, 18)`: its wall planes stand 0.28
        // (within the radius 2 — apart), 5.28, 7.11 and 12.11 (beside) from the axis.
        (
            Family {
                name: "bar cutting the bore, box apart on the bar",
                build: |m| {
                    let t = turned_box(m, 37, [16, 18], 5.0, 10.0);
                    (plate_with_bore(m), bar_and(m, [11.0, 13.0], [5.0, 15.0], t))
                },
                gate: Ok((3, 1, 0, 0)),
                want: three(bored, on_bar_bore.0, on_bar_bore.1),
            },
            Built(1, on_bar_bore.0 - on_bar_bore.1),
        ),
        // The bar `x ∈ [12, 14]` touches the bore along its line `x = 12`; `A − B` is the plate
        // touching itself there, as with the bar alone. The box as above.
        (
            Family {
                name: "bar tangent to the bore, box apart on the bar",
                build: |m| {
                    let t = turned_box(m, 37, [16, 18], 5.0, 10.0);
                    (plate_with_bore(m), bar_and(m, [12.0, 14.0], [5.0, 15.0], t))
                },
                gate: Ok((3, 1, 0, 0)),
                want: [
                    Built(1, bored + on_bar_tangent.0 - on_bar_tangent.1),
                    Refused(SelfTouchingResult),
                    Built(1, on_bar_tangent.1),
                ],
            },
            Built(1, on_bar_tangent.0 - on_bar_tangent.1),
        ),
        // The same tool over a pin of radius 2 standing on the plate; the box as above.
        (
            Family {
                name: "bar cutting the pin, box apart on the bar",
                build: |m| {
                    let t = turned_box(m, 37, [16, 18], 5.0, 10.0);
                    (
                        plate_with_pin(m, [10.0, 10.0, 10.0], 2.0, 3.0),
                        bar_and(m, [11.0, 13.0], [5.0, 15.0], t),
                    )
                },
                gate: Ok((3, 1, 0, 0)),
                want: three(pinned, on_bar_pin.0, on_bar_pin.1),
            },
            Built(1, on_bar_pin.0 - on_bar_pin.1),
        ),
        // A bar through the axis, `x ∈ [10, 15]`, and the pivoted box: the wall plane turned about
        // `(10, 12)` meets `x = 10` on the bar's ruling there. Its wall planes stand 1.20 (apart),
        // 3.80, 9.60 and 14.60 (beside) from the axis.
        (
            Family {
                name: "box wall plane through the bar's ruling",
                build: |m| {
                    let t = pivoted_box(m);
                    (plate_with_bore(m), bar_and(m, [10.0, 15.0], [5.0, 15.0], t))
                },
                gate: Ok((3, 1, 0, 0)),
                want: three(bored, on_ruling.0, on_ruling.1),
            },
            Built(1, on_ruling.0 - on_ruling.1),
        ),
        // The apart box on the bored plate's side, the bar alone on the other: each class's
        // arrangement sections the other operand. The box at `(16, 18)` as above.
        (
            Family {
                name: "box apart on the plate, bar cutting the bore",
                build: |m| {
                    (
                        bored_plate_and_box(m, [16, 18]),
                        cuboid(m, [11.0, 0.0, 5.0], [13.0, 40.0, 15.0]),
                    )
                },
                gate: Ok((3, 1, 0, 0)),
                want: three(bored + 125.0, on_plate.0, on_plate.1),
            },
            Built(1, on_plate.0 - on_plate.1),
        ),
        // The box moved to `(11, 16)`: its faces cross the bar's plane `x = 11` off the lateral's
        // box. Its wall planes stand 1.69 and 1.91 (apart), 6.69 and 6.91 (beside) from the axis.
        // `B − A` is two bodies: the box cuts the bar's upper half through.
        (
            Family {
                name: "box apart on the plate, its faces across the bar's plane",
                build: |m| {
                    (
                        bored_plate_and_box(m, [11, 16]),
                        cuboid(m, [11.0, 0.0, 5.0], [13.0, 40.0, 15.0]),
                    )
                },
                gate: Ok((2, 2, 0, 0)),
                want: three(bored + 125.0, across.0, across.1),
            },
            Built(2, across.0 - across.1),
        ),
        // A cube turned 37 about `x` at `(30, 10, 12)`, reached from the bore by the bar
        // `[11, 28]`. Of its four walls oblique to the bore's axis, two meet the axis at heights
        // 7.85 and 8.87 (apart); the other two stand above the lateral's box, their lowest height
        // over it 13.50 and 13.62 (beside). Its walls `x = 27.5, 32.5` have names.
        (
            Family {
                name: "oblique walls apart, bar cutting the bore",
                build: |m| {
                    let c = cuboid(m, [-2.5; 3], [2.5; 3]);
                    let c = turn(m, c, Axis::X, 37, [0, 0, 0]);
                    let c = shift(m, c, [30, 10, 12]);
                    (plate_with_bore(m), bar_and(m, [11.0, 28.0], [5.0, 15.0], c))
                },
                gate: Ok((2, 2, 0, 0)),
                want: three(bored, oblique.0, oblique.1),
            },
            Built(1, oblique.0 - oblique.1),
        ),
    ];
    for swapped in [false, true] {
        for (f, _) in &fams {
            let mut m = Model::new();
            let (a, b) = (f.build)(&mut m);
            let (a, b) = if swapped { (b, a) } else { (a, b) };
            let setup = plane_index_setup(&m, a, b).expect("the gate passes");
            for (ci, cy) in setup.cyls.iter().enumerate() {
                let apart = crate::planes::CylRecord::Silent(crate::planes::Silence::Apart);
                if !cy.records.iter().any(|&(_, r)| r == apart) {
                    continue;
                }
                let recorded = cy
                    .records
                    .iter()
                    .any(|&(_, r)| r == crate::planes::CylRecord::Cuts)
                    || setup.tangencies.iter().any(|t| t.cyl == ci);
                assert!(
                    recorded,
                    "{} (swapped: {swapped}): the bar leaves a record on cylinder {ci}",
                    f.name
                );
            }
        }
    }
    let ab: Vec<Family> = fams.iter().map(|(f, _)| *f).collect();
    let ba: Vec<Family> = fams
        .iter()
        .map(|(f, cut)| Family {
            want: [f.want[0], *cut, f.want[2]],
            ..*f
        })
        .collect();
    score(&ab, false);
    score(&ba, true);
}

/// **A witness whose name holds a silent class is asked after the ray road declines.** The ray road
/// declines a ring with a pierce corner from a ray plane with a world statement, and `nesting`
/// then skips the named witnesses that would put the same question to it — but `ask` answers a
/// name holding a class the gate passed silent by the ring's box, so that one is no repeat. In this
/// family (the recorded rib wall, a box turned 37 across its plane) twelve questions meet such a
/// witness after the road declined, and each is asked and answers. Read back from the skip's tally
/// (`nesting_probe::DECLINES`):
/// `silent_after` counts them, `named_after` the witnesses the skip should have passed over.
#[test]
fn a_witness_naming_a_silent_class_is_asked_after_the_road_declines() {
    let rows = crate::ledger::owned(|| {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let mut m = Model::new();
            let p = plate_with_bore(&mut m);
            let rib = cuboid(&mut m, [9.0, 30.0, 5.0], [11.0, 40.0, 15.0]);
            let a = op(&mut m, BoolKind::Fuse, p, rib);
            let b = turned_box(&mut m, 37, [10, 25], 5.0, 10.0);
            m.rebuild_adjacency();
            crate::boolean(&mut m, kind, a, b).expect("the turned box meets the ribbed plate");
        }
        crate::nesting::nesting_probe::DECLINES.mine()
    });
    assert!(
        rows.iter().any(|&(_, _, silent_after)| silent_after > 0),
        "no witness naming a silent class was asked after a declining road: {rows:?}"
    );
    for (declines, named_after, _) in rows {
        assert!(declines > 0);
        assert_eq!(
            named_after, 0,
            "a repeating witness was asked after its road declined"
        );
    }
}

/// The half of a `z` cylinder of radius 2 standing on `(10, 10)` over `z ∈ [0, 10]` that lies on
/// `x ≤ 10` — its lateral has no whole rim, so a wall plane cutting it is sectioned. `20π`.
fn half_bar(m: &mut Model) -> Handle<Solid> {
    let c = z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0);
    let cut = cuboid(m, [10.0, 0.0, -1.0], [20.0, 20.0, 11.0]);
    op(m, BoolKind::Cut, c, cut)
}

/// A `5×5` box over `z`, turned 37° about `z` through the half bar's axis `(10, 10)`. Before the
/// turn it spans `x ∈ [8.5, 13.5]` and `y ∈ [10 + y0, 15 + y0]`, so its wall `x = 8.5` turns into
/// a plane 1.5 from the axis on the bar's round side, crossing its arc twice and missing its flat
/// face (which reaches only 1.2 along that plane's normal); the other three walls stand 3.5, `y0`
/// and `y0 + 5` from the axis.
fn box_turned_about_the_bar(m: &mut Model, y0: f64, z: [f64; 2]) -> Handle<Solid> {
    let b = cuboid(m, [8.5, 10.0 + y0, z[0]], [13.5, 15.0 + y0, z[1]]);
    turn(m, b, Axis::Z, 37, [10, 10, 0])
}

/// **A section whose every corner is a pierce point is asked by a pierce corner.** A wall plane
/// of a turned box crosses a half bar's round side twice: on that plane the bar's section is the
/// strip between two rulings, closed by the two caps, and all four corners are pierce points of
/// a class with no world statement — no three-plane name, no rational coordinate, no centre. Its
/// nesting against the plane's other rings is asked from those corners, by rays along the cap's
/// line on the wall plane, the corner placed on that line by the cylinder.
///
/// The gate: the wall 1.5 from the axis sections the bar, the other three walls stand beside it —
/// `(3, 0, 1)`, every family. The answers, by derived volume, in both operand orders:
/// - a box above the bar or beside it on the floor: the Cut and the Common build. The Fuse is two
///   bodies, and the depth of one in the other is a question in space (`assembly::grouping`) whose
///   rays from the box's corners run on the wall plane, where the bar's faces cannot be read yet —
///   refused by name.
/// - the bar inside a named block `[0,20]² × [−2,12]` with the box on the block: on the wall plane
///   the bar's section lies inside the block's, and its corners answer «inside». The Fuse is the
///   block and box, `5600 + 125`; the bar less the block is empty; the block less the bar keeps a
///   cavity, whose depth in space is asked by a ray that meets a face the wall plane bounds or
///   names a corner of, and cannot show it misses — refused by the gate's name.
#[test]
fn a_section_whose_every_corner_is_a_pierce_point_is_asked_by_one() {
    use RejectReason::{CylinderGateUndecided, NoClearRay};
    let pi = std::f64::consts::PI;
    let bar = 20.0 * pi;
    // Each family with the other order's Cut.
    let fams = [
        (
            Family {
                name: "a half bar, a turned box above it, a wall plane through the bar's round side only",
                build: |m| (half_bar(m), box_turned_about_the_bar(m, 6.0, [20.0, 25.0])),
                gate: Ok((3, 0, 1, 0)),
                want: [Refused(NoClearRay), Built(1, bar), Built(0, 0.0)],
            },
            Built(1, 125.0),
        ),
        (
            Family {
                name: "a half bar, a turned box beside it on the floor, a wall plane through the bar's round side only",
                build: |m| (half_bar(m), box_turned_about_the_bar(m, 6.0, [0.0, 5.0])),
                gate: Ok((3, 0, 1, 0)),
                want: [Refused(NoClearRay), Built(1, bar), Built(0, 0.0)],
            },
            Built(1, 125.0),
        ),
        (
            Family {
                name: "a half bar inside a named block, the turned box on the block",
                build: |m| {
                    let b = half_bar(m);
                    let block = cuboid(m, [0.0, 0.0, -2.0], [20.0, 20.0, 12.0]);
                    let t = box_turned_about_the_bar(m, 4.0, [12.0, 17.0]);
                    (b, op(m, BoolKind::Fuse, block, t))
                },
                gate: Ok((3, 0, 1, 0)),
                want: [Built(1, 5725.0), Built(0, 0.0), Built(1, bar)],
            },
            Refused(CylinderGateUndecided),
        ),
    ];
    let ab: Vec<Family> = fams.iter().map(|(f, _)| *f).collect();
    let ba: Vec<Family> = fams
        .iter()
        .map(|(f, cut)| Family {
            want: [f.want[0], *cut, f.want[2]],
            ..*f
        })
        .collect();
    score(&ab, false);
    score(&ba, true);
}

/// A `5×5` box over `z`, `[x0, x0 + 5] × [10 + y0, 15 + y0]` before it is turned `deg` about `z`
/// through `(10, 10)` — the axis of [`half_bar`] and [`whole_bar`]. A turn about a cylinder's own
/// axis keeps each wall's distance from it, so the wall `x = 8` (or `x = 12`) of radius-2 bars
/// touches the bar along one ruling after the turn as before it.
fn box_turned_about_the_axis(
    m: &mut Model,
    deg: i128,
    x0: f64,
    y0: f64,
    z: [f64; 2],
) -> Handle<Solid> {
    let b = cuboid(m, [x0, 10.0 + y0, z[0]], [x0 + 5.0, 15.0 + y0, z[1]]);
    turn(m, b, Axis::Z, deg, [10, 10, 0])
}

/// A `z` cylinder of radius 2 standing on `(10, 10)` over `z ∈ [0, 10]` — whole rims. `40π`.
fn whole_bar(m: &mut Model) -> Handle<Solid> {
    z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0)
}

/// The triangle `(0, 0), (20, 20), (0, 20)` extruded over `z` — a named face on `x = y`.
fn prism_on_the_diagonal(m: &mut Model, z: [f64; 2]) -> Handle<Solid> {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let tri = crate::Profile2d::polygon(vec![p2(0.0, 0.0), p2(20.0, 20.0), p2(0.0, 20.0)])
        .expect("a triangle");
    let frame = crate::SketchFrame::world(m, Axis::Z);
    let Ok(crate::OpOutput::Extrude { solid, .. }) = crate::apply(
        m,
        &crate::Operation::Extrude {
            frame,
            profile: tri,
            dist: z[1] - z[0],
        },
    ) else {
        panic!("the triangle extrudes")
    };
    m.rebuild_adjacency();
    shift(m, solid, [0, 0, z[0] as i128])
}

/// **A wall turned about the cylinder's own axis is read before its turn.** A box turned 37° about
/// a bar's axis has a wall `x = 8` (`x = 12`) that touched the radius-2 bar along one ruling before
/// the turn, and the turn leaves the bar in place, so it touches it still; the interval road cannot
/// prove a distance equal to the radius, and the plane carried back before the turn does
/// (`Model::cylinder_before_chain`). The class grazes the bar (`Silence::Grazes`, the fourth count)
/// whatever its rims.
///
/// The answers, by derived volume, in both operand orders:
/// - a half bar, the box above it — the grazing wall with a wall within the radius whose line with
///   it misses the surface, or with a wall through the axis meeting it on the ruling that wall draws
///   (recorded, `GrazeLine`): the Cut and the Common build. The Fuse is two bodies whose depth is a
///   question in space (`assembly::grouping`), refused by name as for the half bar's other boxes
///   (todo «3차원 덩어리 판정의 광선»).
/// - the half bar inside a named block `[0,20]² × [−2,12]`, the through-axis box on the block: the
///   bar's section on the two walls answers the Fuse and the Common; the block less the bar asks the
///   same question in space and is refused by name.
/// - a whole bar, the same boxes and block: the grazing wall is the fourth count (the filled reading
///   would also take it), the through-axis wall filled (`Apart`), and everything builds but the block
///   less the bar, refused as above.
/// - a named prism on `x = y` across the half bar and a box turned 45° whose wall `x = 8` touches the
///   bar on the ruling `x = y` cuts: a record's line with the grazing class on a lateral face that is
///   not a sectioned class's is refused at the gate. No face of either class reaches that line here
///   (the box's wall stops short of the tangent point), and without the refusal both orders build
///   their derived volumes — `Common = 12π`, the 135° sector of the half disk on `y ≥ x` over the
///   overlap `z ∈ [2, 10]`. The refusal is wider than the harm it stands for, a face on the line
///   meeting a ruling end in one point under two names; that shape is refused earlier today (a named
///   face carrying a corner a nameless class names), so the narrower question cannot be measured
///   yet. The box's corners stay off `x = 10`: a turn of 45° that puts one there ties it with the
///   half bar's flat face in the winding's extreme-point search, a coincidence of the fixture.
#[test]
fn a_wall_turned_about_the_axis_grazing_the_cylinder_is_read_before_its_turn() {
    use RejectReason::{CylinderGateUndecided, NoClearRay};
    let pi = std::f64::consts::PI;
    let (half, whole) = (20.0 * pi, 40.0 * pi);
    // Each family with the other order's Fuse and Cut.
    let fams = [
        (
            Family {
                name: "a half bar, a box above turned about its axis: a grazing wall, a wall within the radius",
                build: |m| {
                    (
                        half_bar(m),
                        box_turned_about_the_axis(m, 37, 8.0, 1.0, [20.0, 25.0]),
                    )
                },
                gate: Ok((2, 0, 1, 1)),
                want: [
                    Refused(CylinderGateUndecided),
                    Built(1, half),
                    Built(0, 0.0),
                ],
            },
            [Refused(NoClearRay), Built(1, 125.0)],
            0,
        ),
        (
            Family {
                name: "a half bar, a box above turned about its axis: a grazing wall, a wall through the axis",
                build: |m| {
                    (
                        half_bar(m),
                        box_turned_about_the_axis(m, 37, 8.0, 0.0, [20.0, 25.0]),
                    )
                },
                gate: Ok((2, 0, 1, 1)),
                want: [
                    Refused(CylinderGateUndecided),
                    Built(1, half),
                    Built(0, 0.0),
                ],
            },
            [Refused(NoClearRay), Built(1, 125.0)],
            1,
        ),
        (
            Family {
                name: "a half bar inside a named block, the through-axis box on the block",
                build: |m| {
                    let b = half_bar(m);
                    let block = cuboid(m, [0.0, 0.0, -2.0], [20.0, 20.0, 12.0]);
                    let t = box_turned_about_the_axis(m, 37, 8.0, 0.0, [12.0, 17.0]);
                    (b, op(m, BoolKind::Fuse, block, t))
                },
                gate: Ok((2, 0, 1, 1)),
                want: [Built(1, 5725.0), Built(0, 0.0), Built(1, half)],
            },
            [Built(1, 5725.0), Refused(CylinderGateUndecided)],
            1,
        ),
        (
            Family {
                name: "a whole bar, a box above turned about its axis: a grazing wall, a wall through the axis",
                build: |m| {
                    (
                        whole_bar(m),
                        box_turned_about_the_axis(m, 37, 12.0, 0.0, [20.0, 25.0]),
                    )
                },
                gate: Ok((2, 1, 0, 1)),
                want: [Built(2, whole + 125.0), Built(1, whole), Built(0, 0.0)],
            },
            [Built(2, whole + 125.0), Built(1, 125.0)],
            0,
        ),
        (
            Family {
                name: "a whole bar, a box above turned about its axis: a grazing wall, a wall within the radius",
                build: |m| {
                    (
                        whole_bar(m),
                        box_turned_about_the_axis(m, 37, 12.0, 1.0, [20.0, 25.0]),
                    )
                },
                gate: Ok((2, 1, 0, 1)),
                want: [Built(2, whole + 125.0), Built(1, whole), Built(0, 0.0)],
            },
            [Built(2, whole + 125.0), Built(1, 125.0)],
            0,
        ),
        (
            Family {
                name: "a whole bar, a box above turned about its axis: a grazing wall, the rest beside",
                build: |m| {
                    (
                        whole_bar(m),
                        box_turned_about_the_axis(m, 37, 12.0, 3.0, [20.0, 25.0]),
                    )
                },
                gate: Ok((3, 0, 0, 1)),
                want: [Built(2, whole + 125.0), Built(1, whole), Built(0, 0.0)],
            },
            [Built(2, whole + 125.0), Built(1, 125.0)],
            0,
        ),
        (
            Family {
                name: "a whole bar inside a named block, the through-axis box on the block",
                build: |m| {
                    let b = whole_bar(m);
                    let block = cuboid(m, [0.0, 0.0, -2.0], [20.0, 20.0, 12.0]);
                    let t = box_turned_about_the_axis(m, 37, 12.0, 0.0, [12.0, 17.0]);
                    (b, op(m, BoolKind::Fuse, block, t))
                },
                gate: Ok((2, 1, 0, 1)),
                want: [Built(1, 5725.0), Built(0, 0.0), Built(1, whole)],
            },
            [Built(1, 5725.0), Refused(CylinderGateUndecided)],
            0,
        ),
        (
            Family {
                name: "a half bar, a named prism on x = y across it, a box turned 45 grazing it on that ruling",
                build: |m| {
                    let b = half_bar(m);
                    let p = prism_on_the_diagonal(m, [2.0, 22.0]);
                    let t = box_turned_about_the_axis(m, 45, 8.0, 4.0, [20.0, 25.0]);
                    (b, op(m, BoolKind::Fuse, p, t))
                },
                gate: Err(CylinderGateUndecided),
                want: [Refused(CylinderGateUndecided); 3],
            },
            [Refused(CylinderGateUndecided); 2],
            0,
        ),
    ];
    // The line a sectioned class and the grazing class share is recorded once per cylinder, and
    // only there.
    for (f, _, lines) in &fams {
        if f.gate.is_err() {
            continue;
        }
        for swapped in [false, true] {
            let mut m = Model::new();
            let (a, b) = (f.build)(&mut m);
            let (a, b) = if swapped { (b, a) } else { (a, b) };
            let setup = plane_index_setup(&m, a, b).expect("the gate passes");
            let recorded: usize = setup.cyls.iter().map(|cy| cy.grazed.len()).sum();
            assert_eq!(
                recorded, *lines,
                "{} (swapped: {swapped}): the shared lines",
                f.name
            );
        }
    }
    let ab: Vec<Family> = fams.iter().map(|(f, _, _)| *f).collect();
    let ba: Vec<Family> = fams
        .iter()
        .map(|(f, [fuse, cut], _)| Family {
            want: [*fuse, *cut, f.want[2]],
            ..*f
        })
        .collect();
    score(&ab, false);
    score(&ba, true);
}
