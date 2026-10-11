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

/// What one kind answers.
#[derive(Clone, Copy, Debug)]
enum Want {
    /// This many solids of this total volume, `validate` clean.
    Built(usize, f64),
    /// Refused by this name.
    Refused(RejectReason),
    /// Refused, by a name the scoreboard does not pin yet: the gate passes the family and a road
    /// that reads the gate's recorded side has not been opened.
    RefusedUnpinned,
}
use Want::{Built, Refused, RefusedUnpinned};

struct Family {
    name: &'static str,
    build: fn(&mut Model) -> (Handle<Solid>, Handle<Solid>),
    /// The gate passes, with this many nameless classes silent about every cylinder, or refuses.
    gate: Result<usize, RejectReason>,
    /// Fuse, Cut, Common.
    want: [Want; 3],
}

#[test]
fn nameless_classes_beside_world_cylinders() {
    use RejectReason::CylinderGateUndecided as Undecided;
    let pi = std::f64::consts::PI;
    let bored = 16000.0 - 40.0 * pi;
    let refused = [Refused(Undecided); 3];
    let unpinned = [RefusedUnpinned; 3];
    let fams = [
        // The control: a quarter turn keeps every name, and the box sits on the plate far from
        // the bore — it overlaps the plate over `z ∈ [5, 10]`, 125.
        Family {
            name: "F0 box turned 90",
            build: |m| (plate_with_bore(m), turned_box(m, 90, [30, 30], 5.0, 10.0)),
            gate: Ok(0),
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
            gate: Ok(4),
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
            gate: Err(Undecided),
            want: refused,
        },
        // The same box lifted off the plate: two bodies under Fuse, the plate under Cut.
        Family {
            name: "F1b box turned 37, above the plate",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [30, 30], 20.0, 5.0)),
            gate: Err(Undecided),
            want: refused,
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
            gate: Ok(4),
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
            gate: Err(Undecided),
            want: refused,
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
            gate: Ok(4),
            want: refused,
        },
        // The box crosses the bore: a true irrational meet.
        Family {
            name: "F1c box turned 37 across the bore",
            build: |m| (plate_with_bore(m), turned_box(m, 37, [10, 10], 5.0, 10.0)),
            gate: Err(Undecided),
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
            gate: Ok(4),
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
            gate: Err(Undecided),
            want: refused,
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
            gate: Ok(4),
            want: unpinned,
        },
        // The playground's example 003 shape: four fillets, four bores, the box in the middle.
        Family {
            name: "F16 four fillets and four bores, box turned 37 in the middle",
            build: |m| (rounded_plate(m, 4, 4), turned_box(m, 37, [0, 0], 6.0, 10.0)),
            gate: Ok(4),
            want: unpinned,
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
            gate: Err(Undecided),
            want: refused,
        },
        // A wall plane of the box cuts a fillet, whose lateral has ruling edges.
        Family {
            name: "F8 one-fillet plate, a wall plane through the fillet",
            build: |m| {
                (
                    rounded_plate(m, 1, 0),
                    turned_box(m, 37, [10, 16], 5.0, 10.0),
                )
            },
            gate: Err(Undecided),
            want: refused,
        },
        Family {
            name: "V1 four fillets, a wall plane through a fillet",
            build: |m| {
                (
                    rounded_plate(m, 4, 0),
                    turned_box(m, 37, [-20, 0], 6.0, 10.0),
                )
            },
            gate: Err(Undecided),
            want: refused,
        },
        // A rib on the plate whose wall `x = 9` runs within the bore's radius (recorded), and a
        // turned box straddling that wall's plane far from the bore.
        Family {
            name: "F14 recorded rib wall, box turned 37 across its plane",
            build: |m| {
                let p = plate_with_bore(m);
                let rib = cuboid(m, [9.0, 30.0, 5.0], [11.0, 40.0, 15.0]);
                let a = op(m, BoolKind::Fuse, p, rib);
                (a, turned_box(m, 37, [10, 25], 5.0, 10.0))
            },
            gate: Ok(4),
            want: unpinned,
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
                    let n = setup.cyls.first().map_or(0, |cy| cy.silent.len());
                    assert!(
                        setup.cyls.iter().all(|cy| cy.silent.len() == n),
                        "{}: every cylinder hears the same silence",
                        f.name
                    );
                    n
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
                (Err(BoolError::Rejected { .. }), RefusedUnpinned) => {}
                (got, _) => panic!("{} {kind:?}: {got:?}, wanted {want:?}", f.name),
            }
        }
    }
}
