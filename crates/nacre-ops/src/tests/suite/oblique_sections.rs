//! **Named plane classes oblique to a cylinder** — the scoreboard.
//!
//! A gusset's slanted face lies on a plane oblique to a bore's axis: the plane meets the cylinder
//! in an ellipse. Each family pairs such a class with world cylinders (a bore, two bores, a tube, a
//! pin, a fillet) and pins two answers: what the gate says, read from [`plane_index_setup`]
//! directly — how many (cylinder, class) pairs it records as the class's section of the cylinder —
//! and what the boolean answers, by kind and in both operand orders ([`every_kind`]). Volumes are
//! derived, not copied.
//!
//! The geometry, shared by the families: the bored plate of [`plate_with_bore`] (`[0,40]²×[0,10]`,
//! a bore of radius 2 about `(10, 10)`), and a gusset — a triangle in the `yz` plane extruded one
//! unit along `x` — whose slanted face lies on `z = y/4 + 5/2`. That plane meets the bore at
//! `z ∈ [4.5, 5.5]`, strictly between the plate's faces, while the gusset stands at `x ∈ [30, 31]`,
//! far from every lateral face.

use super::*;
use Want::{Built, Refused};

type Sketch = [(f64, f64); 3];

/// The gusset on the slanted plane: `(y, z)` corners `(14, 6)`, `(34, 6)`, `(34, 11)` — area 50.
/// In the plate it sinks to the top face `z = 10`, which it crosses at `y = 30`: 48 below, 2 above.
const SUNK: Sketch = [(14.0, 6.0), (34.0, 6.0), (34.0, 11.0)];

/// [`SUNK`] slid up its own plane past the plate's `y = 40` wall — the same class, apart from the
/// plate.
const APART: Sketch = [(54.0, 16.0), (74.0, 16.0), (74.0, 21.0)];

fn cylinder(
    m: &mut Model,
    at: [f64; 3],
    axis: [f64; 3],
    seam: [f64; 3],
    r: f64,
    h: f64,
) -> Handle<Solid> {
    let s = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array(at),
        Vector3::from_array(axis),
        Vector3::from_array(seam),
        r,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    s
}

fn z_cylinder(m: &mut Model, at: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
    cylinder(m, at, [0.0, 0.0, 1.0], [0.0, -1.0, 0.0], r, h)
}

fn op(m: &mut Model, kind: BoolKind, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
    let out = boolean(m, kind, a, b).expect("the fixture's boolean builds")[0];
    m.rebuild_adjacency();
    out
}

/// A triangle in the `yz` plane, `(y, z)` corners counter-clockwise, extruded over
/// `x ∈ [x0, x0 + 1]`.
fn gusset(m: &mut Model, tri: Sketch, x0: f64) -> Handle<Solid> {
    let p2 = |(y, z): (f64, f64)| nacre_math::Point2::from_array([y, z]);
    let edges = vec![
        line(p2(tri[0]), p2(tri[1])),
        line(p2(tri[1]), p2(tri[2])),
        line(p2(tri[2]), p2(tri[0])),
    ];
    let profile = stated(edges).unwrap().remove(0);
    let frame = SketchFrame::world(m, Axis::X);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    )
    .expect("the gusset extrudes") else {
        unreachable!()
    };
    let q = |x: f64| nacre_exact::Rat::from_decimal(x).expect("a short decimal");
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid,
            isometry: nacre_exact::Isometry::translation([q(x0), q(0.0), q(0.0)]),
        },
    )
    .expect("the gusset moves") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// Whether the gate recorded `r` as a class's section of a cylinder drawn on the class's own
/// arrangement. No record states one: the gate refuses every named class oblique to a cylinder
/// whose plane meets a lateral face.
fn states_a_section(_r: crate::planes::CylRecord) -> bool {
    false
}

/// **The gate's verdict, in both operand orders** — `Ok` with the number of (cylinder, class)
/// pairs it records as a section, or the refusal's name.
fn gate(
    build: fn(&mut Model) -> (Handle<Solid>, Handle<Solid>),
    want: Result<usize, RejectReason>,
) {
    for swap in [false, true] {
        let mut m = Model::new();
        let (a, b) = build(&mut m);
        let (x, y) = if swap { (b, a) } else { (a, b) };
        let got = plane_index_setup(&m, x, y)
            .map(|setup| {
                setup
                    .cyls
                    .iter()
                    .map(|cy| {
                        cy.records
                            .iter()
                            .filter(|&&(_, r)| states_a_section(r))
                            .count()
                    })
                    .sum::<usize>()
            })
            .map_err(|e| match e {
                BoolError::Rejected { reason, .. } => reason,
                other => panic!("swap={swap}: the gate failed outside its vocabulary: {other:?}"),
            });
        assert_eq!(got, want, "swap={swap}: the gate");
    }
}

fn pi() -> f64 {
    std::f64::consts::PI
}

/// The bored plate: `16000 − 40π`.
fn bored() -> f64 {
    16000.0 - 40.0 * pi()
}

const OBLIQUE: Want = Refused(RejectReason::ObliqueCylinderCut);
const REFUSED: [[Want; 3]; 2] = [[OBLIQUE; 3]; 2];

/// The sunk gusset beside the bored plate: the class's section of the bore is a whole ellipse
/// between the plate's faces, and the gate refuses the pair.
#[test]
fn a_gusset_sunk_into_a_bored_plate() {
    let build = |m: &mut Model| (plate_with_bore(m), gusset(m, SUNK, 30.0));
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// The same plane, the gusset clear of the plate.
#[test]
fn a_gusset_apart_from_a_bored_plate() {
    let build = |m: &mut Model| (plate_with_bore(m), gusset(m, APART, 30.0));
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// Two bores, `(10, 10)` and `(10, 20)`, both met by the plane between the plate's faces (the
/// second at `z ∈ [7, 8]`): two sections of parallel cylinders on one class.
#[test]
fn a_gusset_sunk_into_a_plate_with_two_bores() {
    let build = |m: &mut Model| {
        let p = plate_with_bore(m);
        let second = z_cylinder(m, [10.0, 20.0, -1.0], 2.0, 12.0);
        (op(m, BoolKind::Cut, p, second), gusset(m, SUNK, 30.0))
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// A tube — radius 2 about `(10, 10)` less radius 1, `z ∈ [0, 10]`, volume `30π` — beside the
/// gusset: one section inside the other on the class.
#[test]
fn a_gusset_beside_a_tube() {
    let build = |m: &mut Model| {
        let outer = z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0);
        let inner = z_cylinder(m, [10.0, 10.0, -1.0], 1.0, 12.0);
        (op(m, BoolKind::Cut, outer, inner), gusset(m, SUNK, 30.0))
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// A pin — radius 2 about `(10, 10)`, `z ∈ [0, 10]`, volume `40π` — beside the gusset.
#[test]
fn a_gusset_beside_a_pin() {
    let build = |m: &mut Model| {
        (
            z_cylinder(m, [10.0, 10.0, 0.0], 2.0, 10.0),
            gusset(m, SUNK, 30.0),
        )
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// The control: a gusset on `z = y/4 − 5`, which meets the bore's cylinder below the plate, so the
/// plane misses every lateral face. Corners `(24, 1)`, `(36, 1)`, `(36, 4)` — area 18, wholly
/// inside the plate.
#[test]
fn a_gusset_whose_plane_misses_the_bore_builds() {
    let build = |m: &mut Model| {
        (
            plate_with_bore(m),
            gusset(m, [(24.0, 1.0), (36.0, 1.0), (36.0, 4.0)], 30.0),
        )
    };
    gate(build, Ok(0));
    every_kind(
        build,
        [
            [Built(1, bored()), Built(1, bored() - 18.0), Built(1, 18.0)],
            [Built(1, bored()), Built(0, 0.0), Built(1, 18.0)],
        ],
    );
}

/// The bored plate with a pin along `x` fused on its top — radius 1 about `y = 32, z = 10.5`,
/// `x ∈ [2, 20]` — whose axis is parallel to the slanted plane, which runs through it: the class
/// carries the pin's rulings beside the bore's section. Above the top the pin keeps the disk less
/// the segment below `z = 10` (distance ½ from the centre): `2π/3 + √3/4` per unit, over 18.
#[test]
fn a_gusset_sunk_into_a_bored_plate_with_a_pin_along_its_plane() {
    let build = |m: &mut Model| {
        let p = plate_with_bore(m);
        let pin = cylinder(
            m,
            [2.0, 32.0, 10.5],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            1.0,
            18.0,
        );
        (op(m, BoolKind::Fuse, p, pin), gusset(m, SUNK, 30.0))
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// A steeper plane, `z = 3y − 25`, meets the bore across both of the plate's faces: the section is
/// arcs of the ellipse, not a whole one. Corners `(10, 5)`, `(12, 5)`, `(12, 11)`.
#[test]
fn a_gusset_whose_plane_crosses_a_bores_rims_is_refused() {
    let build = |m: &mut Model| {
        (
            plate_with_bore(m),
            gusset(m, [(10.0, 5.0), (12.0, 5.0), (12.0, 11.0)], 30.0),
        )
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// The fillet plate (one fillet of radius 5 about `(40, −20)`, `z ∈ [0, 12]`) and a gusset at
/// `x ∈ [0, 1]` on `z = y/4 + 11`, which meets the fillet's quarter at `z ∈ [4.75, 6]`: the
/// fillet's face is no band — rulings bound it — so the section is an arc of the ellipse.
#[test]
fn a_gusset_whose_plane_cuts_a_fillet_is_refused() {
    let build = |m: &mut Model| {
        (
            rounded_plate(m, 1, 0),
            gusset(m, [(-4.0, 10.0), (8.0, 10.0), (8.0, 13.0)], 0.0),
        )
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// The sunk gusset fused with a block `[11, 31]×[9, 15]×[4, 6.5]` whose wall `y = 9` cuts the bore:
/// that wall's trace on the slanted class crosses the bore's section there.
#[test]
fn a_gusset_whose_other_faces_cut_the_bore_is_refused() {
    let build = |m: &mut Model| {
        let g = gusset(m, SUNK, 30.0);
        let block = crate::fixtures::cuboid(
            m,
            Point3::from_array([11.0, 9.0, 4.0]),
            Point3::from_array([31.0, 15.0, 6.5]),
        );
        m.rebuild_adjacency();
        (plate_with_bore(m), op(m, BoolKind::Fuse, g, block))
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// The bored plate with a second bore along `y` — radius 1 about `x = 20, z = 5` — which the
/// slanted plane meets at `y ∈ [6, 14]`, between the plate's walls: the class would carry sections
/// of two cylinders whose axes are not parallel.
#[test]
fn a_gusset_whose_plane_sections_two_crossing_bores_is_refused() {
    let build = |m: &mut Model| {
        let p = plate_with_bore(m);
        let across = cylinder(
            m,
            [20.0, -1.0, 5.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            1.0,
            42.0,
        );
        (op(m, BoolKind::Cut, p, across), gusset(m, SUNK, 30.0))
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// A gusset straddling the bore — `x ∈ [9.5, 10.5]`, corners `(6, 4)`, `(14, 4)`, `(14, 6)` on the
/// same plane: its slanted face crosses the lateral, the cut the chart does not take.
#[test]
fn a_gusset_across_the_bore_is_refused() {
    let build = |m: &mut Model| {
        (
            plate_with_bore(m),
            gusset(m, [(6.0, 4.0), (14.0, 4.0), (14.0, 6.0)], 9.5),
        )
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}
