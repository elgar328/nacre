//! **Named plane classes oblique to a cylinder** — the scoreboard.
//!
//! A gusset's slanted face lies on a plane oblique to a bore's axis: the plane meets the cylinder
//! in an ellipse. Each family pairs such a class with world cylinders (a bore, two bores, a tube, a
//! pin, a fillet, a hole with filleted corners) and pins two answers: what the gate says, read from
//! [`plane_index_setup`] directly — how many (cylinder, class) pairs it records as the class's
//! section of the cylinder — and what the boolean answers, by kind and in both operand orders
//! ([`every_kind`]). Volumes are derived, not copied. The section is a whole ellipse or arcs of one
//! — cut by a rim, by a fillet's rulings, by another face of the class — and either way it closes
//! the class's own arrangement; what stays refused is a face of the class reaching a lateral face,
//! and a class sectioning two cylinders whose axes are not parallel.
//!
//! The geometry, shared by the families: the bored plate of [`plate_with_bore`] (`[0,40]²×[0,10]`,
//! a bore of radius 2 about `(10, 10)`), and a gusset — a triangle in the `yz` plane extruded one
//! unit along `x` — whose slanted face lies on `z = y/4 + 5/2`. That plane meets the bore at
//! `z ∈ [4.5, 5.5]`, strictly between the plate's faces, while the gusset stands at `x ∈ [30, 31]`,
//! far from every lateral face.

use super::*;
use Want::{Built, Refused};

pub(crate) type Sketch = [(f64, f64); 3];

/// The gusset on the slanted plane: `(y, z)` corners `(14, 6)`, `(34, 6)`, `(34, 11)` — area 50.
/// In the plate it sinks to the top face `z = 10`, which it crosses at `y = 30`: 48 below, 2 above.
pub(crate) const SUNK: Sketch = [(14.0, 6.0), (34.0, 6.0), (34.0, 11.0)];

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

/// The plate `[0,40]²×[0,10]` with a square hole `[15, 25]²` whose corners are filleted, radius 2
/// — four quarter cylinders bounded by rulings, each meeting a flat wall of the hole tangentially.
fn filleted_hole_plate(m: &mut Model) -> Handle<Solid> {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let edges = vec![
        line(p2(0.0, 0.0), p2(40.0, 0.0)),
        line(p2(40.0, 0.0), p2(40.0, 40.0)),
        line(p2(40.0, 40.0), p2(0.0, 40.0)),
        line(p2(0.0, 40.0), p2(0.0, 0.0)),
        line(p2(17.0, 15.0), p2(23.0, 15.0)),
        arc_turns(p2(23.0, 17.0), p2(23.0, 15.0), 1),
        line(p2(25.0, 17.0), p2(25.0, 23.0)),
        arc_turns(p2(23.0, 23.0), p2(25.0, 23.0), 1),
        line(p2(23.0, 25.0), p2(17.0, 25.0)),
        arc_turns(p2(17.0, 23.0), p2(17.0, 25.0), 1),
        line(p2(15.0, 23.0), p2(15.0, 17.0)),
        arc_turns(p2(17.0, 17.0), p2(15.0, 17.0), 1),
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
    .expect("the plate extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// A prism over the world-`XY` polygon `pts` (counter-clockwise), `z ∈ [z0, z0 + h]`.
fn block_over(m: &mut Model, pts: &[(f64, f64)], z0: f64, h: f64) -> Handle<Solid> {
    let p2 = |(x, y): (f64, f64)| nacre_math::Point2::from_array([x, y]);
    let edges = (0..pts.len())
        .map(|i| line(p2(pts[i]), p2(pts[(i + 1) % pts.len()])))
        .collect();
    let profile = stated(edges).unwrap().remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("the block extrudes") else {
        unreachable!()
    };
    let q = |x: f64| nacre_exact::Rat::from_decimal(x).expect("a short decimal");
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid,
            isometry: nacre_exact::Isometry::translation([q(0.0), q(0.0), q(z0)]),
        },
    )
    .expect("the block moves") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

fn op(m: &mut Model, kind: BoolKind, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
    let out = boolean(m, kind, a, b).expect("the fixture's boolean builds")[0];
    m.rebuild_adjacency();
    out
}

/// A triangle in the `yz` plane, `(y, z)` corners counter-clockwise, extruded over
/// `x ∈ [x0, x0 + 1]`.
pub(crate) fn gusset(m: &mut Model, tri: Sketch, x0: f64) -> Handle<Solid> {
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

/// [`gusset`] with its coordinates permuted (`X = z`, `Y = x`, `Z = y`): a triangle in the `zx`
/// plane, `(z, x)` corners counter-clockwise, extruded over `y ∈ [y0, y0 + 1]`.
fn gusset_zx(m: &mut Model, tri: Sketch, y0: f64) -> Handle<Solid> {
    let p2 = |(u, v): (f64, f64)| nacre_math::Point2::from_array([u, v]);
    let edges = vec![
        line(p2(tri[0]), p2(tri[1])),
        line(p2(tri[1]), p2(tri[2])),
        line(p2(tri[2]), p2(tri[0])),
    ];
    let profile = stated(edges).unwrap().remove(0);
    let frame = SketchFrame::world(m, Axis::Y);
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
            isometry: nacre_exact::Isometry::translation([q(0.0), q(y0), q(0.0)]),
        },
    )
    .expect("the gusset moves") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// Whether the gate recorded `r` as a class's section of a cylinder drawn on the class's own
/// arrangement.
fn states_a_section(r: crate::planes::CylRecord) -> bool {
    r == crate::planes::CylRecord::Sections
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

/// What a gusset of [`SUNK`] makes with an operand `a` it sinks into by 48 of its 50 ([`sunk_by`]).
fn sunk(a: f64) -> [[Want; 3]; 2] {
    sunk_by(a, 48.0, 2.0)
}

/// What a gusset makes with an operand `a` it sinks into by `below`, standing `above` out of it:
/// Fuse adds the part above, Cut takes the part below, Common is that part — and the other way
/// round, the gusset less the operand is the part above.
fn sunk_by(a: f64, below: f64, above: f64) -> [[Want; 3]; 2] {
    [
        [Built(1, a + above), Built(1, a - below), Built(1, below)],
        [Built(1, a + above), Built(1, above), Built(1, below)],
    ]
}

/// What a gusset of area 50 apart from an operand `a` makes with it: two bodies, then each alone.
fn apart(a: f64) -> [[Want; 3]; 2] {
    [
        [Built(2, a + 50.0), Built(1, a), Built(0, 0.0)],
        [Built(2, a + 50.0), Built(1, 50.0), Built(0, 0.0)],
    ]
}

/// The sunk gusset beside the bored plate: the class's section of the bore is a whole ellipse
/// between the plate's faces, drawn on the class's own arrangement.
#[test]
fn a_gusset_sunk_into_a_bored_plate() {
    let build = |m: &mut Model| (plate_with_bore(m), gusset(m, SUNK, 30.0));
    gate(build, Ok(1));
    every_kind(build, sunk(bored()));
}

/// The same plane, the gusset clear of the plate.
#[test]
fn a_gusset_apart_from_a_bored_plate() {
    let build = |m: &mut Model| (plate_with_bore(m), gusset(m, APART, 30.0));
    gate(build, Ok(1));
    every_kind(build, apart(bored()));
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
    gate(build, Ok(2));
    every_kind(build, sunk(16000.0 - 80.0 * pi()));
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
    gate(build, Ok(2));
    every_kind(build, apart(30.0 * pi()));
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
    gate(build, Ok(1));
    every_kind(build, apart(40.0 * pi()));
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
    gate(build, Ok(1));
    every_kind(
        build,
        sunk(bored() + 18.0 * (2.0 * pi() / 3.0 + 3f64.sqrt() / 4.0)),
    );
}

/// A steeper plane, `z = 3y − 25`, meets the bore across both of the plate's faces (at `y = 25/3`
/// and `35/3`, within the bore's `y ∈ [8, 12]`): the section is two arcs of the ellipse, cut at the
/// rims, and the plate's two faces carry the chords between their ends. Corners `(10, 5)`,
/// `(12, 5)`, `(12, 11)` — area 6, of which the corner above `z = 10` (from `y = 35/3`) is `1/6`.
#[test]
fn a_gusset_whose_plane_crosses_a_bores_rims() {
    let build = |m: &mut Model| {
        (
            plate_with_bore(m),
            gusset(m, [(10.0, 5.0), (12.0, 5.0), (12.0, 11.0)], 30.0),
        )
    };
    gate(build, Ok(1));
    every_kind(build, sunk_by(bored(), 6.0 - 1.0 / 6.0, 1.0 / 6.0));
}

/// A plane `z = y/4 − 2` meets the bore's cylinder at `z ∈ [0, 1]`: it touches the plate's bottom
/// rim at the one point `(10, 8, 0)`, where the section and the bottom face's hole meet with no
/// crossing to name the point by. Refused at the gate, though the gusset — corners `(16, 2)`,
/// `(28, 2)`, `(28, 5)` — stands far from the bore.
#[test]
fn a_gusset_whose_plane_touches_a_bores_rim_is_refused() {
    let build = |m: &mut Model| {
        (
            plate_with_bore(m),
            gusset(m, [(16.0, 2.0), (28.0, 2.0), (28.0, 5.0)], 30.0),
        )
    };
    gate(build, Err(RejectReason::ObliqueCylinderCut));
    every_kind(build, REFUSED);
}

/// [`a_gusset_whose_other_faces_cut_the_bore`] with its coordinates permuted — `X = z`, `Y = x`,
/// `Z = y` — so the bore runs along `x`: the plate `[0,10]×[0,40]²` less a bore of radius 2 about
/// `(Y, Z) = (10, 10)`, the gusset on `X = Z/4 + 5/2` at `Y ∈ [30, 31]`, the block
/// `[4, 6.5]×[11, 31]×[9, 15]`. The volumes are the unpermuted row's. What the permutation moves is
/// the ring's lowest corner: the block's walls `Y = 11` and `Z = 9` cut an L into the section, and
/// the ring's lexicographically smallest node is the L's reflex corner, where its turn is the
/// wrong winding — so the ring is read at the section's extreme of `x`, which lies inside an arc
/// and needs the section's own coordinate gradient there (`g*/(n·m)`, the axis having an `x`
/// component), not the circle's.
#[test]
fn a_gusset_whose_other_faces_cut_a_bore_along_x() {
    let build = |m: &mut Model| {
        let plate = crate::fixtures::cuboid(
            m,
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([10.0, 40.0, 40.0]),
        );
        let bore = cylinder(
            m,
            [-1.0, 10.0, 10.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            2.0,
            12.0,
        );
        let plate = op(m, BoolKind::Cut, plate, bore);
        let g = gusset_zx(m, SUNK, 30.0);
        let block = crate::fixtures::cuboid(
            m,
            Point3::from_array([4.0, 11.0, 9.0]),
            Point3::from_array([6.5, 31.0, 15.0]),
        );
        m.rebuild_adjacency();
        (plate, op(m, BoolKind::Fuse, g, block))
    };
    gate(build, Ok(1));
    let outside = 2.0 + 2.5 * (pi() - 1.0);
    let inside = 349.875 - outside;
    every_kind(
        build,
        [
            [
                Built(1, bored() + outside),
                Built(1, bored() - inside),
                Built(1, inside),
            ],
            [
                Built(1, bored() + outside),
                Built(2, outside),
                Built(1, inside),
            ],
        ],
    );
}

/// The steep gusset of [`a_gusset_whose_plane_crosses_a_bores_rims`] fused with a block whose
/// wall `x + 2y = 137/4` cuts the bore near the top: over the polygon `(65/4, 9)`, `(31, 9)`,
/// `(31, 15)`, `(17/4, 15)` (area `124.5`), `z ∈ [8, 9.5]`. The wall's trace on the slanted class
/// `z = 3y − 25` crosses the bore's section at `z ≈ 9.27`, where the section is far above its
/// centre `z = 5` and the trace leans along it: there the arc bulges to the side of the trace its
/// radial says — `d · (N − c)⊥` — and not the side `d · (N − c)` says, which reads the height too.
/// The block meets the gusset in `9/8` (`y ∈ [11, 12]` above `z = 8` under the slant), so the
/// operand is `6 + 186.75 − 9/8`; outside the plate stand the gusset's `1/6` above the top and the
/// block's part in the bore — the disk of radius 2 beyond the wall, at distance `d = 4.25/√5` from
/// the axis, a segment `4·acos(d/2) − d·√(4 − d²)`, over 1.5.
#[test]
fn a_gusset_whose_other_face_leans_along_the_section_high_above_its_centre() {
    let build = |m: &mut Model| {
        let g = gusset(m, [(10.0, 5.0), (12.0, 5.0), (12.0, 11.0)], 30.0);
        let block = block_over(
            m,
            &[(16.25, 9.0), (31.0, 9.0), (31.0, 15.0), (4.25, 15.0)],
            8.0,
            1.5,
        );
        (plate_with_bore(m), op(m, BoolKind::Fuse, g, block))
    };
    gate(build, Ok(1));
    let d = 4.25 / 5f64.sqrt();
    let segment = 4.0 * (d / 2.0).acos() - d * (4.0 - d * d).sqrt();
    let operand = 6.0 + 124.5 * 1.5 - 1.125;
    let outside = 1.0 / 6.0 + 1.5 * segment;
    let inside = operand - outside;
    every_kind(
        build,
        [
            [
                Built(1, bored() + outside),
                Built(1, bored() - inside),
                Built(1, inside),
            ],
            [
                Built(1, bored() + outside),
                Built(2, outside),
                Built(1, inside),
            ],
        ],
    );
}

/// The plate with the filleted square hole and a gusset at `x ∈ [30, 31]` on `z = (y − 17)/4`, a
/// plane through the corners `(15, 17, 0)` and `(25, 17, 0)` where two fillets meet the hole's
/// walls `x = 15` and `x = 25` on the plate's bottom: the section of each of those fillets would
/// end at a point four surfaces share, which the trace does not name. Refused as the oblique
/// coverage's end, though the gusset — corners `(21, 1)`, `(33, 1)`, `(33, 4)` — stands far from
/// the hole.
#[test]
fn a_gusset_whose_plane_runs_through_a_fillets_corner_is_refused() {
    let build = |m: &mut Model| {
        (
            filleted_hole_plate(m),
            gusset(m, [(21.0, 1.0), (33.0, 1.0), (33.0, 4.0)], 30.0),
        )
    };
    gate(build, Ok(4));
    every_kind(build, REFUSED);
}

/// The fillet plate (one fillet of radius 5 about `(40, −20)`, `z ∈ [0, 12]`, volume
/// `54000 − 12(25 − 25π/4)`) and a gusset at `x ∈ [0, 1]` on `z = y/4 + 11`, which meets the
/// fillet's quarter at `z ∈ [4.75, 6]`: rulings bound the fillet's face, so the section is an arc
/// of the ellipse, ending on them. Corners `(−4, 10)`, `(8, 10)`, `(8, 13)` — area 18, of which the
/// corner above the plate's top `z = 12` (from `y = 4`) is 2.
#[test]
fn a_gusset_whose_plane_cuts_a_fillet() {
    let build = |m: &mut Model| {
        (
            rounded_plate(m, 1, 0),
            gusset(m, [(-4.0, 10.0), (8.0, 10.0), (8.0, 13.0)], 0.0),
        )
    };
    gate(build, Ok(1));
    let plate = 54000.0 - 12.0 * (25.0 - 6.25 * pi());
    every_kind(build, sunk_by(plate, 16.0, 2.0));
}

/// The fillet plate and a steep gusset at `x ∈ [0, 1]` on `z = 3y + 69`, which crosses the
/// fillet's bottom arc at `(44, −23, 0)` and its ruling at `(45, −20)` at `z = 9`: the section's arc
/// on the quarter runs from a crossing of the face's own arc to a crossing of its ruling. Corners
/// `(−22, 3)`, `(−18, 3)`, `(−18, 15)` — area 24, of which the corner above the top `z = 12` (from
/// `y = −19`) is `3/2`.
#[test]
fn a_gusset_whose_plane_crosses_a_fillets_arc() {
    let build = |m: &mut Model| {
        (
            rounded_plate(m, 1, 0),
            gusset(m, [(-22.0, 3.0), (-18.0, 3.0), (-18.0, 15.0)], 0.0),
        )
    };
    gate(build, Ok(1));
    let plate = 54000.0 - 12.0 * (25.0 - 6.25 * pi());
    every_kind(build, sunk_by(plate, 22.5, 1.5));
}

/// The plate `[0,40]²×[0,10]` with a square hole `[15, 25]²` whose four corners are filleted
/// (radius 2), and the sunk gusset beside it: the slanted plane meets the four fillets' quarters at
/// `z ∈ [6.25, 6.75]` and `[8.25, 8.75]`, each in an arc of its own ellipse ending on the fillet's
/// rulings, where the hole's flat walls take the section over. The hole is `100 − 4(4 − π)`.
#[test]
fn a_gusset_sunk_into_a_plate_with_a_filleted_hole() {
    let build = |m: &mut Model| (filleted_hole_plate(m), gusset(m, SUNK, 30.0));
    gate(build, Ok(4));
    every_kind(build, sunk(16000.0 - 10.0 * (84.0 + 4.0 * pi())));
}

/// The sunk gusset fused with a block `[11, 31]×[9, 15]×[4, 6.5]` whose walls `x = 11` and
/// `y = 9` cut the bore: their traces on the slanted class cut the bore's section into arcs. The
/// block (300) and the gusset (50) overlap in `1/8` (`y ∈ [14, 15]` under the slant), so the
/// operand is `349.875`; outside the plate stand the gusset's 2 above the top and the block's part
/// in the bore — the disk of radius 2 about `(10, 10)` beyond `x = 11` and `y = 9`, `π − 1`, over
/// 2.5 — two bodies the other way round.
#[test]
fn a_gusset_whose_other_faces_cut_the_bore() {
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
    gate(build, Ok(1));
    let outside = 2.0 + 2.5 * (pi() - 1.0);
    let inside = 349.875 - outside;
    every_kind(
        build,
        [
            [
                Built(1, bored() + outside),
                Built(1, bored() - inside),
                Built(1, inside),
            ],
            [
                Built(1, bored() + outside),
                Built(2, outside),
                Built(1, inside),
            ],
        ],
    );
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
