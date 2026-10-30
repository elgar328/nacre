//! **A cylinder over a plate's edge** — a pin or boss whose rim crosses the edge of the face it
//! stands on. Where the rim meets the edge's wall the corner lies on two planes only (the face and
//! the overhang's underside are one plane) and on the cylinder: a pierce corner a loop visits
//! between two faces of one plane, and a rim the corners split into arcs. Each body is taken
//! into a second boolean, every kind, both orders, against derived volumes; the bodies are built by
//! the same booleans a user script runs.
//!
//! And the property the families rest on: every loop that visits a vertex names it alike
//! ([`one_name_per_vertex`]).

use super::*;
use nacre_exact::{Isometry, Rat};

use Want::{Built, Refused};

fn cuboid(m: &mut Model, min: [f64; 3], max: [f64; 3]) -> Handle<Solid> {
    let s = crate::fixtures::cuboid(m, Point3::from_array(min), Point3::from_array(max));
    m.rebuild_adjacency();
    s
}

/// A cylinder along the world axis `dir` of radius `r` standing on `at` with height `h`.
fn cylinder(m: &mut Model, at: [f64; 3], dir: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
    let ref_dir = if dir[2] != 0.0 {
        [0.0, -1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let s = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array(at),
        Vector3::from_array(dir),
        Vector3::from_array(ref_dir),
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

/// A prism over the world triangle `tri`, from `z0` up `h`.
fn prism(m: &mut Model, tri: &[[f64; 2]; 3], z0: i128, h: f64) -> Handle<Solid> {
    let profile = Profile2d::polygon(tri.iter().map(|q| p2(q[0], q[1])).collect()).unwrap();
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("the prism extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let out = transform(
        m,
        solid,
        &Isometry::translation([0, 0, z0].map(Rat::from_int)),
    )
    .expect("the prism moves");
    m.rebuild_adjacency();
    out
}

const PI: f64 = std::f64::consts::PI;

/// The area of the disk of radius `r` beyond a chord at distance `d` from its centre.
fn segment(r: f64, d: f64) -> f64 {
    r * r * (d / r).acos() - d * (r * r - d * d).sqrt()
}

/// `∫ √(r² − w²) dw` from `0` to `w`.
fn quarter_integral(r: f64, w: f64) -> f64 {
    (w * (r * r - w * w).max(0.0).sqrt() + r * r * (w / r).asin()) / 2.0
}

/// The area of the triangle `tri` inside the box `[lo, hi]` — clipped edge by edge.
fn triangle_in_box(tri: &[[f64; 2]; 3], lo: [f64; 2], hi: [f64; 2]) -> f64 {
    let clip = |p: &[[f64; 2]], n: [f64; 2], c: f64| {
        let f = |q: [f64; 2]| n[0] * q[0] + n[1] * q[1] - c;
        let mut out = Vec::new();
        for i in 0..p.len() {
            let (a, b) = (p[i], p[(i + 1) % p.len()]);
            let (fa, fb) = (f(a), f(b));
            if fa <= 0.0 {
                out.push(a);
            }
            if (fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0) {
                let t = fa / (fa - fb);
                out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            }
        }
        out
    };
    let p = clip(tri, [1.0, 0.0], hi[0]);
    let p = clip(&p, [-1.0, 0.0], -lo[0]);
    let p = clip(&p, [0.0, 1.0], hi[1]);
    let p = clip(&p, [0.0, -1.0], -lo[1]);
    let n = p.len();
    (0..n)
        .map(|i| p[i][0] * p[(i + 1) % n][1] - p[i][1] * p[(i + 1) % n][0])
        .sum::<f64>()
        .abs()
        / 2.0
}

/// The cells of a second operand of volume `vb` meeting a body of volume `va` in `common`, every
/// kind as built then swapped, one body each.
fn cells(va: f64, vb: f64, common: f64) -> [[Want; 3]; 2] {
    [
        [
            Built(1, va + vb - common),
            Built(1, va - common),
            Built(1, common),
        ],
        [
            Built(1, va + vb - common),
            Built(1, vb - common),
            Built(1, common),
        ],
    ]
}

/// A `5×5×10` box standing on `z = 5` far from every edge: half of it, `125`, inside the plate.
fn far_box(m: &mut Model) -> Handle<Solid> {
    cuboid(m, [22.5, 5.5, 5.0], [27.5, 10.5, 15.0])
}

/// A pin of radius 2 and height 3 meets a box standing past a line across it — the line at
/// `edge`, the axis at `axis`, both measured across the pin — in the part of its disk beyond.
fn beyond(edge: f64, axis: f64) -> f64 {
    3.0 * segment(2.0, edge - axis)
}

/// The pin's corner on the edge's wall: its definition and the planes of the faces at it — the
/// witness that a family reaches a pierce corner on two planes.
fn assert_two_plane_pierce(m: &Model, s: Handle<Solid>, at: [f64; 3]) {
    let mut found = None;
    for &f in &m.shell(m.solid(s).outer).faces {
        let face = m.face(f);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                for v in m.edge(he.edge).vertices {
                    let p = m.vertex_point(v).as_array();
                    if (0..3).all(|i| (p[i] - at[i]).abs() < 1e-9) {
                        found = Some(v);
                    }
                }
            }
        }
    }
    let v = found.expect("a vertex at the corner");
    assert!(
        matches!(m.vertex(v), Vertex::Pierce { .. }),
        "the corner is a pierce point: {:?}",
        m.vertex(v)
    );
    let mut planes: Vec<Handle<nacre_topo::Surface>> = Vec::new();
    for &f in &m.shell(m.solid(s).outer).faces {
        let face = m.face(f);
        let at_v = std::iter::once(&face.outer)
            .chain(face.inner.iter())
            .flat_map(|lp| lp.half_edges.iter())
            .any(|he| m.edge(he.edge).vertices.contains(&v));
        if at_v && matches!(m.surface(face.surface), nacre_topo::Surface::Plane { .. }) {
            let name = m.world_plane_name(face.surface);
            if !planes.iter().any(|&q| m.world_plane_name(q) == name) {
                planes.push(face.surface);
            }
        }
    }
    assert_eq!(planes.len(), 2, "two planes at the corner");
}

#[test]
fn a_pin_over_a_plates_edge_takes_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let pin = |m: &mut Model| plate_with_pin(m, [39.0, 20.0, 10.0], 2.0, 3.0);
    let mut m = Model::new();
    let s = pin(&mut m);
    assert_two_plane_pierce(&m, s, [40.0, 20.0 + 3f64.sqrt(), 10.0]);
    every_kind(|m| (pin(m), far_box(m)), cells(pinned, 250.0, 125.0));
    // A box standing beyond the edge meets the overhang only.
    every_kind(
        |m| (pin(m), cuboid(m, [40.0, 17.0, 8.0], [45.0, 23.0, 14.0])),
        cells(pinned, 180.0, beyond(40.0, 39.0)),
    );
    every_kind(
        |m| (pin(m), cuboid(m, [40.5, 17.0, 8.0], [45.0, 23.0, 14.0])),
        cells(pinned, 162.0, beyond(40.5, 39.0)),
    );
}

/// Under the plate, the pin hangs from its bottom face.
#[test]
fn a_pin_under_a_plates_edge_takes_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let pin = |m: &mut Model| plate_with_pin(m, [39.0, 20.0, -3.0], 2.0, 3.0);
    every_kind(|m| (pin(m), far_box(m)), cells(pinned, 250.0, 125.0));
    every_kind(
        |m| (pin(m), cuboid(m, [40.0, 17.0, -4.0], [45.0, 23.0, 2.0])),
        cells(pinned, 180.0, beyond(40.0, 39.0)),
    );
}

/// Over the corner, the rim crosses both edges.
#[test]
fn a_pin_over_a_plates_corner_takes_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let pin = |m: &mut Model| plate_with_pin(m, [39.0, 39.0, 10.0], 2.0, 3.0);
    every_kind(|m| (pin(m), far_box(m)), cells(pinned, 250.0, 125.0));
    every_kind(
        |m| (pin(m), cuboid(m, [40.0, 36.0, 8.0], [45.0, 42.0, 14.0])),
        cells(pinned, 180.0, beyond(40.0, 39.0)),
    );
}

/// The axis on the edge, half the pin beyond it; and the axis past it, most of the pin beyond.
#[test]
fn a_pin_half_and_mostly_over_the_edge_takes_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let half = |m: &mut Model| plate_with_pin(m, [40.0, 20.0, 10.0], 2.0, 3.0);
    every_kind(|m| (half(m), far_box(m)), cells(pinned, 250.0, 125.0));
    every_kind(
        |m| (half(m), cuboid(m, [40.0, 17.0, 8.0], [45.0, 23.0, 14.0])),
        cells(pinned, 180.0, 6.0 * PI),
    );
    let most = |m: &mut Model| plate_with_pin(m, [41.0, 20.0, 10.0], 2.0, 3.0);
    every_kind(|m| (most(m), far_box(m)), cells(pinned, 250.0, 125.0));
    every_kind(
        |m| (most(m), cuboid(m, [40.0, 17.0, 8.0], [45.0, 23.0, 14.0])),
        cells(pinned, 180.0, 3.0 * (4.0 * PI - segment(2.0, 1.0))),
    );
}

/// On a step: the plate with `[0,36]×[0,40]×[10,12]` fused on, `18880`, and the pin standing on
/// the step over its edge `x = 36`.
#[test]
fn a_pin_over_a_steps_edge_takes_a_second_boolean() {
    let stepped = 18880.0 + 12.0 * PI;
    let pin = |m: &mut Model| {
        let p = cuboid(m, [0.0; 3], [40.0, 40.0, 10.0]);
        let step = cuboid(m, [0.0, 0.0, 10.0], [36.0, 40.0, 12.0]);
        let p = op(m, BoolKind::Fuse, p, step);
        let c = cylinder(m, [35.0, 20.0, 12.0], [0.0, 0.0, 1.0], 2.0, 3.0);
        op(m, BoolKind::Fuse, p, c)
    };
    // The far box stands where the step is: `5×5×7` of it inside.
    every_kind(|m| (pin(m), far_box(m)), cells(stepped, 250.0, 175.0));
    every_kind(
        |m| (pin(m), cuboid(m, [36.0, 17.0, 11.0], [41.0, 23.0, 16.0])),
        cells(stepped, 150.0, beyond(36.0, 35.0)),
    );
}

/// A pin lying on the wall `x = 40`, its axis along `x` at `z = 9`, its rim crossing the top edge.
/// The box standing on the plate's top beyond the wall meets the plate along the edge line only,
/// so their union is refused by name.
#[test]
fn a_pin_on_a_wall_over_its_top_edge_takes_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let pin = |m: &mut Model| {
        let p = cuboid(m, [0.0; 3], [40.0, 40.0, 10.0]);
        let c = cylinder(m, [40.0, 20.0, 9.0], [1.0, 0.0, 0.0], 2.0, 3.0);
        op(m, BoolKind::Fuse, p, c)
    };
    every_kind(|m| (pin(m), far_box(m)), cells(pinned, 250.0, 125.0));
    let touch = beyond(10.0, 9.0);
    let mut want = cells(pinned, 120.0, touch);
    want[0][0] = Refused(RejectReason::NonManifoldResultEdge);
    want[1][0] = Refused(RejectReason::NonManifoldResultEdge);
    every_kind(
        |m| (pin(m), cuboid(m, [40.0, 17.0, 10.0], [45.0, 23.0, 14.0])),
        want,
    );
}

/// The pinned plate with a bore of radius 1 cut far from the pin — the cut is itself a second
/// boolean, and the body carries two cylinders.
#[test]
fn a_pin_over_the_edge_beside_a_bore_takes_a_second_boolean() {
    let body = 16000.0 + 12.0 * PI - 10.0 * PI;
    let pin = |m: &mut Model| {
        let p = plate_with_pin(m, [39.0, 20.0, 10.0], 2.0, 3.0);
        let b = cylinder(m, [10.0, 10.0, -1.0], [0.0, 0.0, 1.0], 1.0, 12.0);
        op(m, BoolKind::Cut, p, b)
    };
    every_kind(|m| (pin(m), far_box(m)), cells(body, 250.0, 125.0));
    every_kind(
        |m| (pin(m), cuboid(m, [40.0, 17.0, 8.0], [45.0, 23.0, 14.0])),
        cells(body, 180.0, beyond(40.0, 39.0)),
    );
}

/// A box fused on the pin's cap and past its rim: the box's wall `x = 21` meets the rim on the
/// cap's plane, where the cap and the box's underside are one class — the same corner on two
/// planes, at the top rim. The box `[21,25]×[18,22]×[13,15]` stands on the cap, `32`, beside
/// nothing else.
#[test]
fn a_box_on_a_pins_cap_past_its_rim_takes_a_second_boolean() {
    let body = 16000.0 + 12.0 * PI + 32.0;
    let capped = |m: &mut Model| {
        let p = plate_with_pin(m, [20.0, 20.0, 10.0], 2.0, 3.0);
        let b = cuboid(m, [21.0, 18.0, 13.0], [25.0, 22.0, 15.0]);
        op(m, BoolKind::Fuse, p, b)
    };
    let mut m = Model::new();
    let s = capped(&mut m);
    assert_two_plane_pierce(&m, s, [21.0, 20.0 + 3f64.sqrt(), 13.0]);
    every_kind(|m| (capped(m), far_box(m)), cells(body, 250.0, 125.0));
    // Under the box past the rim: `[23,25]×[19,21]×[13,14]` of it.
    every_kind(
        |m| (capped(m), cuboid(m, [23.0, 19.0, 12.0], [26.0, 21.0, 14.0])),
        cells(body, 12.0, 4.0),
    );
}

/// The control beside the edge: a pin tangent to the edge `x = 40` from inside keeps its rim one
/// closed edge — no corner on the wall.
#[test]
fn a_pin_tangent_to_the_edge_takes_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let pin = |m: &mut Model| plate_with_pin(m, [38.0, 20.0, 10.0], 2.0, 3.0);
    every_kind(|m| (pin(m), far_box(m)), cells(pinned, 250.0, 125.0));
    // The box `[39,45]×[17,23]×[8,14]`: `[39,40]×[17,23]×[8,10]` of the plate and the pin's disk
    // past `x = 39`.
    every_kind(
        |m| (pin(m), cuboid(m, [39.0, 17.0, 8.0], [45.0, 23.0, 14.0])),
        cells(pinned, 216.0, 12.0 + beyond(39.0, 38.0)),
    );
}

/// A pin of radius 5 at `(37, 36)`, its rim through the plate's corner `(40, 40, 10)` — a corner on
/// three planes and the cylinder — and over both edges beside it.
#[test]
fn a_rim_through_a_plates_corner_takes_a_second_boolean() {
    let pinned = 16000.0 + 75.0 * PI;
    let pin = |m: &mut Model| plate_with_pin(m, [37.0, 36.0, 10.0], 5.0, 3.0);
    every_kind(|m| (pin(m), far_box(m)), cells(pinned, 250.0, 125.0));
    // The box over the corner: `[38,40]²` of the plate under `z = 10`, and the pin's disk in
    // `x, y ≥ 38` — `u = x − 37` from `1` to `√21`, where the circle meets `y = 38`.
    let s21 = 21f64.sqrt();
    let corner = (quarter_integral(5.0, s21) - quarter_integral(5.0, 1.0)) - 2.0 * (s21 - 1.0);
    every_kind(
        |m| (pin(m), cuboid(m, [38.0, 38.0, 8.0], [45.0, 45.0, 14.0])),
        cells(pinned, 294.0, 8.0 + 3.0 * corner),
    );
}

/// The rounded plate of [`rounded_plate`] (one fillet of radius 5 about `(40, −20)`) with a prism
/// on its top whose corner is the fillet's tangent point `(40, −25, 12)` — four planes and the
/// cylinder there: the top, the plate's wall `y = −25`, and the prism's two walls.
fn fillet_plate_with_prism(
    tri: &[[f64; 2]; 3],
    area: f64,
) -> (impl Fn(&mut Model) -> Handle<Solid>, f64) {
    let tri = *tri;
    let plate = 54000.0 - 12.0 * (25.0 - 25.0 * PI / 4.0);
    (
        move |m: &mut Model| {
            let p = rounded_plate(m, 1, 0);
            m.rebuild_adjacency();
            let t = prism(m, &tri, 12, 2.0);
            op(m, BoolKind::Fuse, p, t)
        },
        plate + 2.0 * area,
    )
}

/// The box across the tangent point, `[38,42]×[−27,−21]×[11,14]`: the plate's part under
/// `z = 12` — `[38,40]×[−25,−21]` and the fillet's, `v = y + 20 ∈ [−5, −1]`, `u ∈ [0, min(2,
/// √(25 − v²))]` — and the prism's part above.
fn across_the_tangent_point(tri: &[[f64; 2]; 3]) -> f64 {
    let s21 = 21f64.sqrt();
    let fillet = 2.0 * (s21 - 1.0) + (quarter_integral(5.0, 5.0) - quarter_integral(5.0, s21));
    8.0 + fillet + 2.0 * triangle_in_box(tri, [38.0, -27.0], [42.0, -21.0])
}

#[test]
fn four_planes_at_a_fillets_tangent_point_take_a_second_boolean() {
    for (tri, area, box_less_body) in [
        // Two oblique walls.
        ([[40.0, -25.0], [38.0, -21.0], [36.0, -24.0]], 7.0, 1),
        // The wall `x = 40` and an oblique one, `y = x − 60`: above the plate that wall cuts the
        // box's corner `(38, −21)` off — a triangle of area ½, 2 high, a body of its own in what
        // the body leaves of the box.
        ([[40.0, -25.0], [40.0, -20.0], [36.0, -24.0]], 10.0, 2),
    ] {
        let (body, v) = fillet_plate_with_prism(&tri, area);
        every_kind(
            |m| (body(m), cuboid(m, [-10.0, -5.0, 5.0], [-5.0, 0.0, 15.0])),
            cells(v, 250.0, 175.0),
        );
        let mut want = cells(v, 72.0, across_the_tangent_point(&tri));
        if let Built(_, vol) = want[1][1] {
            want[1][1] = Built(box_less_body, vol);
        }
        every_kind(
            |m| (body(m), cuboid(m, [38.0, -27.0, 11.0], [42.0, -21.0, 14.0])),
            want,
        );
    }
}

/// The controls: a pin inside the outline, and a bore through the edge (a notch, whose wall's
/// loop turns onto a ruling at the corner).
#[test]
fn an_inner_pin_and_a_notch_take_a_second_boolean() {
    let pinned = 16000.0 + 12.0 * PI;
    let inner = |m: &mut Model| plate_with_pin(m, [10.0, 10.0, 10.0], 2.0, 3.0);
    every_kind(|m| (inner(m), far_box(m)), cells(pinned, 250.0, 125.0));
    every_kind(
        |m| (inner(m), cuboid(m, [10.0, 7.0, 8.0], [15.0, 13.0, 14.0])),
        cells(pinned, 180.0, 60.0 + 6.0 * PI),
    );
    let s3 = 3f64.sqrt();
    let notched = 16000.0 - 10.0 * (8.0 * PI / 3.0 + s3);
    let notch = |m: &mut Model| {
        let p = cuboid(m, [0.0; 3], [40.0, 40.0, 10.0]);
        let c = cylinder(m, [39.0, 20.0, -1.0], [0.0, 0.0, 1.0], 2.0, 12.0);
        op(m, BoolKind::Cut, p, c)
    };
    every_kind(|m| (notch(m), far_box(m)), cells(notched, 250.0, 125.0));
    // The box `[38,45]×[17,23]×[5,15]` holds `[38,40]×[17,23]×[5,10]` of the plate less the
    // bore's strip `|x − 39| ≤ 1`, `2√3 + 4π/3` deep 5 — which cuts that block in two.
    let common = 60.0 - 5.0 * (4.0 * PI / 3.0 + 2.0 * s3);
    let mut want = cells(notched, 420.0, common);
    want[0][2] = Built(2, common);
    want[1][2] = Built(2, common);
    every_kind(
        |m| (notch(m), cuboid(m, [38.0, 17.0, 5.0], [45.0, 23.0, 15.0])),
        want,
    );
}

/// **Every loop that visits a vertex names it alike** — the plate faces', the wall's and the
/// lateral's, on bodies whose corners lie on a cylinder and two, three or four planes: a pin over
/// an edge and over a corner, a rim through the plate's corner, a fillet's tangent point under a
/// prism, and a cap pinched where a wedge's apex touches its rim. A pierce corner a loop passes
/// between two faces of one plane has no three-plane name at all; one on three planes and the
/// cylinder would take a three-plane name on the visits between plane faces and its pierce name
/// on the lateral's.
#[test]
fn one_name_per_vertex() {
    type Build = fn(&mut Model) -> Handle<Solid>;
    let bodies: Vec<(&str, Build)> = vec![
        ("pin over the edge", |m| {
            plate_with_pin(m, [39.0, 20.0, 10.0], 2.0, 3.0)
        }),
        ("pin over the corner", |m| {
            plate_with_pin(m, [39.0, 39.0, 10.0], 2.0, 3.0)
        }),
        ("rim through the plate's corner", |m| {
            plate_with_pin(m, [37.0, 36.0, 10.0], 5.0, 3.0)
        }),
        ("four planes at the tangent point", |m| {
            let p = rounded_plate(m, 1, 0);
            m.rebuild_adjacency();
            let t = prism(m, &[[40.0, -25.0], [38.0, -21.0], [36.0, -24.0]], 12, 2.0);
            op(m, BoolKind::Fuse, p, t)
        }),
        ("pinched cap", |m| {
            let w =
                crate::fixtures::wedge_on_a_rim(m, [[1.0, 0.0], [0.5, 0.3], [0.5, -0.3]], 0.5, 2.0);
            m.rebuild_adjacency();
            w
        }),
    ];
    for (name, build) in bodies {
        let mut m = Model::new();
        let s = build(&mut m);
        let names = vertex_names(&m, s);
        assert!(
            names
                .values()
                .any(|ns| ns.iter().any(|n| combinatorics::pierce_name(*n).is_some())),
            "{name}: the body has pierce corners"
        );
        for (v, ns) in &names {
            assert!(
                ns.windows(2).all(|w| w[0] == w[1]),
                "{name}: vertex {v:?} at {:?} is named {ns:?}",
                m.vertex_point(*v).as_array()
            );
        }
    }
}

/// Every name the tracer's loops give each vertex of `s`, one entry per visit — plane faces'
/// outer and hole loops ([`combinatorics::face_vertex_triples`], [`combinatorics::hole_rings`]) and
/// laterals' cycles ([`combinatorics::lateral_cycles`]), set up as the tracer is.
fn vertex_names(
    m: &Model,
    s: Handle<Solid>,
) -> HashMap<Handle<Vertex>, Vec<combinatorics::NodeId>> {
    let faces_tab = collect_planes(m, s).expect("the face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(m, s, &surf_ix).expect("edge incidence");
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = m.world_cylinder_def(surf).expect("a world cylinder");
            let Some(cache) = m.surface_cache(surf).cylinder() else {
                unreachable!()
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                realized: *cache,
                owner: crate::planes::SolidSide::A,
                shared: Vec::new(),
                grazed: Vec::new(),
                records: Vec::new(),
                boxes: Vec::new(),
            }
        })
        .collect();
    let jd = crate::planes::test_judge(&planes);
    let mut out: HashMap<Handle<Vertex>, Vec<combinatorics::NodeId>> = HashMap::new();
    for &fh in &m.shell(m.solid(s).outer).faces {
        let fp = surf_ix[&fh];
        let face = m.face(fh);
        let loops: Vec<&nacre_topo::Loop> = std::iter::once(&face.outer)
            .chain(face.inner.iter())
            .collect();
        let rings = if matches!(plane_ix[fp], crate::planes::ClassIx::Cyl(_)) {
            combinatorics::lateral_cycles(m, fh, fp, &inc, &jd, &plane_ix, &cyls)
        } else {
            combinatorics::face_vertex_triples(m, fh, fp, &inc, &jd, &plane_ix, &cyls).and_then(
                |outer| {
                    let mut all = vec![outer];
                    all.extend(combinatorics::hole_rings(
                        m, fh, fp, &inc, &jd, &plane_ix, &cyls,
                    )?);
                    Ok(all)
                },
            )
        }
        .unwrap_or_else(|e| panic!("a loop of face {fh:?} could not be named: {e:?}"));
        for (lp, ring) in loops.iter().zip(&rings) {
            let Some(nr) = ring.poly() else { continue };
            for (he, &n) in lp.half_edges.iter().zip(&nr.triples) {
                out.entry(crate::he_start(m, *he)).or_default().push(n);
            }
        }
    }
    out
}
