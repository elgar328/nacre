//! **A box joined to a cylinder along its tangent line** — the scoreboard of booleans on the union
//! `U` of the box `[−1, 1] × [0, 3] × [0, 1]` and the unit rod about the `z` axis over `z ∈ cz`.
//! The box's faces `x = ±1` touch the rod along the lines `(±1, 0)`, its edges `x = ±1, y = 0`
//! lie on those lines, and its face `y = 0` runs through the axis — three planes of one solid on
//! one line of its own lateral. Where the rod stands on the box's height the face `y = 0` is buried
//! and `U` is a slot's rounded end (it builds, and is no row here); staggered, a strip of `y = 0`
//! stays and holds the line:
//!
//! - `U9`, the rod over `z ∈ [1/2, 3/2]`: nine faces. The rod's lateral climbs each line from the
//!   half-disk under it, on the side away from the box, to the arc over the box's top on the box's
//!   side — the boundary runs along a line another plane (`x = ±1`) holds and leaves it on the
//!   other side;
//! - `U10`, the rod over `z ∈ [1/4, 3/4]`, inside the box's height: ten faces. The lateral is the
//!   half away from the box, and each line is its edge.
//!
//! Each also mirrored, the box over `y ∈ [−3, 0]`: its face `y = 0` then looks along `+y`, which on
//! the seed plane `y = 0` is the class's stored normal reversed (`frame_sign = −1`).
//!
//! Each against four tools, as Fuse, both Cuts and Common: a box far off; the box
//! `[1/2, 3/2] × [−1/2, 1/2] × [0, 2]` across the line `(1, 0)`; the slot
//! `[−2, 2] × [−1/4, 1/4] × [1/5, 9/10]` across both lines; and the rod's own cylinder over
//! `z ∈ [3/5, 6/5]`. The volumes are closed forms (`volumes`), and OCCT 7.9 builds every boolean
//! of the unmirrored box and slot with them (`the_junctions_closed_forms_are_occts`). Its face
//! counts are these but in two places: each union with the box across has two more, as the
//! unstaggered slot's union does (which builds here too); and `U9`'s box across less `U`, and their
//! common part, have one more.
//! In those two the box's top keeps a sliver between the rod's arc and the wall tangent to it,
//! whose mesh the plane faces' chord budget holds (`mesh_covers_faces`).
//!
//! The rod's own cylinder is refused: the tangent fold finds a tangent piece overlapping a station
//! without matching it (`cyl_chart::stations_of_lines`, todo «접선에서만 만나는 덩이 — 마디가
//! 어긋나는 칼날»).

use super::scoreboard::{OPS, Want, check, counted, refused, run};
use super::*;
use std::f64::consts::PI;

/// The rod's span: `U9`'s or `U10`'s.
const U9: [f64; 2] = [0.5, 1.5];
const U10: [f64; 2] = [0.25, 0.75];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tool {
    Far,
    Across,
    Slot,
    Coaxial,
}

const TOOLS: [Tool; 4] = [Tool::Far, Tool::Across, Tool::Slot, Tool::Coaxial];

/// The unit rod about the `z` axis over `z ∈ [z0, z1]`.
fn rod(m: &mut Model, z0: f64, z1: f64) -> Handle<Solid> {
    let frame = datum_frame(
        m,
        SketchPlane::from_axes(
            Point3::from_array([0.0, 0.0, z0]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
        ),
    );
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile: circle_profile([0.0, 0.0], 1.0),
            dist: z1 - z0,
        },
    )
    .expect("the rod extrudes") else {
        unreachable!()
    };
    solid
}

/// `U`: the box on `y ∈ [0, 3]` (or `[−3, 0]`, `mirrored`) fused with the rod over `cz`.
fn junction(m: &mut Model, cz: [f64; 2], mirrored: bool) -> Handle<Solid> {
    let (y0, y1) = if mirrored { (-3.0, 0.0) } else { (0.0, 3.0) };
    let block = crate::fixtures::cuboid(m, [-1.0, y0, 0.0], [1.0, y1, 1.0]);
    let r = rod(m, cz[0], cz[1]);
    let u = crate::boolean(m, BoolKind::Fuse, block, r).expect("the junction builds");
    assert_eq!(u.len(), 1, "the junction is one body");
    u[0]
}

fn tool(m: &mut Model, t: Tool) -> Handle<Solid> {
    match t {
        Tool::Far => crate::fixtures::cuboid(m, [5.0, 5.0, 0.0], [6.0, 6.0, 1.0]),
        Tool::Across => crate::fixtures::cuboid(m, [0.5, -0.5, 0.0], [1.5, 0.5, 2.0]),
        Tool::Slot => crate::fixtures::cuboid(m, [-2.0, -0.25, 0.2], [2.0, 0.25, 0.9]),
        Tool::Coaxial => rod(m, 0.6, 1.2),
    }
}

/// `∫ₐᵇ √(1 − y²) dy`.
fn chord_area(a: f64, b: f64) -> f64 {
    let f = |y: f64| 0.5 * (y * (1.0 - y * y).sqrt() + y.asin());
    f(b) - f(a)
}

/// The area a tool's section takes from `U`'s section at a height where the box is (`box_here`)
/// and the rod is (`rod_here`). Every tool is symmetric about `y = 0`, so the mirror takes the same.
fn common_section(t: Tool, box_here: bool, rod_here: bool) -> f64 {
    match (t, box_here, rod_here) {
        (_, false, false) | (Tool::Far, _, _) => 0.0,
        // `x ∈ [1/2, 1]` of the box's half, and the rod's half beyond the box out to `x = 1/2`.
        (Tool::Across, true, true) => 0.25 + (chord_area(-0.5, 0.0) - 0.25),
        (Tool::Across, true, false) => 0.25,
        (Tool::Across, false, true) => chord_area(-0.5, 0.5) - 0.5,
        // The strip `|y| ≤ 1/4` across: the box's half whole, the rod's half beyond it.
        (Tool::Slot, true, true) => 0.5 + 2.0 * chord_area(-0.25, 0.0),
        (Tool::Slot, true, false) => 0.5,
        (Tool::Slot, false, true) => 2.0 * chord_area(-0.25, 0.25),
        (Tool::Coaxial, _, true) => PI,
        // The box's half of the disk.
        (Tool::Coaxial, true, false) => PI / 2.0,
    }
}

/// The tool's own height span.
fn tool_span(t: Tool) -> [f64; 2] {
    match t {
        Tool::Far => [0.0, 1.0],
        Tool::Across => [0.0, 2.0],
        Tool::Slot => [0.2, 0.9],
        Tool::Coaxial => [0.6, 1.2],
    }
}

/// `(U, tool, U ∩ tool)` in closed form, slab by slab between the heights where something starts or
/// stops.
fn volumes(cz: [f64; 2], t: Tool) -> (f64, f64, f64) {
    let ts = tool_span(t);
    let mut zs = vec![0.0, 1.0, cz[0], cz[1], ts[0], ts[1]];
    zs.sort_by(f64::total_cmp);
    zs.dedup();
    let (mut u, mut c) = (6.0, 0.0);
    for w in zs.windows(2) {
        let (z0, z1) = (w[0], w[1]);
        let mid = 0.5 * (z0 + z1);
        let box_here = (0.0..=1.0).contains(&mid);
        let rod_here = (cz[0]..=cz[1]).contains(&mid);
        // The rod's half beyond the box, or all of it above the box.
        if rod_here {
            u += (z1 - z0) * if box_here { PI / 2.0 } else { PI };
        }
        if (ts[0]..=ts[1]).contains(&mid) {
            c += (z1 - z0) * common_section(t, box_here, rod_here);
        }
    }
    let tv = match t {
        Tool::Far => 1.0,
        Tool::Across => 2.0,
        Tool::Slot => 4.0 * 0.5 * 0.7,
        Tool::Coaxial => 0.6 * PI,
    };
    (u, tv, c)
}

fn wants(cz: [f64; 2], t: Tool) -> [Want; 4] {
    let (u, tv, c) = volumes(cz, t);
    // Faces per body: U9's, U10's.
    let f = |u9: usize, u10: usize| if cz == U9 { u9 } else { u10 };
    match t {
        Tool::Coaxial => core::array::from_fn(|_| refused(RejectReason::CylinderStagesDisagree)),
        // Far apart: the union is the two, each cut leaves its first operand, nothing in common.
        Tool::Far => [
            counted(&[tv, u], &[6, f(9, 10)]),
            counted(&[u], &[f(9, 10)]),
            counted(&[tv], &[6]),
            counted(&[], &[]),
        ],
        Tool::Across => [
            counted(&[u + tv - c], &[f(15, 16)]),
            counted(&[u - c], &[f(12, 13)]),
            counted(&[tv - c], &[f(13, 14)]),
            counted(&[c], &[f(10, 11)]),
        ],
        // In `U10` the slot parts the rod's half beyond `y = −1/4` from the rest: its segment of
        // the disk over the rod's half height, four faces.
        Tool::Slot if cz == U10 => {
            let part = chord_area(-1.0, -0.25);
            [
                counted(&[u + tv - c], &[17]),
                counted(&[part, u - c - part], &[4, 10]),
                counted(&[tv - c], &[15]),
                counted(&[c], &[12]),
            ]
        }
        Tool::Slot => [
            counted(&[u + tv - c], &[17]),
            counted(&[u - c], &[13]),
            counted(&[tv - c], &[14]),
            counted(&[c], &[10]),
        ],
    }
}

#[test]
fn a_box_joined_to_a_cylinder_along_its_tangent_line() {
    for cz in [U9, U10] {
        for mirrored in [false, true] {
            for t in TOOLS {
                for (op, want) in OPS.into_iter().zip(wants(cz, t)) {
                    let mut m = Model::new();
                    let u = junction(&mut m, cz, mirrored);
                    let b = tool(&mut m, t);
                    let got = run(&mut m, op, u, b);
                    let label = format!("rod {cz:?}, mirrored = {mirrored}, {t:?}, {op:?}");
                    check(&label, &m, got, &want);
                }
            }
        }
    }
}

/// The closed forms against the identities they must satisfy and the values OCCT 7.9 reports for
/// the unmirrored `U ∩ tool`.
#[test]
fn the_junctions_closed_forms_are_occts() {
    let occt = [
        (U9, Tool::Across, 0.592459),
        (U9, Tool::Slot, 0.547897),
        (U10, Tool::Across, 0.364153),
        (U10, Tool::Slot, 0.597371),
    ];
    for (cz, t, common) in occt {
        let (_, _, c) = volumes(cz, t);
        assert!(
            (c - common).abs() < 1e-6,
            "{cz:?} {t:?}: {c} vs OCCT {common}"
        );
    }
    let (u9, ..) = volumes(U9, Tool::Far);
    assert!((u9 - (6.0 + 0.75 * PI)).abs() < 1e-12);
    let (u10, ..) = volumes(U10, Tool::Far);
    assert!((u10 - (6.0 + 0.25 * PI)).abs() < 1e-12);
}
