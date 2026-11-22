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
//! Each against eleven tools, as Fuse, both Cuts and Common: a box far off; the box
//! `[1/2, 3/2] × [−1/2, 1/2] × [0, 2]` across the line `(1, 0)`; the slot
//! `[−2, 2] × [−1/4, 1/4] × [1/5, 9/10]` across both lines; and the rod's own cylinder over eight
//! spans of `z` (`COAXIAL`) — inside the box's height or past it, ending on `U`'s own heights
//! or between them. The coaxial rods are also run as the first operand of the union and the
//! common (`a_coaxial_rod_answers_as_the_first_operand_too`). The volumes are closed forms
//! (`volumes`, `rod_bodies`), and OCCT 7.9 builds every boolean of the unmirrored box and slot
//! with them, and the common of every unmirrored coaxial row
//! (`the_junctions_closed_forms_are_occts`). Its face counts for the box and the slot are these but
//! in two places: each union with the box across has two more, as the unstaggered slot's union
//! does (which builds here too); and `U9`'s box across less `U`, and their common part, have one
//! more. In those two the box's top keeps a sliver between the rod's arc and the wall tangent to
//! it, whose mesh the plane faces' chord budget holds (`mesh_covers_faces`).
//!
//! **A box rounded on one corner against its fillet's cylinder** (`FILLET`) is the same line
//! held by the flats alone: the square `[0, 4]²` two tall with its `(4, 0)` corner rounded by the
//! unit circle about `(3, 1)`, against that circle's rod over six spans — no plane runs through
//! the axis there.
//!
//! **The rod's cylinder notched** (`NOTCHED`) is the same cylinder less its side `x < −1/2`, over
//! the coaxial spans, as the four booleans (`a_notched_rod_on_the_junctions_cylinder`) and as the
//! first operand. Each of its ends is an arc, not a whole rim, and holds the line `(1, 0)` inside
//! it; the notch takes the line `(−1, 0)` away. OCCT 7.9 agrees with these volumes on every
//! unmirrored notched boolean (where they are empty it builds no body) but one: `U9` less the rod
//! over `[1/2, 1]`, where the rod's top above the box and the box beside the void meet along the
//! arc `z = 1` on the box's side and the material runs round through the notch's flat. OCCT builds
//! that as one solid whose arc four faces use; here it is `NonManifoldResultEdge`.
//!
//! OCCT's face counts for the coaxial rows are these in the cuts but four (one more each); in the
//! unions they are these in one row and up to four more in the rest, in the commons these in two
//! rows and up to two more in the rest.
//!
//! A rim of the coaxial rod, or an end arc of the notched one, meets a tangent line inside a piece
//! the tangent wall's class states — the wall's own edge, or the plane through the axis's segment
//! along it — and the cap that bears it touches that class's plane at the point alone. The class
//! cuts the piece there all the same: a lateral face's boundary meeting the line is a vertex of it
//! (`Trace::tangent_points`).

use super::scoreboard::{OPS, Op, Want, check, counted, refused, run};
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
    /// The rod's own cylinder over `z ∈ [z0, z1]`, in tenths.
    Coaxial(i8, i8),
    /// The same, less its side `x < −1/2`.
    Notched(i8, i8),
}

const TOOLS: [Tool; 3] = [Tool::Far, Tool::Across, Tool::Slot];

/// A height in tenths, as a coaxial span states it.
fn tenths(t: i8) -> f64 {
    f64::from(t) / 10.0
}

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
        Tool::Coaxial(z0, z1) => rod(m, tenths(z0), tenths(z1)),
        Tool::Notched(z0, z1) => {
            let r = rod(m, tenths(z0), tenths(z1));
            let side = crate::fixtures::cuboid(m, [-2.0, -2.0, -1.0], [-0.5, 2.0, 3.0]);
            let n = crate::boolean(m, BoolKind::Cut, r, side).expect("the notched rod builds");
            assert_eq!(n.len(), 1, "the notched rod is one body");
            n[0]
        }
    }
}

/// `∫ₐᵇ √(1 − y²) dy`.
fn chord_area(a: f64, b: f64) -> f64 {
    let f = |y: f64| 0.5 * (y * (1.0 - y * y).sqrt() + y.asin());
    f(b) - f(a)
}

/// The area of a rod's section on `U`'s cylinder: the unit disk, or the disk less its segment
/// `x < −1/2`. Both are symmetric about `y = 0`.
fn section(t: Tool) -> f64 {
    match t {
        Tool::Coaxial(..) => PI,
        Tool::Notched(..) => PI - 2.0 * chord_area(-1.0, -0.5),
        _ => unreachable!("{t:?} is no rod on U's cylinder"),
    }
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
        (Tool::Coaxial(..) | Tool::Notched(..), _, true) => section(t),
        // The box's half of the section.
        (Tool::Coaxial(..) | Tool::Notched(..), true, false) => section(t) / 2.0,
    }
}

/// The tool's own height span.
fn tool_span(t: Tool) -> [f64; 2] {
    match t {
        Tool::Far => [0.0, 1.0],
        Tool::Across => [0.0, 2.0],
        Tool::Slot => [0.2, 0.9],
        Tool::Coaxial(z0, z1) | Tool::Notched(z0, z1) => [tenths(z0), tenths(z1)],
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
        Tool::Coaxial(z0, z1) | Tool::Notched(z0, z1) => section(t) * (tenths(z1) - tenths(z0)),
    };
    (u, tv, c)
}

fn wants(cz: [f64; 2], t: Tool) -> [Want; 4] {
    let (u, tv, c) = volumes(cz, t);
    // Faces per body: U9's, U10's.
    let f = |u9: usize, u10: usize| if cz == U9 { u9 } else { u10 };
    match t {
        Tool::Coaxial(..) | Tool::Notched(..) => {
            let row = rod_row(cz, t);
            core::array::from_fn(|i| answer(rod_bodies(cz, t, OPS[i]), row[i]))
        }
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

/// What one coaxial boolean answers: each body's face count, the bodies in the order of the
/// closed form's volumes ([`rod_bodies`], [`fillet_bodies`]), or the refusal it names.
#[derive(Clone, Copy, Debug)]
enum Answer {
    Faces(&'static [usize]),
    Refused(RejectReason),
}

const NON_MANIFOLD: Answer = Answer::Refused(RejectReason::NonManifoldResultEdge);

const fn f(faces: &'static [usize]) -> Answer {
    Answer::Faces(faces)
}

fn answer(volumes: Vec<f64>, a: Answer) -> Want {
    match a {
        Answer::Faces(faces) => counted(&volumes, faces),
        Answer::Refused(reason) => refused(reason),
    }
}

/// A rod's span `[z0, z1]`, in tenths.
type Span = (i8, i8);

/// A coaxial row's six answers ([`COAXIAL`], [`NOTCHED`]).
type CoaxialRow = [Answer; 6];

/// The coaxial rows: per span (tenths), for `U9` then `U10`, the answers as Fuse, `U` less the
/// rod, the rod less `U`, Common, and then the union and the common with the rod first. The
/// mirrored junction answers as the unmirrored one.
const COAXIAL: [(Span, [CoaxialRow; 2]); 8] = [
    (
        (6, 12),
        [
            [f(&[9]), f(&[3, 10]), f(&[]), f(&[3]), f(&[9]), f(&[3])],
            [f(&[9]), f(&[10]), f(&[5]), f(&[5]), f(&[9]), f(&[5])],
        ],
    ),
    (
        (4, 12),
        [
            [f(&[9]), f(&[3, 8]), f(&[4]), f(&[5]), f(&[9]), f(&[5])],
            [f(&[9]), f(&[10]), f(&[5]), f(&[5]), f(&[9]), f(&[5])],
        ],
    ),
    (
        (5, 10),
        [
            [f(&[9]), f(&[3, 8]), f(&[]), f(&[3]), f(&[9]), f(&[3])],
            [f(&[8]), f(&[10]), f(&[4]), f(&[5]), f(&[8]), f(&[5])],
        ],
    ),
    (
        (2, 8),
        [
            [f(&[9]), f(&[11]), f(&[4]), f(&[5]), f(&[9]), f(&[5])],
            [f(&[10]), f(&[10]), f(&[4, 4]), f(&[7]), f(&[10]), f(&[7])],
        ],
    ),
    (
        (6, 9),
        [
            [f(&[9]), f(&[13]), f(&[]), f(&[3]), f(&[9]), f(&[3])],
            [f(&[10]), f(&[12]), f(&[4]), f(&[5]), f(&[10]), f(&[5])],
        ],
    ),
    (
        (0, 20),
        [
            [f(&[7]), f(&[6]), f(&[4, 3]), f(&[5]), f(&[7]), f(&[5])],
            [f(&[7]), f(&[6]), f(&[4, 5]), f(&[7]), f(&[7]), f(&[7])],
        ],
    ),
    (
        (3, 5),
        [
            [f(&[9]), f(&[11]), f(&[4]), f(&[4]), f(&[9]), f(&[4])],
            [f(&[10]), f(&[14]), f(&[]), f(&[3]), f(&[10]), f(&[3])],
        ],
    ),
    (
        (10, 14),
        [
            [f(&[9]), f(&[3, 8]), f(&[]), f(&[3]), f(&[9]), f(&[3])],
            [f(&[13]), f(&[10]), f(&[3]), f(&[]), f(&[13]), f(&[])],
        ],
    ),
];

/// The notched rows: as [`COAXIAL`], with the rod's cylinder less its side `x < −1/2`.
const NOTCHED: [(Span, [CoaxialRow; 2]); 8] = [
    (
        (6, 12),
        [
            [f(&[9]), f(&[13]), f(&[]), f(&[4]), f(&[9]), f(&[4])],
            [f(&[12]), f(&[13]), f(&[6]), f(&[6]), f(&[12]), f(&[6])],
        ],
    ),
    (
        (4, 12),
        [
            [f(&[11]), f(&[13]), f(&[5]), f(&[6]), f(&[11]), f(&[6])],
            [f(&[12]), f(&[13]), f(&[6]), f(&[6]), f(&[12]), f(&[6])],
        ],
    ),
    (
        (5, 10),
        [
            [f(&[9]), NON_MANIFOLD, f(&[]), f(&[4]), f(&[9]), f(&[4])],
            [f(&[11]), f(&[13]), f(&[5]), f(&[6]), f(&[11]), f(&[6])],
        ],
    ),
    (
        (2, 8),
        [
            [f(&[11]), f(&[13]), f(&[5]), f(&[6]), f(&[11]), f(&[6])],
            [f(&[14]), f(&[14]), f(&[5, 5]), f(&[8]), f(&[14]), f(&[8])],
        ],
    ),
    (
        (6, 9),
        [
            [f(&[9]), f(&[13]), f(&[]), f(&[4]), f(&[9]), f(&[4])],
            [f(&[12]), f(&[14]), f(&[5]), f(&[6]), f(&[12]), f(&[6])],
        ],
    ),
    (
        (0, 20),
        [
            [f(&[12]), f(&[11]), f(&[5, 4]), f(&[6]), f(&[12]), f(&[6])],
            [f(&[13]), f(&[12]), f(&[5, 6]), f(&[8]), f(&[13]), f(&[8])],
        ],
    ),
    (
        (3, 5),
        [
            [f(&[11]), f(&[12]), f(&[5]), f(&[5]), f(&[11]), f(&[5])],
            [f(&[10]), f(&[14]), f(&[]), f(&[4]), f(&[10]), f(&[4])],
        ],
    ),
    (
        (10, 14),
        [
            [f(&[9]), f(&[11]), f(&[]), f(&[4]), f(&[9]), f(&[4])],
            [f(&[14]), f(&[10]), f(&[4]), f(&[]), f(&[14]), f(&[])],
        ],
    ),
];

/// The row of rod `t` on `U`'s cylinder, for `U9` or `U10`.
fn rod_row(cz: [f64; 2], t: Tool) -> CoaxialRow {
    let (table, span) = match t {
        Tool::Coaxial(z0, z1) => (&COAXIAL, (z0, z1)),
        Tool::Notched(z0, z1) => (&NOTCHED, (z0, z1)),
        _ => unreachable!("{t:?} is no rod on U's cylinder"),
    };
    let (_, rows) = table
        .iter()
        .find(|(s, _)| *s == span)
        .expect("a coaxial span");
    rows[usize::from(cz == U10)]
}

/// The pieces over `[lo, hi]` of a rod of section `area` on `U`'s cylinder that `U` leaves whole:
/// half the section over the box's height, all of it above.
fn beside_u(lo: f64, hi: f64, area: f64) -> f64 {
    if hi <= lo {
        return 0.0;
    }
    let half = (hi.min(1.0) - lo.max(0.0)).max(0.0);
    area * (hi - lo) - area / 2.0 * half
}

/// Each body's volume, ascending, of the boolean `op` of `U` and the rod `t` on its cylinder — the
/// closed forms of [`volumes`], parted where the result is. `U` less the rod parts the rod's top
/// beyond the cut from the rest only where the cut takes `U`'s whole disk: the whole rod, when the
/// cut reaches the box's top from below or stands above it (`1 ≤ z1 < cz₁`); the notched rod leaves
/// `U`'s side `x < −1/2`, which holds the top on. The rod less `U` is the stretch below `U`'s whole
/// disk and the stretch above it.
fn rod_bodies(cz: [f64; 2], t: Tool, op: Op) -> Vec<f64> {
    let (u, tv, c) = volumes(cz, t);
    let [z0, z1] = tool_span(t);
    let area = section(t);
    let takes_the_disk = matches!(t, Tool::Coaxial(..));
    let mut out = match op {
        Op::Fuse => vec![u + tv - c],
        Op::Common => vec![c],
        Op::CutAb if takes_the_disk && cz[1] > 1.0 && (1.0..cz[1]).contains(&z1) => {
            let top = PI * (cz[1] - z1);
            vec![top, u - c - top]
        }
        Op::CutAb => vec![u - c],
        Op::CutBa => vec![
            beside_u(z0, z1.min(cz[0]), area),
            beside_u(z0.max(cz[1]), z1, area),
        ],
    };
    out.retain(|&v| v > 1e-12);
    out.sort_by(f64::total_cmp);
    out
}

#[test]
fn a_box_joined_to_a_cylinder_along_its_tangent_line() {
    for cz in [U9, U10] {
        for mirrored in [false, true] {
            let coaxial = COAXIAL.iter().map(|&((z0, z1), _)| Tool::Coaxial(z0, z1));
            for t in TOOLS.into_iter().chain(coaxial) {
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

/// **The rod's cylinder notched** — its side `x < −1/2` cut away, over the coaxial spans — against
/// both junctions as the four booleans: the rows state the answers ([`NOTCHED`]).
#[test]
fn a_notched_rod_on_the_junctions_cylinder() {
    for cz in [U9, U10] {
        for mirrored in [false, true] {
            for &((z0, z1), _) in &NOTCHED {
                let t = Tool::Notched(z0, z1);
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

/// **A coaxial rod as the first operand**, whole or notched — the union and the common with the
/// rod first come out as the rows state ([`COAXIAL`], [`NOTCHED`]), with the same closed forms.
#[test]
fn a_coaxial_rod_answers_as_the_first_operand_too() {
    let whole = COAXIAL.iter().map(|&((z0, z1), _)| Tool::Coaxial(z0, z1));
    let notched = NOTCHED.iter().map(|&((z0, z1), _)| Tool::Notched(z0, z1));
    let rods: Vec<Tool> = whole.chain(notched).collect();
    for cz in [U9, U10] {
        for mirrored in [false, true] {
            for &t in &rods {
                let row = rod_row(cz, t);
                for (kind, op, a) in [
                    (BoolKind::Fuse, Op::Fuse, row[4]),
                    (BoolKind::Common, Op::Common, row[5]),
                ] {
                    let mut m = Model::new();
                    let u = junction(&mut m, cz, mirrored);
                    let b = tool(&mut m, t);
                    let got = crate::boolean(&mut m, kind, b, u);
                    let label = format!("rod {cz:?}, mirrored = {mirrored}, {t:?} first, {kind:?}");
                    check(&label, &m, got, &answer(rod_bodies(cz, t, op), a));
                }
            }
        }
    }
}

/// The fillet rows: per span (tenths) of the fillet's rod, the answers as Fuse, the rounded box
/// less the rod, the rod less the box, Common.
const FILLET: [(Span, [Answer; 4]); 6] = [
    ((6, 14), [f(&[7]), f(&[11]), f(&[]), f(&[3])]),
    ((-4, 12), [f(&[8]), f(&[9]), f(&[3]), f(&[3])]),
    ((12, 28), [f(&[8]), f(&[9]), f(&[3]), f(&[3])]),
    ((0, 20), [f(&[7]), f(&[7]), f(&[]), f(&[3])]),
    ((-4, 28), [f(&[9]), f(&[7]), f(&[3, 3]), f(&[3])]),
    ((0, 12), [f(&[7]), f(&[9]), f(&[]), f(&[3])]),
];

/// The square `[0, 4]²`, two tall, its `(4, 0)` corner rounded by the unit circle about `(3, 1)`.
fn rounded_box(m: &mut Model) -> Handle<Solid> {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let profile = stated(vec![
        line(p(0.0, 0.0), p(3.0, 0.0)),
        arc_turns(p(3.0, 1.0), p(3.0, 0.0), 1),
        line(p(4.0, 1.0), p(4.0, 4.0)),
        line(p(4.0, 4.0), p(0.0, 4.0)),
        line(p(0.0, 4.0), p(0.0, 0.0)),
    ])
    .expect("the rounded square")
    .remove(0);
    let frame = SketchFrame::world(m, nacre_exact::Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    )
    .expect("the rounded box extrudes") else {
        unreachable!()
    };
    solid
}

/// The fillet's rod over `z ∈ [z0, z1]`.
fn fillet_rod(m: &mut Model, z0: f64, z1: f64) -> Handle<Solid> {
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
            profile: circle_profile([3.0, 1.0], 1.0),
            dist: z1 - z0,
        },
    )
    .expect("the fillet's rod extrudes") else {
        unreachable!()
    };
    solid
}

/// Each body's volume, ascending, of the fillet boolean `op`: the rounded box is `2·(15 + π/4)`, its
/// section holds the whole disk, and the rod less the box is the stretch below it and the stretch
/// above it.
fn fillet_bodies(span: Span, op: Op) -> Vec<f64> {
    let [z0, z1] = [tenths(span.0), tenths(span.1)];
    let plate = 2.0 * (15.0 + PI / 4.0);
    let rod = PI * (z1 - z0);
    let c = PI * (z1.min(2.0) - z0.max(0.0)).max(0.0);
    let mut out = match op {
        Op::Fuse => vec![plate + rod - c],
        Op::CutAb => vec![plate - c],
        Op::CutBa => vec![
            PI * (z1.min(0.0) - z0).max(0.0),
            PI * (z1 - z0.max(2.0)).max(0.0),
        ],
        Op::Common => vec![c],
    };
    out.retain(|&v| v > 1e-12);
    out.sort_by(f64::total_cmp);
    out
}

#[test]
fn a_box_rounded_on_one_corner_against_its_fillets_cylinder() {
    for &(span, row) in &FILLET {
        for (op, a) in OPS.into_iter().zip(row) {
            let mut m = Model::new();
            let b0 = rounded_box(&mut m);
            let rod = fillet_rod(&mut m, tenths(span.0), tenths(span.1));
            let got = run(&mut m, op, b0, rod);
            check(
                &format!("fillet {span:?}, {op:?}"),
                &m,
                got,
                &answer(fillet_bodies(span, op), a),
            );
        }
    }
}

/// The closed forms against the identities they must satisfy and the values OCCT 7.9 reports for
/// the unmirrored `U ∩ tool`.
#[test]
#[expect(
    clippy::approx_constant,
    reason = "OCCT prints six places: `U9`'s common with the rod over [1/2, 1] is π/2 to them"
)]
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
    // The unmirrored commons of the rods on `U`'s cylinder, whole and notched: `U9`'s, `U10`'s
    // (OCCT builds no body for `U10`'s rod above the box).
    let rods = [
        (Tool::Coaxial(6, 12), 1.88496, 0.863938),
        (Tool::Coaxial(4, 12), 2.35619, 1.49226),
        (Tool::Coaxial(5, 10), 1.5708, 1.1781),
        (Tool::Coaxial(2, 8), 1.41372, 1.72788),
        (Tool::Coaxial(6, 9), 0.942478, 0.706858),
        (Tool::Coaxial(0, 20), 3.92699, 2.35619),
        (Tool::Coaxial(3, 5), 0.314159, 0.628319),
        (Tool::Coaxial(10, 14), 1.25664, 0.0),
        (Tool::Notched(6, 12), 1.51644, 0.695037),
        (Tool::Notched(4, 12), 1.89556, 1.20052),
        (Tool::Notched(5, 10), 1.2637, 0.947778),
        (Tool::Notched(2, 8), 1.13733, 1.39007),
        (Tool::Notched(6, 9), 0.758222, 0.568667),
        (Tool::Notched(0, 20), 3.15926, 1.89556),
        (Tool::Notched(3, 5), 0.252741, 0.505482),
        (Tool::Notched(10, 14), 1.01096, 0.0),
    ];
    for (t, c9, c10) in rods {
        for (cz, common) in [(U9, c9), (U10, c10)] {
            let (u, tv, c) = volumes(cz, t);
            assert!(
                (c - common).abs() < 1e-5,
                "{cz:?} {t:?}: {c} vs OCCT {common}"
            );
            // The bodies' parts make up the wholes.
            let sum = |op| rod_bodies(cz, t, op).iter().sum::<f64>();
            assert!((sum(Op::CutAb) - (u - c)).abs() < 1e-12, "{cz:?} {t:?}");
            assert!((sum(Op::CutBa) - (tv - c)).abs() < 1e-12, "{cz:?} {t:?}");
        }
    }
    let (u9, ..) = volumes(U9, Tool::Far);
    assert!((u9 - (6.0 + 0.75 * PI)).abs() < 1e-12);
    let (u10, ..) = volumes(U10, Tool::Far);
    assert!((u10 - (6.0 + 0.25 * PI)).abs() < 1e-12);
}
