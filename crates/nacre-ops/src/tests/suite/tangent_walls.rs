//! **Walls tangent to a rod, where the material falls apart along the line of contact** — the unit
//! rod `z ∈ [0, 5]`, and the bore it leaves in the block `[−2, 2]² × [0, 5]`, against boxes over
//! `z ∈ [1/2, 7/2]` whose walls touch the lateral along a ruling (`x = ±1` or `y = ±1`), each picture
//! also moved whole by the plane similarities [`TURNED`] and [`ROOT_2`] — the rod's circle and seam,
//! its notch, the box, its window and the block together — so the walls stand off the axes and the
//! radius is stated by its square; each case is a test in each frame (`in_each_frame!`). Where a
//! boolean keeps only the two wedges between the lateral and a wall, the result is two bodies meeting
//! along that line, and the wall's class draws the line as a knife (`WorkingCyl::parting`) so its face
//! parts there (`per_class::knives_as_edges`):
//!
//! - **a box inside the bore's height** — the class of the box's wall holds two rings, the bore's
//!   section and the box's face inside it, and the knife runs from one onto the other and back; the
//!   stretch along the face parts it;
//! - **a windowed box** — the box with a window `[1/2, 3/2] × [−w, w] × [3/2, 5/2]` (`w` two fifths of
//!   the box's half width) through its wall `x = 1` across the line: each of the knife's two stretches
//!   joins the face's outer ring to the window's, and only the two together part the face;
//! - **a half notch** — the rod with `y < 0` taken over `z ∈ [1, 2]`: there the lateral ends on the
//!   line, its edge held by the notch's wall `y = 0`, and the knife over that stretch is the plane
//!   segment of the lateral's half and the wall's (`crossings::fold_onto_rulings`), which joins the
//!   ruling knife's stretches below and above it.
//!
//! **Derived, not read back.** Volumes are closed form ([`hole_runs`](super::hole_runs)' disk–rectangle
//! areas; a window's, the difference of two boxes'). A similarity keeps a picture's bodies and scales
//! its volumes by [`Similar::scale2`], so every frame is held to the stated frame's answers. In the
//! stated frame the body counts come from a voxel model of the same solids outside the kernel
//! (components of the kept set on a grid `0.01` across and `0.05` up), which agrees with every
//! multi-body result the kernel builds in the sweep these cells were taken from (308 of them), and
//! OCCT agrees with the volume and area of every cell, as built, to the six digits it prints.

use super::hole_runs::{BELOW, Body, Rect, SAME, Similar, WHOLE, area};
use super::one_sided_notch::{four, prism, scored};
use super::scoreboard::{OPS, run};
use super::*;

/// The rod as it is.
const PLAIN: Body = Body {
    cuts: &[],
    kept: [&[WHOLE], &[WHOLE], &[WHOLE], &[WHOLE]],
};

/// The height the notch takes the half `y < 0` over — the second of the slabs a [`Body`] keeps.
const NOTCH_Z: [f64; 2] = [1.0, 2.0];

/// The rod with the half `y < 0` taken over [`NOTCH_Z`].
const HALF: Body = Body {
    cuts: &[(BELOW, NOTCH_Z[0], NOTCH_Z[1] - NOTCH_Z[0])],
    kept: [&[WHOLE], &[[-2.0, 0.0, 2.0, 2.0]], &[WHOLE], &[WHOLE]],
};

/// The box's span of height, inside the rod's.
const SPAN: [f64; 2] = [0.5, 3.5];

/// Walls `x = ±1` with the line inside their faces, and walls `y = ±1` likewise.
const ACROSS_X: Rect = [-1.0, -0.25, 1.0, 0.25];
const ACROSS_Y: Rect = [-0.5, -1.0, 0.5, 1.0];
/// Walls `x = ±1` the whole block wide.
const WIDE_X: Rect = [-1.0, -2.0, 1.0, 2.0];
/// All four walls touching the rod.
const SNUG: Rect = [-1.0, -1.0, 1.0, 1.0];

/// The window through the wall `x = 1` of a box over `rect`, and the part of it inside the box.
fn window(rect: Rect) -> (Rect, Rect) {
    let w = 0.4 * rect[1].abs().min(rect[3].abs());
    ([0.5, -w, 1.5, w], [0.5, -w, rect[2], w])
}
const WINDOW_Z: [f64; 2] = [1.5, 2.5];

/// The picture turned by the 3-4-5 rotation: the tangent rulings and the seam `(3/5, 4/5)` stand off
/// the axes.
const TURNED: Similar = Similar { a: 3, b: 4, k: 5 };

/// The picture turned by 45° and scaled by `√2`: the rod's `r² = 2` is stated by its square, and each
/// tangent wall (`x + y = 2`, say) has a rational equation but an irrational unit normal.
const ROOT_2: Similar = Similar { a: 1, b: 1, k: 1 };

/// A `#[test]` of each case in each frame — `case::stated`, `case::turned_3_4_5` ([`TURNED`]) and
/// `case::scaled_by_root_2` ([`ROOT_2`]) — so a frame that fails names itself and the others still
/// run.
macro_rules! in_each_frame {
    ($($case:ident),* $(,)?) => {$(
        mod $case {
            #[test]
            fn stated() {
                super::$case(super::SAME)
            }
            #[test]
            fn turned_3_4_5() {
                super::$case(super::TURNED)
            }
            #[test]
            fn scaled_by_root_2() {
                super::$case(super::ROOT_2)
            }
        }
    )*};
}

/// **Every boolean of `body` (or its bore) with the box over `rect × SPAN`, windowed or not, the
/// picture moved by `frame`** — each at its closed-form volume, `validate` and the mesh (`scored`).
/// `bodies` is the model's count per boolean, in the order of [`OPS`]; `x` is a body touching itself
/// along the frame's image of `(0, −1)` above and below the notch, which the boolean refuses
/// (`SelfTouchingResult`) with a witness on that line there — over the notch the slab under it fills
/// the line.
fn against(frame: Similar, body: &Body, bore: bool, rect: Rect, windowed: bool, bodies: &str) {
    let dz = |z: [f64; 2]| z[1] - z[0];
    let va = body.within(bore, [-9.0, -9.0, 9.0, 9.0], [-9.0, 9.0]);
    let (hole, inside) = window(rect);
    let (mut vb, mut common) = (area(rect) * dz(SPAN), body.within(bore, rect, SPAN));
    if windowed {
        vb -= area(inside) * dz(WINDOW_Z);
        common -= body.within(bore, inside, WINDOW_Z);
    }
    let s = frame.scale2();
    let mut got = String::new();
    for (op, want) in OPS.into_iter().zip(four(va * s, vb * s, common * s)) {
        let mut m = Model::new();
        let a = body.image(&mut m, bore, frame);
        let mut b = prism(&mut m, &frame.corners(rect), SPAN[0], dz(SPAN));
        if windowed {
            let w = prism(&mut m, &frame.corners(hole), WINDOW_Z[0], dz(WINDOW_Z));
            let out = crate::boolean(&mut m, BoolKind::Cut, b, w).expect("the window");
            assert_eq!(out.len(), 1, "the box stays one body");
            b = out[0];
        }
        let label = format!("{frame:?}: bore = {bore}, {rect:?}, windowed = {windowed}, {op:?}");
        if bodies.as_bytes().get(got.len()) == Some(&b'x') {
            let Err(BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                at: Some(crate::RejectWhere::Point(p)),
            }) = run(&mut m, op, a, b)
            else {
                panic!("{label}: not refused as touching itself");
            };
            let [x, y, z] = p.as_array();
            let [wx, wy] = frame.point([0.0, -1.0]).map(nacre_exact::Rat::to_f64);
            assert!(
                (x - wx).abs() < 1e-9 && (y - wy).abs() < 1e-9,
                "{label}: the witness {p:?}"
            );
            let over_notch = NOTCH_Z[0] < z && z < NOTCH_Z[1];
            assert!(
                (SPAN[0]..=SPAN[1]).contains(&z) && !over_notch,
                "{label}: the witness {p:?}"
            );
            got.push('x');
            continue;
        }
        let n =
            scored(&label, &mut m, op, a, b, want).unwrap_or_else(|| panic!("{label}: refused"));
        got.push_str(&n.to_string());
    }
    assert_eq!(
        got, bodies,
        "{frame:?}: bore = {bore}, {rect:?}, windowed = {windowed}"
    );
}

/// A box whose walls touch the rod, less the rod, is four slivers; the bore and the box share the
/// same four.
fn the_wedges_beside_a_tangent_wall_part_its_face(frame: Similar) {
    against(frame, &PLAIN, false, ACROSS_X, false, "1141");
    against(frame, &PLAIN, false, ACROSS_Y, false, "1141");
    against(frame, &PLAIN, true, ACROSS_X, false, "1114");
    against(frame, &PLAIN, true, ACROSS_Y, false, "1114");
}

/// The window's two sides of the line are two bodies: the knife's stretches below and above the
/// window part the wall's face only together.
fn a_window_across_the_line_parts_the_face_on_both_sides_of_it(frame: Similar) {
    against(frame, &PLAIN, false, WIDE_X, true, "1321");
    against(frame, &PLAIN, false, ACROSS_X, true, "1241");
    against(frame, &PLAIN, true, WIDE_X, true, "1112");
    against(frame, &PLAIN, true, ACROSS_X, true, "1114");
}

/// The notch's wall ends the lateral on the line over `z ∈ [1, 2]`, and the wall `x = 1` parts there
/// as well as above and below: two bodies, or three where the walls `y = ±1/4` cut the slivers off.
fn a_notch_ending_the_lateral_on_the_line_parts_the_face_there_too(frame: Similar) {
    against(frame, &HALF, false, WIDE_X, false, "1221");
    against(frame, &HALF, false, ACROSS_X, false, "1131");
    against(frame, &HALF, true, WIDE_X, false, "1112");
    against(frame, &HALF, true, ACROSS_X, false, "1113");
}

/// With the walls `y = ±1` touching too, the slab under the notch joins the two slivers by
/// `y = −1`, which then touch each other along `(0, −1)` above and below it: one body touching
/// itself, refused there; the other booleans build.
fn a_notch_that_joins_the_slivers_by_a_tangent_wall_is_refused_on_that_line(frame: Similar) {
    against(frame, &HALF, false, SNUG, false, "12x1");
    against(frame, &HALF, true, SNUG, false, "111x");
}

in_each_frame!(
    the_wedges_beside_a_tangent_wall_part_its_face,
    a_window_across_the_line_parts_the_face_on_both_sides_of_it,
    a_notch_ending_the_lateral_on_the_line_parts_the_face_there_too,
    a_notch_that_joins_the_slivers_by_a_tangent_wall_is_refused_on_that_line,
);
