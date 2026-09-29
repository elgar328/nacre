//! **One solid's edge on a ruling of the other's cylinder.**
//!
//! Two planes along a cylinder's axis meet in a line along it, so where a prism's corner lies on
//! the cylinder, the prism's edge there lies on a ruling — the whole line, not a point. Two
//! families put it there:
//!
//! - **two secant walls**: both of `B`'s walls at the corner cut `A`'s lateral, one entering `A`;
//! - **a tangent wall and a secant one — the keyhole**: a box as wide as the cylinder's diameter,
//!   its side walls tangent, its corners on the circle. The commonest shape of the two.
//!
//! Plus the shapes that must stay refused by name: both walls entering `A` (an inward wedge), and
//! the corner touching `A` from outside (a line contact).
//!
//! Every family runs over both frames (the world and [`pythagorean_frame`]), the polygon as
//! written and rotated by one vertex (which reorders the plane classes), and four placements
//! of `B` along `A`'s axis ([`HEIGHTS`]).

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, Profile2d, RejectClass, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// Where `B` stands along `A`'s axis, `(lift, height)` — `A` spans `0…2` on the frame. On `A`'s
/// base cap; inside `A`'s span; past `A`'s top; past both caps.
const HEIGHTS: [(f64, f64); 4] = [(0.0, 1.0), (0.5, 1.0), (0.5, 2.0), (-0.5, 3.0)];

/// Every boolean of `A` (extruded over `a`) against the prism over each of `quads`, over both
/// frames, both vertex orders and every [`HEIGHTS`] placement, handed to `check` with a label.
/// Returns how many booleans ran.
fn every_placement(
    a: fn() -> Profile2d,
    quads: &[(&str, &[[f64; 2]])],
    mut check: impl FnMut(&str, &mut Model, Result<Vec<Handle<Solid>>, BoolError>),
) -> usize {
    let mut ran = 0;
    for tilted in [false, true] {
        for &(what, quad) in quads {
            for rotated in [false, true] {
                let mut quad = quad.to_vec();
                if rotated {
                    quad.rotate_left(1);
                }
                for (lift, height) in HEIGHTS {
                    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                        for swapped in [false, true] {
                            let (mut m, x, y) = a_solid_and_a_placed_prism(
                                |m| {
                                    if tilted {
                                        pythagorean_frame(m, Point3::from_array([0.0; 3]))
                                    } else {
                                        nacre_ops::SketchFrame::world(m, nacre_exact::Axis::Z)
                                    }
                                },
                                a(),
                                &quad,
                                lift,
                                height,
                            );
                            let (x, y) = if swapped { (y, x) } else { (x, y) };
                            let at = format!(
                                "{what}, tilted = {tilted}, rotated = {rotated}, \
                                 lift = {lift}, height = {height}, {kind:?} swapped = {swapped}"
                            );
                            let got = boolean(&mut m, kind, x, y);
                            check(&at, &mut m, got);
                            ran += 1;
                        }
                    }
                }
            }
        }
    }
    ran
}

/// Today's answer: refused, as a coverage limit. The class is asserted, not the variant —
/// reason names are engine vocabulary; a refusal that turns `Impossible` or `SuspectedDefect`, or
/// a placement that starts to build, is news here and in the todo item this pins.
fn refused_today(at: &str, _: &mut Model, got: Result<Vec<Handle<Solid>>, BoolError>) {
    match got {
        Err(BoolError::Rejected { reason, .. }) => assert_eq!(
            reason.class(),
            RejectClass::NotSupported,
            "{at}: refused as {reason:?}"
        ),
        other => panic!("{at}: the wall moved — {other:?}"),
    }
}

/// ★ **Two secant walls, one entering `A` — the corner on the rim's seam (`(1, 0)`, `θ = 0` on
/// both frames) and off it (`(0.8, −0.6)`).** Every result would be a manifold solid; every
/// boolean is refused today. A corner `0.001` off the lateral builds.
#[test]
fn a_corner_entering_a_lateral_is_refused_today() {
    let on_the_seam: &[[f64; 2]] = &[[-1.0, 1.6], [1.0, 0.0], [3.0, 0.0], [3.0, 3.6]];
    let off_the_seam: &[[f64; 2]] = &[[-1.0, 1.6], [0.8, -0.6], [3.0, -0.6], [3.0, 3.6]];
    let ran = every_placement(
        unit_disk,
        &[("on the seam", on_the_seam), ("off the seam", off_the_seam)],
        refused_today,
    );
    assert_eq!(ran, 2 * 2 * 2 * 4 * 6, "the family");
    let off: &[[f64; 2]] = &[[-1.0, 1.6], [1.0, 0.001], [3.0, 0.001], [3.0, 3.6]];
    every_placement(unit_disk, &[("0.001 off", off)], |at, _, got| {
        got.unwrap_or_else(|e| panic!("{at}: {e:?}"));
    });
}

/// ★ **The keyhole**: the box `[−1, 1] × [0, 3]` against the unit cylinder — its walls `x = ±1`
/// tangent, its wall `y = 0` through the axis, its corners `(±1, 0)` on the circle. Every result
/// would be a manifold solid; every boolean is refused today.
#[test]
fn a_keyhole_is_refused_today() {
    let keyhole: &[[f64; 2]] = &[[-1.0, 0.0], [1.0, 0.0], [1.0, 3.0], [-1.0, 3.0]];
    let ran = every_placement(unit_disk, &[("the keyhole", keyhole)], refused_today);
    assert_eq!(ran, 2 * 2 * 4 * 6, "the family");
}

/// ★ **A tangent wall on the edge of a flat through `A`'s own axis.** `A` is the half cylinder
/// over [`upper_half_disk`], so its flat `y = 0` meets its lateral in the ruling `(1, 0)` — `A`'s
/// own edge. `B`'s wall `x = 1` is tangent to `A`'s cylinder along that edge, either running
/// across it or ending on it. Refused today — and where `B` reaches past both of `A`'s caps, the
/// corners on the line are `A`'s own (no point of `B` lies on the lateral) and the refusal comes
/// later, from the arrangement's backstops: `RingOrientation` (the wall across the edge, each
/// operand order with `B` first) and `LabelConflict` (the wall ending on it, rotated, `B` first)
/// — `SuspectedDefect`, 18 booleans, pinned here by count.
#[test]
fn a_tangent_wall_on_a_flats_edge_is_refused_today() {
    let across: &[[f64; 2]] = &[[-3.0, -2.0], [1.0, -2.0], [1.0, 0.5], [-3.0, 0.5]];
    let ending: &[[f64; 2]] = &[[-3.0, 0.0], [1.0, 0.0], [1.0, 0.5], [-3.0, 0.5]];
    let mut defects = 0;
    let ran = every_placement(
        upper_half_disk,
        &[("across the edge", across), ("ending on the edge", ending)],
        |at, m, got| match &got {
            Err(BoolError::Rejected { reason, .. })
                if reason.class() == RejectClass::SuspectedDefect =>
            {
                assert!(at.contains("lift = -0.5"), "{at}: {reason:?}");
                defects += 1;
            }
            _ => refused_today(at, m, got),
        },
    );
    assert_eq!(ran, 2 * 2 * 2 * 4 * 6, "the family");
    assert_eq!(defects, 18, "the backstop refusals past both caps");
}

/// ★ **Both walls entering `A` — an inward wedge**, its apex on the lateral at the seam and off
/// it. Refused today, and whatever this becomes it must not be a solid `validate` rejects.
#[test]
fn an_inward_wedge_on_a_lateral_is_refused_today() {
    let on_the_seam: &[[f64; 2]] = &[[1.0, 0.0], [-3.0, 3.0], [-3.0, -3.0]];
    let off_the_seam: &[[f64; 2]] = &[[0.8, -0.6], [-3.0, 3.0], [-3.0, -3.0]];
    let ran = every_placement(
        unit_disk,
        &[("on the seam", on_the_seam), ("off the seam", off_the_seam)],
        refused_today,
    );
    assert_eq!(ran, 2 * 2 * 2 * 4 * 6, "the family");
}

/// ★ **The corner touching `A` from outside — a line contact.** Common would be empty, Cut the
/// operand, Fuse two bodies touching along the line. Refused today.
#[test]
fn a_line_contact_on_a_lateral_is_refused_today() {
    let touching: &[[f64; 2]] = &[[1.0, 0.0], [3.0, -2.0], [4.0, 0.0], [3.0, 2.0]];
    let ran = every_placement(unit_disk, &[("a line contact", touching)], refused_today);
    assert_eq!(ran, 2 * 2 * 4 * 6, "the family");
}
