//! **A plane that *carries* an arrangement line is a fourth plane through every point on it.**
//!
//! Four planes through one point have several valid names, and the engine folds them onto one —
//! otherwise two classes name the same point differently and the seam table, which welds by name,
//! sees two vertices at one coordinate and refuses (`SeamAlias`, a `SuspectedDefect`: the kernel
//! blaming itself).
//!
//! The fold happened for concurrencies at *input* vertices and for the ones where two planes
//! **cross** a line at a point. It did not happen for a plane that **contains** the line: such a
//! plane is in that line's own direction family, so it never becomes a handle and never orders
//! against anything — there was nothing to notice. `Aliases::wall` already recorded it (folding two
//! walls onto one line *is* the statement that a second plane carries it), and the report now reads
//! that back out.
//!
//! ★ The shapes here are the ones a user hit in the playground, reduced. What made them worth a
//! file is that the refused model is **valid** — the prism stops at `y = 0.9` and the wall it is
//! coplanar with lives at `y ∈ [1,2]`, so nothing touches anything.

use nacre_math::Point2;
use nacre_ops::{BoolKind, Operation, Profile2d, RejectReason, SketchFrame, apply, boolean};
use nacre_scalar::Axis;
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn prism(m: &mut Model, axis: Axis, pts: &[[f64; 2]], h: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, axis),
            profile,
            dist: h,
        },
    )
    .expect("extrude") else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// A staircase block whose vertical step stands at `x = step`, and a triangular prism on the ZX
/// plane whose **apex edge** runs at `x = apex`, extruded `dist` along +y.
///
/// The apex edge is the meet of the prism's two slanted faces, so when `apex == step` those two
/// planes and the staircase's step plane all carry **one line** — three planes sharing a line, the
/// configuration this file is about.
fn stair_and_prism(step: f64, apex: f64, dist: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = prism(
        &mut m,
        Axis::Z,
        &[
            [0.0, 0.0],
            [0.0, 2.0],
            [step, 2.0],
            [step, 1.0],
            [0.5, 1.0],
            [0.5, 0.5],
            [2.0, 0.5],
            [2.0, 0.0],
        ],
        1.0,
    );
    // ZX: sketch +u is +z, +v is +x.
    let b = prism(
        &mut m,
        Axis::Y,
        &[[0.5, apex], [0.0, 2.0], [1.0, 2.0]],
        dist,
    );
    (m, a, b)
}

/// ★ Every case in this file validates. It did not always: at `apex = 0.7` a vertex used to land
/// `5.55e-17` off one of its defining planes with a recorded tolerance of `0`, so this helper had a
/// second copy that skipped `validate`. That was a **separate** defect — the tolerance measured the
/// triple the arrangement computed the point with, while the definition named the triple the result
/// re-derived — and it is fixed, so the exception is gone with it. Leaving it would have made this
/// file quiet the next time the same thing broke.
fn fused_volume(step: f64, apex: f64, dist: f64) -> f64 {
    let (mut m, a, b) = stair_and_prism(step, apex, dist);
    let got = boolean(&mut m, BoolKind::Fuse, a, b).expect("the fuse builds");
    m.rebuild_adjacency();
    assert_eq!(got.len(), 1, "one body");
    assert!(nacre_validate::validate(&m).is_empty());
    nacre_props::mass_props(&m, got[0]).expect("props").volume
}

/// ★ **The case the report is for.** The prism's apex edge is coplanar with the staircase's step
/// wall, and the edge pierces the staircase's `y = 0.5` face — four planes at that point, three of
/// them sharing the apex line.
///
/// Before the report this was `SeamAlias`: one class named the point with the folded wall
/// (`[3,6,11]`), another with the two slanted planes (`[3,11,12]`), and the seam table saw two
/// vertices at one coordinate. **The body is valid** — the prism ends at `y = 0.9` and the wall it
/// shares a plane with is at `y ∈ [1,2]`, so they never meet.
#[test]
fn a_prism_coplanar_with_a_far_away_wall_still_builds() {
    let v = fused_volume(1.0, 1.0, 0.9);
    // The staircase is 1×(2·1 + 0.5·1 + 0.5·1) = 1.75 by area × height 1; the prism adds the part
    // of its triangle that is outside the staircase. Locked as a number so a silent change speaks.
    assert!(v > 1.75, "the fuse must add material: volume {v}");
    assert!(
        (v - fused_volume(1.05, 1.0, 0.9)).abs() > 1e-12,
        "the two steps are different shapes; if this is equal the fixture stopped varying"
    );
}

/// The negative control the fix is paired with: nudge the step off the apex's plane and nothing is
/// concurrent. It built before and must keep building.
#[test]
fn a_step_off_the_apex_plane_is_unchanged() {
    let v = fused_volume(1.05, 1.0, 0.9);
    assert!(v > 1.75, "volume {v}");
}

/// The same coplanarity, but the prism stops **before** the face its apex would pierce — so there
/// is no arrangement vertex on the shared line and nothing to name twice. Built before, builds now.
#[test]
fn a_prism_that_stops_short_of_the_crossing_is_unchanged() {
    let v = fused_volume(1.0, 1.0, 0.4);
    assert!(v > 1.75, "volume {v}");
}

/// ★★ **The lock that matters most.** Put the apex on `x = 0.5` instead and the staircase really
/// does have a face there covering the point: the fused surface touches itself along the apex line.
/// That body cannot be built, and this fix must not have widened its way into accepting it —
/// opening a valid model is worth nothing if it also lets an impossible one through.
#[test]
fn a_genuine_self_touch_on_a_shared_line_is_still_refused() {
    let (mut m, a, b) = stair_and_prism(1.0, 0.5, 0.9);
    assert_eq!(
        boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err(),
        nacre_ops::BoolError::Rejected {
            reason: RejectReason::SelfTouchingResult
        },
        "the surface meets itself along the apex line"
    );
}

/// And an apex on no plane of the staircase at all — the ordinary case, which must be untouched by
/// any of this.
#[test]
fn an_apex_on_no_shared_plane_is_unchanged() {
    let v = fused_volume(1.0, 0.7, 0.9);
    assert!(v > 1.75, "volume {v}");
}
