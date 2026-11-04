//! **A solid whose surface touches itself is not a solid.**
//!
//! A boolean can return a body that is connected, has the right volume and passes every topology
//! count while two pieces of its own boundary lie on top of each other. Nothing downstream can
//! represent that body — moving it, cutting it again, meshing it all assume an embedded surface —
//! so the kernel refuses it by name instead of handing it back.
//!
//! The everyday way in is a dimension that makes two surfaces meet exactly: a pocket whose wedge
//! tip lands on the far wall. Parasolid fails the same bodies (`PK_FACE_state_bad_face_face_c`);
//! OCCT accepts them, and this kernel follows Parasolid. The cost is that a plain modelling move is
//! refused — and the reject says only *what* is wrong, since which coincidence was unintended is
//! design intent the kernel cannot see.
//!
//! The propositions here are paired: every fixture that must reject sits next to one that must
//! build, because a check that rejects the wedge by rejecting *everything* would pass half of them.
//! The lateral arms' ruling pairs — an edge along a ruling inside the solid's own lateral, beside
//! two bodies touching on that line — run over both frames with the other ruling families
//! (`coverage/edge_on_a_ruling.rs`, `a_body_touching_itself_along_a_ruling_is_refused`).

use crate::common::{datum_frame, rot30, xf};
use nacre_exact::{Angle, Axis, Isometry, Rat};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, Operation, Profile2d, RejectClass, RejectReason, SketchFrame, SketchPlane,
    apply, boolean,
};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn prism(m: &mut Model, pts: &[[f64; 2]], h: f64) -> Handle<Solid> {
    let frame = SketchFrame::world(m, Axis::Z);
    prism_on(m, frame, pts, h)
}

/// The prism over `pts` on `frame`, `h` along its normal (a negative `h` sweeps against it).
fn prism_on(m: &mut Model, frame: SketchFrame, pts: &[[f64; 2]], h: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    extrude_on(m, frame, profile, h)
}

/// `profile` on `frame`, swept `h` along its normal.
fn extrude_on(m: &mut Model, frame: SketchFrame, profile: Profile2d, h: f64) -> Handle<Solid> {
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("extrude") else {
        panic!()
    };
    solid
}

/// The frame through the origin spanned by `u` and `v` (rational decimals, orthonormal).
fn frame_on(m: &mut Model, u: [f64; 3], v: [f64; 3]) -> SketchFrame {
    datum_frame(
        m,
        SketchPlane::from_axes(
            Point3::from_array([0.0; 3]),
            Vector3::from_array(u),
            Vector3::from_array(v),
        ),
    )
}

/// **A frame slanted against all three world axes** — `u = (0.6, 0.64, 0.48)`,
/// `v = (−0.8, 0.48, 0.36)`, normal `(0, −0.6, 0.8)`: every edge of a prism on it climbs in `z`, so
/// a prism over the first quadrant's corner meets the plane `z = 0` at that one corner.
fn slanted(m: &mut Model) -> SketchFrame {
    frame_on(m, [0.6, 0.64, 0.48], [-0.8, 0.48, 0.36])
}

/// **A frame tipped about the world `x` axis** — `u = (1, 0, 0)`, `v = (0, 0.8, 0.6)`: a prism
/// over the band `y ≥ 0` meets the plane `z = 0` along its bottom edge on the `x` axis.
fn tipped(m: &mut Model) -> SketchFrame {
    frame_on(m, [1.0, 0.0, 0.0], [0.0, 0.8, 0.6])
}

/// `frame` moved `dist` along its normal.
fn offset_of(m: &mut Model, frame: SketchFrame, dist: f64) -> SketchFrame {
    let Ok(nacre_ops::OpOutput::DatumPlane { frame, .. }) = apply(
        m,
        &Operation::DatumPlane {
            def: nacre_ops::DatumDef::Offset { frame, dist },
        },
    ) else {
        panic!("offsetting a frame")
    };
    frame
}

/// A box over `rect` on `frame`, from 5 below it up to it, stated by a sweep against the normal.
fn block_below(m: &mut Model, frame: SketchFrame, rect: [[f64; 2]; 4]) -> Handle<Solid> {
    prism_on(m, frame, &rect, -5.0)
}

/// The world frame or the 3-4-5 one through the origin, with its own axes in the world —
/// `u`, `v` and the normal `n = u × v`.
fn base_frame(m: &mut Model, tilted: bool) -> (SketchFrame, [[f64; 3]; 3]) {
    if tilted {
        (
            crate::common::pythagorean_frame(m, Point3::from_array([0.0; 3])),
            [[0.6, 0.8, 0.0], [-0.48, 0.36, 0.8], [0.64, -0.48, 0.6]],
        )
    } else {
        (
            SketchFrame::world(m, Axis::Z),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
    }
}

/// A unit cube, and a wedge sunk into it whose sharp tip reaches `tip_x`. At `tip_x == 1.0` the tip
/// lands exactly on the cube's far wall and the cut's cavity touches the outer shell along a line;
/// short of that it is an ordinary enclosed cavity.
fn cube_and_wedge(tip_x: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [1.0, 1.0, 1.0]);
    let b = prism(&mut m, &[[tip_x, 0.5], [0.1, 0.1], [0.1, 0.9]], 0.6);
    let b = xf(
        &mut m,
        b,
        Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::new(2, 10).unwrap()]),
    );
    (m, a, b)
}

fn expect_self_touch(m: &mut Model, a: Handle<Solid>, b: Handle<Solid>, what: &str) {
    let live = m.live_solids().to_vec();
    let err = boolean(m, BoolKind::Cut, a, b).expect_err(what);
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "{what}: got {err:?}"
    );
    assert_eq!(
        RejectReason::SelfTouchingResult.class(),
        RejectClass::Impossible,
        "a body cannot be made valid by a later capability"
    );
    // This reject is raised on a fully assembled result, so it could hand back a model with its
    // operands consumed. It must not.
    assert_eq!(m.live_solids(), live, "a reject retired the operands");
}

/// **The case the reject is for.** The wedge's tip lands on `x = 1`; the cut leaves one body
/// whose cavity meets its outer shell along the line `(1, 0.5, 0.2)–(1, 0.5, 0.8)`.
///
/// Without this check the kernel would return it: one solid, volume 0.784, `validate` silent — the
/// contact is invisible to every count, because the wall is *not split* at the line and so every
/// edge is used exactly twice (measured).
#[test]
fn a_wedge_whose_tip_reaches_the_wall_is_not_a_solid() {
    let (mut m, a, b) = cube_and_wedge(1.0);
    expect_self_touch(&mut m, a, b, "the wedge tip on the wall self-touches");
}

/// **The same contact, taken all the way through — and the endpoints stop being the story.**
///
/// Extrude the same wedge through the whole cube instead of stopping inside it. The contact line
/// now spans the wall's full height, so **both of its ends land on the wall's ring** rather than
/// inside it, and a test that asks about endpoints reads "on the boundary, not inside" twice and
/// lets the body through. It would: `Ok(n=1)`, volume `0.64` exactly, `validate` silent (measured)
/// — the wall is not split at the line, so every edge is used twice and no count sees anything.
///
/// The wall's ring is a hexagon visiting both contact points **non-adjacently**; the contact line
/// is its **chord**. Nothing about that is exotic — a through slot is a more ordinary modelling
/// move than a blind pocket, which is why the chord is checked rather than left as a documented
/// blind spot.
#[test]
fn a_wedge_cut_through_the_whole_block_is_not_a_solid() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [1.0, 1.0, 1.0]);
    // Not lifted, and as tall as the cube: the cut goes clean through.
    let b = prism(&mut m, &[[1.0, 0.5], [0.1, 0.1], [0.1, 0.9]], 1.0);
    let live = m.live_solids().to_vec();
    assert!(
        matches!(
            boolean(&mut m, BoolKind::Cut, a, b).unwrap_err(),
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "the through cut leaves the wall touching the pocket along a chord"
    );
    assert_eq!(m.live_solids(), live, "a reject retired the operands");
}

/// The same through cut with the tip off the wall: an ordinary through slot, and it must build.
/// Paired with the test above so the reject is known to be about the coincidence rather than about
/// cutting all the way through.
#[test]
fn a_through_cut_that_stops_short_of_the_wall_builds() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [1.0, 1.0, 1.0]);
    let b = prism(&mut m, &[[0.999, 0.5], [0.1, 0.1], [0.1, 0.9]], 1.0);
    let got = boolean(&mut m, BoolKind::Cut, a, b).expect("the nudged through cut builds");
    assert_eq!(got.len(), 1);
    assert_eq!(
        m.solid(got[0]).cavities.len(),
        0,
        "a through slot opens both caps, so there is no enclosed cavity"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **The same body, turned.** If the check leant on axis-aligned coordinates it would pass here
/// and quietly stop working for every real model; the judgement is on plane triples, so it does
/// not. The rotation is one this kernel's sweep already exercises.
#[test]
fn the_same_contact_is_found_after_a_rotation() {
    let (mut m, a, b) = cube_and_wedge(1.0);
    let turn = || {
        Isometry::rotation(nacre_exact::Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        })
    };
    let a = xf(&mut m, a, turn());
    let b = xf(&mut m, b, turn());
    expect_self_touch(&mut m, a, b, "a rotated self-touch is still a self-touch");
}

/// The neighbouring case: take the tip off the wall and the same cut is an ordinary pocket. It is
/// here so the reject above is known to be about the coincidence and not about wedge-shaped cuts.
#[test]
fn a_wedge_that_stops_short_of_the_wall_builds() {
    let (mut m, a, b) = cube_and_wedge(0.999);
    let got =
        boolean(&mut m, BoolKind::Cut, a, b).expect("the cut with the tip off the wall builds");
    assert_eq!(got.len(), 1);
    assert_eq!(m.solid(got[0]).cavities.len(), 1, "an enclosed cavity");
    assert!(nacre_props::mass_props(&m, got[0]).expect("props").volume > 0.0);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **Two bodies touching *each other* are not a self-touch.** The proposition here is about one
/// solid's surface meeting itself; two boxes sharing exactly an edge is the other situation
/// entirely — nothing is joined, so the answer is the two bodies that went in.
///
/// The failure worth catching is the two questions **collapsing into one**: a self-touch check that
/// fired on this pair would be answering "is there a contact" when the question is "whose surface
/// is it". That failure shows as an error where there should be none, which is a stronger statement
/// than reject-name equality.
#[test]
fn a_contact_between_two_bodies_keeps_its_own_name() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [1.0, 1.0, 1.0]);
    let b = nacre_ops::fixtures::cuboid(&mut m, [1.0, 1.0, 0.0], [2.0, 2.0, 1.0]);
    let out = boolean(&mut m, BoolKind::Fuse, a, b)
        .expect("two bodies meeting each other are not one surface touching itself");
    assert_eq!(out.len(), 2);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **Two bodies touching at a point or along a line are two bodies, in either order** — the
/// twins of the self-touch pairs below. A slanted unit cube stands on the block's top corner by its
/// own lowest corner; a tipped bar lies on the block's top by its bottom edge. Each fuse is two
/// solids, `validate` clean, whatever operand comes first and whether the scene is turned.
///
/// The corner and the edge each lie on a plane of the other solid that their own name does not
/// carry — four planes through the corner, a pencil of three through the edge's line — so the
/// grouping's ray asks a point on a plane through it named without that plane. Ordering along that
/// ray takes a plane of the name that cuts the line ([`pin_on_line`]'s rule); one of the name's
/// planes holds the ray and orders nothing, and reading it would read the point as on the ring.
#[test]
fn two_bodies_touching_at_a_corner_or_an_edge_build_in_either_order() {
    type Case = (
        &'static str,
        [[f64; 2]; 4],
        fn(&mut Model) -> SketchFrame,
        [[f64; 2]; 4],
        f64,
    );
    let cases: [Case; 3] = [
        (
            "a corner on the block's top",
            [[-10.0, -8.0], [3.0, -8.0], [3.0, 4.0], [-10.0, 4.0]],
            slanted,
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            780.0 + 1.0,
        ),
        (
            "a corner on the block's corner",
            [[-10.0, -8.0], [0.0, -8.0], [0.0, 0.0], [-10.0, 0.0]],
            slanted,
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            400.0 + 1.0,
        ),
        (
            "an edge resting on the block's top",
            [[-10.0, -8.0], [6.0, -8.0], [6.0, 4.0], [-10.0, 4.0]],
            tipped,
            [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
            960.0 + 2.0,
        ),
    ];
    for (what, rect, frame, bar, want) in cases {
        for turned in [false, true] {
            for swapped in [false, true] {
                let at = format!("{what}, turned = {turned}, swapped = {swapped}");
                let mut m = Model::new();
                let world = SketchFrame::world(&m, Axis::Z);
                let mut a = block_below(&mut m, world, rect);
                let f = frame(&mut m);
                let mut b = prism_on(&mut m, f, &bar, 1.0);
                if turned {
                    a = xf(&mut m, a, rot30());
                    b = xf(&mut m, b, rot30());
                }
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                let got =
                    boolean(&mut m, BoolKind::Fuse, x, y).unwrap_or_else(|e| panic!("{at}: {e:?}"));
                assert_eq!(got.len(), 2, "{at}: two bodies");
                assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                let v: f64 = got
                    .iter()
                    .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                    .sum();
                assert!((v - want).abs() < 1e-9, "{at}: volume {v}, want {want}");
            }
        }
    }
}

/// **A body resting on itself along an edge is refused, whichever operand comes first.** A tipped
/// bar lies on the block's top by its bottom edge, and a leg of the same bar dips into the block
/// further on, so the fuse is one body whose surface touches itself along that edge; the twin
/// without the leg is two bodies (`two_bodies_touching_at_a_corner_or_an_edge_build_in_either_order`).
///
/// Three planes hold the edge's line — the bar's two and the block's top — and the ends' names
/// carry two of them, which two decided by class order. The faces asked are the ones on whose
/// plane both ends lie, measured exactly, so the block's top is asked in both orders.
#[test]
fn a_bar_resting_on_its_own_block_by_an_edge_is_refused() {
    let legged: &[[f64; 2]] = &[
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 0.5],
        [4.0, 0.5],
        [4.0, -3.0],
        [5.0, -3.0],
        [5.0, 1.0],
        [0.0, 1.0],
    ];
    for turned in [false, true] {
        for swapped in [false, true] {
            let at = format!("turned = {turned}, swapped = {swapped}");
            let mut m = Model::new();
            let world = SketchFrame::world(&m, Axis::Z);
            let mut a = block_below(
                &mut m,
                world,
                [[-10.0, -8.0], [6.0, -8.0], [6.0, 4.0], [-10.0, 4.0]],
            );
            let f = tipped(&mut m);
            let mut b = prism_on(&mut m, f, legged, 1.0);
            if turned {
                a = xf(&mut m, a, rot30());
                b = xf(&mut m, b, rot30());
            }
            let (x, y) = if swapped { (b, a) } else { (a, b) };
            match boolean(&mut m, BoolKind::Fuse, x, y) {
                Err(BoolError::Rejected {
                    reason: RejectReason::SelfTouchingResult,
                    ..
                }) => {}
                other => panic!("{at}: {other:?}"),
            }
        }
    }
}

/// **A body standing on itself at one point is refused, whichever operand comes first.** A
/// slanted prism's lowest corner rests on the block's top — inside the face, or on the face's edge
/// — and a long leg of the same prism sinks into the block further off, so the fuse is one body
/// whose surface touches itself at that corner. No edge of the body lies on the block's top there,
/// so only the point arm asks; on the face's edge the arrangement splits one of the two faces at
/// the corner and not the other, and the edge-use guard would name an open shell further on. The
/// twins without the leg are two bodies
/// (`two_bodies_touching_at_a_corner_or_an_edge_build_in_either_order`).
#[test]
fn a_body_standing_on_itself_at_a_corner_is_refused() {
    let u_shape: &[[f64; 2]] = &[
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 4.0],
        [-8.0, 4.0],
        [-8.0, 1.0],
        [0.0, 1.0],
    ];
    for (what, rect) in [
        (
            "inside the block's top",
            [[-10.0, -8.0], [3.0, -8.0], [3.0, 4.0], [-10.0, 4.0]],
        ),
        (
            "on the block's top edge",
            [[-10.0, -8.0], [0.0, -8.0], [0.0, 4.0], [-10.0, 4.0]],
        ),
    ] {
        for turned in [false, true] {
            for swapped in [false, true] {
                let at = format!("{what}, turned = {turned}, swapped = {swapped}");
                let mut m = Model::new();
                let world = SketchFrame::world(&m, Axis::Z);
                let mut a = block_below(&mut m, world, rect);
                let f = slanted(&mut m);
                let mut b = prism_on(&mut m, f, u_shape, 1.0);
                if turned {
                    a = xf(&mut m, a, rot30());
                    b = xf(&mut m, b, rot30());
                }
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                match boolean(&mut m, BoolKind::Fuse, x, y) {
                    Err(BoolError::Rejected {
                        reason: RejectReason::SelfTouchingResult,
                        ..
                    }) => {}
                    other => panic!("{at}: {other:?}"),
                }
            }
        }
    }
}

/// **A body standing on itself at one point of a face with an arc is refused**, in either operand
/// order — the notched twin of `a_body_standing_on_itself_at_a_corner_is_refused`: a cylinder
/// notches the block's top at its far corner, so the top face's outline is a mixed ring. No edge of
/// the body lies on that face at the corner, so only the point arm asks, and the corner — a
/// three-plane point with rational coordinates — is read against that ring by
/// `rational_point_in_ring`. The notch sits where no plane of the prism meets its side, so the gate
/// passes the pair. The twin without the leg is two bodies.
#[test]
fn a_body_standing_on_its_own_notched_block_at_a_corner_is_refused() {
    let square: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let u_shape: &[[f64; 2]] = &[
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 4.0],
        [-8.0, 4.0],
        [-8.0, 1.0],
        [0.0, 1.0],
    ];
    let notch = std::f64::consts::PI * 0.25 / 4.0 * 5.0; // a quarter disk of radius 0.5, 5 deep
    for (what, top, joined) in [("alone", square, false), ("with a leg", u_shape, true)] {
        for swapped in [false, true] {
            let at = format!("{what}, swapped = {swapped}");
            let mut m = Model::new();
            let world = SketchFrame::world(&m, Axis::Z);
            let block = block_below(
                &mut m,
                world,
                [[-10.0, -8.0], [3.0, -8.0], [3.0, 4.0], [-10.0, 4.0]],
            );
            let floor = offset_of(&mut m, world, -6.0);
            let disk = nacre_ops::from_paths(vec![
                nacre_ops::Ring2d::circle(Point2::from_array([3.0, 4.0]), 0.5).expect("a circle"),
            ])
            .expect("a disk")
            .remove(0);
            let cutter = extrude_on(&mut m, floor, disk, 7.0);
            let a = boolean(&mut m, BoolKind::Cut, block, cutter).expect("the notch")[0];
            let f = slanted(&mut m);
            let b = prism_on(&mut m, f, top, 1.0);
            let (x, y) = if swapped { (b, a) } else { (a, b) };
            match (joined, boolean(&mut m, BoolKind::Fuse, x, y)) {
                (
                    true,
                    Err(BoolError::Rejected {
                        reason: RejectReason::SelfTouchingResult,
                        ..
                    }),
                ) => {}
                (false, Ok(solids)) => {
                    assert_eq!(solids.len(), 2, "{at}: two bodies");
                    assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                    let vol: f64 = solids
                        .iter()
                        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                        .sum();
                    let want = 780.0 - notch + 1.0;
                    assert!((vol - want).abs() < 1e-6, "{at}: volume {vol}, want {want}");
                }
                (_, other) => panic!("{at}: {other:?}"),
            }
        }
    }
}

/// **A bar resting on a notched face of its own block is refused** — the face's outline carries the
/// arc a cylinder cut at its far corner; the bar's bottom edge lies inside the face, so the edge
/// test reads it against that ring (its arc and pierce corners) and the point arm reads the edge's
/// ends on the rational road. The cylinder sits where no slanted plane of the bar meets its side,
/// so the gate passes the pair. The twin without the leg is two bodies. On the world frame and a
/// 3-4-5 one, in either operand order.
#[test]
fn a_bar_resting_on_its_own_notched_block_is_refused() {
    let bar: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]];
    let legged: &[[f64; 2]] = &[
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 0.5],
        [4.0, 0.5],
        [4.0, -3.0],
        [5.0, -3.0],
        [5.0, 1.0],
        [0.0, 1.0],
    ];
    let notch = std::f64::consts::PI * 0.25 / 2.0 * 5.0; // half a disk of radius 0.5, 5 deep
    for (what, top, joined) in [("alone", bar, false), ("with a leg", legged, true)] {
        for tilted in [false, true] {
            for swapped in [false, true] {
                let at = format!("{what}, tilted = {tilted}, swapped = {swapped}");
                let mut m = Model::new();
                let (base, [u, v, n]) = base_frame(&mut m, tilted);
                let block = block_below(
                    &mut m,
                    base,
                    [[-10.0, -11.0], [6.0, -11.0], [6.0, 4.0], [-10.0, 4.0]],
                );
                let floor = offset_of(&mut m, base, -6.0);
                let disk = nacre_ops::from_paths(vec![
                    nacre_ops::Ring2d::circle(Point2::from_array([-8.0, -11.0]), 0.5)
                        .expect("a circle"),
                ])
                .expect("a disk")
                .remove(0);
                let cutter = extrude_on(&mut m, floor, disk, 7.0);
                let a = boolean(&mut m, BoolKind::Cut, block, cutter).expect("the notch")[0];
                // The bar's frame, tipped about the base's `u` towards its normal.
                let tip: [f64; 3] = core::array::from_fn(|k| 0.8 * v[k] + 0.6 * n[k]);
                let f = frame_on(&mut m, u, tip);
                let b = prism_on(&mut m, f, top, 1.0);
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                let got = boolean(&mut m, BoolKind::Fuse, x, y);
                match (joined, got) {
                    (
                        true,
                        Err(BoolError::Rejected {
                            reason: RejectReason::SelfTouchingResult,
                            ..
                        }),
                    ) => {}
                    (false, Ok(solids)) => {
                        assert_eq!(solids.len(), 2, "{at}: two bodies");
                        assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                        let vol: f64 = solids
                            .iter()
                            .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                            .sum();
                        let want = 1200.0 - notch + 2.0;
                        assert!((vol - want).abs() < 1e-6, "{at}: volume {vol}, want {want}");
                    }
                    (_, other) => panic!("{at}: {other:?}"),
                }
            }
        }
    }
}

/// **A turned prism standing on its own bored block at one point is refused.** The prism is
/// turned 30° about the vertical through its corner, so its planes have no world statement, and
/// its corner, resting on the block's top, is named by those planes alone: no rational coordinate
/// says whether it is in the bore's disk. The gate proved the bore's cylinder beside every one of
/// them, and a point on a plane the cylinder lies wholly beside is off the disk — so the top face
/// reads the corner as inside and the point arm refuses the body. The twin without the leg is two
/// bodies. Either operand order.
#[test]
fn a_turned_prism_standing_on_its_own_bored_block_is_refused() {
    let leg: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let u_shape: &[[f64; 2]] = &[
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 4.0],
        [-8.0, 4.0],
        [-8.0, 1.0],
        [0.0, 1.0],
    ];
    let bore = std::f64::consts::PI * 0.25 * 5.0;
    for (what, top, joined) in [("alone", leg, false), ("with a leg", u_shape, true)] {
        for swapped in [false, true] {
            let at = format!("{what}, swapped = {swapped}");
            let mut m = Model::new();
            let world = SketchFrame::world(&m, Axis::Z);
            let block = block_below(
                &mut m,
                world,
                [[-14.0, -12.0], [8.0, -12.0], [8.0, 8.0], [-14.0, 8.0]],
            );
            let floor = offset_of(&mut m, world, -6.0);
            let disk = nacre_ops::from_paths(vec![
                nacre_ops::Ring2d::circle(Point2::from_array([-12.0, 6.0]), 0.5).expect("a circle"),
            ])
            .expect("a disk")
            .remove(0);
            let cutter = extrude_on(&mut m, floor, disk, 7.0);
            let a = boolean(&mut m, BoolKind::Cut, block, cutter).expect("the bore")[0];
            let f = slanted(&mut m);
            let b = prism_on(&mut m, f, top, 1.0);
            let b = xf(&mut m, b, crate::common::rot_iso(Axis::Z, 30));
            let (x, y) = if swapped { (b, a) } else { (a, b) };
            match (joined, boolean(&mut m, BoolKind::Fuse, x, y)) {
                (
                    true,
                    Err(BoolError::Rejected {
                        reason: RejectReason::SelfTouchingResult,
                        ..
                    }),
                ) => {}
                (false, Ok(solids)) => {
                    assert_eq!(solids.len(), 2, "{at}: two bodies");
                    assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                    let vol: f64 = solids
                        .iter()
                        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                        .sum();
                    let want = 2200.0 - bore + 1.0;
                    assert!((vol - want).abs() < 1e-6, "{at}: volume {vol}, want {want}");
                }
                (_, other) => panic!("{at}: {other:?}"),
            }
        }
    }
}

/// **A pin standing in a blind round hole, its top flush with the plate's top, is one body.**
/// The pin's top edges lie on the plate's top plane but over the hole, where that face has no
/// material: the face's bounds are its square and the hole's circle, and both are read — the
/// circle as the two points where the edge's line crosses it. Beside it, a bar standing in the
/// hole that crosses its wall builds as it does without the circle. On the world frame and on a
/// 3-4-5 frame, in either operand order.
#[test]
fn a_pin_flush_in_a_round_hole_is_one_body() {
    let pin: &[[f64; 2]] = &[[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]];
    let bar: &[[f64; 2]] = &[[-0.2, -0.2], [0.98, -0.2], [0.98, 0.2], [-0.2, 0.2]];
    let quarter = std::f64::consts::PI / 2.0; // the hole: radius 1, depth 0.5
    for (what, top, want) in [
        ("a pin flush in the hole", pin, 16.0 - quarter + 0.08),
        ("a bar across the hole's wall", bar, 14.665204),
    ] {
        for tilted in [false, true] {
            for swapped in [false, true] {
                let at = format!("{what}, tilted = {tilted}, swapped = {swapped}");
                let mut m = Model::new();
                let base = if tilted {
                    crate::common::pythagorean_frame(&mut m, Point3::from_array([0.0; 3]))
                } else {
                    SketchFrame::world(&m, Axis::Z)
                };
                let plate = prism_on(
                    &mut m,
                    base,
                    &[[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0], [-2.0, 2.0]],
                    1.0,
                );
                let floor = offset_of(&mut m, base, 0.5);
                let hole = extrude_on(&mut m, floor, crate::common::unit_disk(), 1.0);
                let cup = boolean(&mut m, BoolKind::Cut, plate, hole).expect("the cup")[0];
                let floor = offset_of(&mut m, base, 0.5);
                let b = prism_on(&mut m, floor, top, 0.5);
                let (x, y) = if swapped { (b, cup) } else { (cup, b) };
                let got =
                    boolean(&mut m, BoolKind::Fuse, x, y).unwrap_or_else(|e| panic!("{at}: {e:?}"));
                assert_eq!(got.len(), 1, "{at}: one body");
                assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                let v = nacre_props::mass_props(&m, got[0]).expect("props").volume;
                assert!((v - want).abs() < 1e-6, "{at}: volume {v}, want {want}");
                let t = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
                    .unwrap_or_else(|e| panic!("{at}: the mesher refuses it: {e:?}"));
                assert_eq!(crate::common::non_watertight(&t), 0, "{at}: the mesh leaks");
            }
        }
    }
}

/// **An *area* contact opens the wall instead of leaving a zero-thickness one.**
///
/// This check looks for an edge inside a face, which is a line contact; two faces of one solid
/// overlapping *in area* could in principle escape it, since the boundaries might coincide with no
/// edge strictly interior to anything. So the question is put to the kernel rather than reasoned
/// about: give the wedge a **flat** tip lying exactly on the far wall.
///
/// The answer is that the situation does not arise. The cut comes back as one body with **no
/// cavity** — the pocket opens through the wall, which is the regularized difference — and the
/// volume is exactly the set difference. There is no zero-thickness wall to detect. If that ever
/// changes, this test says so before a silently self-touching body ships.
///
/// Scope, honestly: one construction, not a proof about every area contact.
#[test]
fn an_area_contact_opens_the_wall_rather_than_touching_it() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [1.0, 1.0, 1.0]);
    // The right edge of this quadrilateral lies on x = 1, the cube's wall, over y in [0.4, 0.6].
    let b = prism(
        &mut m,
        &[[1.0, 0.4], [1.0, 0.6], [0.1, 0.9], [0.1, 0.1]],
        0.6,
    );
    let b = xf(
        &mut m,
        b,
        Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::new(2, 10).unwrap()]),
    );
    let got = boolean(&mut m, BoolKind::Cut, a, b).expect("the flat-tipped cut builds");
    assert_eq!(got.len(), 1);
    assert_eq!(
        m.solid(got[0]).cavities.len(),
        0,
        "the pocket opened through the wall, so there is no enclosed cavity to touch it"
    );
    let v = nacre_props::mass_props(&m, got[0]).expect("props").volume;
    assert!(
        (v - 0.73).abs() < 1e-9,
        "volume {v}, expected the set difference"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// The plainest thing a boolean does, next to the rejects: if this moved, the check is reading
/// ordinary tangential contact — which every fuse of two touching boxes has — as a self-touch.
#[test]
fn an_ordinary_fuse_is_untouched() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [2.0, 2.0, 2.0]);
    let b = nacre_ops::fixtures::cuboid(&mut m, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
    let got = boolean(&mut m, BoolKind::Fuse, a, b).expect("the fuse builds");
    assert_eq!(got.len(), 1);
    let v = nacre_props::mass_props(&m, got[0]).expect("props").volume;
    assert!((v - (8.0 + 8.0 - 1.0)).abs() < 1e-9, "volume {v}");
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **A face with a hole, so the even-odd rule is exercised and not merely present.**
///
/// The judge folds *every* ring of a face into one bag and reads even-odd, because a point inside a
/// hole has crossed the boundary twice. The rule needs a fixture rather than a corpus that happens
/// to cover it: the census's edges reach the exact test against 214 faces, none of them holed
/// (measured).
///
/// A plate with a square hole, and a wedge cut whose tip lands exactly on the hole's wall: the
/// contact line lies on a plane the plate carries, and the point is on the **hole's** ring, not the
/// outer one.
#[test]
fn a_contact_on_a_holes_wall_is_seen() {
    let mut m = Model::new();
    // Plate 3×3×1 with a 1×1 hole in the middle.
    let pts =
        |v: &[[f64; 2]]| -> Vec<Point2> { v.iter().map(|&p| Point2::from_array(p)).collect() };
    let profile = Profile2d::with_holes(
        pts(&[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]]),
        vec![pts(&[[1.0, 1.0], [1.0, 2.0], [2.0, 2.0], [2.0, 1.0]])],
    )
    .expect("a plate with a square hole");
    let frame = SketchFrame::world(&m, Axis::Z);
    let nacre_ops::OpOutput::Extrude { solid: a, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    )
    .expect("extrude") else {
        panic!()
    };
    // A wedge in the plate's material whose tip lands on the hole's wall x = 1.
    let b = prism(&mut m, &[[1.0, 1.5], [0.2, 1.1], [0.2, 1.9]], 1.0);
    assert!(
        matches!(
            boolean(&mut m, BoolKind::Cut, a, b).unwrap_err(),
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "the cut's tip touches the hole's wall"
    );
}

// ---------------------------------------------------------------------------------------------
// Chords across a face whose ring carries an arc or a pierce corner.
// ---------------------------------------------------------------------------------------------

/// How the block's top face is trimmed.
#[derive(Clone, Copy, Debug)]
enum Trim {
    /// A plain rectangle `[-10, 6] × [-11, 4]`.
    Plain,
    /// The rectangle with a half-disk notch (r 0.5) cut at `(-5, -11)` through its height.
    Notch,
    /// The same half disk cut only 4 deep — a blind groove, so the front wall `y = −11` stays one
    /// face with the groove's two rulings in its ring.
    Groove,
    /// The rectangle `[-10, 6] × [-11, top]` with every corner filleted (r 0.5).
    Fillet { top: f64 },
}

/// A block 6 deep under the base frame's plane, trimmed by `trim`.
fn trimmed_block(m: &mut Model, base: SketchFrame, trim: Trim) -> Handle<Solid> {
    match trim {
        Trim::Plain | Trim::Notch | Trim::Groove => {
            let block = prism_on(
                m,
                base,
                &[[-10.0, -11.0], [6.0, -11.0], [6.0, 4.0], [-10.0, 4.0]],
                -6.0,
            );
            let depth = match trim {
                Trim::Notch => 7.0,
                Trim::Groove => 4.0,
                Trim::Plain => return block,
                Trim::Fillet { .. } => unreachable!("the fillet has its own arm"),
            };
            let floor = offset_of(m, base, -depth);
            let disk = nacre_ops::from_paths(vec![
                nacre_ops::Ring2d::circle(Point2::from_array([-5.0, -11.0]), 0.5)
                    .expect("a circle"),
            ])
            .expect("a disk")
            .remove(0);
            let cutter = extrude_on(m, floor, disk, depth + 1.0);
            boolean(m, BoolKind::Cut, block, cutter).expect("the notch")[0]
        }
        Trim::Fillet { top } => {
            use crate::stated::{arc_turns, line, stated};
            let p = |x: f64, y: f64| Point2::from_array([x, y]);
            let r = 0.5;
            let (x0, x1, y0, y1) = (-10.0, 6.0, -11.0, top);
            let edges = vec![
                line(p(x0 + r, y0), p(x1 - r, y0)),
                arc_turns(p(x1 - r, y0 + r), p(x1 - r, y0), 1),
                line(p(x1, y0 + r), p(x1, y1 - r)),
                arc_turns(p(x1 - r, y1 - r), p(x1, y1 - r), 1),
                line(p(x1 - r, y1), p(x0 + r, y1)),
                arc_turns(p(x0 + r, y1 - r), p(x0 + r, y1), 1),
                line(p(x0, y1 - r), p(x0, y0 + r)),
                arc_turns(p(x0 + r, y0 + r), p(x0, y0 + r), 1),
            ];
            let profile = stated(edges).expect("the filleted rectangle").remove(0);
            extrude_on(m, base, profile, -6.0)
        }
    }
}

/// `body` fused with a step of height 3 over each of `rects`.
fn with_steps(
    m: &mut Model,
    base: SketchFrame,
    mut body: Handle<Solid>,
    rects: &[[[f64; 2]; 2]],
) -> Handle<Solid> {
    for &[[x0, y0], [x1, y1]] in rects {
        let step = prism_on(m, base, &[[x0, y0], [x1, y0], [x1, y1], [x0, y1]], 3.0);
        body = boolean(m, BoolKind::Fuse, body, step).expect("a step")[0];
    }
    body
}

/// `block` fused onto a plate `[-12, 8] × [-16, 6]` 2 thick under it (from 8 below the base frame
/// to 6 below), and the plate's underside frame, for what stands on the plate beside the block.
fn on_a_plate(
    m: &mut Model,
    base: SketchFrame,
    block: Handle<Solid>,
) -> (SketchFrame, Handle<Solid>) {
    let below = offset_of(m, base, -8.0);
    let plate = prism_on(
        m,
        below,
        &[[-12.0, -16.0], [8.0, -16.0], [8.0, 6.0], [-12.0, 6.0]],
        2.0,
    );
    let body = boolean(m, BoolKind::Fuse, block, plate).expect("the plate")[0];
    (below, body)
}

/// The bar of [`a_bar_resting_on_its_own_notched_block_is_refused`]'s frame: a unit-square
/// section over `[x0, x1]` on the frame tipped about the base's `u` towards its normal, its bottom
/// edge on the base plane along the base's `x` axis when `lift` is 0 (and `lift` up the tipped
/// frame's `v` otherwise).
fn tipped_bar(m: &mut Model, axes: [[f64; 3]; 3], [x0, x1]: [f64; 2], lift: f64) -> Handle<Solid> {
    let [u, v, n] = axes;
    let tip: [f64; 3] = core::array::from_fn(|k| 0.8 * v[k] + 0.6 * n[k]);
    let f = frame_on(m, u, tip);
    prism_on(
        m,
        f,
        &[[x0, lift], [x1, lift], [x1, lift + 1.0], [x0, lift + 1.0]],
        1.0,
    )
}

/// A wedge along the world `x` over `[-8, 4]` whose section in `(y, z)` is the triangle
/// `(0, 0), (0.8, 0.6), (0.6, 0.8)` raised `lift`: both of its faces through the bottom edge
/// lean towards `+y`, so on the far side of the edge both stand above the plane `z = lift`.
fn leaning_wedge(m: &mut Model, lift: f64) -> Handle<Solid> {
    let yz = frame_on(m, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let start = offset_of(m, yz, -8.0);
    prism_on(
        m,
        start,
        &[[0.0, lift], [0.8, 0.6 + lift], [0.6, 0.8 + lift]],
        12.0,
    )
}

/// What a fuse answers.
#[derive(Clone, Copy, Debug)]
enum Want {
    /// This many solids of this total volume, `validate` clean.
    Built(usize, f64),
    /// Refused by the self-touch check, its witness the chord `(−7, 0, 0)`–`(3, 0, 0)` on the
    /// world frame (on the 3-4-5 frame the reason alone is asserted).
    Touches,
    /// Refused by this name.
    Refused(RejectReason),
}

/// Fuses `a` and `b` and holds the answer against `want`.
fn expect_fuse(
    m: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    want: Want,
    world: bool,
    at: &str,
) {
    let got = boolean(m, BoolKind::Fuse, a, b);
    match (want, got) {
        (Want::Built(n, vol), Ok(solids)) => {
            assert_eq!(solids.len(), n, "{at}");
            assert!(nacre_validate::validate(m).is_empty(), "{at}");
            let got: f64 = solids
                .iter()
                .map(|&s| nacre_props::mass_props(m, s).expect("props").volume)
                .sum();
            assert!((got - vol).abs() < 1e-6, "{at}: volume {got}, want {vol}");
        }
        (
            Want::Touches,
            Err(BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                at: witness,
            }),
        ) => {
            if world {
                let Some(nacre_ops::RejectWhere::Segment(ends)) = witness else {
                    panic!("{at}: witness {witness:?}");
                };
                let mut xs = ends.map(|p| p.as_array());
                xs.sort_by(|p, q| p[0].total_cmp(&q[0]));
                let want = [[-7.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
                for (g, w) in xs.iter().zip(&want) {
                    assert!(
                        (0..3).all(|k| (g[k] - w[k]).abs() < 1e-9),
                        "{at}: witness {xs:?}"
                    );
                }
            }
        }
        (Want::Refused(r), Err(BoolError::Rejected { reason, .. })) if reason == r => {}
        (want, got) => panic!("{at}: want {want:?}, got {got:?}"),
    }
}

/// **A bar lying along a chord of its own block's top face is refused, whatever trims the
/// face** — the self-touch edge test against a face whose ring a cylinder cut. The bar is tipped
/// about the base's `x` axis so its bottom edge lies on the top face along `y = 0`, and two steps
/// on the block bury its two ends, so the body is one and its surface meets itself along the chord
/// from `x = −7` to `x = 3`. The top face is plain, or carries a half-disk notch cut in its far
/// side, or has filleted corners — the last two put an arc and pierce corners in its ring, which
/// the edge test reads with the whole cylinder table (`segment_meets_face`). Every one is refused
/// with the chord as its witness; the bar raised off the face is the twin that builds. On the
/// world frame and a 3-4-5 one, in either operand order.
///
/// Around it: the chord left open at one end (one step only) is the same refusal; a tipped bar
/// with no steps, whose edge line runs on to meet the fillet circles, is two bodies touching along
/// that line — its slanted faces stand apart from the fillets and their classes draw the fillets'
/// sections, arcs of ellipses — and, turned 30° about `z`, never reaches the check
/// (`CylinderGateUndecided`: a turned plane has no world statement); a wedge whose two faces through
/// the edge both lean away from the fillet circles its line crosses off their arcs passes the
/// gate, the circles count no crossing there (`arc_crossings`), and it is refused like the bar
/// (its raised twin builds); and a 3-4-5 column standing on a plate against the notched
/// block's front wall — its corner edge on the wall, a line parallel to the cylinder's axis — is
/// refused by the edge test, or, its edge on the notch's own ruling, by the non-manifold edge
/// check.
#[test]
fn a_chord_across_a_face_trimmed_by_a_cylinder_is_refused() {
    let full = [[[3.0, -11.0], [6.0, 4.0]], [[-10.0, -11.0], [-7.0, 4.0]]];
    let partial = [[[3.0, -8.0], [6.0, 4.0]], [[-10.0, -8.0], [-7.0, 4.0]]];
    let notch = std::f64::consts::PI * 0.25 / 2.0 * 6.0;
    let fillets = |depth: f64| 4.0 * (0.25 - std::f64::consts::PI * 0.25 / 4.0) * depth;
    for tilted in [false, true] {
        for (trim, steps, vol) in [
            (Trim::Plain, &full[..], 1720.0),
            (Trim::Notch, &full[..], 1720.0 - notch),
            (
                Trim::Fillet { top: 7.0 },
                &partial[..],
                1728.0 - fillets(6.0) + 216.0 + 10.0,
            ),
        ] {
            for lift in [0.0, 0.5] {
                for swapped in [false, true] {
                    let at =
                        format!("{trim:?}, tilted = {tilted}, lift = {lift}, swapped = {swapped}");
                    let want = if lift == 0.0 {
                        Want::Touches
                    } else {
                        Want::Built(1, vol)
                    };
                    let mut m = Model::new();
                    let (base, axes) = base_frame(&mut m, tilted);
                    let block = trimmed_block(&mut m, base, trim);
                    let a = with_steps(&mut m, base, block, steps);
                    let b = tipped_bar(&mut m, axes, [-8.0, 4.0], lift);
                    let (x, y) = if swapped { (b, a) } else { (a, b) };
                    expect_fuse(&mut m, x, y, want, !tilted, &at);
                }
            }
        }
    }
    for trim in [Trim::Plain, Trim::Notch] {
        let want = Want::Refused(RejectReason::SelfTouchingResult);
        let mut m = Model::new();
        let (base, axes) = base_frame(&mut m, false);
        let block = trimmed_block(&mut m, base, trim);
        let a = with_steps(&mut m, base, block, &full[..1]);
        let b = tipped_bar(&mut m, axes, [-12.0, 4.0], 0.0);
        expect_fuse(&mut m, a, b, want, true, &format!("one step, {trim:?}"));
    }
    // The block 16 × 11.25 × 6 less its fillets, and the bar 12 × 1 × 1 beside it on the top
    // face's line `y = 0`: two bodies touching along a line, as on the plain block.
    for (turned, want) in [
        (false, Want::Built(2, 1080.0 - fillets(6.0) + 12.0)),
        (true, Want::Refused(RejectReason::CylinderGateUndecided)),
    ] {
        let mut m = Model::new();
        let (base, axes) = base_frame(&mut m, false);
        let a = trimmed_block(&mut m, base, Trim::Fillet { top: 0.25 });
        let mut b = tipped_bar(&mut m, axes, [-8.0, 4.0], 0.0);
        if turned {
            b = xf(&mut m, b, rot30());
        }
        let at = format!("an edge line across a fillet arc, turned = {turned}");
        expect_fuse(&mut m, a, b, want, true, &at);
    }
    for lift in [0.0, 0.5] {
        let mut m = Model::new();
        let (base, _) = base_frame(&mut m, false);
        let block = trimmed_block(&mut m, base, Trim::Fillet { top: 0.75 });
        let a = with_steps(
            &mut m,
            base,
            block,
            &[[[3.0, -8.0], [4.5, 1.0]], [[-8.5, -8.0], [-7.0, 1.0]]],
        );
        let b = leaning_wedge(&mut m, lift);
        // The block 16 × 11.75 × 6 less its fillets, two steps 1.5 × 9 × 3, and the wedge's
        // section (0.14) over the 10 between them.
        let vol = 1128.0 - fillets(6.0) + 81.0 + 1.4;
        let at = format!(
            "a wedge whose edge line crosses the fillet circles off their arcs, lift = {lift}"
        );
        let want = if lift == 0.0 {
            Want::Touches
        } else {
            Want::Built(1, vol)
        };
        expect_fuse(&mut m, a, b, want, true, &at);
    }
    for (x0, reason) in [
        (-2.0, RejectReason::SelfTouchingResult),
        (-4.5, RejectReason::NonManifoldResultEdge),
    ] {
        let mut m = Model::new();
        let (base, _) = base_frame(&mut m, false);
        let block = trimmed_block(&mut m, base, Trim::Notch);
        let (below, a) = on_a_plate(&mut m, base, block);
        // A 3-4-5 square column, its corner edge on the front wall `y = −11` at `x0`.
        let ring = [
            [x0, -11.0],
            [x0 - 1.6, -12.2],
            [x0 - 0.4, -13.8],
            [x0 + 1.2, -12.6],
        ];
        let b = prism_on(&mut m, below, &ring, 8.0);
        let at = format!("a column against the front wall at x = {x0}");
        expect_fuse(&mut m, a, b, Want::Refused(reason), true, &at);
    }
}

/// **An edge at a groove's mouth builds** — the chord test (`segment_meets_face`) reading a
/// crossing on a ruling and an arc's departure.
///
/// A column stands on the plate in front of the groove's mouth, narrower than the mouth, its back
/// face on the front wall's plane `y = −11`; its top-back edge lies on a plane face of the same
/// body, so the edge test reads that face's ring along the edge's line:
/// - **against the wall, top at `z = −3`** — the groove is blind (`Trim::Groove`, `z ∈ [−4, 0]`)
///   and the column stands against the wall below it: the line crosses the groove's two rulings
///   in the wall's ring (`crossing_on_ruling`);
/// - **a bracket** — a column off the wall whose arm reaches the mouth at `z ∈ [−3.5, −3]`, touching
///   the wall nowhere: the same crossings with the wall's ring untouched by the bracket;
/// - **flush with the top, at the through notch** — the line `y = −11` on the top face holds both
///   ends of the notch's arc, which leaves it between them (`arc_departure`).
///
/// Each on the world frame and a 3-4-5 one, in either operand order.
///
/// The groove is blind for the first two because a through groove splits the front wall into two
/// faces either side of it, and an edge in front of the mouth overlaps neither face's box — the
/// sieve drops it before the ring is read.
///
/// Measured by breaking each arm in the chord test alone: a ruling crossing left unnamed would
/// refuse the column and the bracket (`RingNaming`); named at the other ruling's root, the column
/// would be refused as touching itself; an arc's departure left unread would refuse the flush
/// column (`PierceVertexUnnamed`), and turned the other way, it would be refused as touching
/// itself. The rows hold only the answer: a sieve that stopped asking these faces would leave them
/// green and the arms unread, so a change to the sieve measures this again by breaking the arms.
///
/// No row here touches: for the line to cross a ruling, the edge's faces must be one ⊥ to the axis
/// and one on the wall's plane facing the other way, and that face fits only where the wall is not,
/// so the edge never runs over the wall's inside. The wedge below has both faces slanted to the
/// axis: it rests on the wall by its edge across the groove and stays a second body.
#[test]
fn an_edge_at_a_grooves_mouth_builds() {
    let disk = std::f64::consts::PI * 0.25 / 2.0;
    let mouth = [[-5.3, -12.0], [-4.7, -12.0], [-4.7, -11.0], [-5.3, -11.0]];
    // The block 16 × 15 × 6 less the groove's half disk over its depth, the plate 20 × 22 × 2, and
    // what stands on it above the plate's top.
    for (row, trim, vol) in [
        (
            "against the wall",
            Trim::Groove,
            1440.0 - disk * 4.0 + 880.0 + 0.6 * 3.0,
        ),
        (
            "a bracket",
            Trim::Groove,
            1440.0 - disk * 4.0 + 880.0 + 0.6 * 3.0 + 0.6 * 2.0 * 0.5 - 0.6 * 0.5,
        ),
        (
            "flush with the top",
            Trim::Notch,
            1440.0 - disk * 6.0 + 880.0 + 0.6 * 6.0,
        ),
    ] {
        for tilted in [false, true] {
            for swapped in [false, true] {
                let mut m = Model::new();
                let (base, _) = base_frame(&mut m, tilted);
                let block = trimmed_block(&mut m, base, trim);
                let (below, a) = on_a_plate(&mut m, base, block);
                let b = match row {
                    "against the wall" => prism_on(&mut m, below, &mouth, 5.0),
                    "a bracket" => {
                        let off = [[-5.3, -13.0], [-4.7, -13.0], [-4.7, -12.0], [-5.3, -12.0]];
                        let column = prism_on(&mut m, below, &off, 5.0);
                        let level = offset_of(&mut m, base, -3.5);
                        let reach = [[-5.3, -13.0], [-4.7, -13.0], [-4.7, -11.0], [-5.3, -11.0]];
                        let arm = prism_on(&mut m, level, &reach, 0.5);
                        boolean(&mut m, BoolKind::Fuse, column, arm).expect("bracket")[0]
                    }
                    "flush with the top" => prism_on(&mut m, below, &mouth, 8.0),
                    _ => unreachable!("a row with no column"),
                };
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                let at = format!("{row}, tilted = {tilted}, swapped = {swapped}");
                expect_fuse(&mut m, x, y, Want::Built(1, vol), !tilted, &at);
            }
        }
    }
    // A wedge along the wall whose edge lies on it at `z = −3` across the groove, both faces
    // leaning off it: slanted to the axis through the groove's side, and touching the block along
    // that edge only — two bodies.
    let mut m = Model::new();
    let (base, _) = base_frame(&mut m, false);
    let block = trimmed_block(&mut m, base, Trim::Groove);
    let (_, a) = on_a_plate(&mut m, base, block);
    let yz = frame_on(&mut m, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let start = offset_of(&mut m, yz, -8.0);
    let wedge = prism_on(
        &mut m,
        start,
        &[[-11.0, -3.0], [-11.8, -2.4], [-11.6, -2.2]],
        6.0,
    );
    let vol = |m: &Model, s| nacre_props::mass_props(m, s).expect("props").volume;
    let want = Want::Built(2, vol(&m, a) + vol(&m, wedge));
    expect_fuse(&mut m, a, wedge, want, true, "a wedge leaning off the wall");
}

/// **A body whose edge touches its own cap's rim at one point is refused** — the edge test's
/// touch: a line tangent to a curved bound meets the face's closed region at a point no vertex
/// names.
///
/// A plate carries a boss of radius 5 up to 4 above it, and a column on the plate carries an arm
/// over the boss's side of it whose bottom face lies on the cap's plane and whose edge there runs
/// along the line tangent to the rim at `(3, 4)` (or at `(1.4, 4.8)`). The arm and the boss meet at
/// that one point and the body is joined through the plate: its surface touches itself there. The
/// cap is a whole disk (a circle bound), or the boss is trimmed to `x ≥ −3` and its top is a ring of
/// an arc and a chord. Each row's twin, the arm raised 0.5, is one body. The control: the boss
/// trimmed to `x ≤ 2`, whose arc stops short of the tangent point — the line touches the circle but
/// not the face, and the body builds. On the world frame and a 3-4-5 one, in either operand order.
#[test]
fn a_body_touching_its_own_cap_rim_at_a_point_is_refused() {
    // The arm's floor over the tangent line through `at`, 2 deep beyond it, and a column under it.
    type Arm = (&'static [[f64; 2]], [[f64; 2]; 4]);
    let tangent_at_3_4: Arm = (
        &[[7.0, 1.0], [8.2, 2.6], [0.2, 8.6], [-1.0, 7.0]],
        [[6.6, 2.0], [7.4, 2.0], [7.4, 2.6], [6.6, 2.6]],
    );
    let tangent_at_1_4_4_8: Arm = (
        &[[6.2, 3.4], [6.76, 5.32], [-2.84, 8.12], [-3.4, 6.2]],
        [[1.6, 5.8], [2.4, 5.8], [2.4, 6.4], [1.6, 6.4]],
    );
    // `(row, trim: the half-plane cut off the boss, arm, refused)`.
    type Row = (&'static str, Option<[[f64; 2]; 4]>, Arm, bool);
    let rows: [Row; 4] = [
        ("a whole cap, tangent at (3, 4)", None, tangent_at_3_4, true),
        (
            "a whole cap, tangent at (1.4, 4.8)",
            None,
            tangent_at_1_4_4_8,
            true,
        ),
        (
            "a D cap (x ≥ −3), tangent at (3, 4)",
            Some([[-6.0, -6.0], [-3.0, -6.0], [-3.0, 6.0], [-6.0, 6.0]]),
            tangent_at_3_4,
            true,
        ),
        (
            "a D cap (x ≤ 2), the tangent point off its arc",
            Some([[2.0, -6.0], [6.0, -6.0], [6.0, 6.0], [2.0, 6.0]]),
            tangent_at_3_4,
            false,
        ),
    ];
    for (row, trim, (arm, column), refused) in rows {
        for tilted in [false, true] {
            for lift in [0.0, 0.5] {
                for swapped in [false, true] {
                    let at =
                        format!("{row}, tilted = {tilted}, lift = {lift}, swapped = {swapped}");
                    let mut m = Model::new();
                    let (base, _) = base_frame(&mut m, tilted);
                    let floor = offset_of(&mut m, base, -2.0);
                    let plate = prism_on(
                        &mut m,
                        floor,
                        &[[-10.0, -10.0], [12.0, -10.0], [12.0, 12.0], [-10.0, 12.0]],
                        2.0,
                    );
                    let disk = nacre_ops::from_paths(vec![
                        nacre_ops::Ring2d::circle(Point2::from_array([0.0, 0.0]), 5.0)
                            .expect("a circle"),
                    ])
                    .expect("a disk")
                    .remove(0);
                    let boss = extrude_on(&mut m, base, disk, 4.0);
                    let mut a = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss")[0];
                    if let Some(half) = trim {
                        let cutter = prism_on(&mut m, base, &half, 5.0);
                        a = boolean(&mut m, BoolKind::Cut, a, cutter).expect("the D")[0];
                    }
                    let level = offset_of(&mut m, base, 4.0 + lift);
                    let floor_of_arm = prism_on(&mut m, level, arm, 2.0);
                    let post = prism_on(&mut m, base, &column, 4.0 + lift);
                    let b = boolean(&mut m, BoolKind::Fuse, floor_of_arm, post).expect("the arm");
                    assert_eq!(b.len(), 1, "{at}: the arm is one body");
                    let b = b[0];
                    let vol = |m: &Model, s| nacre_props::mass_props(m, s).expect("props").volume;
                    let want = vol(&m, a) + vol(&m, b);
                    let (x, y) = if swapped { (b, a) } else { (a, b) };
                    match (
                        refused && lift == 0.0,
                        boolean(&mut m, BoolKind::Fuse, x, y),
                    ) {
                        (
                            true,
                            Err(BoolError::Rejected {
                                reason: RejectReason::SelfTouchingResult,
                                at: witness,
                            }),
                        ) => {
                            if !tilted {
                                // The witness is the arm's edge on the tangent line.
                                let Some(nacre_ops::RejectWhere::Segment(ends)) = witness else {
                                    panic!("{at}: witness {witness:?}");
                                };
                                let line = |p: Point3| {
                                    let [x, y, z] = p.as_array();
                                    let [a, b] = [arm[0], arm[3]];
                                    let cross =
                                        (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]);
                                    cross.abs() < 1e-9 && (z - 4.0).abs() < 1e-9
                                };
                                assert!(ends.iter().all(|&p| line(p)), "{at}: witness {ends:?}");
                            }
                        }
                        (false, Ok(solids)) => {
                            assert_eq!(solids.len(), 1, "{at}: one body");
                            assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                            let got: f64 = solids.iter().map(|&s| vol(&m, s)).sum();
                            assert!((got - want).abs() < 1e-6, "{at}: volume {got}, want {want}");
                        }
                        (_, other) => panic!("{at}: {other:?}"),
                    }
                }
            }
        }
    }
}

/// **A body whose vertex on a cylinder's surface presses a face of its own is refused** — the point
/// arm asks every vertex, whatever kind its name is.
///
/// A plate carries a boss of radius 5 up to 3 above it, and a post on the plate carries a slab
/// whose underside is the plane `z = 10`. A pyramid stands on a column of its own on the plate, its
/// apex up against the slab's underside: on the boss's cylinder's surface above the boss at
/// `(3, 4)` or `(−3, 4)`, where the vertex is represented by its pierce name, or off it at
/// `(3.25, 4)`, where it keeps its three-plane name. (Every coordinate the fixture adds is exact in
/// `f64`: a datum offset off the shortest decimal is a plane with a wide name.) The body is joined through the plate, so its
/// surface touches itself at the apex either way. Lowered 0.5, each is one body. In either operand
/// order.
#[test]
fn a_body_whose_vertex_on_a_cylinder_surface_touches_its_own_face_is_refused() {
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    for (ax, ay) in [(3.0, 4.0), (-3.0, 4.0), (3.25, 4.0)] {
        for drop in [0.0, 0.5] {
            for swapped in [false, true] {
                let at = format!("apex at ({ax}, {ay}), dropped {drop}, swapped = {swapped}");
                let mut m = Model::new();
                let world = SketchFrame::world(&m, Axis::Z);
                let floor = offset_of(&mut m, world, -2.0);
                let plate = prism_on(&mut m, floor, &rect(-10.0, -10.0, 12.0, 12.0), 2.0);
                let disk = nacre_ops::from_paths(vec![
                    nacre_ops::Ring2d::circle(Point2::from_array([0.0, 0.0]), 5.0)
                        .expect("a circle"),
                ])
                .expect("a disk")
                .remove(0);
                let boss = extrude_on(&mut m, world, disk, 3.0);
                let mut a = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss")[0];
                let level = offset_of(&mut m, world, 10.0);
                let slab = prism_on(&mut m, level, &rect(-6.0, -2.0, 8.0, 8.0), 1.0);
                let post = prism_on(&mut m, world, &rect(7.0, 7.0, 8.0, 8.0), 10.0);
                let slab = boolean(&mut m, BoolKind::Fuse, slab, post).expect("the slab");
                assert_eq!(slab.len(), 1, "{at}: the slab on its post is one body");
                a = boolean(&mut m, BoolKind::Fuse, a, slab[0]).expect("the slab")[0];
                // The pyramid: two wedges' common part, apex `(ax, ay, 10 − drop)`, base 3 below.
                let (top, bot) = (10.0 - drop, 7.0 - drop);
                let yz = frame_on(&mut m, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
                let at_x = offset_of(&mut m, yz, ax - 3.5);
                let along_x = prism_on(
                    &mut m,
                    at_x,
                    &[[ay, top], [ay - 2.0, bot], [ay + 2.0, bot]],
                    7.0,
                );
                // This datum's canonical frame reads the sketch as `x = −s`, `z = t`.
                let zx = frame_on(&mut m, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
                let at_y = offset_of(&mut m, zx, ay - 3.5);
                let along_y = prism_on(
                    &mut m,
                    at_y,
                    &[[-ax, top], [-(ax + 2.0), bot], [-(ax - 2.0), bot]],
                    7.0,
                );
                let pyramid =
                    boolean(&mut m, BoolKind::Common, along_x, along_y).expect("the pyramid");
                assert_eq!(pyramid.len(), 1, "{at}: one pyramid");
                let column = prism_on(&mut m, world, &rect(ax - 1.0, ay + 1.0, ax, ay + 2.0), bot);
                let b = boolean(&mut m, BoolKind::Fuse, pyramid[0], column).expect("the column");
                assert_eq!(b.len(), 1, "{at}: the pyramid on its column is one body");
                let b = b[0];
                let vol = |m: &Model, s| nacre_props::mass_props(m, s).expect("props").volume;
                let want = vol(&m, a) + vol(&m, b);
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                match (drop == 0.0, boolean(&mut m, BoolKind::Fuse, x, y)) {
                    (
                        true,
                        Err(BoolError::Rejected {
                            reason: RejectReason::SelfTouchingResult,
                            at: Some(nacre_ops::RejectWhere::Point(p)),
                        }),
                    ) => {
                        let [x, y, z] = p.as_array();
                        assert!(
                            (x - ax).abs() < 1e-9
                                && (y - ay).abs() < 1e-9
                                && (z - 10.0).abs() < 1e-9,
                            "{at}: witness {p:?}"
                        );
                    }
                    (false, Ok(solids)) => {
                        assert_eq!(solids.len(), 1, "{at}: one body");
                        assert!(nacre_validate::validate(&m).is_empty(), "{at}");
                        let got: f64 = solids.iter().map(|&s| vol(&m, s)).sum();
                        assert!((got - want).abs() < 1e-6, "{at}: volume {got}, want {want}");
                    }
                    (_, other) => panic!("{at}: {other:?}"),
                }
            }
        }
    }
}

/// **A body resting a corner on its own rim is a pinch** — refused as `NonManifoldVertex`, the name
/// the same contact gets between flat faces.
///
/// A box's bottom face lies on a cap's plane and its corner on the cap's rim, or on a bored plate's
/// top with its corner on the hole's rim from inside; or an arm's bottom face lies on the cap's
/// plane with its corner on the rim and one side along the rim's tangent there. A column of its own
/// on the plate joins the box (or the arm) to the body elsewhere. The two faces on that plane touch at the corner only, and the plane's
/// arrangement cuts both rings there: the corner becomes a vertex of the cap (or the plate's top),
/// of the cylinder's lateral and of the box at once, so it bounds every face it touches and no
/// face-reading check can see it — it is where two fans of faces meet at one vertex. Each row's
/// twins: raised 0.5, one body; with no column, two bodies. On the world frame and a 3-4-5 one, in
/// either operand order.
#[test]
fn a_body_resting_a_corner_on_its_own_rim_is_a_pinch() {
    type Quad = [[f64; 2]; 4];
    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Quad {
        [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    }
    fn disk() -> Profile2d {
        nacre_ops::from_paths(vec![
            nacre_ops::Ring2d::circle(Point2::from_array([0.0, 0.0]), 5.0).expect("a circle"),
        ])
        .expect("a disk")
        .remove(0)
    }
    fn vol(m: &Model, s: Handle<Solid>) -> f64 {
        nacre_props::mass_props(m, s).expect("props").volume
    }
    // A plate with a boss of radius 5 up to 4 above it, and a box over `top` from the cap's plane
    // up 2 (raised `lift`), carried by a column over `post` when `joined`.
    fn on_a_boss(
        m: &mut Model,
        base: SketchFrame,
        top: Quad,
        post: Quad,
        lift: f64,
        joined: bool,
    ) -> (Handle<Solid>, Handle<Solid>) {
        let floor = offset_of(m, base, -2.0);
        let plate = prism_on(m, floor, &rect(-10.0, -10.0, 12.0, 12.0), 2.0);
        let boss = extrude_on(m, base, disk(), 4.0);
        let a = boolean(m, BoolKind::Fuse, plate, boss).expect("the boss")[0];
        let level = offset_of(m, base, 4.0 + lift);
        let mut b = prism_on(m, level, &top, 2.0);
        if joined {
            let column = prism_on(m, base, &post, 4.0 + lift);
            b = boolean(m, BoolKind::Fuse, b, column).expect("the arm")[0];
        }
        (a, b)
    }
    // A plate 5 thick with a blind hole of radius 5, 3 deep, and a box over `top` on the plate's
    // top (raised `lift`), carried by a column from the hole's floor when `joined`.
    fn over_a_hole(
        m: &mut Model,
        base: SketchFrame,
        top: Quad,
        post: Quad,
        lift: f64,
        joined: bool,
    ) -> (Handle<Solid>, Handle<Solid>) {
        let floor = offset_of(m, base, -5.0);
        let plate = prism_on(m, floor, &rect(-10.0, -10.0, 12.0, 12.0), 5.0);
        let deep = offset_of(m, base, -3.0);
        let cutter = extrude_on(m, deep, disk(), 4.0);
        let a = boolean(m, BoolKind::Cut, plate, cutter).expect("the hole")[0];
        let level = if lift == 0.0 {
            base
        } else {
            offset_of(m, base, lift)
        };
        let mut b = prism_on(m, level, &top, 2.0);
        if joined {
            let column = prism_on(m, deep, &post, 3.0 + lift);
            b = boolean(m, BoolKind::Fuse, b, column).expect("the box")[0];
        }
        (a, b)
    }
    type Build =
        fn(&mut Model, SketchFrame, Quad, Quad, f64, bool) -> (Handle<Solid>, Handle<Solid>);
    let rows: [(&str, Build, Quad, Quad, [f64; 3]); 4] = [
        (
            "a box corner on the boss's rim at (3, 4)",
            on_a_boss,
            rect(3.0, 4.0, 10.0, 10.0),
            rect(8.0, 8.0, 10.0, 10.0),
            [3.0, 4.0, 4.0],
        ),
        (
            "a box corner on the boss's rim at (4, 3)",
            on_a_boss,
            rect(4.0, 3.0, 10.0, 10.0),
            rect(8.0, 8.0, 10.0, 10.0),
            [4.0, 3.0, 4.0],
        ),
        (
            "an arm's corner on the rim's tangent point (3, 4), one side along the tangent",
            on_a_boss,
            [[3.0, 4.0], [7.0, 1.0], [8.2, 2.6], [4.2, 5.6]],
            rect(6.6, 2.0, 7.4, 2.6),
            [3.0, 4.0, 4.0],
        ),
        (
            "a box corner on the hole's rim at (3, 4), from inside",
            over_a_hole,
            rect(-1.0, 0.0, 3.0, 4.0),
            rect(0.0, 1.0, 1.0, 2.0),
            [3.0, 4.0, 0.0],
        ),
    ];
    let mut wrong = Vec::new();
    for (row, build, top, post, corner) in rows {
        for tilted in [false, true] {
            for (lift, joined) in [(0.0, true), (0.5, true), (0.0, false)] {
                for swapped in [false, true] {
                    let at = format!(
                        "{row}, tilted = {tilted}, lift = {lift}, joined = {joined}, \
                         swapped = {swapped}"
                    );
                    let mut m = Model::new();
                    let (base, _) = base_frame(&mut m, tilted);
                    let (a, b) = build(&mut m, base, top, post, lift, joined);
                    let want = vol(&m, a) + vol(&m, b);
                    let (x, y) = if swapped { (b, a) } else { (a, b) };
                    let got = boolean(&mut m, BoolKind::Fuse, x, y);
                    let pinch = lift == 0.0 && joined;
                    match (pinch, got) {
                        (
                            true,
                            Err(BoolError::Rejected {
                                reason: RejectReason::NonManifoldVertex,
                                at: Some(nacre_ops::RejectWhere::Point(p)),
                            }),
                        ) => {
                            let off = p
                                .as_array()
                                .iter()
                                .zip(corner)
                                .any(|(g, w)| (g - w).abs() > 1e-9);
                            if !tilted && off {
                                wrong.push(format!("{at}: witness {p:?}"));
                            }
                        }
                        (false, Ok(solids)) => {
                            let n = if joined { 1 } else { 2 };
                            let got: f64 = solids.iter().map(|&s| vol(&m, s)).sum();
                            if solids.len() != n
                                || !nacre_validate::validate(&m).is_empty()
                                || (got - want).abs() > 1e-6
                            {
                                wrong.push(format!(
                                    "{at}: {} bodies, volume {got} (want {want})",
                                    solids.len()
                                ));
                            }
                        }
                        (_, other) => wrong.push(format!("{at}: {other:?}")),
                    }
                }
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} cells:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// **A body touching its own lateral face at one point is refused** — the lateral arms: a vertex on
/// the side of the boss, and an edge tangent to it, of a part that joins the body through the
/// plate.
///
/// The plate `[−10, 12]² × [−2, 0]` carries a boss of radius 5 about `z` up to `z = 4`. A slanted
/// prism's corner rests on the boss's side at `(3, 4, 2)`, every edge from it leaving the cylinder
/// outward, and a column carries it to the plate; and a wedge's apex edge, along `(0, 0.8, 0.6)`
/// through `(5, 0, 2)`, is tangent to the side, the wedge carried by its own column. Each touches
/// the lateral at that point only and the body is one through the plate, so its surface touches
/// itself there. Each row's twin, the part moved off the side by 0.5, is one body. In either
/// operand order.
#[test]
fn a_body_touching_its_own_lateral_at_a_point_is_refused() {
    let stated = |m: &mut Model, o: [f64; 3], u: [f64; 3], v: [f64; 3]| {
        datum_frame(
            m,
            nacre_ops::SketchPlane::from_axes(
                Point3::from_array(o),
                Vector3::from_array(u),
                Vector3::from_array(v),
            ),
        )
    };
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    for row in ["a corner on the side", "an edge tangent to the side"] {
        for lift in [0.0, 0.5] {
            for swapped in [false, true] {
                let at = format!("{row}, lift = {lift}, swapped = {swapped}");
                let mut m = Model::new();
                let world = SketchFrame::world(&m, Axis::Z);
                let floor = offset_of(&mut m, world, -2.0);
                let plate = prism_on(&mut m, floor, &rect(-10.0, -10.0, 12.0, 12.0), 2.0);
                let disk = nacre_ops::from_paths(vec![
                    nacre_ops::Ring2d::circle(Point2::from_array([0.0, 0.0]), 5.0)
                        .expect("a circle"),
                ])
                .expect("a disk")
                .remove(0);
                let boss = extrude_on(&mut m, world, disk, 4.0);
                let a = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss")[0];
                let (part, column) = if row == "a corner on the side" {
                    let p = [3.0 + 0.6 * lift, 4.0 + 0.8 * lift, 2.0];
                    let f = stated(&mut m, p, [0.8, -0.6, 0.0], [0.36, 0.48, -0.8]);
                    let w = prism_on(&mut m, f, &[[0.0, 0.0], [1.0, 2.0], [-1.0, 2.0]], 2.0);
                    (w, prism_on(&mut m, world, &rect(3.6, 4.9, 4.4, 5.7), 1.6))
                } else {
                    let f = stated(
                        &mut m,
                        [5.0 + lift, -1.6, 0.4],
                        [1.0, 0.0, 0.0],
                        [0.0, 0.6, -0.8],
                    );
                    let w = prism_on(&mut m, f, &[[0.0, 0.0], [1.0, 0.5], [1.0, -0.5]], 4.0);
                    let col = rect(5.6 + lift, -2.2, 7.0 + lift, -0.8);
                    (w, prism_on(&mut m, world, &col, 1.0))
                };
                let b = boolean(&mut m, BoolKind::Fuse, part, column).expect("the part");
                assert_eq!(b.len(), 1, "{at}: the part is one body");
                let b = b[0];
                let vol = |m: &Model, s| nacre_props::mass_props(m, s).expect("props").volume;
                let want = if lift == 0.0 {
                    Want::Touches
                } else {
                    Want::Built(1, vol(&m, a) + vol(&m, b))
                };
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                expect_fuse(&mut m, x, y, want, false, &at);
            }
        }
    }
}

/// **A part on a plane tangent to the boss is asked where it meets the side** — the tangent
/// hand-off: a point of the lateral a tangency row's stretch holds is the tangency verdict's, and
/// any other point the part meets the side at is a touch.
///
/// The plate and boss of [`a_body_touching_its_own_lateral_at_a_point_is_refused`] (radius 5 about
/// `z`, up to `z = 4`), and a triangular prism with its face on the tangent plane `x = 5`, extruded
/// outward to `x = 6`, carried to the plate by a column or standing apart. Four ways its face meets
/// the tangent line `x = 5, y = 0`:
/// - **a corner on the line** at `(5, 0, 2)`, the face beside it: the face meets the line at that
///   point only, so no row speaks, and the joined body touches its own side there;
/// - **a face along the line** over `2 ≤ z ≤ 3.5`: a line contact inside the lateral's span, which
///   the tangency verdict judges, and the joined body pinches along it;
/// - **a corner on the top rim** at `(5, 0, 4)`, the face above it: the lateral ends there, so the
///   face meets no stretch of it;
/// - **an edge along the top rim's tangent**, the face above it: its edge meets the rim at
///   `(5, 0, 4)` only.
///
/// Joined, each touches itself; moved off the side by 0.5, or standing apart, it builds — one body
/// or two. In either operand order.
#[test]
fn a_part_on_a_tangent_plane_is_asked_where_it_meets_the_side() {
    let stated = |m: &mut Model, o: [f64; 3], u: [f64; 3], v: [f64; 3]| {
        datum_frame(
            m,
            nacre_ops::SketchPlane::from_axes(
                Point3::from_array(o),
                Vector3::from_array(u),
                Vector3::from_array(v),
            ),
        )
    };
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    // Each row: the face's frame (origin on the line, `u`, `v` — the normal `u × v` is `+x`), its
    // triangle, and the column `[x0, y0, x1, y1]` and its height.
    #[allow(clippy::type_complexity)]
    let rows: [(&str, [[f64; 3]; 3], [[f64; 2]; 3], [f64; 4], f64); 4] = [
        (
            "a corner on the line",
            [[5.0, 0.0, 2.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 0.0], [2.0, -1.0], [2.0, 1.0]],
            [5.3, 1.5, 5.9, 1.9],
            2.6,
        ),
        (
            "a face along the line",
            [[5.0, 0.0, 2.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0], [1.5, -1.0], [1.5, 1.0]],
            [5.3, 0.3, 5.9, 0.6],
            3.2,
        ),
        (
            "a corner on the top rim",
            [[5.0, 0.0, 4.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0], [1.0, -0.5], [1.0, 0.5]],
            [5.3, 0.1, 5.9, 0.3],
            4.8,
        ),
        (
            "an edge along the top rim's tangent",
            [[5.0, 0.0, 4.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
            [[0.0, -1.0], [0.0, 1.0], [0.5, 0.0]],
            [5.3, 0.3, 5.9, 0.6],
            4.1,
        ),
    ];
    for (row, [o, u, v], tri, col, h) in rows {
        for joined in [true, false] {
            for lift in [0.0, 0.5] {
                for swapped in [false, true] {
                    let at =
                        format!("{row}, joined = {joined}, lift = {lift}, swapped = {swapped}");
                    let mut m = Model::new();
                    let world = SketchFrame::world(&m, Axis::Z);
                    let floor = offset_of(&mut m, world, -2.0);
                    let plate = prism_on(&mut m, floor, &rect(-10.0, -10.0, 12.0, 12.0), 2.0);
                    let disk = nacre_ops::from_paths(vec![
                        nacre_ops::Ring2d::circle(Point2::from_array([0.0, 0.0]), 5.0)
                            .expect("a circle"),
                    ])
                    .expect("a disk")
                    .remove(0);
                    let boss = extrude_on(&mut m, world, disk, 4.0);
                    let a = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss")[0];
                    let f = stated(&mut m, [o[0] + lift, o[1], o[2]], u, v);
                    let part = prism_on(&mut m, f, &tri, 1.0);
                    let b = if joined {
                        let c = rect(col[0] + lift, col[1], col[2] + lift, col[3]);
                        let column = prism_on(&mut m, world, &c, h);
                        let b = boolean(&mut m, BoolKind::Fuse, part, column).expect("the part");
                        assert_eq!(b.len(), 1, "{at}: the part is one body");
                        b[0]
                    } else {
                        part
                    };
                    let vol = |m: &Model, s| nacre_props::mass_props(m, s).expect("props").volume;
                    let want = match (joined, lift == 0.0) {
                        (true, true) => Want::Touches,
                        (true, false) => Want::Built(1, vol(&m, a) + vol(&m, b)),
                        (false, _) => Want::Built(2, vol(&m, a) + vol(&m, b)),
                    };
                    let (x, y) = if swapped { (b, a) } else { (a, b) };
                    expect_fuse(&mut m, x, y, want, false, &at);
                }
            }
        }
    }
}
