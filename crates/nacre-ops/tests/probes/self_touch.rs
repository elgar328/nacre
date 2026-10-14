//! **A solid whose surface touches itself is not a solid.**
//!
//! A boolean can return a body that is connected, has the right volume and passes every topology
//! count while two pieces of its own boundary lie on top of each other. Nothing downstream can
//! represent that body — moving it, cutting it again, meshing it all assume an embedded surface —
//! so the kernel refuses it by name instead of handing it back.
//!
//! The everyday way in is a dimension that makes two surfaces meet exactly: a pocket whose wedge
//! tip lands on the far wall. Parasolid fails the same bodies (`PK_FACE_state_bad_face_face_c`);
//! OCCT accepts them, and this kernel follows Parasolid. The cost is that a plain modelling move
//! is now refused — and the reject says only *what* is wrong, since which coincidence was
//! unintended is design intent the kernel cannot see.
//!
//! ★ The propositions here are paired: every fixture that must reject sits next to one that must
//! build, because a check that rejects the wedge by rejecting *everything* would pass half of them.
//! The curved arm's pairs — an edge along a ruling inside the solid's own lateral, beside two
//! bodies touching on that line — run over both frames with the other ruling families
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
    m.rebuild_adjacency();
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

/// ★ **A frame slanted against all three world axes** — `u = (0.6, 0.64, 0.48)`,
/// `v = (−0.8, 0.48, 0.36)`, normal `(0, −0.6, 0.8)`: every edge of a prism on it climbs in `z`, so
/// a prism over the first quadrant's corner meets the plane `z = 0` at that one corner.
fn slanted(m: &mut Model) -> SketchFrame {
    frame_on(m, [0.6, 0.64, 0.48], [-0.8, 0.48, 0.36])
}

/// ★ **A frame tipped about the world `x` axis** — `u = (1, 0, 0)`, `v = (0, 0.8, 0.6)`: a prism
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
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
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
        "a body cannot be made valid by a later milestone"
    );
    // ★ This is the kernel's first reject raised on a fully assembled result, so it is the first
    // one that could have handed back a model with its operands consumed. It must not.
    assert_eq!(m.live_solids(), live, "a reject retired the operands");
}

/// ★ **The case the reject is for.** The wedge's tip lands on `x = 1`; the cut leaves one body
/// whose cavity meets its outer shell along the line `(1, 0.5, 0.2)–(1, 0.5, 0.8)`.
///
/// Before this check the kernel returned it: one solid, volume 0.784, `validate` silent — the
/// contact is invisible to every count, because the wall is *not split* at the line and so every
/// edge is still used exactly twice.
#[test]
fn a_wedge_whose_tip_reaches_the_wall_is_not_a_solid() {
    let (mut m, a, b) = cube_and_wedge(1.0);
    expect_self_touch(&mut m, a, b, "the wedge tip on the wall self-touches");
}

/// ★★ **The same contact, taken all the way through — and the endpoints stop being the story.**
///
/// Extrude the same wedge through the whole cube instead of stopping inside it. The contact line
/// now spans the wall's full height, so **both of its ends land on the wall's ring** rather than
/// inside it, and a test that asks about endpoints reads "on the boundary, not inside" twice and
/// lets the body through. It did: `Ok(n=1)`, volume `0.64` exactly, `validate` silent — the wall is
/// not split at the line, so every edge is still used twice and no count sees anything.
///
/// The wall's ring is a hexagon visiting both contact points **non-adjacently**; the contact line is
/// its **chord**. Nothing about that is exotic — a through slot is a more ordinary modelling move
/// than a blind pocket, which is why the blind spot was worth closing rather than documenting.
#[test]
fn a_wedge_cut_through_the_whole_block_is_not_a_solid() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
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
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let b = prism(&mut m, &[[0.999, 0.5], [0.1, 0.1], [0.1, 0.9]], 1.0);
    let got = boolean(&mut m, BoolKind::Cut, a, b).expect("the nudged through cut builds");
    m.rebuild_adjacency();
    assert_eq!(got.len(), 1);
    assert_eq!(
        m.solid(got[0]).cavities.len(),
        0,
        "a through slot opens both caps, so there is no enclosed cavity"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// ★ **The same body, turned.** If the check leant on axis-aligned coordinates it would pass here
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
    m.rebuild_adjacency();
    assert_eq!(got.len(), 1);
    assert_eq!(m.solid(got[0]).cavities.len(), 1, "an enclosed cavity");
    assert!(nacre_props::mass_props(&m, got[0]).expect("props").volume > 0.0);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// ★ **Two bodies touching *each other* are not a self-touch.** The proposition here is about one
/// solid's surface meeting itself; two boxes sharing exactly an edge is the other situation
/// entirely — nothing is joined, so the answer is the two bodies that went in.
///
/// The failure worth catching is the two questions **collapsing into one**: a self-touch check that
/// fired on this pair would be answering "is there a contact" when the question is "whose surface
/// is it". That failure now shows as an error where there should be none, which is a stronger
/// statement than the reject-name equality this used to make.
#[test]
fn a_contact_between_two_bodies_keeps_its_own_name() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, b)
        .expect("two bodies meeting each other are not one surface touching itself");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 2);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// ★ **Two bodies touching at a point or along a line are two bodies, in either order** — the
/// twins of the self-touch pairs below. A slanted unit cube stands on the block's top corner by its
/// own lowest corner; a tipped bar lies on the block's top by its bottom edge. Each fuse is two
/// solids, `validate` clean, whatever operand comes first and whether the scene is turned.
///
/// The corner and the edge each lie on a plane of the other solid that their own name does not
/// carry — four planes through the corner, a pencil of three through the edge's line — so the
/// grouping's ray asks a point on a plane through it named without that plane. Ordering along that
/// ray takes a plane of the name that cuts the line ([`pin_on_line`]'s rule); one of the name's
/// planes holds the ray and orders nothing, and reading it read the point as on the ring.
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
                m.rebuild_adjacency();
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

/// ★ **A body resting on itself along an edge is refused, whichever operand comes first.** A tipped
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

/// ★ **A body standing on itself at one point is refused, whichever operand comes first.** A
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

/// ★ **A bar resting on a notched face of its own block is refused** — the face's outline carries
/// the arc a cylinder cut at its far corner, so the edge test does not read it (a mixed ring); the
/// bar's bottom edge has its ends inside the face, and the point arm reads them on the rational
/// road. The cylinder sits where no slanted plane of the bar meets its side, so the gate passes
/// the pair. The twin without the leg is two bodies. On the world frame and a 3-4-5 one, in either
/// operand order.
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
                m.rebuild_adjacency();
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
                        m.rebuild_adjacency();
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

/// ★ **A turned prism standing on its own bored block at one point is refused.** The prism is
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
            m.rebuild_adjacency();
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
                    m.rebuild_adjacency();
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

/// ★ **A pin standing in a blind round hole, its top flush with the plate's top, is one body.**
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
                m.rebuild_adjacency();
                let floor = offset_of(&mut m, base, 0.5);
                let b = prism_on(&mut m, floor, top, 0.5);
                let (x, y) = if swapped { (b, cup) } else { (cup, b) };
                let got =
                    boolean(&mut m, BoolKind::Fuse, x, y).unwrap_or_else(|e| panic!("{at}: {e:?}"));
                m.rebuild_adjacency();
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

/// ★ **An *area* contact opens the wall instead of leaving a zero-thickness one.**
///
/// This check looks for an edge inside a face, which is a line contact; two faces of one solid
/// overlapping *in area* could in principle escape it, since the boundaries might coincide with no
/// edge strictly interior to anything. So the question was put to the kernel rather than reasoned
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
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
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
    m.rebuild_adjacency();
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
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, 1.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    m.rebuild_adjacency();
    let got = boolean(&mut m, BoolKind::Fuse, a, b).expect("the fuse builds");
    assert_eq!(got.len(), 1);
    let v = nacre_props::mass_props(&m, got[0]).expect("props").volume;
    assert!((v - (8.0 + 8.0 - 1.0)).abs() < 1e-9, "volume {v}");
    assert!(nacre_validate::validate(&m).is_empty());
}

/// ★ **A face with a hole, so the even-odd rule is exercised and not merely present.**
///
/// The judge folds *every* ring of a face into one bag and reads even-odd, because a point inside a
/// hole has crossed the boundary twice. Writing it as "inside the outer ring **and** outside each
/// hole" was one of the two shapes that got the earlier drafts wrong, so the rule needs a fixture
/// rather than a corpus that happens to cover it: the census's edges reach the exact test against
/// 214 faces, none of them holed (measured).
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
    m.rebuild_adjacency();
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
