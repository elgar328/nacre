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

use nacre_exact::{Angle, Axis, Isometry, Rat};
use nacre_math::{Point2, Point3};
use nacre_ops::{
    BoolError, BoolKind, Operation, Profile2d, RejectClass, RejectReason, SketchFrame, apply,
    boolean,
};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn prism(m: &mut Model, pts: &[[f64; 2]], h: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
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

fn moved(m: &mut Model, s: Handle<Solid>, isometry: Isometry) -> Handle<Solid> {
    let nacre_ops::OpOutput::Transform { solid } =
        apply(m, &Operation::Transform { solid: s, isometry }).expect("transform")
    else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// A unit cube, and a wedge sunk into it whose sharp tip reaches `tip_x`. At `tip_x == 1.0` the tip
/// lands exactly on the cube's far wall and the cut's cavity touches the outer shell along a line;
/// short of that it is an ordinary enclosed cavity.
fn cube_and_wedge(tip_x: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let b = prism(&mut m, &[[tip_x, 0.5], [0.1, 0.1], [0.1, 0.9]], 0.6);
    let b = moved(
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
    let a = m.add_cuboid(
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
    let a = m.add_cuboid(
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
    let a = moved(&mut m, a, turn());
    let b = moved(&mut m, b, turn());
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
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
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
    let b = moved(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = m.add_cuboid(
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
/// rather than a corpus that happens to cover it: 52 of the corpus's candidate faces are holed
/// (measured), but nothing there would fail if the hole handling broke.
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
