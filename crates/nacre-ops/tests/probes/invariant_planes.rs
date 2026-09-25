//! **A motion that fixes a plane restates nothing — the plane keeps its world statement.**
//!
//! The invariant-plane restatement: a rigid motion whose rotation axis is parallel
//! to a plane's normal, and whose translation slides within it, maps the plane onto itself as
//! a set. Such a plane's image is the source statement verbatim, so the transform pushes the
//! same statement and interns back onto the **source handle** — the road `Copy` already takes
//! for the identity motion — instead of minting a moved twin that carries a motion node.
//!
//! What that buys, each locked below: identity heals (the turned cap *is* the seed plane,
//! one handle), the cap never gains a history (a second turn re-qualifies), and its sketch frame is
//! the world one (`face_sketch_frame`; the contract sweep checks it bit for bit).

use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::Point2;
use nacre_ops::{OpOutput, Operation, Profile2d, SketchFrame, apply, face_sketch_frame};
use nacre_store::Handle;
use nacre_topo::{Model, PlanePoints, Solid, Surface};

fn block(m: &mut Model) -> Handle<Solid> {
    let world = SketchFrame::world(m, Axis::Z);
    let profile = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([2.0, 0.0]),
        Point2::from_array([2.0, 2.0]),
        Point2::from_array([0.0, 2.0]),
    ])
    .expect("profile");
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: world,
            profile,
            dist: 1.0,
        },
    )
    .expect("the block") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

fn rot_z30(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).expect("angle"),
            }),
        },
    )
    .expect("the turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// The z-cap surfaces of a solid, told apart from the walls by their truth's points (all
/// sharing one z) — not by handle order, which this test is about.
fn caps_and_walls(m: &Model, s: Handle<Solid>) -> (Vec<Handle<Surface>>, Vec<Handle<Surface>>) {
    let (mut caps, mut walls) = (vec![], vec![]);
    for &fh in &m.shell(m.solid(s).outer).faces {
        let surf = m.face(fh).surface;
        let Surface::Plane {
            points: PlanePoints::Known(p),
            ..
        } = m.surface(surf)
        else {
            panic!("a block face states known points");
        };
        if p[0][2] == p[1][2] && p[1][2] == p[2][2] {
            caps.push(surf);
        } else {
            walls.push(surf);
        }
    }
    (caps, walls)
}

#[test]
fn a_turned_cap_is_the_seed_plane_and_never_gains_a_history() {
    let mut m = Model::new();
    let s = block(&mut m);
    let (src_caps, _) = caps_and_walls(&m, s);

    let s = rot_z30(&mut m, s);
    let (caps, walls) = caps_and_walls(&m, s);
    assert_eq!(caps.len(), 2, "a block has two z-caps");

    // Identity heals: the turned bottom cap *is* the world XY seed — one plane, one handle —
    // and the top cap interned back onto its own extrude-minted source.
    let seed = m.world_plane(Axis::Z);
    assert!(
        caps.contains(&seed),
        "the turned bottom cap must intern back onto the seed: caps {caps:?}, seed {seed:?}"
    );
    for c in &caps {
        assert!(
            src_caps.contains(c),
            "a turned cap must be its source handle, not a moved twin: {c:?} vs {src_caps:?}"
        );
        assert!(
            matches!(m.surface(*c), Surface::Plane { motion: None, .. }),
            "a fixed plane records no motion"
        );
    }
    // The walls genuinely moved: each records the chain. This is what keeps the lock honest —
    // remove the restatement and the caps join them (new handles, `Some` motions).
    for w in &walls {
        assert!(
            matches!(
                m.surface(*w),
                Surface::Plane {
                    motion: Some(_),
                    ..
                }
            ),
            "a turned wall records its motion"
        );
    }

    // Both cap faces answer their world frame (the contract sweep verifies the returned frame
    // realizes bit-identical to the pad's; here the lock is that the answer exists at all).
    for &fh in &m.shell(m.solid(s).outer).faces {
        if caps.contains(&m.face(fh).surface) {
            face_sketch_frame(&m, fh).expect("a fixed cap hosts a sketch frame");
        }
    }

    // A second turn re-qualifies: the cap stayed world-stated, so the restatement's "source
    // carries no motion" condition holds again and the caps remain history-free at depth 2.
    let s = rot_z30(&mut m, s);
    let (caps2, _) = caps_and_walls(&m, s);
    assert!(caps2.contains(&seed), "still the seed after a second turn");
    for c in &caps2 {
        assert!(
            matches!(m.surface(*c), Surface::Plane { motion: None, .. }),
            "a twice-turned cap still records no motion"
        );
    }

    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the model stays valid"
    );
}
