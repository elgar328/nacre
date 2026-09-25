//! **Every live planar face records its exact points.**
//!
//! No producer makes a point-less plane — `from_axes` lifts its axes, a circle prism states
//! its caps, an overflowing exact move records a node. The points live in the surface's truth
//! variant, so the type guarantees it; this battery is the producer-path smoke test — one model
//! per producer path — and the cache/truth agreement sweep.
//!
//! Cylinder *lateral* surfaces are not planes: their truth is a `CylinderDef`, so the sweep only
//! checks that their cache and truth agree on the kind.

use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::{Point3, Vector3};
use nacre_ops::SketchFrame;
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};
use nacre_topo::Model;

use crate::fixtures::{datum_frame, p2};

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![p2(a, a), p2(b, a), p2(b, b), p2(a, b)]).unwrap()
}

fn extrude(m: &mut Model, plane: SketchPlane, profile: Profile2d, dist: f64) {
    // The plane is stated first — a sketch names a plane the model holds.
    let frame = datum_frame(m, plane);
    let out = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("extrude");
    assert!(matches!(out, OpOutput::Extrude { .. }));
    m.rebuild_adjacency();
}

/// Every live planar face of `m` records its three exact points — asserted through the truth
/// store. The type guarantees the points (a `Plane` truth *is* three points); what the sweep
/// adds is the cache/truth variant agreement, which is not type-carried.
fn assert_all_planes_record_points(m: &Model, what: &str) {
    let mut seen = 0usize;
    for &s in m.live_solids() {
        for &sh in std::iter::once(&m.solid(s).outer).chain(m.solid(s).cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let surf = m.face(fh).surface;
                match (m.surface_cache(surf), m.surface(surf)) {
                    (nacre_geom::Surface::Plane(_), nacre_topo::Surface::Plane { .. }) => {
                        seen += 1;
                    }
                    (nacre_geom::Surface::Cylinder(_), nacre_topo::Surface::Cylinder { .. }) => {
                        // a lateral: its truth is a `CylinderDef`, not three points
                    }
                    (cache, truth) => panic!(
                        "{what}: cache and truth disagree about what this surface is: \
                         {cache:?} vs {truth:?}"
                    ),
                }
            }
        }
    }
    assert!(seen > 0, "{what}: the sweep saw no planar faces at all");
}

/// One model per producer path; the sweep is the assertion.
#[test]
fn every_live_planar_face_records_its_points() {
    // ① Plain construction: cuboid + world extrude.
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    m.rebuild_adjacency();
    extrude(&mut m, SketchPlane::world_xy(), square(10.0, 12.0), 1.5);
    assert_all_planes_record_points(&m, "construction");

    // ② The frame roads: a named tilted plane, and an axes-only tilted frame.
    let named = SketchPlane::from_origin_normal(
        Point3::from_array([0.25, -0.5, 1.5]),
        Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
    )
    .unwrap();
    extrude(&mut m, named, square(0.5, 2.5), 1.1);
    let axes_only = SketchPlane::from_axes(
        Point3::from_array([20.0, 0.0, 0.0]),
        Vector3::from_array([0.7123456789012345, 0.5876543210987654, 0.4098765432101234]),
        Vector3::from_array([-0.5876543210987654, 0.7123456789012345, 0.1234567890123456]),
    );
    // (The def's presence and wideness are pinned by the terminal lock in `nacre-ops`'s own
    // tests; here the sweep below is the assertion.)
    extrude(&mut m, axes_only, square(0.1, 1.6), 0.7);
    assert_all_planes_record_points(&m, "frame roads");

    // ③ Face features: pad and pocket share `extrude_and_boolean`.
    let top = *m
        .shell(m.solid(base).outer)
        .faces
        .iter()
        .find(|&&fh| {
            let s = m.face(fh).surface;
            m.surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, -1.0])
        })
        .expect("the cuboid top");
    let out = apply(
        &mut m,
        &Operation::PadOnFace {
            face: top,
            profile: square(1.0, 2.0),
            dist: 0.5,
        },
    )
    .expect("pad");
    let OpOutput::PadOnFace { solid: padded, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert_all_planes_record_points(&m, "pad");

    // ④ Motions: an inexact rotation (recorded node), a translation, a mirror.
    let turned = {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: padded,
                isometry: Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(30)).expect("angle"),
                }),
            },
        )
        .expect("rotate") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid
    };
    let OpOutput::Mirror { solid: mirrored } = apply(
        &mut m,
        &Operation::Mirror {
            solid: turned,
            axis: Axis::X,
            offset: Rat::new(-11, 2).unwrap(),
        },
    )
    .expect("mirror") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert_all_planes_record_points(&m, "motions");

    // ⑤ A boolean mints no surface, so its result inherits the record — but the sweep is what
    // says so, not the argument.
    let mut mb = Model::new();
    let a = mb.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = mb.add_cuboid(
        Point3::from_array([1.0, 1.0, 1.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    mb.rebuild_adjacency();
    nacre_ops::boolean(&mut mb, BoolKind::Fuse, a, b).expect("fuse");
    mb.rebuild_adjacency();
    assert_all_planes_record_points(&mb, "boolean");

    // ⑥ Cylinder caps — the lateral face is the skipped population, visibly.
    let mut mc = Model::new();
    nacre_ops::fixtures::cylinder_with_seam(
        &mut mc,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        1.0,
        2.0,
    );
    mc.rebuild_adjacency();
    assert_all_planes_record_points(&mc, "cylinder caps");

    let _ = mirrored;
}

/// **The world sugar names the seeded plane, and a sketch on it sits on that very handle.**
///
/// The frame **names** the seed up front, so "the base cap is the seed" is true by
/// construction and asserting it alone would be a tautology.
///
/// What still has teeth is the sugar: `SketchFrame::world` must name the seed the model planted,
/// not some other plane on the same geometry — and a sketch placed through it must land on that
/// handle. If the seeding or the arbitrary-axis convention moved, this is what notices.
/// ★★ A `z = 0` sketch's base cap **is** the seeded XY plane — the operation and the
/// pre-seeded vocabulary meet at one handle. And replay determinism holds with seeds included:
/// the same log twice gives the same handles.
#[test]
fn the_world_sugar_names_the_seed_a_sketch_then_sits_on() {
    let log = [Operation::Extrude {
        frame: SketchFrame::world(&Model::new(), Axis::Z),
        profile: square(0.0, 2.0),
        dist: 1.0,
    }];
    let m1 = nacre_ops::replay(&log).unwrap();
    let m2 = nacre_ops::replay(&log).unwrap();
    let cap_surface = |m: &Model| {
        let s = m.live_solids()[0];
        // faces[0] of the extrude output is the base cap; find via the z=0 name instead of
        // output plumbing.
        m.shell(m.solid(s).outer)
            .faces
            .iter()
            .map(|&fh| m.face(fh).surface)
            .find(|su| su == &m.world_plane(nacre_exact::Axis::Z))
            .expect("the z = 0 base cap must intern onto the seed")
    };
    let (c1, c2) = (cap_surface(&m1), cap_surface(&m2));
    assert_eq!(c1, m1.world_plane(nacre_exact::Axis::Z));
    assert_eq!(c1, c2, "replay is deterministic down to seeded handles");
}
