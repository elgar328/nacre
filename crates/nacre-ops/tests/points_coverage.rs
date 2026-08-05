//! **The S6b precondition, asserted: every live planar face records its exact points.**
//!
//! S6a drained the point-less-plane producers one by one — `from_axes` lifts its axes,
//! `add_cylinder` records its caps, an overflowing exact move records a node — and this is the
//! gate that kept them drained: one model per producer path, and a sweep at the end that
//! refused any live planar face without recorded points. S6b then absorbed the points into the
//! surface's truth variant, so the type now guarantees what the sweep used to check — the
//! battery stays as the producer-path smoke test and the cache/truth agreement sweep.
//!
//! Cylinder *lateral* surfaces are the deliberate exception: a curved surface's truth arrives
//! with M6, and booleans already reject it honestly (`CylinderFace`).

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_topo::Model;

fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![p2(a, a), p2(b, a), p2(b, b), p2(a, b)]).unwrap()
}

fn extrude(m: &mut Model, plane: SketchPlane, profile: Profile2d, dist: f64) {
    let out = apply(
        m,
        &Operation::Extrude {
            plane,
            profile,
            dist,
        },
    )
    .expect("extrude");
    assert!(matches!(out, OpOutput::Extrude { .. }));
    m.rebuild_adjacency();
}

/// Every live planar face of `m` records its three exact points — asserted through the truth
/// store. Since S6b the type guarantees it (a `Plane` truth *is* three points), so the sweep is
/// a retrospective record of what this gate had to check while the point-less population still
/// existed — plus the cache/truth variant agreement, which is not type-carried.
fn assert_all_planes_record_points(m: &Model, what: &str) {
    let mut seen = 0usize;
    for &s in &m.live_solids {
        for &sh in std::iter::once(&m.solids.get(s).outer).chain(m.solids.get(s).cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                let surf = m.faces.get(fh).surface;
                match (m.surface(surf), m.surface_truth(surf)) {
                    (nacre_geom::Surface::Plane(_), nacre_topo::SurfaceTruth::Plane { .. }) => {
                        seen += 1;
                    }
                    (
                        nacre_geom::Surface::Cylinder(_),
                        nacre_topo::SurfaceTruth::Cylinder { .. },
                    ) => {
                        // curved truth is M6's; booleans reject it honestly today
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

    // ② The frame roads: a named tilted plane, and an axes-only tilted frame (S6a's opening).
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
        .shells
        .get(m.solids.get(base).outer)
        .faces
        .iter()
        .find(|&&fh| {
            let s = m.faces.get(fh).surface;
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
                    point: [Rat::from_int(0); 3],
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

    // ⑥ Cylinder caps (S6a) — the lateral face is the skipped population, visibly.
    let mut mc = Model::new();
    mc.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    mc.rebuild_adjacency();
    assert_all_planes_record_points(&mc, "cylinder caps");

    let _ = mirrored;
}

/// ★★ S9: a `z = 0` sketch's base cap **is** the seeded XY plane — the operation and the
/// pre-seeded vocabulary meet at one handle. And replay determinism holds with seeds included:
/// the same log twice gives the same handles.
#[test]
fn a_world_sketch_base_cap_is_the_seeded_plane() {
    let log = [Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: square(0.0, 2.0),
        dist: 1.0,
    }];
    let m1 = nacre_ops::replay(&log).unwrap();
    let m2 = nacre_ops::replay(&log).unwrap();
    let cap_surface = |m: &Model| {
        let s = m.live_solids[0];
        // faces[0] of the extrude output is the base cap; find via the z=0 name instead of
        // output plumbing.
        m.shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .map(|&fh| m.faces.get(fh).surface)
            .find(|su| su == &m.world_plane(nacre_scalar::Axis::Z))
            .expect("the z = 0 base cap must intern onto the seed")
    };
    let (c1, c2) = (cap_surface(&m1), cap_surface(&m2));
    assert_eq!(c1, m1.world_plane(nacre_scalar::Axis::Z));
    assert_eq!(c1, c2, "replay is deterministic down to seeded handles");
}
