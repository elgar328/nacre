//! **The S6b precondition, asserted: every live planar face records its exact points.**
//!
//! S6a drained the point-less-plane producers one by one — `from_axes` lifts its axes,
//! `add_cylinder` records its caps, an overflowing exact move records a node — and this is the
//! gate that keeps them drained: one model per producer path, and a sweep at the end that
//! refuses any live planar face without a `surface_points` entry. When `Surface::Plane` absorbs
//! its points as a variant (S6b), a failure here is the population that would have made that
//! change refuse to compile honestly.
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

/// S6b commit-1 mirror invariant: for every live face's surface, the truth store says exactly
/// what the (transitional) side tables say — same existence, same motion, same points. This is
/// the bridge that lets consumers move from the tables to the truth one commit at a time; it
/// dies with the tables.
fn assert_truth_mirrors_the_tables(m: &Model, what: &str) {
    use nacre_topo::{PlanePoints, SurfaceDef, SurfaceTruth};
    for &s in &m.live_solids {
        for &sh in std::iter::once(&m.solids.get(s).outer).chain(m.solids.get(s).cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                let surf = m.faces.get(fh).surface;
                let def = m.surface_defs.get(&surf).copied();
                let pts = m.surface_points.get(&surf).copied();
                match (m.surface_truth(surf), m.surface(surf)) {
                    (Some(SurfaceTruth::Cylinder { motion }), _) => {
                        assert!(
                            matches!(m.surface(surf), nacre_geom::Surface::Cylinder(_)),
                            "{what}: a Cylinder truth on a non-cylinder cache"
                        );
                        let want = match def {
                            Some(SurfaceDef::Moved { motion }) => Some(motion),
                            _ => None,
                        };
                        assert_eq!(*motion, want, "{what}: cylinder motion mirror");
                    }
                    (Some(SurfaceTruth::Plane { points, motion }), _) => {
                        let PlanePoints::Known(p) = points;
                        assert_eq!(Some(*p), pts, "{what}: plane points mirror");
                        let want = match def {
                            Some(SurfaceDef::Moved { motion }) => Some(motion),
                            _ => None,
                        };
                        assert_eq!(*motion, want, "{what}: plane motion mirror");
                        assert!(
                            !matches!(def, Some(SurfaceDef::Inexact)),
                            "{what}: an Inexact def with a truth entry"
                        );
                    }
                    (None, nacre_geom::Surface::Plane(_)) => {
                        assert!(
                            pts.is_none() || matches!(def, Some(SurfaceDef::Inexact)),
                            "{what}: a point-bearing, non-Inexact plane without truth"
                        );
                    }
                    (None, nacre_geom::Surface::Cylinder(_)) => {
                        panic!("{what}: a cylinder push always records truth");
                    }
                }
            }
        }
    }
}

/// Every live planar face of `m` records its three exact points.
fn assert_all_planes_record_points(m: &Model, what: &str) {
    let (mut with, mut without) = (0usize, 0usize);
    for &s in &m.live_solids {
        for &sh in std::iter::once(&m.solids.get(s).outer).chain(m.solids.get(s).cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                let surf = m.faces.get(fh).surface;
                if !matches!(m.surface(surf), nacre_geom::Surface::Plane(_)) {
                    continue; // curved truth is M6's; booleans reject it honestly today
                }
                if m.surface_points.contains_key(&surf) {
                    with += 1;
                } else {
                    without += 1;
                }
            }
        }
    }
    assert!(with > 0, "{what}: the sweep saw no planar faces at all");
    assert_eq!(
        without,
        0,
        "{what}: {without} of {} live planar faces record no exact points",
        with + without
    );
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
    assert_truth_mirrors_the_tables(&m, "construction");

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
    assert_truth_mirrors_the_tables(&m, "frame roads");

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
    assert_truth_mirrors_the_tables(&m, "pad");

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
    assert_truth_mirrors_the_tables(&m, "motions");

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
    assert_truth_mirrors_the_tables(&mb, "boolean");

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
    assert_truth_mirrors_the_tables(&mc, "cylinder caps");

    let _ = mirrored;
}
