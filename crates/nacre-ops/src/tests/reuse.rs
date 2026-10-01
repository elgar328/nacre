use super::*;

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: crate::SketchPlane) -> crate::SketchFrame {
    match crate::apply(
        m,
        &crate::Operation::DatumPlane {
            def: crate::DatumDef::Stated(plane),
        },
    ) {
        Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

use crate::{OpOutput, Operation, apply};
use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::{Point2, Point3, Vector3};

fn square(a: f64, b: f64) -> crate::Profile2d {
    crate::Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
    .unwrap()
}

/// ★ **Witnesses are solved from the definition, not lifted from the cache.** A box stated at
/// `0.3` has corners `f64` cannot hold; the rounded cache lifted back to a rational is a
/// different point, claimed with tol 0. The definition road solves
/// the corner from its three plane names and states the rounding it carries: nonzero here,
/// exactly zero for an integer box.
#[test]
fn witnesses_are_solved_from_the_definition_and_carry_their_rounding() {
    let mut m = Model::new();
    let unit = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let dec = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.3, -1.0, 0.3]),
        Point3::from_array([0.7, 2.0, 0.7]),
    );
    m.rebuild_adjacency();

    let unit_pts = solid_points(&m, unit).expect("an integer box answers");
    assert_eq!(unit_pts.len(), 8);
    assert!(
        unit_pts.iter().all(|p| p.tol() == [0.0; 3]),
        "an integer corner is an f64: tol 0"
    );

    let dec_pts = solid_points(&m, dec).expect("a decimal box answers");
    assert_eq!(dec_pts.len(), 8);
    for p in &dec_pts {
        // The coordinate is the definition's nearest f64 — for a constructed box the same
        // value the cache holds, so what the road changes is the *tolerance*, not the point.
        let c = p.coord();
        assert!(c.iter().all(|x| [0.3, 0.7, -1.0, 2.0].contains(x)), "{c:?}");
        assert!(
            p.tol().iter().all(|&t| t > 0.0),
            "0.3/0.7 are not f64 — the witness must say so on every axis: {:?}",
            p.tol()
        );
    }
}

/// ★ **A boolean's result answers too**: its vertices' definition has a rational base.
#[test]
fn a_boolean_result_answers_from_its_definition() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0; 3]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0; 3]),
        Point3::from_array([3.0; 3]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean::boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 1);
    let pts =
        solid_points(&m, out[0]).expect("a result's vertices are three-plane meets with names");
    assert!(pts.len() >= 8, "{}", pts.len());
    assert!(pts.iter().all(|p| p.tol() == [0.0; 3]), "integer corners");
}

/// The def road's answer on one solid: `Some` with every point realized finite, or `None`.
fn assert_answers(m: &Model, s: Handle<Solid>, want_some: bool, what: &str) {
    match solid_points(m, s) {
        Some(pts) => {
            assert!(
                want_some,
                "{what}: expected a decline, got {} points",
                pts.len()
            );
            assert!(!pts.is_empty(), "{what}: answered with no points");
            for p in &pts {
                assert!(
                    p.coord().iter().all(|c| c.is_finite()),
                    "{what}: a realized coordinate is not finite"
                );
            }
        }
        None => assert!(!want_some, "{what}: expected an answer, got a decline"),
    }
}

/// ★★★ **The lock on the population that reaches this module.** A decimal-framed
/// prism's constructed corners solve from their carriers' in-frame names, whose rational
/// Cramer overflows on every one (measured 8/8 — while all eight points fit `Rat`), and
/// [`solid_points`] gives the whole solid up on the first failure — so a framed operand
/// would cost a boolean its entire class reuse.
#[test]
fn a_framed_prisms_corners_solve_for_reuse() {
    let mut m = Model::new();
    let plane = crate::SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let frame = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid: prism, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: crate::Profile2d::polygon(vec![
                Point2::from_array([0.1111111111111111, 0.1234567890123456]),
                Point2::from_array([4.123456789012345, 0.2345678901234567]),
                Point2::from_array([3.9876543210987654, 3.1234567890123459]),
                Point2::from_array([0.2222222222222222, 2.765432109876543]),
            ])
            .unwrap(),
            dist: 2.5,
        },
    )
    .expect("the framed prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    // ① The solve: all eight corners — the narrow Cramer alone answered none of them.
    let pts = solid_points(&m, prism).expect("a framed prism's corners solve");
    assert_eq!(pts.len(), 8, "eight corners, none given up");

    // ② One point, one realization road: every replayed coordinate IS a stored vertex
    // coordinate, bit for bit. A ulp here would mean construction and replay realize one
    // rational through two roads — a finding, not a tolerance.
    let mut stored: std::collections::HashSet<[u64; 3]> = std::collections::HashSet::new();
    let sol = m.solid(prism);
    for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
        for &fh in &m.shell(sh).faces {
            for &he in &m.face(fh).outer.half_edges {
                stored.insert(
                    m.vertex_point(he_start(&m, he))
                        .as_array()
                        .map(f64::to_bits),
                );
            }
        }
    }
    assert_eq!(stored.len(), 8, "a prism has eight distinct corners");
    for p in &pts {
        assert!(
            stored.contains(&p.coord().map(f64::to_bits)),
            "a replayed corner {:?} is not any stored coordinate",
            p.coord()
        );
    }

    // ③ The gate bites in the consultation direction that matters: a class owned by the world
    // cuboid asks for the *prism's* points, so it can leave `Arrange`. (The other direction —
    // prism-owned classes consulting the cuboid — proves nothing here.)
    let cub = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-2.0, -2.0, -2.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    m.rebuild_adjacency();
    let setup = crate::arrangement::plane_index_setup(&m, prism, cub).expect("plane setup");
    let plans = class_plans(
        &m,
        ClassReuse::Proved,
        BoolKind::Fuse,
        prism,
        cub,
        &setup.geom,
        &setup.class_owner,
    );
    let opened = plans
        .iter()
        .zip(&setup.class_owner)
        .filter(|(p, o)| **o == Some(SolidSide::B) && **p != ClassPlan::Arrange)
        .count();
    assert!(
        opened > 0,
        "no cuboid-owned class left Arrange — the reopened road was not exercised"
    );
}

/// ★★ **Which populations the def road answers for**, pinned per producer. Constructed and
/// turned solids answer; so does a boolean's result — the road does not read the cache, and
/// a result's corner is a three-plane meet with names like any other, so ③
/// answers. The one recorded decline is the mixed-frame population (④): there is no
/// sketch-frame base vertex to answer it through, and the def road
/// declines honestly — reuse falls back to Arrange, which is slower and never wrong.
#[test]
fn the_def_road_answers_for_the_populations_it_can_name() {
    let mut m = Model::new();
    // ① A constructed box, decimal-friendly and decimal-unfriendly corners.
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 2.0, 3.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.1, 5.0, 0.3]),
        Point3::from_array([1.7, 6.9, 2.2]),
    );
    m.rebuild_adjacency();
    assert_answers(&m, a, true, "constructed (integer corners)");
    assert_answers(&m, b, true, "constructed (decimal corners)");

    // ② Turned: an inexact rotation records a chain.
    let turned = {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: b,
                isometry: Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
                }),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        solid
    };
    m.rebuild_adjacency();
    assert_answers(&m, turned, true, "moved (rotated cuboid)");

    // ③ A fuse's result: its definition has a rational base, so it answers (reading the
    //    cache instead, a measured point "has no rational base").
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 2.5, 3.5]),
    );
    m.rebuild_adjacency();
    let fused = {
        let OpOutput::Boolean { solids } = apply(
            &mut m,
            &Operation::Boolean {
                kind: crate::BoolKind::Fuse,
                a,
                b: c,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        solids[0]
    };
    m.rebuild_adjacency();
    assert_answers(&m, fused, true, "discovered (a fused result)");

    // ④ The recorded decline, pinned: a tilted-frame prism's base ring sits under a
    //    world-stated cap (mixed frames), and no rational pullback exists — so the def road
    //    declines.
    let tilted = crate::SketchPlane::from_origin_normal(
        Point3::from_array([0.25, -0.5, 1.5]),
        Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
    )
    .unwrap();
    let prism = {
        let __frame0 = datum_frame(&mut m, tilted);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __frame0,
                profile: square(0.5, 2.5),
                dist: 1.1,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        solid
    };
    m.rebuild_adjacency();
    assert_answers(
        &m,
        prism,
        false,
        "mixed frames (a tilted prism's base ring)",
    );
}
