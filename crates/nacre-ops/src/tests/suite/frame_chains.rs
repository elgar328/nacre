//! **Statements behind a frame node, in the world** — a cylinder drawn on a datum whose axes are
//! irrational as a frame but whose normal has a rational length, and what its chain's fold carries.
//!
//! The datum is the one an application's `plane({ origin, xPoint, yHint })` takes: normal
//! `(2, 3, 6)` (length `7`), `+u` along `(−3, 2, 0)` (length `√13`). Its map is irrational, and yet
//! a circle at its origin extrudes to a cylinder whose origin, axis and seam direction are rational
//! in the world, and whose far cap `n·p = n·o + 4·7` is too. Volumes are derived, not copied.

use super::*;
use nacre_exact::{Isometry, Rat};

const NORMAL_SEVEN: Axes = ([-3.0, 2.0, 0.0], [-12.0, -18.0, 13.0]);

/// A circle of radius `1` at `(cu, cv)` on the datum through `origin` with axes `x`, `y`,
/// extruded `4` along the datum's normal.
fn datum_cylinder(
    m: &mut Model,
    origin: [f64; 3],
    (x, y): Axes,
    centre: [f64; 2],
) -> Handle<Solid> {
    let frame = datum_frame(m, origin, (x, y));
    let ring =
        crate::Ring2d::circle(nacre_math::Point2::from_array(centre), 1.0).expect("a circle");
    extrude(
        m,
        frame,
        Profile2d::from_normalized_rings(ring, vec![]),
        4.0,
    )
}

fn datum_frame(m: &mut Model, origin: [f64; 3], (x, y): Axes) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(SketchPlane::from_axes(
                Point3::from_array(origin),
                Vector3::from_array(x),
                Vector3::from_array(y),
            )),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("the datum: {other:?}"),
    }
}

fn extrude(m: &mut Model, frame: SketchFrame, profile: Profile2d, dist: f64) -> Handle<Solid> {
    match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    ) {
        Ok(OpOutput::Extrude { solid, .. }) => solid,
        other => panic!("the extrude: {other:?}"),
    }
}

fn lateral(m: &Model, s: Handle<Solid>) -> Handle<Surface> {
    m.shell(m.solid(s).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|&h| matches!(m.surface(h), nacre_topo::Surface::Cylinder { .. }))
        .expect("the lateral")
}

fn volume(m: &Model, solids: &[Handle<Solid>]) -> f64 {
    solids
        .iter()
        .map(|&s| nacre_props::mass_props(m, s).expect("a volume").volume)
        .sum()
}

const PI4: f64 = 4.0 * std::f64::consts::PI;

/// A move applied to a solid after its frame.
type Move = Box<dyn Fn(&mut Model, Handle<Solid>) -> Handle<Solid>>;

/// A datum's axes, as [`datum_frame`] takes them.
type Axes = ([f64; 3], [f64; 3]);

/// The three answers against a box far away: the cylinder and the box, the cylinder, nothing.
fn against_a_far_box(what: &str, make: impl Fn(&mut Model) -> Handle<Solid>) {
    for (kind, want, pieces) in [
        (BoolKind::Fuse, PI4 + 8.0, 2),
        (BoolKind::Cut, PI4, 1),
        (BoolKind::Common, 0.0, 0),
    ] {
        let mut m = Model::new();
        let c = make(&mut m);
        let b = crate::fixtures::cuboid(&mut m, [50.0, 50.0, 50.0], [52.0, 52.0, 52.0]);
        let got = boolean(&mut m, kind, c, b).unwrap_or_else(|e| panic!("{what} {kind:?}: {e:?}"));
        assert_eq!(got.len(), pieces, "{what} {kind:?}");
        let v = volume(&m, &got);
        assert!(
            (v - want).abs() < 1e-9,
            "{what} {kind:?}: volume {v}, want {want}"
        );
        assert!(nacre_validate::validate(&m).is_empty(), "{what} {kind:?}");
    }
}

/// **A cylinder drawn at the origin of a datum whose normal has a rational length has a world
/// statement, and booleans like any other.** The datum's frame is irrational (`|u| = √13`); the
/// chain's fold carries the cylinder's origin, axis and seam direction and both caps to the
/// rationals all the same, at the world origin and away from it.
#[test]
fn a_cylinder_at_a_rational_normal_datums_origin_booleans_in_the_world() {
    for origin in [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]] {
        let mut m = Model::new();
        let c = datum_cylinder(&mut m, origin, NORMAL_SEVEN, [0.0, 0.0]);
        let world = m
            .world_cylinder_def(lateral(&m, c))
            .expect("a world statement");
        let r = |x: i128| Rat::from_int(x);
        let along = [r(2), r(3), r(6)];
        assert!(
            nacre_exact::parallel_rat(&world.dir(), &along),
            "{origin:?}: the axis runs along the datum's normal"
        );
        against_a_far_box(&format!("{origin:?}"), |m| {
            datum_cylinder(m, origin, NORMAL_SEVEN, [0.0, 0.0])
        });
    }
}

/// **The fold composes with what moves the body afterwards** — a translation and a quarter turn
/// after the frame each keep the world statement and the answers.
#[test]
fn a_frame_cylinder_moved_afterwards_keeps_its_world_statement() {
    let r = |n: i128, d: i128| Rat::new(n, d).expect("a rational");
    let moves: [(&str, Move); 2] = [
        (
            "translated",
            Box::new(move |m, s| {
                crate::transform::transform(
                    m,
                    s,
                    &Isometry::translation([r(3, 2), r(0, 1), r(-1, 1)]),
                )
                .expect("moves")
            }),
        ),
        (
            "a quarter turn",
            Box::new(move |m, s| {
                crate::transform::transform(
                    m,
                    s,
                    &Isometry::rotation(nacre_exact::Rotation {
                        axis: Axis::X,
                        pivot: [r(0, 1), r(1, 1), r(0, 1)],
                        angle: nacre_exact::Angle::from_deg(r(90, 1)).expect("an angle"),
                    }),
                )
                .expect("turns")
            }),
        ),
    ];
    for (what, moved) in &moves {
        let mut m = Model::new();
        let c = datum_cylinder(&mut m, [1.0, 2.0, 3.0], NORMAL_SEVEN, [0.0, 0.0]);
        let c = moved(&mut m, c);
        assert!(
            m.world_cylinder_def(lateral(&m, c)).is_ok(),
            "{what}: a world statement"
        );
        against_a_far_box(what, |m| {
            let c = datum_cylinder(m, [1.0, 2.0, 3.0], NORMAL_SEVEN, [0.0, 0.0]);
            moved(m, c)
        });
    }
}

/// **A mirror after the frame reverses the chain's handedness, and the planes' world facing
/// follows it.** A block on the datum (the mirror takes planar solids), mirrored in `z = −1/2`:
/// every face whose plane has a world name faces, by that name and its sense, the way its own
/// loop winds — read off the vertices' realized coordinates, which go through the replay and not
/// through the fold. And the block booleans like any other.
#[test]
fn a_mirrored_frame_block_faces_the_way_its_loops_wind() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let block = |m: &mut Model| -> Handle<Solid> {
        let frame = datum_frame(m, [1.0, 2.0, 3.0], NORMAL_SEVEN);
        let b = extrude(
            m,
            frame,
            Profile2d::polygon(vec![p(0.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(0.0, 2.0)])
                .expect("a rectangle"),
            4.0,
        );
        match apply(
            m,
            &Operation::Mirror {
                solid: b,
                axis: Axis::Z,
                offset: Rat::new(-1, 2).expect("a rational"),
            },
        ) {
            Ok(OpOutput::Mirror { solid }) => solid,
            other => panic!("the mirror: {other:?}"),
        }
    };
    let mut m = Model::new();
    let b = block(&mut m);
    let mut named = 0;
    for &f in &m.shell(m.solid(b).outer).faces {
        let face = m.face(f);
        let Ok(name) = m.world_plane_name(face.surface) else {
            continue;
        };
        if m.surface_motion(face.surface)
            .is_some_and(|leaf| m.chain_through_frame(leaf))
        {
            named += 1;
        }
        let sense = m
            .world_plane_name_sense(face.surface)
            .expect("a world name faces some way");
        let pts: Vec<Point3> = face
            .outer
            .half_edges
            .iter()
            .map(|&he| m.vertex_point(m.he_start(he)))
            .collect();
        let mut newell = [0.0; 3];
        for k in 0..pts.len() {
            let (a, c) = (pts[k].as_array(), pts[(k + 1) % pts.len()].as_array());
            newell[0] += (a[1] - c[1]) * (a[2] + c[2]);
            newell[1] += (a[2] - c[2]) * (a[0] + c[0]);
            newell[2] += (a[0] - c[0]) * (a[1] + c[1]);
        }
        let n = name.narrow().expect("a narrow name");
        let along: f64 = (0..3).map(|k| n[k].to_f64() * newell[k]).sum();
        let want = f64::from(sense.sign() * face.orientation.sign());
        assert!(
            along * want > 0.0,
            "a face's world name and sense against its loop: {along} want sign {want}"
        );
    }
    assert_eq!(
        named, 2,
        "the far cap and the wall through the origin have world names through the frame"
    );
    let c = crate::fixtures::cuboid(&mut m, [50.0, 50.0, 50.0], [52.0, 52.0, 52.0]);
    let got = boolean(&mut m, BoolKind::Fuse, b, c).expect("the fuse");
    let v = volume(&m, &got);
    assert!((v - 40.0).abs() < 1e-9, "volume {v}");
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **A block on the same datum halves the cylinder** — the block's wall through the datum's
/// origin is the plane `û·x = 0` there, rational in the world, and it cuts the cylinder through
/// its axis: the common part is half the cylinder, `2π`, and the cut keeps the other half.
#[test]
fn a_block_on_the_same_datum_halves_the_cylinder() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    for (kind, want) in [(BoolKind::Common, PI4 / 2.0), (BoolKind::Cut, PI4 / 2.0)] {
        let mut m = Model::new();
        let c = datum_cylinder(&mut m, [0.0, 0.0, 0.0], NORMAL_SEVEN, [0.0, 0.0]);
        let frame = datum_frame(&mut m, [0.0, 0.0, 0.0], NORMAL_SEVEN);
        let block = extrude(
            &mut m,
            frame,
            Profile2d::polygon(vec![p(0.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(0.0, 2.0)])
                .expect("a rectangle"),
            4.0,
        );
        match boolean(&mut m, kind, c, block) {
            Ok(got) => {
                let v = volume(&m, &got);
                assert!((v - want).abs() < 1e-9, "{kind:?}: volume {v}, want {want}");
                assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
            }
            Err(e) => panic!("{kind:?}: {e:?}"),
        }
    }
}

/// **What the fold cannot carry is refused by its cause.** A circle off the datum's origin puts the
/// cylinder's axis one along `û` — at `u/√13` — from it, irrational; a datum through `(1, 2, 3)`
/// whose normal is `(1, 1, 0)` keeps the cylinder rational but puts its far cap at
/// `x + y = 3 + 4√2`. Neither has a world statement where it matters, and the gate says so by name.
#[test]
fn an_irrational_statement_behind_a_frame_is_refused_by_its_cause() {
    let sqrt2 = ([-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let cases: [(&str, [f64; 3], Axes, [f64; 2]); 2] = [
        ("off the origin", [1.0, 2.0, 3.0], NORMAL_SEVEN, [1.0, 0.0]),
        ("a √2 normal", [1.0, 2.0, 3.0], sqrt2, [0.0, 0.0]),
    ];
    for (what, origin, axes, centre) in cases {
        let mut m = Model::new();
        let c = datum_cylinder(&mut m, origin, axes, centre);
        if what == "off the origin" {
            assert_eq!(
                m.world_cylinder_def(lateral(&m, c)),
                Err(nacre_topo::WorldFail::Unstated),
                "{what}: the axis is off the rationals"
            );
        } else {
            assert!(
                m.world_cylinder_def(lateral(&m, c)).is_ok(),
                "{what}: the cylinder itself is rational"
            );
            let unstated = m
                .shell(m.solid(c).outer)
                .faces
                .iter()
                .map(|&f| m.face(f).surface)
                .filter(|&h| matches!(m.surface(h), nacre_topo::Surface::Plane { .. }))
                .filter(|&h| m.world_plane_name(h) == Err(nacre_topo::WorldFail::Unstated))
                .count();
            assert_eq!(
                unstated, 1,
                "{what}: the far cap, and only it, has no world name"
            );
        }
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let mut m = Model::new();
            let c = datum_cylinder(&mut m, origin, axes, centre);
            let b = crate::fixtures::cuboid(&mut m, [50.0, 50.0, 50.0], [52.0, 52.0, 52.0]);
            assert_eq!(
                boolean(&mut m, kind, c, b),
                Err(BoolError::Rejected {
                    reason: RejectReason::CylinderGateUndecided,
                    at: None
                }),
                "{what} {kind:?}"
            );
        }
    }
}

/// **The cylinder drills a world plate** — the way an application uses it: a plate `[−5,5]²×[1,2]`
/// that the cylinder crosses whole (its tilted caps reach `z = 0 ± 0.52` and `24/7 ± 0.52`). The
/// plate's faces meet the lateral in ellipses; a horizontal section of a radius-`1` cylinder whose
/// axis leans to `n_z/|n| = 6/7` is an ellipse of area `π·7/6`, so the common part is `7π/6`, the
/// drilled plate `100 − 7π/6`, the union `100 + 4π − 7π/6`.
#[test]
fn a_frame_cylinder_drills_a_world_plate() {
    let through = 7.0 * std::f64::consts::PI / 6.0;
    for (kind, want) in [
        (BoolKind::Common, through),
        (BoolKind::Cut, 100.0 - through),
        (BoolKind::Fuse, 100.0 + PI4 - through),
    ] {
        let mut m = Model::new();
        let c = datum_cylinder(&mut m, [0.0, 0.0, 0.0], NORMAL_SEVEN, [0.0, 0.0]);
        let plate = crate::fixtures::cuboid(&mut m, [-5.0, -5.0, 1.0], [5.0, 5.0, 2.0]);
        let (a, b) = if kind == BoolKind::Cut {
            (plate, c)
        } else {
            (c, plate)
        };
        let got = boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
        let v = volume(&m, &got);
        assert!((v - want).abs() < 1e-9, "{kind:?}: volume {v}, want {want}");
        assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
    }
}

/// **A wide name behind a frame is asked its geometry before its width.** A pad on the datum whose
/// corners are long decimals of mixed scale gives one wall a `Wide` name; that wall misses the
/// datum's origin, so its world equation `n·p = n·t + δ·√13` is irrational — `Unstated`, not the
/// width of a name that has no rational image to be too wide for.
#[test]
fn a_wide_wall_behind_a_frame_is_unstated_for_its_geometry() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let frame = datum_frame(&mut m, [1.0, 2.0, 3.0], NORMAL_SEVEN);
    let pad = extrude(
        &mut m,
        frame,
        Profile2d::polygon(vec![
            p(1.234567890123457e-7, -9.876543210987654),
            p(19.87654321098765, -1.23456789012345e-7),
            p(13.141592653589793, 12.718281828459045),
            p(2.718281828459045e-7, 7.071067811865476),
        ])
        .expect("a quadrilateral"),
        4.0,
    );
    let wide: Vec<_> = m
        .shell(m.solid(pad).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .filter(|&h| m.surface_name(h).is_some_and(|n| n.narrow().is_none()))
        .collect();
    assert_eq!(wide.len(), 1, "premise: one wall's name is wide");
    assert!(
        m.surface_motion(wide[0])
            .is_some_and(|leaf| m.chain_through_frame(leaf)),
        "premise: the wall is behind the datum's frame"
    );
    assert_eq!(
        m.world_plane_name(wide[0]),
        Err(nacre_topo::WorldFail::Unstated)
    );
}

/// **A chain holds a frame whether or not it folds.** A block on the datum turned `37°` about `z`
/// has no fold (the turn is off the quarters); [`Model::chain_through_frame`] answers what walking
/// each face's chain finds — the walls and the far cap behind the datum's frame, the bottom cap
/// (the datum plane itself, stated in the world) not.
#[test]
fn a_turned_frame_block_is_through_a_frame_without_a_fold() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let frame = datum_frame(&mut m, [1.0, 2.0, 3.0], NORMAL_SEVEN);
    let block = extrude(
        &mut m,
        frame,
        Profile2d::polygon(vec![p(0.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(0.0, 2.0)])
            .expect("a rectangle"),
        4.0,
    );
    let turned = crate::transform::transform(
        &mut m,
        block,
        &Isometry::rotation(nacre_exact::Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: nacre_exact::Angle::from_deg(Rat::from_int(37)).expect("an angle"),
        }),
    )
    .expect("turns");
    let mut framed = 0;
    for &f in &m.shell(m.solid(turned).outer).faces {
        let Some(leaf) = m.surface_motion(m.face(f).surface) else {
            continue;
        };
        let mut walked = false;
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            walked |= matches!(m.motion(h).motion, nacre_topo::Motion::Frame { .. });
            cur = m.motion(h).parent;
        }
        assert_eq!(m.chain_through_frame(leaf), walked);
        if walked {
            framed += 1;
            assert_eq!(
                m.chain_dir_along(leaf, [Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]),
                Err(nacre_topo::WorldFail::Unstated),
                "premise: the chain has no fold (with one, ŵ's line would be rational)"
            );
        }
    }
    assert_eq!(framed, 5, "the four walls and the far cap");
}

/// **A cylinder's statement asks the geometry of all its parts before any width.** A circle at the
/// origin of the `√2` datum (`|u| = √2`, `|v| = 1`) whose seam starts at `(3/5, 4/5)` has a seam
/// reference mixing `û` and `v̂` in irrational ratio; two steps whose offsets cannot compose in
/// `Rat` (`1/10²⁰ + 1/(10²⁰ + 1)`) lose the origin. The seam is irrational whatever the origin's
/// width: `Unstated`, and the gate refuses by that name.
#[test]
fn a_cylinder_with_a_lost_origin_is_unstated_for_its_seam() {
    let r = |n: i128, d: i128| Rat::new(n, d).expect("a rational");
    let e20 = 10_i128.pow(20);
    let make = |m: &mut Model, steps: &[Rat]| {
        let frame = datum_frame(m, [0.0; 3], ([-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]));
        let ring = crate::Ring2d::new(
            vec![[r(3, 5), r(4, 5)]],
            vec![nacre_geom::mixed::Edge2d::Arc {
                center: [r(0, 1), r(0, 1)],
                r2: r(1, 1),
                ccw: true,
            }],
        )
        .expect("a circle");
        let mut c = extrude(
            m,
            frame,
            Profile2d::from_normalized_rings(ring, vec![]),
            4.0,
        );
        for &s in steps {
            c = crate::transform::transform(m, c, &Isometry::translation([s, r(0, 1), r(0, 1)]))
                .expect("moves");
        }
        c
    };
    for steps in [vec![r(1, e20)], vec![r(1, e20), r(1, e20 + 1)]] {
        let mut m = Model::new();
        let c = make(&mut m, &steps);
        assert_eq!(
            m.world_cylinder_def(lateral(&m, c)),
            Err(nacre_topo::WorldFail::Unstated),
            "{} steps",
            steps.len()
        );
        let b = crate::fixtures::cuboid(&mut m, [50.0, 50.0, 50.0], [52.0, 52.0, 52.0]);
        assert_eq!(
            boolean(&mut m, BoolKind::Fuse, c, b),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                at: None
            }),
            "{} steps",
            steps.len()
        );
    }
}
