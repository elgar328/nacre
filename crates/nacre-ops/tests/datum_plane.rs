//! **A datum is a plane the model holds because someone named it.**
//!
//! Before this, a log could name exactly two kinds of plane: the three world planes `Model::new`
//! seeds (handles 0/1/2) and the plane of a face it had already built (`face_sketch_frame`).
//! Every other plane lived *inside* `Operation::Extrude` as a value, which is why that variant is
//! the last one still carrying a `SketchPlane` — and why the vocabulary swap that replaces it
//! needs this operation first.
//!
//! What the tests here pin is not that the call works but that the **arena answers the same way it
//! always did**: the same plane keeps one handle, the stored truth is the caller's own points, and
//! the whole thing replays.

use nacre_geom::Surface;
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{DatumDef, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, replay};
use nacre_scalar::Axis;
use nacre_topo::{FramePlacement, Model, PlanePoints, SurfaceTruth};

fn datum(plane: SketchPlane) -> Operation {
    Operation::DatumPlane {
        def: DatumDef::Stated(plane),
    }
}

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
    .expect("a square is a fair profile")
}

/// ★ **Stating a world plane returns the seed, and the arena does not grow.**
///
/// This is the lock on the cache-normal convention. `Model::new` seeds the world planes facing
/// `−axis` (S9 measured that `+axis` flipped 781 stored cap normals for nothing), `extrude` pushes
/// its base cap as `−plane.normal()`, and `PlaneGeom::frame_sign` reads a stored normal against a
/// root face's outward — which for a base cap is `−N`. A datum that pushed `+normal` would still
/// intern, but it would come back `flipped`, and every face built on it afterwards would have to
/// be spelled the other way round to compensate. `flipped == false` here says the datum joined the
/// convention rather than starting a second one.
#[test]
fn stating_a_world_plane_returns_the_seed_it_already_is() {
    for (axis, plane) in [
        (Axis::Z, SketchPlane::world_xy()),
        (Axis::X, SketchPlane::world_yz()),
        (Axis::Y, SketchPlane::world_zx()),
    ] {
        let mut m = Model::new();
        let before = m.surface_count();
        let seed = m.world_plane(axis);
        let OpOutput::DatumPlane { plane: h, .. } = apply(&mut m, &datum(plane)).expect("stated")
        else {
            unreachable!()
        };
        assert_eq!(h, seed, "{axis:?}: a world plane is the seeded one");
        assert_eq!(
            m.surface_count(),
            before,
            "{axis:?}: interning answered, so the arena must not grow"
        );
        // The convention check, stated where it can fail: the seed's cache and the datum's face
        // the same way, so nothing downstream has to compensate.
        let Surface::Plane(p) = m.surface(h) else {
            unreachable!()
        };
        let axis_v = match axis {
            Axis::X => [1.0, 0.0, 0.0],
            Axis::Y => [0.0, 1.0, 0.0],
            Axis::Z => [0.0, 0.0, 1.0],
        };
        assert!(
            p.normal().dot(Vector3::from_array(axis_v)) < 0.0,
            "{axis:?}: the stored cache faces −axis, as the seeds and base caps do"
        );
    }
}

/// The stored **truth** is the caller's own points, in the world, with no motion — asserted as
/// truth rather than inferred from coordinates.
#[test]
fn a_tilted_datum_records_the_points_the_caller_wrote() {
    let mut m = Model::new();
    let origin = Point3::from_array([1.0, 2.0, 3.0]);
    let sp = SketchPlane::from_origin_normal(origin, Vector3::from_array([1.0, 1.0, 1.0]))
        .expect("a nonzero normal names a plane");
    let OpOutput::DatumPlane { plane, frame } = apply(&mut m, &datum(sp)).expect("stated") else {
        unreachable!()
    };

    let SurfaceTruth::Plane { points, motion } = m.surface_truth(plane) else {
        unreachable!("a datum is a plane")
    };
    assert_eq!(*motion, None, "a world statement records no motion");
    let PlanePoints::Known(pts) = points;
    // The first point is the sketch origin, exactly — the caller's decimals, lifted once.
    assert_eq!(
        pts[0].map(|r| r.to_f64()),
        [1.0, 2.0, 3.0],
        "points[0] is the stated origin"
    );
    // Every recorded point satisfies the plane it names: `n·(p − o) == 0` in exact rationals.
    let zero = nacre_scalar::Rat::from_int(0);
    let n = [1i128, 1, 1].map(nacre_scalar::Rat::from_int);
    for p in pts.iter().skip(1) {
        let dot = (0..3).fold(zero, |acc, k| {
            let d = p[k].checked_sub(pts[0][k]).expect("small integers");
            acc.checked_add(d.checked_mul(n[k]).expect("small integers"))
                .expect("small integers")
        });
        assert_eq!(dot, zero, "a point off its plane");
    }

    // ★ The frame carries the caller's statement, not a derivation — the ZX convention is the
    // case that forbids deriving, so `Named` is unconditional.
    assert_eq!(frame.plane(), plane);
    assert!(!frame.flip(), "flip is the consumer's to measure");
    let FramePlacement::Named { origin: o, .. } = frame.placement() else {
        unreachable!("a stated datum names its placement")
    };
    assert_eq!(o.map(|r| r.to_f64()), [1.0, 2.0, 3.0]);
}

/// A plane with no exact statement is a **named reject**, not a panic and not a silent f64 plane.
#[test]
fn a_plane_without_an_exact_form_is_rejected_by_name() {
    let mut m = Model::new();
    // An origin past the decimal window (`~1e38`) leaves `def` empty while the f64 axes are still
    // perfectly ordinary — the exact shape `PlaneWithoutExactForm` is for.
    let sp = SketchPlane::from_axes(
        Point3::from_array([1e39, 0.0, 0.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
    );
    assert!(matches!(
        apply(&mut m, &datum(sp)),
        Err(OpError::PlaneWithoutExactForm)
    ));
    assert_eq!(
        m.surface_count(),
        3,
        "a rejected datum leaves the seeds alone"
    );
}

/// ★★ **The seam this whole stage exists to open**: a plane stated up front is the plane the
/// extrude then builds on — one handle, and the arena no larger for having said it twice.
///
/// This is the direct evidence that the vocabulary swap ahead (`Extrude` naming a `SketchFrame`
/// instead of carrying a `SketchPlane`) is not a change of population. `tests/plane_anchor.rs`
/// measures the other half — that the *cache* the datum leaves behind does not move the result.
#[test]
fn a_datum_is_the_plane_the_extrude_then_builds_on() {
    let sp = SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5]));

    let mut plain = Model::new();
    let OpOutput::Extrude { faces, .. } = apply(
        &mut plain,
        &Operation::Extrude {
            plane: sp,
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("plain extrude") else {
        unreachable!()
    };
    let base_cap = plain.faces.get(faces[0]).surface;

    let mut stated = Model::new();
    let OpOutput::DatumPlane { plane: h, .. } = apply(&mut stated, &datum(sp)).expect("stated")
    else {
        unreachable!()
    };
    let n_after_datum = stated.surface_count();
    let OpOutput::Extrude { faces, .. } = apply(
        &mut stated,
        &Operation::Extrude {
            plane: sp,
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("extrude on the stated plane") else {
        unreachable!()
    };
    assert_eq!(
        stated.faces.get(faces[0]).surface,
        h,
        "the base cap is the datum, not a second plane on the same geometry"
    );
    assert_eq!(
        base_cap.index(),
        h.index(),
        "and it lands at the index the plain extrude gave it"
    );
    assert_eq!(
        stated.surface_count() - n_after_datum,
        plain.surface_count() - 3 - 1,
        "the extrude minted the same surfaces either way, minus the one already stated"
    );
}

/// A log carrying a datum replays, and an orphan datum is not a validation failure.
///
/// A datum nothing sketches on is reachable from no face — exactly like the three seeds, which
/// have been permanent orphans since `Model::new` began planting them. `validate` counts the live
/// reachable set, and neither `nacre-step` nor `nacre-tess` walks the surface store at all, so an
/// unused datum cannot leak into output either.
#[test]
fn a_log_with_a_datum_replays_and_validates() {
    let sp = SketchPlane::from_origin_normal(
        Point3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .expect("a plane");
    let log = vec![
        datum(sp),
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    ];

    let mut scratch = Model::new();
    for op in &log {
        apply(&mut scratch, op).expect("each step is legal");
    }
    scratch.rebuild_adjacency();

    let replayed = replay(&log).expect("a datum log replays");
    assert_eq!(replayed.surface_count(), scratch.surface_count());
    assert_eq!(replayed.vertices.len(), scratch.vertices.len());
    assert_eq!(replayed.faces.len(), scratch.faces.len());
    assert_eq!(
        replayed.live_solids.to_vec(),
        scratch.live_solids.to_vec(),
        "replay reproduces the live set, indices included"
    );
    assert!(nacre_validate::validate(&replayed).is_empty());

    // The orphan on its own.
    let only_datum = replay(&[datum(sp)]).expect("a datum alone replays");
    assert_eq!(only_datum.surface_count(), 4, "three seeds and the datum");
    assert!(
        nacre_validate::validate(&only_datum).is_empty(),
        "a plane nothing sits on is not a violation — the seeds are orphans too"
    );
}
