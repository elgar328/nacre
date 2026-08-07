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

// ─────────────────────────────────────────────────────────────────────────────
// `Offset` — d away, said inside the plane's own frame
// ─────────────────────────────────────────────────────────────────────────────

fn offset(frame: nacre_ops::SketchFrame, dist: f64) -> Operation {
    Operation::DatumPlane {
        def: DatumDef::Offset { frame, dist },
    }
}

/// A box, and the frame of its top cap — a face frame, so `+d` means outward.
fn box_top_frame(m: &mut Model) -> nacre_ops::SketchFrame {
    let OpOutput::Extrude { faces, .. } = apply(
        m,
        &Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("a box") else {
        unreachable!()
    };
    nacre_ops::face_sketch_frame(m, faces[1]).expect("the top cap has a frame")
}

/// ★★ **An offset of a world-liftable plane is said in the world, so it is the same handle as
/// every other statement of that plane.**
///
/// This is the normalization that keeps "same plane, same handle" true. `push_plane` keys on
/// `(name, motion)`, so stating `z = 1` under a frame node would file it away from the `z = 1` a
/// box's cap already occupies — one geometry, two handles, and flush contact (which is decided by
/// comparing handles) quietly stops recognizing itself.
#[test]
fn an_offset_of_a_world_plane_is_the_plane_the_world_already_names() {
    let mut m = Model::new();
    let frame = box_top_frame(&mut m); // the cap on z = 1
    let cap_plane = frame.plane();

    // Down by one: z = 0, which the box's own base cap already occupies (the seeded XY plane).
    let OpOutput::DatumPlane { plane, .. } = apply(&mut m, &offset(frame, -1.0)).expect("offset")
    else {
        unreachable!()
    };
    assert_eq!(
        plane,
        m.world_plane(Axis::Z),
        "z = 0 is the seeded XY plane, however it was reached"
    );
    let SurfaceTruth::Plane { motion, .. } = m.surface_truth(plane) else {
        unreachable!()
    };
    assert_eq!(*motion, None, "a world-liftable offset makes no frame node");

    // Up by one from the cap: a genuinely new plane at z = 2, still stated in the world.
    let before = m.surface_count();
    let OpOutput::DatumPlane { plane: up, .. } =
        apply(&mut m, &offset(frame, 1.0)).expect("offset")
    else {
        unreachable!()
    };
    assert_eq!(m.surface_count(), before + 1, "z = 2 is new");
    assert_ne!(up, cap_plane);
    let SurfaceTruth::Plane { points, motion } = m.surface_truth(up) else {
        unreachable!()
    };
    assert_eq!(*motion, None);
    let PlanePoints::Known(pts) = points;
    for p in pts {
        assert_eq!(p[2].to_f64(), 2.0, "every recorded point is on z = 2");
    }
}

/// **Zero is a named reject** — in *both* roads, which is the half that is easy to miss.
///
/// On a world-liftable frame a zero offset would simply intern back to its own plane and look
/// harmless. On a genuinely tilted frame it would not: the base plane is rational in the world, so
/// the node-omission road does not apply, and the offset would be filed under
/// `([0,0,1,0], Some(node))` while the plane itself sits at `(name_world, None)` — the same
/// geometry, two handles. Rejecting zero everywhere is what keeps the two roads honest with each
/// other.
#[test]
fn a_zero_offset_is_rejected_on_both_roads() {
    let mut m = Model::new();
    let world = box_top_frame(&mut m);
    assert!(matches!(
        apply(&mut m, &offset(world, 0.0)),
        Err(OpError::ZeroOffset)
    ));

    let tilted = tilted_face_frame(&mut m);
    let before = m.surface_count();
    assert!(matches!(
        apply(&mut m, &offset(tilted, 0.0)),
        Err(OpError::ZeroOffset)
    ));
    assert_eq!(m.surface_count(), before, "a rejected offset mints nothing");
}

/// A prism on a tilted plane, and the frame of the face the sketch sat on.
fn tilted_face_frame(m: &mut Model) -> nacre_ops::SketchFrame {
    let sp = SketchPlane::from_origin_normal(
        Point3::from_array([10.0, 10.0, 10.0]),
        Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .expect("a tilted plane");
    let OpOutput::Extrude { faces, .. } = apply(
        m,
        &Operation::Extrude {
            plane: sp,
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("a tilted prism") else {
        unreachable!()
    };
    nacre_ops::face_sketch_frame(m, faces[1]).expect("the far cap has a frame")
}

/// ★★ **On a tilted face the offset is stated in the frame, exactly** — the point of the whole
/// variant. `p + d·n̂` is irrational in the world because `n̂` carries a square root; `w = d` in the
/// frame is three written decimals.
#[test]
fn a_tilted_offset_is_exact_inside_the_frame() {
    let mut m = Model::new();
    let frame = tilted_face_frame(&mut m);
    let OpOutput::DatumPlane { plane, frame: out } =
        apply(&mut m, &offset(frame, 0.5)).expect("offset")
    else {
        unreachable!()
    };

    let SurfaceTruth::Plane { points, motion } = m.surface_truth(plane) else {
        unreachable!()
    };
    assert!(
        motion.is_some(),
        "a tilted offset is written under a frame node"
    );
    let PlanePoints::Known(pts) = points;
    let half = nacre_scalar::Rat::from_decimal(0.5).unwrap();
    let (zero, one) = (
        nacre_scalar::Rat::from_int(0),
        nacre_scalar::Rat::from_int(1),
    );
    assert_eq!(
        *pts,
        [[zero, zero, half], [one, zero, half], [zero, one, half]],
        "the frame coordinates are exactly (0,0,d),(1,0,d),(0,1,d)"
    );

    // The realized geometry is d away from the base plane, measured against the base's own cache.
    let nacre_geom::Surface::Plane(base_pl) = m.surface(frame.plane()) else {
        unreachable!()
    };
    let nacre_geom::Surface::Plane(off_pl) = m.surface(plane) else {
        unreachable!()
    };
    let gap = base_pl.distance(off_pl.origin());
    println!(
        "stat tilted_offset requested=0.5 realized={gap:.17e} dev={:.3e}",
        (gap - 0.5).abs()
    );
    assert!(
        (gap - 0.5).abs() < 1e-12,
        "the offset lands where it was asked: {gap}"
    );

    // ★ The offset plane's own frame is the base frame slid along ŵ: its coefficients in the base
    // frame are `[0,0,1,−d]`, and the arbitrary-axis rule's vertical branch (`ŷ × ẑ = x̂`) hands
    // back the base `û`. Without that, a sketch on an offset plane would find its `(0,0)`
    // somewhere else entirely.
    assert_eq!(out.plane(), plane);
    assert!(matches!(out.placement(), FramePlacement::Canonical));
}

/// ★★ **Two positive controls for the normalization.**
///
/// (i) `flip` is folded into the sign, so the same `d` through frames facing opposite ways names
/// the two *different* planes either side — proof the fold is load-bearing and not decoration.
/// (ii) Two different frames of one plane, meaning the same side, name **one** handle — proof the
/// identity is keyed on `(plane, signed distance)` and not on which frame happened to be held.
#[test]
fn an_offsets_identity_is_the_plane_and_the_signed_distance() {
    let mut m = Model::new();
    let frame = box_top_frame(&mut m);
    let plane = frame.plane();

    // (i) opposite sides are different planes.
    let up = match apply(&mut m, &offset(frame, 1.0)).expect("up") {
        OpOutput::DatumPlane { plane, .. } => plane,
        _ => unreachable!(),
    };
    let down = match apply(&mut m, &offset(frame, -1.0)).expect("down") {
        OpOutput::DatumPlane { plane, .. } => plane,
        _ => unreachable!(),
    };
    assert_ne!(up, down, "±d must name the two sides, not one plane");

    // (ii) a differently-placed frame on the same plane, same side, same distance.
    let named = nacre_ops::SketchFrame::named(
        &m,
        plane,
        Point3::from_array([2.0, 2.0, 1.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .expect("a named frame on the cap");
    assert_ne!(
        named.placement(),
        frame.placement(),
        "the control is vacuous unless the two frames really differ"
    );
    let before = m.surface_count();
    let again = match apply(&mut m, &offset(named, 1.0)).expect("up again") {
        OpOutput::DatumPlane { plane, .. } => plane,
        _ => unreachable!(),
    };
    assert_eq!(
        again, up,
        "one plane, one handle — an origin and a +u do not move a parallel plane"
    );
    assert_eq!(
        m.surface_count(),
        before,
        "and nothing was minted to say it twice"
    );
}

/// A log whose offset names a plane by handle replays — and a handle past the end is a named
/// reject in both profiles, not a panic.
#[test]
fn an_offset_log_replays_and_a_stale_handle_is_rejected_by_name() {
    let mut m = Model::new();
    let frame = box_top_frame(&mut m);
    let up = offset(frame, 1.0);
    apply(&mut m, &up).expect("offset");
    m.rebuild_adjacency();

    let log = vec![
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
        up,
    ];
    let replayed = replay(&log).expect("an offset log replays");
    assert_eq!(replayed.surface_count(), m.surface_count());
    assert_eq!(replayed.live_solids.to_vec(), m.live_solids.to_vec());
    assert!(nacre_validate::validate(&replayed).is_empty());

    // The same operation with nothing built before it names a surface that does not exist.
    let alone = replay(&log[1..]);
    match alone {
        Err(OpError::LogHandleOutOfRange { cell, .. }) => {
            assert_eq!(cell, nacre_ops::LogCell::Surface)
        }
        other => panic!(
            "expected a named reject, got {:?}",
            other.map(|_| "a model")
        ),
    }
}

/// ★★ **The `flip` fold, proved load-bearing.**
///
/// Two boxes stacked on `z = 1` give two faces on **one plane** whose outward normals oppose: the
/// lower box's top cap and the upper box's base cap. `face_sketch_frame` measures `flip` against
/// each face's outward, so the two frames name the same `Handle<Surface>` and disagree about which
/// way `ŵ` runs. The *same* `+d` through them must therefore name the two different planes either
/// side — which is what "`+dist` runs along the frame's `ŵ`" promises a caller, and what folding
/// `flip` into the sign is for.
///
/// Without the fold both would build in the canonical frame and land on the same plane, and a
/// caller asking for "0.5 outward" from the underside would get 0.5 the other way. The pairing
/// with [`an_offsets_identity_is_the_plane_and_the_signed_distance`] is deliberate: that one shows
/// what must *not* distinguish two frames (origin, `+u`), this one shows what must.
#[test]
fn the_flip_fold_decides_which_side_and_it_bites() {
    let mut m = Model::new();
    let lower = box_top_frame(&mut m); // outward +z on z = 1
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            plane: SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 1.0])),
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("a box stacked on the first") else {
        unreachable!()
    };
    let upper = nacre_ops::face_sketch_frame(&m, faces[0]).expect("its base cap"); // outward −z

    assert_eq!(
        lower.plane(),
        upper.plane(),
        "the two faces share one plane handle, or the control proves nothing"
    );
    assert_ne!(
        lower.flip(),
        upper.flip(),
        "their outward normals oppose, so the measured flips must differ"
    );

    let side = |f: nacre_ops::SketchFrame, m: &mut Model| match apply(m, &offset(f, 0.5)) {
        Ok(OpOutput::DatumPlane { plane, .. }) => plane,
        other => panic!("offset failed: {other:?}"),
    };
    let a = side(lower, &mut m);
    let b = side(upper, &mut m);
    assert_ne!(
        a, b,
        "the same +d through opposed frames must name the two sides — the fold is what does it"
    );

    // And the geometry says which is which: +d from the +z face is above, from the −z face below.
    let z_of = |h| {
        let nacre_geom::Surface::Plane(p) = m.surface(h) else {
            unreachable!()
        };
        p.origin().as_array()[2]
    };
    assert!(z_of(a) > 1.0, "outward from the top cap is up: {}", z_of(a));
    assert!(
        z_of(b) < 1.0,
        "outward from the base cap is down: {}",
        z_of(b)
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Which way a frame faces — the half a plane handle cannot carry
// ─────────────────────────────────────────────────────────────────────────────

/// Realize a frame's world basis the way every consumer does, through a prism it builds.
///
/// `frame_world_basis` is crate-private, so this reads the direction back the way an operation
/// would: extrude a unit square on the frame and look at which side the far cap landed on.
fn frame_sweep_direction(m: &mut Model, frame: nacre_ops::SketchFrame) -> [f64; 3] {
    let before = m.surface_count();
    let OpOutput::DatumPlane { plane: up, .. } = apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Offset { frame, dist: 1.0 },
        },
    )
    .expect("a unit offset names the side the frame faces") else {
        unreachable!()
    };
    let _ = before;
    let nacre_geom::Surface::Plane(base) = m.surface(frame.plane()) else {
        unreachable!()
    };
    let nacre_geom::Surface::Plane(off) = m.surface(up) else {
        unreachable!()
    };
    let (a, b) = (base.origin().as_array(), off.origin().as_array());
    core::array::from_fn(|k| b[k] - a[k])
}

/// ★★★ **The measurement this whole stage rests on.**
///
/// A plane's canonical name has no direction and planes intern, so stating `z = 0` facing `+ẑ`
/// and stating it facing `−ẑ` produce **one handle**. If the frame handed back carried
/// `flip: false` in both cases, the second caller would be answered "up" when they said "down" —
/// and an operation that took a frame where it used to take a plane would sweep the wrong way.
///
/// So `DatumPlane` measures `flip` against the normal the caller's own point order fixes. Same
/// handle, opposite flips, and each frame faces the way its caller stated.
#[test]
fn two_statements_of_one_plane_share_a_handle_and_keep_their_directions() {
    let mut m = Model::new();
    let up = SketchPlane::through_points(
        Point3::origin(),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    )
    .expect("z = 0 facing +z");
    let down = SketchPlane::through_points(
        Point3::origin(),
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
    )
    .expect("z = 0 facing −z");
    assert!(
        up.normal().as_array()[2] > 0.0 && down.normal().as_array()[2] < 0.0,
        "the fixture is vacuous unless the two statements really oppose"
    );

    let state = |m: &mut Model, sp: SketchPlane| match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(sp),
        },
    ) {
        Ok(OpOutput::DatumPlane { plane, frame }) => (plane, frame),
        other => panic!("datum failed: {other:?}"),
    };
    let (hu, fu) = state(&mut m, up);
    let (hd, fd) = state(&mut m, down);

    assert_eq!(
        hu, hd,
        "one geometric plane, one handle — interning answered"
    );
    assert_ne!(
        fu.flip(),
        fd.flip(),
        "and the two frames must disagree about which way ŵ runs, or the direction was lost"
    );

    // Read the directions back through an operation, not through a private helper.
    let (du, dd) = (
        frame_sweep_direction(&mut m, fu),
        frame_sweep_direction(&mut m, fd),
    );
    assert!(du[2] > 0.0, "the +z statement's frame faces up: {du:?}");
    assert!(dd[2] < 0.0, "the −z statement's frame faces down: {dd:?}");
}

/// The seed sugar and a datum of the same world plane **mean** the same frame.
///
/// ★ They are not the same *value*, and that is worth knowing rather than asserting away: the
/// sugar spells XY and YZ `Canonical` (derived) while a datum always spells `Named` (the caller's
/// own `points[0]` and `points[1] − points[0]`). Two spellings of one frame — which is harmless
/// here because a world plane's basis lifts exactly, so **neither road builds a motion node** and
/// the arena cannot tell them apart. On a *tilted* plane the two spellings would be two nodes,
/// which is correct: there they really are two different sketch coordinate systems.
///
/// What must agree is the meaning: same plane, same `flip`, same realized `ŵ`. The axes are
/// pinned separately by [`the_zx_sugar_is_not_the_canonical_frame`].
#[test]
fn the_world_sugar_means_what_a_datum_of_the_same_plane_means() {
    for (axis, stated) in [
        (Axis::Z, SketchPlane::world_xy()),
        (Axis::X, SketchPlane::world_yz()),
        (Axis::Y, SketchPlane::world_zx()),
    ] {
        let mut m = Model::new();
        let sugar = nacre_ops::SketchFrame::world(&m, axis);
        assert_eq!(sugar.plane(), m.world_plane(axis));
        assert!(!sugar.flip(), "{axis:?}: the sugar states no flip");

        let OpOutput::DatumPlane { frame, .. } = apply(
            &mut m,
            &Operation::DatumPlane {
                def: DatumDef::Stated(stated),
            },
        )
        .expect("stating a world plane") else {
            unreachable!()
        };
        assert_eq!(frame.plane(), sugar.plane(), "{axis:?}: same plane");
        assert_eq!(frame.flip(), sugar.flip(), "{axis:?}: same sense");

        let (a, b) = (
            frame_sweep_direction(&mut m, sugar),
            frame_sweep_direction(&mut m, frame),
        );
        for k in 0..3 {
            assert!(
                (a[k] - b[k]).abs() < 1e-12,
                "{axis:?}: the two spellings must face the same way — {a:?} vs {b:?}"
            );
        }
    }
}

/// ★★ **The ZX trap, pinned as a positive control.**
///
/// `SketchFrame::canonical` on the ZX seed gives `+u = −x̂` (the arbitrary-axis derivation) while
/// the convention gives `+u = +ẑ`. They are *different frames on one plane*, and the sugar exists
/// precisely because of that. Without this assertion nothing would record that the two differ,
/// and a later simplification ("the sugar is just `canonical`") would pass every other test.
#[test]
fn the_zx_sugar_is_not_the_canonical_frame() {
    let m = Model::new();
    let zx = m.world_plane(Axis::Y);
    let sugar = nacre_ops::SketchFrame::world(&m, Axis::Y);
    let derived = nacre_ops::SketchFrame::canonical(zx);
    assert_eq!(sugar.plane(), derived.plane(), "same plane");
    assert_ne!(
        sugar.placement(),
        derived.placement(),
        "ZX is the axis where the convention cannot be derived — if these ever agree, either the \
         seeding or the arbitrary-axis rule moved, and every ZX sketch turned with it"
    );

    // XY and YZ are the other half of the claim: there, derived *is* the convention.
    for axis in [Axis::Z, Axis::X] {
        assert_eq!(
            nacre_ops::SketchFrame::world(&m, axis),
            nacre_ops::SketchFrame::canonical(m.world_plane(axis)),
            "{axis:?}: derived and stated agree here, so the sugar must not invent a placement"
        );
    }
}

/// ★ The world frames face `+axis`, and that is a measurement.
///
/// Each seed's canonical `ŵ` realizes to `+axis` because of how `Model::new` writes the seeds'
/// points and how the canonical name normalizes — not because anything derives it. Re-seed
/// differently and `SketchFrame::world` turns silently; this is what stops that.
#[test]
fn a_world_frame_faces_its_axis() {
    for (axis, expect) in [
        (Axis::Z, [0.0, 0.0, 1.0]),
        (Axis::X, [1.0, 0.0, 0.0]),
        (Axis::Y, [0.0, 1.0, 0.0]),
    ] {
        let mut m = Model::new();
        let f = nacre_ops::SketchFrame::world(&m, axis);
        let d = frame_sweep_direction(&mut m, f);
        for k in 0..3 {
            assert!(
                (d[k] - expect[k]).abs() < 1e-12,
                "{axis:?}: a unit offset along the frame should land at {expect:?}, got {d:?}"
            );
        }
    }
}

/// The flip a datum measures and the sign `Offset` folds must agree, or an offset from a
/// "downward" plane would climb.
#[test]
fn the_datums_flip_and_the_offsets_fold_agree() {
    let mut m = Model::new();
    let down = SketchPlane::from_origin_normal(
        Point3::from_array([0.0, 0.0, 2.0]),
        Vector3::from_array([0.0, 0.0, -1.0]),
    )
    .expect("z = 2 facing down");
    let OpOutput::DatumPlane { frame, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(down),
        },
    )
    .expect("stated") else {
        unreachable!()
    };
    assert!(
        frame.flip(),
        "the canonical ŵ of z = 2 is +ẑ, so a −ẑ statement flips"
    );
    let d = frame_sweep_direction(&mut m, frame);
    assert!(
        d[2] < 0.0,
        "+dist from a downward-facing frame must go down: {d:?}"
    );
}
