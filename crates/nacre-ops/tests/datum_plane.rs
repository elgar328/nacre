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
use nacre_ops::SketchFrame;
use nacre_ops::{DatumDef, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, replay};
use nacre_scalar::Axis;
use nacre_topo::{FramePlacement, Model, PlanePoints, SurfaceTruth};

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

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
/// its base cap as `−plane.normal()`, and `WorkingPlane::frame_sign` reads a stored normal against a
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
    let PlanePoints::Known(pts) = points else {
        panic!("a stated datum records its points by value, not by handle")
    };
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

/// ★★★★★ **A plane the model cannot name cannot host a sketch — and that is the wall the last
/// migration row actually stands behind.**
///
/// `docs/truth-and-cache.md` recorded S5(ii)-2 as judging-layer work. It is not. A datum through
/// vertices in mixed frames has no exact name, and from there:
///
/// ```text
/// no name → frame_chain declines → no SketchFrame → never a base cap → never in a plane table
/// ```
///
/// so `WorkingPlaneDef::Through` has no producer and cannot have one until a nameless plane can be
/// framed. That was read out of the code; this runs it.
///
/// ★ **Not a backstop for a future stage.** `Model::push_plane` is `pub` and leaves a nameless
/// plane behind for a collinear triple (`intern_plane` builds no key), and `SketchFrame::canonical`
/// is `pub` and checks nothing — so the pair is reachable from outside this crate today and had no
/// test. Production is safe because every push passes a non-collinearity gate first (`PlaneDef`'s
/// constructors call `plane_name_exact`; the prism's caps and walls gate on
/// `Plane::from_point_normal`/`through_points` and `Profile2d::check`), which is what
/// `every_plane_still_has_a_name` observes — but that is *we do not*, not *it cannot*.
///
/// ★★ **The control is half the test.** An unregistered plane's name is *derivable* and merely
/// absent from the side table; a mixed-frame `Through` plane's name **does not exist**. Different
/// facts — but indistinguishable *to a frame*, because `frame_chain`'s first line and
/// `SketchFrame::named` both read that one table and there is no other route. Pushing the same
/// points through `push_plane`, where the only thing that changes is that the name gets recorded,
/// is what says so. A seed plane would not: it would answer "it extruded because it is world XY".
#[test]
fn a_plane_with_no_name_cannot_host_a_sketch() {
    let r = nacre_scalar::Rat::from_int;
    let pts = [[r(0), r(0), r(3)], [r(1), r(0), r(3)], [r(0), r(1), r(3)]];
    let origin = Point3::from_array([0.0, 0.0, 3.0]);
    // `−normal`, the convention every datum and base cap keeps.
    let cache = nacre_geom::Plane::from_point_normal(origin, Vector3::from_array([0.0, 0.0, -1.0]))
        .expect("z = 3 is a plane");
    let sketch = |frame| Operation::Extrude {
        frame,
        profile: square(0.0, 2.0),
        dist: 1.0,
    };

    let mut m = Model::new();
    let nameless = m.push_plane_unregistered(cache, pts);
    assert!(
        !m.surface_name.contains_key(&nameless),
        "the instrument did not produce a nameless plane"
    );

    let before = m.live_solids.clone();
    assert_eq!(
        apply(&mut m, &sketch(SketchFrame::canonical(nameless))).unwrap_err(),
        OpError::PlaneWithoutExactForm,
        "a sketch on a plane with no name must be a named reject"
    );
    assert_eq!(
        m.live_solids, before,
        "a refusal must leave the model alone"
    );

    // ★ The other door onto the same plane answers the same way. One failure, one description —
    // two descriptions of one thing is the shape this kernel has been bitten by repeatedly.
    assert_eq!(
        SketchFrame::named(&m, nameless, origin, Vector3::from_array([1.0, 0.0, 0.0])).unwrap_err(),
        OpError::PlaneWithoutExactForm,
        "the two doors onto a nameless plane describe the failure differently"
    );

    // ★★ The control: same cache, same points, and the name recorded.
    let mut named = Model::new();
    let (h, _) = named.push_plane(cache, pts, None);
    assert!(
        named.surface_name.contains_key(&h),
        "push_plane derives the name from the points"
    );
    apply(&mut named, &sketch(SketchFrame::canonical(h)))
        .expect("the same plane, named, hosts the same sketch");
}

/// ★★ **The seam this stage opened, now closed: an extrude's base cap *is* the frame's plane.**
///
/// Before S5(i)-b this compared two roads — "state the plane first" against "let the extrude
/// create it" — and measured that they met at one handle. That comparison no longer exists,
/// because an extrude names its plane and there is no road where the plane is absent. What
/// remains, and is still not free, is the claim underneath it: the operation **reuses** the
/// handle it was given rather than pushing a second plane on the same geometry, and stating that
/// plane again afterwards mints nothing.
///
/// (The other half — that which point anchored the plane's cache does not move the result — is
/// `tests/plane_anchor.rs`, where both arms are now production roads.)
#[test]
fn an_extrudes_base_cap_is_the_frame_it_was_given() {
    let sp = SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5]));
    let mut m = Model::new();

    let OpOutput::DatumPlane { plane: h, frame } = apply(&mut m, &datum(sp)).expect("stated")
    else {
        unreachable!()
    };
    let after_datum = m.surface_count();

    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("extrude on the stated plane") else {
        unreachable!()
    };
    assert_eq!(
        m.faces.get(faces[0]).surface,
        h,
        "the base cap is the plane the frame named, not a second one on the same geometry"
    );
    // ★ Three, not six. The prism wants a base cap, a top cap and four walls; the base is the
    // handle it was given, and two of the walls (`x = 0`, `y = 0`) are the **seeded** YZ and ZX
    // planes. Interning answers for those the same way it answers for the base — which is the
    // point: a plane is a thing the model has, not a thing each operation restates.
    assert_eq!(
        m.surface_count() - after_datum,
        3,
        "top cap + two new walls; the base and the two axis walls were already there"
    );

    // Saying the same plane again is interning's job, not the arena's.
    let before = m.surface_count();
    let OpOutput::DatumPlane { plane: again, .. } = apply(&mut m, &datum(sp)).expect("again")
    else {
        unreachable!()
    };
    assert_eq!(again, h);
    assert_eq!(
        m.surface_count(),
        before,
        "one plane, one handle, no growth"
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
    // The log is assembled against the model it is applied to — a frame names a handle, and a
    // handle only means something in its own arena until `replay` re-anchors it.
    let mut scratch = Model::new();
    let log = vec![
        datum(sp),
        Operation::Extrude {
            frame: SketchFrame::world(&scratch, Axis::Z),
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    ];
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
            frame: SketchFrame::world(m, Axis::Z),
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
    let PlanePoints::Known(pts) = points else {
        panic!("a stated datum records its points by value, not by handle")
    };
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
    let __g209 = datum_frame(m, sp);
    let OpOutput::Extrude { faces, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __g209,
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
    let PlanePoints::Known(pts) = points else {
        panic!("a stated datum records its points by value, not by handle")
    };
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
            frame: SketchFrame::world(&m, Axis::Z),
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
    let __f106 = datum_frame(
        &mut m,
        SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 1.0])),
    );
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f106,
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

// =================================================================================================
// A datum that names vertices instead of coordinates (S5(ii), stage 1)
// =================================================================================================

/// Three corners of a **tilted** prism whose vertices are discovered — the population where the
/// coordinate road is measured to produce a different plane (`tests/point_width.rs`).
fn tilted_prism_with_pocket() -> Model {
    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let mut m = Model::new();
    let plane = SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let frame = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(0.1111111111111111, 0.1234567890123456),
                p2(4.123456789012345, 0.2345678901234567),
                p2(3.9876543210987654, 3.1234567890123459),
                p2(0.2222222222222222, 2.765432109876543),
            ])
            .unwrap(),
            dist: 2.5,
        },
    )
    .expect("the base prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let wall = *m
        .shells
        .get(m.solids.get(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            let s = m.faces.get(f).surface;
            m.surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
        })
        .expect("the tilted-wall population");
    let sp = nacre_ops::face_plane(&m, wall).expect("planar");
    let d = nacre_props::face_props(&m, wall).unwrap().centroid - sp.origin();
    let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
    apply(
        &mut m,
        &Operation::PocketOnFace {
            face: wall,
            profile: Profile2d::polygon(vec![
                p2(cu - 0.3, cv - 0.3),
                p2(cu + 0.3, cv - 0.3),
                p2(cu + 0.3, cv + 0.3),
                p2(cu - 0.3, cv + 0.3),
            ])
            .unwrap(),
            dist: 0.4,
        },
    )
    .expect("the pocket");
    m.rebuild_adjacency();
    m
}

/// Every live vertex of `m`, in walk order.
fn live_verts(m: &Model) -> Vec<nacre_store::Handle<nacre_topo::Vertex>> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &s in &m.live_solids {
        let sol = m.solids.get(s);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                let f = m.faces.get(fh);
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for &he in &lp.half_edges {
                        let vh = m.he_start(he);
                        if seen.insert(vh) {
                            out.push(vh);
                        }
                    }
                }
            }
        }
    }
    out
}

/// Three vertices that share one frame, solve, **and name a plane the model does not already
/// hold**.
///
/// ★★★ That last clause is not fussiness — without it the fixture measures nothing. Three corners
/// of one face are coplanar with that face, so `push_plane_through` interns onto the existing
/// handle and the stored truth stays `Known`: the datum works, but the storing path never runs.
/// (Which is also the design saying what it says — *"interning is deterministic, the first pusher
/// wins"* — so the fixture has to reach across faces to see the new variant at all.)
fn three_solvable(m: &Model) -> [nacre_store::Handle<nacre_topo::Vertex>; 3] {
    let mut ok = Vec::new();
    for vh in live_verts(m) {
        if m.through_points_rat([vh, vh, vh]).is_some() {
            ok.push(vh);
        }
    }
    for i in 0..ok.len() {
        for j in (i + 1)..ok.len() {
            for k in (j + 1)..ok.len() {
                let t = [ok[i], ok[j], ok[k]];
                let mut sorted = t;
                sorted.sort_by_key(|v| v.index());
                let Some(name) = m.plane_name_through(sorted) else {
                    continue;
                };
                let mut n = 0u32;
                let mut already = false;
                while let Some(h) = m.surface_handle_at(n) {
                    n += 1;
                    already |= m.surface_name.get(&h) == Some(&name);
                }
                if !already {
                    return t;
                }
            }
        }
    }
    panic!("no solvable non-collinear triple in this fixture")
}

/// ★★★★★ **The capability, closed.** The datum built by *naming* three vertices is the plane
/// actually through them — and it is **not** the plane their coordinates produce.
///
/// `tests/point_width.rs` measured the second half of that: on this population, all 220 triples
/// spelled in coordinates name a different plane. This is the other side — the same three
/// vertices, named, landing on the exact one.
#[test]
fn a_datum_through_vertices_is_the_plane_the_coordinates_miss() {
    let mut m = tilted_prism_with_pocket();
    let vs = three_solvable(&m);
    let mut sorted = vs;
    sorted.sort_by_key(|v| v.index());
    let exact = m.plane_name_through(sorted).expect("nameable");

    let OpOutput::DatumPlane { plane, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("a datum through three solvable vertices") else {
        unreachable!()
    };
    assert_eq!(
        m.surface_name.get(&plane),
        Some(&exact),
        "the datum must be the exact plane through those vertices"
    );

    // ★ The negative control that gives the assertion teeth: the same three vertices *spelled in
    // coordinates* land somewhere else. Without this the test would pass on a kernel that quietly
    // rounded, since both roads would then agree.
    let lift = |vh| {
        let a = m
            .vertex_point(vh)
            .as_array()
            .map(nacre_scalar::Rat::from_decimal);
        [a[0].unwrap(), a[1].unwrap(), a[2].unwrap()]
    };
    let rounded =
        nacre_scalar::plane_name_exact(lift(vs[0]), lift(vs[1]), lift(vs[2])).expect("nameable");
    assert_ne!(
        exact, rounded,
        "this fixture cannot show the gap — the coordinate road happened to agree"
    );
}

/// ★★★ **One plane, two orders, opposite sides.** Sorting keeps identity; the caller's order
/// keeps direction. Asserting only "one handle" would measure half of it and pass on a kernel
/// where the caller cannot choose a side at all.
#[test]
fn reversing_the_vertex_order_keeps_the_handle_and_flips_the_frame() {
    let mut m = tilted_prism_with_pocket();
    let [a, b, c] = three_solvable(&m);

    let datum = |m: &mut Model, vs: [nacre_store::Handle<nacre_topo::Vertex>; 3]| {
        let OpOutput::DatumPlane { plane, frame } = apply(
            m,
            &Operation::DatumPlane {
                def: DatumDef::ThroughVertices(vs),
            },
        )
        .expect("datum") else {
            unreachable!()
        };
        (plane, frame)
    };
    let (p1, f1) = datum(&mut m, [a, b, c]);
    let (p2, f2) = datum(&mut m, [a, c, b]);
    assert_eq!(p1, p2, "the same three vertices are one plane");
    assert_ne!(
        f1.flip(),
        f2.flip(),
        "the caller's order is the only place direction can live"
    );
}

/// ★★ **Each refusal has its own name, and each name has a fixture that fires it.** A cause with
/// no fixture is a cause that was never measured — and `VerticesInMixedFrames` in particular is
/// how the next stage will count what it is opening.
///
/// ★ `VertexPointTooWide` is deliberately absent: reaching it needs a vertex whose carriers are
/// all named yet whose meeting point leaves `Rat`, and this kernel's corpus produces no wide
/// carrier at all (measured `wide_carrier = 0`). Recorded as unfired rather than faked — the
/// stage that removes the requirement removes the variant with it.
#[test]
fn a_datum_through_vertices_refuses_by_cause() {
    let mut m = tilted_prism_with_pocket();
    let [a, b, _] = three_solvable(&m);
    let before = m.live_solids.clone();
    let through = |m: &mut Model, vs: [nacre_store::Handle<nacre_topo::Vertex>; 3]| {
        apply(
            m,
            &Operation::DatumPlane {
                def: DatumDef::ThroughVertices(vs),
            },
        )
        .unwrap_err()
    };

    assert_eq!(through(&mut m, [a, b, a]), OpError::DuplicateVertex);
    assert_eq!(
        m.live_solids, before,
        "a refusal must leave the model alone"
    );

    // ★ A cylinder's seam vertex is the reachable `VertexNotThreePlane` case: its pair pins a
    // curve, so no three planes name it and no exact coordinate follows from its definition. The
    // box beside it supplies the other two — a cylinder alone has only the two seam vertices.
    let mut cy = Model::new();
    cy.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    cy.add_cuboid(
        Point3::from_array([4.0, 0.0, 0.0]),
        Point3::from_array([5.0, 1.0, 1.0]),
    );
    cy.rebuild_adjacency();
    let seam = live_verts(&cy)
        .into_iter()
        .find(|v| matches!(cy.vertices.get(*v).def, nacre_topo::VertexDef::OnSeam(_)))
        .expect("a cylinder has seam vertices");
    let corners: Vec<_> = live_verts(&cy)
        .into_iter()
        .filter(|v| {
            matches!(
                cy.vertices.get(*v).def,
                nacre_topo::VertexDef::ThreePlane(_)
            )
        })
        .collect();
    assert_eq!(
        through(&mut cy, [seam, corners[0], corners[1]]),
        OpError::VertexNotThreePlane
    );

    // ★★ Two boxes stacked share a vertical line, so three of their corners are collinear — the
    // one shape in this vocabulary that puts three vertices on a line.
    let mut st = Model::new();
    st.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    st.add_cuboid(
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    st.rebuild_adjacency();
    let on_axis: Vec<_> = live_verts(&st)
        .into_iter()
        .filter(|v| {
            let c = st.vertex_point(*v).as_array();
            c[0] == 0.0 && c[1] == 0.0
        })
        .collect();
    assert!(on_axis.len() >= 3, "the stack must share a corner line");
    assert_eq!(
        through(&mut st, [on_axis[0], on_axis[1], on_axis[2]]),
        OpError::CollinearVertices
    );

    // ★ The pure-vertices-in-different-frames shape is **no longer on this list** — open item
    // 16's first wall opened it, and `a_nameless_datum_hosts_a_sketch_end_to_end` is where it
    // now lives as an acceptance. What stays refused under `VerticesInMixedFrames` is a single
    // vertex whose own carriers straddle (`a_datum_on_straddling_carriers_has_no_name`).
}

/// ★★★★★ **The first wall of open item 16, opened end to end**: a datum through three vertices
/// that are each exact in their own frame — but not in each other's — gets a *judged* frame,
/// hosts a sketch, and the whole thing replays.
///
/// The population S5(ii)-1 refused with `VerticesInMixedFrames` splits: this (the caller's
/// vertices differ) is now accepted; a single straddling vertex still refuses. The plane has
/// **no name** — its exact world coefficients are irrational — so it interns by statement, and
/// the frame is derived from the defining points as intervals at a fixed rung.
#[test]
fn a_nameless_datum_hosts_a_sketch_end_to_end() {
    let mut mx = Model::new();
    let fixed = mx.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let spun = mx.add_cuboid(
        Point3::from_array([4.0, 0.0, 0.0]),
        Point3::from_array([5.0, 1.0, 1.0]),
    );
    mx.rebuild_adjacency();
    let OpOutput::Transform { solid: spun } = apply(
        &mut mx,
        &Operation::Transform {
            solid: spun,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    mx.rebuild_adjacency();
    let of = |m: &Model, s| {
        let mut out = Vec::new();
        for &fh in &m.shells.get(m.solids.get(s).outer).faces {
            for lp in std::iter::once(&m.faces.get(fh).outer).chain(m.faces.get(fh).inner.iter()) {
                for &he in &lp.half_edges {
                    out.push(m.he_start(he));
                }
            }
        }
        out
    };
    let (still, turned) = (of(&mx, fixed), of(&mx, spun));
    let vs = [still[0], still[1], turned[0]];

    let OpOutput::DatumPlane { plane, frame } = apply(
        &mut mx,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("the pure-mixed datum is accepted now") else {
        unreachable!()
    };
    assert!(
        !mx.surface_name.contains_key(&plane),
        "its exact world coefficients are irrational — a name here would be an invention"
    );
    assert!(
        matches!(
            mx.surface_truth(plane),
            SurfaceTruth::Plane {
                points: PlanePoints::Through(_),
                motion: None,
            }
        ),
        "the truth is the statement: handles, and no motion of its own"
    );

    // ★ Statement interning, end to end: the same three vertices again are the same handle.
    let OpOutput::DatumPlane { plane: again, .. } = apply(
        &mut mx,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("restating is fine") else {
        unreachable!()
    };
    assert_eq!(plane, again, "one statement, one handle");

    // The sketch: an extrude on the judged frame, validated.
    let OpOutput::Extrude { .. } = apply(
        &mut mx,
        &Operation::Extrude {
            frame,
            profile: square(0.0, 1.0),
            dist: 0.5,
        },
    )
    .expect("a sketch on the judged frame") else {
        unreachable!()
    };
    mx.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&mx).is_empty(),
        "the prism on a nameless datum must be a valid closed b-rep"
    );

    // ★ And the same session rebuilt from scratch reaches the same arena — the judged road is
    // deterministic. The scratch model mints its **own** handles (using `mx`'s here would trip
    // the cross-store guard, and rightly — that is the guard working); determinism is asserted
    // on the indices and the arena count.
    let mut scratch = Model::new();
    let f2 = scratch.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let s2 = scratch.add_cuboid(
        Point3::from_array([4.0, 0.0, 0.0]),
        Point3::from_array([5.0, 1.0, 1.0]),
    );
    scratch.rebuild_adjacency();
    let OpOutput::Transform { solid: s2 } = apply(
        &mut scratch,
        &Operation::Transform {
            solid: s2,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    scratch.rebuild_adjacency();
    let (still2, turned2) = (of(&scratch, f2), of(&scratch, s2));
    let vs2 = [still2[0], still2[1], turned2[0]];
    assert_eq!(
        [vs2[0].index(), vs2[1].index(), vs2[2].index()],
        [vs[0].index(), vs[1].index(), vs[2].index()],
        "the scratch construction mints the same indices — the premise of comparing arenas"
    );
    let OpOutput::DatumPlane {
        plane: plane2,
        frame: frame2,
    } = apply(
        &mut scratch,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs2),
        },
    )
    .expect("the same statement is accepted again")
    else {
        unreachable!()
    };
    assert_eq!(
        plane2.index(),
        plane.index(),
        "same statement, same arena slot"
    );
    apply(
        &mut scratch,
        &Operation::Extrude {
            frame: frame2,
            profile: square(0.0, 1.0),
            dist: 0.5,
        },
    )
    .expect("the same sketch extrudes again");
    assert_eq!(
        scratch.surface_count(),
        mx.surface_count(),
        "same statements, same arena — minus nothing"
    );
}

/// ★★★★ **The second wall, heterogeneous half: the nameless plane's face sits in a judgment
/// table and a boolean runs over it.**
///
/// The prism's base cap *is* the nameless datum surface, so `collect_planes` meets a `Through`
/// truth with no rational solve and builds the witness triangle from the vertices' own chains —
/// three exact points whose chains differ, which the judging layer never forbade. Every wall of
/// the prism also realizes through the judged frame node (`Motion::Frame` → the nameless plane →
/// `FrameThrough`), so this exercises the recursion end to end without a separate probe.
///
/// `Common` against a box that swallows the prism is the crisp oracle: the result is the prism
/// itself, volume `1 × 1 × 0.5` **exactly as stated in the sketch frame** — a rigid frame does
/// not change volume, whatever its coefficients are.
#[test]
fn a_prism_on_a_nameless_datum_survives_a_boolean() {
    let mut m = Model::new();
    let fixed = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let spun = m.add_cuboid(
        Point3::from_array([4.0, 0.0, 0.0]),
        Point3::from_array([5.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let OpOutput::Transform { solid: spun } = apply(
        &mut m,
        &Operation::Transform {
            solid: spun,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let corners = |m: &Model, s| {
        let mut out = Vec::new();
        for &fh in &m.shells.get(m.solids.get(s).outer).faces {
            for lp in std::iter::once(&m.faces.get(fh).outer).chain(m.faces.get(fh).inner.iter()) {
                for &he in &lp.half_edges {
                    out.push(m.he_start(he));
                }
            }
        }
        out
    };
    let (still, turned) = (corners(&m, fixed), corners(&m, spun));
    let OpOutput::DatumPlane { plane, frame } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices([still[0], still[1], turned[0]]),
        },
    )
    .expect("the pure-mixed datum") else {
        unreachable!()
    };
    assert!(!m.surface_name.contains_key(&plane), "nameless, still");
    let OpOutput::Extrude { solid: prism, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square(0.0, 1.0),
            dist: 0.5,
        },
    )
    .expect("prism on the judged frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    // A box that swallows the prism wherever the judged frame put it.
    let block = m.add_cuboid(
        Point3::from_array([-12.0, -12.0, -12.0]),
        Point3::from_array([12.0, 12.0, 12.0]),
    );
    m.rebuild_adjacency();
    let out = nacre_ops::boolean(&mut m, nacre_ops::BoolKind::Common, block, prism)
        .expect("the boolean over a nameless class must answer, not panic");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 1, "block ∩ prism is the prism");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "and it is a valid closed b-rep"
    );
    let vol = nacre_props::mass_props(&m, out[0])
        .expect("measurable")
        .volume;
    assert!(
        (vol - 0.5).abs() < 1e-9,
        "1 × 1 × 0.5 wherever the frame sits — got {vol}"
    );
}

/// ★★★ **The mixed-frame cause that is not the caller's doing — one solid, one vertex, carriers in
/// two frames.**
///
/// `a_datum_through_vertices_refuses_by_cause` fires `VerticesInMixedFrames` by combining two
/// *solids*, which a caller chooses to do. The population that actually dominates is the other
/// one: a boolean between a turned operand and a still one hands back **one** solid whose corners
/// are the meeting of an unmoved wall and two turned ones, so a single vertex straddles. Measured,
/// **12 of that solid's 20 vertices** are in that state (`tests/point_width.rs`), and until this
/// there was no fixture for it — the label was only ever fired by the caller-error shape.
///
/// ★ **The control is the minimal difference**: three vertices that solve, and the same triple with
/// one of them swapped for a straddler. Nothing else about the model changes, so the swap is the
/// whole of the cause.
#[test]
fn a_datum_on_straddling_carriers_has_no_name() {
    use nacre_topo::VertexDef;

    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([3.3, 3.3, -1.0]),
        Point3::from_array([7.7, 7.7, 11.0]),
    );
    let OpOutput::Transform { solid: b } = apply(
        &mut m,
        &Operation::Transform {
            solid: b,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let cut = nacre_ops::boolean(&mut m, nacre_ops::BoolKind::Cut, a, b).expect("cut");
    m.rebuild_adjacency();
    // ★ **The boolean's own result** — the whole point is that the straddling is something the
    // *kernel* produced, not something a caller reached across two solids to make. (The untouched
    // second cuboid is still live; walking every live solid would readmit exactly the caller-error
    // shape this test is not about.) The cut severs a corner, so the result is more than one piece
    // — that count is recorded rather than asserted, since nothing here depends on it.
    let mine: std::collections::HashSet<_> = {
        let mut out = std::collections::HashSet::new();
        for &sh in cut
            .iter()
            .flat_map(|&sh| {
                let s = m.solids.get(sh);
                std::iter::once(s.outer).chain(s.cavities.iter().copied())
            })
            .collect::<Vec<_>>()
            .iter()
        {
            for &fh in &m.shells.get(sh).faces {
                let f = m.faces.get(fh);
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for &he in &lp.half_edges {
                        out.insert(m.he_start(he));
                    }
                }
            }
        }
        out
    };

    // Split them by whether their own three carriers share a motion.
    let (mut pure, mut straddling) = (Vec::new(), Vec::new());
    for vh in mine {
        let VertexDef::ThreePlane(tri) = m.vertices.get(vh).def else {
            continue;
        };
        let mine = m.plane_motion(tri[0]);
        if tri.iter().any(|h| m.plane_motion(*h) != mine) {
            straddling.push(vh);
        } else {
            pure.push(vh);
        }
    }
    assert!(
        !straddling.is_empty(),
        "a turned-against-still cut must leave carriers straddling — the fixture stopped measuring"
    );
    assert!(pure.len() >= 3, "and some corners must still be pure");

    // A triple that does solve, so the swap below has something to be the difference from.
    let solvable = (0..pure.len())
        .flat_map(|i| ((i + 1)..pure.len()).map(move |j| (i, j)))
        .flat_map(|(i, j)| ((j + 1)..pure.len()).map(move |k| [i, j, k]))
        .map(|ix| ix.map(|i| pure[i]))
        .find(|t| m.through_points_rat(*t).is_some())
        .expect("three pure vertices that share one frame");

    // ★★★★★ **16-2: the straddling vertex is accepted** — its definition (the meet of its three
    // carriers) is complete even though its coordinate exists in no frame at all. This block used
    // to assert `VerticesInMixedFrames`; the acceptance below is the same fixture with the wall
    // moved, end to end: datum → sketch → extrude → validate → deterministic rebuild.
    let swapped = [straddling[0], solvable[1], solvable[2]];
    assert!(
        m.through_points_rat(swapped).is_none(),
        "a straddling vertex has no rational coordinate in any single frame"
    );
    let OpOutput::DatumPlane { plane, frame } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(swapped),
        },
    )
    .expect("a straddling-vertex datum is accepted now (16-2)") else {
        unreachable!()
    };
    assert!(
        !m.surface_name.contains_key(&plane),
        "its exact world coefficients are irrational — nameless, by statement"
    );
    // Restating is the same handle — statement interning holds for meets too.
    let OpOutput::DatumPlane { plane: again, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(swapped),
        },
    )
    .expect("restating") else {
        unreachable!()
    };
    assert_eq!(plane, again, "one statement, one handle");

    let OpOutput::Extrude { solid: prism, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square(0.0, 1.0),
            dist: 0.5,
        },
    )
    .expect("a sketch on the meet-judged frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the prism on a straddle datum must be a valid closed b-rep"
    );

    // ★★★ **The honest boundary until 16-3**: this prism's base cap is a plane defined through
    // an implicit point, which the judging table cannot seat yet — a boolean touching it must
    // refuse **by its own name**, not by a name about rational width.
    let block = m.add_cuboid(
        Point3::from_array([-30.0, -30.0, -30.0]),
        Point3::from_array([30.0, 30.0, 30.0]),
    );
    m.rebuild_adjacency();
    let live_before = m.live_solids.clone();
    match nacre_ops::boolean(&mut m, nacre_ops::BoolKind::Common, block, prism) {
        Err(nacre_ops::BoolError::Unsupported {
            reason: nacre_ops::RejectReason::ImplicitPlaneUnsupported,
        }) => {}
        other => panic!("expected the implicit-plane refusal, got {other:?}"),
    }
    assert_eq!(m.live_solids, live_before, "a refusal leaves the live set");

    // Deterministic rebuild: the same construction accepts the same statements into the same
    // arena slots — the meet road is judged once at a fixed rung, so replay holds.
    let mut scratch = Model::new();
    let a2 = scratch.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    let b2 = scratch.add_cuboid(
        Point3::from_array([3.3, 3.3, -1.0]),
        Point3::from_array([7.7, 7.7, 11.0]),
    );
    let OpOutput::Transform { solid: b2 } = apply(
        &mut scratch,
        &Operation::Transform {
            solid: b2,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    scratch.rebuild_adjacency();
    nacre_ops::boolean(&mut scratch, nacre_ops::BoolKind::Cut, a2, b2).expect("cut");
    scratch.rebuild_adjacency();
    let swapped2 = swapped.map(|v| {
        scratch
            .vertex_handle_at(v.index())
            .expect("same construction, same indices")
    });
    let OpOutput::DatumPlane { plane: plane2, .. } = apply(
        &mut scratch,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(swapped2),
        },
    )
    .expect("the same statement is accepted again") else {
        unreachable!()
    };
    assert_eq!(
        plane2.index(),
        plane.index(),
        "same statement, same arena slot"
    );
}

/// ★★★★ **A solid built on a vertex-named datum moves exactly once.**
///
/// Two failures are possible and they point opposite ways, so both are measured:
/// **zero** times, if the mover took the no-node path — the f64 cache would travel while the
/// truth kept pointing at vertices that did not move, leaving the two describing different
/// planes; and **twice**, if the definition were transported *and* a node recorded.
///
/// ★ The translation is far from the origin on purpose. Near zero, "twice as far" and "not at
/// all" both look like zero, and the fixture would pass whatever the code did.
#[test]
fn a_prism_on_a_vertex_named_datum_moves_exactly_once() {
    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let mut m = tilted_prism_with_pocket();
    let vs = three_solvable(&m);
    let OpOutput::DatumPlane { plane, frame } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("datum") else {
        unreachable!()
    };
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(-0.2, -0.2),
                p2(0.2, -0.2),
                p2(0.2, 0.2),
                p2(-0.2, 0.2),
            ])
            .unwrap(),
            dist: 0.5,
        },
    )
    .expect("a prism raised on the vertex-named datum") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    let base_before = nacre_ops::face_plane(&m, base_cap(&m, solid, plane)).expect("planar");
    let shift = [7.0, 11.0, 13.0]; // far from the origin, so 0× and 2× are distinguishable
    let OpOutput::Transform { solid: moved } = apply(
        &mut m,
        &Operation::Transform {
            solid,
            isometry: nacre_scalar::Isometry::translation(
                shift.map(|x| nacre_scalar::Rat::from_decimal(x).unwrap()),
            ),
        },
    )
    .expect("moving a solid built on a vertex-named datum") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    let cap = m
        .shells
        .get(m.solids.get(moved).outer)
        .faces
        .iter()
        .copied()
        .find(|&f| {
            matches!(
                m.surface_truth(m.faces.get(f).surface),
                SurfaceTruth::Plane {
                    points: PlanePoints::Through(_),
                    ..
                }
            )
        })
        .expect("the moved base cap still points at vertices");
    let after = nacre_ops::face_plane(&m, cap).expect("planar");

    // The plane's own origin travelled once along the shift.
    let d = after.origin() - base_before.origin();
    let along = [d.as_array()[0], d.as_array()[1], d.as_array()[2]];
    for (got, want) in along.iter().zip(shift) {
        assert!(
            (got - want).abs() < 1e-9,
            "the base cap moved {along:?}, expected {shift:?} — 0× means the truth stayed \
             behind, 2× means it was transported and re-moved"
        );
    }

    // ★ And the truth still agrees with the cache: the definition names the *same* vertices plus
    // a node, so the plane derived from the definition must contain the cache's own origin.
    let SurfaceTruth::Plane {
        points: PlanePoints::Through(named),
        motion,
    } = m.surface_truth(m.faces.get(cap).surface)
    else {
        unreachable!()
    };
    assert_eq!(*named, {
        let mut s = vs;
        s.sort_by_key(|v| v.index());
        s
    });
    assert!(
        motion.is_some(),
        "a Through plane must take the recorded-node path, never the transporting one"
    );
}

/// The prism's base cap — the face whose surface is the datum handle.
fn base_cap(
    m: &Model,
    solid: nacre_store::Handle<nacre_topo::Solid>,
    plane: nacre_store::Handle<Surface>,
) -> nacre_store::Handle<nacre_topo::Face> {
    m.shells
        .get(m.solids.get(solid).outer)
        .faces
        .iter()
        .copied()
        .find(|&f| m.faces.get(f).surface == plane)
        .expect("the base cap is the datum's own plane")
}

/// ★★ **`Copy` must not refuse it either.** `defs_are_remappable` asks whether every *vertex*
/// definition names surfaces this walk will remap; a `Through` plane is the first reference
/// pointing the other way, at vertices that need not belong to the solid at all. It is not
/// remapped — the node carries the motion — but "not needed" is an argument, and this is the
/// observation.
#[test]
fn a_solid_on_a_vertex_named_datum_can_be_copied() {
    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let mut m = tilted_prism_with_pocket();
    let vs = three_solvable(&m);
    let OpOutput::DatumPlane { frame, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("datum") else {
        unreachable!()
    };
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(-0.2, -0.2),
                p2(0.2, -0.2),
                p2(0.2, 0.2),
                p2(-0.2, 0.2),
            ])
            .unwrap(),
            dist: 0.5,
        },
    )
    .expect("prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let out = apply(&mut m, &Operation::Copy { solid });
    assert!(
        matches!(out, Ok(OpOutput::Copy { .. })),
        "copy refused a solid whose base cap points at vertices: {out:?}"
    );
}

/// ★★★★ **The invariant S2 established, re-scoped for open item 16.** What S2 drained was the
/// *record-less* plane — no points, no truth to derive anything from. A nameless `Through`
/// statement (mixed frames — its exact world coefficients are irrational) is not that: the truth
/// is complete, only a rational description of it does not exist, and rule 6 wrote that
/// population into the design from the start. So the invariant is now: **a plane without a name
/// is exactly a `Through` statement, and everything else keeps its name.** This model holds no
/// nameless planes at all (the datum here is rational-closure), so the sweep also proves the
/// named road did not lose anyone.
#[test]
fn every_plane_still_has_a_name() {
    let mut m = tilted_prism_with_pocket();
    let vs = three_solvable(&m);
    apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("datum");
    // ★ Walked through `surface_handle_at`, which is the only door out of S1's seal — the store
    // itself is private, so this sweep sees exactly what a caller can.
    let mut planes = 0;
    let mut i = 0u32;
    while let Some(h) = m.surface_handle_at(i) {
        i += 1;
        if let SurfaceTruth::Plane { points, .. } = m.surface_truth(h) {
            planes += 1;
            assert!(
                m.surface_name.contains_key(&h) || matches!(points, PlanePoints::Through(_)),
                "a plane with no name that is not a Through statement — \
                 the record-less population S2 drained is back"
            );
        }
    }
    assert!(planes > 0, "the sweep found no planes to check");
}

/// ★★★★★ **A datum through vertices must not collide with a plane it is not.**
///
/// A prism raised on a tilted frame states its far cap **inside that frame**, where it is
/// `w = dist` — canonical name `[0,0,1,−dist]`, which is *letter for letter* the name of the world
/// plane `z = dist`. What keeps those apart is the motion in the interning key
/// (`SurfaceKey = (name, motion)`): the far cap is filed under `(that name, Some(node))`.
///
/// A datum through the far cap's corners solves them **in the same frame**, so it derives the same
/// name — and if it files that under `None`, it is claiming to be a world plane. The model here
/// holds `z = dist` as a box top, so the claim is resolved by interning handing back **the box's
/// face**: a plane somewhere else entirely, which the caller then sketches on.
///
/// The assertion is the harm, not the mechanism: *this datum is not the box top*. That stays true
/// however the cause is later described.
#[test]
fn a_datum_through_frame_local_vertices_is_not_a_world_plane() {
    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let mut m = Model::new();

    // A plane whose normal has `n·n = 3` — not a perfect square, so `exact_frame` declines and the
    // prism is written against a recorded frame node (scalar: "it is `3` for `[1,1,1]`").
    let tilted = SketchPlane::from_origin_normal(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .expect("a tilted plane");
    let frame = datum_frame(&mut m, tilted);
    let dist = 2.0;
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(0.0, 0.0),
                p2(1.0, 0.0),
                p2(1.0, 1.0),
                p2(0.0, 1.0),
            ])
            .unwrap(),
            dist,
        },
    )
    .expect("a prism on the tilted frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    // The world plane `z = dist`, as an ordinary box's top face.
    let boxy = m.add_cuboid(
        Point3::from_array([10.0, 10.0, 0.0]),
        Point3::from_array([11.0, 11.0, dist]),
    );
    m.rebuild_adjacency();
    let box_top = m
        .shells
        .get(m.solids.get(boxy).outer)
        .faces
        .iter()
        .map(|&f| m.faces.get(f).surface)
        .find(|s| {
            m.surface_name
                .get(s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, -dist])
        })
        .expect("the box top names z = dist");

    // The far cap's corners: every one of their carriers is written in the prism's frame, so the
    // three share a motion and the datum is buildable. (Base-cap corners would not: the base cap
    // is the stated plane with no motion, so those triples are mixed-frame and rejected.)
    let far_cap = m
        .shells
        .get(m.solids.get(solid).outer)
        .faces
        .iter()
        .map(|&f| m.faces.get(f).surface)
        .find(|s| {
            matches!(
                m.surface_truth(*s),
                SurfaceTruth::Plane {
                    motion: Some(_),
                    ..
                }
            ) && m
                .surface_name
                .get(s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, -dist])
        })
        .expect("the far cap is `w = dist` in the frame — the premise of this test");
    assert_ne!(
        far_cap, box_top,
        "the far cap and the box top share a name and are kept apart by the motion — \
         if they were already one handle this test would be vacuous"
    );

    let corners: Vec<_> = live_verts(&m)
        .into_iter()
        .filter(|v| match m.vertices.get(*v).def {
            nacre_topo::VertexDef::ThreePlane(tri) => tri.contains(&far_cap),
            _ => false,
        })
        .collect();
    assert!(corners.len() >= 3, "the far cap has corners to name");

    let OpOutput::DatumPlane { plane, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices([corners[0], corners[1], corners[2]]),
        },
    )
    .expect("three far-cap corners share one frame, so this datum is buildable") else {
        unreachable!()
    };

    // ★ The harm, stated as the harm.
    assert_ne!(
        plane, box_top,
        "the datum through the prism's far-cap corners came back as the box's top face — \
         a plane in a different place, which the caller would now sketch on"
    );
    // ★ And it *is* that far cap: the same three points, so the same plane.
    assert_eq!(
        plane, far_cap,
        "a datum through three corners of the far cap is the far cap"
    );

    // ★★ The truth says which frame its points are written in — the disambiguator that keeps the
    // two same-named planes apart.
    let SurfaceTruth::Plane { motion, .. } = m.surface_truth(plane) else {
        unreachable!()
    };
    assert!(
        motion.is_some(),
        "a datum whose vertices live in a frame must record that frame"
    );

    // ★★★ And the f64 cache is a *world* description of that same plane: put the cache's own
    // anchor back through the definition and the residual must vanish. A cache built from frame
    // coordinates would land far off, which is the other half of the defect.
    let sp = nacre_ops::face_plane(&m, {
        *m.shells
            .get(m.solids.get(solid).outer)
            .faces
            .iter()
            .find(|&&f| m.faces.get(f).surface == plane)
            .expect("the far cap is a face of this prism")
    })
    .expect("planar");
    // ★ The statement is **"the named vertices lie on the cached plane"**, not "the cache's
    // origin is one of them" — `face_plane`'s origin is the canonical foot of perpendicular, so
    // comparing positions would pass or fail for reasons that have nothing to do with the defect.
    let n = sp.normal().as_array();
    let o = sp.origin().as_array();
    for (i, v) in corners.iter().take(3).enumerate() {
        let p = m.vertex_point(*v).as_array();
        let d = (0..3).map(|k| n[k] * (p[k] - o[k])).sum::<f64>().abs();
        assert!(
            d < 1e-9,
            "named vertex {i} sits {d} off the plane's own f64 cache — \
             a frame coordinate was used as a world point"
        );
    }
}
