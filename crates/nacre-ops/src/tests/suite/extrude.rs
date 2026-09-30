//! Extrusion: profiles, outward normals, swept prisms and the surfaces a prism reuses.

use super::*;

#[test]
fn square_extrudes_to_a_cube() {
    let m = replay(&[extrude_log_op(square(), 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertex_count(), 8);
    assert_eq!(m.edge_count(), 12);
    assert_eq!(m.face_count(), 6);
    assert_eq!(m.solid_count(), 1);

    let mut got: Vec<[f64; 3]> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .map(|h| (h, m.vertex(h)))
        .map(|(vh, _)| m.vertex_point(vh).as_array())
        .collect();
    let mut want = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let key = |p: &[f64; 3]| (p[0].to_bits(), p[1].to_bits(), p[2].to_bits());
    got.sort_by_key(key);
    want.sort_by_key(key);
    assert_eq!(got, want);
}

#[test]
fn triangle_extrudes_to_a_prism() {
    let tri = Profile2d::polygon(vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(1.0, 1.5)]).unwrap();
    let m = replay(&[extrude_log_op(tri, 3.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertex_count(), 6);
    assert_eq!(m.edge_count(), 9);
    assert_eq!(m.face_count(), 5);
}

#[test]
fn pentagon_extrudes_clean() {
    let m = replay(&[extrude_log_op(regular_ngon(5, 2.0), 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertex_count(), 10);
    assert_eq!(m.face_count(), 7);
}

#[test]
fn concave_l_profile_is_valid() {
    // An L-shape (a reflex vertex) — a simple concave hexagon.
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertex_count(), 12);
    assert_eq!(m.face_count(), 8);
}

/// `FaceInfo::n_out` is documented as the single source of "outward". Two
/// independent sources say which way that is: the ring, which winds CCW about the
/// outward normal, and the b-rep's own `Surface` plus `Orientation`. They must agree
/// on every face of every solid.
///
/// `n_out` reads the second source (the stored orientation), so what this pins is
/// the first: the witness triangle `outer_tri` picks must span the ring's winding.
/// The turn at the first non-collinear corner is the winding **only when that corner
/// is convex**, and the four fixtures below would pass that reading by accident —
/// none starts its cap loop one vertex before a reflex corner. `rotated_l_prism`
/// does, and it is the same solid.
///
/// ★ And a loop's **arcs** wind it too: the three-quarter disk's caps are three points — the
/// centre and the arc's two ends — whose chord triangle turns the other way round from the face,
/// the long arc enclosing the region on the chords' far side. Both caps walk the arc, in
/// opposite directions, in the world frame and on a tilted one.
#[test]
fn outward_normals_agree_with_their_orientation() {
    let mut cube = Model::new();
    let c = cube.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let (ml, sl) = l_prism();
    let (mu, su) = u_prism();
    let (mr, sr) = rotated_l_prism();
    let sector = |tilted: bool| {
        let mut m = Model::new();
        let frame = if tilted {
            pythagorean_frame(&mut m, Point3::from_array([0.0; 3]))
        } else {
            SketchFrame::world(&m, nacre_exact::Axis::Z)
        };
        let profile = stated(vec![
            line(p2(0.0, 0.0), p2(1.0, 0.0)),
            arc_turns(p2(0.0, 0.0), p2(1.0, 0.0), 3),
            line(p2(0.0, -1.0), p2(0.0, 0.0)),
        ])
        .unwrap()
        .remove(0);
        let Ok(OpOutput::Extrude { solid, .. }) = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 1.0,
            },
        ) else {
            panic!("the three-quarter disk extrudes")
        };
        m.rebuild_adjacency();
        (m, solid)
    };
    let (mq, sq) = sector(false);
    let (mt, st) = sector(true);
    for (name, m, s) in [
        ("cube", &cube, c),
        ("l_prism", &ml, sl),
        ("u_prism", &mu, su),
        ("rotated_l_prism", &mr, sr),
        ("three-quarter disk", &mq, sq),
        ("three-quarter disk, tilted", &mt, st),
    ] {
        for row in &collect_planes(m, s).unwrap() {
            // A lateral has no one normal to agree with.
            let crate::planes::FaceRow::Plane(pi) = row else {
                continue;
            };
            let (tri, _) = crate::planes::outer_tri(m, m.face(pi.face.expect("a model face")))
                .expect("a corner");
            let n_out = pi.plane.normal() * f64::from(pi.orient_sign);
            let cos = (tri[1] - tri[0])
                .cross(tri[2] - tri[0])
                .normalize()
                .expect("a widest corner spans area")
                .dot(n_out);
            assert!(
                cos > 0.5,
                "{name}: the witness triangle does not span its face's stated outward"
            );
        }
    }
}

/// `pt3_base_collinear` is exact on the pre-rotation rational bases: three genuinely
/// collinear points stay collinear under rotation (→ skipped), and a real sliver (one point
/// off the line) is never falsely called collinear (→ its crossing is kept, no silent-wrong).
#[test]
fn pt3_base_collinear_exact() {
    use nacre_exact::{Angle, Axis, Rat};
    let ang = Angle::from_deg(Rat::from_int(37)).unwrap();
    let piv = [Rat::from_int(2), Rat::from_int(-1), Rat::from_int(0)];
    let rp = |x: i128, y: i128, z: i128| {
        WitnessPoint::at([Rat::from_int(x), Rat::from_int(y), Rat::from_int(z)]).rotate_about(
            Axis::Z,
            ang,
            piv,
        )
    };
    // (0,0,0), (2,4,6), (1,2,3): all on the line t·(1,2,3) → collinear.
    assert!(pt3_base_collinear(&rp(0, 0, 0), &rp(2, 4, 6), &rp(1, 2, 3)));
    // (1,2,4) is off that line (z), a real nonzero-area triangle → not collinear.
    assert!(!pt3_base_collinear(
        &rp(0, 0, 0),
        &rp(2, 4, 6),
        &rp(1, 2, 4)
    ));
}

/// Explicit sharing: a prism built with a shared base-cap surface reuses that `Surface` handle for
/// its flush cap, with the orientation its caller states — here the one a pad on a `+z` face states
/// (the face's own, reversed), so the materialized outward normal is `−sweep`.
#[test]
fn build_prism_base_cap_reuses_shared_surface() {
    let mut m = Model::new();
    // A face-plane surface with outward normal +z (as a face on the base solid).
    let r = nacre_exact::Rat::from_int;
    // `z = 0` interns onto the seeded world plane, which faces `−z`: `flipped` says the handle
    // faces the other way from this statement, so a face on it facing `+z` is `Reversed`.
    let (sf, flipped) = m.push_plane(
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap(),
        [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
        None,
        nacre_topo::Orientation::Forward,
    );
    let base_pts = [
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    ];
    let (_prism, faces) = build_prism(
        &mut m,
        swept_world(base_pts.to_vec(), Vector3::from_array([0.0, 0.0, 1.0])),
        vec![],
        Vector3::from_array([0.0, 0.0, 1.0]),
        Some((
            sf,
            // A pad's base cap: the face's orientation on the surface, reversed.
            if flipped {
                nacre_topo::Orientation::Forward
            } else {
                nacre_topo::Orientation::Reversed
            },
        )),
        None,
    )
    .unwrap();
    let cap = m.face(faces[0]); // base cap is pushed first
    // Shared handle (was a fresh push before overhaul #3).
    assert_eq!(cap.surface, sf, "base cap reuses the shared surface handle");
    // The stated orientation: materialized outward normal is −sweep (−z).
    let nacre_geom::Surface::Plane(p) = m.surface_cache(cap.surface) else {
        unreachable!()
    };
    let sign = f64::from(cap.orientation.sign());
    let materialized = p.normal() * sign;
    assert!(
        (materialized - Vector3::from_array([0.0, 0.0, -1.0])).norm() < 1e-12,
        "materialized cap normal stays −z, got {materialized:?}"
    );
}

/// ★★★★★ **A plane's record and its motion are one statement.**
///
/// No motion means *"these points speak about the world"*, a motion means *"about the
/// pre-motion frame"* — so a base cap that takes the caller's world triple must carry no
/// motion, and one that falls back to the prism's own ring must take the frame's motion with
/// it. There is nothing else to keep in step: the plane's canonical name is derived from
/// whichever triple is recorded.
///
/// ★★ **Two halves, and they cannot come apart.** Coefficients supplied beside the points and
/// chosen by a **separate** `or_else` would record world coefficients beside the prism's own
/// frame ring whenever the points overflow and the coefficients do not — one plane stated two
/// ways, which an agreement filter cannot catch either, since `c · p` overflows at exactly those
/// widths. With no coefficient parameter the pairing is structural; this pins the choice that is
/// left.
#[test]
fn a_prisms_base_cap_records_the_frame_its_def_names() {
    let r = nacre_exact::Rat::from_int;
    let base_pts = [
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    ];
    let prism = |pts: Option<[[nacre_exact::Rat; 3]; 3]>| {
        let mut m = Model::new();
        let (_prism, faces) = build_prism(
            &mut m,
            swept_world(base_pts.to_vec(), Vector3::from_array([0.0, 0.0, 1.0])),
            vec![],
            Vector3::from_array([0.0, 0.0, 1.0]),
            None,
            pts,
        )
        .unwrap();
        let surf = m.face(faces[0]).surface; // base cap is pushed first
        (m, surf)
    };

    // ★ The caller's triple, when there is one — here a plane nowhere near the prism, so a
    // record that ignored it would be visibly different rather than coincidentally equal.
    let far_pts = [[r(0), r(0), r(3)], [r(1), r(0), r(3)], [r(0), r(1), r(3)]];
    let (m, surf) = prism(Some(far_pts));
    assert_eq!(
        m.surface(surf),
        &nacre_topo::Surface::Plane {
            points: nacre_topo::PlanePoints::Known(far_pts),
            motion: None,
            // The triple spans `+z`; a base cap faces against the sweep.
            sense: nacre_topo::Orientation::Reversed,
        },
        "the caller's triple was not the one recorded (or gained a motion)"
    );
    assert_eq!(
        m.surface_name.get(&surf),
        Some(&nacre_exact::PlaneName::Narrow([r(0), r(0), r(1), r(-3)])),
        "the name was not derived from the triple that was recorded"
    );

    // ★ And with no caller statement, the ring answers — its own exact cap triple, and the
    // name derived from it (`z = 0`, visibly different from the caller's `z = 3` above).
    // There is no f64 fallback, so there is no point-less cap to pin: the negative cases live
    // at the operation as named rejects.
    let (m, surf) = prism(None);
    assert!(
        matches!(m.surface(surf), nacre_topo::Surface::Plane { .. }),
        "the ring's triple is recorded"
    );
    assert_eq!(
        m.surface_name.get(&surf),
        Some(&nacre_exact::PlaneName::Narrow([r(0), r(0), r(1), r(0)])),
        "the name is derived from the ring's own plane"
    );
}

/// **A shared `Surface` handle is one class on its own** — `plane_classes` merges by handle
/// before it asks the judge anything, so two rows on one surface are one class even when their
/// witnesses describe two different planes (which no producer makes — interning is what shared
/// the handle — so a hand-built table is the only way to separate the two questions).
#[test]
fn plane_classes_merge_a_shared_handle_before_judging() {
    let mut m = Model::new();
    let r = nacre_exact::Rat::from_int;
    let shared = m.push_plane_unregistered(
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0])).unwrap(),
        [[r(0); 3], [r(0), r(1), r(0)], [r(0), r(0), r(1)]],
        nacre_topo::Orientation::Forward,
    );
    let fh = m.push_face_unchecked(Face {
        surface: shared,
        outer: Loop { half_edges: vec![] },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    let plane_x0 =
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0])).unwrap();
    let plane_z0 =
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap();
    // Each witness is three NON-collinear points of its own plane, so the judge can tell the two
    // rows apart and the merge below rests on the handle alone.
    let mk = |plane, tri: [Point3; 3]| {
        crate::planes::FaceRow::Plane(FaceInfo {
            // Unmoved and hand-built: nothing to record, and the base frame is unused anyway.
            base_rat: None,
            world: None,
            name: None,
            motion: None,
            surf: shared,
            face: Some(fh),
            plane,
            // Unread: this table only ever reaches `Judge::planes_coplanar`, which decides on the
            // witnesses' definitions (`tri_pt3`) — the rows have no world name.
            orient_sign: 1,
            tri_pt3: tri.map(|p| {
                nacre_judge::WitnessPoint::at_nearest(
                    p.as_array()
                        .map(|x| nacre_exact::Rat::try_from_f64(x).expect("exact")),
                )
            }),
            rotated: false,
        })
    };
    let p = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
    let planes = vec![
        mk(plane_x0, [p(0., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]), // in x = 0
        mk(plane_z0, [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.)]), // in z = 0
    ];
    // The definitions say these really are two different planes…
    assert!(!crate::planes::test_judge(&planes).planes_coplanar(0, 1));
    // …and the shared handle alone makes them one class.
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&planes));
    assert_eq!(canon[0], canon[1], "one handle, one class");
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn prop_regular_ngon_on_xy_is_clean(
        n in 3usize..8,
        r in 0.5f64..10.0,
        dist in 0.1f64..10.0,
    ) {
        let m = replay(&[extrude_log_op(regular_ngon(n, r), dist)]).unwrap();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        prop_assert_eq!(m.vertex_count(), 2 * n);
        prop_assert_eq!(m.face_count(), n + 2);
    }

    #[test]
    fn prop_ngon_on_arbitrary_plane_is_clean(
        n in 3usize..8,
        nx in -1.0f64..1.0,
        ny in -1.0f64..1.0,
        nz in -1.0f64..1.0,
        dist in 0.1f64..10.0,
    ) {
        let normal = Vector3::from_array([nx, ny, nz]);
        prop_assume!(normal.norm() > 0.1); // skip near-zero normals
        let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
        let __plane = plane;
        let mut __scratch205 = Model::new();
        let __g204 = datum_frame(&mut __scratch205, __plane);
        let m = replay(&[
            Operation::DatumPlane { def: DatumDef::Stated(__plane) },
            Operation::Extrude {
                frame: __g204,
            profile: regular_ngon(n, 2.0),
            dist,
        }])
        .unwrap();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        prop_assert_eq!(m.face_count(), n + 2);
    }

    /// A blind pocket on a randomly-slanted face: the arrangement must give a valid solid of the
    /// right volume or reject honestly — **never panic**. Drives general (non-axis) plane normals
    /// through the dir-sign guard and the `angular_order`/`turn_at` consumers that read its zeros.
    ///
    /// Locks a `D = 0` panic on general normals: a triple naming one geometric plane through
    /// two coincident faces (symmetric normals like `(1,1,1)` clear it, `(0.446, 0.737, 0.990)`
    /// does not). Naming every plane by its class makes those two faces one index, so the
    /// degenerate triple cannot form.
    #[test]
    fn pocket_on_a_random_slanted_face_is_valid_or_rejects(
        nx in -1.0f64..1.0,
        ny in -1.0f64..1.0,
        nz in 0.2f64..1.0, // keep the normal clear of the sketch's degenerate zero
    ) {
        let normal = Vector3::from_array([nx, ny, nz]);
        prop_assume!(normal.norm() > 0.3);
        let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
        let mut m = Model::new();
        let big = Profile2d::polygon(vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)]).unwrap();
        let frame = datum_frame(&mut m, plane);
        let OpOutput::Extrude { faces, .. } =
            apply(&mut m, &Operation::Extrude { frame, profile: big, dist: 2.0 }).unwrap()
        else { unreachable!() };
        match apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)) {
            Ok(OpOutput::PocketOnFace { solid, .. }) => {
                m.rebuild_adjacency();
                prop_assert!(nacre_validate::validate(&m).is_empty());
                let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
                prop_assert!((vol - (8.0 - 0.16 * 0.5)).abs() <= 1e-9 * 8.0, "volume {}", vol);
            }
            Ok(_) => prop_assert!(false, "unexpected op output"),
            Err(_) => {} // an honest reject is acceptable; a panic is not (and would fail the test)
        }
    }

    /// A random boss on a random box stays a valid b-rep (any interior
    /// profile, any positive height).
    #[test]
    fn prop_pad_stays_valid(
        sx in 0.5f64..5.0,
        sy in 0.5f64..5.0,
        sz in 0.5f64..5.0,
        h in 0.05f64..0.15,
        dist in 0.1f64..5.0,
    ) {
        let rect = Profile2d::polygon(vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)]).unwrap();
        let mut m = Model::new();
        let __w5 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
            frame: __w5,
            profile: rect,
            dist: sz,
        }).unwrap() else { unreachable!() };
        let hole = Profile2d::polygon(vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)]).unwrap();
        apply(&mut m, &Operation::PadOnFace { face: faces[1], profile: hole, dist }).unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A random blind pocket on a random box stays valid. `dist ≤ 0.8 < sz`
    /// keeps the pocket from punching through the box (height `sz ≥ 1`).
    #[test]
    fn prop_pocket_stays_valid(
        sx in 0.5f64..5.0,
        sy in 0.5f64..5.0,
        sz in 1.0f64..5.0,
        h in 0.05f64..0.15,
        dist in 0.1f64..0.8,
    ) {
        let rect = Profile2d::polygon(vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)]).unwrap();
        let mut m = Model::new();
        let __w4 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
            frame: __w4,
            profile: rect,
            dist: sz,
        }).unwrap() else { unreachable!() };
        let hole = Profile2d::polygon(vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)]).unwrap();
        apply(&mut m, &Operation::PocketOnFace { face: faces[1], profile: hole, dist }).unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
    }
}
