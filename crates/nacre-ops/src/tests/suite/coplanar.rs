//! Coplanar contact: coincident merges, boss fuses, plane classes and `unify_coplanar_faces`.

use super::*;

#[test]
fn fuse_a_boss_onto_a_non_convex_solid() {
    // A contained boss on the top of an L-prism (non-convex kept `a`); the contained-coplanar
    // Fuse admits it. Volume 3 (L) + 0.4²·0.5 = 3.08.
    let (mut m, l) = l_prism(); // L footprint area 3, height 1
    let boss = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 1.0]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, l, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 3.08).abs() < 1e-12, "volume {vol}");
    // The boss top cap sits on the z = 1.5 plane, its outward normal +z.
    assert!(has_face_on_plane(
        &m,
        r,
        Point3::from_array([0.5, 0.5, 1.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    ));
}

#[test]
fn fuse_a_non_convex_profile_boss() {
    // An L-shaped boss (non-convex cutter `b`) on a cube top; the contained-coplanar Fuse
    // carries the L footprint as a hole.
    // Volume 1 (cube) + 0.12 (L area) · 0.4 = 1.048.
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let l_base: Vec<Point3> = [
        [0.3, 0.3],
        [0.7, 0.3],
        [0.7, 0.5],
        [0.5, 0.5],
        [0.5, 0.7],
        [0.3, 0.7],
    ]
    .iter()
    .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
    .collect();
    let (boss, _) = build_prism(
        &mut m,
        swept_world(l_base, Vector3::from_array([0.0, 0.0, 0.4])),
        vec![],
        Vector3::from_array([0.0, 0.0, 1.0]),
        None,
        None,
    )
    .unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, cube, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.048).abs() < 1e-12, "volume {vol}");
    // The L boss top cap sits on the z = 1.4 plane, its outward normal +z.
    assert!(has_face_on_plane(
        &m,
        r,
        Point3::from_array([0.4, 0.4, 1.4]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    ));
}

#[test]
fn cut_by_an_overhanging_boss_carrying_a_pin_owes_a_notch() {
    // Same seating, but the tool carries a pin reaching below the contact plane, so the plane
    // no longer separates the solids and the cut owes a real notch (1 − 0.2·0.2·0.5 = 0.98).
    //
    // This used to be an honest reject: the tool's z=1 cap is an annulus-like face whose
    // *inner* edge rides the pin's walls, and reading its occupancy off the ring's flank put
    // the material on the wrong side, so the class would not label. With the side read from
    // the ring's travel instead (`arrangement::run_body_above`), the notch comes out at the
    // hand-computed volume with a clean model.
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let pin = m.add_cuboid(
        Point3::from_array([0.55, 0.55, 0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    m.rebuild_adjacency();
    let tool = boolean_one(&mut m, BoolKind::Fuse, block, pin).unwrap();
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, tool).unwrap().volume - 1.02).abs() < 1e-12,
        "the pinned tool itself"
    );
    let notched = boolean_one(&mut m, BoolKind::Cut, base, tool).expect("the notch is buildable");
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, notched).unwrap().volume - 0.98).abs() < 1e-12,
        "the notch the tool owes: {}",
        nacre_props::mass_props(&m, notched).unwrap().volume
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "a notched result is still a clean model"
    );
}

// Plane-class canonicalization — coplanar walls of the two operands fold into one line.
#[test]
fn plane_classes_merge_a_shared_wall() {
    // Two unit cubes side by side share the plane x=1 (a's +x wall, b's -x wall — the same
    // plane, opposite normals). `plane_classes` must merge those two into one line class and
    // keep the far walls (a's x=0, b's x=2) distinct.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let planes_a = collect_planes(&m, a).unwrap();
    let na = planes_a.len();
    let mut planes = planes_a;
    planes.extend(collect_planes(&m, b).unwrap());
    // Find a plane by outward-normal x-sign and its x coordinate, within an index range.
    let find = |rng: std::ops::Range<usize>, nx: f64, x: f64| -> usize {
        rng.clone()
            .find(|&i| {
                let n = planes[i].plane().n_out.as_array();
                n[0] * nx > 0.5 && (planes[i].plane().tri[0].as_array()[0] - x).abs() < 1e-9
            })
            .expect("plane")
    };
    let a_xp = find(0..na, 1.0, 1.0); // a's +x wall at x=1
    let b_xm = find(na..planes.len(), -1.0, 1.0); // b's -x wall at x=1
    let a_xm = find(0..na, -1.0, 0.0); // a's -x wall at x=0
    let b_xp = find(na..planes.len(), 1.0, 2.0); // b's +x wall at x=2
    let canon = plane_classes(&crate::planes::test_judge(&planes));
    assert_eq!(canon[a_xp], canon[b_xm], "shared x=1 wall is one class");
    assert_ne!(canon[a_xm], canon[b_xp], "far walls stay distinct");
    assert_ne!(canon[a_xp], canon[a_xm], "x=1 and x=0 are different lines");
    // Every b face but its +x wall is coplanar with an a face, so 12 planes fold to 7 classes.
    let distinct: std::collections::HashSet<usize> = canon.iter().copied().collect();
    assert_eq!(distinct.len(), na + 1, "only b's far wall is a new class");
    // The class root is the smallest index in the class (deterministic canon).
    assert_eq!(canon[a_xp], a_xp.min(b_xm));
}

/// ★★★★ **Two faces of one judged surface are one class.** A nameless plane has no
/// name to intern classes by, so its class merging rests on the probes: both faces carry the
/// *same statement* → the same frame probes → identical chains → `shared_base` cancels the
/// motion and the exact predicate answers a **proved zero**. This is the argument turned into
/// a run: a straddle-datum prism is severed through its cap, leaving two faces on the one
/// judged surface, and `plane_classes` must fold them into one line.
#[test]
fn two_faces_of_one_judged_surface_are_one_class() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([3.3, 3.3, -1.0]),
        Point3::from_array([7.7, 7.7, 11.0]),
    );
    m.rebuild_adjacency();
    let OpOutput::Transform { solid: b } = crate::apply(
        &mut m,
        &Operation::Transform {
            solid: b,
            isometry: nacre_exact::Isometry::rotation(nacre_exact::Rotation {
                axis: nacre_exact::Axis::Z,
                pivot: [nacre_exact::Rat::from_int(0); 3],
                angle: nacre_exact::Angle::from_deg(nacre_exact::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
    m.rebuild_adjacency();

    // A straddling vertex of the cut, plus two pure corners that share its plane sanely.
    let mut straddle = None;
    let mut pure = Vec::new();
    for &s in &cut {
        for &fh in &m.shell(m.solid(s).outer).faces {
            for lp in std::iter::once(&m.face(fh).outer).chain(m.face(fh).inner.iter()) {
                for &he in &lp.half_edges {
                    let vh = m.he_start(he);
                    let nacre_topo::Vertex::ThreePlane(_) = *m.vertex(vh) else {
                        continue;
                    };
                    // ★ The kernel says which corners it can place in one frame; comparing the
                    // three carriers' motions here would be `vertex_meet`'s rule written twice,
                    // and since the invariant-plane restatement that copy calls a turned solid's
                    // corner straddling when the chain-fixes licence places it.
                    if m.vertex_meet(vh).is_none() {
                        straddle.get_or_insert(vh);
                    } else if !pure.contains(&vh) {
                        pure.push(vh);
                    }
                }
            }
        }
    }
    let straddle = straddle.expect("the cut leaves straddling corners");
    let vs = [straddle, pure[0], pure[1]];
    let OpOutput::DatumPlane { plane, frame } = crate::apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("the straddle datum") else {
        unreachable!()
    };
    assert!(!m.surface_name.contains_key(&plane), "nameless");
    let square = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([1.0, 0.0]),
        Point2::from_array([1.0, 1.0]),
        Point2::from_array([0.0, 1.0]),
    ])
    .unwrap();
    let OpOutput::Extrude { solid: prism, .. } = crate::apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square,
            dist: 0.5,
        },
    )
    .expect("prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    // Sever the prism through the middle with a thin box crossing the whole height: two
    // pieces, each with a base-cap face on the SAME judged surface.
    let (o, u, _v, w) = crate::rotated_vertex::frame_world_basis(
        &m,
        plane,
        &nacre_topo::FramePlacement::Canonical,
        frame.flip(),
    )
    .expect("the judged frame realizes");
    // A knife centred over the prism: origin + 0.5·û ± thickness, spanning v and w amply.
    let centre = Point3::from_array([
        o[0] + 0.5 * u[0] + 0.5 * w[0] * 0.0,
        o[1] + 0.5 * u[1],
        o[2] + 0.5 * u[2],
    ]);
    let knife = m.add_cuboid(
        centre + Vector3::from_array([-0.1, -5.0, -5.0]),
        centre + Vector3::from_array([0.1, 5.0, 5.0]),
    );
    m.rebuild_adjacency();
    let pieces = crate::boolean(&mut m, BoolKind::Cut, prism, knife).expect("sever");
    m.rebuild_adjacency();
    assert!(pieces.len() >= 2, "the knife must sever the prism");

    // The judgment table over the two pieces: their base-cap faces share the judged surface
    // and must fold into one class.
    let planes_a = collect_planes(&m, pieces[0]).unwrap();
    let na = planes_a.len();
    let mut planes = planes_a;
    planes.extend(collect_planes(&m, pieces[1]).unwrap());
    let on_datum: Vec<usize> = (0..planes.len())
        .filter(|&i| planes[i].surf() == plane)
        .collect();
    assert!(
        on_datum.len() >= 2
            && on_datum.iter().any(|&i| i < na)
            && on_datum.iter().any(|&i| i >= na),
        "both pieces must carry a face on the judged surface"
    );
    let canon = plane_classes(&crate::planes::test_judge(&planes));
    let first = canon[on_datum[0]];
    for &i in &on_datum[1..] {
        assert_eq!(
            canon[i], first,
            "two faces of one judged surface must be one class — same statement, same probes"
        );
    }
}

#[test]
fn overhang_boss_with_a_non_convex_footprint() {
    // An L-shaped (non-convex) boss footprint overhanging a cube edge, built exactly. Volume =
    // cube 1.0 + L-prism (area 0.9·0.2 + 0.3·0.2 = 0.24) · height 0.4 = 1.096.
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let l_base: Vec<Point3> = [
        [0.3, 0.3],
        [1.2, 0.3],
        [1.2, 0.5],
        [0.6, 0.5],
        [0.6, 0.7],
        [0.3, 0.7],
    ]
    .iter()
    .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
    .collect();
    let (l_tool, _) = build_prism(
        &mut m,
        swept_world(l_base, Vector3::from_array([0.0, 0.0, 0.4])),
        vec![],
        Vector3::from_array([0.0, 0.0, 1.0]),
        None,
        None,
    )
    .unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, cube, l_tool).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.096).abs() < 1e-12, "volume {vol}");
}

// ---- `unify_coplanar_faces`: the general (interface-free) coplanar merge, on hand-built
// `LocalFace` lists the coincident goldens above never reach — chains, opposite normals,
// holes, seam edges, and the asymmetric T-junction the global dissolve exists to prevent. ----

/// An axis-aligned `FaceInfo` at `d` along its normal, with a **non-degenerate `tri`** whose
/// right-hand normal is `n_out`. The merge reads more than `n_out` now — `loop_winding` and
/// `point_in_ring` name their arguments by plane and evaluate exact predicates on `tri` — so a
/// dummy triangle would make those answers meaningless.
fn mk_axis_plane(m: &mut Model, axis: usize, d: f64, positive: bool) -> WorkingPlane {
    let mut n = [0.0; 3];
    n[axis] = if positive { 1.0 } else { -1.0 };
    let normal = Vector3::from_array(n);
    let mut at = [0.0; 3];
    at[axis] = d;
    let origin = Point3::from_array(at);
    let plane = Plane::from_point_normal(origin, normal).unwrap();
    // Unregistered on purpose: these fixtures push one geometric plane as *two* handles
    // (`positive` both ways), which interning would collapse. The truth (the same tri
    // computed below, lifted) is still stated — nothing point-less enters the arena.
    let lift = |p: Point3| {
        p.as_array()
            .map(|x| nacre_exact::Rat::from_decimal(x).unwrap())
    };
    let (ti, tj) = ((axis + 1) % 3, (axis + 2) % 3);
    let (ti, tj) = if positive { (ti, tj) } else { (tj, ti) };
    let stepr = |k: usize| {
        let mut q = at;
        q[k] += 1.0;
        Point3::from_array(q)
    };
    let surf = m.push_plane_unregistered(
        plane,
        [lift(origin), lift(stepr(ti)), lift(stepr(tj))],
        nacre_topo::Orientation::Forward,
    );
    let face = m.push_face_unchecked(Face {
        surface: surf,
        outer: Loop { half_edges: vec![] },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    // Two in-plane directions whose cross product is `+normal`, so `tri` winds outward.
    let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
    let (i, j) = if positive { (i, j) } else { (j, i) };
    let step = |k: usize| {
        let mut q = at;
        q[k] += 1.0;
        Point3::from_array(q)
    };
    let _ = face;
    let tri = [origin, step(i), step(j)];
    WorkingPlane {
        // A hand-built table has no recorded coefficients; the composed-rotation route
        // declines and the fixture takes the same escalating path it always did.
        base_rat: None,
        world_rat: None,
        name_ints: None,
        base: crate::planes::BaseFrame::none(),
        surf,
        plane,
        tri,
        tri_pt3: tri.map(|p| {
            nacre_judge::WitnessPoint::at_nearest(
                p.as_array()
                    .map(|x| nacre_exact::Rat::try_from_f64(x).expect("exact")),
            )
        }),
        rotated: false,
        frame_sign: 1, // `plane` is built from `normal`, so the two agree
        exact_coeffs: WorkingPlane::reconcile(&plane, tri, false).0,
        exact_normal: WorkingPlane::reconcile(&plane, tri, false).1,
    }
}

#[test]
fn unify_merges_a_coplanar_chain() {
    // Three unit squares on z=0 (+z), tiled in x, each sharing a vertical edge with the
    // next. One plane class, one normal ⇒ all fuse into a single face; the four
    // straight-angle mid-edge vertices dissolve, leaving one 4-corner rectangle.
    //
    // Named the way the arrangement names things — every vertex is the meeting of three
    // planes — because the straight-angle test reads those triples.
    let mut m = Model::new();
    let p = vec![
        mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0, the shared class
        mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
        mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
        mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
        mk_axis_plane(&mut m, 0, 3.0, true),  // 4: x=3
        mk_axis_plane(&mut m, 1, 0.0, false), // 5: y=0
        mk_axis_plane(&mut m, 1, 1.0, true),  // 6: y=1
    ];
    let _canon: Vec<usize> = (0..p.len()).collect();
    let v = |x: usize, y: usize| NodeId::three_planes(Canon3::three([0, x, y])); // class, x-plane, y-plane
    let (c00, c10, c20, c30) = (v(1, 5), v(2, 5), v(3, 5), v(4, 5));
    let (c01, c11, c21, c31) = (v(1, 6), v(2, 6), v(3, 6), v(4, 6));
    let faces = vec![
        face(0, vec![c00, c10, c11, c01], vec![]),
        face(0, vec![c10, c20, c21, c11], vec![]),
        face(0, vec![c20, c30, c31, c21], vec![]),
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p), &[]).unwrap();
    assert_eq!(out.len(), 1, "three coplanar faces fuse into one");
    let l = &out[0].outer.expect_ring();
    assert_eq!(l.len(), 4, "straight-angle mid vertices dissolved: {l:?}");
    for c in [c00, c30, c31, c01] {
        assert!(l.contains(&c), "corner kept");
    }
    for c in [c10, c20, c11, c21] {
        assert!(!l.contains(&c), "mid vertex dropped");
    }
}

#[test]
fn an_overhang_fuse_keeps_the_two_z1_caps_separate() {
    // An overhanging boss splits `z = 1` between two coplanar faces with **opposite** outward
    // normals — the base's exposed top (`+z`) and the boss underside (`-z`). They must not be
    // fused into one face: their `flip` differs, so `unify`'s `(plane_idx, flip)` group key
    // keeps them apart.
    //
    // The invariant is exercised here on the production path, in the default `cargo test` run:
    // `overhang_fuse_then_cut_matches_occt`
    // proves it against OCCT but is `#[ignore]`, so this hand-computed volume is the non-ignored
    // guard. A wrong merge collapses the topology — the volume shifts or `validate` speaks.
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!(
        (vol - 1.5).abs() < 1e-12,
        "base 1 + boss 0.5, no overlap: {vol}"
    );
}

#[test]
fn a_hole_filled_by_two_faces_still_merges() {
    // What matters is that the merge is not special-cased to "a hole filled by exactly one
    // neighbour". A [0,3]² face with a [1,2]² hole, and that hole filled by **two** pieces split
    // at x=1.5: every ring edge between them is carried in both directions, so erasing interior
    // boundary leaves only the outer square — one face, no hole, whatever the filling is cut
    // into. This is the case that separates a general rule from a bespoke one.
    let mut m = Model::new();
    let p = vec![
        mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
        mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
        mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
        mk_axis_plane(&mut m, 0, 1.5, true),  // 3: x=1.5, where the filling is split
        mk_axis_plane(&mut m, 0, 2.0, true),  // 4: x=2
        mk_axis_plane(&mut m, 0, 3.0, true),  // 5: x=3
        mk_axis_plane(&mut m, 1, 0.0, false), // 6: y=0
        mk_axis_plane(&mut m, 1, 1.0, true),  // 7: y=1
        mk_axis_plane(&mut m, 1, 2.0, true),  // 8: y=2
        mk_axis_plane(&mut m, 1, 3.0, true),  // 9: y=3
    ];
    let _canon: Vec<usize> = (0..p.len()).collect();
    let v = |x: usize, y: usize| NodeId::three_planes(Canon3::three([0, x, y]));
    let (o00, o30, o33, o03) = (v(1, 6), v(5, 6), v(5, 9), v(1, 9));
    let (h11, h12, h22, h21) = (v(2, 7), v(2, 8), v(4, 8), v(4, 7));
    let (m12, m11) = (v(3, 8), v(3, 7)); // the split points on the hole's top and bottom
    let faces = vec![
        // Outer square with the hole, wound the way `emit_faces` states it: outer CCW, hole CW.
        face(
            0,
            vec![o00, o30, o33, o03],
            vec![vec![h11, h12, m12, h22, h21, m11]],
        ),
        face(0, vec![h11, m11, m12, h12], vec![]), // left filler
        face(0, vec![m11, h21, h22, m12], vec![]), // right filler
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p), &[]).unwrap();
    assert_eq!(out.len(), 1, "the hole is filled, so one face remains");
    assert!(out[0].inner.is_empty(), "and it has no hole left");
    assert_eq!(out[0].outer.expect_ring().len(), 4, "just the outer square");
    for c in [o00, o30, o33, o03] {
        assert!(out[0].outer.expect_ring().contains(&c), "outer corner kept");
    }
}

#[test]
fn unify_keeps_a_vertex_that_is_a_corner_elsewhere() {
    // F0,F1 on z=0 merge; their shared-edge endpoint (1,0,0) is a straight angle on the
    // merged face but a real corner on a perpendicular face G (plane y=0). Global degree 3
    // ⇒ it is NOT dissolved — a per-face local rule would have, opening a T-junction. Its
    // twin (1,1,0), on the merged face only, IS dissolved.
    let mut m = Model::new();
    let p = vec![
        mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
        mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
        mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
        mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
        mk_axis_plane(&mut m, 1, 0.0, false), // 4: y=0, the perpendicular face's plane
        mk_axis_plane(&mut m, 1, 1.0, true),  // 5: y=1
        mk_axis_plane(&mut m, 2, 1.0, true),  // 6: z=1
    ];
    let _canon: Vec<usize> = (0..p.len()).collect();
    let (v000, v100, v200) = (
        NodeId::three_planes(Canon3::three([0, 1, 4])),
        NodeId::three_planes(Canon3::three([0, 2, 4])),
        NodeId::three_planes(Canon3::three([0, 3, 4])),
    );
    let (v010, v110, v210) = (
        NodeId::three_planes(Canon3::three([0, 1, 5])),
        NodeId::three_planes(Canon3::three([0, 2, 5])),
        NodeId::three_planes(Canon3::three([0, 3, 5])),
    );
    let (v101, v201) = (
        NodeId::three_planes(Canon3::three([2, 4, 6])),
        NodeId::three_planes(Canon3::three([3, 4, 6])),
    );
    let faces = vec![
        face(0, vec![v000, v100, v110, v010], vec![]),
        face(0, vec![v100, v200, v210, v110], vec![]),
        face(4, vec![v200, v100, v101, v201], vec![]), // perpendicular, not coplanar
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p), &[]).unwrap();
    assert_eq!(out.len(), 2, "z=0 pair merges; G stays");
    let merged = out.iter().find(|lf| lf.surf.plane() == 0).unwrap();
    assert!(
        merged.outer.expect_ring().contains(&v100),
        "corner-elsewhere vertex kept (no T-junction)"
    );
    assert!(
        !merged.outer.expect_ring().contains(&v110),
        "pure straight-angle vertex dropped"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Diagonal corner overlaps (clean seam): fuse/cut volumes match the
    /// independent AABB formula (not nacre's own common).
    #[test]
    fn fuse_cut_diagonal_boxes_match_aabb(
        amin in prop::array::uniform3(-3.0f64..3.0),
        aext in prop::array::uniform3(1.0f64..3.0),
        t in prop::array::uniform3(0.15f64..0.6),
        s in prop::array::uniform3(0.3f64..2.0),
    ) {
        let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
        let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
        let bmax: [f64; 3] = std::array::from_fn(|i| amax[i] + s[i]);
        let ov: f64 = (0..3).map(|i| amax[i] - bmin[i]).product();
        let va: f64 = aext.iter().product();
        let vb: f64 = (0..3).map(|i| bmax[i] - bmin[i]).product();

        let build = || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
            let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
            (m, a, b)
        };

        let (mut m1, a1, b1) = build();
        let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1);
        prop_assume!(rf.is_ok()); // skip rare coplanar/degenerate configs
        let rf = rf.unwrap();
        m1.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m1).is_empty());
        let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
        prop_assert!((vf - (va + vb - ov)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

        let (mut m2, a2, b2) = build();
        let rc = boolean_one(&mut m2, BoolKind::Cut, a2, b2).unwrap();
        m2.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m2).is_empty());
        let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;
        prop_assert!((vc - (va - ov)).abs() <= 1e-9 * va, "cut {vc}");
    }

    /// Matched-footprint stacked boxes (coincident z-interface): fuse volume
    /// is the sum, cut is A, common is empty.
    #[test]
    fn stacked_boxes_merge_volumes(
        x0 in -3.0f64..3.0,
        y0 in -3.0f64..3.0,
        dx in 0.5f64..3.0,
        dy in 0.5f64..3.0,
        z0 in -3.0f64..3.0,
        h1 in 0.5f64..3.0,
        h2 in 0.5f64..3.0,
    ) {
        let (x1, y1) = (x0 + dx, y0 + dy);
        let (zm, z1) = (z0 + h1, z0 + h1 + h2);
        let build = || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([x0, y0, z0]), Point3::from_array([x1, y1, zm]));
            let b = m.add_cuboid(Point3::from_array([x0, y0, zm]), Point3::from_array([x1, y1, z1]));
            (m, a, b)
        };
        let (va, vb) = (dx * dy * h1, dx * dy * h2);

        let (mut m1, a1, b1) = build();
        // Not `prop_assume!`: a matched-footprint stack is squarely in coverage whatever the
        // dimensions are, so a reject here is a defect, not an uninteresting sample. Assuming
        // it away is how this property went on passing while the kernel aborted on 2% of the
        // space and rejected 95% of it (family #3).
        let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1)
            .expect("stacked boxes fuse at any dimensions");
        m1.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m1).is_empty());
        let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
        prop_assert!((vf - (va + vb)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

        let (mut m2, a2, b2) = build();
        // The stack shares only its interface plane ⇒ no volume in common, at any dimensions.
        prop_assert!(boolean(&mut m2, BoolKind::Common, a2, b2).unwrap().is_empty());

        let (mut m3, a3, b3) = build();
        let rc = boolean_one(&mut m3, BoolKind::Cut, a3, b3).unwrap();
        let vc = nacre_props::mass_props(&m3, rc).unwrap().volume;
        prop_assert!((vc - va).abs() <= 1e-9 * va, "cut {vc}");
    }
}
