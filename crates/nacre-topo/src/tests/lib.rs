use super::*;
use nacre_math::Vector3;
use proptest::prelude::*;

/// ★★★ **The defect the rational coefficients exist to remove, pinned from both sides.**
///
/// Two boxes meet on the plane `x = 3` with faces of different size. `Plane` stores an
/// un-normalized normal whose length follows that size, so `d = −raw·origin` is a differently
/// rounded product on each side and the two coefficient vectors are **not exactly
/// proportional** — the f64 test says "different planes" about one plane. Measured across the
/// census, 18 pairs are merged only because a second test looks at the faces' coordinates
/// instead.
///
/// The rational coefficients are built from the corners the caller wrote and canonicalized, so
/// they have no scale to disagree about and come out **equal**.
///
/// Both halves are load-bearing. If the first assertion ever fails the f64 defect was fixed
/// somewhere else and this test should be re-read, not deleted; if the second fails the
/// rational path stopped reaching these surfaces.
#[test]
fn two_faces_of_one_plane_disagree_in_f64_and_agree_in_the_rationals() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 2.2, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([3.0, 0.0, 0.0]),
        Point3::from_array([5.0, 13.2, 1.0]),
    );
    // Each solid's face on x = 3: a's outward +X, b's outward −X.
    let face_on_x3 = |s: Handle<Solid>, want_x: f64| -> Handle<Surface> {
        let shell = m.solid(s).outer;
        *m.shell(shell)
            .faces
            .iter()
            .map(|&fh| &m.face(fh).surface)
            .find(|&&sh| match m.surface_cache(sh) {
                nacre_geom::Surface::Plane(p) => {
                    let [a, b, c, d] = p.coefficients();
                    b == 0.0 && c == 0.0 && a != 0.0 && (-d / a - want_x).abs() < 1e-12
                }
                nacre_geom::Surface::Cylinder(_) => false,
            })
            .expect("a face on x = 3")
    };
    let (sa, sb) = (face_on_x3(a, 3.0), face_on_x3(b, 3.0));

    // ★ **They are one handle now** — that is what the rational coefficients bought.
    assert_eq!(sa, sb, "one plane, one surface");
    assert!(m.surface_name.contains_key(&sa), "and it is recorded");

    // The f64 defect that made this necessary, shown on the planes themselves rather than
    // through the model, since the model no longer holds two of them. `Plane` keeps an
    // un-normalized normal whose length follows the face's size, so `d = −raw·origin` is a
    // differently rounded product on each side and the two vectors are not exactly
    // proportional — the f64 test says "different planes" about one plane.
    let wall = |dy: f64| {
        Plane::through_points(
            Point3::from_array([3.0, 0.0, 0.0]),
            Point3::from_array([3.0, dy, 0.0]),
            Point3::from_array([3.0, 0.0, 1.0]),
        )
        .expect("non-degenerate")
    };
    let (pa, pb) = (wall(2.2), wall(13.2));
    assert!(
        !nacre_geom::intersect::planes_coplanar(&pa, &pb),
        "f64 coefficients of one plane at two face sizes: {:?} vs {:?}",
        pa.coefficients(),
        pb.coefficients()
    );
    let rat = |dy: f64| {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).expect("decimal");
        nacre_exact::plane_through_points(
            [r(3.0), r(0.0), r(0.0)],
            [r(3.0), r(dy), r(0.0)],
            [r(3.0), r(0.0), r(1.0)],
        )
    };
    assert_eq!(
        rat(2.2),
        rat(13.2),
        "the rationals have no scale to disagree about"
    );
    assert!(rat(2.2).is_some());
}

/// ★★★ **The third door: a corner whose carriers came by different roads still solves.**
///
/// Built synthetically — three planes, two of them carrying **different** translation chains —
/// because waiting for a boolean to produce the shape would measure "what passed" rather than
/// "what was fixed". The two older doors want one frame (all-world, or one shared chain) and
/// neither holds here; the world road transports each carrier's *name* and meets the three in
/// the world, so the answer carries no leaf and the caller replays nothing.
///
/// Negative controls, both load-bearing: a **rotated** carrier has no rational world name and
/// still declines, and a corner the shared-chain door can answer keeps answering **through
/// that door** (`leaf = Some`) — which is what says the new road is a fallback and not a
/// re-spelling of what already worked.
#[test]
fn a_corner_of_two_translation_chains_solves_in_the_world() {
    use nacre_exact::{Angle, Axis, Rat};
    let mut m = Model::new();
    let r = Rat::from_int;
    // Three axis planes through the origin, pushed as their own statements.
    let plane = |m: &mut Model, n: [f64; 3], pts: [[i128; 3]; 3]| {
        m.push_plane(
            nacre_geom::Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array(n),
            )
            .expect("unit normal"),
            pts.map(|p| p.map(r)),
            None,
        )
        .0
    };
    let px = plane(&mut m, [1.0, 0.0, 0.0], [[0, 0, 0], [0, 1, 0], [0, 0, 1]]);
    let py = plane(&mut m, [0.0, 1.0, 0.0], [[0, 0, 0], [1, 0, 0], [0, 0, 1]]);
    let pz = plane(&mut m, [0.0, 0.0, 1.0], [[0, 0, 0], [1, 0, 0], [0, 1, 0]]);
    // Two different translation chains, and a third carrier that stays in the world.
    let t =
        |m: &mut Model, o: [i128; 3]| m.push_motion(Motion::Translate { offset: o.map(r) }, None);
    let t1 = t(&mut m, [2, 0, 0]);
    let t2 = t(&mut m, [0, 3, 0]);
    // ★ The cache is the **realized world** surface and the truth (points + leaf) stays in the
    // pre-motion frame — that is what `transform` writes, so the fixture writes it too. Giving
    // a moved plane its pre-motion cache would model a state the kernel never builds.
    let moved = |m: &mut Model,
                 src: Handle<Surface>,
                 leaf: Handle<MotionNode>,
                 n: [f64; 3],
                 at: [f64; 3]| {
        let Surface::Plane {
            points: PlanePoints::Known(pts),
            ..
        } = m.surface(src).clone()
        else {
            unreachable!("an axis plane states points")
        };
        let world = nacre_geom::Plane::from_point_normal(
            Point3::from_array(at),
            nacre_math::Vector3::from_array(n),
        )
        .expect("unit normal");
        m.push_plane(world, pts, Some(leaf)).0
    };
    // x = 0 moved by (2,0,0) → the world plane x = 2; y = 0 moved by (0,3,0) → y = 3.
    let mx = moved(&mut m, px, t1, [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
    let my = moved(&mut m, py, t2, [0.0, 1.0, 0.0], [0.0, 3.0, 0.0]);
    let v = m.push_vertex(
        Vertex::ThreePlane([mx, my, pz]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.0, 3.0, 0.0]),
        },
    );
    let (p, frame) = m
        .vertex_meet(v)
        .expect("two translation chains still state the world");
    assert_eq!(frame, None, "a world answer carries no chain to replay");
    assert_eq!(
        p.narrow().map(|c| c.map(|x| x.to_f64())),
        Some([2.0, 3.0, 0.0]),
        "the corner is where the two moved planes and the world plane meet"
    );

    // ① A rotated carrier has no rational world name — the road declines, honestly. A quarter
    // turn carries x = 0 to y = 0, so the cache is exact, and pairing it with the *translated*
    // x = 2 wall and z = 0 keeps the triple a proper three-plane point: what declines here is
    // the rotation, not a degeneracy.
    let spin = m.push_motion(
        Motion::Rotate {
            axis: Axis::Z,
            pivot: [r(0); 3],
            angle: Angle::from_deg(r(90)).expect("angle"),
        },
        None,
    );
    let turned = moved(&mut m, px, spin, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
    let v_turned = m.push_vertex(
        Vertex::ThreePlane([turned, mx, pz]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.0, 0.0, 0.0]),
        },
    );
    assert!(
        m.vertex_meet(v_turned).is_none(),
        "a rotated carrier has no world name to meet with"
    );

    // ② The shared-chain door still answers through itself: same leaf on both moved carriers.
    let mx2 = moved(&mut m, px, t1, [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
    let my2 = moved(&mut m, py, t1, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
    let v_shared = m.push_vertex(
        Vertex::ThreePlane([mx2, my2, pz]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.0, 0.0, 0.0]),
        },
    );
    let (_, frame) = m.vertex_meet(v_shared).expect("one chain, one frame");
    assert_eq!(
        frame,
        Some(t1),
        "a corner the shared-chain door answers keeps taking that door"
    );
}

/// ★ **Two doors, one fact.** `through_meets` used to solve each vertex inline; that body is
/// now `vertex_meet`, and the three-vertex door is its caller. The two must agree point for
/// point — otherwise the extraction quietly created a second spelling of the solve, which is
/// this repo's dominant defect shape.
///
/// The frame comes back too, and on an unmoved box it is the world (`None`). That is the bit a
/// consumer comparing against world coordinates has to demand: `through_meets` only asks the
/// three to *agree*, which would pass a solid whose points are all stated pre-motion.
#[test]
fn one_vertex_and_three_vertices_solve_the_same_meet() {
    let m = build([0.0; 3], [2.0, 3.0, 5.0]);
    let vs: Vec<Handle<Vertex>> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .collect();
    assert!(vs.len() >= 3, "a box has corners");
    let tri = [vs[0], vs[1], vs[2]];
    let together = m
        .through_meets(tri)
        .expect("a box's corners share the world");
    for (i, v) in tri.iter().enumerate() {
        let (alone, frame) = m.vertex_meet(*v).expect("a corner is a three-plane point");
        assert_eq!(
            alone.narrow(),
            together[i].narrow(),
            "vertex {i} solves differently through the two doors"
        );
        assert_eq!(
            frame, None,
            "an unmoved box states its corners in the world"
        );
    }
}

fn build(min: [f64; 3], max: [f64; 3]) -> Model {
    let mut m = Model::new();
    m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
    m.rebuild_adjacency();
    m
}

fn he_start(m: &Model, he: HalfEdge) -> Handle<Vertex> {
    let [a, b] = m.edge(he.edge).vertices;
    if he.forward { a } else { b }
}
fn he_end(m: &Model, he: HalfEdge) -> Handle<Vertex> {
    let [a, b] = m.edge(he.edge).vertices;
    if he.forward { b } else { a }
}
fn face_plane_normal(m: &Model, f: &Face) -> Vector3 {
    match m.surface_cache(f.surface) {
        nacre_geom::Surface::Plane(p) => p.normal(),
        // Planar-only helper: callers filter to plane faces (caps), never cylinders.
        nacre_geom::Surface::Cylinder(_) => {
            unreachable!("face_plane_normal called on a curved face")
        }
    }
}
fn face_centroid(m: &Model, f: &Face) -> Point3 {
    let pts: Vec<Point3> = f
        .outer
        .half_edges
        .iter()
        .map(|he| m.vertex_point(he_start(m, *he)))
        .collect();
    Point3::centroid(&pts).unwrap()
}

// --- golden ---

#[test]
fn cuboid_counts() {
    let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    assert_eq!(m.vertex_count(), 8);
    assert_eq!(m.edge_count(), 12);
    assert_eq!(m.face_count(), 6);
    assert_eq!(m.shell_count(), 1);
    assert_eq!(m.solid_count(), 1);
    // 6, composed with the seeds: the origin box's bottom/left/front
    // faces intern onto the three seeded world planes (same name, same handle), so the
    // arena holds 3 seeds + 3 fresh (top/back/right). Seeding adds nothing here precisely
    // because the seeds are these planes.
    assert_eq!(m.surfaces.len(), 6);
    assert_eq!(
        m.edge_cache.len(),
        m.edge_count(),
        "the curve cache stays index-parallel"
    );
}

#[test]
fn cuboid_corner_points() {
    let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
    let expected = vec![
        [-2.0, 1.0, 0.0],
        [3.0, 1.0, 0.0],
        [3.0, 4.0, 0.0],
        [-2.0, 4.0, 0.0],
        [-2.0, 1.0, 10.0],
        [3.0, 1.0, 10.0],
        [3.0, 4.0, 10.0],
        [-2.0, 4.0, 10.0],
    ];
    let got: Vec<[f64; 3]> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .map(|h| (h, m.vertex(h)))
        .map(|(vh, _)| m.vertex_point(vh).as_array())
        .collect();
    assert_eq!(got, expected);
}

#[test]
fn face_normals_point_outward() {
    let m = build([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
    let center = Point3::origin().lerp(Point3::from_array([2.0, 3.0, 4.0]), 0.5);
    let mut i = 0u32;
    while let Some(h_) = m.face_handle_at(i) {
        i += 1;
        let f = m.face(h_);
        let outward = face_plane_normal(&m, f).dot(face_centroid(&m, f) - center);
        assert!(outward > 0.0);
    }
}

#[test]
fn every_edge_used_twice_opposite() {
    let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    assert_eq!(m.adj.edge_uses.len(), 12);
    for uses in m.adj.edge_uses.values() {
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0].1, uses[1].1);
    }
}

#[test]
fn every_vertex_incident_to_three_edges() {
    let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    assert_eq!(m.adj.vertex_edges.len(), 8);
    for edges in m.adj.vertex_edges.values() {
        assert_eq!(edges.len(), 3);
    }
}

#[test]
fn outer_loops_are_closed() {
    let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut i = 0u32;
    while let Some(h_) = m.face_handle_at(i) {
        i += 1;
        let f = m.face(h_);
        let hes = &f.outer.half_edges;
        assert_eq!(hes.len(), 4);
        for i in 0..hes.len() {
            assert_eq!(he_end(&m, hes[i]), he_start(&m, hes[(i + 1) % hes.len()]));
        }
    }
}

/// ★★★★ **The write doors bite** — the negative control for [`Model::push_face`] and
/// [`Model::push_shell`].
///
/// ★ It is measured that *every face the product builds closes its loop*: one temporary
/// assertion, the whole suite, 1324 tests, a single red — and that one was a fixture whose
/// own comment called it a franken-face. That is a **different proposition** from *the
/// assertion refuses a face that does not close*. The first says the population is clean;
/// the second says the door has teeth. Only this test says the second, and without it the
/// doors could assert nothing at all and every green would still be green.
///
/// Each case violates **exactly one** assertion and the panic **message is checked**, because
/// asking only "did it panic" lets a case go green for the wrong reason: a face
/// cloned out of a *throwaway* model makes walking its loop hit
/// `Store`'s cross-model handle guard before ever reaching the door's own assertion. The face
/// comes from the very model it is pushed into, and only the out-of-bounds handles are
/// strangers — those are read with `.index()`, which no guard sees. The assertions are
/// `debug_assert`, so this test is `cfg(debug_assertions)` — the shape the foreign-handle
/// lock above already uses. Panic output is left unsuppressed on purpose: swapping the
/// panic hook is global state, and the suite runs its tests in parallel.
#[test]
#[cfg(debug_assertions)]
fn the_write_doors_refuse_what_their_invariants_forbid() {
    // The panic message is checked, not just the panic: the claim above is that each case
    // trips *its own* assertion, and only the message can say which one bit.
    fn refuses(what: &str, expect: &str, call: impl FnOnce()) {
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(call));
        let e = r
            .err()
            .unwrap_or_else(|| panic!("{what}: the door let an invalid cell into the arena"));
        let msg = e
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| e.downcast_ref::<&str>().copied())
            .unwrap_or("<non-string panic>");
        assert!(
            msg.contains(expect),
            "{what}: tripped a different assertion — wanted {expect:?}, got {msg:?}"
        );
    }
    let cube = || build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let sample = |m: &Model| {
        m.face(m.face_handle_at(0).expect("a cuboid has faces"))
            .clone()
    };

    // A bigger model mints handles this one does not hold. It is the only way to name an
    // out-of-bounds cell: `handle_at` answers `None` past the end, by design.
    let mut big = cube();
    let _ = big.add_cuboid(
        Point3::from_array([2.0, 2.0, 2.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let stranger_surface = big
        .surface_handle_at((big.surface_count() - 1) as u32)
        .expect("the bigger model has surfaces");
    let stranger_face = big
        .face_handle_at((big.face_count() - 1) as u32)
        .expect("the bigger model has faces");

    // (1) a face names a surface the arena does not hold
    let mut m = cube();
    let mut f = sample(&m);
    f.surface = stranger_surface;
    refuses(
        "push_face / surface in bounds",
        "a face names a surface the arena does not hold",
        || {
            m.push_face(f);
        },
    );

    // (2) a face's outer loop has no half-edges
    let mut m = cube();
    let mut f = sample(&m);
    f.outer.half_edges.clear();
    refuses(
        "push_face / outer loop non-empty",
        "a face's outer loop has no half-edges",
        || {
            m.push_face(f);
        },
    );

    // (3) a face's inner loop has no half-edges — the outer one is left intact so this
    //     case cannot be carried by (2).
    let mut m = cube();
    let mut f = sample(&m);
    f.inner.push(Loop {
        half_edges: Vec::new(),
    });
    refuses(
        "push_face / inner loop non-empty",
        "a face's inner loop has no half-edges",
        || {
            m.push_face(f);
        },
    );

    // (4) a face loop does not close: flipping one half-edge swaps its end for its start,
    //     which is precisely what walking the loop is there to catch.
    let mut m = cube();
    let mut f = sample(&m);
    f.outer.half_edges[0].forward = !f.outer.half_edges[0].forward;
    refuses(
        "push_face / loop closes",
        "a face loop does not close",
        || {
            m.push_face(f);
        },
    );

    // (5) a shell with no faces bounds nothing
    let mut m = cube();
    refuses(
        "push_shell / non-empty",
        "a shell with no faces bounds nothing",
        || {
            m.push_shell(Shell { faces: Vec::new() });
        },
    );

    // (6) a shell names a face the arena does not hold
    let mut m = cube();
    refuses(
        "push_shell / faces in bounds",
        "a shell names a face the arena does not hold",
        || {
            m.push_shell(Shell {
                faces: vec![stranger_face],
            });
        },
    );
}

/// ★★★ **A foreign handle does not read a cache** — the guard [`Store::get`] owns, kept on
/// the index-parallel cache reads too ([`Model::debug_guard`]).
///
/// ★ This is a **regression test with a measured history**: `Model::surface` used to be a
/// `Store::get` and had the guard for free; the arena flip made it a `Vec` index and dropped
/// it, and `vertex_point`/`edge_curve` had never had it — a handle from another model reached
/// all three and answered with the wrong cell. The guard is `cfg(debug_assertions)`, so the
/// test is too. Panic output is left unsuppressed on purpose: swapping the panic hook is
/// global state, and the suite runs tests in parallel.
#[test]
#[cfg(debug_assertions)]
fn a_foreign_handle_cannot_read_a_cache() {
    let mut a = Model::new();
    let mut b = Model::new();
    for m in [&mut a, &mut b] {
        let _ = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
    }
    let surf = b.surface_handle_at(4).expect("b has surfaces");
    let vert = b.vertex_handle_at(3).expect("b has vertices");
    let edge = b.edges.handle_at(3).expect("b has edges");

    let refuses = |what: &str, call: &dyn Fn()| {
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(call));
        assert!(
            r.is_err(),
            "{what}: a handle minted by another model must not read this model's cell"
        );
    };
    refuses("surface_cache", &|| {
        let _ = a.surface_cache(surf);
    });
    refuses("surface", &|| {
        let _ = a.surface(surf);
    });
    refuses("vertex_point", &|| {
        let _ = a.vertex_point(vert);
    });
    refuses("vertex_cache", &|| {
        let _ = a.vertex_cache(vert);
    });
    refuses("edge_curve", &|| {
        let _ = a.edge_curve(edge);
    });
}

/// ★★★★ **What the flip bought: the surface cache is writable, and the truth is not.**
///
/// Before the arena held the truth, a refinement pass had nowhere to write — the realization
/// lived in a `Store`, which is append-only and sealed, and the only mutable copy was the
/// *truth*. This test is the warrant, and it could not have been written then: it rewrites a
/// cache entry in place and reads it back through [`Model::surface_cache`], while the truth the same
/// handle names is untouched.
///
/// ★ It writes through the private field on purpose — that is the door the refinement pass
/// (`realize_surface`) will take, from inside this crate. Outside, there is no door at all,
/// which is the other half of the warrant and what [`Model::surface_cache`]'s `compile_fail` pins.
#[test]
fn the_surface_cache_is_writable_and_the_truth_is_not() {
    let mut m = Model::new();
    let h = m.world_plane(nacre_exact::Axis::Z);
    let truth_before = m.surface(h).clone();
    let cache_before = m.surface_cache(h).clone();

    // A different plane in the same slot — what a refinement at a higher precision does in
    // kind, if not in size.
    let refined = nacre_geom::Surface::Plane(
        Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.25]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .expect("a unit normal names a plane"),
    );
    m.surface_cache[h.index() as usize] = SurfaceCache {
        realized: refined.clone(),
    };

    assert_eq!(
        *m.surface_cache(h),
        refined,
        "the cache took the new realization"
    );
    assert_ne!(
        *m.surface_cache(h),
        cache_before,
        "and it is not the old one"
    );
    assert_eq!(
        *m.surface(h),
        truth_before,
        "the truth a handle names is untouched — it is the arena entry, and the arena is sealed"
    );
    assert_eq!(
        m.surface_cache.len(),
        m.surface_count(),
        "index-parallel, still"
    );
}

/// ★★ The «discard and regenerate» warrant: the edge-curve cache rebuilt from the
/// carriers and endpoints is bit-identical to the one `push_edge` filled eagerly — proof
/// that nothing in it was truth.
#[test]
fn edge_cache_discard_and_regenerate_bit_identical() {
    let mut m = Model::new();
    m.add_cuboid(
        Point3::from_array([-2.0, 1.0, 0.0]),
        Point3::from_array([3.0, 4.0, 10.0]),
    );
    m.add_cylinder(
        Point3::from_array([8.0, 0.0, 0.0]),
        Vector3::from_array([0.3, -0.4, 1.0]),
        1.25,
        2.5,
    );
    let snapshot = |m: &Model| -> Vec<Curve> {
        (0..m.edge_count() as u32)
            .filter_map(|i| m.edge_handle_at(i))
            .map(|h| (h, m.edge(h)))
            .map(|(eh, _)| m.edge_curve(eh).clone())
            .collect()
    };
    let before = snapshot(&m);
    m.rebuild_edge_cache();
    // ⚠ True of a model nothing has refined. [`Model::refine_vertex_cache`] moves coordinates,
    // and a rebuild after *that* is not a no-op — which is the whole reason the paying caller
    // re-derives. This asserts the derivation is stable, not that rebuilding is always free.
    assert_eq!(
        before,
        snapshot(&m),
        "regeneration must reproduce the cache bit for bit"
    );
}

#[test]
fn edge_endpoints_lie_on_their_curve() {
    let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let e = m.edge(eh);
        let curve = m.edge_curve(eh);
        let [a, b] = e.vertices;
        assert!(curve.contains(m.vertex_point(a), 1e-9));
        assert!(curve.contains(m.vertex_point(b), 1e-9));
    }
}

#[test]
fn euler_poincare_holds() {
    let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let (v, e, f) = (
        m.vertex_count() as i64,
        m.edge_count() as i64,
        m.face_count() as i64,
    );
    // V − E + F = 2(S − G) + L_i, with S=1, G=0, L_i=0. Formal validate:
    // nacre-validate (next unit).
    assert_eq!(v - e + f, 2);
}

// --- cylinder (seam b-rep) --- (validate-clean lives in nacre-validate)

fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
    let mut m = Model::new();
    m.add_cylinder(Point3::from_array(base), Vector3::from_array(axis), r, h);
    m.rebuild_adjacency();
    m
}

#[test]
fn cylinder_counts_and_euler() {
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    assert_eq!(m.vertex_count(), 2);
    assert_eq!(m.edge_count(), 3);
    assert_eq!(m.face_count(), 3);
    assert_eq!(m.shell_count(), 1);
    assert_eq!(m.solid_count(), 1);
    // Euler χ = V − E + F = 2 (one shell, genus 0, no inner loops).
    assert_eq!(
        m.vertex_count() as i64 - m.edge_count() as i64 + m.face_count() as i64,
        2
    );
}

#[test]
fn cylinder_geometry() {
    // +Z axis, r=2, h=5. Seam direction is X (least-aligned axis of +Z).
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    let pts: Vec<[f64; 3]> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .map(|h| (h, m.vertex(h)))
        .map(|(vh, _)| m.vertex_point(vh).as_array())
        .collect();
    // Seam direction for +Z is any_perpendicular([0,0,1]) = X×Z = [0,-1,0], so
    // the seam vertices sit at radius 2 along −Y, at z=0 and z=5.
    assert_eq!(pts, vec![[0.0, -2.0, 0.0], [0.0, -2.0, 5.0]]);
    // Two rim circles carry a Circle; the straight seam a Line.
    let mut circles = 0;
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        if let Curve::Circle(c) = m.edge_curve(eh) {
            assert_eq!(c.radius(), 2.0);
            circles += 1;
        }
    }
    assert_eq!(circles, 2);
}

#[test]
fn cylinder_caps_point_outward() {
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    // Planar caps only (the lateral cylindrical face has no single normal).
    let mut caps = 0;
    let mut i = 0u32;
    while let Some(h_) = m.face_handle_at(i) {
        i += 1;
        let f = m.face(h_);
        if matches!(m.surface_cache(f.surface), nacre_geom::Surface::Plane(_)) {
            let n = face_plane_normal(&m, f).as_array();
            // Bottom cap → −Z, top cap → +Z (outward along the axis).
            assert!(n == [0.0, 0.0, -1.0] || n == [0.0, 0.0, 1.0]);
            caps += 1;
        }
    }
    assert_eq!(caps, 2);
}

/// ★★ **A seam vertex states its two carriers** — `OnSeam([lateral, its own cap])`.
/// The pair pins the rim circle; the cylinder's `ref_dir` truth pins the point
/// (see `Vertex::OnSeam`). The lateral surface must be the cylinder, the other its cap.
#[test]
fn a_seam_vertex_states_its_rim_carriers() {
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    let defs: Vec<_> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .map(|h| *m.vertex(h))
        .collect();
    assert_eq!(
        defs.len(),
        2,
        "a cylinder has exactly its two seam vertices"
    );
    for (i, d) in defs.iter().enumerate() {
        let Vertex::OnSeam([a, b]) = d else {
            panic!("a seam vertex carries OnSeam, got {d:?}")
        };
        assert!(
            matches!(m.surface_cache(*a), nacre_geom::Surface::Cylinder(_)),
            "first carrier is the lateral cylinder"
        );
        let nacre_geom::Surface::Plane(p) = m.surface_cache(*b) else {
            panic!("second carrier is the cap plane")
        };
        // Bottom vertex names the bottom cap (through z = 0), top the top cap (z = 5).
        let z = if i == 0 { 0.0 } else { 5.0 };
        assert_eq!(p.distance(Point3::from_array([0.0, 0.0, z])), 0.0);
    }
}

#[test]
fn cylinder_seam_edge_is_self_adjacent() {
    // The novel topology: the seam edge is used twice by the SAME lateral
    // face with opposite orientation (a valid non-manifold-looking but
    // manifold seam). Every edge is still used exactly twice, opposite.
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    let mut seam_uses = None;
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let uses = &m.adj.edge_uses[&eh];
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0].1, uses[1].1); // opposite orientation
        // The seam's discriminator IS the invariant: self-adjacent carriers.
        if m.edge(eh).surfaces[0] == m.edge(eh).surfaces[1] {
            seam_uses = Some(uses.clone());
        }
    }
    let uses = seam_uses.expect("a seam line edge exists");
    assert_eq!(uses[0].0, uses[1].0); // both uses are the same (lateral) face
}

// --- proptest ---

fn box_strategy() -> impl Strategy<Value = ([f64; 3], [f64; 3])> {
    (
        prop::array::uniform3(-1e3f64..1e3),
        prop::array::uniform3(1e-2f64..1e3),
    )
        .prop_map(|(min, ext)| {
            let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
            (min, max)
        })
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn prop_cuboid_structural((min, max) in box_strategy()) {
        let m = build(min, max);
        prop_assert_eq!(m.vertex_count(), 8);
        prop_assert_eq!(m.edge_count(), 12);
        prop_assert_eq!(m.face_count(), 6);
        prop_assert_eq!(m.adj.edge_uses.len(), 12);
        for uses in m.adj.edge_uses.values() {
            prop_assert_eq!(uses.len(), 2);
            prop_assert_ne!(uses[0].1, uses[1].1);
        }
        prop_assert_eq!(m.adj.vertex_edges.len(), 8);
        for edges in m.adj.vertex_edges.values() {
            prop_assert_eq!(edges.len(), 3);
        }
    }

    #[test]
    fn prop_cuboid_outward_normals((min, max) in box_strategy()) {
        let m = build(min, max);
        let center = Point3::from_array(min).lerp(Point3::from_array(max), 0.5);
        let mut i = 0u32;
        while let Some(h_) = m.face_handle_at(i) {
            i += 1;
            let f = m.face(h_);
            prop_assert!(face_plane_normal(&m, f).dot(face_centroid(&m, f) - center) > 0.0);
        }
    }

    #[test]
    fn prop_cuboid_endpoints_on_curves((min, max) in box_strategy()) {
        let m = build(min, max);
        let scale = 1e-6 * (max.iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        let mut i = 0u32;
        while let Some(eh) = m.edge_handle_at(i) {
            i += 1;
            let e = m.edge(eh);
            let curve = m.edge_curve(eh);
            let [a, b] = e.vertices;
            prop_assert!(curve.contains(m.vertex_point(a), scale));
            prop_assert!(curve.contains(m.vertex_point(b), scale));
        }
    }

    #[test]
    fn prop_cylinder_structural(
        base in prop::array::uniform3(-1e3f64..1e3),
        axis in prop::array::uniform3(-1.0f64..1.0),
        r in 0.5f64..10.0,
        h in 0.1f64..10.0,
    ) {
        let axis = Vector3::from_array(axis);
        prop_assume!(axis.norm() > 0.1); // skip near-zero axes
        let mut m = Model::new();
        m.add_cylinder(Point3::from_array(base), axis, r, h);
        m.rebuild_adjacency();
        prop_assert_eq!(m.vertex_count(), 2);
        prop_assert_eq!(m.edge_count(), 3);
        prop_assert_eq!(m.face_count(), 3);
        // Every edge used exactly twice, opposite orientation (seam included).
        for uses in m.adj.edge_uses.values() {
            prop_assert_eq!(uses.len(), 2);
            prop_assert_ne!(uses[0].1, uses[1].1);
        }
    }

    /// The same structure over the population that **actually stressed the exact
    /// arithmetic**: an axis a hair off `ẑ`, whose tiny components carry a full f64's worth
    /// of decimal digits and so lift to rationals with ~10²⁰ denominators. Squaring one of
    /// those leaves `i128`, which is what the parallelism test used to refuse — and
    /// `add_cylinder` read that refusal as "degenerate cylinder" and panicked.
    ///
    /// ★ The sibling above cannot stand in for this: its uniform axis reaches this family
    /// about **0.01%** of the time (measured), which is why the defect sat green for a
    /// milestone and then surfaced from one unlucky seed. Here it is ~80%.
    #[test]
    fn prop_cylinder_near_axis_aligned_is_built_not_refused(
        u in -1.0f64..1.0,
        v in -1.0f64..1.0,
        k in 1i32..9,
        j in 1i32..9,
        r in 0.5f64..10.0,
        h in 0.1f64..10.0,
    ) {
        let axis = Vector3::from_array([u * 10f64.powi(-k), v * 10f64.powi(-j), 1.0]);
        let mut m = Model::new();
        m.add_cylinder(Point3::origin(), axis, r, h);
        m.rebuild_adjacency();
        prop_assert_eq!(m.vertex_count(), 2);
        prop_assert_eq!(m.edge_count(), 3);
        prop_assert_eq!(m.face_count(), 3);
        for uses in m.adj.edge_uses.values() {
            prop_assert_eq!(uses.len(), 2);
            prop_assert_ne!(uses[0].1, uses[1].1);
        }
    }
}

// --- reversed_shell (M5 containment building block) ---

#[test]
fn orientation_flip_is_involution() {
    assert_eq!(Orientation::Forward.flipped(), Orientation::Reversed);
    assert_eq!(Orientation::Reversed.flipped(), Orientation::Forward);
    assert_eq!(
        Orientation::Forward.flipped().flipped(),
        Orientation::Forward
    );
}

#[test]
fn reversed_shell_toggles_orientation_and_reverses_loops() {
    let mut m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let outer = m.solid(m.live_solids[0]).outer;
    let f0 = m.shell(outer).faces[0];
    let orig = m.face(f0).clone();

    let rev_shell = m.reversed_shell(outer);
    // Fresh cells (not reused faces), fresh shell.
    assert_ne!(rev_shell, outer);
    let rf0 = m.shell(rev_shell).faces[0];
    assert_ne!(rf0, f0);
    let rev = m.face(rf0);

    assert_eq!(rev.surface, orig.surface); // surface reused
    assert_eq!(rev.orientation, orig.orientation.flipped());
    let n = orig.outer.half_edges.len();
    assert_eq!(rev.outer.half_edges.len(), n);
    // Reversed winding: he[i] mirrors orig[n-1-i] with the edge reused and
    // the traversal direction flipped.
    for i in 0..n {
        let o = orig.outer.half_edges[n - 1 - i];
        let r = rev.outer.half_edges[i];
        assert_eq!(r.edge, o.edge);
        assert_ne!(r.forward, o.forward);
    }
}

#[test]
fn reversed_shell_is_a_valid_manifold() {
    // Reversing every face's winding preserves the b-rep manifold: each edge
    // is still used by exactly two faces with opposed half-edges. The
    // reversed shell reuses the cube's edges, but the source solid is
    // superseded, so `Adjacency` (reachable-scoped) counts only the reversed
    // faces — no 4-use false positive.
    let mut m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let outer = m.solid(m.live_solids[0]).outer;
    let rev = m.reversed_shell(outer);
    let s = m.push_solid(Solid {
        outer: rev,
        cavities: vec![],
    });
    m.live_solids.retain(|&h| h == s); // supersede the original cube
    m.rebuild_adjacency();
    assert_eq!(m.adj.edge_uses.len(), 12);
    for uses in m.adj.edge_uses.values() {
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0].1, uses[1].1); // opposite forward
    }
}

/// ★★★★★ **A plane is named by its points, and by nothing else.**
///
/// The three assertions are the three ways this can go wrong. A name has to come out for an
/// ordinary plane; it has to be the plane the points are **on**, not some other; and it has to
/// come out even where the derivation's narrow route cannot reach — that last one is what a
/// producer used to lose a name to, and losing a name closes the exact road for everything
/// built on that plane.
///
/// ★ There is no "wrong name" case left to test. Coefficients are no longer something a
/// producer can hand in, so a plane cannot be stated twice — the state the old agreement filter
/// watched for is now unspellable.
#[test]
fn a_plane_is_named_by_its_points() {
    use nacre_exact::Rat;
    let r = Rat::from_int;
    let pl = |z: f64| {
        nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, z]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap()
    };

    let mut m = Model::new();
    let pts = [[r(0), r(0), r(0)], [r(1), r(0), r(0)], [r(0), r(1), r(0)]];
    let (ok, _) = m.push_plane(pl(0.0), pts, None);
    assert_eq!(
        m.surface_name.get(&ok),
        Some(&nacre_exact::PlaneName::Narrow([r(0), r(0), r(1), r(0)])),
        "the plane z = 0 was not named, or was named as something else"
    );

    // ★★★★ **Wide enough that the narrow derivation gives up.** Decimal arithmetic reduces to
    // denominators that are powers of two on one coordinate and powers of five on another, and
    // a triple of *coprime* denominators needs their product to state its plane — which
    // `(b − a) × (c − a)` then needs squared.
    //
    // ★ An earlier spelling here used points chosen to overflow the old **agreement check**
    // (`c · p`) and asserted the derivation gave up on them too. It does not: those are
    // different products, and the derivation went through. Two propositions, one fixture.
    let q = |n: i128, d: i128| nacre_exact::Rat::new(n, d).unwrap();
    let wide_pts = [
        [q(1, 1 << 53), r(0), r(0)],
        [r(0), q(1, 5i128.pow(23)), r(0)],
        [r(0), r(0), q(1, (1 << 40) * 5i128.pow(11))],
    ];
    assert_eq!(
        nacre_exact::plane_through_points(wide_pts[0], wide_pts[1], wide_pts[2]),
        None,
        "the narrow route was expected to overflow here — the case has stopped being the case"
    );
    let (wide, _) = m.push_plane(pl(2.0), wide_pts, None);
    // ★ The stored value carries the proposition — `Narrow` says "fits i128" directly.
    // (This used to compare the global counter before/after, which races against other
    // tests pushing wide planes in parallel; the value cannot.)
    let name = *m
        .surface_name
        .get(&wide)
        .expect("a plane the narrow route cannot reach went unnamed")
        .narrow()
        .expect("this fixture's canonical answer is small — it must be stored Narrow");
    // ★★ And it names *these* points' plane — exact rationals, no tolerance.
    for p in &wide_pts {
        let mut acc = name[3];
        for k in 0..3 {
            acc = acc
                .checked_add(name[k].checked_mul(p[k]).expect("no overflow"))
                .expect("no overflow");
        }
        assert_eq!(
            acc,
            r(0),
            "a recorded point is off the name derived from it"
        );
    }

    // ★ (A plane with no points cannot be pushed at all any more — the "no points, no
    // name" arm retired with the old API; the type is the assertion now.)
}

/// ★★★★★ **A wide name interns — and opens no shortcut**.
///
/// The fixture is a triple whose **canonical answer** exceeds `i128` (cross-product terms
/// multiply two ~2^90 coprime numerators — the same triple `nacre-exact` locks as `Wide`).
/// Without a name two statements of such a plane would be two handles;
/// it interns like any other.
///
/// ★ Note the population: an axis-aligned cuboid — even on seventeen-digit corners — is NOT
/// wide, because its canonical answers are tiny (`[10^13, 0, 0, −c]`); only the *narrow
/// route's intermediates* overflow there, and the big route names those.
/// ★★ Measured: today's sketch→extrude walls cap out
/// around ~115 bits (the decimal window bounds the products), so a wide name currently
/// arises only from hand-stated triples like this one — datum planes are the
/// production source.
///
/// What a wide name still does not do: `narrow()` is `None`, so the **narrow shortcuts**
/// (`base_rat`, integer Shewchuk, `Isometry` transport) decline exactly as they do on a
/// missing name. A wide name **does** host a sketch frame — through the
/// arbitrary-precision road (`FramePlacement::Canonical`, locked in nacre-ops) — so what
/// this pins is the shortcut boundary, not a frame one.
#[test]
fn a_wide_plane_interns_but_opens_no_narrow_shortcut() {
    let q = |n: i128, d: i128| nacre_exact::Rat::new(n, d).unwrap();
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let a = [q(big1, 3), q(big2, 7), q(0, 1)];
    let b = [q(-big2, 5), q(big1, 11), q(0, 1)];
    let c = [q(1, 13), q(1, 17), q(1, 19)];
    // Fixture qualification: the narrow route gives up on these points.
    assert_eq!(
        nacre_exact::plane_through_points(a, b, c),
        None,
        "the narrow route was expected to overflow here — the case has stopped being the case"
    );

    let pl = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let mut m = Model::new();
    let before = WIDE_PLANES.load(std::sync::atomic::Ordering::Relaxed);
    let (first, _) = m.push_plane(pl, [a, b, c], None);
    assert!(
        WIDE_PLANES.load(std::sync::atomic::Ordering::Relaxed) > before,
        "an answer past i128 must be counted as wide"
    );
    let name = m
        .surface_name
        .get(&first)
        .expect("a plane too wide for i128 must still be named");
    assert!(
        name.narrow().is_none(),
        "a wide name must not open the narrow shortcuts (base_rat, Shewchuk)"
    );

    // ★ The same plane stated again — permuted, even — is the same handle now.
    let (second, _) = m.push_plane(pl, [a, b, c], None);
    assert_eq!(first, second, "one wide plane, one handle");
    let (permuted, _) = m.push_plane(pl, [b, c, a], None);
    assert_eq!(
        first, permuted,
        "two spellings of one wide plane must intern"
    );
}

/// ★★ **A statement the name key cannot hold interns by the statement**.
///
/// The fixture reaches namelessness through `OnSeam` carriers — the cheapest population
/// `plane_name_through` declines inside this crate. Semantically an ops producer would
/// refuse this particular datum by cause; the door's contract is narrower ("store the
/// statement, once") and holds for every nameless reason identically, which is what is
/// pinned here. The real mixed-frame population is exercised end-to-end in `nacre-ops`.
#[test]
fn a_nameless_through_statement_interns_by_its_statement() {
    let mut m = Model::new();
    m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    m.add_cuboid(
        Point3::from_array([4.0, 0.0, 0.0]),
        Point3::from_array([5.0, 1.0, 1.0]),
    );
    let mut vs: Vec<Handle<Vertex>> = Vec::new();
    let mut i = 0u32;
    while let Some(h) = m.vertex_handle_at(i) {
        i += 1;
        if matches!(*m.vertex(h), Vertex::OnSeam(_)) && vs.len() < 2 {
            vs.push(h);
        } else if matches!(*m.vertex(h), Vertex::ThreePlane(_)) && vs.len() == 2 {
            vs.push(h);
            break;
        }
    }
    let mut triple: [Handle<Vertex>; 3] = [vs[0], vs[1], vs[2]];
    triple.sort_by_key(|v| v.index());
    assert!(
        m.plane_name_through(triple).is_none(),
        "the fixture must be nameless, or this test measures the name road"
    );

    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.5]),
        Vector3::from_array([0.0, 0.0, -1.0]),
    )
    .unwrap();
    let before = m.surface_count();
    let (h, flipped) = m.push_plane_through(cache, triple, None);
    assert!(!flipped);
    assert!(
        !m.surface_name.contains_key(&h),
        "a nameless statement must not invent a name"
    );
    assert_eq!(m.surface_count(), before + 1, "stored once");

    // The same statement again — one handle; and a cache built facing the other way is the
    // same plane with `flipped` reported, exactly as the name road reports it.
    let (again, flipped_same) = m.push_plane_through(cache, triple, None);
    assert_eq!(h, again, "one statement, one handle");
    assert!(!flipped_same);
    let reversed = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let (still, flipped_now) = m.push_plane_through(reversed, triple, None);
    assert_eq!(h, still, "direction is not part of the statement");
    assert!(
        flipped_now,
        "but the survivor's other-way cache is reported"
    );
    assert_eq!(m.surface_count(), before + 1, "and nothing new was stored");

    // A different motion is a different statement — the same rule the name key keeps.
    let node = m.push_motion(
        Motion::Translate {
            offset: [
                nacre_exact::Rat::from_int(1),
                nacre_exact::Rat::from_int(0),
                nacre_exact::Rat::from_int(0),
            ],
        },
        None,
    );
    let (moved, _) = m.push_plane_through(cache, triple, Some(node));
    assert_ne!(h, moved, "the motion belongs in the statement key");
}

/// ★★★ **The three world planes are born with the model** — deterministic handles,
/// canonical truth, −axis caches — and `Default` is the same seeded model (the unseeded
/// back door is closed).
#[test]
fn a_new_model_carries_the_three_world_planes() {
    let r = nacre_exact::Rat::from_int;
    for m in [Model::new(), Model::default()] {
        assert_eq!(m.surface_count(), 3);
        let want = [
            (
                nacre_exact::Axis::Z,
                [r(0), r(0), r(1), r(0)],
                [0.0, 0.0, -1.0],
            ),
            (
                nacre_exact::Axis::X,
                [r(1), r(0), r(0), r(0)],
                [-1.0, 0.0, 0.0],
            ),
            (
                nacre_exact::Axis::Y,
                [r(0), r(1), r(0), r(0)],
                [0.0, -1.0, 0.0],
            ),
        ];
        for (i, (axis, name, cache_n)) in want.into_iter().enumerate() {
            let h = m.world_plane(axis);
            assert_eq!(h.index() as usize, i, "deterministic seed handles");
            assert_eq!(
                m.surface_name.get(&h),
                Some(&nacre_exact::PlaneName::Narrow(name))
            );
            assert!(matches!(m.surface(h), Surface::Plane { motion: None, .. }));
            let nacre_geom::Surface::Plane(pl) = m.surface_cache(h) else {
                panic!("a seed is a plane")
            };
            assert_eq!(
                pl.normal().as_array(),
                cache_n,
                "the cache points down the −axis (the base-cap sense)"
            );
        }
    }
}

/// ★★★★★ **The anchor derivation's proposition, at its smallest**: a plane's realized cache
/// does not depend on **which point of it the producer anchored at**.
///
/// One truth, two producers: each builds its own `Plane` at a different point of that same
/// plane, facing the same way. Their values differ — asserted first, on what the *producers*
/// built, because the door now overwrites what it is handed and reading it back would compare
/// the derivation against itself.
///
/// ⚠★★★ **The stronger claim is not made, deliberately** — one
/// plane through two *different* triples demanding one cache. Interning folds two such
/// statements into **one handle and one truth** (the first pusher's), so no model can hold
/// them at once — the kernel never needed that property, and buying it meant a canonical row,
/// which is measured and rejected (32 census result rows moved, an exact-coefficient road
/// lost, `n·n` overflowing for 8 planes). What survives is the reachable half.
///
/// ★ This lock watches the **door**; `plane_anchor.rs`'s
/// `which_surface_caches_two_anchors_leave_disagreeing` watches the **pipeline**. Fix the door
/// while an operation still routes around it and this one is green while that one is red.
#[test]
fn a_derived_plane_cache_does_not_depend_on_which_anchor_states_it() {
    let r = nacre_exact::Rat::from_int;
    let f = |q: [nacre_exact::Rat; 3]| {
        Point3::from_array([q[0].to_f64(), q[1].to_f64(), q[2].to_f64()])
    };
    // 2x + 3y + 6z = 6, stated through its three axis intercepts.
    let pts = [[r(3), r(0), r(0)], [r(0), r(2), r(0)], [r(0), r(0), r(1)]];
    let normal = Vector3::from_array([2.0, 3.0, 6.0]);
    // Two anchors on that plane: the truth's own first point, and a point that is none of the
    // three (2·1.5 + 3·1 + 6·0 = 6). The second is the one a producer picks by accident.
    let anchors = [f(pts[0]), Point3::from_array([1.5, 1.0, 0.0])];
    let mut stated = Vec::new();
    let mut held = Vec::new();
    for a in anchors {
        let mut m = Model::new();
        let cache = nacre_geom::Plane::from_point_normal(a, normal).expect("a nonzero normal");
        let (h, _) = m.push_plane(cache, pts, None);
        stated.push(cache);
        held.push(m.surface_cache(h).clone());
    }
    // ⚠★★★ **Compared as whole caches, not by `coefficients()`.** On this plane both anchors
    // give `d = −6` *exactly*, so the two rows are identical while the producers genuinely
    // disagree — about the anchor, which is the only thing at stake here. A control built on
    // `coefficients()` fires on a fixture that is perfectly good; the sibling file names the
    // same blindness (`plane_anchor.rs`'s `an_axis_aligned_plane_is_anchor_blind`).
    let whole = |p: nacre_geom::Plane| nacre_geom::Surface::Plane(p);
    assert_ne!(
        surface_bits(&whole(stated[0])),
        surface_bits(&whole(stated[1])),
        "the two producers already agree — this fixture measures nothing"
    );
    assert_eq!(
        surface_bits(&held[0]),
        surface_bits(&held[1]),
        "the door kept an anchor the producer chose"
    );
    let nacre_geom::Surface::Plane(d) = &held[0] else {
        panic!("a plane is stored as a plane")
    };
    assert_eq!(
        d.origin().as_array(),
        f(pts[0]).as_array(),
        "the anchor is the truth's first point, realized"
    );
}

/// **A seed's anchor is its truth's first point**, which for the world planes is the origin —
/// and the seeds are the most-interned planes there are, so this is where a derivation that
/// moved anchors would be felt first.
///
/// ⚠ This does not assert that the derived **normal** matches the stored one or carries
/// no negative zero: with the row copied verbatim both are **tautologies**, and a lock
/// that cannot fail is worse than no lock. It aims at the one thing the derivation decides.
#[test]
fn a_derived_seed_plane_anchors_at_its_truths_first_point() {
    let m = Model::new();
    for axis in [
        nacre_exact::Axis::Z,
        nacre_exact::Axis::X,
        nacre_exact::Axis::Y,
    ] {
        let h = m.world_plane(axis);
        let derived = m
            .derive_surface_cache(h)
            .expect("a seed is named and unmoved");
        let nacre_geom::Surface::Plane(d) = &derived else {
            panic!("a seed is a plane")
        };
        assert_eq!(
            d.origin().as_array(),
            [0.0; 3],
            "a seed's truth is `[[0,0,0], u, v]`, so its anchor is the origin"
        );
        assert_eq!(
            surface_bits(&derived),
            surface_bits(m.surface_cache(h)),
            "the seeds were already anchored there — this derivation must not move them"
        );
    }
}

/// A cylinder's **position** descends exactly — `origin` from the statement's rational point
/// and `radius` from its correctly-rounded `√r²`. Its two directions do not: `from_axis`
/// normalizes both, which is the same `√` wall a plane's unit normal sits behind.
#[test]
fn a_derived_cylinder_restates_its_own_axis_and_radius() {
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    let mut seen = 0;
    let mut i = 0u32;
    while let Some(sh) = m.surface_handle_at(i) {
        i += 1;
        if !matches!(m.surface(sh), Surface::Cylinder { .. }) {
            continue;
        }
        seen += 1;
        let derived = m
            .derive_surface_cache(sh)
            .expect("an unmoved cylinder derives");
        assert_eq!(
            surface_bits(&derived),
            surface_bits(m.surface_cache(sh)),
            "an axis-aligned cylinder's statement realizes to the cache it was pushed with"
        );
    }
    assert_eq!(seen, 1, "one lateral surface");
}

/// ★★ The seeds are the interning survivors — an origin cuboid's bottom/left/front
/// faces carry the seed handles, so "the world plane" and "that face's plane" are one
/// surface, stated once.
#[test]
fn an_origin_cuboids_axis_faces_intern_onto_the_seeds() {
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let face_surfaces: Vec<_> = m
        .shell(m.solid(s).outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .collect();
    for axis in [
        nacre_exact::Axis::Z,
        nacre_exact::Axis::X,
        nacre_exact::Axis::Y,
    ] {
        assert!(
            face_surfaces.contains(&m.world_plane(axis)),
            "the {axis:?}-normal face at 0 must be the seed itself"
        );
    }
    // And nothing pointless was minted: 3 seeds + the 3 off-origin faces.
    assert_eq!(m.surface_count(), 6);
}

/// ★★★ **A cylinder's caps record their three exact points**, as `add_cuboid`'s faces do.
/// The direct evidence that the record is a real name and not a dead entry: a cap that
/// shares a plane with a box face **interns to the same handle**, which no point-less
/// surface could ever do.
/// `surface_handle_at` gives back the handle the arena already issued, and nothing more.
///
/// The pair to this is [`Model::surface_cache`]'s `compile_fail` lock, which still refuses to open
/// the store: this adds a *name* for an index round-trip that `iter`-style access could
/// already express, not a new power. Past the end is `None`, because the question it answers
/// is existence.
#[test]
fn surface_handle_at_is_the_handle_the_arena_issued() {
    let mut m = Model::new();
    for axis in [
        nacre_exact::Axis::Z,
        nacre_exact::Axis::X,
        nacre_exact::Axis::Y,
    ] {
        let h = m.world_plane(axis);
        assert_eq!(
            m.surface_handle_at(h.index()),
            Some(h),
            "a seeded plane's index names it back"
        );
    }
    let n = m.surface_count() as u32;
    assert_eq!(m.surface_handle_at(n), None, "past the end is None");
    assert_eq!(m.surface_handle_at(u32::MAX), None);

    // A surface pushed after the seeds is reachable by its own index too.
    let cuboid = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let s = m.face(m.shell(m.solid(cuboid).outer).faces[0]).surface;
    assert_eq!(m.surface_handle_at(s.index()), Some(s));
}

#[test]
fn a_cylinders_caps_record_points_and_intern_with_a_coplanar_face() {
    let mut m = Model::new();
    let cuboid = m.add_cuboid(
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([-1.0, -1.0, 2.0]),
    );
    let cyl = m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    m.rebuild_adjacency();
    // Every planar face of the cylinder carries points and a derived name.
    let planar: Vec<_> = m
        .shell(m.solid(cyl).outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .filter(|&s| matches!(m.surface_cache(s), nacre_geom::Surface::Plane(_)))
        .collect();
    assert_eq!(planar.len(), 2, "two caps");
    for s in &planar {
        assert!(
            matches!(m.surface(*s), Surface::Plane { .. }),
            "a cap without truth"
        );
        assert!(m.surface_name.contains_key(s), "a cap without a name");
    }
    // The top cap lies on `z = 2`, the same plane as the box's top face — one plane, one
    // handle, across two producers.
    let box_top = m
        .shell(m.solid(cuboid).outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .find(|s| {
            m.surface_name
                .get(s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, -2.0])
        })
        .expect("the box top names z = 2");
    assert!(
        planar.contains(&box_top),
        "the coplanar cap and box face must intern to one handle"
    );
}
/// ★★★★★ **A plane can point at three vertices, and it is the same plane as the values.**
///
/// The corners of a unit box name three coordinate planes each, so a datum through them is
/// derivable — and the plane through `(1,0,0)`, `(0,1,0)`, `(0,0,1)` has a name like any
/// other. That the name comes out at all is what makes every road below (frames, far caps,
/// predicates) work unchanged: they read the name, and the name arrives at push.
#[test]
fn a_plane_stated_through_vertices_is_named_like_any_other() {
    let (mut m, vs) = box_corner_vertices();
    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([1.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .unwrap();
    let (h, _) = m.push_plane_through(cache, vs, None);
    let name = m
        .surface_name
        .get(&h)
        .expect("a Through plane must be named");
    assert_eq!(
        name.narrow().map(|c| c.map(|r| r.to_f64())),
        Some([1.0, 1.0, 1.0, -1.0]),
        "x + y + z = 1 is the plane through those three corners"
    );
    assert!(matches!(
        m.surface(h),
        Surface::Plane {
            points: PlanePoints::Through(_),
            ..
        }
    ));
}

/// ★★ **Interning does not care which way the plane was stated.** A `Known` push of the same
/// plane returns the handle the `Through` push made — one plane, one handle, which is the
/// property `PlanePoints` must not break by adding a variant.
#[test]
fn a_through_plane_and_a_known_plane_that_are_one_plane_share_a_handle() {
    let (mut m, vs) = box_corner_vertices();
    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([1.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .unwrap();
    let (through, _) = m.push_plane_through(cache, vs, None);
    let r = |n: i128, d: i128| nacre_exact::Rat::new(n, d).unwrap();
    let (known, _) = m.push_plane(
        cache,
        [
            [r(1, 1), r(0, 1), r(0, 1)],
            [r(0, 1), r(1, 1), r(0, 1)],
            [r(0, 1), r(0, 1), r(1, 1)],
        ],
        None,
    );
    assert_eq!(through, known, "one plane, two statements, one handle");
    // ★ And the first statement wins, as interning is documented to: the truth still points
    // at vertices. Which variant a plane ends up with is arena history, not a promise.
    assert!(matches!(
        m.surface(known),
        Surface::Plane {
            points: PlanePoints::Through(_),
            ..
        }
    ));
}

/// ★★★ **The reference goes backwards in time only** — append-only storage. A datum is pushed after the
/// vertices it names, and those vertices name surfaces pushed before *them* — so the walk
/// plane → vertex → surface strictly descends in index and cannot cycle.
#[test]
fn a_through_plane_points_only_at_older_cells() {
    let (mut m, vs) = box_corner_vertices();
    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([1.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .unwrap();
    let (h, _) = m.push_plane_through(cache, vs, None);
    let Surface::Plane {
        points: PlanePoints::Through(named),
        ..
    } = m.surface(h)
    else {
        unreachable!()
    };
    for v in named {
        assert!(v.index() < m.vertex_count() as u32);
        let Vertex::ThreePlane(tri) = *m.vertex(*v) else {
            unreachable!()
        };
        for s in tri {
            assert!(
                s.index() < h.index(),
                "a vertex's carrier must be older than the datum that names the vertex"
            );
        }
    }
}

/// Three corners of the unit box, sorted — each the meeting of three coordinate planes.
fn box_corner_vertices() -> (Model, [Handle<Vertex>; 3]) {
    let mut m = Model::new();
    m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let mut want = Vec::new();
    for i in 0..m.vertex_count() as u32 {
        let vh = m.vertices.handle_at(i).unwrap();
        let c = m.vertex_point(vh).as_array();
        if c.iter().filter(|x| **x == 1.0).count() == 1 && c.iter().all(|x| *x == 0.0 || *x == 1.0)
        {
            want.push(vh);
        }
    }
    want.sort_by_key(|v| v.index());
    let vs = [want[0], want[1], want[2]];
    (m, vs)
}
/// ★★★★ **Both counters see the second producer too.**
///
/// `push_plane_through` began as a copy of `push_plane`'s tail and the copy silently dropped
/// [`WIDE_PLANES`] and [`SEEDED_HITS`] — the census bridges `stat wide_planes` and
/// `stat seeded_hits` went blind on a road that *feeds* them — a datum plane's name can be
/// `Wide`. Sharing one interning tail is
/// the fix; this is the observation that it worked, and it is two assertions because a
/// counter that cannot move is indistinguishable from a population that never arrives.
#[test]
fn a_through_plane_is_counted_by_both_bridges() {
    use std::sync::atomic::Ordering::Relaxed;

    // ① A `Through` datum landing on a **seed** — three corners of a unit box on `z = 0`.
    let (mut m, vs) = seed_plane_corners();
    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let seeded_before = SEEDED_HITS.load(Relaxed);
    let (h, _) = m.push_plane_through(cache, vs, None);
    assert!(
        (h.index() as usize) < 3,
        "three corners of the box's z = 0 face are the world XY seed"
    );
    assert!(
        SEEDED_HITS.load(Relaxed) > seeded_before,
        "a Through statement that interns onto a seed must be counted like any other"
    );

    // ② A `Through` datum whose name needs the wide vessel. The carriers are stated with
    // coprime ~2^90 numerators, so the canonical coefficients leave `i128` with no content to
    // divide out — the same construction `nacre-exact` uses to reach that arm.
    let wide_before = WIDE_PLANES.load(Relaxed);
    let (mut m2, vs2) = wide_named_corner();
    let cache2 = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .unwrap();
    let (h2, _) = m2.push_plane_through(cache2, vs2, None);
    assert!(
        m2.surface_name
            .get(&h2)
            .is_some_and(|n| n.narrow().is_none()),
        "this fixture must actually produce a wide name, or ② measures nothing"
    );
    assert!(
        WIDE_PLANES.load(Relaxed) > wide_before,
        "a Through statement with a wide name must reach the wide-name bridge"
    );
}

/// Three corners of the unit box that lie on `z = 0` — the world XY seed.
fn seed_plane_corners() -> (Model, [Handle<Vertex>; 3]) {
    let mut m = Model::new();
    m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let mut want = Vec::new();
    for i in 0..m.vertex_count() as u32 {
        let vh = m.vertices.handle_at(i).unwrap();
        if m.vertex_point(vh).as_array()[2] == 0.0 {
            want.push(vh);
        }
    }
    want.sort_by_key(|v| v.index());
    let vs = [want[0], want[1], want[2]];
    (m, vs)
}

/// A model whose three planes meet at a corner and whose *plane through* those corners needs
/// the wide vessel. Built through the test door so the fixture is exactly three planes.
fn wide_named_corner() -> (Model, [Handle<Vertex>; 3]) {
    let r = |n: i128, d: i128| nacre_exact::Rat::new(n, d).unwrap();
    let (b1, b2) = ((1i128 << 90) + 1, (1i128 << 90) + 3);
    let mut m = Model::new();
    // Nine planes: three per vertex, each triple meeting at a point whose coordinates carry
    // the coprime numerators, so the plane through the three points is wide.
    let mut vs = Vec::new();
    for (k, off) in [(0i128, 0i128), (1, 1), (2, 3)].iter().enumerate() {
        let _ = k;
        let (a, c) = (b1 + off.0, b2 + off.1);
        let tri = [
            axis_plane_at(&mut m, 0, r(1, a)),
            axis_plane_at(&mut m, 1, r(1, c)),
            axis_plane_at(&mut m, 2, r(off.0 + 1, b1)),
        ];
        let coord = Point3::from_array([
            1.0 / a as f64,
            1.0 / c as f64,
            (off.0 + 1) as f64 / b1 as f64,
        ]);
        vs.push(m.push_vertex(Vertex::ThreePlane(tri), PointCache::Unrealized { coord }));
    }
    vs.sort_by_key(|v| v.index());
    let out = [vs[0], vs[1], vs[2]];
    (m, out)
}

/// The plane `x_axis = value`, pushed with its exact triple.
fn axis_plane_at(m: &mut Model, axis: usize, value: nacre_exact::Rat) -> Handle<Surface> {
    let z = nacre_exact::Rat::from_int(0);
    let one = nacre_exact::Rat::from_int(1);
    let mut pts = [[z; 3]; 3];
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    for (i, p) in pts.iter_mut().enumerate() {
        p[axis] = value;
        if i == 1 {
            p[u] = one;
        }
        if i == 2 {
            p[v] = one;
        }
    }
    let mut n = [0.0; 3];
    n[axis] = 1.0;
    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array(core::array::from_fn(|k| {
            if k == axis { value.to_f64() } else { 0.0 }
        })),
        nacre_math::Vector3::from_array(n),
    )
    .unwrap();
    m.push_plane(cache, pts, None).0
}

/// **A swap renames the root, because it renames the line.**
///
/// The expected values are written out from the rule's own sentence rather than produced by
/// calling the constructor a second time — an oracle built from the function under test
/// measures the function against itself.
#[test]
fn a_swapped_pair_renames_the_same_root() {
    let (a, b) = (2usize, 9usize);
    assert_eq!(
        QuadRoot::canonical([b, a], QuadRoot::Lo),
        ([a, b], QuadRoot::Hi)
    );
    assert_eq!(
        QuadRoot::canonical([b, a], QuadRoot::Hi),
        ([a, b], QuadRoot::Lo)
    );
}

/// **A pair already in order keeps its root** — the half its sibling above cannot see.
///
/// ★ A constructor that toggled *unconditionally* satisfies
/// [`a_swapped_pair_renames_the_same_root`] on its own: both spellings flip and the equality
/// survives the double negation. Only this test separates "toggle on swap" from "toggle".
#[test]
fn an_ordered_pair_keeps_its_root() {
    let (a, b) = (2usize, 9usize);
    assert_eq!(
        QuadRoot::canonical([a, b], QuadRoot::Lo),
        ([a, b], QuadRoot::Lo)
    );
    assert_eq!(
        QuadRoot::canonical([a, b], QuadRoot::Hi),
        ([a, b], QuadRoot::Hi)
    );
}

/// **A tangency keeps its name through a swap** — its two roots coincide, so the swap sends
/// the one point to itself.
///
/// ★ This is the case `Double` exists for: with a tangency stored as `Lo`, a remap that
/// toggles unconditionally brings the same point back named `Hi` — one point, two
/// names, which is exactly what canonicalization exists to prevent.
#[test]
fn a_tangency_keeps_its_name_through_a_swap() {
    let (a, b) = (2usize, 9usize);
    assert_eq!(
        QuadRoot::canonical([b, a], QuadRoot::Double),
        ([a, b], QuadRoot::Double)
    );
    assert_eq!(
        QuadRoot::canonical([a, b], QuadRoot::Double),
        ([a, b], QuadRoot::Double)
    );
}

/// ★★ **The geometry the three locks above are about.** Without it they are statements about
/// a convention someone chose; with it they are statements about where the points are.
///
/// `z = 0` and `x = 0` meet on the y-axis, and the cylinder about the x-axis with `r = 2` is
/// crossed there at `(0, ∓2, 0)`. Canonical normals give `ℓ = (0,0,1) × (1,0,0) = +y`, so
/// solved in that order the smaller parameter is `y = −2`. Solved the other way `ℓ` is `−y`,
/// and the algebra `lo′ = −hi`, `hi′ = −lo` says `s′[1]` names the *same point* as `s[0]`.
#[test]
fn swapping_the_planes_trades_the_two_roots_names() {
    use nacre_exact::Rat;
    use nacre_exact::quad::{CylinderMeet, branch_point_f64, plane_plane_cylinder};
    let r = Rat::from_int;
    let z0 = [r(0), r(0), r(1), r(0)];
    let x0 = [r(1), r(0), r(0), r(0)];
    let (origin, dir) = ([r(0), r(0), r(0)], [r(1), r(0), r(0)]);
    let solve = |p1: &[Rat; 4], p2: &[Rat; 4]| {
        let meet = plane_plane_cylinder(p1, p2, &origin, &dir, &nacre_exact::BigRat::from(r(4))); // radius 2, as r²
        match meet {
            Some(CylinderMeet::Pair { line, s }) => (line, s),
            other => panic!("the y-axis crosses this cylinder twice: {other:?}"),
        }
    };
    let (l_ab, s_ab) = solve(&z0, &x0);
    let (l_ba, s_ba) = solve(&x0, &z0);
    let near = |p: [f64; 3], q: [f64; 3]| (0..3).all(|i| (p[i] - q[i]).abs() < 1e-12);
    let at = |l: &nacre_exact::quad::MeetLine, s| branch_point_f64(l, s);
    assert!(near(at(&l_ab, &s_ab[0]), [0.0, -2.0, 0.0]), "lo is y = −2");
    assert!(near(at(&l_ab, &s_ab[1]), [0.0, 2.0, 0.0]), "hi is y = +2");
    assert!(
        near(at(&l_ba, &s_ba[1]), [0.0, -2.0, 0.0]),
        "hi′ is the old lo"
    );
    assert!(
        near(at(&l_ba, &s_ba[0]), [0.0, 2.0, 0.0]),
        "lo′ is the old hi"
    );
}

/// ★★★★ **`supersede_live` keeps the survivors in their original order — and nothing else in
/// this tree would notice if it did not.**
///
/// `nacre_step::to_step` exports the live set *in order*, the census reads the **arena** (it
/// never sees live order), and there is no golden STEP text anywhere. So a permutation here
/// would travel all the way out to the exported file unseen. The order is therefore locked at
/// the door itself, and again on the export side (`nacre-step`'s
/// `superseding_a_solid_leaves_the_export_order_alone`).
///
/// ⚠ The oracle is the **construction order**, not a second call to the same machinery —
/// comparing this against a hand-written `retain` would be one implementation checking itself.
#[test]
fn supersede_live_preserves_the_order_of_the_survivors() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
    let c = m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
    assert_eq!(
        m.live_solids(),
        [a, b, c].as_slice(),
        "the fixture's premise"
    );

    m.supersede_live(&[b]);
    assert_eq!(
        m.live_solids(),
        [a, c].as_slice(),
        "the middle solid goes and the order stays"
    );

    // It drops sets, not only singletons — that is what the 19 `retain` call sites ask for.
    m.supersede_live(&[a, c]);
    assert!(m.live_solids().is_empty());
}

/// `restore_live` is the rollback half: what a rejected operation puts back.
#[test]
fn restore_live_puts_the_snapshot_back() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
    let snapshot = m.live_solids().to_vec();

    m.supersede_live(&[a]);
    assert_eq!(m.live_solids(), [b].as_slice());

    m.restore_live(snapshot);
    assert_eq!(m.live_solids(), [a, b].as_slice(), "order comes back too");
}
