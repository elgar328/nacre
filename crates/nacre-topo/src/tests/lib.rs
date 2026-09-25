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

    // The f64 defect this guards against, shown on the planes themselves rather than
    // through the model, since the model holds only one of them. `Plane` keeps an
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
/// Controls, both load-bearing: a carrier turned a **quarter** folds and meets in the world while
/// one turned **37°** has no rational world name and still declines, and a corner the
/// shared-chain door can answer keeps answering **through that door** (`leaf = Some`) — which
/// is what says the new road is a fallback and not a re-spelling of what already worked.
#[test]
fn a_corner_of_two_translation_chains_solves_in_the_world() {
    use nacre_exact::{Angle, Axis, Rat};
    let mut m = Model::new();
    let r = Rat::from_int;
    // Three axis planes through the origin, pushed as their own statements.
    let plane = |m: &mut Model, n: [f64; 3], pts: [[i128; 3]; 3], sense: Orientation| {
        m.push_plane(
            nacre_geom::Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array(n),
            )
            .expect("unit normal"),
            pts.map(|p| p.map(r)),
            None,
            sense,
        )
        .0
    };
    let fwd = Orientation::Forward;
    let px = plane(
        &mut m,
        [1.0, 0.0, 0.0],
        [[0, 0, 0], [0, 1, 0], [0, 0, 1]],
        fwd,
    );
    // `x̂ × ẑ = −ŷ`: these points span the other way from the `+y` cache.
    let py = plane(
        &mut m,
        [0.0, 1.0, 0.0],
        [[0, 0, 0], [1, 0, 0], [0, 0, 1]],
        fwd.flipped(),
    );
    let pz = plane(
        &mut m,
        [0.0, 0.0, 1.0],
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
        fwd,
    );
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
                 at: [f64; 3],
                 sense: Orientation| {
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
        m.push_plane(world, pts, Some(leaf), sense).0
    };
    // x = 0 moved by (2,0,0) → the world plane x = 2; y = 0 moved by (0,3,0) → y = 3.
    let mx = moved(&mut m, px, t1, [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], fwd);
    // `py` interned onto the seeded y = 0, whose stored points `(0, ẑ, x̂)` span `+y` — the sense
    // is against the points the arena holds, not the ones this fixture wrote.
    let my = moved(&mut m, py, t2, [0.0, 1.0, 0.0], [0.0, 3.0, 0.0], fwd);
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

    // ① A quarter turn folds, so a turned carrier has a world name and the road meets it: x = 0
    // turned a quarter about z is y = 0, which meets the translated x = 2 wall and z = 0 at
    // (2, 0, 0).
    let quarter = m.push_motion(
        Motion::Rotate {
            axis: Axis::Z,
            pivot: [r(0); 3],
            angle: Angle::from_deg(r(90)).expect("angle"),
        },
        None,
    );
    let turned = moved(&mut m, px, quarter, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], fwd);
    let v_turned = m.push_vertex(
        Vertex::ThreePlane([turned, mx, pz]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.0, 0.0, 0.0]),
        },
    );
    let (p, frame) = m
        .vertex_meet(v_turned)
        .expect("a quarter turn folds into the world");
    assert_eq!(frame, None, "a folded carrier answers in the world");
    assert_eq!(
        p.narrow().map(|c| c.map(|x| x.to_f64())),
        Some([2.0, 0.0, 0.0]),
        "the turned wall meets the moved wall and the floor where the turn put it"
    );

    // ①′ The control: a turn off the quarters has irrational cos and sin, does not fold, and the
    // same corner under it still declines — what says the answer above is the fold's.
    let spin = m.push_motion(
        Motion::Rotate {
            axis: Axis::Z,
            pivot: [r(0); 3],
            angle: Angle::from_deg(r(37)).expect("angle"),
        },
        None,
    );
    let (c37, s37) = (37f64.to_radians().cos(), 37f64.to_radians().sin());
    let spun = moved(&mut m, px, spin, [c37, s37, 0.0], [0.0, 0.0, 0.0], fwd);
    let v_spun = m.push_vertex(
        Vertex::ThreePlane([spun, mx, pz]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.0, -2.0 * c37 / s37, 0.0]),
        },
    );
    assert!(
        m.vertex_meet(v_spun).is_none(),
        "a turn the rationals cannot state has no world name to meet with"
    );

    // ② The shared-chain door still answers through itself: same leaf on both moved carriers.
    let mx2 = moved(&mut m, px, t1, [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], fwd);
    let my2 = moved(&mut m, py, t1, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], fwd);
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

/// ★ **Two doors, one fact.** The three-vertex door `through_meets` calls `vertex_meet` per
/// vertex; the two must agree point for point — otherwise there is a second spelling of the
/// solve, which is this repo's dominant defect shape.
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
/// ★ `Model::surface`, `vertex_point` and `edge_curve` index a `Vec`, so without the guard a
/// handle from another model reaches all three and answers with the wrong cell. The guard is
/// `cfg(debug_assertions)`, so the
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

/// ★★★★ **The surface cache is writable, and the truth is not.**
///
/// The arena holds the truth, so the realization beside it can be rewritten in place — a
/// refinement pass has somewhere to write. This test is the warrant: it rewrites a cache entry in
/// place and reads it back through [`Model::surface_cache`], while the truth the same handle names
/// is untouched.
///
/// ★ It writes through the private field on purpose — that is the door a surface refinement
/// pass inside this crate would take. Outside, there is no door at all,
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
/// come out even where the derivation's narrow route cannot reach — that last one is where a
/// producer would lose a name, and losing a name closes the exact road for everything built on
/// that plane.
///
/// ★ There is no "wrong name" case to test. Coefficients are not something a producer can hand
/// in, so a plane cannot be stated twice — a disagreement is unspellable.
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
    let (ok, _) = m.push_plane(pl(0.0), pts, None, Orientation::Forward);
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
    // ★ Points chosen to overflow an **agreement check** (`c · p`) do not make the derivation
    // give up: those are different products, and the derivation goes through.
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
    let (wide, _) = m.push_plane(pl(2.0), wide_pts, None, Orientation::Forward);
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

    // ★ (A plane with no points cannot be pushed at all; the type is the assertion.)
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
    let (first, _) = m.push_plane(pl, [a, b, c], None, Orientation::Forward);
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
    let (second, _) = m.push_plane(pl, [a, b, c], None, Orientation::Forward);
    assert_eq!(first, second, "one wide plane, one handle");
    let (permuted, _) = m.push_plane(pl, [b, c, a], None, Orientation::Forward);
    assert_eq!(
        first, permuted,
        "two spellings of one wide plane must intern"
    );
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
        let (h, _) = m.push_plane(cache, pts, None, Orientation::Forward);
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
    // The lateral alone, stated through the kernel's own door: `+z` from the origin, `r = 2`.
    let mut m = Model::new();
    let r = |x: i128| nacre_exact::Rat::from_int(x);
    m.push_cylinder(
        nacre_geom::Cylinder::from_axis(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            2.0,
        )
        .expect("a cylinder"),
        CylinderDef::new(
            [r(0), r(0), r(0)],
            [r(0), r(0), r(1)],
            [r(1), r(0), r(0)],
            nacre_exact::BigRat::square_of(r(2)),
        )
        .expect("a statement"),
        None,
    );
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
    let (h, _) = m.push_plane_through(cache, vs, None, Orientation::Forward);
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
    let (through, _) = m.push_plane_through(cache, vs, None, Orientation::Forward);
    let r = |n: i128, d: i128| nacre_exact::Rat::new(n, d).unwrap();
    let (known, _) = m.push_plane(
        cache,
        [
            [r(1, 1), r(0, 1), r(0, 1)],
            [r(0, 1), r(1, 1), r(0, 1)],
            [r(0, 1), r(0, 1), r(1, 1)],
        ],
        None,
        Orientation::Forward,
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
    let (h, _) = m.push_plane_through(cache, vs, None, Orientation::Forward);
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
    let (h, _) = m.push_plane_through(cache, vs, None, Orientation::Forward);
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
    let (h2, _) = m2.push_plane_through(cache2, vs2, None, Orientation::Forward);
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
    m.push_plane(cache, pts, None, Orientation::Forward).0
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

/// ★ **A plane's cache faces the way its truth says** — the push turns a cache stated the other
/// way, exactly (both normals negated, the anchor kept), and leaves an agreeing one bit for bit.
///
/// Through the planting door, because the named door asserts the two agree before anything is
/// stored — so a disagreeing cache can only arrive where nothing asserts, which is exactly where
/// the turn has to hold on its own.
#[test]
fn a_plane_cache_follows_the_sense_its_truth_states() {
    let r = nacre_exact::Rat::from_int;
    let pts = [[r(0), r(0), r(2)], [r(1), r(0), r(2)], [r(0), r(1), r(2)]]; // spans +z
    let up = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 2.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 3.0]),
    )
    .expect("a plane");
    let bits = |p: &nacre_geom::Plane| {
        let c = p.coefficients();
        (c.map(f64::to_bits), p.origin().as_array().map(f64::to_bits))
    };

    let mut m = Model::new();
    // Agreeing: `Forward` against a `+z` cache — nothing moves.
    let h = m.push_plane_unregistered(up, pts, Orientation::Forward);
    let nacre_geom::Surface::Plane(kept) = *m.surface_cache(h) else {
        unreachable!()
    };
    assert_eq!(
        bits(&kept),
        bits(&up),
        "an agreeing cache is kept bit for bit"
    );

    // Disagreeing: `Reversed` against the same `+z` cache — the cache turns to `−z`.
    let h = m.push_plane_unregistered(up, pts, Orientation::Reversed);
    let nacre_geom::Surface::Plane(turned) = *m.surface_cache(h) else {
        unreachable!()
    };
    assert_eq!(turned, up.reversed(), "the cache follows the truth's sense");
    assert_eq!(turned.origin(), up.origin(), "only the sense moved");
    assert!(!m.align_cache_sense(h), "and asking again turns nothing");
}

/// **Which way a plane's own name faces is read off its truth** — `sense` against its points'
/// turn — and it is the answer the aligned cache gives: the three seeds (cache down each axis, name
/// up it), a `Known` plane stated either way round and with either sense, a `Through` datum, and a
/// wide name. The oracle is the cache normal against the name, an independent `f64` reading of
/// two normals of one plane.
#[test]
fn a_planes_name_sense_is_read_from_its_truth() {
    use nacre_exact::Rat;
    let expected = |m: &Model, h: Handle<Surface>| {
        let nacre_geom::Surface::Plane(p) = m.surface_cache(h) else {
            unreachable!()
        };
        let c = m.surface_name.get(&h).expect("a named plane").coeff_ints();
        let dot: f64 = (0..3)
            .map(|k| p.normal().as_array()[k] * c[k].to_string().parse::<f64>().expect("finite"))
            .sum();
        if dot > 0.0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    };
    let mut m = Model::new();
    for i in 0..3 {
        let h = m.surface_handle_at(i).expect("a seed");
        assert_eq!(
            m.plane_name_sense(h),
            Some(Orientation::Reversed),
            "seed {i}"
        );
        assert_eq!(m.plane_name_sense(h), Some(expected(&m, h)), "seed {i}");
    }
    let r = Rat::from_int;
    let turning_up = |z: i128| [[r(0), r(0), r(z)], [r(1), r(0), r(z)], [r(0), r(1), r(z)]];
    let turning_down = |z: i128| [[r(0), r(0), r(z)], [r(0), r(1), r(z)], [r(1), r(0), r(z)]];
    let cache = |z: f64, nz: f64| {
        nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, z]),
            Vector3::from_array([0.0, 0.0, nz]),
        )
        .unwrap()
    };
    for (pts, sense, nz, want) in [
        (
            turning_up(3),
            Orientation::Forward,
            1.0,
            Orientation::Forward,
        ),
        (
            turning_down(5),
            Orientation::Forward,
            -1.0,
            Orientation::Reversed,
        ),
        (
            turning_up(7),
            Orientation::Reversed,
            -1.0,
            Orientation::Reversed,
        ),
        (
            turning_down(9),
            Orientation::Reversed,
            1.0,
            Orientation::Forward,
        ),
    ] {
        let z = pts[0][2].to_f64();
        let (h, flipped) = m.push_plane(cache(z, nz), pts, None, sense);
        assert!(!flipped);
        assert_eq!(m.plane_name_sense(h), Some(want), "z = {z}");
        assert_eq!(m.plane_name_sense(h), Some(expected(&m, h)), "z = {z}");
    }
    // A wide name, cache from the points' own `f64` turn.
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let (b1, b2) = ((1i128 << 90) + 1, (1i128 << 90) + 3);
    let pts = [
        [q(b1, 3), q(b2, 7), q(0, 1)],
        [q(-b2, 5), q(b1, 11), q(0, 1)],
        [q(1, 13), q(1, 17), q(1, 19)],
    ];
    let f = |p: [Rat; 3]| Point3::from_array(p.map(|x| x.to_f64()));
    let wide_cache =
        nacre_geom::Plane::through_points(f(pts[0]), f(pts[1]), f(pts[2])).expect("a plane");
    let (h, _) = m.push_plane(wide_cache, pts, None, Orientation::Forward);
    assert!(m.surface_name.get(&h).is_some_and(|n| n.narrow().is_none()));
    assert_eq!(m.plane_name_sense(h), Some(expected(&m, h)), "wide");
    // A `Through` datum over three box corners.
    let (mut m, vs) = box_corner_vertices();
    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .unwrap();
    let (h, _) = m.push_plane_through(cache, vs, None, Orientation::Forward);
    assert_eq!(m.plane_name_sense(h), Some(expected(&m, h)), "through");
}

/// **Which way a plane's world name faces follows its chain** — the name's own sense, times the
/// chain's determinant (a reflection reverses the points' turn), times the sign canonicalizing the
/// carried name may take. Planes stated both ways round and with both senses, under reflections,
/// quarter turns and a translation. The oracle is independent of that algebra: the plane's points
/// carried to the world ([`Model::chain_point_rat`]) and crossed in `f64`, times the sense,
/// against the world name.
#[test]
fn a_planes_world_name_sense_follows_its_chain() {
    use nacre_exact::{Angle, Axis, Rat};
    let r = Rat::from_int;
    let mut m = Model::new();
    let motions = [
        Motion::Mirror {
            axis: Axis::Z,
            offset: r(0),
        },
        Motion::Mirror {
            axis: Axis::X,
            offset: r(1),
        },
        Motion::Rotate {
            axis: Axis::X,
            pivot: [r(0); 3],
            angle: Angle::from_deg(r(90)).expect("angle"),
        },
        Motion::Rotate {
            axis: Axis::Y,
            pivot: [r(1), r(0), r(0)],
            angle: Angle::from_deg(r(270)).expect("angle"),
        },
        Motion::Translate {
            offset: [r(1), r(-2), r(3)],
        },
    ];
    let statements = [
        [[r(0), r(0), r(3)], [r(1), r(0), r(3)], [r(0), r(1), r(3)]],
        [[r(0), r(0), r(3)], [r(0), r(1), r(3)], [r(1), r(0), r(3)]],
        // Tilted: normal (1, 2, 3) and its reverse.
        [[r(6), r(0), r(0)], [r(0), r(3), r(0)], [r(0), r(0), r(2)]],
        [[r(6), r(0), r(0)], [r(0), r(0), r(2)], [r(0), r(3), r(0)]],
    ];
    let mut checked = 0;
    for motion in motions {
        let leaf = m.push_motion(motion, None);
        for pts in statements {
            for sense in [Orientation::Forward, Orientation::Reversed] {
                let world: Vec<Point3> = pts
                    .iter()
                    .map(|&p| {
                        let q = m.chain_point_rat(leaf, p).expect("a folding chain");
                        Point3::from_array(q.map(|x| x.to_f64()))
                    })
                    .collect();
                let turn = (world[1] - world[0]).cross(world[2] - world[0]);
                let facing = turn * f64::from(sense.sign());
                let cache =
                    nacre_geom::Plane::from_point_normal(world[0], facing).expect("a plane");
                let (h, flipped) = m.push_plane(cache, pts, Some(leaf), sense);
                let name = m
                    .world_plane_name(h)
                    .expect("a folding chain names the world");
                let c = name.narrow().expect("narrow");
                let dot: f64 = (0..3).map(|k| facing.as_array()[k] * c[k].to_f64()).sum();
                // An interning hit keeps the survivor's sense; the new statement faces the other
                // way from it exactly when `flipped`.
                let along = (dot > 0.0) != flipped;
                let want = if along {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                assert_eq!(
                    m.world_plane_name_sense(h),
                    Some(want),
                    "{motion:?} {pts:?} {sense:?}"
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 40);
    // Unmoved, the world name is the name.
    let (h, _) = m.push_plane(
        nacre_geom::Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 5.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap(),
        [[r(0), r(0), r(5)], [r(1), r(0), r(5)], [r(0), r(1), r(5)]],
        None,
        Orientation::Forward,
    );
    assert_eq!(m.world_plane_name_sense(h), m.plane_name_sense(h));
}

/// **A chain pushed node by node folds to its nodes applied root first** — the one lock on the
/// order `push_motion` composes in (the census has no chain with two non-identity linear parts).
/// The oracle is the existing per-motion arithmetic, stepped; a turn off the quarters is in the
/// mix, and both sides must then answer `None`.
mod chain_folds {
    use super::*;
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};

    fn r(n: i128, d: i128) -> Rat {
        Rat::new(n, d).expect("a rational")
    }

    fn small() -> impl Strategy<Value = Rat> {
        (-12i128..=12, 1i128..=6).prop_map(|(n, d)| r(n, d))
    }

    fn point() -> impl Strategy<Value = [Rat; 3]> {
        proptest::array::uniform3(small())
    }

    fn axis_of(k: u8) -> Axis {
        [Axis::X, Axis::Y, Axis::Z][k as usize % 3]
    }

    fn motion() -> impl Strategy<Value = Motion> {
        prop_oneof![
            4 => point().prop_map(|offset| Motion::Translate { offset }),
            4 => (0u8..3, point(), 0i128..4).prop_map(|(a, pivot, q)| Motion::Rotate {
                axis: axis_of(a),
                pivot,
                angle: Angle::from_deg(Rat::from_int(90 * q)).expect("a quarter"),
            }),
            3 => (0u8..3, small()).prop_map(|(a, offset)| Motion::Mirror { axis: axis_of(a), offset }),
            1 => (0u8..3, point()).prop_map(|(a, pivot)| Motion::Rotate {
                axis: axis_of(a),
                pivot,
                angle: Angle::from_deg(Rat::from_int(37)).expect("an angle"),
            }),
        ]
    }

    fn step_point(m: Motion, p: [Rat; 3]) -> Option<[Rat; 3]> {
        match m {
            Motion::Translate { offset } => Isometry::translation(offset).point_rat(p),
            Motion::Rotate { axis, pivot, angle } => {
                Isometry::rotation(Rotation { axis, pivot, angle }).point_rat(p)
            }
            Motion::Mirror { axis, offset } => nacre_exact::mirror_point_rat(p, axis, offset),
            Motion::Frame { .. } => None,
        }
    }

    fn step_plane(m: Motion, c: [Rat; 4]) -> Option<[Rat; 4]> {
        match m {
            Motion::Translate { offset } => Isometry::translation(offset).plane_coeffs(c),
            Motion::Rotate { axis, pivot, angle } => {
                Isometry::rotation(Rotation { axis, pivot, angle }).plane_coeffs(c)
            }
            Motion::Mirror { axis, offset } => nacre_exact::mirror_plane_coeffs(c, axis, offset),
            Motion::Frame { .. } => None,
        }
    }

    proptest! {
        #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
            proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
        ))]
        #[test]
        fn a_pushed_chain_folds_to_its_nodes_root_first(
            chain in proptest::collection::vec(motion(), 1..6),
            p in point(),
            c in proptest::array::uniform4(small()),
        ) {
            let mut m = Model::new();
            let mut leaf = None;
            for &node in &chain {
                leaf = Some(m.push_motion(node, leaf));
            }
            let leaf = leaf.expect("a nonempty chain");
            let stepped = chain.iter().try_fold(p, |q, &n| step_point(n, q));
            prop_assert_eq!(m.chain_point_rat(leaf, p), stepped);
            prop_assume!(c[..3].iter().any(|x| *x != Rat::from_int(0)));
            let stepped_plane = chain.iter().try_fold(c, |q, &n| step_plane(n, q));
            prop_assert_eq!(m.chain_plane_coeffs(leaf, c), stepped_plane);
        }
    }
}

/// ★★ **An edge's kind is the truth's, not the cache's tolerance.** Two planes a hair off the
/// relations a plane can have with a `+z` cylinder — one tilted `1e-10` from ⊥ the axis, one
/// tilted `1e-10` from ∥ it — are both **oblique** in their exact statements: the section is an
/// ellipse, which no edge here can carry. Read off `f64` normals under a tolerance, the first
/// would pass for a circle and the second for a ruling line.
#[test]
fn a_plane_a_hair_off_either_relation_is_oblique() {
    let mut m = Model::new();
    let r = |x: i128| Rat::from_int(x);
    let d = |x: f64| Rat::from_decimal(x).expect("a short decimal");
    let cyl = m.push_cylinder(
        Cylinder::from_axis(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            1.0,
        )
        .expect("a cylinder"),
        CylinderDef::new(
            [r(0), r(0), r(0)],
            [r(0), r(0), r(1)],
            [r(1), r(0), r(0)],
            nacre_exact::BigRat::square_of(r(1)),
        )
        .expect("a statement"),
        None,
    );
    for (name, pts) in [
        (
            "a hair off ⊥",
            [
                [r(0), r(0), r(0)],
                [r(1), r(0), d(1e-10)],
                [r(0), r(1), r(0)],
            ],
        ),
        (
            "a hair off ∥",
            [
                [r(0), r(0), r(0)],
                [r(1), r(0), r(0)],
                [r(0), d(1e-10), r(1)],
            ],
        ),
    ] {
        let f = |p: [Rat; 3]| Point3::from_array(p.map(|x| x.to_f64()));
        let cache = Plane::through_points(f(pts[0]), f(pts[1]), f(pts[2])).expect("a plane");
        let (plane, _) = m.push_plane(cache, pts, None, Orientation::Forward);
        let v = m.push_vertex(
            Vertex::OnSeam([cyl, plane]),
            PointCache::Unrealized {
                coord: Point3::from_array([1.0, 0.0, 0.0]),
            },
        );
        assert_eq!(
            m.derive_edge_curve([plane, cyl], [v, v]),
            Err(EdgeDecline::Oblique),
            "{name}: an oblique section is no circle and no ruling"
        );
    }
}

/// ★★ **A plane the cylinder's turn leaves in place is read back through that turn.** An `x`
/// cylinder turned 30° about `z` has an irrational world axis, so no world statement; the world
/// plane `z = 1` is square to the turn's axis. Carried forward, the plane reaches the world and the
/// cylinder's axis stops at the turn — but the turn fixes `ẑ`, so the plane's normal carried
/// **back** through it is `ẑ` again, square to the axis `x̂`: the section is a ruling, a line.
#[test]
fn a_plane_the_turn_fixes_meets_a_turned_cylinder_along_a_ruling() {
    let mut m = Model::new();
    let r = |x: i128| Rat::from_int(x);
    let turn = m.push_motion(
        Motion::Rotate {
            axis: Axis::Z,
            pivot: [r(0), r(0), r(0)],
            angle: Angle::from_deg(r(30)).expect("an angle"),
        },
        None,
    );
    let (c, s) = (0.75f64.sqrt(), 0.5);
    let cyl = m.push_cylinder(
        Cylinder::from_axis(
            Point3::origin(),
            Vector3::from_array([c, s, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            1.0,
        )
        .expect("a cylinder"),
        CylinderDef::new(
            [r(0), r(0), r(0)],
            [r(1), r(0), r(0)],
            [r(0), r(0), r(1)],
            nacre_exact::BigRat::square_of(r(1)),
        )
        .expect("a statement"),
        Some(turn),
    );
    let plane_pts = [[r(0), r(0), r(1)], [r(1), r(0), r(1)], [r(0), r(1), r(1)]];
    let cache = Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .expect("a plane");
    let (plane, _) = m.push_plane(cache, plane_pts, None, Orientation::Forward);
    let mut at = |x: f64| {
        m.push_vertex(
            Vertex::OnSeam([cyl, plane]),
            PointCache::Unrealized {
                coord: Point3::from_array([x * c, x * s, 1.0]),
            },
        )
    };
    let (a, b) = (at(0.0), at(2.0));
    assert!(
        matches!(
            m.derive_edge_curve([plane, cyl], [a, b]),
            Ok(Curve::Line(_))
        ),
        "the plane square to the turn's axis meets the turned `x` cylinder along a ruling"
    );
}
