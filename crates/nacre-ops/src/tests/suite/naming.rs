//! What a result is called: vertex plane triples, which piece a void belongs to, replay
//! determinism,
//! plane classes and interning, canonical triples.

use super::*;

/// **`Origin` no longer tells result faces apart.** The arrangement names every vertex
/// it emits by the three planes meeting there, so an operand corner the cut never touched
/// comes back as `Discovered`, exactly like a seam vertex. Nothing carries over as
/// `Constructed` (the arrangement builds no vertex from an original handle).
///
/// This is a contract, not a curiosity: `pipeline.rs`'s island test selected a face by
/// "all its vertices are `Discovered`", which was unique under the old engine and is
/// now true of *every* face. It flipped the wrong face and only the last assertion
/// noticed. Selecting a face by provenance is what this locks out.
///
/// The subject is `Cut(l, stub)` — the **holed** result, so `face_half_edges` walks
/// `inner` rings too (`count_discovered` walks only `outer` and would miss them).
///
/// **Unrotated only.** Rotating a boolean result re-marks these vertices `Rotated` over
/// a `Discovered` base — that is `transform_rotate_boolean_result_keeps_discovered_base`,
/// and this lock must not be read as contradicting it.
#[test]
fn an_unrotated_boolean_names_every_vertex_by_its_plane_triple() {
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();

    let mut seen = std::collections::HashSet::new();
    let mut holed = 0;
    for sh in solid_shell_handles(&m, r) {
        for &fh in &m.shell(sh).faces {
            let face = m.face(fh);
            holed += usize::from(!face.inner.is_empty());
            for he in face_half_edges(face) {
                for vh in m.edge(he.edge).vertices {
                    if !seen.insert(vh) {
                        continue;
                    }
                    assert!(
                        matches!(*m.vertex(vh), Vertex::ThreePlane(_))
                            && matches!(m.vertex_cache(vh), nacre_topo::PointCache::Bounded { .. }),
                        "vertex {:?} is {:?}, not a realized plane triple",
                        m.vertex_point(vh).as_array(),
                        *m.vertex(vh)
                    );
                }
            }
        }
    }
    assert_eq!(holed, 1, "the blind dimple leaves exactly one holed face");
    assert_eq!(seen.len(), 20, "the L's 12 corners + the dimple's 8");
}

/// A sever that also leaves a surviving cavity: a hollow box whose void sits to one side,
/// cut by a slab that severs it without touching the void. The x<2 piece keeps the void as a
/// cavity, the x>2 piece is solid — two outward shells *and* one inward. `point_in_component`
/// assigns the void to the x<2 piece that nests it (a containment test), rather than rejecting.
#[test]
fn severed_with_cavity_assigns_the_void_to_its_piece() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    // Void near the x-low side (1×2×2 = 4), clear of the x=2 cut.
    let inner = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 2.5, 2.5]),
    );
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(m.solid(hollow).cavities.len(), 1);
    // A slab spanning full y,z, thin in x at x∈[2,2.2] — severs into x<2 (holds the void, vol
    // 2·3·3 − 4 = 14) and x>2 (solid, vol 0.8·3·3 = 7.2).
    let slab = m.add_cuboid(
        Point3::from_array([2.0, -1.0, -1.0]),
        Point3::from_array([2.2, 4.0, 4.0]),
    );
    let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
    assert_eq!(
        solids.len(),
        2,
        "the slab severs the hollow box into two pieces"
    );
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // Exactly one piece owns the void; volumes match the hand calculation.
    let with_cav: Vec<_> = solids
        .iter()
        .filter(|&&s| !m.solid(s).cavities.is_empty())
        .collect();
    assert_eq!(
        with_cav.len(),
        1,
        "the void is assigned to exactly one piece"
    );
    let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
    let hollow_piece = *with_cav[0];
    assert!(
        (vol(hollow_piece) - 14.0).abs() < 1e-9,
        "hollow piece {}",
        vol(hollow_piece)
    );
    let total: f64 = solids.iter().map(|&s| vol(s)).sum();
    assert!((total - 21.2).abs() < 1e-9, "total {total}");
}

/// The adjacent case: a cut that passes *through* the void opens it — the void wall becomes
/// exterior boundary, so no cavity survives. Handled by the plain sever path (no cavity to
/// assign), not the containment code, but pinned so a regression there is caught.
#[test]
fn a_cut_through_the_void_leaves_no_cavity() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3])); // void 1³
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(m.solid(hollow).cavities.len(), 1);
    // Slab x∈[1.4,1.6] passes through the void (x∈[1,2]) → severs AND opens the void.
    let slab = m.add_cuboid(
        Point3::from_array([1.4, -1.0, -1.0]),
        Point3::from_array([1.6, 4.0, 4.0]),
    );
    let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // The void is opened, so neither piece keeps a cavity; material = 26 − (1.8 − 0.2) = 24.4.
    let total_cavities: usize = solids.iter().map(|&s| m.solid(s).cavities.len()).sum();
    assert_eq!(
        total_cavities, 0,
        "the cut opened the void — no surviving cavity"
    );
    let total: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((total - 24.4).abs() < 1e-9, "total {total}");
}

/// Nested cavities: a hollow box A ([0,6]³ − [1,5]³ void) with a smaller hollow box B
/// ([2,4]³ − [2.5,3.5]³ void) floating inside A's void. `Fuse(A,B)` is one arrangement with
/// four components (two materials, two voids); B's void is contained by **both** A's outer
/// shell and B's own, so the containment assignment must pick the **innermost** (B), not A.
/// The result is two solids, each keeping its own void (A: 216−64 = 152, B: 8−1 = 7).
#[test]
fn a_void_nested_in_a_floating_island_goes_to_the_inner_solid() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([6.0; 3]));
    let void = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([5.0; 3]));
    let a = boolean_one(&mut m, BoolKind::Cut, big, void).unwrap();
    m.rebuild_adjacency();
    let bbig = m.add_cuboid(Point3::from_array([2.0; 3]), Point3::from_array([4.0; 3]));
    let bvoid = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
    let b = boolean_one(&mut m, BoolKind::Cut, bbig, bvoid).unwrap();
    m.rebuild_adjacency();
    let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // Each solid keeps exactly one void — B's void was assigned to B (innermost), not A.
    let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
    for &s in &solids {
        assert_eq!(
            m.solid(s).cavities.len(),
            1,
            "each piece keeps its own void"
        );
    }
    let mut vols: Vec<f64> = solids.iter().map(|&s| vol(s)).collect();
    vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
    assert!((vols[0] - 7.0).abs() < 1e-9, "inner {}", vols[0]);
    assert!((vols[1] - 152.0).abs() < 1e-9, "outer {}", vols[1]);
}

#[test]
fn replay_is_deterministic() {
    let log = vec![extrude_log_op(square(), 1.0)];
    let m1 = replay(&log).unwrap();
    let m2 = replay(&log).unwrap();
    let pts = |m: &Model| {
        (0..m.vertex_count() as u32)
            .filter_map(|i| m.vertex_handle_at(i))
            .map(|h| (h, m.vertex(h)))
            .map(|(vh, _)| m.vertex_point(vh).as_array())
            .collect::<Vec<_>>()
    };
    assert_eq!(pts(&m1), pts(&m2));
    assert_eq!(m1.edge_count(), m2.edge_count());
    assert_eq!(m1.face_count(), m2.face_count());
}

#[test]
fn two_extrudes_make_two_solids() {
    let far = SketchPlane::world_xy().with_origin(Point3::from_array([5.0, 0.0, 0.0]));
    // ★ The second plane is not a seed, so the log has to state it — and a datum's handle is
    // not known until the datum runs. So the log is assembled against a **scratch model built
    // the same way**: a frame that names surface *N* there names surface *N* in the replay,
    // because `replay` re-anchors indices. That is what a recording
    // session does, and R is what makes it sound.
    let mut scratch = Model::new();
    let first = extrude_op(&scratch, square(), 1.0);
    apply(&mut scratch, &first).unwrap();
    let far_frame = datum_frame(&mut scratch, far);
    let log = vec![
        first,
        Operation::DatumPlane {
            def: DatumDef::Stated(far),
        },
        Operation::Extrude {
            frame: far_frame,
            profile: square(),
            dist: 1.0,
        },
    ];
    let m = replay(&log).unwrap();
    assert_eq!(m.solid_count(), 2);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **Chaining onto a fused boss.** The fuse leaves the base's `z=1` face a *ring* — a face with
/// a hole where the boss sits — and the second boolean cuts through both. Every plane class the
/// cut opens then meets that ring along the **hole's own edge**, which is the case that used to
/// label inconsistently: the ring's neighbouring vertices there point *into* the hole, so
/// reading the occupied side off a flank put the material on the wrong side of `W`. The side
/// now comes from the ring's travel ([`arrangement::run_body_above`]), and the run leaves as its own
/// homogeneous segment rather than being swallowed by the straddling stretch beside it.
///
/// Hand volume: `1 + 0.5·0.5·1` fused, less the cutter's `0.2·0.2` column over `z ∈ [0.5, 2]`
/// — `1.25 − 0.06 = 1.19`. The same shape is scored against OCCT by
/// `boss_fuse_then_cut_matches_occt`, but that oracle is `#[ignore]`d, so this is the copy that
/// runs on every `cargo test`.
#[test]
fn a_boss_fused_then_cut_through() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.25, 0.25, 1.0]),
        Point3::from_array([0.75, 0.75, 2.0]),
    );
    let bossed = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, bossed).unwrap().volume - 1.25).abs() < 1e-12,
        "the fused boss itself"
    );
    let cutter = m.add_cuboid(
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, bossed, cutter).expect("the chained cut");
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, r).unwrap().volume - 1.19).abs() < 1e-12,
        "base + boss less the drilled column: {}",
        nacre_props::mass_props(&m, r).unwrap().volume
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "a chained result is still a clean model"
    );
}

/// One geometric plane is one class **whatever the two faces' sizes**.
///
/// ★ **The reason changed, and that is the news.** This used to assert that the coefficient
/// test *could not* prove these coplanar — two walls of one plane at different face sizes have
/// un-normalized 4-vectors that are not exactly proportional — and that the coordinate branch
/// was what earned the merge. Surfaces are interned on their canonical rational coefficients
/// now, so the two walls are handed **one handle**, and the merge is a handle comparison
/// before any geometry is asked. The f64 non-proportionality is still real and still pinned,
/// on the planes themselves, in `nacre_topo`'s
/// `two_faces_of_one_plane_disagree_in_f64_and_agree_in_the_rationals`.
#[test]
fn one_plane_is_one_class_whatever_the_face_size() {
    let (dx, dy) = (1.628165457453874f64, 0.5f64);
    let (z0, h1, h2) = (0.11200046228159026f64, 0.5f64, 2.07926124157585f64);
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, z0]),
        Point3::from_array([dx, dy, z0 + h1]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, z0 + h1]),
        Point3::from_array([dx, dy, z0 + h1 + h2]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: _planes,
        plane_ix,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    // The two `+X` walls: same plane x = dx, different face sizes (heights h1 vs h2). Two
    // *faces*, so this searches the face table — the plane table holds one entry for both,
    // which is the property under test.
    let x_walls: Vec<usize> = (0..faces_tab.len())
        .filter(|&i| {
            faces_tab[i].plane().n_out.as_array() == [1.0, 0.0, 0.0]
                && (faces_tab[i].plane().tri[0].as_array()[0] - dx).abs() < 1e-12
        })
        .collect();
    assert_eq!(x_walls.len(), 2, "one wall from each box: {x_walls:?}");
    let (i, j) = (x_walls[0], x_walls[1]);
    assert_eq!(
        faces_tab[i].surf(),
        faces_tab[j].surf(),
        "two faces, one surface — interning collapsed them at construction"
    );
    assert!(
        crate::planes::test_judge(&faces_tab).planes_coplanar(i, j),
        "and the geometry agrees, so nothing rests on the handle alone"
    );
    assert_eq!(plane_ix[i], plane_ix[j], "so they are one plane-table row");
}

/// ★★★ **Solving a vertex's three planes lands on its coordinate.**
///
/// This is the premise the point types rest on: a point *is* the
/// meeting of three surfaces, and the stored coordinate is a rounded answer to that question.
/// Stage 0 measured it by reading the census; this asserts it on data the kernel itself built,
/// which is a different claim — the definitions have to be *right*, not merely present.
///
/// ★ **Split by origin, because one row proves nothing.** A `Discovered` vertex's coordinate
/// was produced by solving exactly this triple, so agreement there is an identity. The rows
/// that carry weight are `Constructed` and `Moved`, where the coordinate came from somewhere
/// else entirely — construction arithmetic, or a motion replayed on a base point.
#[test]
fn a_vertex_definition_solves_to_its_own_coordinate() {
    // ★ **Three solids, because one would not exercise three kinds of vertex.** The first
    // draft of this test used a fused-then-turned-then-cut solid and reported
    // `Constructed (0,0) Discovered (24,24) Moved (0,0)` with a worst error of exactly zero —
    // a boolean recomputes every vertex it emits, so the only row present was the tautological
    // one. A plain box keeps its constructed corners; a turned box keeps them as `Moved`.
    let mut m = Model::new();
    let plain = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let to_turn = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let OpOutput::Transform { solid: turned } = apply(
        &mut m,
        &Operation::Transform {
            solid: to_turn,
            isometry: nacre_exact::Isometry::rotation(nacre_exact::Rotation {
                axis: nacre_exact::Axis::Z,
                pivot: [nacre_exact::Rat::from_int(0); 3],
                angle: nacre_exact::Angle::from_deg(nacre_exact::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!("transform yields Transform output")
    };
    m.rebuild_adjacency();
    let a = m.add_cuboid(
        Point3::from_array([10.0, 0.0, 0.0]),
        Point3::from_array([12.0, 3.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([11.0, 1.0, 0.5]),
        Point3::from_array([14.0, 2.0, 2.5]),
    );
    m.rebuild_adjacency();
    let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse")[0];
    m.rebuild_adjacency();

    // Rows by **producer family**: a plainly built box, a rotated one, a boolean result.
    let mut counts = [(0usize, 0usize); 3]; // (with a three-plane definition, total)
    let mut worst = [0.0f64; 3];
    let mut diam = 0.0f64;
    for (kind, solid) in [plain, turned, fused].into_iter().enumerate() {
        let shell = m.solid(solid).outer;
        let mut seen: Vec<Handle<Vertex>> = Vec::new();
        for &fh in &m.shell(shell).faces {
            let face = m.face(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in m.edge(he.edge).vertices.iter() {
                        if !seen.contains(&vh) {
                            seen.push(vh);
                        }
                    }
                }
            }
        }
        for &vh in &seen {
            let coord = m.vertex_point(vh);
            for c in coord.as_array() {
                diam = diam.max(c.abs());
            }
            counts[kind].1 += 1;
            let Vertex::ThreePlane(planes) = *m.vertex(vh) else {
                continue; // a seam vertex names a curve, not a point — no solve to check
            };
            counts[kind].0 += 1;
            let coeffs = planes.map(|s| match m.surface_cache(s) {
                nacre_geom::Surface::Plane(p) => p.coefficients(),
                nacre_geom::Surface::Cylinder(_) => panic!("ThreePlane named a cylinder"),
            });
            let solved = solve_three_planes(coeffs).expect("three planes meeting at a point");
            let d = solved
                .iter()
                .zip(coord.as_array())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            worst[kind] = worst[kind].max(d);
        }
    }

    let limit = diam * 2f64.powi(-40);
    eprintln!(
        "[definition coverage] plain {:?} worst {:e} | turned {:?} worst {:e} | \
             fused {:?} worst {:e} | limit {:e}",
        counts[0], worst[0], counts[1], worst[1], counts[2], worst[2], limit
    );
    for (i, name) in ["plain", "turned", "fused"].iter().enumerate() {
        assert!(
            counts[i].1 > 0,
            "{name} is not exercised — the row proves nothing"
        );
        assert_eq!(
            counts[i].0, counts[i].1,
            "{name} vertices without a definition: {:?}",
            counts[i]
        );
        assert!(
            worst[i] <= limit,
            "{name}: a definition solved {:e} away from its own coordinate (limit {limit:e})",
            worst[i]
        );
    }
    // ★ The **fused** row is the tautological one — a boolean vertex's coordinate *came
    // from* solving this very triple, so a zero there says nothing. `plain` and `turned`
    // are the claims.
    assert_eq!(worst[2], 0.0, "a boolean vertex is its own solve, exactly");
}

/// The negative control for the assertion above: point a definition at the wrong plane and the
/// solve must land somewhere else. Without this, an agreement test passes on any model whose
/// planes happen to be near each other.
#[test]
fn a_wrong_plane_in_a_definition_is_caught() {
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    let shell = m.solid(s).outer;
    let surfaces: Vec<_> = m
        .shell(shell)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .collect();
    // Take a real corner's triple and swap one plane for another face of the same box.
    let corner = m
        .shell(shell)
        .faces
        .first()
        .map(|&fh| {
            let face = m.face(fh);
            m.edge(face.outer.half_edges[0].edge).vertices
        })
        .expect("a face with a loop")[0];
    let Vertex::ThreePlane(mut planes) = *m.vertex(corner) else {
        panic!("a constructed corner has a definition")
    };
    let good = solve_three_planes(planes.map(|h| match m.surface_cache(h) {
        nacre_geom::Surface::Plane(p) => p.coefficients(),
        nacre_geom::Surface::Cylinder(_) => unreachable!(),
    }))
    .expect("meets at a point");
    // Any face not already in the triple. Every one of them moves the point (or leaves the
    // three not meeting at all, which is just as good a refutation).
    let other = *surfaces
        .iter()
        .find(|h| !planes.contains(h))
        .expect("a fourth face");
    planes[0] = other;
    let bad = solve_three_planes(planes.map(|h| match m.surface_cache(h) {
        nacre_geom::Surface::Plane(p) => p.coefficients(),
        nacre_geom::Surface::Cylinder(_) => unreachable!(),
    }));
    let moved = bad.is_none_or(|bad| {
        good.iter()
            .zip(bad)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max)
            > 0.1
    });
    assert!(
        moved,
        "swapping a plane must move the solved point, or the assertion above is vacuous"
    );
}

/// Three planes by Cramer, or `None` when they do not meet in a point.
fn solve_three_planes(p: [[f64; 4]; 3]) -> Option<[f64; 3]> {
    let n = |i: usize| [p[i][0], p[i][1], p[i][2]];
    let (a, b, c) = (n(0), n(1), n(2));
    let det3 = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det3([a, b, c]);
    if d == 0.0 {
        return None;
    }
    let rhs = [-p[0][3], -p[1][3], -p[2][3]];
    Some(core::array::from_fn(|j| {
        let mut m = [a, b, c];
        for (row, r) in m.iter_mut().zip(rhs) {
            row[j] = r;
        }
        det3(m) / d
    }))
}

/// ★ **A collinear midpoint is dissolved at construction, so the prism it builds is its
/// clean twin's, bit for bit — and every corner keeps a three-plane definition.**
///
/// Undissolved, this profile builds *seven* faces whose split bottom edge's two walls intern to
/// one surface handle — a vertex named `[S, S, cap]`, a line and not a point, `definition:
/// None`. The constructor deletes the flat corner (a lossless normalization: the shape is
/// identical), so that population cannot reach the topology at all.
#[test]
fn a_collinear_midpoint_profile_builds_its_clean_twin_bit_for_bit() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let build = |profile: Profile2d| {
        let mut m = Model::new();
        let __w3 = SketchFrame::world(&m, Axis::Z);
        apply(
            &mut m,
            &Operation::Extrude {
                frame: __w3,
                profile,
                dist: 1.0,
            },
        )
        .expect("extrude");
        m
    };
    // A unit square whose bottom edge carries a redundant midpoint — and the square itself.
    let split = build(
        Profile2d::polygon(vec![
            p(0.0, 0.0),
            p(0.5, 0.0), // collinear with its neighbours — dissolved at construction
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
        ])
        .unwrap(),
    );
    let clean = build(
        Profile2d::polygon(vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]).unwrap(),
    );
    // Bit-for-bit the same model: same counts, same coordinates, same surface wiring.
    assert_eq!(split.face_count(), clean.face_count(), "6 faces, not 7");
    assert_eq!(split.vertex_count(), clean.vertex_count());
    for i in 0..split.vertex_count() as u32 {
        let (sv, cv) = (
            split.vertex_handle_at(i).expect("in range"),
            clean.vertex_handle_at(i).expect("in range"),
        );
        assert_eq!(
            split.vertex_point(sv).as_array(),
            clean.vertex_point(cv).as_array(),
            "coordinates"
        );
    }
    for i in 0..split.face_count() as u32 {
        let s = split.face(split.face_handle_at(i).expect("in range"));
        let c = clean.face(clean.face_handle_at(i).expect("in range"));
        assert_eq!(s.surface, c.surface, "surface wiring");
    }
    // And the corner population is whole: every vertex holds a three-plane definition.
    let mut i = 0u32;
    while let Some(h_) = split.vertex_handle_at(i) {
        i += 1;
        let v = split.vertex(h_);
        assert!(
            matches!(*v, Vertex::ThreePlane(_)),
            "a corner without a three-plane definition survived: {v:?}"
        );
    }
}

/// ★ **Two *separated* collinear walls still intern to one surface** — the re-pin of what
/// `a_collinear_profile_vertex_gives_its_two_walls_one_surface` used to hold.
///
/// The dissolve pass only deletes flat corners (adjacent same-plane walls); two edges of a
/// notched profile lying on one line are legitimate geometry, their walls are two statements
/// of one plane, and interning makes them one handle.
/// Adjacent walls can no longer collide, so this is where "same plane = same handle" stays
/// pinned.
#[test]
fn two_separated_collinear_walls_intern_to_one_surface() {
    let mut m = Model::new();
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    // A right-edge notch: the two vertical segments at `x = 4` are collinear, not adjacent.
    let profile = Profile2d::polygon(vec![
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 1.0),
        p(3.0, 1.0),
        p(3.0, 2.0),
        p(4.0, 2.0),
        p(4.0, 3.0),
        p(0.0, 3.0),
    ])
    .unwrap();
    // Self-qualification: the dissolve pass must have left all eight corners standing.
    assert_eq!(profile.outer().vertices().len(), 8, "no corner is flat");
    let __w2 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w2,
            profile,
            dist: 1.0,
        },
    )
    .expect("extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    let shell = m.solid(solid).outer;
    let surfaces: Vec<_> = m
        .shell(shell)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .collect();
    // Eight profile points ⇒ eight wall quads, plus two caps.
    assert_eq!(surfaces.len(), 10, "eight walls and two caps");
    let mut distinct = surfaces.clone();
    distinct.sort_unstable_by_key(|h| h.index());
    distinct.dedup();
    assert_eq!(
        distinct.len(),
        9,
        "the notch's two x=4 walls are one plane, so one handle: {surfaces:?}"
    );
}

/// A vertex where one plane is split between two faces is named by **the planes that touch it**,
/// not by the loop's two neighbours.
///
/// The overhang chain puts the base's exposed top and the cantilever's underside on one plane
/// (`z = 1`, opposite normals — `unify` rightly keeps them apart, `canon` rightly calls them one
/// class). At `(1, 0.25, 1)` the side wall's loop runs straight through their shared line, so the
/// old rule named that vertex with the same plane twice: a triple defining no point, which the
/// exact predicates — whose precondition is `D ≠ 0` — aborted on. Measured as the only
/// path a degenerate triple reached them (16 arrivals in the OCCT suite, now 0).
#[test]
fn a_vertex_is_named_by_the_planes_that_touch_it() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let cutter = m.add_cuboid(
        Point3::from_array([1.1, 0.35, 0.5]),
        Point3::from_array([1.4, 0.65, 2.5]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, overhung, cutter).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let _ = &faces_tab;
    let mut checked = 0usize;
    for sh in solid_shell_handles(&m, overhung) {
        for &fh in &m.shell(sh).faces {
            let p = surf_ix[&fh];
            let tris = combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                .unwrap()
                .poly()
                .expect("a poly outer")
                .triples
                .clone();
            let tris = names(&tris);
            for t in &tris {
                // The triple is already dense plane ids: distinct means three real planes.
                assert!(
                    t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                    "vertex triple {t:?} names one plane twice"
                );
                // A name that denotes three distinct classes must denote a real point.
                assert!(
                    three_planes(
                        &planes[t[0]].plane,
                        &planes[t[1]].plane,
                        &planes[t[2]].plane
                    )
                    .is_some(),
                    "triple {t:?} defines no point"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "the chained operand has vertices to name");
}

/// **Dense plane ids are order-isomorphic to the sparse roots.** `dense_planes` ranks the class
/// roots, so any comparison, sort or lex-min over plane indices reads the same either way.
///
/// This is a **migration gate, not a permanent invariant**: it exists so the claim is measured
/// before the split rides on it, and it retires with `canon` — its subject, not its coverage,
/// is what goes away.
#[test]
fn dense_plane_ids_are_monotone_in_canon() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    // An overhanging boss splits `z = 1` between two faces, so classes really do merge and the
    // ranking really does compress — without that the map is the identity and proves nothing.
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let probe = m.add_cuboid(
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    // Rebuild the pieces `dense_planes` consumes, so this locks its contract without needing
    // `canon` to escape `plane_index_setup`. `plane_classes` is the same union-find the setup
    // runs; `dense_planes` the same ranking.
    let mut faces = collect_planes(&m, chained).unwrap();
    faces.extend(collect_planes(&m, probe).unwrap());
    let canon = plane_classes(&crate::planes::test_judge(&faces));
    let (geom, plane_ix, _cyls) = dense_planes(&faces, &canon);
    assert!(
        canon.iter().enumerate().any(|(i, &c)| c != i),
        "fixture has no split plane — the invariant would be vacuous"
    );
    assert!(
        geom.len() < canon.len(),
        "the ranking must actually compress"
    );
    for i in 0..canon.len() {
        for j in 0..canon.len() {
            assert_eq!(
                canon[i].cmp(&canon[j]),
                plane_ix[i].plane().cmp(&plane_ix[j].plane()),
                "faces {i}/{j}: canon {}/{} vs dense {}/{}",
                canon[i],
                canon[j],
                plane_ix[i].plane(),
                plane_ix[j].plane()
            );
        }
    }
}

/// **A plane triple is always in class form.** `planes` is a per-face table, so the same
/// `usize` could mean "face" or "plane"; producers settle it by emitting class roots, and a
/// consumer's raw `==` then means "same plane". Four silent-wrong bugs on this branch came from
/// the two meanings meeting in one comparison, so the invariant is asserted, not assumed.
#[test]
fn plane_triples_are_always_canon() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    // An *overhanging* boss splits `z = 1` between two faces with opposite normals (the base's
    // exposed top and the boss underside) — the shape that makes "face index" and "plane index"
    // differ at all. A boss sitting wholly inside the top merges into one holed face instead.
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let probe = m.add_cuboid(
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, chained, probe).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    // The fixture must actually merge two faces into one plane, or this proves nothing.
    assert!(
        planes.len() < faces_tab.len(),
        "fixture has no split plane — the invariant would be vacuous"
    );
    // A producer hands out dense plane ids (`loop_triples` maps face indices through
    // `plane_ix`), so "every element is a plane, not a face" is now the type, not a runtime
    // check. What remains testable is that the ids are in range and sorted-distinct.
    let mut checked = 0usize;
    for sh in solid_shell_handles(&m, chained) {
        for &fh in &m.shell(sh).faces {
            let p = surf_ix[&fh];
            let mut rings = vec![
                combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                    .unwrap()
                    .poly()
                    .expect("a poly outer")
                    .triples
                    .clone(),
            ];
            rings.extend(
                combinatorics::hole_rings(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                    .unwrap()
                    .into_iter()
                    .filter_map(|lr| lr.poly().map(|nr| nr.triples.clone())),
            );
            for t in rings.iter().flat_map(|r| names(r)) {
                for &k in &t {
                    assert!(
                        k < planes.len(),
                        "triple {t:?} names {k}, out of the plane table"
                    );
                }
                assert!(
                    t[0] < t[1] && t[1] < t[2],
                    "triple {t:?} is not three distinct planes in sorted order"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "the chained operand has vertices to name");
}

/// **A point reads zero on each of its three defining planes.** `Judge::orient3d`'s on-plane
/// shortcut is a raw `==` against the triple, so a vertex on the query plane must name it by the
/// same id the query uses. The face/plane split makes that automatic — a plane has exactly one
/// id now, so the old failure (a vertex named by face 6 of the `z = 1` class invisible to a
/// query about face 1 of it) cannot be expressed. What is left to check is the identity itself.
#[test]
fn a_vertex_on_the_cut_plane_reads_zero_whichever_face_names_it() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let probe = m.add_cuboid(
        Point3::from_array([1.1, 0.35, 0.5]),
        Point3::from_array([1.4, 0.65, 2.5]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, chained, probe).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    assert!(
        planes.len() < faces_tab.len(),
        "fixture has no split plane — the sibling faces this used to distinguish"
    );
    let mut on_plane = 0usize;
    for sh in solid_shell_handles(&m, chained) {
        for &fh in &m.shell(sh).faces {
            let p = surf_ix[&fh];
            let tris = combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                .unwrap()
                .poly()
                .expect("a poly outer")
                .triples
                .clone();
            let tris = names(&tris);
            for t in &tris {
                // The vertex lies on exactly its three defining planes; each must read 0.
                for &q in t {
                    assert_eq!(
                        combinatorics::side_of(
                            &jd,
                            &[],
                            combinatorics::NodeId::three_planes(combinatorics::Canon3::three(*t)),
                            q
                        ),
                        Some(0),
                        "vertex {t:?} lies on plane {q} but does not read 0"
                    );
                    on_plane += 1;
                }
            }
        }
    }
    assert!(on_plane > 0, "some vertex lies on some queried plane");
}
