use super::*;
use nacre_exact::Angle;

/// ★★ **A moved pierce vertex names the crossing it moved to.** A definition's root is an
/// order along `ℓ = n₁ × n₂`, and a restatement may spell a moved plane with the opposite
/// normal — one reversal, trading `Lo` and `Hi`. Found by the commutation oracle the moment
/// the cache started reading the definition: under a 90° turn a boss's four pierce
/// vertices realized to the *other* crossing, half a boss away, and the stored `f64` had hidden
/// the wrong label since it was written. Locked on three quadrantal turns; the reflection half
/// of the rule is paid at the same site but a mirrored cylinder has no exact world statement
/// yet, so its pierce vertices keep the construction's figure and cannot be asked here.
#[test]
fn a_moved_pierce_vertex_names_the_crossing_it_moved_to() {
    use nacre_exact::{Isometry, Rotation};
    for (what, axis, deg) in [
        ("turned 90° about x", Axis::X, 90),
        ("turned 90° about z", Axis::Z, 90),
        ("turned 270° about y", Axis::Y, 270),
    ] {
        let iso = Isometry::rotation(Rotation {
            axis,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        });
        let motion = Xform::Rigid(&iso);
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([4.0, 4.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        let fused = crate::boolean::boolean(&mut m, crate::BoolKind::Fuse, plate, boss)
            .expect("a boss fuses onto its plate");
        m.rebuild_adjacency();
        let pierce_points = |m: &Model, s: Handle<Solid>| -> Vec<Point3> {
            let mut out = Vec::new();
            for &fh in &m.shell(m.solid(s).outer).faces {
                for he in &m.face(fh).outer.half_edges {
                    for vh in m.edge(he.edge).vertices {
                        if matches!(*m.vertex(vh), Vertex::Pierce { .. }) {
                            assert!(
                                matches!(m.vertex_cache(vh), PointCache::Bounded { .. }),
                                "{what}: a pierce vertex realizes: {:?}",
                                m.vertex_cache(vh)
                            );
                            out.push(m.vertex_point(vh));
                        }
                    }
                }
            }
            out
        };
        let images: Vec<Point3> = pierce_points(&m, fused[0])
            .into_iter()
            .map(|p| motion.point(p))
            .collect();
        assert!(!images.is_empty(), "{what}: the boss pierces its plate");
        let moved = transform_solid(&mut m, fused[0], &motion).expect("a quadrantal turn");
        m.rebuild_adjacency();
        let after = pierce_points(&m, moved);
        assert_eq!(after.len(), images.len(), "{what}");
        for p in &after {
            // The other crossing is half a boss away; the right one is within rounding.
            let nearest = images
                .iter()
                .map(|q| (*p - *q).norm())
                .fold(f64::INFINITY, f64::min);
            assert!(
                nearest < 1e-9,
                "{what}: a moved pierce vertex realized {nearest:e} from every image — the other crossing"
            );
        }
    }
}

/// **"The statements could carry this step" never excuses a chain from recording it.**
///
/// A motion every statement can carry needs no node *of its own* — but if the datum already has
/// a history, its statement is the plane before the chain and stays verbatim, so the chain has to
/// keep reproducing the motion; a chain missing a link describes the datum as it was before that
/// link. Silently.
///
/// Asserted on the rule rather than on a symptom, for a translation and for a reflection.
#[test]
fn an_exact_motion_is_still_recorded_once_there_is_a_history() {
    let mut m = Model::new();
    let root = m.push_motion(
        Motion::Rotate {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        },
        None,
    );
    let place = Isometry::translation([Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)]);
    let flip = Xform::mirror(Axis::X, Rat::from_int(0));

    for (what, motion, kind) in [
        ("translation", &Xform::Rigid(&place), "Translate"),
        ("reflection", &flip, "Mirror"),
    ] {
        // No history: an exact motion records nothing, and the coordinates stay the truth.
        assert_eq!(
            chain_motion(&mut m, None, motion, Carry::Full),
            None,
            "an exact {what} over no history records nothing"
        );
        // With a history: recorded anyway, and the new leaf hangs off the old one.
        let leaf = chain_motion(&mut m, Some(root), motion, Carry::Full)
            .unwrap_or_else(|| panic!("a {what} over a history is always recorded"));
        assert_ne!(leaf, root);
        assert_eq!(m.motion(leaf).parent, Some(root));
        // And it records **this** motion — a chain that keeps the wrong kind of node
        // reproduces the wrong datum just as silently as one that keeps none.
        let got = match m.motion(leaf).motion {
            Motion::Rotate { .. } => "Rotate",
            Motion::Translate { .. } => "Translate",
            Motion::Mirror { .. } => "Mirror",
            Motion::Frame { .. } => "Frame",
        };
        assert_eq!(got, kind, "a {what} records a {kind} node");
    }
}
/// A unit cube whose top cap carries a **deep** statement of `z = 1` — a triple with a `5⁴²`
/// denominator — returned with that denominator's `deep`, the triple and its surface.
fn deep_topped_cube(
    m: &mut Model,
) -> (
    Handle<Solid>,
    Rat,
    [[Rat; 3]; 3],
    nacre_store::Handle<nacre_topo::Surface>,
) {
    // ★★★ **State the deep plane first and let the cuboid intern onto it.** The three points
    // below name `z = 1` — the same plane the cuboid's top cap names — so the box's extrude
    // finds this handle by canonical name and the cap carries *this* statement. No test-only door
    // and no mutation: the arena holds the truth, and interning is the production road onto it.
    // (`z = 1` rather than `z = 0` because the world seeds already state the three origin
    // planes, and interning would hand back a seed's shallow triple.)
    let deep = Rat::new(1, 5i128.pow(42)).unwrap();
    let pts = [
        [deep, Rat::from_int(0), Rat::from_int(1)],
        [Rat::from_int(1), Rat::from_int(0), Rat::from_int(1)],
        [Rat::from_int(0), Rat::from_int(1), Rat::from_int(1)],
    ];
    let (deep_surf, _) = m.push_plane(
        Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .expect("a unit normal names a plane"),
        pts,
        None,
        nacre_topo::Orientation::Forward,
    );
    let s = crate::fixtures::cuboid(
        m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    assert!(
        m.shell(m.solid(s).outer)
            .faces
            .iter()
            .any(|&f| m.face(f).surface == deep_surf),
        "the cuboid's top cap must intern onto the deep statement"
    );
    (s, deep, pts, deep_surf)
}

/// ★★★ **A move whose rational point transport would overflow records a node instead of
/// dropping the points.** One surface's stored triple has a `5⁴²` denominator, so `q + t` with
/// `t = 2⁻³⁰` needs `lcm(5⁴², 2³⁰) ≈ 2.4e38 > i128`. Taking the no-node path would push the moved
/// surface point-less; `carry_of` puts the whole solid on the recorded path, and the original
/// triple survives verbatim as the pre-motion truth.
#[test]
fn an_overflowing_exact_move_records_a_node_and_keeps_the_points() {
    let mut m = Model::new();
    let (s, deep, pts, _) = deep_topped_cube(&mut m);
    // Fixture qualification: the rational side overflows.
    let t = Rat::new(1, 1 << 30).unwrap();
    assert!(deep.checked_add(t).is_none(), "the transport must overflow");

    // The z component keeps the deep z-plane off the invariant branch (a purely-x
    // translation *fixes* it, and a fixed plane is restated with no transport at all —
    // this test is about the transport overflowing).
    let iso = Isometry::translation([t, Rat::from_int(0), Rat::from_int(1)]);
    let moved = transform_solid(&mut m, s, &Xform::Rigid(&iso)).unwrap();
    m.rebuild_adjacency();
    let moved_surf = m
        .shell(m.solid(moved).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|&s2| {
            matches!(
                m.surface(s2),
                nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Known(p),
                    ..
                } if *p == pts
            )
        })
        .expect("the deep triple must survive the move verbatim (pre-motion truth)");
    assert!(
        matches!(
            m.surface(moved_surf),
            nacre_topo::Surface::Plane {
                motion: Some(_),
                ..
            }
        ),
        "the overflow must force the recorded path, not drop the points"
    );
}

/// ★★★ **Where one rigid motion and its two halves store differently, the world is still the
/// same.** `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)` holds as a statement about the
/// world; how it is stored follows what each step can carry. On the deep-topped cube with
/// `R = rz90` about the origin and `t = (2⁻³⁰, 0, 1)`:
///
/// - **at once**, the turn takes the deep coordinate from `x` to `y`, where `t` adds `0`, so every
///   statement moves exactly and the whole motion is carried — no node anywhere;
/// - **in two**, `R` fixes the top plane (its statement stays verbatim, `deep` still in `x`) and
///   is carried, and then `t` adds `2⁻³⁰` to `deep`, which leaves `i128` — so `T` is recorded on
///   every face.
///
/// The two models must still name every face's plane the same in the world and realize every
/// vertex to the same bits.
#[test]
fn an_overflow_splits_the_representation_not_the_world() {
    use nacre_exact::Rotation;
    let turn = Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
    };
    let shift = [
        Rat::new(1, 1 << 30).unwrap(),
        Rat::from_int(0),
        Rat::from_int(1),
    ];
    let said = |m: &Model, s: Handle<Solid>| -> (Vec<String>, Vec<[u64; 3]>, bool) {
        let faces = m.shell(m.solid(s).outer).faces.clone();
        let mut names: Vec<String> = faces
            .iter()
            .map(|&f| format!("{:?}", m.world_plane_name(m.face(f).surface)))
            .collect();
        names.sort();
        let mut bits: Vec<[u64; 3]> = faces
            .iter()
            .flat_map(|&f| m.face(f).outer.half_edges.clone())
            .map(|he| m.vertex_point(m.he_start(he)).as_array().map(f64::to_bits))
            .collect();
        bits.sort_unstable();
        bits.dedup();
        let recorded = faces
            .iter()
            .any(|&f| m.plane_motion(m.face(f).surface).is_some());
        (names, bits, recorded)
    };

    let mut at_once = Model::new();
    let (s, ..) = deep_topped_cube(&mut at_once);
    let rigid = Isometry::rigid(turn, shift);
    let a = transform_solid(&mut at_once, s, &Xform::Rigid(&rigid)).unwrap();
    at_once.rebuild_adjacency();

    let mut in_two = Model::new();
    let (s, ..) = deep_topped_cube(&mut in_two);
    let r = Isometry::rotation(turn);
    let t = Isometry::translation(shift);
    let b = transform_solid(&mut in_two, s, &Xform::Rigid(&r)).unwrap();
    in_two.rebuild_adjacency();
    let b = transform_solid(&mut in_two, b, &Xform::Rigid(&t)).unwrap();
    in_two.rebuild_adjacency();

    let (names_a, bits_a, recorded_a) = said(&at_once, a);
    let (names_b, bits_b, recorded_b) = said(&in_two, b);
    assert!(
        !recorded_a,
        "at once, every statement moves exactly: nothing is recorded"
    );
    assert!(
        recorded_b,
        "in two, the translation overflows the verbatim top: it is recorded"
    );
    assert!(
        names_a.iter().all(|n| n != "None"),
        "every face names its world plane: {names_a:?}"
    );
    assert_eq!(names_a, names_b, "one world, two representations");
    assert_eq!(bits_a, bits_b, "and one realization");
}

/// ★★ **A moved cylinder records its history instead of silently degrading.**
/// The truth variant has a motion slot of its own, and this is it working.
/// ★ Positive control: the remap gate **bites** — a walkable solid whose vertex
/// definition names a foreign surface is rejected with `OriginNotOnSolid`.
/// Without it, "the gate fired zero times across the
/// suite" would be indistinguishable from "the gate checks nothing".
#[test]
fn a_foreign_definition_is_rejected() {
    let mut m = Model::new();
    // A plane of this solid's own...
    let (own, _) = m.push_plane(
        Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 9.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap(),
        [
            [Rat::from_int(0), Rat::from_int(0), Rat::from_int(9)],
            [Rat::from_int(1), Rat::from_int(0), Rat::from_int(9)],
            [Rat::from_int(0), Rat::from_int(1), Rat::from_int(9)],
        ],
        None,
        nacre_topo::Orientation::Forward,
    );
    // ...and a vertex whose definition names the world seeds — surfaces this solid's face
    // set does not contain.
    let foreign = Vertex::ThreePlane([
        m.world_plane(Axis::Z),
        m.world_plane(Axis::X),
        m.world_plane(Axis::Y),
    ]);
    let mk_v = |m: &mut Model, x: f64| {
        m.push_vertex(
            foreign,
            PointCache::Unrealized {
                coord: Point3::from_array([x, 0.0, 9.0]),
            },
        )
    };
    let v0 = mk_v(&mut m, 0.0);
    let v1 = mk_v(&mut m, 1.0);
    let e = m
        .push_edge([own, m.world_plane(Axis::Z)], [v0, v1])
        .unwrap();
    let f = m.push_face_unchecked(Face {
        surface: own,
        outer: Loop {
            half_edges: vec![HalfEdge {
                edge: e,
                forward: true,
            }],
        },
        inner: vec![],
        orientation: nacre_topo::Orientation::Forward,
    });
    let sh = m.push_shell_unchecked(Shell { faces: vec![f] });
    let s = m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    let zero = Isometry::translation([Rat::from_int(0); 3]);
    assert_eq!(
        transform(&mut m, s, &zero),
        Err(OpError::OriginNotOnSolid),
        "a foreign definition must be rejected before the walk"
    );
}

/// ★★ A moved cylinder's seam vertices carry `OnSeam` re-pointed at the **twin's own**
/// surfaces — the carrier pair moves with the solid, like every other definition.
#[test]
fn a_moved_cylinders_seam_defs_repoint_to_the_twin() {
    let mut m = Model::new();
    let s = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([1.0, 2.0, 0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        1.5,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    let iso = Isometry::rotation(nacre_exact::Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: nacre_exact::Angle::from_deg(Rat::from_int(31)).expect("angle"),
    });
    let out = transform(&mut m, s, &iso).expect("rotate the cylinder");
    let mut twin_surfs = std::collections::HashSet::new();
    for &fh in &m.shell(m.solid(out).outer).faces {
        twin_surfs.insert(m.face(fh).surface);
    }
    let mut seams = 0;
    for &fh in &m.shell(m.solid(out).outer).faces {
        let face = m.face(fh);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                for &vh in &m.edge(he.edge).vertices {
                    if let Vertex::OnSeam(pair) = *m.vertex(vh) {
                        seams += 1;
                        assert!(
                            pair.iter().all(|c| twin_surfs.contains(c)),
                            "a moved seam def must name the twin's own surfaces"
                        );
                    }
                }
            }
        }
    }
    assert!(seams > 0, "the sweep saw no seam vertices at all");
}

#[test]
fn a_rotated_cylinder_records_its_motion() {
    let mut m = Model::new();
    let s = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        1.0,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let iso = Isometry::rotation(nacre_exact::Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: nacre_exact::Angle::from_deg(Rat::from_int(31)).expect("angle"),
    });
    let turned = transform_solid(&mut m, s, &Xform::Rigid(&iso)).unwrap();
    m.rebuild_adjacency();
    let lateral = m
        .shell(m.solid(turned).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|&su| matches!(m.surface_cache(su), nacre_geom::Surface::Cylinder(_)))
        .expect("a cylinder keeps its lateral face");
    // An inexact turn records a node, and the def is carried **verbatim** — the recorded
    // node states its cylinder before the motion (the plane rule).
    match m.surface(lateral) {
        nacre_topo::Surface::Cylinder {
            def,
            motion: Some(_),
        } => {
            let z = Rat::from_int(0);
            assert_eq!(def.origin(), [z; 3], "pre-motion statement, verbatim");
            assert_eq!(def.dir(), [z, z, Rat::from_int(1)]);
        }
        other => panic!("a moved cylinder's truth must carry the motion, got {other:?}"),
    }
    // ★ **`validate` on the *moved* cylinder population.** The net's rule that reads a cap's
    // rim circle against its plane reads two caches that a motion carries together. If they
    // ever stopped travelling together this is where it would show.
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
}

/// ★ The remap of a `Pierce` definition is not `.map(remap)` — pass 1 issues new
/// surface handles in face-traversal order, so the two planes' handle order can invert,
/// and re-sorting flips the canonical line direction ℓ = n₁×n₂, so the root must toggle
/// with the swap or `Lo` silently names the other point.
///
/// The fixture makes the swap *actually happen*: the pierce planes are the cylinder's
/// bottom cap (a fresh handle after a z-translation) and the world x = 0 seed (invariant
/// under that translation — it keeps handle 1), so `[cap(0), x0(1)]` remaps to
/// `[N, 1] → sorted [1, N]` — swapped. Geometry agrees with the toggle: with planes
/// `[cap, x0]` the canonical ℓ is +y and the y = −2 point is `Lo`; with `[x0, cap′]` ℓ
/// is −y and that same point is `Hi`.
#[test]
fn a_pierce_definition_swap_toggles_its_root() {
    use nacre_topo::QuadRoot;
    let mut m = Model::new();
    let s = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    )
    .solid;
    m.rebuild_adjacency();
    let shell = m.solid(s).outer;
    let faces = m.shell(shell).faces.clone();
    let lateral = faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|&su| matches!(m.surface_cache(su), nacre_geom::Surface::Cylinder(_)))
        .expect("lateral");
    let bottom = m.world_plane(Axis::Z); // the z = 0 cap interned onto the world seed
    let x0 = m.world_plane(Axis::X);
    assert!(bottom.index() < x0.index(), "the fixture's premise");
    // The two pierce points of {z = 0} ∧ {x = 0} against the cylinder: (0, ∓2, 0).
    // Canonical normals (0,0,1) × (1,0,0) = +y, so y = −2 is the smaller parameter: Lo.
    let v_lo = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, x0],
            cylinder: lateral,
            root: QuadRoot::Lo,
        },
        PointCache::Unrealized {
            coord: Point3::from_array([0.0, -2.0, 0.0]),
        },
    );
    let v_hi = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, x0],
            cylinder: lateral,
            root: QuadRoot::Hi,
        },
        PointCache::Unrealized {
            coord: Point3::from_array([0.0, 2.0, 0.0]),
        },
    );
    // Wire the vertices into the solid (a franken-face on the x = 0 seed): transform
    // remaps only what its face walk reaches, and `defs_are_remappable` requires every
    // carrier among the face surfaces.
    let edge = m
        .push_edge([bottom, x0], [v_lo, v_hi])
        .expect("a line through distinct endpoints");
    let franken = m.push_face(Face {
        surface: x0,
        outer: Loop {
            // A closed degenerate loop: v_lo -> v_hi -> v_lo. One half-edge would not
            // close, and `Model::push_face` now says so — the scaffold this test needs is a
            // face that *names the carriers*, never a malformed one.
            half_edges: vec![
                HalfEdge {
                    edge,
                    forward: true,
                },
                HalfEdge {
                    edge,
                    forward: false,
                },
            ],
        },
        inner: vec![],
        orientation: nacre_topo::Orientation::Forward,
    });
    let mut new_faces = faces.clone();
    new_faces.push(franken);
    let sh = m.push_shell(Shell { faces: new_faces });
    let franken_solid = m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    m.restore_live(vec![franken_solid]);

    let iso = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
    let moved = transform_solid(&mut m, franken_solid, &Xform::Rigid(&iso)).expect("moves");

    // Find the two pierce vertices of the moved solid and read their defs.
    let mut seen = Vec::new();
    let solid = m.solid(moved).clone();
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &fh in &m.shell(sh).faces {
            let face = m.face(fh).clone();
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in m.edge(he.edge).vertices.iter() {
                        if let Vertex::Pierce { planes, root, .. } = *m.vertex(vh) {
                            seen.push((m.vertex_point(vh).as_array(), planes, root));
                        }
                    }
                }
            }
        }
    }
    seen.sort_by(|a, b| a.0[1].partial_cmp(&b.0[1]).unwrap());
    seen.dedup_by_key(|e| e.0[1] as i64);
    assert_eq!(seen.len(), 2, "both pierce vertices survive the move");
    for (p, planes, root) in &seen {
        assert!(
            planes[0].index() < planes[1].index(),
            "stored order stays ascending"
        );
        // x = 0 kept its seed handle (invariant restatement); the cap moved to a fresh
        // one — so the pair swapped, and the root must have toggled with it.
        assert_eq!(planes[0], x0, "the surviving seed now sorts first");
        let want = if p[1] < 0.0 {
            QuadRoot::Hi
        } else {
            QuadRoot::Lo
        };
        assert_eq!(
            *root, want,
            "at y = {}: ℓ flipped to −y, so the root names the same point only if it \
                 toggled",
            p[1]
        );
    }
}

#[test]
fn an_exact_turn_transports_a_cylinders_truth_instead_of_recording() {
    let mut m = Model::new();
    let s = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.5, -1.25, 2.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        1.5,
        2.5,
    )
    .solid;
    m.rebuild_adjacency();
    let iso = Isometry::rotation(nacre_exact::Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: nacre_exact::Angle::from_deg(Rat::from_int(90)).expect("angle"),
    });
    let turned = transform_solid(&mut m, s, &Xform::Rigid(&iso)).unwrap();
    m.rebuild_adjacency();
    let lateral = m
        .shell(m.solid(turned).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|&su| matches!(m.surface_cache(su), nacre_geom::Surface::Cylinder(_)))
        .expect("a cylinder keeps its lateral face");
    // A 90°-family turn is exact: nothing is recorded, and the def rides the very
    // transport the probe checked — origin through `point_rat`, directions through
    // `dir_rat` (the pivot cancels), radius invariant. (x, y) ↦ (−y, x).
    let d = |x: f64| Rat::from_decimal(x).expect("decimal");
    match m.surface(lateral) {
        nacre_topo::Surface::Cylinder { def, motion: None } => {
            assert_eq!(def.origin(), [d(1.25), d(0.5), d(2.0)]);
            assert_eq!(def.dir(), [d(0.0), d(0.0), d(1.0)]);
            assert_eq!(
                def.ref_dir(),
                [d(1.0), d(0.0), d(0.0)],
                "seam turned with it"
            );
            assert_eq!(
                *def.r2(),
                nacre_exact::BigRat::from(d(2.25)),
                "the squared radius is rigid-invariant"
            );
        }
        other => panic!("an exact turn must transport, not record — got {other:?}"),
    }
}
