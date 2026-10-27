use super::*;

/// The canonical form is the sorted triple, and the constructor is the only way to get one —
/// so the same three planes named in any order are **one** name.
///
/// The expected value is written out rather than derived, because a test that re-derives it
/// through the constructor would agree with the constructor however the constructor behaved.
#[test]
fn three_planes_names_one_vertex_however_it_is_spelled() {
    let canonical = NodeId(NodeKind::ThreePlane([2, 5, 9]));
    for spelling in [
        [2, 5, 9],
        [2, 9, 5],
        [5, 2, 9],
        [5, 9, 2],
        [9, 2, 5],
        [9, 5, 2],
    ] {
        assert_eq!(
            NodeId::three_planes(Canon3::three(spelling)),
            canonical,
            "{spelling:?} names the same vertex as [2, 5, 9]"
        );
    }
    assert_eq!(
        three_plane_name(canonical),
        Some([2, 5, 9]),
        "and the door agrees"
    );
}

/// **`Ord` is the bare triple's lexicographic order** — the proposition the whole migration
/// to this type stands on. Rules read it that would answer differently if it moved:
/// `Aliases::union_point`'s "smallest name wins", `merge_component`'s sorted ring starts, and
/// `reuse::canonical`'s ring rotation.
///
/// Each expected sign is written by hand from the pair, not computed from either operand, so
/// the test cannot agree with a wrong implementation by sharing its derivation.
#[test]
fn a_name_orders_like_the_triple_it_is() {
    use std::cmp::Ordering::{Equal, Greater, Less};
    for (a, b, want) in [
        ([0, 1, 2], [0, 1, 3], Less),    // last component decides
        ([0, 1, 3], [0, 1, 2], Greater), // …and antisymmetrically
        ([0, 2, 9], [0, 3, 4], Less),    // middle decides before last
        ([1, 0, 0], [0, 9, 9], Greater), // first decides before middle
        ([4, 4, 4], [4, 4, 4], Equal),
    ] {
        assert_eq!(
            NodeId(NodeKind::ThreePlane(a)).cmp(&NodeId(NodeKind::ThreePlane(b))),
            want,
            "{a:?} vs {b:?}"
        );
    }
}

/// **A pierce point named either way round is one name** — the [`NodeId::three_planes`]
/// proposition for the second variant, where the canonicalization is bigger than a sort.
///
/// The expected values are written from the rule's sentence, not by calling the constructor a
/// second time. The rule itself is [`crate::combinatorics::QuadRoot::canonical`]'s, and its own
/// locks (`swapping_the_planes_trades_the_two_roots_names` below) are what say the *geometry*
/// agrees; this only says the name reads it.
#[test]
fn a_pierce_is_one_name_however_the_pair_is_handed_in() {
    use crate::combinatorics::QuadRoot::{Double, Hi, Lo};
    for (root, other) in [(Lo, Hi), (Hi, Lo), (Double, Double)] {
        assert_eq!(
            NodeId::pierce(9, 2, 7, root),
            NodeId::pierce(2, 9, 7, other),
            "{root:?} handed in descending is {other:?} in stored order"
        );
    }
    // And the two roots of one pair stay **two** names — a canonicalization that collapsed
    // them would make this pass by making everything equal.
    assert_ne!(NodeId::pierce(2, 9, 7, Lo), NodeId::pierce(2, 9, 7, Hi));
    // A different cylinder through the same two planes is a different point.
    assert_ne!(NodeId::pierce(2, 9, 7, Lo), NodeId::pierce(2, 9, 8, Lo));
}

/// **The door tells the two variants apart, and `Ord` places them deterministically.**
///
/// The order itself is a *choice* — derived `Ord` puts every `ThreePlane` before every `Pierce`, so
/// "the smallest name wins" gains a systematic lean toward three-plane names at the six rules that
/// read it. It is written down here rather than left to be discovered.
#[test]
fn the_two_variants_are_told_apart_and_ordered() {
    let three = NodeId::three_planes(Canon3::three([9, 2, 5]));
    let pierce = NodeId::pierce(2, 9, 7, crate::combinatorics::QuadRoot::Lo);
    assert_eq!(three_plane_name(three), Some([2, 5, 9]));
    assert_eq!(
        three_plane_name(pierce),
        None,
        "a pierce point has no triple"
    );
    assert!(
        three < pierce,
        "declaration order: ThreePlane before Pierce"
    );
    // The probe helper drops exactly the pierce node, and keeps the ring's order otherwise.
    assert_eq!(
        three_plane_probes([pierce, three, pierce]),
        vec![[2, 5, 9]],
        "a probe list may lose a member; that is its licence"
    );
}

/// **The chart's shape, stated from the rule rather than from itself.**
///
/// One normal carries every clause: `n = 3e19 * (1, 2, 3)` is **wide** (a component of `9e19`,
/// whose square leaves `i128`, so an unrescaled `e2` is `None`) and **tilted** (`n_k != 0` for
/// the chosen `k`, without which two of the four clauses below cannot fail).
///
/// **The expected axes are written out by hand from `e1 = ê_k × n`, `e2 = n × e1`, never by calling
/// `of_normal` a second time.** Comparing `of_normal(K·n)` against `of_normal(n)` would let a chart
/// that ignores magnitude entirely (dropping a coordinate) pass on both sides while regressing the
/// census by nine rows (measured). An oracle has to come from the inputs.
///
/// Cross-check on `e2`: `n × (ê₀ × n) = |n|²ê₀ − n₀n = 14(1,0,0) − (1,2,3) = (13, −2, −3)`,
/// the same vector the cross product gives — two derivations, one answer.
#[test]
fn the_chart_is_an_orthogonal_in_plane_frame_along_the_normals_own_directions() {
    use nacre_exact::{Orient, Rat, dot_sign_rat, parallel_rat};
    let q = Rat::from_int;
    let k = 30_000_000_000_000_000_000i128; // 3e19
    let n = [q(k), q(2 * k), q(3 * k)];
    let chart = Chart2dRat::of_normal(&n)
        .expect("a wide normal is a rescale away from a chart, not a refusal");
    let (e1, e2) = chart.axes();
    for (got, want, name) in [
        (e1, [q(0), q(-3), q(2)], "e1 = ê₀ × n"),
        (e2, [q(13), q(-2), q(-3)], "e2 = n × e1"),
    ] {
        assert!(
            parallel_rat(got, &want) && dot_sign_rat(got, &want) == Orient::Positive,
            "{name}: {got:?} is not along {want:?}"
        );
        assert_eq!(
            dot_sign_rat(got, &n),
            Orient::Zero,
            "{name} is a direction *in* the plane, and `ring_interior_candidates` walks it as one"
        );
    }
    assert_eq!(
        dot_sign_rat(e1, e2),
        Orient::Zero,
        "the frame is orthogonal (not orthonormal), which is what makes `axes`' sentence — \
             the parity walks its ray along e1 — true of the world and not just of the chart"
    );
}

/// **A direction sign reads the truth, never the plane cache.** A prism over
/// `(0,0)·(1, 10⁻¹⁷)·(1,1)·(0,1)`: its wall from the origin runs along `(1, 10⁻¹⁷)`, so the
/// plane's world name is `(10⁻¹⁷·k, −k, 0, 0)` canonicalized — a first nonzero component of
/// `10⁻¹⁷` against the largest. The class's cache is then replaced by a rounded image of the same
/// plane with that component's sign flipped: an error of `2·10⁻¹⁷` against a unit normal, inside
/// one ulp of the large component, so nothing reading the cache as a rounding could object.
///
/// The class's facing ([`stored_coeffs_rat`]) and the canonical → outward sign
/// ([`Judge::outward_fix`](crate::tolerant::Judge::outward_fix)) come out the same: they read the world name and its sense. Comparing
/// the name's first nonzero component with the cache's reads the flipped sign, which the contrast
/// below states.
#[test]
fn a_direction_sign_does_not_read_the_plane_cache() {
    use crate::{OpOutput, Operation, SketchFrame, apply, from_rings};
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, nacre_exact::Axis::Z);
    let profile = from_rings(vec![vec![
        p2(0.0, 0.0),
        p2(1.0, 1e-17),
        p2(1.0, 1.0),
        p2(0.0, 1.0),
    ]])
    .expect("a quadrilateral")
    .remove(0);
    let Ok(OpOutput::Extrude { solid, .. }) = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    ) else {
        panic!("the prism")
    };
    m.rebuild_adjacency();
    let faces = crate::planes::collect_planes(&m, solid).expect("the face table");
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces));
    let (mut planes, _, _) = crate::planes::dense_planes(&faces, &canon);
    let zero = nacre_exact::Rat::from_int(0);
    let world = |k: usize| crate::planes::test_judge(&planes).world_rat(k);
    let c = (0..planes.len())
        .find(|&k| world(k).is_some_and(|w| w[0] != zero && w[1] != zero))
        .expect("the slanted wall");
    let w = world(c).expect("a world name");
    assert!(
        w[0].to_f64().abs() < 1e-15 * w[1].to_f64().abs(),
        "the name's first nonzero component is tiny: {w:?}"
    );
    let answers = |planes: &[WorkingPlane]| {
        let jd = crate::planes::test_judge(planes);
        (jd.stored_world_rat(c), jd.outward_fix(c))
    };
    let truth = answers(&planes);
    assert!(truth.0.is_some() && truth.1.is_some(), "both answer");

    let n = planes[c].plane.normal().as_array();
    assert!(
        n[0] != 0.0 && n[0].abs() < 1e-15,
        "the cache carries the tiny component: {n:?}"
    );
    let lying = nacre_geom::Plane::from_point_normal(
        planes[c].plane.origin(),
        nacre_math::Vector3::from_array([-n[0], n[1], n[2]]),
    )
    .expect("a plane");
    planes[c].plane = lying;
    // The contrast: with the lying cache planted, the first nonzero component of the name, compared
    // with the cache's, reads the other way — a direction sign computed that way would flip.
    assert_ne!(
        (w[0] > zero) == (lying.normal().as_array()[0] > 0.0),
        (w[0] > zero) == (n[0] > 0.0),
        "the planted cache flips the first-component reading"
    );
    assert_eq!(answers(&planes), truth, "the direction signs read no cache");
}

/// **An arc's span holds the seam point like any other** — `arc_span` against an oracle derived
/// apart from it: a point `P` of a circle is on the counter-clockwise arc `A → B` iff it lies to the
/// right of the chord, `orient2d(A, B, P) < 0`, viewed with the axis toward the viewer.
///
/// The circle is radius 5 about `+z` through the origin, and its twelve points with integer
/// coordinates (`(±5, 0)`, `(0, ±5)`, `(±3, ±4)`, `(±4, ±3)`) are every arc's ends and every
/// root. The seam is put on one of them three ways — and, the fourth, between them — so ends and
/// roots on the seam, arcs across it and arcs beside it all occur: every ordered pair of distinct
/// ends against every point, `AtLo`/`AtHi` at the ends. And every point as a **closed** arc's one
/// node (`ArcEnds::Closed`): every other point is inside, the node is `AtBothEnds` — while the
/// same two ends given as an *open* arc are a zero-span arc, which is refused.
#[test]
fn an_arc_span_reads_the_seam_point_as_any_other() {
    use nacre_exact::quad::{CylinderMeet, QuadVal, plane_plane_cylinder};
    use nacre_exact::{BigRat, Orient, Rat};
    let ri = Rat::from_int;
    let mut pts: Vec<(i128, i128)> = vec![(5, 0), (0, 5), (-5, 0), (0, -5)];
    for (a, b) in [(3, 4), (4, 3)] {
        for (sa, sb) in [(1, 1), (1, -1), (-1, 1), (-1, -1)] {
            pts.push((sa * a, sb * b));
        }
    }
    let r2 = BigRat::from(ri(25));
    for seam in [[1, 0], [0, -1], [3, 4], [1, 2]] {
        let def = nacre_topo::CylinderDef::new(
            [ri(0); 3],
            [ri(0), ri(0), ri(1)],
            [ri(seam[0]), ri(seam[1]), ri(0)],
            r2.clone(),
        )
        .expect("a cylinder");
        // The point `(x, y)` as the meet of `x = px` and `z = 0` and the root on its side of `y`.
        let at = |(px, py): (i128, i128)| {
            let meet = plane_plane_cylinder(
                &[ri(1), ri(0), ri(0), ri(-px)],
                &[ri(0), ri(0), ri(1), ri(0)],
                &def.origin(),
                &def.dir(),
                def.r2(),
            )
            .expect("fits");
            let (line, roots) = match meet {
                CylinderMeet::Pair { line, s } => (line, s.to_vec()),
                CylinderMeet::Tangent { line, s } => (line, vec![QuadVal::from_rat(s)]),
                other => panic!("({px}, {py}) is on the circle: {other:?}"),
            };
            let y_of = |s: &QuadVal| {
                QuadVal::from_rat(line.base()[1])
                    .checked_add(&s.checked_mul_rat(line.dir()[1]).unwrap())
                    .unwrap()
            };
            let s = *roots
                .iter()
                .find(|s| {
                    y_of(s)
                        .checked_sub(&QuadVal::from_rat(ri(py)))
                        .unwrap()
                        .sign()
                        == Orient::Zero
                })
                .expect("the root at the point");
            (line, s)
        };
        let mut decided = 0usize;
        for &a in &pts {
            for &b in &pts {
                if a == b {
                    continue;
                }
                for &p in &pts {
                    let want = if p == a {
                        ArcSpan::AtLo
                    } else if p == b {
                        ArcSpan::AtHi
                    } else {
                        let o = (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0);
                        if o < 0 {
                            ArcSpan::Inside
                        } else {
                            ArcSpan::Outside
                        }
                    };
                    let got = arc_span(&def, ArcEnds::Open(&at(a), &at(b)), &at(p));
                    assert_eq!(
                        got,
                        Some(want),
                        "seam {seam:?}: arc {a:?} → {b:?}, point {p:?}"
                    );
                    decided += 1;
                }
            }
        }
        for &a in &pts {
            assert_eq!(
                arc_span(&def, ArcEnds::Open(&at(a), &at(a)), &at(a)),
                None,
                "seam {seam:?}: an open arc from {a:?} to itself spans nothing"
            );
            for &p in &pts {
                let want = if p == a {
                    ArcSpan::AtBothEnds
                } else {
                    ArcSpan::Inside
                };
                let got = arc_span(&def, ArcEnds::Closed(&at(a)), &at(p));
                assert_eq!(
                    got,
                    Some(want),
                    "seam {seam:?}: closed arc at {a:?}, point {p:?}"
                );
                decided += 1;
            }
        }
        assert_eq!(decided, 12 * 11 * 12 + 12 * 12, "seam {seam:?}");
    }
}

/// **A node's interval holds the node, whichever kind of name it is** — [`node_point_hp`] against
/// [`node_point_approx`], over every ring node of every face a boolean emits. The fixture is a
/// plate with a bore fused with a wall whose plane runs inside the bore's radius, so the wall
/// class's rings carry the bore's pierce points beside three-plane corners, and the caps carry
/// arcs.
///
/// Bracketed, not contained: the `f64` is a rounding of the true point, and the interval is far
/// narrower than an ulp, so the true point is inside the interval while the `f64` may sit an ulp
/// outside it. The assertion reads both bounds of a window around the `f64`.
#[test]
fn a_nodes_interval_holds_it_for_either_kind_of_name() {
    use crate::arrangement::{PlaneSetup, plane_index_setup};
    use nacre_math::Point3;
    let mut m = Model::new();
    let plate = crate::tests::plate_with_bore(&mut m);
    let wall = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([9.0, 0.0, 5.0]),
        Point3::from_array([11.0, 40.0, 15.0]),
    );
    m.rebuild_adjacency();
    let setup = plane_index_setup(&m, plate, wall).expect("the gate passes");
    let PlaneSetup {
        planes: faces_tab,
        geom,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        class_owner,
        n_a,
        standard,
        notes,
        cyls,
        ..
    } = &setup;
    let jd = Judge::new(geom, *standard, notes);
    let trace_in = trace_input(
        &m,
        [(plate, inc_a), (wall, inc_b)],
        surf_ix,
        faces_tab.len(),
        &jd,
        plane_ix,
        cyls,
    );
    let (faces, _, _) = crate::arrangement::trace_result_faces_full_for_test(
        &m,
        crate::BoolKind::Fuse,
        plate,
        wall,
        &jd,
        faces_tab,
        plane_ix,
        cyls,
        *n_a,
        class_owner,
        &trace_in,
    )
    .expect("the fuse traces");
    let prec = jd.standard.prec;
    let mut kinds = [0usize; 2];
    for f in &faces {
        for &n in f.outer.ring().into_iter().flat_map(|r| &r.nodes) {
            let x = node_point_hp(&jd, cyls, n).expect("an interval for every emitted node");
            let p = node_point_approx(&jd, cyls, n).expect("a realized point");
            for k in 0..3 {
                let at = |d: f64| {
                    nacre_exact::HpBounded::of_rat(
                        nacre_exact::Rat::try_from_f64(p[k] + d).unwrap(),
                        prec,
                    )
                };
                assert_eq!(
                    x[k].sub(&at(1e-9), prec).sign(),
                    Some(false),
                    "{n:?} axis {k}"
                );
                assert_eq!(
                    x[k].sub(&at(-1e-9), prec).sign(),
                    Some(true),
                    "{n:?} axis {k}"
                );
            }
            kinds[usize::from(pierce_name(n).is_some())] += 1;
        }
    }
    assert!(
        kinds[0] > 0 && kinds[1] > 0,
        "both kinds of name are read: {kinds:?}"
    );
}

/// A radius-5 cylinder on the `z` axis (`z ∈ [−1, 11]`) and a `5×5` box turned 37° about `z` at the
/// axis — every wall plane cuts the cylinder and none has a world name. The table only, without the
/// gate (which would refuse this pair), with the cylinder's def; `wall(d)` finds the wall `n·x = d`
/// (`n = (cos 37°, sin 37°, 0)`) and `top` the cap `z = 11`.
fn nameless_walls_through_a_cylinder() -> (
    crate::arrangement::PlaneSetup,
    nacre_topo::CylinderDef,
    impl Fn(&crate::arrangement::PlaneSetup, f64) -> usize,
    impl Fn(&crate::arrangement::PlaneSetup) -> usize,
) {
    use crate::arrangement::plane_index_setup_inner;
    use nacre_exact::{Angle, Isometry, Rat, Rotation};
    use nacre_math::{Point3, Vector3};
    let mut m = Model::new();
    let cyl = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        5.0,
        12.0,
    )
    .solid;
    let block = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-2.5, -2.5, 0.0]),
        Point3::from_array([2.5, 2.5, 10.0]),
    );
    m.rebuild_adjacency();
    let block = crate::transform::transform(
        &mut m,
        block,
        &Isometry::rotation(Rotation {
            axis: nacre_exact::Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(37)).expect("an angle"),
        }),
    )
    .expect("the block turns");
    m.rebuild_adjacency();
    let (setup, cyl_surfs) = plane_index_setup_inner(&m, cyl, block).expect("the table builds");
    let nacre_topo::Surface::Cylinder { def, .. } = m.surface(cyl_surfs[0]) else {
        panic!("the cylinder's surface");
    };
    let def = def.clone();
    let on = |setup: &crate::arrangement::PlaneSetup, wc: usize, p: [f64; 3]| {
        setup.geom[wc]
            .plane
            .signed_distance(nacre_math::Point3::from_array(p))
            .abs()
            < 1e-9
    };
    let wall = move |setup: &crate::arrangement::PlaneSetup, d: f64| {
        let (c, s) = 37f64.to_radians().sin_cos();
        let foot = [d * s, d * c, 11.0];
        let along = [foot[0] - c, foot[1] + s, 3.0];
        (0..setup.geom.len())
            .find(|&wc| {
                setup.geom[wc].world.is_none() && on(setup, wc, foot) && on(setup, wc, along)
            })
            .expect("a wall with no world name")
    };
    let top = move |setup: &crate::arrangement::PlaneSetup| {
        (0..setup.geom.len())
            .find(|&wc| {
                setup.geom[wc].world.is_some()
                    && on(setup, wc, [0.0, 0.0, 11.0])
                    && on(setup, wc, [1.0, 2.0, 11.0])
            })
            .expect("the top cap z = 11")
    };
    (setup, def, wall, top)
}

/// **A pierce point on a plane with no world statement is an interval box, and the questions read
/// it.** A `5×5` box turned 37° about `z` at the axis of a radius-5 cylinder: every wall plane cuts
/// the cylinder, and none has a world name. The wall `n·x = 2.5` (`n = (cos 37°, sin 37°, 0)`) meets
/// the top cap `z = 11` in two points `2.5·n ± √18.75·n⊥`; the meet of that pair is
/// [`PiercePoint::Bounded`], its two roots are those two points (in some order), the point is on its
/// own planes by name, and three points of the top circle are in one cyclic order.
#[test]
fn a_pierce_point_on_a_nameless_plane_is_a_box_the_questions_read() {
    let (setup, def, wall_at, top_at) = nameless_walls_through_a_cylinder();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let def = &def;
    use nacre_exact::Rat;
    let k = 0; // the one cylinder class
    let (c, s) = 37f64.to_radians().sin_cos();
    let (n, nt) = ([s, c, 0.0], [-c, s, 0.0]);
    let foot = [2.5 * n[0], 2.5 * n[1], 11.0];
    let (wall, top) = (wall_at(&setup, 2.5), top_at(&setup));
    let t = 18.75f64.sqrt();
    let want = [
        [foot[0] + t * nt[0], foot[1] + t * nt[1], 11.0],
        [foot[0] - t * nt[0], foot[1] - t * nt[1], 11.0],
    ];
    let holds = |b: &[nacre_exact::HpBounded; 3], p: [f64; 3]| {
        (0..3).all(|i| {
            let lo = nacre_exact::HpBounded::of_rat(
                Rat::try_from_f64(p[i] - 1e-9).expect("finite"),
                128,
            );
            let hi = nacre_exact::HpBounded::of_rat(
                Rat::try_from_f64(p[i] + 1e-9).expect("finite"),
                128,
            );
            b[i].sub(&lo, 128).sign() == Some(true) && hi.sub(&b[i], 128).sign() == Some(true)
        })
    };
    let mut seen = [false; 2];
    let mut names = Vec::new();
    for root in [
        crate::combinatorics::QuadRoot::Lo,
        crate::combinatorics::QuadRoot::Hi,
    ] {
        let name = NodeId::pierce(wall, top, k, root);
        let meet = pierce_meet(&jd, k, def, name).expect("a bounded meet");
        assert!(meet.exact().is_none(), "no exact meet without a world name");
        let at = meet.point_hp(&jd).expect("its box");
        let hit: Vec<usize> = (0..2).filter(|&i| holds(&at, want[i])).collect();
        assert_eq!(
            hit.len(),
            1,
            "{root:?}'s box holds exactly one of the two points"
        );
        seen[hit[0]] = true;
        assert_eq!(meet.side(&jd, wall), Some(0), "on its own wall by name");
        assert_eq!(meet.side(&jd, top), Some(0), "on its own cap by name");
        names.push(name);
    }
    assert!(seen[0] && seen[1], "the two roots are the two points");
    // A third point of the top circle — the opposite wall's root — sits in exactly one of the two
    // arcs the first two cut the circle into.
    let opposite = wall_at(&setup, -2.5);
    let third = NodeId::pierce(opposite, top, k, crate::combinatorics::QuadRoot::Lo);
    let one = between_on_circle_hp(&jd, k, def, names[0], names[1], third).expect("decided");
    let other = between_on_circle_hp(&jd, k, def, names[1], names[0], third).expect("decided");
    assert_ne!(
        one, other,
        "a point of the circle is in exactly one of the two arcs"
    );
    // The closed arc one of them cuts — its whole turn back to it — holds every other point of
    // the circle and not its node.
    for (lo, n, want) in [
        (names[0], third, true),
        (names[0], names[1], true),
        (names[0], names[0], false),
    ] {
        assert_eq!(
            between_on_circle_hp(&jd, k, def, lo, lo, n),
            Some(want),
            "closed arc at {lo:?}, point {n:?}"
        );
    }
}

/// **An arc's crossings with a nameless plane are counted from boxes, both of them.** On the top
/// circle of the same cylinder, the wall `n·x = −2.5` cuts the circle at two points; of the two arcs
/// between them, the one on the far side holds both points where the wall `n·x = 2.5` cuts the
/// circle, and the other holds none — the count the ray road's walk reads where an arc's two ends
/// lie on one side of the plane.
#[test]
fn an_arc_crossed_twice_by_a_nameless_plane_counts_two() {
    let (setup, def, wall_at, top_at) = nameless_walls_through_a_cylinder();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (k, top) = (0, top_at(&setup));
    let (near, far) = (wall_at(&setup, 2.5), wall_at(&setup, -2.5));
    let ends = [
        NodeId::pierce(far, top, k, crate::combinatorics::QuadRoot::Lo),
        NodeId::pierce(far, top, k, crate::combinatorics::QuadRoot::Hi),
    ];
    let counts: Vec<usize> = [true, false]
        .into_iter()
        .map(|ccw| {
            arc_crossings(&jd, &def, top, near, k, ccw, ends, None)
                .expect("decided")
                .len()
        })
        .collect();
    let mut sorted = counts.clone();
    sorted.sort_unstable();
    assert_eq!(
        sorted,
        [0, 2],
        "one arc holds both crossings, the other none: {counts:?}"
    );
}

/// Keys with many ties, each tagged by where it started — the tag is what shows stability.
fn tagged(keys: &[u8]) -> Vec<(u8, usize)> {
    keys.iter().copied().zip(0..).collect()
}

proptest::proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// On a comparison that never fails, `try_sort_by` is the standard stable sort — ties in the
    /// order they came, on random, already sorted and reversed input, short and past the
    /// insertion length.
    #[test]
    fn a_sort_that_cannot_fail_is_the_standard_stable_sort(
        keys in proptest::collection::vec(0u8..6, 0..90),
        shape in 0u8..3,
    ) {
        let mut keys = keys;
        match shape {
            1 => keys.sort_unstable(),
            2 => keys.sort_unstable_by(|a, b| b.cmp(a)),
            _ => {}
        }
        let mut want = tagged(&keys);
        want.sort_by_key(|a| a.0);
        let mut got = tagged(&keys);
        try_sort_by(&mut got, |a, b| Ok::<_, ()>(a.0.cmp(&b.0))).unwrap();
        proptest::prop_assert_eq!(got, want);
    }
}

/// A comparison that fails for every pair one element takes part in — the shape a
/// cylinder-pinned point with no readable side has — reads as `Equal` nowhere: the sort stops
/// with the failure, does not panic, and leaves the slice a permutation of its input (a caller
/// names the failure from what the slice holds).
#[test]
fn a_failing_comparison_stops_the_sort_and_keeps_every_element() {
    for seed in 0..200u64 {
        let mut s = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let mut next = || {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            s >> 33
        };
        let v: Vec<(u64, bool)> = (0..40).map(|_| (next() % 1000, next() % 10 == 0)).collect();
        if !v.iter().any(|x| x.1) {
            continue;
        }
        let mut got = v.clone();
        let out = try_sort_by(&mut got, |a, b| {
            if a.1 || b.1 {
                Err(())
            } else {
                Ok(a.0.cmp(&b.0))
            }
        });
        assert_eq!(
            out,
            Err(()),
            "seed {seed}: a failed comparison must come back"
        );
        let (mut a, mut b) = (v, got);
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "seed {seed}: the slice must still hold its input");
    }
}

/// The same, failing **inside a merge**: each half holds one marked element and the comparison
/// fails only where two marked ones meet, so both halves sort and the merge is what stops. (With
/// failures spread as above, the insertions of the halves meet them first and no merge is ever
/// interrupted — a merge written straight into the slice would pass that test.)
#[test]
fn a_comparison_failing_mid_merge_keeps_every_element() {
    // Evens then odds, `20` and `21` marked: the halves sort without meeting a second mark, and
    // the merge, walking both, reaches the two marks together.
    let v: Vec<(u32, bool)> = (0..20)
        .map(|k| (2 * k, 2 * k == 20))
        .chain((0..20).map(|k| (2 * k + 1, 2 * k + 1 == 21)))
        .collect();
    let mut got = v.clone();
    let out = try_sort_by(&mut got, |a, b| {
        if a.1 && b.1 {
            Err(())
        } else {
            Ok(a.0.cmp(&b.0))
        }
    });
    assert_eq!(out, Err(()), "the two marks meet only in the merge");
    let (mut a, mut b) = (v, got);
    a.sort_unstable();
    b.sort_unstable();
    assert_eq!(
        a, b,
        "a merge that stops must leave the slice holding its input"
    );
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
/// A constructor that toggled *unconditionally* satisfies
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
/// This is the case `Double` exists for: with a tangency stored as `Lo`, a remap that
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

/// **The geometry the three locks above are about.** Without it they are statements about
/// a convention someone chose; with it they are statements about where the points are.
///
/// `z = 0` and `x = 0` meet on the y-axis, and the cylinder about the x-axis with `r = 2` is
/// crossed there at `(0, ∓2, 0)`. Canonical normals give `ℓ = (0,0,1) × (1,0,0) = +y`, so
/// solved in that order the smaller parameter is `y = −2`. Solved the other way `ℓ` is `−y`,
/// and the algebra `lo′ = −hi`, `hi′ = −lo` says `s′[1]` names the *same point* as `s[0]`.
#[test]
fn swapping_the_planes_trades_the_two_roots_names() {
    use nacre_exact::Rat;
    use nacre_exact::quad::{CylinderMeet, branch_point_approx, plane_plane_cylinder};
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
    let at = |l: &nacre_exact::quad::MeetLine, s| branch_point_approx(l, s);
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
