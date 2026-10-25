//! Corners where two arcs of different circles meet: `Vertex::CylinderPair`.

use super::*;
use nacre_exact::{Angle, Isometry, Rat, Rotation};
use nacre_topo::PointCache;

/// The leaf: two quarter arcs of radius 5 between `(0, 0)` and `(5, 5)`, about `(5, 0)` and
/// `(0, 5)`, extruded `h` along `+z`. `from_far` starts the pen at `(5, 5)` instead, which pushes
/// the `(0, 5)` arc's wall first.
fn leaf(m: &mut Model, h: f64, from_far: bool) -> Handle<Solid> {
    let edges = if from_far {
        vec![
            arc_turns(p2(0.0, 5.0), p2(5.0, 5.0), -1),
            arc_turns(p2(5.0, 0.0), p2(0.0, 0.0), -1),
        ]
    } else {
        vec![
            arc_turns(p2(5.0, 0.0), p2(0.0, 0.0), -1),
            arc_turns(p2(0.0, 5.0), p2(5.0, 5.0), -1),
        ]
    };
    let profile = stated(edges).expect("a leaf").remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("a leaf extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// Every vertex of `s`, once, in loop order.
fn vertices_of(m: &Model, s: Handle<Solid>) -> Vec<Handle<Vertex>> {
    let mut out = Vec::new();
    for &fh in &m.shell(m.solid(s).outer).faces {
        let f = m.face(fh);
        for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
            for he in &lp.half_edges {
                for vh in m.edge(he.edge).vertices {
                    if !out.contains(&vh) {
                        out.push(vh);
                    }
                }
            }
        }
    }
    out
}

/// The centre `(x, y)` of a cylinder the leaf's frame states (axis along `z`).
fn centre_xy(m: &Model, cyl: Handle<Surface>) -> [f64; 2] {
    let nacre_topo::Surface::Cylinder { def, .. } = m.surface(cyl) else {
        panic!("a cylinder")
    };
    let o = def.origin();
    [o[0].to_f64(), o[1].to_f64()]
}

/// **A leaf builds, and its four corners are points on two cylinders.** The prism is valid, its
/// volume is `25(π/2 − 1)·h`, and every corner is realized exactly at the point the sketch stated.
/// The ranks are counted by hand about `cylinders[0]`, whose `ref_dir` is the frame's `x̂`: about
/// `(5, 0)` the corner `(5, 5)` sits at `θ = π/2` and `(0, 0)` at `π`; about `(0, 5)`, `(5, 5)` is the
/// seam point (`θ = 0`) and `(0, 0)` is at `3π/2`. Drawn from either end, so each cylinder is the
/// one the rank is read about once — the seam point is ranked in the second.
#[test]
fn a_leaf_builds_and_ranks_its_corners_about_the_lower_cylinder() {
    let h = 2.0;
    for (from_far, anchor) in [(false, [5.0, 0.0]), (true, [0.0, 5.0])] {
        let mut m = Model::new();
        let s = leaf(&mut m, h, from_far);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let want = 25.0 * (std::f64::consts::FRAC_PI_2 - 1.0) * h;
        let got = nacre_props::mass_props(&m, s).expect("props").volume;
        assert!((got - want).abs() < 1e-9, "volume {got} vs {want}");
        let vs = vertices_of(&m, s);
        assert_eq!(vs.len(), 4);
        for vh in vs {
            let Vertex::CylinderPair {
                cylinders, rank, ..
            } = *m.vertex(vh)
            else {
                panic!("an arc–arc corner is a cylinder pair: {:?}", m.vertex(vh))
            };
            assert_eq!(
                centre_xy(&m, cylinders[0]),
                anchor,
                "the lower handle's cylinder"
            );
            let PointCache::Realized { coord, .. } = m.vertex_cache(vh) else {
                panic!("{:?}", m.vertex_cache(vh))
            };
            let p = coord.as_array();
            let want_rank = match [p[0], p[1]] {
                [5.0, 5.0] => 0,
                [0.0, 0.0] => 1,
                other => panic!("a corner at {other:?}"),
            };
            assert_eq!(rank, want_rank, "{p:?} about {anchor:?}");
            assert!(p[2] == 0.0 || p[2] == h, "{p:?}");
            let again = crate::realize_vertex(&m, vh, crate::Precision::NearestF64)
                .expect("realized")
                .to_f64();
            assert_eq!(
                again.map(|(x, _)| x),
                Some(p),
                "the realization is the cache"
            );
        }
        super::mesh_covers_faces("the leaf", &m, &[s]);
        assert!(
            nacre_step::to_step(&m, fixtures::STAMP)
                .expect("the leaf exports")
                .contains("MANIFOLD_SOLID_BREP")
        );
        assert_eq!(crate::refine_caches(&mut m), crate::RefineReport::default());
    }
}

/// **Two tangent circles meet once.** An S: a bump over `x ∈ [0, 2]` and a dip under
/// `x ∈ [2, 4]`, the two half circles touching at `(2, 0)`, closed by three straight sides
/// tangent at the ends. The joint is `rank: 0`; the bump and the dip cancel, so the volume is the
/// `4 × 3` rectangle's.
#[test]
fn an_s_joint_of_two_tangent_circles_is_one_point() {
    let h = 1.0;
    let profile = stated(vec![
        arc_turns(p2(1.0, 0.0), p2(0.0, 0.0), -2),
        arc_turns(p2(3.0, 0.0), p2(2.0, 0.0), 2),
        line(p2(4.0, 0.0), p2(4.0, -3.0)),
        line(p2(4.0, -3.0), p2(0.0, -3.0)),
        line(p2(0.0, -3.0), p2(0.0, 0.0)),
    ])
    .expect("an S")
    .remove(0);
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: s, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("the S extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let got = nacre_props::mass_props(&m, s).expect("props").volume;
    assert!((got - 12.0 * h).abs() < 1e-9, "volume {got}");
    let joints: Vec<_> = vertices_of(&m, s)
        .into_iter()
        .filter(|&vh| matches!(m.vertex(vh), Vertex::CylinderPair { .. }))
        .collect();
    assert_eq!(joints.len(), 2, "the joint, on both caps");
    for vh in joints {
        assert!(matches!(m.vertex(vh), Vertex::CylinderPair { rank: 0, .. }));
        let PointCache::Realized { coord, .. } = m.vertex_cache(vh) else {
            panic!("{:?}", m.vertex_cache(vh))
        };
        assert_eq!([coord[0], coord[1]], [2.0, 0.0]);
    }
    super::mesh_covers_faces("the S", &m, &[s]);
}

/// **The rank is read from the cylinder's own statement.** The leaf's two corners about the
/// `(5, 0)` cylinder, with its seam turned four ways: the order follows the seam, and the corner
/// on the seam ranks first.
#[test]
fn the_rank_follows_the_seam_of_the_first_cylinder() {
    let r = Rat::from_int;
    let z = [r(0), r(0), r(1)];
    let cyl = |o: [i128; 3], e: [i128; 3]| {
        nacre_topo::CylinderDef::new(o.map(r), z, e.map(r), nacre_exact::BigRat::from(r(25)))
            .expect("a cylinder")
    };
    let other = cyl([0, 5, 0], [1, 0, 0]);
    let corner = |p: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)| -> [Rat; 2] {
        let t = p.1.as_rat().expect("rational corners");
        let q = p.0.point_at(t).expect("no overflow");
        [q[0], q[1]]
    };
    let (near, far) = ([r(0), r(0)], [r(5), r(5)]);
    // θ about (5, 0): (5, 5) is +y, (0, 0) is −x.
    for (seam, want) in [
        ([1, 0, 0], [far, near]),  // seam +x: π/2, then π
        ([-1, 0, 0], [near, far]), // seam −x: (0, 0) on it
        ([0, 1, 0], [far, near]),  // seam +y: (5, 5) on it
        ([0, -1, 0], [far, near]), // seam −y: π, then 3π/2
    ] {
        // The rulings' feet are on `z = 0`, the plane across the axis through `(5, 0, 0)`.
        let got = crate::realize::pair_rulings(&cyl([5, 0, 0], seam), &other).expect("no overflow");
        assert_eq!(
            got.iter().map(corner).collect::<Vec<_>>(),
            want.to_vec(),
            "seam {seam:?}"
        );
    }
}

/// **A moved leaf keeps its corners' names, and the realization finds them where the motion put
/// them.** Moved, turned a quarter (both fold, so the corners stay exact), turned 30° (no fold: the
/// met road) about `y` and about the caps' normal, and moved twice — a turn about the caps' normal
/// fixes the caps, which keep their statement while the walls take the node, so a second motion
/// puts the caps and the walls on different chains. Each time valid, every corner still the same
/// rank, and the cache the moved point. Mirroring is refused for a solid with cylinders
/// (`mirroring_a_cylinder_is_declined`), so no reflection reaches a rank.
#[test]
fn a_moved_leaf_keeps_its_corner_names() {
    let r = |n: i128, d: i128| Rat::new(n, d).expect("a rational");
    let turn = |axis, deg| {
        Isometry::rotation(Rotation {
            axis,
            pivot: [r(1, 2), r(-3, 1), r(0, 1)],
            angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
        })
    };
    let up = Isometry::translation([r(0, 1), r(0, 1), r(5, 2)]);
    for (name, moves) in [
        (
            "moved",
            vec![Isometry::translation([r(3, 2), r(-7, 1), r(4, 1)])],
        ),
        ("a quarter about z", vec![turn(Axis::Z, 90)]),
        ("a quarter about x", vec![turn(Axis::X, 90)]),
        ("30° about y", vec![turn(Axis::Y, 30)]),
        ("30° about z", vec![turn(Axis::Z, 30)]),
        ("30° about z, then up", vec![turn(Axis::Z, 30), up]),
        (
            "30° about z, then 30° about x",
            vec![turn(Axis::Z, 30), turn(Axis::X, 30)],
        ),
        (
            "30° about z, then a quarter about x",
            vec![turn(Axis::Z, 30), turn(Axis::X, 90)],
        ),
    ] {
        let mut m = Model::new();
        let mut s = leaf(&mut m, 2.0, false);
        let before: Vec<(u8, Point3)> = vertices_of(&m, s)
            .into_iter()
            .map(|vh| match *m.vertex(vh) {
                Vertex::CylinderPair { rank, .. } => (rank, m.vertex_point(vh)),
                ref other => panic!("{other:?}"),
            })
            .collect();
        for iso in &moves {
            s = fixtures::xf(&mut m, s, *iso);
        }
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{name}: {:?}",
            nacre_validate::validate(&m)
        );
        let after = vertices_of(&m, s);
        assert_eq!(after.len(), before.len(), "{name}");
        for (vh, (rank, p)) in after.into_iter().zip(before) {
            assert!(
                matches!(*m.vertex(vh), Vertex::CylinderPair { rank: k, .. } if k == rank),
                "{name}: {:?}",
                m.vertex(vh)
            );
            let want = Point3::from_array(
                moves
                    .iter()
                    .fold(p.as_array(), |q, iso| iso.point_approx(q)),
            );
            let got = m.vertex_point(vh);
            assert!(
                (got - want).norm() < 1e-12,
                "{name}: corner at {got:?}, the motion puts it at {want:?}"
            );
            let again = crate::realize_vertex(&m, vh, crate::Precision::NearestF64)
                .unwrap_or_else(|e| panic!("{name}: {e:?}"))
                .to_f64();
            assert_eq!(again.map(|(x, _)| x), Some(got.as_array()), "{name}");
        }
    }
}

/// **A leaf stands on every face of a turned box.** The box is turned 30° about `x`, so its
/// faces sit on the turn's node; a sketch on one of them is a frame node of its own, the leaf's
/// walls stand on it, and the base cap is the box's face — a different chain. The rank is the
/// walls' alone (their shared ruling), so each corner is named, and the realization meets the
/// ruling with the face in the world.
#[test]
fn a_leaf_stands_on_every_face_of_a_turned_box() {
    let r = Rat::from_int;
    let mut m = Model::new();
    let block = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-10.0, -10.0, -10.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    let block = fixtures::xf(
        &mut m,
        block,
        Isometry::rotation(Rotation {
            axis: Axis::X,
            pivot: [r(0), r(0), r(0)],
            angle: Angle::from_deg(r(30)).expect("an angle"),
        }),
    );
    m.rebuild_adjacency();
    let faces = m.shell(m.solid(block).outer).faces.clone();
    assert_eq!(faces.len(), 6);
    let profile = stated(vec![
        arc_turns(p2(5.0, 0.0), p2(0.0, 0.0), -1),
        arc_turns(p2(0.0, 5.0), p2(5.0, 5.0), -1),
    ])
    .expect("a leaf")
    .remove(0);
    for fh in faces {
        let frame = crate::face_sketch_frame(&m, fh).expect("a planar live face");
        let OpOutput::Extrude { solid: leaf, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: profile.clone(),
                dist: 2.0,
            },
        )
        .unwrap_or_else(|e| panic!("a leaf on {fh:?}: {e:?}")) else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let want = 25.0 * (std::f64::consts::FRAC_PI_2 - 1.0) * 2.0;
        let got = nacre_props::mass_props(&m, leaf).expect("props").volume;
        assert!((got - want).abs() < 1e-9, "{fh:?}: volume {got} vs {want}");
        for vh in vertices_of(&m, leaf) {
            assert!(matches!(m.vertex(vh), Vertex::CylinderPair { .. }));
            let again = crate::realize_vertex(&m, vh, crate::Precision::NearestF64)
                .unwrap_or_else(|e| panic!("{fh:?}: {e:?}"))
                .to_f64()
                .map(|(x, _)| x);
            assert_eq!(again, Some(m.vertex_point(vh).as_array()), "{fh:?}");
        }
    }
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
}

/// **A boolean refuses an operand with an arc–arc corner by that name**, whatever the other
/// operand is: one crossing the leaf's corner, one far away — the name is the operand's, and the
/// operands stay live.
#[test]
fn a_boolean_refuses_a_two_cylinder_edge_at_the_entry() {
    for (lo, hi) in [
        ([2.0, 2.0, -1.0], [8.0, 8.0, 1.0]),
        ([20.0, 20.0, 0.0], [21.0, 21.0, 1.0]),
    ] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let mut m = Model::new();
            let l = leaf(&mut m, 2.0, false);
            let b = crate::fixtures::cuboid(&mut m, Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            let live = m.live_solids().to_vec();
            for (x, y) in [(l, b), (b, l)] {
                let out = crate::boolean(&mut m, kind, x, y);
                assert!(
                    matches!(
                        out,
                        Err(BoolError::Rejected {
                            reason: RejectReason::CylinderPairEdge,
                            at: None
                        })
                    ),
                    "{kind:?}: {out:?}"
                );
                assert_eq!(m.live_solids(), live.as_slice(), "the operands stay");
            }
        }
    }
}

/// **A plane along the axis names no point on the ruling, and says so.** A cylinder-pair definition
/// handed a box's side wall instead of a cap: every carrier is stated in the world, so the exact
/// road sees the plane run along the axis — no precision would place the point, so the answer is
/// `NoCurvedPoint`, not `Undecided` (which a cache reads as «ask with more bits»).
#[test]
fn a_cylinder_pair_on_a_plane_along_the_axis_names_no_point() {
    let mut m = Model::new();
    let s = leaf(&mut m, 2.0, false);
    let Vertex::CylinderPair { cylinders, .. } = *m.vertex(vertices_of(&m, s)[0]) else {
        panic!("a cylinder pair")
    };
    let block = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([10.0, 0.0, 0.0]),
        Point3::from_array([11.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let wall = m
        .shell(m.solid(block).outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .find(|&h| {
            m.surface_name(h)
                .and_then(|n| n.narrow().copied())
                .is_some_and(|c| c[2] == Rat::from_int(0))
        })
        .expect("a side wall");
    let vh = crate::realize::push_vertex_realized(
        &mut m,
        Vertex::CylinderPair {
            plane: wall,
            cylinders,
            rank: 0,
        },
        PointCache::Unrealized {
            coord: Point3::from_array([0.0; 3]),
        },
        crate::realize::ChainLink::Fresh,
    );
    assert_eq!(
        crate::realize_vertex(&m, vh, crate::Precision::NearestF64).err(),
        Some(crate::RealizeError::NoCurvedPoint)
    );
    assert!(matches!(m.vertex_cache(vh), PointCache::Unrealized { .. }));
}
