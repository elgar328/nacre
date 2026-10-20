//! Loop nesting, ring winding and ray casting on planar faces; voids and islands.

use super::*;

/// A slab over the pocketed cube, its underside at height `z0`. The rectangle is
/// asymmetric so that the cube's four vertical edges, which pierce the underside at
/// `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`, miss its fan diagonals; a square slab has all
/// four apexes degenerate at once.
///
/// `z0 = 0.3` runs below the pocket floor, so the slab's underside meets only the
/// cube's outer walls. `z0 = 0.7` runs between the floor and the lid and meets the
/// pocket walls as well — two loops, nested.
fn pocket_and_slab(z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, pc) = pocketed_cube();
    let slab = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.2, -0.25, z0]),
        Point3::from_array([1.3, 1.2, 1.5]),
    );
    (m, slab, pc)
}

/// Nested loops: a pocket sliced between its floor and its lid is enough to reach them. The
/// slab's underside carries the cube's cross-section as one loop and the pocket's inside
/// it — a loop within a loop, which the pairwise guard over-rejected.
///
/// `Cut` opens, both orders: `A ∩ B = B ∩ {z ≥ 0.7}` = `0.30 − 0.048 = 0.252`, slab `1.74`,
/// pocketed cube `0.92`, so `Cut(slab,pc) = 1.488` and `Cut(pc,slab) = 0.668`. `Cut(pc,slab)`
/// is the one that makes the cube cross-section an *island with a hole* (the slab is B, its
/// underside dropped, so no region survives to own the loops); `Cut(slab,pc)` keeps the slab
/// region and hangs the pocket loop in it as a hole beside an island.
///
/// `Fuse` seals the pocket (`[0.3,0.7]² × [0.5,0.7]`, capped by the slab at `z = 0.7`) into
/// an enclosed cavity — a second shell. The assembly builds it: the union material is
/// `2.408` and the result carries one cavity of volume `0.032`, two shells, `validate`
/// clean (the void's inward orientation). Both orders (Fuse is commutative).
#[test]
fn a_slab_between_the_lid_and_the_floor_nests_two_loops() {
    for (swap, expect) in [(false, 1.488), (true, 0.668)] {
        let (mut m, slab, pc) = pocket_and_slab(0.7);
        let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
        let r = boolean_one(&mut m, BoolKind::Cut, x, y).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "Cut swap={swap}: {vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - expect).abs() < 1e-9,
            "Cut swap={swap}: {} vs {expect}",
            props.volume
        );
    }
    // Fuse seals the pocket into a cavity — a second shell.
    for swap in [false, true] {
        let (mut m, slab, pc) = pocket_and_slab(0.7);
        let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
        let r = boolean_one(&mut m, BoolKind::Fuse, x, y).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "Fuse swap={swap}: {vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 2.408).abs() < 1e-9,
            "Fuse swap={swap}: {}",
            props.volume
        );
        assert_eq!(m.solid(r).cavities.len(), 1, "Fuse swap={swap}");
        assert_eq!(m.reachable().shells.len(), 2, "Fuse swap={swap}");
    }
}

/// A ring's nodes as coordinates — each triple is three planes, so its point is their meet.
fn ring_points(planes: &[WorkingPlane], ring: &[[usize; 3]]) -> Vec<[f64; 3]> {
    ring.iter()
        .map(|t| {
            three_planes(
                &planes[t[0]].plane,
                &planes[t[1]].plane,
                &planes[t[2]].plane,
            )
            .unwrap()
            .as_array()
        })
        .collect()
}

#[test]
fn a_hole_winds_clockwise_and_an_island_counter_clockwise() {
    // The two rings of a holed face are stored with opposite windings — that is what makes one
    // a hole and the other its outer boundary — and `loop_winding` must read exactly that.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap()
        )
        .unwrap(),
        1
    );
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &hole).unwrap()
        )
        .unwrap(),
        -1
    );

    // Nothing but the ring's direction went into that. Reversing it by hand agrees.
    for (name, ring, want) in [("outer", &outer, -1i8), ("hole", &hole, 1)] {
        let mut reversed = ring.clone();
        reversed.reverse();
        assert_eq!(
            combinatorics::loop_winding(
                &crate::planes::test_judge(&planes),
                NO_CYLS,
                p,
                &combinatorics::ring_from_names(p, &reversed).unwrap()
            )
            .unwrap(),
            want,
            "{name} reversed"
        );
    }
}

/// **A ring pinched at its extreme node winds by its outermost passage.** Two rings on the plane
/// `z = 1`, each passing `P = (2, 2)` — their lexicographically smallest point, no arc anywhere —
/// twice, built straight from names so the shape is exactly this:
///
/// - two clockwise lobes touching at `P` (the triangles `P (4,4) (4,3)` and `P (4,1) (4,0)`) — two
///   holes meeting at a point, the contour of the face around them: `−1`;
/// - a counter-clockwise triangle `P (4,0) (4,4)` with a clockwise hole `P (3,2.5) (3,1.5)` inside
///   it, touching its outer boundary at `P` — one ring of an outer boundary and a hole: `+1`.
///
/// Each is asked from all six starting nodes and reversed, so neither the start index nor which of
/// the two passages comes first decides it. The arrangement delivers rings of the first kind — a
/// plate whose hole touches a boss's footprint at one corner gives one hole ring around both
/// triangles ([`a_plate_pinched_where_a_hole_meets_a_boss_meshes`]); here the ring is built from
/// names so the rule's many-passage reading is held on exactly this shape.
#[test]
fn a_ring_pinched_at_its_extreme_node_winds_by_its_outermost_passage() {
    let mut m = Model::new();
    let at = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
    let wall = |m: &mut Model, a: (f64, f64), b: (f64, f64)| {
        plane_through(m, [at(a.0, a.1, 0.0), at(b.0, b.1, 0.0), at(a.0, a.1, 1.0)])
    };
    let top = plane_through(
        &mut m,
        [at(0.0, 0.0, 1.0), at(1.0, 0.0, 1.0), at(0.0, 1.0, 1.0)],
    );
    let planes = vec![
        top,                                  // 0: z = 1
        wall(&mut m, (2.0, 2.0), (4.0, 3.0)), // 1: through P along (2, 1)
        wall(&mut m, (2.0, 2.0), (4.0, 4.0)), // 2: through P along (1, 1)
        wall(&mut m, (2.0, 2.0), (4.0, 1.0)), // 3: through P along (2, −1)
        wall(&mut m, (2.0, 2.0), (4.0, 0.0)), // 4: through P along (1, −1)
        wall(&mut m, (4.0, 0.0), (4.0, 1.0)), // 5: x = 4
        wall(&mut m, (3.0, 0.0), (3.0, 1.0)), // 6: x = 3
    ];
    let jd = crate::planes::test_judge(&planes);
    let named = |ring: [[usize; 3]; 6]| -> Vec<[usize; 3]> {
        ring.map(|mut t| {
            t.sort_unstable();
            t
        })
        .to_vec()
    };
    // `P` is named by the two walls the ring arrives and leaves along at each passage.
    let two_holes = named([
        [0, 4, 2],
        [0, 2, 5],
        [0, 5, 1],
        [0, 1, 3],
        [0, 3, 5],
        [0, 5, 4],
    ]);
    let outer_and_hole = named([
        [0, 3, 4],
        [0, 4, 5],
        [0, 5, 2],
        [0, 2, 1],
        [0, 1, 6],
        [0, 6, 3],
    ]);
    for (name, ring, want) in [
        ("two holes", two_holes, -1i8),
        ("outer and hole", outer_and_hole, 1),
    ] {
        for reversed in [false, true] {
            for start in 0..ring.len() {
                let mut r = ring.clone();
                if reversed {
                    r.reverse();
                }
                r.rotate_left(start);
                let edges = combinatorics::ring_from_names(0, &r).expect("the names chain");
                let got = combinatorics::loop_winding(&jd, NO_CYLS, 0, &edges);
                let want = if reversed { -want } else { want };
                assert_eq!(
                    got.ok(),
                    Some(want),
                    "{name}, reversed = {reversed}, start = {start}"
                );
            }
        }
    }
}

/// A triangular prism over `tri` standing on the plane `z`, `height` tall.
fn triangle_prism(m: &mut Model, tri: [[f64; 2]; 3], z: f64, height: f64) -> Handle<Solid> {
    let frame = datum_frame(
        m,
        crate::SketchPlane::from_axes(
            Point3::from_array([0.0, 0.0, z]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
        ),
    );
    let profile = Profile2d::polygon(
        tri.iter()
            .map(|&q| nacre_math::Point2::from_array(q))
            .collect(),
    )
    .expect("a triangle");
    let Ok(crate::OpOutput::Extrude { solid, .. }) = crate::apply(
        m,
        &crate::Operation::Extrude {
            frame,
            profile,
            dist: height,
        },
    ) else {
        panic!("extruding a triangle")
    };
    m.rebuild_adjacency();
    solid
}

/// Whether some face of `solid` has a loop that passes the vertex at `p` twice.
fn a_loop_passes_twice(m: &Model, solid: Handle<Solid>, p: [f64; 3]) -> bool {
    let shell = m.solid(solid).outer;
    m.shell(shell).faces.iter().any(|&fh| {
        let f = m.face(fh);
        std::iter::once(&f.outer).chain(f.inner.iter()).any(|l| {
            let at_p: Vec<_> = l
                .half_edges
                .iter()
                .map(|he| {
                    let [a, b] = m.edge(he.edge).vertices;
                    if he.forward { a } else { b }
                })
                .filter(|&v| m.vertex_point(v).as_array() == p)
                .collect();
            at_p.len() == 2 && at_p[0] == at_p[1]
        })
    })
}

/// **A plate pinched where a hole meets a boss builds, and its top meshes.** A `6 × 6 × 1`
/// plate, a triangular hole `(2,2) (4,3) (4,4)` through it, and a triangular boss of height 1 on
/// its top with a corner at `(2, 2)` too. The arrangement names `(2, 2, 1)` once, and the top
/// face's boundary passes it twice — with the boss inside the plate, one hole ring around both
/// triangles (one of its two sectors at `(2, 2)` is 278°); with a boss corner on the plate's edge
/// at `(4, 0)`, the outer ring, which passes `(4, 0)` twice as well. Every operation order and
/// hole height is asked; the mesh must cover each face and enclose the solid (`mesh_covers_faces`),
/// as the boolean's own exit census also asserts. Boss first with its corner on the edge makes the
/// hole's boolean take a pinched face as its input — see
/// [`a_pinched_body_is_an_operand_like_any_other`].
#[test]
fn a_plate_pinched_where_a_hole_meets_a_boss_meshes() {
    let hole = [[2.0, 2.0], [4.0, 3.0], [4.0, 4.0]];
    let inside = [[2.0, 2.0], [4.0, 0.5], [4.0, 1.0]];
    let on_edge = [[2.0, 2.0], [4.0, 0.0], [4.0, 1.0]];
    for (name, boss, hole_floor, hole_height, boss_first) in [
        ("boss inside, hole through", inside, -1.0, 3.0, false),
        ("boss inside, hole flush", inside, 0.0, 1.0, false),
        ("boss inside, boss first", inside, -1.0, 3.0, true),
        ("boss on the edge, hole through", on_edge, -1.0, 3.0, false),
        ("boss on the edge, hole flush", on_edge, 0.0, 1.0, false),
        ("boss on the edge, boss first", on_edge, -1.0, 3.0, true),
    ] {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([6.0, 6.0, 1.0]),
        );
        let one = |out: Vec<Handle<Solid>>| {
            assert_eq!(out.len(), 1, "{name}: one body");
            out[0]
        };
        let out = if boss_first {
            let b = triangle_prism(&mut m, boss, 1.0, 1.0);
            let bossed = one(crate::boolean(&mut m, crate::BoolKind::Fuse, plate, b).expect(name));
            m.rebuild_adjacency();
            let h = triangle_prism(&mut m, hole, hole_floor, hole_height);
            one(crate::boolean(&mut m, crate::BoolKind::Cut, bossed, h).expect(name))
        } else {
            let h = triangle_prism(&mut m, hole, hole_floor, hole_height);
            let holed = one(crate::boolean(&mut m, crate::BoolKind::Cut, plate, h).expect(name));
            m.rebuild_adjacency();
            let b = triangle_prism(&mut m, boss, 1.0, 1.0);
            one(crate::boolean(&mut m, crate::BoolKind::Fuse, holed, b).expect(name))
        };
        m.rebuild_adjacency();
        assert!(
            a_loop_passes_twice(&m, out, [2.0, 2.0, 1.0]),
            "{name}: the top's boundary passes (2, 2, 1) twice, as one vertex"
        );
        mesh_covers_faces(name, &m, &[out]);
    }
}

/// **A pinched body is an operand like any other.** The plate with a boss corner on its edge
/// (`P1`: the top face's one ring passes `(4, 0, 1)` twice) and the same plate with a bore through
/// it (`P3`, so a circle's chord stations share the top face's line with the pinch), each against
/// boxes placed where the pinch matters: far from it, around it, with a wall `x = 4` through it,
/// wholly below the top (`z ≤ 0.5`) — there the top's class is out of the box's reach, so the
/// pinched face passes through unarranged and the debug build compares the two routes' faces —
/// and standing on the boss's face `x = 4`, touching it, where the pinched face's cells merge back
/// around the pinch (on the top, and on the underside of a plate whose boss hangs below it).
/// Every kind, both orders: the result validates, its mesh covers it, and its volume is the one
/// written by hand from `A`, `B` and `A ∩ B` (`Fuse = A + B − I`, `Common = I`, `A − B = A − I`,
/// `B − A = B − I`).
///
/// Every line through the pinch is a class here — the boss's walls `x = 4` and `x + y = 4` and
/// the plate's `y = 0` — and on each the two passes of the ring leave one name twice at one
/// station (a run's end and a touch, or two runs meeting end to start). That is the reading
/// [`crate::arrangement`]'s line sweep exists for.
#[test]
fn a_pinched_body_is_an_operand_like_any_other() {
    use crate::BoolKind::{Common, Cut, Fuse};
    let pt = |a: [f64; 3]| Point3::from_array(a);
    let quarter_bore = std::f64::consts::PI / 4.0;
    #[derive(Clone, Copy, PartialEq)]
    enum Body {
        Top,
        Bored,
        Under,
        Holed,
    }
    use Body::{Bored, Holed, Top, Under};
    // The plate `[0, 6]² × [0, 1]` with the triangle `(2,2) (4,0) (4,1)` of area 1 standing on it,
    // height 1 — or hanging under it (`Under`, so the pinched face faces `−z`) — and, for `P3`
    // (`Bored`), a bore of radius ½ at `(1.5, 4)` through the plate alone. `Holed` pinches a hole
    // ring instead: the boss `(2,2) (4,0.5) (4,1)` of area ½ stands inside the plate and the hole
    // `(2,2) (4,3) (4,4)` of area 1 meets it at `(2, 2)`, so the top's one hole ring passes there
    // twice.
    let pinched = |m: &mut Model, body: Body| -> Handle<Solid> {
        let plate = crate::fixtures::cuboid(m, pt([0.0, 0.0, 0.0]), pt([6.0, 6.0, 1.0]));
        let floor = if body == Under { -1.0 } else { 1.0 };
        let corner = if body == Holed { 0.5 } else { 0.0 };
        let boss = triangle_prism(m, [[2.0, 2.0], [4.0, corner], [4.0, 1.0]], floor, 1.0);
        let mut s = crate::boolean(m, Fuse, plate, boss).expect("the boss")[0];
        m.rebuild_adjacency();
        if body == Bored {
            let bore = crate::fixtures::cylinder(
                m,
                pt([1.5, 4.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                3.0,
            )
            .solid;
            s = crate::boolean(m, Cut, s, bore).expect("the bore")[0];
            m.rebuild_adjacency();
        }
        if body == Holed {
            let hole = triangle_prism(m, [[2.0, 2.0], [4.0, 3.0], [4.0, 4.0]], -1.0, 3.0);
            s = crate::boolean(m, Cut, s, hole).expect("the hole")[0];
            m.rebuild_adjacency();
        }
        let pinch = match body {
            Under => [4.0, 0.0, 0.0],
            Holed => [2.0, 2.0, 1.0],
            Top | Bored => [4.0, 0.0, 1.0],
        };
        assert!(
            a_loop_passes_twice(m, s, pinch),
            "the operand's pinched face passes {pinch:?} twice"
        );
        s
    };
    // (name, body, box, A, B, A ∩ B)
    let rows = [
        // The box's corner of the plate: 0.25 × 1.
        (
            "far",
            Top,
            [[5.0, 5.0, -1.0], [5.5, 5.5, 3.0]],
            37.0,
            1.0,
            0.25,
        ),
        // Plate 1 × ½ × ½, and the boss over x ∈ [3.5, 4], y ≤ ½: ∫(x − 3.5) dx = 1/8, times ½.
        (
            "around",
            Top,
            [[3.5, -0.5, 0.5], [4.5, 0.5, 1.5]],
            37.0,
            1.0,
            0.3125,
        ),
        // Plate 1 × ½ × ½; the boss has no volume beyond x = 4.
        (
            "wall x = 4",
            Top,
            [[4.0, -1.0, 0.5], [5.0, 0.5, 1.5]],
            37.0,
            1.5,
            0.25,
        ),
        // Plate 0.25 × ½.
        (
            "below the top",
            Top,
            [[5.0, 5.0, -1.0], [5.5, 5.5, 0.5]],
            37.0,
            0.375,
            0.125,
        ),
        (
            "bored, far",
            Bored,
            [[5.0, 5.0, -1.0], [5.5, 5.5, 3.0]],
            37.0 - quarter_bore,
            1.0,
            0.25,
        ),
        (
            "bored, around",
            Bored,
            [[3.5, -0.5, 0.5], [4.5, 0.5, 1.5]],
            37.0 - quarter_bore,
            1.0,
            0.3125,
        ),
        // A box sharing the boss's face `x = 4` and standing on the top: they touch, nothing
        // overlaps. The top's cells on `z = 1` merge back around the pinch.
        (
            "on the boss's face, corner at the pinch",
            Top,
            [[4.0, 0.0, 1.0], [5.0, 1.0, 2.0]],
            37.0,
            1.0,
            0.0,
        ),
        (
            "on the boss's face, off the pinch",
            Top,
            [[4.0, 0.5, 1.0], [5.0, 1.5, 2.0]],
            37.0,
            1.0,
            0.0,
        ),
        // The same under the plate, where the pinched face faces `−z`.
        (
            "under, on the boss's face, corner at the pinch",
            Under,
            [[4.0, 0.0, -1.0], [5.0, 1.0, 0.0]],
            37.0,
            1.0,
            0.0,
        ),
        (
            "under, on the boss's face, off the pinch",
            Under,
            [[4.0, 0.5, -1.0], [5.0, 1.5, 0.0]],
            37.0,
            1.0,
            0.0,
        ),
        // The pinched hole ring merged back: a box standing on the top away from it, and one
        // against the boss's face `x = 4`.
        (
            "holed, on the top",
            Holed,
            [[0.5, 0.5, 1.0], [1.5, 1.5, 2.0]],
            35.5,
            1.0,
            0.0,
        ),
        (
            "holed, on the boss's face",
            Holed,
            [[4.0, 0.5, 1.0], [5.0, 1.0, 2.0]],
            35.5,
            0.5,
            0.0,
        ),
    ];
    for (name, body, [lo, hi], a_vol, b_vol, i_vol) in rows {
        for (kind, swapped, want) in [
            (Fuse, false, a_vol + b_vol - i_vol),
            (Fuse, true, a_vol + b_vol - i_vol),
            (Common, false, i_vol),
            (Common, true, i_vol),
            (Cut, false, a_vol - i_vol),
            (Cut, true, b_vol - i_vol),
        ] {
            let what = format!("{name}, {kind:?}, swapped = {swapped}");
            let mut m = Model::new();
            let a = pinched(&mut m, body);
            let b = crate::fixtures::cuboid(&mut m, pt(lo), pt(hi));
            m.rebuild_adjacency();
            let (x, y) = if swapped { (b, a) } else { (a, b) };
            let out =
                crate::boolean(&mut m, kind, x, y).unwrap_or_else(|e| panic!("{what}: {e:?}"));
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{what}: {vs:?}");
            mesh_covers_faces(&what, &m, &out);
            let v: f64 = out
                .iter()
                .map(|&s| {
                    nacre_props::mass_props(&m, s)
                        .expect("mass properties")
                        .volume
                })
                .sum();
            assert!((v - want).abs() < 1e-9, "{what}: volume {v}, want {want}");
        }
    }
}

/// **A face pinched at its lowest point meshes on a mirrored chart.** A triangular plate
/// `(0,0) (4,4) (−2,4)` over `z ∈ [0, 1]` with a triangular boss `(0,0) (1,2) (−0.5,2)` under it
/// over `z ∈ [−1, 0]`, sharing the plate's corner: the plate's underside is the triangle less the
/// boss's footprint, so its outer ring passes `(0, 0, 0)` twice, and that corner is the ring's
/// lowest point. The underside faces `−z`, so its loop is clockwise in the chart the mesher
/// projects it to and the chart must swap its axes — a decision `ring_orientation` makes at that
/// very point. Which pass comes first is not this test's to vary: the arrangement starts this
/// loop at the pinched corner with the reflex pass first (four plate and boss shapes, measured), so
/// a reading of the first pass alone would mesh it too. That the answer does not depend on the
/// order is held where the rings are built by hand
/// (`a_ring_pinched_at_its_extreme_point_winds_by_its_first_ray`, every start); this holds the
/// chart's road. (A hole through the plate with its apex at the corner is not this shape: the two
/// wedges of plate left on either side would meet along the corner's vertical edge, and the
/// boolean refuses that as `NonManifoldResultEdge`.)
#[test]
fn a_face_pinched_at_its_lowest_point_meshes_on_a_mirrored_chart() {
    let mut m = Model::new();
    let plate = triangle_prism(&mut m, [[0.0, 0.0], [4.0, 4.0], [-2.0, 4.0]], 0.0, 1.0);
    let boss = triangle_prism(&mut m, [[0.0, 0.0], [1.0, 2.0], [-0.5, 2.0]], -1.0, 1.0);
    let out = crate::boolean(&mut m, crate::BoolKind::Fuse, plate, boss).expect("the union");
    assert_eq!(out.len(), 1, "one body");
    m.rebuild_adjacency();
    assert!(
        a_loop_passes_twice(&m, out[0], [0.0, 0.0, 0.0]),
        "the underside passes its corner twice, as one vertex"
    );
    mesh_covers_faces("triangular plate pinched at its corner", &m, &out);
}

/// **An antenna at the extreme node cancels — its tip or its base.** Same plane table as
/// [`a_ring_pinched_at_its_extreme_node_winds_by_its_outermost_passage`]: a clockwise triangle
/// `P (4,4) (4,3)` with an antenna out of `P` to `(4,1)` and back — the antenna's two departures
/// at `P` come on different passages (out on one, back on the next), and once they cancel the
/// triangle's corner reads `−1`. A ring that is nothing but antennas — `(4,0) → (4,1) → P → (4,1)
/// → (4,0)` — leaves no departure at `P` at all: no side, `StraightAngle`.
#[test]
fn an_antenna_at_the_extreme_node_cancels() {
    let mut m = Model::new();
    let at = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
    let wall = |m: &mut Model, a: (f64, f64), b: (f64, f64)| {
        plane_through(m, [at(a.0, a.1, 0.0), at(b.0, b.1, 0.0), at(a.0, a.1, 1.0)])
    };
    let top = plane_through(
        &mut m,
        [at(0.0, 0.0, 1.0), at(1.0, 0.0, 1.0), at(0.0, 1.0, 1.0)],
    );
    let planes = vec![
        top,                                  // 0: z = 1
        wall(&mut m, (2.0, 2.0), (4.0, 3.0)), // 1: through P along (2, 1)
        wall(&mut m, (2.0, 2.0), (4.0, 4.0)), // 2: through P along (1, 1)
        wall(&mut m, (2.0, 2.0), (4.0, 1.0)), // 3: through P along (2, −1)
        wall(&mut m, (2.0, 2.0), (4.0, 0.0)), // 4: through P along (1, −1)
        wall(&mut m, (4.0, 0.0), (4.0, 1.0)), // 5: x = 4
    ];
    let jd = crate::planes::test_judge(&planes);
    let sorted = |ring: &[[usize; 3]]| -> Vec<[usize; 3]> {
        ring.iter()
            .map(|t| {
                let mut t = *t;
                t.sort_unstable();
                t
            })
            .collect()
    };
    let wind = |r: &[[usize; 3]]| {
        let edges = combinatorics::ring_from_names(0, r).expect("the names chain");
        combinatorics::loop_winding(&jd, NO_CYLS, 0, &edges)
    };
    let based = sorted(&[[0, 3, 2], [0, 2, 5], [0, 5, 1], [0, 1, 3], [0, 3, 5]]);
    for start in 0..based.len() {
        let mut r = based.clone();
        r.rotate_left(start);
        assert_eq!(
            wind(&r).ok(),
            Some(-1),
            "antenna based at P, start = {start}"
        );
    }
    let only = sorted(&[[0, 5, 4], [0, 5, 3], [0, 3, 1], [0, 3, 5]]);
    for start in 0..only.len() {
        let mut r = only.clone();
        r.rotate_left(start);
        assert!(
            matches!(
                wind(&r),
                Err(BoolError::Rejected {
                    reason: RejectReason::StraightAngle,
                    ..
                })
            ),
            "antennas only, start = {start}: {:?}",
            wind(&r)
        );
    }
}

#[test]
fn a_reflex_node_turns_against_its_ring() {
    // The L cap's outer ring is a hexagon with exactly one reflex corner, at `(1, 1)`. A
    // convex node turns with the ring and the reflex one turns against it, so the turn signs
    // are *not* all equal — which is why the winding cannot be read off an arbitrary node.
    let (planes, p, outer, _) = holed_face_rings("dimple");
    let pts = ring_points(&planes, &outer);
    let reflex = pts
        .iter()
        .position(|q| near(Point3::from_array(*q), [1.0, 1.0, 1.0]))
        .expect("the reflex node");
    assert_eq!(
        combinatorics::turn_at(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap(),
            reflex
        )
        .unwrap(),
        -1
    );
    let turns: Vec<i8> = (0..outer.len())
        .map(|i| {
            combinatorics::turn_at(
                &crate::planes::test_judge(&planes),
                NO_CYLS,
                p,
                &combinatorics::ring_from_names(p, &outer).unwrap(),
                i,
            )
            .unwrap()
        })
        .collect();
    assert_eq!(
        turns.iter().filter(|&&t| t == -1).count(),
        1,
        "one reflex corner: {turns:?}"
    );

    // The node `loop_winding` lands on is the lexicographically least — a hull vertex, where
    // the turn *is* the winding. The test finds it by reading coordinates; `loop_winding`
    // finds it with an exact predicate.
    let lo = (0..pts.len())
        .min_by(|&i, &j| pts[i].partial_cmp(&pts[j]).unwrap())
        .unwrap();
    assert_ne!(lo, reflex);
    assert_eq!(
        combinatorics::turn_at(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap(),
            lo
        )
        .unwrap(),
        1
    );

    // The teeth. A ring is a cycle, so its winding cannot depend on where the walk began. Start it
    // at the reflex node and a `turn_at(ring[0])` implementation reads the reflex sign. Measured:
    // without this rotation, such an implementation passes every assertion above.
    let mut rotated = outer.clone();
    rotated.rotate_left(reflex);
    assert_eq!(
        combinatorics::turn_at(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &rotated).unwrap(),
            0
        )
        .unwrap(),
        -1
    );
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &rotated).unwrap()
        )
        .unwrap(),
        1
    );
}

#[test]
fn a_loop_is_inside_the_face_it_was_found_on() {
    // A hole ring never touches its face's outer ring: every node is *strictly* inside, so
    // `point_in_ring` must say so for all of them. The outer ring is the L's cap, a hexagon
    // with a reflex corner, so this is not a convex test.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    for t in &hole {
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                &[],
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
    }
}

#[test]
fn square_and_l_holes_are_inside_their_faces_too() {
    // Two hole shapes that must not move: a square one and an L-shaped one. Both sit
    // strictly inside the same reflex hexagon, and **every** clear ray agrees — the parity
    // cannot depend on which ray was cast, which is a second machine for free.
    for which in ["dimple", "ell"] {
        let (planes, p, outer, hole) = holed_face_rings(which);
        for t in &hole {
            let rays = combinatorics::every_ray(
                &crate::planes::test_judge(&planes),
                &[],
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap(),
            )
            .unwrap();
            assert!(!rays.is_empty(), "{which}: no clear ray");
            assert!(rays.iter().all(|&x| x), "{which}: {rays:?}");
        }
    }
}

#[test]
fn a_loop_is_placed_by_where_it_is_not_by_how_it_winds() {
    // Containment is about **where** a ring is, never which way it runs. A hole ring is stored
    // clockwise about the face normal and its outer ring counter-clockwise, and neither
    // direction may enter the answer: reversing either must change nothing.
    //
    // And containment is **not symmetric** — the classic way to get this wrong is a test that
    // only ever asks it one way round.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    let (mut rev_outer, mut rev_hole) = (outer.clone(), hole.clone());
    rev_outer.reverse();
    rev_hole.reverse();
    for t in &hole {
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                &[],
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                &[],
                p,
                *t,
                &combinatorics::ring_from_names(p, &rev_outer).unwrap()
            )
            .unwrap(),
            "reversing the outer ring must not move the hole"
        );
    }
    for t in &outer {
        assert!(
            !combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                &[],
                p,
                *t,
                &combinatorics::ring_from_names(p, &hole).unwrap()
            )
            .unwrap()
        );
        assert!(
            !combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                &[],
                p,
                *t,
                &combinatorics::ring_from_names(p, &rev_hole).unwrap()
            )
            .unwrap(),
            "nor may reversing the hole swallow the outer ring"
        );
    }
}

#[test]
fn every_clear_ray_agrees() {
    // The ring is simple, so the parity cannot depend on the ray. `every_ray` returns one
    // answer per usable candidate and they must be unanimous; a disagreement means the ray
    // choice leaked into the result.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    let rays = combinatorics::every_ray(
        &crate::planes::test_judge(&planes),
        &[],
        p,
        hole[0],
        &combinatorics::ring_from_names(p, &outer).unwrap(),
    )
    .unwrap();
    assert!(!rays.is_empty(), "at least one candidate is clear");
    assert!(rays.iter().all(|&x| x), "and they agree: inside — {rays:?}");
}

#[test]
fn a_ring_inside_a_ring_is_what_nesting_looks_like() {
    // Nested loops on one face have no operand in the suite that produces them — a polyhedral
    // torus would.
    // The detector can still be aimed at real geometry: a holed face *is* a ring inside a ring,
    // which fires exactly the condition the nesting brick asks about. It is the detector under
    // test, not the fixture.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    assert!(
        combinatorics::point_in_ring(
            &crate::planes::test_judge(&planes),
            &[],
            p,
            hole[0],
            &combinatorics::ring_from_names(p, &outer).unwrap()
        )
        .unwrap()
    );
    assert!(
        !combinatorics::point_in_ring(
            &crate::planes::test_judge(&planes),
            &[],
            p,
            outer[0],
            &combinatorics::ring_from_names(p, &hole).unwrap()
        )
        .unwrap()
    );
}

/// A 40×40 plate with the given outline and a through pocket at `[20,30]×[10,20]`, as the
/// holed cap's two rings. The pocket's four wall **planes** are what the fixtures below aim
/// rays along; the outline decides which of them a ring node sits on.
fn grazed_plate(outline: &[[f64; 2]]) -> (Vec<WorkingPlane>, usize, Ring, Ring) {
    let profile = Profile2d::polygon(outline.iter().map(|&p| p2(p[0], p[1])).collect()).unwrap();
    let mut m = replay(&[extrude_log_op(profile, 12.0)]).unwrap();
    let plate = m.live_solids()[0];
    // Through, so the cap really is holed rather than dimpled.
    let pocket = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([20.0, 10.0, -1.0]),
        Point3::from_array([30.0, 20.0, 13.0]),
    );
    holed_face_rings_of(m, plate, pocket)
}

/// The ring node whose point is `(x, y)` — fixtures are written in coordinates and the engine
/// answers in plane triples, so this is where the two meet.
fn node_at(planes: &[WorkingPlane], ring: &Ring, x: f64, y: f64) -> [usize; 3] {
    let pts = ring_points(planes, ring);
    let i = pts
        .iter()
        .position(|p| (p[0] - x).abs() < 1e-9 && (p[1] - y).abs() < 1e-9)
        .unwrap_or_else(|| panic!("no node at ({x}, {y}) among {pts:?}"));
    ring[i]
}

/// **A ray that grazes a corner is still a ray** — the four ways a ring can meet the ray's
/// line, each with an answer known by hand.
///
/// A ring node *on* the line makes the parity ambiguous; throwing such a candidate away, and
/// both of a vertex's candidates with it, sends the whole question back as `NoClearRay`, and a
/// band of rotation angles dies of it. The rule that resolves it is the one
/// `trace_transversal_face` uses: look at the node's two off-line neighbours —
/// **opposite sides is a crossing, equal sides a touch**.
///
/// So both fixtures block **both** candidates. That matters: with one candidate left clear the
/// answer comes out anyway and the test would be green with or without the rule, measuring
/// nothing. Here `point_in_ring` is `NoClearRay` without the rule and the truth with it.
///
/// The outer outlines are chosen so that no vertex sits mid-run on a straight edge —
/// `Profile2d` dissolves those at construction — and the ring lengths are asserted so that a
/// dissolve fails loudly instead of quietly weakening the fixture.
#[test]
fn a_ray_that_grazes_a_corner_still_answers() {
    // The pocket's own corner `(20,10)`: its two candidate lines are `y=10` and `x=20`.
    //
    // `y=10` carries the run — the outline's `(5,10)→(0,10)` edge lies *on* it — and the run's
    // flanks `(5,20)` and `(0,0)` are on opposite sides, so the boundary crosses there.
    // `x=20` is grazed by the single corner `(20,30)`, whose flanks `(40,40)` and `(0,40)` are
    // also opposite. One crossing each way: the pocket corner is inside the plate.
    let a = [
        [0.0, 0.0],
        [40.0, 0.0],
        [40.0, 40.0],
        [20.0, 30.0],
        [0.0, 40.0],
        [0.0, 20.0],
        [5.0, 20.0],
        [5.0, 10.0],
        [0.0, 10.0],
    ];
    let (planes, p, outer, hole) = grazed_plate(&a);
    assert_eq!(outer.len(), 9, "every corner of the outline survived");
    let jd = crate::planes::test_judge(&planes);
    let outer_edges = combinatorics::ring_from_names(p, &outer).unwrap();
    let hole_edges = combinatorics::ring_from_names(p, &hole).unwrap();

    // Run-crossing and isolated-crossing, both saying "inside".
    let v = node_at(&planes, &hole, 20.0, 10.0);
    assert!(
        combinatorics::point_in_ring(&jd, &[], p, v, &outer_edges).unwrap(),
        "the pocket's corner is inside the plate"
    );
    let rays = combinatorics::every_ray(&jd, &[], p, v, &outer_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(rays.iter().all(|&x| x), "and unanimous — {rays:?}");

    // Negative control, and the touch rule's own lock: from the outline's `(0,10)` the
    // `y=10` line grazes the pocket's `(20,10)→(30,10)` edge, whose flanks `(30,20)` and
    // `(20,20)` are on the *same* side. Nothing crossed, so `(0,10)` is outside the pocket —
    // and counting that touch as a crossing would make this ray disagree with the `x=0` one.
    let w = node_at(&planes, &outer, 0.0, 10.0);
    assert!(
        !combinatorics::point_in_ring(&jd, &[], p, w, &hole_edges).unwrap(),
        "a plate corner is not inside the pocket"
    );
    let rays = combinatorics::every_ray(&jd, &[], p, w, &hole_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(!rays.iter().any(|&x| x), "and unanimous — {rays:?}");
}

/// The isolated **touch** — the branch fixture A never reaches.
///
/// The outline's bottom notch rises to `(12,10)` and turns straight back down, so both its
/// neighbours `(8,0)` and `(16,0)` are below `y=10`: the ring touched the ray's line without
/// crossing it. Counting it as a crossing flips the parity of the `−x` ray alone, and the two
/// directions of one candidate would then contradict each other.
#[test]
fn a_ring_that_touches_the_ray_and_turns_back_crossed_nothing() {
    let b = [
        [0.0, 0.0],
        [8.0, 0.0],
        [12.0, 10.0], // touches y=10 and turns back — the branch under test
        [16.0, 0.0],
        [40.0, 0.0],
        [40.0, 40.0],
        [20.0, 30.0], // blocks the x=20 candidate, so neither is clear
        [0.0, 40.0],
        [0.0, 20.0],
        [5.0, 10.0], // crosses y=10
        [0.0, 5.0],
    ];
    let (planes, p, outer, hole) = grazed_plate(&b);
    assert_eq!(outer.len(), 11, "every corner of the outline survived");
    let jd = crate::planes::test_judge(&planes);
    let outer_edges = combinatorics::ring_from_names(p, &outer).unwrap();
    let v = node_at(&planes, &hole, 20.0, 10.0);
    assert!(
        combinatorics::point_in_ring(&jd, &[], p, v, &outer_edges).unwrap(),
        "the pocket's corner is inside the plate"
    );
    let rays = combinatorics::every_ray(&jd, &[], p, v, &outer_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(rays.iter().all(|&x| x), "and unanimous — {rays:?}");
}
