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

    // ★ The teeth. A ring is a cycle, so its winding cannot depend on where the walk began.
    // Start it at the reflex node and a `turn_at(ring[0])` implementation reads the reflex
    // sign — the exact fault `outer_tri` shipped. Measured: without this rotation, such an
    // implementation passes every assertion above.
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
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
    }
}

#[test]
fn the_older_loops_are_inside_their_faces_too() {
    // Two hole shapes that must not move: a square one and an L-shaped one. Both sit
    // strictly inside the same reflex hexagon, and **every** clear ray agrees — the parity
    // cannot depend on which ray was cast, which is a second machine for free.
    for which in ["dimple", "ell"] {
        let (planes, p, outer, hole) = holed_face_rings(which);
        for t in &hole {
            let rays = combinatorics::every_ray(
                &crate::planes::test_judge(&planes),
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
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
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
                p,
                *t,
                &combinatorics::ring_from_names(p, &hole).unwrap()
            )
            .unwrap()
        );
        assert!(
            !combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
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
            p,
            hole[0],
            &combinatorics::ring_from_names(p, &outer).unwrap()
        )
        .unwrap()
    );
    assert!(
        !combinatorics::point_in_ring(
            &crate::planes::test_judge(&planes),
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

/// ★ **A ray that grazes a corner is still a ray** — the four ways a ring can meet the ray's
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
        combinatorics::point_in_ring(&jd, p, v, &outer_edges).unwrap(),
        "the pocket's corner is inside the plate"
    );
    let rays = combinatorics::every_ray(&jd, p, v, &outer_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(rays.iter().all(|&x| x), "and unanimous — {rays:?}");

    // ★ Negative control, and the touch rule's own lock: from the outline's `(0,10)` the
    // `y=10` line grazes the pocket's `(20,10)→(30,10)` edge, whose flanks `(30,20)` and
    // `(20,20)` are on the *same* side. Nothing crossed, so `(0,10)` is outside the pocket —
    // and counting that touch as a crossing would make this ray disagree with the `x=0` one.
    let w = node_at(&planes, &outer, 0.0, 10.0);
    assert!(
        !combinatorics::point_in_ring(&jd, p, w, &hole_edges).unwrap(),
        "a plate corner is not inside the pocket"
    );
    let rays = combinatorics::every_ray(&jd, p, w, &hole_edges).unwrap();
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
        combinatorics::point_in_ring(&jd, p, v, &outer_edges).unwrap(),
        "the pocket's corner is inside the plate"
    );
    let rays = combinatorics::every_ray(&jd, p, v, &outer_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(rays.iter().all(|&x| x), "and unanimous — {rays:?}");
}
