use super::*;
use nacre_exact::{Angle, Axis, Rat};

/// ★★★★★ **The same circle and the same ruling, and the arc decides.**
///
/// The unit-ish circle of radius `6/5` about the origin in `z = 0`, and a cylinder along `x`
/// through the origin with radius `1` — so that plane cuts it in the two rulings `y = ±1`, and
/// the extent every piece is measured along is `e = n × m = ŷ`.
///
/// The whole circle reaches `y ∈ [−6/5, 6/5]` and holds both rulings. The **left** arc — from
/// `(−24/25, 18/25)` to `(−24/25, −18/25)`, counter-clockwise about `+ẑ`, so it runs through
/// `−x̂` and holds neither `±ŷ` — reaches only `[−18/25, 18/25]` and holds neither. A **top**
/// arc holds `+ŷ`, so it reaches `+6/5` and holds the `+1` ruling but not the `−1`.
///
/// Asking the circle instead would refuse a plate whose
/// corner fillets never come near the drill's rulings.
#[test]
fn the_arc_decides_which_ruling_a_circle_holds() {
    use nacre_exact::{MeetPoint, cylinder_ruling_reached_extent};
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let coeffs = [z(0), z(0), z(1), z(0)]; // the plane `z = 0`
    let axis = [z(0), z(0), z(1)]; // the circle's carrier
    let (o, m, r) = (
        [z(0); 3],
        [z(1), z(0), z(0)],
        nacre_exact::BigRat::from(z(1)),
    ); // the ruling cylinder, along `x`
    let e = [z(0), z(1), z(0)]; // `n × m`
    let centre = [z(0); 3];
    let rho2 = nacre_exact::BigRat::from(q(36, 25)); // ρ = 6/5, stated as its square
    let holds = |arc: Option<&RimArc>, side: i8| {
        let (lo, hi) = arc_ends_along(&centre, &rho2, &axis, arc, &e).expect("an extent");
        cylinder_ruling_reached_extent(
            &coeffs,
            &nacre_exact::StripReach {
                lo: (&MeetPoint::Narrow(lo.0), &lo.1),
                hi: Some((&MeetPoint::Narrow(hi.0), &hi.1)),
            },
            &o,
            &m,
            &r,
            side,
            // Production's question: a crossing, not a touch. Every case below is
            // clear of the boundary either way, which is why the boundary has its own lock.
            false,
        )
    };
    // The whole circle spans `±6/5` and holds both.
    for side in [-1i8, 1] {
        assert!(holds(None, side), "the whole circle, side {side}");
    }
    // The left arc spans `±18/25` and holds neither — the arc's own reach, not its circle's.
    let left = RimArc {
        from: [q(-24, 25), q(18, 25), z(0)],
        to: [q(-24, 25), q(-18, 25), z(0)],
    };
    for side in [-1i8, 1] {
        assert!(!holds(Some(&left), side), "the left arc, side {side}");
    }
    // A top arc holds `+ŷ`, so it reaches `+6/5` — and still never the other ruling.
    let top = RimArc {
        from: [q(24, 25), q(18, 25), z(0)],
        to: [q(-24, 25), q(18, 25), z(0)],
    };
    assert!(holds(Some(&top), 1), "the top arc reaches the `+1` ruling");
    assert!(!holds(Some(&top), -1), "and never the `−1` one");
}

/// ★★★★★ **An arc's reach, by hand, and which arc it is.**
///
/// The unit circle about the origin in the plane `z = 0`, axis `+ẑ`, and the **first quadrant**
/// as `from = (1,0) → to = (0,1)` (counter-clockwise about that axis). Its points are
/// `(cos θ, sin θ)` for `θ ∈ [0, π/2]`, so along `d = (1,1,0)` it reaches
/// `cos θ + sin θ ∈ [1, √2]` — **one end rational, the other a radical**, which is the shape
/// no symmetric margin can state and the reason this rule exists.
///
/// ★ The last block is the lock that matters most: swapping the ends names the **complementary**
/// arc, whose reach is genuinely different. Reading the face's normal instead of the carrier's
/// axis would make exactly that swap.
#[test]
fn an_arcs_reach_is_its_own_and_not_its_complements() {
    let z = |v: i128| Rat::from_int(v);
    let one = nacre_exact::BigRat::from(z(1)); // r² = 1
    let axis = [z(0), z(0), z(1)];
    let quadrant = RimArc {
        from: [z(1), z(0), z(0)],
        to: [z(0), z(1), z(0)],
    };
    // Along `+x̂`: the peak is on the arc (it is `from`), so that end carries the radical and
    // the other stops at an end of the arc — `[0, 1]`, the quadrant's own span in `x`.
    assert_eq!(
        arc_extent(Some(&quadrant), &one, &axis, &[z(1), z(0), z(0)]),
        Some((z(0), z(0), z(0), z(1)))
    );
    // Along the diagonal: `[1, √2]`. `ρ² = r²|d⊥|² = 2`, and the low end is the rational `1`.
    assert_eq!(
        arc_extent(Some(&quadrant), &one, &axis, &[z(1), z(1), z(0)]),
        Some((z(1), z(0), z(0), z(2)))
    );
    // A whole circle reaches alike both ways, and `d ∥ axis` has no radial term at all.
    assert_eq!(
        arc_extent(None, &one, &axis, &[z(1), z(0), z(0)]),
        Some((z(0), z(1), z(0), z(1)))
    );
    assert_eq!(
        arc_extent(Some(&quadrant), &one, &axis, &axis),
        Some((z(0), z(0), z(0), z(0)))
    );
    // ★ The complement is a different set: `[−1, 1]`, not `[0, 1]`.
    let rest = RimArc {
        from: quadrant.to,
        to: quadrant.from,
    };
    assert_eq!(
        arc_extent(Some(&rest), &one, &axis, &[z(1), z(0), z(0)]),
        Some((z(0), z(1), z(0), z(1)))
    );
}

/// **The pierce tolerance is a measurement, and each of its terms can move.**
///
/// Hand geometry — no class indices, so nothing here is copied from an engine run: the turned
/// boss's second crossing, where the meet line of `y = 0` and `x = 4` pierces a cylinder of
/// radius `0.5` about the `+X` axis through `(y, z) = (0.25, 2)`. The true point is
/// `(4, 0, 2 − √3/4)` (irrational on purpose — every term sees real rounding, not exact zeros).
///
/// ★ The second probe moves **along the meet line**: both planes and the line stay at zero, so
/// only the cylinder-surface term can see it — which is how this asserts that term is
/// *exercised*, not merely present ([[instrument-decides-the-answer]]'s negative control).
#[test]
fn pierce_vertex_tol_measures_and_each_term_moves() {
    let pa = Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .unwrap();
    let pb = Plane::from_point_normal(
        Point3::from_array([4.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
    )
    .unwrap();
    let cyl = nacre_geom::Cylinder::from_axis(
        Point3::from_array([0.0, 0.25, 2.0]),
        nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
        0.5,
    )
    .unwrap();
    let s = 2.0 - 3.0f64.sqrt() / 4.0;
    let p = Point3::from_array([4.0, 0.0, s]);
    assert!(
        pierce_vertex_tol(p, &pa, &pb, &cyl) < 1e-12,
        "the true crossing measures at rounding scale: {}",
        pierce_vertex_tol(p, &pa, &pb, &cyl)
    );
    // Along the meet line: planes and line stay zero, the cylinder term alone answers.
    let along = Point3::from_array([4.0, 0.0, s - 1e-6]);
    assert!(
        pierce_vertex_tol(along, &pa, &pb, &cyl) > 1e-7,
        "the cylinder-surface term is exercised"
    );
    // Off a plane: the instrument moves by the full offset.
    let off = Point3::from_array([4.0, 1e-6, s]);
    assert!(
        pierce_vertex_tol(off, &pa, &pb, &cyl) > 0.9e-6,
        "a plane term is exercised"
    );
}

/// **A loop with points that do not turn still states its own outward direction.**
///
/// Two faces sharing an edge list the same vertices along it, so a pad that covers part of a
/// face leaves the neighbouring walls with points strung along one line. Those points cannot
/// be removed — that would leave a T-vertex — so the population is every partially-covering
/// pad, and `outer_tri` has to survive it.
///
/// ★ **Turned coordinates are the whole difficulty.** Axis-aligned collinear points cancel to
/// exactly `0.0`, which a "first corner that is not flat" rule skips correctly. Turned, the
/// cancellation leaves ~2⁻⁵³ — not zero, so that rule takes it, and the direction that comes
/// back is rounding: measured 90° off its own plane, which is `orient_sign`, the winding of the
/// witness triangle every predicate borrows, and which side of the plane holds material.
///
/// Asserted on `collect_planes`' own output rather than through `debug_assert`, so it is a
/// measurement in release too. `tests/probes/collinear_loop_points.rs` sweeps the same proposition
/// from outside the crate over a band of angles.
#[test]
fn a_loop_whose_points_do_not_all_turn_still_faces_the_right_way() {
    use crate::{OpOutput, Operation, Profile2d, SketchFrame, apply};
    use nacre_exact::{Isometry, Rotation};
    use nacre_math::Point2;

    let rect = |a: f64, b: f64, c: f64, d: f64| {
        Profile2d::polygon(vec![
            Point2::from_array([a, b]),
            Point2::from_array([c, b]),
            Point2::from_array([c, d]),
            Point2::from_array([a, d]),
        ])
        .expect("a rectangle is a fair profile")
    };

    let mut m = Model::new();
    let world = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: world,
            profile: rect(0.0, 0.0, 2.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("the block") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let OpOutput::Transform { solid } = apply(
        &mut m,
        &Operation::Transform {
            solid,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).expect("angle"),
            }),
        },
    )
    .expect("the turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // A pad covering part of the face: the walls around it inherit the split edge.
    let face = m.shell(m.solid(solid).outer).faces[0];
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face,
            profile: rect(0.25, 0.25, 1.25, 1.25),
            dist: 1.0,
        },
    )
    .expect("the pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    let faces = collect_planes(&m, solid).expect("the planes of a padded block");
    // The population has to be present, or this asserts over nothing.
    let straight = faces
        .iter()
        .filter(|f| {
            let n = m
                .face(f.face().expect("a model face"))
                .outer
                .half_edges
                .len();
            n > 4
        })
        .count();
    assert!(
        straight > 0,
        "no face carries a split edge — the fixture stopped reaching the population"
    );
    // The lock is on the *winding* side now: `n_out` is read off the stored
    // orientation, so `plane.normal()·n_out` is ±1 by construction and asserts
    // nothing. What the widest-corner rule still owns is the triangle — revert
    // `outer_tri` to "first non-flat corner" and the cross below is rounding
    // noise again, failing this in release.
    for row in &faces {
        let f = row.plane();
        let cos = (f.tri[1] - f.tri[0])
            .cross(f.tri[2] - f.tri[0])
            .normalize()
            .expect("a widest corner spans area")
            .dot(f.n_out);
        assert!(
            cos > 0.5,
            "a face's outer triangle is {cos:.3} of the way to perpendicular against its \
                 stated outward"
        );
    }
}

/// **A judgement's headroom is relative to its model, not carved out of a shared ceiling.**
///
/// The two budgets answer different questions — how deep the model is, and how thin a witness
/// it may still judge — and sharing one absolute number silently couples them: a solid turned
/// a thousand times would leave a hard judgement no room, so the same sliver would be judged
/// in a fresh model and abandoned in a turned one. Turning a model must not change what
/// counts as judgeable, so this pins the headroom to a constant *above the model's own
/// precision*, at both ends of the depth range.
#[test]
fn the_climbing_headroom_survives_a_deep_model() {
    let deg = Angle::from_deg(Rat::from_int(37)).expect("angle");
    let mut p = WitnessPoint::at([Rat::from_int(1), Rat::from_int(2), Rat::from_int(3)]);
    let mut seen = Vec::new();
    for turn in 0..=200 {
        if turn == 0 || turn == 20 || turn == 200 {
            let j = standard_for_points(std::slice::from_ref(&p));
            assert_eq!(
                j.cap,
                j.prec + CLIMB_HEADROOM,
                "turn {turn}: headroom is not the model's own precision plus a constant"
            );
            seen.push(j.prec);
        }
        p = p.rotate(Axis::Z, deg);
    }
    // …and the precision really did grow with the history, or the test above is vacuous.
    assert!(
        seen[0] < seen[2],
        "precision did not grow with the rotation history: {seen:?}"
    );
}

/// **Plane classes merge only what interning merged** — over two tables whose different walls
/// share a rounded image: the plane cache's coefficients bit for bit (`rounded_twin_walls`), and
/// the face corners' caches on one `f64` plane (`rounded_corner_walls`).
///
/// The oracle is interning, not the name the judge reads: every plane here is unmoved, and an
/// unmoved plane's canonical name is its interning key, so two handles are two planes and the
/// class count must be the surface count. Either cache, admitted as evidence, merges the two walls
/// and comes up one short.
#[test]
fn plane_classes_merge_only_what_interning_merged() {
    for (case, (m, b, c)) in [
        ("twin caches", crate::tests::rounded_twin_walls()),
        ("corner caches", crate::tests::rounded_corner_walls()),
    ] {
        for (x, y) in [(b, c), (c, b)] {
            let mut faces = collect_planes(&m, x).unwrap();
            faces.extend(collect_planes(&m, y).unwrap());
            let surfaces: std::collections::HashSet<_> =
                faces.iter().map(|f| f.plane().surf).collect();
            assert!(
                surfaces.iter().all(|&s| m.plane_motion(s).is_none()),
                "{case}: the oracle needs unmoved planes"
            );
            let mut roots = plane_classes(&test_judge(&faces));
            roots.sort_unstable();
            roots.dedup();
            assert_eq!(roots.len(), surfaces.len(), "{case}: classes vs surfaces");
        }
    }
}
