use super::*;
use nacre_math::Vector3;
use nacre_predicates::orient2d;
use proptest::prelude::*;

/// The half-open rule's whole table: a step straddles iff exactly one end is strictly
/// above (an end on the ray is not), a straddling step crosses iff the probe is left of it
/// going up or right of it going down, and a straddling step through the probe is `None`.
#[test]
fn the_ray_crossing_rule_is_half_open() {
    use Orient::{Negative as N, Positive as P, Zero as Z};
    for ya in [P, N, Z] {
        for yb in [P, N, Z] {
            let straddles = (ya == P) != (yb == P);
            assert_eq!(
                ray_straddle(ya, yb),
                straddles.then_some(yb == P),
                "{ya:?} {yb:?}"
            );
            for side in [P, N, Z] {
                let want = match (straddles, side) {
                    (false, _) => Some(false),
                    (true, Z) => None,
                    (true, s) => Some((s == P) == (yb == P)),
                };
                assert_eq!(
                    ray_step_crossing(ya, yb, side),
                    want,
                    "{ya:?} {yb:?} {side:?}"
                );
            }
        }
    }
    // The corner on the ray, spelled out: the step leaving it upward counts when the corner
    // is right of the probe (the probe left of the rising step) and not otherwise; the step
    // leaving it downward never counts, nor does a step along the ray — so two steps
    // leaving a corner to opposite sides count once and two leaving to the same side count
    // twice or not at all.
    let up = |side| ray_step_crossing(Z, P, side);
    assert_eq!((up(P), up(N), up(Z)), (Some(true), Some(false), None));
    for side in [P, N, Z] {
        assert_eq!(
            ray_step_crossing(Z, N, side),
            Some(false),
            "downward {side:?}"
        );
        assert_eq!(
            ray_step_crossing(N, Z, side),
            Some(false),
            "arriving from below {side:?}"
        );
        assert_eq!(
            ray_step_crossing(Z, Z, side),
            Some(false),
            "along the ray {side:?}"
        );
    }
}

fn plane(origin: [f64; 3], normal: [f64; 3]) -> Plane {
    Plane::from_point_normal(Point3::from_array(origin), Vector3::from_array(normal)).unwrap()
}

// --- proptest ---

fn vec3() -> impl Strategy<Value = Vector3> {
    prop::array::uniform3(-1.0f64..1.0).prop_map(Vector3::from_array)
}

// --- three_planes ---

#[test]
fn three_planes_axis_gives_unit_vertex() {
    // x = 1, y = 1, z = 1 ⇒ (1, 1, 1).
    let v = three_planes(
        &plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        &plane([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]),
        &plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
    )
    .unwrap();
    assert_eq!(v.as_array(), [1.0, 1.0, 1.0]);
}

#[test]
fn three_planes_parallel_pair_is_none() {
    // x = 0 and x = 1 are parallel ⇒ no vertex, whatever the third plane.
    assert!(
        three_planes(
            &plane([0.0; 3], [1.0, 0.0, 0.0]),
            &plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            &plane([0.0; 3], [0.0, 1.0, 0.0]),
        )
        .is_none()
    );
}

/// Well-conditioned unit normals: three vectors whose triple product is well
/// clear of zero (so the vertex is well-conditioned).
fn three_unit_normals() -> impl Strategy<Value = (Vector3, Vector3, Vector3)> {
    (vec3(), vec3(), vec3()).prop_filter_map("zero/near-coplanar normals", |(a, b, c)| {
        let n1 = a.normalize()?;
        let n2 = b.normalize()?;
        let n3 = c.normalize()?;
        (n1.dot(n2.cross(n3)).abs() >= 0.1).then_some((n1, n2, n3))
    })
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Three well-conditioned planes through a target point recover it.
    #[test]
    fn three_planes_recovers_constructed_vertex(
        p in prop::array::uniform3(-100.0f64..100.0),
        (n1, n2, n3) in three_unit_normals(),
    ) {
        let p = Point3::from_array(p);
        let v = three_planes(
            &Plane::from_point_normal(p, n1).unwrap(),
            &Plane::from_point_normal(p, n2).unwrap(),
            &Plane::from_point_normal(p, n3).unwrap(),
        )
        .unwrap();
        let mag = p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max);
        prop_assert!(v.distance(p) <= 1e-6 * (mag + 1.0));
    }

}

// --- ring_self_intersection_rat ---

/// Lift literal coordinates into the truth the predicates judge. `unwrap` is fine here: every
/// fixture coordinate is a short decimal, comfortably inside the window.
fn ring(pts: &[[f64; 2]]) -> Vec<[Rat; 2]> {
    pts.iter().map(|p| [d(p[0]), d(p[1])]).collect()
}

fn d(x: f64) -> Rat {
    Rat::from_decimal(x).unwrap()
}

/// The net that catches the likeliest way to get this wrong: edges `0` and `n-1` share a point
/// too, so a non-cyclic adjacency test reports *every* ring. Convex, non-convex, and a flat
/// (collinear) corner all have to pass.
#[test]
fn simple_polygons_have_no_self_intersection() {
    let square = ring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection_rat(&square), None);
    let triangle = ring(&[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection_rat(&triangle), None);
    // Reflex corner (an L).
    let l = ring(&[
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 2.0],
        [2.0, 2.0],
        [2.0, 4.0],
        [0.0, 4.0],
    ]);
    assert_eq!(ring_self_intersection_rat(&l), None);
    // A flat corner: `(2,0)` sits mid-run on a straight edge. Collinear but not overlapping.
    let flat = ring(&[[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection_rat(&flat), None);
}

#[test]
fn a_bowtie_crosses_itself() {
    let bowtie = ring(&[[0.0, 0.0], [4.0, 4.0], [4.0, 0.0], [0.0, 4.0]]);
    // Edges 0 (`(0,0)→(4,4)`) and 2 (`(4,0)→(0,4)`) are the crossing pair.
    assert_eq!(ring_self_intersection_rat(&bowtie), Some((0, 2)));
}

/// A touch is as fatal as a crossing: the ring has no strict inside at the pinch point.
#[test]
fn a_pinch_touching_without_crossing_is_rejected() {
    let pinch = ring(&[
        [0.0, 0.0],
        [2.0, 2.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [2.0, 2.0],
        [0.0, 4.0],
    ]);
    assert!(ring_self_intersection_rat(&pinch).is_some());
}

/// A spike shorter than the edge it retraces, and one longer. Neither *isolates* the two-sided
/// adjacency test — with four or more points a long spike puts the ring's start point in the
/// interior of the doubling-back edge, so the non-adjacent rule catches it first. The case that
/// needs both directions is the collinear triple below.
#[test]
fn a_spike_doubling_back_is_rejected_either_length() {
    let short = ring(&[[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 4.0]]);
    assert!(ring_self_intersection_rat(&short).is_some());
    let long = ring(&[[0.0, 0.0], [4.0, 0.0], [-2.0, 0.0], [0.0, 4.0]]);
    assert!(ring_self_intersection_rat(&long).is_some());
}

#[test]
fn a_repeated_consecutive_point_is_a_zero_length_edge() {
    let dup = ring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection_rat(&dup), Some((1, 1)));
    // Also across the wrap: the last point repeats the first.
    let closed = ring(&[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0], [0.0, 0.0]]);
    assert_eq!(ring_self_intersection_rat(&closed), Some((3, 3)));
}

/// A ring of collinear points encloses nothing, and a signed area cannot orient it. It falls
/// out of the same rule — the return leg always overlaps the outbound one.
///
/// **This is the net for the two-sided adjacency test.** A brute force over every ring of 3–5
/// points on a 4×4 grid found the one-sided and two-sided predicates disagreeing on 88 rings,
/// *all* of them collinear triples like this one: at each corner the shared point sits at an
/// end of the retraced span, so only the `u ∈ [s, v]` direction fires.
#[test]
fn a_zero_area_ring_is_rejected() {
    let flat = ring(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]);
    assert!(ring_self_intersection_rat(&flat).is_some());
}

// --- the rational sign ---

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// The sign primitive agrees with the f64 predicate wherever both are exact — integer
    /// coordinates are exact in both worlds, so a disagreement is a bug in one of them.
    #[test]
    fn the_rational_orientation_agrees_with_shewchuk_on_exact_input(
        xs in prop::array::uniform6(-64i32..64),
    ) {
        let f = |i: usize| [xs[2 * i] as f64, xs[2 * i + 1] as f64];
        let r = |i: usize| [Rat::from_int(xs[2 * i] as i128), Rat::from_int(xs[2 * i + 1] as i128)];
        let sign = orient2d(f(0), f(1), f(2));
        let sign = (sign > 0.0) as i8 - (sign < 0.0) as i8;
        prop_assert_eq!(sign, orient2d_rat(r(0), r(1), r(2)));
    }
}

/// **The ring predicates judge the rational truth because the two worlds disagree — here is the
/// disagreement.**
///
/// `(0, 0.1) → (0.1, 0.2) → (0.2, 0.3)` is collinear in the decimals the author wrote
/// (slope one), but not in the binary values the f64s hold: `0.2` is exactly `2·0.1bin`,
/// while `0.3bin ≠ 3·0.1bin`, so Shewchuk's exact sign of the *binary* points is nonzero.
/// The truth is what was written, which is why `Profile2d`
/// judges the rational side of this fork.
#[test]
fn a_decimal_collinearity_the_binary_points_do_not_have() {
    let (a, b, c) = ([0.0, 0.1], [0.1, 0.2], [0.2, 0.3]);
    assert_ne!(orient2d(a, b, c), 0.0, "binary: a hair off the line");
    let lift = |p: [f64; 2]| [d(p[0]), d(p[1])];
    assert_eq!(orient2d_rat(lift(a), lift(b), lift(c)), 0, "decimal: on it");
}

/// A flat corner is strictly mid-run; everything `check` must still see is not one.
#[test]
fn a_flat_corner_is_strictly_mid_run() {
    let p = |x: f64, y: f64| [d(x), d(y)];
    assert!(
        flat_corner(p(0.0, 0.0), p(2.0, 0.0), p(4.0, 0.0)),
        "mid-run"
    );
    assert!(
        !flat_corner(p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0)),
        "a real corner"
    );
    // A repeated point is a zero-length edge to report, not a corner to erase.
    assert!(
        !flat_corner(p(0.0, 0.0), p(4.0, 0.0), p(4.0, 0.0)),
        "repeated point"
    );
    // A spike's tip is collinear but not between its neighbours — a self-intersection to report.
    assert!(
        !flat_corner(p(0.0, 0.0), p(4.0, 0.0), p(2.0, 0.0)),
        "a spike's tip"
    );
    // Decided on the decimal truth, where `0.1·0.2·0.3` is a line.
    assert!(
        flat_corner(p(0.0, 0.1), p(0.1, 0.2), p(0.2, 0.3)),
        "decimal collinear"
    );
}
