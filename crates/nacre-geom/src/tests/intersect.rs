use super::*;
use nacre_math::Vector3;
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

// --- golden ---

#[test]
fn plane_plane_axis_planes_give_z_axis() {
    // x = 0 (normal +x) ∩ y = 0 (normal +y) = the z-axis.
    let line = plane_plane(
        &plane([0.0; 3], [1.0, 0.0, 0.0]),
        &plane([0.0; 3], [0.0, 1.0, 0.0]),
    )
    .unwrap();
    let dir = line.direction().as_array();
    assert!(dir[0].abs() < 1e-15 && dir[1].abs() < 1e-15 && dir[2].abs() > 1.0 - 1e-15);
    let o = line.origin().as_array();
    assert!(o[0].abs() < 1e-15 && o[1].abs() < 1e-15);
}

#[test]
fn plane_plane_parallel_is_none() {
    // Same normal, different offset — parallel, never meet.
    assert!(
        plane_plane(
            &plane([0.0; 3], [0.0, 0.0, 1.0]),
            &plane([0.0, 0.0, 3.0], [0.0, 0.0, 1.0])
        )
        .is_none()
    );
    // Anti-parallel normals are parallel too.
    assert!(
        plane_plane(
            &plane([0.0; 3], [0.0, 0.0, 1.0]),
            &plane([0.0, 0.0, 3.0], [0.0, 0.0, -1.0])
        )
        .is_none()
    );
}

// --- proptest ---

fn coord() -> impl Strategy<Value = f64> {
    -100.0f64..100.0
}

fn vec3() -> impl Strategy<Value = Vector3> {
    prop::array::uniform3(-1.0f64..1.0).prop_map(Vector3::from_array)
}

/// Two planes whose normals are well-separated (so their line is
/// well-conditioned).
fn two_planes() -> impl Strategy<Value = (Plane, Plane)> {
    (
        prop::array::uniform3(coord()),
        vec3(),
        prop::array::uniform3(coord()),
        vec3(),
    )
        .prop_filter_map("zero/near-parallel normals", |(o1, n1, o2, n2)| {
            let a = Plane::from_point_normal(Point3::from_array(o1), n1)?;
            let b = Plane::from_point_normal(Point3::from_array(o2), n2)?;
            (a.normal().cross(b.normal()).norm() >= 0.1).then_some((a, b))
        })
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// The returned line lies on both planes and runs perpendicular to both
    /// normals.
    #[test]
    fn plane_plane_line_lies_on_both((a, b) in two_planes(), t in -100.0f64..100.0) {
        let line = plane_plane(&a, &b).unwrap();
        let dir = line.direction();
        prop_assert!(a.normal().dot(dir).abs() <= 1e-9);
        prop_assert!(b.normal().dot(dir).abs() <= 1e-9);
        let p = line.point_at(t);
        let mag = p.as_array().iter().map(|v| v.abs()).fold(0.0, f64::max);
        let scale = 1e-9 * (mag + 1.0);
        prop_assert!(a.distance(p) <= scale);
        prop_assert!(b.distance(p) <= scale);
    }
}

// --- three_planes / three_plane_orient3d ---

/// The sign of an f64 (`+1`/`-1`/`0`) — not `f64::signum`, which maps `0.0`
/// to `+1.0`; a coplanar `orient3d` (exactly `0.0`) must read as `0`.
fn sign_f64(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

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

fn ivec3() -> impl Strategy<Value = [i64; 3]> {
    prop::array::uniform3(-30i64..=30)
}

/// Pinned against `plane_plane` itself: the sign must agree with the dot of
/// that function's actual line direction and `c`'s normal. The seam ordering
/// in `nacre-ops` assumes exactly this coupling.
#[test]
fn three_plane_cmp_coord_orders_two_meets() {
    // The unit cube's corners `(0,0,0)` and `(1,1,0)`, each as a triple of its faces.
    let px = |d: f64| {
        Plane::from_point_normal(
            Point3::from_array([d, 0.0, 0.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
        )
        .unwrap()
    };
    let py = |d: f64| {
        Plane::from_point_normal(
            Point3::from_array([0.0, d, 0.0]),
            Vector3::from_array([0.0, 1.0, 0.0]),
        )
        .unwrap()
    };
    let pz =
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap();
    let (x0, x1, y0, y1) = (px(0.0), px(1.0), py(0.0), py(1.0));
    let a = [&x0, &y0, &pz];
    let b = [&x1, &y1, &pz];
    assert_eq!(three_plane_cmp_coord(a, b, 0), -1);
    assert_eq!(three_plane_cmp_coord(b, a, 0), 1);
    assert_eq!(three_plane_cmp_coord(a, b, 1), -1);
    assert_eq!(three_plane_cmp_coord(a, b, 2), 0); // both on z = 0

    // A `Plane`'s normal is unit, so the coefficients carry a `d` that is not an
    // integer here — the predicate is scale-invariant and never divides.
    let tilt = Plane::through_points(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let c = [&x0, &y0, &tilt];
    assert_eq!(three_plane_cmp_coord(c, c, 2), 0);
    assert_eq!(three_plane_cmp_coord(a, c, 2), -1); // (0,0,0) below (0,0,1)
}

#[test]
fn plane_pair_dir_sign_agrees_with_plane_plane() {
    let cases = [
        // x=0 ∩ y=0 is the z axis, direction (1,0,0)×(0,1,0) = +z.
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 1),
        // Flip c's normal ⇒ flip the sign.
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0], -1),
        // Swap a and b ⇒ the line reverses ⇒ flip the sign.
        ([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], -1),
        // c parallel to the line ⇒ normals coplanar ⇒ 0.
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], 0),
    ];
    for (na, nb, nc, expect) in cases {
        let (a, b, c) = (
            plane([0.0; 3], na),
            plane([0.0; 3], nb),
            plane([0.0; 3], nc),
        );
        assert_eq!(
            plane_pair_dir_sign(&a, &b, &c),
            expect,
            "{na:?} {nb:?} {nc:?}"
        );
        let l = plane_plane(&a, &b).expect("a and b are not parallel here");
        let dot = l.direction().dot(c.normal());
        let dot_sign = if dot > 1e-12 {
            1
        } else if dot < -1e-12 {
            -1
        } else {
            0
        };
        assert_eq!(dot_sign, expect, "plane_plane's direction disagrees");
    }
}

/// The two anti-parallel side walls of a `(1,1,1)`-slanted prism, as they actually came off
/// `build_prism`. Their `raw` rows sum to an exact zero in the 2×2 minors,
/// so `det = 0` against any third plane — but the `sqrt`-rounded **unit** rows do not, and the
/// old `normal()`-based predicate returned `±1`, admitting a triple its consumer rejects as
/// `D = 0`. Pinned with the real coefficients so a regression to `normal()` fails with the
/// reason attached. `raw` is set verbatim via `from_point_normal` (which stores its argument as
/// `raw`); the origin is irrelevant here (the predicate reads only the normal rows).
#[test]
fn parallel_planes_stay_degenerate_after_normalization() {
    let p = |raw: [f64; 3]| {
        Plane::from_point_normal(Point3::origin(), Vector3::from_array(raw)).unwrap()
    };
    let a = p([
        -2.220446049250313e-16,
        -2.828427124746191,
        2.8284271247461907,
    ]);
    let b = p([0.0, 2.8284271247461907, -2.8284271247461907]);
    // A third, independent plane (the prism's tilted cap direction).
    let c = p([
        -0.5773502691896258,
        -0.5773502691896258,
        -0.5773502691896258,
    ]);
    assert_eq!(
        plane_pair_dir_sign(&a, &b, &c),
        0,
        "anti-parallel walls have no well-conditioned common line — det is exactly 0"
    );
    // The bug this pins: the unit-normal determinant does NOT vanish.
    let unit_det = nacre_predicates::det3_sign([
        a.normal().as_array(),
        b.normal().as_array(),
        c.normal().as_array(),
    ]);
    assert_ne!(
        unit_det, 0,
        "sanity: normalized normals round the exact zero away — the reason we read coefficients"
    );
}

#[test]
fn planes_coplanar_names_the_same_plane_regardless_of_scale_or_direction() {
    // z = 0, built three ways: two positive normals of different magnitude and
    // one opposite normal. All name the same plane.
    let a = Plane::through_points(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    )
    .unwrap();
    let bigger = Plane::through_points(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 0.0, 0.0]),
        Point3::from_array([0.0, 3.0, 0.0]),
    )
    .unwrap();
    let opposite = Plane::through_points(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
    )
    .unwrap();
    assert!(planes_coplanar(&a, &bigger));
    assert!(planes_coplanar(&a, &opposite));
    // Parallel but offset (z = 1) and non-parallel (y = 0): distinct planes.
    assert!(!planes_coplanar(
        &a,
        &plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0])
    ));
    assert!(!planes_coplanar(&a, &plane([0.0; 3], [0.0, 1.0, 0.0])));
}

/// `+1` means `V` is on the triangle's right-hand-normal side — pinned here
/// because the whole seam-ordering algebra in `nacre-ops` hangs off it, and
/// the doc comment said the opposite until this test existed.
#[test]
fn three_plane_orient3d_is_positive_on_the_rh_normal_side() {
    let at = |n: [f64; 3], d: f64| {
        Plane::from_point_normal(
            Point3::from_array([n[0] * d, n[1] * d, n[2] * d]),
            Vector3::from_array(n),
        )
        .unwrap()
    };
    // x=1, y=1, z=1 meet at V = (1,1,1).
    let (px, py, pz) = (
        at([1.0, 0.0, 0.0], 1.0),
        at([0.0, 1.0, 0.0], 1.0),
        at([0.0, 0.0, 1.0], 1.0),
    );
    // Triangle in the z=0 plane, CCW seen from +z ⇒ RH normal is +z. V is above it.
    let (q, r, s) = (
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    );
    assert_eq!(three_plane_orient3d(&px, &py, &pz, q, r, s), 1);
    // Swapping two triangle points flips its RH normal, hence the sign.
    assert_eq!(three_plane_orient3d(&px, &py, &pz, r, q, s), -1);
    // A plane's *stored* normal is not the triangle's RH normal: `pz` above has
    // normal +z, but the triangle (q, s, r) spans the same plane with RH normal
    // −z. Reading the convention off `Plane::normal()` would invert the answer.
    assert_eq!(three_plane_orient3d(&px, &py, &pz, q, s, r), -1);
}

/// `plane_side` is the all-explicit twin, and shares the convention exactly: put the
/// implicit point's coordinates in and the two agree, sign for sign. Pinning that
/// here is what lets a caller mix them without thinking.
#[test]
fn plane_side_shares_three_plane_orient3d_s_convention() {
    // The same triangle in `z = 0`, RH normal `+z`.
    let (q, r, s) = (
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    );
    let v = Point3::from_array([1.0, 1.0, 1.0]); // above it
    assert_eq!(plane_side([q, r, s], v), 1);
    assert_eq!(plane_side([r, q, s], v), -1); // flip the triangle, flip the sign
    assert_eq!(
        plane_side([q, r, s], Point3::from_array([1.0, 1.0, -1.0])),
        -1
    );
    // On the plane is exactly zero, however far from the triangle itself.
    assert_eq!(
        plane_side([q, r, s], Point3::from_array([9.0, -4.0, 0.0])),
        0
    );
    assert_eq!(plane_side([q, r, s], q), 0);
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// The two orient3d handoffs are one predicate seen from two sides: a plane
    /// triple's meet, fed to `plane_side` as coordinates, gives the same sign the
    /// implicit form gives without ever building it.
    #[test]
    fn prop_plane_side_agrees_with_three_plane_orient3d(
        v in prop::array::uniform3(-20.0f64..20.0),
        t in prop::array::uniform3(prop::array::uniform3(-20.0f64..20.0)),
    ) {
        let tri = t.map(Point3::from_array);
        let e1 = tri[1] - tri[0];
        let e2 = tri[2] - tri[0];
        prop_assume!(e1.cross(e2).norm() > 1e-6);
        // Three axis planes meeting exactly at `v`.
        let at = |n: [f64; 3]| {
            Plane::from_point_normal(Point3::from_array(v), Vector3::from_array(n)).unwrap()
        };
        let (px, py, pz) = (
            at([1.0, 0.0, 0.0]),
            at([0.0, 1.0, 0.0]),
            at([0.0, 0.0, 1.0]),
        );
        prop_assert_eq!(
            plane_side(tri, Point3::from_array(v)),
            three_plane_orient3d(&px, &py, &pz, tri[0], tri[1], tri[2])
        );
    }
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

    /// End-to-end wiring: `three_plane_orient3d` (via `coefficients()`) agrees
    /// with `orient3d` evaluated at the constructed vertex. Because `Plane`
    /// normalizes, the implicit vertex is `p + O(1e-16)`, not exactly `p`; so
    /// the config is restricted to well-conditioned normals and a non-coplanar
    /// (p, q, r, s) — there the true `orient3d` is a nonzero integer, which the
    /// tiny vertex perturbation cannot flip. This exercises the new geom code
    /// (coefficient extraction + assembly), not `indirect_orient3d` itself.
    #[test]
    fn three_plane_orient3d_matches_materialized(
        p in ivec3(),
        normals in prop::array::uniform3(ivec3()),
        q in ivec3(),
        r in ivec3(),
        s in ivec3(),
    ) {
        let pf = Point3::from_array(p.map(|v| v as f64));
        let planes: Option<Vec<Plane>> = normals
            .iter()
            .map(|n| {
                Plane::from_point_normal(pf, Vector3::from_array(n.map(|v| v as f64)))
            })
            .collect();
        prop_assume!(planes.is_some()); // reject a zero integer normal
        let planes = planes.unwrap();
        let u: Vec<_> = planes.iter().map(|pl| pl.normal()).collect();
        prop_assume!(u[0].dot(u[1].cross(u[2])).abs() >= 0.1); // well-conditioned

        let qf = Point3::from_array(q.map(|v| v as f64));
        let rf = Point3::from_array(r.map(|v| v as f64));
        let sf = Point3::from_array(s.map(|v| v as f64));
        let expected = sign_f64(nacre_predicates::orient3d(
            pf.as_array(),
            qf.as_array(),
            rf.as_array(),
            sf.as_array(),
        ));
        prop_assume!(expected != 0); // coplanar: perturbation could flip it

        prop_assert_eq!(
            three_plane_orient3d(&planes[0], &planes[1], &planes[2], qf, rf, sf),
            expected
        );
    }
}

#[test]
fn ray_and_segment_face_cross_handoff() {
    // CCW triangle in z = 0, right-hand normal +z.
    let tri = [
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    ];
    // Ray from below, straight up through the interior ⇒ forward Cross(+1).
    assert_eq!(
        ray_face_cross(
            Point3::from_array([0.25, 0.25, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            tri,
        ),
        RayCross::Cross(1)
    );
    // Ray pointing away from the triangle ⇒ Miss.
    assert_eq!(
        ray_face_cross(
            Point3::from_array([0.25, 0.25, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            tri,
        ),
        RayCross::Miss
    );
    // Segment straddling the plane through the interior ⇒ Cross(+1).
    assert_eq!(
        segment_face_cross(
            Point3::from_array([0.25, 0.25, -1.0]),
            Point3::from_array([0.25, 0.25, 1.0]),
            tri,
        ),
        SegCross::Cross(1)
    );
}

// --- ring_self_intersection ---

fn ring(pts: &[[f64; 2]]) -> Vec<Point2> {
    pts.iter().map(|p| Point2::from_array(*p)).collect()
}

/// The net that catches the likeliest way to get this wrong: edges `0` and `n-1` share a point
/// too, so a non-cyclic adjacency test reports *every* ring. Convex, non-convex, and a flat
/// (collinear) corner all have to pass.
#[test]
fn simple_polygons_have_no_self_intersection() {
    let square = ring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection(&square), None);
    let triangle = ring(&[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection(&triangle), None);
    // Reflex corner (an L).
    let l = ring(&[
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 2.0],
        [2.0, 2.0],
        [2.0, 4.0],
        [0.0, 4.0],
    ]);
    assert_eq!(ring_self_intersection(&l), None);
    // A flat corner: `(2,0)` sits mid-run on a straight edge. Collinear but not overlapping.
    let flat = ring(&[[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection(&flat), None);
}

#[test]
fn a_bowtie_crosses_itself() {
    let bowtie = ring(&[[0.0, 0.0], [4.0, 4.0], [4.0, 0.0], [0.0, 4.0]]);
    // Edges 0 (`(0,0)→(4,4)`) and 2 (`(4,0)→(0,4)`) are the crossing pair.
    assert_eq!(ring_self_intersection(&bowtie), Some((0, 2)));
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
    assert!(ring_self_intersection(&pinch).is_some());
}

/// A spike shorter than the edge it retraces, and one longer. Neither *isolates* the two-sided
/// adjacency test — with four or more points a long spike puts the ring's start point in the
/// interior of the doubling-back edge, so the non-adjacent rule catches it first. The case that
/// needs both directions is the collinear triple below.
#[test]
fn a_spike_doubling_back_is_rejected_either_length() {
    let short = ring(&[[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 4.0]]);
    assert!(ring_self_intersection(&short).is_some());
    let long = ring(&[[0.0, 0.0], [4.0, 0.0], [-2.0, 0.0], [0.0, 4.0]]);
    assert!(ring_self_intersection(&long).is_some());
}

#[test]
fn a_repeated_consecutive_point_is_a_zero_length_edge() {
    let dup = ring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection(&dup), Some((1, 1)));
    // Also across the wrap: the last point repeats the first.
    let closed = ring(&[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0], [0.0, 0.0]]);
    assert_eq!(ring_self_intersection(&closed), Some((3, 3)));
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
    assert!(ring_self_intersection(&flat).is_some());
}

// --- the rational twins ---

/// Lift literal coordinates into the truth the twins judge. `unwrap` is fine here: every
/// fixture coordinate is a short decimal, comfortably inside the window.
fn rring(pts: &[[f64; 2]]) -> Vec<[Rat; 2]> {
    pts.iter().map(|p| [d(p[0]), d(p[1])]).collect()
}

fn d(x: f64) -> Rat {
    Rat::from_decimal(x).unwrap()
}

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

    /// The twins are the same walkers. On exact (integer) coordinates every branch condition
    /// evaluates identically in f64 and in `Rat`, so all three must return the very same
    /// values — indices included — on any input, simple or not.
    #[test]
    fn the_rational_walkers_answer_what_the_f64_walkers_do(
        pts in prop::collection::vec(prop::array::uniform2(-8i32..8), 3..7),
        probe in prop::array::uniform2(-8i32..8),
    ) {
        let f: Vec<Point2> = pts.iter().map(|p| Point2::from_array([p[0] as f64, p[1] as f64])).collect();
        let r: Vec<[Rat; 2]> = pts.iter().map(|p| [Rat::from_int(p[0] as i128), Rat::from_int(p[1] as i128)]).collect();
        prop_assert_eq!(ring_self_intersection(&f), ring_self_intersection_rat(&r));
        let (fa, fb) = f.split_at(f.len() / 2);
        let (ra, rb) = r.split_at(r.len() / 2);
        prop_assert_eq!(rings_cross(fa, fb), rings_cross_rat(ra, rb));
        let fp = Point2::from_array([probe[0] as f64, probe[1] as f64]);
        let rp = [Rat::from_int(probe[0] as i128), Rat::from_int(probe[1] as i128)];
        prop_assert_eq!(point_in_ring_2d(fp, &f), point_in_ring_2d_rat(rp, &r));
    }
}

/// The fixture gallery above, re-judged on the truth — same verdicts, same indices.
#[test]
fn the_rational_self_intersection_matches_the_gallery() {
    let flat = rring(&[[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    assert_eq!(
        ring_self_intersection_rat(&flat),
        None,
        "a flat corner is legal"
    );
    let bowtie = rring(&[[0.0, 0.0], [4.0, 4.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection_rat(&bowtie), Some((0, 2)));
    let dup = rring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(ring_self_intersection_rat(&dup), Some((1, 1)));
    let spike = rring(&[[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 4.0]]);
    assert!(ring_self_intersection_rat(&spike).is_some());
}

/// ★★ **The twins exist because the two worlds disagree — here is the disagreement.**
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

/// The dissolve pass: flat corners go, everything `check` must still see survives.
#[test]
fn dissolving_flat_corners_keeps_every_reportable_defect() {
    // A flat corner dissolves, leaving the plain square.
    let flat = rring(&[[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    let square = rring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
    assert_eq!(drop_collinear_midpoints(flat), square);
    // A clean ring is untouched.
    assert_eq!(drop_collinear_midpoints(square.clone()), square);
    // A repeated point is NOT a flat corner — it must survive to be reported as a
    // zero-length edge, not silently erased as if the author had drawn it once.
    let dup = rring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
    assert_eq!(drop_collinear_midpoints(dup.clone()), dup);
    // A spike's tip is collinear but not between its neighbours — it survives for the
    // self-intersection report.
    let spike = rring(&[[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 4.0]]);
    assert_eq!(drop_collinear_midpoints(spike.clone()), spike);
    // Four points on one line: removing one midpoint makes the next one flat — the scan
    // repeats to a fixpoint.
    let run = rring(&[
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [3.0, 0.0],
        [3.0, 3.0],
        [0.0, 3.0],
    ]);
    let clean = rring(&[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]]);
    assert_eq!(drop_collinear_midpoints(run), clean);
    // A fully-collinear ring collapses below three points and is returned for `check` to
    // reject as degenerate.
    let line = rring(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]);
    assert_eq!(drop_collinear_midpoints(line).len(), 2);
}
