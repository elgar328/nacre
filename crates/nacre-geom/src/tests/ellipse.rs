use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn unit_circle() -> Ellipse {
    Ellipse::from_center_normal(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        1.0,
    )
    .unwrap()
}

// --- golden ---

/// `angle_of` is `point_at`'s inverse, on all four quadrants and on the seam — and it
/// answers for a point off the circle by its in-plane direction (a scaled point maps to the
/// same angle), which is how the arc consumers ask it about realized endpoints.
#[test]
fn angle_of_inverts_point_at() {
    let c = Ellipse::from_center_normal(
        Point3::from_array([1.0, -2.0, 0.5]),
        Vector3::from_array([0.0, 3.0, 4.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.5,
    )
    .unwrap();
    for k in 0..8 {
        let theta = std::f64::consts::TAU * (k as f64) / 8.0;
        let back = c.angle_of(c.point_at(theta));
        let diff = (back - theta).abs();
        assert!(
            diff < 1e-12 || (std::f64::consts::TAU - diff) < 1e-12,
            "theta {theta} came back as {back}"
        );
    }
    // Off-circle: the direction speaks, not the distance.
    let p = c.center() + (c.point_at(1.0) - c.center()) * 0.25;
    assert!((c.angle_of(p) - 1.0).abs() < 1e-12);
    // The seam itself is angle zero, not 2π.
    assert!(c.angle_of(c.point_at(0.0)) < 1e-15);
}

#[test]
fn unit_circle_basics() {
    let c = unit_circle();
    // sin_cos(0) = (0, 1) exactly, so these are exact.
    assert_eq!(c.point_at(0.0).as_array(), [1.0, 0.0, 0.0]);
    assert_eq!(c.derivative(0.0).as_array(), [0.0, 1.0, 0.0]);
    // Distances are exact (same bits as `f64::sqrt`).
    assert_eq!(c.distance(Point3::origin()), 1.0); // the axis point, = radius
    assert_eq!(c.distance(Point3::from_array([1.0, 0.0, 0.0])), 0.0); // on the circle
    assert_eq!(c.distance(Point3::from_array([2.0, 0.0, 0.0])), 1.0); // radially outside
    assert_eq!(
        c.distance(Point3::from_array([0.0, 0.0, 5.0])),
        26.0_f64.sqrt()
    ); // up the axis
}

#[test]
fn ref_dir_is_orthogonalized_and_normalized() {
    // A ref_dir with a normal component and non-unit length is projected and
    // normalized: normal +Z, ref_dir (2, 0, 7) → angle-0 axis +X.
    let c = Ellipse::from_center_normal(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([2.0, 0.0, 7.0]),
        3.0,
    )
    .unwrap();
    assert_eq!(c.ref_dir().as_array(), [1.0, 0.0, 0.0]);
    assert_eq!(c.point_at(0.0).as_array(), [3.0, 0.0, 0.0]);
}

#[test]
fn degenerate_constructions_return_none() {
    let c = Point3::origin();
    let z = Vector3::from_array([0.0, 0.0, 1.0]);
    let x = Vector3::from_array([1.0, 0.0, 0.0]);
    assert!(Ellipse::from_center_normal(c, Vector3::zero(), x, 1.0).is_none()); // zero normal
    assert!(Ellipse::from_center_normal(c, z, z, 1.0).is_none()); // ref_dir ∥ normal
    assert!(Ellipse::from_center_normal(c, z, x, 0.0).is_none()); // zero radius
    assert!(Ellipse::from_center_normal(c, z, x, -2.0).is_none()); // negative radius
}

// --- proptest ---

fn pt3() -> impl Strategy<Value = Point3> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Point3::from_array)
}
fn vec3() -> impl Strategy<Value = Vector3> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Vector3::from_array)
}

/// A well-conditioned circle: normal and ref_dir are kept clearly non-parallel
/// so the in-plane frame is stable (mirrors `plane_and_points`).
fn circle() -> impl Strategy<Value = Ellipse> {
    (pt3(), vec3(), vec3(), 0.5f64..10.0).prop_filter_map(
        "degenerate circle frame",
        |(center, normal, ref_dir, radius)| {
            let separated = matches!(
                (normal.normalize(), ref_dir.normalize()),
                (Some(n), Some(r)) if n.cross(r).norm() >= 1e-3
            );
            if separated {
                Ellipse::from_center_normal(center, normal, ref_dir, radius)
            } else {
                None
            }
        },
    )
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn derivative_matches_central_difference(c in circle(), theta in -10.0f64..10.0) {
        // The design-mandated numerical-differentiation check. The derivative is
        // center-independent (impl uses only ref_dir/binormal/radius), so evaluate
        // it on an origin-centered copy — differencing `point_at` at a large center
        // would just add catastrophic cancellation from the huge shared term.
        let c0 =
            Ellipse::from_center_normal(Point3::origin(), c.axis(), c.ref_dir(), c.radius())
                .unwrap();
        let h = 1e-6;
        let central = (c0.point_at(theta + h) - c0.point_at(theta - h)) / (2.0 * h);
        let tol = 1e-6 * (c0.radius() + 1.0);
        prop_assert!((c0.derivative(theta) - central).norm() <= tol);
    }

    #[test]
    fn point_on_circle_has_zero_distance(c in circle(), theta in -10.0f64..10.0) {
        let p = c.point_at(theta);
        let scale = 1e-6 * (p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!(c.distance(p) <= scale);
    }

    #[test]
    fn derivative_norm_is_radius(c in circle(), theta in -10.0f64..10.0) {
        prop_assert!(approx_eq(c.derivative(theta).norm(), c.radius(), EPS, EPS));
    }

    #[test]
    fn derivative_is_perpendicular_to_radius(c in circle(), theta in -10.0f64..10.0) {
        // Velocity ⟂ radius vector — an analytic cross-check independent of the
        // finite difference above.
        let radial = c.point_at(theta) - c.center();
        let scale = 1e-6 * (c.radius() * c.radius() + 1.0);
        prop_assert!(radial.dot(c.derivative(theta)).abs() <= scale);
    }
}

// --- a tilted curve: the section of a cylinder by a plane oblique to its axis ---

/// A cylinder's rim (centre near the origin — the maths is translation-free, and `translated`
/// is locked on its own) tilted onto a plane whose normal stays clear of the cylinder's axis
/// direction's perpendicular (`|n·axis| ≥ 0.2`, so the semi-major is at most five radii).
fn tilted() -> impl Strategy<Value = Ellipse> {
    (
        prop::array::uniform3(-3.0f64..3.0),
        prop::array::uniform3(-1.0f64..1.0),
        prop::array::uniform3(-1.0f64..1.0),
        prop::array::uniform3(-1.0f64..1.0),
        0.5f64..5.0,
    )
        .prop_filter_map(
            "a degenerate frame or a plane too close to the axis",
            |(c, a, r, n, radius)| {
                let (a, r, n) = (
                    Vector3::from_array(a).normalize()?,
                    Vector3::from_array(r),
                    Vector3::from_array(n).normalize()?,
                );
                if a.cross(r).norm() < 1e-2 || n.dot(a).abs() < 0.2 || n.cross(a).norm() < 1e-3 {
                    return None;
                }
                Ellipse::from_center_normal(Point3::from_array(c), a, r, radius)?.tilted_into(n)
            },
        )
}

/// The curve's own plane, from three of its points — an oracle that owes nothing to
/// `plane_normal`.
fn plane_of(e: Ellipse) -> Vector3 {
    let (p0, p1, p2) = (e.point_at(0.3), e.point_at(2.1), e.point_at(4.4));
    (p1 - p0).cross(p2 - p0).normalize().unwrap()
}

/// The arc from `t0` over `dt`, as a fine polygon closed by its chord: `(area, centroid)` in its
/// plane — the numbers `segment_area` and `segment_centroid` claim, measured.
fn segment_by_polygon(e: Ellipse, t0: f64, dt: f64) -> (f64, Point3) {
    let n = 4000;
    let pts: Vec<Point3> = (0..=n)
        .map(|i| e.point_at(t0 + dt * i as f64 / n as f64))
        .collect();
    let normal = plane_of(e);
    let base = pts[0];
    let (mut area, mut moment) = (0.0, Vector3::zero());
    for w in pts[1..].windows(2) {
        let (u, v) = (w[0] - base, w[1] - base);
        let a = 0.5 * u.cross(v).dot(normal);
        area += a;
        moment += (u + v) * (a / 3.0);
    }
    (area, base + moment * (1.0 / area))
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]

    /// Every point lies on its cylinder (`radius` off the axis) and on its plane, and `angle_of`
    /// reads back the angle it was made at — the tilt is along the axis, so θ is the cylinder's.
    #[test]
    fn a_tilted_point_lies_on_its_cylinder_and_plane(e in tilted(), theta in 0.0f64..std::f64::consts::TAU) {
        let p = e.point_at(theta);
        let w = p - e.center();
        let off_axis = (w - e.axis() * w.dot(e.axis())).norm();
        prop_assert!((off_axis - e.radius()).abs() <= 1e-9 * (e.radius() + 1.0));
        prop_assert!(w.dot(e.plane_normal()).abs() <= 1e-9 * (e.semi_major() + 1.0));
        prop_assert!(plane_of(e).cross(e.plane_normal()).norm() <= 1e-9);
        prop_assert!(e.plane_normal().dot(e.axis()) > 0.0, "the plane normal on the axis's side");
        prop_assert!((e.angle_of(p) - theta).abs() <= 1e-9);
        prop_assert!(e.distance(p) <= 1e-9 * (e.semi_major() + 1.0));
        prop_assert!((e.z_at(theta) - w.dot(e.axis())).abs() <= 1e-9 * (e.semi_major() + 1.0));
    }

    #[test]
    fn a_tilted_derivative_matches_central_difference(e in tilted(), theta in 0.0f64..std::f64::consts::TAU) {
        let h = 1e-6;
        let central = (e.point_at(theta + h) - e.point_at(theta - h)) / (2.0 * h);
        prop_assert!((e.derivative(theta) - central).norm() <= 1e-6 * (e.semi_major() + 1.0));
    }

    /// The segment's area and centroid, and the whole curve's area, against a fine polygon.
    #[test]
    fn a_tilted_segment_is_the_polygons(e in tilted(), t0 in 0.0f64..std::f64::consts::TAU, dt in 0.2f64..6.2) {
        let (area, centroid) = segment_by_polygon(e, t0, dt);
        let scale = e.radius() * e.semi_major();
        prop_assert!((e.segment_area(dt) - area).abs() <= 1e-5 * scale, "{} vs {}", e.segment_area(dt), area);
        prop_assert!((e.segment_centroid(t0, dt) - centroid).norm() <= 1e-4 * e.semi_major());
        let (whole, at) = segment_by_polygon(e, 0.0, std::f64::consts::TAU);
        prop_assert!((e.area() - whole).abs() <= 1e-5 * scale);
        prop_assert!((at - e.center()).norm() <= 1e-6 * e.semi_major());
    }

    /// The half-extent bounds every sampled point and is reached; the semi-axes are the farthest
    /// and nearest points from the centre; the largest radius of curvature is at least the
    /// curvature radius anywhere.
    #[test]
    fn a_tilted_curves_extent_axes_and_curvature(e in tilted()) {
        let n = 2000;
        let mut reach = [0.0f64; 3];
        let (mut far, mut near) = (0.0f64, f64::MAX);
        let mut curve_radius = 0.0f64;
        for i in 0..n {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            let d = e.point_at(t) - e.center();
            for (k, x) in d.as_array().into_iter().enumerate() {
                reach[k] = reach[k].max(x.abs());
            }
            far = far.max(d.norm());
            near = near.min(d.norm());
            let (d1, h) = (e.derivative(t), 1e-5);
            let d2 = (e.derivative(t + h) - e.derivative(t - h)) / (2.0 * h);
            curve_radius = curve_radius.max(d1.norm().powi(3) / d1.cross(d2).norm());
        }
        let half = e.half_extent();
        for k in 0..3 {
            prop_assert!(half[k] + 1e-9 >= reach[k] && half[k] - reach[k] <= 1e-4 * e.semi_major());
        }
        prop_assert!((e.semi_major() - far).abs() <= 1e-5 * far);
        prop_assert!((e.radius() - near).abs() <= 1e-5 * far);
        prop_assert!(e.max_curvature_radius() + 1e-6 >= curve_radius * (1.0 - 1e-4));
    }

    /// **What a STEP file writes is this curve, the same way round.** The placement's
    /// `center + semi_1 cos u·ref_dir + semi_2 sin u·(normal × ref_dir)` is the point at `θ* + u`,
    /// `θ*` where the major direction points — so the file's parameter runs as θ does.
    #[test]
    fn a_tilted_curve_is_written_as_itself(e in tilted(), u in 0.0f64..std::f64::consts::TAU) {
        let w = e.written();
        prop_assert!(w.normal.dot(w.ref_dir).abs() <= 1e-12);
        prop_assert!(w.semi_1 >= w.semi_2);
        let y = w.normal.cross(w.ref_dir);
        let file = w.center + w.ref_dir * (w.semi_1 * u.cos()) + y * (w.semi_2 * u.sin());
        let [t0, t1] = e.axial_trig();
        let ours = e.point_at(t1.atan2(t0) + u);
        prop_assert!((file - ours).norm() <= 1e-9 * (e.semi_major() + 1.0), "{file:?} vs {ours:?}");
    }
}

/// A plane across the axis has no tilt to give — `n·axis` is all of `n`, and the tilt comes out
/// `[0, 0]` — and a plane along it meets no curve; both refuse rather than make a curve.
#[test]
fn a_plane_across_or_along_the_axis_tilts_nothing() {
    let c = unit_circle();
    assert!(
        c.tilted_into(Vector3::from_array([0.0, 0.0, 1.0]))
            .is_none()
    );
    assert!(
        c.tilted_into(Vector3::from_array([1.0, 0.0, 0.0]))
            .is_none()
    );
    let e = c
        .tilted_into(Vector3::from_array([1.0, 0.0, 1.0]).normalize().unwrap())
        .expect("a 45° plane tilts the circle");
    assert!(!e.is_circle() && c.is_circle());
    // `z = −x` on the plane `x + z = 0`: the tilt is `[−1, 0]`, the semi-major `√2`.
    assert_eq!(e.axial_trig(), [-1.0, 0.0]);
    assert!((e.semi_major() - 2f64.sqrt()).abs() < 1e-15);
}

/// **A circle's numbers are the circle formulas', bit for bit** — every method keeps the circle's
/// arithmetic where there is no tilt, which is what keeps every circle a consumer reads unmoved.
#[test]
fn a_circle_is_the_circle_formulas_to_the_bit() {
    let c = Ellipse::from_center_normal(
        Point3::from_array([1.5, -2.25, 0.75]),
        Vector3::from_array([0.0, 3.0, 4.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.5,
    )
    .unwrap();
    let (n, x, r, o) = (c.axis(), c.ref_dir(), c.radius(), c.center());
    let y = n.cross(x);
    for theta in [0.0, 0.7, 2.0, 3.9, 5.5] {
        let (sin, cos) = f64::sin_cos(theta);
        assert_eq!(c.point_at(theta), o + x * (r * cos) + y * (r * sin));
        assert_eq!(c.derivative(theta), x * (-r * sin) + y * (r * cos));
        assert_eq!(c.segment_area(theta), 0.5 * r * r * (theta - theta.sin()));
        assert_eq!(c.z_at(theta), 0.0);
    }
    assert_eq!(c.plane_normal(), n);
    assert_eq!(c.area(), std::f64::consts::PI * r * r);
    assert_eq!(c.max_curvature_radius(), r);
    let na = n.as_array();
    assert_eq!(
        c.half_extent(),
        std::array::from_fn(|i| r * (1.0 - na[i] * na[i]).max(0.0).sqrt())
    );
    let w = c.written();
    assert_eq!(
        (w.center, w.normal, w.ref_dir, w.semi_1, w.semi_2),
        (o, n, x, r, r)
    );
}
