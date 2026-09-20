use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::from_array([x, y, z])
}

/// A clamped uniform knot vector for `num_cp` control points of `degree`.
fn clamped_uniform_knots(degree: usize, num_cp: usize) -> Vec<f64> {
    let interior = num_cp - degree - 1;
    let mut k = vec![0.0; degree + 1];
    for i in 1..=interior {
        k.push(i as f64 / (interior + 1) as f64);
    }
    k.extend(std::iter::repeat_n(1.0, degree + 1));
    k
}

// --- golden ---

#[test]
fn clamped_endpoints_interpolate() {
    let c = NurbsCurve::new(
        2,
        vec![p3(0.0, 0.0, 0.0), p3(1.0, 2.0, 0.0), p3(3.0, 0.0, 0.0)],
        vec![1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    )
    .unwrap();
    assert_eq!(c.point_at(0.0).as_array(), [0.0, 0.0, 0.0]);
    assert_eq!(c.point_at(1.0).as_array(), [3.0, 0.0, 0.0]);
}

#[test]
fn degree_one_is_a_polyline() {
    let c = NurbsCurve::new(
        1,
        vec![p3(0.0, 0.0, 0.0), p3(2.0, 0.0, 0.0)],
        vec![1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .unwrap();
    assert_eq!(c.point_at(0.5).as_array(), [1.0, 0.0, 0.0]);
}

#[test]
fn rational_quadratic_is_a_circular_arc() {
    // Standard exact quarter-circle: corner weight cos(45°) = √2/2.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let c = NurbsCurve::new(
        2,
        vec![p3(1.0, 0.0, 0.0), p3(1.0, 1.0, 0.0), p3(0.0, 1.0, 0.0)],
        vec![1.0, w, 1.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    )
    .unwrap();
    let mid = c.point_at(0.5);
    let r = (mid - Point3::origin()).norm();
    assert!(approx_eq(r, 1.0, EPS, EPS), "off the unit circle: r={r}");
}

#[test]
fn degenerate_constructions_return_none() {
    let cp = vec![p3(0.0, 0.0, 0.0), p3(1.0, 0.0, 0.0), p3(2.0, 0.0, 0.0)];
    let ok = (
        2,
        cp.clone(),
        vec![1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    );
    assert!(NurbsCurve::new(ok.0, ok.1.clone(), ok.2.clone(), ok.3.clone()).is_some());
    // degree 0
    assert!(
        NurbsCurve::new(0, cp.clone(), vec![1.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]).is_none()
    );
    // too few control points for degree 3
    assert!(NurbsCurve::new(3, cp.clone(), vec![1.0, 1.0, 1.0], ok.3.clone()).is_none());
    // weight count mismatch
    assert!(NurbsCurve::new(2, cp.clone(), vec![1.0, 1.0], ok.3.clone()).is_none());
    // non-positive weight
    assert!(NurbsCurve::new(2, cp.clone(), vec![1.0, 0.0, 1.0], ok.3.clone()).is_none());
    // wrong knot count
    assert!(
        NurbsCurve::new(2, cp.clone(), vec![1.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]).is_none()
    );
    // decreasing knots
    assert!(
        NurbsCurve::new(
            2,
            cp,
            vec![1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 0.5, 1.0]
        )
        .is_none()
    );
}

#[test]
fn degree_one_derivative_is_constant() {
    // A line C(u) = (1−u)P0 + u·P1 over [0,1] has C'(u) = P1 − P0.
    let c = NurbsCurve::new(
        1,
        vec![p3(1.0, 2.0, 3.0), p3(4.0, 6.0, 8.0)],
        vec![1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .unwrap();
    assert_eq!(c.derivative(0.5).as_array(), [3.0, 4.0, 5.0]);
}

// --- surface golden ---

/// A bilinear (degree 1×1) surface over the unit square.
fn bilinear(p00: Point3, p10: Point3, p01: Point3, p11: Point3) -> NurbsSurface {
    NurbsSurface::new(
        1,
        1,
        vec![vec![p00, p01], vec![p10, p11]], // [i][j]: i→u, j→v
        vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .unwrap()
}

#[test]
fn surface_clamped_corners_interpolate() {
    let s = bilinear(
        p3(0.0, 0.0, 0.0),
        p3(2.0, 0.0, 1.0),
        p3(0.0, 3.0, 2.0),
        p3(2.0, 3.0, 5.0),
    );
    assert_eq!(s.point_at(0.0, 0.0).as_array(), [0.0, 0.0, 0.0]); // P[0][0]
    assert_eq!(s.point_at(1.0, 1.0).as_array(), [2.0, 3.0, 5.0]); // P[1][1]
}

#[test]
fn bilinear_midpoint_is_corner_average() {
    let s = bilinear(
        p3(0.0, 0.0, 0.0),
        p3(2.0, 0.0, 0.0),
        p3(0.0, 2.0, 0.0),
        p3(2.0, 2.0, 4.0),
    );
    assert_eq!(s.point_at(0.5, 0.5).as_array(), [1.0, 1.0, 1.0]); // average of the four
}

#[test]
fn rational_surface_is_a_quarter_cylinder() {
    // u: exact quarter-circle (weight √2/2 on the middle column); v: linear
    // extrusion 0→h. Asymmetric 3×2 grid catches u/v and i/j swaps.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let h = 4.0;
    let s = NurbsSurface::new(
        2,
        1,
        vec![
            vec![p3(1.0, 0.0, 0.0), p3(1.0, 0.0, h)],
            vec![p3(1.0, 1.0, 0.0), p3(1.0, 1.0, h)],
            vec![p3(0.0, 1.0, 0.0), p3(0.0, 1.0, h)],
        ],
        vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .unwrap();
    for &v in &[0.0, 0.3, 0.7, 1.0] {
        let p = s.point_at(0.5, v).as_array();
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        assert!(
            approx_eq(r, 1.0, EPS, EPS),
            "off the unit cylinder at v={v}: r={r}"
        );
        assert!(approx_eq(p[2], v * h, EPS, EPS), "wrong height at v={v}");
    }
}

#[test]
fn surface_degenerate_constructions_return_none() {
    let ok_cp = vec![
        vec![p3(0.0, 0.0, 0.0), p3(0.0, 1.0, 0.0)],
        vec![p3(1.0, 0.0, 0.0), p3(1.0, 1.0, 0.0)],
    ];
    let ok_w = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
    let ku = vec![0.0, 0.0, 1.0, 1.0];
    let mk = |du, dv, cp, w, ku: Vec<f64>, kv: Vec<f64>| NurbsSurface::new(du, dv, cp, w, ku, kv);
    assert!(mk(1, 1, ok_cp.clone(), ok_w.clone(), ku.clone(), ku.clone()).is_some());
    assert!(mk(0, 1, ok_cp.clone(), ok_w.clone(), ku.clone(), ku.clone()).is_none()); // degree_u 0
    // non-rectangular grid
    let ragged = vec![
        vec![p3(0.0, 0.0, 0.0), p3(0.0, 1.0, 0.0)],
        vec![p3(1.0, 0.0, 0.0)],
    ];
    assert!(mk(1, 1, ragged, ok_w.clone(), ku.clone(), ku.clone()).is_none());
    // weight grid mismatch
    assert!(
        mk(
            1,
            1,
            ok_cp.clone(),
            vec![vec![1.0, 1.0]],
            ku.clone(),
            ku.clone()
        )
        .is_none()
    );
    // non-positive weight
    assert!(
        mk(
            1,
            1,
            ok_cp.clone(),
            vec![vec![1.0, 0.0], vec![1.0, 1.0]],
            ku.clone(),
            ku.clone()
        )
        .is_none()
    );
    // wrong knot count
    assert!(
        mk(
            1,
            1,
            ok_cp.clone(),
            ok_w.clone(),
            vec![0.0, 0.0, 1.0],
            ku.clone()
        )
        .is_none()
    );
    // decreasing knots
    assert!(mk(1, 1, ok_cp, ok_w, vec![0.0, 1.0, 0.0, 1.0], ku).is_none());
}

#[test]
fn bilinear_partials_are_known() {
    // S(u,v)=(1−u)(1−v)P00+u(1−v)P10+(1−u)v P01+uv P11 over unit square.
    // ∂u(.5,.5)=½(P10−P00)+½(P11−P01); ∂v(.5,.5)=½(P01−P00)+½(P11−P10).
    let s = bilinear(
        p3(0.0, 0.0, 0.0),
        p3(1.0, 0.0, 0.0),
        p3(0.0, 1.0, 0.0),
        p3(1.0, 1.0, 2.0),
    );
    assert_eq!(s.du(0.5, 0.5).as_array(), [1.0, 0.0, 1.0]);
    assert_eq!(s.dv(0.5, 0.5).as_array(), [0.0, 1.0, 1.0]);
}

// --- proptest ---

fn nurbs() -> impl Strategy<Value = NurbsCurve> {
    (2usize..=3, 4usize..=7).prop_flat_map(|(degree, m)| {
        let pts = prop::collection::vec(prop::array::uniform3(-10.0f64..10.0), m);
        let ws = prop::collection::vec(0.5f64..2.0, m);
        (Just(degree), pts, ws).prop_map(|(degree, pts, ws)| {
            let cp: Vec<Point3> = pts.into_iter().map(Point3::from_array).collect();
            let knots = clamped_uniform_knots(degree, cp.len());
            NurbsCurve::new(degree, cp, ws, knots).unwrap()
        })
    })
}

fn nurbs_surface() -> impl Strategy<Value = NurbsSurface> {
    (2usize..=3, 2usize..=3, 3usize..=5, 3usize..=5).prop_flat_map(|(du, dv, nu, nv)| {
        let nu = nu.max(du + 1);
        let nv = nv.max(dv + 1);
        let pts = prop::collection::vec(
            prop::collection::vec(prop::array::uniform3(-10.0f64..10.0), nv),
            nu,
        );
        let ws = prop::collection::vec(prop::collection::vec(0.5f64..2.0, nv), nu);
        (Just((du, dv)), pts, ws).prop_map(move |((du, dv), pts, ws)| {
            let cp: Vec<Vec<Point3>> = pts
                .into_iter()
                .map(|row| row.into_iter().map(Point3::from_array).collect())
                .collect();
            let ku = clamped_uniform_knots(du, nu);
            let kv = clamped_uniform_knots(dv, nv);
            NurbsSurface::new(du, dv, cp, ws, ku, kv).unwrap()
        })
    })
}

/// **Is `t` far enough from every interior knot for a central difference to be second-order?**
///
/// ★★★ **The oracle, not the surface, is what a knot breaks.** A central difference is
/// second-order only where the function is C²; at a simple interior knot a degree-`p` spline is
/// only C^(p−1), so a sample within `h` of one differences *across* the break and the error
/// degrades to O(h) — larger than these tests' tolerance, while `derivative`/`du`/`dv`
/// themselves are perfectly correct there (the curve is still C¹).
///
/// Measured, and the reason this is written down rather than absorbed into a looser tolerance:
/// proptest found `v = 0.66666679` against a knot at `2/3` — `1.2e-7` away with `h = 1e-6`. The
/// window is `2h` wide out of a unit domain, so roughly one sample in `10⁵` lands in it; that is
/// the shape a fixed corpus never finds and a proptest eventually does. Loosening the tolerance
/// instead would blind the check everywhere to buy nothing here.
///
/// ★ Both derivative oracles read this. The curve one had the same latent flaw and had simply
/// not been hit — one rule, one place.
fn clear_of_knots(t: f64, knots: &[f64], degree: usize, h: f64) -> bool {
    knots[degree + 1..knots.len() - degree - 1]
        .iter()
        .all(|k| (t - k).abs() > 2.0 * h)
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn derivative_matches_central_difference(c in nurbs(), s in 0.05f64..0.95) {
        let (lo, hi) = c.domain();
        let u = lo + s * (hi - lo);
        let h = 1e-6;
        prop_assume!(clear_of_knots(u, c.knots(), c.degree(), h));
        let central = (c.point_at(u + h) - c.point_at(u - h)) / (2.0 * h);
        let extent: f64 = c
            .control_points()
            .iter()
            .map(|p| (*p - Point3::origin()).norm())
            .fold(0.0, f64::max);
        let tol = 1e-6 * (extent + 1.0);
        prop_assert!((c.derivative(u) - central).norm() <= tol);
    }

    #[test]
    fn point_lies_in_control_aabb(c in nurbs(), s in 0.0f64..1.0) {
        let (lo, hi) = c.domain();
        let p = c.point_at(lo + s * (hi - lo)).as_array();
        for (axis, &pa) in p.iter().enumerate() {
            let mn = c.control_points().iter().map(|q| q.as_array()[axis]).fold(f64::MAX, f64::min);
            let mx = c.control_points().iter().map(|q| q.as_array()[axis]).fold(f64::MIN, f64::max);
            prop_assert!(pa >= mn - 1e-9 && pa <= mx + 1e-9);
        }
    }

    #[test]
    fn clamps_outside_domain(c in nurbs()) {
        let (lo, hi) = c.domain();
        prop_assert_eq!(c.point_at(lo - 5.0).as_array(), c.point_at(lo).as_array());
        prop_assert_eq!(c.point_at(hi + 5.0).as_array(), c.point_at(hi).as_array());
    }

    #[test]
    fn surface_partials_match_central_differences(
        s in nurbs_surface(),
        su in 0.05f64..0.95,
        sv in 0.05f64..0.95,
    ) {
        let ((u0, u1), (v0, v1)) = s.domain();
        let (u, v) = (u0 + su * (u1 - u0), v0 + sv * (v1 - v0));
        let h = 1e-6;
        prop_assume!(clear_of_knots(u, s.knots_u(), s.degree_u(), h));
        prop_assume!(clear_of_knots(v, s.knots_v(), s.degree_v(), h));
        let cu = (s.point_at(u + h, v) - s.point_at(u - h, v)) / (2.0 * h);
        let cv = (s.point_at(u, v + h) - s.point_at(u, v - h)) / (2.0 * h);
        let extent: f64 = s.control_points().iter().flatten()
            .map(|p| (*p - Point3::origin()).norm()).fold(0.0, f64::max);
        let tol = 1e-6 * (extent + 1.0);
        prop_assert!((s.du(u, v) - cu).norm() <= tol);
        prop_assert!((s.dv(u, v) - cv).norm() <= tol);
    }

    #[test]
    fn surface_normal_is_unit_and_perpendicular(
        s in nurbs_surface(),
        su in 0.1f64..0.9,
        sv in 0.1f64..0.9,
    ) {
        let ((u0, u1), (v0, v1)) = s.domain();
        let (u, v) = (u0 + su * (u1 - u0), v0 + sv * (v1 - v0));
        let (du, dv) = (s.du(u, v), s.dv(u, v));
        // Skip near-degenerate parameterizations (parallel/zero partials).
        prop_assume!(du.cross(dv).norm() > 1e-6 * (du.norm() * dv.norm() + 1.0));
        let n = s.normal_at(u, v).unwrap();
        prop_assert!(approx_eq(n.norm(), 1.0, EPS, EPS));
        let scale = du.norm().max(dv.norm()) + 1.0;
        prop_assert!(n.dot(du).abs() <= 1e-6 * scale);
        prop_assert!(n.dot(dv).abs() <= 1e-6 * scale);
    }

    #[test]
    fn surface_point_lies_in_control_aabb(s in nurbs_surface(), su in 0.0f64..1.0, sv in 0.0f64..1.0) {
        let ((u0, u1), (v0, v1)) = s.domain();
        let p = s.point_at(u0 + su * (u1 - u0), v0 + sv * (v1 - v0)).as_array();
        for (axis, &pa) in p.iter().enumerate() {
            let coords = || s.control_points().iter().flatten().map(|q| q.as_array()[axis]);
            let mn = coords().fold(f64::MAX, f64::min);
            let mx = coords().fold(f64::MIN, f64::max);
            prop_assert!(pa >= mn - 1e-9 && pa <= mx + 1e-9);
        }
    }
}
