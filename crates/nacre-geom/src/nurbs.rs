//! NURBS curves — rational B-spline evaluation (design §3).
//!
//! Implemented from *The NURBS Book* (Piegl & Tiller): FindSpan (A2.1),
//! BasisFuns (A2.2), DersBasisFuns (A2.3), and the rational curve point /
//! derivative (A4.1 / A4.2). The basis helpers are free functions keyed on
//! `(degree, knots)` so a future `NurbsSurface` reuses them in u and v.

use nacre_math::{Point3, Vector3};

/// A rational B-spline (NURBS) curve: control points with positive weights over
/// a non-decreasing knot vector.
///
/// Invariant (constructor-enforced): `degree ≥ 1`; at least `degree + 1` control
/// points; one positive weight per control point; `knots.len() == n + degree + 2`
/// (`n` = last control-point index); knots non-decreasing. **Not `Copy`** (heap
/// control points), same rationale as the `Curve`/`Surface` enums.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsCurve {
    degree: usize,
    control_points: Vec<Point3>,
    weights: Vec<f64>,
    knots: Vec<f64>,
}

impl NurbsCurve {
    /// Build a curve, validating all invariants. `None` if `degree` is 0, there
    /// are fewer than `degree + 1` control points, the weight count differs, any
    /// weight is not positive (rejects `0`, negatives, `NaN`), the knot count is
    /// not `control_points.len() + degree + 1`, or the knots decrease.
    pub fn new(
        degree: usize,
        control_points: Vec<Point3>,
        weights: Vec<f64>,
        knots: Vec<f64>,
    ) -> Option<NurbsCurve> {
        if degree < 1 || control_points.len() < degree + 1 {
            return None;
        }
        if weights.len() != control_points.len() {
            return None;
        }
        if !weights.iter().all(|&w| w > 0.0) {
            return None;
        }
        if knots.len() != control_points.len() + degree + 1 {
            return None;
        }
        if !knots.windows(2).all(|w| w[1] >= w[0]) {
            return None;
        }
        Some(NurbsCurve {
            degree,
            control_points,
            weights,
            knots,
        })
    }

    /// The polynomial degree.
    #[inline]
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// The control points.
    #[inline]
    pub fn control_points(&self) -> &[Point3] {
        &self.control_points
    }

    /// The per-control-point weights (all positive).
    #[inline]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// The knot vector (non-decreasing).
    #[inline]
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// The parameter interval `[knots[degree], knots[n+1]]` the curve spans.
    #[inline]
    pub fn domain(&self) -> (f64, f64) {
        let n = self.control_points.len() - 1;
        (self.knots[self.degree], self.knots[n + 1])
    }

    /// The point at parameter `u` (clamped to the domain) — the rational curve
    /// point `Σ Nᵢ(u)·wᵢ·Pᵢ / Σ Nᵢ(u)·wᵢ` (A4.1).
    pub fn point_at(&self, u: f64) -> Point3 {
        let u = self.clamp(u);
        let p = self.degree;
        let span = find_span(p, &self.knots, u);
        let n = basis_funs(p, &self.knots, span, u);

        let mut num = [0.0; 3];
        let mut den = 0.0;
        for (j, &nj) in n.iter().enumerate() {
            let i = span - p + j;
            let cp = self.control_points[i].as_array();
            let nw = nj * self.weights[i];
            num[0] += nw * cp[0];
            num[1] += nw * cp[1];
            num[2] += nw * cp[2];
            den += nw;
        }
        Point3::from_array([num[0] / den, num[1] / den, num[2] / den])
    }

    /// The first parametric derivative `dC/du` at `u` (clamped to the domain),
    /// via the homogeneous derivative projected back (A4.2, `k = 1`).
    pub fn derivative(&self, u: f64) -> Vector3 {
        let u = self.clamp(u);
        let p = self.degree;
        let span = find_span(p, &self.knots, u);
        let ders = ders_basis_funs(p, &self.knots, span, u, 1);

        // Homogeneous numerator A and weight w, orders 0 and 1.
        let mut a = [[0.0; 3]; 2];
        let mut wd = [0.0; 2];
        for (j, (&d0, &d1)) in ders[0].iter().zip(&ders[1]).enumerate() {
            let i = span - p + j;
            let w = self.weights[i];
            let cp = self.control_points[i].as_array();
            let (nw0, nw1) = (d0 * w, d1 * w);
            a[0][0] += nw0 * cp[0];
            a[0][1] += nw0 * cp[1];
            a[0][2] += nw0 * cp[2];
            a[1][0] += nw1 * cp[0];
            a[1][1] += nw1 * cp[1];
            a[1][2] += nw1 * cp[2];
            wd[0] += nw0;
            wd[1] += nw1;
        }
        // C0 = A0 / w0; C1 = (A1 − w1·C0) / w0.
        let c0 = [a[0][0] / wd[0], a[0][1] / wd[0], a[0][2] / wd[0]];
        Vector3::from_array([
            (a[1][0] - wd[1] * c0[0]) / wd[0],
            (a[1][1] - wd[1] * c0[1]) / wd[0],
            (a[1][2] - wd[1] * c0[2]) / wd[0],
        ])
    }

    #[inline]
    fn clamp(&self, u: f64) -> f64 {
        let (lo, hi) = self.domain();
        u.clamp(lo, hi)
    }
}

/// The knot span index containing `u` (*The NURBS Book* A2.1). `u` must lie in
/// the domain `[knots[degree], knots[n+1]]`.
fn find_span(degree: usize, knots: &[f64], u: f64) -> usize {
    let n = knots.len() - degree - 2; // last control-point index
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[degree] {
        return degree;
    }
    let (mut lo, mut hi) = (degree, n + 1);
    let mut mid = (lo + hi) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

/// The `degree + 1` nonzero basis functions at `u` in `span` (A2.2).
fn basis_funs(degree: usize, knots: &[f64], span: usize, u: f64) -> Vec<f64> {
    let mut n = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Basis functions and their derivatives up to order `n` at `u` in `span`
/// (A2.3). Row `k` holds the `k`-th derivatives of the `degree + 1` nonzero
/// basis functions.
fn ders_basis_funs(degree: usize, knots: &[f64], span: usize, u: f64, n: usize) -> Vec<Vec<f64>> {
    let p = degree;
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r]; // lower triangle: knot differences
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp; // upper triangle: basis values
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    let mut ders = vec![vec![0.0; p + 1]; n + 1];
    for (j, d) in ders[0].iter_mut().enumerate() {
        *d = ndu[j][p]; // 0-th derivatives = basis functions
    }

    // Compute the derivatives (A2.3 inner loop) via the two-row `a` scratch.
    for r in 0..=p {
        let mut a = vec![vec![0.0; p + 1]; 2];
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=n {
            let mut dd = 0.0;
            let rk = r as isize - k as isize;
            let pk = p as isize - k as isize;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[(pk + 1) as usize][rk as usize];
                dd = a[s2][0] * ndu[rk as usize][pk as usize];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk { k - 1 } else { p - r };
            for j in j1..=j2 {
                a[s2][j] =
                    (a[s1][j] - a[s1][j - 1]) / ndu[(pk + 1) as usize][(rk + j as isize) as usize];
                dd += a[s2][j] * ndu[(rk + j as isize) as usize][pk as usize];
            }
            if r <= pk as usize {
                a[s2][k] = -a[s1][k - 1] / ndu[(pk + 1) as usize][r];
                dd += a[s2][k] * ndu[r][pk as usize];
            }
            ders[k][r] = dd;
            std::mem::swap(&mut s1, &mut s2);
        }
    }

    // Multiply by the factorial factors p, p(p−1), … .
    let mut acc = p;
    for (k, ders_k) in ders.iter_mut().enumerate().skip(1) {
        for d in ders_k.iter_mut() {
            *d *= acc as f64;
        }
        acc *= p.saturating_sub(k);
    }
    ders
}

#[cfg(test)]
mod tests {
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

    proptest! {
        #[test]
        fn derivative_matches_central_difference(c in nurbs(), s in 0.05f64..0.95) {
            let (lo, hi) = c.domain();
            let u = lo + s * (hi - lo);
            let h = 1e-6;
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
    }
}
