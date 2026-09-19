//! NURBS curves — rational B-spline evaluation.
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

/// A rectangular grid of control points, `[i][j]` with `i` along u, `j` along v.
type ControlGrid = Vec<Vec<Point3>>;

/// A rational B-spline (NURBS) **surface**: a tensor product of the curve
/// construction in u and v. Same invariants as [`NurbsCurve`] in
/// each direction, over a rectangular control grid. **Not `Copy`** (heap grid).
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsSurface {
    degree_u: usize,
    degree_v: usize,
    control_points: ControlGrid,
    weights: Vec<Vec<f64>>,
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
}

impl NurbsSurface {
    /// Build a surface, validating all invariants. `None` if either degree is 0,
    /// the control grid is not rectangular or a dimension is below `degree + 1`,
    /// the weight grid does not match, any weight is not positive (rejects `0`,
    /// negatives, `NaN`), a knot count is wrong, or a knot vector decreases.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        control_points: ControlGrid,
        weights: Vec<Vec<f64>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
    ) -> Option<NurbsSurface> {
        if degree_u < 1 || degree_v < 1 {
            return None;
        }
        let n_u = control_points.len();
        if n_u < degree_u + 1 {
            return None;
        }
        let n_v = control_points[0].len();
        if n_v < degree_v + 1 {
            return None;
        }
        // Rectangular grid, matching weight grid, positive weights.
        if control_points.iter().any(|row| row.len() != n_v) {
            return None;
        }
        if weights.len() != n_u || weights.iter().any(|row| row.len() != n_v) {
            return None;
        }
        if !weights.iter().flatten().all(|&w| w > 0.0) {
            return None;
        }
        if knots_u.len() != n_u + degree_u + 1 || knots_v.len() != n_v + degree_v + 1 {
            return None;
        }
        if !knots_u.windows(2).all(|w| w[1] >= w[0]) || !knots_v.windows(2).all(|w| w[1] >= w[0]) {
            return None;
        }
        Some(NurbsSurface {
            degree_u,
            degree_v,
            control_points,
            weights,
            knots_u,
            knots_v,
        })
    }

    /// The u degree.
    #[inline]
    pub fn degree_u(&self) -> usize {
        self.degree_u
    }

    /// The v degree.
    #[inline]
    pub fn degree_v(&self) -> usize {
        self.degree_v
    }

    /// The control grid (`[i][j]`, i along u, j along v).
    #[inline]
    pub fn control_points(&self) -> &[Vec<Point3>] {
        &self.control_points
    }

    /// The weight grid (all positive).
    #[inline]
    pub fn weights(&self) -> &[Vec<f64>] {
        &self.weights
    }

    /// The u knot vector.
    #[inline]
    pub fn knots_u(&self) -> &[f64] {
        &self.knots_u
    }

    /// The v knot vector.
    #[inline]
    pub fn knots_v(&self) -> &[f64] {
        &self.knots_v
    }

    /// The parameter rectangle `((u0, u1), (v0, v1))` the surface spans.
    #[inline]
    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        let n_u = self.control_points.len() - 1;
        let n_v = self.control_points[0].len() - 1;
        (
            (self.knots_u[self.degree_u], self.knots_u[n_u + 1]),
            (self.knots_v[self.degree_v], self.knots_v[n_v + 1]),
        )
    }

    /// The point at `(u, v)` (clamped to the domain) — the rational tensor-product
    /// surface point (A4.3).
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let (u, v) = self.clamp(u, v);
        let (su, nu, sv, nv) = self.spans_and_bases(u, v);
        let (pu, pv) = (self.degree_u, self.degree_v);

        let mut num = [0.0; 3];
        let mut den = 0.0;
        for (k, &nuk) in nu.iter().enumerate() {
            let i = su - pu + k;
            for (l, &nvl) in nv.iter().enumerate() {
                let j = sv - pv + l;
                let cp = self.control_points[i][j].as_array();
                let nw = nuk * nvl * self.weights[i][j];
                num[0] += nw * cp[0];
                num[1] += nw * cp[1];
                num[2] += nw * cp[2];
                den += nw;
            }
        }
        Point3::from_array([num[0] / den, num[1] / den, num[2] / den])
    }

    /// The first partial derivative `∂S/∂u` at `(u, v)` (A4.4).
    pub fn du(&self, u: f64, v: f64) -> Vector3 {
        self.partial(u, v, 1, 0)
    }

    /// The first partial derivative `∂S/∂v` at `(u, v)` (A4.4).
    pub fn dv(&self, u: f64, v: f64) -> Vector3 {
        self.partial(u, v, 0, 1)
    }

    /// The unit surface normal `normalize(∂S/∂u × ∂S/∂v)` at `(u, v)`; `None` at a
    /// degenerate point (parallel or zero partials).
    pub fn normal_at(&self, u: f64, v: f64) -> Option<Vector3> {
        self.du(u, v).cross(self.dv(u, v)).normalize()
    }

    /// The `(a, b)`-order rational partial derivative (A4.4), for `a + b == 1`.
    fn partial(&self, u: f64, v: f64, a: usize, b: usize) -> Vector3 {
        let (u, v) = self.clamp(u, v);
        let (pu, pv) = (self.degree_u, self.degree_v);
        let su = find_span(pu, &self.knots_u, u);
        let sv = find_span(pv, &self.knots_v, v);
        let du = ders_basis_funs(pu, &self.knots_u, su, u, a);
        let dv = ders_basis_funs(pv, &self.knots_v, sv, v, b);

        // Homogeneous numerator/weight at orders (0,0) and (a,b).
        let (mut a00, mut a_ab) = ([0.0; 3], [0.0; 3]);
        let (mut w00, mut w_ab) = (0.0, 0.0);
        for (k, &d0uk) in du[0].iter().enumerate() {
            let i = su - pu + k;
            let dauk = du[a][k];
            for (l, &d0vl) in dv[0].iter().enumerate() {
                let j = sv - pv + l;
                let cp = self.control_points[i][j].as_array();
                let w = self.weights[i][j];
                let f00 = d0uk * d0vl * w;
                let fab = dauk * dv[b][l] * w;
                a00[0] += f00 * cp[0];
                a00[1] += f00 * cp[1];
                a00[2] += f00 * cp[2];
                a_ab[0] += fab * cp[0];
                a_ab[1] += fab * cp[1];
                a_ab[2] += fab * cp[2];
                w00 += f00;
                w_ab += fab;
            }
        }
        // S_ab = (A_ab − W_ab·S00) / W00, with S00 = A00 / W00.
        let s00 = [a00[0] / w00, a00[1] / w00, a00[2] / w00];
        Vector3::from_array([
            (a_ab[0] - w_ab * s00[0]) / w00,
            (a_ab[1] - w_ab * s00[1]) / w00,
            (a_ab[2] - w_ab * s00[2]) / w00,
        ])
    }

    fn spans_and_bases(&self, u: f64, v: f64) -> (usize, Vec<f64>, usize, Vec<f64>) {
        let su = find_span(self.degree_u, &self.knots_u, u);
        let nu = basis_funs(self.degree_u, &self.knots_u, su, u);
        let sv = find_span(self.degree_v, &self.knots_v, v);
        let nv = basis_funs(self.degree_v, &self.knots_v, sv, v);
        (su, nu, sv, nv)
    }

    #[inline]
    fn clamp(&self, u: f64, v: f64) -> (f64, f64) {
        let ((u0, u1), (v0, v1)) = self.domain();
        (u.clamp(u0, u1), v.clamp(v0, v1))
    }
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
        let mk =
            |du, dv, cp, w, ku: Vec<f64>, kv: Vec<f64>| NurbsSurface::new(du, dv, cp, w, ku, kv);
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
}
