//! `Circle` — a full circle in 3-D: the periodic carrier for circular edges.

use nacre_math::{Point3, Vector3};

/// A full circle, stored as a `center`, a **unit** `normal` (the plane it lies
/// in), a **unit** `ref_dir` in that plane (the angle-0 direction), and a
/// positive `radius`.
///
/// Invariant: `normal` and `ref_dir` are unit and mutually orthogonal, and
/// `radius > 0` — the constructor normalizes, orthogonalizes, and rejects a zero
/// axis or non-positive radius. Only `(center, normal, ref_dir, radius)` are
/// stored; the second in-plane axis (`normal × ref_dir`) is derived on demand so
/// there is one fewer orthonormality invariant to keep.
///
/// A `Circle` is the *full, unbounded* curve, parameterized by angle θ (radians)
/// from `ref_dir`. An edge trims it to an arc via its endpoint vertices — a
/// closed edge (`bounds: None`) is the whole circle — the same carrier-vs-trim
/// split as [`Line`](crate::Line) + `Edge` (design §3). This matches the STEP
/// `CIRCLE` (its `axis2_placement_3d` is exactly `center`/`normal`/`ref_dir`).
///
/// Same `PartialEq` / no-`Eq` / no-`Hash` rationale as [`Line`](crate::Line):
/// exact `==` is for tests only; on-circle queries use [`Circle::distance`] /
/// [`Circle::contains`] with a caller-supplied tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle {
    center: Point3,
    normal: Vector3,
    ref_dir: Vector3,
    radius: f64,
}

impl Circle {
    /// A circle in the plane through `center` perpendicular to `normal`, with the
    /// angle-0 direction taken from `ref_dir` (its component in the plane, via
    /// Gram–Schmidt) and the given `radius`.
    ///
    /// Returns `None` if `normal` is zero, `ref_dir` is parallel to `normal`
    /// (no in-plane component), or `radius` is not positive (rejects `0`,
    /// negatives, and `NaN`). Named `from_*` (not `new`) as a fallible
    /// constructor, mirroring [`Plane::from_point_normal`](crate::Plane).
    pub fn from_center_normal(
        center: Point3,
        normal: Vector3,
        ref_dir: Vector3,
        radius: f64,
    ) -> Option<Circle> {
        let normal = normal.normalize()?;
        // Strip the normal component so the angle-0 axis lies in the plane.
        let ref_dir = (ref_dir - normal * ref_dir.dot(normal)).normalize()?;
        (radius > 0.0).then_some(Circle {
            center,
            normal,
            ref_dir,
            radius,
        })
    }

    /// The center.
    #[inline]
    pub fn center(&self) -> Point3 {
        self.center
    }

    /// The unit normal of the circle's plane.
    #[inline]
    pub fn normal(&self) -> Vector3 {
        self.normal
    }

    /// The unit angle-0 direction (in-plane).
    #[inline]
    pub fn ref_dir(&self) -> Vector3 {
        self.ref_dir
    }

    /// The radius.
    #[inline]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The circle translated by `offset` — the center shifts, the plane
    /// (normal, ref_dir) and radius are unchanged.
    #[inline]
    pub fn translated(self, offset: Vector3) -> Circle {
        Circle {
            center: self.center + offset,
            normal: self.normal,
            ref_dir: self.ref_dir,
            radius: self.radius,
        }
    }

    /// The second in-plane axis (angle 90°): `normal × ref_dir` (unit).
    #[inline]
    fn binormal(self) -> Vector3 {
        self.normal.cross(self.ref_dir)
    }

    /// The point at angle `theta` (radians):
    /// `center + r·(cos θ · ref_dir + sin θ · (normal × ref_dir))`.
    #[inline]
    pub fn point_at(self, theta: f64) -> Point3 {
        let (sin, cos) = theta.sin_cos();
        self.center + self.ref_dir * (self.radius * cos) + self.binormal() * (self.radius * sin)
    }

    /// The angle θ ∈ [0, 2π) of `p`'s direction from the center — the inverse of
    /// [`point_at`](Circle::point_at)'s parameterization (`atan2` of the projections onto the
    /// two in-plane axes). `p` need not lie on the circle: the answer is the angle of its
    /// in-plane direction, and a point on the axis (both projections zero) returns `0.0` by
    /// `atan2`'s own convention rather than refusing — callers pass realized edge endpoints,
    /// which the store guarantees off-axis.
    ///
    /// ★ This is what turns the arc convention (M6-2b: an edge's stored `[from, to]` order is
    /// CCW about the axis) into numbers: consumers take `Δθ = (θ_to − θ_from).rem_euclid(τ)`
    /// and never re-derive direction from anything else.
    #[inline]
    pub fn angle_of(self, p: Point3) -> f64 {
        let w = p - self.center;
        let theta = w.dot(self.binormal()).atan2(w.dot(self.ref_dir));
        theta.rem_euclid(std::f64::consts::TAU)
    }

    /// The parametric derivative `dP/dθ` at `theta`:
    /// `r·(−sin θ · ref_dir + cos θ · (normal × ref_dir))`.
    ///
    /// Magnitude is `radius` (**not** unit) and never zero — it is the velocity
    /// along the curve, the quantity a finite difference of [`point_at`] returns.
    /// (Term follows *The NURBS Book*, the M3 reference.)
    ///
    /// [`point_at`]: Circle::point_at
    #[inline]
    pub fn derivative(self, theta: f64) -> Vector3 {
        let (sin, cos) = theta.sin_cos();
        self.ref_dir * (-self.radius * sin) + self.binormal() * (self.radius * cos)
    }

    /// Distance from `p` to the nearest point on the circle — the tolerance-free
    /// residual (the primitive, mirroring [`Line::distance`](crate::Line)).
    ///
    /// Axial/radial decomposition `√(a² + (ρ − r)²)` with `a = (p − c)·n` and
    /// `ρ = ‖(p − c) − a·n‖`. Closed-form and safe on the axis (`ρ = 0`), where
    /// the nearest-point direction is undefined but the distance is not.
    #[inline]
    pub fn distance(self, p: Point3) -> f64 {
        let w = p - self.center;
        let axial = w.dot(self.normal);
        let radial = (w - self.normal * axial).norm();
        let dr = radial - self.radius;
        (axial * axial + dr * dr).sqrt()
    }

    /// Whether `p` lies on the circle within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx_eq;
    use proptest::prelude::*;

    const EPS: f64 = 1e-9;

    fn unit_circle() -> Circle {
        Circle::from_center_normal(
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
        let c = Circle::from_center_normal(
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
        let c = Circle::from_center_normal(
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
        assert!(Circle::from_center_normal(c, Vector3::zero(), x, 1.0).is_none()); // zero normal
        assert!(Circle::from_center_normal(c, z, z, 1.0).is_none()); // ref_dir ∥ normal
        assert!(Circle::from_center_normal(c, z, x, 0.0).is_none()); // zero radius
        assert!(Circle::from_center_normal(c, z, x, -2.0).is_none()); // negative radius
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
    fn circle() -> impl Strategy<Value = Circle> {
        (pt3(), vec3(), vec3(), 0.5f64..10.0).prop_filter_map(
            "degenerate circle frame",
            |(center, normal, ref_dir, radius)| {
                let separated = matches!(
                    (normal.normalize(), ref_dir.normalize()),
                    (Some(n), Some(r)) if n.cross(r).norm() >= 1e-3
                );
                if separated {
                    Circle::from_center_normal(center, normal, ref_dir, radius)
                } else {
                    None
                }
            },
        )
    }

    proptest! {
        #[test]
        fn derivative_matches_central_difference(c in circle(), theta in -10.0f64..10.0) {
            // The design-mandated numerical-differentiation check. The derivative is
            // center-independent (impl uses only ref_dir/binormal/radius), so evaluate
            // it on an origin-centered copy — differencing `point_at` at a large center
            // would just add catastrophic cancellation from the huge shared term.
            let c0 =
                Circle::from_center_normal(Point3::origin(), c.normal(), c.ref_dir(), c.radius())
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
}
