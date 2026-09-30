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
/// closed edge is the whole circle — its two endpoint vertices are the same seam
/// vertex (`[v, v]`), the same carrier-vs-trim split as [`Line`](crate::Line) +
/// `Edge`. This matches the STEP `CIRCLE` (its `axis2_placement_3d`
/// is exactly `center`/`normal`/`ref_dir`).
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
    /// ★ This is what turns the arc convention (an edge's stored `[from, to]` order is
    /// CCW about the axis) into numbers: consumers take `Δθ = (θ_to − θ_from).rem_euclid(τ)`
    /// and never re-derive direction from anything else.
    #[inline]
    pub fn angle_of(self, p: Point3) -> f64 {
        let w = p - self.center;
        let theta = w.dot(self.binormal()).atan2(w.dot(self.ref_dir));
        theta.rem_euclid(std::f64::consts::TAU)
    }

    /// The area of the **circular segment** spanned by a CCW angle `dtheta` ∈ [0, 2π] — the
    /// region between the chord and the arc: `r²(Δθ − sin Δθ)/2`. This is the one spelling of
    /// the number every arc consumer adds to a chord polygon (props' integrals, validate's
    /// winding witness, the boolean's `outer_tri`); the sign — does this traversal add or remove the bulge — is the
    /// caller's, read from the `[from, to]`-CCW convention.
    #[inline]
    pub fn segment_area(self, dtheta: f64) -> f64 {
        0.5 * self.radius * self.radius * (dtheta - dtheta.sin())
    }

    /// The centroid of that circular segment, for the segment starting at angle `theta0` and
    /// spanning CCW `dtheta`: on the bisector, `4r·sin³(Δθ/2) / (3(Δθ − sin Δθ))` from the
    /// center. A degenerate span (`Δθ − sin Δθ` ≈ 0) returns the chord midpoint's direction at
    /// distance `r` (the limit), rather than dividing by zero — its weight is zero anyway.
    #[inline]
    pub fn segment_centroid(self, theta0: f64, dtheta: f64) -> Point3 {
        let denom = 3.0 * (dtheta - dtheta.sin());
        let dist = if denom > 0.0 {
            4.0 * self.radius * (0.5 * dtheta).sin().powi(3) / denom
        } else {
            self.radius
        };
        let mid = self.point_at(theta0 + 0.5 * dtheta);
        let dir = (mid - self.center) * (1.0 / self.radius);
        self.center + dir * dist
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
#[path = "tests/circle.rs"]
mod tests;
