//! `Ellipse` — the periodic carrier of a rim, the curve a plane meets a cylinder in.

use nacre_math::{Point3, Vector3};

/// **A rim's carrier — the curve a plane meets a cylinder in**, parameterized by the cylinder's
/// own angle θ (radians from `ref_dir` about `axis`): `P(θ) = center + r(cos θ·ref_dir + sin θ·
/// (axis × ref_dir)) + z(θ)·axis`, with `z(θ) = t₀ cos θ + t₁ sin θ`. A plane across the axis
/// leaves `z ≡ 0` and the curve is the circle of radius `r` about `center`; a plane oblique to it
/// tilts the circle into the ellipse it is the axial projection of — the same θ, so an arc's ends,
/// its order about the axis and the lateral chart's coordinate are one reading.
///
/// Stored: a `center` (where the axis meets the plane), a **unit** `axis`, a **unit** `ref_dir`
/// across it, the cylinder's `radius`, and the `tilt` `[t₀, t₁]` — `None` for a circle. Which one
/// a rim is comes from the exact relation of its plane to the cylinder's axis, not from the
/// numbers here, so the circle's maths runs on the circle's numbers, operation for operation, in
/// every method below.
///
/// Invariant: `axis` and `ref_dir` are unit and mutually orthogonal, `radius > 0`, and a `tilt`
/// is finite and not `[0, 0]` — [`Ellipse::from_center_normal`] normalizes, orthogonalizes, and
/// rejects a zero axis or non-positive radius; [`Ellipse::from_unit_frame`] takes a frame that
/// already is one; [`Ellipse::tilted_into`] tilts one onto a plane. Only those are stored; the
/// second in-plane axis (`axis × ref_dir`) is derived on demand so there is one fewer
/// orthonormality invariant to keep.
///
/// An `Ellipse` is the *full, unbounded* curve. An edge trims it to an arc via its endpoint
/// vertices — a closed edge is the whole turn — its two endpoint vertices are the same seam vertex
/// (`[v, v]`), the same carrier-vs-trim split as [`Line`](crate::Line) + `Edge`. A circle matches
/// the STEP `CIRCLE` (its `axis2_placement_3d` is exactly `center`/`axis`/`ref_dir`), an ellipse
/// the STEP `ELLIPSE` [`Ellipse::written`] gives.
///
/// Same `PartialEq` / no-`Eq` / no-`Hash` rationale as [`Line`](crate::Line):
/// exact `==` is for tests only; on-curve queries use [`Ellipse::distance`] /
/// [`Ellipse::contains`] with a caller-supplied tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ellipse {
    center: Point3,
    axis: Vector3,
    ref_dir: Vector3,
    radius: f64,
    tilt: Option<[f64; 2]>,
}

/// **The numbers a STEP file writes for an [`Ellipse`]** — a placement and two semi-axes, the
/// first along `ref_dir`. For a circle: its `center`, `axis`, `ref_dir`, and its radius twice. For
/// an ellipse: its `center`, its plane's unit normal on the side of the axis, the unit major
/// direction, the semi-major and the semi-minor (the cylinder's radius). Either way `ref_dir`'s
/// angle 0 runs counter-clockwise about `normal` the way θ runs about the axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Written {
    pub center: Point3,
    pub normal: Vector3,
    pub ref_dir: Vector3,
    pub semi_1: f64,
    pub semi_2: f64,
}

impl Ellipse {
    /// A circle in the plane through `center` perpendicular to `normal`, with the
    /// angle-0 direction taken from `ref_dir` (its component in the plane, via
    /// Gram–Schmidt) and the given `radius`.
    ///
    /// Returns `None` if `normal` is zero, `ref_dir` is parallel to `normal`
    /// (no in-plane component), or `radius` is not positive (rejects `0`,
    /// negatives, and `NaN`). Named `from_*` (not `new`) as a fallible
    /// constructor, mirroring [`Plane::from_point_normal`](crate::Plane). A frame that is
    /// already a correctly rounded unit pair takes [`Ellipse::from_unit_frame`] — normalizing it
    /// again moves its bits.
    pub fn from_center_normal(
        center: Point3,
        normal: Vector3,
        ref_dir: Vector3,
        radius: f64,
    ) -> Option<Ellipse> {
        let normal = normal.normalize()?;
        // Strip the normal component so the angle-0 axis lies in the plane.
        let ref_dir = (ref_dir - normal * ref_dir.dot(normal)).normalize()?;
        (radius > 0.0).then_some(Ellipse {
            center,
            axis: normal,
            ref_dir,
            radius,
            tilt: None,
        })
    }

    /// A circle whose `normal` and `ref_dir` **are already** the unit, mutually perpendicular
    /// vectors to keep, stored bit for bit — the frame a realization rounded correctly (a
    /// cylinder cache's axis and seam direction), which [`Ellipse::from_center_normal`] would round
    /// again: it normalizes the one and strips the other in `f64`, and a component whose truth is
    /// `0` comes back as a residue. The mirror of
    /// [`Cylinder::from_unit_frame`](crate::Cylinder::from_unit_frame). `None` for a radius that is
    /// not positive.
    #[inline]
    pub fn from_unit_frame(
        center: Point3,
        normal: Vector3,
        ref_dir: Vector3,
        radius: f64,
    ) -> Option<Ellipse> {
        debug_assert!(
            (normal.dot(normal) - 1.0).abs() < 1e-12
                && (ref_dir.dot(ref_dir) - 1.0).abs() < 1e-12
                && normal.dot(ref_dir).abs() < 1e-12,
            "a unit ref_dir across a unit normal: {normal:?} {ref_dir:?}"
        );
        (radius > 0.0).then_some(Ellipse {
            center,
            axis: normal,
            ref_dir,
            radius,
            tilt: None,
        })
    }

    /// **This circle tilted onto the plane through its `center` with unit `normal`** — the curve
    /// the plane meets the circle's cylinder in, at the same θ: `t = −r(n·ref_dir, n·(axis ×
    /// ref_dir)) / (n·axis)`. `None` where the plane holds the axis direction (`n·axis` is zero —
    /// no curve), where the tilt is not finite, or where it comes out `[0, 0]` or leaves a major
    /// direction that does not normalize (a plane the numbers cannot tell from one across the axis
    /// — the curve has no major direction to write), so every tilted curve has one.
    pub fn tilted_into(self, normal: Vector3) -> Option<Ellipse> {
        let na = normal.dot(self.axis);
        if na == 0.0 {
            return None;
        }
        let k = -self.radius / na;
        let tilt = [
            k * normal.dot(self.ref_dir),
            k * normal.dot(self.binormal()),
        ];
        if !tilt.iter().all(|t| t.is_finite()) || tilt == [0.0, 0.0] {
            return None;
        }
        let curve = Ellipse {
            tilt: Some(tilt),
            ..self
        };
        curve.major_direction().map(|_| curve)
    }

    /// The unit major direction of a tilted curve — the point at θ* = atan2(t₁, t₀), where
    /// `|P − c|² = r² + (t₀ cos θ + t₁ sin θ)²` peaks: `r(t₀ û + t₁ ŵ) + (t₀² + t₁²)·axis`. `None`
    /// for a circle (every direction is one) or a vector that does not normalize.
    fn major_direction(&self) -> Option<Vector3> {
        let [t0, t1] = self.tilt?;
        ((self.ref_dir * t0 + self.binormal() * t1) * self.radius + self.axis * (t0 * t0 + t1 * t1))
            .normalize()
    }

    /// The center.
    #[inline]
    pub fn center(&self) -> Point3 {
        self.center
    }

    /// The cylinder's unit axis — the direction θ turns about.
    #[inline]
    pub fn axis(&self) -> Vector3 {
        self.axis
    }

    /// The unit normal of the curve's plane, on the side of the axis (`plane_normal · axis > 0`)
    /// — the axis itself for a circle; for an ellipse `r·axis − t₀·ref_dir − t₁·(axis × ref_dir)`,
    /// normalized (the cross of the two conjugate semi-diameters at θ = 0 and 90°).
    #[inline]
    pub fn plane_normal(&self) -> Vector3 {
        match self.tilt {
            None => self.axis,
            Some([t0, t1]) => {
                let n = self.axis * self.radius - self.ref_dir * t0 - self.binormal() * t1;
                n * (1.0 / n.norm())
            }
        }
    }

    /// The unit angle-0 direction, across the axis.
    #[inline]
    pub fn ref_dir(&self) -> Vector3 {
        self.ref_dir
    }

    /// The cylinder's radius — the circle's, the ellipse's semi-minor.
    #[inline]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Whether this is a circle — the plane across the axis, as the exact relation decided.
    #[inline]
    pub fn is_circle(&self) -> bool {
        self.tilt.is_none()
    }

    /// `[t₀, t₁]`: the curve's height along the axis above `center` is `t₀ cos θ + t₁ sin θ` —
    /// `[0, 0]` for a circle.
    #[inline]
    pub fn axial_trig(&self) -> [f64; 2] {
        self.tilt.unwrap_or([0.0, 0.0])
    }

    /// The curve's height along the axis above `center` at `theta`.
    #[inline]
    pub fn z_at(self, theta: f64) -> f64 {
        match self.tilt {
            None => 0.0,
            Some([t0, t1]) => {
                let (sin, cos) = theta.sin_cos();
                t0 * cos + t1 * sin
            }
        }
    }

    /// The semi-major — the radius for a circle, `√(r² + t₀² + t₁²)` for an ellipse.
    #[inline]
    pub fn semi_major(&self) -> f64 {
        match self.tilt {
            None => self.radius,
            Some([t0, t1]) => (self.radius * self.radius + t0 * t0 + t1 * t1).sqrt(),
        }
    }

    /// The area the whole curve bounds in its plane — `πr²` for a circle, `πr·a` for an ellipse
    /// (`a` the semi-major).
    #[inline]
    pub fn area(&self) -> f64 {
        match self.tilt {
            None => std::f64::consts::PI * self.radius * self.radius,
            Some(_) => std::f64::consts::PI * self.radius * self.semi_major(),
        }
    }

    /// **How the curve bends over `[theta0, theta1]`** (`theta0 ≤ theta1`, unwrapped) — `(ρ, R)`:
    /// the largest radius a chord there sags by and the most the tangent turns per radian of θ
    /// there. `(r, 1)` for a circle everywhere; for an ellipse `(a, a/b)` over a range holding an
    /// end of the major axis (`a` the [`semi_major`](Ellipse::semi_major), `b` the cylinder's
    /// radius), less away from it — down to `(b, b/a)` over one holding an end of the minor axis,
    /// where the curve turns slower than its circle.
    ///
    /// The curve is the circle of radius `r` lifted linearly along the axis, θ kept, so a chord
    /// from θ₀ to θ₁ is parallel to the tangent at `(θ₀ + θ₁)/2` (a linear map keeps the circle's
    /// parallels). With `|P′(θ)|² = r²·q`, `q = 1 + (A/r)²·sin²(θ − θ*)` (`A = √(t₀² + t₁²)`, `θ*`
    /// where the major axis points), the tangent turns at `a/(r·q)` and the chord over a step `h`
    /// stands off the arc by `a(1 − cos(h/2))/√q` at its middle. Both are largest where `q` is
    /// least on the range — `1` where it holds `θ*` or `θ* + π`, else at its nearer end. So two
    /// consecutive chords inside the range turn by at most `R` times their two half-steps, and
    /// each sags at most `ρ(1 − cos(h/2))`.
    #[inline]
    pub fn pace_over(&self, theta0: f64, theta1: f64) -> (f64, f64) {
        match self.tilt {
            None => (self.radius, 1.0),
            Some([t0, t1]) => {
                let pi = std::f64::consts::PI;
                let major = t1.atan2(t0);
                let (lo, hi) = (theta0 - major, theta1 - major);
                let least = if (lo / pi).ceil() * pi <= hi {
                    0.0
                } else {
                    lo.sin().powi(2).min(hi.sin().powi(2))
                };
                let lift = (t0 * t0 + t1 * t1) / (self.radius * self.radius);
                let q = 1.0 + lift * least;
                let a = self.semi_major();
                (a / q.sqrt(), a / self.radius / q)
            }
        }
    }

    /// The half-extent of the whole curve along each world axis — `r·√(1 − n_i²)` for a circle
    /// (`n` its axis), `√(U_i² + V_i²)` for an ellipse with conjugate semi-diameters `U`, `V` at
    /// θ = 0 and 90°.
    #[inline]
    pub fn half_extent(&self) -> [f64; 3] {
        match self.tilt {
            None => {
                let n = self.axis.as_array();
                let r = self.radius;
                std::array::from_fn(|i| r * (1.0 - n[i] * n[i]).max(0.0).sqrt())
            }
            Some([t0, t1]) => {
                let u = (self.ref_dir * self.radius + self.axis * t0).as_array();
                let v = (self.binormal() * self.radius + self.axis * t1).as_array();
                std::array::from_fn(|i| (u[i] * u[i] + v[i] * v[i]).sqrt())
            }
        }
    }

    /// The numbers a STEP file writes for this curve ([`Written`]).
    pub fn written(&self) -> Written {
        match self.major_direction() {
            None => Written {
                center: self.center,
                normal: self.axis,
                ref_dir: self.ref_dir,
                semi_1: self.radius,
                semi_2: self.radius,
            },
            Some(major) => Written {
                center: self.center,
                normal: self.plane_normal(),
                ref_dir: major,
                semi_1: self.semi_major(),
                semi_2: self.radius,
            },
        }
    }

    /// The curve translated by `offset` — the center shifts; the axis, `ref_dir`, radius and tilt
    /// are unchanged.
    #[inline]
    pub fn translated(self, offset: Vector3) -> Ellipse {
        Ellipse {
            center: self.center + offset,
            ..self
        }
    }

    /// The second axis across the cylinder (angle 90°): `axis × ref_dir` (unit).
    #[inline]
    fn binormal(self) -> Vector3 {
        self.axis.cross(self.ref_dir)
    }

    /// The point at angle `theta` (radians):
    /// `center + r·(cos θ · ref_dir + sin θ · (axis × ref_dir)) + z(θ)·axis`.
    #[inline]
    pub fn point_at(self, theta: f64) -> Point3 {
        let (sin, cos) = theta.sin_cos();
        let round = self.center
            + self.ref_dir * (self.radius * cos)
            + self.binormal() * (self.radius * sin);
        match self.tilt {
            None => round,
            Some([t0, t1]) => round + self.axis * (t0 * cos + t1 * sin),
        }
    }

    /// **The farthest the curve from `theta0` to `theta1` stands off the chord between its two
    /// ends** — the distance from the chord's line of the curve at 65 even angles over the range,
    /// both ends and the middle among them (the middle is where a chord's sag is greatest —
    /// [`pace_over`](Ellipse::pace_over) — so the reading is that sag, not a sampled guess under
    /// it). A chord of no length reads the distance from its one point. Test-only (`test-util`):
    /// the bound [`pace_over`](Ellipse::pace_over) gives is held to it here, and a mesh's
    /// polylines are held to the chord budget with it.
    #[cfg(any(test, feature = "test-util"))]
    pub fn chord_standoff(self, theta0: f64, theta1: f64) -> f64 {
        let (a, b) = (self.point_at(theta0), self.point_at(theta1));
        let along = (b - a).normalize();
        (0..=64)
            .map(|i| {
                let p = self.point_at(theta0 + (theta1 - theta0) * f64::from(i) / 64.0) - a;
                match along {
                    Some(d) => (p - d * p.dot(d)).norm(),
                    None => p.norm(),
                }
            })
            .fold(0.0f64, f64::max)
    }

    /// The angle θ ∈ [0, 2π) of `p`'s direction from the axis — the inverse of
    /// [`point_at`](Ellipse::point_at)'s parameterization (`atan2` of the projections onto the
    /// two axes across the cylinder; the part along the axis drops out, so a tilt does not enter).
    /// `p` need not lie on the curve: the answer is the angle of its direction off the axis, and a
    /// point on the axis (both projections zero) returns `0.0` by `atan2`'s own convention rather
    /// than refusing — callers pass realized edge endpoints, which the store guarantees off-axis.
    ///
    /// This is what turns the arc convention (an edge's stored `[from, to]` order is
    /// CCW about the axis) into numbers: consumers take `Δθ = (θ_to − θ_from).rem_euclid(τ)`
    /// for two distinct ends — a closed edge is the whole turn, `τ` (`Model::derive_edge_curve`)
    /// — and never re-derive direction from anything else.
    #[inline]
    pub fn angle_of(self, p: Point3) -> f64 {
        let w = p - self.center;
        let theta = w.dot(self.binormal()).atan2(w.dot(self.ref_dir));
        theta.rem_euclid(std::f64::consts::TAU)
    }

    /// The area of the **segment** spanned by a CCW angle `dtheta` ∈ [0, 2π] — the region of the
    /// curve's plane between the chord and the arc: `r²(Δθ − sin Δθ)/2` for a circle; an ellipse is
    /// the circle's affine image along the axis, so its segment is the circle's times
    /// `|U × V| / r² = a/r` (`a` the semi-major). This is the one spelling of the number every arc
    /// consumer adds to a chord polygon (props' integrals, validate's winding witness, the
    /// boolean's `outer_tri`), about [`Ellipse::plane_normal`]; the sign — does this traversal add
    /// or remove the bulge — is the caller's, read from the `[from, to]`-CCW convention.
    #[inline]
    pub fn segment_area(self, dtheta: f64) -> f64 {
        match self.tilt {
            None => 0.5 * self.radius * self.radius * (dtheta - dtheta.sin()),
            Some(_) => 0.5 * self.radius * self.semi_major() * (dtheta - dtheta.sin()),
        }
    }

    /// The centroid of that segment, for the segment starting at angle `theta0` and spanning CCW
    /// `dtheta`. For a circle: on the bisector, `4r·sin³(Δθ/2) / (3(Δθ − sin Δθ))` from the center;
    /// a degenerate span (`Δθ − sin Δθ` ≈ 0) returns the chord midpoint's direction at distance
    /// `r` (the limit), rather than dividing by zero — its weight is zero anyway. For an ellipse:
    /// the circle's segment centroid carried along the axis onto the plane (a centroid follows an
    /// affine map).
    #[inline]
    pub fn segment_centroid(self, theta0: f64, dtheta: f64) -> Point3 {
        let denom = 3.0 * (dtheta - dtheta.sin());
        let dist = if denom > 0.0 {
            4.0 * self.radius * (0.5 * dtheta).sin().powi(3) / denom
        } else {
            self.radius
        };
        let flat = Ellipse { tilt: None, ..self };
        let mid = flat.point_at(theta0 + 0.5 * dtheta);
        let dir = (mid - self.center) * (1.0 / self.radius);
        let at = self.center + dir * dist;
        match self.tilt {
            None => at,
            Some([t0, t1]) => {
                let d = at - self.center;
                let z = (t0 * d.dot(self.ref_dir) + t1 * d.dot(self.binormal())) / self.radius;
                at + self.axis * z
            }
        }
    }

    /// The parametric derivative `dP/dθ` at `theta`:
    /// `r·(−sin θ · ref_dir + cos θ · (axis × ref_dir)) + z'(θ)·axis`.
    ///
    /// Never zero — it is the velocity along the curve, the quantity a finite difference of
    /// [`point_at`] returns (its magnitude is `radius` for a circle, **not** unit). (Term follows
    /// *The NURBS Book*.)
    ///
    /// [`point_at`]: Ellipse::point_at
    #[inline]
    pub fn derivative(self, theta: f64) -> Vector3 {
        let (sin, cos) = theta.sin_cos();
        let round = self.ref_dir * (-self.radius * sin) + self.binormal() * (self.radius * cos);
        match self.tilt {
            None => round,
            Some([t0, t1]) => round + self.axis * (t1 * cos - t0 * sin),
        }
    }

    /// Distance from `p` to the curve — the tolerance-free residual (the primitive, mirroring
    /// [`Line::distance`](crate::Line)).
    ///
    /// A circle's is exact: the axial/radial decomposition `√(a² + (ρ − r)²)` with `a = (p − c)·n`
    /// and `ρ = ‖(p − c) − a·n‖`, closed-form and safe on the axis (`ρ = 0`), where the
    /// nearest-point direction is undefined but the distance is not. An ellipse's nearest point
    /// solves a quartic, so this is the distance to the curve's point at `p`'s own angle
    /// ([`Ellipse::angle_of`]) — never less than the nearest distance, and equal to it for a point
    /// on the curve.
    #[inline]
    pub fn distance(self, p: Point3) -> f64 {
        match self.tilt {
            None => {
                let w = p - self.center;
                let axial = w.dot(self.axis);
                let radial = (w - self.axis * axial).norm();
                let dr = radial - self.radius;
                (axial * axial + dr * dr).sqrt()
            }
            Some(_) => (p - self.point_at(self.angle_of(p))).norm(),
        }
    }

    /// Whether `p` lies on the curve within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }
}

#[cfg(test)]
#[path = "tests/ellipse.rs"]
mod tests;
