//! Closed-form intersections among analytic surfaces (design §3, §8 M5).
//!
//! Robustness-sensitive intersection code is isolated in this one module
//! (design §3). Everything here is the **construction (coordinate)** side of the
//! precision split (design §3): the returned `Line`/`Point3` are f64 *caches* —
//! the defining surfaces are the truth. So near-degenerate inputs are gated by a
//! conditioning threshold (there is no exact answer to cache), whereas exact
//! **sign** decisions live in `nacre-predicates`.

use crate::{Line, Plane};
use nacre_math::Point3;

/// `sin²θ` below which two plane normals count as parallel. Unit normals make
/// `‖n1 × n2‖² = sin²θ ∈ [0, 1]`, so this absolute cutoff is scale-free.
const PARALLEL_EPS: f64 = 1e-16;

/// The line where two planes meet, or `None` if they are parallel (or
/// coincident).
///
/// Closed form — no SSI march or spline cache (design §8 M5). Direction is
/// `n1 × n2`; the base point is the point of the line closest to the origin.
/// Returns a bare [`Line`] (the closed-form truth), not a cached
/// `Curve::Intersection` (that variant, holding `Handle<Surface>`, is for the
/// marched intersections of M7).
///
/// The base point is `p0 = (h1·(n2×d) + h2·(d×n1)) / (d·d)` with `d = n1×n2` and
/// `hi = ni·originᵢ`; one verifies `n1·p0 = h1` and `n2·p0 = h2` (using
/// `n1·(n2×d) = ‖d‖²` and `n1·(d×n1) = 0`), so `p0` lies on both planes, and
/// `p0 ⟂ d` makes it the closest to the origin.
pub fn plane_plane(a: &Plane, b: &Plane) -> Option<Line> {
    let n1 = a.normal();
    let n2 = b.normal();
    let d = n1.cross(n2);
    let dd = d.norm_squared();
    if dd <= PARALLEL_EPS {
        return None; // parallel or coincident — no unique line
    }
    // Plane i is `nᵢ·X = hᵢ` with `hᵢ = nᵢ·originᵢ`.
    let h1 = n1.dot(a.origin() - Point3::origin());
    let h2 = n2.dot(b.origin() - Point3::origin());
    let base = Point3::origin() + (h1 * n2.cross(d) + h2 * d.cross(n1)) / dd;
    // `dd > 0` guarantees `d` is nonzero, so this is always `Some`; keep the
    // `Option` rather than unwrapping (defensive).
    Line::from_point_direction(base, d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;
    use proptest::prelude::*;

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
}
