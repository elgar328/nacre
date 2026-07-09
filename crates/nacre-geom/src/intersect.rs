//! Closed-form intersections among analytic surfaces (design §3, §8 M5).
//!
//! Robustness-sensitive intersection code is isolated in this one module
//! (design §3). Everything here is the **construction (coordinate)** side of the
//! precision split (design §3): the returned `Line`/`Point3` are f64 *caches* —
//! the defining surfaces are the truth. So near-degenerate inputs are gated by a
//! conditioning threshold (there is no exact answer to cache), whereas exact
//! **sign** decisions live in `nacre-predicates`.

use crate::{Line, Plane};
use nacre_math::{Point3, Vector3};

/// The exact ray/segment vs triangle crossing outcomes, re-exported so callers
/// (`nacre-ops`) reach them through geom without depending on `nacre-predicates`
/// directly (the same `topo → geom → predicates` layering as the surface
/// handoffs).
pub use nacre_predicates::{RayCross, SegCross};

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

/// `|det|` below which three unit normals are too close to a common plane to
/// yield a usable vertex coordinate. `det = n1·(n2×n3) ∈ [−1, 1]` for unit
/// normals, so this is scale-free; and since the coordinate error is roughly
/// `ε_f64 / |det|`, the cutoff bounds that error (~`1e-7` here) rather than
/// being an arbitrary constant.
const COPLANAR_DET_EPS: f64 = 1e-9;

/// The point where three planes meet, or `None` if they do not meet in a
/// well-conditioned point (two parallel, or the three normals near-coplanar).
///
/// Closed-form Cramer: `P = (h1·(n2×n3) + h2·(n3×n1) + h3·(n1×n2)) / det`,
/// `det = n1·(n2×n3)` (the scalar triple product = determinant of rows
/// `[n1, n2, n3]`), `hᵢ = nᵢ·originᵢ`.
///
/// The f64 coordinate is a **cache** — the three planes are the truth (design
/// §3/§4). Its residual to the planes is the `Origin::Discovered` tolerance,
/// which the caller (M5-c) measures when it forms the vertex (closed-form
/// residual is recomputable, unlike an iterative solve's, so it is not returned
/// here). The gate is the conditioning threshold [`COPLANAR_DET_EPS`], not the
/// exact `det3_sign`: an exact-nonzero determinant can still round its f64
/// counterpart to ~0 and blow the coordinate up, and this function's output is
/// an inexact coordinate, so an exact existence test would be the wrong gate.
/// (The exact predicate is used where a *decision* must be exact — M5-c's
/// combinatorial "do these meet".)
pub fn three_planes(a: &Plane, b: &Plane, c: &Plane) -> Option<Point3> {
    let n1 = a.normal();
    let n2 = b.normal();
    let n3 = c.normal();
    let det = n1.dot(n2.cross(n3));
    if det.abs() <= COPLANAR_DET_EPS {
        return None; // parallel or near-coplanar normals — no usable vertex
    }
    let h1 = n1.dot(a.origin() - Point3::origin());
    let h2 = n2.dot(b.origin() - Point3::origin());
    let h3 = n3.dot(c.origin() - Point3::origin());
    let num = h1 * n2.cross(n3) + h2 * n3.cross(n1) + h3 * n1.cross(n2);
    Some(Point3::origin() + num / det)
}

/// The exact sign of `orient3d(V, q, r, s)`, where `V` is the implicit point at
/// which planes `a, b, c` meet and `q, r, s` are explicit points — the
/// geom→predicates handoff (design §9).
///
/// **Sign convention.** `orient3d(a,b,c,d) = det[a−d, b−d, c−d]`, so with `V`
/// first this is `(V − s) · ((q − s) × (r − s))`, and `(q−s)×(r−s) = (r−q)×(s−q)`
/// is the right-hand normal of the triangle `(q, r, s)`. Therefore
///
/// > `+1` ⇔ `V` lies on the **right-hand-normal side** of the triangle `(q, r, s)`,
/// > `-1` on the other side, `0` on its plane.
///
/// Callers that pass a face's outer-CCW triple get "positive = outside" for free,
/// because an outer-CCW loop's RH normal is the outward normal. Do not substitute
/// a `Plane`'s stored normal for the triangle's — the two need not agree in sign.
///
/// The three planes' [`coefficients`](Plane::coefficients) fill the rows of a
/// [`nacre_predicates::ThreePlane`], keeping the predicate dependency inside
/// geom (so `topo → geom → predicates`, and `topo` need not know predicates).
/// `V` is never materialized — the sign is exact even though `V`'s coordinates
/// are generally irrational (design §3, Attene 2020).
///
/// **Precondition (inherited from [`nacre_predicates::indirect_orient3d`]):**
/// `a, b, c` must meet in a single point (`det` of their normals ≠ 0), and
/// `q, r, s` must be non-collinear. A degenerate input yields `0`, which in a
/// release build is indistinguishable from a true coplanar `0` (the debug
/// assertion is compiled out) — callers must pass a valid vertex, as M5-c does.
pub fn three_plane_orient3d(
    a: &Plane,
    b: &Plane,
    c: &Plane,
    q: Point3,
    r: Point3,
    s: Point3,
) -> i8 {
    let tp = nacre_predicates::ThreePlane([a.coefficients(), b.coefficients(), c.coefficients()]);
    nacre_predicates::indirect_orient3d(&tp, q.as_array(), r.as_array(), s.as_array())
}

/// The exact sign of `(nA × nB) · nC` — how the line `A ∩ B` runs relative to
/// plane `C`'s normal.
///
/// The line is oriented the way [`plane_plane`] orients it (direction `nA × nB`),
/// so a caller may use this to reason about order along that line without ever
/// building it. `0` iff the three normals are coplanar, i.e. the planes have no
/// well-conditioned common point — the same condition that makes
/// [`three_planes`] return `None`, but decided exactly.
///
/// Exact: `det3_sign` of the three normals as rows, since
/// `det[nA; nB; nC] = nA · (nB × nC) = (nA × nB) · nC`. Not an f64 dot product —
/// the sign of a decision, so it belongs to the predicate side of the precision
/// split (design §3).
pub fn plane_pair_dir_sign(a: &Plane, b: &Plane, c: &Plane) -> i8 {
    nacre_predicates::det3_sign([
        a.normal().as_array(),
        b.normal().as_array(),
        c.normal().as_array(),
    ])
}

/// Exact forward-ray/triangle crossing for kernel types — the geom→predicates
/// handoff for [`nacre_predicates::ray_triangle_cross`] (design §8 M5-d
/// point-in-polyhedron). `d` is the ray direction; `tri` is a single triangle
/// (a caller triangulates a face into these). See [`RayCross`] for the outcome.
pub fn ray_face_cross(p: Point3, d: Vector3, tri: [Point3; 3]) -> RayCross {
    nacre_predicates::ray_triangle_cross(
        p.as_array(),
        d.as_array(),
        tri[0].as_array(),
        tri[1].as_array(),
        tri[2].as_array(),
    )
}

/// Exact segment/triangle crossing for kernel types — the handoff for
/// [`nacre_predicates::segment_triangle_cross`]. See [`SegCross`].
pub fn segment_face_cross(a: Point3, b: Point3, tri: [Point3; 3]) -> SegCross {
    nacre_predicates::segment_triangle_cross(
        a.as_array(),
        b.as_array(),
        tri[0].as_array(),
        tri[1].as_array(),
        tri[2].as_array(),
    )
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

    proptest! {
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
}
