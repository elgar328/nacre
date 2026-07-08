//! Exact mass properties (volume, surface area) of a nacre solid.
//!
//! Computed **analytically on the exact geometry** via the divergence theorem —
//! not from a tessellation, so the numbers match an analytic oracle (OCCT) to
//! near machine precision rather than to mesh resolution. This crate is a
//! read-only analysis over the topology [`Model`]; it lives above `nacre-topo`
//! (not inside it) so the topology crate stays a truth-only aggregate, mirroring
//! `nacre-validate` / `nacre-tess` / `nacre-step`.
//!
//! Consumers: user "part volume / surface area" queries, the M5 boolean
//! volume-conservation invariants (proptest), and the OCCT oracle diff
//! (`nacre-oracle`).
//!
//! Precondition: the solid is a valid closed, outward-oriented b-rep (what the
//! producers emit and `validate` accepts). Coordinates are the cache side of the
//! truth/cache split (design.md §0), so f64 arithmetic here is appropriate.

use nacre_geom::{Curve, Surface};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Face, Loop, Model, Orientation, Solid};

use std::f64::consts::PI;

/// Volume and surface area of one solid. (Centroid/inertia are deferred until
/// asymmetric solids make them worth the extra oracle plumbing — M4 face ops.)
///
/// No `PartialEq`: the fields are `f64`, and production geometry never compares
/// coordinates with a bare `==` (tests use a tolerance).
#[derive(Clone, Copy, Debug)]
pub struct MassProps {
    pub volume: f64,
    pub area: f64,
}

/// A shape this analysis does not yet cover.
///
/// There is deliberately **no** `UnsupportedSurface`: the surface `match` is
/// exhaustive over [`Surface`], so a new variant (NURBS, sphere…) is a compile
/// error here until it is handled — the same "new variant forces handling"
/// idiom `Surface::distance` and `to_step` use.
#[derive(Debug)]
pub enum MassError {
    /// A planar face bounded by something other than a straight polygon or a
    /// single full circle (e.g. a future line+arc mix). None occur today.
    UnsupportedBoundary,
    /// A **curved** face carries an inner loop (a hole). Planar holes are
    /// supported; a hole in a cylindrical face has no producer yet, so it is
    /// rejected rather than mis-integrated.
    UnsupportedInnerLoop,
    /// The solid has inner cavity shells; hollow-solid mass is deferred to the
    /// boolean milestones (M5+). No producer creates cavities today.
    HasCavities,
}

/// Exact mass properties of one solid via the divergence theorem.
///
/// `V = (1/3) ∮ (r − R)·n̂ dA`, area `= Σ face area`, summed face by face. The
/// local reference point `R` (a vertex of the solid) is subtracted before the
/// flux dot product: the closed-surface identity `∮ n̂ dA = 0` makes `V`
/// independent of `R`, while choosing `R` near the solid removes the
/// catastrophic cancellation a far-from-origin placement would otherwise cause.
pub fn mass_props(model: &Model, solid: Handle<Solid>) -> Result<MassProps, MassError> {
    let solid = model.solids.get(solid);
    if !solid.cavities.is_empty() {
        return Err(MassError::HasCavities);
    }
    let shell = model.shells.get(solid.outer);

    // R = the first vertex of the first face; any vertex on the solid works.
    let first_face = model.faces.get(shell.faces[0]);
    let reference = model
        .vertices
        .get(he_start(model, first_face.outer.half_edges[0])?)
        .point;

    let mut volume_flux = 0.0;
    let mut area = 0.0;
    for &face in &shell.faces {
        let (a, flux) = face_contribution(model, model.faces.get(face), reference)?;
        area += a;
        volume_flux += flux;
    }
    Ok(MassProps {
        volume: volume_flux / 3.0,
        area,
    })
}

/// `(area, raw volume flux)` of one face — flux is `∮ (r − R)·n̂ dA` without the
/// global `1/3`. Outward normal = surface natural normal × orientation sign.
fn face_contribution(
    model: &Model,
    face: &Face,
    reference: Point3,
) -> Result<(f64, f64), MassError> {
    let sign = match face.orientation {
        Orientation::Forward => 1.0,
        Orientation::Reversed => -1.0,
    };
    match model.surfaces.get(face.surface) {
        Surface::Plane(plane) => {
            // Outer boundary, minus each inner loop (a hole): area and first
            // moment are additive, so both the area and the flux subtract the
            // hole's contribution (divergence theorem, n̂ constant on a plane).
            let normal = plane.normal() * sign;
            let (mut area, centroid) = planar_face(model, &face.outer)?;
            let mut flux = normal.dot(centroid - reference) * area;
            for hole in &face.inner {
                let (a_in, c_in) = planar_face(model, hole)?;
                area -= a_in;
                flux -= normal.dot(c_in - reference) * a_in;
            }
            Ok((area, flux))
        }
        Surface::Cylinder(cyl) => {
            // No producer puts a hole in a curved face yet (imprint is planar).
            if !face.inner.is_empty() {
                return Err(MassError::UnsupportedInnerLoop);
            }
            // Full 2π lateral band (seam model A, design §9). Area = 2πr·h; the
            // flux integral ∮(r−R)·n̂ over the full band is 2πr²·h (the axial and
            // off-axis terms cancel over the closed angle), so it needs no R.
            let radius = cyl.radius();
            let axis = cyl.axis().direction();
            let height = axial_span(model, &face.outer, axis)?;
            let area = 2.0 * PI * radius * height;
            let flux = sign * 2.0 * PI * radius * radius * height;
            Ok((area, flux))
        }
    }
}

/// Area and area-weighted centroid of a planar face: a straight polygon, or a
/// single full circle (a cylinder cap).
fn planar_face(model: &Model, outer: &Loop) -> Result<(f64, Point3), MassError> {
    let all_lines = outer
        .half_edges
        .iter()
        .all(|he| matches!(edge_curve(model, *he), Curve::Line(_)));

    if all_lines {
        let points = loop_points(model, outer)?;
        Ok(polygon_area_centroid(&points))
    } else if outer.half_edges.len() == 1 {
        match edge_curve(model, outer.half_edges[0]) {
            Curve::Circle(circle) => {
                let area = PI * circle.radius() * circle.radius();
                Ok((area, circle.center()))
            }
            _ => Err(MassError::UnsupportedBoundary),
        }
    } else {
        Err(MassError::UnsupportedBoundary)
    }
}

/// Area and area-weighted centroid of a simple (convex or concave) planar
/// polygon, via a signed triangle fan from the first vertex.
fn polygon_area_centroid(points: &[Point3]) -> (f64, Point3) {
    let base = points[0];
    // Area vector A_vec = ½ Σ (vᵢ − v₀) × (vᵢ₊₁ − v₀); |A_vec| is the true area.
    let mut area_vec = Vector3::zero();
    for pair in points[1..].windows(2) {
        area_vec += (pair[0] - base).cross(pair[1] - base);
    }
    let unit = area_vec.normalize().unwrap_or(Vector3::zero());

    // Area-weighted centroid: signed triangle areas (about `unit`) weight each
    // triangle centroid. Σ signed-2area = |A_vec|, so this is exact for concave
    // faces too — a plain vertex average would be wrong.
    let mut weighted = Vector3::zero();
    let mut weight = 0.0;
    for pair in points[1..].windows(2) {
        let tri = (pair[0] - base).cross(pair[1] - base);
        let signed = tri.dot(unit);
        let centroid_rel = ((pair[0] - base) + (pair[1] - base)) * (1.0 / 3.0);
        weighted += centroid_rel * signed;
        weight += signed;
    }
    let area = 0.5 * area_vec.norm();
    let centroid = base + weighted * (1.0 / weight);
    (area, centroid)
}

/// Axial extent of a face's loop: the span of its vertices projected on `axis`
/// (the two seam vertices for a cylinder band → the exact height).
fn axial_span(model: &Model, outer: &Loop, axis: Vector3) -> Result<f64, MassError> {
    let points = loop_points(model, outer)?;
    // Project relative to the first point to keep the numbers small.
    let origin = points[0];
    let mut min = 0.0;
    let mut max = 0.0;
    for &p in &points[1..] {
        let t = (p - origin).dot(axis);
        min = f64::min(min, t);
        max = f64::max(max, t);
    }
    Ok(max - min)
}

/// The ordered start vertices of a loop's half-edges (the boundary polyline).
fn loop_points(model: &Model, outer: &Loop) -> Result<Vec<Point3>, MassError> {
    outer
        .half_edges
        .iter()
        .map(|&he| Ok(model.vertices.get(he_start(model, he)?).point))
        .collect()
}

/// The start vertex of a half-edge (`bounds[0]` if forward, else `bounds[1]`).
/// A loop edge with no endpoints (`bounds: None`) is not a valid solid boundary.
fn he_start(
    model: &Model,
    he: nacre_topo::HalfEdge,
) -> Result<Handle<nacre_topo::Vertex>, MassError> {
    let bounds = model
        .edges
        .get(he.edge)
        .bounds
        .ok_or(MassError::UnsupportedBoundary)?;
    Ok(if he.forward { bounds[0] } else { bounds[1] })
}

/// The `Curve` carried by a half-edge's edge.
fn edge_curve(model: &Model, he: nacre_topo::HalfEdge) -> &Curve {
    model.curves.get(model.edges.get(he.edge).curve)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::{Point2, Vector3};
    use nacre_ops::{OpOutput, Operation, Profile2d, SketchPlane, apply};

    /// Relative-or-absolute comparison sized for accumulated f64 error.
    fn close(got: f64, expected: f64) -> bool {
        (got - expected).abs() <= 1e-9 * expected.abs().max(1.0)
    }

    // --- golden ---

    #[test]
    fn cube_mass_matches_analytic() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let props = mass_props(&m, s).unwrap();
        assert!(close(props.volume, 24.0), "volume {}", props.volume); // 2·3·4
        assert!(close(props.area, 52.0), "area {}", props.area); // 2(6+8+12)
    }

    #[test]
    fn cylinder_mass_matches_analytic() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        let props = mass_props(&m, s).unwrap();
        assert!(close(props.volume, 20.0 * PI), "volume {}", props.volume); // πr²h
        assert!(close(props.area, 28.0 * PI), "area {}", props.area); // 2πr² + 2πrh
    }

    #[test]
    fn concave_extrude_mass() {
        // An L-shaped profile: a 2×2 square with the top-right 1×1 corner removed.
        // Concave, so the area-weighted centroid differs from the vertex average —
        // a vertex-average bug in the flux would fail the volume assertion.
        let pts = [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let profile = Profile2d {
            points: pts.iter().map(|&p| Point2::from_array(p)).collect(),
        };
        let dist = 3.0;
        let mut m = Model::new();
        let op = Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist,
        };
        let OpOutput::Extrude { solid: s, .. } = apply(&mut m, &op).unwrap() else {
            unreachable!("extrude yields Extrude output");
        };

        // Independent expected base area via the 2D shoelace (not mass_props).
        let a_l = shoelace(&pts);
        assert!((a_l - 3.0).abs() < 1e-12); // 4 − 1

        let props = mass_props(&m, s).unwrap();
        assert!(close(props.volume, a_l * dist), "volume {}", props.volume);
    }

    /// Build a `size` cube, mass its solid, then imprint a `2·hw` square hole in
    /// its top face and mass the result. Returns `(before, after)`.
    fn cube_then_imprint(size: f64, hw: f64) -> (MassProps, MassProps) {
        let sq = |s: f64| Profile2d {
            points: [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        };
        let mut m = Model::new();
        let OpOutput::Extrude { solid, faces } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq(size),
                dist: size,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let before = mass_props(&m, solid).unwrap();
        let hole = Profile2d {
            points: [[-hw, -hw], [hw, -hw], [hw, hw], [-hw, hw]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        };
        let OpOutput::ImprintSketch { solid: new, .. } = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: faces[1],
                profile: hole,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        (before, mass_props(&m, new).unwrap())
    }

    #[test]
    fn imprint_preserves_mass() {
        // Imprint is a coplanar subdivision: the hole cut from the outer face is
        // filled by the region face, so volume and area are unchanged. This
        // exercises the inner-loop subtraction (the outer face now has a hole).
        let (before, after) = cube_then_imprint(1.0, 0.2);
        assert!(close(before.volume, 1.0));
        assert!(close(before.area, 6.0));
        assert!(close(after.volume, before.volume), "vol {}", after.volume);
        assert!(close(after.area, before.area), "area {}", after.area);
    }

    /// Pad a `2·hw` square boss of height `dist` on a `size` cube's top face,
    /// returning the padded solid's mass.
    fn cube_then_pad(size: f64, hw: f64, dist: f64) -> MassProps {
        let sq = |s: f64| Profile2d {
            points: [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        };
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq(size),
                dist: size,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let boss = Profile2d {
            points: [[-hw, -hw], [hw, -hw], [hw, hw], [-hw, hw]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        };
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: faces[1],
                profile: boss,
                dist,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        mass_props(&m, solid).unwrap()
    }

    #[test]
    fn pad_boss_mass() {
        // Unit cube + a 0.4-square boss (area 0.16, perimeter 1.6) of height 0.5.
        // Volume = 1 + 0.16·0.5 = 1.08; area = 6 + 1.6·0.5 = 6.8. This is the
        // *absolute* test of the inner-loop subtraction (the hole is not filled).
        let m = cube_then_pad(1.0, 0.2, 0.5);
        assert!(close(m.volume, 1.08), "vol {}", m.volume);
        assert!(close(m.area, 6.8), "area {}", m.area);
    }

    /// Carve a `2·hw` square pocket of depth `dist` into a `size` cube's top face.
    fn cube_then_pocket(size: f64, hw: f64, dist: f64) -> MassProps {
        let sq = |s: f64| Profile2d {
            points: [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        };
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq(size),
                dist: size,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let pocket = Profile2d {
            points: [[-hw, -hw], [hw, -hw], [hw, hw], [-hw, hw]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        };
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: pocket,
                dist,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        mass_props(&m, solid).unwrap()
    }

    #[test]
    fn pocket_mass() {
        // Unit cube − a 0.4-square pocket of depth 0.5. Volume = 1 − 0.16·0.5 =
        // 0.92 (the inward walls contribute negatively); area = 6 + 1.6·0.5 = 6.8
        // (same as the boss). Exercises the hole subtraction with inward walls.
        let m = cube_then_pocket(1.0, 0.2, 0.5);
        assert!(close(m.volume, 0.92), "vol {}", m.volume);
        assert!(close(m.area, 6.8), "area {}", m.area);
    }

    /// Unsigned area of a 2D polygon (independent check for the concave test).
    fn shoelace(pts: &[[f64; 2]]) -> f64 {
        let mut two_area = 0.0;
        for i in 0..pts.len() {
            let a = pts[i];
            let b = pts[(i + 1) % pts.len()];
            two_area += a[0] * b[1] - b[0] * a[1];
        }
        0.5 * two_area.abs()
    }

    // --- proptest ---

    use proptest::prelude::*;

    fn wide_point() -> impl Strategy<Value = Point3> {
        prop::array::uniform3(-1e6f64..1e6).prop_map(Point3::from_array)
    }

    proptest! {
        /// Random box, possibly far from the origin — exercises the R
        /// cancellation guard as well as the polygon path.
        #[test]
        fn prop_cuboid_mass(
            min in wide_point(),
            a in 0.1f64..1e3, b in 0.1f64..1e3, c in 0.1f64..1e3,
        ) {
            let max = min + Vector3::from_array([a, b, c]);
            let mut m = Model::new();
            let s = m.add_cuboid(min, max);
            let props = mass_props(&m, s).unwrap();
            prop_assert!(close(props.volume, a * b * c), "vol {} vs {}", props.volume, a*b*c);
            prop_assert!(close(props.area, 2.0 * (a*b + b*c + c*a)), "area {}", props.area);
        }

        /// Random cylinder with an arbitrary axis direction — exercises the
        /// non-axis-aligned frame, the lateral sign, and the 2π band height.
        #[test]
        fn prop_cylinder_mass(
            base in wide_point(),
            dir in prop::array::uniform3(-1e3f64..1e3)
                .prop_map(Vector3::from_array)
                .prop_filter("axis must be clearly nonzero", |v| v.norm() > 1e-3),
            r in 0.5f64..10.0, h in 0.5f64..1e3,
        ) {
            let mut m = Model::new();
            let s = m.add_cylinder(base, dir, r, h);
            let props = mass_props(&m, s).unwrap();
            prop_assert!(close(props.volume, PI * r * r * h), "vol {} vs {}", props.volume, PI*r*r*h);
            prop_assert!(close(props.area, 2.0 * PI * r * (r + h)), "area {}", props.area);
        }

        /// Imprinting a hole preserves volume and area for any (interior) hole
        /// size. `hw ≤ 0.3` keeps the hole's circumradius `hw√2 ≈ 0.42` below the
        /// top face's inradius `size/2 ≥ 0.5`, so the profile stays interior.
        #[test]
        fn prop_imprint_preserves_mass(size in 1.0f64..5.0, hw in 0.05f64..0.3) {
            let (before, after) = cube_then_imprint(size, hw);
            prop_assert!(close(after.volume, before.volume), "vol {} vs {}", after.volume, before.volume);
            prop_assert!(close(after.area, before.area), "area {} vs {}", after.area, before.area);
        }

        /// A boss adds `A_p·dist` of volume and `P·dist` of surface (the top hole
        /// area cancels the boss cap). For a `2·hw` square: `A_p = 4hw²`, `P = 8hw`.
        /// This absolutely exercises the inner-loop subtraction (uncancelled hole).
        #[test]
        fn prop_pad_volume(size in 1.0f64..5.0, hw in 0.05f64..0.3, dist in 0.1f64..5.0) {
            let m = cube_then_pad(size, hw, dist);
            let vol = size * size * size + 4.0 * hw * hw * dist;
            let area = 6.0 * size * size + 8.0 * hw * dist;
            prop_assert!(close(m.volume, vol), "vol {} vs {}", m.volume, vol);
            prop_assert!(close(m.area, area), "area {} vs {}", m.area, area);
        }

        /// A pocket *removes* `A_p·dist` of volume (inward walls) while adding the
        /// same `P·dist` of surface as a boss. `dist ≤ 0.9 < size` keeps it blind.
        #[test]
        fn prop_pocket_volume(size in 1.0f64..5.0, hw in 0.05f64..0.3, dist in 0.1f64..0.9) {
            let m = cube_then_pocket(size, hw, dist);
            let vol = size * size * size - 4.0 * hw * hw * dist;
            let area = 6.0 * size * size + 8.0 * hw * dist;
            prop_assert!(close(m.volume, vol), "vol {} vs {}", m.volume, vol);
            prop_assert!(close(m.area, area), "area {} vs {}", m.area, area);
        }
    }
}
