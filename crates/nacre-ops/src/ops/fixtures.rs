//! **Solids built through the product's own operations** — test-only (`test-util`).
//!
//! ★ A fixture here states nothing the product does not: it is a sequence of calls an application
//! makes (a datum plane, a circle sketched on it, an extrude). That is the difference from a
//! constructor that writes cells directly — such a door skips the realization funnel, and a
//! population built through it can hold what no operation produces (a tilted cap lifted from
//! realized `f64`, whose rim is an ellipse while its edge cache says circle).

use super::*;

/// A cylinder built by [`cylinder`] or [`cylinder_with_seam`], with its faces **named** — the
/// extrude's `[base, top, sides…]` order is a fact of `OpOutput::Extrude`, and a caller reading
/// a face by index would be guessing it.
#[derive(Clone, Copy, Debug)]
pub struct CylinderSolid {
    pub solid: Handle<Solid>,
    /// The cap on the sketch plane, facing against the axis.
    pub base: Handle<Face>,
    /// The far cap, facing along the axis.
    pub top: Handle<Face>,
    pub lateral: Handle<Face>,
}

/// A closed cylinder from `base` along `axis` for `height`, of `radius`: the plane through `base`
/// with normal `axis` stated as a datum ([`SketchPlane::from_origin_normal`]), a whole circle
/// about its origin, extruded. The seam sits where the product puts it — the frame's `+x̂`.
///
/// Panics on a statement the product refuses (a zero axis, a value outside the decimal window):
/// a fixture's input is the test's own.
pub fn cylinder(
    m: &mut Model,
    base: Point3,
    axis: Vector3,
    radius: f64,
    height: f64,
) -> CylinderSolid {
    let plane = SketchPlane::from_origin_normal(base, axis).expect("a nonzero cylinder axis");
    extruded_circle(m, plane, radius, height)
}

/// [`cylinder`] with the seam written out — for a test whose proposition involves **where the
/// seam is**. The sketch frame is `+x̂ = seam`, `+ŷ = axis × seam`
/// ([`SketchPlane::from_axes`]), so the circle's seam lands on the `seam` side of the axis and
/// the frame's normal is `axis`. `seam` must be perpendicular to `axis`; `axis.any_perpendicular()`
/// is the usual choice. Exact when the two and their cross product are the decimals written (an
/// axis-aligned pair always is).
pub fn cylinder_with_seam(
    m: &mut Model,
    base: Point3,
    axis: Vector3,
    seam: Vector3,
    radius: f64,
    height: f64,
) -> CylinderSolid {
    extruded_circle(
        m,
        SketchPlane::from_axes(base, seam, axis.cross(seam)),
        radius,
        height,
    )
}

fn extruded_circle(m: &mut Model, plane: SketchPlane, radius: f64, height: f64) -> CylinderSolid {
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating the cylinder's base plane: {other:?}"),
    };
    let ring = Ring2d::circle(Point2::from_array([0.0, 0.0]), radius).expect("a stated radius");
    let profile = Profile2d::from_normalized_rings(ring, vec![]);
    match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: height,
        },
    ) {
        Ok(OpOutput::Extrude { solid, faces }) => {
            let [base, top, lateral] = faces[..] else {
                panic!("a circle extrudes to three faces, got {}", faces.len())
            };
            CylinderSolid {
                solid,
                base,
                top,
                lateral,
            }
        }
        other => panic!("extruding the cylinder's circle: {other:?}"),
    }
}
