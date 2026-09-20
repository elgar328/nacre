use super::*;

use nacre_ops::{DatumDef, OpOutput, Operation, SketchFrame, SketchPlane, apply};

use nacre_exact::Axis;

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

use nacre_math::{Point3, Vector3};

use nacre_ops::boolean;

use nacre_ops::{BoolError, BoolKind};

use nacre_props::mass_props;

use std::f64::consts::PI;

/// Test shim: a boolean whose result is exactly one solid (a boolean may return several).
fn boolean_one(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    let solids = boolean(model, kind, a, b)?;
    assert_eq!(
        solids.len(),
        1,
        "boolean_one: expected one solid, got {}",
        solids.len()
    );
    Ok(solids[0])
}

fn approx(a: f64, b: f64) -> bool {
    let diff = (a - b).abs();
    diff <= 1e-6 || diff <= 1e-4 * a.abs().max(b.abs())
}

/// A concave L-prism (volume 3) plus a box. Rebuilt per test because `boolean`
/// supersedes its operands.
fn l_prism_and_box(lo: [f64; 3], hi: [f64; 3]) -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
    (m, l, bx)
}

/// The concave L-prism alone (volume 3).
fn l_prism() -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let l = extrude(
        &mut m,
        SketchPlane::world_xy(),
        &[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ],
    );
    (m, l)
}

/// Extrude a closed profile 1.0 along `plane`'s normal.
fn extrude(m: &mut Model, plane: SketchPlane, pts: &[[f64; 2]]) -> Handle<Solid> {
    extrude_dist(m, plane, pts, 1.0)
}

/// Extrude a closed profile `dist` along `plane`'s normal.
fn extrude_dist(m: &mut Model, plane: SketchPlane, pts: &[[f64; 2]], dist: f64) -> Handle<Solid> {
    use nacre_math::Point2;
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};
    let profile = Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).unwrap();
    let __f108 = datum_frame(m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __f108,
            profile,
            dist,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output");
    };
    solid
}

mod booleans;
mod nonconvex;
mod pockets_and_sweeps;
mod props;
mod rotated;
