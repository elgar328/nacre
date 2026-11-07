//! **Results with elliptical edges, read by OCCT** — a cylinder cut by planes oblique to its axis.
//! Each result is exported (its rims as `ELLIPSE`s), read back with `checkshape` on (a faulty read
//! is refused, `OracleError::Faulty`), and its volume and face count compared with nacre's.

use super::*;

/// A rod `r = 2` about `z` through `(0, 3)`, `z ∈ [0, 10]`, less the prism over the plane
/// `z = 5 + k·x` (its `zx` section between that plane and `z = 12` over `x ∈ [−3, 3]`, extruded over
/// `y ∈ [0, 6]`): the rod's end cut obliquely, one whole elliptical rim.
fn cut_rod(k: f64) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let rod = nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, 3.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.0,
        10.0,
    )
    .solid;
    let pts = [
        (5.0 - 3.0 * k, -3.0),
        (5.0 + 3.0 * k, 3.0),
        (12.0, 3.0),
        (12.0, -3.0),
    ]
    .map(|(u, v)| nacre_math::Point2::from_array([u, v]));
    let ring = nacre_ops::Ring2d::polygon_decimal(pts.to_vec()).expect("a quadrilateral");
    let profile = nacre_ops::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = SketchFrame::world(&m, Axis::Y);
    let Ok(OpOutput::Extrude { solid: prism, .. }) = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 6.0,
        },
    ) else {
        panic!("the prism extrudes")
    };
    let out = boolean(&mut m, BoolKind::Cut, rod, prism).expect("the rod cuts");
    (m, out[0])
}

/// A cylinder `r = 1/2`, length 2, on the rational frame of normal `(0.64, −0.48, 0.6)` at
/// `(1, 1, 0)`, against the box `[0, 2]³` — every plane of the box oblique to its axis.
fn tilted_against_box(kind: BoolKind) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [0.0; 3], [2.0; 3]);
    let frame = datum_frame(
        &mut m,
        SketchPlane::from_axes(
            Point3::from_array([1.0, 1.0, 0.0]),
            Vector3::from_array([0.6, 0.8, 0.0]),
            Vector3::from_array([-0.48, 0.36, 0.8]),
        ),
    );
    let ring = nacre_ops::Ring2d::circle(nacre_math::Point2::from_array([0.0, 0.0]), 0.5)
        .expect("a circle");
    let profile = nacre_ops::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let Ok(OpOutput::Extrude { solid: cyl, .. }) = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    ) else {
        panic!("the tilted cylinder extrudes")
    };
    let out = boolean(&mut m, kind, a, cyl).expect("the boolean builds");
    assert_eq!(out.len(), 1, "{kind:?}: one body");
    (m, out[0])
}

fn agrees(m: &Model, s: Handle<Solid>, what: &str) {
    let occt = crate::occt_props_of(m).unwrap_or_else(|e| panic!("{what}: {e:?}"));
    let ours = mass_props(m, s).expect("props");
    assert!(
        approx(occt.volume, ours.volume),
        "{what}: OCCT volume {} vs nacre {}",
        occt.volume,
        ours.volume
    );
    let sol = m.solid(s);
    let faces: usize = std::iter::once(sol.outer)
        .chain(sol.cavities.iter().copied())
        .map(|sh| m.shell(sh).faces.len())
        .sum();
    assert_eq!(occt.faces, faces, "{what}: faces");
}

/// The rod's oblique end, at two slopes: OCCT reads the `ELLIPSE` rim and its cap as a valid solid
/// of nacre's volume (`20π`) and faces.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_rods_oblique_end_reads_as_a_valid_solid() {
    for k in [0.5, 2.0] {
        let (m, s) = cut_rod(k);
        agrees(&m, s, &format!("k = {k}"));
        assert!((mass_props(&m, s).unwrap().volume - 20.0 * PI).abs() < 1e-9);
    }
}

/// The tilted cylinder against the box, the one-body results: the union, the box less the
/// cylinder, and the common part — the volume the set identities alone do not fix
/// (`every_kind_with_an_oblique_cylinder_builds`).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_tilted_cylinder_against_a_box_reads_as_valid_solids() {
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        let (m, s) = tilted_against_box(kind);
        agrees(&m, s, &format!("{kind:?}"));
    }
}
