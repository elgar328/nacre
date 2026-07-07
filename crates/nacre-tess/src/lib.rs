//! Mesh export for the nacre kernel.
//!
//! **M1 bootstrap.** Until the provenance-tagged tessellation layer exists
//! (design §5, M3), this crate offers only [`to_obj`]: it triangulates a
//! model's planar faces directly (fan over each face loop) and writes Wavefront
//! OBJ. It reads vertex coordinates and face connectivity from `nacre-topo`
//! only — it does not sample surfaces, so it is correct for the flat, convex
//! faces of M1 (e.g. a cuboid) but not curved geometry. In M3 this grows to
//! consume a real `Tessellation`.

use nacre_topo::Model;
use std::fmt::Write;

/// Export a model to Wavefront OBJ text (vertices shared; each face fan-
/// triangulated). Bootstrap — see the crate docs for the scope.
///
/// Assumes every edge is bounded (M1); a closed edge (`bounds: None`, M3) would
/// panic. Faces are emitted in their loop winding, which for M1's outward-wound
/// `Orientation::Forward` faces yields outward-facing triangles.
pub fn to_obj(model: &Model) -> String {
    let mut out = String::new();
    // Writing to a String is infallible; unwrap keeps `unused_must_use` quiet.
    writeln!(
        out,
        "# nacre OBJ export (bootstrap: direct planar triangulation)"
    )
    .unwrap();

    // Vertices in store order → OBJ indices 1..=n (matches each handle's index).
    for (_, v) in model.vertices.iter() {
        let [x, y, z] = v.point.as_array();
        writeln!(out, "v {} {} {}", x, y, z).unwrap();
    }

    for (_, face) in model.faces.iter() {
        // Ordered boundary as 1-based OBJ vertex indices.
        let boundary: Vec<u32> = face
            .outer
            .half_edges
            .iter()
            .map(|he| {
                let bounds = model
                    .edges
                    .get(he.edge)
                    .bounds
                    .expect("M1: bounded edges only (closed edges arrive in M3)");
                let start = if he.forward { bounds[0] } else { bounds[1] };
                start.index() + 1
            })
            .collect();
        // Fan triangulation — valid for convex faces (M1 faces are convex quads).
        for i in 1..boundary.len().saturating_sub(1) {
            writeln!(out, "f {} {} {}", boundary[0], boundary[i], boundary[i + 1]).unwrap();
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Point3;
    use proptest::prelude::*;
    use std::collections::HashSet;

    fn cube(min: [f64; 3], max: [f64; 3]) -> Model {
        let mut m = Model::new();
        m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
        m
    }

    fn v_lines(obj: &str) -> Vec<&str> {
        obj.lines().filter(|l| l.starts_with("v ")).collect()
    }
    fn f_lines(obj: &str) -> Vec<&str> {
        obj.lines().filter(|l| l.starts_with("f ")).collect()
    }
    fn parse3(rest: &str) -> Vec<f64> {
        rest.split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect()
    }

    #[test]
    fn unit_cube_obj_shape() {
        let obj = to_obj(&cube([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]));
        assert!(obj.lines().next().unwrap().starts_with('#'));
        assert_eq!(v_lines(&obj).len(), 8);

        let faces = f_lines(&obj);
        assert_eq!(faces.len(), 12); // 6 quads × 2 triangles
        let mut used = HashSet::new();
        for f in faces {
            let idx = parse3(&f[2..]);
            assert_eq!(idx.len(), 3);
            for &i in &idx {
                assert!((1.0..=8.0).contains(&i));
                used.insert(i as u32);
            }
        }
        assert_eq!(used.len(), 8); // every vertex referenced
    }

    #[test]
    fn vertices_round_trip() {
        let obj = to_obj(&cube([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]));
        let got: Vec<[f64; 3]> = v_lines(&obj)
            .into_iter()
            .map(|l| {
                let c = parse3(&l[2..]);
                [c[0], c[1], c[2]]
            })
            .collect();
        let expected = vec![
            [-2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 4.0, 0.0],
            [-2.0, 4.0, 0.0],
            [-2.0, 1.0, 10.0],
            [3.0, 1.0, 10.0],
            [3.0, 4.0, 10.0],
            [-2.0, 4.0, 10.0],
        ];
        assert_eq!(got, expected);
    }

    proptest! {
        #[test]
        fn arbitrary_box_obj_shape(
            min in prop::array::uniform3(-1e3f64..1e3),
            ext in prop::array::uniform3(1e-2f64..1e3),
        ) {
            let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
            let obj = to_obj(&cube(min, max));
            prop_assert_eq!(v_lines(&obj).len(), 8);
            let faces = f_lines(&obj);
            prop_assert_eq!(faces.len(), 12);
            for f in faces {
                let idx = parse3(&f[2..]);
                prop_assert_eq!(idx.len(), 3);
                for &i in &idx {
                    prop_assert!((1.0..=8.0).contains(&i));
                }
            }
        }
    }
}
