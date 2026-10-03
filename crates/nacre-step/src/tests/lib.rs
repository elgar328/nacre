use super::*;
use nacre_math::{Point3, Vector3};
use nacre_topo::Solid;
use step_io::read;
use step_io::scene::geometry::{CurveKind, SurfaceKind};

fn cuboid(min: [f64; 3], max: [f64; 3]) -> Model {
    let mut m = Model::new();
    nacre_ops::fixtures::cuboid(&mut m, Point3::from_array(min), Point3::from_array(max));
    m
}

#[test]
fn cube_round_trips_through_step_io_reader() {
    let text = to_step(&cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).expect("export");
    let (model, report) = read(text.as_bytes()).expect("re-read");

    // No dropped/orphan entities — the structural lint (step-loupe's report).
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
    // A solid-bearing part promotes to ABSR, not a plain SR.
    assert_eq!(
        model.advanced_brep_shape_representation_arena.items.len(),
        1
    );
    assert_eq!(model.shape_representation_arena.items.len(), 0);

    let scene = model.scene();
    let solids: Vec<_> = scene.all_solids().collect();
    assert_eq!(solids.len(), 1);
    let faces: Vec<_> = solids[0].faces().collect();
    assert_eq!(faces.len(), 6);
    for face in &faces {
        assert!(matches!(face.surface().kind(), SurfaceKind::Plane(_)));
        let bounds: Vec<_> = face.bounds().collect();
        assert_eq!(bounds.len(), 1);
        let edges: Vec<_> = bounds[0].oriented_edges().collect();
        assert_eq!(edges.len(), 4);
        for (edge, _forward) in edges {
            assert!(matches!(edge.curve().kind(), CurveKind::Line(_)));
        }
    }
}

#[test]
fn output_has_expected_entities_and_ap242e2_schema() {
    let text = to_step(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])).expect("export");
    for needle in [
        "MANIFOLD_SOLID_BREP",
        "CLOSED_SHELL",
        "ADVANCED_FACE",
        "PLANE",
        "LINE",
        "CARTESIAN_POINT",
    ] {
        assert!(text.contains(needle), "missing {needle}");
    }
    // AP242 Edition 2 schema token.
    assert!(text.contains("442 3 1 4"), "not AP242 Ed2");
}

#[test]
fn asymmetric_box_round_trips() {
    let text = to_step(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])).expect("export");
    let (_, report) = read(text.as_bytes()).expect("re-read");
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
}

#[test]
fn single_solid_export_isolates_one_solid() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([2.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );

    // The whole live model exports both solids.
    let both = to_step(&m).expect("export both");
    let (model_both, _) = read(both.as_bytes()).expect("re-read");
    assert_eq!(model_both.scene().all_solids().count(), 2);

    // A single-solid export isolates exactly one solid with its six faces.
    for h in [a, b] {
        let text = to_step_solid(&m, h).expect("export one");
        let (model, report) = read(text.as_bytes()).expect("re-read");
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
        let scene = model.scene();
        let solids: Vec<_> = scene.all_solids().collect();
        assert_eq!(solids.len(), 1);
        assert_eq!(solids[0].faces().count(), 6);
    }
}

/// ★★★★ **The live set's order reaches the exported file, so superseding may not permute it.**
///
/// `to_step` walks `Model::live_solids` in order. Nothing else in the tree would catch a
/// permutation — the census reads the arena (never live order) and there is no golden STEP
/// text — so the export side locks it here, beside the door's own lock in `nacre-topo`.
///
/// ⚠ The oracle is the **construction order and distinctive coordinates**, not a second call
/// to the same door: checking `supersede_live` against a hand-written `retain` would be one
/// implementation checking itself. The superseded solid's coordinate must also be *absent*,
/// which is what proves the drop happened at all.
#[test]
fn superseding_a_solid_leaves_the_export_order_alone() {
    let mut m = Model::new();
    nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([111.0; 3]),
        Point3::from_array([112.0; 3]),
    );
    let middle = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([777.0; 3]),
        Point3::from_array([778.0; 3]),
    );
    nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([333.0; 3]),
        Point3::from_array([334.0; 3]),
    );

    m.supersede_live(&[middle]);
    let text = to_step(&m).expect("export the two survivors");

    assert!(
        !text.contains("777."),
        "the superseded solid still exported"
    );
    let first = text.find("111.").expect("first survivor's coordinate");
    let third = text.find("333.").expect("second survivor's coordinate");
    assert!(
        first < third,
        "the export reordered the live set: 111. at {first}, 333. at {third}"
    );
}

#[test]
fn hollow_solid_round_trips_as_brep_with_voids() {
    // A 4-cube with a concentric 2-cube void, built with the cavity
    // producer's primitive (`reversed_shell`): the inner shell reversed
    // inward. Exports as a BREP_WITH_VOIDS and reads back with one cavity.
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0; 3]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0; 3]),
        Point3::from_array([3.0; 3]),
    );
    let b_outer = m.solid(b).outer;
    let void = m.reversed_shell(b_outer);
    let a_outer = m.solid(a).outer;
    let hollow = m.push_solid(Solid {
        outer: a_outer,
        cavities: vec![void],
    });
    m.restore_live(vec![hollow]); // supersede the two source cubes

    let text = to_step(&m).expect("export hollow");
    assert!(text.contains("BREP_WITH_VOIDS"), "no void entity in output");

    let (model, report) = read(text.as_bytes()).expect("re-read");
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
    let scene = model.scene();
    let solids: Vec<_> = scene.all_solids().collect();
    assert_eq!(solids.len(), 1);
    assert_eq!(solids[0].faces().count(), 6); // outer shell
    let voids = solids[0].voids();
    assert_eq!(voids.len(), 1);
    assert_eq!(voids[0].len(), 6); // one cavity shell, 6 faces
}

/// ★★ **A stated rim's centre reaches the file as the truth's nearest `f64`.** The edge cache
/// holds it (nacre-ops' `a_stated_rims_centre_is_the_truths_nearest`); this asks the next layer —
/// what each `CIRCLE`'s placement in the written file says, read back. Stated tilted axes off the
/// origin; the oracle is `base` and `base + h·û`, computed here in `Rat` from the fixture's
/// decimals and a unit axis rational by hand.
///
/// ★ The sweep also counts the rims whose centre the caches' own `f64` meet (`line_plane` of the
/// cylinder cache's axis and the cap cache) gets wrong, and requires one — a sweep where that
/// meet happens to be right everywhere cannot tell the two roads apart.
#[test]
fn a_stated_rims_centre_reaches_the_file_as_the_truths_nearest() {
    use nacre_exact::Rat;
    use step_io::generated::model::{Axis2PlacementRef, CartesianPointRef};
    let r = |n: i128, d: i128| Rat::new(n, d).expect("a rational");
    let axes: [([f64; 3], [Rat; 3]); 3] = [
        ([3.0, 4.0, 0.0], [r(3, 5), r(4, 5), r(0, 1)]),
        ([0.0, 3.0, 4.0], [r(0, 1), r(3, 5), r(4, 5)]),
        ([12.0, 0.0, 5.0], [r(12, 13), r(0, 1), r(5, 13)]),
    ];
    let bases: [[f64; 3]; 2] = [[0.1, 0.7, 0.3], [1.5, -2.25, 12.75]];
    let height = 1.1;
    let mut caches_miss = 0;
    for (axis, unit) in axes {
        for base in bases {
            let mut m = Model::new();
            nacre_ops::fixtures::cylinder(
                &mut m,
                Point3::from_array(base),
                Vector3::from_array(axis),
                0.3,
                height,
            );
            let b = base.map(|x| Rat::from_decimal(x).expect("a decimal"));
            let h = Rat::from_decimal(height).expect("a decimal");
            let top: [Rat; 3] = core::array::from_fn(|k| {
                b[k].checked_add(h.checked_mul(unit[k]).unwrap()).unwrap()
            });
            let mut want: Vec<[u64; 3]> = [b, top]
                .iter()
                .map(|p| p.map(|x| x.to_f64().to_bits()))
                .collect();
            want.sort();

            for eh in (0..m.edge_count() as u32).filter_map(|i| m.edge_handle_at(i)) {
                let nacre_geom::Curve::Circle(c) = m.edge_curve(eh) else {
                    continue;
                };
                let [s0, s1] = m.edge(eh).surfaces;
                let (cyl, cap) = match m.surface_cache(s0) {
                    nacre_geom::Surface::Cylinder(_) => (s0, s1),
                    nacre_geom::Surface::Plane(_) => (s1, s0),
                };
                let (nacre_geom::Surface::Cylinder(cc), nacre_geom::Surface::Plane(pc)) =
                    (m.surface_cache(cyl), m.surface_cache(cap))
                else {
                    unreachable!("a rim is a cylinder against a plane")
                };
                let meet = nacre_geom::intersect::line_plane(&cc.axis(), pc)
                    .expect("the axis crosses the cap");
                if meet.as_array().map(f64::to_bits) != c.center().as_array().map(f64::to_bits) {
                    caches_miss += 1;
                }
            }

            let text = to_step(&m).expect("export");
            let (model, report) = read(text.as_bytes()).expect("re-read");
            assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
            let mut got: Vec<[u64; 3]> = model
                .circle_arena
                .items
                .iter()
                .map(|c| {
                    let Axis2PlacementRef::Axis2Placement3d(p) = c.position else {
                        panic!("a 3-D circle has a 3-D placement")
                    };
                    let CartesianPointRef::CartesianPoint(o) =
                        model.axis2_placement3d_arena.get(p.0).location
                    else {
                        panic!("a placement's location is a cartesian point")
                    };
                    let xyz = &model.cartesian_point_arena.get(o.0).coordinates;
                    [xyz[0].to_bits(), xyz[1].to_bits(), xyz[2].to_bits()]
                })
                .collect();
            got.sort();
            got.dedup();
            assert_eq!(
                got, want,
                "axis {axis:?} base {base:?}: the file's rim centres"
            );
        }
    }
    assert!(
        caches_miss > 0,
        "the caches' f64 meet is right on every rim here — the lock cannot tell the roads apart"
    );
}

#[test]
fn cylinder_round_trips_through_step_io_reader() {
    let mut m = Model::new();
    nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    );
    let text = to_step(&m).expect("export cylinder");

    for needle in ["CYLINDRICAL_SURFACE", "CIRCLE", "442 3 1 4"] {
        assert!(text.contains(needle), "missing {needle}");
    }

    let (model, report) = read(text.as_bytes()).expect("re-read");
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);

    let scene = model.scene();
    let solids: Vec<_> = scene.all_solids().collect();
    assert_eq!(solids.len(), 1);
    let faces: Vec<_> = solids[0].faces().collect();
    assert_eq!(faces.len(), 3);

    let (mut cylindrical, mut planes) = (0, 0);
    for face in &faces {
        match face.surface().kind() {
            SurfaceKind::Cylindrical(_) => cylindrical += 1,
            SurfaceKind::Plane(_) => planes += 1,
            other => panic!("unexpected surface kind: {other:?}"),
        }
        for bound in face.bounds() {
            for (edge, _forward) in bound.oriented_edges() {
                assert!(matches!(
                    edge.curve().kind(),
                    CurveKind::Line(_) | CurveKind::Circle(_)
                ));
            }
        }
    }
    assert_eq!((cylindrical, planes), (1, 2));

    // The lateral (cylindrical) face's loop reuses the seam edge → 4 oriented edges.
    let lateral = faces
        .iter()
        .find(|f| matches!(f.surface().kind(), SurfaceKind::Cylindrical(_)))
        .unwrap();
    let edges: Vec<_> = lateral.bounds().next().unwrap().oriented_edges().collect();
    assert_eq!(edges.len(), 4);
}
