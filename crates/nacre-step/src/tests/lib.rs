use super::*;
use nacre_math::{Point3, Vector3};
use nacre_topo::Solid;
use step_io::read;
use step_io::scene::geometry::{CurveKind, SurfaceKind};

/// The header time stamp these tests write — fixed, so a file is a function of its model.
const STAMP: &str = "2026-10-05T00:00:00Z";

fn cuboid(min: [f64; 3], max: [f64; 3]) -> Model {
    let mut m = Model::new();
    nacre_ops::fixtures::cuboid(&mut m, Point3::from_array(min), Point3::from_array(max));
    m
}

#[test]
fn cube_round_trips_through_step_io_reader() {
    let text = to_step(&cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), STAMP).expect("export");
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
    let text = to_step(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]), STAMP).expect("export");
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
    let text = to_step(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]), STAMP).expect("export");
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
    let both = to_step(&m, STAMP).expect("export both");
    let (model_both, _) = read(both.as_bytes()).expect("re-read");
    assert_eq!(model_both.scene().all_solids().count(), 2);

    // A single-solid export isolates exactly one solid with its six faces.
    for h in [a, b] {
        let text = to_step_solids(&m, &[h], STAMP).expect("export one");
        let (model, report) = read(text.as_bytes()).expect("re-read");
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
        let scene = model.scene();
        let solids: Vec<_> = scene.all_solids().collect();
        assert_eq!(solids.len(), 1);
        assert_eq!(solids[0].faces().count(), 6);
    }

    // A chosen pair is both, and an empty choice is a file with no solid in it.
    let pair = to_step_solids(&m, &[a, b], STAMP).expect("export the pair");
    let (model_pair, _) = read(pair.as_bytes()).expect("re-read");
    assert_eq!(model_pair.scene().all_solids().count(), 2);
    let none = to_step_solids(&m, &[], STAMP).expect("export nothing");
    let (model_none, report) = read(none.as_bytes()).expect("re-read");
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
    assert_eq!(model_none.scene().all_solids().count(), 0);
}

/// **The header's time stamp is the caller's, verbatim — so a file is a function of its model
/// and that stamp.** The kernel reads no clock: `wasm32-unknown-unknown` has none (asking the
/// system for the time there panics), and a stamp taken from the clock would make two exports of
/// one model differ. Two stamps prove the field is the argument rather than a coincidence of the
/// moment, and two models built the same way give the same bytes, header included.
#[test]
fn the_time_stamp_is_the_callers_and_the_file_repeats() {
    let lo = [-2.0, 1.0, 0.0];
    let hi = [3.0, 4.0, 10.0];
    let m = cuboid(lo, hi);
    let text = to_step(&m, STAMP).expect("export");
    assert!(
        text.contains(&format!("'{STAMP}'")),
        "the stamp is not in the header"
    );
    let other = "1999-12-31T23:59:59Z";
    let text_other = to_step(&m, other).expect("export");
    assert!(
        text_other.contains(&format!("'{other}'")),
        "the stamp is not the argument"
    );
    assert_eq!(text, to_step(&cuboid(lo, hi), STAMP).expect("export"));
}

/// **The live set's order reaches the exported file, so superseding may not permute it.**
///
/// `to_step` walks `Model::live_solids` in order. Nothing else in the tree would catch a
/// permutation — the census reads the arena (never live order) and there is no golden STEP
/// text — so the export side locks it here, beside the door's own lock in `nacre-topo`.
///
/// The oracle is the **construction order and distinctive coordinates**, not a second call
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
    let text = to_step(&m, STAMP).expect("export the two survivors");

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

/// **A superseded solid named to `to_step_solids` is written as the stores hold it.** The live set
/// is the current model — `to_step` leaves the superseded box out (the test above) — but the arena
/// never edits a cell, so the box is still a whole shape of the model's history, and a caller who
/// names it gets that shape: its coordinate, its six faces, a file that reads back clean.
#[test]
fn a_superseded_solid_named_to_the_export_is_written_whole() {
    let mut m = Model::new();
    nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([111.0; 3]),
        Point3::from_array([112.0; 3]),
    );
    let gone = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([777.0; 3]),
        Point3::from_array([778.0; 3]),
    );
    m.supersede_live(&[gone]);

    let text = to_step_solids(&m, &[gone], STAMP).expect("export the superseded box");
    assert!(
        text.contains("777."),
        "the superseded box's coordinate is missing"
    );
    assert!(!text.contains("111."), "the live box was written too");
    let (model, report) = read(text.as_bytes()).expect("re-read");
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
    let scene = model.scene();
    let solids: Vec<_> = scene.all_solids().collect();
    assert_eq!(solids.len(), 1);
    assert_eq!(solids[0].faces().count(), 6);
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

    let text = to_step(&m, STAMP).expect("export hollow");
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

/// **A stated rim's centre reaches the file as the truth's nearest `f64`.** The edge cache
/// holds it (nacre-ops' `a_stated_rims_centre_is_the_truths_nearest`); this asks the next layer —
/// what each `CIRCLE`'s placement in the written file says, read back. Stated tilted axes off the
/// origin; the oracle is `base` and `base + h·û`, computed here in `Rat` from the fixture's
/// decimals and a unit axis rational by hand.
///
/// The sweep also counts the rims whose centre the caches' own `f64` meet (`line_plane` of the
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
                    nacre_topo::SurfaceCache::Cylinder { .. } => (s0, s1),
                    nacre_topo::SurfaceCache::Plane { .. } => (s1, s0),
                };
                let (Some(cc), Some(pc)) = (
                    m.surface_cache(cyl).cylinder(),
                    m.surface_cache(cap).plane(),
                ) else {
                    unreachable!("a rim is a cylinder against a plane")
                };
                let meet = nacre_geom::intersect::line_plane(&cc.axis(), pc)
                    .expect("the axis crosses the cap");
                if meet.as_array().map(f64::to_bits) != c.center().as_array().map(f64::to_bits) {
                    caches_miss += 1;
                }
            }

            let text = to_step(&m, STAMP).expect("export");
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

/// **A straight edge's direction reaches the file as its edge cache's.** The cache holds the
/// truth's direction (nacre-ops' census line lock); this asks the written `LINE`s, read back: each
/// `DIRECTION` is a cached direction bit for bit and each `VECTOR` magnitude is 1. Boxes turned by
/// whole degrees that are not quarters, requiring one edge whose endpoint-derived direction
/// differs from its cache — where writing the endpoints' difference would put another direction
/// in the file.
#[test]
fn a_straight_edges_direction_reaches_the_file_as_its_cache() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    use step_io::generated::model::{DirectionRef, VectorRef};
    let bits = |v: [f64; 3]| v.map(f64::to_bits);
    let mut endpoints_differ = 0;
    for (name, axis, deg) in [("37° about y", Axis::Y, 37), ("23° about x", Axis::X, 23)] {
        let mut m = Model::new();
        let s = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.1, 0.7, 0.3]),
            Point3::from_array([1.2, 1.5, 2.4]),
        );
        let turn = Isometry::rotation(Rotation {
            axis,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).expect("a whole-degree angle"),
        });
        let Ok(nacre_ops::OpOutput::Transform { solid: s }) = nacre_ops::apply(
            &mut m,
            &nacre_ops::Operation::Transform {
                solid: s,
                isometry: turn,
            },
        ) else {
            panic!("{name}: the turn applies")
        };
        m.rebuild_adjacency();
        let reach = m.reachable();
        let mut want: Vec<[u64; 3]> = Vec::new();
        for e in (0..m.edge_count() as u32).filter_map(|i| m.edge_handle_at(i)) {
            if !reach.edges.contains(&e) {
                continue;
            }
            let nacre_geom::Curve::Line(l) = m.edge_curve(e) else {
                panic!("{name}: a box's edges are straight")
            };
            let cached = l.direction().as_array();
            let [v0, v1] = m.edge(e).vertices;
            let d = m.vertex_point(v1) - m.vertex_point(v0);
            if bits(d.normalize().expect("an edge has length").as_array()) != bits(cached) {
                endpoints_differ += 1;
            }
            want.push(bits(cached));
        }
        want.sort();

        let text = to_step_solids(&m, &[s], STAMP).expect("export");
        let (model, report) = read(text.as_bytes()).expect("re-read");
        assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
        let mut got: Vec<[u64; 3]> = model
            .line_arena
            .items
            .iter()
            .map(|line| {
                let VectorRef::Vector(v) = &line.dir else {
                    panic!("a plain VECTOR")
                };
                let vector = model.vector_arena.get(v.0);
                assert_eq!(
                    vector.magnitude.to_bits(),
                    1.0_f64.to_bits(),
                    "{name}: magnitude"
                );
                let DirectionRef::Direction(dr) = &vector.orientation else {
                    panic!("a plain DIRECTION")
                };
                let r = &model.direction_arena.get(dr.0).direction_ratios;
                bits([r[0], r[1], r[2]])
            })
            .collect();
        got.sort();
        assert_eq!(got.len(), 12, "{name}: a box's twelve LINEs");
        assert_eq!(got, want, "{name}: the file's directions");
    }
    assert!(
        endpoints_differ > 0,
        "every edge here has endpoints that already give its direction — the lock cannot tell the \
         roads apart"
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
    let text = to_step(&m, STAMP).expect("export cylinder");

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

    // The lateral (cylindrical) face is bounded by its two rims and no seam: two bounds of one
    // circle each.
    let lateral = faces
        .iter()
        .find(|f| matches!(f.surface().kind(), SurfaceKind::Cylindrical(_)))
        .unwrap();
    let bounds: Vec<usize> = lateral
        .bounds()
        .map(|b| b.oriented_edges().count())
        .collect();
    assert_eq!(bounds, vec![1, 1]);
}

/// **A leaf round-trips**: two quarter arcs of radius 5 between `(0, 0)` and `(5, 5)`, extruded.
/// Its two risers are edges between two cylindrical faces — lines, the ruling the two arcs' walls
/// share — and the reader takes the file whole: one solid, two caps, two cylindrical walls, every
/// edge a line or a circle.
#[test]
fn a_leaf_with_two_cylinders_on_one_edge_round_trips() {
    let r = nacre_exact::Rat::from_int;
    let arc = |c: [i128; 2]| nacre_ops::Edge2d::Arc {
        center: c.map(r),
        r2: r(25),
        ccw: false,
    };
    let ring = nacre_ops::Ring2d::new(
        vec![[r(0), r(0)], [r(5), r(5)]],
        vec![arc([5, 0]), arc([0, 5])],
    )
    .expect("a leaf");
    let profile = nacre_ops::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let mut m = Model::new();
    let frame = nacre_ops::SketchFrame::world(&m, nacre_exact::Axis::Z);
    nacre_ops::apply(
        &mut m,
        &nacre_ops::Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    )
    .expect("the leaf extrudes");
    let text = to_step(&m, STAMP).expect("export the leaf");

    let (model, report) = read(text.as_bytes()).expect("re-read");
    assert!(report.dropped.is_empty(), "drops: {:?}", report.dropped);
    let scene = model.scene();
    let solids: Vec<_> = scene.all_solids().collect();
    assert_eq!(solids.len(), 1);
    let faces: Vec<_> = solids[0].faces().collect();
    let (mut cylindrical, mut planes, mut lines) = (0, 0, 0);
    for face in &faces {
        match face.surface().kind() {
            SurfaceKind::Cylindrical(_) => cylindrical += 1,
            SurfaceKind::Plane(_) => planes += 1,
            other => panic!("unexpected surface kind: {other:?}"),
        }
        for bound in face.bounds() {
            for (edge, _forward) in bound.oriented_edges() {
                match edge.curve().kind() {
                    CurveKind::Line(_) => lines += 1,
                    CurveKind::Circle(_) => {}
                    other => panic!("unexpected curve kind: {other:?}"),
                }
            }
        }
    }
    assert_eq!((cylindrical, planes), (2, 2));
    // Each riser bounds both walls.
    assert_eq!(lines, 4, "two risers, each used by two faces");
}

/// An output that records the length of every write it receives.
struct Recording {
    bytes: Vec<u8>,
    writes: Vec<usize>,
}

impl std::io::Write for Recording {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.bytes.extend_from_slice(buf);
        self.writes.push(buf.len());
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// **`write_step` streams: the file reaches the output in pieces, not as one text at the end.**
/// That is what it is for — a caller writing to a file never holds the whole file. Thirty boxes make
/// a file of 185 KB, past the writer's 64 KiB chunk twice, so it arrives in at least three writes
/// (two chunks and the close); an adapter that gathered the text and handed it over once would show
/// one. Joined, the pieces are the file [`to_step`] returns.
#[test]
fn write_step_streams_the_file_in_pieces() {
    let mut m = Model::new();
    for i in 0..30 {
        let x = f64::from(2 * i);
        nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([x, 0.0, 0.0]),
            Point3::from_array([x + 1.0, 1.0, 1.0]),
        );
    }
    let out = Recording {
        bytes: Vec::new(),
        writes: Vec::new(),
    };
    let out = write_step(out, &m, m.live_solids(), STAMP).expect("export");
    let text = to_step(&m, STAMP).expect("export as text");

    assert!(
        text.len() > 2 * 64 * 1024,
        "the file is {} bytes",
        text.len()
    );
    assert!(out.writes.len() >= 3, "writes: {:?}", out.writes);
    let largest = out.writes.iter().max().copied().unwrap_or(0);
    assert!(largest < text.len(), "one write of the whole file");
    assert_eq!(out.bytes, text.as_bytes());
}

/// An output that refuses every write.
#[derive(Debug)]
struct Refusing;

impl std::io::Write for Refusing {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("refused"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// **A failing output is `StepError::Io`, carrying the output's own error** — not the writer's
/// message flattened into `Backend`, so the caller can ask what went wrong (a full disk, a closed
/// pipe). A small file reaches the output only when it is finished, so that is where it fails.
#[test]
fn a_failing_output_is_an_io_error() {
    let m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let err = write_step(Refusing, &m, m.live_solids(), STAMP).expect_err("the output refuses");
    match err {
        StepError::Io(e) => assert_eq!(e.to_string(), "refused"),
        other => panic!("not an io error: {other:?}"),
    }
}

/// **Parts are named `Part 1`, `Part 2`, … in the order given, and the header names the kernel
/// and its version.** A name a viewer shows as it is: with no name, FreeCAD shows `Shapes`,
/// step-loupe the definition's id and OCCT nothing, and several parts cannot be told apart. The
/// kernel holds no names, so a single solid written alone is always `Part 1`.
#[test]
fn parts_are_numbered_and_the_header_names_the_kernel() {
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

    let both = to_step(&m, STAMP).expect("export both");
    let first = both
        .find("PRODUCT('Part 1','Part 1'")
        .expect("the first part's name");
    let second = both
        .find("PRODUCT('Part 2','Part 2'")
        .expect("the second part's name");
    assert!(first < second);
    let system = format!("'nacre {}'", env!("CARGO_PKG_VERSION"));
    assert!(both.contains(&system), "the header does not name {system}");

    for h in [a, b] {
        let one = to_step_solids(&m, &[h], STAMP).expect("export one");
        assert!(one.contains("PRODUCT('Part 1','Part 1'"));
        assert!(!one.contains("'Part 2'"));
    }
}
