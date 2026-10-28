//! Does the facade actually stand on its own?
//!
//! Every path here goes through `nacre::`. A missing re-export is not a subtle
//! degradation — this file stops compiling, which is the point: a list of `pub use`
//! lines proves nothing by itself.

use nacre::prelude::*;
use nacre_exact::Axis;
use nacre_ops::SketchFrame;

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
    .unwrap()
}

fn extrude(m: &mut Model, profile: Profile2d, dist: f64) -> Handle<Solid> {
    match apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
            profile,
            dist,
        },
    ) {
        Ok(OpOutput::Extrude { solid, .. }) => solid,
        other => panic!("extrude: {other:?}"),
    }
}

/// Sketch → extrude → **transform** → boolean → validate → mass properties → STEP,
/// with `use nacre::prelude::*` as the only import.
///
/// The `Transform` step is not decoration. `Operation::Transform` carries a
/// `exact::Isometry`, so a facade that forgot to re-export `nacre-exact` would leave
/// a consumer unable to *name* the value this operation needs — the one leak of a lower
/// layer into `ops`' public API. Exercising it turns that claim into a compile-time
/// fact. `bounds` is here for the same reason: the newest query layer must be reachable.
#[test]
fn the_whole_pipeline_runs_through_the_facade_alone() {
    let mut model = Model::new();

    // Two unit-ish plates: one at the origin, one slid over to overlap it by half.
    let a = extrude(&mut model, square(0.0, 4.0), 2.0);
    let b = extrude(&mut model, square(0.0, 4.0), 2.0);
    let two = Rat::from_int(2);
    let zero = Rat::from_int(0);
    let b = match apply(
        &mut model,
        &Operation::Transform {
            solid: b,
            isometry: Isometry::translation([two, zero, zero]),
        },
    ) {
        Ok(OpOutput::Transform { solid }) => solid,
        other => panic!("transform: {other:?}"),
    };

    let joined = boolean(&mut model, BoolKind::Fuse, a, b).unwrap();
    assert_eq!(joined.len(), 1);
    model.rebuild_adjacency();

    let issues = validate(&model);
    assert!(issues.is_empty(), "{issues:?}");

    // The union spans x in [0, 6], y in [0, 4], z in [0, 2]: volume 48, and a bounding
    // box the shifted copy widened.
    let props = mass_props(&model, joined[0]).unwrap();
    assert!((props.volume - 48.0).abs() < 1e-9, "{}", props.volume);
    let (lo, hi) = bounds(&model, joined[0]).unwrap();
    assert!(
        (lo[0] - 0.0).abs() < 1e-12 && (hi[0] - 6.0).abs() < 1e-12,
        "{lo:?} {hi:?}"
    );
    assert!((centroid(&model, joined[0]).unwrap()[0] - 3.0).abs() < 1e-9);

    let step = nacre::step::to_step(&model, "2026-10-05T00:00:00Z").unwrap();
    assert!(
        step.starts_with("ISO-10303-21;"),
        "{}",
        &step[..40.min(step.len())]
    );
}

/// The sketch front door and the face queries, also facade-only. `from_rings` takes closed
/// rings as written; `face_props` is what lets a caller *name* a face by its geometry rather
/// than by an index.
#[test]
fn the_sketch_front_door_and_face_queries_are_reachable() {
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let profiles = from_rings(vec![vec![
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(0.0, 4.0),
    ]])
    .unwrap();
    assert_eq!(profiles.len(), 1);

    let mut model = Model::new();
    let solid = extrude(&mut model, profiles.into_iter().next().unwrap(), 1.0);
    model.rebuild_adjacency();

    // Pick the face looking straight up, the way a script names one.
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let top = *model
        .shell(model.solid(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            face_props(&model, f)
                .unwrap()
                .normal
                .is_some_and(|n| n.dot(up) > 0.99)
        })
        .expect("a prism has a top face");
    assert!((face_props(&model, top).unwrap().area - 16.0).abs() < 1e-12);

    // And the sketch plane a pad/pocket would place a profile in.
    let plane = face_plane(&model, top).unwrap();
    assert!(
        (plane.origin()[2] - 1.0).abs() < 1e-12,
        "{:?}",
        plane.origin()
    );

    // A malformed sketch is refused by name, not silently built.
    let bowtie = from_rings(vec![vec![
        p(0.0, 0.0),
        p(4.0, 4.0),
        p(4.0, 0.0),
        p(0.0, 4.0),
    ]]);
    assert!(matches!(
        bowtie,
        Err(SketchError::RingSelfIntersects { .. })
    ));
}

/// The topology layer's test-only primitive, reached through the facade. It is enabled
/// by this crate's `[dev-dependencies]` rather than its own `test-util` feature, so this
/// runs in a plain `cargo test` instead of only under an extra flag.
#[test]
fn the_test_util_primitive_is_forwarded() {
    let mut model = Model::new();
    let cube = nacre_ops::fixtures::cuboid(
        &mut model,
        Point3::origin(),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    model.rebuild_adjacency();
    assert!((mass_props(&model, cube).unwrap().volume - 24.0).abs() < 1e-12);

    // Tessellation and OBJ live one module away; check they are wired at all.
    let mesh = tessellate(&model, &TessConfig::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
    let _: Tessellation = mesh;
    assert!(mesh.to_obj(&model).contains("v "));
}

/// **A manifest invariant, guarded where it can actually be checked.**
///
/// A consumer turns rayon off with `nacre = { default-features = false }` — which only
/// works because *this* crate takes `nacre-ops` without its defaults and re-enables
/// them through its own `parallel` feature. Drop that one attribute and rayon becomes
/// unremovable: `cargo tree -p nacre --no-default-features` still lists it (measured).
///
/// The property lives in `Cargo.toml`, not in code, and the pre-commit hook does not
/// build with feature flags — so it is asserted against the manifest itself. Clumsy,
/// but a comment nobody re-runs is worse: this is the failure mode most likely to slip
/// back in, and it breaks the wasm consumer silently.
#[test]
fn the_ops_dependency_stays_default_free_so_consumers_can_drop_rayon() {
    let manifest = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml"))
        .expect("the facade's own manifest");
    let line = manifest
        .lines()
        .find(|l| l.starts_with("nacre-ops ="))
        .expect("nacre-ops must be a dependency");
    assert!(
        line.contains("default-features = false"),
        "nacre-ops must be taken without default features, or `nacre \
         --no-default-features` still links rayon and a wasm consumer is stuck: {line}"
    );
}

/// **Can a caller outside the kernel actually name three vertices?**
///
/// The datum that closes the coordinate gap is only useful if the vertices can be pointed at from
/// the facade, and today there is no query for that — the caller walks
/// `solids → shells → faces → loops → half_edges` and calls `he_start`. Every step of that is
/// public, so it *works*; this test is what says so, and it is also the evidence for whether a
/// derived query is worth adding. Read it: if the walk below looks like something a script author
/// should never write, that is the argument, made of code rather than of guesswork.
#[test]
fn a_caller_can_name_vertices_and_build_a_datum_through_them() {
    use nacre::ops::{DatumDef, OpOutput, Operation, apply};
    use nacre::prelude::*;

    let mut m = nacre::topo::Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 4.0]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, 2.0]),
        Point3::from_array([3.0, 3.0, 5.0]),
    );
    nacre::ops::boolean(&mut m, nacre::ops::BoolKind::Cut, a, b).expect("a pocketed box");
    m.rebuild_adjacency();

    // The walk. This is the whole vocabulary a caller has for "give me a corner".
    let mut corners = Vec::new();
    for &s in m.live_solids() {
        let sol = m.solid(s);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let f = m.face(fh);
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for &he in &lp.half_edges {
                        let v = m.he_start(he);
                        if !corners.contains(&v) {
                            corners.push(v);
                        }
                    }
                }
            }
        }
    }
    assert!(corners.len() >= 3, "a cut box has corners to name");

    // Pick a triple that names a plane, then build the datum and raise something on it.
    let mut triple = None;
    'search: for i in 0..corners.len() {
        for j in (i + 1)..corners.len() {
            for k in (j + 1)..corners.len() {
                let t = [corners[i], corners[j], corners[k]];
                let mut sorted = t;
                sorted.sort_by_key(|v| v.index());
                if m.plane_name_through(sorted).is_some() {
                    triple = Some(t);
                    break 'search;
                }
            }
        }
    }
    let triple = triple.expect("some triple of a cut box's corners names a plane");

    let OpOutput::DatumPlane { frame, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(triple),
        },
    )
    .expect("a datum through named vertices, from the facade alone") else {
        unreachable!()
    };
    let out = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: nacre::ops::Profile2d::polygon(vec![
                Point2::from_array([-0.3, -0.3]),
                Point2::from_array([0.3, -0.3]),
                Point2::from_array([0.3, 0.3]),
                Point2::from_array([-0.3, 0.3]),
            ])
            .unwrap(),
            dist: 0.5,
        },
    );
    assert!(
        matches!(out, Ok(OpOutput::Extrude { .. })),
        "the datum a caller named must be sketchable: {out:?}"
    );
}

/// Every public error type of the facade is an `Error` a `Box<dyn Error>` — or an
/// `anyhow::Error` — takes: `Send + Sync + 'static`. A listed type that loses a bound stops this
/// file compiling; a new error type joins the list by hand.
#[test]
fn every_error_type_is_a_std_error() {
    fn is_error<T: std::error::Error + Send + Sync + 'static>() {}
    is_error::<nacre::ops::BoolError>();
    is_error::<nacre::ops::OpError>();
    is_error::<nacre::ops::SketchError>();
    is_error::<nacre::ops::RealizeError>();
    is_error::<nacre::ops::CacheDecline>();
    is_error::<nacre::step::StepError>();
    is_error::<nacre::tess::TessError>();
    is_error::<nacre::props::PropsError>();
    is_error::<nacre::topo::EdgeDecline>();
    is_error::<nacre::geom::mixed::Undecidable>();
}

fn disk(m: &mut Model) -> Handle<Solid> {
    let circle = Ring2d::circle(Point2::from_array([0.0, 0.0]), 1.0).unwrap();
    extrude(m, from_paths(vec![circle]).unwrap().remove(0), 1.0)
}

/// What the doors actually return, read as text: a refused boolean says its class and reason,
/// `OpError::Boolean` says the same line and adds no source, and a failed write hands its
/// `io::Error` on as the source rather than repeating it.
#[test]
fn errors_from_the_doors_say_one_line_and_chain_their_cause() {
    use std::error::Error;

    // Two identical cylinders share one lateral surface — a boolean the kernel refuses by name.
    let mut m = Model::new();
    let (a, b) = (disk(&mut m), disk(&mut m));
    let refused = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();
    let BoolError::Rejected { reason, .. } = refused else {
        panic!("{refused:?}")
    };
    assert_eq!(reason, RejectReason::CylinderPairContact);
    let line = refused.to_string();
    assert!(
        line.starts_with(&format!(
            "boolean rejected, not supported: {}",
            reason.as_str()
        )),
        "{line}"
    );

    let (c, d) = (disk(&mut m), disk(&mut m));
    let wrapped = apply(
        &mut m,
        &Operation::Boolean {
            kind: BoolKind::Fuse,
            a: c,
            b: d,
        },
    )
    .unwrap_err();
    assert!(matches!(wrapped, OpError::Boolean(_)), "{wrapped:?}");
    assert_eq!(wrapped.to_string(), line);
    assert!(wrapped.source().is_none());

    struct Refuses;
    impl std::io::Write for Refuses {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("disk full"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let solids: Vec<_> = m.live_solids().to_vec();
    let Err(failed) = nacre::step::write_step(Refuses, &m, &solids, "") else {
        panic!("a writer that refuses every byte cannot take a file")
    };
    assert!(
        matches!(failed, nacre::step::StepError::Io(_)),
        "{failed:?}"
    );
    assert_eq!(failed.to_string(), "writing the STEP output failed");
    let cause = failed.source().expect("the io::Error is the source");
    assert_eq!(
        cause
            .downcast_ref::<std::io::Error>()
            .map(|e| e.to_string()),
        Some("disk full".to_string())
    );

    let no_radius = Ring2d::circle(Point2::from_array([0.0, 0.0]), 0.0).unwrap_err();
    assert_eq!(
        no_radius.to_string(),
        "the radius 0 at (0, 0) is not positive"
    );
}

/// `?` carries every door's error into one `Box<dyn Error>`, and the caller can still take it
/// back out by type.
#[test]
fn question_mark_carries_kernel_errors_into_a_box() {
    fn bored(m: &mut Model) -> Result<f64, Box<dyn std::error::Error>> {
        let circle = Ring2d::circle(Point2::from_array([2.0, 2.0]), 1.0)?;
        let plate = extrude(m, square(0.0, 4.0), 1.0);
        let bore = extrude(m, from_paths(vec![circle])?.remove(0), 1.0);
        let result = boolean(m, BoolKind::Cut, plate, bore)?;
        Ok(mass_props(m, result[0])?.volume)
    }
    let volume = bored(&mut Model::new()).unwrap();
    assert!(
        (volume - (16.0 - std::f64::consts::PI)).abs() < 1e-9,
        "{volume}"
    );

    fn fused(m: &mut Model) -> Result<(), Box<dyn std::error::Error>> {
        let (a, b) = (disk(m), disk(m));
        boolean(m, BoolKind::Fuse, a, b)?;
        Ok(())
    }
    let boxed = fused(&mut Model::new()).unwrap_err();
    assert!(
        matches!(
            boxed.downcast_ref::<BoolError>(),
            Some(BoolError::Rejected { .. })
        ),
        "{boxed}"
    );
}

/// The std convention for every variant's line: something to say, on one line, starting
/// lowercase, with no closing full stop. The `match` in each `Display` is what makes a new
/// variant choose its words; this list only reads how they are written.
#[test]
fn every_error_line_follows_the_std_convention() {
    use nacre::ops::{CacheDecline, LogCell, ProfileRing, RealizeError, RejectWhere};
    use nacre::props::PropsError;
    use nacre::step::StepError;
    use nacre::tess::TessError;
    use nacre::topo::EdgeDecline;

    let p = Point3::from_array([1.0, 2.0, 3.0]);
    let rejected = |at| BoolError::Rejected {
        reason: RejectReason::CylinderPairContact,
        at,
    };
    let (outer, hole) = (ProfileRing::Outer, ProfileRing::Hole(0));
    let xy = [1.0, 2.0];
    let lines: Vec<String> = [
        rejected(None).to_string(),
        rejected(Some(RejectWhere::Point(p))).to_string(),
        rejected(Some(RejectWhere::Segment([p, p]))).to_string(),
        BoolError::InputNotLive.to_string(),
    ]
    .into_iter()
    .chain(
        [
            OpError::DegenerateProfile,
            OpError::ProfileOutsideDecimalWindow {
                ring: hole,
                point: 2,
            },
            OpError::SelfIntersectingProfile {
                ring: outer,
                edges: (0, 2),
            },
            OpError::ZeroLengthProfileEdge {
                ring: outer,
                edge: 1,
            },
            OpError::ProfileRingsMeet { a: outer, b: hole },
            OpError::HoleNotInsideOuter { hole: 0 },
            OpError::NestedHole { outer: 0, inner: 1 },
            OpError::ProfileUndecidable,
            OpError::ArcSweepNotQuarterTurn,
            OpError::ZeroDistance,
            OpError::DistOutsideDecimalWindow,
            OpError::PlaneWithoutExactForm,
            OpError::FrameOutsideDecimalWindow,
            OpError::OriginNotOnPlane,
            OpError::RefDirParallelToNormal,
            OpError::DegenerateGeometry,
            OpError::NonPlanarFace,
            OpError::FaceNotInLiveSolid,
            OpError::Boolean(BoolError::InputNotLive),
            OpError::SolidNotLive,
            OpError::OriginNotOnSolid,
            OpError::MirrorNotPlanar,
            OpError::ZeroOffset,
            OpError::DuplicateVertex,
            OpError::CollinearVertices,
            OpError::VertexNotThreePlane,
            OpError::VerticesInMixedFrames,
            OpError::ThroughFrameUndecided,
            OpError::LogHandleOutOfRange {
                cell: LogCell::Surface,
                index: 7,
            },
        ]
        .iter()
        .map(ToString::to_string),
    )
    .chain(
        [
            SketchError::DegenerateRing { ring: 0 },
            SketchError::RingsMeet { a: 0, b: 1 },
            SketchError::UnevenRing {
                vertices: 3,
                edges: 2,
            },
            SketchError::ZeroLengthEdge { edge: 1 },
            SketchError::RingSelfIntersects {
                ring: 0,
                at: [xy, xy],
            },
            SketchError::OutsideDecimalWindow { at: xy },
            SketchError::ArcEndOffCircle {
                center: xy,
                start: xy,
                end: xy,
            },
            SketchError::ArcRadiusMismatch {
                center: xy,
                start: xy,
                stated: 1.5,
            },
            SketchError::ZeroLengthArc { at: xy },
            SketchError::ArcTurnsOutOfRange { turns: 4 },
            SketchError::NonPositiveRadius {
                center: xy,
                radius: -1.0,
            },
            SketchError::Undecidable,
        ]
        .iter()
        .map(ToString::to_string),
    )
    .chain(
        [
            RealizeError::NoMeet,
            RealizeError::NoPoint,
            RealizeError::NoMotionChain,
            RealizeError::NoCurvedPoint,
            RealizeError::Undecided,
            RealizeError::Unrepresentable,
        ]
        .iter()
        .map(ToString::to_string),
    )
    .chain([
        CacheDecline::CostCap.to_string(),
        CacheDecline::Cannot(RealizeError::NoPoint).to_string(),
        StepError::Backend("a header string is too long".to_string()).to_string(),
        StepError::Io(std::io::Error::other("disk full")).to_string(),
    ])
    .chain(
        [
            TessError::DegenerateRing,
            TessError::SelfTouchingBoundary,
            TessError::HoleWinding,
            TessError::OverBudget,
        ]
        .iter()
        .map(ToString::to_string),
    )
    .chain(
        [
            PropsError::UnsupportedBoundary,
            PropsError::UnsupportedInnerLoop,
        ]
        .iter()
        .map(ToString::to_string),
    )
    .chain(
        [
            EdgeDecline::Coincident,
            EdgeDecline::Oblique,
            EdgeDecline::Unstated,
            EdgeDecline::TwoCylinders,
            EdgeDecline::Degenerate,
        ]
        .iter()
        .map(ToString::to_string),
    )
    .chain([nacre::geom::mixed::Undecidable.to_string()])
    .collect();

    let broken: Vec<&String> = lines
        .iter()
        .filter(|l| {
            l.is_empty()
                || l.contains('\n')
                || l.ends_with('.')
                || l.starts_with(|c: char| c.is_uppercase())
        })
        .collect();
    assert!(broken.is_empty(), "{broken:#?}");
}
