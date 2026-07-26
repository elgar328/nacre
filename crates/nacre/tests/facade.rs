//! Does the facade actually stand on its own?
//!
//! Every path here goes through `nacre::`. A missing re-export is not a subtle
//! degradation — this file stops compiling, which is the point: a list of `pub use`
//! lines proves nothing by itself.

use nacre::prelude::*;

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
}

fn extrude(m: &mut Model, profile: Profile2d, dist: f64) -> Handle<Solid> {
    match apply(
        m,
        &Operation::Extrude {
            plane: SketchPlane::world_xy(),
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
/// `scalar::Isometry`, so a facade that forgot to re-export `nacre-scalar` would leave
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

    let step = nacre::step::to_step(&model).unwrap();
    assert!(
        step.starts_with("ISO-10303-21;"),
        "{}",
        &step[..40.min(step.len())]
    );
}

/// The sketch front door and the face queries, also facade-only. `from_edges` takes
/// loose segments in any order; `face_props` is what lets a caller *name* a face by its
/// geometry rather than by an index.
#[test]
fn the_sketch_front_door_and_face_queries_are_reachable() {
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let profiles = from_edges(vec![
        Edge2d::line(p(0.0, 0.0), p(4.0, 0.0)),
        Edge2d::line(p(4.0, 0.0), p(4.0, 4.0)),
        Edge2d::line(p(4.0, 4.0), p(0.0, 4.0)),
        Edge2d::line(p(0.0, 4.0), p(0.0, 0.0)),
    ])
    .unwrap();
    assert_eq!(profiles.len(), 1);

    let mut model = Model::new();
    let solid = extrude(&mut model, profiles.into_iter().next().unwrap(), 1.0);
    model.rebuild_adjacency();

    // Pick the face looking straight up, the way a script names one.
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let top = *model
        .shells
        .get(model.solids.get(solid).outer)
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
    assert!((plane.origin[2] - 1.0).abs() < 1e-12, "{:?}", plane.origin);

    // A malformed sketch is refused by name, not silently built.
    let bowtie = from_edges(vec![
        Edge2d::line(p(0.0, 0.0), p(4.0, 4.0)),
        Edge2d::line(p(4.0, 4.0), p(4.0, 0.0)),
        Edge2d::line(p(4.0, 0.0), p(0.0, 4.0)),
        Edge2d::line(p(0.0, 4.0), p(0.0, 0.0)),
    ]);
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
    let cube = model.add_cuboid(Point3::origin(), Point3::from_array([2.0, 3.0, 4.0]));
    model.rebuild_adjacency();
    assert!((mass_props(&model, cube).unwrap().volume - 24.0).abs() < 1e-12);

    // Tessellation and OBJ live one module away; check they are wired at all.
    let mesh = tessellate(&model, &TessConfig::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
    let _: Tessellation = mesh;
    assert!(nacre::tess::to_obj(&model).unwrap().contains("v "));
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
