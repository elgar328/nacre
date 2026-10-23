//! **Every plane faces the way its truth says** — the sense lock over the populations where a
//! sense rule is easiest to get wrong, gathered on purpose: reflections carried into points and
//! recorded as nodes, turns, frames, statements of handles in both orders, and an offset under a
//! reflection (on a plain block and on one whose `Through` face makes every motion a node), and a
//! nameless statement still, moved and mirrored. The census holds the same lock over its corpus;
//! its reflections are few and its nameless planes none.

use super::*;
use nacre_exact::{Angle, Isometry, Rat, Rotation};
use nacre_math::Vector3;

fn world_box(m: &mut Model) -> Handle<Solid> {
    let frame = crate::SketchFrame::world(m, Axis::Z);
    match crate::apply(
        m,
        &crate::Operation::Extrude {
            frame,
            profile: square(),
            dist: 1.0,
        },
    ) {
        Ok(crate::OpOutput::Extrude { solid, .. }) => solid,
        other => panic!("a box: {other:?}"),
    }
}

fn moved(m: &mut Model, solid: Handle<Solid>, isometry: Isometry) -> Handle<Solid> {
    match crate::apply(m, &crate::Operation::Transform { solid, isometry }) {
        Ok(crate::OpOutput::Transform { solid }) => solid,
        other => panic!("a transform: {other:?}"),
    }
}

fn mirrored(m: &mut Model, solid: Handle<Solid>) -> Handle<Solid> {
    match crate::apply(
        m,
        &crate::Operation::Mirror {
            solid,
            axis: Axis::X,
            offset: Rat::from_int(0),
        },
    ) {
        Ok(crate::OpOutput::Mirror { solid }) => solid,
        other => panic!("a mirror: {other:?}"),
    }
}

fn turned_30(m: &mut Model, solid: Handle<Solid>) -> Handle<Solid> {
    moved(
        m,
        solid,
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(30)).expect("an angle"),
        }),
    )
}

/// A prism on `frame` — how a datum plane comes to be a face's surface, which is what the
/// lock walks.
fn prism_on(m: &mut Model, frame: crate::SketchFrame) -> Handle<Solid> {
    match crate::apply(
        m,
        &crate::Operation::Extrude {
            frame,
            profile: square(),
            dist: 0.5,
        },
    ) {
        Ok(crate::OpOutput::Extrude { solid, .. }) => solid,
        other => panic!("a prism on a datum: {other:?}"),
    }
}

fn datum(m: &mut Model, def: crate::DatumDef) -> crate::SketchFrame {
    match crate::apply(m, &crate::Operation::DatumPlane { def }) {
        Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("a datum: {other:?}"),
    }
}

/// Three consecutive corners of a solid's top cap.
fn top_corners(m: &Model, solid: Handle<Solid>) -> [Handle<nacre_topo::Vertex>; 3] {
    let top = m.shell(m.solid(solid).outer).faces[1];
    let he = &m.face(top).outer.half_edges;
    [0, 1, 2].map(|k| crate::he_start(m, he[k]))
}

/// The top cap's plane of `solid`, found by where its cache stands.
fn top_plane(m: &Model, solid: Handle<Solid>, z: f64) -> Handle<nacre_topo::Surface> {
    m.shell(m.solid(solid).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|&s| {
            matches!(m.surface_cache(s).plane(), Some(p)
                if (p.normal().as_array()[2].abs() - 1.0).abs() < 1e-12
                    && (p.origin().as_array()[2] - z).abs() < 1e-12)
        })
        .expect("a cap at that height")
}

#[test]
fn every_plane_faces_the_way_its_truth_says() {
    let mut cases: Vec<(&str, Model)> = Vec::new();

    // Turned: a rotation recorded as a node.
    let mut m = Model::new();
    let b = world_box(&mut m);
    turned_30(&mut m, b);
    cases.push(("turned", m));

    // Reflected, carried into the points: a fresh box's mirror is exact.
    let mut m = Model::new();
    let b = world_box(&mut m);
    mirrored(&mut m, b);
    cases.push(("mirror carried", m));

    // Reflected, recorded: a rounding translation makes a history, and the mirror follows it.
    let mut m = Model::new();
    let b = world_box(&mut m);
    let b = moved(
        &mut m,
        b,
        Isometry::translation([
            Rat::from_int(0),
            Rat::from_int(0),
            Rat::from_decimal(0.1).expect("0.1"),
        ]),
    );
    let b = mirrored(&mut m, b);
    // …and an offset of its top cap, which rides that reflection.
    let top = top_plane(&m, b, 1.1);
    let off = datum(
        &mut m,
        crate::DatumDef::Offset {
            frame: crate::SketchFrame::canonical(top),
            dist: 1.0,
        },
    );
    prism_on(&mut m, off);
    cases.push(("mirror recorded + offset", m));

    // The same on a block whose every motion is a node (its `Through` face): the mirror is
    // recorded behind the lift whatever the lift is, and the offset rides it.
    let mut m = Model::new();
    let b = super::fixtures::through_mirrored_block(&mut m, 1.0, Rat::new(1, 10).expect("1/10"));
    let top = top_plane(&m, b, 1.1);
    let off = datum(
        &mut m,
        crate::DatumDef::Offset {
            frame: crate::SketchFrame::canonical(top),
            dist: 1.0,
        },
    );
    prism_on(&mut m, off);
    cases.push(("through: mirror recorded + offset", m));

    // A frame: a prism on a tilted plane rides a frame node; its mirror adds a reflection.
    let mut m = Model::new();
    let tilted = crate::SketchPlane::from_origin_normal(
        Point3::origin(),
        Vector3::from_array([1.0, 2.0, 2.0]),
    )
    .expect("a tilted plane");
    let frame = datum_frame(&mut m, tilted);
    let p = prism_on(&mut m, frame);
    mirrored(&mut m, p);
    cases.push(("frame + mirror", m));

    // Statements of handles, in an even and an odd order of the same three corners.
    for (name, order) in [("through even", [0, 1, 2]), ("through odd", [1, 0, 2])] {
        let mut m = Model::new();
        let b = world_box(&mut m);
        let c = top_corners(&m, b);
        let f = datum(
            &mut m,
            crate::DatumDef::ThroughVertices(order.map(|k| c[k])),
        );
        prism_on(&mut m, f);
        cases.push((name, m));
    }

    // Statements of handles under a turn: the corners share the rotated box's chain.
    let mut m = Model::new();
    let b = world_box(&mut m);
    let b = turned_30(&mut m, b);
    let c = top_corners(&m, b);
    let f = datum(&mut m, crate::DatumDef::ThroughVertices([c[1], c[0], c[2]]));
    prism_on(&mut m, f);
    cases.push(("through turned", m));

    // A nameless statement — vertices exact in their own frames and not in each other's — and the
    // same plane moved and mirrored, which leaves the push door nothing to turn its cache by.
    let nameless = |m: &mut Model| -> Handle<Solid> {
        let still = world_box(m);
        let spun = world_box(m);
        let spun = moved(
            m,
            spun,
            Isometry::translation([Rat::from_int(4), Rat::from_int(0), Rat::from_int(0)]),
        );
        let spun = moved(
            m,
            spun,
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(37)).expect("an angle"),
            }),
        );
        let top = top_corners(m, still);
        let base = m.shell(m.solid(spun).outer).faces[0];
        let low = crate::he_start(m, m.face(base).outer.half_edges[0]);
        let f = datum(m, crate::DatumDef::ThroughVertices([top[0], top[1], low]));
        prism_on(m, f)
    };
    let mut m = Model::new();
    nameless(&mut m);
    cases.push(("nameless", m));
    let mut m = Model::new();
    let p = nameless(&mut m);
    let p = moved(
        &mut m,
        p,
        Isometry::translation([
            Rat::new(1, 3).expect("1/3"),
            Rat::from_int(2),
            Rat::from_int(0),
        ]),
    );
    moved(
        &mut m,
        p,
        Isometry::rotation(Rotation {
            axis: Axis::X,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(20)).expect("an angle"),
        }),
    );
    cases.push(("nameless moved", m));
    let mut m = Model::new();
    let p = nameless(&mut m);
    mirrored(&mut m, p);
    cases.push(("nameless mirrored", m));
    // A nameless statement with a straddling vertex — a corner of two boxes turned about different
    // axes and fused — whose witness is its judged frame's probes, not its own points: still, moved
    // and mirrored.
    let straddling = |m: &mut Model| -> Handle<Solid> {
        let a = crate::fixtures::cuboid(
            m,
            nacre_math::Point3::from_array([0.0, 0.0, 0.0]),
            nacre_math::Point3::from_array([2.0, 2.0, 2.0]),
        );
        let b = crate::fixtures::cuboid(
            m,
            nacre_math::Point3::from_array([1.0, 1.0, 1.0]),
            nacre_math::Point3::from_array([3.0, 3.0, 3.0]),
        );
        m.rebuild_adjacency();
        let turned = |deg: i128, axis| {
            Isometry::rotation(Rotation {
                axis,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
            })
        };
        let a = moved(m, a, turned(37, Axis::Z));
        let b = moved(m, b, turned(20, Axis::X));
        let Ok(crate::OpOutput::Boolean { solids }) = crate::apply(
            m,
            &crate::Operation::Boolean {
                kind: crate::BoolKind::Fuse,
                a,
                b,
            },
        ) else {
            panic!("the fuse")
        };
        m.rebuild_adjacency();
        let straddler = m
            .reachable_from(&solids)
            .vertices
            .into_iter()
            .filter(|&v| {
                matches!(m.vertex(v), nacre_topo::Vertex::ThreePlane(t)
                    if t.iter().all(|h| m.surface_name(*h).is_some()))
                    && m.vertex_meet(v).is_none()
            })
            .min_by_key(|v| v.index())
            .expect("a fused corner straddles");
        let still = world_box(m);
        let still = moved(
            m,
            still,
            Isometry::translation([Rat::from_int(10), Rat::from_int(0), Rat::from_int(0)]),
        );
        let top = top_corners(m, still);
        let f = datum(
            m,
            crate::DatumDef::ThroughVertices([straddler, top[0], top[1]]),
        );
        prism_on(m, f)
    };
    let mut m = Model::new();
    straddling(&mut m);
    cases.push(("nameless straddling", m));
    let mut m = Model::new();
    let p = straddling(&mut m);
    let p = moved(
        &mut m,
        p,
        Isometry::translation([
            Rat::new(1, 3).expect("1/3"),
            Rat::from_int(2),
            Rat::from_int(0),
        ]),
    );
    moved(
        &mut m,
        p,
        Isometry::rotation(Rotation {
            axis: Axis::X,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(20)).expect("an angle"),
        }),
    );
    cases.push(("nameless straddling moved", m));
    let mut m = Model::new();
    let p = straddling(&mut m);
    mirrored(&mut m, p);
    cases.push(("nameless straddling mirrored", m));
    // The premise: each holds a live face on a nameless `Through` plane, the moved and mirrored ones
    // give it a motion of its own, and the straddling ones have a defining vertex in two frames.
    for (name, m) in &cases[cases.len() - 6..] {
        let moved_case = name.ends_with("moved") || name.ends_with("mirrored");
        let straddle_case = name.contains("straddling");
        let straddles =
            |vs: &[Handle<nacre_topo::Vertex>; 3]| vs.iter().any(|&v| m.vertex_meet(v).is_none());
        let found = m.live_solids().iter().any(|&s| {
            m.shell(m.solid(s).outer).faces.iter().any(|&f| {
                let h = m.face(f).surface;
                matches!(
                    m.surface(h),
                    nacre_topo::Surface::Plane {
                        points: nacre_topo::PlanePoints::Through(vs),
                        motion,
                        ..
                    } if motion.is_some() == moved_case && straddles(vs) == straddle_case
                ) && m.surface_name(h).is_none()
            })
        });
        assert!(
            found,
            "{name}: no live face on a nameless plane of that motion"
        );
        // And each such plane's cache is realized from its truth — the funnel's own answer asked
        // again, as the census holds every other plane to it.
        for &s in m.live_solids() {
            for &f in &m.shell(m.solid(s).outer).faces {
                let h = m.face(f).surface;
                let nameless = matches!(
                    m.surface(h),
                    nacre_topo::Surface::Plane {
                        points: nacre_topo::PlanePoints::Through(_),
                        ..
                    }
                ) && m.surface_name(h).is_none();
                if !nameless {
                    continue;
                }
                assert_eq!(
                    m.surface_cache_standing(h),
                    nacre_topo::CacheStanding::Realized,
                    "{name}: a nameless plane's cache is its truth's realization"
                );
                assert_eq!(
                    crate::realize_surface_cache(m, h).map(|(c, _)| c).as_ref(),
                    Some(&m.surface_cache(h).to_geom()),
                    "{name}: and the funnel's answer, asked again"
                );
            }
        }
    }

    let (mut mirrored_seen, mut turned_seen) = (0, 0);
    for (name, m) in &cases {
        let a = crate::audit_plane_senses(m);
        assert!(a.disagree.is_empty(), "{name}: opposed {:?}", a.disagree);
        assert_eq!(
            a.unmeasured, 0,
            "{name}: the lock could not look at every plane"
        );
        assert!(a.agree > 0, "{name}: nothing was judged");
        if name.starts_with("nameless") {
            assert!(a.nameless > 0, "{name}: no nameless plane was read");
        }
        if name.starts_with("through: mirror") {
            assert!(
                a.mirrored > 0,
                "{name}: no plane rode the recorded reflection"
            );
        }
        mirrored_seen += a.mirrored;
        turned_seen += a.turned;
    }
    // The populations the fixtures are gathered for are really there.
    assert!(mirrored_seen > 0, "no plane rode a reflection");
    assert!(turned_seen > 0, "no plane rode a turn or a frame");
}
