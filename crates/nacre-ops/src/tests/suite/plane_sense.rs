//! **Every plane faces the way its truth says** — the sense lock over the populations where a
//! sense rule is easiest to get wrong, gathered on purpose: reflections carried into points and
//! recorded as nodes, turns, frames, statements of handles in both orders, and an offset under a
//! reflection (on a plain block and on one whose `Through` face makes every motion a node). The
//! census holds the same lock over its corpus; its reflections are few.

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
            matches!(m.surface_cache(s), nacre_geom::Surface::Plane(p)
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

    let (mut mirrored_seen, mut turned_seen) = (0, 0);
    for (name, m) in &cases {
        let a = crate::audit_plane_senses(m);
        assert!(a.disagree.is_empty(), "{name}: opposed {:?}", a.disagree);
        assert_eq!(
            a.unmeasured, 0,
            "{name}: the lock could not look at every plane"
        );
        assert!(a.agree > 0, "{name}: nothing was judged");
        if name.starts_with("through: mirror") {
            assert!(
                a.mirrored > 0,
                "{name}: no plane rode the recorded reflection"
            );
        }
        mirrored_seen += a.mirrored;
        turned_seen += a.turned;
    }
    // The populations the fixtures were gathered for are really there.
    assert!(mirrored_seen > 0, "no plane rode a reflection");
    assert!(turned_seen > 0, "no plane rode a turn or a frame");
}
