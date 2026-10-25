//! Corners where a plane's line crosses a cylinder: `Vertex::Pierce`, named by its place about the
//! cylinder and realized along its ruling.

use super::*;
use nacre_exact::Rat;

/// `p` turned `deg` degrees about the world `axis` through the origin, in `f64` — the oracle a
/// realized point is held to, never the road that realizes it.
fn turned_f64(p: [f64; 3], axis: Axis, deg: f64) -> [f64; 3] {
    let (s, c) = deg.to_radians().sin_cos();
    let [x, y, z] = p;
    match axis {
        Axis::X => [x, c * y - s * z, s * y + c * z],
        Axis::Y => [c * x + s * z, y, -s * x + c * z],
        Axis::Z => [c * x - s * y, s * x + c * y, z],
    }
}

/// The pierce corners a solid reaches, each once.
fn pierce_corners(m: &Model, s: Handle<Solid>) -> Vec<Handle<Vertex>> {
    let mut out = Vec::new();
    for &f in &m.shell(m.solid(s).outer).faces {
        let face = m.face(f);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                for v in m.edge(he.edge).vertices {
                    if matches!(m.vertex(v), Vertex::Pierce { .. }) && !out.contains(&v) {
                        out.push(v);
                    }
                }
            }
        }
    }
    out
}

/// **A turned pierce corner is realized where the turn took it.** Each shape is built twice: once
/// unmoved, where its pierce corners stand on world statements, and once turned off the quarters
/// (and on, after such a turn). Every pierce corner of the turned solid realizes, nearest-`f64`,
/// within `1e-9` of an unmoved corner turned in `f64` — the other crossing of the same line is
/// over a unit away. The shapes are the ones whose corners had no road before: a windowed boss,
/// whose walls and cylinder share the turn's node or whose wall the turn fixes, and a lateral cut
/// at both rims turned twice, whose wall the first turn fixes and the second moves on a node of
/// its own — the tail of the cylinder's chain.
#[test]
fn a_turned_pierce_corner_is_realized_where_the_turn_took_it() {
    type Build = fn(&mut Model) -> Handle<Solid>;
    type Shape = (&'static str, Build, &'static [(Axis, i128)]);
    let shapes: [Shape; 7] = [
        (
            "window",
            |m| crate::fixtures::windowed_boss(m, 0.6, -1.0),
            &[(Axis::Y, 37)],
        ),
        (
            "window twice",
            |m| crate::fixtures::windowed_boss(m, 0.6, -1.0),
            &[(Axis::Y, 37), (Axis::X, 23)],
        ),
        (
            "window, then a quarter about z",
            |m| crate::fixtures::windowed_boss(m, 0.6, -1.0),
            &[(Axis::Y, 37), (Axis::Z, 90)],
        ),
        (
            "window, then a quarter about x",
            |m| crate::fixtures::windowed_boss(m, 0.6, -1.0),
            &[(Axis::Y, 37), (Axis::X, 90)],
        ),
        (
            "both rims fused",
            |m| crate::fixtures::cut_at_both_rims(m, BoolKind::Fuse, -2.0),
            &[(Axis::X, 37), (Axis::Y, 23)],
        ),
        (
            "both rims cut",
            |m| crate::fixtures::cut_at_both_rims(m, BoolKind::Cut, -2.0),
            &[(Axis::X, 37), (Axis::Y, 23)],
        ),
        (
            "both rims cut, once",
            |m| crate::fixtures::cut_at_both_rims(m, BoolKind::Cut, -2.0),
            &[(Axis::X, 37)],
        ),
    ];
    for (name, build, turns) in shapes {
        let mut m0 = Model::new();
        let s0 = build(&mut m0);
        m0.rebuild_adjacency();
        let want: Vec<[f64; 3]> = pierce_corners(&m0, s0)
            .into_iter()
            .map(|v| {
                turns
                    .iter()
                    .fold(m0.vertex_point(v).as_array(), |p, &(axis, deg)| {
                        turned_f64(p, axis, deg as f64)
                    })
            })
            .collect();
        assert!(!want.is_empty(), "{name}: the shape has pierce corners");
        let mut m = Model::new();
        let mut s = build(&mut m);
        for &(axis, deg) in turns {
            s = transform(&mut m, s, &rot_iso(axis, deg)).expect("a turn");
        }
        m.rebuild_adjacency();
        let got = pierce_corners(&m, s);
        assert_eq!(got.len(), want.len(), "{name}");
        for v in got {
            let p = crate::realize_vertex(&m, v, crate::Precision::NearestF64)
                .unwrap_or_else(|e| panic!("{name}: {:?} {e:?}", m.vertex(v)))
                .to_f64()
                .map(|(x, _)| x)
                .expect("a finite point");
            assert!(
                want.iter()
                    .any(|w| (0..3).all(|k| (w[k] - p[k]).abs() < 1e-9)),
                "{name}: {:?} realizes at {p:?}, where no corner was turned to",
                m.vertex(v)
            );
            assert!(
                matches!(m.vertex_cache(v), nacre_topo::PointCache::Realized { .. }),
                "{name}: {:?}",
                m.vertex_cache(v)
            );
        }
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{name}: {:?}",
            nacre_validate::validate(&m)
        );
    }
}

/// A slot `2r` wide between `y0` and `y1` on the frame's `y` axis: two straight walls `x = ±r`
/// tangent to two half cylinders. Its corners sit at angle `0` and `π` about the arcs' `ref_dir`
/// (the frame's `x̂`) — tangent corners on the seam ray.
fn upright_slot(r: f64, y0: f64, y1: f64) -> Profile2d {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    stated(vec![
        line(p2(r, y0), p2(r, y1)),
        arc_turns(p2(0.0, y1), p2(r, y1), 2),
        line(p2(-r, y1), p2(-r, y0)),
        arc_turns(p2(0.0, y0), p2(-r, y0), 2),
    ])
    .expect("a slot")
    .remove(0)
}

fn prism_on(m: &mut Model, frame: SketchFrame, profile: Profile2d, h: f64) -> Handle<Solid> {
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("a prism") else {
        unreachable!()
    };
    solid
}

/// **A turned slot is realized, and the export door has nothing left to raise.** An upright slot —
/// its tangent corners on the seam ray, rank 0 — turned 30° about its own axis (the caps keep their
/// statements, the walls and cylinders take the node) and about `x` (everything moves): every
/// pierce corner realizes, `refine_caches` leaves no vertex unrealized, and the solid is valid.
#[test]
fn a_turned_slot_is_realized_and_leaves_nothing_for_the_export() {
    for axis in [Axis::Z, Axis::X] {
        let mut m = Model::new();
        let frame = SketchFrame::world(&m, Axis::Z);
        let s = prism_on(&mut m, frame, upright_slot(2.0, 0.0, 10.0), 3.0);
        let s = transform(&mut m, s, &rot_iso(axis, 30)).expect("a turn");
        m.rebuild_adjacency();
        let corners = pierce_corners(&m, s);
        assert_eq!(corners.len(), 8, "{axis:?}");
        for v in &corners {
            assert!(
                matches!(m.vertex(*v), Vertex::Pierce { rank: 0, .. }),
                "{axis:?}: a tangent corner is its line's one crossing: {:?}",
                m.vertex(*v)
            );
            assert!(
                crate::realize_vertex(&m, *v, crate::Precision::NearestF64).is_ok(),
                "{axis:?}: {:?}",
                m.vertex(*v)
            );
        }
        let report = crate::refine_caches_of(&mut m, &[s]);
        assert_eq!(report.vertices.left_unrealized, 0, "{axis:?}: {report:?}");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{axis:?}: {:?}",
            nacre_validate::validate(&m)
        );
    }
}

/// **A slot stands on every face of a turned box.** The box is turned 30° about `x`, so its
/// faces sit on the turn's node; a sketch on one of them is a frame node of its own, the slot's
/// walls and cylinders stand on it, and the base cap is the box's face — a different chain. The
/// rank is the wall's and the cylinder's alone (their shared ruling), so each corner is named,
/// and the realization meets the ruling with the face in the world.
#[test]
fn a_slot_stands_on_every_face_of_a_turned_box() {
    let r = Rat::from_int;
    let mut m = Model::new();
    let block = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-10.0, -10.0, -10.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    let block = fixtures::xf(
        &mut m,
        block,
        nacre_exact::Isometry::rotation(nacre_exact::Rotation {
            axis: Axis::X,
            pivot: [r(0), r(0), r(0)],
            angle: nacre_exact::Angle::from_deg(r(30)).expect("an angle"),
        }),
    );
    m.rebuild_adjacency();
    let faces = m.shell(m.solid(block).outer).faces.clone();
    assert_eq!(faces.len(), 6);
    for fh in faces {
        let frame = crate::face_sketch_frame(&m, fh).expect("a planar live face");
        let slot = prism_on(&mut m, frame, upright_slot(1.0, -3.0, 3.0), 2.0);
        m.rebuild_adjacency();
        let want = (2.0 * 6.0 + std::f64::consts::PI) * 2.0;
        let got = nacre_props::mass_props(&m, slot).expect("props").volume;
        assert!((got - want).abs() < 1e-9, "{fh:?}: volume {got} vs {want}");
        for v in pierce_corners(&m, slot) {
            let again = crate::realize_vertex(&m, v, crate::Precision::NearestF64)
                .unwrap_or_else(|e| panic!("{fh:?}: {e:?}"))
                .to_f64()
                .map(|(x, _)| x);
            assert_eq!(again, Some(m.vertex_point(v).as_array()), "{fh:?}");
        }
    }
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
}
