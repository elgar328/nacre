//! S8 커밋-1 잠금: **모서리의 담체는 인접성과 일치한다** — 모든 생산 경로에서, 살아 있는
//! 모든 모서리에 대해, `edge.surfaces` (생산자가 실어 준 쌍) == 그 모서리를 쓰는 두 면의
//! surface 다중집합 (인접성이 관측한 쌍). seam 은 같은 면이 두 번 쓰므로 [cyl, cyl] 로
//! 일치한다.
//!
//! 담체는 파생이 아니라 진술이므로(4평면 동시성 — `Ring.walls` doc), 이 잠금은 진술과
//! 관측이 갈라지는 순간을 생산 경로별로 잡는 반증 장치다. 특히 boolean 은 정리 패스
//! (`unify_coplanar_faces`)가 면-surface 를 재배선한 **뒤의** 결과를 검사한다.

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_topo::Model;

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
    .expect("a square is a fair profile")
}

/// Every live edge's stated carriers vs the multiset of face surfaces adjacency sees using it.
fn assert_carriers_agree(m: &Model, what: &str) {
    let mut seen = 0usize;
    for &s in &m.live_solids {
        for &fh in &m.shells.get(m.solids.get(s).outer).faces {
            let f = m.faces.get(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    seen += 1;
                    let e = m.edges.get(he.edge);
                    // This use's face surface must be one of the stated carriers.
                    assert!(
                        e.surfaces.contains(&f.surface),
                        "{what}: a face rides an edge whose carriers do not name it"
                    );
                }
            }
        }
    }
    // And globally: collect uses per edge, compare multisets.
    let mut uses: std::collections::HashMap<_, Vec<_>> = std::collections::HashMap::new();
    for &s in &m.live_solids {
        for &fh in &m.shells.get(m.solids.get(s).outer).faces {
            let f = m.faces.get(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    uses.entry(he.edge).or_default().push(f.surface);
                }
            }
        }
    }
    for (eh, mut face_surfs) in uses {
        assert_eq!(
            face_surfs.len(),
            2,
            "{what}: a live manifold edge has exactly two uses"
        );
        face_surfs.sort_by_key(|h| h.index());
        assert_eq!(
            face_surfs[..],
            m.edges.get(eh).surfaces[..],
            "{what}: stated carriers != observed adjacency"
        );
    }
    assert!(seen > 0, "{what}: the sweep saw no half-edges at all");
}

#[test]
fn edge_carriers_agree_with_adjacency() {
    // ① Plain construction: cuboid + cylinder.
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    m.add_cylinder(
        Point3::from_array([10.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.5,
        2.0,
    );
    m.rebuild_adjacency();
    assert_carriers_agree(&m, "construction");

    // ② A tilted extrude (the frame road) — walls/caps carriers from `sweep_ring`.
    let tilted = SketchPlane::from_origin_normal(
        Point3::from_array([0.25, -0.5, 1.5]),
        Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
    )
    .unwrap();
    apply(
        &mut m,
        &Operation::Extrude {
            plane: tilted,
            profile: square(0.5, 2.5),
            dist: 1.1,
        },
    )
    .expect("tilted extrude");
    m.rebuild_adjacency();
    assert_carriers_agree(&m, "tilted extrude");

    // ③ Booleans, including populations the cleaning pass actually rewrites:
    //    a coplanar fuse (shared plane, unify_coplanar_faces merges), then a severing cut.
    let other = {
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: square(2.0, 6.0),
                dist: 1.0,
            },
        )
        .expect("coplanar operand") else {
            unreachable!()
        };
        solid
    };
    let fused = {
        let OpOutput::Boolean { solids } = apply(
            &mut m,
            &Operation::Boolean {
                kind: BoolKind::Fuse,
                a: base,
                b: other,
            },
        )
        .expect("coplanar fuse") else {
            unreachable!()
        };
        solids[0]
    };
    m.rebuild_adjacency();
    assert_carriers_agree(&m, "coplanar fuse (cleaned)");

    let cutter = {
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: square(2.5, 3.5),
                dist: 5.0,
            },
        )
        .expect("cutter") else {
            unreachable!()
        };
        solid
    };
    apply(
        &mut m,
        &Operation::Boolean {
            kind: BoolKind::Cut,
            a: fused,
            b: cutter,
        },
    )
    .expect("cut");
    m.rebuild_adjacency();
    assert_carriers_agree(&m, "cut");

    // ④ Motions: a rotated copy (surf_map re-pointing) — cuboid family and cylinder family.
    let mut m2 = Model::new();
    let box2 = m2.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 2.0, 3.0]),
    );
    m2.add_cylinder(
        Point3::from_array([5.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    m2.rebuild_adjacency();
    apply(
        &mut m2,
        &Operation::Transform {
            solid: box2,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
            }),
        },
    )
    .expect("rotate");
    m2.rebuild_adjacency();
    assert_carriers_agree(&m2, "transform");
}
