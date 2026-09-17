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
use nacre_ops::{DatumDef, SketchFrame};

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}
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
        for &fh in &m.shell(m.solid(s).outer).faces {
            let f = m.face(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    seen += 1;
                    let e = m.edge(he.edge);
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
        for &fh in &m.shell(m.solid(s).outer).faces {
            let f = m.face(fh);
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
            m.edge(eh).surfaces[..],
            "{what}: stated carriers != observed adjacency"
        );
    }
    assert!(seen > 0, "{what}: the sweep saw no half-edges at all");
}

/// S8 커밋-2 잠금: **파생 곡선 == 저장 곡선** — 담체·끝점에서 `derive_edge_curve` 로 다시
/// 이끌어낸 곡선이, 생산자가 저장한 곡선과 일치한다.
///
/// * 직선: **원점 비트 동일** + 방향 상대 편차 ≤ 1e-12. 첫-구성 직선은 방향까지 비트
///   동일하다(파생이 생산자와 같은 표현식·같은 끝점 순서 — 카운터가 센다). **이동된**
///   직선은 다르다 — 실측 반박: transform pass 2 는 방향 벡터를 직접 회전하고 파생은
///   이동된 끝점 차의 재정규화라 마지막 ulp 가 갈린다. 직선 기하는 어디서도 좌표로
///   관측되지 않으므로(S8 조사) ulp 는 무해하고, 여기서 수치로 못박는다.
/// * 원: 상대 편차 ≤ 1e-12 (코퍼스-수치 게이트; 최대 편차를 찍는다 — 관문 규칙).
struct DeriveStats {
    max_circle_dev: f64,
    max_line_dir_dev: f64,
    lines_bit_identical: usize,
    lines_total: usize,
}

fn assert_derived_matches_stored(m: &Model, what: &str, st: &mut DeriveStats) {
    use nacre_geom::Curve;

    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(b.abs()).max(1.0);
    let mut seen = std::collections::HashSet::new();
    for &s in &m.live_solids {
        for &fh in &m.shell(m.solid(s).outer).faces {
            let f = m.face(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    if !seen.insert(he.edge) {
                        continue;
                    }
                    let e = m.edge(he.edge);
                    let derived = m
                        .derive_edge_curve(e.surfaces, e.vertices)
                        .expect("a live edge's curve must derive");
                    match (m.edge_curve(he.edge), &derived) {
                        (Curve::Line(stored), Curve::Line(d)) => {
                            st.lines_total += 1;
                            assert_eq!(
                                stored.origin(),
                                d.origin(),
                                "{what}: a derived line's origin differs from the stored one"
                            );
                            let mut dev = 0.0f64;
                            for k in 0..3 {
                                dev = dev.max(rel(
                                    stored.direction().as_array()[k],
                                    d.direction().as_array()[k],
                                ));
                            }
                            assert!(
                                dev <= 1e-12,
                                "{what}: a derived line's direction deviates {dev:e}"
                            );
                            st.max_line_dir_dev = st.max_line_dir_dev.max(dev);
                            if dev == 0.0 {
                                st.lines_bit_identical += 1;
                            }
                        }
                        (Curve::Circle(stored), Curve::Circle(d)) => {
                            let mut dev: f64 = rel(stored.radius(), d.radius());
                            for k in 0..3 {
                                dev = dev
                                    .max(rel(
                                        stored.center().as_array()[k],
                                        d.center().as_array()[k],
                                    ))
                                    .max(rel(
                                        stored.normal().as_array()[k],
                                        d.normal().as_array()[k],
                                    ))
                                    .max(rel(
                                        stored.ref_dir().as_array()[k],
                                        d.ref_dir().as_array()[k],
                                    ));
                            }
                            assert!(
                                dev <= 1e-12,
                                "{what}: a derived circle deviates {dev:e} from the stored one"
                            );
                            st.max_circle_dev = st.max_circle_dev.max(dev);
                        }
                        (stored, d) => {
                            panic!("{what}: curve kind changed: {stored:?} vs {d:?}")
                        }
                    }
                }
            }
        }
    }
    assert!(st.lines_total > 0, "{what}: the sweep derived no lines");
}

#[test]
fn derived_curves_match_stored() {
    let mut st = DeriveStats {
        max_circle_dev: 0.0,
        max_line_dir_dev: 0.0,
        lines_bit_identical: 0,
        lines_total: 0,
    };

    // ① Construction, including a tilted (irrational-normalization) cylinder axis.
    let mut m = Model::new();
    m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    m.add_cylinder(
        Point3::from_array([10.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.5,
        2.0,
    );
    let tilted_cyl = m.add_cylinder(
        Point3::from_array([20.0, 1.0, 0.5]),
        Vector3::from_array([0.3, -0.4, 1.0]),
        0.7,
        3.0,
    );
    m.rebuild_adjacency();
    assert_derived_matches_stored(&m, "construction", &mut st);

    // ② A tilted extrude and a boolean (Discovered edges).
    let tilted = SketchPlane::from_origin_normal(
        Point3::from_array([0.25, -0.5, 1.5]),
        Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
    )
    .unwrap();
    let prism = {
        let __frame1 = datum_frame(&mut m, tilted);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __frame1,
                profile: square(0.5, 2.5),
                dist: 1.1,
            },
        )
        .expect("tilted extrude") else {
            unreachable!()
        };
        solid
    };
    let cutter = {
        let __w2 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __w2,
                profile: square(1.0, 2.0),
                dist: 4.0,
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
            a: prism,
            b: cutter,
        },
    )
    .expect("cut");
    m.rebuild_adjacency();
    assert_derived_matches_stored(&m, "extrude+cut", &mut st);

    // ③ ★ Moved cylinders — the one commit where pass 2 (`transform_curve`) still exists, so
    //    «derive from moved caches» vs «pass 2's directly-transformed circle» can be compared.
    //    (A mirrored cylinder has no population: `MirrorNotPlanar` rejects it.)
    let moved = {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: tilted_cyl,
                isometry: Isometry::translation([
                    Rat::from_decimal(1.25).unwrap(),
                    Rat::from_decimal(-0.5).unwrap(),
                    Rat::from_decimal(2.0).unwrap(),
                ]),
            },
        )
        .expect("translate cylinder") else {
            unreachable!()
        };
        solid
    };
    apply(
        &mut m,
        &Operation::Transform {
            solid: moved,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
            }),
        },
    )
    .expect("rotate cylinder");
    m.rebuild_adjacency();
    assert_derived_matches_stored(&m, "moved cylinders", &mut st);

    println!(
        "stat derived_curves lines {}/{} bit-identical, max line-dir dev {:e}, max circle dev {:e}",
        st.lines_bit_identical, st.lines_total, st.max_line_dir_dev, st.max_circle_dev
    );
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
    let __frame0 = datum_frame(&mut m, tilted);
    apply(
        &mut m,
        &Operation::Extrude {
            frame: __frame0,
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
        let __w1 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __w1,
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
        let __w0 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __w0,
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
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
            }),
        },
    )
    .expect("rotate");
    m2.rebuild_adjacency();
    assert_carriers_agree(&m2, "transform");
}
