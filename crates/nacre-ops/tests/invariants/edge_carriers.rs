//! Lock: **an edge's carriers agree with adjacency** — on every production path, for every live
//! edge, `edge.surfaces` (the pair its producer stated) equals the multiset of surfaces of the two
//! faces that use the edge (the pair adjacency observes).
//!
//! Carriers are a statement, not a derivation (four-plane concurrency — see `Ring.walls`' doc),
//! so this lock is the falsifier that catches, per production path, the moment statement and
//! observation part. For a boolean in particular it checks the result **after** the cleaning pass
//! (`unify_coplanar_faces`) has rewired the faces' surfaces.

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::SketchFrame;
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};

use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_topo::Model;

use crate::fixtures::datum_frame;

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
    for &s in m.live_solids() {
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
    for &s in m.live_solids() {
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

/// Lock: **derived curve == stored curve** — the curve re-derived from the carriers and the end
/// points by `derive_edge_curve`, with this caller giving nothing, agrees with the one the producer
/// stored (a line on a plane without a world name reads the direction its push realized, which the
/// model keeps for the carrier pair; where nothing was kept it is re-derived along its end points —
/// the ulps between the two are the bound below).
///
/// * Line: **origin bit-identical**, direction within a relative 1e-12. A first-construction line
///   is bit-identical in direction too (the derivation uses the producer's expression and end-point
///   order — the counter counts them). A **moved** line is not: transform's second pass rotates
///   the direction vector itself, while the derivation re-normalizes the difference of the moved
///   end points, and the last ulp differs. A line's geometry is never observed through coordinates,
///   so the ulp is harmless, and it is pinned here as a number.
/// * Circle: within a relative 1e-12 (a corpus-number gate; the largest deviation is printed).
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
    for &s in m.live_solids() {
        for &fh in &m.shell(m.solid(s).outer).faces {
            let f = m.face(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    if !seen.insert(he.edge) {
                        continue;
                    }
                    let e = m.edge(he.edge);
                    let derived = m
                        .derive_edge_curve(e.surfaces, e.vertices, |_| nacre_topo::EdgeGiven::NONE)
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
    nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([10.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        1.5,
        2.0,
    );
    let tilted_cyl = nacre_ops::fixtures::cylinder(
        &mut m,
        Point3::from_array([20.0, 1.0, 0.5]),
        Vector3::from_array([0.3, -0.4, 1.0]),
        0.7,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    assert_derived_matches_stored(&m, "construction", &mut st);

    // ② A tilted extrude and a boolean (arrangement-made edges).
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

    // ③ ★ Moved cylinders: «derive from moved caches» against the stored circle.
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
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([10.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
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
    let box2 = nacre_ops::fixtures::cuboid(
        &mut m2,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 2.0, 3.0]),
    );
    nacre_ops::fixtures::cylinder_with_seam(
        &mut m2,
        Point3::from_array([5.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
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
