//! **A rod's end cut by a plane oblique to its axis** — the one shape whose rim is a whole ellipse,
//! so every consumer of an elliptical edge (mesh, measures, `validate`, STEP) is locked on it.
//!
//! The rod: radius 2 about `z` through `(0, 3)`, seam `+x`, `z ∈ [0, 10]`. The cut: a prism over
//! the plane `z = 5 + k·x` — its section in `zx` the region between that plane and `z = Z` over
//! `x ∈ [−3, 3]`, extruded over `y ∈ [0, 6]` — for `k = 1/2`, and `k = 2` (the cap 63° off square,
//! past where a winding read along the axis would call it misoriented).
//!
//! Derived, not copied: the rod's top rim is `z(θ) = 5 + 2k cos θ`, so the part below the plane has
//! volume `πr²·5 = 20π` (the part above too), its lateral area `2πr·5 = 20π`, the elliptical cap's
//! `πr·a` with `a = r√(1 + k²)`, and centroid `x̄ = ∬x·z_top / V = k/5`, `z̄ = ∬z_top²/2 / V =
//! 2.5 + k²/10`; the prism is `6·6·(Z − 5)`.

use super::*;

const PI: f64 = std::f64::consts::PI;
const SLOPES: [f64; 2] = [0.5, 2.0];

/// The prism's roof: clear of the rod's top and of the plane's highest point over `x ∈ [−3, 3]`.
fn roof(k: f64) -> f64 {
    12.0f64.max(5.0 + 3.0 * k + 1.0)
}

fn rod(m: &mut Model) -> Handle<Solid> {
    crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([0.0, 3.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.0,
        10.0,
    )
    .solid
}

fn prism(m: &mut Model, k: f64) -> Handle<Solid> {
    let top = roof(k);
    // The world `zx` frame: `u = z`, `v = x`, extruded along `+y`.
    let pts = [
        (5.0 - 3.0 * k, -3.0),
        (5.0 + 3.0 * k, 3.0),
        (top, 3.0),
        (top, -3.0),
    ]
    .map(|(u, v)| nacre_math::Point2::from_array([u, v]));
    let ring = crate::Ring2d::polygon_decimal(pts.to_vec()).expect("a quadrilateral");
    let profile = crate::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = SketchFrame::world(m, Axis::Y);
    let Ok(OpOutput::Extrude { solid, .. }) = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 6.0,
        },
    ) else {
        panic!("the prism extrudes")
    };
    solid
}

/// The rod less the prism: the rod cut obliquely at `z = 5 + k·x`.
fn cut_rod(k: f64) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let (r, p) = (rod(&mut m), prism(&mut m, k));
    let out = crate::boolean(&mut m, BoolKind::Cut, r, p).expect("the rod cuts");
    assert_eq!(out.len(), 1);
    (m, out[0])
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs().max(1.0)
}

/// Every kind in both orders builds, at the derived volumes.
#[test]
fn a_rods_end_cut_obliquely_builds_every_kind() {
    for k in SLOPES {
        let p = 36.0 * (roof(k) - 5.0);
        every_kind(
            |m: &mut Model| (rod(m), prism(m, k)),
            [
                [
                    Want::Built(1, 20.0 * PI + p),
                    Want::Built(1, 20.0 * PI),
                    Want::Built(1, 20.0 * PI),
                ],
                [
                    Want::Built(1, 20.0 * PI + p),
                    Want::Built(1, p - 20.0 * PI),
                    Want::Built(1, 20.0 * PI),
                ],
            ],
        );
    }
}

/// The cut rod is a valid solid, its top rim one whole ellipse on the oblique plane, and that
/// rim's seam vertex where the generator through `(2, 3)` meets the plane, `(2, 3, 5 + 2k)`.
#[test]
fn a_rods_oblique_end_is_one_whole_ellipse() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let violations = nacre_validate::validate(&m);
        assert!(violations.is_empty(), "{k}: {violations:?}");
        let ellipses: Vec<_> = m
            .reachable_from(&[s])
            .edges
            .iter()
            .filter_map(|&e| match m.edge_curve(e) {
                nacre_geom::Curve::Ellipse(c) if !c.is_circle() => Some((e, *c)),
                _ => None,
            })
            .collect();
        assert_eq!(ellipses.len(), 1, "{k}: one elliptical rim");
        let [v0, v1] = m.edge(ellipses[0].0).vertices;
        assert_eq!(v0, v1, "{k}: a whole rim");
        assert!(
            matches!(m.vertex(v0), nacre_topo::Vertex::OnSeam(_)),
            "{k}: its vertex is the seam's"
        );
        let p = m.vertex_point(v0).as_array();
        assert!(
            close(p[0], 2.0, 1e-15) && close(p[1], 3.0, 1e-15) && close(p[2], 5.0 + 2.0 * k, 1e-15),
            "{k}: {p:?}"
        );
        assert!((ellipses[0].1.semi_major() - 2.0 * (1.0 + k * k).sqrt()).abs() < 1e-12);
    }
}

/// Volume, area and centroid against the derived values, and a box that holds the shape.
#[test]
fn a_rods_oblique_end_measures_as_derived() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let mp = nacre_props::mass_props(&m, s).expect("props");
        assert!(
            close(mp.volume, 20.0 * PI, 1e-12),
            "{k}: volume {}",
            mp.volume
        );
        let area = 4.0 * PI + 20.0 * PI + 4.0 * (1.0 + k * k).sqrt() * PI;
        assert!(
            close(mp.area, area, 1e-12),
            "{k}: area {} vs {area}",
            mp.area
        );
        let c = nacre_props::centroid(&m, s).expect("centroid").as_array();
        let (x, z) = (k / 5.0, 2.5 + k * k / 10.0);
        assert!(
            close(c[0], x, 1e-12) && close(c[1], 3.0, 1e-12) && close(c[2], z, 1e-12),
            "{k}: {c:?}"
        );
        let (lo, hi) = nacre_props::bounds(&m, s).expect("bounds");
        let (lo, hi) = (lo.as_array(), hi.as_array());
        let top = 5.0 + 2.0 * k;
        assert!(
            lo[0] <= -2.0 && hi[0] >= 2.0 && lo[2] <= 0.0 && hi[2] >= top,
            "{lo:?} {hi:?}"
        );
        assert!(
            hi[2] <= top + 1e-12 && lo[1] >= 1.0 - 1e-12,
            "{lo:?} {hi:?}"
        );
    }
}

/// The mesh: no lateral triangle leans off the cylinder, and the triangles enclose the solid's
/// volume to the mesh's resolution.
#[test]
fn a_rods_oblique_end_meshes() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let cfg = nacre_tess::TessConfig::default();
        let t = nacre_tess::tessellate(&m, &cfg).expect("the mesh");
        let leaning = nacre_tess::leaning_laterals(&m, &t, cfg.max_angle_deg);
        assert!(leaning.is_empty(), "{leaning:?}");
        let mut volume = 0.0;
        for f in &m.reachable_from(&[s]).faces {
            for &tri in &t.by_face[f] {
                let [a, b, c] = t.triangles.get(tri).vertices.map(|v| t.vertices.get(v).pos);
                let o = Point3::origin();
                volume += (a - o).dot((b - o).cross(c - o)) / 6.0;
            }
        }
        assert!(close(volume, 20.0 * PI, 2e-3), "{k}: mesh volume {volume}");
    }
}

/// STEP writes the rim as an `ELLIPSE` and proves an uncertainty for it.
#[test]
fn a_rods_oblique_end_writes_an_ellipse() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let text = nacre_step::to_step(&m, "2000-01-01T00:00:00").expect("STEP");
        assert!(text.contains("ELLIPSE("), "{k}: no ELLIPSE written");
        match nacre_step::uncertainty(&m, &[s]) {
            nacre_step::Uncertainty::Proven(d) => assert!(d < 1e-12, "{k}: uncertainty {d}"),
            other => panic!("{k}: not proven: {other:?}"),
        }
    }
}

/// **A result with an elliptical edge is refused as an operand, at the door** — every kind, both
/// orders, against a box beside it: `ObliqueCylinderCut`, raised by the boolean's entrance (the
/// stages behind it read a curve on a cylinder as a circle across its axis), not by any of them.
#[test]
fn a_cut_rod_fed_back_is_refused_at_the_door() {
    use crate::error::raise_probe::RAISED;
    crate::ledger::owned(|| {
        for swap in [false, true] {
            for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                let (mut m, s) = cut_rod(0.5);
                let block = crate::fixtures::cuboid(&mut m, [1.0, 1.0, -1.0], [4.0, 5.0, 2.0]);
                let (x, y) = if swap { (block, s) } else { (s, block) };
                let r0 = RAISED.len();
                match crate::boolean(&mut m, kind, x, y) {
                    Err(BoolError::Rejected {
                        reason: RejectReason::ObliqueCylinderCut,
                        ..
                    }) => {}
                    other => panic!("{kind:?} swap={swap}: {other:?}"),
                }
                let sites: Vec<&str> = RAISED.mine_since(r0).into_iter().map(|(_, f)| f).collect();
                assert!(
                    sites
                        .iter()
                        .all(|f| f.ends_with("nacre-ops/src/boolean.rs")),
                    "{kind:?} swap={swap}: raised at {sites:?}"
                );
            }
        }
    });
}

/// **A moved cut rod keeps its edges' kinds** — translated by a third and turned a quarter about
/// `z`, it validates, keeps one elliptical rim (re-derived from its moved carriers, the kind theirs),
/// and its volume.
#[test]
fn a_moved_cut_rod_keeps_its_elliptical_rim() {
    use nacre_exact::{Isometry, Rat};
    let third = Rat::new(1, 3).expect("a third");
    let zero = Rat::from_int(0);
    let (mut m, mut s) = cut_rod(0.5);
    for iso in [
        Isometry::translation([third, zero, third]),
        super::fixtures::rot_iso(Axis::Z, 90),
    ] {
        s = super::fixtures::xf(&mut m, s, iso);
    }
    let violations = nacre_validate::validate(&m);
    assert!(violations.is_empty(), "{violations:?}");
    let ellipses = m
        .reachable_from(&[s])
        .edges
        .iter()
        .filter(|&&e| matches!(m.edge_curve(e), nacre_geom::Curve::Ellipse(c) if !c.is_circle()))
        .count();
    assert_eq!(ellipses, 1);
    let v = nacre_props::mass_props(&m, s).expect("props").volume;
    assert!(close(v, 20.0 * PI, 1e-12), "{v}");
}
