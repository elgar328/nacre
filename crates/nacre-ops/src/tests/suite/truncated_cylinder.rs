//! **A cylinder cut by a plane oblique to its axis, built by hand — until a boolean builds one.**
//!
//! The consumers of an elliptical edge (mesh, measures, `validate`, STEP) are locked here on the
//! one shape that is all ellipse: the product's cylinder (radius 2 about `z`, seam `+x`, `z ∈ [0,
//! 10]`) with its top cap replaced by the plane `z = 5 + k·x` — `k = 1/2`, and `k = 2` (the cap
//! 63° off square, past where a winding read along the axis would call it misoriented). A class's stopper refuses every
//! boolean result that keeps such an edge, so no operation makes this shape yet; it is the one
//! exception to «fixtures go through the product's operations» (`crate::fixtures`), and it is
//! assembled the product's way — every cell copied from the product's cylinder but the top rim,
//! whose vertex and edge go through the realization funnel (`push_vertex_realized`,
//! `push_edge_realized`) exactly as the assembly pushes a rim. Retired when booleans build these
//! shapes and the locks move onto them.
//!
//! Derived, not copied: the top rim is `z(θ) = 5 + 2k cos θ` (`x = 2 cos θ`), so the solid's volume
//! is `πr²·5 = 20π`, its lateral area `2πr·5 = 20π`, the elliptical cap's `πr·a` with `a = r√(1 +
//! k²)`, and the centroid `x̄ = ∬x·z_top / V = 4πk / 20π = k/5`, `z̄ = ∬z_top²/2 / V = (100π +
//! 4πk²) / 40π = 2.5 + k²/10`.

use super::*;
use crate::realize::{ChainLink, PlaneMemo, push_edge_realized, push_vertex_realized};
use nacre_topo::{
    Face, HalfEdge, Loop, Orientation, PointCache, Shell, Solid as TopoSolid, Vertex,
};

const PI: f64 = std::f64::consts::PI;
const SLOPES: [f64; 2] = [0.5, 2.0];

fn truncated(k: f64) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let c = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.0,
        10.0,
    );
    let plane = match crate::apply(
        &mut m,
        &crate::Operation::DatumPlane {
            def: crate::DatumDef::Stated(
                crate::SketchPlane::from_origin_normal(
                    Point3::from_array([0.0, 0.0, 5.0]),
                    Vector3::from_array([-k, 0.0, 1.0]),
                )
                .expect("a plane"),
            ),
        },
    ) {
        Ok(crate::OpOutput::DatumPlane { plane, .. }) => plane,
        other => panic!("the oblique cap's datum: {other:?}"),
    };
    let (base, top, lateral) = (c.base, m.face(c.top).clone(), m.face(c.lateral).clone());
    let lat = lateral.surface;
    let v = push_vertex_realized(
        &mut m,
        Vertex::OnSeam([lat, plane]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.0, 0.0, 5.0 + 2.0 * k]),
        },
        ChainLink::Fresh,
    );
    let e = push_edge_realized(&mut m, [lat, plane], [v, v], &mut PlaneMemo::default())
        .expect("an elliptical rim derives its curve");
    let old_rim = top.outer.half_edges[0];
    let swap = |lp: &Loop| Loop {
        half_edges: lp
            .half_edges
            .iter()
            .map(|he| HalfEdge {
                edge: if he.edge == old_rim.edge { e } else { he.edge },
                forward: he.forward,
            })
            .collect(),
    };
    let up = m
        .surface_cache(plane)
        .plane()
        .expect("a plane")
        .normal()
        .as_array()[2]
        > 0.0;
    let new_top = m.push_face(Face {
        surface: plane,
        outer: swap(&top.outer),
        inner: Vec::new(),
        orientation: if up {
            Orientation::Forward
        } else {
            Orientation::Reversed
        },
    });
    let new_lateral = m.push_face(Face {
        surface: lat,
        outer: swap(&lateral.outer),
        inner: lateral.inner.iter().map(swap).collect(),
        orientation: lateral.orientation,
    });
    let shell = m.push_shell(Shell {
        faces: vec![base, new_top, new_lateral],
    });
    let solid = m.push_solid(TopoSolid {
        outer: shell,
        cavities: Vec::new(),
    });
    m.supersede_live(&[c.solid]);
    (m, solid)
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs().max(1.0)
}

/// The shape is a valid solid, its top rim an ellipse (no circle) on the oblique plane, and its
/// seam vertex where the generator through `(2, 0)` meets the plane, `(2, 0, 5 + 2k)`.
#[test]
fn a_hand_built_oblique_cut_cylinder_is_a_valid_solid() {
    for k in SLOPES {
        let (m, s) = truncated(k);
        let violations = nacre_validate::validate(&m);
        assert!(violations.is_empty(), "{violations:?}");
        let ellipses: Vec<_> = m
            .reachable_from(&[s])
            .edges
            .iter()
            .filter_map(|&e| match m.edge_curve(e) {
                nacre_geom::Curve::Ellipse(c) if !c.is_circle() => Some((e, *c)),
                _ => None,
            })
            .collect();
        assert_eq!(ellipses.len(), 1, "one elliptical rim");
        let [v, _] = m.edge(ellipses[0].0).vertices;
        let p = m.vertex_point(v).as_array();
        assert!(
            close(p[0], 2.0, 1e-15) && p[1].abs() < 1e-15 && close(p[2], 5.0 + 2.0 * k, 1e-15),
            "{p:?}"
        );
        assert!((ellipses[0].1.semi_major() - 2.0 * (1.0 + k * k).sqrt()).abs() < 1e-12);
    }
}

/// Volume, area and centroid against the derived values; the faces' own areas and centroids too.
#[test]
fn a_hand_built_oblique_cut_cylinder_measures_as_derived() {
    for k in SLOPES {
        let (m, s) = truncated(k);
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
            close(c[0], x, 1e-12) && c[1].abs() < 1e-12 && close(c[2], z, 1e-12),
            "{k}: {c:?}"
        );
        let (lo, hi) = nacre_props::bounds(&m, s).expect("bounds");
        let (lo, hi) = (lo.as_array(), hi.as_array());
        // The box holds the shape and is no looser than the whole ellipse's reach.
        let top = 5.0 + 2.0 * k;
        assert!(
            lo[0] <= -2.0 && hi[0] >= 2.0 && lo[2] <= 0.0 && hi[2] >= top,
            "{lo:?} {hi:?}"
        );
        assert!(
            hi[2] <= top + 1e-12 && lo[1] >= -2.0 - 1e-12,
            "{lo:?} {hi:?}"
        );
    }
}

/// The mesh: every face meshes, no lateral triangle leans off the cylinder, and the triangles'
/// enclosed volume is the solid's to the mesh's resolution.
#[test]
fn a_hand_built_oblique_cut_cylinder_meshes() {
    for k in SLOPES {
        let (m, s) = truncated(k);
        let cfg = nacre_tess::TessConfig::default();
        let t = nacre_tess::tessellate(&m, &cfg).expect("the mesh");
        let leaning = nacre_tess::leaning_laterals(&m, &t, cfg.max_angle_deg);
        assert!(leaning.is_empty(), "{leaning:?}");
        let faces = m.reachable_from(&[s]).faces;
        let mut volume = 0.0;
        for f in &faces {
            for &tri in &t.by_face[f] {
                let [a, b, c] = t.triangles.get(tri).vertices.map(|v| t.vertices.get(v).pos);
                volume += (a - Point3::origin())
                    .dot((b - Point3::origin()).cross(c - Point3::origin()))
                    / 6.0;
            }
        }
        assert!(close(volume, 20.0 * PI, 2e-3), "{k}: mesh volume {volume}");
    }
}

/// STEP writes the rim as an `ELLIPSE` and proves an uncertainty for it.
#[test]
fn a_hand_built_oblique_cut_cylinder_writes_an_ellipse() {
    for k in SLOPES {
        let (m, s) = truncated(k);
        let text = nacre_step::to_step(&m, "2000-01-01T00:00:00").expect("STEP");
        assert!(text.contains("ELLIPSE("), "no ELLIPSE written");
        match nacre_step::uncertainty(&m, &[s]) {
            nacre_step::Uncertainty::Proven(d) => assert!(d < 1e-12, "uncertainty {d}"),
            other => panic!("not proven: {other:?}"),
        }
    }
}
