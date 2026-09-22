//! **A realized vertex sits on the surfaces it is defined by.**
//!
//! The coordinate is a cache; the truth is "where these three surfaces meet". A realized cache is
//! the definition's nearest `f64`, bit for bit, and `nacre-validate` holds it within the
//! construction epsilon (`EPS_CONSTRUCTED`) of those surfaces. So the distance is measured against
//! the surfaces the vertex is *defined* by — the result's own, which the assembly **re-names** it
//! in
//! wherever four planes concur — not against the triple the arrangement computed the point with.
//!
//! On an axis-aligned plane the two agree exactly. A rotated plane passes within an ulp or two, so
//! a figure measured against the computing triple bounds a distance nobody checks while the checker
//! measures a different one: 454 of 91,394 result vertices in the corpus, short by at most 2.1e-14
//! (about two ulps at these coordinates), and a `replay` proptest run drew one.
//!
//! ★ The proptest seed that found it is kept, but a seed is a lucky draw. What this file asserts is
//! the **proposition**, on shapes chosen because they must exercise it.

use nacre_exact::Axis;
use nacre_math::{Point2, Point3};
use nacre_ops::{
    BoolKind, Operation, Precision, Profile2d, SketchFrame, apply, boolean, realize_vertex,
};
use nacre_store::Handle;
use nacre_topo::{Model, PointCache, Solid};
use nacre_validate::EPS_CONSTRUCTED;

fn prism(m: &mut Model, pts: &[[f64; 2]], h: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
            profile,
            dist: h,
        },
    )
    .expect("extrude") else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// **Every realized vertex is its definition's realization, bit for bit, and sits within
/// `EPS_CONSTRUCTED` of every surface it is defined by** — the proposition `nacre-validate`
/// enforces, asserted here directly so the failure names this file rather than arriving as a
/// generic "model is invalid".
///
/// Returns how many vertices were **realized from their definition**, so a fixture that stopped
/// producing any cannot pass by measuring nothing.
///
/// ⚠ Only realized (`Bounded`) vertices are walked — the cache claims nothing about the rest — and
/// the callers' `> 0` is what keeps the skip from swallowing the whole fixture.
fn every_vertex_matches_its_definition(m: &Model) -> usize {
    let mut measured = 0;
    let mut i = 0u32;
    while let Some(vh) = m.vertex_handle_at(i) {
        i += 1;
        let v = m.vertex(vh);
        // What the cache knows, and the figure it is held to: a realized coordinate to the
        // realization bit for bit, and to its carriers within the construction epsilon
        // `nacre-validate` applies to it (its bound speaks of the coordinate, not of the cached
        // carriers). A coordinate with no realization behind it is skipped — the cache claims
        // nothing about it, so there is nothing here to hold it to.
        let tol = match *m.vertex_cache(vh) {
            PointCache::Unrealized { .. } | PointCache::Ceiling { .. } => continue,
            PointCache::Bounded { coord, .. } => {
                let (realized, _) = realize_vertex(m, vh, Precision::NearestF64)
                    .expect("a Bounded vertex realizes")
                    .to_f64()
                    .expect("decided");
                assert_eq!(
                    coord.as_array(),
                    realized,
                    "vertex {}: the cache is the realization",
                    vh.index()
                );
                EPS_CONSTRUCTED
            }
        };
        measured += 1;
        for sh in v.carriers() {
            let residual = m.surface_cache(sh).distance(m.vertex_point(vh));
            assert!(
                residual <= tol,
                "vertex {} at {:?} is {residual:e} from surface {} but is held to {tol:e}",
                vh.index(),
                m.vertex_point(vh).as_array(),
                sh.index(),
            );
        }
    }
    measured
}

/// ★ **The shape that re-names a vertex onto a rotated plane.**
///
/// Chosen by measurement, not by reasoning: two boxes with one rotated 30/45/60° read like the
/// right thing but never make a vertex whose definition is re-derived onto another plane, so they
/// pass either way. This does. It is the staircase and prism of `concurrent_line.rs`: the prism's
/// apex edge pierces the staircase's step face, four planes meet at that point, and the assembly
/// re-names the vertex.
///
/// ★★ Measured against the computing triple instead, this goes red. A lock that stays green either
/// way measures nothing.
#[test]
fn a_pierced_step_realizes_on_its_defining_planes() {
    for apex in [0.7_f64, 1.05, 1.2] {
        let mut m = Model::new();
        let a = prism(
            &mut m,
            &[
                [0.0, 0.0],
                [0.0, 2.0],
                [1.0, 2.0],
                [1.0, 1.0],
                [0.5, 1.0],
                [0.5, 0.5],
                [2.0, 0.5],
                [2.0, 0.0],
            ],
            1.0,
        );
        // ZX: sketch +u is +z, +v is +x — the prism's apex edge runs along y.
        let profile = Profile2d::polygon(
            [[0.5, apex], [0.0, 2.0], [1.0, 2.0]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        )
        .expect("profile");
        let frame = SketchFrame::world(&m, Axis::Y);
        let nacre_ops::OpOutput::Extrude { solid: b, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 0.9,
            },
        )
        .expect("extrude") else {
            panic!()
        };
        m.rebuild_adjacency();
        let got = boolean(&mut m, BoolKind::Fuse, a, b).expect("the fuse builds");
        m.rebuild_adjacency();
        assert!(!got.is_empty());
        let n = every_vertex_matches_its_definition(&m);
        assert!(n > 0, "apex {apex}: no vertex was realized");
        assert!(nacre_validate::validate(&m).is_empty(), "apex {apex}");
    }
}

/// The negative control: axis-aligned, where the computing triple and the re-named one agree
/// exactly. If this fails, the ordinary case broke rather than the rotated one.
#[test]
fn an_axis_aligned_fuse_is_unchanged() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, 1.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    m.rebuild_adjacency();
    let got = boolean(&mut m, BoolKind::Fuse, a, b).expect("the fuse builds");
    m.rebuild_adjacency();
    assert_eq!(got.len(), 1);
    assert!(every_vertex_matches_its_definition(&m) > 0);
    assert!(nacre_validate::validate(&m).is_empty());
}
