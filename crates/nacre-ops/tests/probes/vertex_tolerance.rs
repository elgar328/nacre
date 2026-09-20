//! **A vertex's recorded tolerance measures the planes that vertex is defined by.**
//!
//! The coordinate is a cache; the truth is "where these three surfaces meet". The tolerance is the
//! bridge — it says how far the cache sits from the truth — and `nacre-validate` holds the model to
//! it. So the number has to be measured against the surfaces the vertex is *defined* by, and for a
//! while it was not: the arrangement computed the point from one plane triple and measured the
//! tolerance against that triple, then the assembly **re-named** the vertex in the result's own
//! surfaces — a different triple wherever four planes concur — and carried the old figure across.
//!
//! The carried figure was fine while nothing was rotated: an axis-aligned plane passes exactly
//! through the points it defines, so re-naming changed nothing to measure. A rotated plane passes
//! within an ulp or two instead, and then the tolerance bounded a distance nobody would check while
//! the checker measured a different one. 454 of 91,394 result vertices in the corpus were short
//! (by at most 2.1e-14, about two ulps at these coordinates) — and a `replay` proptest run
//! eventually drew one and failed.
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

/// **Every vertex of the model sits within its own recorded tolerance of every surface it is
/// defined by** — the proposition `nacre-validate` enforces, asserted here directly so the failure
/// names this file rather than arriving as a generic "model is invalid".
///
/// Returns how many vertices were **realized from their definition**, so a fixture that stopped
/// producing any cannot pass by measuring nothing.
///
/// ⚠ It used to count "vertices carrying a recorded tolerance". The cache stores none now (cell
/// 54), so the population this walks is the realized one and everything else is skipped — the
/// callers' `> 0` is what keeps the skip from swallowing the whole fixture.
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
                "vertex {} at {:?} is {residual:e} from surface {} but records tol {tol:e}",
                vh.index(),
                m.vertex_point(vh).as_array(),
                sh.index(),
            );
        }
    }
    measured
}

/// ★ **The shape that actually produces a short tolerance.**
///
/// Chosen by measurement, not by reasoning: my first attempt here was two boxes with one rotated
/// 30/45/60°, which reads like the right thing and **passes without the fix** — it never makes a
/// vertex whose definition is re-derived onto a plane the tolerance never measured. This one does.
/// It is the staircase and prism of `concurrent_line.rs`: the prism's apex edge pierces the
/// staircase's step face, four planes meet at that point, and the assembly re-names the vertex.
///
/// ★★ Verified by removing the fix and watching this go red. A lock that stays green either way
/// measures nothing, which is how the first version of this file was written.
#[test]
fn a_pierced_step_records_tolerances_that_hold() {
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
        assert!(n > 0, "apex {apex}: no vertex carried a measured tolerance");
        assert!(nacre_validate::validate(&m).is_empty(), "apex {apex}");
    }
}

/// The negative control: axis-aligned, where the carried figure was always right. If this ever
/// fails, the fix broke the ordinary case rather than the rotated one.
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
