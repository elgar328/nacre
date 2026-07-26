//! Boolean-algebra invariants over **grid-sampled** solids (property tests).
//!
//! Random *float* boxes almost never land on an exact coordinate coincidence, so they only explore
//! the generic space the engine already handles; the silent-wrongs live at coincidences (shared
//! planes, coincident corners). So this samples small **integer** coordinates, where coplanar
//! contact / corner-edge-face touch / containment / separation all occur with positive probability.
//!
//! A reject is never a silent-wrong (those are `Ok(invalid)` or `Ok(wrong volume)`), so the net is:
//! on `Ok`, the result must `validate` clean and satisfy the algebra; on `Err`, skip. `family #3`
//! (never assume rejects away) is honored by the explicit `.expect(Ok)` success pins below — if the
//! engine ever starts rejecting the core cases, those fail rather than the net going vacuously green.

#![allow(unused_imports)]
use crate::common::*;
use nacre_math::{Point2, Point3};
use nacre_ops::{BoolKind, OpOutput, Profile2d, apply, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};
use proptest::prelude::*;

/// An L-profile: the rectangle `[0,w]×[0,h]` minus the top-right notch `[nx,w]×[ny,h]` (a reflex,
/// non-convex operand). Area = `w·h − (w−nx)·(h−ny)`.
fn l_profile(w: i64, h: i64, nx: i64, ny: i64) -> Profile2d {
    let p = |x: i64, y: i64| Point2::from_array([x as f64, y as f64]);
    Profile2d::polygon(vec![
        p(0, 0),
        p(w, 0),
        p(w, ny),
        p(nx, ny),
        p(nx, h),
        p(0, h),
    ])
}

/// A cuboid from integer grid coordinates.
fn grid_box(m: &mut Model, lo: [i64; 3], ext: [i64; 3]) -> Handle<Solid> {
    let f = |v: i64| v as f64;
    m.add_cuboid(
        Point3::from_array([f(lo[0]), f(lo[1]), f(lo[2])]),
        Point3::from_array([f(lo[0] + ext[0]), f(lo[1] + ext[1]), f(lo[2] + ext[2])]),
    )
}

fn box_vol(ext: [i64; 3]) -> f64 {
    (ext[0] * ext[1] * ext[2]) as f64
}

/// AABB overlap volume — the independent oracle (pure integer arithmetic).
fn aabb_overlap(alo: [i64; 3], aext: [i64; 3], blo: [i64; 3], bext: [i64; 3]) -> f64 {
    let mut p: i64 = 1;
    for i in 0..3 {
        let ov = (alo[i] + aext[i]).min(blo[i] + bext[i]) - alo[i].max(blo[i]);
        p *= ov.max(0);
    }
    p as f64
}

fn total_vol(m: &Model, solids: &[Handle<Solid>]) -> f64 {
    solids
        .iter()
        .map(|&s| nacre_props::mass_props(m, s).unwrap().volume)
        .sum()
}

/// Run one boolean; on `Ok`, require the result to be valid and return its total volume; on `Err`,
/// return `None` (a reject is never a silent-wrong). Panics (a proptest failure) on an invalid Ok.
fn run(m: &mut Model, kind: BoolKind, a: Handle<Solid>, b: Handle<Solid>) -> Option<f64> {
    match boolean(m, kind, a, b) {
        Ok(solids) => {
            m.rebuild_adjacency();
            let issues = nacre_validate::validate(m);
            assert!(
                issues.is_empty(),
                "boolean {kind:?} returned an invalid solid: {issues:?}"
            );
            Some(total_vol(m, &solids))
        }
        Err(_) => None,
    }
}

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (a.abs() + b.abs() + 1.0)
}

proptest! {
    /// vol(Fuse)=vA+vB−ovlp, vol(Common)=ovlp, vol(Cut)=vA−ovlp against the independent AABB oracle,
    /// over every grid overlap style (separate / contained / partial / face-contact); corner/edge
    /// contact is non-manifold → `Err` → skipped.
    #[test]
    fn inclusion_exclusion(
        alo in prop::array::uniform3(0i64..5), aext in prop::array::uniform3(1i64..5),
        blo in prop::array::uniform3(0i64..5), bext in prop::array::uniform3(1i64..5),
    ) {
        let build = || {
            let mut m = Model::new();
            let a = grid_box(&mut m, alo, aext);
            let b = grid_box(&mut m, blo, bext);
            (m, a, b)
        };
        let (va, vb) = (box_vol(aext), box_vol(bext));
        let ovlp = aabb_overlap(alo, aext, blo, bext);

        let (mut mf, af, bf) = build();
        if let Some(vf) = run(&mut mf, BoolKind::Fuse, af, bf) {
            prop_assert!(approx(vf, va + vb - ovlp), "fuse {vf} vs {}", va + vb - ovlp);
        }
        let (mut mc, ac, bc) = build();
        if let Some(vc) = run(&mut mc, BoolKind::Common, ac, bc) {
            prop_assert!(approx(vc, ovlp), "common {vc} vs {ovlp}");
        }
        let (mut mk, ak, bk) = build();
        if let Some(vk) = run(&mut mk, BoolKind::Cut, ak, bk) {
            prop_assert!(approx(vk, va - ovlp), "cut {vk} vs {}", va - ovlp);
        }
    }

    /// vol(A) = vol(A∩B) + vol(A∖B), when both are computable.
    #[test]
    fn cut_identity(
        alo in prop::array::uniform3(0i64..5), aext in prop::array::uniform3(1i64..5),
        blo in prop::array::uniform3(0i64..5), bext in prop::array::uniform3(1i64..5),
    ) {
        let build = || {
            let mut m = Model::new();
            let a = grid_box(&mut m, alo, aext);
            let b = grid_box(&mut m, blo, bext);
            (m, a, b)
        };
        let va = box_vol(aext);
        let (mut mc, ac, bc) = build();
        let (mut mk, ak, bk) = build();
        if let (Some(vc), Some(vk)) = (
            run(&mut mc, BoolKind::Common, ac, bc),
            run(&mut mk, BoolKind::Cut, ak, bk),
        ) {
            prop_assert!(approx(va, vc + vk), "A {va} vs common {vc} + cut {vk}");
        }
    }

    /// Fuse and Common commute: `f(a,b)` and `f(b,a)` agree on **both** the Ok/Err verdict and, when
    /// Ok, the volume. An order-dependent volume is a silent-wrong; an order-dependent *verdict* is a
    /// robustness defect (the arrangement is a function of the geometry, not the operand order). A
    /// measure-zero contact — two boxes meeting only at an edge or corner — is the sharp case: it is
    /// empty either way, and both orders must now agree on that (they once didn't, when a pinched
    /// unbounded contour was rejected or accepted depending on its ring's start index; see
    /// `loop_winding`).
    #[test]
    fn fuse_common_commute(
        alo in prop::array::uniform3(0i64..5), aext in prop::array::uniform3(1i64..5),
        blo in prop::array::uniform3(0i64..5), bext in prop::array::uniform3(1i64..5),
    ) {
        let build = || {
            let mut m = Model::new();
            let a = grid_box(&mut m, alo, aext);
            let b = grid_box(&mut m, blo, bext);
            (m, a, b)
        };
        for kind in [BoolKind::Fuse, BoolKind::Common] {
            let (mut m1, a1, b1) = build();
            let (mut m2, a2, b2) = build();
            let (x, y) = (run(&mut m1, kind, a1, b1), run(&mut m2, kind, b2, a2));
            prop_assert_eq!(
                x.is_some(), y.is_some(),
                "{:?} order-dependent verdict: {:?} vs {:?}", kind, x, y
            );
            if let (Some(x), Some(y)) = (x, y) {
                prop_assert!(approx(x, y), "{:?} order-dependent volume {x} vs {y}", kind);
            }
        }
    }
}

proptest! {
    /// Non-convex operand: an extruded L-prism cut by a grid box still satisfies
    /// `vol(L) = vol(L∩box) + vol(L∖box)` — the reflex corner and holed footprint don't break the
    /// difference/intersection split. `vol(L)` is the independent profile-area × depth.
    #[test]
    fn l_prism_cut_identity(
        w in 2i64..6, h in 2i64..6, d in 1i64..4,
        nxr in 0i64..5, nyr in 0i64..5,
        blo in prop::array::uniform3(0i64..6), bext in prop::array::uniform3(1i64..6),
    ) {
        // Notch strictly inside the rectangle — derived (not `prop_assume`d) so no samples reject.
        let nx = 1 + nxr % (w - 1);
        let ny = 1 + nyr % (h - 1);
        let vl = (w * h - (w - nx) * (h - ny)) as f64 * d as f64;
        let build = || {
            let mut m = Model::new();
            let OpOutput::Extrude { solid: a, .. } =
                apply(&mut m, &extrude_op(l_profile(w, h, nx, ny), d as f64)).unwrap()
            else {
                unreachable!()
            };
            let b = grid_box(&mut m, blo, bext);
            (m, a, b)
        };
        let (mut mc, ac, bc) = build();
        let (mut mk, ak, bk) = build();
        if let (Some(vc), Some(vk)) = (
            run(&mut mc, BoolKind::Common, ac, bc),
            run(&mut mk, BoolKind::Cut, ak, bk),
        ) {
            prop_assert!(approx(vl, vc + vk), "L {vl} vs common {vc} + cut {vk}");
        }
    }
}

/// Idempotence — capability probe first (identical operands are a full coplanar coincidence, not a
/// measure-zero one). If this succeeds, the property holds for any grid box; if the engine cannot
/// do identical operands, this test documents the gap rather than a proptest going permanently red.
#[test]
fn idempotence_of_identical_boxes() {
    let build = || {
        let mut m = Model::new();
        let a = grid_box(&mut m, [0, 0, 0], [2, 3, 4]);
        let b = grid_box(&mut m, [0, 0, 0], [2, 3, 4]);
        (m, a, b)
    };
    let va = box_vol([2, 3, 4]);
    let (mut mf, af, bf) = build();
    let vf = run(&mut mf, BoolKind::Fuse, af, bf).expect("Fuse(A, A) is in coverage");
    assert!(approx(vf, va), "Fuse(A,A) volume {vf} vs {va}");
    let (mut mc, ac, bc) = build();
    let vc = run(&mut mc, BoolKind::Common, ac, bc).expect("Common(A, A) is in coverage");
    assert!(approx(vc, va), "Common(A,A) volume {vc} vs {va}");
}

/// Vacuousness pins (family #3): representative in-coverage grid configs must succeed, so the
/// property net above cannot go quietly green by the engine rejecting the core cases.
#[test]
fn core_configs_succeed() {
    // partial overlap, containment, face-contact stack.
    for (alo, aext, blo, bext) in [
        ([0, 0, 0], [3, 3, 3], [1, 1, 1], [3, 3, 3]), // partial
        ([0, 0, 0], [4, 4, 4], [1, 1, 1], [2, 2, 2]), // B contained in A
        ([0, 0, 0], [2, 2, 2], [0, 0, 2], [2, 2, 2]), // face-contact stack (z)
    ] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let mut m = Model::new();
            let a = grid_box(&mut m, alo, aext);
            let b = grid_box(&mut m, blo, bext);
            run(&mut m, kind, a, b).unwrap_or_else(|| {
                panic!("{kind:?} on a core config must succeed: {alo:?} {blo:?}")
            });
        }
    }
}
