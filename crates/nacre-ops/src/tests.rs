use super::*;

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: crate::SketchPlane) -> crate::SketchFrame {
    match crate::apply(
        m,
        &crate::Operation::DatumPlane {
            def: crate::DatumDef::Stated(plane),
        },
    ) {
        Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}
use crate::tolerant::Judge;
use crate::transform::transform;
use crate::{boolean::*, ops::*, planes::*};
use nacre_cip::WitnessPoint;
use nacre_geom::intersect::{planes_coplanar, three_planes};
use nacre_geom::{Plane, Surface};
use nacre_scalar::Axis;
use nacre_topo::{Loop, Orientation, VertexDef};
use proptest::prelude::*;
use std::collections::HashMap;

/// Test shim: a boolean whose result is exactly one solid. Most tests operate on a single
/// body; this asserts that and returns the lone handle, so call sites read as before while
/// `boolean` itself returns the full `Vec` (cell 0.4 multi-solid).
fn boolean_one(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    let solids = boolean(model, kind, a, b)?;
    assert_eq!(
        solids.len(),
        1,
        "boolean_one: expected one solid, got {}",
        solids.len()
    );
    Ok(solids[0])
}

fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

/// Is there an outer-shell face on the plane through `pt` with normal `n`, oriented that way?
/// The production path names a cap by the *face* that made it (`find_face_coplanar_with`); a
/// test that wants to say "a face sits on z = 1.5 facing +z" has no such face in hand, and
/// asserting geometry from coordinates is exactly what a test may do.
fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
    let Some(target) = Plane::from_point_normal(pt, n) else {
        return false;
    };
    let shell = m.solids.get(solid).outer;
    m.shells.get(shell).faces.iter().any(|&fh| {
        let f = m.faces.get(fh);
        let Surface::Plane(plane) = m.surface(f.surface) else {
            return false;
        };
        let sign = f64::from(f.orientation.sign());
        planes_coplanar(plane, &target) && (plane.normal() * sign).dot(n) > 0.0
    })
}

fn square() -> Profile2d {
    Profile2d::polygon(vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)]).unwrap()
}

/// A world-XY extrude for a log that is **replayed** rather than applied to a live model.
///
/// ★ Its frame names a plane in a *throwaway* model, and that is sound for one reason: a
/// log's handles are index vocabulary, and `replay` re-anchors them onto the model it builds
/// (`docs/design.md` §2). The world planes are seeded at fixed indices, so `world(Axis::Z)`
/// names the same plane in every model. Do **not** `apply` one of these to a live model —
/// that is the cross-model misuse `Store::get`'s debug guard exists to catch.
fn extrude_log_op(profile: Profile2d, dist: f64) -> Operation {
    extrude_op(&Model::new(), profile, dist)
}

fn extrude_op(m: &Model, profile: Profile2d, dist: f64) -> Operation {
    Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
        profile,
        dist,
    }
}

fn regular_ngon(n: usize, r: f64) -> Profile2d {
    let points = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            p2(r * a.cos(), r * a.sin())
        })
        .collect();
    Profile2d::polygon(points).unwrap()
}

#[test]
fn square_extrudes_to_a_cube() {
    let m = replay(&[extrude_log_op(square(), 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertices.len(), 8);
    assert_eq!(m.edges.len(), 12);
    assert_eq!(m.faces.len(), 6);
    assert_eq!(m.solids.len(), 1);

    let mut got: Vec<[f64; 3]> = m
        .vertices
        .iter()
        .map(|(vh, _)| m.vertex_point(vh).as_array())
        .collect();
    let mut want = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let key = |p: &[f64; 3]| (p[0].to_bits(), p[1].to_bits(), p[2].to_bits());
    got.sort_by_key(key);
    want.sort_by_key(key);
    assert_eq!(got, want);
}

#[test]
fn triangle_extrudes_to_a_prism() {
    let tri = Profile2d::polygon(vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(1.0, 1.5)]).unwrap();
    let m = replay(&[extrude_log_op(tri, 3.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertices.len(), 6);
    assert_eq!(m.edges.len(), 9);
    assert_eq!(m.faces.len(), 5);
}

#[test]
fn pentagon_extrudes_clean() {
    let m = replay(&[extrude_log_op(regular_ngon(5, 2.0), 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertices.len(), 10);
    assert_eq!(m.faces.len(), 7);
}

#[test]
fn concave_l_profile_is_valid() {
    // An L-shape (a reflex vertex) — a simple concave hexagon.
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.vertices.len(), 12);
    assert_eq!(m.faces.len(), 8);
}

/// The L-prism: profile `[(0,0),(2,0),(2,1),(1,1),(1,2),(0,2)]` extruded to
/// z ∈ [0,1]. Material = bottom bar (x∈[0,2],y∈[0,1]) ∪ left bar (x∈[0,1],
/// y∈[1,2]); the notch (x∈[1,2],y∈[1,2]) is empty.
fn l_prism() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    let s = m.live_solids[0];
    (m, s)
}

/// The same L-prism, its profile started one vertex earlier so the reflex corner
/// `(1,1)` lands at index 1 of the cap's loop. Geometrically identical.
fn rotated_l_prism() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
        p2(0.0, 0.0),
        p2(2.0, 0.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    let s = m.live_solids[0];
    (m, s)
}

/// `FaceInfo::n_out` is documented as the single source of "outward". Two
/// independent sources say which way that is: the ring, which winds CCW about the
/// outward normal, and the b-rep's own `Surface` plus `Orientation`. They must agree
/// on every face of every solid.
///
/// `n_out` reads the second source now (the stored orientation), so what this
/// pins is the first: the witness triangle `outer_tri` picks must span the ring's
/// winding. `outer_tri` used to read the turn at the first non-collinear corner,
/// which is the winding **only when that corner is convex**. The four fixtures
/// below were safe by accident — none starts its cap loop one vertex before a
/// reflex corner. `rotated_l_prism` does, and it is the same solid.
#[test]
fn outward_normals_agree_with_their_orientation() {
    let mut cube = Model::new();
    let c = cube.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let (ml, sl) = l_prism();
    let (mu, su) = u_prism();
    let (mr, sr) = rotated_l_prism();
    for (name, m, s) in [
        ("cube", &cube, c),
        ("l_prism", &ml, sl),
        ("u_prism", &mu, su),
        ("rotated_l_prism", &mr, sr),
    ] {
        for pi in &collect_planes(m, s).unwrap() {
            let pi = pi.plane();
            let cos = (pi.tri[1] - pi.tri[0])
                .cross(pi.tri[2] - pi.tri[0])
                .normalize()
                .expect("a widest corner spans area")
                .dot(pi.n_out);
            assert!(
                cos > 0.5,
                "{name}: the witness triangle does not span its face's stated outward"
            );
        }
    }
}

/// `pt3_base_collinear` is exact on the pre-rotation rational bases: three genuinely
/// collinear points stay collinear under rotation (→ skipped), and a real sliver (one point
/// off the line) is never falsely called collinear (→ its crossing is kept, no silent-wrong).
#[test]
fn pt3_base_collinear_exact() {
    use nacre_scalar::{Angle, Axis, Rat};
    let ang = Angle::from_deg(Rat::from_int(37)).unwrap();
    let piv = [Rat::from_int(2), Rat::from_int(-1), Rat::from_int(0)];
    let rp = |x: i128, y: i128, z: i128| {
        WitnessPoint::at([Rat::from_int(x), Rat::from_int(y), Rat::from_int(z)]).rotate_about(
            Axis::Z,
            ang,
            piv,
        )
    };
    // (0,0,0), (2,4,6), (1,2,3): all on the line t·(1,2,3) → collinear.
    assert!(pt3_base_collinear(&rp(0, 0, 0), &rp(2, 4, 6), &rp(1, 2, 3)));
    // (1,2,4) is off that line (z), a real nonzero-area triangle → not collinear.
    assert!(!pt3_base_collinear(
        &rp(0, 0, 0),
        &rp(2, 4, 6),
        &rp(1, 2, 4)
    ));
}

/// The L-prism with a `[0.1,0.9]³` box strictly inside its bottom bar
/// (non-coplanar coordinates ⇒ no shared face planes). `V_L = 3`, `V_box =
/// 0.512`.
fn l_and_inner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
    (m, l, bx)
}

/// The L-prism with a box biting its convex corner `(2, 0)` — the first
/// non-convex *overlap* (a real single-chord seam), M5-d2. The box spans
/// `x∈[1.3,2.4]`, `y∈[-0.3,0.4]`, `z∈[0.2,1.4]`: it straddles the corner in x
/// and y, and its z-range pokes above the L (`z=1`) while its floor `z=0.2`
/// sits inside — so every crossing edge is a clean straddle (no edge tunnels
/// fully through the other) and no box face is coplanar with an L face. The
/// span is deliberately asymmetric so no seam point lands on a face centre
/// (where both fan diagonals cross and every apex would graze).
/// Overlap = `x∈[1.3,2]·y∈[0,0.4]·z∈[0.2,1]` = `0.224`;
/// `V_L=3`, `V_box=1.1·0.7·1.2=0.924`.
fn l_and_corner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(
        Point3::from_array([1.3, -0.3, 0.2]),
        Point3::from_array([2.4, 0.4, 1.4]),
    );
    (m, l, bx)
}

/// The L with a box straddling its reflex corner (1,1): a *single* chord with one
/// reflex bend. The box vertex `(1.6,1.6,·)` sits in the L's notch — inside the
/// convex hull, outside the L — exactly where a convex half-space test would
/// misclassify it `Inside`. Overlap = `xy(1.0 − notch 0.36) · z(0.8)` = `0.512`.
fn l_and_reflex_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(
        Point3::from_array([0.6, 0.6, 0.2]),
        Point3::from_array([1.6, 1.6, 1.4]),
    );
    (m, l, bx)
}

// --- Rotated booleans go live (overhaul 3d-i, `ROTATED_UNSUPPORTED` retired) ---
// A boolean commutes with a rigid motion, so rotating both operands by the same
// irrational-angle isometry must give the rigid image of the unrotated result — identical
// volume, solid count, and cavity count, and still valid. These are the first live proof
// that the CIP-wired machinery (arrangement, seam, in/out, outer/cavity — 3a–3c-vi) is
// sound end-to-end on rotated (rounded-irrational) geometry.

/// A cutter's convex corner landing exactly on the target's concave corner — three planes
/// (x=1,y=1,z=1) meeting at one point (1,1,1), six faces there — makes a **non-manifold pinch**
/// (two face-fans touching at the point). No valid 2-manifold solid exists, so `boolean` rejects
/// with the clear `NON_MANIFOLD_VERTEX` reason (not the incidental `EULER_PARITY`), leaving the
/// live model untouched. **Rotation-independent**: both the axis-aligned and the rotated framings
/// (exact and CIP-kernel paths) hit the same pinch and reject. R = a cube minus a far-corner
/// octant; C = a cube whose +corner is that removed octant's inner corner. (Coplanar contact away
/// from a corner works — [`a_rotated_boolean_result_can_be_cut_again`].)
#[test]
fn a_corner_coincident_cut_is_rejected_not_silently_wrong() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    // `rotate`: None = axis-aligned (exact path); Some = the result and cutter tilted (CIP path).
    let run = |rotate: bool| {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        // C's +corner is (1,1,1) = R's concave corner → three shared planes meet there.
        let c = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([1.0; 3]),
        );
        m.rebuild_adjacency();
        let (r, c) = if rotate {
            let r = transform(&mut m, r, &iso).unwrap();
            m.rebuild_adjacency();
            let c = transform(&mut m, c, &iso).unwrap();
            m.rebuild_adjacency();
            (r, c)
        } else {
            (r, c)
        };
        let live = m.live_solids.clone();
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, r, c),
            RejectReason::NonManifoldVertex,
        );
        assert_eq!(m.live_solids, live, "reject must not mutate the live set");
    };
    run(false); // axis-aligned
    run(true); // rotated
}

/// Two cubes touching only at the corner (1,1,1): nothing is joined, so the Fuse comes back as
/// the two bodies it was handed. ★ This used to be `NON_MANIFOLD_VERTEX` — true of the single
/// welded body the reconstruction built then, and beside the point once the pieces are minted
/// per solid. The pinch reject is still there for a body that touches *itself*
/// (`tests/contact_separates.rs`).
#[test]
fn two_cubes_touching_at_a_corner_fuse_to_two_bodies() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("a point contact separates");
    assert_eq!(out.len(), 2);
    m.rebuild_adjacency();
    assert_eq!(nacre_validate::validate(&m), Vec::new());
}

/// The kernel's validator catches the pinch too (it previously only exposed it indirectly as
/// `EulerParity`). Bypass `boolean`'s reject via `arrangement::boolean` to obtain a malformed
/// solid, then confirm `validate` reports a `NonManifoldVertex`.
///
/// ★ The shape has to be one that **cannot** part: two bodies touching at a corner now come
/// back as two bodies and there is no malformed solid to obtain from them. So A and B meet only
/// at `(2,2,1)` while two bridges run around the contact and join them elsewhere — the material
/// loops, the pinch is real, and the vertex check (which lives in `boolean`, not in the
/// arrangement) is the one this bypasses.
#[test]
fn validate_reports_the_non_manifold_pinch() {
    let mut m = Model::new();
    let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
        let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        s
    };
    let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let g1 = cub(&mut m, [1.0, 0.3, 0.2], [4.5, 1.3, 0.8]);
    let g2 = cub(&mut m, [3.5, 0.3, 0.2], [4.5, 3.8, 1.6]);
    let b = cub(&mut m, [2.0, 2.0, 1.0], [4.0, 4.0, 2.0]);
    let t1 = boolean(&mut m, BoolKind::Fuse, a, g1).expect("a and g1 overlap");
    m.rebuild_adjacency();
    let t2 = boolean(&mut m, BoolKind::Fuse, t1[0], g2).expect("g1 and g2 overlap");
    m.rebuild_adjacency();
    crate::arrangement::boolean(&mut m, BoolKind::Fuse, t2[0], b).unwrap();

    m.rebuild_adjacency();
    let issues = nacre_validate::validate(&m);
    assert!(
        issues
            .iter()
            .any(|v| matches!(v, nacre_validate::Violation::NonManifoldVertex { .. })),
        "validate must flag the pinch: {issues:?}"
    );
}

/// Predicates over a rotated result's witness planes are rotation-invariant against the same
/// result unrotated — the provenance witness (with outward winding) defines the exact plane,
/// so `Judge::orient3d` agrees on every definite triple. A regression guard on the witness itself.
#[test]
fn rotated_result_witness_predicates_are_invariant() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let table = |m: &Model, s: Handle<Solid>| {
        let f = collect_planes(m, s).unwrap();
        let c = plane_classes(&crate::planes::test_judge(&f));
        dense_planes(&f, &c).0
    };
    let pu = table(&m, r);
    let r2 = transform(&mut m, r, &iso).unwrap();
    m.rebuild_adjacency();
    let pr = table(&m, r2);
    assert_eq!(pu.len(), pr.len(), "rotation preserves the plane count");
    let n = pu.len();
    let indep = |p: &[WorkingPlane], a: usize, b: usize, c: usize| {
        let nrm = |k: usize| p[k].plane.normal();
        nrm(a).dot(nrm(b).cross(nrm(c))).abs() > 0.3
    };
    let mut disagree = 0;
    for p in 0..n {
        for q in (p + 1)..n {
            for rr in (q + 1)..n {
                if !indep(&pu, p, q, rr) {
                    continue;
                }
                for j in 0..n {
                    if j == p || j == q || j == rr {
                        continue;
                    }
                    let su = crate::planes::test_judge(&pu).orient3d(p, q, rr, j);
                    let sr = crate::planes::test_judge(&pr).orient3d(p, q, rr, j);
                    if su != 0 && sr != 0 && su != sr {
                        eprintln!("DISAGREE orient3d ({p},{q},{rr},{j}): u={su} r={sr}");
                        disagree += 1;
                    }
                }
            }
        }
    }
    assert_eq!(
        disagree, 0,
        "{disagree} predicate disagreements (witness wrong)"
    );
}

/// Rotating a boolean *result* and feeding it back into a boolean (was `ROTATED_UNSUPPORTED`):
/// `collect_planes` now witnesses each rotated seam face's plane through provenance — its
/// plane is `R(π)` for an operand plane `π`, recovered from the operand face still on `π` and
/// rotated by the face's own chain. A boolean commutes with a rigid motion, so the rotated
/// chain's result matches the unrotated chain's (volume, solid count) and stays valid.
#[test]
fn a_rotated_boolean_result_can_be_cut_again() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    // Chain: R = Cut(A, B) removes a far-corner octant; then Cut(R, C) removes a near one.
    let build = |m: &mut Model| -> (Handle<Solid>, Handle<Solid>) {
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let r = boolean_one(m, BoolKind::Cut, a, b).unwrap();
        // A clean slab that severs R at x = 0.5 — no plane of C coincides with any of R's
        // (avoids the separate rotated-coplanar-contact gap; isolates the witness).
        let c = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([0.5, 4.0, 4.0]),
        );
        (r, c)
    };
    // Unrotated reference (reuse already works when nothing is rotated).
    let mut m0 = Model::new();
    let (r0, c0) = build(&mut m0);
    m0.rebuild_adjacency();
    let ref_out = boolean(&mut m0, BoolKind::Cut, r0, c0).unwrap();
    m0.rebuild_adjacency();
    let ref_vol: f64 = ref_out
        .iter()
        .map(|&s| nacre_props::mass_props(&m0, s).unwrap().volume)
        .sum();
    // Rotated: turn the *result* R (and C) by the same isometry, then reuse R.
    let mut m = Model::new();
    let (r, c) = build(&mut m);
    m.rebuild_adjacency();
    let r = transform(&mut m, r, &iso).unwrap();
    m.rebuild_adjacency();
    let c = transform(&mut m, c, &iso).unwrap();
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Cut, r, c).unwrap();
    m.rebuild_adjacency();
    let issues = nacre_validate::validate(&m);
    assert!(
        issues.is_empty(),
        "rotated-result reuse must be valid: {issues:?}"
    );
    let vol: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert_eq!(
        out.len(),
        ref_out.len(),
        "solid count invariant under rotation"
    );
    assert!(
        (vol - ref_vol).abs() < 1e-6,
        "rotated reuse volume {vol} vs unrotated {ref_vol}"
    );
}

/// Result-reuse rotation stress: build R with a first boolean, then feed R into a second
/// boolean with a fresh cutter C — once unrotated, once with R and C rotated by the same
/// isometry. A boolean commutes with a rigid motion, so the rotated reuse must equal the
/// unrotated one (volume, solid count, cavity count) or be an honest reject — never silently
/// wrong. This is the invariant on the newly-enabled rotated-*Discovered* geometry (every
/// vertex of a boolean result is `Discovered`, so its rotation exercises the provenance
/// witness on every face). `#[ignore]`: rotated booleans escalate to astro-float (~1–3 s).
#[test]
#[ignore = "slow: rotated result-reuse booleans (run with --ignored)"]
fn rotated_result_reuse_stress() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    #[derive(PartialEq, Debug)]
    enum Out {
        Rej,
        Ok(f64, usize, usize),
    }
    let rot = |axis: Axis, deg: i128, piv: [i128; 3]| {
        Isometry::rotation(Rotation {
            axis,
            point: [
                Rat::from_int(piv[0]),
                Rat::from_int(piv[1]),
                Rat::from_int(piv[2]),
            ],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    };
    // (first kind, second kind, |m| -> (a, b, c)). R = kind1(a, b); out = kind2(R, c).
    type Build = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>, Handle<Solid>)>;
    let cuboid = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
        m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
    };
    let fixtures: Vec<(&str, BoolKind, BoolKind, Build)> = vec![
        (
            "corner_then_slab",
            BoolKind::Cut,
            BoolKind::Cut,
            Box::new(move |m: &mut Model| {
                let a = cuboid(m, [0.0; 3], [2.0; 3]);
                let b = cuboid(m, [1.0; 3], [3.0; 3]);
                let c = cuboid(m, [-1.0, -1.0, -1.0], [0.5, 4.0, 4.0]);
                (a, b, c)
            }),
        ),
        (
            "fuse_then_bite",
            BoolKind::Fuse,
            BoolKind::Cut,
            Box::new(move |m: &mut Model| {
                let a = cuboid(m, [0.0; 3], [2.0, 1.0, 1.0]);
                let b = cuboid(m, [0.0, 0.0, 0.0], [1.0, 2.0, 1.0]);
                let c = cuboid(m, [1.3, 1.3, -1.0], [3.0, 3.0, 2.0]);
                (a, b, c)
            }),
        ),
        (
            "cut_then_fuse",
            BoolKind::Cut,
            BoolKind::Fuse,
            Box::new(move |m: &mut Model| {
                let a = cuboid(m, [0.0; 3], [2.0; 3]);
                let b = cuboid(m, [1.3, 1.3, 1.3], [3.0, 3.0, 3.0]);
                let c = cuboid(m, [-0.7, 0.4, 0.4], [0.3, 1.4, 1.4]);
                (a, b, c)
            }),
        ),
    ];
    let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
        ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
        ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
        (
            "Z50>Y37",
            vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::Y, 37, [0, 0, 1])],
        ),
    ];
    let run = |k1: BoolKind, k2: BoolKind, build: &Build, isos: &[Isometry]| -> Out {
        let mut m = Model::new();
        let (a, b, c) = build(&mut m);
        m.rebuild_adjacency();
        let Ok(r) = boolean(&mut m, k1, a, b) else {
            return Out::Rej;
        };
        assert_eq!(r.len(), 1, "first boolean is a single solid");
        let mut r = r[0];
        let mut c = c;
        m.rebuild_adjacency();
        for iso in isos {
            r = transform(&mut m, r, iso).unwrap();
            m.rebuild_adjacency();
            c = transform(&mut m, c, iso).unwrap();
            m.rebuild_adjacency();
        }
        match boolean(&mut m, k2, r, c) {
            Ok(solids) => {
                m.rebuild_adjacency();
                assert!(
                    nacre_validate::validate(&m).is_empty(),
                    "INVALID rotated reuse result"
                );
                let vol: f64 = solids
                    .iter()
                    .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                    .sum();
                let cav: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
                Out::Ok(vol, solids.len(), cav)
            }
            Err(_) => Out::Rej,
        }
    };
    let vclose =
        |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());
    let (mut success, mut reject, mut silent) = (0, 0, 0);
    for (fname, k1, k2, build) in &fixtures {
        let base = run(*k1, *k2, build, &[]);
        for (rname, isos) in &isos_list {
            let r = run(*k1, *k2, build, isos);
            match (&base, &r) {
                (_, Out::Rej) => reject += 1,
                (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2))
                    if vclose(*v1, *v2) && s1 == s2 && c1 == c2 =>
                {
                    success += 1
                }
                _ => {
                    silent += 1;
                    eprintln!("SILENT-WRONG {fname} {rname} base={base:?} rot={r:?}");
                }
            }
        }
    }
    eprintln!("REUSE STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent}");
    assert_eq!(
        silent, 0,
        "a rotated result-reuse boolean was silently wrong"
    );
    assert!(
        success >= 1,
        "at least one rotated reuse must actually succeed"
    );
}

/// Adversarial rotation stress (overhaul 3d-iii): many fixtures × kinds × rotations
/// (single-axis, and Euler chains reaching arbitrary orientation) confirm the DNA
/// invariant — a rotated boolean is *never silently wrong*: its result either equals the
/// unrotated one (a boolean commutes with a rigid motion, so volume/solid-count/cavity-count
/// are invariant) or is an honest reject. `#[ignore]`: each rotated boolean escalates its
/// CIP predicates to astro-float and costs ~0.5–2.5 s, so this runs on demand, not per commit
/// (the invariance regression guard is the fast `rotated_*` tests above).
#[test]
#[ignore = "slow: rotated booleans ~2s each (run with --ignored)"]
fn rotation_invariance_stress() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    #[derive(PartialEq, Debug)]
    enum Out {
        Rej,
        Ok(f64, usize, usize),
    }
    let rot = |axis: Axis, deg: i128, piv: [i128; 3]| {
        Isometry::rotation(Rotation {
            axis,
            point: [
                Rat::from_int(piv[0]),
                Rat::from_int(piv[1]),
                Rat::from_int(piv[2]),
            ],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    };
    let run = |build: &dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>),
               kind: BoolKind,
               isos: &[Isometry]|
     -> Out {
        let (mut m, mut a, mut b) = build();
        for iso in isos {
            a = transform(&mut m, a, iso).unwrap();
            m.rebuild_adjacency();
            b = transform(&mut m, b, iso).unwrap();
            m.rebuild_adjacency();
        }
        match boolean(&mut m, kind, a, b) {
            Ok(solids) => {
                m.rebuild_adjacency();
                assert!(
                    nacre_validate::validate(&m).is_empty(),
                    "INVALID rotated result"
                );
                let vol: f64 = solids
                    .iter()
                    .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                    .sum();
                let cav: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
                Out::Ok(vol, solids.len(), cav)
            }
            Err(_) => Out::Rej,
        }
    };
    let vclose =
        |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());
    let matches = |base: &Out, r: &Out| match (base, r) {
        (Out::Rej, Out::Rej) => true,
        (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2)) => vclose(*v1, *v2) && s1 == s2 && c1 == c2,
        _ => false,
    };
    let cube = |lo: [f64; 3], hi: [f64; 3]| {
        move || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
            (m, a, b)
        }
    };
    type Build = Box<dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>)>;
    let fixtures: Vec<(&str, Build)> = vec![
        ("corner", Box::new(l_and_corner_box)),
        (
            "rod",
            Box::new(|| {
                let (m, l, r) = l_and_rod();
                (m, r, l) // sever = Cut(rod, l): a=rod, b=l
            }),
        ),
        ("inner", Box::new(l_and_inner_box)),
        ("cube_corner", Box::new(cube([0.5; 3], [1.5; 3]))),
    ];
    let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
        ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
        ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
        (
            "Z50>X37",
            vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::X, 37, [0, 0, 1])],
        ),
        // A three-axis Euler chain reaches an arbitrary orientation (axes are X/Y/Z only).
        (
            "Z30>X30>Y73",
            vec![
                rot(Axis::Z, 30, [1, 1, 0]),
                rot(Axis::X, 30, [0, 0, 0]),
                rot(Axis::Y, 73, [0, 2, 0]),
            ],
        ),
    ];
    let (mut success, mut reject, mut silent, mut skipped) = (0, 0, 0, 0);
    for (fname, build) in &fixtures {
        for kind in [BoolKind::Cut, BoolKind::Fuse, BoolKind::Common] {
            let base = run(build.as_ref(), kind, &[]);
            for (rname, isos) in &isos_list {
                if matches!(base, Out::Rej) {
                    skipped += 1;
                    continue;
                }
                let r = run(build.as_ref(), kind, isos);
                if matches!(r, Out::Rej) {
                    reject += 1; // an honest reject is acceptable, not a counterexample
                } else if matches(&base, &r) {
                    success += 1;
                } else {
                    silent += 1;
                    eprintln!("SILENT-WRONG {fname} {kind:?} {rname} base={base:?} rot={r:?}");
                }
            }
        }
    }
    eprintln!(
        "ROTATION STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent} skipped(base_rej)={skipped}"
    );
    assert_eq!(
        silent, 0,
        "a rotated boolean was silently wrong (valid but != unrotated)"
    );
}

/// A signature that changes if `Store::push` order (hence handle identity) changes:
/// vertex points in handle order — exactly what `assemble_fuse_cut` assigns by first
/// appearance across `faces` — plus edge/face/solid counts and sorted volumes.
#[cfg(feature = "parallel")]
fn model_sig(m: &Model, solids: &[Handle<Solid>]) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    let _ = write!(s, "S{}", solids.len());
    for (vh, _) in m.vertices.iter() {
        let p = m.vertex_point(vh).as_array();
        let _ = write!(
            s,
            "|{:x},{:x},{:x}",
            p[0].to_bits(),
            p[1].to_bits(),
            p[2].to_bits()
        );
    }
    let _ = write!(s, "|E{}F{}", m.edges.len(), m.faces.len());
    let mut vols: Vec<u64> = solids
        .iter()
        .map(|&sh| nacre_props::mass_props(m, sh).unwrap().volume.to_bits())
        .collect();
    vols.sort_unstable();
    let _ = write!(s, "|V{vols:?}");
    s
}

/// The parallel boolean must be bit-identical regardless of rayon thread count — replay
/// determinism (DNA) requires thread-order independence. Each fixture's result under a
/// 1-thread pool must equal the default many-thread result, run repeatedly so scheduling
/// jitter would show.
///
/// **★ Two things this test has to keep honest about itself.**
///
/// First, it passed for the whole stretch when there was *no* parallelism — the cutover
/// took the old engine's `par_iter` calls with it, and a one-thread pool is trivially
/// equal to a many-thread one when nothing forks. So the fixtures must have real
/// parallel width: the fin fold below reaches into the dozens of plane classes, where
/// the two-ngon fuse has barely a dozen.
///
/// Second, `model_sig` compares coordinates, counts and volumes — **not the report**.
/// `Notes::sorted` is a stable sort by site, so any two entries sharing a site keep
/// their *arrival* order, and under `parallel` that is the schedule's to decide. So the
/// signature includes `BoolReport`, which is the only thing that would catch it.
#[cfg(feature = "parallel")]
#[test]
fn parallel_boolean_is_thread_order_independent() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    // A hub with fins arrayed around it — every fin is turned by an angle with no exact
    // f64, so every judgement is on the toleranced path and the report is non-empty.
    let fin_fold = |n: i128| -> String {
        let mut m = replay(&[extrude_log_op(
            Profile2d::polygon(vec![
                p2(-3.0, -3.0),
                p2(3.0, -3.0),
                p2(3.0, 3.0),
                p2(-3.0, 3.0),
            ])
            .unwrap(),
            2.0,
        )])
        .unwrap();
        let mut acc = m.live_solids[0];
        let mut sig = String::new();
        for i in 0..n {
            let fin = {
                let __w7 = SketchFrame::world(&m, Axis::Z);
                let out = ops::apply(
                    &mut m,
                    &Operation::Extrude {
                        frame: __w7,
                        profile: Profile2d::polygon(vec![
                            p2(2.0, -0.4),
                            p2(8.0, -0.4),
                            p2(8.0, 0.4),
                            p2(2.0, 0.4),
                        ])
                        .unwrap(),
                        dist: 1.0,
                    },
                )
                .unwrap();
                m.rebuild_adjacency();
                match out {
                    OpOutput::Extrude { solid, .. } => solid,
                    o => panic!("{o:?}"),
                }
            };
            let fin = transform(
                &mut m,
                fin,
                &Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::new(360 * i, n).unwrap()).unwrap(),
                }),
            )
            .unwrap();
            m.rebuild_adjacency();
            let (solids, report) = boolean_with_report(&mut m, BoolKind::Fuse, acc, fin).unwrap();
            m.rebuild_adjacency();
            acc = solids[0];
            // **The report travels in the signature**, not just the geometry. Measured on
            // this fixture: up to 841 coincidences per boolean, so `loosest`'s `max_by_key`
            // is choosing among hundreds of candidates — which is the tie-break that a
            // schedule could otherwise decide.
            assert!(
                i == 0 || report.coincidences > 0,
                "the report is empty, so comparing it proves nothing"
            );
            sig.push_str(&format!("{report:?}"));
        }
        sig.push_str(&model_sig(&m, &[acc]));
        sig
    };

    let build = || {
        let ngon = |n: usize, r: f64, cx: f64, cy: f64| {
            Profile2d::polygon(
                (0..n)
                    .map(|i| {
                        let ang = std::f64::consts::TAU * (i as f64) / (n as f64);
                        p2(cx + r * ang.cos(), cy + r * ang.sin())
                    })
                    .collect(),
            )
            .unwrap()
        };
        let mut m = replay(&[
            extrude_log_op(ngon(16, 2.0, 0.0, 0.0), 3.0),
            extrude_log_op(ngon(16, 2.0, 2.5, 0.5), 3.0),
        ])
        .unwrap();
        let a = m.live_solids[0];
        let b = m.live_solids[1];
        let up = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
        let b = transform(&mut m, b, &up).unwrap();
        m.rebuild_adjacency();
        let tilt = rot_iso(Axis::X, 30);
        let a = transform(&mut m, a, &tilt).unwrap();
        m.rebuild_adjacency();
        let b = transform(&mut m, b, &tilt).unwrap();
        m.rebuild_adjacency();
        (m, a, b)
    };
    let two_ngons = || {
        let (mut m, a, b) = build();
        let (solids, report) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        format!("{report:?}{}", model_sig(&m, &solids))
    };
    let pool1 = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let seven_fins = || fin_fold(7);
    // **The only fixture here whose alias table is not empty.** `Aliases::record` returns
    // immediately below four planes, so every other model leaves the table at zero and
    // never exercises the snapshot-and-absorb a parallel round is built on. Measured on
    // this one: 36 aliases, settled over two rounds.
    let four_plane = || {
        let (mut m, target, bar) = four_plane_model(0.2, 45);
        let (solids, report) = boolean_with_report(&mut m, BoolKind::Cut, target, bar).unwrap();
        m.rebuild_adjacency();
        format!("{report:?}{}", model_sig(&m, &solids))
    };
    for (name, run) in [
        (
            "two rotated ngons",
            &two_ngons as &(dyn Fn() -> String + Sync),
        ),
        (
            "a seven-fin fold",
            &seven_fins as &(dyn Fn() -> String + Sync),
        ),
        (
            "a four-plane concurrency",
            &four_plane as &(dyn Fn() -> String + Sync),
        ),
    ] {
        let reference = pool1.install(run);
        // A signature that came back empty would make this vacuous whatever the schedule.
        assert!(!reference.is_empty(), "{name}: nothing to compare");
        for _ in 0..4 {
            assert_eq!(
                run(),
                reference,
                "{name}: the result depends on thread order"
            );
        }
    }
}

/// A U-prism: a bottom bar `y∈[0,1]` with two prongs rising from it. The prong
/// tops sit at *different* heights (y=2.3 and y=2.0) on purpose — level tops
/// would be coplanar faces, which the pre-cutover `has_coplanar_pair` door guard
/// rejected before the seam machinery ran. That guard is gone; the staggering stays
/// as this fixture's pinned shape. Area 3 + 1 + 1.3, extruded 1.0 ⇒ volume 5.3.
fn u_prism() -> (Model, Handle<Solid>) {
    let u = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(3.0, 0.0),
        p2(3.0, 2.3),
        p2(2.0, 2.3),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(u, 1.0)]).unwrap();
    let s = m.live_solids[0];
    (m, s)
}

/// The U with a slab shearing off both prong tops. The slab overhangs the U in
/// x and z, so **every slab edge lies outside the U** (pierces nothing) and every
/// U edge either straddles cleanly (one crossing of the slab's `y=1.5` face) or
/// misses. No edge threads the other solid, so the arcs are the whole story.
fn u_and_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, u) = u_prism();
    let slab = m.add_cuboid(
        Point3::from_array([-0.5, 1.5, -0.5]),
        Point3::from_array([3.5, 2.5, 1.5]),
    );
    (m, u, slab)
}

#[test]
fn u_prism_is_valid() {
    // Pin the fixture itself: a mistyped profile could still trip `multichord`
    // below, for the wrong reason.
    let (m, u) = u_prism();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, u).unwrap().volume;
    assert!((vol - 5.3).abs() < 1e-9, "volume {vol}");
}

// ---- arrangement: seam segment gathering (M5-d3 cell 3b) ----

fn near(a: Point3, b: [f64; 3]) -> bool {
    (a - Point3::from_array(b)).norm() < 1e-9
}

proptest! {}

/// A big cube whose `y=0, z=0` edge is crossed **twice** by the seam: the notch
/// spans `x∈[3,7]` and hangs below both `y=0` and `z=0`, so that one edge enters
/// and leaves it. Extents are asymmetric so no crossing lands on a face centre or a
/// fan diagonal — a debt to the fan, which cell (5a) deleted. Kept: one variable at
/// a time.
///
/// `edge_seam` maps an edge to *one* seam triple. This input is what that map
/// cannot represent — and `boolean` rejects it today (see
/// `an_edge_crossed_twice_is_rejected`).
fn cube_and_notch() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let y = m.add_cuboid(
        Point3::from_array([3.0, -1.0, -1.0]),
        Point3::from_array([7.0, 1.4, 1.2]),
    );
    (m, a, y)
}

/// A thin rod skewering the L's bottom bar in `z`, both ends outside. Each of its four
/// vertical edges pierces the L's two caps, so the caps take a closed seam loop and the
/// rod's walls take two chords apiece — and each wall's vertical edges are crossed
/// **twice**, leaving runs with no vertex at all.
fn l_and_rod() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let rod = m.add_cuboid(
        Point3::from_array([0.3, 0.3, -0.5]),
        Point3::from_array([0.5, 0.6, 1.5]),
    );
    (m, l, rod)
}

/// A slab over the pocketed cube, its underside at height `z0`. The rectangle is
/// asymmetric so that the cube's four vertical edges, which pierce the underside at
/// `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`, miss its fan diagonals; a square slab has all
/// four apexes degenerate at once.
///
/// `z0 = 0.3` runs below the pocket floor, so the slab's underside meets only the
/// cube's outer walls. `z0 = 0.7` runs between the floor and the lid and meets the
/// pocket walls as well — two loops, nested.
fn pocket_and_slab(z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, pc) = pocketed_cube();
    let slab = m.add_cuboid(
        Point3::from_array([-0.2, -0.25, z0]),
        Point3::from_array([1.3, 1.2, 1.5]),
    );
    (m, slab, pc)
}

/// Nested loops, at last (cell 3f-7). Cell 3f-3 argued nesting was reachable and reached
/// for a polyhedral torus; a pocket sliced between its floor and its lid is enough. The
/// slab's underside carries the cube's cross-section as one loop and the pocket's inside
/// it — a loop within a loop, which the pairwise guard over-rejected.
///
/// `Cut` opens, both orders: `A ∩ B = B ∩ {z ≥ 0.7}` = `0.30 − 0.048 = 0.252`, slab `1.74`,
/// pocketed cube `0.92`, so `Cut(slab,pc) = 1.488` and `Cut(pc,slab) = 0.668`. `Cut(pc,slab)`
/// is the one that makes the cube cross-section an *island with a hole* (the slab is B, its
/// underside dropped, so no region survives to own the loops); `Cut(slab,pc)` keeps the slab
/// region and hangs the pocket loop in it as a hole beside an island.
///
/// `Fuse` seals the pocket (`[0.3,0.7]² × [0.5,0.7]`, capped by the slab at `z = 0.7`) into
/// an enclosed cavity — a second shell. Cell (5c) assembles it: the union material is
/// `2.408` and the result carries one cavity of volume `0.032`, two shells, `validate`
/// clean (the void's inward orientation). Both orders (Fuse is commutative).
#[test]
fn a_slab_between_the_lid_and_the_floor_nests_two_loops() {
    for (swap, expect) in [(false, 1.488), (true, 0.668)] {
        let (mut m, slab, pc) = pocket_and_slab(0.7);
        let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
        let r = boolean_one(&mut m, BoolKind::Cut, x, y).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "Cut swap={swap}: {vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - expect).abs() < 1e-9,
            "Cut swap={swap}: {} vs {expect}",
            props.volume
        );
    }
    // Fuse seals the pocket into a cavity — a second shell, assembled by cell (5c).
    for swap in [false, true] {
        let (mut m, slab, pc) = pocket_and_slab(0.7);
        let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
        let r = boolean_one(&mut m, BoolKind::Fuse, x, y).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "Fuse swap={swap}: {vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 2.408).abs() < 1e-9,
            "Fuse swap={swap}: {}",
            props.volume
        );
        assert_eq!(m.solids.get(r).cavities.len(), 1, "Fuse swap={swap}");
        assert_eq!(m.reachable().shells.len(), 2, "Fuse swap={swap}");
    }
}

/// The L with a box biting its reflex corner and poking out the top. The box top
/// (z=1.2) clears the L's z=1 **deliberately**: sunk inside the L's slab, the L's
/// vertical edges at (2,1) and (1,1) would pierce the box's bottom *and* top face,
/// which `pierced_multi` used to reject. Cell 3e-3 supports it; the fixture keeps its
/// clearance so that it goes on testing one thing.
fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.2]),
        Point3::from_array([2.5, 1.5, 1.2]),
    );
    (m, l, bx)
}

/// The L-prism and an L-shaped bar lying in its notch, biting two convex corners of
/// the L's top face. The bar spans `z ∈ [0.5, 1.5]`, so its body clears the cap.
///
/// Each bite crosses **two different** edges of the cap, which is exactly why no edge
/// is pierced twice — the bar takes corners, not edges. Two chords, no closed loop.
fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bar = Profile2d::polygon(vec![
        p2(1.8, 0.8),
        p2(2.1, 0.8),
        p2(2.1, 2.1),
        p2(0.8, 2.1),
        p2(0.8, 1.8),
        p2(1.8, 1.8),
    ])
    .unwrap();
    let __w6 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: b, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w6,
            profile: bar,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    (m, l, b)
}

/// The L-prism with an **L-shaped** stub standing wholly inside its top face,
/// `z ∈ [0.5, 1.5]`. The seam on the cap is a closed loop with a reflex node — the
/// suite's first non-convex inner loop, and the shape a winding must be read from.
///
/// Its coordinates dodge the cap's fan diagonals from `(0,0)` (`y = x`, `y = x/2`,
/// `y = 2x`), which `segment_crosses_face` would graze.
fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let ell = Profile2d::polygon(vec![
        p2(0.2, 0.25),
        p2(0.85, 0.25),
        p2(0.85, 0.4),
        p2(0.35, 0.4), // reflex
        p2(0.35, 0.9),
        p2(0.2, 0.9),
    ])
    .unwrap();
    let __f103 = datum_frame(
        &mut m,
        SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
    );
    let OpOutput::Extrude { solid: stub, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f103,
            profile: ell,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    (m, l, stub)
}

/// A ring's nodes as coordinates — each triple is three planes, so its point is their meet.
fn ring_points(planes: &[WorkingPlane], ring: &[[usize; 3]]) -> Vec<[f64; 3]> {
    ring.iter()
        .map(|t| {
            three_planes(
                &planes[t[0]].plane,
                &planes[t[1]].plane,
                &planes[t[2]].plane,
            )
            .unwrap()
            .as_array()
        })
        .collect()
}

#[test]
fn a_hole_winds_clockwise_and_an_island_counter_clockwise() {
    // The two rings of a holed face are stored with opposite windings — that is what makes one
    // a hole and the other its outer boundary — and `loop_winding` must read exactly that.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap()
        )
        .unwrap(),
        1
    );
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            p,
            &combinatorics::ring_from_names(p, &hole).unwrap()
        )
        .unwrap(),
        -1
    );

    // Nothing but the ring's direction went into that. Reversing it by hand agrees.
    for (name, ring, want) in [("outer", &outer, -1i8), ("hole", &hole, 1)] {
        let mut reversed = ring.clone();
        reversed.reverse();
        assert_eq!(
            combinatorics::loop_winding(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &reversed).unwrap()
            )
            .unwrap(),
            want,
            "{name} reversed"
        );
    }
}

#[test]
fn a_reflex_node_turns_against_its_ring() {
    // The L cap's outer ring is a hexagon with exactly one reflex corner, at `(1, 1)`. A
    // convex node turns with the ring and the reflex one turns against it, so the turn signs
    // are *not* all equal — which is why the winding cannot be read off an arbitrary node.
    let (planes, p, outer, _) = holed_face_rings("dimple");
    let pts = ring_points(&planes, &outer);
    let reflex = pts
        .iter()
        .position(|q| near(Point3::from_array(*q), [1.0, 1.0, 1.0]))
        .expect("the reflex node");
    assert_eq!(
        combinatorics::turn_at(
            &crate::planes::test_judge(&planes),
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap(),
            reflex
        )
        .unwrap(),
        -1
    );
    let turns: Vec<i8> = (0..outer.len())
        .map(|i| {
            combinatorics::turn_at(
                &crate::planes::test_judge(&planes),
                p,
                &combinatorics::ring_from_names(p, &outer).unwrap(),
                i,
            )
            .unwrap()
        })
        .collect();
    assert_eq!(
        turns.iter().filter(|&&t| t == -1).count(),
        1,
        "one reflex corner: {turns:?}"
    );

    // The node `loop_winding` lands on is the lexicographically least — a hull vertex, where
    // the turn *is* the winding. The test finds it by reading coordinates; `loop_winding`
    // finds it with an exact predicate.
    let lo = (0..pts.len())
        .min_by(|&i, &j| pts[i].partial_cmp(&pts[j]).unwrap())
        .unwrap();
    assert_ne!(lo, reflex);
    assert_eq!(
        combinatorics::turn_at(
            &crate::planes::test_judge(&planes),
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap(),
            lo
        )
        .unwrap(),
        1
    );

    // ★ The teeth. A ring is a cycle, so its winding cannot depend on where the walk began.
    // Start it at the reflex node and a `turn_at(ring[0])` implementation reads the reflex
    // sign — the exact fault `outer_tri` shipped. Measured: without this rotation, such an
    // implementation passes every assertion above.
    let mut rotated = outer.clone();
    rotated.rotate_left(reflex);
    assert_eq!(
        combinatorics::turn_at(
            &crate::planes::test_judge(&planes),
            p,
            &combinatorics::ring_from_names(p, &rotated).unwrap(),
            0
        )
        .unwrap(),
        -1
    );
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            p,
            &combinatorics::ring_from_names(p, &rotated).unwrap()
        )
        .unwrap(),
        1
    );
}

/// The L-prism and a П-shaped staple straddling the L's reflex corner. The profile
/// lives in the **XZ** sketch plane and extrudes along `−y`, so the L's cap (`z = 1`)
/// is *parallel* to the extrusion axis and the staple's section there falls into two
/// pieces: one wholly inside the cap, one wrapping the corner `(1,1)`.
///
/// That parallelism is the whole point. A prism cut by a plane **perpendicular** to
/// its axis meets a face in the profile, which is connected — so every component of
/// `profile ∩ f` reaches `∂f`, and a face can never carry both an arc and a loop. Every
/// earlier attempt at such a fixture died on that.
///
/// Leg bottoms sit at `z = 0.5` and `z = 0.45`: two coplanar faces of *one* operand
/// tripped the pre-cutover door guard, exactly as `u_prism`'s staggered prongs avoid.
/// And the legs span `y ∈ [0.65, 1.3]`, not `[0.7, 1.3]`, because `(1.4, 0.7)` lies on
/// the cap's fan diagonal `y = x/2` and `segment_crosses_face` would graze it.
fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let staple = Profile2d::polygon(vec![
        p2(0.1, 0.5),
        p2(0.6, 0.5),
        p2(0.6, 1.3),
        p2(0.8, 1.3),
        p2(0.8, 0.45),
        p2(1.4, 0.45),
        p2(1.4, 1.5),
        p2(0.1, 1.5),
    ])
    .unwrap();
    let __frame0 = datum_frame(
        &mut m,
        SketchPlane::from_origin_normal(
            Point3::from_array([0.0, 1.3, 0.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
        )
        .expect("a unit normal"),
    );
    let OpOutput::Extrude { solid: st, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __frame0,
            profile: staple,
            dist: 0.65,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    (m, l, st)
}

/// One face ring, as plane triples.
type Ring = Vec<[usize; 3]>;

/// A **holed reflex face** from the live engine, as plane triples: `(planes, p, outer, hole)`.
///
/// `Cut(L-prism, stub)` leaves the L's top cap carrying a hole — `"dimple"` a square one,
/// `"ell"` an L-shaped one. Both rings come from [`combinatorics::face_vertex_triples`] and
/// [`combinatorics::hole_rings`], which the boolean itself uses, so the fixture exercises only code
/// the kernel runs.
///
/// This replaces a helper that built its rings from `seam_paths_on`/`orient_seam_loop` — the
/// retired seam engine. The *properties* below are about `point_in_ring`/`every_ray`, which are
/// live and load-bearing (`nest_cells` picks a hole's host with them, `unify_coplanar_faces`
/// groups by them), so they had to be re-homed rather than deleted with their old fixture.
fn holed_face_rings(which: &str) -> (Vec<WorkingPlane>, usize, Ring, Ring) {
    let (m, l, stub) = if which == "dimple" {
        l_and_dimple()
    } else {
        l_and_ell_stub()
    };
    let (planes, p, outer, hole) = holed_face_rings_of(m, l, stub);
    assert_eq!(outer.len(), 6, "{which}: the L's cap is a reflex hexagon");
    (planes, p, outer, hole)
}

/// Cut `stub` out of `l` and hand back the first holed face's two rings, with the plane table
/// they are named in. The ring-length assertion belongs to the caller: a fixture built to put
/// a specific corner on a specific line **must** say how many corners it expects, because
/// `Profile2d` dissolves a vertex that sits mid-run on a straight edge and a silently
/// dissolved corner is a fixture that measures nothing.
fn holed_face_rings_of(
    mut m: Model,
    l: Handle<Solid>,
    stub: Handle<Solid>,
) -> (Vec<WorkingPlane>, usize, Ring, Ring) {
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).expect("the cut");
    m.rebuild_adjacency();
    let faces_tab = collect_planes(&m, r).unwrap();
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        surf_ix.insert(pi.face().expect("a real face table"), i);
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, _cyls) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(&m, r, &surf_ix).unwrap();
    for &fh in &m.shells.get(m.solids.get(r).outer).faces {
        let fp = surf_ix[&fh];
        let holes = combinatorics::hole_rings(
            &m,
            fh,
            fp,
            &inc,
            &crate::planes::test_judge(&planes),
            &plane_ix,
        )
        .unwrap();
        if let Some(hole) = holes.into_iter().next() {
            let outer = combinatorics::face_vertex_triples(
                &m,
                fh,
                fp,
                &inc,
                &crate::planes::test_judge(&planes),
                &plane_ix,
            )
            .unwrap();
            return (
                planes,
                plane_ix[fp].plane(),
                outer.poly().expect("a poly outer").triples.clone(),
                hole.poly().expect("a poly hole").triples.clone(),
            );
        }
    }
    panic!("no holed face");
}

#[test]
fn a_loop_is_inside_the_face_it_was_found_on() {
    // A hole ring never touches its face's outer ring: every node is *strictly* inside, so
    // `point_in_ring` must say so for all of them. The outer ring is the L's cap, a hexagon
    // with a reflex corner, so this is not a convex test.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    for t in &hole {
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
    }
}

#[test]
fn the_older_loops_are_inside_their_faces_too() {
    // Two hole shapes that must not move: a square one and an L-shaped one. Both sit
    // strictly inside the same reflex hexagon, and **every** clear ray agrees — the parity
    // cannot depend on which ray was cast, which is a second machine for free.
    for which in ["dimple", "ell"] {
        let (planes, p, outer, hole) = holed_face_rings(which);
        for t in &hole {
            let rays = combinatorics::every_ray(
                &crate::planes::test_judge(&planes),
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap(),
            )
            .unwrap();
            assert!(!rays.is_empty(), "{which}: no clear ray");
            assert!(rays.iter().all(|&x| x), "{which}: {rays:?}");
        }
    }
}

#[test]
fn a_loop_is_placed_by_where_it_is_not_by_how_it_winds() {
    // Containment is about **where** a ring is, never which way it runs. A hole ring is stored
    // clockwise about the face normal and its outer ring counter-clockwise, and neither
    // direction may enter the answer: reversing either must change nothing.
    //
    // And containment is **not symmetric** — the classic way to get this wrong is a test that
    // only ever asks it one way round.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    let (mut rev_outer, mut rev_hole) = (outer.clone(), hole.clone());
    rev_outer.reverse();
    rev_hole.reverse();
    for t in &hole {
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                *t,
                &combinatorics::ring_from_names(p, &outer).unwrap()
            )
            .unwrap()
        );
        assert!(
            combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                *t,
                &combinatorics::ring_from_names(p, &rev_outer).unwrap()
            )
            .unwrap(),
            "reversing the outer ring must not move the hole"
        );
    }
    for t in &outer {
        assert!(
            !combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                *t,
                &combinatorics::ring_from_names(p, &hole).unwrap()
            )
            .unwrap()
        );
        assert!(
            !combinatorics::point_in_ring(
                &crate::planes::test_judge(&planes),
                p,
                *t,
                &combinatorics::ring_from_names(p, &rev_hole).unwrap()
            )
            .unwrap(),
            "nor may reversing the hole swallow the outer ring"
        );
    }
}

#[test]
fn every_clear_ray_agrees() {
    // The ring is simple, so the parity cannot depend on the ray. `every_ray` returns one
    // answer per usable candidate and they must be unanimous; a disagreement means the ray
    // choice leaked into the result.
    //
    // (Its ancestor also pinned that *half* the candidates were unusable — that count came
    // from the retired staple fixture, whose loop and arc shared a plane. The holed L cap has
    // no such sharing, so only the unanimity survives the move.)
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    let rays = combinatorics::every_ray(
        &crate::planes::test_judge(&planes),
        p,
        hole[0],
        &combinatorics::ring_from_names(p, &outer).unwrap(),
    )
    .unwrap();
    assert!(!rays.is_empty(), "at least one candidate is clear");
    assert!(rays.iter().all(|&x| x), "and they agree: inside — {rays:?}");
}

#[test]
fn a_ring_inside_a_ring_is_what_nesting_looks_like() {
    // `nested_loops` has no operand in the suite that produces it — a polyhedral torus would.
    // The detector can still be aimed at real geometry: a holed face *is* a ring inside a ring,
    // which fires exactly the condition the nesting brick asks about. It is the detector under
    // test, not the fixture.
    let (planes, p, outer, hole) = holed_face_rings("dimple");
    assert!(
        combinatorics::point_in_ring(
            &crate::planes::test_judge(&planes),
            p,
            hole[0],
            &combinatorics::ring_from_names(p, &outer).unwrap()
        )
        .unwrap()
    );
    assert!(
        !combinatorics::point_in_ring(
            &crate::planes::test_judge(&planes),
            p,
            outer[0],
            &combinatorics::ring_from_names(p, &hole).unwrap()
        )
        .unwrap()
    );
}

/// A 40×40 plate with the given outline and a through pocket at `[20,30]×[10,20]`, as the
/// holed cap's two rings. The pocket's four wall **planes** are what the fixtures below aim
/// rays along; the outline decides which of them a ring node sits on.
fn grazed_plate(outline: &[[f64; 2]]) -> (Vec<WorkingPlane>, usize, Ring, Ring) {
    let profile = Profile2d::polygon(outline.iter().map(|&p| p2(p[0], p[1])).collect()).unwrap();
    let mut m = replay(&[extrude_log_op(profile, 12.0)]).unwrap();
    let plate = m.live_solids[0];
    // Through, so the cap really is holed rather than dimpled.
    let pocket = m.add_cuboid(
        Point3::from_array([20.0, 10.0, -1.0]),
        Point3::from_array([30.0, 20.0, 13.0]),
    );
    holed_face_rings_of(m, plate, pocket)
}

/// The ring node whose point is `(x, y)` — fixtures are written in coordinates and the engine
/// answers in plane triples, so this is where the two meet.
fn node_at(planes: &[WorkingPlane], ring: &Ring, x: f64, y: f64) -> [usize; 3] {
    let pts = ring_points(planes, ring);
    let i = pts
        .iter()
        .position(|p| (p[0] - x).abs() < 1e-9 && (p[1] - y).abs() < 1e-9)
        .unwrap_or_else(|| panic!("no node at ({x}, {y}) among {pts:?}"));
    ring[i]
}

/// ★ **A ray that grazes a corner is still a ray** — the four ways a ring can meet the ray's
/// line, each with an answer known by hand.
///
/// A ring node *on* the line used to make the parity ambiguous, so the candidate was thrown
/// away; with both of a vertex's candidates thrown away the whole question came back
/// `NO_CLEAR_RAY`, and a band of rotation angles died of it. The rule that resolves it is the
/// one `trace_transversal_face` has always used: look at the node's two off-line neighbours —
/// **opposite sides is a crossing, equal sides a touch**.
///
/// So both fixtures block **both** candidates. That matters: with one candidate left clear the
/// answer comes out anyway and the test would be green before the fix as well as after,
/// measuring nothing. Here `point_in_ring` is `NO_CLEAR_RAY` before and the truth after.
///
/// The outer outlines are chosen so that no vertex sits mid-run on a straight edge —
/// `Profile2d` dissolves those at construction — and the ring lengths are asserted so that a
/// dissolve fails loudly instead of quietly weakening the fixture.
#[test]
fn a_ray_that_grazes_a_corner_still_answers() {
    // The pocket's own corner `(20,10)`: its two candidate lines are `y=10` and `x=20`.
    //
    // `y=10` carries the run — the outline's `(5,10)→(0,10)` edge lies *on* it — and the run's
    // flanks `(5,20)` and `(0,0)` are on opposite sides, so the boundary crosses there.
    // `x=20` is grazed by the single corner `(20,30)`, whose flanks `(40,40)` and `(0,40)` are
    // also opposite. One crossing each way: the pocket corner is inside the plate.
    let a = [
        [0.0, 0.0],
        [40.0, 0.0],
        [40.0, 40.0],
        [20.0, 30.0],
        [0.0, 40.0],
        [0.0, 20.0],
        [5.0, 20.0],
        [5.0, 10.0],
        [0.0, 10.0],
    ];
    let (planes, p, outer, hole) = grazed_plate(&a);
    assert_eq!(outer.len(), 9, "every corner of the outline survived");
    let jd = crate::planes::test_judge(&planes);
    let outer_edges = combinatorics::ring_from_names(p, &outer).unwrap();
    let hole_edges = combinatorics::ring_from_names(p, &hole).unwrap();

    // Run-crossing and isolated-crossing, both saying "inside".
    let v = node_at(&planes, &hole, 20.0, 10.0);
    assert!(
        combinatorics::point_in_ring(&jd, p, v, &outer_edges).unwrap(),
        "the pocket's corner is inside the plate"
    );
    let rays = combinatorics::every_ray(&jd, p, v, &outer_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(rays.iter().all(|&x| x), "and unanimous — {rays:?}");

    // ★ Negative control, and the touch rule's own lock: from the outline's `(0,10)` the
    // `y=10` line grazes the pocket's `(20,10)→(30,10)` edge, whose flanks `(30,20)` and
    // `(20,20)` are on the *same* side. Nothing crossed, so `(0,10)` is outside the pocket —
    // and counting that touch as a crossing would make this ray disagree with the `x=0` one.
    let w = node_at(&planes, &outer, 0.0, 10.0);
    assert!(
        !combinatorics::point_in_ring(&jd, p, w, &hole_edges).unwrap(),
        "a plate corner is not inside the pocket"
    );
    let rays = combinatorics::every_ray(&jd, p, w, &hole_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(!rays.iter().any(|&x| x), "and unanimous — {rays:?}");
}

/// The isolated **touch** — the branch fixture A never reaches.
///
/// The outline's bottom notch rises to `(12,10)` and turns straight back down, so both its
/// neighbours `(8,0)` and `(16,0)` are below `y=10`: the ring touched the ray's line without
/// crossing it. Counting it as a crossing flips the parity of the `−x` ray alone, and the two
/// directions of one candidate would then contradict each other.
#[test]
fn a_ring_that_touches_the_ray_and_turns_back_crossed_nothing() {
    let b = [
        [0.0, 0.0],
        [8.0, 0.0],
        [12.0, 10.0], // touches y=10 and turns back — the branch under test
        [16.0, 0.0],
        [40.0, 0.0],
        [40.0, 40.0],
        [20.0, 30.0], // blocks the x=20 candidate, so neither is clear
        [0.0, 40.0],
        [0.0, 20.0],
        [5.0, 10.0], // crosses y=10
        [0.0, 5.0],
    ];
    let (planes, p, outer, hole) = grazed_plate(&b);
    assert_eq!(outer.len(), 11, "every corner of the outline survived");
    let jd = crate::planes::test_judge(&planes);
    let outer_edges = combinatorics::ring_from_names(p, &outer).unwrap();
    let v = node_at(&planes, &hole, 20.0, 10.0);
    assert!(
        combinatorics::point_in_ring(&jd, p, v, &outer_edges).unwrap(),
        "the pocket's corner is inside the plate"
    );
    let rays = combinatorics::every_ray(&jd, p, v, &outer_edges).unwrap();
    assert_eq!(rays.len(), 4, "both candidates are usable: {rays:?}");
    assert!(rays.iter().all(|&x| x), "and unanimous — {rays:?}");
}

/// The L with a stub rising out of its top face, footprint strictly inside that
/// face. Unlike the rod of `drill_through_the_l`, the stub enters the L from within,
/// so each of its vertical edges crosses exactly one face and its bottom ring stays
/// inside — one chord, no tunnel, and the seam loop is the whole story.
fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let stub = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.5]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    (m, l, stub)
}

/// **`Origin` no longer tells result faces apart.** The arrangement names every vertex
/// it emits by the three planes meeting there, so an operand corner the cut never touched
/// comes back as `Discovered`, exactly like a seam vertex. Nothing carries over as
/// `Constructed` (the arrangement builds no vertex from an original handle).
///
/// This is a contract, not a curiosity: `pipeline.rs`'s island test selected a face by
/// "all its vertices are `Discovered`", which was unique under the old engine and is
/// now true of *every* face. It flipped the wrong face and only the last assertion
/// noticed. Selecting a face by provenance is what this locks out.
///
/// The subject is `Cut(l, stub)` — the **holed** result, so `face_half_edges` walks
/// `inner` rings too (`count_discovered` walks only `outer` and would miss them).
///
/// **Unrotated only.** Rotating a boolean result re-marks these vertices `Rotated` over
/// a `Discovered` base — that is `transform_rotate_boolean_result_keeps_discovered_base`,
/// and this lock must not be read as contradicting it.
#[test]
fn an_unrotated_boolean_names_every_vertex_by_its_plane_triple() {
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();

    let mut seen = std::collections::HashSet::new();
    let mut holed = 0;
    for sh in solid_shell_handles(&m, r) {
        for &fh in &m.shells.get(sh).faces {
            let face = m.faces.get(fh);
            holed += usize::from(!face.inner.is_empty());
            for he in face_half_edges(face) {
                for vh in m.edges.get(he.edge).vertices {
                    if !seen.insert(vh) {
                        continue;
                    }
                    assert!(
                        matches!(m.vertices.get(vh).def, VertexDef::ThreePlane(_))
                            && m.vertex_tol(vh).is_some(),
                        "vertex {:?} is {:?}, not a measured plane triple",
                        m.vertex_point(vh).as_array(),
                        m.vertices.get(vh).def
                    );
                }
            }
        }
    }
    assert_eq!(holed, 1, "the blind dimple leaves exactly one holed face");
    assert_eq!(seen.len(), 20, "the L's 12 corners + the dimple's 8");
}

/// The unit cube with a 0.4-square pocket, 0.5 deep, in its top face: the void
/// is `[0.3,0.7]² × [0.5,1]` and the solid measures `1 − 0.16·0.5 = 0.92`. Its
/// lid is the only face in the suite that carries an inner loop.
fn pocketed_cube() -> (Model, Handle<Solid>) {
    let (mut m, top) = cube_with_top();
    let OpOutput::PocketOnFace { solid, .. } =
        apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap()
    else {
        unreachable!()
    };
    (m, solid)
}

/// A sever that also leaves a surviving cavity: a hollow box whose void sits to one side,
/// cut by a slab that severs it without touching the void. The x<2 piece keeps the void as a
/// cavity, the x>2 piece is solid — two outward shells *and* one inward. `point_in_component`
/// assigns the void to the x<2 piece that nests it (a containment test), rather than rejecting.
#[test]
fn severed_with_cavity_assigns_the_void_to_its_piece() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    // Void near the x-low side (1×2×2 = 4), clear of the x=2 cut.
    let inner = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 2.5, 2.5]),
    );
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(m.solids.get(hollow).cavities.len(), 1);
    // A slab spanning full y,z, thin in x at x∈[2,2.2] — severs into x<2 (holds the void, vol
    // 2·3·3 − 4 = 14) and x>2 (solid, vol 0.8·3·3 = 7.2).
    let slab = m.add_cuboid(
        Point3::from_array([2.0, -1.0, -1.0]),
        Point3::from_array([2.2, 4.0, 4.0]),
    );
    let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
    assert_eq!(
        solids.len(),
        2,
        "the slab severs the hollow box into two pieces"
    );
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // Exactly one piece owns the void; volumes match the hand calculation.
    let with_cav: Vec<_> = solids
        .iter()
        .filter(|&&s| !m.solids.get(s).cavities.is_empty())
        .collect();
    assert_eq!(
        with_cav.len(),
        1,
        "the void is assigned to exactly one piece"
    );
    let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
    let hollow_piece = *with_cav[0];
    assert!(
        (vol(hollow_piece) - 14.0).abs() < 1e-9,
        "hollow piece {}",
        vol(hollow_piece)
    );
    let total: f64 = solids.iter().map(|&s| vol(s)).sum();
    assert!((total - 21.2).abs() < 1e-9, "total {total}");
}

/// The adjacent case: a cut that passes *through* the void opens it — the void wall becomes
/// exterior boundary, so no cavity survives. Handled by the plain sever path (no cavity to
/// assign), not the containment code, but pinned so a regression there is caught.
#[test]
fn a_cut_through_the_void_leaves_no_cavity() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3])); // void 1³
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(m.solids.get(hollow).cavities.len(), 1);
    // Slab x∈[1.4,1.6] passes through the void (x∈[1,2]) → severs AND opens the void.
    let slab = m.add_cuboid(
        Point3::from_array([1.4, -1.0, -1.0]),
        Point3::from_array([1.6, 4.0, 4.0]),
    );
    let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // The void is opened, so neither piece keeps a cavity; material = 26 − (1.8 − 0.2) = 24.4.
    let total_cavities: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
    assert_eq!(
        total_cavities, 0,
        "the cut opened the void — no surviving cavity"
    );
    let total: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((total - 24.4).abs() < 1e-9, "total {total}");
}

/// Nested cavities: a hollow box A ([0,6]³ − [1,5]³ void) with a smaller hollow box B
/// ([2,4]³ − [2.5,3.5]³ void) floating inside A's void. `Fuse(A,B)` is one arrangement with
/// four components (two materials, two voids); B's void is contained by **both** A's outer
/// shell and B's own, so the containment assignment must pick the **innermost** (B), not A.
/// The result is two solids, each keeping its own void (A: 216−64 = 152, B: 8−1 = 7).
#[test]
fn a_void_nested_in_a_floating_island_goes_to_the_inner_solid() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([6.0; 3]));
    let void = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([5.0; 3]));
    let a = boolean_one(&mut m, BoolKind::Cut, big, void).unwrap();
    m.rebuild_adjacency();
    let bbig = m.add_cuboid(Point3::from_array([2.0; 3]), Point3::from_array([4.0; 3]));
    let bvoid = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
    let b = boolean_one(&mut m, BoolKind::Cut, bbig, bvoid).unwrap();
    m.rebuild_adjacency();
    let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // Each solid keeps exactly one void — B's void was assigned to B (innermost), not A.
    let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
    for &s in &solids {
        assert_eq!(
            m.solids.get(s).cavities.len(),
            1,
            "each piece keeps its own void"
        );
    }
    let mut vols: Vec<f64> = solids.iter().map(|&s| vol(s)).collect();
    vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
    assert!((vols[0] - 7.0).abs() < 1e-9, "inner {}", vols[0]);
    assert!((vols[1] - 152.0).abs() < 1e-9, "outer {}", vols[1]);
}

#[test]
fn replay_is_deterministic() {
    let log = vec![extrude_log_op(square(), 1.0)];
    let m1 = replay(&log).unwrap();
    let m2 = replay(&log).unwrap();
    let pts = |m: &Model| {
        m.vertices
            .iter()
            .map(|(vh, _)| m.vertex_point(vh).as_array())
            .collect::<Vec<_>>()
    };
    assert_eq!(pts(&m1), pts(&m2));
    assert_eq!(m1.edges.len(), m2.edges.len());
    assert_eq!(m1.faces.len(), m2.faces.len());
}

#[test]
fn two_extrudes_make_two_solids() {
    let far = SketchPlane::world_xy().with_origin(Point3::from_array([5.0, 0.0, 0.0]));
    // ★ The second plane is not a seed, so the log has to state it — and a datum's handle is
    // not known until the datum runs. So the log is assembled against a **scratch model built
    // the same way**: a frame that names surface *N* there names surface *N* in the replay,
    // because `replay` re-anchors indices (`docs/design.md` §2). That is what a recording
    // session does, and R is what makes it sound.
    let mut scratch = Model::new();
    let first = extrude_op(&scratch, square(), 1.0);
    apply(&mut scratch, &first).unwrap();
    let far_frame = datum_frame(&mut scratch, far);
    let log = vec![
        first,
        Operation::DatumPlane {
            def: DatumDef::Stated(far),
        },
        Operation::Extrude {
            frame: far_frame,
            profile: square(),
            dist: 1.0,
        },
    ];
    let m = replay(&log).unwrap();
    assert_eq!(m.solids.len(), 2);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Extrude a unit cube and return `(model, top face handle)`.
fn cube_with_top() -> (Model, Handle<Face>) {
    let mut m = Model::new();
    let op = extrude_op(&m, square(), 1.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
        unreachable!()
    };
    let top = faces[1]; // base, top, sides…
    (m, top)
}

/// The `0.4` square on `[0.3, 0.7]²` of the unit cube's lid.
///
/// ★ **On a lid these are world coordinates.** The sketch origin is the world origin projected
/// onto the plane, and the arbitrary-axis convention gives `n = ẑ` the axes `u = +x̂, v = +ŷ`,
/// so a frame point `(a, b)` is world `(a, b, 1)`.
fn small_square() -> Profile2d {
    Profile2d::polygon(vec![p2(0.3, 0.7), p2(0.3, 0.3), p2(0.7, 0.3), p2(0.7, 0.7)]).unwrap()
}

/// **Chaining onto a fused boss.** The fuse leaves the base's `z=1` face a *ring* — a face with
/// a hole where the boss sits — and the second boolean cuts through both. Every plane class the
/// cut opens then meets that ring along the **hole's own edge**, which is the case that used to
/// label inconsistently: the ring's neighbouring vertices there point *into* the hole, so
/// reading the occupied side off a flank put the material on the wrong side of `W`. The side
/// now comes from the ring's travel ([`arrangement::run_body_above`]), and the run leaves as its own
/// homogeneous segment rather than being swallowed by the straddling stretch beside it.
///
/// Hand volume: `1 + 0.5·0.5·1` fused, less the cutter's `0.2·0.2` column over `z ∈ [0.5, 2]`
/// — `1.25 − 0.06 = 1.19`. The same shape is scored against OCCT by
/// `boss_fuse_then_cut_matches_occt`, but that oracle is `#[ignore]`d, so this is the copy that
/// runs on every `cargo test`.
#[test]
fn a_boss_fused_then_cut_through() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.25, 0.25, 1.0]),
        Point3::from_array([0.75, 0.75, 2.0]),
    );
    let bossed = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, bossed).unwrap().volume - 1.25).abs() < 1e-12,
        "the fused boss itself"
    );
    let cutter = m.add_cuboid(
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, bossed, cutter).expect("the chained cut");
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, r).unwrap().volume - 1.19).abs() < 1e-12,
        "base + boss less the drilled column: {}",
        nacre_props::mass_props(&m, r).unwrap().volume
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "a chained result is still a clean model"
    );
}

/// An exact world-frame `Swept` from decimal f64 points — the truth-stating successor of the
/// retired `Swept::along` (S6b): the fallback is gone, so a test states its rings the way a
/// producer does. The f64 base is kept as handed in (decimals realize back bit-identically);
/// the top is the realization of the exact sum.
fn swept_world(base: Vec<Point3>, sweep: Vector3) -> crate::exact::Swept {
    let lift =
        |p: [f64; 3]| p.map(|x| nacre_scalar::Rat::from_decimal(x).expect("decimal fixture"));
    let sv = lift(sweep.as_array());
    let rb: Vec<[nacre_scalar::Rat; 3]> = base.iter().map(|p| lift(p.as_array())).collect();
    let rt: Vec<[nacre_scalar::Rat; 3]> = rb
        .iter()
        .map(|b| core::array::from_fn(|i| b[i].checked_add(sv[i]).expect("fixture widths")))
        .collect();
    let top = crate::exact::realize(&rt);
    crate::exact::Swept {
        base,
        top,
        exact: crate::exact::SweptRat {
            base: rb,
            top: rt,
            motion: None,
        },
    }
}

/// explicit sharing (overhaul #3): a prism built with a shared base-cap
/// surface reuses that `Surface` handle for its flush cap, and reconciles the
/// cap's face orientation so the materialized outward normal stays `−sweep`.
#[test]
fn build_prism_base_cap_reuses_shared_surface() {
    let mut m = Model::new();
    // A face-plane surface with outward normal +z (as a face on the base solid).
    let r = nacre_scalar::Rat::from_int;
    let (sf, _) = m.push_plane(
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap(),
        [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
        None,
    );
    let base_pts = [
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    ];
    let (_prism, faces) = build_prism(
        &mut m,
        swept_world(base_pts.to_vec(), Vector3::from_array([0.0, 0.0, 1.0])),
        vec![],
        Vector3::from_array([0.0, 0.0, 1.0]),
        Some(sf),
        None,
    )
    .unwrap();
    let cap = m.faces.get(faces[0]); // base cap is pushed first
    // Shared handle (was a fresh push before overhaul #3).
    assert_eq!(cap.surface, sf, "base cap reuses the shared surface handle");
    // Orientation reconciled: materialized outward normal is −sweep (−z).
    let Surface::Plane(p) = m.surface(cap.surface) else {
        unreachable!()
    };
    let sign = f64::from(cap.orientation.sign());
    let materialized = p.normal() * sign;
    assert!(
        (materialized - Vector3::from_array([0.0, 0.0, -1.0])).norm() < 1e-12,
        "materialized cap normal stays −z, got {materialized:?}"
    );
}

/// ★★★★★ **A plane's record and its motion are one statement.**
///
/// `Constructed` means *"these points speak about the world"*, `Moved` means *"about the
/// pre-motion frame"* — so a base cap that takes the caller's world triple must be
/// `Constructed`, and one that falls back to the prism's own ring must take the frame's `def`
/// with it. There is nothing else to keep in step: the plane's canonical name is derived from
/// whichever triple is recorded.
///
/// ★★ **It used to be three halves and they could come apart.** Coefficients were supplied
/// beside the points and chosen by a **separate** `or_else`, so a caller whose points
/// overflowed while their coefficients did not got world coefficients recorded beside the
/// prism's own frame ring — one plane stated two ways. The agreement filter could not catch it
/// either, since `c · p` overflows at exactly those widths. Removing the coefficient parameter
/// is what made the pairing structural; this pins what is left of the choice.
#[test]
fn a_prisms_base_cap_records_the_frame_its_def_names() {
    let r = nacre_scalar::Rat::from_int;
    let base_pts = [
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    ];
    let prism = |pts: Option<[[nacre_scalar::Rat; 3]; 3]>| {
        let mut m = Model::new();
        let (_prism, faces) = build_prism(
            &mut m,
            swept_world(base_pts.to_vec(), Vector3::from_array([0.0, 0.0, 1.0])),
            vec![],
            Vector3::from_array([0.0, 0.0, 1.0]),
            None,
            pts,
        )
        .unwrap();
        let surf = m.faces.get(faces[0]).surface; // base cap is pushed first
        (m, surf)
    };

    // ★ The caller's triple, when there is one — here a plane nowhere near the prism, so a
    // record that ignored it would be visibly different rather than coincidentally equal.
    let far_pts = [[r(0), r(0), r(3)], [r(1), r(0), r(3)], [r(0), r(1), r(3)]];
    let (m, surf) = prism(Some(far_pts));
    assert_eq!(
        m.surface_truth(surf),
        &nacre_topo::SurfaceTruth::Plane {
            points: nacre_topo::PlanePoints::Known(far_pts),
            motion: None,
        },
        "the caller's triple was not the one recorded (or gained a motion)"
    );
    assert_eq!(
        m.surface_name.get(&surf),
        Some(&nacre_scalar::PlaneName::Narrow([r(0), r(0), r(1), r(-3)])),
        "the name was not derived from the triple that was recorded"
    );

    // ★ And with no caller statement, the ring answers — its own exact cap triple, and the
    // name derived from it (`z = 0`, visibly different from the caller's `z = 3` above).
    // This used to pin the *absence* of both — the f64 fallback's point-less cap — and the
    // fallback is gone (S6b): the negative pins live at the operation as named rejects now.
    let (m, surf) = prism(None);
    assert!(
        matches!(
            m.surface_truth(surf),
            nacre_topo::SurfaceTruth::Plane { .. }
        ),
        "the ring's triple is recorded"
    );
    assert_eq!(
        m.surface_name.get(&surf),
        Some(&nacre_scalar::PlaneName::Narrow([r(0), r(0), r(1), r(0)])),
        "the name is derived from the ring's own plane"
    );
}

/// The handle branch of `shares_or_coplanar` is load-bearing: a shared
/// `Surface` handle reports coplanar even when the stored `plane` values are
/// *not* geometrically coplanar (so the fallback would not fire). This is the
/// path a referenced coplanar contact takes; on axis-aligned M5 it is redundant
/// with the geometric test, but the branch must work for rotated frames.
#[test]
fn shares_or_coplanar_uses_the_handle_branch() {
    let mut m = Model::new();
    let r = nacre_scalar::Rat::from_int;
    let shared = m.push_plane_unregistered(
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0])).unwrap(),
        [[r(0); 3], [r(0), r(1), r(0)], [r(0), r(0), r(1)]],
    );
    let fh = m.faces.push(Face {
        surface: shared,
        outer: Loop { half_edges: vec![] },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    let plane_x0 =
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0])).unwrap();
    let plane_z0 =
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap();
    // Each `tri` is three NON-collinear points of its own plane. A degenerate `tri` (three
    // equal points) would make every `orient3d` vanish, so the coordinate branch would report
    // coplanar and this test would pass without the handle branch ever mattering.
    let mk = |plane, tri: [Point3; 3]| {
        crate::planes::FaceRow::Plane(FaceInfo {
            // Unmoved and hand-built: nothing to record, and the base frame is unused anyway.
            base_rat: None,
            name: None,
            motion: None,
            surf: shared,
            face: Some(fh),
            plane,
            tri,
            n_out: Vector3::from_array([0.0; 3]),
            // Unread: this table only ever reaches `Judge::planes_coplanar`, which decides on `tri`.
            orient_sign: 1,
            tri_pt3: tri.map(|p| nacre_cip::WitnessPoint::exact(p.as_array()).expect("exact")),
            rotated: false,
        })
    };
    let p = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
    let planes = vec![
        mk(plane_x0, [p(0., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]), // in x = 0
        mk(plane_z0, [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.)]), // in z = 0
    ];
    // Neither fallback fires: the coefficients are not proportional, and the coordinates say
    // these really are two different planes.
    assert!(!planes_coplanar(
        &planes[0].plane().plane,
        &planes[1].plane().plane
    ));
    assert!(!crate::planes::test_judge(&planes).planes_coplanar(0, 1));
    // The shared handle alone makes them coplanar-by-reference.
    assert!(shares_or_coplanar(
        &crate::planes::test_judge(&planes),
        0,
        1
    ));
}

fn pocket_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
    Operation::PocketOnFace {
        face,
        profile,
        dist,
    }
}

proptest! {
    #[test]
    fn prop_regular_ngon_on_xy_is_clean(
        n in 3usize..8,
        r in 0.5f64..10.0,
        dist in 0.1f64..10.0,
    ) {
        let m = replay(&[extrude_log_op(regular_ngon(n, r), dist)]).unwrap();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        prop_assert_eq!(m.vertices.len(), 2 * n);
        prop_assert_eq!(m.faces.len(), n + 2);
    }

    #[test]
    fn prop_ngon_on_arbitrary_plane_is_clean(
        n in 3usize..8,
        nx in -1.0f64..1.0,
        ny in -1.0f64..1.0,
        nz in -1.0f64..1.0,
        dist in 0.1f64..10.0,
    ) {
        let normal = Vector3::from_array([nx, ny, nz]);
        prop_assume!(normal.norm() > 0.1); // skip near-zero normals
        let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
        let __plane = plane;
        let mut __scratch205 = Model::new();
        let __g204 = datum_frame(&mut __scratch205, __plane);
        let m = replay(&[
            Operation::DatumPlane { def: DatumDef::Stated(__plane) },
            Operation::Extrude {
                frame: __g204,
            profile: regular_ngon(n, 2.0),
            dist,
        }])
        .unwrap();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        prop_assert_eq!(m.faces.len(), n + 2);
    }

    /// A blind pocket on a randomly-slanted face: the arrangement must give a valid solid of the
    /// right volume or reject honestly — **never panic**. Drives general (non-axis) plane normals
    /// through the dir-sign guard and the `angular_order`/`turn_at` consumers that read its zeros.
    ///
    /// Was `#[ignore]`d for a residual `D = 0` panic on general normals: a triple naming one
    /// geometric plane through two coincident faces. Symmetric normals like `(1,1,1)` cleared
    /// it, `(0.446, 0.737, 0.990)` did not. **Un-ignored 2026-07-22** — naming every plane by
    /// its class made those two faces one index, so the degenerate triple can no longer form.
    #[test]
    fn pocket_on_a_random_slanted_face_is_valid_or_rejects(
        nx in -1.0f64..1.0,
        ny in -1.0f64..1.0,
        nz in 0.2f64..1.0, // keep the normal clear of the sketch's degenerate zero
    ) {
        let normal = Vector3::from_array([nx, ny, nz]);
        prop_assume!(normal.norm() > 0.3);
        let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
        let mut m = Model::new();
        let big = Profile2d::polygon(vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)]).unwrap();
        let frame = datum_frame(&mut m, plane);
        let OpOutput::Extrude { faces, .. } =
            apply(&mut m, &Operation::Extrude { frame, profile: big, dist: 2.0 }).unwrap()
        else { unreachable!() };
        match apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)) {
            Ok(OpOutput::PocketOnFace { solid, .. }) => {
                m.rebuild_adjacency();
                prop_assert!(nacre_validate::validate(&m).is_empty());
                let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
                prop_assert!((vol - (8.0 - 0.16 * 0.5)).abs() <= 1e-9 * 8.0, "volume {}", vol);
            }
            Ok(_) => prop_assert!(false, "unexpected op output"),
            Err(_) => {} // an honest reject is acceptable; a panic is not (and would fail the test)
        }
    }

    /// A random boss on a random box stays a valid b-rep (any interior
    /// profile, any positive height).
    #[test]
    fn prop_pad_stays_valid(
        sx in 0.5f64..5.0,
        sy in 0.5f64..5.0,
        sz in 0.5f64..5.0,
        h in 0.05f64..0.15,
        dist in 0.1f64..5.0,
    ) {
        let rect = Profile2d::polygon(vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)]).unwrap();
        let mut m = Model::new();
        let __w5 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
            frame: __w5,
            profile: rect,
            dist: sz,
        }).unwrap() else { unreachable!() };
        let hole = Profile2d::polygon(vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)]).unwrap();
        apply(&mut m, &Operation::PadOnFace { face: faces[1], profile: hole, dist }).unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A random blind pocket on a random box stays valid. `dist ≤ 0.8 < sz`
    /// keeps the pocket from punching through the box (height `sz ≥ 1`).
    #[test]
    fn prop_pocket_stays_valid(
        sx in 0.5f64..5.0,
        sy in 0.5f64..5.0,
        sz in 1.0f64..5.0,
        h in 0.05f64..0.15,
        dist in 0.1f64..0.8,
    ) {
        let rect = Profile2d::polygon(vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)]).unwrap();
        let mut m = Model::new();
        let __w4 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
            frame: __w4,
            profile: rect,
            dist: sz,
        }).unwrap() else { unreachable!() };
        let hole = Profile2d::polygon(vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)]).unwrap();
        apply(&mut m, &Operation::PocketOnFace { face: faces[1], profile: hole, dist }).unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
    }
}

// ---- boolean API (M5-c3 commit 1) ----

fn two_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    (m, a, b)
}

/// Outer A = [0,3]³ (volume 27) with inner B = [1,2]³ (volume 1) strictly
/// inside it — the containment fixture (returns `(model, outer, inner)`).
fn nested_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    (m, a, b)
}

/// The reported four-plane model: a unit cube with a block fused on its top, and a bar spun
/// `deg`° about an axis in the block's `x = 0.5` plane. At 45° with `half_z == 0.2` the bar's
/// half-width and its pivot-to-bottom offset are equal, so its bottom corner edge lands in that
/// plane and the y-planes cutting the edge become four-plane vertices.
fn four_plane_model(half_z: f64, deg: i128) -> (Model, Handle<Solid>, Handle<Solid>) {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let target = boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses on")[0];
    m.rebuild_adjacency();
    let bar = m.add_cuboid(
        Point3::from_array([0.3, -0.5, 1.0 - half_z]),
        Point3::from_array([0.7, 1.5, 1.0 + half_z]),
    );
    m.rebuild_adjacency();
    let bar = transform(
        &mut m,
        bar,
        &Isometry::rotation(Rotation {
            axis: Axis::Y,
            point: [Rat::new(1, 2).unwrap(), Rat::from_int(0), Rat::from_int(1)],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        }),
    )
    .unwrap();
    m.rebuild_adjacency();
    (m, target, bar)
}

/// **What the corpus actually contains in the way of concurrent vertices — and whether the
/// engine's rule for noticing them is complete.**
///
/// The four-plane work rests on one premise: a point's identity can be made a function of the
/// *set* of planes through it, because **every producer that meets the point derives the same
/// set**. Until now that was checked on a single model. This measures it against ground truth
/// (every plane asked, not the rule asking itself) over the whole fixture corpus.
///
/// Two things are asserted, and the second is the load-bearing one:
///
/// 1. **Exactly four.** No corpus point has five or more planes through it. That matters
///    because it is what makes both discovery rules complete: the trace learns `{wc} ∪ t`, and
///    with `|S| = 4` that *is* `S`. A five-plane point would leave it one short — so if this
///    ever fires, the identity rule needs the full set from somewhere else, and the message
///    says which model found it.
/// 2. **The trace's rule reproduces ground truth.** Wherever the trace would notice a
///    concurrency (a name on class `wc` that does not mention `wc`), the set it would record
///    equals the set every plane agrees on.
///
/// The count is printed rather than pinned: this is a *measurement*, and a number here would
/// only pin today's fixture list.
#[test]
fn concurrent_vertices_are_four_planes_and_the_trace_sees_all_of_them() {
    let mut boxed: Vec<(String, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
    macro_rules! fixture {
        ($name:ident) => {{
            let (m, a, b) = $name();
            boxed.push((stringify!($name).to_string(), m, a, b));
        }};
    }
    fixture!(two_boxes);
    fixture!(nested_boxes);
    fixture!(cube_and_notch);
    fixture!(stacked_cubes);
    fixture!(l_and_corner_box);
    fixture!(l_and_reflex_box);
    fixture!(l_and_inner_box);
    fixture!(l_and_popup_box);
    fixture!(l_and_notch_bar);
    fixture!(l_and_ell_stub);
    fixture!(l_and_staple);
    fixture!(l_and_dimple);
    fixture!(l_and_rod);
    fixture!(u_and_slab);
    // ...and the same shapes tilted, since a rotated operand is where concurrencies actually
    // turn up: an axis-aligned corpus would measure the easy half and call it the whole.
    macro_rules! tilted {
        ($name:ident) => {{
            let (mut m, a, b) = $name();
            let iso = rot_iso(nacre_scalar::Axis::Z, 30);
            let a = transform(&mut m, a, &iso).unwrap();
            m.rebuild_adjacency();
            let b = transform(&mut m, b, &iso).unwrap();
            m.rebuild_adjacency();
            boxed.push((format!("{} (tilted)", stringify!($name)), m, a, b));
        }};
    }
    tilted!(two_boxes);
    tilted!(nested_boxes);
    tilted!(cube_and_notch);
    tilted!(stacked_cubes);
    tilted!(l_and_corner_box);
    tilted!(l_and_reflex_box);
    tilted!(l_and_inner_box);
    tilted!(l_and_popup_box);
    tilted!(l_and_notch_bar);
    tilted!(l_and_ell_stub);
    tilted!(l_and_staple);
    tilted!(l_and_dimple);
    tilted!(l_and_rod);
    tilted!(u_and_slab);

    // ★ And the models that actually have one. Without these the sweep asserts nothing: the
    // corpus above turns out to carry **no** concurrent vertex at all, so it can measure the
    // blast radius of the coming stages but not the discovery rule. These are the reported
    // model (a bar spun 45° whose bottom corner edge lands in the block's x = 0.5 plane) and
    // its variants; `rejects.rs` pins the same shape from outside.
    for (label, half_z, deg) in [
        ("four_plane (the reported model)", 0.2, 45),
        ("four_plane at 44 deg (a near miss)", 0.2, 44),
        ("four_plane, taller bar (a near miss)", 0.3, 45),
    ] {
        let (m, t, bar) = four_plane_model(half_z, deg);
        boxed.push((label.to_string(), m, t, bar));
    }

    let (mut with_any, mut total, mut trace_rule_checked) = (0usize, 0usize, 0usize);
    let mut lines_found = 0usize;
    for (name, m, a, b) in &boxed {
        let found = arrangement::concurrency_audit(m, *a, *b).unwrap();
        if !found.is_empty() {
            with_any += 1;
        }
        total += found.len();
        for c in &found {
            lines_found += c.lines.len();
            assert_eq!(
                c.planes.len(),
                4,
                "{name}: a {}-plane point at {:?} — the trace learns only `{{wc}} ∪ t`, \
                     which is four names at most, so this one would be discovered incomplete",
                c.planes.len(),
                c.planes
            );
            if !c.triple.contains(&c.wc) {
                trace_rule_checked += 1;
                let mut derived = c.triple.to_vec();
                derived.push(c.wc);
                derived.sort_unstable();
                assert_eq!(
                    derived, c.planes,
                    "{name}: on class {} the trace would record {derived:?} for the point \
                         named {:?}, but every plane says {:?}",
                    c.wc, c.triple, c.planes
                );
            }
        }
    }
    println!(
        "concurrency audit: {with_any}/{} fixtures carry a concurrent vertex, {total} in total",
        boxed.len()
    );
    // A sweep that quietly stops finding anything reads as agreement, so say what was actually
    // exercised — and fail if the load-bearing assertion never ran.
    assert!(
        trace_rule_checked > 0,
        "no observation reached the trace's discovery condition, so the rule was not measured"
    );
    // The same discovery must also yield the *line* aliases — three planes sharing a line show
    // up as a sub-triple of `S` that names no point. In this corpus every concurrency comes
    // from exactly that (a tool edge lying in a target plane), so each one carries one.
    assert!(
        lines_found > 0,
        "no line-sharing triple was derived, yet every concurrency here comes from one"
    );
    println!("  and {lines_found} carried a line-sharing triple");
    println!("  of which {trace_rule_checked} exercised the trace's discovery rule");
}

/// **Every producer states its side in the label frame** — the invariant family #2 restored,
/// swept over the whole two-solid corpus (prints the interesting classes with `--nocapture`).
///
/// The arrangement states its cell labels as `[*_above, *_below]` about one direction per plane
/// class: the class root's **stored surface normal** (`Seated{body_above}` and `emit_faces`'
/// `flip` are written against it). `combinatorics::side_of` answers in the root's **outward** frame
/// instead, and the two are opposite exactly when the root face is `Reversed`
/// (`orient_sign == -1`) — which no `add_cuboid` face ever is, but a face an earlier boolean
/// re-emitted flipped is. `graze_above` read `side_of` raw, so on a pocket wall it flipped the
/// wrong label bit. It no longer reads a point's side at all — [`arrangement::run_body_above`] derives
/// the occupied side from the ring's travel, and the frame term cancels there because
/// `order_along`'s direction and the label frame are defined by the same stored normal — but
/// this class stays the corpus's only crossed-frame witness, so it is what would catch a
/// producer that regresses to a raw `side_of`.
///
/// What this pins, measured before the fix (2026-07-22):
/// - the pocket fixture has 5 `orient_sign == -1` classes carrying seated *and* graze segments
///   (4 walls + the floor); every other fixture has **no** `Reversed` root at all, which is why
///   the whole corpus passed with the frames crossed and why converting cannot regress it;
/// - the four wall classes stopped at `loop_orient_mismatch`; the floor class did **not** — its
///   four rim edges all carry a graze, so seated and graze were wrong *together*, consistently,
///   and the label survived verification while being inverted (a silent wrong, not a reject);
/// - the hand-derived sides on those classes, so a re-crossed frame fails here first.
#[test]
fn every_producer_states_its_side_in_the_label_frame() {
    let mut boxed: Vec<(&str, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
    macro_rules! fixture {
        ($name:ident) => {{
            let (m, a, b) = $name();
            boxed.push((stringify!($name), m, a, b));
        }};
    }
    fixture!(two_boxes);
    fixture!(nested_boxes);
    fixture!(cube_and_notch);
    fixture!(stacked_cubes);
    fixture!(l_and_corner_box);
    fixture!(l_and_reflex_box);
    fixture!(l_and_inner_box);
    fixture!(l_and_popup_box);
    fixture!(l_and_notch_bar);
    fixture!(l_and_ell_stub);
    fixture!(l_and_staple);
    fixture!(l_and_dimple);
    fixture!(l_and_rod);
    fixture!(u_and_slab);
    {
        // The pocket family: `pocketed_cube` is itself a boolean result, so its pocket walls
        // are `Reversed` faces. Box coordinates are `pocket_corner_cut`'s.
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.85, 0.85, 0.85]),
            Point3::from_array([1.15, 1.15, 1.15]),
        );
        boxed.push(("pocket_corner_cut", m, pc, bx));
    }

    let mut reversed_with_graze: Vec<String> = Vec::new();
    let mut reversed_seated_only: Vec<String> = Vec::new();
    for (name, m, a, b) in &boxed {
        let audits = arrangement::frame_audit(m, BoolKind::Cut, *a, *b).unwrap();
        for au in &audits {
            let interesting = au.orient_sign < 0 || au.failed_at.is_some();
            if !interesting {
                continue;
            }
            let where_ = format!(
                "{name}: wc={} pt={:?} n={:?} orient_sign={} seated={:?} graze={:?} trans={} \
                     declined={:?} failed_at={:?}",
                au.wc,
                au.root_point,
                au.root_normal,
                au.orient_sign,
                au.seated,
                au.grazes,
                au.transversals,
                au.declined,
                au.failed_at,
            );
            println!("{where_}");
            if au.orient_sign < 0 {
                if au.grazes.is_empty() {
                    if !au.seated.is_empty() {
                        reversed_seated_only.push(where_);
                    }
                } else {
                    reversed_with_graze.push(where_);
                }
            }
        }
    }
    println!("--- reversed-root classes carrying a graze (the set the fix moves) ---");
    for r in &reversed_with_graze {
        println!("  {r}");
    }
    println!("--- reversed-root classes with seated but no graze (the alternative's risk) ---");
    for r in &reversed_seated_only {
        println!("  {r}");
    }
    // The pocket fixture must keep supplying such classes, or this test has stopped exercising
    // the crossed-frame configuration and would pass vacuously.
    assert_eq!(
        reversed_with_graze.len(),
        5,
        "the pocket's 4 walls + floor are the corpus's only reversed-root classes with a graze"
    );
    assert!(
        reversed_with_graze.iter().all(|r| r.starts_with("pocket")),
        "no other fixture may have one: {reversed_with_graze:?}"
    );

    // Hand-derived sides on the pocket wall class `x = 0.7` (root = the pocket's +x wall, whose
    // outward normal points into the void, so the stored normal `+x` makes "above" the material
    // side `x > 0.7`). The wall is seated with its body above; the two side walls and the floor
    // graze it from `x < 0.7`, i.e. below. Crossed frames invert the grazes.
    //
    // The **box's top face** grazes it too, from `x > 0.7`: the pocket's opening makes that face
    // a notched region whose edge rides this plane with the material outside the pocket. That is
    // a run whose flanks *differ*, which the engine used to read as a straddling transversal —
    // this class is the corpus's only `Reversed` root, so it is also the only place the frame
    // handling of `arrangement::run_body_above` is exercised against a crossed frame: the three
    // `false` entries below are the pre-existing answers, unchanged by the new rule.
    let (m, pc, bx) = boxed
        .iter()
        .find_map(|(n, m, a, b)| (*n == "pocket_corner_cut").then_some((m, *a, *b)))
        .unwrap();
    let wall = arrangement::frame_audit(m, BoolKind::Cut, pc, bx)
        .unwrap()
        .into_iter()
        .find(|au| au.root_point == [0.7, 0.7, 1.0] && au.root_normal == [1.0, 0.0, 0.0])
        .expect("the x=0.7 pocket wall class");
    assert_eq!(wall.orient_sign, -1, "the pocket wall is a Reversed face");
    assert_eq!(
        wall.seated,
        vec![true; 4],
        "body above = material at x > 0.7"
    );
    assert_eq!(
        wall.grazes,
        vec![true, false, false, false],
        "the top face grazes from x > 0.7 (its pocket-opening edge, material outside); \
             the two side walls and the floor graze from x < 0.7"
    );
    assert_eq!(wall.failed_at, None, "the class labels consistently");
}

/// Lower corner of a solid's outer-shell vertex bounding box (for translation
/// tests: a rigid move shifts it by exactly the offset).
fn bbox_lo(m: &Model, s: Handle<Solid>) -> [f64; 3] {
    let mut lo = [f64::INFINITY; 3];
    let sh = m.solids.get(s).outer;
    for &fh in &m.shells.get(sh).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            {
                for vh in m.edges.get(he.edge).vertices {
                    let p = m.vertex_point(vh).as_array();
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k]);
                    }
                }
            }
        }
    }
    lo
}

fn count_discovered(m: &Model, s: Handle<Solid>) -> usize {
    let mut seen = std::collections::HashSet::new();
    let mut n = 0;
    let sh = m.solids.get(s).outer;
    for &fh in &m.shells.get(sh).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            {
                for vh in m.edges.get(he.edge).vertices {
                    if seen.insert(vh) && m.vertex_tol(vh).is_some() {
                        n += 1;
                    }
                }
            }
        }
    }
    n
}

fn test_iso() -> (nacre_scalar::Isometry, [f64; 3]) {
    use nacre_scalar::Rat;
    (
        nacre_scalar::Isometry::translation([
            Rat::new(7, 2).unwrap(),
            Rat::from_int(-4),
            Rat::from_int(11),
        ]),
        [3.5, -4.0, 11.0],
    )
}

/// A rational translation supersedes a cuboid: rigid, so volume/area are
/// invariant and the bounding box shifts by exactly the offset; validate/tess/
/// STEP all accept the moved solid, and the input drops from `live_solids`.
#[test]
fn transform_translate_cuboid() {
    let (iso, off) = test_iso();
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap();
    let lo0 = bbox_lo(&m, c);

    let c2 = transform(&mut m, c, &iso).unwrap();
    m.rebuild_adjacency();

    assert_eq!(m.live_solids, vec![c2], "input superseded");
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap();
    assert!(
        (after.volume - before.volume).abs() < 1e-12,
        "volume invariant"
    );
    assert!((after.area - before.area).abs() < 1e-12, "area invariant");
    let lo1 = bbox_lo(&m, c2);
    for k in 0..3 {
        assert!(
            (lo1[k] - (lo0[k] + off[k])).abs() < 1e-12,
            "bbox shifted by offset"
        );
    }
    assert!(nacre_tess::to_obj(&m).is_ok(), "moved solid tessellates");
    assert!(
        nacre_step::to_step(&m)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP"),
        "moved solid exports to STEP"
    );
}

/// Transforming a boolean *result* (which carries `Discovered` seam vertices) does not
/// **downgrade** those vertices to `Constructed` — the failure this guards.
///
/// A motion that records a forest node supersedes the origin with `Moved { base, motion }`,
/// and the definition is preserved *through the base*: the base is the seam vertex, still in
/// the arena with its `ThreePlane` definition, and the node says how it moved. A motion that
/// records nothing instead remaps the definition's planes in place. Either way the truth
/// survives; only its spelling depends on whether the motion was worth recording.
#[test]
fn transform_translate_preserves_discovered_definition() {
    let (iso, _) = test_iso();
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    let before = nacre_props::mass_props(&m, r).unwrap().volume;
    let disc = count_discovered(&m, r);
    assert!(
        disc > 0,
        "the Cut result must have Discovered seam vertices"
    );

    let r2 = transform(&mut m, r, &iso).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    // Every seam vertex still names its three-plane definition — directly, or through the
    // base of the motion that moved it. None fell back to `Constructed`.
    let named = boundary_verts(&m, r2)
        .into_iter()
        // A measured tolerance survives an exact move and is dropped by a recorded one
        // (S7's tol rule), so what "still named" means here is the definition: every seam
        // vertex names three planes, whichever road it travelled.
        .filter(|&vh| matches!(m.vertices.get(vh).def, VertexDef::ThreePlane(_)))
        .count();
    assert_eq!(named, disc, "seam definitions preserved");
    let after = nacre_props::mass_props(&m, r2).unwrap().volume;
    assert!((after - before).abs() < 1e-12, "volume invariant");
}

/// Replay determinism (DNA 3): the same construction + transform reproduces the
/// same geometry and the same handle down to the index.
///
/// The `assert_eq!` below compares handles minted by two *different* `Model`s, which is
/// legal only because `Handle`'s equality is its index — the very premise `replay` now
/// relies on. `tests/replay.rs` measures that premise directly instead of assuming it.
#[test]
fn transform_is_deterministic() {
    let (iso, _) = test_iso();
    let build = || {
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c2 = transform(&mut m, c, &iso).unwrap();
        (bbox_lo(&m, c2), c2)
    };
    assert_eq!(build(), build(), "same ops → same geometry and handle");
}

/// The bit-identity guard, on the one producer that can break it.
///
/// A definition and its cached coordinate must agree **exactly**, and they do only because the
/// replay performs the very same float operations the producer did. A reflection is the step
/// where that is easiest to lose — `WitnessPoint::mirror` must walk the same `2c − x` that
/// `AxisMirror::point` just walked — and a chain that ends in one is what this checks.
#[test]
fn a_mirrored_rotated_vertex_reconstructs_from_its_definition() {
    use nacre_scalar::{Axis, Rat};
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let r = transform(&mut m, c, &rot30()).unwrap();
    m.rebuild_adjacency();
    let mirrored = crate::transform::mirror(&mut m, r, Axis::X, Rat::from_int(0)).unwrap();
    m.rebuild_adjacency();

    let shell = m.solids.get(mirrored).outer;
    let mut checked = 0;
    for &fh in &m.shells.get(shell).faces.clone() {
        for he in &m.faces.get(fh).outer.half_edges.clone() {
            for vh in m.edges.get(he.edge).vertices.iter() {
                // The image keeps its definition, and the definition reproduces the
                // coordinate: solve the corner's three planes in the frame their names are
                // stated in, replay the chain the *faces* record (S7 — the vertex has no
                // motion of its own), and the answer is the stored coordinate bit for bit.
                let VertexDef::ThreePlane(tri) = m.vertices.get(*vh).def else {
                    unreachable!("a cuboid corner is a three-plane point")
                };
                let motion_of = |h| match m.surface_truth(h) {
                    nacre_topo::SurfaceTruth::Plane { motion, .. } => *motion,
                    nacre_topo::SurfaceTruth::Cylinder { motion, .. } => *motion,
                };
                // The caps were fixed by the Z turn and by the X mirror (normal ⊥ both),
                // so they are world-stated; the moved carriers share the one leaf whose
                // chain the fixed caps provably survive — the reconstruction below runs
                // through the walls' recorded [rotate, mirror] history unchanged.
                let leaves: Vec<_> = tri.iter().filter_map(|&h| motion_of(h)).collect();
                let rotation = *leaves
                    .first()
                    .expect("a mirrored image's corner names a moved carrier");
                assert!(
                    leaves.iter().all(|&l| l == rotation),
                    "the corner's moved planes share one motion leaf"
                );
                let mut coeffs = [[Rat::from_int(0); 4]; 3];
                for (o, h) in coeffs.iter_mut().zip(tri) {
                    *o = *m.surface_name.get(&h).unwrap().narrow().unwrap();
                }
                let base = nacre_scalar::three_planes_rat(coeffs)
                    .expect("three distinct planes meet in a point");
                let replayed = crate::rotated_vertex::replay_chain_coord(
                    &m,
                    [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
                    rotation,
                )
                .expect("the root coordinate lifts to an exact rational");
                assert_eq!(
                    replayed,
                    m.vertex_point(*vh).as_array(),
                    "definition reproduces the stored coordinate bit for bit"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "walked no vertices");
}

/// Replay determinism (DNA 3) for the one additive operation: a copy reproduces the same
/// geometry *and* the same handle index, and leaves the same live set behind it.
///
/// The `assert_eq!` below compares handles minted by two *different* `Model`s, which is
/// legal only because `Handle`'s equality is its index — the very premise `replay` now
/// relies on. `tests/replay.rs` measures that premise directly instead of assuming it.
#[test]
fn copy_is_deterministic() {
    let build = || {
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let twin = crate::transform::copy(&mut m, c).unwrap();
        (bbox_lo(&m, twin), twin, m.live_solids.clone())
    };
    assert_eq!(build(), build(), "same ops → same geometry and handles");
}

/// A genuinely tilted rigid rotation: 30° about Z through the rational axis
/// point (1,1,0). Non-90° and non-axis-aligned, so it exercises the Rotated
/// origin and the boolean reject guard (unlike the 90° family, which stays exact).
fn rot30() -> nacre_scalar::Isometry {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    })
}

/// Distinct outer-shell vertex points of a solid (dedup by handle).
fn outer_points(m: &Model, s: Handle<Solid>) -> Vec<[f64; 3]> {
    let mut seen = std::collections::HashSet::new();
    let mut pts = Vec::new();
    let sh = m.solids.get(s).outer;
    for &fh in &m.shells.get(sh).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            {
                for vh in m.edges.get(he.edge).vertices {
                    if seen.insert(vh) {
                        pts.push(m.vertex_point(vh).as_array());
                    }
                }
            }
        }
    }
    pts
}

/// A non-90° rotation genuinely tilts the solid: rigid (volume/area invariant),
/// validate/tess/STEP clean, a known corner lands at its exact rotated image, the
/// faces record their motion (`solid_is_rotated`), and a boolean against it now
/// runs (a *mixed*-rotation cut: rotated `c2` minus an axis-aligned `d` it contains, so
/// `d` becomes a cavity — overhaul 3d-i retired the `ROTATED_UNSUPPORTED` entry guard).
#[test]
fn transform_rotate_cuboid_tilts_and_cuts() {
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap();
    let c2 = transform(&mut m, c, &rot30()).unwrap();
    m.rebuild_adjacency();

    assert_eq!(m.live_solids, vec![c2], "input superseded");
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap();
    assert!(
        (after.volume - before.volume).abs() < 1e-9,
        "volume invariant"
    );
    assert!((after.area - before.area).abs() < 1e-9, "area invariant");
    assert!(nacre_tess::to_obj(&m).is_ok(), "rotated solid tessellates");
    assert!(
        nacre_step::to_step(&m)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP"),
        "rotated solid exports to STEP"
    );

    assert!(solid_is_rotated(&m, c2), "the faces record their rotation");

    // Corner (0,0,0) rotates about pivot (1,1) by 30°: dx=dy=-1, so
    // x' = 1 - cos30 + sin30, y' = 1 - sin30 - cos30, z' = 0.
    let (c30, s30) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let want = [1.0 - c30 + s30, 1.0 - s30 - c30, 0.0];
    let pts = outer_points(&m, c2);
    assert!(
        pts.iter()
            .any(|p| p.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9)),
        "corner (0,0,0) rotated to its exact image {want:?}; got {pts:?}"
    );

    // A mixed-rotation cut now runs: the axis-aligned `d` sits inside the rotated `c2`, so
    // `Cut(c2, d)` leaves `c2` with `d` carved out as a cavity (volume 24 − 1 = 23).
    let d = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "mixed cut is valid"
    );
    assert_eq!(
        m.solids.get(r).cavities.len(),
        1,
        "the contained box is a cavity"
    );
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 23.0).abs() < 1e-9, "volume {vol}");
}

/// A 90° rotation about Z is axis-aligned and exact: the solid stays
/// `Constructed` (`solid_is_rotated` false), volume is exact, and boolean is
/// still allowed — a following cut succeeds and validates.
#[test]
fn transform_rotate_90_is_exact_and_allows_boolean() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let rot90 = Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
    });
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c2 = transform(&mut m, c, &rot90).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    assert!(!solid_is_rotated(&m, c2), "90° stays exact (Constructed)");
    assert_eq!(
        nacre_props::mass_props(&m, c2).unwrap().volume,
        24.0,
        "exact volume"
    );

    // c rotated 90° about origin occupies x∈[-3,0], y∈[0,2], z∈[0,4].
    // Cut with d = [-1,0.5,1]-[0.5,1.5,2]: overlap volume 1 → 24 − 1 = 23.
    let d = m.add_cuboid(
        Point3::from_array([-1.0, 0.5, 1.0]),
        Point3::from_array([0.5, 1.5, 2.0]),
    );
    let r = boolean(&mut m, BoolKind::Cut, c2, d).expect("exact rotation → boolean allowed");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r[0]).unwrap().volume;
    assert!((vol - 23.0).abs() < 1e-9, "cut volume {vol}");
}

/// Rotating a boolean *result* marks its `Discovered` seam vertices `Rotated`
/// over the pre-rotation vertex as base: validate stays clean, volume is
/// invariant, and at least one Rotated base is a Discovered vertex (the seam
/// definition is preserved through the rotation, not downgraded).
#[test]
fn transform_rotate_boolean_result_keeps_discovered_base() {
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    assert!(
        count_discovered(&m, r) > 0,
        "Cut result has Discovered seams"
    );
    let before = nacre_props::mass_props(&m, r).unwrap().volume;

    let r2 = transform(&mut m, r, &rot30()).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    assert!(
        solid_is_rotated(&m, r2),
        "rotated result carries Rotated origin"
    );
    let after = nacre_props::mass_props(&m, r2).unwrap().volume;
    assert!((after - before).abs() < 1e-9, "volume invariant");

    // A seam vertex of the *moved* result still names three planes, and those planes are
    // the moved ones (S7: the vertex follows its faces' motion — there is no base vertex
    // left to chase).
    let sh = m.solids.get(r2).outer;
    let own_surfaces: std::collections::HashSet<_> = m
        .shells
        .get(sh)
        .faces
        .iter()
        .map(|&fh| m.faces.get(fh).surface)
        .collect();
    let mut found_moved_carrier = false;
    for &fh in &m.shells.get(sh).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            for vh in m.edges.get(he.edge).vertices {
                let VertexDef::ThreePlane(tri) = m.vertices.get(vh).def else {
                    continue;
                };
                // The carriers are the result's own planes — never a superseded
                // pre-move twin. (A restated cap's handle *equals* its pre-move
                // handle by intern-back, which is why this is a subset check on the
                // live faces' surfaces rather than "all carriers moved".)
                for h in tri {
                    assert!(
                        own_surfaces.contains(&h),
                        "a vertex carrier is not one of the solid's own planes"
                    );
                }
                if tri.iter().any(|&h| {
                    !matches!(
                        m.surface_truth(h),
                        nacre_topo::SurfaceTruth::Plane { motion: None, .. }
                    )
                }) {
                    found_moved_carrier = true;
                }
            }
        }
    }
    assert!(
        found_moved_carrier,
        "a moved result's vertices name its moved planes"
    );
}

/// Replay determinism (DNA 3): the same construction + rotation reproduces the
/// same geometry and the same handle.
///
/// The `assert_eq!` below compares handles minted by two *different* `Model`s, which is
/// legal only because `Handle`'s equality is its index — the very premise `replay` now
/// relies on. `tests/replay.rs` measures that premise directly instead of assuming it.
#[test]
fn transform_rotate_is_deterministic() {
    let build = || {
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c2 = transform(&mut m, c, &rot30()).unwrap();
        (bbox_lo(&m, c2), c2)
    };
    assert_eq!(build(), build(), "same ops → same geometry and handle");
}

/// A rotation `Transform` flows through `apply` and marks the result Rotated.
#[test]
fn transform_rotate_op_applies() {
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let out = apply(
        &mut m,
        &Operation::Transform {
            solid: c,
            isometry: rot30(),
        },
    )
    .unwrap();
    let OpOutput::Transform { solid } = out else {
        panic!("expected Transform output, got {out:?}");
    };
    assert_eq!(m.live_solids, vec![solid]);
    assert!(solid_is_rotated(&m, solid));
}

fn boundary_verts(m: &Model, s: Handle<Solid>) -> Vec<Handle<Vertex>> {
    let mut seen = std::collections::HashSet::new();
    let mut vs = Vec::new();
    let sh = m.solids.get(s).outer;
    for &fh in &m.shells.get(sh).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            {
                for vh in m.edges.get(he.edge).vertices {
                    if seen.insert(vh) {
                        vs.push(vh);
                    }
                }
            }
        }
    }
    vs
}

/// Whether any wall of `s` is a **rotated image** — what `planes::solid_is_rotated` used to
/// ask of the vertices, now asked of the surfaces that actually record it.
fn solid_is_rotated(m: &Model, s: Handle<Solid>) -> bool {
    m.shells.get(m.solids.get(s).outer).faces.iter().any(|&fh| {
        matches!(
            m.surface_truth(m.faces.get(fh).surface),
            nacre_topo::SurfaceTruth::Plane {
                motion: Some(_),
                ..
            } | nacre_topo::SurfaceTruth::Cylinder {
                motion: Some(_),
                ..
            }
        )
    })
}

fn rot_iso(axis: nacre_scalar::Axis, deg: i128) -> nacre_scalar::Isometry {
    use nacre_scalar::{Angle, Isometry, Rat, Rotation as SRot};
    Isometry::rotation(SRot {
        axis,
        point: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
    })
}

/// Chain: `(node_count, axes-root-to-leaf)` of the **fullest** face history of `s` — the
/// walls', which go through every motion. It used to read `faces[0]`, but that is a cap,
/// and a cap a rotation fixes is restated world-side now (no motion, or a shorter chain
/// begun by a later non-fixing motion) — the forest's story lives on the faces that
/// genuinely moved. `None` when no face carries a rotation.
///
/// ★ S7 dropped a third field, `base_is_rotated` — "does this vertex's base vertex itself
/// carry a motion?", the one-hop invariant that kept a replay from applying the same motion
/// twice. There is no base vertex any more (a moved corner is the intersection of its moved
/// planes), so double application is **unrepresentable** rather than merely untrue: the type
/// absorbed the invariant, and a probe field that could only ever read `false` would be
/// theatre.
fn forest_probe(m: &Model, s: Handle<Solid>) -> Option<(usize, Vec<nacre_scalar::Axis>)> {
    let mut best: Option<Vec<nacre_scalar::Axis>> = None;
    for &fh in &m.shells.get(m.solids.get(s).outer).faces {
        let &nacre_topo::SurfaceTruth::Plane {
            motion: Some(rotation),
            ..
        } = m.surface_truth(m.faces.get(fh).surface)
        else {
            continue;
        };
        let mut axes = Vec::new();
        let mut cur = Some(rotation);
        while let Some(h) = cur {
            let n = m.motion(h);
            if let nacre_topo::Motion::Rotate { axis, .. } = n.motion {
                axes.push(axis);
            }
            cur = n.parent;
        }
        axes.reverse();
        if best.as_ref().is_none_or(|b| axes.len() > b.len()) {
            best = Some(axes);
        }
    }
    best.map(|axes| (axes.len(), axes))
}

/// **A chained boolean must not lose exactness.**
///
/// A rotated solid's face coordinates are rounded, so its planes are truthful only through an
/// exact *definition* (`FaceInfo::tri_pt3` built from a rotation history). A boolean's *result*
/// is just as rotated as its operands — but the result carries no rotation provenance, so
/// `collect_planes` describes every one of its faces by `WitnessPoint::exact` of the rounded triangle
/// and the kernel starts treating a rounded copy as the truth. That is what makes one wall
/// become two plane classes on the next operation.
///
/// The invariant: **every face of a boolean between rotated operands is described by a
/// rotation definition, not by its rounded coordinates.**
#[test]
fn a_chained_boolean_keeps_its_faces_exact() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.0]),
        Point3::from_array([1.0, 1.0, 3.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.5, -0.2, 1.0]),
        Point3::from_array([4.0, 0.2, 3.0]),
    );
    m.rebuild_adjacency();
    let a = transform(&mut m, a, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let b = transform(&mut m, b, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    // The operands hold up: a rotated solid's faces do carry definitions. The caps'
    // *truth* is world-stated since the invariant-plane restatement, but the judging
    // table re-chains a chain-fixed plane (the mirror takes the strongest description),
    // so the whole table reads rotated — exactly as it did when the twin carried the
    // motion itself.
    for (name, s) in [("operand a", a), ("operand b", b)] {
        let planes = crate::planes::collect_planes(&m, s).expect("planes");
        assert!(
            planes.iter().all(|f| f.plane().rotated),
            "{name}: a rotated operand's faces must be described by their rotation"
        );
    }
    let r = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse")[0];
    m.rebuild_adjacency();
    let planes = crate::planes::collect_planes(&m, r).expect("planes");
    let described = planes.iter().filter(|f| f.plane().rotated).count();
    assert_eq!(
        described,
        planes.len(),
        "the result of a rotated boolean is rotated too: {described}/{} faces carry a \
             definition, the rest are rounded coordinates declared exact",
        planes.len()
    );
}

/// **The same motion, applied twice, is the same node.**
///
/// The motion handle is the canonical name of "which motion", and judgments use it to decide
/// whether a whole judgement can be answered exactly in the pre-motion frame. Without
/// interning, turning two solids by the same 30° makes two nodes, their shared motion stops
/// cancelling, and rotating a model turns its exact questions into assumed ones — which is
/// what `a_shared_rotation_still_assumes_nothing` measures. (The identity it replaced was a
/// 64-bit hash of the chain's contents, where a collision would have answered a *different*
/// question with full confidence.)
#[test]
fn the_same_motion_applied_twice_is_one_node() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([2.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let a = transform(&mut m, a, &rot_iso(Axis::X, 30)).unwrap();
    m.rebuild_adjacency();
    let b = transform(&mut m, b, &rot_iso(Axis::X, 30)).unwrap();
    m.rebuild_adjacency();
    fn leaf(m: &Model, s: Handle<Solid>) -> Handle<nacre_topo::MotionNode> {
        let sh = m.solids.get(s).outer;
        let fh = m.shells.get(sh).faces[0];
        match m.surface_truth(m.faces.get(fh).surface) {
            nacre_topo::SurfaceTruth::Plane {
                motion: Some(motion),
                ..
            } => *motion,
            other => panic!("a rotated solid's faces record their motion, got {other:?}"),
        }
    }
    assert_eq!(leaf(&m, a), leaf(&m, b), "one motion, one node");
    // …and a *different* motion is a different node, or the identity would be worthless.
    let c = m.add_cuboid(
        Point3::from_array([5.0, 0.0, 0.0]),
        Point3::from_array([6.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let c = transform(&mut m, c, &rot_iso(Axis::X, 31)).unwrap();
    m.rebuild_adjacency();
    assert_ne!(
        leaf(&m, a),
        leaf(&m, c),
        "different motions must not share a node"
    );
}

/// **A boolean result rotated again continues its history — per wall.**
///
/// One solid does not have one rotation history. A result's vertices are all `Discovered`, so
/// asking them "what rotation is this solid at" answers `None` and the next rotation would
/// start a fresh root — replaying a pre-first-rotation witness through only the *second*
/// rotation, which is a plane that does not exist. And its walls can come from operands
/// rotated by different angles, so there is no single answer to give. Each surface therefore
/// chains from its own leaf.
#[test]
fn a_rerotated_boolean_result_continues_each_walls_history() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    // Different angles, so the two operands' walls carry genuinely different histories.
    let a = transform(&mut m, a, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let b = transform(&mut m, b, &rot_iso(Axis::Z, 50)).unwrap();
    m.rebuild_adjacency();
    let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap()[0];
    m.rebuild_adjacency();
    let r = transform(&mut m, r, &rot_iso(Axis::Z, 20)).unwrap();
    m.rebuild_adjacency();

    let mut leaves = std::collections::HashSet::new();
    let mut caps = 0usize;
    let sh = m.solids.get(r).outer;
    for &fh in &m.shells.get(sh).faces {
        let s = m.faces.get(fh).surface;
        match m.surface_truth(s) {
            &nacre_topo::SurfaceTruth::Plane {
                motion: Some(rotation),
                ..
            } => {
                assert_eq!(
                    crate::rotated_vertex::motion_chain(&m, rotation)
                        .expect("an axis-aligned history holds no frame")
                        .len(),
                    2,
                    "both rotations, once each"
                );
                leaves.insert(rotation);
            }
            // Every rotation in this chain is about Z, so the z-caps are fixed by all
            // of them and stay restated — the only faces allowed to carry nothing.
            nacre_topo::SurfaceTruth::Plane { motion: None, .. } => {
                let n = m.surface_name.get(&s).unwrap().narrow().unwrap();
                assert!(
                    n[0] == nacre_scalar::Rat::from_int(0)
                        && n[1] == nacre_scalar::Rat::from_int(0),
                    "only a Z-fixed cap may go without a history here"
                );
                caps += 1;
            }
            other => panic!("unexpected truth {other:?}"),
        }
    }
    assert_eq!(leaves.len(), 2, "the two operands' histories stay apart");
    assert!(caps > 0, "the restated caps are present in the result");
}

/// **A copy of a rotated solid is still exactly defined.**
///
/// `copy` is `transform` under a *zero* translation, so a surface rule that reads "any
/// translation makes a rotated plane inexpressible" swallows it — and then a copy of a rotated
/// solid cannot take part in a boolean at all, though its geometry is bit-identical to the
/// original's. The identity is not a translation.
#[test]
fn a_copy_of_a_rotated_solid_answers_like_the_original() {
    use nacre_scalar::Axis;
    let build = |use_copy: bool| {
        let mut m = Model::new();
        let hub = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, 0.0]),
            Point3::from_array([1.0, 1.0, 3.0]),
        );
        let fin = m.add_cuboid(
            Point3::from_array([0.5, -0.2, 1.0]),
            Point3::from_array([4.0, 0.2, 3.0]),
        );
        m.rebuild_adjacency();
        let mut fin = transform(&mut m, fin, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        if use_copy {
            fin = crate::transform::copy(&mut m, fin).unwrap();
            m.rebuild_adjacency();
        }
        let r = boolean(&mut m, BoolKind::Fuse, hub, fin).expect("fuse");
        nacre_props::mass_props(&m, r[0]).expect("props").volume
    };
    // Bit-identical, not merely close: the copy walks the same definitions.
    assert_eq!(build(true), build(false));
}

fn translate_iso(off: [i128; 3]) -> nacre_scalar::Isometry {
    use nacre_scalar::{Isometry, Rat};
    Isometry::translation([
        Rat::from_int(off[0]),
        Rat::from_int(off[1]),
        Rat::from_int(off[2]),
    ])
}

fn rigid_iso(axis: nacre_scalar::Axis, deg: i128, off: [i128; 3]) -> nacre_scalar::Isometry {
    use nacre_scalar::{Angle, Isometry, Rat, Rotation as SRot};
    Isometry::rigid(
        SRot {
            axis,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        },
        [
            Rat::from_int(off[0]),
            Rat::from_int(off[1]),
            Rat::from_int(off[2]),
        ],
    )
}

/// A boolean produces `Discovered` seam vertices with tol 0 (exact axis-aligned
/// intersections). Before exact quadrantal realization, rotating them exactly 90°
/// left an ~8e-17 f64 residual that exceeded tol 0 → `VertexOffSurface`. Now the
/// rotation is exact, so the residual stays 0 and validate is clean — both for a
/// pure 90° rotation and for a rigid 90°+translation (the offset cancels in
/// vertex−plane, so it does not reintroduce a residual).
#[test]
fn boolean_result_rotated_90_validates() {
    use nacre_scalar::Axis;
    for iso in [rot_iso(Axis::Z, 90), rigid_iso(Axis::Z, 90, [5, -3, 2])] {
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        let before = nacre_props::mass_props(&m, r).unwrap().volume;
        let r2 = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(
            vs.is_empty(),
            "exact 90° realization → validate clean: {vs:?}"
        );
        let after = nacre_props::mass_props(&m, r2).unwrap().volume;
        assert!((after - before).abs() < 1e-12, "volume invariant");
    }
}

/// A 90°-family rotation lands a cuboid's vertices exactly on the axis-aligned grid
/// (no ~6e-17 spurious offset): the corner (2,3,4) rotated 90° about Z maps to
/// exactly (-3,2,4).
#[test]
fn rotate_90_lands_vertices_exactly() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c2 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
    m.rebuild_adjacency();
    let pts = outer_points(&m, c2);
    // (2,3,4) about Z by 90°: (x,y)→(-y,x) → (-3,2,4). Bit-exact.
    assert!(
        pts.iter().any(|p| *p == [-3.0, 2.0, 4.0]),
        "corner lands exactly on the grid; got {pts:?}"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Re-rotating about the same axis chains a second forest node onto the first
/// (this cell does not bundle): the leaf's parent is the earlier rotation, `base`
/// stays the Constructed root, and the solid remains a rigid (volume/area-invariant)
/// `Rotated` solid that validate/tess/STEP accept and a boolean now runs against.
#[test]
fn rerotate_same_axis_chains() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap();
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::Z, 20)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c2),
        Some((2, vec![Axis::Z, Axis::Z])),
        "two Z nodes chained, base = the Constructed root"
    );
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap();
    assert!(
        (after.volume - before.volume).abs() < 1e-9,
        "volume invariant"
    );
    assert!((after.area - before.area).abs() < 1e-9, "area invariant");
    assert!(nacre_tess::to_obj(&m).is_ok());
    assert!(
        nacre_step::to_step(&m)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP")
    );
    assert!(solid_is_rotated(&m, c2), "re-rotated solid stays Rotated");
    // A boolean against the chain-rotated solid runs (guard retired, 3d-i): the axis-aligned
    // `d` inside the re-rotated `c2` is carved out, and the result is a valid solid.
    let d = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "chain-rotated cut is valid"
    );
}

/// Re-rotating about a different axis chains a node whose parent is the first
/// rotation (root → Z → X), `base` still the root.
#[test]
fn rerotate_different_axis_chains() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap().volume;
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c2),
        Some((2, vec![Axis::Z, Axis::X])),
        "chain root→Z→X"
    );
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap().volume;
    assert!((after - before).abs() < 1e-9, "volume invariant");
}

/// An **exact** (90°-family) rotation applied to an already-rotated solid is still
/// recorded as a chain node — the composite is inexact (an ancestor is), so the
/// forest must stay complete (1b silently dropped it). The solid stays Rotated and
/// boolean-rejected.
#[test]
fn rerotate_exact_after_inexact_records_node() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 90)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c2),
        Some((2, vec![Axis::Z, Axis::X])),
        "the exact 90°X is recorded as a chain node, not dropped"
    );
    assert!(
        solid_is_rotated(&m, c2),
        "composite is inexact → still Rotated"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A translation between two same-axis rotations forces a chain (this cell never
/// bundles anyway): the forest records both rotations, `base` stays the root, and
/// the result is rigid and valid — sound with no adjacency guard (each rotation is
/// its own node).
#[test]
fn rerotate_across_translation_chains() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap().volume;
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &translate_iso([5, -3, 2])).unwrap();
    m.rebuild_adjacency();
    let c3 = transform(&mut m, c2, &rot_iso(Axis::Z, 20)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c3),
        Some((2, vec![Axis::Z, Axis::Z])),
        "both rotations recorded; the intervening translation is not in the forest"
    );
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c3).unwrap().volume;
    assert!((after - before).abs() < 1e-9, "volume invariant");
}

/// A three-axis chain records three nodes root→Z→X→Y.
#[test]
fn rerotate_deep_chain() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
    m.rebuild_adjacency();
    let c3 = transform(&mut m, c2, &rot_iso(Axis::Y, 20)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c3),
        Some((3, vec![Axis::Z, Axis::X, Axis::Y])),
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A fresh rotation of a Constructed solid is unchanged from 1b (B0): an inexact
/// angle records a single root node; a 90°-family angle stays Constructed (no node,
/// boolean allowed). Guards that the B0/B1 split preserves fresh-rotation behavior.
#[test]
fn fresh_rotation_of_constructed_unchanged() {
    use nacre_scalar::Axis;
    // inexact → one root node, base = the cuboid's Constructed vertices.
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    assert_eq!(forest_probe(&m, c1), Some((1, vec![Axis::Z])));

    // exact 90° → Constructed (no node), boolean allowed.
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
    m.rebuild_adjacency();
    assert!(!solid_is_rotated(&m, c1), "fresh 90° stays Constructed");
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(forest_probe(&m, c1), None, "no rotation node");
}

/// Replay determinism (DNA 3): a re-rotation sequence reproduces the same forest
/// and handles.
///
/// The `assert_eq!` below compares handles minted by two *different* `Model`s, which is
/// legal only because `Handle`'s equality is its index — the very premise `replay` now
/// relies on. `tests/replay.rs` measures that premise directly instead of assuming it.
#[test]
fn rerotate_is_deterministic() {
    use nacre_scalar::Axis;
    let build = || {
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
        (bbox_lo(&m, c2), c2)
    };
    assert_eq!(build(), build(), "same ops → same geometry and handle");
}

// ---- coincident-coplanar merge (M5-c5) ----

fn stacked_cubes() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    (m, a, b)
}

#[test]
fn fuse_a_boss_onto_a_non_convex_solid() {
    // A contained boss on the top of an L-prism (non-convex kept `a`). The convexity gate
    // used to decline this to the seam path, which rejected the seamless contact; the
    // contained-coplanar Fuse now admits it. Volume 3 (L) + 0.4²·0.5 = 3.08.
    let (mut m, l) = l_prism(); // L footprint area 3, height 1
    let boss = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 1.0]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, l, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 3.08).abs() < 1e-12, "volume {vol}");
    // The boss top cap sits on the z = 1.5 plane, its outward normal +z.
    assert!(has_face_on_plane(
        &m,
        r,
        Point3::from_array([0.5, 0.5, 1.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    ));
}

#[test]
fn fuse_a_non_convex_profile_boss() {
    // An L-shaped boss (non-convex cutter `b`) on a cube top. The gate used to decline the
    // non-convex prism; the contained-coplanar Fuse now carries the L footprint as a hole.
    // Volume 1 (cube) + 0.12 (L area) · 0.4 = 1.048.
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let l_base: Vec<Point3> = [
        [0.3, 0.3],
        [0.7, 0.3],
        [0.7, 0.5],
        [0.5, 0.5],
        [0.5, 0.7],
        [0.3, 0.7],
    ]
    .iter()
    .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
    .collect();
    let (boss, _) = build_prism(
        &mut m,
        swept_world(l_base, Vector3::from_array([0.0, 0.0, 0.4])),
        vec![],
        Vector3::from_array([0.0, 0.0, 1.0]),
        None,
        None,
    )
    .unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, cube, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.048).abs() < 1e-12, "volume {vol}");
    // The L boss top cap sits on the z = 1.4 plane, its outward normal +z.
    assert!(has_face_on_plane(
        &m,
        r,
        Point3::from_array([0.4, 0.4, 1.4]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    ));
}

#[test]
fn cut_by_an_overhanging_boss_carrying_a_pin_owes_a_notch() {
    // Same seating, but the tool carries a pin reaching below the contact plane, so the plane
    // no longer separates the solids and the cut owes a real notch (1 − 0.2·0.2·0.5 = 0.98).
    //
    // This used to be an honest reject: the tool's z=1 cap is an annulus-like face whose
    // *inner* edge rides the pin's walls, and reading its occupancy off the ring's flank put
    // the material on the wrong side, so the class would not label. With the side read from
    // the ring's travel instead (`arrangement::run_body_above`), the notch comes out at the
    // hand-computed volume with a clean model.
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let pin = m.add_cuboid(
        Point3::from_array([0.55, 0.55, 0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    m.rebuild_adjacency();
    let tool = boolean_one(&mut m, BoolKind::Fuse, block, pin).unwrap();
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, tool).unwrap().volume - 1.02).abs() < 1e-12,
        "the pinned tool itself"
    );
    let notched = boolean_one(&mut m, BoolKind::Cut, base, tool).expect("the notch is buildable");
    m.rebuild_adjacency();
    assert!(
        (nacre_props::mass_props(&m, notched).unwrap().volume - 0.98).abs() < 1e-12,
        "the notch the tool owes: {}",
        nacre_props::mass_props(&m, notched).unwrap().volume
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "a notched result is still a clean model"
    );
}

// A1: plane-class canonicalization — coplanar walls of the two operands fold into one line.
#[test]
fn plane_classes_merge_a_shared_wall() {
    // Two unit cubes side by side share the plane x=1 (a's +x wall, b's -x wall — the same
    // plane, opposite normals). `plane_classes` must merge those two into one line class and
    // keep the far walls (a's x=0, b's x=2) distinct.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let planes_a = collect_planes(&m, a).unwrap();
    let na = planes_a.len();
    let mut planes = planes_a;
    planes.extend(collect_planes(&m, b).unwrap());
    // Find a plane by outward-normal x-sign and its x coordinate, within an index range.
    let find = |rng: std::ops::Range<usize>, nx: f64, x: f64| -> usize {
        rng.clone()
            .find(|&i| {
                let n = planes[i].plane().n_out.as_array();
                n[0] * nx > 0.5 && (planes[i].plane().tri[0].as_array()[0] - x).abs() < 1e-9
            })
            .expect("plane")
    };
    let a_xp = find(0..na, 1.0, 1.0); // a's +x wall at x=1
    let b_xm = find(na..planes.len(), -1.0, 1.0); // b's -x wall at x=1
    let a_xm = find(0..na, -1.0, 0.0); // a's -x wall at x=0
    let b_xp = find(na..planes.len(), 1.0, 2.0); // b's +x wall at x=2
    let canon = plane_classes(&crate::planes::test_judge(&planes));
    assert_eq!(canon[a_xp], canon[b_xm], "shared x=1 wall is one class");
    assert_ne!(canon[a_xm], canon[b_xp], "far walls stay distinct");
    assert_ne!(canon[a_xp], canon[a_xm], "x=1 and x=0 are different lines");
    // Every b face but its +x wall is coplanar with an a face, so 12 planes fold to 7 classes.
    let distinct: std::collections::HashSet<usize> = canon.iter().copied().collect();
    assert_eq!(distinct.len(), na + 1, "only b's far wall is a new class");
    // The class root is the smallest index in the class (deterministic canon).
    assert_eq!(canon[a_xp], a_xp.min(b_xm));
}

/// ★★★★ **Two faces of one judged surface are one class** (16-3). A nameless plane has no
/// name to intern classes by, so its class merging rests on the probes: both faces carry the
/// *same statement* → the same frame probes → identical chains → `shared_base` cancels the
/// motion and the exact predicate answers a **proved zero**. This is the argument turned into
/// a run: a straddle-datum prism is severed through its cap, leaving two faces on the one
/// judged surface, and `plane_classes` must fold them into one line.
#[test]
fn two_faces_of_one_judged_surface_are_one_class() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([3.3, 3.3, -1.0]),
        Point3::from_array([7.7, 7.7, 11.0]),
    );
    m.rebuild_adjacency();
    let OpOutput::Transform { solid: b } = crate::apply(
        &mut m,
        &Operation::Transform {
            solid: b,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: nacre_scalar::Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
    m.rebuild_adjacency();

    // A straddling vertex of the cut, plus two pure corners that share its plane sanely.
    let mut straddle = None;
    let mut pure = Vec::new();
    for &s in &cut {
        for &fh in &m.shells.get(m.solids.get(s).outer).faces {
            for lp in std::iter::once(&m.faces.get(fh).outer).chain(m.faces.get(fh).inner.iter()) {
                for &he in &lp.half_edges {
                    let vh = m.he_start(he);
                    let nacre_topo::VertexDef::ThreePlane(tri) = m.vertices.get(vh).def else {
                        continue;
                    };
                    let mine = m.plane_motion(tri[0]);
                    if tri.iter().any(|h| m.plane_motion(*h) != mine) {
                        straddle.get_or_insert(vh);
                    } else if !pure.contains(&vh) {
                        pure.push(vh);
                    }
                }
            }
        }
    }
    let straddle = straddle.expect("the cut leaves straddling corners");
    let vs = [straddle, pure[0], pure[1]];
    let OpOutput::DatumPlane { plane, frame } = crate::apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        },
    )
    .expect("the straddle datum") else {
        unreachable!()
    };
    assert!(!m.surface_name.contains_key(&plane), "nameless");
    let square = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([1.0, 0.0]),
        Point2::from_array([1.0, 1.0]),
        Point2::from_array([0.0, 1.0]),
    ])
    .unwrap();
    let OpOutput::Extrude { solid: prism, .. } = crate::apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square,
            dist: 0.5,
        },
    )
    .expect("prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    // Sever the prism through the middle with a thin box crossing the whole height: two
    // pieces, each with a base-cap face on the SAME judged surface.
    let (o, u, _v, w) = crate::rotated_vertex::frame_world_basis(
        &m,
        plane,
        &nacre_topo::FramePlacement::Canonical,
        frame.flip(),
    )
    .expect("the judged frame realizes");
    // A knife centred over the prism: origin + 0.5·û ± thickness, spanning v and w amply.
    let centre = Point3::from_array([
        o[0] + 0.5 * u[0] + 0.5 * w[0] * 0.0,
        o[1] + 0.5 * u[1],
        o[2] + 0.5 * u[2],
    ]);
    let knife = m.add_cuboid(
        centre + Vector3::from_array([-0.1, -5.0, -5.0]),
        centre + Vector3::from_array([0.1, 5.0, 5.0]),
    );
    m.rebuild_adjacency();
    let pieces = crate::boolean(&mut m, BoolKind::Cut, prism, knife).expect("sever");
    m.rebuild_adjacency();
    assert!(pieces.len() >= 2, "the knife must sever the prism");

    // The judgment table over the two pieces: their base-cap faces share the judged surface
    // and must fold into one class.
    let planes_a = collect_planes(&m, pieces[0]).unwrap();
    let na = planes_a.len();
    let mut planes = planes_a;
    planes.extend(collect_planes(&m, pieces[1]).unwrap());
    let on_datum: Vec<usize> = (0..planes.len())
        .filter(|&i| planes[i].surf() == plane)
        .collect();
    assert!(
        on_datum.len() >= 2
            && on_datum.iter().any(|&i| i < na)
            && on_datum.iter().any(|&i| i >= na),
        "both pieces must carry a face on the judged surface"
    );
    let canon = plane_classes(&crate::planes::test_judge(&planes));
    let first = canon[on_datum[0]];
    for &i in &on_datum[1..] {
        assert_eq!(
            canon[i], first,
            "two faces of one judged surface must be one class — same statement, same probes"
        );
    }
}

#[test]
fn overhang_boss_with_a_non_convex_footprint() {
    // An L-shaped (non-convex) boss footprint overhanging a cube edge. The old convexity-gated
    // overhang detector declined this, so it used to be an honest reject; the F2 dispatch
    // collapse hands it to the unified coplanar driver, which builds it exactly. Volume =
    // cube 1.0 + L-prism (area 0.9·0.2 + 0.3·0.2 = 0.24) · height 0.4 = 1.096.
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let l_base: Vec<Point3> = [
        [0.3, 0.3],
        [1.2, 0.3],
        [1.2, 0.5],
        [0.6, 0.5],
        [0.6, 0.7],
        [0.3, 0.7],
    ]
    .iter()
    .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
    .collect();
    let (l_tool, _) = build_prism(
        &mut m,
        swept_world(l_base, Vector3::from_array([0.0, 0.0, 0.4])),
        vec![],
        Vector3::from_array([0.0, 0.0, 1.0]),
        None,
        None,
    )
    .unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, cube, l_tool).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.096).abs() < 1e-12, "volume {vol}");
}

// ---- `unify_coplanar_faces`: the general (interface-free) coplanar merge, on hand-built
// `LocalFace` lists the coincident goldens above never reach — chains, opposite normals,
// holes, seam edges, and the asymmetric T-junction the global dissolve exists to prevent. ----

/// An axis-aligned `FaceInfo` at `d` along its normal, with a **non-degenerate `tri`** whose
/// right-hand normal is `n_out`. The merge reads more than `n_out` now — `loop_winding` and
/// `point_in_ring` name their arguments by plane and evaluate exact predicates on `tri` — so a
/// dummy triangle would make those answers meaningless.
fn mk_axis_plane(m: &mut Model, axis: usize, d: f64, positive: bool) -> WorkingPlane {
    let mut n = [0.0; 3];
    n[axis] = if positive { 1.0 } else { -1.0 };
    let normal = Vector3::from_array(n);
    let mut at = [0.0; 3];
    at[axis] = d;
    let origin = Point3::from_array(at);
    let plane = Plane::from_point_normal(origin, normal).unwrap();
    // Unregistered on purpose: these fixtures push one geometric plane as *two* handles
    // (`positive` both ways), which interning would collapse. The truth (the same tri
    // computed below, lifted) is still stated — nothing point-less enters the arena.
    let lift = |p: Point3| {
        p.as_array()
            .map(|x| nacre_scalar::Rat::from_decimal(x).unwrap())
    };
    let (ti, tj) = ((axis + 1) % 3, (axis + 2) % 3);
    let (ti, tj) = if positive { (ti, tj) } else { (tj, ti) };
    let stepr = |k: usize| {
        let mut q = at;
        q[k] += 1.0;
        Point3::from_array(q)
    };
    let surf = m.push_plane_unregistered(plane, [lift(origin), lift(stepr(ti)), lift(stepr(tj))]);
    let face = m.faces.push(Face {
        surface: surf,
        outer: Loop { half_edges: vec![] },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    // Two in-plane directions whose cross product is `+normal`, so `tri` winds outward.
    let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
    let (i, j) = if positive { (i, j) } else { (j, i) };
    let step = |k: usize| {
        let mut q = at;
        q[k] += 1.0;
        Point3::from_array(q)
    };
    let _ = face;
    let tri = [origin, step(i), step(j)];
    WorkingPlane {
        // A hand-built table has no recorded coefficients; the composed-rotation route
        // declines and the fixture takes the same escalating path it always did.
        base_rat: None,
        name_ints: None,
        base: crate::planes::BaseFrame::none(),
        surf,
        plane,
        tri,
        tri_pt3: tri.map(|p| nacre_cip::WitnessPoint::exact(p.as_array()).expect("exact")),
        rotated: false,
        frame_sign: 1, // `plane` is built from `normal`, so the two agree
        exact_coeffs: WorkingPlane::reconcile(&plane, tri, false).0,
        exact_normal: WorkingPlane::reconcile(&plane, tri, false).1,
    }
}

/// A `FaceInfo` for the `unify_coplanar_faces` tests, which read none of its geometry; the rest is a
/// valid-but-unreferenced dummy (`surf`/`face`/`plane` are never dereferenced there).
fn face(plane_idx: usize, nodes: Vec<Node>, inner: Vec<Vec<Node>>) -> LocalFace {
    LocalFace {
        surf: crate::planes::ClassIx::Plane(plane_idx),
        outer: crate::boolean::Bound::Ring(crate::boolean::Ring::from_clean_names(
            plane_idx, nodes,
        )),
        inner: inner
            .into_iter()
            .map(|r| {
                crate::boolean::Bound::Ring(crate::boolean::Ring::from_clean_names(plane_idx, r))
            })
            .collect(),
        flip: false,
    }
}

#[test]
fn unify_merges_a_coplanar_chain() {
    // Three unit squares on z=0 (+z), tiled in x, each sharing a vertical edge with the
    // next. One plane class, one normal ⇒ all fuse into a single face; the four
    // straight-angle mid-edge vertices dissolve, leaving one 4-corner rectangle.
    //
    // Named the way the arrangement names things — every vertex is the meeting of three
    // planes — because the straight-angle test reads those triples. (The old `Orig` fixture
    // exercised a path the engine stopped producing when it went all-`Seam`.)
    let mut m = Model::new();
    let p = vec![
        mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0, the shared class
        mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
        mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
        mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
        mk_axis_plane(&mut m, 0, 3.0, true),  // 4: x=3
        mk_axis_plane(&mut m, 1, 0.0, false), // 5: y=0
        mk_axis_plane(&mut m, 1, 1.0, true),  // 6: y=1
    ];
    let _canon: Vec<usize> = (0..p.len()).collect();
    let v = |x: usize, y: usize| Node::Seam([0, x, y]); // sorted: class, x-plane, y-plane
    let (c00, c10, c20, c30) = (v(1, 5), v(2, 5), v(3, 5), v(4, 5));
    let (c01, c11, c21, c31) = (v(1, 6), v(2, 6), v(3, 6), v(4, 6));
    let faces = vec![
        face(0, vec![c00, c10, c11, c01], vec![]),
        face(0, vec![c10, c20, c21, c11], vec![]),
        face(0, vec![c20, c30, c31, c21], vec![]),
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p)).unwrap();
    assert_eq!(out.len(), 1, "three coplanar faces fuse into one");
    let l = &out[0].outer.expect_ring();
    assert_eq!(l.len(), 4, "straight-angle mid vertices dissolved: {l:?}");
    for c in [c00, c30, c31, c01] {
        assert!(l.contains(&c), "corner kept");
    }
    for c in [c10, c20, c11, c21] {
        assert!(!l.contains(&c), "mid vertex dropped");
    }
}

#[test]
fn an_overhang_fuse_keeps_the_two_z1_caps_separate() {
    // An overhanging boss splits `z = 1` between two coplanar faces with **opposite** outward
    // normals — the base's exposed top (`+z`) and the boss underside (`-z`). They must not be
    // fused into one face: their `flip` differs, so `unify`'s `(plane_idx, flip)` group key
    // keeps them apart.
    //
    // This replaces two retired tests (`unify_keeps_opposite_normal_coplanar`,
    // `unify_keeps_holed_faces`) that built `Node::Orig` faces to exercise a passthrough the
    // arrangement never triggers — it emits all-`Seam`. The real invariant is exercised here on
    // the production path, in the default `cargo test` run: `overhang_fuse_then_cut_matches_occt`
    // proves it against OCCT but is `#[ignore]`, so this hand-computed volume is the non-ignored
    // guard. A wrong merge collapses the topology — the volume shifts or `validate` speaks.
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!(
        (vol - 1.5).abs() < 1e-12,
        "base 1 + boss 0.5, no overlap: {vol}"
    );
}

#[test]
fn a_hole_filled_by_two_faces_still_merges() {
    // Replaces `unify_skips_seam_shared_edges`, whose premise is gone twice over: the `Orig`
    // gate it pinned was removed when the engine went all-`Seam`, and `splice_along`, whose
    // panic it guarded against, no longer exists.
    //
    // What matters now is that the merge is not special-cased to "a hole filled by exactly one
    // neighbour". A [0,3]² face with a [1,2]² hole, and that hole filled by **two** pieces split
    // at x=1.5: every ring edge between them is carried in both directions, so erasing interior
    // boundary leaves only the outer square — one face, no hole, whatever the filling is cut
    // into. This is the case that separates a general rule from a bespoke one.
    let mut m = Model::new();
    let p = vec![
        mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
        mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
        mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
        mk_axis_plane(&mut m, 0, 1.5, true),  // 3: x=1.5, where the filling is split
        mk_axis_plane(&mut m, 0, 2.0, true),  // 4: x=2
        mk_axis_plane(&mut m, 0, 3.0, true),  // 5: x=3
        mk_axis_plane(&mut m, 1, 0.0, false), // 6: y=0
        mk_axis_plane(&mut m, 1, 1.0, true),  // 7: y=1
        mk_axis_plane(&mut m, 1, 2.0, true),  // 8: y=2
        mk_axis_plane(&mut m, 1, 3.0, true),  // 9: y=3
    ];
    let _canon: Vec<usize> = (0..p.len()).collect();
    let v = |x: usize, y: usize| Node::Seam([0, x, y]);
    let (o00, o30, o33, o03) = (v(1, 6), v(5, 6), v(5, 9), v(1, 9));
    let (h11, h12, h22, h21) = (v(2, 7), v(2, 8), v(4, 8), v(4, 7));
    let (m12, m11) = (v(3, 8), v(3, 7)); // the split points on the hole's top and bottom
    let faces = vec![
        // Outer square with the hole, wound the way `emit_faces` states it: outer CCW, hole CW.
        face(
            0,
            vec![o00, o30, o33, o03],
            vec![vec![h11, h12, m12, h22, h21, m11]],
        ),
        face(0, vec![h11, m11, m12, h12], vec![]), // left filler
        face(0, vec![m11, h21, h22, m12], vec![]), // right filler
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p)).unwrap();
    assert_eq!(out.len(), 1, "the hole is filled, so one face remains");
    assert!(out[0].inner.is_empty(), "and it has no hole left");
    assert_eq!(out[0].outer.expect_ring().len(), 4, "just the outer square");
    for c in [o00, o30, o33, o03] {
        assert!(out[0].outer.expect_ring().contains(&c), "outer corner kept");
    }
}

#[test]
fn unify_keeps_a_vertex_that_is_a_corner_elsewhere() {
    // F0,F1 on z=0 merge; their shared-edge endpoint (1,0,0) is a straight angle on the
    // merged face but a real corner on a perpendicular face G (plane y=0). Global degree 3
    // ⇒ it is NOT dissolved — a per-face local rule would have, opening a T-junction. Its
    // twin (1,1,0), on the merged face only, IS dissolved.
    let mut m = Model::new();
    let p = vec![
        mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
        mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
        mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
        mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
        mk_axis_plane(&mut m, 1, 0.0, false), // 4: y=0, the perpendicular face's plane
        mk_axis_plane(&mut m, 1, 1.0, true),  // 5: y=1
        mk_axis_plane(&mut m, 2, 1.0, true),  // 6: z=1
    ];
    let _canon: Vec<usize> = (0..p.len()).collect();
    let (v000, v100, v200) = (
        Node::Seam([0, 1, 4]),
        Node::Seam([0, 2, 4]),
        Node::Seam([0, 3, 4]),
    );
    let (v010, v110, v210) = (
        Node::Seam([0, 1, 5]),
        Node::Seam([0, 2, 5]),
        Node::Seam([0, 3, 5]),
    );
    let (v101, v201) = (Node::Seam([2, 4, 6]), Node::Seam([3, 4, 6]));
    let faces = vec![
        face(0, vec![v000, v100, v110, v010], vec![]),
        face(0, vec![v100, v200, v210, v110], vec![]),
        face(4, vec![v200, v100, v101, v201], vec![]), // perpendicular, not coplanar
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p)).unwrap();
    assert_eq!(out.len(), 2, "z=0 pair merges; G stays");
    let merged = out.iter().find(|lf| lf.surf.plane() == 0).unwrap();
    assert!(
        merged.outer.expect_ring().contains(&v100),
        "corner-elsewhere vertex kept (no T-junction)"
    );
    assert!(
        !merged.outer.expect_ring().contains(&v110),
        "pure straight-angle vertex dropped"
    );
}

proptest! {
    /// Diagonal corner overlaps (clean seam): fuse/cut volumes match the
    /// independent AABB formula (not nacre's own common).
    #[test]
    fn fuse_cut_diagonal_boxes_match_aabb(
        amin in prop::array::uniform3(-3.0f64..3.0),
        aext in prop::array::uniform3(1.0f64..3.0),
        t in prop::array::uniform3(0.15f64..0.6),
        s in prop::array::uniform3(0.3f64..2.0),
    ) {
        let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
        let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
        let bmax: [f64; 3] = std::array::from_fn(|i| amax[i] + s[i]);
        let ov: f64 = (0..3).map(|i| amax[i] - bmin[i]).product();
        let va: f64 = aext.iter().product();
        let vb: f64 = (0..3).map(|i| bmax[i] - bmin[i]).product();

        let build = || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
            let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
            (m, a, b)
        };

        let (mut m1, a1, b1) = build();
        let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1);
        prop_assume!(rf.is_ok()); // skip rare coplanar/degenerate configs
        let rf = rf.unwrap();
        m1.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m1).is_empty());
        let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
        prop_assert!((vf - (va + vb - ov)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

        let (mut m2, a2, b2) = build();
        let rc = boolean_one(&mut m2, BoolKind::Cut, a2, b2).unwrap();
        m2.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m2).is_empty());
        let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;
        prop_assert!((vc - (va - ov)).abs() <= 1e-9 * va, "cut {vc}");
    }

    /// Matched-footprint stacked boxes (coincident z-interface): fuse volume
    /// is the sum, cut is A, common is empty.
    #[test]
    fn stacked_boxes_merge_volumes(
        x0 in -3.0f64..3.0,
        y0 in -3.0f64..3.0,
        dx in 0.5f64..3.0,
        dy in 0.5f64..3.0,
        z0 in -3.0f64..3.0,
        h1 in 0.5f64..3.0,
        h2 in 0.5f64..3.0,
    ) {
        let (x1, y1) = (x0 + dx, y0 + dy);
        let (zm, z1) = (z0 + h1, z0 + h1 + h2);
        let build = || {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([x0, y0, z0]), Point3::from_array([x1, y1, zm]));
            let b = m.add_cuboid(Point3::from_array([x0, y0, zm]), Point3::from_array([x1, y1, z1]));
            (m, a, b)
        };
        let (va, vb) = (dx * dy * h1, dx * dy * h2);

        let (mut m1, a1, b1) = build();
        // Not `prop_assume!`: a matched-footprint stack is squarely in coverage whatever the
        // dimensions are, so a reject here is a defect, not an uninteresting sample. Assuming
        // it away is how this property went on passing while the kernel aborted on 2% of the
        // space and rejected 95% of it (family #3).
        let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1)
            .expect("stacked boxes fuse at any dimensions");
        m1.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m1).is_empty());
        let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
        prop_assert!((vf - (va + vb)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

        let (mut m2, a2, b2) = build();
        // The stack shares only its interface plane ⇒ no volume in common, at any dimensions.
        prop_assert!(boolean(&mut m2, BoolKind::Common, a2, b2).unwrap().is_empty());

        let (mut m3, a3, b3) = build();
        let rc = boolean_one(&mut m3, BoolKind::Cut, a3, b3).unwrap();
        let vc = nacre_props::mass_props(&m3, rc).unwrap().volume;
        prop_assert!((vc - va).abs() <= 1e-9 * va, "cut {vc}");
    }
}

/// One geometric plane is one class **whatever the two faces' sizes**.
///
/// ★ **The reason changed, and that is the news.** This used to assert that the coefficient
/// test *could not* prove these coplanar — two walls of one plane at different face sizes have
/// un-normalized 4-vectors that are not exactly proportional — and that the coordinate branch
/// was what earned the merge. Surfaces are interned on their canonical rational coefficients
/// now, so the two walls are handed **one handle**, and the merge is a handle comparison
/// before any geometry is asked. The f64 non-proportionality is still real and still pinned,
/// on the planes themselves, in `nacre_topo`'s
/// `two_faces_of_one_plane_disagree_in_f64_and_agree_in_the_rationals`.
#[test]
fn one_plane_is_one_class_whatever_the_face_size() {
    let (dx, dy) = (1.628165457453874f64, 0.5f64);
    let (z0, h1, h2) = (0.11200046228159026f64, 0.5f64, 2.07926124157585f64);
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, z0]),
        Point3::from_array([dx, dy, z0 + h1]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, z0 + h1]),
        Point3::from_array([dx, dy, z0 + h1 + h2]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: _planes,
        plane_ix,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    // The two `+X` walls: same plane x = dx, different face sizes (heights h1 vs h2). Two
    // *faces*, so this searches the face table — the plane table holds one entry for both,
    // which is the property under test.
    let x_walls: Vec<usize> = (0..faces_tab.len())
        .filter(|&i| {
            faces_tab[i].plane().n_out.as_array() == [1.0, 0.0, 0.0]
                && (faces_tab[i].plane().tri[0].as_array()[0] - dx).abs() < 1e-12
        })
        .collect();
    assert_eq!(x_walls.len(), 2, "one wall from each box: {x_walls:?}");
    let (i, j) = (x_walls[0], x_walls[1]);
    assert_eq!(
        faces_tab[i].surf(),
        faces_tab[j].surf(),
        "two faces, one surface — interning collapsed them at construction"
    );
    assert!(
        crate::planes::test_judge(&faces_tab).planes_coplanar(i, j),
        "and the geometry agrees, so nothing rests on the handle alone"
    );
    assert_eq!(plane_ix[i], plane_ix[j], "so they are one plane-table row");
}

/// ★★★ **Solving a vertex's three planes lands on its coordinate.**
///
/// This is the premise the whole of `docs/truth-and-cache.md` rests on: a point *is* the
/// meeting of three surfaces, and the stored coordinate is a rounded answer to that question.
/// Stage 0 measured it by reading the census; this asserts it on data the kernel itself built,
/// which is a different claim — the definitions have to be *right*, not merely present.
///
/// ★ **Split by origin, because one row proves nothing.** A `Discovered` vertex's coordinate
/// was produced by solving exactly this triple, so agreement there is an identity. The rows
/// that carry weight are `Constructed` and `Moved`, where the coordinate came from somewhere
/// else entirely — construction arithmetic, or a motion replayed on a base point.
#[test]
fn a_vertex_definition_solves_to_its_own_coordinate() {
    // ★ **Three solids, because one would not exercise three kinds of vertex.** The first
    // draft of this test used a fused-then-turned-then-cut solid and reported
    // `Constructed (0,0) Discovered (24,24) Moved (0,0)` with a worst error of exactly zero —
    // a boolean recomputes every vertex it emits, so the only row present was the tautological
    // one. A plain box keeps its constructed corners; a turned box keeps them as `Moved`.
    let mut m = Model::new();
    let plain = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let to_turn = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let OpOutput::Transform { solid: turned } = apply(
        &mut m,
        &Operation::Transform {
            solid: to_turn,
            isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                axis: nacre_scalar::Axis::Z,
                point: [nacre_scalar::Rat::from_int(0); 3],
                angle: nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!("transform yields Transform output")
    };
    m.rebuild_adjacency();
    let a = m.add_cuboid(
        Point3::from_array([10.0, 0.0, 0.0]),
        Point3::from_array([12.0, 3.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([11.0, 1.0, 0.5]),
        Point3::from_array([14.0, 2.0, 2.5]),
    );
    m.rebuild_adjacency();
    let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse")[0];
    m.rebuild_adjacency();

    // Rows by **producer family** (S7: `Origin`'s three tags are gone, and the families
    // they used to stand in for are exactly these three fixtures — a plainly built box, a
    // rotated one, a boolean result).
    let mut counts = [(0usize, 0usize); 3]; // (with a three-plane definition, total)
    let mut worst = [0.0f64; 3];
    let mut diam = 0.0f64;
    for (kind, solid) in [plain, turned, fused].into_iter().enumerate() {
        let shell = m.solids.get(solid).outer;
        let mut seen: Vec<Handle<Vertex>> = Vec::new();
        for &fh in &m.shells.get(shell).faces {
            let face = m.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in m.edges.get(he.edge).vertices.iter() {
                        if !seen.contains(&vh) {
                            seen.push(vh);
                        }
                    }
                }
            }
        }
        for &vh in &seen {
            let coord = m.vertex_point(vh);
            for c in coord.as_array() {
                diam = diam.max(c.abs());
            }
            counts[kind].1 += 1;
            let VertexDef::ThreePlane(planes) = m.vertices.get(vh).def else {
                continue; // a seam vertex names a curve, not a point — no solve to check
            };
            counts[kind].0 += 1;
            let coeffs = planes.map(|s| match m.surface(s) {
                nacre_geom::Surface::Plane(p) => p.coefficients(),
                nacre_geom::Surface::Cylinder(_) => panic!("ThreePlane named a cylinder"),
            });
            let solved = solve_three_planes(coeffs).expect("three planes meeting at a point");
            let d = solved
                .iter()
                .zip(coord.as_array())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            worst[kind] = worst[kind].max(d);
        }
    }

    let limit = diam * 2f64.powi(-40);
    eprintln!(
        "[definition coverage] plain {:?} worst {:e} | turned {:?} worst {:e} | \
             fused {:?} worst {:e} | limit {:e}",
        counts[0], worst[0], counts[1], worst[1], counts[2], worst[2], limit
    );
    for (i, name) in ["plain", "turned", "fused"].iter().enumerate() {
        assert!(
            counts[i].1 > 0,
            "{name} is not exercised — the row proves nothing"
        );
        assert_eq!(
            counts[i].0, counts[i].1,
            "{name} vertices without a definition: {:?}",
            counts[i]
        );
        assert!(
            worst[i] <= limit,
            "{name}: a definition solved {:e} away from its own coordinate (limit {limit:e})",
            worst[i]
        );
    }
    // ★ The **fused** row is the tautological one — a boolean vertex's coordinate *came
    // from* solving this very triple, so a zero there says nothing. `plain` and `turned`
    // are the claims. (Same argument as before S7; the row moved from a tag to a producer.)
    assert_eq!(worst[2], 0.0, "a boolean vertex is its own solve, exactly");
}

/// The negative control for the assertion above: point a definition at the wrong plane and the
/// solve must land somewhere else. Without this, an agreement test passes on any model whose
/// planes happen to be near each other.
#[test]
fn a_wrong_plane_in_a_definition_is_caught() {
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    let shell = m.solids.get(s).outer;
    let surfaces: Vec<_> = m
        .shells
        .get(shell)
        .faces
        .iter()
        .map(|&fh| m.faces.get(fh).surface)
        .collect();
    // Take a real corner's triple and swap one plane for another face of the same box.
    let corner = m
        .shells
        .get(shell)
        .faces
        .first()
        .map(|&fh| {
            let face = m.faces.get(fh);
            m.edges.get(face.outer.half_edges[0].edge).vertices
        })
        .expect("a face with a loop")[0];
    let VertexDef::ThreePlane(mut planes) = m.vertices.get(corner).def else {
        panic!("a constructed corner has a definition")
    };
    let good = solve_three_planes(planes.map(|h| match m.surface(h) {
        nacre_geom::Surface::Plane(p) => p.coefficients(),
        nacre_geom::Surface::Cylinder(_) => unreachable!(),
    }))
    .expect("meets at a point");
    // Any face not already in the triple. Every one of them moves the point (or leaves the
    // three not meeting at all, which is just as good a refutation).
    let other = *surfaces
        .iter()
        .find(|h| !planes.contains(h))
        .expect("a fourth face");
    planes[0] = other;
    let bad = solve_three_planes(planes.map(|h| match m.surface(h) {
        nacre_geom::Surface::Plane(p) => p.coefficients(),
        nacre_geom::Surface::Cylinder(_) => unreachable!(),
    }));
    let moved = bad.is_none_or(|bad| {
        good.iter()
            .zip(bad)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max)
            > 0.1
    });
    assert!(
        moved,
        "swapping a plane must move the solved point, or the assertion above is vacuous"
    );
}

/// Three planes by Cramer, or `None` when they do not meet in a point.
fn solve_three_planes(p: [[f64; 4]; 3]) -> Option<[f64; 3]> {
    let n = |i: usize| [p[i][0], p[i][1], p[i][2]];
    let (a, b, c) = (n(0), n(1), n(2));
    let det3 = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det3([a, b, c]);
    if d == 0.0 {
        return None;
    }
    let rhs = [-p[0][3], -p[1][3], -p[2][3]];
    Some(core::array::from_fn(|j| {
        let mut m = [a, b, c];
        for (row, r) in m.iter_mut().zip(rhs) {
            row[j] = r;
        }
        det3(m) / d
    }))
}

/// ★ **A collinear midpoint is dissolved at construction, so the prism it builds is its
/// clean twin's, bit for bit — and every corner keeps a three-plane definition** (S3; the
/// last piece of `docs/truth-and-cache.md` Q2 ②).
///
/// Before S3 this profile built *seven* faces whose split bottom edge's two walls interned to
/// one surface handle — a vertex named `[S, S, cap]`, a line and not a point, `definition:
/// None`. The constructor now deletes the flat corner (a lossless normalization: the shape is
/// identical), so that population cannot reach the topology at all.
#[test]
fn a_collinear_midpoint_profile_builds_its_clean_twin_bit_for_bit() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let build = |profile: Profile2d| {
        let mut m = Model::new();
        let __w3 = SketchFrame::world(&m, Axis::Z);
        apply(
            &mut m,
            &Operation::Extrude {
                frame: __w3,
                profile,
                dist: 1.0,
            },
        )
        .expect("extrude");
        m
    };
    // A unit square whose bottom edge carries a redundant midpoint — and the square itself.
    let split = build(
        Profile2d::polygon(vec![
            p(0.0, 0.0),
            p(0.5, 0.0), // collinear with its neighbours — dissolved at construction
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
        ])
        .unwrap(),
    );
    let clean = build(
        Profile2d::polygon(vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]).unwrap(),
    );
    // Bit-for-bit the same model: same counts, same coordinates, same surface wiring.
    assert_eq!(split.faces.len(), clean.faces.len(), "6 faces, not 7");
    assert_eq!(split.vertices.len(), clean.vertices.len());
    for ((sv, _), (cv, _)) in split.vertices.iter().zip(clean.vertices.iter()) {
        assert_eq!(
            split.vertex_point(sv).as_array(),
            clean.vertex_point(cv).as_array(),
            "coordinates"
        );
    }
    for ((_, s), (_, c)) in split.faces.iter().zip(clean.faces.iter()) {
        assert_eq!(s.surface, c.surface, "surface wiring");
    }
    // And the corner population is whole: every vertex holds a three-plane definition.
    for (_, v) in split.vertices.iter() {
        assert!(
            matches!(v.def, VertexDef::ThreePlane(_)),
            "a corner without a three-plane definition survived: {v:?}"
        );
    }
}

/// ★ **Two *separated* collinear walls still intern to one surface** — the re-pin of what
/// `a_collinear_profile_vertex_gives_its_two_walls_one_surface` used to hold.
///
/// The dissolve pass only deletes flat corners (adjacent same-plane walls); two edges of a
/// notched profile lying on one line are legitimate geometry, their walls are two statements
/// of one plane, and interning makes them one handle (`docs/truth-and-cache.md`, 「남은 것」 3).
/// Adjacent walls can no longer collide, so this is where "same plane = same handle" stays
/// pinned.
#[test]
fn two_separated_collinear_walls_intern_to_one_surface() {
    let mut m = Model::new();
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    // A right-edge notch: the two vertical segments at `x = 4` are collinear, not adjacent.
    let profile = Profile2d::polygon(vec![
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 1.0),
        p(3.0, 1.0),
        p(3.0, 2.0),
        p(4.0, 2.0),
        p(4.0, 3.0),
        p(0.0, 3.0),
    ])
    .unwrap();
    // Self-qualification: the dissolve pass must have left all eight corners standing.
    assert_eq!(profile.outer().points().len(), 8, "no corner is flat");
    let __w2 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w2,
            profile,
            dist: 1.0,
        },
    )
    .expect("extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    let shell = m.solids.get(solid).outer;
    let surfaces: Vec<_> = m
        .shells
        .get(shell)
        .faces
        .iter()
        .map(|&fh| m.faces.get(fh).surface)
        .collect();
    // Eight profile points ⇒ eight wall quads, plus two caps.
    assert_eq!(surfaces.len(), 10, "eight walls and two caps");
    let mut distinct = surfaces.clone();
    distinct.sort_unstable_by_key(|h| h.index());
    distinct.dedup();
    assert_eq!(
        distinct.len(),
        9,
        "the notch's two x=4 walls are one plane, so one handle: {surfaces:?}"
    );
}

/// A vertex where one plane is split between two faces is named by **the planes that touch it**,
/// not by the loop's two neighbours.
///
/// The overhang chain puts the base's exposed top and the cantilever's underside on one plane
/// (`z = 1`, opposite normals — `unify` rightly keeps them apart, `canon` rightly calls them one
/// class). At `(1, 0.25, 1)` the side wall's loop runs straight through their shared line, so the
/// old rule named that vertex with the same plane twice: a triple defining no point, which the
/// exact predicates — whose precondition is `D ≠ 0` — aborted on. Measured 2026-07-22 as the only
/// path a degenerate triple reached them (16 arrivals in the OCCT suite, now 0).
#[test]
fn a_vertex_is_named_by_the_planes_that_touch_it() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let cutter = m.add_cuboid(
        Point3::from_array([1.1, 0.35, 0.5]),
        Point3::from_array([1.4, 0.65, 2.5]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, overhung, cutter).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let _ = &faces_tab;
    let mut checked = 0usize;
    for sh in solid_shell_handles(&m, overhung) {
        for &fh in &m.shells.get(sh).faces {
            let p = surf_ix[&fh];
            let tris = combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix)
                .unwrap()
                .poly()
                .expect("a poly outer")
                .triples
                .clone();
            for t in &tris {
                // The triple is already dense plane ids: distinct means three real planes.
                assert!(
                    t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                    "vertex triple {t:?} names one plane twice"
                );
                // A name that denotes three distinct classes must denote a real point.
                assert!(
                    three_planes(
                        &planes[t[0]].plane,
                        &planes[t[1]].plane,
                        &planes[t[2]].plane
                    )
                    .is_some(),
                    "triple {t:?} defines no point"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "the chained operand has vertices to name");
}

/// **Dense plane ids are order-isomorphic to the sparse roots.** `dense_planes` ranks the class
/// roots, so any comparison, sort or lex-min over plane indices reads the same either way.
///
/// This is a **migration gate, not a permanent invariant**: it exists so the claim is measured
/// before the split rides on it, and it retires with `canon` — its subject, not its coverage,
/// is what goes away.
#[test]
fn dense_plane_ids_are_monotone_in_canon() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    // An overhanging boss splits `z = 1` between two faces, so classes really do merge and the
    // ranking really does compress — without that the map is the identity and proves nothing.
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let probe = m.add_cuboid(
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    // Rebuild the pieces `dense_planes` consumes, so this locks its contract without needing
    // `canon` to escape `plane_index_setup`. `plane_classes` is the same union-find the setup
    // runs; `dense_planes` the same ranking.
    let mut faces = collect_planes(&m, chained).unwrap();
    faces.extend(collect_planes(&m, probe).unwrap());
    let canon = plane_classes(&crate::planes::test_judge(&faces));
    let (geom, plane_ix, _cyls) = dense_planes(&faces, &canon);
    assert!(
        canon.iter().enumerate().any(|(i, &c)| c != i),
        "fixture has no split plane — the invariant would be vacuous"
    );
    assert!(
        geom.len() < canon.len(),
        "the ranking must actually compress"
    );
    for i in 0..canon.len() {
        for j in 0..canon.len() {
            assert_eq!(
                canon[i].cmp(&canon[j]),
                plane_ix[i].plane().cmp(&plane_ix[j].plane()),
                "faces {i}/{j}: canon {}/{} vs dense {}/{}",
                canon[i],
                canon[j],
                plane_ix[i].plane(),
                plane_ix[j].plane()
            );
        }
    }
}

/// **A plane triple is always in class form.** `planes` is a per-face table, so the same
/// `usize` could mean "face" or "plane"; producers settle it by emitting class roots, and a
/// consumer's raw `==` then means "same plane". Four silent-wrong bugs on this branch came from
/// the two meanings meeting in one comparison, so the invariant is asserted, not assumed.
#[test]
fn plane_triples_are_always_canon() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    // An *overhanging* boss splits `z = 1` between two faces with opposite normals (the base's
    // exposed top and the boss underside) — the shape that makes "face index" and "plane index"
    // differ at all. A boss sitting wholly inside the top merges into one holed face instead.
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let probe = m.add_cuboid(
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, chained, probe).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    // The fixture must actually merge two faces into one plane, or this proves nothing.
    assert!(
        planes.len() < faces_tab.len(),
        "fixture has no split plane — the invariant would be vacuous"
    );
    // A producer hands out dense plane ids (`loop_triples` maps face indices through
    // `plane_ix`), so "every element is a plane, not a face" is now the type, not a runtime
    // check. What remains testable is that the ids are in range and sorted-distinct.
    let mut checked = 0usize;
    for sh in solid_shell_handles(&m, chained) {
        for &fh in &m.shells.get(sh).faces {
            let p = surf_ix[&fh];
            let mut rings = vec![
                combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix)
                    .unwrap()
                    .poly()
                    .expect("a poly outer")
                    .triples
                    .clone(),
            ];
            rings.extend(
                combinatorics::hole_rings(&m, fh, p, &inc_a, &jd, &plane_ix)
                    .unwrap()
                    .into_iter()
                    .filter_map(|lr| lr.poly().map(|nr| nr.triples.clone())),
            );
            for t in rings.iter().flatten() {
                for &k in t {
                    assert!(
                        k < planes.len(),
                        "triple {t:?} names {k}, out of the plane table"
                    );
                }
                assert!(
                    t[0] < t[1] && t[1] < t[2],
                    "triple {t:?} is not three distinct planes in sorted order"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "the chained operand has vertices to name");
}

/// **A point reads zero on each of its three defining planes.** `Judge::orient3d`'s on-plane
/// shortcut is a raw `==` against the triple, so a vertex on the query plane must name it by the
/// same id the query uses. The face/plane split makes that automatic — a plane has exactly one
/// id now, so the old failure (a vertex named by face 6 of the `z = 1` class invisible to a
/// query about face 1 of it) cannot be expressed. What is left to check is the identity itself.
#[test]
fn a_vertex_on_the_cut_plane_reads_zero_whichever_face_names_it() {
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    let probe = m.add_cuboid(
        Point3::from_array([1.1, 0.35, 0.5]),
        Point3::from_array([1.4, 0.65, 2.5]),
    );
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, chained, probe).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    assert!(
        planes.len() < faces_tab.len(),
        "fixture has no split plane — the sibling faces this used to distinguish"
    );
    let mut on_plane = 0usize;
    for sh in solid_shell_handles(&m, chained) {
        for &fh in &m.shells.get(sh).faces {
            let p = surf_ix[&fh];
            let tris = combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix)
                .unwrap()
                .poly()
                .expect("a poly outer")
                .triples
                .clone();
            for t in &tris {
                // The vertex lies on exactly its three defining planes; each must read 0.
                for &q in t {
                    assert_eq!(
                        combinatorics::side_of(&jd, *t, q),
                        0,
                        "vertex {t:?} lies on plane {q} but does not read 0"
                    );
                    on_plane += 1;
                }
            }
        }
    }
    assert!(on_plane > 0, "some vertex lies on some queried plane");
}

// ---- boolean Common algorithm (M5-c3 commit 2) ----

#[test]
fn common_rejects_non_planar_input() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let cyl = m.add_cylinder(
        Point3::from_array([1.0, 1.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    assert_rejects(
        || boolean_one(&mut m, BoolKind::Common, a, cyl),
        RejectReason::SeatedCylinderCap,
    );
}

proptest! {
    /// Overlapping axis-aligned boxes: the intersection volume equals the
    /// independent AABB-overlap product (mixed A/B axis-aligned vertices).
    #[test]
    fn common_axis_boxes_volume_matches_aabb_overlap(
        amin in prop::array::uniform3(-5.0f64..5.0),
        aext in prop::array::uniform3(1.0f64..4.0),
        t in prop::array::uniform3(0.05f64..0.7),
        bext in prop::array::uniform3(1.0f64..4.0),
    ) {
        let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
        let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
        let bmax: [f64; 3] = std::array::from_fn(|i| bmin[i] + bext[i]);
        let expected: f64 = (0..3)
            .map(|i| (amax[i].min(bmax[i]) - bmin[i]).max(0.0))
            .product();
        prop_assume!(expected > 1e-3);

        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
        let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
        let res = boolean_one(&mut m, BoolKind::Common, a, b);
        prop_assume!(res.is_ok()); // skip rare coplanar/degenerate configs
        let r = res.unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        prop_assert!((vol - expected).abs() <= 1e-9 * expected.max(1.0), "{vol} vs {expected}");
    }

    /// A tilted square prism (oblique planes) intersected with a big enclosing
    /// box is the prism — exercises non-axis-aligned face normals (the in/out
    /// sign and CCW ordering) with an independent oracle (the prism's own mass).
    #[test]
    fn common_tilted_prism_with_enclosing_box_is_the_prism(
        nx in -0.5f64..0.5,
        ny in -0.5f64..0.5,
    ) {
        let plane = SketchPlane::from_origin_normal(
            Point3::origin(),
            Vector3::from_array([nx, ny, 1.0]),
        )
        .unwrap();
        let __plane = plane;
        let mut __scratch203 = Model::new();
        let __g202 = datum_frame(&mut __scratch203, __plane);
        let mut m = replay(&[
            Operation::DatumPlane { def: DatumDef::Stated(__plane) },
            Operation::Extrude {
                frame: __g202,
            profile: square(),
            dist: 1.0,
        }])
        .unwrap();
        let prism = *m.live_solids.first().unwrap();
        let vol_prism = nacre_props::mass_props(&m, prism).unwrap().volume;
        let c = m.add_cuboid(Point3::from_array([-10.0; 3]), Point3::from_array([10.0; 3]));
        let res = boolean_one(&mut m, BoolKind::Common, prism, c);
        prop_assume!(res.is_ok());
        let r = res.unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        let vol_r = nacre_props::mass_props(&m, r).unwrap().volume;
        prop_assert!(
            (vol_r - vol_prism).abs() <= 1e-9 * vol_prism.max(1.0),
            "{vol_r} vs {vol_prism}"
        );
    }
}
/// ★★★★ **Padding the same footprint twice must not leave zero-area faces.**
///
/// It did, and `validate` reported nothing. Padding `1.1` and then `6.6` on the top of a
/// cuboid left two faces of area `2.2e-16` on the plane `z = 2.1`, whose long edges sat one
/// ulp apart (`-0.6` against `-0.6000000000000001`).
///
/// The ulp came from the sketch frame's **origin**. `face_frame` takes it from the face's area
/// centroid, and `ring_area_centroid` was rounding twice more than it needed to, so the first
/// pad's top face reported its centre as `-1.11e-16` instead of `0`. The second pad then placed
/// the same profile one ulp away from where the first had placed it, and the kernel — correctly
/// — built the one-ulp-wide faces that answer describes.
///
/// ★ `SketchPlane::exact()` does not catch this: it checks only that the axes are orthonormal,
/// never the origin.
#[test]
fn padding_one_footprint_twice_leaves_no_zero_area_face() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let rect = || {
        Profile2d::polygon(vec![p(-1.0, -0.6), p(1.0, -0.6), p(1.0, 0.6), p(-1.0, 0.6)]).unwrap()
    };
    let mut m = Model::new();
    let mut solid = m.add_cuboid(
        Point3::from_array([-2.0, -2.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    m.rebuild_adjacency();
    // The topmost face pointing up. Everything here is axis-aligned, so this is unambiguous.
    let top = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        let shell = m.solids.get(s).outer;
        *m.shells
            .get(shell)
            .faces
            .iter()
            .max_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    crate::ops::face_plane(m, f)
                        .ok()
                        .filter(|sp| sp.normal().as_array()[2] > 0.5)
                        .map(|sp| sp.origin.as_array()[2])
                        .unwrap_or(f64::NEG_INFINITY)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a face")
    };
    for dist in [1.1, 6.6] {
        let face = top(&m, solid);
        let OpOutput::PadOnFace { solid: out, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face,
                profile: rect(),
                dist,
            },
        )
        .expect("pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid = out;
    }
    let shell = m.solids.get(solid).outer;
    let degenerate: Vec<_> = m
        .shells
        .get(shell)
        .faces
        .iter()
        .filter_map(|&f| {
            let area = nacre_props::face_props(&m, f).ok()?.area;
            (area < 1e-9).then_some((f, area))
        })
        .collect();
    assert!(
        degenerate.is_empty(),
        "zero-area faces survived: {degenerate:?}"
    );
    // A plain box with one rib on top: 6 + 5 walls/cap, no leftovers from the seam.
    assert_eq!(m.shells.get(shell).faces.len(), 11);
}
/// ★★★★★ **The target: two ways of reaching one height land on one plane, far from the
/// origin.**
///
/// The frame's origin is where the drift used to enter — `face_frame` took it from the face's
/// area centroid, computed in f64 from the face's own vertices, and `exact.rs` lifted that as
/// truth. Padding the same footprint `1.1` then `6.6` put the second profile an ulp from the
/// first and left faces of area `2.2e-16`.
///
/// Placed **far from the world origin**, because that is where the projected origin is least
/// like the old centroid — if anything about the new rule were fragile with distance, a
/// hundred units of it would show here.
#[test]
fn two_routes_to_one_height_share_a_plane_far_from_the_origin() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    // On a lid, frame coordinates are world x and y — the origin is the world origin projected
    // onto the plane and the axes are `u = +x̂`, `v = +ŷ`. So this ring is world
    // x ∈ [99, 101], y ∈ [99.4, 100.6].
    let rect = || {
        Profile2d::polygon(vec![
            p(99.0, 100.6),
            p(99.0, 99.4),
            p(101.0, 99.4),
            p(101.0, 100.6),
        ])
        .unwrap()
    };
    let top = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        let shell = m.solids.get(s).outer;
        *m.shells
            .get(shell)
            .faces
            .iter()
            .max_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    crate::ops::face_plane(m, f)
                        .ok()
                        .filter(|sp| sp.normal().as_array()[2] > 0.5)
                        .map(|sp| sp.origin.as_array()[2])
                        .unwrap_or(f64::NEG_INFINITY)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a face")
    };
    let build = |dists: &[f64]| -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let mut solid = m.add_cuboid(
            Point3::from_array([98.0, 98.0, 0.0]),
            Point3::from_array([102.0, 102.0, 1.0]),
        );
        m.rebuild_adjacency();
        for &dist in dists {
            let face = top(&m, solid);
            let OpOutput::PadOnFace { solid: out, .. } = apply(
                &mut m,
                &Operation::PadOnFace {
                    face,
                    profile: rect(),
                    dist,
                },
            )
            .expect("pad") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            solid = out;
        }
        (m, solid)
    };
    let (m1, one) = build(&[7.7]);
    let (m2, two) = build(&[1.1, 6.6]);
    // The two boss tops are the same plane, to the bit.
    let z = |m: &Model, s| {
        crate::ops::face_plane(m, top(m, s))
            .unwrap()
            .origin
            .as_array()[2]
    };
    assert_eq!(
        z(&m1, one),
        z(&m2, two),
        "7.7 and 1.1+6.6 disagree on the cap plane"
    );
    // And the two-step route left nothing degenerate behind.
    let shell = m2.solids.get(two).outer;
    let degenerate: Vec<_> = m2
        .shells
        .get(shell)
        .faces
        .iter()
        .filter_map(|&f| {
            let area = nacre_props::face_props(&m2, f).ok()?.area;
            (area < 1e-9).then_some((f, area))
        })
        .collect();
    assert!(
        degenerate.is_empty(),
        "zero-area faces survived: {degenerate:?}"
    );
    assert_eq!(
        m1.shells.get(m1.solids.get(one).outer).faces.len(),
        m2.shells.get(shell).faces.len(),
        "the two routes did not build the same solid"
    );
}

/// ★★★ **One plane, one origin** — even when two faces of it were made by different operations.
///
/// This is what the area centroid could not promise: it was a property of the *face*, so a
/// boolean that reshaped one face moved its sketch origin away from its coplanar neighbour's.
/// The projection is a property of the plane, and surfaces are interned, so the two cannot
/// disagree. Only the origin is asserted — the axes follow the face's `Orientation`, which is
/// today's behaviour and a separate question.
#[test]
fn two_faces_of_one_plane_share_a_sketch_origin() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    // A notch out of one end: the lid `z = 1` becomes two faces of the same plane.
    let cutter = m.add_cuboid(
        Point3::from_array([0.8, -1.0, 0.5]),
        Point3::from_array([1.2, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, a, cutter).expect("cut");
    m.rebuild_adjacency();
    let lids: Vec<_> = m
        .shells
        .get(m.solids.get(r).outer)
        .faces
        .iter()
        .filter(|&&f| {
            crate::ops::face_plane(&m, f)
                .is_ok_and(|sp| sp.normal().as_array()[2] > 0.5 && sp.origin.as_array()[2] == 1.0)
        })
        .copied()
        .collect();
    assert_eq!(lids.len(), 2, "the notch should leave two lid faces");
    let o = |f| crate::ops::face_plane(&m, f).unwrap().origin.as_array();
    assert_eq!(
        o(lids[0]),
        o(lids[1]),
        "coplanar faces disagree on the origin"
    );
    assert_eq!(
        o(lids[0]),
        [0.0, 0.0, 1.0],
        "the lid's origin is the world origin projected"
    );
}
/// ★★★★★ **The convention, face by face.** This table *is* the rule — every property below is a
/// consequence of it, and pinning the consequences without pinning the table would let a
/// different rule that happens to satisfy them slip in.
#[test]
fn the_six_axis_directions_get_the_frames_the_convention_names() {
    let v = |a: [f64; 3]| Vector3::from_array(a);
    for (n, u, w) in [
        // The lid: this is the row that must equal `SketchPlane::world_xy`.
        ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        // Every wall: `v` is +ẑ.
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        ([-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    ] {
        let (gu, gv) = crate::ops::frame_axes(v(n)).expect("a unit normal has frame axes");
        assert_eq!(gu.as_array(), u, "u for normal {n:?}");
        assert_eq!(gv.as_array(), w, "v for normal {n:?}");
        // ★ And the axes stay exactly representable, so the rational construction path still
        // fires — losing that would drop every axis-aligned model to f64 silently.
        let plane = SketchPlane::from_axes(Point3::origin(), gu, gv);
        assert!(plane.exact().is_some(), "exact path lost for normal {n:?}");
    }
}

/// The lid's frame and [`SketchPlane::world_xy`] name the same plane, so they must name it the
/// same way. They did not: `any_perpendicular` gave the lid `u = −ŷ, v = +x̂`, ninety degrees
/// round, and the kernel carried both spellings at once.
#[test]
fn a_lid_gets_the_same_frame_as_the_world_xy_plane() {
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let (u, v) = crate::ops::frame_axes(up).unwrap();
    let w = SketchPlane::world_xy();
    assert_eq!(u.as_array(), w.x_axis.as_array());
    assert_eq!(v.as_array(), w.y_axis.as_array());
}

/// **On anything but a horizontal face, `v` points up.** `u = ẑ × n` is horizontal, so
/// `v·ẑ = 1 − n_z² > 0` whenever `n` is not vertical — the reason a sketch on a wall has "up"
/// where a person expects it. Checked on tilts the axis-aligned table cannot reach.
#[test]
fn every_non_horizontal_face_has_its_v_pointing_up() {
    for n in [
        [1.0, 1.0, 0.0],
        [0.6, 0.0, 0.8],
        [-0.3, 0.5, -0.81],
        [0.0, 1.0, 0.001],
        [7.0, -13.0, 5.0],
    ] {
        let n = Vector3::from_array(n).normalize().unwrap();
        let (u, v) = crate::ops::frame_axes(n).unwrap();
        assert_eq!(u.as_array()[2], 0.0, "u must be horizontal for {n:?}");
        assert!(v.as_array()[2] > 0.0, "v points down for {n:?}");
    }
}

/// `u ⊥ n`, `|u| = 1`, and `(u, v, n)` right-handed — for the tilted normals too, where the
/// table above says nothing.
#[test]
fn the_reference_axis_is_a_unit_normal_perpendicular() {
    for n in [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 1.0, 1.0],
        [-2.0, 0.5, 3.25],
        [1e-9, 0.0, 1.0],
    ] {
        let n = Vector3::from_array(n).normalize().unwrap();
        let (u, v) = crate::ops::frame_axes(n).unwrap();
        assert!((u.norm() - 1.0).abs() < 1e-15, "|u| for {n:?}");
        assert!(u.dot(n).abs() < 1e-15, "u·n for {n:?}");
        assert!((u.cross(v) - n).norm() < 1e-15, "handedness for {n:?}");
    }
}

/// ★★ **The jump at the poles is intended, not a bug to be fixed later.**
///
/// No continuous tangent frame exists on the sphere, so some set of normals must jump; this
/// convention spends that budget on the two poles and nowhere else. A normal a billionth off
/// vertical takes the other branch and lands ninety degrees away — pinned here so the next
/// reader can see it was chosen.
#[test]
fn the_frame_jumps_at_the_poles_and_that_is_the_deal() {
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let tilted = Vector3::from_array([1e-9, 0.0, 1.0]).normalize().unwrap();
    assert_eq!(
        crate::ops::frame_axes(up).unwrap().0.as_array(),
        [1.0, 0.0, 0.0]
    );
    assert_eq!(
        crate::ops::frame_axes(tilted).unwrap().0.as_array(),
        [0.0, 1.0, 0.0]
    );
}

/// No `-0.0` reaches a caller. It compares equal to `0.0` and lifts to the same rational, so
/// this is presentation only — but a frame printed as `[-0.0, 1.0, 0.0]` reads like a defect.
#[test]
fn no_axis_component_is_negative_zero() {
    for n in [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ] {
        let n = Vector3::from_array(n);
        let (u, v) = crate::ops::frame_axes(n).unwrap();
        for c in u.as_array().iter().chain(v.as_array().iter()) {
            assert!(
                !(*c == 0.0 && c.is_sign_negative()),
                "negative zero in {n:?}'s frame"
            );
        }
    }
}

/// The zero vector has no frame — the only `None`.
#[test]
fn the_zero_vector_has_no_frame_axes() {
    assert!(crate::ops::frame_axes(Vector3::zero()).is_none());
}
/// ★★★★★ **A boss on a tilted face, and another beside it — `pad` used to lose the second one.**
///
/// The prism's far cap and the first boss's cap are one plane, and the boolean says so: its
/// classes are decided with evidence, and on a tilted face that evidence is a composed-rotation
/// proof or a coincidence within the limit, never a handle match — the cap's surface has no
/// rational coefficients to intern by, so it is minted fresh.
///
/// `find_face_coplanar_with` then had to guess which surface the class had collapsed to, from
/// handles and an exact `plane_side` on f64 points. Both miss: the survivor carries the *other*
/// operand's surface, and the two f64 planes sit `1.8e-15` apart. `pad` turned "cap not found"
/// into a hard error and **threw away a correct solid** — the volume was already right.
///
/// Now it asks the boolean instead.
#[test]
fn a_second_boss_on_a_tilted_face_keeps_its_cap() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let mut s = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    m.rebuild_adjacency();
    // Two turns, so the face's normal is off every world axis and its frame is not exact.
    for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::rotation(Rotation {
                    axis,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                }),
            },
        )
        .expect("turn") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
    }
    // Where the original +z went, so the face can be picked without a heuristic tie.
    let (sy, cy) = (53f64).to_radians().sin_cos();
    let (sz, cz) = (17f64).to_radians().sin_cos();
    let up = Vector3::from_array([cz * sy, sz * sy, cy]);
    // The **original** tilted face, not a boss raised on it: lowest along `up` among the faces
    // that point that way. Taking the highest would stack the second boss on the first.
    let facing = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        *m.shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .filter(|&&f| crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99))
            .min_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a face along up")
    };
    // The frame is a function of the plane, so it survives the face being reshaped by the first
    // pad — both columns are placed from one reading.
    let (cu, cv) = {
        let f = facing(&m, s);
        let sp = crate::ops::face_plane(&m, f).expect("planar");
        let d = nacre_props::face_props(&m, f).unwrap().centroid - sp.origin;
        (d.dot(sp.x_axis), d.dot(sp.y_axis))
    };
    let before = nacre_props::mass_props(&m, s).unwrap().volume;
    let mut caps = Vec::new();
    for (lo, hi) in [(-1.0f64, -0.4f64), (0.4, 1.0)] {
        let f = facing(&m, s);
        let profile = Profile2d::polygon(vec![
            p(cu + lo, cv - 0.5),
            p(cu + hi, cv - 0.5),
            p(cu + hi, cv + 0.5),
            p(cu + lo, cv + 0.5),
        ])
        .unwrap();
        let OpOutput::PadOnFace { solid, top_face } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: f,
                profile,
                dist: 7.7,
            },
        )
        .expect("a boss on a tilted face keeps its cap") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
        caps.push(top_face);
    }
    // Each column is 0.6 × 1.0 × 7.7 = 4.62 of material.
    let after = nacre_props::mass_props(&m, s).unwrap().volume;
    assert!(
        (after - before - 2.0 * 4.62).abs() < 1e-9,
        "volume {before} -> {after}"
    );
    // ★ And the handle it returned is the boss top, not some other face that happened to pass:
    // area 0.6, outward along the face's normal.
    for cap in caps {
        let props = nacre_props::face_props(&m, cap).expect("a planar cap");
        assert!((props.area - 0.6).abs() < 1e-9, "cap area {}", props.area);
        assert!(
            props.normal.is_some_and(|n| n.dot(up) > 0.99),
            "cap faces the wrong way"
        );
    }
    assert!(nacre_validate::validate(&m).is_empty());
}

// ---- tilted-face sketches: what the kernel does today (2b's corpus) ----
//
// ★★★ These pin **today's** behaviour, not a wish. The suite had almost none of them, so the
// stage that makes tilted frames exact would otherwise be built with nothing to measure against.
// Each says which part 2b is expected to change and which part must not move.

/// A profile edge running `(1, 2)` sweeps a wall whose normal is `(2, 1, 0)` — and **no
/// rational-degree rotation reaches it** (the angle is `atan(1/2)`). That is the case
/// `docs/truth-and-cache.md` names for `Motion::Frame`.
///
/// What holds today and must keep holding:
///
/// * the wall's **world coefficients are rational** — `[2, 1, 0, −10]`, so a frame built on it
///   has exact data to derive from;
/// * its sketch frame follows the arbitrary-axis convention — origin at the world origin's
///   projection `(4, 2, 0) = (10/5)·(2,1,0)`, `u = ẑ × n` normalized, `v = +ẑ`;
/// * a pad on it works — **through `Motion::Frame`**, since stage 2 (an earlier line here
///   said "through the f64 path", which stopped being true then).
///
/// What the final assertion pins is narrower than the test's old name suggests: the
/// *reported* `SketchPlane`'s world axes are irrational, so `exact()` is `None` — true
/// before S4 and after, because it is about the world lift, not about the frame road the
/// operation actually takes. ★ The axes must **not** move — this plane's own frame is the
/// world, so `ẑ × n` is the same vector before and after.
#[test]
fn a_sketch_on_a_prism_side_wall_takes_the_f64_path_today() {
    let (m, wall) = prism_with_a_slanted_wall();
    let sp = crate::ops::face_plane(&m, wall).expect("planar");
    let c = m
        .surface_name
        .get(&m.faces.get(wall).surface)
        .expect("a world-frame wall has a name")
        .narrow()
        .expect("a world-frame wall's name is narrow");
    assert_eq!(c.map(|r| r.to_f64()), [2.0, 1.0, 0.0, -10.0]);
    assert_eq!(
        sp.origin.as_array(),
        [4.0, 2.0, 0.0],
        "the projected origin"
    );
    assert_eq!(
        sp.y_axis.as_array(),
        [0.0, 0.0, 1.0],
        "v points up on a wall"
    );
    assert!(
        (sp.x_axis - Vector3::from_array([-1.0, 2.0, 0.0]).normalize().unwrap()).norm() < 1e-15,
        "u is ẑ × n normalized, got {:?}",
        sp.x_axis.as_array()
    );
    // ★ The gap 2b closes: the frame is not exact, so nothing built here records coefficients.
    assert!(
        sp.exact().is_none(),
        "a tilted frame has no rational form today"
    );
}

/// The same wall, actually used: a boss on it comes out right through the f64 path. Pinned so
/// that making the frame exact cannot change the **answer**, only how it is recorded.
#[test]
fn a_boss_on_a_slanted_wall_is_correct_today() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let profile = centred_on(&m, wall, 0.5);
    let OpOutput::PadOnFace { solid, top_face } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 1.0,
        },
    )
    .expect("pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // Base prism 15 × 3 = 45, plus a 1 × 1 × 1 boss.
    let props = nacre_props::mass_props(&m, solid).unwrap();
    assert!(
        (props.volume - 46.0).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!(
        (nacre_props::face_props(&m, top_face).unwrap().area - 1.0).abs() < 1e-9,
        "the boss top is 1 × 1"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// ★★★ **A sketch on a wall raised from a sketch on a wall** — the nesting `Motion::Frame` was
/// redesigned for. The second wall's *world* normal is irrational, so its frame cannot be
/// written down by naming a normal; only by naming the plane.
///
/// It works today, through f64. Pinned because nesting is where a frame that names its plane
/// by handle must terminate its recursion.
#[test]
fn a_sketch_on_a_wall_raised_from_a_slanted_wall_works_today() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let profile = centred_on(&m, wall, 0.5);
    let OpOutput::PadOnFace { solid, top_face } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 1.0,
        },
    )
    .expect("pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // A side face of that boss: not the cap, not on the original wall's plane.
    let cap_n = nacre_props::face_props(&m, top_face)
        .unwrap()
        .normal
        .unwrap();
    let side = *m
        .shells
        .get(m.solids.get(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            f != top_face
                && nacre_props::face_props(&m, f).is_ok_and(|p| {
                    (p.area - 1.0).abs() < 1e-9
                        && p.normal.is_some_and(|n| n.dot(cap_n).abs() < 0.5)
                })
        })
        .expect("a boss side face");
    let profile = centred_on(&m, side, 0.2);
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: side,
            profile,
            dist: 0.5,
        },
    )
    .expect("nested pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let props = nacre_props::mass_props(&m, solid).unwrap();
    // The nested boss is 0.4 × 0.4 × 0.5 = 0.08.
    assert!(
        (props.volume - 46.08).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A pocket on the slanted wall — the other face-based operation, so the sweep runs inward.
#[test]
fn a_pocket_in_a_slanted_wall_is_correct_today() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let profile = centred_on(&m, wall, 0.5);
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: wall,
            profile,
            dist: 0.5,
        },
    )
    .expect("pocket") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let props = nacre_props::mass_props(&m, solid).unwrap();
    assert!(
        (props.volume - 44.5).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A pentagonal prism whose fourth wall is slanted, and that wall's handle. Footprint area 15
/// (a 4×4 square less the 1×2 triangle the slant cuts off), swept 3.
fn prism_with_a_slanted_wall() -> (Model, Handle<Face>) {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let __w1 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w1,
            profile: Profile2d::polygon(vec![
                p(0.0, 0.0),
                p(4.0, 0.0),
                p(4.0, 2.0),
                p(3.0, 4.0),
                p(0.0, 4.0),
            ])
            .unwrap(),
            dist: 3.0,
        },
    )
    .expect("extrude") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let wall = *m
        .shells
        .get(m.solids.get(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            crate::ops::face_plane(&m, f).is_ok_and(|sp| {
                let n = sp.normal().as_array();
                n[0].abs() > 0.1 && n[1].abs() > 0.1 && n[2].abs() < 1e-12
            })
        })
        .expect("a slanted wall");
    (m, wall)
}

/// A `2·half` square centred on `face`, in that face's own sketch frame. The frame is a
/// function of the plane, so reading it once is enough even if the face is reshaped later.
/// ★★★★★ **The payoff, and its limit — both measured.**
///
/// Two bosses of the same height on one tilted face used to be two plane records that agreed
/// only if their f64 coefficients happened to. Written in the plane's own frame they are both
/// `w = 7.7`, and `SurfaceKey` is `(coefficients, motion)` — so they are **one
/// `Handle<Surface>` at construction**, before anything is compared. And every face of the
/// result states itself exactly, where before a tilted sketch recorded nothing at all.
///
/// ★★★ **What this does *not* buy, stated plainly**: `7.7` against `1.1 + 6.6` — the target
/// the plan named. Stacking sketches the second boss on the **first boss's cap**, which is a
/// different plane and therefore a different frame, so the two caps come out `w = 7.7` and
/// `w = 6.6` — two exact descriptions of one plane that `SurfaceKey` cannot equate. That is
/// not a regression (before this they had no descriptions at all, and the merge still happens
/// through the judge), but the plan's headline claim only holds for sketches sharing a frame,
/// which is what this pins instead.
#[test]
fn two_bosses_on_one_tilted_face_share_a_cap_plane_by_name() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let __w0 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w0,
            profile: Profile2d::polygon(vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)])
                .unwrap(),
            dist: 3.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let mut s = solid;
    for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::rotation(Rotation {
                    axis,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                }),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
    }
    let (sy, cy) = (53f64).to_radians().sin_cos();
    let (sz, cz) = (17f64).to_radians().sin_cos();
    let up = Vector3::from_array([cz * sy, sz * sy, cy]);
    let facing = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        *m.shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .filter(|&&f| crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99))
            .min_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .unwrap()
    };
    let mut caps = Vec::new();
    for (lo, hi) in [(-1.0f64, -0.4f64), (0.4, 1.0)] {
        let f = facing(&m, s);
        let sp = crate::ops::face_plane(&m, f).unwrap();
        let d = nacre_props::face_props(&m, f).unwrap().centroid - sp.origin;
        let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
        let profile = Profile2d::polygon(vec![
            p(cu + lo, cv - 0.5),
            p(cu + hi, cv - 0.5),
            p(cu + hi, cv + 0.5),
            p(cu + lo, cv + 0.5),
        ])
        .unwrap();
        let OpOutput::PadOnFace { solid, top_face } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: f,
                profile,
                dist: 7.7,
            },
        )
        .expect("boss") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
        let su = m.faces.get(top_face).surface;
        assert_eq!(
            m.surface_name
                .get(&su)
                .and_then(|n| n.narrow())
                .map(|c| c.map(|r| r.to_f64())),
            Some([0.0, 0.0, 10.0, -77.0]),
            "★ a cap raised in a frame records `w = 7.7` there, exactly"
        );
        caps.push(su);
    }
    assert_eq!(caps[0], caps[1], "one plane, one handle, by name");
    // Rule audit: how many of the result's face surfaces state themselves exactly.
    let sh = m.solids.get(s).outer;
    let (mut with, mut without) = (0, 0);
    for &f in &m.shells.get(sh).faces {
        let su = m.faces.get(f).surface;
        if m.surface_name.contains_key(&su) {
            with += 1
        } else {
            without += 1
        }
    }
    assert_eq!(
        (with, without),
        (16, 0),
        "★ every face of a twice-turned, twice-bossed result states itself exactly"
    );
}

/// ★★★★★ **A plane states itself exactly, and its f64 axes are the realization of that.**
///
/// This is the whole point of closing the struct: `(1, 1, 1)` is coefficients `[1, 1, 1, 0]`,
/// three integers, while the unit axes derived from it square to `0.9999999999999999…`. The
/// old API stored only the axes and threw the normal away, so nothing exact survived the door.
#[test]
fn a_named_plane_records_what_its_caller_stated() {
    // The canonical name is *derived* from the definition's points now; reading it back is
    // how the old coefficient assertions keep their meaning.
    let f = |d: &PlaneDef| {
        let p = d.points();
        nacre_scalar::plane_name_exact(p[0], p[1], p[2])
            .expect("a definition names a plane")
            .narrow()
            .expect("these fixtures are narrow")
            .map(|r| r.to_f64())
    };
    // The three world planes, with the axes the script layer documents.
    for (p, want, u) in [
        (
            SketchPlane::world_xy(),
            [0.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ),
        (
            SketchPlane::world_yz(),
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ),
        // ★ `ẑ × n` would give `−x̂` here; a named plane says `+u = ẑ` instead.
        (
            SketchPlane::world_zx(),
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ),
    ] {
        let d = p.def.expect("a world plane states itself");
        assert_eq!(f(&d), want);
        assert_eq!(d.ref_dir().map(|r| r.to_f64()), u);
        assert_eq!(p.x_axis().as_array(), u, "the axis follows the definition");
    }
    // A tilted normal the caller wrote: exact coefficients, though its axes never can be.
    let tilt =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    assert_eq!(f(&tilt.def.unwrap()), [1.0, 1.0, 1.0, 0.0]);
    assert!(
        tilt.exact().is_none(),
        "★ the axes still have no exact form — that is what the definition exists to replace"
    );
    // Three written points: the plane is exact and `+u` runs toward `x_point`.
    let tp = SketchPlane::through_points(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([1.0, 2.0, 0.0]),
        Point3::from_array([1.0, 0.0, 3.0]),
    )
    .unwrap();
    let d = tp.def.unwrap();
    assert_eq!(f(&d), [1.0, 0.0, 0.0, -1.0], "the plane x = 1");
    assert_eq!(
        d.ref_dir().map(|r| r.to_f64()),
        [0.0, 2.0, 0.0],
        "x_point − origin"
    );
    assert_eq!(d.origin().map(|r| r.to_f64()), [1.0, 0.0, 0.0]);
    // Moving the sketch origin keeps the plane and moves only `(0, 0)`.
    let moved = tp.with_origin(Point3::from_array([1.0, 5.0, 5.0]));
    let m = moved.def.unwrap();
    assert_eq!(f(&m), f(&d), "same plane");
    assert_eq!(m.origin().map(|r| r.to_f64()), [1.0, 5.0, 5.0]);
    // ★★★★★ **The invariant that used to tie the two halves together — the origin is *on*
    // the plane — is structural now: the origin IS `points[0]`, so there are no halves to
    // disagree (the failure this guards against cost a boolean 0.04 of volume once). The
    // loop keeps the check as a derivation audit: substituting the origin into the *derived*
    // name must still give zero, or the derivation itself is wrong.
    for p in [
        SketchPlane::world_xy(),
        SketchPlane::world_yz(),
        SketchPlane::world_zx(),
        SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
        SketchPlane::world_zx().with_origin(Point3::from_array([10.0, 20.0, 5.0])),
        SketchPlane::from_origin_normal(
            Point3::from_array([0.0, 1.3, 0.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
        )
        .unwrap(),
        tilt,
        tp,
        moved,
    ] {
        let d = p.def.expect("stated");
        let pts = d.points();
        let c = nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2])
            .expect("a definition names a plane")
            .narrow()
            .copied()
            .expect("these fixtures are narrow");
        let mut s = c[3];
        for (ck, ok) in c.iter().zip(d.origin()) {
            s = s.checked_add(ck.checked_mul(ok).unwrap()).unwrap();
        }
        assert_eq!(
            s,
            nacre_scalar::Rat::from_int(0),
            "the origin must lie on the plane its points name: {:?} vs {:?}",
            c.map(|r| r.to_f64()),
            d.origin().map(|r| r.to_f64())
        );
    }
    // ★ And the axes-only route records the axes' decimal truth (S6a — this used to assert
    // `def.is_none()`, "honest about having no definition"; the honest statement now is the
    // definition itself, `[o, o + x, o + y]`).
    let axes = SketchPlane::from_axes(
        Point3::from_array([1.0, 2.0, 3.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
    );
    let r = |v: [f64; 3]| v.map(|x| nacre_scalar::Rat::from_decimal(x).unwrap());
    assert_eq!(
        axes.def.expect("axes state their truth").points(),
        [r([1.0, 2.0, 3.0]), r([2.0, 2.0, 3.0]), r([1.0, 3.0, 3.0])],
        "o, o + x, o + y"
    );
    // A degenerate pair still names nothing.
    assert!(
        SketchPlane::from_axes(
            Point3::origin(),
            Vector3::from_array([1.0, 0.0, 0.0]),
            Vector3::from_array([2.0, 0.0, 0.0]),
        )
        .def
        .is_none(),
        "parallel axes name no plane"
    );
}

/// ★★★★★ **A prism raised on a named tilted plane states every one of its faces.**
///
/// In world coordinates that plane's axes are irrational, so the whole prism used to drop to
/// f64 and record nothing. Two things fixed it: the base cap **is** the plane the caller
/// named, so it states itself in the world; and the walls and far cap are built **inside that
/// plane's frame**, where the axes are `x̂`/`ŷ` and the profile's own decimals are the truth.
///
/// ★★★ **The base cap stays in the world on purpose.** Writing it as `[0,0,1,0]` in this
/// prism's frame would be a second exact description of one plane under a different
/// `SurfaceKey` — the duplication this work exists to remove. Stated in the world it is
/// `Constructed`, its judgment stays exact, and two extrudes share it whatever frames they chose.
#[test]
fn a_prism_on_a_named_tilted_plane_states_all_of_its_faces() {
    let plane =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    assert!(plane.exact().is_none(), "the axes have no exact form");
    let mut m = Model::new();
    let __f102 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f102,
            profile: square(),
            dist: 1.0,
        },
    )
    .expect("extrude on a tilted plane") else {
        unreachable!()
    };
    let coeffs = |f: Handle<Face>, m: &Model| {
        m.surface_name
            .get(&m.faces.get(f).surface)
            .and_then(|n| n.narrow())
            .map(|c| c.map(|r| r.to_f64()))
    };
    assert_eq!(
        coeffs(faces[0], &m),
        Some([1.0, 1.0, 1.0, 0.0]),
        "★ the base cap is the caller's plane, in the world"
    );
    assert!(
        matches!(
            m.surface_truth(m.faces.get(faces[0]).surface),
            nacre_topo::SurfaceTruth::Plane { motion: None, .. }
        ),
        "★ and it carries no motion, so its judgment stays exact"
    );
    assert_eq!(
        coeffs(faces[1], &m),
        Some([0.0, 0.0, 1.0, -1.0]),
        "★ the far cap is `w = dist` in the frame"
    );
    // ★ Every face now states itself — that is the whole measurement.
    for &f in &faces {
        assert!(coeffs(f, &m).is_some(), "a face with no exact plane");
    }

    // ★★★★★ **Two extrudes on one named plane put their far caps on one handle** — by name,
    // at construction, with no f64 comparison. That is what the frame buys over the f64 path,
    // where the two would agree only if their rounded coefficients happened to.
    let cap_of = |m: &mut Model, d: f64| -> Handle<nacre_geom::Surface> {
        let __g201 = datum_frame(m, plane);
        let OpOutput::Extrude { faces, .. } = apply(
            m,
            &Operation::Extrude {
                frame: __g201,
                profile: square(),
                dist: d,
            },
        )
        .expect("extrude") else {
            unreachable!()
        };
        m.faces.get(faces[1]).surface
    };
    let a = cap_of(&mut m, 2.5);
    let b = cap_of(&mut m, 2.5);
    assert_eq!(a, b, "one height on one plane is one plane");
    assert_ne!(a, cap_of(&mut m, 2.6), "and a different height is not");
}

/// ★★★★★ **S4's end-to-end lock: the f64-fallback chain is cut on the `n·n`-overflow
/// population** — the measured 1.6%, and the one today's producers actually reach.
///
/// A prism raised on a fully tilted, exactly-orthonormal decimal frame from a 16-digit
/// profile has walls whose names run to ~110 bits: **narrow names whose squared lengths
/// (~2^220) overflow `i128`**, so `plane_frame_default`/`plane_frame_named` hard-declined
/// them and a pad on such a wall fell to `Swept::along` — every new face point-less, the
/// chain that kept the `Inexact` state alive. Now the frame realizes through the wide
/// road and **every face of the result states its exact points**.
///
/// ★ Probed while building this fixture: a sketch→extrude wall's canonical name caps out
/// around ~115 bits (the profile's decimal window bounds the products), so **`Wide`-named
/// faces do not arise from today's construction route at all** — the `Wide` opening is
/// locked at the basis level (`a_wide_plane_hosts_a_canonical_frame`) and the end-to-end
/// chain here on the population that exists. The fixture qualifies itself rather than
/// assuming its population (the S2 census lesson).
#[test]
fn a_pad_on_a_wall_with_overflowing_squares_takes_the_exact_road() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    // A fully tilted, exactly-orthonormal decimal frame: u·u = v·v = 1 and u·v = 0 hold in
    // the lifted rationals, and the normal u×v = (0.64, −0.48, 0.6) is not axis-aligned —
    // so the walls' names mix all three coordinates.
    let plane = crate::ops::SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let __f101 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f101,
            profile: Profile2d::polygon(vec![
                p(0.1111111111111111, 0.1234567890123456),
                p(4.123456789012345, 0.2345678901234567),
                p(3.9876543210987654, 3.1234567890123459),
                p(0.2222222222222222, 2.765432109876543),
            ])
            .unwrap(),
            dist: 2.5,
        },
    )
    .expect("extrude on a Pythagorean frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // Fixture qualification: a wall whose name is narrow but whose squared lengths are not
    // — the exact population the narrow frame derivation hard-declines.
    let shell = m.solids.get(solid).outer;
    let wall = *m
        .shells
        .get(shell)
        .faces
        .iter()
        .find(|&&f| {
            let s = m.faces.get(f).surface;
            m.surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
        })
        .expect("an nn-overflow wall — retune the fixture constants if this fails");
    let profile = centred_on(&m, wall, 0.3);
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 0.4,
        },
    )
    .expect("S4: a pad on a wide-named wall must build") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // The chain is cut: every face of the result states its exact points.
    let mut missing = 0;
    let mut total = 0;
    for &f in &m.shells.get(m.solids.get(solid).outer).faces {
        total += 1;
        if !matches!(
            m.surface_truth(m.faces.get(f).surface),
            nacre_topo::SurfaceTruth::Plane { .. }
        ) {
            missing += 1;
        }
    }
    assert_eq!(
        missing, 0,
        "{missing} of {total} faces carry no exact points — the f64 fallback fired"
    );
}

/// ★★★★★ S6a terminal lock: **an axes-only tilted frame extrudes on the exact road.**
///
/// The sibling of `a_pad_on_a_wall_with_overflowing_squares_takes_the_exact_road`, for the
/// population that S4 could not reach: a `from_axes` plane whose axes never lift to exact
/// orthonormal rationals (a rotated frame — dev-log's "#28: the caller who only has axes",
/// the first refutation of `Inexact`'s removal). Before S6a its `def` was `None` by design
/// and the whole prism fell silently to f64 — every surface point-less, undemotable under
/// any later motion. Now the axes' decimal truth is the definition, the canonical name is
/// `Wide` (asserted — the full-width crosses exceed `i128`), `WideFrame::named_of` realizes
/// the frame, and **every face of the prism records its exact points**.
#[test]
fn a_prism_on_an_axes_only_tilted_frame_takes_the_exact_road() {
    let mut m = Model::new();
    let plane = crate::ops::SketchPlane::from_axes(
        Point3::from_array([0.2547863291057384, -0.5123456789012345, 1.5432109876543211]),
        Vector3::from_array([0.7123456789012345, 0.5876543210987654, 0.4098765432101234]),
        Vector3::from_array([-0.5876543210987654, 0.7123456789012345, 0.1234567890123456]),
    );
    // Fixture qualification: no world lift (the frame road is the only exact road), the
    // definition exists, and its name is genuinely Wide.
    assert!(
        plane.exact().is_none(),
        "the axes must not lift orthonormal"
    );
    let d = plane.def.expect("S6a: axes state their decimal truth");
    let pts = d.points();
    assert!(
        nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2])
            .expect("a plane")
            .narrow()
            .is_none(),
        "the fixture was chosen to have a Wide name — retune the axes if this fails"
    );
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let __f100 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f100,
            profile: Profile2d::polygon(vec![
                p(0.1234567890123456, 0.2345678901234567),
                p(2.765432109876543, 0.3456789012345678),
                p(2.543210987654321, 1.9876543210987654),
                p(0.3456789012345678, 1.8765432109876543),
            ])
            .unwrap(),
            dist: 1.3,
        },
    )
    .expect("S6a: an axes-only tilted extrude must build exactly") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // The end-to-end claim: nothing fell to f64 — every face states its exact points.
    let mut missing = 0;
    let mut total = 0;
    for &f in &m.shells.get(m.solids.get(solid).outer).faces {
        total += 1;
        if !matches!(
            m.surface_truth(m.faces.get(f).surface),
            nacre_topo::SurfaceTruth::Plane { .. }
        ) {
            missing += 1;
        }
    }
    assert_eq!(
        missing, 0,
        "{missing} of {total} faces carry no exact points — the f64 fallback fired"
    );
    // Geometry check: the realized frame is orthonormal, so the prism's volume is the
    // profile's own area times the sweep — independent of the tilt.
    let props = nacre_props::mass_props(&m, solid).unwrap();
    let ring = [
        [0.1234567890123456, 0.2345678901234567],
        [2.765432109876543, 0.3456789012345678],
        [2.543210987654321, 1.9876543210987654],
        [0.3456789012345678, 1.8765432109876543],
    ];
    let mut area2 = 0.0f64;
    for i in 0..4 {
        let (a, b) = (ring[i], ring[(i + 1) % 4]);
        area2 += a[0] * b[1] - b[0] * a[1];
    }
    let want = (area2 / 2.0).abs() * 1.3;
    assert!(
        (props.volume - want).abs() < 1e-9,
        "volume {} vs analytic {want}",
        props.volume
    );
}

/// ★★★ S6b: **what the f64 fallback used to build silently is a named reject now.** A plane
/// with no exact statement — axes outside the decimal window — and a sweep distance the
/// window cannot hold each get their own name at the operation's door. The prisms these
/// used to build recorded no exact points and could not survive a motion; the reject is the
/// honest form of the same fact.
#[test]
fn a_prism_the_exact_arithmetic_cannot_state_is_refused_by_name() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let square =
        || Profile2d::polygon(vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]).unwrap();
    // ① Axes outside the decimal window: no def, no world lift — no exact form at all.
    let far = crate::ops::SketchPlane::from_axes(
        Point3::origin(),
        Vector3::from_array([1e300, 0.0, 0.0]),
        Vector3::from_array([0.0, 1e300, 0.0]),
    );
    assert!(
        far.exact().is_none() && far.def.is_none(),
        "fixture: no exact statement"
    );
    // A plane with no exact statement cannot even be stated as a datum, which is where the
    // rejection now lands — one step earlier than it used to, and by the same name.
    assert_eq!(
        apply(
            &mut Model::new(),
            &Operation::DatumPlane {
                def: DatumDef::Stated(far)
            },
        ),
        Err(OpError::PlaneWithoutExactForm)
    );
    // ② A sweep distance outside the window — the profile-coordinate rule's sibling.
    assert_eq!(
        apply(
            &mut Model::new(),
            &Operation::Extrude {
                frame: SketchFrame::world(&Model::new(), Axis::Z),
                profile: square(),
                dist: 1e300,
            },
        ),
        Err(OpError::DistOutsideDecimalWindow)
    );
}

fn centred_on(m: &Model, face: Handle<Face>, half: f64) -> Profile2d {
    let sp = crate::ops::face_plane(m, face).expect("planar");
    let d = nacre_props::face_props(m, face).unwrap().centroid - sp.origin;
    let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    Profile2d::polygon(vec![
        p(cu - half, cv - half),
        p(cu + half, cv - half),
        p(cu + half, cv + half),
        p(cu - half, cv + half),
    ])
    .unwrap()
}

/// ★★ S9, the node-omission normalization's ground: **a seeded world plane's canonical frame
/// realizes bit-exactly on the unit axes** — origin `(0,0,0)`, every axis component `0.0` or
/// `±1.0`, no rounding anywhere. Such a frame's axes lift to exact orthonormal rationals, so
/// the operation's `exact()` gate elides the node (rule 207's normalization) and loses
/// nothing: the frame the node would state is the world statement already made.
///
/// The expected axes are the **arbitrary-axis convention's**, not `axis_plane`'s script
/// triples: for the ZX plane the convention's `+u` is `ẑ × ŷ = −x̂` where the script's is
/// `+ẑ` — the documented case a caller states through a `Named` placement instead. The seed
/// *points* carry the script triple; the *canonical frame* is the plane's own.
#[test]
fn a_seeded_planes_canonical_frame_is_the_world_basis_exactly() {
    use nacre_scalar::Axis;
    let m = Model::new();
    // (axis, expected û, v̂, ŵ) — ŵ is the +axis (the canonical name's own sign).
    let want = [
        (Axis::Z, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        (Axis::X, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        (Axis::Y, [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
    ];
    for (axis, u, v, w) in want {
        let h = m.world_plane(axis);
        let (o, ru, rv, rw) = crate::rotated_vertex::frame_world_basis(
            &m,
            h,
            &nacre_topo::FramePlacement::Canonical,
            false,
        )
        .expect("a seed always carries a name");
        assert_eq!(
            (o, ru, rv, rw),
            ([0.0; 3], u, v, w),
            "{axis:?}: not the world basis"
        );
    }
}

// ---------------------------------------------------------------------------
// The rational road (M6-2a C4a): `combinatorics::point_in_faces_rat`.
//
// Its production consumer is the cylinder band (C4b), but the road answers a question about
// **planar** solids, so it is locked here against one — a non-convex L-prism, where a convex
// test would be wrong and the notch is the counterexample the plan's uniform-slab theorem was
// rewritten around.
// ---------------------------------------------------------------------------

/// A solid's faces as the rational road takes them: `(plane class, rings of node triples)`.
fn component_triples(
    m: &Model,
    s: Handle<Solid>,
    setup: &PlaneSetup,
    jd: &Judge<'_, WorkingPlane>,
) -> combinatorics::ComponentTriples {
    let mut out = combinatorics::ComponentTriples::new();
    for sh in solid_shell_handles(m, s) {
        for &fh in &m.shells.get(sh).faces {
            let fp = setup.surf_ix[&fh];
            let cls = setup.plane_ix[fp].plane();
            let outer =
                combinatorics::face_vertex_triples(m, fh, fp, &setup.inc_a, jd, &setup.plane_ix)
                    .expect("a planar operand names its outer loop");
            let mut rings = vec![outer.poly().expect("a poly outer").triples.clone()];
            for h in combinatorics::hole_rings(m, fh, fp, &setup.inc_a, jd, &setup.plane_ix)
                .expect("a planar operand names its holes")
            {
                rings.push(h.poly().expect("a poly hole").triples.clone());
            }
            out.push((cls, rings));
        }
    }
    out
}

/// The L-prism's own profile as the differential's **oracle**: a 2D crossing parity on the ring
/// the fixture was written from, plus the extrusion's z range. It reads the test's inputs, never
/// the model the road reads — the two must not share a derivation.
fn l_prism_oracle(p: [f64; 3]) -> bool {
    let ring = [
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ];
    let inside_2d = matches!(
        nacre_geom::intersect::point_in_ring_2d(p2(p[0], p[1]), &ring),
        nacre_geom::intersect::RingSide::Inside
    );
    inside_2d && p[2] > 0.0 && p[2] < 1.0
}

fn rat3(p: [f64; 3]) -> [Rat; 3] {
    p.map(|x| Rat::from_decimal(x).expect("a fixture coordinate is a decimal"))
}

/// A grid of witnesses over the L-prism and its notch, road against oracle — **and both ray
/// directions agree**, which is the road's own claim: parity does not depend on where you look
/// from. The notch is why: it is outside the material and inside the bounding box, so a convex
/// or "z-range" answer would be wrong on a quarter of the grid.
#[test]
fn the_rational_road_agrees_with_the_profile_over_the_l_notch() {
    let (m, s) = l_prism();
    let setup = plane_index_setup(&m, s, s).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let faces = component_triples(&m, s, &setup, &jd);

    let coords = [0.25, 0.75, 1.25, 1.75, 2.25];
    let mut inside = 0;
    let mut answers = 0;
    let mut silent: Vec<[f64; 3]> = Vec::new();
    let mut abstained: Vec<([f64; 3], [i128; 3])> = Vec::new();
    for x in coords {
        for y in coords {
            let p = [x, y, 0.5];
            let want = l_prism_oracle(p);
            let mut answered = 0;
            for dir in [[0, 0, 1], [1, 0, 0], [0, 1, 0], [1, 1, 1]] {
                let d = dir.map(Rat::from_int);
                let got = combinatorics::point_in_faces_rat(&jd, &rat3(p), &d, &faces)
                    .unwrap_or_else(|e| panic!("{p:?} along {dir:?}: {e:?}"));
                // An abstention is a legal answer (a grazing ray), so it is not asserted
                // against here — it is collected and named below instead.
                let Some(got) = got else {
                    abstained.push((p, dir));
                    continue;
                };
                assert_eq!(got, want, "point {p:?} along {dir:?}");
                answered += 1;
            }
            if answered == 0 {
                silent.push(p);
            }
            answers += answered;
            inside += usize::from(want);
        }
    }
    // The instrument moved: the grid is not all-outside or all-inside — the notch column is
    // genuinely void while the bar beside it is material.
    assert_eq!(inside, 12, "12 of 25 grid points are material");
    // ★ **Every abstention is named, not counted.** The three are the diagonal ray from a
    // point 0.25 short of a vertical edge, which meets that edge exactly — at the reflex corner
    // (1,1) and at the convex corners (1,2) and (2,1). That the grazing detector fires there,
    // and only there, is this grid's negative control: a road that never abstained and a road
    // that abstained everywhere would both fail this line.
    assert_eq!(
        abstained,
        vec![
            ([0.75, 0.75, 0.5], [1, 1, 1]),
            ([0.75, 1.75, 0.5], [1, 1, 1]),
            ([1.75, 0.75, 0.5], [1, 1, 1]),
        ],
        "the grazing rays are exactly the diagonals through the vertical edges"
    );
    // And no point was left undecided: a grazed direction is repaired by the others, which is
    // the abstention protocol's whole claim.
    assert!(silent.is_empty(), "points no direction decided: {silent:?}");
    assert_eq!(answers, 25 * 4 - 3, "every other query decided");
}

/// The three non-generic rays, each abstaining rather than guessing: a ray that **starts on** the
/// solid's surface, a ray that **runs inside** a face's plane, and a ray that **grazes an edge**.
/// Abstention is a fact about the ray — the same point answers along another direction, which the
/// test shows in the same breath.
#[test]
fn a_non_generic_ray_abstains_and_another_direction_answers() {
    let (m, s) = l_prism();
    let setup = plane_index_setup(&m, s, s).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let faces = component_triples(&m, s, &setup, &jd);
    let ask = |p: [f64; 3], dir: [i128; 3]| {
        combinatorics::point_in_faces_rat(&jd, &rat3(p), &dir.map(Rat::from_int), &faces).unwrap()
    };

    // ① The origin is on the bottom cap, inside its material: "inside" has no answer there.
    assert_eq!(ask([0.5, 0.5, 0.0], [0, 0, 1]), None);
    // ② The same origin, a direction lying in that cap's plane: no crossing parity at all.
    assert_eq!(ask([0.5, 0.5, 0.0], [1, 0, 0]), None);
    // ③ A grazing ray: from (0.5,0.5,0.5) toward (2,1,0.5), the convex corner where the x=2 and
    //    y=1 walls meet — the hit lands on both faces' boundary.
    assert_eq!(ask([0.5, 0.5, 0.5], [3, 1, 0]), None);
    // …and the same point, off the corner, answers.
    assert_eq!(ask([0.5, 0.5, 0.5], [1, 0, 0]), Some(true));
}

/// A class with no exact description is **not** an abstention: no direction repairs it, so the
/// road raises `WitnessNotRational` rather than letting a substrate limit wear "no clear ray"'s
/// clothes. Measured on a rotated L-prism, whose classes carry realized coefficients only.
#[test]
fn a_rotated_class_refuses_the_rational_road_by_name() {
    use nacre_scalar::{Angle, Isometry, Rotation};
    let (mut m, s) = l_prism();
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    let t = transform(&mut m, s, &iso).unwrap();
    let setup = plane_index_setup(&m, t, t).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    assert!(
        setup.geom.iter().any(|p| p.rotated),
        "the fixture must actually be rotated"
    );
    let faces = component_triples(&m, t, &setup, &jd);
    let err = combinatorics::point_in_faces_rat(
        &jd,
        &rat3([0.5, 0.5, 0.5]),
        &[0, 0, 1].map(Rat::from_int),
        &faces,
    )
    .expect_err("a rotated class has no exact description to divide by");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::WitnessNotRational,
                ..
            }
        ),
        "{err:?}"
    );
}
