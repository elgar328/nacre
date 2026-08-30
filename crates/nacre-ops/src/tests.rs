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
use crate::combinatorics::NodeId;
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

/// ★ **A fixture with no cylinders, said as a fact rather than left as a hole.** The table is how
/// a `NodeId::Branch` reaches its definition, so an all-plane fixture has nothing to put in it —
/// and a bare `&[]` at a call site reads like something forgotten.
const NO_CYLS: &[crate::planes::WorkingCyl] = &[];

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
pub(crate) fn extrude_log_op(profile: Profile2d, dist: f64) -> Operation {
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
            NO_CYLS,
            p,
            &combinatorics::ring_from_names(p, &outer).unwrap()
        )
        .unwrap(),
        1
    );
    assert_eq!(
        combinatorics::loop_winding(
            &crate::planes::test_judge(&planes),
            NO_CYLS,
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
                NO_CYLS,
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
            NO_CYLS,
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
                NO_CYLS,
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
            NO_CYLS,
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
            NO_CYLS,
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
            NO_CYLS,
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

/// A ring of vertex names as the plane triples these fixtures assert about.
fn names(ring: &[combinatorics::NodeId]) -> Ring {
    ring.iter()
        .map(|&n| combinatorics::three_plane_name(n).expect("a three-plane node"))
        .collect()
}

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
            &[],
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
                &[],
            )
            .unwrap();
            return (
                planes,
                plane_ix[fp].plane(),
                names(&outer.poly().expect("a poly outer").triples),
                names(&hole.poly().expect("a poly hole").triples),
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
            world_rat: None,
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
            let triple = combinatorics::three_plane_name(c.triple).expect("a three-plane node");
            if !triple.contains(&c.wc) {
                trace_rule_checked += 1;
                let mut derived = triple.to_vec();
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

/// Does any surface of this solid record a motion? **The premise every fixture below asserts
/// first**: a translation whose `f64` landing happens to be exact records nothing, and a fixture
/// that silently took that road would prove the opposite of what it claims (measured — a boss
/// fixture did exactly that while the road under test was still shut).
fn carries_motion(m: &Model, s: Handle<Solid>) -> bool {
    crate::planes::solid_shell_handles(m, s)
        .into_iter()
        .flat_map(|sh| m.shells.get(sh).faces.clone())
        .any(|fh| m.plane_motion(m.faces.get(fh).surface).is_some())
}

/// A cylinder tool translated onto a plate, then cutting it — **the hole-pattern idiom**, and
/// the smallest shape the moved-cylinder road serves. The offset is non-dyadic on purpose, so
/// the move records a chain (asserted) and the boolean has to carry the statement out to the
/// world itself.
#[test]
fn a_translated_tool_cuts() {
    use nacre_scalar::Rat;
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = m.add_cylinder(
        Point3::from_array([7.3, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.1,
        30.0,
    );
    m.rebuild_adjacency();
    let tool = transform(
        &mut m,
        tool,
        &nacre_scalar::Isometry::translation([
            Rat::try_from_f64(10.7).unwrap(),
            Rat::from_int(0),
            Rat::from_int(0),
        ]),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(carries_motion(&m, tool), "the move records a chain");
    let out = boolean(&mut m, BoolKind::Cut, plate, tool).expect("the moved tool cuts");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    let want = 40.0 * 40.0 * 10.0 - std::f64::consts::PI * 2.1 * 2.1 * 10.0;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
}

/// The same total move **split in two** — the chain is walked and folded, not read one node
/// deep. Same solid, so the same volume: the two roads must agree to the tolerance the oracle
/// is stated at.
#[test]
fn a_chained_translation_folds() {
    use nacre_scalar::Rat;
    let shift = |m: &mut Model, s, x: f64| {
        let s = transform(
            m,
            s,
            &nacre_scalar::Isometry::translation([
                Rat::try_from_f64(x).unwrap(),
                Rat::from_int(0),
                Rat::from_int(0),
            ]),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    };
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = m.add_cylinder(
        Point3::from_array([7.3, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.1,
        30.0,
    );
    m.rebuild_adjacency();
    let tool = shift(&mut m, tool, 5.2);
    let tool = shift(&mut m, tool, 5.5);
    assert!(carries_motion(&m, tool), "the moves record a chain");
    let out = boolean(&mut m, BoolKind::Cut, plate, tool).expect("the chained tool cuts");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    let want = 40.0 * 40.0 * 10.0 - std::f64::consts::PI * 2.1 * 2.1 * 10.0;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
}

/// A **bored body** translated and fused onto its twin: the plane side of the same fact — the
/// walls perpendicular to the move carry the chain, and the cylinder roads need their world
/// coefficients. The offset moves all three axes, so no wall of either body is shared (a shared
/// wall is the contact family, another cell).
#[test]
fn a_translated_bored_body_fuses() {
    use nacre_scalar::Rat;
    let bored = |m: &mut Model| {
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 10.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([7.3, 7.3, -5.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.1,
            30.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    };
    let mut m = Model::new();
    let a = bored(&mut m);
    let b = bored(&mut m);
    let b = transform(
        &mut m,
        b,
        &nacre_scalar::Isometry::translation(
            [25.7, 3.3, 2.0].map(|c| Rat::try_from_f64(c).unwrap()),
        ),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(carries_motion(&m, b), "the move records a chain");
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved body fuses");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // Two bored plates, their overlap counted once, and the part of b's bore that a fills back.
    let pi = std::f64::consts::PI;
    let bore = pi * 2.1 * 2.1;
    let want = 2.0 * (16000.0 - bore * 10.0) - (14.3 * 36.7 * 8.0) + bore * 8.0;
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
}

/// ★★ **Two bored plates meeting face to face fuse.** The shared wall is interior, so the result
/// keeps no face on it — and the corners that sat on it have to *dissolve*, which they can only do
/// once each cell's bottom (and top) become **one** face. That merge used to be skipped outright
/// whenever a member carried a circle hole, so the corners stayed corners, kept naming the dropped
/// wall, and the assembly refused the whole boolean (`VertexNamesAbsentSurface` at the corner —
/// measured). The merge now carries the circles through, and this pins the result the same way the
/// pass's own charter should have: right volume, clean `validate`, watertight mesh.
///
/// ★ **The merged face holds two circle holes** (one bore per cell), so the owner-assignment loop
/// actually runs rather than falling out on a single candidate.
#[test]
fn two_bored_plates_fuse_face_to_face() {
    let cell = |m: &mut Model, x0: f64| {
        let plate = m.add_cuboid(
            Point3::from_array([x0, 0.0, 0.0]),
            Point3::from_array([x0 + 20.0, 20.0, 10.0]),
        );
        m.rebuild_adjacency();
        let bore = m.add_cylinder(
            Point3::from_array([x0 + 14.0, 14.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            12.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    };
    let mut m = Model::new();
    let a = cell(&mut m, 0.0);
    let b = cell(&mut m, 20.0);
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the plates fuse");
    assert_eq!(out.len(), 1, "one solid");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let want = 2.0 * (20.0 * 20.0 * 10.0 - std::f64::consts::PI * 4.0 * 10.0);
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
    // The corner on the dropped wall is gone: no result vertex sits at (20, 20, 0).
    let reach = m.reachable();
    assert!(
        !reach.vertices.iter().any(|&vh| {
            let p = m.vertex_point(vh).as_array();
            (p[0] - 20.0).abs() < 1e-9 && (p[1] - 20.0).abs() < 1e-9 && p[2].abs() < 1e-9
        }),
        "the straight-angle corner on the shared wall dissolved"
    );
    // And the two bores' rims both survived as holes of the merged caps.
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    assert_eq!(
        uses.values().filter(|&&n| n != 2).count(),
        0,
        "the mesh is watertight"
    );
}

/// ★★★ **A 2×2 grid: an array fused in x, then moved in y and fused again.** The second
/// generation is what makes this hard — the first fuse leaves a body whose surfaces come from
/// *two* provenances (the original cell, and the one that moved), so moving that result puts
/// carriers with different chains on one corner. The exact-corner road used to want one shared
/// frame and declined, and the population gate's face test then could not judge the face at all
/// (`CylinderGateUndecided`, measured on the user's script). Each carrier states the world, so
/// the corner solves there.
///
/// The oracle is **proportionality**: n cells of one part must weigh exactly n times one cell.
/// That is what says nothing was lost or double-counted at the joins.
///
/// ★★★★★ **The cell's features are not decoration — they are what makes this fixture go red.**
/// A plain plate, or a plate with a bore or two, fuses into a 2×2 grid *without* the world road:
/// its second-generation corners never land where the population gate has to judge them. The
/// first fixture written for this lock was exactly that, and it passed with the road switched
/// off — measuring nothing. These three features (one pocket, two through-bores, taken from the
/// user's own part) are the smallest set measured to refuse `CylinderGateUndecided` with the road
/// off and build with it on. The relation is **not monotone** — the user's cell cut down to two
/// pockets and four bores *builds*, while one of those pockets with two of the other bores (the
/// pair here) refuses — so treat the numbers as pinned: changing one moves the fixture out of the
/// population it exists to hold.
#[test]
fn a_two_by_two_grid_fuses() {
    use nacre_scalar::Rat;
    let cell = |m: &mut Model, at: [f64; 3]| {
        let plate = m.add_cuboid(
            Point3::from_array(at),
            Point3::from_array([at[0] + 86.0, at[1] + 86.0, at[2] + 71.5]),
        );
        m.rebuild_adjacency();
        let pocket = m.add_cuboid(
            Point3::from_array([at[0] + 1.0, at[1] + 38.6, at[2]]),
            Point3::from_array([at[0] + 33.2, at[1] + 58.6, at[2] + 68.5]),
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
        m.rebuild_adjacency();
        let bore = m.add_cylinder(
            Point3::from_array([at[0] + 17.1, at[1] + 8.9, at[2]]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.22,
            143.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore).expect("bore 1")[0];
        m.rebuild_adjacency();
        let bore2 = m.add_cylinder(
            Point3::from_array([at[0] + 46.75, at[1] + 37.35, at[2]]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.34,
            143.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore2).expect("bore 2")[0];
        m.rebuild_adjacency();
        out
    };
    let shift = |m: &mut Model, s, o: [f64; 3]| {
        let s = transform(
            m,
            s,
            &nacre_scalar::Isometry::translation(o.map(|c| Rat::try_from_f64(c).unwrap())),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    };
    // One cell, for the oracle.
    let one = {
        let mut m = Model::new();
        let c = cell(&mut m, [0.0; 3]);
        nacre_props::mass_props(&m, c).unwrap().volume
    };
    // The grid: (uc ∪ uc→x) ∪ (that ∪ →y).
    let mut m = Model::new();
    let a = cell(&mut m, [0.0; 3]);
    let b = cell(&mut m, [0.0; 3]);
    let b = shift(&mut m, b, [86.0, 0.0, 0.0]);
    let row = boolean(&mut m, BoolKind::Fuse, a, b).expect("the x pair fuses")[0];
    m.rebuild_adjacency();
    let c = cell(&mut m, [0.0; 3]);
    let d = cell(&mut m, [0.0; 3]);
    let d = shift(&mut m, d, [86.0, 0.0, 0.0]);
    let row2 = boolean(&mut m, BoolKind::Fuse, c, d).expect("the second row fuses")[0];
    m.rebuild_adjacency();
    // ★ The second generation: a row that already mixes two provenances, moved again.
    let row2 = shift(&mut m, row2, [0.0, 86.0, 0.0]);
    assert!(carries_motion(&m, row2), "the move records a chain");
    let out = boolean(&mut m, BoolKind::Fuse, row, row2).expect("the grid fuses");
    assert_eq!(out.len(), 1, "one solid");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!(
        (v - 4.0 * one).abs() <= 1e-9 * one,
        "{v} vs {} (4 x {one})",
        4.0 * one
    );
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    assert_eq!(
        uses.values().filter(|&&n| n != 2).count(),
        0,
        "the mesh is watertight"
    );

    // ★★ **The corners this opened are a datum's carriers too, so the ledger says so here.**
    // The grid's vertices split: some answer in the world (their carriers each state it), some in
    // a shared frame — both must be present, or the fixture is not the mixed body it claims to be.
    // A datum through three of the world-answered ones then gets a **name**, which is the road
    // `ThroughStatement` takes when the three meets agree on a frame. Before this cell the whole
    // grid could not be built, so this capability has no earlier behaviour to change — it is new,
    // and cheap to hold here on a model that is already standing.
    let mut world: Vec<Handle<nacre_topo::Vertex>> = Vec::new();
    let mut framed = 0usize;
    {
        let sol = m.solids.get(out[0]).clone();
        let mut seen: Vec<Handle<nacre_topo::Vertex>> = Vec::new();
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                for he in &m.faces.get(fh).outer.half_edges {
                    let vh = m.he_start(*he);
                    if seen.contains(&vh) {
                        continue;
                    }
                    seen.push(vh);
                    match m.vertex_meet(vh) {
                        Some((_, None)) => world.push(vh),
                        Some((_, Some(_))) => framed += 1,
                        None => {}
                    }
                }
            }
        }
    }
    world.sort_by_key(|v| v.index());
    assert!(
        world.len() >= 3 && framed > 0,
        "a two-generation grid answers some corners in the world and some in a frame \
         (world {}, framed {framed})",
        world.len()
    );
    let Ok(OpOutput::DatumPlane { plane, .. }) = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices([world[0], world[1], world[2]]),
        },
    ) else {
        panic!("a datum through three world-answered corners")
    };
    assert!(
        m.surface_name.contains_key(&plane),
        "the datum is named, not judged"
    );
}

/// ★★★★ **A part does not care which side of itself a boss sits on.**
///
/// With the axis exactly on a wall, the `+x` and `+y` walls built and the `-x` and `-y` walls
/// answered `RingOrientation` — a **SuspectedDefect** name, on an input a script writes by
/// accident. The cause was one cross product reading two spellings of the same fact: the chord's
/// sense is made against **stored** frames, and the class's *canonical* name opposes it on half
/// the classes, so the product was right only where the two agreed. On the other half the chord's
/// sense flipped, one cell of six (the overhang below the plate) came back wound the other way,
/// and the walk found two outer contours where a one-component class has one.
///
/// ★ **The lock the code asked for.** That site's own comment said the `wc` half "is still
/// lock-invisible … so that half stands on the convention argument, recorded rather than
/// assumed". This is the fixture that sees it: the four walls have to weigh the same, and the
/// value is hand-derived rather than copied from a run.
#[test]
fn a_boss_on_any_wall_weighs_the_same() {
    // plate 4×4×2 = 32; the boss (r = 0.5, h = 4) runs clear through it, half of it inside the
    // material, so the shared part is half a cylinder of the plate's height.
    let plate = 32.0;
    let boss = std::f64::consts::PI * 0.25 * 4.0;
    let shared = 0.5 * std::f64::consts::PI * 0.25 * 2.0;
    let run = |base: [f64; 3], kind: BoolKind| -> (f64, usize) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
        m.rebuild_adjacency();
        let out = boolean(&mut m, kind, a, b).expect("a boss on a wall builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight"
        );
        let sol = m.solids.get(out[0]).clone();
        let faces = std::iter::once(&sol.outer)
            .chain(sol.cavities.iter())
            .map(|&sh| m.shells.get(sh).faces.len())
            .sum();
        (nacre_props::mass_props(&m, out[0]).unwrap().volume, faces)
    };
    let mut fuse_faces: Vec<usize> = Vec::new();
    // ★★ **All four walls, all three kinds.** `+y`'s cut and common used to answer `MissingSeam`
    // and stood here as an exception: the cap's wrap arc had to be split at its rim's seam
    // vertex, and that vertex existed only when a *band* happened to claim the same rim — which
    // on this wall nothing does. Registering a rim from the arc that needs it (`wrapping_rim`)
    // retired the exception; the two rows are exact against the same closed forms as the others.
    for base in [
        [0.0, 2.0, -1.0], // -x wall
        [4.0, 2.0, -1.0], // +x
        [2.0, 0.0, -1.0], // -y
        [2.0, 4.0, -1.0], // +y
    ] {
        let (v, f) = run(base, BoolKind::Fuse);
        assert!(
            (v - (plate + boss - shared)).abs() < 1e-12,
            "fuse {base:?}: {v}"
        );
        fuse_faces.push(f);
        let (v, _) = run(base, BoolKind::Cut);
        assert!((v - (plate - shared)).abs() < 1e-12, "cut {base:?}: {v}");
        let (v, _) = run(base, BoolKind::Common);
        assert!((v - shared).abs() < 1e-12, "common {base:?}: {v}");
    }
    // ★★ **The mirror invariant a volume cannot state.** Four congruent solids must also be built
    // out of the same number of faces: a volume can come out right while a face is split or merged
    // away, and a face wound the wrong way is exactly the shape of what went wrong here.
    assert!(
        fuse_faces.windows(2).all(|w| w[0] == w[1]),
        "the four walls give congruent solids, so their face counts must agree: {fuse_faces:?}"
    );
    // ★ The control that says the defect needed the cylinder: the same straddle with a box tool
    // built on every wall before this fix and must go on doing so.
    for (lo, hi) in [
        ([-0.5, 1.5, -1.0], [0.5, 2.5, 3.0]),
        ([3.5, 1.5, -1.0], [4.5, 2.5, 3.0]),
        ([1.5, -0.5, -1.0], [2.5, 0.5, 3.0]),
        ([1.5, 3.5, -1.0], [2.5, 4.5, 3.0]),
    ] {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the box straddle builds");
        m.rebuild_adjacency();
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - 35.0).abs() < 1e-12, "box {lo:?}: {v}");
    }
}

/// ★★★★★ **The eye that was missing: does the mesh actually cover the face it approximates?**
///
/// A boolean result was once shipped whose lateral mesh lost **16% of its area** — four triangles
/// spanned half a turn of the cylinder as flat chords, and the boss rendered as two cones. Every
/// standing instrument was green at the time: `validate` clean, the mesh watertight, the **exact**
/// volume right (`mass_props` integrates the b-rep, not the mesh), the face counts as predicted.
/// Nothing asked the one question that mattered — whether the triangles lie on the surface they
/// claim — so nothing answered it.
///
/// Two readings, and **neither is derived from the other**: a face's triangles against
/// [`nacre_props::face_props`] (an exact boundary integral) and the whole solid's triangles against
/// [`nacre_props::mass_props`] (the divergence theorem on the b-rep).
///
/// **The budget is derived, not chosen.** A circle sampled into `n` chords is approximated by the
/// inscribed regular `n`-gon, whose area is `sinc(2π/n) = 1 − (2π/n)²/6` of the true one; a
/// cylinder's lateral loses the arc-versus-chord ratio `sinc(π/n) = 1 − (π/n)²/6`, four times
/// less. Both have `n ≥ 360°/Δθ = 180` (`circle_segments`' angular term holds whatever the
/// radius), so the worst is the disk's **2.0e-4**, and a **relative 1e-3** leaves five times that
/// while sitting two hundred times below the defect this exists to catch. The error takes
/// **either sign**: a plate whose circular bite is inscribed comes out slightly *larger*.
#[test]
fn the_mesh_covers_the_faces_it_approximates() {
    // 5 × the derived worst case (2.0e-4, a full disk at the angular budget).
    const BUDGET: f64 = 1e-3;
    let check = |name: &str, m: &Model, solids: &[Handle<Solid>]| {
        let mesh = nacre_tess::tessellate(m, &nacre_tess::TessConfig::default()).expect("tess");
        for &s in solids {
            let sol = m.solids.get(s).clone();
            let mut mesh_volume = 0.0f64;
            for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
                for fh in m.shells.get(sh).faces.clone() {
                    let exact = nacre_props::face_props(m, fh)
                        .unwrap_or_else(|e| panic!("{name}: face_props {e:?}"))
                        .area;
                    let mut area = 0.0f64;
                    for &th in mesh.by_face.get(&fh).map(|v| v.as_slice()).unwrap_or(&[]) {
                        let t = mesh.triangles.get(th);
                        let p: Vec<_> = t
                            .vertices
                            .iter()
                            .map(|&h| mesh.vertices.get(h).pos)
                            .collect();
                        area += (p[1] - p[0]).cross(p[2] - p[0]).norm() * 0.5;
                        // The signed volume of the tetrahedron on the origin; summed over an
                        // outward-oriented closed mesh it is the volume that mesh encloses.
                        let (a, b, c) = (p[0].as_array(), p[1].as_array(), p[2].as_array());
                        mesh_volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
                            - a[1] * (b[0] * c[2] - b[2] * c[0])
                            + a[2] * (b[0] * c[1] - b[1] * c[0]))
                            / 6.0;
                    }
                    assert!(
                        (area - exact).abs() <= BUDGET * exact,
                        "{name}: a face's mesh area {area} is not its own {exact}"
                    );
                }
            }
            let exact = nacre_props::mass_props(m, s).expect("mass_props").volume;
            assert!(
                (mesh_volume - exact).abs() <= BUDGET * exact,
                "{name}: the mesh encloses {mesh_volume}, the solid is {exact}"
            );
        }
    };
    let plate = |m: &mut Model| {
        m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        )
    };
    // The plate alone is the instrument's own control: nothing is curved, so both readings are
    // exact to rounding and a failure here would be the instrument's, not the mesh's.
    {
        let mut m = Model::new();
        let a = plate(&mut m);
        m.rebuild_adjacency();
        check("plate", &m, &[a]);
    }
    for (name, base, h) in [
        ("wall -x", [0.0, 2.0, -1.0], 4.0),
        ("wall +x", [4.0, 2.0, -1.0], 4.0),
        ("wall -y", [2.0, 0.0, -1.0], 4.0),
        ("wall +y", [2.0, 4.0, -1.0], 4.0),
        ("corner", [4.0, 4.0, -1.0], 4.0),
        ("through", [2.0, 2.0, -1.0], 4.0),
        ("on top", [2.0, 2.0, 2.0], 1.0),
        ("flush", [2.0, 2.0, 0.0], 3.0),
        // The half-height wall boss: its upper cap sits inside the plate (D2b opened it), and
        // its mirror with the lower cap inside — one lateral face each since D4, whose chain
        // rim passes the seam on a **wrap arc** here (the seam is at −y, outside the plate),
        // where the boss on the x = 40 wall in `bands` runs *along* the seam ruling.
        ("half wall", [2.0, 0.0, -1.0], 2.0),
        ("half wall, cap below", [2.0, 0.0, 1.0], 2.0),
    ] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let mut m = Model::new();
            let a = plate(&mut m);
            let b = m.add_cylinder(
                Point3::from_array(base),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                h,
            );
            m.rebuild_adjacency();
            // A population this ladder does not build yet is not this test's business; what it
            // *does* build has to be drawn correctly.
            let Ok(out) = boolean(&mut m, kind, a, b) else {
                continue;
            };
            if out.is_empty() {
                continue;
            }
            m.rebuild_adjacency();
            check(&format!("{name} {kind:?}"), &m, &out);
        }
    }
}

/// ★★★ **A boss standing on a wall has ONE lateral face, not three.**
///
/// The band pass emits a lateral surface in as many pieces as the arrangement cut it into: the
/// band under the plate, the half-band beside it, the band above. The two circles between those
/// pieces bound **nothing** — the surface runs smooth across them — so drawn they are a line
/// ringing a boss that has none, which is what a user reported seeing. `unify_curved_faces` erases
/// them, and what is left is a band with one notch punched out of its side.
///
/// The negative controls are the point of the test: a boss standing on the middle of the plate
/// really *is* two lateral faces (the plate interrupts it), and a bore really is one. Both must
/// come out with their face counts untouched, or the pass is erasing boundaries that exist.
#[test]
fn a_boss_on_a_wall_has_one_lateral_face() {
    let run = |base: [f64; 3], h: f64, kind: BoolKind| -> (f64, usize, usize, Vec<usize>) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            h,
        );
        m.rebuild_adjacency();
        let out = boolean(&mut m, kind, a, b).expect("the boss builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "validate {base:?}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight {base:?}"
        );
        let sol = m.solids.get(out[0]).clone();
        let (mut total, mut lateral, mut holes) = (0usize, 0usize, Vec::new());
        for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
            for fh in m.shells.get(sh).faces.clone() {
                total += 1;
                let f = m.faces.get(fh);
                if matches!(m.surface(f.surface), nacre_geom::Surface::Cylinder(_)) {
                    lateral += 1;
                    holes.push(f.inner.len());
                }
            }
        }
        (
            nacre_props::mass_props(&m, out[0]).unwrap().volume,
            total,
            lateral,
            holes,
        )
    };
    // The plate is 32; the boss (r = 0.5, h = 4) runs clear through it with half its section
    // buried, so a wall boss weighs `32 + π/4·4 − ½·π/4·2` and a corner one buries a quarter.
    let quarter = std::f64::consts::PI * 0.25;
    // ★★ **The notch is spelled one of two ways, and which one is not a choice.** Where it meets
    // the seam generator the outer walk bridges it in (`band_loop`) and the face carries **no**
    // inner loop; where it misses the seam it stays an honest hole. `ref_dir` is a function of the
    // axis alone, so four congruent walls stand in four different relations to θ = 0 — and that
    // asymmetry has to show up **nowhere else**, which is what the rest of this loop pins: one
    // lateral face and exactly two faces fewer, on every one of them.
    for (name, base, faces, notch_holes) in [
        ("-x", [0.0, 2.0, -1.0], 10, 0),
        ("+x", [4.0, 2.0, -1.0], 10, 0),
        ("-y", [2.0, 0.0, -1.0], 10, 1),
        ("+y", [2.0, 4.0, -1.0], 10, 0),
    ] {
        let (v, total, lateral, holes) = run(base, 4.0, BoolKind::Fuse);
        assert!(
            (v - (32.0 + quarter * 4.0 - 0.5 * quarter * 2.0)).abs() < 1e-12,
            "{name}: {v}"
        );
        assert_eq!(lateral, 1, "{name}: one lateral face");
        assert_eq!(holes, vec![notch_holes], "{name}: the notch's spelling");
        assert_eq!(
            total, faces,
            "{name}: exactly two faces fewer than the three pieces"
        );
    }
    // The corner: two walls cut the cylinder, so the surviving panel's four corners name **two
    // different** wall classes — and the merge is the same one.
    let (v, total, lateral, holes) = run([4.0, 4.0, -1.0], 4.0, BoolKind::Fuse);
    assert!(
        (v - (32.0 + quarter * 4.0 - 0.25 * quarter * 2.0)).abs() < 1e-12,
        "corner: {v}"
    );
    assert_eq!((lateral, holes, total), (1, vec![0], 9), "corner");
    // ★★ **The negative controls — seams that are real must survive untouched.**
    for (name, base, h, kind, lateral, total) in [
        // The plate genuinely interrupts this one: two lateral faces, and neither has a hole.
        ("through boss", [2.0, 2.0, -1.0], 4.0, BoolKind::Fuse, 2, 10),
        ("bore", [2.0, 2.0, -1.0], 4.0, BoolKind::Cut, 1, 7),
        ("boss on top", [2.0, 2.0, 2.0], 1.0, BoolKind::Fuse, 1, 8),
        ("flush boss", [2.0, 2.0, 0.0], 3.0, BoolKind::Fuse, 1, 8),
    ] {
        let (_, t, l, holes) = run(base, h, kind);
        assert_eq!((l, t), (lateral, total), "{name} must not move");
        assert!(holes.iter().all(|&n| n == 0), "{name}: no notch");
    }
}

/// A plate, one cylinder op, then a second — the population the corpus did not have.
pub(crate) fn chained(
    first: ([f64; 3], f64, BoolKind),
    second: ([f64; 3], f64, BoolKind),
) -> (Model, Result<Vec<Handle<Solid>>, BoolError>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array(first.0), up, 0.5, first.1);
    m.rebuild_adjacency();
    let out = boolean(&mut m, first.2, plate, a).expect("the first op builds");
    assert_eq!(out.len(), 1, "the first op is one solid");
    m.rebuild_adjacency();
    let b = m.add_cylinder(Point3::from_array(second.0), up, 0.5, second.1);
    m.rebuild_adjacency();
    let r = boolean(&mut m, second.2, out[0], b);
    (m, r)
}

/// ★★★★ **Two cylinder operations in a row — measured for the first time.**
///
/// The corpus had **no chained-cylinder fixture at all**: every boss and bore stood on a plain
/// plate. So "a plate takes a second cylinder" was neither locked nor known, and it turned out to
/// be *half* true — this is the half that always built, kept as its own regression guard.
///
/// A bore or a boss standing on the plate's face leaves the operand's rings plane-named (a hole is
/// a full circle, the one curved loop the tracer's road speaks), so the next cylinder is business
/// as usual. A boss standing on a **wall** did not, for six wall-names running; it builds now, and
/// [`a_chained_cylinder_bounded_by_the_first_builds`] holds that half.
///
/// The volumes are derived, not copied: the plate is 96, a bore through it removes `π/4·2`, a boss
/// on top adds `π/4·1`, and a wall boss adds `π/4·4` with half its section buried (`−½·π/4·2`).
#[test]
fn chained_cylinder_operations_that_build_today_still_build() {
    let quarter = std::f64::consts::PI * 0.25;
    let bore = ([2.0, 2.0, -1.0], 4.0, BoolKind::Cut);
    let top_boss = ([2.0, 2.0, 2.0], 1.0, BoolKind::Fuse);
    for (name, first, second, volume) in [
        (
            "bore then bore",
            bore,
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            96.0 - 2.0 * quarter * 2.0,
        ),
        (
            "bore then top boss",
            bore,
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            96.0 - quarter * 2.0 + quarter,
        ),
        (
            "top boss then top boss",
            top_boss,
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            96.0 + 2.0 * quarter,
        ),
        (
            "top boss then wall boss",
            top_boss,
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            96.0 + quarter + (quarter * 4.0 - 0.5 * quarter * 2.0),
        ),
        (
            "bore then wall boss",
            bore,
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            96.0 - quarter * 2.0 + (quarter * 4.0 - 0.5 * quarter * 2.0),
        ),
    ] {
        let (mut m, r) = chained(first, second);
        let out = r.unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(out.len(), 1, "{name}: one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "{name}: validate");
        let v = nacre_props::mass_props(&m, out[0]).expect("mass").volume;
        assert!((v - volume).abs() < 1e-9, "{name}: {v} vs {volume}");
    }
    // ★ A wall boss as a *middle* operation used to be the one that could not follow. It builds
    // now; its volumes live in `a_chained_cylinder_bounded_by_the_first_builds`, beside the
    // mechanism that opened it.
}

/// ★★★★★ **The gate decides a boss-seated wall, and records the pair — and nothing else in the
/// crate would notice if it stopped.**
///
/// The clearance test's whole content is a `bool`, and on a chained operand its answer travels
/// exactly one place: a `false` sends the seated pair to the `d = 0` record-and-pass arm, which
/// puts it in `crossings`, which is what makes the tracer's ruling and chord arms fire at all. A
/// wrong `true` would refuse nothing, panic nowhere and change no output — today's population
/// declines a step later either way (`BranchNode`, the lock below) — so **the suite would stay
/// green with the branch corner's answer thrown away entirely**. Measured, by throwing it away:
/// 316 green. This is the lock that sees it.
///
/// ★★ **What it sees is that the corner is *readable*, not what it says** — recorded rather than
/// claimed. Forcing every branch corner to one strip side leaves even this lock green, because the
/// verdict here is settled on the *other* separating axis: the wall's plate corners already sit
/// inside the boss's axial span, so `along` fails whatever the strip half answers. The corner's
/// side becomes decisive only for a face that clears along the axis and has to be judged across
/// it, and no fixture builds one yet. Removing the read itself is what turns this red
/// (`CurvedOperandBoundary`, measured).
#[test]
fn the_gate_records_a_wall_the_boss_is_seated_on() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
    m.rebuild_adjacency();
    let b = m.add_cylinder(Point3::from_array([6.0, 2.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    // ★ The gate *deciding* is half the claim — it used to stop at the first corner the boss made.
    let setup = crate::planes::plane_index_setup(&m, first, b).expect("the gate decides");
    // ★★ And recording is the other half. The boss sits **on** `y = 0`, so those faces are not
    // clear of it; that is what the record-and-pass arm is for, and a missing entry here means the
    // clearance test answered "clear" about a wall a cylinder is standing in.
    //
    // ★ The pair is checked by **what it means**, not by its indices: the arm records exactly the
    // classes whose plane the cylinder's axis lies *on*, so that is the assertion — renumbering
    // the classes cannot make it pass or fail for the wrong reason.
    assert_eq!(setup.crossings.len(), 1, "one seated pair, not more");
    for &(c, ci) in &setup.crossings {
        let coeffs = setup.geom[c].world_rat.expect("a named wall class");
        assert_eq!(
            nacre_scalar::point_plane_clearance_rat(
                &coeffs,
                &setup.cyls[ci].def.origin(),
                nacre_scalar::Rat::from_int(0)
            ),
            nacre_scalar::Orient::Zero,
            "the recorded pair is a cylinder seated exactly on that class's plane"
        );
    }
}

/// ★★★★★ **The rule that replaced the two plane fences is *differenced* against them, not argued
/// equal to them.**
///
/// They are the same predicate wherever the fence is well posed — the fence plane is the class that
/// pins one end, so it crosses the line exactly at that end and "the side the far endpoint is on" is
/// the half-line from there. That argument has a hole: where the far endpoint sits *on* the fence
/// plane, `want` is `0` and every nonzero side reads as outside. So the arc split computes both
/// verdicts for every crossing it examines and counts the disagreements, and this reads the counter.
///
/// ★★★★★ **Three negative controls, and one of them says the green is thinner than it looks.**
/// ☑ Inverting `closed_contains`' verdict inside the probe: **red, 8 of 8** — both roads are really
/// read. ☑ Letting the fence side check only *one* of its two fences: **still green**. ☑ Letting it
/// check **neither**, so it answers "inside" unconditionally: **still green** — which says that on
/// this fixture *no crossing is out of extent*, so the two predicates agree by both saying yes. The
/// case that would separate them — a segment ending strictly inside the disk, whose line's far
/// crossing lies past its end — is what a **chained** operand produces, and that population is still
/// behind the gate. So: this locks that the roads are wired to the same question, and the census is
/// what locks the answers.
///
/// ★ The disagreement half is asserted **globally**, not as a delta: any fixture anywhere in this
/// binary that does reach an out-of-extent crossing has to agree too, whatever order it runs in.
///
/// ★ It expires with the gate: the fences need both endpoints' rational coordinates, which is the
/// demand the next rung removes.
#[test]
fn the_extent_rule_agrees_with_the_fences_it_replaced() {
    use crate::arrangement::extent_probe::{ASKED, DISAGREED};
    use core::sync::atomic::Ordering::Relaxed;
    let before = (ASKED.load(Relaxed), DISAGREED.load(Relaxed));
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let boss = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
    let after = (ASKED.load(Relaxed), DISAGREED.load(Relaxed));
    assert!(
        after.0 > before.0,
        "the probe never ran, so it measured nothing: asked {} -> {}",
        before.0,
        after.0
    );
    assert_eq!(
        after.1, 0,
        "the line-order rule and the plane fences disagreed on {} of the {} crossings this binary \
         has examined",
        after.1, after.0
    );
}

/// ★★★★★ **The order rule may read the retired ruler backwards, but it may never reshuffle it.**
///
/// The arc split used to sort its points by their parameter on `branch_meet`'s canonical meet line;
/// it asks [`combinatorics::order_located`] now. The two do **not** agree pointwise — the ruler's
/// direction comes from the classes' rational coefficients and the rule's axis sign from the judge's
/// stored planes, and those two spellings name the same plane without naming the same side. A
/// wholesale reversal is harmless: the pieces come out in the comparator's own ascending order and
/// `forward` is taken with that same comparator, so the two cancel. A **partial** disagreement would
/// not cancel — it would put sub-segments between the wrong pairs of points — and the bit census
/// cannot see it, because it sorts before it hashes.
///
/// So the proposition is per segment: **the same order, or exactly its reverse, and never in
/// between.** ☑ Flipping one comparison inside the probe turns this red.
#[test]
fn the_order_rule_never_reshuffles_the_ruler_it_replaced() {
    use crate::arrangement::order_probe::{
        EQUALITY_DISAGREED, REVERSED, SCRAMBLED, SEGMENTS, WC_ABOVE_WALL, WC_BELOW_WALL,
    };
    use core::sync::atomic::Ordering::Relaxed;
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let boss = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
    assert!(
        SEGMENTS.load(Relaxed) > 0,
        "the probe never ran, so it measured nothing"
    );
    assert_eq!(
        SCRAMBLED.load(Relaxed),
        0,
        "{} of the {} segments this binary split came out neither in the ruler's order nor exactly \
         reversed ({} were reversed)",
        SCRAMBLED.load(Relaxed),
        SEGMENTS.load(Relaxed),
        REVERSED.load(Relaxed)
    );
    assert_eq!(
        EQUALITY_DISAGREED.load(Relaxed),
        0,
        "the two roads disagreed {} times about whether two split points are the same place — the \
         `CoincidentNodes` refusal is reachable from one and not the other",
        EQUALITY_DISAGREED.load(Relaxed)
    );
    // ★★★ **Both plane orders are exercised, and the invariant above is why neither is "the right
    // one".** ☑ Measured over the whole binary: **36 of 256** segments have `wc > wall`, where the
    // sorted pair and the call-order pair are opposite calls. Swapping the pair flips *every*
    // comparison on a segment, so it can only turn "same" into "reversed" — which the assertion
    // above already says is harmless. The sorted pair is chosen for agreeing with
    // `NodeId::Branch`'s convention, not because a fixture prefers it.
    //
    // ★ That count is **recorded, not asserted**: these counters accumulate across the binary, and
    // this test cannot know what has run before it. Only the two "never" claims above are safe to
    // assert from here. The relation is read off `wc` and `wall`, which nothing in this cell moves.
    let _ = (WC_ABOVE_WALL.load(Relaxed), WC_BELOW_WALL.load(Relaxed));
}

/// ★★★★★ **A holed lateral's ruling comes out in three pieces, and the middle one grazes.**
///
/// A wall through the boss's axis meets its lateral in two rulings, and where the boss is buried in
/// the plate the lateral does not *cross* that wall — it ends at it. So the mark is
/// `Transversal · Graze · Transversal` along the axis, not one full-height crossing.
///
/// ★★ **Held at the trace, because there is no result to open**: the operation this exercises still
/// refuses further down (the assembly's own incompleteness), so the face count that would show the
/// ghost wall face gone cannot be read. The trace's own answer is what survives, the way
/// `arrangement`'s `sides == [-1, 1]` lock already holds the hole-free case.
///
/// ★ The claim is over **every** carved ruling this binary produces, in whatever order the tests
/// ran — all of them come from the wall-boss family.
#[test]
fn a_holed_laterals_ruling_grazes_where_the_hole_is() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
    m.rebuild_adjacency();
    let b = m.add_cylinder(Point3::from_array([6.0, 2.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    let _ = boolean(&mut m, BoolKind::Cut, first, b);
    // ★ The ledger is the binary's: since E2-2 a panel's ruling (one graze) and a chain's (a
    // transversal then a graze) are recorded too. Read this fixture's boss — origin `(2, 0, −1)`,
    // rulings spanning `t ∈ [0, 4]` — and nothing else.
    let carved: Vec<Vec<crate::arrangement::SegKind>> = crate::arrangement::ruling_probe::CARVED
        .lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        .filter(|c| c.origin == [2.0, 0.0, -1.0] && c.span == [0.0, 4.0])
        .map(|c| c.kinds.clone())
        .collect();
    assert!(
        !carved.is_empty(),
        "the probe never ran, so it measured nothing"
    );
    for kinds in &carved {
        assert!(
            matches!(
                kinds[..],
                [
                    crate::arrangement::SegKind::Transversal { .. },
                    crate::arrangement::SegKind::Graze { .. },
                    crate::arrangement::SegKind::Transversal { .. }
                ]
            ),
            "a carved ruling came out as {kinds:?}"
        );
    }
    // ★ Both rulings are carved, not just one: the hole has a vertical edge on each.
    assert!(carved.len() >= 2, "only {} ruling carved", carved.len());
    // ★★★★★ **And which side the face occupies, against the fixture's own geometry.** The buried
    // half of the boss is the `y > 0` one, so along the hole's vertical edges the lateral survives
    // at `y < 0` — the stored-normal side exactly when that normal points at `−y`. ☑ Flipping any
    // factor of the derivation turns this red; nothing downstream does, yet.
    let sides: Vec<(bool, f64)> = crate::arrangement::ruling_probe::GRAZE_SIDE
        .lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        // ★ Only the **hole's** runs: the same boss's Cut result (the census's) puts a *panel* on
        // the same wall and span, and its face lies on the other side — measured (`body_above`
        // false against the hole's true), which is exactly what a filter by origin alone let in.
        .filter(|g| {
            g.kind == combinatorics::CycleKind::Hole
                && g.origin == [2.0, 0.0, -1.0]
                && g.span[1] - g.span[0] == 2.0
        })
        .map(|g| (g.body_above, g.ny))
        .collect();
    assert!(!sides.is_empty(), "the side probe never ran");
    for &(body_above, ny) in &sides {
        assert!(
            ny.abs() > 0.5,
            "the wall's stored normal is not axis-aligned: {ny}"
        );
        assert_eq!(
            body_above,
            ny < 0.0,
            "a hole's ruling graze states the wrong side of the wall: {body_above} against ny={ny}"
        );
    }
}

/// ★★★★★ **Which of the two arcs is the hole — measured, not argued.**
///
/// A plane cuts a circle in two points, and the whole difficulty of a lateral face's hole is which
/// of the two arcs it is: the hole's two ruling edges lie on **one** plane here, so that plane's
/// sides cannot tell them apart. `cycle_on_class` answers from the ring's own winding (material is
/// on the left of the ring's travel), and gets there through three signs — the face's
/// `orient_sign`, the class's `frame_sign`, and whether the class's stored normal points along the
/// axis. Reversed, the answer is not approximate but **exactly opposite**: the graze would sit on
/// the band and the crossing inside the hole.
///
/// ★★ **Nothing downstream notices yet** — these fixtures stop at `loop_winding` before
/// `label_cells` could judge a label — so the arc is realized and judged here instead. The kernel
/// reads no coordinate to choose it; this reads one afterwards.
///
/// ☑ **Which factors this actually sees.** Turning each one off in turn: reversing either arm's
/// winding is **red**, dropping `plus_t_is_above` is **red**, and dropping the root restatement
/// (the ruling named for a different ⊥ plane) is **red**. Dropping `frame_sign` or the face's
/// `orient_sign` is **green** — both are `+1` everywhere in today's holed population (an
/// `add_cuboid` face is never `Reversed`, and a boss's lateral faces outward), so this lock cannot
/// see them and their reasons stand on the derivation alone. A **bore** with a hole in its lateral
/// is what would exercise them.
///
/// ★ The claim is over **every** arc this binary carves, in whatever order the tests ran: all of
/// them come from the wall-boss family, where the plate stands on `y > 0` and the boss is centred
/// on the wall `y = 0`, so the buried half — the hole — is the `y > 0` one. A fixture that buries
/// a lateral somewhere else would need its own reading of "which side is the hole", and this
/// assertion is where that would show up.
#[test]
fn a_holes_arcs_run_the_way_the_hole_lies() {
    // ★ Two second operands, because they reach **different arms** of the walk. A cylinder rising
    // from `z = −1` puts its cap on `z = 0`, which is the hole's own low rim — an on-line *run*.
    // One rising from `z = 1` puts a cap **strictly inside** the hole's `z ∈ (0, 2)`, where the
    // ring crosses the class on its two ruling edges instead — the *crossing* arm, whose extra
    // step is restating the ruling's root for a different ⊥ plane.
    for base in [-1.0, 1.0] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([12.0, 4.0, 2.0]),
        );
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
        m.rebuild_adjacency();
        let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
        m.rebuild_adjacency();
        let b = m.add_cylinder(Point3::from_array([6.0, 2.0, base]), up, 0.5, 4.0);
        m.rebuild_adjacency();
        let _ = boolean(&mut m, BoolKind::Cut, first, b);
    }
    // ★ The ledger is the binary's, not this test's: every fixture that carves a cycle writes to
    // it, and since E2-2 a panel's and a chain's arcs (which face any way) are carved too. Read
    // only this fixture's boss (`(2, 0, −1)`) and only its **hole** — the sentence is about holes.
    let mids: Vec<[f64; 3]> = crate::arrangement::arc_probe::MIDS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        .filter(|m| m.kind == combinatorics::CycleKind::Hole && m.origin == [2.0, 0.0, -1.0])
        .map(|m| m.dir)
        .collect();
    assert!(
        !mids.is_empty(),
        "the probe never ran, so it measured nothing"
    );
    assert!(
        mids.iter().all(|d| d[1] > 0.0),
        "a hole's arc was stated running the wrong way round: midpoints {mids:?}"
    );
}

/// ★★★★★ **A boolean's own result takes a second cylinder — the whole road, end to end.**
///
/// A boss standing on a **wall** was the one chained operand that did not build. It buries half the
/// boss's lateral in the plate, so that face comes back as a band with a **hole**, and every layer
/// downstream had a sentence that was false about it. The wall moved six times as those were
/// closed one at a time — `TraceDeclined { BranchNode }`, `WitnessNotRational`,
/// [`RejectReason::CurvedStraightRun`], [`RejectReason::OpenResultShell`],
/// [`RejectReason::LabelConflict`] when one of two cancelling falsehoods was fixed, then
/// `OpenResultShell` again with the arrangement whole — and this lock was the map of that walk.
/// It is now the map of the far side.
///
/// ★★★★★ **The last false sentence was that a label answers two questions.** A label says where
/// *material* is; the panel road also read it as saying whether the lateral face is **there**. On
/// the first fusion the two sectors of the holed rim carry different labels (the buried one is
/// inside the *other* solid) and the road happened to be right; on a second operation the plate is
/// already own material, the two sectors' labels are **literally identical**, and both survived —
/// filling the hole in and leaving its four boundary edges claimed once each. The trace had said
/// which sector was a face all along, per arc, since the hole taught it to mark by angular extent:
/// `Graze` where the face ends, `Transversal` where it runs through. `bands::face_spans` asks it.
///
/// ☑ Measured across these four fixtures: **four** sectors dropped for existence, exactly one per
/// operation — the buried half — and the volumes below are what says that count is right.
///
/// ★ The volumes are derived, not copied: the plate is 96, the wall boss adds `π/4·4` with half its
/// section buried (`−½·π/4·2`), a through bore removes `π/4·2`, a boss on top adds `π/4·1`, and a
/// bore *on the wall* removes only the half-section the plate holds (`−½·π/4·2`).
#[test]
fn a_chained_cylinder_bounded_by_the_first_builds() {
    let q = std::f64::consts::PI * 0.25;
    let wall_boss = ([2.0, 0.0, -1.0], 4.0, BoolKind::Fuse);
    for (name, second, volume) in [
        (
            "bore",
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            96.0 + 3.0 * q - 2.0 * q,
        ),
        (
            "boss on top",
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            96.0 + 3.0 * q + q,
        ),
        (
            "another wall boss",
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            96.0 + 3.0 * q + 3.0 * q,
        ),
        (
            "a bore on the same wall",
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Cut),
            96.0 + 3.0 * q - q,
        ),
    ] {
        let (mut m, r) = chained(wall_boss, second);
        let out = r.unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(out.len(), 1, "{name}: one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "{name}: validate");
        let v = nacre_props::mass_props(&m, out[0]).expect("mass").volume;
        assert!((v - volume).abs() < 1e-9, "{name}: {v} vs {volume}");
    }
    // ★★★ **And that the existence gate is what did it.** The band road's `panel_probe` used to
    // record it here; since D3 the chart's census holds it where the cells are read:
    // `exist_marks_false` (a sector dropped because the face is not there — the four above, held
    // as growth ≥ 4 in `cyl_chart::tests::the_cells_read_their_chamber_from_the_horizontal_lines`)
    // and `arcs_no_mark`/`arcs_multi_mark` (every cut end read carries exactly one lateral mark of
    // its own solid — the doc that says why the `Seated` skip is load-bearing lives on `face_spans`).
}

/// ★★★★★ **A boolean's result carries names the tracer's plane-only road cannot read — and it
/// says so instead of aborting.**
///
/// A boss standing on a wall leaves the plate's caps bitten by an **arc** and the wall split in two
/// by the boss's **rulings**, so those faces' rings run along a cylinder. [`combinatorics`]'s ring
/// naming asks two plane-only questions of every edge — the carried wall class and the vertex's
/// three-plane name — and both go through `ClassIx::plane`, which **panics** on a cylinder. That
/// accessor is right to: for its forty-odd other callers a cylinder there is an upstream filter
/// bug. This is the caller whose input is an *operand*, so this is where the filter belongs.
///
/// ★ **Measured before the filter existed: this very call panicked** ("a plane-only path got
/// cylinder class 0"). The population gate refuses such an operand long before the tracer sees it,
/// so nothing in production reaches this today — which is exactly why the lock calls the function
/// **directly**, the same way the ring-naming goldens next door do. The gate opens in the cell that
/// builds the road; this is the net that has to be under it first.
///
/// The untouched walls are the negative control: they are still named, so the filter is a filter
/// and not a blanket refusal.
#[test]
fn an_operand_bounded_by_a_cylinder_is_named_in_class_space() {
    // ★★ **The corner boss is here because the wall boss is symmetric.** Its axis lies *on* the
    // wall, so the two rulings sit either side of that plane and a root read the wrong way round
    // lands on a mirror-image point — a lock that only ever saw this fixture could pass with the
    // restatement inverted. The corner's two walls break that symmetry.
    let (mut ups, mut sides): (Vec<bool>, Vec<i8>) = (Vec::new(), Vec::new());
    for at in [[2.0, 0.0, -1.0], [12.0, 0.0, -1.0]] {
        let (u, s) = named_in_class_space(at);
        ups.extend(u);
        sides.extend(s);
    }
    // ★★ **The direction assertions inside are only worth their ink if the fixtures move them** —
    // a boss whose rulings all climbed, or all sat one side, would let a constant pass for a
    // derivation. It takes **both** fixtures to move both bits, and that is the geometry rather
    // than a gap: the wall boss is symmetric about its wall, so its two rulings straddle that
    // plane and `side` takes both values; the corner boss keeps **one** ruling per wall (the
    // other is buried), and those two are measured against *different* planes, so nothing says
    // they should oppose. Measured — the corner run alone gives `[1, 1]`.
    assert!(
        ups.contains(&true) && ups.contains(&false),
        "the fixtures exercise both senses of `up`: {ups:?}"
    );
    assert!(
        sides.contains(&1) && sides.contains(&-1),
        "the fixtures exercise both sides of the axis plane: {sides:?}"
    );
}

fn named_in_class_space(at: [f64; 3]) -> (Vec<bool>, Vec<i8>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array(at),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
    m.rebuild_adjacency();
    let r = out[0];

    let faces_tab = collect_planes(&m, r).expect("the result's face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(&m, r, &surf_ix).expect("edge incidence");
    // The cylinder table the tracer would hold, built the way `cylinder_gate` builds it — the gate
    // itself cannot be asked here, because refusing this very operand is its job.
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = crate::planes::world_cylinder_def(&m, surf).expect("a world cylinder");
            let nacre_geom::Surface::Cylinder(cache) = m.surface(surf) else {
                unreachable!("a cylinder class names a cylinder")
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                cache: *cache,
            }
        })
        .collect();

    let jd = crate::planes::test_judge(&planes);
    let (mut curved_walls, mut branch_names) = (0usize, 0usize);
    let (mut ups, mut sides): (Vec<bool>, Vec<i8>) = (Vec::new(), Vec::new());
    for &fh in &m.shells.get(m.solids.get(r).outer).faces {
        let fp = surf_ix[&fh];
        if matches!(plane_ix[fp], crate::planes::ClassIx::Cyl(_)) {
            continue; // the tracer skips a lateral face; its loops are the rims
        }
        let ring = combinatorics::face_vertex_triples(&m, fh, fp, &inc, &jd, &plane_ix, &cyls)
            .unwrap_or_else(|e| panic!("face {:?}: {e:?}", fh.index()));
        let Some(nr) = ring.poly() else { continue };
        let hes = &m.faces.get(fh).outer.half_edges;
        let n = hes.len();
        assert_eq!(
            nr.triples.len(),
            n,
            "face {:?}: one name per corner",
            fh.index()
        );
        let ends = |k: usize| m.edges.get(hes[k].edge).vertices;
        // The corner the walk stands on when it starts edge `k` — the vertex edge `k-1` and edge
        // `k` share. Read off the ring rather than off either edge's stored pair, because a
        // ruling's pair carries no order (`edge_for` keys it unordered) and this is the walk.
        let corner = |k: usize| {
            let (a, b) = (ends((k + n - 1) % n), ends(k));
            *a.iter()
                .find(|x| b.contains(x))
                .expect("consecutive edges share a corner")
        };
        for (i, &node) in nr.triples.iter().enumerate() {
            // The carrier and the corner are independent facts, and both are asserted: a curved
            // wall must be an arc or a ruling of the cylinder its far face is on, never a plane.
            let straight = matches!(m.edge_curve(hes[i].edge), nacre_geom::Curve::Line(_));
            match nr.walls[i] {
                crate::boolean::Wall::Plane(_) => {}
                crate::boolean::Wall::Arc { .. } => {
                    assert!(
                        !straight,
                        "face {:?} edge {i}: an arc carrier on a straight curve",
                        fh.index()
                    );
                    curved_walls += 1;
                }
                // ★★★★★ **The direction bits, against the walk and against the geometry.** The
                // carrier alone says *which surface*; `up` and `side` say *which way* and *which
                // of the two rulings*, and a consumer that reads them wrong builds a face wound
                // backwards or seated on the far side of the cylinder. Both are measured off
                // realized coordinates on purpose — the code derives them without any (`side`
                // exactly through `quad::plane_side`, `up` from the cutting planes' axial
                // parameters), so the coordinate is a genuinely second road to the same bit.
                crate::boolean::Wall::Ruling { cyl: k, side, up } => {
                    assert!(
                        straight,
                        "face {:?} edge {i}: a ruling carrier on a curved curve",
                        fh.index()
                    );
                    let crate::planes::ClassIx::Plane(near) = plane_ix[fp] else {
                        unreachable!("a lateral face was skipped above")
                    };
                    let axis = cyls[k].cache.axis();
                    let axial = |p: Point3| (p - axis.origin()).dot(axis.direction());
                    assert_eq!(
                        up,
                        axial(m.vertex_point(corner((i + 1) % n)))
                            > axial(m.vertex_point(corner(i))),
                        "face {:?} edge {i}: `up` disagrees with the walk",
                        fh.index()
                    );
                    // `side` is the sign of `(x − o) · (m × n̂)`. ★★ **`n̂` has to be the class's
                    // `world_rat` normal, and no other spelling of the class's plane will do** —
                    // measured, by writing `plane.normal()` here first and watching it come out
                    // opposed (`[0,-1,0]` against `[0,1,0]`, and `frame_sign` is `+1`, so that is
                    // not the reconciliation either). A class has no outward normal to agree on:
                    // `world_rat` is the plane's *name*, whose sign is its own, while `plane` is a
                    // stored surface's. `side` is a **label** telling the two rulings apart, so it
                    // is well defined exactly as long as every road spells `n̂` the one way.
                    // Sharing that input leaves the two roads independent where it counts — the
                    // code decides in exact `quad::plane_side` on the branch meet, this in `f64`
                    // on the realized vertex.
                    let wr = combinatorics::class_coeffs_rat(&jd, near).expect("a named class");
                    let n_hat =
                        Vector3::from_array([wr[0].to_f64(), wr[1].to_f64(), wr[2].to_f64()]);
                    let cross = axis.direction().cross(n_hat);
                    let s = (m.vertex_point(corner(i)) - axis.origin()).dot(cross);
                    assert_eq!(
                        side,
                        if s > 0.0 { 1 } else { -1 },
                        "face {:?} edge {i}: `side` disagrees with the geometry ({s})",
                        fh.index()
                    );
                    ups.push(up);
                    sides.push(side);
                    curved_walls += 1;
                }
            }
            let Some((_, cyl, _)) = combinatorics::branch_name(node) else {
                continue;
            };
            branch_names += 1;
            // ★★★★★ **Two roads, one point.** The name is this boolean's classes; the coordinate is
            // what the *previous* boolean realized from *its* classes. Realizing the one and
            // measuring it against the other is what says the restatement — the pair's order and
            // each normal's sign — came out right: get either wrong and the name designates the
            // **other root**, a visibly different point on the far ruling.
            let got = combinatorics::branch_point(&jd, cyl, &cyls[cyl].def, node)
                .expect("the name realizes");
            let want = m.vertex_point(corner(i)).as_array();
            let d: f64 = (0..3)
                .map(|k| (got[k] - want[k]).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(
                d < 1e-9,
                "face {:?} corner {i}: the name realizes at {got:?}, the vertex is at {want:?}",
                fh.index()
            );
        }
    }
    // Four faces run along the boss — the two plate caps it bit an arc out of, and the two halves
    // its rulings split the wall into — and each contributes two curved edges and two branch
    // corners. The boss's own caps are full circles (the one curved loop this road already spoke)
    // and the three untouched walls are plain.
    assert_eq!(curved_walls, 4, "edges riding the boss");
    assert_eq!(branch_names, 8, "corners named as branch points");
    (ups, sides)
}

/// ★★★★★ **A cylinder-pinned end is ordered, not refused — and the sense is the geometry's.**
///
/// `edge_dir` is the one place a direction is made, and its cylinder arm used to say `RingNaming`
/// by name: a branch point has no third plane, and the integer predicate wants one. Measured
/// before this landed — every one of the 40 pins below came back refused. Now they order through
/// the `a + b√c` tower, and this fixes what the answer must be.
///
/// ★★★ **The oracle is a second road, and only the half that matters is second.** The direction
/// `n_p × n_q` is read from the same raw coefficients the code reads — deliberately, the way the
/// ruling lock next door shares `world_rat`: a *label*'s reference frame has to be one spelling or
/// the two roads are not comparing the same thing. What is independent is the part under test —
/// the **order of the two points on that line**: the code decides it exactly (a rational meet
/// against a branch root, or two roots against each other), the oracle realizes both points in
/// `f64` and subtracts. Measured `|t|` from 1.5 to 456, so nothing here is decided in the noise.
///
/// ★★ **The population is one-sided and that is a fact, not a blind instrument.** All 40 order
/// `-1`: a face's outer ring travels counter-clockwise about that face's own outward normal, which
/// is the direction `order_along` sorts by, so agreement is structural — two boss positions, one
/// of them with its axis reversed, and three `Cut` notches do not move it. What moves is **swapping the two ends**,
/// and the oracle swaps with it, so both signs are measured against something that could disagree.
#[test]
fn a_cylinder_pinned_end_orders_through_the_tower() {
    let mut pins = 0usize;
    for (at, dir, kind) in [
        ([2.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Fuse),
        ([12.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Fuse),
        ([2.0, 0.0, 3.0], [0.0, 0.0, -1.0], BoolKind::Fuse),
        ([2.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Cut),
        ([12.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Cut),
        ([6.0, 4.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Cut),
    ] {
        pins += pinned_ends_ordered(at, dir, kind);
    }
    // ★ The count is the lock on the *population*: let the fixtures stop producing branch-pinned
    // ends and every assertion below would pass vacuously. 40 while a ring split at a seam joint
    // was declined; **48** since E1 ③ names it — the `[6, 4, −1]` boss sits on the +y wall, so
    // its bite on the plate's caps wraps the seam (θ = 0 is at −y), and those rings' ends join.
    assert_eq!(pins, 48, "cylinder-pinned ring ends exercised");
}

fn pinned_ends_ordered(at: [f64; 3], dir: [f64; 3], kind: BoolKind) -> usize {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(Point3::from_array(at), Vector3::from_array(dir), 0.5, 4.0);
    m.rebuild_adjacency();
    let out = boolean(&mut m, kind, plate, boss).expect("the boss builds");
    m.rebuild_adjacency();
    let r = out[0];
    let faces_tab = collect_planes(&m, r).expect("the result's face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(&m, r, &surf_ix).expect("edge incidence");
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = crate::planes::world_cylinder_def(&m, surf).expect("a world cylinder");
            let nacre_geom::Surface::Cylinder(cache) = m.surface(surf) else {
                unreachable!()
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                cache: *cache,
            }
        })
        .collect();
    let jd = crate::planes::test_judge(&planes);
    let mut pins = 0usize;
    for &fh in &m.shells.get(m.solids.get(r).outer).faces {
        let fp = surf_ix[&fh];
        let crate::planes::ClassIx::Plane(p) = plane_ix[fp] else {
            continue;
        };
        let Ok(ring) = combinatorics::face_vertex_triples(&m, fh, fp, &inc, &jd, &plane_ix, &cyls)
        else {
            continue;
        };
        let Some(nr) = ring.poly() else { continue };
        let n = nr.triples.len();
        for i in 0..n {
            let crate::boolean::Wall::Plane(q) = nr.walls[i] else {
                continue;
            };
            let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
            let pin = |x: combinatorics::NodeId| match combinatorics::three_plane_name(x) {
                Some(t) => t
                    .iter()
                    .copied()
                    .find(|&c| c != p && c != q)
                    .map(combinatorics::EndPin::Class),
                None => Some(combinatorics::EndPin::Cylinder),
            };
            let (Some(pa), Some(pb)) = (pin(a), pin(b)) else {
                continue;
            };
            if matches!(pa, combinatorics::EndPin::Class(_))
                && matches!(pb, combinatorics::EndPin::Class(_))
            {
                continue;
            }
            pins += 1;
            // The f64 road: realize both points and dot the difference with `n_p × n_q` taken
            // from the raw coefficients — the same direction `plane_pair_dir_sign` reads.
            let xyz = |x: combinatorics::NodeId| -> [f64; 3] {
                match combinatorics::branch_name(x) {
                    Some((_, cyl, _)) => {
                        combinatorics::branch_point(&jd, cyl, &cyls[cyl].def, x).expect("realizes")
                    }
                    None => {
                        let c = combinatorics::node_coords_rat(&jd, x).expect("coords");
                        [c[0].to_f64(), c[1].to_f64(), c[2].to_f64()]
                    }
                }
            };
            let nv = |c: usize| {
                let k = jd.planes[c].plane.coefficients();
                Vector3::from_array([k[0], k[1], k[2]])
            };
            let d = nv(p).cross(nv(q));
            let (xa, xb) = (xyz(a), xyz(b));
            let t: f64 = (0..3).map(|k| (xa[k] - xb[k]) * d.as_array()[k]).sum();
            assert!(t.abs() > 1e-6, "the oracle decides in the noise: {t}");
            let want = if t > 0.0 { 1i8 } else { -1 };
            for (x, y, w) in [((a, pa), (b, pb), want), ((b, pb), (a, pa), -want)] {
                assert_eq!(
                    combinatorics::order_pinned(&jd, &cyls, p, q, x, y),
                    Some(w),
                    "p={p} q={q} {x:?} vs {y:?}"
                );
                // ★ `edge_dir` inverts the order to get the travel sense, and it is the carrier
                // it pairs with — not a bare `i8` — that the walk consumes.
                let dir = combinatorics::edge_dir(&jd, &cyls, p, q, x, y)
                    .unwrap_or_else(|e| panic!("p={p} q={q}: {e:?}"));
                assert!(
                    matches!(dir, combinatorics::EdgeDir::Line { carrier, sense }
                             if carrier == q && sense == -w),
                    "p={p} q={q}: {dir:?}"
                );
            }
        }
    }
    pins
}

/// ★★★ **A circle can be an interior boundary, and then the merge erases it.**
///
/// A tool that only *touches* removes nothing, and the planar engine has always said so plainly:
/// a plate cut by a box resting on its top face comes back as the plate — same volume, **six
/// faces**, no imprint left behind (measured, both for face contact and for a box straddling an
/// edge). That measurement is the oracle here; there is no sentence in the design documents that
/// decides it.
///
/// The cylinder twin used to refuse. Its result was already right — exact volume, `validate`
/// clean — but the top plane came back as **two** faces: the plate's top with a circular hole, and
/// the disk filling it. Nothing joined them, because a disk's boundary is a `Bound::Circle` with
/// no nodes at all, so the merge's edge table could not see it; the corners left on that circle
/// went on naming the tool's cylinder after every face of it was gone, and the assembly refused
/// the whole boolean (`VertexNamesAbsentSurface`).
///
/// What connects them is that the other face holds **the same cylinder class** as a hole — and one
/// plane class with one cylinder class names one circle, so that coincidence *is* adjacency.
#[test]
fn a_disk_merges_into_the_face_it_lies_in() {
    let plate_and_boss = |m: &mut Model, base: [f64; 3], h: f64| {
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            h,
        );
        m.rebuild_adjacency();
        (a, b)
    };
    let faces_of = |m: &Model, s: Handle<Solid>| -> usize {
        let sol = m.solids.get(s).clone();
        std::iter::once(&sol.outer)
            .chain(sol.cavities.iter())
            .map(|&sh| m.shells.get(sh).faces.len())
            .sum()
    };
    let run = |base: [f64; 3], h: f64, kind: BoolKind| -> (f64, usize) {
        let mut m = Model::new();
        let (a, b) = plate_and_boss(&mut m, base, h);
        let out = boolean(&mut m, kind, a, b).expect("the boolean builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // ★ The merged face has to survive tessellation too: `validate` reads the topology store,
        // and a boundary this pass rewrote is exactly the kind a mesher can drop a piece of.
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight"
        );
        (
            nacre_props::mass_props(&m, out[0]).unwrap().volume,
            faces_of(&m, out[0]),
        )
    };

    // ① A boss standing on the top face, cut: the plate, untouched — the planar twin's answer.
    assert_eq!(run([2.0, 2.0, 2.0], 1.0, BoolKind::Cut), (32.0, 6));
    // ② A boss sunk to the top face and fused: it adds nothing, so likewise.
    assert_eq!(run([2.0, 2.0, 1.0], 1.0, BoolKind::Fuse), (32.0, 6));
    // ③ A boss through the plate whose cap is flush with the bottom, fused. The bottom is **one**
    //    face (plate bottom + the boss's cap): 6 plate faces, of which the bottom absorbed the
    //    disk, plus the lateral and the boss's top cap = 8. This is the seam the user could see.
    let (v, f) = run([2.0, 2.0, 0.0], 3.0, BoolKind::Fuse);
    assert_eq!(f, 8, "the bottom is one face");
    // The boss spans z ∈ [0, 3] and the plate z ∈ [0, 2], so only its last millimetre of height
    // adds material.
    assert!(
        (v - (32.0 + std::f64::consts::PI * 0.25)).abs() < 1e-12,
        "{v}"
    );

    // ★ Negative controls — a circle that separates *material from void* must survive.
    // ④ A through bore: its circle is shared with the cylinder, which is not in the plane's group,
    //    so nothing joins and the hole stays a hole.
    let (v, f) = run([2.0, 2.0, -1.0], 4.0, BoolKind::Cut);
    assert_eq!(f, 7, "plate faces with two mouths, plus the bore's lateral");
    assert!(
        (v - (32.0 - std::f64::consts::PI * 0.25 * 2.0)).abs() < 1e-12,
        "{v}"
    );
    // ⑤ A boss standing on the top face, **fused**: the contact disk is interior and never was a
    //    face, so this pass has nothing to do and the answer must not move.
    let (v, f) = run([2.0, 2.0, 2.0], 1.0, BoolKind::Fuse);
    assert_eq!(f, 8, "annulus + 5 plate faces + lateral + cap");
    assert!(
        (v - (32.0 + std::f64::consts::PI * 0.25)).abs() < 1e-12,
        "{v}"
    );
}

/// ★ **A wall whose plane passes near a hole, on a body that moved** — the face-level clearance
/// test's own frame question. The infinite plane `y = 47.8` clears the bore at `(17.1, 48.6)`
/// by 0.8 with `r = 2.12`, so the cheap test fails and each face on that class has to answer for
/// itself; the pocket wall that carries it sits at `x ∈ [60.3, 73.2]`, nowhere near the bore, and
/// says so — but only if its corners are read in the same frame as the axis. The user's 12-up
/// array is this shape, and it stopped here after the class descriptions were carried out.
#[test]
fn a_moved_face_answers_the_clearance_test() {
    use nacre_scalar::Rat;
    let cell = |m: &mut Model| {
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([86.0, 86.0, 71.5]),
        );
        let pocket = m.add_cuboid(
            Point3::from_array([60.3, 47.8, 0.0]),
            Point3::from_array([73.2, 64.5, 68.5]),
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
        m.rebuild_adjacency();
        let bore = m.add_cylinder(
            Point3::from_array([17.1, 48.6, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.12,
            80.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    };
    let mut m = Model::new();
    let a = cell(&mut m);
    let b = cell(&mut m);
    let b = transform(
        &mut m,
        b,
        &nacre_scalar::Isometry::translation([
            Rat::try_from_f64(100.3).unwrap(),
            Rat::from_int(0),
            Rat::from_int(0),
        ]),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(carries_motion(&m, b), "the move records a chain");
    // The cells stand clear of each other (a *touching* pair is the contact family, another
    // cell), so this is an ordinary two-body result — what it pins is that the gate **decided**
    // at all, rather than declining because a corner was stated in another frame.
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved cell is judged");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let pi = std::f64::consts::PI;
    let one = 86.0 * 86.0 * 71.5 - 12.9 * 16.7 * 68.5 - pi * 2.12 * 2.12 * 71.5;
    let v: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((v - 2.0 * one).abs() <= 1e-9 * one, "{v} vs {}", 2.0 * one);
}

/// ★ The negative control: a **rotated** cylinder has no exact world description, and the gate
/// says so by its own name rather than measuring across two frames. (A 90°-family turn keeps a
/// datum exact and records nothing, so the angle here is one that does record.)
#[test]
fn a_rotated_cylinder_is_still_undecided() {
    use nacre_scalar::{Angle, Rat, Rotation};
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = m.add_cylinder(
        Point3::from_array([18.0, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.1,
        30.0,
    );
    m.rebuild_adjacency();
    let tool = transform(
        &mut m,
        tool,
        &nacre_scalar::Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(20), Rat::from_int(20), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(31)).unwrap(),
        }),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(carries_motion(&m, tool), "the turn records a chain");
    let live = m.live_solids.clone();
    let err = boolean(&mut m, BoolKind::Cut, plate, tool).expect_err("no world description");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                ..
            }
        ),
        "{err:?}"
    );
    assert_eq!(m.live_solids, live, "the live set survives the refusal");
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
    assert!(
        nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok(),
        "moved solid tessellates"
    );
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
    assert!(
        nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok(),
        "rotated solid tessellates"
    );
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
    assert!(nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok());
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
                    let nacre_topo::VertexDef::ThreePlane(_) = m.vertices.get(vh).def else {
                        continue;
                    };
                    // ★ The kernel says which corners it can place in one frame; comparing the
                    // three carriers' motions here would be `vertex_meet`'s rule written twice,
                    // and since the invariant-plane restatement that copy calls a turned solid's
                    // corner straddling when the chain-fixes licence places it.
                    if m.vertex_meet(vh).is_none() {
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
        world_rat: None,
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
fn face(plane_idx: usize, nodes: Vec<NodeId>, inner: Vec<Vec<NodeId>>) -> LocalFace {
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
    let v = |x: usize, y: usize| NodeId::three_planes([0, x, y]); // class, x-plane, y-plane
    let (c00, c10, c20, c30) = (v(1, 5), v(2, 5), v(3, 5), v(4, 5));
    let (c01, c11, c21, c31) = (v(1, 6), v(2, 6), v(3, 6), v(4, 6));
    let faces = vec![
        face(0, vec![c00, c10, c11, c01], vec![]),
        face(0, vec![c10, c20, c21, c11], vec![]),
        face(0, vec![c20, c30, c31, c21], vec![]),
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p), &[]).unwrap();
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
    // `unify_keeps_holed_faces`) that built pass-through faces (a retired second node variant) to
    // exercise a passthrough the arrangement never triggers — every node it emits is a three-plane
    // one. The real invariant is exercised here on
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
    let v = |x: usize, y: usize| NodeId::three_planes([0, x, y]);
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
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p), &[]).unwrap();
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
        NodeId::three_planes([0, 1, 4]),
        NodeId::three_planes([0, 2, 4]),
        NodeId::three_planes([0, 3, 4]),
    );
    let (v010, v110, v210) = (
        NodeId::three_planes([0, 1, 5]),
        NodeId::three_planes([0, 2, 5]),
        NodeId::three_planes([0, 3, 5]),
    );
    let (v101, v201) = (
        NodeId::three_planes([2, 4, 6]),
        NodeId::three_planes([3, 4, 6]),
    );
    let faces = vec![
        face(0, vec![v000, v100, v110, v010], vec![]),
        face(0, vec![v100, v200, v210, v110], vec![]),
        face(4, vec![v200, v100, v101, v201], vec![]), // perpendicular, not coplanar
    ];
    let out = unify_coplanar_faces(faces, &crate::planes::test_judge(&p), &[]).unwrap();
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
            let tris = combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                .unwrap()
                .poly()
                .expect("a poly outer")
                .triples
                .clone();
            let tris = names(&tris);
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
                combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                    .unwrap()
                    .poly()
                    .expect("a poly outer")
                    .triples
                    .clone(),
            ];
            rings.extend(
                combinatorics::hole_rings(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                    .unwrap()
                    .into_iter()
                    .filter_map(|lr| lr.poly().map(|nr| nr.triples.clone())),
            );
            for t in rings.iter().flat_map(|r| names(r)) {
                for &k in &t {
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
            let tris = combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &jd, &plane_ix, &[])
                .unwrap()
                .poly()
                .expect("a poly outer")
                .triples
                .clone();
            let tris = names(&tris);
            for t in &tris {
                // The vertex lies on exactly its three defining planes; each must read 0.
                for &q in t {
                    assert_eq!(
                        combinatorics::side_of(
                            &jd,
                            &[],
                            combinatorics::NodeId::three_planes(*t),
                            q
                        ),
                        Some(0),
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

/// The oblique twin of the seated `Common` in `bands`: the same box, but the cylinder's axis runs
/// down the body diagonal, so none of the box's planes is either ⊥ or ∥ to it. Every crossing is
/// an ellipse — the population M6-3 opens, and the one this door still names.
#[test]
fn common_rejects_an_oblique_cylinder() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let cyl = m.add_cylinder(
        Point3::from_array([1.0, 1.0, 0.0]),
        Vector3::from_array([1.0, 1.0, 1.0]).normalize().unwrap(),
        0.5,
        2.0,
    );
    assert_rejects(
        || boolean_one(&mut m, BoolKind::Common, a, cyl),
        RejectReason::ObliqueCylinderCut,
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

// ── K2: a cylinder in the operation log ───────────────────────────────────────────────────────
//
// `Model::add_cylinder` is a test convenience behind a feature; the road an application takes is
// `Operation::Cylinder`, which is what these gates measure. What the frame buys is that the
// statement never leaves the rationals: the axis is the frame's unit normal, the seam its `+u`.

/// The lateral surface's exact truth in a model holding exactly one cylinder.
fn lone_cylinder_def(m: &Model) -> nacre_topo::CylinderDef {
    let mut found = None;
    for (h, _) in m.faces.iter() {
        let s = m.faces.get(h).surface;
        if let nacre_topo::SurfaceTruth::Cylinder { def, .. } = m.surface_truth(s) {
            found = Some(def.clone());
        }
    }
    found.expect("a cylinder solid has a lateral face")
}

fn cylinder_op(m: &Model, center: [f64; 2], radius: f64, dist: f64) -> Operation {
    Operation::Cylinder {
        frame: SketchFrame::world(m, Axis::Z),
        center,
        radius,
        dist,
    }
}

/// ★★ **The frame keeps the statement rational.** On world XY the axis is exactly `(0,0,1)` and
/// the seam reference exactly `(1,0,0)` — small integers, not decimals lifted back out of a
/// normalization. The centre and radius are the caller's own written decimals, lifted once.
///
/// This is the whole reason the op exists as it does: `add_cylinder` reaches the same shape by
/// normalizing an f64 axis, and the `1/3`-base gate in `nacre-topo` shows what that costs.
#[test]
fn a_cylinder_op_states_its_axis_and_seam_as_integers() {
    let mut m = Model::new();
    let op = cylinder_op(&m, [0.5, 0.25], 0.1, 2.0);
    let OpOutput::Cylinder { faces, .. } = apply(&mut m, &op).expect("a world-XY cylinder builds")
    else {
        panic!("a cylinder op answers with a cylinder");
    };
    let def = lone_cylinder_def(&m);
    let int = |n: i128| nacre_scalar::Rat::from_int(n);
    let rat = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
    assert_eq!(
        def.dir(),
        [int(0), int(0), int(1)],
        "the frame's unit normal"
    );
    assert_eq!(def.ref_dir(), [int(1), int(0), int(0)], "the frame's +u");
    assert_eq!(
        def.origin(),
        [rat(1, 2), rat(1, 4), int(0)],
        "the centre, placed"
    );
    assert_eq!(def.radius(), rat(1, 10));
    assert_eq!(faces.len(), 3, "lateral, bottom cap, top cap");
}

/// A cylinder in a log replays to the same model — handles and coordinates both. The op is only
/// worth having if the log stays reproducible with it in there.
#[test]
fn a_cylinder_op_replays_deterministically() {
    let log = vec![
        extrude_log_op(square(), 1.0),
        cylinder_op(&Model::new(), [0.5, 0.5], 0.2, 1.0),
    ];
    let (m1, m2) = (replay(&log).unwrap(), replay(&log).unwrap());
    let pts = |m: &Model| {
        m.vertices
            .iter()
            .map(|(vh, _)| m.vertex_point(vh).as_array())
            .collect::<Vec<_>>()
    };
    assert_eq!(pts(&m1), pts(&m2));
    assert_eq!(m1.faces.len(), m2.faces.len());
    assert_eq!(m1.surface_count(), m2.surface_count());
}

/// ★★ **Every refusal is a name, and none of them leaves a cell behind.** The topo entry panics
/// on the first two; an application's numbers are input, so here they are values. The store
/// lengths are the other half — a refusal that had already pushed would shift every later log
/// index (`apply`'s own doc says so).
#[test]
fn a_refused_cylinder_op_is_named_and_leaves_nothing_behind() {
    let mut m = Model::new();
    // A cylinder surface to aim a frame at — the one thing a `SketchFrame` can name that is not
    // a plane.
    let fixture = cylinder_op(&m, [0.0, 0.0], 1.0, 1.0);
    let OpOutput::Cylinder { faces, .. } =
        apply(&mut m, &fixture).expect("the fixture cylinder builds")
    else {
        panic!("a cylinder op answers with a cylinder");
    };
    let on_lateral = SketchFrame::canonical(m.faces.get(faces[0]).surface);

    let before = (
        m.surface_count(),
        m.vertices.iter().count(),
        m.edges.iter().count(),
        m.faces.len(),
        m.live_solids.len(),
    );
    let wide = 1e300;
    let cases: [(Operation, OpError); 6] = [
        (
            cylinder_op(&m, [0.0, 0.0], 0.0, 1.0),
            OpError::NonPositiveRadius,
        ),
        (
            cylinder_op(&m, [0.0, 0.0], -1.0, 1.0),
            OpError::NonPositiveRadius,
        ),
        (
            cylinder_op(&m, [0.0, 0.0], 1.0, 0.0),
            OpError::NonPositiveDistance,
        ),
        (
            cylinder_op(&m, [0.0, 0.0], wide, 1.0),
            OpError::RadiusOutsideDecimalWindow,
        ),
        (
            cylinder_op(&m, [wide, 0.0], 1.0, 1.0),
            OpError::CenterOutsideDecimalWindow,
        ),
        (
            Operation::Cylinder {
                frame: on_lateral,
                center: [0.0, 0.0],
                radius: 1.0,
                dist: 1.0,
            },
            OpError::NonPlanarFace,
        ),
    ];
    for (op, want) in cases {
        assert_eq!(apply(&mut m, &op).unwrap_err(), want, "op: {op:?}");
    }
    assert_eq!(
        (
            m.surface_count(),
            m.vertices.iter().count(),
            m.edges.iter().count(),
            m.faces.len(),
            m.live_solids.len()
        ),
        before,
        "a refusal is decided before anything is pushed"
    );
}

/// ★★★ **What the app will do, done through the log alone**: a plate, a drill standing through
/// it, and a `Cut`.
///
/// The drill **overshoots** the plate, the way a "through all" hole is drawn. A drill that stops
/// exactly on the plate's own faces works too (`bands`' seated-cap block measures that family);
/// this fixture is the overshooting gesture a kit step emits.
#[test]
fn a_logged_cylinder_cuts_a_through_hole() {
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    // The log is assembled against a scratch model built the same way, so its handles are the
    // index vocabulary `replay` re-anchors (the `two_extrudes_make_two_solids` precedent).
    let mut scratch = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&scratch, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: a, .. } = apply(&mut scratch, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let below = Operation::DatumPlane {
        def: DatumDef::Stated(
            crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
        ),
    };
    let OpOutput::DatumPlane { frame, .. } = apply(&mut scratch, &below).unwrap() else {
        panic!("a datum answers with a datum");
    };
    let drill = Operation::Cylinder {
        frame,
        center: [2.0, 2.0],
        radius: 0.5,
        dist: 3.0,
    };
    let OpOutput::Cylinder { solid: b, .. } = apply(&mut scratch, &drill).unwrap() else {
        panic!("a cylinder op answers with a cylinder");
    };
    let log = vec![
        extrude,
        below,
        drill,
        Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    ];
    let m = replay(&log).expect("a logged drill");
    assert_eq!(m.live_solids.len(), 1, "one drilled plate");
    let v = nacre_props::mass_props(&m, m.live_solids[0])
        .expect("props")
        .volume;
    let want = 4.0 * 4.0 * 2.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!(
        (v - want).abs() < 1e-9,
        "volume {v} is not the plate minus the bore {want}"
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
}

/// ★★ **The negative control the hole above needs**: the same drill, moved clear of the plate,
/// removes nothing. Without it, "the volume dropped by πr²t" could be read as "any cylinder in
/// the log makes a hole".
#[test]
fn a_cylinder_clear_of_the_plate_removes_nothing() {
    let mut m = Model::new();
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: plate_h, .. } = apply(&mut m, &extrude).expect("the plate")
    else {
        panic!("an extrude answers with an extrude");
    };
    let before = nacre_props::mass_props(&m, plate_h).expect("props").volume;
    let above = datum_frame(
        &mut m,
        crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 3.0])),
    );
    let clear = Operation::Cylinder {
        frame: above,
        center: [2.0, 2.0],
        radius: 0.5,
        dist: 1.0,
    };
    let OpOutput::Cylinder { solid: drill, .. } = apply(&mut m, &clear).expect("a clear cylinder")
    else {
        panic!("a cylinder op answers with a cylinder");
    };
    let solids = boolean(&mut m, BoolKind::Cut, plate_h, drill).expect("a disjoint cut answers");
    assert_eq!(solids.len(), 1, "cutting away nothing leaves one solid");
    let after = nacre_props::mass_props(&m, solids[0])
        .expect("props")
        .volume;
    assert!(
        (after - before).abs() < 1e-9,
        "a drill clear of the plate removed {} of material",
        before - after
    );
}

/// ★★ **A face's frame faces *out* of its solid, so a cylinder on it stands on top rather than
/// drilling in.** Asked as a placement question, which is where it is decided — the seam
/// vertices sit at the plate's top face and one unit above it, not below.
///
/// This is the honest statement of what the log expresses today: `flip` is *measured* by the
/// consuming operation (`measured_frame`), never stated by a caller, so drilling into a face is a
/// composite that measures `toward` inward — `PocketOnFace`'s sibling, and not yet written.
#[test]
fn a_cylinder_on_a_face_frame_stands_outward() {
    let mut m = Model::new();
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude).expect("the plate builds") else {
        panic!("an extrude answers with an extrude");
    };
    // faces[1] is the top cap; its frame is measured against that face's outward normal.
    let top = face_sketch_frame(&m, faces[1]).expect("a face has a frame");
    let before = m.vertices.iter().count();
    let boss = Operation::Cylinder {
        frame: top,
        center: [2.0, 2.0],
        radius: 0.5,
        dist: 1.0,
    };
    apply(&mut m, &boss).expect("a cylinder on a face frame builds");
    let seam_z: Vec<f64> = m
        .vertices
        .iter()
        .skip(before)
        .filter(|(_, v)| matches!(v.def, VertexDef::OnSeam(_)))
        .map(|(h, _)| m.vertex_point(h).as_array()[2])
        .collect();
    assert_eq!(seam_z.len(), 2, "one seam vertex per rim");
    assert!(
        seam_z.iter().all(|z| *z >= 2.0 - 1e-12),
        "the cylinder stands outward from the face, not into the plate: {seam_z:?}"
    );
    assert!(
        seam_z.iter().any(|z| (*z - 3.0).abs() < 1e-12),
        "and reaches its full height above it: {seam_z:?}"
    );
}

/// ★★ **A tilted frame builds a cylinder and cannot yet cut with it — and both halves are the
/// point.** The frame has no exact rational basis, so the statement is written in the plane's own
/// frame and a motion node carries it out (the road a tilted prism already takes). The gate then
/// declines a *moved* cylinder by name, because its truth is stated before its motion.
///
/// If the first half ever fails, the op stopped taking the frame-node road; if the second starts
/// succeeding, M6-3 landed and this test should be re-read, not deleted.
#[test]
fn a_cylinder_on_a_tilted_frame_is_built_and_honestly_declined() {
    let mut m = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: square(),
        dist: 1.0,
    };
    let OpOutput::Extrude { solid: cube_h, .. } = apply(&mut m, &extrude).expect("the cube builds")
    else {
        panic!("an extrude answers with an extrude");
    };
    let tilted = datum_frame(
        &mut m,
        crate::SketchPlane::from_origin_normal(
            Point3::from_array([0.5, 0.5, -1.0]),
            nacre_math::Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
        )
        .expect("a tilted plane"),
    );
    let op = Operation::Cylinder {
        frame: tilted,
        center: [0.0, 0.0],
        radius: 0.2,
        dist: 4.0,
    };
    let OpOutput::Cylinder {
        solid: drill,
        faces,
    } = apply(&mut m, &op).expect("a tilted cylinder")
    else {
        panic!("a cylinder op answers with a cylinder");
    };
    let lateral = m.faces.get(faces[0]).surface;
    match m.surface_truth(lateral) {
        nacre_topo::SurfaceTruth::Cylinder { motion, .. } => assert!(
            motion.is_some(),
            "a tilted frame states its cylinder inside a motion node"
        ),
        _ => panic!("the first face is the lateral"),
    }
    // ★ The caches are realized from the same *local* rationals the truth is stated in — the
    // prism road's deal (`SweptRat::motion`). `validate` reads truth against cache, so a
    // world/local mix-up here would surface as `CylinderTruthCacheMismatch` rather than as a
    // wrong model much later.
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert!(
        matches!(
            boolean(&mut m, BoolKind::Cut, cube_h, drill),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                ..
            })
        ),
        "a moved cylinder is declined by name, not silently mis-cut"
    );
}

/// ★★★ **A cylinder's two caps face opposite ways — even when their planes already exist.**
///
/// The bottom cap lies on the very plane the frame names, so it always interns; a plane's
/// canonical name has no direction, so the surface handed back can be the one that faces the
/// other way. Measured before the fix: a cylinder on a plate's top face came out with **both
/// caps facing +Z**, because the b-rep body dropped `push_plane`'s `flipped` bit (the bit
/// `add_cuboid` has always honoured). Not a solid at all, and no earlier fixture looked.
///
/// The three cases are the three ways a frame's plane can arrive: a seeded world plane, one a
/// `DatumPlane` stated, and one an earlier face already put there.
#[test]
fn a_cylinders_caps_face_opposite_ways_however_their_planes_arrived() {
    let outward = |m: &Model, fh: Handle<Face>| -> [f64; 3] {
        let f = m.faces.get(fh);
        let Surface::Plane(p) = m.surface(f.surface) else {
            panic!("a cap is planar")
        };
        let s = f.orientation.sign() as f64;
        p.normal().as_array().map(|c| c * s)
    };
    let check = |m: &Model, faces: [Handle<Face>; 3], what: &str| {
        let (bot, top) = (outward(m, faces[1]), outward(m, faces[2]));
        assert!(
            bot[2] < 0.0,
            "{what}: the bottom cap must face down, got {bot:?}"
        );
        assert!(
            top[2] > 0.0,
            "{what}: the top cap must face up, got {top:?}"
        );
        assert!(
            nacre_validate::validate(m).is_empty(),
            "{what}: {:?}",
            nacre_validate::validate(m)
        );
    };

    // (a) a seeded world plane
    let mut m = Model::new();
    let op = cylinder_op(&m, [0.0, 0.0], 1.0, 2.0);
    let OpOutput::Cylinder { faces, .. } = apply(&mut m, &op).unwrap() else {
        panic!("a cylinder op answers with a cylinder");
    };
    check(&m, faces, "world XY");

    // (b) a plane the log stated as a datum, facing +Z
    let mut m = Model::new();
    let below = datum_frame(
        &mut m,
        crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
    );
    let op = Operation::Cylinder {
        frame: below,
        center: [0.0, 0.0],
        radius: 1.0,
        dist: 2.0,
    };
    let OpOutput::Cylinder { faces, .. } = apply(&mut m, &op).unwrap() else {
        panic!("a cylinder op answers with a cylinder");
    };
    check(&m, faces, "stated datum");

    // (c) the plane of a face that already exists — the case that was wrong
    let mut m = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: square(),
        dist: 1.0,
    };
    let OpOutput::Extrude { faces: plate, .. } = apply(&mut m, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let top = face_sketch_frame(&m, plate[1]).expect("a face has a frame");
    let op = Operation::Cylinder {
        frame: top,
        center: [0.5, 0.5],
        radius: 0.2,
        dist: 1.0,
    };
    let OpOutput::Cylinder { faces, .. } = apply(&mut m, &op).unwrap() else {
        panic!("a cylinder op answers with a cylinder");
    };
    check(&m, faces, "an existing face's plane");
}

/// ★★ **A hole must wind the other way round, and now something checks it.**
///
/// `build_prism` is the single place that decides winding — outer CCW about the
/// sweep, holes the opposite — and until now nothing verified that decision:
/// `shell_signed_volume` takes an *unsigned* area and subtracts, so it assumes the
/// convention rather than reading it.
///
/// The defect is built by hand, and the twin's **orientation** is what gets flipped
/// rather than the loop's direction: reversing a rim half-edge would make
/// `NonOpposedEdge` fire first (the rim is shared with the cylinder's wall) and this
/// check's specificity would be lost. Flipping the flag leaves edge traversal alone,
/// so the two violations that come back are both this rule's — and **`cos`'s sign
/// says which loop spoke**: `−1` the outer polygon, `+1` the inner rim.
#[test]
fn a_holes_winding_is_checked_too() {
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let mut scratch = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&scratch, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: a, .. } = apply(&mut scratch, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let below = Operation::DatumPlane {
        def: DatumDef::Stated(
            crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
        ),
    };
    let OpOutput::DatumPlane { frame, .. } = apply(&mut scratch, &below).unwrap() else {
        panic!("a datum answers with a datum");
    };
    let drill = Operation::Cylinder {
        frame,
        center: [2.0, 2.0],
        radius: 0.5,
        dist: 3.0,
    };
    let OpOutput::Cylinder { solid: b, .. } = apply(&mut scratch, &drill).unwrap() else {
        panic!("a cylinder op answers with a cylinder");
    };
    let mut m = replay(&[
        extrude,
        below,
        drill,
        Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    ])
    .expect("a logged drill");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the drilled plate is clean to begin with: {:?}",
        nacre_validate::validate(&m)
    );

    // A drilled face: an outer polygon with one rim hole.
    let solid = m.live_solids[0];
    let shell = m.solids.get(solid).outer;
    let faces = m.shells.get(shell).faces.clone();
    let victim = *faces
        .iter()
        .find(|&&fh| {
            let f = m.faces.get(fh);
            f.inner.len() == 1 && f.inner[0].half_edges.len() == 1
        })
        .expect("a through hole leaves two drilled faces");
    let twin = {
        let f = m.faces.get(victim).clone();
        m.faces.push(nacre_topo::Face {
            orientation: f.orientation.flipped(),
            ..f
        })
    };
    let sh = m.shells.push(nacre_topo::Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(nacre_topo::Solid {
        outer: sh,
        cavities: vec![],
    });
    m.live_solids.retain(|&s| s == replaced);
    m.rebuild_adjacency();

    let vs = nacre_validate::validate(&m);
    let cosines: Vec<f64> = vs
        .iter()
        .filter_map(|v| match v {
            nacre_validate::Violation::FaceMisoriented { face, cos } if *face == twin => Some(*cos),
            _ => None,
        })
        .collect();
    assert_eq!(vs.len(), 2, "both loops of that face should speak: {vs:?}");
    assert!(
        cosines.iter().any(|c| *c < -0.5),
        "the outer polygon must report a flipped flag: {cosines:?}"
    );
    assert!(
        cosines.iter().any(|c| *c > 0.5),
        "the hole must report that it now winds like an outer loop — this is the \
         assertion that says the inner branch ran: {cosines:?}"
    );
}

/// The polygonal twin of `a_holes_winding_is_checked_too`: a prism with a square
/// hole. Same rule, same convention, and — measured across the suite — a real
/// population (444 polygonal hole loops, all of them already honouring it).
///
/// Worth its own fixture because the two shapes reach the rule by different roads:
/// a rim answers from its circle, a polygon from its Newell sum, and only running
/// both says the shared `loop_winding` serves both.
#[test]
fn a_polygonal_holes_winding_is_checked_too() {
    let profile = Profile2d::with_holes(
        vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)],
        vec![vec![p2(1.0, 1.0), p2(3.0, 1.0), p2(3.0, 3.0), p2(1.0, 3.0)]],
    )
    .expect("a square ring");
    let mut m = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile,
        dist: 1.0,
    };
    let OpOutput::Extrude { solid, .. } = apply(&mut m, &extrude).expect("a ring prism") else {
        panic!("an extrude answers with an extrude");
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the ring prism is clean to begin with: {:?}",
        nacre_validate::validate(&m)
    );

    let shell = m.solids.get(solid).outer;
    let faces = m.shells.get(shell).faces.clone();
    let victim = *faces
        .iter()
        .find(|&&fh| {
            let f = m.faces.get(fh);
            f.inner.len() == 1 && f.inner[0].half_edges.len() >= 3
        })
        .expect("a ring prism has two holed caps");
    let twin = {
        let f = m.faces.get(victim).clone();
        m.faces.push(nacre_topo::Face {
            orientation: f.orientation.flipped(),
            ..f
        })
    };
    let sh = m.shells.push(nacre_topo::Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(nacre_topo::Solid {
        outer: sh,
        cavities: vec![],
    });
    m.live_solids.retain(|&s| s == replaced);
    m.rebuild_adjacency();

    let vs = nacre_validate::validate(&m);
    let cosines: Vec<f64> = vs
        .iter()
        .filter_map(|v| match v {
            nacre_validate::Violation::FaceMisoriented { face, cos } if *face == twin => Some(*cos),
            _ => None,
        })
        .collect();
    assert_eq!(vs.len(), 2, "both loops of that face should speak: {vs:?}");
    assert!(
        cosines.iter().any(|c| *c < -0.5) && cosines.iter().any(|c| *c > 0.5),
        "one report per loop, told apart by sign: {cosines:?}"
    );
}

/// ★★ **A boss's wall faces away from its axis; a bore's faces toward it.**
///
/// Same geometry, opposite sense — the material is inside the wall for one and outside it for
/// the other — and the boolean says so in its own words when it mints a curved face: *"a
/// cylinder's stored normal points away from its axis, and that is the face's outward normal
/// for a boss … `flip` (material outside the wall: a hole) is what reverses it."*
///
/// This asks that sentence back from a different layer: `props::face_normal_at` reads the
/// stored `orientation` and turns the radial direction by it. If the two ever disagreed, a
/// viewer would light a bore inside out and nothing else would notice.
#[test]
fn a_bores_wall_faces_its_axis_and_a_bosss_faces_away() {
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let mut scratch = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&scratch, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: a, .. } = apply(&mut scratch, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let below = Operation::DatumPlane {
        def: DatumDef::Stated(
            crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
        ),
    };
    let OpOutput::DatumPlane { frame, .. } = apply(&mut scratch, &below).unwrap() else {
        panic!("a datum answers with a datum");
    };
    let drill = Operation::Cylinder {
        frame,
        center: [2.0, 2.0],
        radius: 0.5,
        dist: 3.0,
    };
    let OpOutput::Cylinder { solid: b, .. } = apply(&mut scratch, &drill).unwrap() else {
        panic!("a cylinder op answers with a cylinder");
    };

    // The axis of both the bore and the free-standing drill, as a line to measure against.
    let axis_at = |z: f64| Point3::from_array([2.0, 2.0, z]);
    let radial_sense = |m: &Model, solid: Handle<Solid>| -> f64 {
        let shell = m.solids.get(solid).outer;
        let wall = *m
            .shells
            .get(shell)
            .faces
            .iter()
            .find(|&&fh| matches!(m.surface(m.faces.get(fh).surface), Surface::Cylinder(_)))
            .expect("a cylindrical wall");
        // A point on that wall: the seam vertex of one of its rims.
        let he = m.faces.get(wall).outer.half_edges[0];
        let p = m.vertex_point(m.edges.get(he.edge).vertices[0]);
        let n = nacre_props::face_normal_at(m, wall, p).expect("off the axis");
        let out = p - axis_at(p.as_array()[2]);
        n.dot(out)
    };

    // The drill on its own is a boss: its wall faces outward.
    assert!(
        radial_sense(&scratch, b) > 0.0,
        "a free-standing cylinder's wall faces away from its axis"
    );

    // The same cylinder, once it is a hole, faces the other way.
    let m = replay(&[
        extrude,
        below,
        drill,
        Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    ])
    .expect("a logged drill");
    assert!(
        radial_sense(&m, m.live_solids[0]) < 0.0,
        "a bore's wall faces its axis — the material is outside the wall"
    );
}

/// ★★ **A tunnel crossing a bore is not a tunnel touching a bore** — the model a user brought,
/// and the pair rule's real question.
///
/// A plate, a vertical bore through it, and a horizontal tunnel six away with radii summing to
/// four. Nothing touches, and the volume says so exactly. It used to be refused as "two
/// cylinders touch or overlap" because the gate could only measure the distance between
/// *parallel* axes and read its own blind spot as contact.
///
/// ★ The pair rule fires on the **second** boolean — the first result already carries a
/// cylindrical face — so a single cut would not exercise it at all.
///
/// The second half is what keeps the first honest: slide the tunnel onto the bore and the
/// refusal must come back, because two cylinders that really do meet cross in a quartic curve
/// nothing here computes. (Measured with the guard removed: the result kept the *same* volume
/// as the disjoint case — a silent wrong answer, validate clean.)
#[test]
fn a_tunnel_clear_of_a_bore_is_cut_and_one_through_it_is_refused() {
    let build = |x: f64| -> (Model, Handle<Solid>, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([-10.0, -10.0, -1.5]),
            Point3::from_array([10.0, 10.0, 1.5]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([0.0, 0.0, -10.0]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            20.0,
        );
        let tunnel = m.add_cylinder(
            Point3::from_array([x, -25.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
            1.0,
            50.0,
        );
        m.rebuild_adjacency();
        (m, plate, bore, tunnel)
    };

    // Six apart: the two cuts go through, and the volume is the plate minus both voids.
    let (mut m, plate, bore, tunnel) = build(-6.0);
    let drilled = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
    m.rebuild_adjacency();
    let both = boolean(&mut m, BoolKind::Cut, drilled, tunnel).expect("the tunnel cuts too")[0];
    let vol = nacre_props::mass_props(&m, both).expect("props").volume;
    let want = 20.0 * 20.0 * 3.0 - std::f64::consts::PI * 9.0 * 3.0 - std::f64::consts::PI * 20.0;
    assert!(
        (vol - want).abs() < 1e-9,
        "plate minus a bore and a tunnel: {vol} vs {want}"
    );
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );

    // Straight through the bore: still refused, and by the name that now states a fact.
    let (mut m, plate, bore, tunnel) = build(0.0);
    let drilled = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
    m.rebuild_adjacency();
    assert!(
        matches!(
            boolean(&mut m, BoolKind::Cut, drilled, tunnel),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderPairContact,
                ..
            })
        ),
        "a tunnel through the bore is a cylinder pair that meets"
    );
}

/// **The mesh census is live** — the invariant itself is asserted where the fact is made.
///
/// `boolean::tess_census::record` checks every result as it is produced ("undrawable only for a
/// named reason"), because a test reading the vector afterwards sees only the booleans that ran
/// before it. What is left for a test is that the census is **running at all**: a hook that
/// silently stopped recording would take the whole guarantee with it and nothing would go red.
///
/// ★ The count is deliberately not asserted. It is whatever the suite happened to run before this
/// test, and it moves with every fixture added — the population's *size* is pinned where the two
/// known members are, beside their own geometry (`bands::tests`' two tangency fixtures).
#[test]
fn the_mesh_census_is_running() {
    let seen = crate::boolean::tess_census::MESHED
        .lock()
        .expect("the census lock is never held across a panic")
        .clone();
    assert!(!seen.is_empty(), "the census never ran");
    assert!(
        seen.iter().any(|r| r.is_ok()),
        "the census recorded no mesh at all"
    );
}

/// **A ruling carries the label of the cell inside the cylinder** — capability D's vertical answer.
///
/// The chart's horizontal lines have had their answer since the band pass: a ⊥ class's
/// [`crate::arrangement::Label`] holds the material on **both** sides of its plane, which is why
/// the cell reader reads two disk labels and asks them to agree. The vertical lines had none, and
/// this is it — the label of the cell a ruling borders on the axis side of its wall.
///
/// ★★★★★ **The side is derived, and an independent description checks it.**
/// `ruling_interior_is_even` composes three sentences already in the file (`RulingCarrier::side`,
/// the material-on-the-left convention, and `world_rat_sense`'s lift between the rational name and
/// the stored normal) into `side · κ`. The check shares no step with that: **inside the cylinder
/// the lateral's own solid has material and outside it does not**, which is the same content rule
/// `ArcLabels`' doc set its own side by. It is asserted at the record, in `per_class`.
///
/// ★★★★ **And the check is what watches this sign — the volume oracle cannot.** Everywhere else in
/// this ladder a side selector is guarded by the through-boss volume (a global flip passes the
/// relative locks), but that only works for a sign production *reads*. This label is
/// `#[cfg(test)]`, so no volume moves whatever it says. ☑ Flipping `side · κ` turns the check from
/// every agreement into a contradiction — it has eyes on 97% of the population.
///
/// ☑ Measured over the whole binary: **862 ruling pieces, 862 labelled** (`world_rat_sense`
/// declines none) · **838 agree, 0 contradict, 24 blind**. The blind ones are real and expected —
/// a boss surrounded by its own plate has that solid's material on *both* sides of the ruling, so
/// the content cannot tell, and only the derivation speaks there.
#[test]
fn a_ruling_labels_the_cell_inside_the_cylinder() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds");
    let lab = crate::arrangement::ruling_probe::LABELLED
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    let chk = crate::arrangement::ruling_probe::SIDE_CHECK
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    assert!(
        !lab.is_empty(),
        "the probe never ran, so it measured nothing"
    );
    // ★ Universal over every recorded piece, which no interleaving can break — a filtered or
    // parallel run may bring more entries, never different ones.
    assert!(
        lab.iter().all(|&b| b),
        "a ruling piece got no label: {} of {} unlabelled",
        lab.iter().filter(|b| !**b).count(),
        lab.len()
    );
    assert!(
        !chk.is_empty(),
        "the side check never ran beside the labels"
    );
    assert!(
        chk.iter().all(|c| *c != Some(false)),
        "a ruling label contradicted the content check"
    );
    // ★ And the check is not vacuous: it decides for most of the population, not none of it.
    assert!(
        chk.contains(&Some(true)),
        "the content check decided nothing at all"
    );
}

/// **The census's cylinder corpus, in the lib suite** (D3). `tests/census.rs` records these
/// families' digests, but it is an integration test and the lib is compiled without `cfg(test)`
/// there — so the chart's census, which asserts where the facts are made, had never seen them.
/// D2b learned that the hard way: `wal corner-lo` slipped past the shadow comparison. What these
/// lock is only the **outcome** (built, or the refusal's name); the geometry stays the digest's.
/// Fixtures copied from `tests/census.rs` (`rul`, `wal`, `cap`, `ct2`, `trc` families).
///
/// ★ `wal corner-lo` refuses `CylinderGateUndecided` from the chart's emitter: the plate classes'
/// disk cells carry no B material at the (0,0) corner while the caps do (an arrangement label
/// defect, D2b's finding) — the cells' ends disagree and the emitter refuses rather than read
/// either. The census's guard for exactly that shape is `src2_disagree == 0 || emitted.is_err()`.
///
/// ★ `rul flush` (E3-b/c): a straddling boss whose caps sit flush with **both** plate planes. Cut
/// and Common assemble clean and build; the Fuse's two cap chords are each an interior boundary
/// between a half-disk 2-gon and its neighbour that the coplanar merge abstains on (the chord and
/// its arcs collide in the node-pair edge key), and the shipped pair spelled a stated plane-self
/// edge — validate's producer bug. The whole-result check refuses that shape by the cleaning's
/// own name (`CoplanarMerge`, the chord's corner as witness) until the grouping-arm cell teaches
/// the merge to thread mixed rings.
#[test]
fn census_corpus_cylinder_families_build_or_refuse_by_name() {
    fn tr(m: &mut Model, s: Handle<Solid>, v: [f64; 3]) -> Handle<Solid> {
        let s = transform(
            m,
            s,
            &nacre_scalar::Isometry::translation(v.map(|x| Rat::try_from_f64(x).unwrap())),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    }
    fn plate(m: &mut Model, hi: [f64; 3]) -> Handle<Solid> {
        m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array(hi))
    }
    fn cyl(m: &mut Model, base: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
        m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            r,
            h,
        )
    }
    /// A `20 × 20 × 10` cell at `x0`, optionally pocketed, then bored (`ct2`).
    fn ct2cell(m: &mut Model, x0: f64, pocket: bool) -> Handle<Solid> {
        let plate = m.add_cuboid(
            Point3::from_array([x0, 0.0, 0.0]),
            Point3::from_array([x0 + 20.0, 20.0, 10.0]),
        );
        m.rebuild_adjacency();
        let mut out = plate;
        if pocket {
            let p = m.add_cuboid(
                Point3::from_array([x0 + 2.0, 2.0, 4.0]),
                Point3::from_array([x0 + 8.0, 8.0, 10.0]),
            );
            m.rebuild_adjacency();
            out = boolean(m, BoolKind::Cut, out, p).expect("the pocket cuts")[0];
            m.rebuild_adjacency();
        }
        let bore = cyl(m, [x0 + 14.0, 14.0, -1.0], 2.0, 12.0);
        m.rebuild_adjacency();
        out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    }
    type Fx = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);
    /// A family: its name, its two operands, and the outcome per kind (fuse, cut, common).
    type Family = (&'static str, Fx, [Result<usize, RejectReason>; 3]);
    let ok = |n: usize| Ok::<usize, RejectReason>(n);
    let fams: Vec<Family> = vec![
        (
            "rul offmid",
            |m| {
                let a = plate(m, [40.0, 40.0, 20.0]);
                let b = cyl(m, [40.0, 10.0, -10.0], 5.0, 50.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "rul flush",
            |m| {
                let a = plate(m, [40.0, 40.0, 20.0]);
                let b = cyl(m, [40.0, 20.0, 0.0], 5.0, 20.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [Err(RejectReason::CoplanarMerge), ok(1), ok(1)],
        ),
        (
            "wal corner-lo",
            |m| {
                let a = plate(m, [4.0, 4.0, 2.0]);
                let b = cyl(m, [0.0, 0.0, -1.0], 0.5, 4.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [Err(RejectReason::CylinderGateUndecided); 3],
        ),
        (
            "cap sunk",
            |m| {
                let a = plate(m, [4.0, 4.0, 2.0]);
                let b = cyl(m, [2.0, 2.0, 1.0], 0.5, 1.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "trc boss",
            |m| {
                let a = plate(m, [40.0, 40.0, 10.0]);
                let b = cyl(m, [10.0, 20.0, 5.0], 3.0, 8.0);
                m.rebuild_adjacency();
                let b = tr(m, b, [15.3, 0.1, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "trc onaxis",
            |m| {
                let p = plate(m, [40.0, 40.0, 10.0]);
                let bore = cyl(m, [18.0, 20.0, -5.0], 2.1, 30.0);
                m.rebuild_adjacency();
                let holed = boolean(m, BoolKind::Cut, p, bore).expect("the bore cuts")[0];
                m.rebuild_adjacency();
                let twin = cyl(m, [7.3, 20.0, -5.0], 2.1, 30.0);
                m.rebuild_adjacency();
                let twin = tr(m, twin, [10.7, 0.0, 0.0]);
                (holed, twin)
            },
            [Err(RejectReason::CylinderPairContact); 3],
        ),
        (
            "ct2 pocketed",
            |m| {
                let a = ct2cell(m, 0.0, true);
                let b = ct2cell(m, 20.0, true);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "ct2 chain",
            |m| {
                let a = ct2cell(m, 0.0, false);
                let b = ct2cell(m, 20.0, false);
                let ab = boolean(m, BoolKind::Fuse, a, b).expect("the first pair fuses")[0];
                m.rebuild_adjacency();
                let c = ct2cell(m, 40.0, false);
                (ab, c)
            },
            [ok(1), ok(1), ok(0)],
        ),
    ];
    for (name, fx, want) in &fams {
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
        {
            let mut m = Model::new();
            let (a, b) = fx(&mut m);
            match (boolean(&mut m, kind, a, b), want) {
                (Ok(out), Ok(n)) => {
                    assert_eq!(out.len(), *n, "{name} {kind:?}: solids");
                    m.rebuild_adjacency();
                    let issues = nacre_validate::validate(&m);
                    assert!(issues.is_empty(), "{name} {kind:?}: {issues:?}");
                }
                (Err(BoolError::Rejected { reason, .. }), Err(r)) => {
                    assert_eq!(&reason, r, "{name} {kind:?}: the refusal's name");
                }
                (got, _) => panic!("{name} {kind:?}: {got:?}, wanted {want:?}"),
            }
        }
    }
}

/// The `xy` families of the census corpus (rows of bored cells fused in two generations, and the
/// turned negative control), in the lib suite for the same reason as
/// [`census_corpus_cylinder_families_build_or_refuse_by_name`].
#[test]
fn census_corpus_xy_generations_build_or_refuse_by_name() {
    use nacre_scalar::{Angle, Isometry, Rotation};
    fn tr(m: &mut Model, s: Handle<Solid>, v: [f64; 3]) -> Handle<Solid> {
        let s = transform(
            m,
            s,
            &Isometry::translation(v.map(|x| Rat::try_from_f64(x).unwrap())),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    }
    fn cell(m: &mut Model) -> Handle<Solid> {
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([20.0, 20.0, 10.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([6.3, 6.3, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.1,
            12.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    }
    fn row(m: &mut Model, n: usize, step: [f64; 3]) -> Handle<Solid> {
        let mut acc = cell(m);
        for i in 1..n {
            let c = cell(m);
            let c = tr(m, c, step.map(|v| v * i as f64));
            acc = boolean(m, BoolKind::Fuse, acc, c).expect("the row fuses")[0];
            m.rebuild_adjacency();
        }
        acc
    }
    type Fx = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);
    /// A family: its name, its two operands, and the outcome per kind (fuse, cut, common).
    type Family = (&'static str, Fx, [Result<usize, RejectReason>; 3]);
    let ok = |n: usize| Ok::<usize, RejectReason>(n);
    let fams: Vec<Family> = vec![
        (
            "xy grid",
            |m| {
                let a = row(m, 2, [20.0, 0.0, 0.0]);
                let b = row(m, 2, [20.0, 0.0, 0.0]);
                let b = tr(m, b, [0.0, 20.0, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "xy yx",
            |m| {
                let a = row(m, 2, [0.0, 20.0, 0.0]);
                let b = row(m, 2, [0.0, 20.0, 0.0]);
                let b = tr(m, b, [20.0, 0.0, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "xy three",
            |m| {
                let a = row(m, 3, [20.0, 0.0, 0.0]);
                let b = row(m, 3, [20.0, 0.0, 0.0]);
                let b = tr(m, b, [0.0, 20.0, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "xy turned",
            |m| {
                let a = row(m, 2, [20.0, 0.0, 0.0]);
                let b = row(m, 2, [20.0, 0.0, 0.0]);
                let b = transform(
                    m,
                    b,
                    &Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
                    }),
                )
                .unwrap();
                m.rebuild_adjacency();
                (a, b)
            },
            [Err(RejectReason::CylinderGateUndecided); 3],
        ),
    ];
    for (name, fx, want) in &fams {
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
        {
            let mut m = Model::new();
            let (a, b) = fx(&mut m);
            match (boolean(&mut m, kind, a, b), want) {
                (Ok(out), Ok(n)) => {
                    assert_eq!(out.len(), *n, "{name} {kind:?}: solids");
                    m.rebuild_adjacency();
                    let issues = nacre_validate::validate(&m);
                    assert!(issues.is_empty(), "{name} {kind:?}: {issues:?}");
                }
                (Err(BoolError::Rejected { reason, .. }), Err(r)) => {
                    assert_eq!(&reason, r, "{name} {kind:?}: the refusal's name");
                }
                (got, _) => panic!("{name} {kind:?}: {got:?}, wanted {want:?}"),
            }
        }
    }
}

/// ★★★★★ **A tool that reaches the recovered band — the lock the far cube cannot be.** The
/// re-operation census cuts each result with a cube far outside it, so a wrong band would pass
/// it unseen. Here the wall +x boss's result (a band whose notch met the seam and was spliced
/// into the outer walk; E2 reads it back as a hole) is cut by a cuboid whose walls clear the
/// boss's strip and whose ⊥ caps meet the lateral. Volumes derived, not copied: the plate is
/// 4·4·2 = 32; the boss `r = 0.5` on the wall `x = 4` adds its outer half over its height 4
/// (`½·π/4·4 = π/2`) and its inner half outside the plate, below and above (`½·π/4·2 = π/4`) —
/// the inner half inside the plate is plate already — so the fuse is `32 + 3π/4`.
///
/// * A cuboid `[3, 6] × [−1, 5] × [2.5, 4]`: its lower cap crosses the band above the plate,
///   where the circle is whole (the band's rims and holes are read from its cycles, and the
///   notch below carves nothing here), and its upper cap is beyond the band. It removes the boss
///   over `z ∈ [2.5, 3]`, both halves: `π/8`. (A lower cap **on** the plate's top, `z = 2`,
///   refuses `NoClearRay` today — a coplanar-contact containment question, not this road's.)
/// * A cuboid `[3, 6] × [−1, 5] × [0.5, 1.5]`: its caps run **through** the notch, where the
///   band's circle is only the outer half (the Crossing arm) — and through the plate's wall face
///   `x = 4`, whose ring carries the boss's **rulings**; the planar scan's crossing arm has no
///   plane class beside a curved carrier and declines `CurvedRingWall` (E3's wall). Locked as
///   the honest refusal it is today, so the day E3 opens it this test says so.
#[test]
fn a_spliced_band_is_cut_across_its_notch() {
    let pi = std::f64::consts::PI;
    let build = || {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
        m.rebuild_adjacency();
        let v0 = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!(
            (v0 - (32.0 + 3.0 * pi / 4.0)).abs() < 1e-9,
            "the first op's volume: {v0}"
        );
        (m, out[0], v0)
    };
    // (a) caps across the band above the plate and beyond it.
    {
        let (mut m, r0, v0) = build();
        let tool = m.add_cuboid(
            Point3::from_array([3.0, -1.0, 2.5]),
            Point3::from_array([6.0, 5.0, 4.0]),
        );
        m.rebuild_adjacency();
        let r = boolean(&mut m, BoolKind::Cut, r0, tool)
            .unwrap_or_else(|e| panic!("the band is cut above the plate: {e:?}"));
        assert_eq!(r.len(), 1, "one solid");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, r[0]).expect("props").volume;
        assert!(
            (v - (v0 - pi / 8.0)).abs() < 1e-9,
            "{v} vs {}",
            v0 - pi / 8.0
        );
    }
    // (b) caps through the notch: the plate's wall face carries the boss's rulings, and the
    // caps' classes cross it there — the scan names each crossing as a branch node (E3). The
    // tool takes the plate's `x ∈ [3, 4]` slab (4) and the boss's outer half over the slab's
    // height (π/8); its wall `x = 3` clears the boss (1 > r). Beside the volume, the probe: every
    // crossing named in this binary lies on its class, its face plane and the cylinder, on the
    // ruling side the ring carried — the f64 road to the sign `crossing_on_ruling` chose by.
    {
        let (mut m, r0, v0) = build();
        let tool = m.add_cuboid(
            Point3::from_array([3.0, -1.0, 0.5]),
            Point3::from_array([6.0, 5.0, 1.5]),
        );
        m.rebuild_adjacency();
        let r = boolean(&mut m, BoolKind::Cut, r0, tool)
            .unwrap_or_else(|e| panic!("the notch is cut across its rulings: {e:?}"));
        assert_eq!(r.len(), 1, "one solid");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, r[0]).expect("props").volume;
        assert!(
            (v - (v0 - (4.0 + pi / 8.0))).abs() < 1e-9,
            "{v} vs {}",
            v0 - (4.0 + pi / 8.0)
        );
        let hits = crate::arrangement::crossing_probe::HITS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        // Two cap classes × the wall face's two halves, one ruling each.
        assert!(hits.len() >= 4, "the probe saw {} crossings", hits.len());
        for h in &hits {
            assert!(
                h.off.iter().all(|&d| d <= 1e-9),
                "a named crossing sits off its surfaces: {h:?}"
            );
            if !h.arc {
                assert_eq!(
                    h.side_f64, h.side,
                    "a named crossing lies on the other ruling: {h:?}"
                );
            }
        }
    }
}

/// **The boss corpus** — a 4×4×2 plate at the origin and an r = 0.5 boss, `(name, base, height)`,
/// one row per way a boss can sit on a plate: through it, standing on it, flush with its floor, on
/// each wall, at a corner, off the wall's middle, and the half-height variants whose lateral is a
/// chain. Shared by the re-operation census and the crossing census so the two read one corpus.
const BOSS_FAMILIES: [(&str, [f64; 3], f64); 14] = [
    ("through", [2.0, 2.0, -1.0], 4.0),
    ("on top", [2.0, 2.0, 2.0], 1.0),
    ("flush", [2.0, 2.0, 0.0], 3.0),
    ("wall -y", [2.0, 0.0, -1.0], 4.0),
    ("wall +y", [2.0, 4.0, -1.0], 4.0),
    ("wall -x", [0.0, 2.0, -1.0], 4.0),
    ("wall +x", [4.0, 2.0, -1.0], 4.0),
    ("corner", [4.0, 4.0, -1.0], 4.0),
    ("corner-lo", [0.0, 0.0, -1.0], 4.0),
    ("offmid", [4.0, 1.0, -1.0], 5.0),
    ("half wall", [2.0, 0.0, -1.0], 2.0),
    ("half wall, cap below", [2.0, 0.0, 1.0], 2.0),
    ("half +x", [4.0, 2.0, -1.0], 2.0),
    ("half +x, cap below", [4.0, 2.0, 1.0], 2.0),
];

/// The plate and boss of a [`BOSS_FAMILIES`] row, adjacency rebuilt.
fn boss_family(base: [f64; 3], h: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array(base),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        h,
    );
    m.rebuild_adjacency();
    (m, plate, boss)
}

/// A result's lateral faces' boundary cycles by kind — (rims, chains, panels, holes) summed
/// over the laterals; an unnamed cycle panics (the census locks that none is). The tracer's
/// inputs are set up as `pinned_ends_ordered` does.
fn lateral_cycle_census(m: &Model, r: Handle<Solid>) -> [usize; 4] {
    let faces_tab = collect_planes(m, r).expect("the result's face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(m, r, &surf_ix).expect("edge incidence");
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = crate::planes::world_cylinder_def(m, surf).expect("a world cylinder");
            let nacre_geom::Surface::Cylinder(cache) = m.surface(surf) else {
                unreachable!()
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                cache: *cache,
            }
        })
        .collect();
    let jd = crate::planes::test_judge(&planes);
    let mut tally = [0usize; 4];
    for &fh in &m.shells.get(m.solids.get(r).outer).faces {
        let fp = surf_ix[&fh];
        if !matches!(plane_ix[fp], crate::planes::ClassIx::Cyl(_)) {
            continue;
        }
        let cycles = combinatorics::lateral_cycles(m, fh, fp, &inc, &jd, &plane_ix, &cyls)
            .unwrap_or_else(|e| panic!("a lateral's cycles could not be named: {e:?}"));
        for (kind, _) in cycles {
            tally[match kind {
                combinatorics::CycleKind::Rim => 0,
                combinatorics::CycleKind::Chain => 1,
                combinatorics::CycleKind::Panel => 2,
                combinatorics::CycleKind::Hole => 3,
            }] += 1;
        }
    }
    tally
}

/// **The re-operation census** (E1-0): can a boolean's *result* be an operand again? For every
/// family × kind on a 4×4×2 plate with an r = 0.5 boss, the result is cut with a cube that lies
/// far outside it (`[20, 21]³`) — the geometry cannot change, so the second boolean exercises only
/// what an operand needs: every face named as rings, every lateral row stated. What is locked is
/// the **outcome by name** (`Ok` with the volume unchanged, or the refusal's name — a
/// `TraceDeclined` by its `kind`, the face handle being a witness rather than a lock).
///
/// Today's table, measured before anything was built: **Ok 9 · `DegenerateFace` 14 (a two-edge
/// cap: arc + chord, no three loop points to spread) · `CylSpan` 10 (a lateral whose outer loop is
/// not two rims and the seam: a hole spliced into the outer walk, a panel, a chain rim) ·
/// `OuterRing` 5 (a wrap arc split at its `OnSeam` vertex, whose two neighbours are arcs on one
/// cylinder — no third surface names the joint) · first op refused 3 · empty 1**. The rungs that
/// follow change this table one named cause at a time, and each predicts its cells:
///
/// * E1 ②: a loop's corner is its half-edge's start and a circle-bounded face states its
///   triangle from the plane — `DegenerateFace` 14 → 0, every one landing on the lateral row's
///   `CylSpan` (those caps sit beside a panel or chain lateral, which the tracer reaches first);
///   Ok 9 unchanged.
/// * E1 ③: a seam joint is not a corner — the two legs of a wrap arc are one step — `OuterRing`
///   5 → 0, all landing on `CylSpan` (wall +y's notch holds the seam and is spliced into the
///   lateral's outer walk; the half walls' laterals are chains); Ok 9 unchanged. After E1 no
///   row is `OuterRing`/`HoleRing`/`CylFaceHole`/`DegenerateFace`: every ring of a result face
///   has a name, and what remains is the lateral's outer loop (E2).
/// * E2 ②b: the tracer reads a lateral's boundary cycles (`FaceLoops::cycles`) instead of its
///   span — a band's rims and holes, the spliced hole recovered — and declines a panel or a
///   chain rim by name until the chart can hold them: the five Fuse rows whose hole met the seam
///   (wall +y/−x/+x, corner, offmid) become **Ok**; `CylSpan` 29 → 24 (every Cut/Common row is a
///   panel, every half row a chain or a panel); Ok 9 → 14.
/// * E2-2 ①: the circle road reads every cycle (a panel, a chain rim) with an optional outer
///   answer, while the rulings road still states bands only. The reported reject is the lowest
///   class's first declined face, and it moves: a panel's ⊥ classes no longer decline, so the
///   wall/offmid Common rows and every half row surface the **two-edge cap** leaving its chord's
///   class — `CurvedDeparture` 17 (E3-c's population); the Cut rows and corner Common (a three-node
///   cap) stay `CylSpan` on the rulings road: 24 → 7. Ok 14 unchanged.
/// * E2-2 ②: the rulings road sweeps every cycle — the six Cut rows (a panel: the notch's wall)
///   and corner Common become **Ok** with their volumes unchanged; `CylSpan` 7 → 0, Ok 14 → 21.
/// * E3-c: the ring walk names the side an arc departs to, so a half-disk cap is a graze along
///   its chord — every half Fuse and Common row and the wall/offmid Common rows build with their
///   volumes unchanged (13 rows); the half **Cut** rows (a notch with a half-disk ceiling) reach
///   the chart and its emitter refuses a cell it cannot read (`CylinderGateUndecided`, `End::Other`
///   — E2-2's remaining item); Ok 21 → 34, `CurvedDeparture` 17 → 0 (the name is gone).
///
/// ★ The probe beside it counts the loops carrying an `OnSeam` joint at any vertex — outer loops
/// of plane faces and holes of every face (a lateral's outer loop is joined at the seam by
/// construction, so it is not counted). Today: **8 outer, 0 holes** — the population the
/// seam-joint rung reads, and the zero that says a hole never carries one (a hole touching the
/// seam is spliced into the outer walk).
#[test]
fn reop_census_families_reoperate_or_decline_by_name() {
    #[derive(Debug, PartialEq, Eq, Clone, Copy)]
    enum Reop {
        Ok,
        Empty,
        First(RejectReason),
        Rejected(RejectReason),
        Declined(DeclineKind),
    }
    use Reop::*;
    // A lateral's boundary cycles per kind, summed over the result's lateral faces:
    // (rims, chains, panels, holes).
    type Cyc = [usize; 4];
    const RIMS: Cyc = [2, 0, 0, 0];
    const RIMS_HOLE: Cyc = [2, 0, 0, 1];
    const CHAIN: Cyc = [1, 1, 0, 0];
    const PANEL: Cyc = [0, 0, 1, 0];
    const NONE: Cyc = [0; 4];
    // (name, [Fuse, Cut, Common], the cycles each result's laterals carry) — one row per
    // `BOSS_FAMILIES` entry, in its order.
    type Family = (&'static str, [Reop; 3], [Cyc; 3]);
    let families: [Family; 14] = [
        ("through", [Ok, Ok, Ok], [[4, 0, 0, 0], RIMS, RIMS]),
        // A boss standing on the plate: the cut removes nothing and leaves no lateral at all.
        ("on top", [Ok, Ok, Empty], [RIMS, NONE, NONE]),
        ("flush", [Ok, Ok, Ok], [RIMS, RIMS, RIMS]),
        ("wall -y", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("wall +y", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("wall -x", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("wall +x", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("corner", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        (
            "corner-lo",
            [First(RejectReason::CylinderGateUndecided); 3],
            [NONE, NONE, NONE],
        ),
        ("offmid", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        (
            "half wall",
            [Ok, Rejected(RejectReason::CylinderGateUndecided), Ok],
            [CHAIN, PANEL, PANEL],
        ),
        (
            "half wall, cap below",
            [Ok, Rejected(RejectReason::CylinderGateUndecided), Ok],
            [CHAIN, PANEL, PANEL],
        ),
        (
            "half +x",
            [Ok, Rejected(RejectReason::CylinderGateUndecided), Ok],
            [CHAIN, PANEL, PANEL],
        ),
        (
            "half +x, cap below",
            [Ok, Rejected(RejectReason::CylinderGateUndecided), Ok],
            [CHAIN, PANEL, PANEL],
        ),
    ];
    let mut seam_joints = (0usize, 0usize);
    let mut table: Vec<String> = Vec::new();
    let mut mismatches = 0usize;
    let mut cycle_mismatches: Vec<String> = Vec::new();
    for ((name, want, want_cyc), &(fam, base, h)) in families.into_iter().zip(&BOSS_FAMILIES) {
        assert_eq!(name, fam, "the table's rows follow the corpus");
        for ((kind, want), want_cyc) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
            .zip(want_cyc)
        {
            let (mut m, plate, boss) = boss_family(base, h);
            let got = match boolean(&mut m, kind, plate, boss) {
                Err(BoolError::Rejected { reason, .. }) => First(reason),
                Err(e) => panic!("{name} {kind:?}: first op {e:?}"),
                Result::Ok(out) if out.is_empty() => Empty,
                Result::Ok(out) => {
                    m.rebuild_adjacency();
                    // The seam-joint probe, over the result's faces.
                    let solid = m.solids.get(out[0]);
                    for sh in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
                        for &fh in &m.shells.get(sh).faces {
                            let face = m.faces.get(fh);
                            let lateral = matches!(
                                m.surface_truth(face.surface),
                                nacre_topo::SurfaceTruth::Cylinder { .. }
                            );
                            let joint = |lp: &nacre_topo::Loop| {
                                lp.half_edges.len() >= 2
                                    && lp.half_edges.iter().any(|&he| {
                                        matches!(
                                            m.vertices.get(m.he_start(he)).def,
                                            nacre_topo::VertexDef::OnSeam(_)
                                        )
                                    })
                            };
                            if !lateral && joint(&face.outer) {
                                seam_joints.0 += 1;
                            }
                            seam_joints.1 += face.inner.iter().filter(|lp| joint(lp)).count();
                        }
                    }
                    // E2-0: the lateral faces' boundary cycles, named as the tracer will read
                    // them — the outer loop cut at its slits, the spliced hole recovered.
                    let got_cyc = lateral_cycle_census(&m, out[0]);
                    if got_cyc != want_cyc {
                        cycle_mismatches
                            .push(format!("{name} {kind:?}: {got_cyc:?} ← want {want_cyc:?}"));
                    }
                    let v0 = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                    let far =
                        m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
                    m.rebuild_adjacency();
                    match boolean(&mut m, BoolKind::Cut, out[0], far) {
                        Result::Ok(r) => {
                            assert_eq!(r.len(), 1, "{name} {kind:?}: the far cut keeps one solid");
                            m.rebuild_adjacency();
                            let issues = nacre_validate::validate(&m);
                            assert!(issues.is_empty(), "{name} {kind:?}: {issues:?}");
                            let v = nacre_props::mass_props(&m, r[0]).expect("props").volume;
                            assert!(
                                (v - v0).abs() < 1e-9,
                                "{name} {kind:?}: the far cut changed the volume {v0} → {v}"
                            );
                            Ok
                        }
                        Err(BoolError::Rejected {
                            reason: RejectReason::TraceDeclined { kind, .. },
                            ..
                        }) => Declined(kind),
                        Err(BoolError::Rejected { reason, .. }) => Rejected(reason),
                        Err(e) => panic!("{name} {kind:?}: {e:?}"),
                    }
                }
            };
            table.push(format!(
                "{name} {kind:?}: {got:?}{}",
                if got == want { "" } else { "  ← want " }
            ));
            if got != want {
                let last = table.len() - 1;
                table[last].push_str(&format!("{want:?}"));
                mismatches += 1;
            }
        }
    }
    // The whole table at once, so a moved cell is read beside its neighbours.
    assert_eq!(mismatches, 0, "re-operation table:\n{}", table.join("\n"));
    assert!(
        cycle_mismatches.is_empty(),
        "lateral cycles:\n{}",
        cycle_mismatches.join("\n")
    );
    // The distribution the doc states, so a drift in the table is read as a whole.
    let count = |p: fn(&Reop) -> bool| -> usize {
        families
            .iter()
            .flat_map(|(_, w, _)| w.iter())
            .filter(|r| p(r))
            .count()
    };
    assert_eq!(count(|r| *r == Ok), 34, "{table:?}");
    assert_eq!(
        count(|r| *r == Rejected(RejectReason::CylinderGateUndecided)),
        4,
        "a half boss's Cut: the notch's ceiling is a half-disk the chart cannot yet read"
    );
    assert_eq!(
        count(|r| *r == Rejected(RejectReason::DegenerateFace)),
        0,
        "a circle-bounded face is never degenerate"
    );
    assert_eq!(count(|r| *r == Declined(DeclineKind::CylSpan)), 0);
    assert_eq!(
        count(|r| *r == Declined(DeclineKind::OuterRing)),
        0,
        "every ring of a result face has a name"
    );
    assert_eq!(seam_joints, (8, 0), "seam-joint loops (outer, holes)");
}

/// The volume a tool box removes from a family's **first** result, summed from parts: the plate's
/// box ∩ tool, the boss's disk ∩ the tool's footprint times the axial overlap, and the doubly
/// counted disk ∩ plate ∩ tool — `Fuse = plate + boss − both`, `Cut = plate − both`,
/// `Common = both`. A disk clipped by a rectangle is `πr²/2ᵏ` where `k` counts the rectangle's
/// edges through the disk's centre; an edge that cuts the disk anywhere else is not in this corpus
/// and panics rather than approximating.
fn removed_by(kind: BoolKind, base: [f64; 3], h: f64, tool: [[f64; 3]; 2]) -> f64 {
    let seg = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] { [a[0].max(b[0]), a[1].min(b[1])] };
    let len = |a: [f64; 2]| (a[1] - a[0]).max(0.0);
    let plate = [[0.0, 0.0, 0.0], [4.0, 4.0, 2.0]];
    let axis = |b: [[f64; 3]; 2], i: usize| [b[0][i], b[1][i]];
    let box_vol = |a: [[f64; 3]; 2], b: [[f64; 3]; 2]| -> f64 {
        (0..3).map(|i| len(seg(axis(a, i), axis(b, i)))).product()
    };
    let disk_in_rect = |rect: [[f64; 2]; 2]| -> f64 {
        let r = 0.5;
        let mut halvings = 0;
        for (i, c) in [base[0], base[1]].into_iter().enumerate() {
            let (lo, hi) = (rect[i][0], rect[i][1]);
            if hi <= lo {
                return 0.0;
            }
            for (edge, inward) in [(lo, c - lo), (hi, hi - c)] {
                if inward == 0.0 {
                    halvings += 1;
                } else if inward <= -r {
                    return 0.0; // the disk lies wholly beyond this edge
                } else if inward < r {
                    panic!("the tool edge at {edge} cuts the disk off-centre (centre {c})");
                }
            }
        }
        std::f64::consts::PI * r * r / f64::from(1 << halvings)
    };
    let tool_xy = [axis(tool, 0), axis(tool, 1)];
    let plate_xy = [seg(tool_xy[0], [0.0, 4.0]), seg(tool_xy[1], [0.0, 4.0])];
    let boss_z = [base[2], base[2] + h];
    let l_bt = len(seg(boss_z, axis(tool, 2)));
    let l_pbt = len(seg(seg(boss_z, [0.0, 2.0]), axis(tool, 2)));
    let plate_tool = box_vol(plate, tool);
    let boss_tool = l_bt * disk_in_rect(tool_xy);
    let both = l_pbt * disk_in_rect(plate_xy);
    match kind {
        BoolKind::Fuse => plate_tool + boss_tool - both,
        BoolKind::Cut => plate_tool - both,
        BoolKind::Common => both,
    }
}

/// **The crossing census** (E3-0): a boolean's result cut by a tool whose faces actually **cross**
/// the result's rings — where the far cube of the re-operation census could see only whether the
/// rows are stated. Four tools per family × kind: a slab through the plate's middle
/// (`z ∈ [0.5, 1.5]`, its ⊥ caps crossing the wall faces that carry a boss's rulings — E3's wall),
/// a slab above it (`[2.5, 3.5]`, crossing the bands standing on the plate), a slab below
/// (`[−1.5, −0.5]`), and a box whose wall passes **through the boss's axis** (so it crosses the
/// caps' arcs and cuts the lateral along two rulings). What is locked is the outcome by name —
/// `Ok(n)` with the volume equal to the first result's minus [`removed_by`] — and the counts.
///
/// Today's table, measured before anything was built (168 cells): **Ok 36** (every band result ×
/// slab; the mid slab leaves two solids) · **`CurvedRingWall` 11** — the wall-boss Fuse rows × mid
/// slab (6: the ⊥ caps cross the plate's wall face on the boss's rulings) and × through-axis wall
/// (5: the wall crosses the bitten cap's arc) · `CylSpan` 96 (panel and chain laterals, E2-2) ·
/// `NoClearRay` 8 (an interior boss × through-axis wall: the wall halves the cap's circular hole;
/// corner × its coplanar wall) · `BranchVertexUnnamed` 1 (offmid Fuse × top slab) · first op
/// refused 12 · empty 4. The rungs that follow move this table one named cause at a time:
///
/// * E3 ②: the scan names a crossing on a ruling as a branch node, and a cell's nesting reads a
///   mixed ring. The **wall slab** column (added here, 210 cells) is the rung's own population:
///   wall ±x/±y and offmid Fuse × wall slab are **Ok(1)** with their exact volumes. The mid slab
///   crosses the same rulings and then splits the result in two, where the grouping road's
///   plane-walls shim refuses the mixed rings — `CurvedRingWall` 6 → `BranchVertexUnnamed` 5 +
///   `RulingBoundNotYet` 1 (corner: the chart's `End::Other`, E2-2); the through-axis wall's
///   `NoClearRay` 8 → `BranchVertexUnnamed` 6 (the chord's cell has no cylinder for its corners,
///   `coord_key`) + 2 (Common). Ok 49 · `CurvedRingWall` 5 (arc crossings) · `CylSpan` 120 ·
///   `BranchVertexUnnamed` 12 · `RulingBoundNotYet` 2 · `NoClearRay` 2 · first 15 · empty 5.
/// * E2-2 ①: the circle road reads panels and chains; the first reject moves off the ⊥ classes,
///   and the rows whose result has a two-edge cap (wall/offmid Common, every half family) surface
///   `CurvedDeparture` on every tool — 85; `CylSpan` 120 → 35 (the Cut rows and corner Common,
///   still on the rulings road).
/// * E2-2 ②: the rulings road sweeps every cycle, so a Cut result's panel is stated: the six wall
///   Cut rows read like the Fuse rows (top/bottom Ok(1), the **wall slab Ok(1) with its exact
///   volume**, the mid slab on the grouping road, the through-axis wall on the bite's arc) and the
///   corner's coplanar tool reaches the chord's cell; corner Common (a quarter cylinder alone)
///   builds under the top and bottom slabs and finds no clear ray once a slab parts it. `CylSpan`
///   0 · Ok 69 · `BranchVertexUnnamed` 19 · `CurvedRingWall` 10 · `NoClearRay` 4 ·
///   `CylinderGateUndecided` 1 (the chart's walk has an open end at a wall with one ruling).
/// * E3-c: the half-disk caps are traced, so `CurvedDeparture` 85 → 0: the half Fuse rows read
///   like the wall Fuse rows (top/bottom/wall slab Ok(1) with exact volumes), the Common rows and
///   the half Cut rows reach the next walls by name — the coplanar cleaning pass's plane-only ray
///   (`NoClearRay` 4 → 40) and the chart's unreadable cell (`CylinderGateUndecided` 1 → 17); the
///   through-axis wall still crosses the caps' arcs (`CurvedRingWall` 10 → 27, E3-b). Ok 69 → 81.
/// * E3-b: the scan names a crossing on an **arc**, the caps' chords join the overlay (so a
///   through-axis tool wall splits them), and a chord's sense reads its canonical root — the
///   whole through-axis wall column completes its trace and lands on the chart's unreadable cell:
///   `CurvedRingWall` 27 → 0, `CylinderGateUndecided` 17 → 44. The name is gone from the corpus.
#[test]
fn crossing_census_slabs_and_through_axis_walls_by_name() {
    #[derive(Debug, PartialEq, Clone, Copy)]
    enum Cross {
        Ok(usize),
        First,
        Empty,
        Rejected(RejectReason),
        Declined(DeclineKind),
    }
    use Cross::*;
    let mid = [[-1.0, -1.0, 0.5], [5.0, 5.0, 1.5]];
    let top = [[-1.0, -1.0, 2.5], [5.0, 5.0, 3.5]];
    let bottom = [[-1.0, -1.0, -1.5], [5.0, 5.0, -0.5]];
    // The through-axis wall: `y = by` for a boss on an x-wall or in the interior, `x = bx` for one
    // on a y-wall — the wall that is not the plate's own.
    let axis_wall = |b: [f64; 3]| -> [[f64; 3]; 2] {
        if b[1] == 0.0 || b[1] == 4.0 {
            [[b[0], b[1] - 3.0, -2.0], [b[0] + 3.0, b[1] + 3.0, 6.0]]
        } else {
            [[b[0] - 3.0, b[1], -2.0], [b[0] + 3.0, b[1] + 3.0, 6.0]]
        }
    };
    // The wall slab: the mid slab's height, but only past the line one unit inside the wall the
    // boss stands on — its cap classes cross that wall face on the boss's rulings while the plate
    // stays one solid (E3-a's own question, apart from the mid slab's second one: splitting the
    // result). For an interior boss it clears the boss and is a control.
    let wall_slab = |b: [f64; 3]| -> [[f64; 3]; 2] {
        if b[0] == 4.0 {
            [[3.0, -1.0, 0.5], [6.0, 5.0, 1.5]]
        } else if b[0] == 0.0 {
            [[-2.0, -1.0, 0.5], [1.0, 5.0, 1.5]]
        } else if b[1] == 4.0 {
            [[-1.0, 3.0, 0.5], [5.0, 6.0, 1.5]]
        } else if b[1] == 0.0 {
            [[-1.0, -2.0, 0.5], [5.0, 1.0, 1.5]]
        } else {
            [[3.0, -1.0, 0.5], [6.0, 5.0, 1.5]]
        }
    };
    const TOOLS: [&str; 5] = ["mid", "top", "bottom", "axis wall", "wall slab"];
    // Per family, per kind: the five tools' outcomes, in `TOOLS` order.
    use RejectReason::{BranchVertexUnnamed, NoClearRay, RulingBoundNotYet};
    // A band result: the mid slab parts the plate into two solids, the top and bottom slabs cut
    // the standing boss, the wall slab clears the boss. The through-axis wall halves the cap's
    // circular hole: a chord edge with branch ends rides the **plane** (`chord_on_class`), so a
    // cell bounded by the chord alone carries no cylinder for its corners — `coord_key` refuses
    // (`BranchVertexUnnamed`) — and for Common, whose result is the bore alone, the point
    // classification has no clear ray at the hole's diameter.
    const BAND: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)];
    const BAND_COMMON: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Rejected(NoClearRay), Ok(1)];
    // A wall boss fused: the wall slab's caps cross the plate's wall face on the boss's
    // **rulings** — named as branch nodes since E3 — and the cut builds with its exact volume.
    // The mid slab does the same and then splits the result in two, and the grouping road
    // (`Ring::edges`, a legacy plane-walls shim) refuses the mixed rings by name; the through-axis
    // wall crosses the bitten cap's **arc**, which waits on the lateral's ruling sweep.
    const WALL_FUSE: [Cross; 5] = [
        Rejected(BranchVertexUnnamed),
        Ok(1),
        Ok(1),
        Rejected(RejectReason::CylinderGateUndecided),
        Ok(1),
    ];
    // A wall boss cut (its lateral a panel — the notch's wall): the top and bottom slabs miss
    // the result (nothing of the boss stands outside the plate), the wall slab crosses the
    // panel's rulings and builds with its exact volume (E2-2's own population), the mid slab
    // splits the result and the grouping road refuses the mixed rings, and the through-axis wall
    // crosses the bite's arcs.
    const PANEL_CUT: [Cross; 5] = [
        Rejected(BranchVertexUnnamed),
        Ok(1),
        Ok(1),
        Rejected(RejectReason::CylinderGateUndecided),
        Ok(1),
    ];
    // A half boss fused (its lateral a chain, its top cap a half-disk): since E3-c the cap's
    // chord is a graze on the wall class and the result builds — the top and bottom slabs miss
    // the boss, the wall slab crosses the chain's rulings with its exact volume, the mid slab
    // splits the result (grouping road), the through-axis wall crosses the cap's arc.
    const HALF_FUSE: [Cross; 5] = [
        Rejected(BranchVertexUnnamed),
        Ok(1),
        Ok(1),
        Rejected(RejectReason::CylinderGateUndecided),
        Ok(1),
    ];
    // A Common result — the inner half-cylinder alone, two half-disk caps and a flat side — under
    // any slab: the coplanar cleaning pass nests its rings through the plane-only ray road, which
    // has no clear ray from a branch corner (the grouping road's shim, the next wall by name).
    const COMMON: [Cross; 5] = [
        Rejected(NoClearRay),
        Rejected(NoClearRay),
        Rejected(NoClearRay),
        Rejected(RejectReason::CylinderGateUndecided),
        Rejected(NoClearRay),
    ];
    // A half boss's Cut (a notch with a half-disk ceiling): the chart's emitter meets a cell it
    // cannot read (`End::Other`), whatever the tool.
    const HALF_CUT: [Cross; 5] = [
        Rejected(RejectReason::CylinderGateUndecided),
        Rejected(RejectReason::CylinderGateUndecided),
        Rejected(RejectReason::CylinderGateUndecided),
        Rejected(RejectReason::CylinderGateUndecided),
        Rejected(RejectReason::CylinderGateUndecided),
    ];
    let want: [[[Cross; 5]; 3]; 14] = [
        [BAND, BAND, BAND_COMMON], // through
        // on top: the cut leaves the plate alone.
        [BAND, [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)], [Empty; 5]],
        [BAND, BAND, BAND_COMMON],      // flush
        [WALL_FUSE, PANEL_CUT, COMMON], // wall -y
        [WALL_FUSE, PANEL_CUT, COMMON], // wall +y
        [WALL_FUSE, PANEL_CUT, COMMON], // wall -x
        [WALL_FUSE, PANEL_CUT, COMMON], // wall +x
        // corner: either slab cuts the Fuse's circle at four rulings (two walls), and the chart
        // cannot pair a cut end with its rim (`End::Other`, E2-2) — the emitter refuses by name;
        // the through-axis tool's wall `y = 4` is the plate's own wall (coplanar contact), which
        // after E2-2 reaches the chord's cell (`coord_key`) for the Cut and a quarter lateral whose
        // wall has one ruling in the interval for the Common (`CylinderGateUndecided`: the chart's
        // walk has an open end there). The Common is a quarter cylinder alone, so a slab through
        // it leaves two solids and the point classification finds no clear ray.
        [
            [
                Rejected(RulingBoundNotYet),
                Ok(1),
                Ok(1),
                Rejected(BranchVertexUnnamed),
                Rejected(RulingBoundNotYet),
            ],
            [
                Rejected(BranchVertexUnnamed),
                Ok(1),
                Ok(1),
                Rejected(BranchVertexUnnamed),
                Ok(1),
            ],
            [
                Rejected(NoClearRay),
                Ok(1),
                Ok(1),
                Rejected(RejectReason::CylinderGateUndecided),
                Rejected(NoClearRay),
            ],
        ],
        [[First; 5]; 3], // corner-lo: the first op is refused
        // offmid: the top slab's cap at z = 2.5 meets the notch's rulings above the plate — a
        // branch vertex the result cannot name yet (outside this rung).
        [
            [
                Rejected(BranchVertexUnnamed),
                Rejected(BranchVertexUnnamed),
                Ok(1),
                Rejected(RejectReason::CylinderGateUndecided),
                Ok(1),
            ],
            PANEL_CUT,
            COMMON,
        ],
        [HALF_FUSE, HALF_CUT, COMMON], // half wall
        [HALF_FUSE, HALF_CUT, COMMON], // half wall, cap below
        [HALF_FUSE, HALF_CUT, COMMON], // half +x
        [HALF_FUSE, HALF_CUT, COMMON], // half +x, cap below
    ];
    let mut table: Vec<String> = Vec::new();
    let mut mismatches = 0usize;
    let mut tally: Vec<(Cross, usize)> = Vec::new();
    for (&(name, base, h), want) in BOSS_FAMILIES.iter().zip(want) {
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
        {
            let tools = [mid, top, bottom, axis_wall(base), wall_slab(base)];
            for ((tool_name, tool), want) in TOOLS.into_iter().zip(tools).zip(want) {
                let (mut m, plate, boss) = boss_family(base, h);
                let got = match boolean(&mut m, kind, plate, boss) {
                    Err(BoolError::Rejected { .. }) => First,
                    Err(e) => panic!("{name} {kind:?}: first op {e:?}"),
                    Result::Ok(out) if out.is_empty() => Empty,
                    Result::Ok(out) => {
                        m.rebuild_adjacency();
                        let v0 = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                        let t =
                            m.add_cuboid(Point3::from_array(tool[0]), Point3::from_array(tool[1]));
                        m.rebuild_adjacency();
                        match boolean(&mut m, BoolKind::Cut, out[0], t) {
                            Result::Ok(r) => {
                                m.rebuild_adjacency();
                                let issues = nacre_validate::validate(&m);
                                assert!(
                                    issues.is_empty(),
                                    "{name} {kind:?} × {tool_name}: {issues:?}"
                                );
                                let v: f64 = r
                                    .iter()
                                    .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                                    .sum();
                                let expect = v0 - removed_by(kind, base, h, tool);
                                assert!(
                                    (v - expect).abs() < 1e-9,
                                    "{name} {kind:?} × {tool_name}: volume {v} ≠ {v0} − removed = {expect}"
                                );
                                Ok(r.len())
                            }
                            Err(BoolError::Rejected {
                                reason: RejectReason::TraceDeclined { kind, .. },
                                ..
                            }) => Declined(kind),
                            Err(BoolError::Rejected { reason, .. }) => Rejected(reason),
                            Err(e) => panic!("{name} {kind:?} × {tool_name}: {e:?}"),
                        }
                    }
                };
                match tally.iter_mut().find(|(c, _)| *c == got) {
                    Some((_, n)) => *n += 1,
                    None => tally.push((got, 1)),
                }
                let ok = want == got;
                table.push(format!(
                    "{name} {kind:?} × {tool_name}: {got:?}{}",
                    if ok {
                        String::new()
                    } else {
                        format!("  ← want {want:?}")
                    }
                ));
                if !ok {
                    mismatches += 1;
                }
            }
        }
    }
    assert_eq!(
        mismatches,
        0,
        "crossing table:\n{}\n\ntally: {tally:?}",
        table.join("\n")
    );
    // The distribution the doc states, so a moved cell is read as a whole: 210 cells.
    let count = |p: fn(&Cross) -> bool| -> usize {
        want.iter().flatten().flatten().filter(|c| p(c)).count()
    };
    assert_eq!(count(|c| matches!(c, Ok(_))), 86, "{tally:?}");
    assert_eq!(count(|c| *c == Rejected(NoClearRay)), 40);
    assert_eq!(count(|c| *c == Rejected(BranchVertexUnnamed)), 18);
    assert_eq!(count(|c| *c == Rejected(RulingBoundNotYet)), 2);
    assert_eq!(
        count(|c| *c == Rejected(RejectReason::CylinderGateUndecided)),
        44
    );
    assert_eq!(count(|c| *c == First), 15);
    assert_eq!(count(|c| *c == Empty), 5);
}

/// **The ⊥ road runs on a panel and on a chain, and says what the walk-through predicts** (E2-2 ①).
/// Locked through the road's own ledger (`cycle_probe`) because neither face reaches a result yet:
/// the panel is refused one road later (the rulings road still states bands only), the chain's
/// re-operation is refused at its cap's chord (`CurvedDeparture`) — and every class is traced by
/// hand (`trace_every_class`), since the production driver stops at the first decline and, in
/// serial mode, would never reach the classes this reads. Stations are axis parameters from the
/// boss's base (z = −1), so `t = z + 1`.
///
/// Panel (wall +x Cut, cut by the wall slab): at `z = 0.5` the class runs strictly inside the
/// panel's range — outer `Crosses`, the two rulings crossed carve **one** extent (the outer half)
/// away, one span left;
/// at `z = 0` the class is the panel's own lower arc — no outer answer, one graze carved, one
/// span. Chain (half wall Fuse, cut by the far cube): `z = −1` is the whole rim — a graze, nothing
/// carved; `z = 0` is inside — `Crosses`, the plate-bottom arc carves a graze on the inner half,
/// two spans; `z = 1` is the chain's top — no outer answer, the boss-top arc grazes, one span.
#[test]
fn a_cycle_is_carved_on_its_classes() {
    {
        let (mut m, plate, boss) = boss_family([4.0, 2.0, -1.0], 4.0);
        let out = boolean(&mut m, BoolKind::Cut, plate, boss).expect("the notch builds");
        m.rebuild_adjacency();
        let t = m.add_cuboid(
            Point3::from_array([3.0, -1.0, 0.5]),
            Point3::from_array([6.0, 5.0, 1.5]),
        );
        m.rebuild_adjacency();
        crate::arrangement::trace_every_class(&m, out[0], t).expect("the panel's classes trace");
    }
    {
        let (mut m, plate, boss) = boss_family([2.0, 0.0, -1.0], 2.0);
        let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the half boss builds");
        m.rebuild_adjacency();
        let far = m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
        m.rebuild_adjacency();
        crate::arrangement::trace_every_class(&m, out[0], far).expect("the chain's classes trace");
    }
    use crate::arrangement::CylOnClass;
    let hits = crate::arrangement::cycle_probe::HITS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    // `kinds = [rims, chains, panels, holes]`; every recorded entry of the shape at the station
    // must read the same, and at least one must exist. A graze's `body_above` is written in the
    // class's stored frame, so the rim check asks only for the kind of answer.
    let check = |kinds: [usize; 4],
                 t: f64,
                 outer: fn(Option<CylOnClass>) -> bool,
                 carved: usize,
                 spans: usize| {
        let rows: Vec<_> = hits
            .iter()
            .filter(|h| h.kinds == kinds && h.t == t)
            .collect();
        assert!(!rows.is_empty(), "no record for {kinds:?} at t = {t}");
        for h in rows {
            assert!(outer(h.outer), "{kinds:?} at t = {t}: {h:?}");
            assert_eq!(
                (h.carved, h.spans),
                (carved, spans),
                "{kinds:?} at t = {t}: {h:?}"
            );
        }
    };
    let crosses = |o: Option<CylOnClass>| o == Some(CylOnClass::Crosses);
    let grazes = |o: Option<CylOnClass>| matches!(o, Some(CylOnClass::Grazes { .. }));
    let absent = |o: Option<CylOnClass>| o.is_none();
    let panel = [0, 0, 1, 0];
    check(panel, 1.5, crosses, 1, 1);
    check(panel, 1.0, absent, 1, 1);
    let chain = [1, 1, 0, 0];
    check(chain, 0.0, grazes, 0, 1);
    check(chain, 1.0, crosses, 1, 2);
    check(chain, 2.0, absent, 1, 1);
}

/// **The rulings road sweeps a chain rim** (E2-2 ②): on the half wall boss's wall class `y = 0`,
/// each ruling is `Transversal` from the whole rim at `z = −1` up to the plate's bottom `z = 0`
/// (both halves of the boss are there), then a `Graze` along the chain's own ruling edge up to
/// the boss's top `z = 1` — the collinear run whose arcs arrive from the outer half and leave into
/// the inner one toggles the face off above it. A lock of the **rule**, not of a result: every
/// re-operation of a half boss is still refused at its cap's chord (`CurvedDeparture`, E3-c), so
/// the classes are traced by hand and the road's ledger is read for this boss (`origin
/// (2, 0, −1)`, rulings spanning `t ∈ [0, 2]`).
#[test]
fn a_chain_sweeps_its_rulings() {
    let (mut m, plate, boss) = boss_family([2.0, 0.0, -1.0], 2.0);
    let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the half boss builds");
    m.rebuild_adjacency();
    let far = m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
    m.rebuild_adjacency();
    crate::arrangement::trace_every_class(&m, out[0], far).expect("the chain's classes trace");
    let rows: Vec<Vec<crate::arrangement::SegKind>> = crate::arrangement::ruling_probe::CARVED
        .lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        .filter(|c| c.origin == [2.0, 0.0, -1.0] && c.span == [0.0, 2.0])
        .map(|c| c.kinds.clone())
        .collect();
    assert!(rows.len() >= 2, "both rulings sweep: {rows:?}");
    for kinds in &rows {
        assert!(
            matches!(
                kinds[..],
                [
                    crate::arrangement::SegKind::Transversal { .. },
                    crate::arrangement::SegKind::Graze { .. }
                ]
            ),
            "a chain's ruling swept as {kinds:?}"
        );
    }
}

/// **The ring walk cuts a run where an arc departs, and the departure's side flanks the pieces**
/// (E3-c). A hand-built ring on the unit cube's classes: `(0,0,0) → (1,0,0) → (1,1,0) → (1,1,1) →
/// (1,0,1) → (0,0,1)`, read against the class `x = 1` — four nodes on the line between two off it
/// (both on the `x < 1` side). With the edge `(1,1,0) → (1,1,1)` declared a departure to side σ,
/// the line is met in two runs, `[(1,0,0),(1,1,0)]` and `[(1,1,1),(1,0,1)]`: the first is flanked
/// by the off-line node and σ, the second by σ and the off-line node — so both pieces cross
/// exactly when σ is the *other* side, and negating σ swaps the answer. With the edge on the line
/// the four nodes are one run that touches and turns back. The closure supplies σ, so no
/// cylinder is needed: this locks the walk's rule, and `arc_departure_side`'s sign is locked by
/// the re-operation census (a wall boss's plate face is such a run with the arc's own σ).
#[test]
fn the_walk_cuts_a_run_at_a_departure() {
    use crate::combinatorics::{EdgeMeet, Feature, RingWalk, ring_against_plane, side_of};
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
    m.rebuild_adjacency();
    let setup = crate::planes::plane_index_setup(&m, a, b).expect("setup");
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    // A class by the coordinate all three of its triangle's points share.
    let class = |axis: usize, at: f64| {
        (0..setup.geom.len())
            .find(|&c| {
                setup.geom[c]
                    .tri
                    .iter()
                    .all(|p| (p.as_array()[axis] - at).abs() < 1e-12)
            })
            .expect("a class of the unit cube")
    };
    let (x0, x1, y0, y1, z0, z1) = (
        class(0, 0.0),
        class(0, 1.0),
        class(1, 0.0),
        class(1, 1.0),
        class(2, 0.0),
        class(2, 1.0),
    );
    let ring = [
        NodeId::three_planes([x0, y0, z0]),
        NodeId::three_planes([x1, y0, z0]),
        NodeId::three_planes([x1, y1, z0]),
        NodeId::three_planes([x1, y1, z1]),
        NodeId::three_planes([x1, y0, z1]),
        NodeId::three_planes([x0, y0, z1]),
    ];
    let off = side_of(&jd, &[], ring[0], x1).expect("a side");
    assert_ne!(off, 0);
    let runs = |walk: RingWalk| -> Vec<(usize, usize, bool, i8)> {
        let RingWalk::Met(f) = walk else {
            panic!("the ring meets the line")
        };
        f.into_iter()
            .map(|f| match f {
                Feature::Run {
                    first,
                    len,
                    flanks_differ,
                    flank,
                } => (first, len, flanks_differ, flank),
                Feature::Crossing { .. } => panic!("no edge crosses x = 1 strictly"),
            })
            .collect()
    };
    for sigma in [off, -off] {
        let got = runs(ring_against_plane(&jd, &[], &ring, x1, |i| {
            Some(if i == 2 {
                EdgeMeet::Departs(sigma)
            } else {
                EdgeMeet::On
            })
        }));
        let crosses = sigma != off;
        assert_eq!(
            got,
            vec![(1, 2, crosses, off), (3, 2, crosses, sigma)],
            "σ = {sigma}, off-line side {off}"
        );
    }
    let got = runs(ring_against_plane(&jd, &[], &ring, x1, |_| {
        Some(EdgeMeet::On)
    }));
    assert_eq!(got, vec![(1, 4, false, off)]);
}
