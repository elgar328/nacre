//! Wall faces: the gate asks the boundary, not the infinite plane.

use super::*;

//
// ★★ The uniform-slab theorem's premise is about the other operand's **boundary**. Testing the
// wall's infinite *plane* is a cheaper sufficient condition — and it refuses a whole family the
// engine serves: a body standing well clear of a bore whose wall plane, extended, happens to
// pass through it. These lock what the boundary test opens and what it does not.

/// A plate with a bore, fused to a boss standing far away in `x` — whose `y = 12` wall **plane**
/// stands only 2 from the bore's axis (`r = 3`). The boss's face is 20 away; nothing meets.
fn plate_bore_and_boss(hole_y: f64, boss_z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([8.0, hole_y, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        7.0,
    )
    .solid;
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    let boss = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([28.0, 4.0, boss_z0]),
        Point3::from_array([36.0, 12.0, 8.0]),
    );
    m.rebuild_adjacency();
    (m, holed, boss)
}

#[test]
fn a_boss_whose_wall_plane_crosses_a_distant_bore_still_fuses() {
    let (mut m, holed, boss) = plate_bore_and_boss(10.0, 5.0);
    let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("the boss fuses");
    assert_eq!(out.len(), 1, "one body");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★ **The d = 0 member of the family: the wall plane runs exactly through the bore's
/// axis.** The gate passes it the same way (the boss's face clears the
/// footprint along the span), and the rulings road stays **silent** — its trigger is the
/// gate's carried `crossings` record, which is empty for every gate-passed pair. This
/// geometry had no corpus row (its siblings are all d = 2), and it is the population an
/// unconditionally-firing contribution arm broke — measured, 14 arc tests red — so it locks
/// "an empty record changes nothing".
#[test]
fn a_boss_whose_wall_plane_holds_the_bores_axis_still_fuses() {
    let (mut m, holed, boss) = plate_bore_and_boss(12.0, 5.0);
    let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("the boss fuses");
    assert_eq!(out.len(), 1, "one body");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// The d = 0 member with **caps in play**: a cylinder tool standing clear of the plate, its
/// axis exactly on the plate's `x = 40` wall plane and its caps strictly inside the plate's
/// height (no coplanar contact anywhere). Two disjoint bodies is the valid fuse answer, and
/// the chord arm — like the ruling arm — stays silent behind the empty record.
#[test]
fn a_capped_tool_on_the_plates_wall_plane_fuses_as_two_bodies() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let tool = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([40.0, 30.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, plate, tool).expect("a distant fuse");
    assert_eq!(out.len(), 2, "two disjoint bodies");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let mut vols: Vec<f64> = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .collect();
    vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
    let want = [std::f64::consts::PI * 9.0 * 3.0, 4000.0];
    assert!(
        (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
        "{vols:?} vs {want:?}"
    );
}

/// The `Cut` twin — the same wall planes reach the gate whichever way the operation runs.
/// The tool is sunk into the plate (`z` from 3) so it removes material: a body merely *resting*
/// on the top face is a coplanar contact and meets a different guard entirely, which would make
/// this fixture measure that instead of the wall rule.
#[test]
fn the_same_boss_cuts_the_bored_plate() {
    let (mut m, holed, boss) = plate_bore_and_boss(10.0, 3.0);
    let out = crate::boolean(&mut m, BoolKind::Cut, holed, boss).expect("the boss cuts");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 - 8.0 * 8.0 * 2.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **Two bores, one wall plane crossing both.** The rule is per (cylinder, class), so a second
/// cylinder is a second set of questions about the same face — and the face answers for each.
#[test]
fn one_wall_plane_may_cross_two_bores() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let mut holed = plate;
    for x in [8.0, 20.0] {
        let hole = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([x, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            3.0,
            7.0,
        )
        .solid;
        m.rebuild_adjacency();
        holed = crate::boolean(&mut m, BoolKind::Cut, holed, hole).expect("bore")[0];
    }
    let boss = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([28.0, 4.0, 5.0]),
        Point3::from_array([36.0, 12.0, 8.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("two bores and a boss");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - 2.0 * std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **A wall face that really does cross the bore builds.** The plate's own `y = 20`
/// wall would clear, so the tool here is a slab whose face runs right across the hole — the
/// wall rule's true population. The gate records the pair and the tracer cuts the bore's lateral
/// along two rulings
/// (2 from the axis, r = 3) and its caps along the chord: one body, the exact volume.
#[test]
fn a_wall_face_that_really_crosses_the_bore_builds() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        7.0,
    )
    .solid;
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    // A slab covering x ∈ [0, 40]: its y = 12 face passes straight through the bore.
    let slab = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 12.0, 0.0]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect("a real crossing builds");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 1);
    assert!(nacre_validate::validate(&m).is_empty());
    // The slab's box less the bore's `y ≥ 12` segment (d = 2, r = 3), taken from the bored plate.
    let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let pi = std::f64::consts::PI;
    let want = 4000.0 - 45.0 * pi - (1600.0 - 5.0 * seg(2.0, 3.0));
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **Exact tangency, and the operation is what refuses it.** The slab's face stands exactly
/// `r` from the axis, so it touches the bore's lateral along one line. The gate passes
/// that — what convicts this shape is the *verdict*: cutting the slab away leaves the
/// material as the **two wedges** between the parabola and the plane, which meet only on the
/// line, and both are bounded by the same faces so they are one solid. `SelfTouchingResult`.
///
/// ★ Its siblings are the control: the same wall with the boss *inside* the plate builds under
/// `Fuse` and `Common` and only `Cut` pinches, and a boss tangent from **outside** comes back
/// as two valid bodies. So this is not "the gate refuses tangencies"; it is one operation's
/// answer about one shape.
#[test]
fn a_wall_face_tangent_to_the_bore_pinches_under_cut() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        7.0,
    )
    .solid;
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    let slab = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 13.0, 0.0]), // y = 13 is exactly r = 3 from y = 10
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    m.rebuild_adjacency();
    let err = crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect_err("tangency");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "{err:?}"
    );
}

/// ★★★★★ **A blind stud tangent to the wall, and the three answers it must get.** A unit cube
/// and a stud whose axis stands `0.3` from the origin with `r = 0.2`: the wall `x = 0.5` is
/// **exactly** `r` away, which is what a round set of dimensions produces. The three operations
/// differ, and they differ without any of them being named — `keep` is asked about the three
/// regions beside the tangent line and the answers fall out.
///
/// ★ This is not *the user's script*: the app's `cylinder({center})`
/// anchors the cylinder at mid-height, so the user's stud runs `z ∈ [−1, 1]` **through** the
/// cube. That model is the next test; this one — base at
/// `z = 0`, one cap pinched — is its blind neighbour, and the shape OCCT measured.
///
/// Volumes **derived, not copied**: the stud's footprint is `x ∈ [0.1, 0.5] × y ∈ [−0.2, 0.2]`,
/// inside the cube's, and only `z ∈ [0, 0.5]` overlaps — so the shared volume is `π r² · 0.5`.
#[test]
fn a_blind_stud_tangent_to_the_wall_gives_three_answers() {
    let pi = std::f64::consts::PI;
    let (whole, shared) = (pi * 0.04 * 2.0, pi * 0.04 * 0.5);
    let build = || {
        let mut m = Model::new();
        let cube = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([-0.5; 3]),
            Point3::from_array([0.5; 3]),
        );
        let stud = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([0.3, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.2,
            2.0,
        )
        .solid;
        m.rebuild_adjacency();
        (m, cube, stud)
    };
    for (kind, want) in [
        (BoolKind::Fuse, 1.0 + whole - shared),
        (BoolKind::Common, shared),
    ] {
        let (mut m, a, b) = build();
        let out = crate::boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1, "{kind:?}");
        assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
    }
    // Cutting the stud out leaves the material as the **two wedges** beside the tangent line,
    // which meet only on it — one solid whose surface touches itself.
    let (mut m, a, b) = build();
    let err = crate::boolean(&mut m, BoolKind::Cut, a, b).expect_err("the bore pinches");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "{err:?}"
    );
}

/// ★★★★★ **The user's script — `cuboid()` fused with `cylinder({r: 0.2, h: 2, center: [0.3,0,0]})`
/// — and its three answers.** `center` is mid-height, so the stud runs `z ∈ [−1, 1]` through
/// the cube and is tangent to the wall `x = 0.5` along the cube's whole height. Both caps are
/// pinched (the stud's circle touches each square's edge at one point) — and
/// both are bridged and drawn: this is the model the tessellator was opened for.
///
/// Volumes derived: the stud is `π r² · 2`, the part inside the cube `π r² · 1`.
#[test]
fn the_users_through_stud_gives_three_answers() {
    let pi = std::f64::consts::PI;
    let (whole, shared) = (pi * 0.04 * 2.0, pi * 0.04 * 1.0);
    let build = || {
        let mut m = Model::new();
        let cube = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([-0.5; 3]),
            Point3::from_array([0.5; 3]),
        );
        let stud = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([0.3, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.2,
            2.0,
        )
        .solid;
        m.rebuild_adjacency();
        (m, cube, stud)
    };
    for (kind, want) in [
        (BoolKind::Fuse, 1.0 + whole - shared),
        (BoolKind::Common, shared),
    ] {
        let (mut m, a, b) = build();
        let out = crate::boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1, "{kind:?}");
        assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
        if kind == BoolKind::Fuse {
            let faces: usize = std::iter::once(m.solid(out[0]).outer)
                .map(|sh| m.shell(sh).faces.len())
                .sum();
            assert_eq!(
                faces, 10,
                "four walls, two pinched caps, two bands, two stud caps"
            );
            crate::tests::mesh_covers_faces("the user's through stud", &m, &out);
        }
    }
    // Cutting the stud out leaves the two wedges beside the tangent line, meeting only on
    // it, the whole height of the cube — one solid whose surface touches itself.
    let (mut m, a, b) = build();
    let err = crate::boolean(&mut m, BoolKind::Cut, a, b).expect_err("the bore pinches");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "{err:?}"
    );
}

/// ★★★★★ **The user's cross studs — a stud along `z` and one along `y` through the cube,
/// fused in either order.** After the first fuse the `z` stud's remaining lateral faces sit
/// at `|z| ≥ 0.5` and the `y` stud's whole surface within `|z| ≤ 0.2`: the two cylinder
/// classes share no face, and the arrangement builds the result: 14 faces, volume
/// `1 + 0.08π`, meshed. Refusing the pair on the distance between their *axes* — zero, since
/// they cross at the origin — would be a fact about two infinite surfaces, not about any face.
///
/// Volumes: each stud is `0.08π`, half of it inside the cube; the second stud meets the
/// first only inside the cube.
#[test]
fn the_users_cross_studs_build_in_either_order() {
    let pi = std::f64::consts::PI;
    let z_stud = |m: &mut Model| {
        crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([0.0, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.2,
            2.0,
        )
        .solid
    };
    let y_stud = |m: &mut Model| {
        crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([0.0, -1.0, 0.0]),
            Vector3::from_array([0.0, 1.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        )
        .solid
    };
    for z_first in [true, false] {
        let build = || {
            let mut m = Model::new();
            let cube = crate::fixtures::cuboid(
                &mut m,
                Point3::from_array([-0.5; 3]),
                Point3::from_array([0.5; 3]),
            );
            let (first, second) = if z_first {
                (z_stud(&mut m), y_stud(&mut m))
            } else {
                (y_stud(&mut m), z_stud(&mut m))
            };
            m.rebuild_adjacency();
            let studded = crate::boolean(&mut m, BoolKind::Fuse, cube, first)
                .expect("the first stud fuses")[0];
            m.rebuild_adjacency();
            (m, studded, second)
        };
        for (kind, want, faces_want) in [
            (BoolKind::Fuse, 1.0 + 0.08 * pi, Some(14)),
            (BoolKind::Common, 0.04 * pi, Some(3)),
            (BoolKind::Cut, 1.0, None),
        ] {
            let (mut m, a, b) = build();
            let out = crate::boolean(&mut m, kind, a, b)
                .unwrap_or_else(|e| panic!("{kind:?} (z first: {z_first}): {e:?}"));
            m.rebuild_adjacency();
            assert_eq!(out.len(), 1, "{kind:?} (z first: {z_first})");
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{kind:?}: {:?}",
                nacre_validate::validate(&m)
            );
            let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
            assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
            if let Some(n) = faces_want {
                let faces = m.shell(m.solid(out[0]).outer).faces.len();
                assert_eq!(faces, n, "{kind:?} (z first: {z_first})");
            }
            crate::tests::mesh_covers_faces("the user's cross studs", &m, &out);
        }
    }
}

/// **The same two studs without the cube really cross**, and the gate still says so: their
/// faces meet along a quartic curve (M6b). This is the negative control of the face-level
/// clearance — letting perpendicular pairs through blindly makes this fixture come back as
/// two separate, individually valid solids, which is wrong.
#[test]
fn crossing_studs_are_still_a_cylinder_pair_that_meets() {
    let mut m = Model::new();
    let z = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.2,
        2.0,
    )
    .solid;
    let y = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, -1.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let err = crate::boolean(&mut m, BoolKind::Fuse, z, y).expect_err("the studs cross");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::CylinderPairContact,
                ..
            }
        ),
        "{err:?}"
    );
}

/// **A cross above a stud**: a `y` cylinder at `z = 2.5` over a `z` cylinder ending at
/// `z = 2`. Their axes cross (distance zero), so the surface rule refuses, but the upper
/// one's reach along `z` is `[2.3, 2.7]` and the lower one's face spans `[0, 2]` — clear
/// by the face rule: two solids that never touch, and the boolean says so.
#[test]
fn a_cross_above_a_stud_is_two_solids() {
    let mut m = Model::new();
    let low = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0; 3]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.2,
        2.0,
    )
    .solid;
    let high = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.0, -1.0, 2.5]),
        Vector3::from_array([0.0, 1.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, low, high).unwrap_or_else(|e| panic!("{e:?}"));
    m.rebuild_adjacency();
    assert_eq!(out.len(), 2, "two bodies that never touch");
    assert!(nacre_validate::validate(&m).is_empty());
    for s in &out {
        let v = nacre_props::mass_props(&m, *s).unwrap().volume;
        assert!((v - 0.08 * std::f64::consts::PI).abs() < 1e-9, "{v}");
    }
}

/// **An oblique cross above the stud** — the same, with the upper axis `(0, 3, 4)`. The pair
/// rule would clear it too (`lateral_reach`'s `d·m ≠ 0` arm), but the upper cylinder's caps are
/// planes oblique to the lower axis, so the plane–cylinder gate asks first — of the *faces*
/// (`oblique_plane_clears`): does every lateral face of the stud miss that infinite plane?
///
/// Both answers, one axis apart in height. At `z = 3` the lower cap's plane `3y + 4z = 9` stands
/// past the stud's whole reach along its normal (`3·0.2 + 4·2 = 8.6`): cleared, and the fuse is
/// two bodies. At `z = 2.4` the plane `3y + 4z = 6.6` runs through the stud's slab, the face
/// question cannot show a miss, and the pair is refused by its caps — a lock on today's name for
/// the cell that opens the ellipse road. (A 3-4-5 axis, so the upper cylinder is stated in the
/// world: an axis like `(0, 1, 1)` has no rational unit frame, rides a frame node, and the gate
/// refuses it earlier, as having no world statement.)
#[test]
fn an_oblique_cross_above_the_stud_is_decided_by_its_caps_planes() {
    for (z, clears) in [(3.0, true), (2.4, false)] {
        let mut m = Model::new();
        let low = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([0.0; 3]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.2,
            2.0,
        )
        .solid;
        let high = crate::fixtures::cylinder(
            &mut m,
            Point3::from_array([0.0, -1.0, z]),
            Vector3::from_array([0.0, 3.0, 4.0]),
            0.2,
            2.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, low, high);
        if clears {
            let out = out.unwrap_or_else(|e| panic!("z = {z}: {e:?}"));
            m.rebuild_adjacency();
            assert_eq!(out.len(), 2, "z = {z}: two bodies that never touch");
            assert!(nacre_validate::validate(&m).is_empty());
            for s in &out {
                let v = nacre_props::mass_props(&m, *s).unwrap().volume;
                assert!((v - 0.08 * std::f64::consts::PI).abs() < 1e-9, "{v}");
            }
        } else {
            let err = out.expect_err("oblique caps across the stud");
            assert!(
                matches!(
                    err,
                    BoolError::Rejected {
                        reason: RejectReason::ObliqueCylinderCut,
                        ..
                    }
                ),
                "z = {z}: {err:?}"
            );
        }
    }
}

/// ★★★★★ **The bridge pre-pass splits the wall's two straight edges under the through stud.**
///
/// A stud through the cube (`center` anchoring, `z ∈ [−1, 1]`) pinches **both** caps: on
/// each, the stud's circle touches the square's `x = 0.5` edge at one point, and that edge
/// is shared with the `x = 0.5` wall. Before the pre-pass a straight edge carries exactly its
/// two ends, so the touching sample — a vertex of the circle — is a vertex of neither the
/// square ring nor the wall. The pre-pass inserts it into both edges' polylines, so the cap
/// gains its twin and the wall gains the same vertex (no T-vertex, no crack).
///
/// Read through the test-only `bridge_report`, because `tessellate` still refuses these caps
/// (the bridge itself is the next commit) and discards everything on the way out.
#[test]
fn the_bridge_prepass_splits_the_walls_edges_under_both_caps() {
    let mut m = Model::new();
    let cube = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5; 3]),
        Point3::from_array([0.5; 3]),
    );
    // The seam on `−y`, off the touching point `(0.5, 0)` — on `+x` it would *be* that point, a
    // vertex of the b-rep rather than a sample of the circle.
    let stud = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.3, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.2,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, cube, stud).expect("the through stud fuses");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    let (report, t) = nacre_tess::bridge_report(&m, &nacre_tess::TessConfig::default());
    assert!(
        report.declined.is_empty(),
        "declined: {:?}",
        report.declined
    );
    assert_eq!(report.splits.len(), 2, "splits: {:?}", report.splits);
    assert_ne!(
        report.splits[0].edge, report.splits[1].edge,
        "one split per cap"
    );
    for s in &report.splits {
        assert!(
            matches!(m.edge_curve(s.edge), nacre_geom::Curve::Line(_)),
            "the split edge is straight"
        );
        let poly = &t.by_edge[&s.edge];
        assert_eq!(poly.len(), 3, "two ends and the touching sample");
        assert_eq!(poly[s.at], s.vertex);
        assert!(
            matches!(
                t.vertices.get(s.vertex).origin,
                nacre_tess::TessOrigin::OnEdge { .. }
            ),
            "the inserted vertex is the circle's own sample, not a new one"
        );
    }
    // The blind stud (base anchoring) pinches one cap only.
    let mut m = Model::new();
    let cube = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5; 3]),
        Point3::from_array([0.5; 3]),
    );
    let stud = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.3, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.2,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    crate::boolean(&mut m, BoolKind::Fuse, cube, stud).expect("the blind stud fuses");
    m.rebuild_adjacency();
    let (report, _) = nacre_tess::bridge_report(&m, &nacre_tess::TessConfig::default());
    assert!(report.declined.is_empty(), "{:?}", report.declined);
    assert_eq!(report.splits.len(), 1, "{:?}", report.splits);
}

/// ★★★★★ **A tangency from *outside* is two bodies, not a refusal** — the control that says
/// the rule is not "reject every tangency", and the measurement that corrected it.
///
/// The boss's axis stands at `x = −0.5` with `r = 0.5`, so it touches the plate's wall `x = 0`
/// from the void side. `Fuse` keeps the lens (the boss) and the far side (the plate) but not
/// the wedges between — two lumps meeting only on the line — and because they share no face
/// they are **two valid solids**, which is what the kernel returns. `Cut` removes nothing and
/// `Common` is empty.
#[test]
fn a_boss_tangent_from_outside_is_two_bodies() {
    let pi = std::f64::consts::PI;
    let build = || {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([-0.5, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        (m, plate, boss)
    };
    for (kind, bodies, want) in [
        (BoolKind::Fuse, 2, 32.0 + pi),
        (BoolKind::Cut, 1, 32.0),
        (BoolKind::Common, 0, 0.0),
    ] {
        let (mut m, a, b) = build();
        let out = crate::boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
        m.rebuild_adjacency();
        assert_eq!(out.len(), bodies, "{kind:?}");
        assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
        let v: f64 = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
    }
}

/// ★★★★★ **The same two lumps, joined elsewhere — and this one *is* a pinch.** The block has a
/// notch, the boss stands in it tangent to the notch's wall `x = 2` from the void side, and it
/// overlaps the block in `y` (the notch is `0.8` wide, the boss `1.0`). So the lens and the far
/// side are the same body, and the tangent line is where that body's surface meets itself.
///
/// ★ This is why the verdict asks the **grouping** and not only the geometry: without it, the
/// answer here would be the one the fixture above earns, and this solid shipped `Ok` with
/// `validate` clean and a mesh — measured before the check existed.
#[test]
fn a_boss_tangent_in_a_notch_pinches_the_block_it_joins() {
    let mut m = Model::new();
    let outer = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([6.0, 6.0, 2.0]),
    );
    let notch = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([2.0, 2.6, -0.5]),
        Point3::from_array([6.5, 3.4, 2.5]),
    );
    m.rebuild_adjacency();
    let block = crate::boolean(&mut m, BoolKind::Cut, outer, notch).expect("notch")[0];
    m.rebuild_adjacency();
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.5, 3.0, -0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    let err = crate::boolean(&mut m, BoolKind::Fuse, block, boss).expect_err("a joined pinch");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::SelfTouchingResult,
                ..
            }
        ),
        "{err:?}"
    );
}

/// ★★★★★ **A tangent wall through a window in the lateral does not touch it.** The slab's face
/// `x = 1` is tangent to [`crate::fixtures::windowed_boss`]'s cylinder along `x = 1, y = 0`, but
/// only over `z ∈ [1.7, 2.3]`, inside the window where the lateral face is absent — so cutting the
/// boss from the slab leaves the slab's face whole and a bridge of slab between the bore's two
/// sides. One body, its volume the slab less the cylinder's slice less what the window took out of
/// it: the window's cross-section inside the unit circle is `w·√(1−w²) + asin(w) − w`.
///
/// ★ The pair beside it is the pinch: a slab longer than the window touches the lateral above and
/// below it, and the result's surface meets itself along those two stretches.
#[test]
fn a_tangent_wall_through_a_lateral_window_builds() {
    for w in [0.3f64, 0.6] {
        for seam in [1.0, -1.0] {
            let mut m = Model::new();
            let boss = crate::fixtures::windowed_boss(&mut m, w, seam);
            let slab = crate::fixtures::cuboid(
                &mut m,
                Point3::from_array([-2.0, -2.0, 1.7]),
                Point3::from_array([1.0, 2.0, 2.3]),
            );
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, BoolKind::Cut, slab, boss)
                .unwrap_or_else(|e| panic!("w {w} seam {seam}: {e:?}"));
            assert_eq!(out.len(), 1, "w {w} seam {seam}");
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty(), "w {w} seam {seam}");
            let window = w * (1.0 - w * w).sqrt() + w.asin() - w;
            let want = 7.2 - (std::f64::consts::PI - window) * 0.6;
            let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
            assert!((v - want).abs() < 1e-9, "w {w} seam {seam}: {v} vs {want}");
        }
    }
}

#[test]
fn a_tangent_wall_longer_than_the_window_still_pinches() {
    for w in [0.3, 0.6] {
        for seam in [1.0, -1.0] {
            let mut m = Model::new();
            let boss = crate::fixtures::windowed_boss(&mut m, w, seam);
            let slab = crate::fixtures::cuboid(
                &mut m,
                Point3::from_array([-2.0, -2.0, 1.0]),
                Point3::from_array([1.0, 2.0, 3.0]),
            );
            m.rebuild_adjacency();
            let err = crate::boolean(&mut m, BoolKind::Cut, slab, boss)
                .expect_err("the slab touches the lateral above and below the window");
            assert!(
                matches!(
                    err,
                    BoolError::Rejected {
                        reason: RejectReason::SelfTouchingResult,
                        ..
                    }
                ),
                "w {w} seam {seam}: {err:?}"
            );
        }
    }
}

/// Every boolean of two disjoint-but-for-a-line operands, both orders: `Fuse` is the two
/// operands as two bodies, each `Cut` the minuend unchanged, `Common` empty — valid, and the
/// volumes the operands' own, read before the boolean.
fn two_bodies_touching_on_a_line(
    name: &str,
    build: impl Fn() -> (Model, Handle<Solid>, Handle<Solid>),
) {
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        for swapped in [false, true] {
            let (mut m, a, b) = build();
            let (a, b) = if swapped { (b, a) } else { (a, b) };
            let vol = |m: &Model, s| nacre_props::mass_props(m, s).unwrap().volume;
            let (va, vb) = (vol(&m, a), vol(&m, b));
            let out = crate::boolean(&mut m, kind, a, b)
                .unwrap_or_else(|e| panic!("{name}: {kind:?} swapped {swapped}: {e:?}"));
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{name}: {kind:?} swapped {swapped}: {vs:?}");
            let mut got: Vec<f64> = out.iter().map(|&s| vol(&m, s)).collect();
            got.sort_by(f64::total_cmp);
            let mut want = match kind {
                BoolKind::Fuse => vec![va, vb],
                BoolKind::Cut => vec![va],
                BoolKind::Common => vec![],
            };
            want.sort_by(f64::total_cmp);
            assert_eq!(
                got.len(),
                want.len(),
                "{name}: {kind:?} swapped {swapped}: bodies"
            );
            for (g, w) in got.iter().zip(&want) {
                assert!(
                    (g - w).abs() < 1e-9,
                    "{name}: {kind:?} swapped {swapped}: {got:?} vs {want:?}"
                );
            }
        }
    }
}

/// ★ **A face far away on the wall's plane joins nothing.** The boss touches the box's wall
/// `x = 0` from outside along `(0, 0)` — the lens and the far side, two bodies — and it is joined
/// by a bar to a block whose face lies on that same plane at `y ∈ [10, 12]`, nowhere near the box.
/// The union is still the two operands touching along the line. Read by the plane's class, the
/// verdict found the wall's class and the cylinder's in the boss's body and refused
/// `SelfTouchingResult`; the lens and the far side are in one body only when the result is one.
#[test]
fn a_far_face_on_the_walls_plane_joins_nothing() {
    two_bodies_touching_on_a_line("boss, bar and block", || {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([-2.0, -3.0, 0.0]),
            Point3::from_array([0.0, 3.0, 1.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([1.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            1.0,
            1.0,
        )
        .solid;
        let bar = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.5, 0.0, 0.0]),
            Point3::from_array([1.5, 11.0, 1.0]),
        );
        m.rebuild_adjacency();
        let b = crate::boolean(&mut m, BoolKind::Fuse, boss, bar).expect("boss and bar")[0];
        let block = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, 10.0, 0.0]),
            Point3::from_array([2.0, 12.0, 1.0]),
        );
        m.rebuild_adjacency();
        let b = crate::boolean(&mut m, BoolKind::Fuse, b, block).expect("and the block")[0];
        m.rebuild_adjacency();
        (m, a, b)
    });
}

/// ★ **A row through a wall face's notch touches nothing.** The box's face `x = 0` is a U — a slot
/// `y ∈ [−2, 2]`, `z ∈ [1, 3]` cut through the box — and the boss floats in the slot's void, its
/// surface tangent to the plane `x = 0` along `(0, 0)` where the face is not. Over the boss's
/// height the line runs through the notch alone, so the face touches nothing on the lateral and
/// the gate writes no row — though the U's corners lie on both sides of the line. Every boolean
/// builds.
#[test]
fn a_tangency_through_a_wall_faces_notch_touches_nothing() {
    two_bodies_touching_on_a_line("notched box", || {
        let mut m = Model::new();
        let block = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, -3.0, 0.0]),
            Point3::from_array([4.0, 3.0, 2.0]),
        );
        let slot = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([-1.0, -2.0, 1.0]),
            Point3::from_array([5.0, 2.0, 3.0]),
        );
        m.rebuild_adjacency();
        let a = crate::boolean(&mut m, BoolKind::Cut, block, slot).expect("the slot")[0];
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([1.0, 0.0, 1.2]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            1.0,
            0.6,
        )
        .solid;
        m.rebuild_adjacency();
        (m, a, b)
    });
}

/// Where the tangent line meets the wall plane `x = 0`: a slot `z ∈ [1, 3]` cut through the
/// block, which leaves the face a U (`Notch`, the block 2 high) or a face with a rectangular hole
/// (the rest, the block 4 high).
#[derive(Clone, Copy, Debug)]
enum WallGap {
    /// The slot `y ∈ [−2, 2]`: the line runs through the U's gap.
    Notch,
    /// The slot `y ∈ [−2, 2]`: the line runs through the hole.
    Hole,
    /// The slot `y ∈ [0, 2]`: the line runs along the hole's edge `y = 0`, which the slot's wall
    /// `y = 0` holds too.
    HoleAlongLine,
    /// `Hole`, and a step cut out of the block's far end, `x ∈ [4, 6]`, `y ≥ 0`, whose face `y = 0`
    /// holds the line far from it.
    HoleBesideStep,
}

/// **A wide boss tangent to a wall face only where the face is not** — and cutting material
/// beside it. The block `[0, 6] × [−3, 3]` has the gap of [`WallGap`] in its face `x = 0`; the
/// boss (axis `(2.5, 0)` along `+z`, `r = 2.5`) touches the plane `x = 0` along `(0, 0)`, which
/// within the boss's height runs through the gap alone, and it reaches past `|y| = 2` into the
/// block's material on both sides of the slot. Returns the model, the block and the boss.
fn wide_boss_in_a_wall_gap(gap: WallGap) -> (Model, Handle<Solid>, Handle<Solid>) {
    let (top, z0, h) = match gap {
        WallGap::Notch => (2.0, 1.2, 0.6),
        _ => (4.0, 1.5, 1.0),
    };
    let slot_y0 = match gap {
        WallGap::HoleAlongLine => 0.0,
        _ => -2.0,
    };
    let mut m = Model::new();
    let block = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, -3.0, 0.0]),
        Point3::from_array([6.0, 3.0, top]),
    );
    let slot = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, slot_y0, 1.0]),
        Point3::from_array([7.0, 2.0, 3.0]),
    );
    m.rebuild_adjacency();
    let mut a = crate::boolean(&mut m, BoolKind::Cut, block, slot).expect("the slot")[0];
    if let WallGap::HoleBesideStep = gap {
        let step = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([4.0, 0.0, -1.0]),
            Point3::from_array([7.0, 4.0, 5.0]),
        );
        m.rebuild_adjacency();
        a = crate::boolean(&mut m, BoolKind::Cut, a, step).expect("the step")[0];
    }
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.5, 0.0, z0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.5,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    (m, a, b)
}

/// ★★ **A tangency through a wall face's gap writes no row, so the cut beside it builds.** The
/// line `(0, 0)` touches the plane `x = 0` only inside the slot over the boss's height, so no
/// wall face touches it there; the boss still overlaps the block in the two circular segments past
/// `|y| = 2`. The U face (and the face around the hole) has corners on both sides of the line, so
/// a row read off the outer loop alone convicts `block − boss` of a pinch it does not have
/// (`SelfTouchingResult`, measured).
///
/// The oracle is the inputs': the overlap is two segments of the boss's circle cut by the chords
/// `y = ±2`, `2h(r²·acos(d/r) − d·√(r² − d²))` with `d = 2`.
#[test]
fn a_tangency_through_a_wall_faces_gap_beside_material_builds() {
    let (r, d): (f64, f64) = (2.5, 2.0);
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    for gap in [WallGap::Notch, WallGap::Hole] {
        let h = match gap {
            WallGap::Notch => 0.6,
            _ => 1.0,
        };
        let overlap = 2.0 * segment * h;
        for (kind, swapped) in [
            (BoolKind::Fuse, false),
            (BoolKind::Fuse, true),
            (BoolKind::Cut, false),
            (BoolKind::Cut, true),
        ] {
            let (mut m, block, boss) = wide_boss_in_a_wall_gap(gap);
            let vol = |m: &Model, s| nacre_props::mass_props(m, s).unwrap().volume;
            let (vblock, vboss) = (vol(&m, block), vol(&m, boss));
            let (a, b) = if swapped {
                (boss, block)
            } else {
                (block, boss)
            };
            let out = crate::boolean(&mut m, kind, a, b)
                .unwrap_or_else(|e| panic!("{gap:?} {kind:?} swapped {swapped}: {e:?}"));
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{gap:?} {kind:?} swapped {swapped}: {vs:?}");
            assert_eq!(out.len(), 1, "{gap:?} {kind:?} swapped {swapped}");
            let want = match (kind, swapped) {
                (BoolKind::Fuse, _) => vblock + vboss - overlap,
                (_, false) => vblock - overlap,
                (_, true) => vboss - overlap,
            };
            let got = vol(&m, out[0]);
            assert!(
                (got - want).abs() < 1e-9,
                "{gap:?} {kind:?} swapped {swapped}: {got} vs {want}"
            );
        }
    }
}

/// **The pinch formula, as a truth table** — four configurations × three operations × both
/// owner orders, checked against the geometry by hand rather than against an engine run.
///
/// Read the rows as: the cylinder's side of the wall plane is (or is not) the wall face's
/// material side, and the cylinder keeps its material inside (a boss) or outside (a bore).
/// ★ The third row's `Common` is the one that says this is not "a rule about `Cut`": a bore
/// intersected with the wall's solid leaves the two wedges alone. ★ The lens and the far side
/// come apart in one cell only — the union with a boss touching the wall from outside, in either
/// owner order — which is what lets the verdict read that shape as the two operands' own.
#[test]
fn the_pinch_formula_is_a_truth_table() {
    use crate::assembly::{Pinch, pinch};
    use crate::planes::SolidSide::{A, B};
    let (w, lf) = (Some(Pinch::Wedges), Some(Pinch::LensAndFar));
    // (lens_in_wall_solid, cyl_orient, [Fuse, Cut, Common] with the wall on A, … on B).
    for (lens, orient, on_a, on_b) in [
        (true, 1i8, [None, w, None], [None, None, None]), // an inside boss: only A − B pinches
        (false, 1, [lf, None, None], [lf, None, None]),   // an outside boss: only Fuse, lens + far
        (true, -1, [None, None, w], [None, None, w]),     // a bore on the material side: Common
        (false, -1, [None, None, None], [None, w, None]), // a bore on the void side: B − A only
    ] {
        for (i, kind) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .enumerate()
        {
            for (side, want) in [(A, on_a[i]), (B, on_b[i])] {
                assert_eq!(
                    pinch(kind, side, lens, orient),
                    want,
                    "wall on {side:?}, lens {lens}, orient {orient}, {kind:?}"
                );
            }
        }
    }
}

/// ★★★★★ **The disk on a wall plane is read, and it clears.**
///
/// The plate carries a **crosswise** bore, so the plane `x = 30` holds that bore's circular cap
/// face: an outer loop of one arc and one seam vertex. That plane is also parallel to the
/// vertical drill's axis and passes within `r` of it, so the face test is reached. Answering
/// «cannot read this shape» there would stop the boolean — a sound barrier
/// (a rule reading vertices alone would find "every vertex on one side" true
/// of a **single point**) but one that turns away the true answer along with the false one.
///
/// ☑ The disk's own numbers: centre `(30, 10, 5)` radius `2`, against a strip centred on
/// `y = 4` of half-width `√(9 − 1) ≈ 2.83`. Six apart, so it clears by more than its radius —
/// and now says so. The volume is the oracle that the answer is not merely *an* answer.
#[test]
fn a_disk_face_on_a_wall_plane_is_read_and_clears() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 10.0]),
    );
    // A crosswise bore along +X, ending inside the plate at x = 30: its cap lies on x = 30.
    let cross = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([-1.0, 10.0, 5.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, -1.0]),
        2.0,
        31.0,
    )
    .solid;
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, cross).expect("crosswise bore")[0];
    // A vertical drill whose axis stands 1 from the plane x = 30 — inside its radius 3 — and
    // ★ **6 from the crosswise bore's axis**, clear of the radius sum 5, so the pair rule is
    // not what answers here. (At y = 10 the two axes actually meet and this fixture would be
    // measuring `CylinderPairContact` instead — the adjacent proposition.)
    let drill = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([31.0, 4.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        12.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, bored, drill).expect("the drill cuts");
    assert_eq!(out.len(), 1, "one body");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    // 40 × 20 × 10, less the crosswise bore (r 2, thirty long inside) and the drill (r 3,
    // through the ten of thickness). Their axes pass six apart, clear of the radius sum.
    let want = 8000.0 - std::f64::consts::PI * (4.0 * 30.0 + 9.0 * 10.0);
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★★ **The same boss's `Common` is two thin pieces, and their depths are decided.** What the
/// boss and the block share is the two circular segments past `|y| = 2` (and, with the slot's
/// edge on the tangent line, the boss's half disk `y < 0` beside one of them) — two components,
/// and the depth of each is asked of the other by rays (`assembly::grouping`). The rays from the
/// chords' middles run in the plane through the boss's axis and its seam, so every arc root they
/// meet sits on the seam: those are ordered as the cyclic order's first point
/// (`nacre_exact::quad::SeamOrder::seam_first`), and a root at the probe on the circle off the
/// arc is no boundary of that ring. Read as ties instead, every probe abstained and the boolean
/// was refused `NoClearRay` here, in both orders (measured).
///
/// The oracle is the inputs': a segment is `h(r²·acos(d/r) − d·√(r² − d²))` with `d = 2`, the half
/// disk `h·πr²/2`.
#[test]
fn a_wide_boss_in_a_wall_gap_common_builds() {
    let (r, d): (f64, f64) = (2.5, 2.0);
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let half_disk = std::f64::consts::PI * r * r / 2.0;
    for gap in [WallGap::Notch, WallGap::Hole, WallGap::HoleAlongLine] {
        let (h, near) = match gap {
            WallGap::Notch => (0.6, segment),
            WallGap::HoleAlongLine => (1.0, half_disk),
            _ => (1.0, segment),
        };
        let mut want = [h * segment, h * near];
        want.sort_by(f64::total_cmp);
        for swapped in [false, true] {
            let (mut m, block, boss) = wide_boss_in_a_wall_gap(gap);
            let (a, b) = if swapped {
                (boss, block)
            } else {
                (block, boss)
            };
            let out = crate::boolean(&mut m, BoolKind::Common, a, b)
                .unwrap_or_else(|e| panic!("{gap:?} swapped {swapped}: {e:?}"));
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{gap:?} swapped {swapped}: {vs:?}");
            let mut got: Vec<f64> = out
                .iter()
                .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                .collect();
            got.sort_by(f64::total_cmp);
            assert_eq!(got.len(), 2, "{gap:?} swapped {swapped}: {got:?}");
            for (g, w) in got.iter().zip(&want) {
                assert!(
                    (g - w).abs() < 1e-9,
                    "{gap:?} swapped {swapped}: {got:?} vs {want:?}"
                );
            }
        }
    }
}

/// ★★ **A hole whose edge lies on the tangent line hands the line to the structure.** The slot
/// `y ∈ [0, 2]` leaves the face `x = 0` a hole whose edge `y = 0` is the tangent line over the
/// boss's height, and the slot's own wall `y = 0` holds the line too. The face ends on the line
/// there rather than crossing it, so the row says so (`runs_through` false), the line is an edge
/// the wall face and the slot's wall both end on, and the shell guard judges the contact. The face
/// has corners on both sides of the line, so a row read off the outer loop alone is no edge, and
/// the gate refuses every operation `TangentLineInAnotherPlane` (measured).
///
/// The oracle is the inputs': the boss meets the block in its half disk `y < 0` and the segment
/// past `y = 2`, `h(πr²/2 + r²·acos(d/r) − d·√(r² − d²))` with `d = 2`.
#[test]
fn a_hole_edge_on_the_tangent_line_hands_the_line_to_the_structure() {
    let (m, block, boss) = wide_boss_in_a_wall_gap(WallGap::HoleAlongLine);
    let rows = crate::arrangement::plane_index_setup(&m, block, boss)
        .expect("the gate passes")
        .tangencies;
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert!(!rows[0].runs_through && rows[0].line_is_an_edge, "{rows:?}");
    let (r, d, h): (f64, f64, f64) = (2.5, 2.0, 1.0);
    let overlap = h
        * (std::f64::consts::PI * r * r / 2.0 + r * r * (d / r).acos()
            - d * (r * r - d * d).sqrt());
    for (kind, swapped) in [
        (BoolKind::Fuse, false),
        (BoolKind::Fuse, true),
        (BoolKind::Cut, false),
        (BoolKind::Cut, true),
    ] {
        let (mut m, block, boss) = wide_boss_in_a_wall_gap(WallGap::HoleAlongLine);
        let vol = |m: &Model, s| nacre_props::mass_props(m, s).unwrap().volume;
        let (vblock, vboss) = (vol(&m, block), vol(&m, boss));
        let (a, b) = if swapped {
            (boss, block)
        } else {
            (block, boss)
        };
        let out = crate::boolean(&mut m, kind, a, b)
            .unwrap_or_else(|e| panic!("{kind:?} swapped {swapped}: {e:?}"));
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} swapped {swapped}: {vs:?}");
        assert_eq!(out.len(), 1, "{kind:?} swapped {swapped}");
        let want = match (kind, swapped) {
            (BoolKind::Fuse, _) => vblock + vboss - overlap,
            (_, false) => vblock - overlap,
            (_, true) => vboss - overlap,
        };
        let got = vol(&m, out[0]);
        assert!(
            (got - want).abs() < 1e-9,
            "{kind:?} swapped {swapped}: {got} vs {want}"
        );
    }
}

/// **A tangent line through a hole, held by a far face's plane, is refused by the rulings road
/// today.** The line `(0, 0)` runs through the hole of `x = 0`, so no wall face touches it and
/// the gate writes no row; but the step's face `y = 0` is a secant on the boss's axis plane, and
/// with a tangent class of the other solid it records the line as a shared ruling. That road does
/// not arrange this shape yet and refuses it by name (`RulingBoundNotYet`) — every operation, either
/// order. The same block with the boss moved off the wall (axis at `x = 2.6`) builds all six, so
/// the refusal is the tangent line's, not the step's. A row read off the outer loop alone has the
/// gate refuse it `TangentLineInAnotherPlane` instead, whose six-region sentence is false here:
/// nothing touches the line.
#[test]
fn a_tangent_line_through_a_hole_held_by_a_far_step_is_refused_today() {
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        for swapped in [false, true] {
            let (mut m, block, boss) = wide_boss_in_a_wall_gap(WallGap::HoleBesideStep);
            let (a, b) = if swapped {
                (boss, block)
            } else {
                (block, boss)
            };
            match crate::boolean(&mut m, kind, a, b) {
                Err(BoolError::Rejected { reason, .. }) => assert_eq!(
                    reason,
                    RejectReason::RulingBoundNotYet,
                    "{kind:?} swapped {swapped}"
                ),
                other => panic!("{kind:?} swapped {swapped}: {other:?}"),
            }
        }
    }
}
