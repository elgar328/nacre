//! Wall faces: the gate asks the boundary, not the infinite plane.

use super::*;

//
// ★★ The uniform-slab theorem's premise is about the other operand's **boundary**. The gate
// used to test the wall's infinite *plane*, which is a cheaper sufficient condition — and it
// refused a whole family the engine serves: a body standing well clear of a bore whose wall
// plane, extended, happens to pass through it. These lock what that opened and what it did not.

/// A plate with a bore, fused to a boss standing far away in `x` — whose `y = 12` wall **plane**
/// stands only 2 from the bore's axis (`r = 3`). The boss's face is 20 away; nothing meets.
fn plate_bore_and_boss(hole_y: f64, boss_z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = m.add_cylinder(
        Point3::from_array([8.0, hole_y, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        7.0,
    );
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    let boss = m.add_cuboid(
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
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let tool = m.add_cylinder(
        Point3::from_array([40.0, 30.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        3.0,
    );
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
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let mut holed = plate;
    for x in [8.0, 20.0] {
        let hole = m.add_cylinder(
            Point3::from_array([x, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        holed = crate::boolean(&mut m, BoolKind::Cut, holed, hole).expect("bore")[0];
    }
    let boss = m.add_cuboid(
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
/// wall rule's true population. It used to be the fence (`WallMeetsLateral`, a
/// reason since retired); the gate records the pair now and the tracer cuts the bore's lateral
/// along two rulings
/// (2 from the axis, r = 3) and its caps along the chord: one body, the exact volume.
#[test]
fn a_wall_face_that_really_crosses_the_bore_builds() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = m.add_cylinder(
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        7.0,
    );
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    // A slab covering x ∈ [0, 40]: its y = 12 face passes straight through the bore.
    let slab = m.add_cuboid(
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
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = m.add_cylinder(
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        7.0,
    );
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    let slab = m.add_cuboid(
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
        let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let stud = m.add_cylinder(
            Point3::from_array([0.3, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
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
        let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let stud = m.add_cylinder(
            Point3::from_array([0.3, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
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
/// classes share no face, and the arrangement builds the result — measured before the gate
/// opened, with the pair let through by hand: 14 faces, volume `1 + 0.08π`, meshed. The
/// gate used to refuse the pair on the distance between their *axes*, zero since they cross
/// at the origin — a fact about two infinite surfaces, not about any face.
///
/// Volumes: each stud is `0.08π`, half of it inside the cube; the second stud meets the
/// first only inside the cube.
#[test]
fn the_users_cross_studs_build_in_either_order() {
    let pi = std::f64::consts::PI;
    let z_stud = |m: &mut Model| {
        m.add_cylinder(
            Point3::from_array([0.0, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        )
    };
    let y_stud = |m: &mut Model| {
        m.add_cylinder(
            Point3::from_array([0.0, -1.0, 0.0]),
            Vector3::from_array([0.0, 1.0, 0.0]),
            0.2,
            2.0,
        )
    };
    for z_first in [true, false] {
        let build = || {
            let mut m = Model::new();
            let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
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
    let z = m.add_cylinder(
        Point3::from_array([0.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    );
    let y = m.add_cylinder(
        Point3::from_array([0.0, -1.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
        0.2,
        2.0,
    );
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
    let low = m.add_cylinder(
        Point3::from_array([0.0; 3]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    );
    let high = m.add_cylinder(
        Point3::from_array([0.0, -1.0, 2.5]),
        Vector3::from_array([0.0, 1.0, 0.0]),
        0.2,
        2.0,
    );
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

/// **An oblique cross above the stud** — the same, with the upper axis `(0, 1, 1)`. The
/// face rule would clear it too (`lateral_reach`'s `d·m ≠ 0` arm), but the pair never
/// reaches the pair loop: the upper cylinder's caps are planes oblique to the lower axis,
/// and the plane–cylinder gate refuses those without asking whether they clear — the same
/// proposition still spelled at surface level there. A lock on today's name, for the cell
/// that opens that arm to flip.
#[test]
fn an_oblique_cross_above_the_stud_is_refused_by_its_caps() {
    let mut m = Model::new();
    let low = m.add_cylinder(
        Point3::from_array([0.0; 3]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    );
    let high = m.add_cylinder(
        Point3::from_array([0.0, -1.0, 3.0]),
        Vector3::from_array([0.0, 1.0, 1.0]),
        0.2,
        2.0,
    );
    m.rebuild_adjacency();
    let err = crate::boolean(&mut m, BoolKind::Fuse, low, high).expect_err("oblique caps");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::ObliqueCylinderCut,
                ..
            }
        ),
        "{err:?}"
    );
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
    let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
    let stud = m.add_cylinder(
        Point3::from_array([0.3, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    );
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
    let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
    let stud = m.add_cylinder(
        Point3::from_array([0.3, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.2,
        2.0,
    );
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
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([-0.5, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
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
    let outer = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([6.0, 6.0, 2.0]),
    );
    let notch = m.add_cuboid(
        Point3::from_array([2.0, 2.6, -0.5]),
        Point3::from_array([6.5, 3.4, 2.5]),
    );
    m.rebuild_adjacency();
    let block = crate::boolean(&mut m, BoolKind::Cut, outer, notch).expect("notch")[0];
    m.rebuild_adjacency();
    let boss = m.add_cylinder(
        Point3::from_array([2.5, 3.0, -0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        3.0,
    );
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

/// **The pinch formula, as a truth table** — four configurations × three operations × both
/// owner orders, checked against the geometry by hand rather than against an engine run.
///
/// Read the rows as: the cylinder's side of the wall plane is (or is not) the wall face's
/// material side, and the cylinder keeps its material inside (a boss) or outside (a bore).
/// ★ The third row's `Common` is the one that says this is not "a rule about `Cut`": a bore
/// intersected with the wall's solid leaves the two wedges alone.
#[test]
fn the_pinch_formula_is_a_truth_table() {
    use crate::assembly::lumps_fall_apart;
    use crate::planes::SolidSide::{A, B};
    // (lens_in_wall_solid, cyl_orient, [Fuse, Cut, Common]) with the wall on side A.
    for (lens, orient, want) in [
        (true, 1i8, [false, true, false]), // an inside boss: only Cut pinches
        (false, 1, [true, false, false]),  // an outside boss: only Fuse splits into lumps
        (true, -1, [false, false, true]),  // a bore on the material side: only Common
        (false, -1, [false, false, false]), // a bore on the void side: never
    ] {
        for (i, kind) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .enumerate()
        {
            assert_eq!(
                lumps_fall_apart(kind, A, lens, orient),
                want[i],
                "wall on A, lens {lens}, orient {orient}, {kind:?}"
            );
        }
    }
    // ★ Swapping the owners is not a symmetry: `Fuse` and `Common` are commutative but `Cut`
    // is not, so the same geometry judged with the wall on `B` reads `A − B` the other way.
    assert!(lumps_fall_apart(BoolKind::Cut, A, true, 1));
    assert!(!lumps_fall_apart(BoolKind::Cut, B, true, 1));
    assert!(lumps_fall_apart(BoolKind::Cut, B, false, -1));
    assert!(!lumps_fall_apart(BoolKind::Cut, A, false, -1));
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
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 10.0]),
    );
    // A crosswise bore along +X, ending inside the plate at x = 30: its cap lies on x = 30.
    let cross = m.add_cylinder(
        Point3::from_array([-1.0, 10.0, 5.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.0,
        31.0,
    );
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, cross).expect("crosswise bore")[0];
    // A vertical drill whose axis stands 1 from the plane x = 30 — inside its radius 3 — and
    // ★ **6 from the crosswise bore's axis**, clear of the radius sum 5, so the pair rule is
    // not what answers here. (At y = 10 the two axes actually meet and this fixture would be
    // measuring `CylinderPairContact` instead — the adjacent proposition.)
    let drill = m.add_cylinder(
        Point3::from_array([31.0, 4.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        12.0,
    );
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
