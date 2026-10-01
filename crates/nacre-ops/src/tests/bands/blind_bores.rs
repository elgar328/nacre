//! A blind bore is usable, not just buildable.

use super::*;

//
// ★★ A cylinder's lateral face touching a ⊥ class at its **rim** contributes a `Graze`, the
// same as any planar wall meeting a class along an edge, and `edge_mask`'s `Graze > Seated`
// precedence exists exactly for the corner where the two disagree — a **blind bore's ceiling**.
// Without the graze, the cap's seated rule flips the wrong label bit, the band's two ends
// contradict each other, and the *next* boolean on that solid comes back
// `CylinderGateUndecided`. A through bore has no ceiling.

/// A plate bored to `hole_h` deep, ready to be operated on again.
fn plate_with_a_bore(hole_h: f64) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([8.0, 10.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        hole_h,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore");
    (m, out[0])
}

/// A boss well clear of the bore (walls 4+ from the axis, r = 3), fused onto the plate.
fn fuse_a_clear_boss(m: &mut Model, s: Handle<Solid>, z0: f64) -> Vec<Handle<Solid>> {
    let boss = crate::fixtures::cuboid(
        m,
        Point3::from_array([4.0, 4.0, z0]),
        Point3::from_array([12.0, 6.0, 8.0]),
    );
    m.rebuild_adjacency();
    crate::boolean(m, BoolKind::Fuse, s, boss).expect("the boss")
}

#[test]
fn a_blind_bore_can_be_fused_onto_afterwards() {
    let (mut m, s) = plate_with_a_bore(3.0);
    let out = fuse_a_clear_boss(&mut m, s, 5.0);
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert_eq!(
        lateral_face_counts(&m, out[0]),
        vec![1],
        "the bore's wall survives"
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - std::f64::consts::PI * 27.0 + 48.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// The `Cut` twin. ★ The tool is sunk into the plate (`z` from 3): a body merely *resting* on
/// the top face is a coplanar contact that meets `StraightAngle`, a different refusal
/// entirely, and this fixture would then be measuring that instead of the bore.
#[test]
fn a_blind_bore_can_be_cut_afterwards() {
    let (mut m, s) = plate_with_a_bore(3.0);
    let tool = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([4.0, 4.0, 3.0]),
        Point3::from_array([12.0, 6.0, 8.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, s, tool).expect("the cut");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - std::f64::consts::PI * 27.0 - 8.0 * 2.0 * 2.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★★ **The script a person actually writes**: a blind hole, then a second hole somewhere
/// else. It goes through the cylinder-pair rule and two lateral surfaces at once — a road the
/// boss fixtures above do not take.
#[test]
fn a_blind_bore_can_be_drilled_beside() {
    let (mut m, s) = plate_with_a_bore(3.0);
    let second = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([20.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        7.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, s, second).expect("the second hole");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert_eq!(
        lateral_face_counts(&m, out[0]),
        vec![1, 1],
        "two bores, one wall each"
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4000.0 - std::f64::consts::PI * 27.0 - std::f64::consts::PI * 45.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★ **The control**: a through bore and an overshooting one already worked, and must keep
/// working. Without this, a change that merely made *any* second boolean succeed would look
/// like the repair.
#[test]
fn a_through_bore_is_unaffected() {
    for h in [5.0, 7.0] {
        let (mut m, s) = plate_with_a_bore(h);
        let out = fuse_a_clear_boss(&mut m, s, 5.0);
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 45.0 + 48.0;
        assert!((v - want).abs() < 1e-9, "h={h}: {v} vs {want}");
    }
}

/// **An enclosed cylindrical void is a cavity, not a second body.** A `[0,4]³` box `Cut` by a
/// radius-`0.5`, height-`2` cylinder wholly inside it: the boolean's two components are the
/// box's shell and the void's, and the void must come back **depth 1** — odd, so a hole.
///
/// ★★★ **This is the population the coordinate probe exists for, and the only one that makes
/// it say `true`.** The void's boundary is two disks and a band, so it carries no vertex at
/// all and `nodes_of` comes back empty.
/// Two *disjoint* bodies exercise the same road, but their answer is "outside" either way, so
/// a road that always said `false` would pass them; here it would make the void a **second
/// solid** and the box's volume whole.
#[test]
fn an_enclosed_cylindrical_void_is_a_cavity() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0; 3]),
    );
    let void = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 2.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, a, void).expect("a void inside a box");
    assert_eq!(out.len(), 1, "one body, hollow — not two");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert_eq!(m.solid(out[0]).cavities.len(), 1, "the void is a cavity");
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 64.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!(
        (v - want).abs() < 1e-9,
        "the box less the void: {v} vs {want}"
    );
    // Two shells, genus 0 each: `χ = V − E + F − L = 2(S − G) = 4`.
    let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
    assert_eq!(
        v_n - e_n + f_n - l_n,
        4,
        "two genus-0 shells: V{v_n} E{e_n} F{f_n} L{l_n}"
    );
}

/// **A bore and a sealed void in one body.** `[0,8]×[0,4]×[0,4]` with a through bore at
/// `(2,2)` and, `4` away, a wholly enclosed cylinder at `(6,2)` — `128 − π − π/2`.
///
/// ★ It is the only fixture where the **coordinate** probe meets a face with a hole: the
/// void has no vertex, so it is probed from a cap centre, and the body it asks about has
/// annular caps. Two cylinders in one boolean also need `cylinders_clear`, which the `4`
/// between the axes supplies.
#[test]
fn a_bore_and_a_sealed_void_live_in_one_body() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([8.0, 4.0, 4.0]),
    );
    let bore = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 2.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        6.0,
    )
    .solid;
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, a, bore).expect("the bore");
    let void = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([6.0, 2.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, bored[0], void).expect("the sealed void");
    assert_eq!(out.len(), 1, "one body");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert_eq!(m.solid(out[0]).cavities.len(), 1, "the void is a cavity");
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 128.0 - std::f64::consts::PI * 0.25 * 4.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!((v - want).abs() < 1e-9, "the box less both: {v} vs {want}");
}

/// **A cavity finds its owner when there is more than one material.** The void fixture above
/// takes a shortcut the code makes explicit — one material owns every cavity, no search — so
/// the search itself is only reached when two materials stand apart. Here the hollow box is
/// fused with a distant second box, and the void must still come back as the **hollow box's**
/// cavity rather than a third body.
///
/// ★ It is the only fixture that runs `first_deciding` over a **coordinate** probe in the
/// cavity-owner search: the void has no vertex, and every other multi-material case in the
/// suite has cavities with corners.
#[test]
fn a_cavity_with_no_vertex_still_finds_its_owner() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0; 3]),
    );
    let void = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 2.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let hollow = crate::boolean(&mut m, BoolKind::Cut, a, void).expect("a void inside a box");
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([10.0, 10.0, 10.0]),
        Point3::from_array([14.0, 14.0, 14.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, hollow[0], b).expect("two bodies apart");
    assert_eq!(out.len(), 2, "the void is a cavity, not a third body");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let mut seen: Vec<(f64, usize)> = out
        .iter()
        .map(|&s| {
            (
                nacre_props::mass_props(&m, s).expect("props").volume,
                m.solid(s).cavities.len(),
            )
        })
        .collect();
    seen.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite"));
    let hollowed = 64.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!(
        (seen[0].0 - hollowed).abs() < 1e-9 && seen[0].1 == 1,
        "the hollow box keeps its one cavity: {seen:?}"
    );
    assert!(
        (seen[1].0 - 64.0).abs() < 1e-9 && seen[1].1 == 0,
        "the second box is solid: {seen:?}"
    );
}

/// **A drill that misses.** `Cut` returns the box untouched and `Fuse` returns two bodies —
/// the population where a cylinder is present but no circle is ever emitted, so the whole
/// curved path has to stay out of the way.
#[test]
fn a_cylinder_that_misses_changes_nothing_and_fuses_apart() {
    let (mut m, a, b) = {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0; 3]),
        );
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([10.0, 10.0, 0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            1.0,
        )
        .solid;
        m.rebuild_adjacency();
        (m, a, b)
    };
    let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that misses");
    assert_eq!(cut.len(), 1);
    let v = nacre_props::mass_props(&m, cut[0]).expect("props").volume;
    assert!((v - 8.0).abs() < 1e-12, "the box is untouched: {v}");

    // ★★ **The fuse this test's name promises.** Two disjoint bodies means `n = 2`, so each is
    // classified by a ray probe — one that must not ask a band face for its plane class (a
    // panic): the ray counts a cylinder's crossings, and the bare cylinder — whose boundary
    // carries **no vertex at all** — is probed from a cap
    // disk's centre instead of from a corner.
    //
    // ★ The values are derived from the inputs, not read back from the result: a `[0,2]³` box
    // is `8`, and a radius-`0.5`, height-`1` cylinder is `π/4`. Two separate bodies, each a
    // sphere topologically (`χ = 2`), and `validate` clean.
    let mut m2 = Model::new();
    let a2 = crate::fixtures::cuboid(
        &mut m2,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0; 3]),
    );
    let b2 = crate::fixtures::cylinder_with_seam(
        &mut m2,
        Point3::from_array([10.0, 10.0, 0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        1.0,
    )
    .solid;
    m2.rebuild_adjacency();
    let out = crate::boolean(&mut m2, BoolKind::Fuse, a2, b2).expect("two bodies apart");
    assert_eq!(out.len(), 2, "a fuse of two disjoint bodies is two bodies");
    assert!(
        nacre_validate::validate(&m2).is_empty(),
        "{:?}",
        nacre_validate::validate(&m2)
    );
    let mut vols: Vec<f64> = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m2, s).expect("props").volume)
        .collect();
    vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
    let want = [std::f64::consts::PI * 0.25, 8.0];
    assert!(
        (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
        "the cylinder and the box, each untouched: {vols:?}"
    );
    for &s in &out {
        let (v_n, e_n, f_n, l_n) = euler_counts(&m2, s);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            2,
            "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
    }
}

/// **A blind hole, end to end.** The bore stops inside the box, so its own cap closes the
/// bottom — a disk face that only the seated-circle path can emit.
#[test]
fn a_blind_hole_is_built_and_measures_what_it_should() {
    let (mut m, a, b) = box_and_drill(-1.0, 2.0); // z ∈ [−1, 1]
    let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("the blind bore cuts");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 8.0 - std::f64::consts::PI * 0.25 * 1.0; // the bore reaches z = 1
    assert!((v - want).abs() < 1e-9, "volume {v} vs {want}");
}

/// ★ **The uniform-slab theorem's counterexample, as a band case.** A cylinder standing in an
/// L-prism's notch is clear of every wall by more than `r` — the refuted "z-range" rule would
/// keep its wall — but it is outside the material, so `Cut` changes nothing and no band
/// survives.
#[test]
fn a_cylinder_in_the_notch_keeps_no_band() {
    let profile = crate::ops::Profile2d::polygon(vec![
        nacre_math::Point2::from_array([0.0, 0.0]),
        nacre_math::Point2::from_array([2.0, 0.0]),
        nacre_math::Point2::from_array([2.0, 1.0]),
        nacre_math::Point2::from_array([1.0, 1.0]),
        nacre_math::Point2::from_array([1.0, 2.0]),
        nacre_math::Point2::from_array([0.0, 2.0]),
    ])
    .unwrap();
    let mut m = crate::ops::replay(&[crate::tests::extrude_log_op(profile, 1.0)]).unwrap();
    let a = m.live_solids()[0];
    // In the notch (x,y ∈ [1,2]²), a slim drill standing clear of both notch walls.
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([1.5, 1.5, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.25,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    let (out, _) = bands(&m, a, b, BoolKind::Cut);
    assert!(
        out.is_empty(),
        "the notch is void — nothing to cut: {out:?}"
    );

    // ★ …and end to end: `Cut` returns the prism untouched. That is the counterexample's real
    // assertion — a band pass answering by z-range would run a wall through the notch, and
    // this volume would move.
    let before = nacre_props::mass_props(&m, a).expect("props").volume;
    let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that changes nothing");
    assert_eq!(cut.len(), 1);
    let after = nacre_props::mass_props(&m, cut[0]).expect("props").volume;
    assert!(
        (after - before).abs() < 1e-12,
        "the notch drill removes nothing: {before} → {after}"
    );
}
