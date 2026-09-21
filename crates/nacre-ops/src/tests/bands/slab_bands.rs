//! Slab bands: what a band knows, where it ends, a drilled box, which bands each boolean keeps.

use super::*;

/// **A through hole.** The drill runs from z=−1 to z=3 through a `[0,2]³` box; `Cut` keeps
/// exactly the wall between the box's two caps, and the material there is **outside** the
/// cylinder, so the face is flipped against the stored (outward) normal.
#[test]
fn a_through_hole_keeps_the_band_between_the_caps_and_flips_it() {
    let (m, a, b) = box_and_drill(-1.0, 4.0);
    let (out, ts) = bands(&m, a, b, BoolKind::Cut);
    assert_eq!(out.len(), 1, "one surviving band: {out:?}");
    assert_eq!(ends(&out[0], &ts), (0.0, 2.0));
    assert!(out[0].flip, "a hole's wall faces its own axis");
    assert!(matches!(out[0].surf, ClassIx::Cyl(0)));
}

/// The same drill, fused instead: the material is **inside** the cylinder exactly where the
/// box is not, so the two protruding stretches survive and the buried one does not — and
/// neither survivor is flipped.
#[test]
fn a_fused_drill_keeps_the_two_protruding_bands_unflipped() {
    let (m, a, b) = box_and_drill(-1.0, 4.0);
    let (out, ts) = bands(&m, a, b, BoolKind::Fuse);
    let mut spans: Vec<(f64, f64)> = out.iter().map(|lf| ends(lf, &ts)).collect();
    spans.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    assert_eq!(spans, vec![(-1.0, 0.0), (2.0, 3.0)], "{out:?}");
    assert!(out.iter().all(|lf| !lf.flip), "a boss's wall faces outward");
}

/// `Common` keeps the buried stretch and nothing else — the mirror of the fuse case, and the
/// one that would pass if the keep rule ignored the operand order.
#[test]
fn a_common_keeps_only_the_buried_band() {
    let (m, a, b) = box_and_drill(-1.0, 4.0);
    let (out, ts) = bands(&m, a, b, BoolKind::Common);
    assert_eq!(out.len(), 1, "{out:?}");
    assert_eq!(ends(&out[0], &ts), (0.0, 2.0));
    assert!(!out[0].flip, "the kept material is inside the wall");
}

/// **A blind hole.** The drill stops at z=1 inside the box: the wall survives from the box's
/// cap down to the drill's own cap, which is where the bottom disk closes it.
#[test]
fn a_blind_hole_stops_at_the_drills_own_cap() {
    let (m, a, b) = box_and_drill(-1.0, 2.0); // z ∈ [−1, 1]
    let (out, ts) = bands(&m, a, b, BoolKind::Cut);
    assert_eq!(out.len(), 1, "{out:?}");
    assert_eq!(ends(&out[0], &ts), (0.0, 1.0));
    assert!(out[0].flip);
}

/// ★★ **The first end-to-end cylinder boolean.** A `[0,2]³` box drilled through by a
/// radius-0.5 bore: seven faces (four walls, two drilled caps, one hole wall), genus 1, and a
/// volume of `8 − π·0.25·2`. The volume is the **winding lock** — props integrates by the
/// divergence theorem, so a hole loop wound the wrong way returns `8 + πr²h` instead.
#[test]
fn a_through_hole_is_built_and_measures_what_it_should() {
    let (mut m, a, b) = box_and_drill(-1.0, 4.0);
    let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("the drill cuts");
    assert_eq!(out.len(), 1, "one body");
    let s = out[0];
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let faces = crate::planes::solid_shell_handles(&m, s)
        .into_iter()
        .map(|sh| m.shell(sh).faces.len())
        .sum::<usize>();
    assert_eq!(faces, 7, "4 walls + 2 drilled caps + the bore's wall");
    // ★ **Genus 1, re-derived rather than quoted.** An earlier plan wrote `V10−E15+F7−L2`
    // from memory; the counts are measured here and what is asserted is the relation
    // (`χ = V − E + F − L = 2(S − G)`, so one shell with one through hole gives `χ = 0`).
    let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
    assert_eq!(
        v_n - e_n + f_n - l_n,
        0,
        "genus 1: V{v_n} E{e_n} F{f_n} L{l_n}"
    );
    let v = nacre_props::mass_props(&m, s)
        .expect("a closed solid has mass props")
        .volume;
    let want = 8.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!(
        (v - want).abs() < 1e-6,
        "volume {v} vs {want} — a bore removes material, it does not add it"
    );
    // The same box without the bore, for the sign of the correction rather than its value.
    let mut m2 = Model::new();
    let plain = m2.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    m2.rebuild_adjacency();
    let plain_v = nacre_props::mass_props(&m2, plain)
        .expect("a box has mass props")
        .volume;
    assert!(v < plain_v, "the bore removes material: {v} vs {plain_v}");
}

/// ★★ **Two bores in one plate — the case the band pass was rebuilt for.**
///
/// The second cut's counterpart already carries a cylinder face (the first bore's wall), and
/// the old witness road could not describe such a body, so this used to decline by name. The
/// band pass now reads the arrangement's own disk labels, and a curved counterpart is no
/// longer a question anyone has to answer.
///
/// ★ It is also where the wall's **own** membership stopped being assumed: the first bore's
/// wall bounds the *plate*, and inside that circle the plate has no material — the opposite of
/// a drill, which fills its own cylinder. Assuming the drill's case (as the pass did while
/// only drills were tested) puts the second bore's `keep` on the wrong chamber.
#[test]
fn a_two_hole_plate_drills_both() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 2.0, 1.0]),
    );
    let d1 = m.add_cylinder(
        Point3::from_array([1.0, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.25,
        3.0,
    );
    let d2 = m.add_cylinder(
        Point3::from_array([3.0, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.25,
        3.0,
    );
    m.rebuild_adjacency();
    let one = crate::boolean(&mut m, BoolKind::Cut, a, d1).expect("first bore");
    let v1 = nacre_props::mass_props(&m, one[0]).expect("props").volume;
    assert!(
        (v1 - (8.0 - std::f64::consts::PI * 0.0625)).abs() < 1e-9,
        "one bore: {v1}"
    );
    let two = crate::boolean(&mut m, BoolKind::Cut, one[0], d2).expect("second bore");
    assert_eq!(two.len(), 1, "one body");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v2 = nacre_props::mass_props(&m, two[0]).expect("props").volume;
    let want = 8.0 - 2.0 * std::f64::consts::PI * 0.0625;
    assert!((v2 - want).abs() < 1e-9, "two bores: {v2} vs {want}");
    // Genus 2: two through holes in one shell.
    let reach = m.reachable();
    let l: i64 = reach
        .faces
        .iter()
        .map(|fh| m.face(*fh).inner.len() as i64)
        .sum();
    let chi = reach.vertices.len() as i64 - reach.edges.len() as i64 + reach.faces.len() as i64 - l;
    assert_eq!(chi, -2, "two handles: χ = 2(1 − 2)");
}

/// **Three bores.** Two was the case the label route was built for; three is the check that
/// nothing in it counts to two — each cut's counterpart carries one more cylinder face than
/// the last.
#[test]
fn three_bores_in_a_row() {
    let mut m = Model::new();
    let mut solid = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([6.0, 2.0, 1.0]),
    );
    for x in [1.0, 3.0, 5.0] {
        let d = m.add_cylinder(
            Point3::from_array([x, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        solid = crate::boolean(&mut m, BoolKind::Cut, solid, d).expect("a bore")[0];
    }
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).expect("props").volume;
    let want = 12.0 - 3.0 * std::f64::consts::PI * 0.0625;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **A drilled plate cut by a body flush with its faces.** `Common` of a drilled plate with a
/// box covering half of it: the two operands share the plate's top and bottom planes, and the
/// bore's rims end on *holed* faces rather than on disks.
#[test]
fn a_drilled_plate_can_be_cut_by_a_body_flush_with_its_faces() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 2.0, 1.0]),
    );
    let d = m.add_cylinder(
        Point3::from_array([1.0, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.25,
        3.0,
    );
    m.rebuild_adjacency();
    let drilled = crate::boolean(&mut m, BoolKind::Cut, plate, d).expect("bore")[0];
    let half = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Common, drilled, half).expect("flush common");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 4.0 - std::f64::consts::PI * 0.0625;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **A cylinder standing far away, whose caps happen to land on the box's own planes.** Two
/// coplanar planes are one class whichever way they sit, so this cylinder's *disks* share a
/// class with the box's faces while touching nothing — the shape that made the retired seated
/// rule mistake a shared class for a contact. `Cut` takes nothing away.
#[test]
fn a_distant_cap_on_the_boxs_own_plane_removes_nothing() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    // Far from the box in x/y, but its caps land exactly on z = 0 and z = 2.
    let b = m.add_cylinder(
        Point3::from_array([10.0, 10.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that misses");
    assert_eq!(out.len(), 1);
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    assert!((v - 8.0).abs() < 1e-12, "the box is untouched: {v}");
}
