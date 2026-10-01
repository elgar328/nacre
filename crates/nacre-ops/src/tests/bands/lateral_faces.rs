//! One surface, several lateral faces: the band belongs to the face; the footprint along the
//! cylinder.

use super::*;

/// The shape at the heart of it: a bored plate whose bore's **middle** is cut away, leaving two
/// disjoint bands on one lateral surface. Returns the model and that solid.
fn plate_with_a_split_bore(cut_x0: f64, cut_x1: f64) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([10.0; 3]),
    );
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([5.0, 5.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        12.0,
    )
    .solid;
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    // Walls 3 from the axis (r = 2), so the wall rule is not what this measures.
    let mid = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([cut_x0, cut_x0, 4.0]),
        Point3::from_array([cut_x1, cut_x1, 6.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, bored, mid).expect("the middle cut");
    (m, out[0])
}

/// The lid tool the locks below cut with.
fn add_lid(m: &mut Model) -> Handle<Solid> {
    let lid = crate::fixtures::cuboid(
        m,
        Point3::from_array([0.0, 0.0, 8.0]),
        Point3::from_array([10.0, 10.0, 12.0]),
    );
    m.rebuild_adjacency();
    lid
}

/// **Making it.** The split bore is an ordinary, correct solid — and nothing measured that
/// until this defect was found, which is why it is its own test: if building it ever breaks,
/// the test below (which *uses* it) must not be the one that goes red.
#[test]
fn a_bore_cut_across_the_middle_leaves_two_bands_on_one_surface() {
    let (m, s) = plate_with_a_split_bore(2.0, 8.0);
    assert_eq!(
        lateral_face_counts(&m, s),
        vec![2],
        "one lateral surface, two faces"
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, s).expect("props").volume;
    // 1000 − bore(π·4·10) − the middle cube outside the bore(6·6·2 − π·4·2)
    let want = 1000.0 - std::f64::consts::PI * 40.0 - (72.0 - std::f64::consts::PI * 8.0);
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★★ **Using it.** The two-banded solid as an operand. Letting the first lateral face speak for
/// the class clips the second band's span away — the result shell comes back open
/// (`OpenResultShell`), naming a symptom of our own omission rather than anything about the input.
#[test]
fn a_solid_with_two_bands_on_one_surface_can_be_cut_again() {
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let lid = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 8.0]),
        Point3::from_array([10.0, 10.0, 12.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, s, lid).expect("the lid cut");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    // ★ The property under repair, counted directly — a volume alone does not say "two bands
    // came out", and this defect is exactly one band going missing.
    assert_eq!(
        lateral_face_counts(&m, out[0]),
        vec![2],
        "both bands survive the cut"
    );
    let p = nacre_props::mass_props(&m, out[0]).expect("props");
    let want = 1000.0
        - std::f64::consts::PI * 40.0
        - (72.0 - std::f64::consts::PI * 8.0)
        - (200.0 - std::f64::consts::PI * 8.0);
    assert!((p.volume - want).abs() < 1e-9, "{} vs {want}", p.volume);
    // ★★ **Area is the sharper oracle for this defect.** A missing band does not change the
    // volume — the boolean simply refuses — but it takes its `2πr·h` out of the surface. The
    // two bands here are `z ∈ [0,4]` and `[6,8]`, so `16π + 8π` of lateral area must be in
    // this number: 640 + 8π once every planar face is counted.
    let want_area = 640.0 + std::f64::consts::PI * 8.0;
    assert!(
        (p.area - want_area).abs() < 1e-9,
        "{} vs {want_area}",
        p.area
    );
}

/// ★ **The negative control.** Move the middle cut into a corner, away from the bore: the wall
/// stays **one** face and the same three steps work whatever the lateral rows say. Without this, a
/// change that merely made *any* third boolean succeed would look like a repair.
#[test]
fn a_middle_cut_that_misses_the_bore_leaves_one_band() {
    let (mut m, s) = plate_with_a_split_bore(0.0, 2.0);
    assert_eq!(lateral_face_counts(&m, s), vec![1], "the wall is untouched");
    let lid = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 8.0]),
        Point3::from_array([10.0, 10.0, 12.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, s, lid).expect("the lid cut");
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 1000.0 - std::f64::consts::PI * 40.0 - 8.0 - (200.0 - std::f64::consts::PI * 8.0);
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **The mechanism, read directly.** Two rows, one per lateral face, with the spans the two
/// bands actually occupy — and the `t` axis here starts at the drill's origin `z = −1`, so the
/// bands `z ∈ [0,4]` and `[6,10]` are `t ∈ [1,5]` and `[7,11]`.
#[test]
fn cyl_rows_gives_one_row_per_lateral_face() {
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let lid = add_lid(&mut m);
    let setup = plane_index_setup(&m, s, lid).unwrap();
    let rows = cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).expect("rows");
    let spans: Vec<[f64; 2]> = rows
        .iter()
        .map(|r| [r.span[0].to_f64(), r.span[1].to_f64()])
        .collect();
    assert_eq!(
        spans,
        vec![[1.0, 5.0], [7.0, 11.0]],
        "one row per face, in t order"
    );
    assert!(rows.iter().all(|r| r.class == 0), "both on the one surface");
}

/// **A class the face rows never name is the stages disagreeing, not an input.** The real setup
/// with its lateral rows relabelled to class 1 leaves class 0 with no row — a wiring failure no
/// boolean reaches, planted so the guard's name is seen to bite.
#[test]
fn a_cylinder_class_no_face_row_names_is_a_defect() {
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let lid = add_lid(&mut m);
    let mut setup = plane_index_setup(&m, s, lid).unwrap();
    for ix in &mut setup.plane_ix {
        if *ix == ClassIx::Cyl(0) {
            *ix = ClassIx::Cyl(1);
        }
    }
    let got = cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).err();
    assert!(
        matches!(
            got,
            Some(BoolError::Rejected {
                reason: RejectReason::CylinderStagesDisagree,
                ..
            })
        ),
        "{got:?}"
    );
}

/// ★★ **The band pass makes nothing in the gap.** This is where the right fix parts company
/// with the plausible one: merging the two faces' spans into a single `min..max` would put a
/// band across `z ∈ [4,6]`, where the solid has no lateral face at all.
#[test]
fn no_band_is_invented_where_the_solid_has_no_face() {
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let lid = add_lid(&mut m);
    let (out, ts) = bands(&m, s, lid, BoolKind::Cut);
    let mut spans: Vec<(f64, f64)> = out.iter().map(|lf| ends(lf, &ts)).collect();
    spans.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (lo, hi) in &spans {
        // Disjoint from the open gap (4, 6): a band either ends at or below 4, or starts at or
        // above 6. ★ Measured: with the spans merged into one `min..max` this fires with
        // "a band 4..6 crosses the gap", which is the whole reason the clause is here.
        assert!(
            *hi <= 4.0 || *lo >= 6.0,
            "a band {lo}..{hi} crosses the gap z 4..6 where there is no face"
        );
    }
    assert!(
        spans.iter().any(|(lo, _)| (*lo - 0.0).abs() < 1e-9),
        "the lower band is emitted: {spans:?}"
    );
    assert!(
        spans.iter().any(|(lo, _)| (*lo - 6.0).abs() < 1e-9),
        "the upper band is emitted too: {spans:?}"
    );
}

//
// ★★ A wall parallel to the axis meets the cylinder in a **rectangle** of the wall's own
// plane — the strip across, a lateral face's span along — so "does this face miss it" is one
// question with two separating axes. Reading only the first refuses every wall that stands
// clear of the cylinder along its length.
//
// ★ With several lateral faces there are several rectangles: the strip is shared, the spans
// are not. That is why the gap between two bands is passable at all, and it is not a special
// case bolted on — a row's span is the existence truth at an uncut end, so no band is ever
// built in a gap and there is no premise there for a wall to break.

/// **The case this opened.** The two-banded bore of `plate_with_a_split_bore` has bands at
/// `z ∈ [0,4]` and `[6,10]` — `t ∈ [1,5]` and `[7,11]`. A tool whose `y = 6` wall stands only 1
/// from the axis (`r = 2`) crosses the strip, and its face runs the plate's whole width, so the
/// axis across cannot clear it. Along the axis it is `t ∈ [5.5, 6.5]` — inside the gap, missing
/// both rectangles.
#[test]
fn a_wall_face_in_the_gap_between_two_bands_clears() {
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let tool = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 4.5]),
        Point3::from_array([10.0, 6.0, 5.5]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, s, tool).expect("the wall clears in t");
    assert_eq!(out.len(), 1, "one body");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert_eq!(
        lateral_face_counts(&m, out[0]),
        vec![2],
        "the tool passes between the bands and leaves both"
    );
    // 1000 − bore(π·4·10) − the middle cut(6·6·2 − π·4·2) − this tool, which meets the solid
    // over 10×6 minus the 6×4 the middle cut already took, one deep.
    let want =
        1000.0 - std::f64::consts::PI * 40.0 - (72.0 - std::f64::consts::PI * 8.0) - (60.0 - 24.0);
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★★ **The two axes are one judgement, measured as one.** The same tool three ways against
/// the same two-banded bore, with only the numbers moved: clearing *either* axis is enough,
/// and clearing neither is the **record**, not a refusal. A rule that
/// had merely gained a second, independent test would pass (a) and (b) too — what this pins
/// is (c), that the two are OR-ed rather than each able to wave a face through on its own
/// terms: a face waved through is an unlisted pair the tracer stays silent for, and the bore
/// would be left uncut where the face crosses it — (c)'s exact volume is what says the pair
/// was recorded.
#[test]
fn either_axis_clears_the_footprint_and_neither_does_not() {
    // (a) Across only: the `y = 6` face sits at `x ∈ [0, 2.5]`, clear of the strip
    //     `x ∈ [3.27, 6.73]`, while its plane still crosses the bore. Along the axis it runs
    //     `t ∈ [4,8]`, straddling both bands. ★ The `x = 2.5` wall must stand clear of the axis
    //     by more than `r`, or *it* becomes the face under test — at `x = 3` it is exactly
    //     tangent, and this arm would then measure the tangency instead of the clearance.
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let across = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 3.0]),
        Point3::from_array([2.5, 6.0, 7.0]),
    );
    m.rebuild_adjacency();
    crate::boolean(&mut m, BoolKind::Cut, s, across).expect("clear across the strip");

    // (b) Along only: the face spans the full width, so only the gap saves it.
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let along = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 4.5]),
        Point3::from_array([10.0, 6.0, 5.5]),
    );
    m.rebuild_adjacency();
    crate::boolean(&mut m, BoolKind::Cut, s, along).expect("clear along the axis");

    // (c) Neither: full width *and* straddling both bands — a genuine crossing of both
    // laterals (`y = 6`, 1 from the axis, r = 2), which the gate
    // records: two ruling pieces per band at `x = 5 ± √3`, the chord on the caps and on
    // the middle cut's ceiling and floor. One body, the exact volume.
    let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
    let neither = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 3.0]),
        Point3::from_array([10.0, 6.0, 7.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, s, neither).expect("crosses both, builds");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 1);
    assert!(nacre_validate::validate(&m).is_empty());
    let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let pi = std::f64::consts::PI;
    // The plate less the bore, less the middle cut (less the bore inside it), less the tool's
    // box (less the bore's `y ≤ 6` part over its height, less the middle cut inside it — which
    // had already lost the bore's `y ≤ 6` part).
    let disk_le6 = 4.0 * pi - seg(1.0, 2.0);
    let want =
        1000.0 - 40.0 * pi - (72.0 - 8.0 * pi) - (240.0 - 4.0 * disk_le6 - (48.0 - 2.0 * disk_le6));
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}
