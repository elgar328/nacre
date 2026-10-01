//! Seated caps: a cap flush on the other body's face.

use super::*;

//
// ★★ These five are the seated population. What actually makes a
// seated circle hard is its boundary meeting the counterpart's — and a boundary is either an
// edge on a plane (parallel to the axis → the wall rule, oblique → asked to miss every lateral
// face, else `ObliqueCylinderCut`) or another cylinder's rim (→ proved apart per
// face pair, else `CylinderPairContact`). ★ The wall rule is not a *fence*:
// it records a crossing or a tangency and the roads behind it answer, so
// what still stands at the end of this block is the oblique cut that does meet a face and the
// cylinder pair whose faces do meet.
//
// ★★★ **Measured against the pre-deletion kernel, all seven of these came back
// `SeatedCylinderCap` — the two fences included.** The rule stood before the wall and the
// curved-depth rules and answered for them: an overhanging boss and a pair of coplanar-capped
// cylinders are refused for reasons that have nothing to do with seating, and the sentence a
// user got named the seating anyway. Deleting it does not only open the five; it lets the two
// that stay shut say what is actually in the way.

/// **(a) A through hole whose two caps are flush with the plate's own faces.** The drill neither
/// overshoots nor stops short: both cap planes coincide with the plate's, and every ⊥ class in
/// the operation carries faces of both operands.
#[test]
fn a_flush_through_drill_bores_the_plate() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let drill = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 2.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a flush through hole");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 32.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **(b) A blind hole seated on the face it is drilled from.** Only the lower cap is flush;
/// the upper one stops inside the plate and closes the bore itself.
#[test]
fn a_blind_drill_seated_on_the_plates_base_bores_it() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let drill = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 2.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        1.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a seated blind hole");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 32.0 - std::f64::consts::PI * 0.25;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **(c) A boss standing on the plate.** `Fuse` with the cylinder's base cap flush on the
/// plate's top face — the seating a person draws first.
#[test]
fn a_boss_seated_on_the_plate_fuses_to_it() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 2.0, 2.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        1.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("a seated boss");
    assert_eq!(out.len(), 1, "one body, not two");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 32.0 + std::f64::consts::PI * 0.25;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **(d) `Common` where the cylinder's caps sit on the box's own faces.** The intersection is
/// the whole cylinder — every one of its faces is seated or shared.
#[test]
fn a_flush_cylinder_meets_the_box_in_itself() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0; 3]),
    );
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([1.0, 1.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        2.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Common, a, b).expect("a flush common");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = std::f64::consts::PI * 0.25 * 2.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **(i) Drilling the floor of a pocket.** The seated face here is not an original operand
/// face at all — it is a *floor the previous boolean made*, so the seating arrives through the
/// arrangement rather than from the modeller. A blind bore down from that floor.
#[test]
fn a_drill_seated_on_a_pocket_floor_bores_it() {
    let mut m = Model::new();
    let block = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([6.0; 3]),
    );
    // Open at the top — a tool that stopped inside would leave a void, not a pocket.
    let pocket_tool = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, 2.0]),
        Point3::from_array([5.0, 5.0, 7.0]),
    );
    m.rebuild_adjacency();
    let pocketed = crate::boolean(&mut m, BoolKind::Cut, block, pocket_tool).expect("pocket")[0];
    // The bore hangs from the pocket floor (z = 2) down into the material below it.
    let drill = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([3.0, 3.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        1.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, pocketed, drill).expect("a floor bore");
    assert_eq!(out.len(), 1);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 216.0 - 64.0 - std::f64::consts::PI * 0.25;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// One rulings-road boolean on the 40×40×20 plate, end to end on the production road:
/// build, retire the operands, validate clean, answer the exact volume, tessellate
/// watertight. The volume oracle is what finally *measures* the ladder's sign roster —
/// a flipped panel winding, disk-side selector, ruling turn or chord sense moves it.
fn through_boss_builds(kind: BoolKind, base: [f64; 3], want: f64) {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 20.0]),
    );
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array(base),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        5.0,
        50.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, kind, plate, boss).expect("the rulings road builds");
    assert_eq!(out.len(), 1, "one solid");
    m.rebuild_adjacency();
    assert_eq!(
        m.live_solids().to_vec(),
        out,
        "the operands retired, the result lives"
    );
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
        .expect("the panels tessellate");
    let mut uses: std::collections::HashMap<(u32, u32), usize> = std::collections::HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    assert_eq!(
        uses.values().filter(|&&n| n != 2).count(),
        0,
        "the mesh is watertight"
    );
}

/// ★ **The through-boss builds** — the rulings ladder's milestone: a boss standing through
/// the plate's wall (axis exactly on the `x = 40` plane). Fuse keeps the outer half
/// (plate + the boss outside it), cut
/// carves the notch, common keeps the inner half-cylinder — each an exact closed form, and
/// the three exercise both complementary θ-sectors.
#[test]
fn a_through_boss_fuses() {
    let pi = std::f64::consts::PI;
    through_boss_builds(BoolKind::Fuse, [40.0, 20.0, -10.0], 32000.0 + 1000.0 * pi);
}

#[test]
fn a_through_boss_cuts_a_notch() {
    let pi = std::f64::consts::PI;
    through_boss_builds(BoolKind::Cut, [40.0, 20.0, -10.0], 32000.0 - 250.0 * pi);
}

#[test]
fn a_through_boss_common_is_the_inner_half() {
    let pi = std::f64::consts::PI;
    through_boss_builds(BoolKind::Common, [40.0, 20.0, -10.0], 250.0 * pi);
}

/// ★ A boss on the plate's **corner** — its axis on *two* wall planes, both recorded, each
/// circle cut by both walls: the rim pairing (`(wall class, root)`) and the θ-panel
/// machinery quarter the lateral. Measured to assemble in the lift probe; the volume pins
/// it (plate + three quarters of the cylinder outside).
#[test]
fn a_corner_boss_fuses() {
    let pi = std::f64::consts::PI;
    through_boss_builds(BoolKind::Fuse, [40.0, 40.0, -10.0], 32000.0 + 1125.0 * pi);
}

/// ★★★★★ **Neither exclusion of this fence stands.** The **offset** crossing
/// (`0 <` distance `< r`) builds exactly.
/// The **tangent** wall (distance exactly `r`) «assembles
/// a volume-correct zero-thickness pinch `validate` cannot see» — and that sentence
/// is *exactly right*: with the gate passing it, `Cut` returns `Ok`, `validate` is clean, and
/// nothing in the kernel sees the contact. So the answer is not to keep refusing at the
/// gate; it is to give that pinch a judge (`assembly::tangency_reject`), which is what this row
/// exercises — the operation decides, and the ones that do not pinch build.
///
/// ★ The **half-height** boss, whose upper cap
/// sits inside the plate's material, has its own row: the chart reads it (a band below the
/// plate, the outer sector beside it), so it builds — [`Self::a_half_height_boss_builds`].
#[test]
fn one_chord_formula_covers_the_offset_wall_and_its_tangent_limit() {
    let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let pi = std::f64::consts::PI;
    // ★★★★★ **One expression, both rows — which is the evidence that the rule is general.**
    // The boss's part inside the plate is the disk less the `x > 40` circular segment, over the
    // plate's height; the tangent row is that same expression at `d = r`, where `seg(5, 5) = 0`
    // and the whole disk is inside. Nothing here is copied from an engine run: the tangent
    // volume is *derived* by walking the offset row's formula to its limit.
    let want = |d: f64| 32000.0 + 1250.0 * pi - (25.0 * pi - seg(d, 5.0)) * 20.0;
    for (base, d) in [([38.0, 20.0, -10.0], 2.0), ([35.0, 20.0, -10.0], 5.0)] {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            5.0,
            50.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
            .unwrap_or_else(|e| panic!("{base:?}: {e:?}"));
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1, "{base:?}");
        assert!(nacre_validate::validate(&m).is_empty(), "{base:?}");
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - want(d)).abs() < 1e-9, "{base:?}: {v} vs {}", want(d));
    }
}

/// ★★★★★ **The half-height boss builds.** A boss
/// through the plate's wall whose upper cap (z = 10) sits inside the plate: below the plate
/// the lateral is a whole band, beside it only the outer sector survives, and the cap is a
/// half-disk. The band road refused this end (`RulingBoundNotYet`); the chart's cells read
/// it off the same labels. Volumes derived, not copied: the boss outside the plate is the
/// outer half-cylinder over the full height (`π·25·20/2 = 250π`) plus the inner half below
/// the plate (`π·25·10/2 = 125π`); the boss inside the plate is the inner half over
/// `z ∈ [0, 10]` (`125π`).
///
/// The lateral is **one** face: the band and the outer panel merge into a band whose
/// upper rim is a **wrapping chain** (the z = 0 inner arc, a ruling, the z = 10 outer arc,
/// a ruling). ★ Both
/// senses: the boss with its cap inside the plate (`z0 = −10`, a `hi` chain) and its mirror
/// with its base inside (`z0 = 10`, a `lo` chain) — the walk's rotation rule is measured
/// on each rather than assumed symmetric. Same volumes by symmetry.
#[test]
fn a_half_height_boss_builds() {
    let pi = std::f64::consts::PI;
    for (z0, kind, want, lateral) in [
        (-10.0, BoolKind::Fuse, 32000.0 + 375.0 * pi, vec![1]),
        (-10.0, BoolKind::Cut, 32000.0 - 125.0 * pi, vec![1]),
        (-10.0, BoolKind::Common, 125.0 * pi, vec![1]),
        (10.0, BoolKind::Fuse, 32000.0 + 375.0 * pi, vec![1]),
        (10.0, BoolKind::Cut, 32000.0 - 125.0 * pi, vec![1]),
        (10.0, BoolKind::Common, 125.0 * pi, vec![1]),
    ] {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([40.0, 20.0, z0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            5.0,
            20.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, kind, plate, boss)
            .unwrap_or_else(|e| panic!("{kind:?} z0 {z0}: the half-height boss builds: {e:?}"));
        assert_eq!(out.len(), 1, "{kind:?}: one solid");
        m.rebuild_adjacency();
        assert_eq!(
            m.live_solids().to_vec(),
            out,
            "{kind:?}: the operands retired"
        );
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{kind:?}: {issues:?}");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - want).abs() <= 1e-9 * want, "{kind:?}: {v} vs {want}");
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            lateral,
            "{kind:?}: lateral faces"
        );
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the half-height boss tessellates");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "{kind:?}: the mesh is watertight"
        );
    }
}

/// **The straddling boss builds.** A boss hanging over the
/// plate's edge — its rim circle cut by the plate top's boundary segment — fuses into one
/// valid solid.
///
/// The volume and the mixed-loop integrals are pinned by
/// [`Self::a_cut_rim_boolean_builds_a_complete_solid`]; here the result answers the two
/// whole-model judges: `validate` (whose winding check reads the arcs as witnesses now) and
/// the tessellation (watertight, arcs sampled as sub-arcs).
#[test]
fn a_boss_overhanging_the_plates_edge_builds() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([4.0, 2.0, 2.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        1.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
    assert_eq!(out.len(), 1, "one fused solid");
    m.rebuild_adjacency();
    assert_eq!(
        m.live_solids().to_vec(),
        out,
        "the operands retired, the result lives"
    );
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
        .expect("the arcs tessellate");
    let mut uses: std::collections::HashMap<(u32, u32), usize> = std::collections::HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let open = uses.values().filter(|&&n| n != 2).count();
    assert_eq!(open, 0, "the mesh is watertight");
}

/// ★★ **A boss standing over a bore, its whole outline *inside* the rim.** No segment crosses
/// the circle, so the circle keeps its closed cell — what the shape really asks is that the
/// **disk cell host a polygon**, which `nest_cells` does.
///
/// ★★★ **The population immediately found a live defect in the nesting predicate.** With the
/// disk allowed to host, the *circle*'s contour came back as a hole of the **square** — the
/// footprint `[7,9]×[9,11]` straddles the axis `(8,10)`, so "the circle's centre is inside the
/// ring" is true in the nesting that does not hold. One witness point says *whether* two
/// disjoint loops nest, never *which way*; the other direction is what settles it
/// (`cell_in_cell`).
#[test]
fn a_segment_inside_the_rim_builds_the_boss_over_the_hole() {
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
    // Footprint `[7,9] × [9,11]` sits wholly inside the rim `(8,10)`, `r = 3`; `z ∈ [5,8]`
    // leaves the bore's span `t ∈ [1,6]` clear along the axis, so the wall rule passes it on.
    let boss = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([7.0, 9.0, 5.0]),
        Point3::from_array([9.0, 11.0, 8.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("a boss over a bore");
    // ★ The boss's whole footprint is inside the rim, so it stands over the **hole** and
    // touches no material: two bodies, not one. That is forced by the *geometry* — the
    // footprint is strictly inside the rim, so there is no material to reach — and not
    // by a refusal: the wall rule
    // serves that population instead of fencing it off.
    assert_eq!(out.len(), 2, "the boss touches nothing");
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
    let want = [
        2.0 * 2.0 * 3.0,
        40.0 * 20.0 * 5.0 - std::f64::consts::PI * 9.0 * 5.0,
    ];
    assert!(
        (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
        "the boss and the bored plate: {vols:?} vs {want:?}"
    );
    // The plate keeps its through hole (χ = 0), the boss is a box (χ = 2).
    let mut chi: Vec<i64> = out
        .iter()
        .map(|&s| {
            let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
            v_n - e_n + f_n - l_n
        })
        .collect();
    chi.sort_unstable();
    assert_eq!(chi, vec![0, 2], "genus 1 and genus 0");
}
