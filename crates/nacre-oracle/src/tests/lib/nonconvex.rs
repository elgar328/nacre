//! Non-convex operands against OCCT: L and U prisms, slots, notches, overhangs, staples, dimples.

use super::*;

/// The box bites the L's convex corner `(2, 0)`: a single chord, monotone bends.
/// The non-convex overlap fixture of M5-d2.
fn l_and_corner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    l_prism_and_box([1.3, -0.3, 0.2], [2.4, 0.4, 1.4])
}

/// The box straddles the L's *reflex* corner `(1, 1)`: one chord with one reflex
/// bend, and that bend projects outside its chord's endpoints.
fn l_and_reflex_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    l_prism_and_box([0.6, 0.6, 0.2], [1.6, 1.6, 1.4])
}

/// The box crosses the reflex corner and pops out the L's top: on its bottom face
/// the seam is a staircase whose two bends turn opposite ways.
/// The reconstructed face is a correct simple polygon.
fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    l_prism_and_box([0.5, 0.5, 0.2], [2.5, 1.5, 1.2])
}

/// nacre's non-convex single-chord overlap vs OCCT, `Cut`. The seam is a real
/// boundary crossing, so this scores the boolean's exact classification and seam
/// reconstruction against an independent kernel.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn nonconvex_overlap_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_corner_box();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "non-convex overlap cut: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// The `Fuse` counterpart — the box protrudes past the L, so the union grows.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn nonconvex_overlap_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_corner_box();
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "non-convex overlap fuse: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// The `Common` counterpart — a non-convex `Common` (the intersection is
/// the corner bite, `0.224`). Same seam as the Cut/Fuse
/// above, only the keep/flip table differs.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn nonconvex_overlap_common_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_corner_box();
    let occt = occt_boolean_of(&m, OcctBool::Common, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Common, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "vol {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs {}",
        nacre.area,
        occt.area
    );
}

/// A slotted bar — a cuboid with a full-width groove — carries two coplanar top strips
/// that share one Surface. Cell coplanar-narrow lets it chain a second boolean; OCCT
/// scores the blind pocket cut into one strip on the exported b-rep.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slotted_bar_pocket_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let bar = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );
    let groove = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, -0.5, 0.5]),
        Point3::from_array([2.0, 1.5, 1.5]),
    );
    let slotted = boolean_one(&mut m, BoolKind::Cut, bar, groove).unwrap();
    m.rebuild_adjacency();
    let pocket = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.2, 0.2, 0.7]),
        Point3::from_array([0.5, 0.5, 1.5]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Cut, slotted, pocket).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, slotted, pocket).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "slotted bar pocket cut: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// A blind pocket cut into a face (cell coplanar-contact-cut): a prism inside the base with
/// its top flush is carved out. OCCT scores the pocket and a transversal Cut chained onto
/// the (non-convex) pocketed solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_cut_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let prism = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.25, 0.25, 0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let occt_pocket = occt_boolean_of(&m, OcctBool::Cut, base, prism).unwrap();
    let pocketed = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, pocketed).unwrap().volume, occt_pocket.volume),
        "pocket cut: {} vs {}",
        mass_props(&m, pocketed).unwrap().volume,
        occt_pocket.volume
    );
    // Chain a transversal cut at a corner, away from the pocket and off its face planes
    // (the pocketed solid is non-convex → seam path).
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.8, 0.8, 0.3]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, pocketed, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, pocketed, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "pocket then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A boss fuses onto a face inside its boundary (cell coplanar-contact-boss): the base's
/// top gains the boss footprint as a hole and the boss rides on it. OCCT scores the fuse
/// and a Cut chained onto the bossed solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn boss_fuse_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.25, 0.25, 1.0]),
        Point3::from_array([0.75, 0.75, 2.0]),
    );
    let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, boss).unwrap();
    let bossed = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, bossed).unwrap().volume, occt_fuse.volume),
        "boss fuse: {} vs {}",
        mass_props(&m, bossed).unwrap().volume,
        occt_fuse.volume
    );
    // Chain a cut through the boss.
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.4, 0.4, 0.5]),
        Point3::from_array([0.6, 0.6, 2.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, bossed, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, bossed, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "boss then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A boss overhangs a single edge of a face (cell coplanar-contact-overhang): part fuses,
/// part cantilevers. OCCT scores the fuse and a transversal Cut chained onto the (non-convex)
/// overhanging solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn overhang_fuse_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.25, 1.0]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, boss).unwrap();
    let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, overhung).unwrap().volume, occt_fuse.volume),
        "overhang fuse: {} vs {}",
        mass_props(&m, overhung).unwrap().volume,
        occt_fuse.volume
    );
    // Chain a transversal cut drilling straight through the cantilever (the overhung solid
    // is non-convex → seam path).
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.1, 0.35, 0.5]),
        Point3::from_array([1.4, 0.65, 2.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, overhung, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, overhung, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "overhang then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A single-edge overhang Cut carves an edge-slot breaking out through a wall (cell
/// coplanar-contact-overhang-cut). OCCT scores the slot and a transversal Cut chained onto
/// the (non-convex) slotted solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn edge_slot_cut_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let prism = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.25, 0.5]),
        Point3::from_array([1.5, 0.75, 1.0]),
    );
    let occt_slot = occt_boolean_of(&m, OcctBool::Cut, base, prism).unwrap();
    let slotted = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, slotted).unwrap().volume, occt_slot.volume),
        "edge slot: {} vs {}",
        mass_props(&m, slotted).unwrap().volume,
        occt_slot.volume
    );
    // Chain a transversal cut drilling through the base away from the slot (the slotted
    // solid is non-convex → seam path).
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.15, 0.8, -0.5]),
        Point3::from_array([0.35, 0.95, 1.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, slotted, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, slotted, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "slot then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A corner-overhanging boss fuses onto a face (cell coplanar-contact-overhang-corner): the
/// boss footprint swallows a base-top corner, crossing two edges. OCCT scores the fuse and a
/// transversal Cut chained onto the (non-convex, L-cantilever) result.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn corner_overhang_fuse_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, boss).unwrap();
    let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, overhung).unwrap().volume, occt_fuse.volume),
        "corner fuse: {} vs {}",
        mass_props(&m, overhung).unwrap().volume,
        occt_fuse.volume
    );
    // Chain a transversal cut drilling through the L cantilever's outer corner (x>1, y>1).
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.1, 1.1, 0.5]),
        Point3::from_array([1.4, 1.4, 2.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, overhung, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, overhung, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "corner then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A spanning-slab boss fuses onto a face (cell coplanar-contact-overhang-multi): the slab
/// crosses the base top, splitting each contact face into two pieces. OCCT scores the fuse and
/// a transversal Cut chained onto the (non-convex, two-cantilever) result.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn spanning_slab_fuse_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let slab = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5, 0.4, 1.0]),
        Point3::from_array([1.5, 0.6, 2.0]),
    );
    let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, slab).unwrap();
    let slabbed = boolean_one(&mut m, BoolKind::Fuse, base, slab).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, slabbed).unwrap().volume, occt_fuse.volume),
        "slab fuse: {} vs {}",
        mass_props(&m, slabbed).unwrap().volume,
        occt_fuse.volume
    );
    // Chain a transversal cut drilling through the x>1 cantilever piece.
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.1, 0.45, 0.5]),
        Point3::from_array([1.4, 0.55, 2.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, slabbed, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, slabbed, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "slab then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A spanning-slab Cut carves a channel breaking out two opposite walls (cell
/// coplanar-contact-overhang-slab-cut). OCCT scores the channel and a transversal Cut chained
/// onto the (non-convex) channelled solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_channel_cut_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let slab = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5, 0.4, 0.5]),
        Point3::from_array([1.5, 0.6, 1.0]),
    );
    let occt_channel = occt_boolean_of(&m, OcctBool::Cut, base, slab).unwrap();
    let channelled = boolean_one(&mut m, BoolKind::Cut, base, slab).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(
            mass_props(&m, channelled).unwrap().volume,
            occt_channel.volume
        ),
        "slab channel: {} vs {}",
        mass_props(&m, channelled).unwrap().volume,
        occt_channel.volume
    );
    // Chain a transversal drill through the base away from the channel (y > 0.6).
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.2, 0.75, -0.5]),
        Point3::from_array([0.4, 0.95, 1.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, channelled, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, channelled, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "channel then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A corner slot Cut breaks out two adjacent walls (cell coplanar-contact-overhang-corner-cut).
/// OCCT scores the corner slot and a transversal Cut chained onto the (non-convex) result.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn corner_slot_cut_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let corner = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.0]),
    );
    let occt_slot = occt_boolean_of(&m, OcctBool::Cut, base, corner).unwrap();
    let slotted = boolean_one(&mut m, BoolKind::Cut, base, corner).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, slotted).unwrap().volume, occt_slot.volume),
        "corner slot: {} vs {}",
        mass_props(&m, slotted).unwrap().volume,
        occt_slot.volume
    );
    // Chain a transversal drill through the base at the opposite (0,0) corner.
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.05, 0.05, -0.5]),
        Point3::from_array([0.25, 0.25, 1.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, slotted, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, slotted, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "corner then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// An L-step Cut breaks out three walls at once — two corners (x=0, x=1) and a fully covered
/// wall (y=1) — via the general N-wall path (cell coplanar-contact-overhang-cut-general). OCCT
/// scores the L rebate and a transversal Cut chained onto the (non-convex) result.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn l_step_cut_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let prism = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.0]),
    );
    let occt_step = occt_boolean_of(&m, OcctBool::Cut, base, prism).unwrap();
    let stepped = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, stepped).unwrap().volume, occt_step.volume),
        "l step: {} vs {}",
        mass_props(&m, stepped).unwrap().volume,
        occt_step.volume
    );
    // Chain a transversal drill through the full-height front strip (y < 0.5).
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.1, 0.1, -0.5]),
        Point3::from_array([0.3, 0.3, 1.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, stepped, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, stepped, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "l step then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// A `Common` of a top-flush overhang pair is the convex overlap R = a ∩ b (cell
/// coplanar-contact-overhang-common). The prism hangs past one base edge; OCCT scores R.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn edge_overhang_common_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let prism = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.3, 0.5, 0.5]),
        Point3::from_array([0.7, 1.5, 1.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Common, base, prism).unwrap();
    let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "edge overhang common: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
}

/// A `Common` whose small box is contained in the base's top face (same-normal coplanar cap) yet
/// pokes out the bottom. OCCT cross-checks the result (vol 0.25).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn contained_common_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let box_ = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.25, 0.25, -0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Common, base, box_).unwrap();
    let r = boolean_one(&mut m, BoolKind::Common, base, box_).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "contained common: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
}

/// A corner overhang `Common`: the prism swallows the base's (1,1) corner, so R meets at a
/// corner column (cell coplanar-contact-overhang-common). OCCT cross-checks that cc topology.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn corner_overhang_common_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let prism = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Common, base, prism).unwrap();
    let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "corner overhang common: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
}

/// A blind pocket carved by a non-convex (L-shaped) cutter — the contained-coplanar Cut now
/// admits non-convex operands. OCCT scores the L pocket.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn non_convex_profile_pocket_matches_occt() {
    use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, -1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 0.5]),
    );
    let l = Profile2d::polygon(
        [
            [-0.3, -0.3],
            [0.3, -0.3],
            [0.3, 0.0],
            [0.0, 0.0],
            [0.0, 0.3],
            [-0.3, 0.3],
        ]
        .iter()
        .map(|&p| nacre_math::Point2::from_array(p))
        .collect(),
    )
    .unwrap();
    let __w7 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: lp, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w7,
            profile: l,
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let occt = occt_boolean_of(&m, OcctBool::Cut, base, lp).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, base, lp).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "non-convex profile pocket: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
}

/// A blind pocket carved into an already-pocketed (non-convex) cube. OCCT scores the second
/// pocket cut on the concave kept solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_into_non_convex_solid_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, pc) = nacre_ops::fixtures::pocketed_cube();
    let corner = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.05, 0.1, 0.6]),
        Point3::from_array([0.25, 0.2, 1.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Cut, pc, corner).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, pc, corner).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "pocket into non-convex solid: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
    // The box stands in material beside the void (`x < 0.3`): `0.92 − 0.2·0.1·0.4`.
    let v = mass_props(&m, r).unwrap().volume;
    assert!(
        approx(v, 0.912),
        "pocket into non-convex solid: {v} vs hand 0.912"
    );
}

/// A boss raised by a non-convex (L-shaped) prism — the contained-coplanar Fuse now admits
/// non-convex operands. The base's top sits at z = 0, the L boss extrudes onto it flush.
/// OCCT scores the L boss.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn non_convex_profile_boss_matches_occt() {
    use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, -1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 0.0]),
    );
    let l = Profile2d::polygon(
        [
            [-0.3, -0.3],
            [0.3, -0.3],
            [0.3, 0.0],
            [0.0, 0.0],
            [0.0, 0.3],
            [-0.3, 0.3],
        ]
        .iter()
        .map(|&p| nacre_math::Point2::from_array(p))
        .collect(),
    )
    .unwrap();
    let __w5 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: lb, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w5,
            profile: l,
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let occt = occt_boolean_of(&m, OcctBool::Fuse, base, lb).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, base, lb).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "non-convex profile boss: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
}

/// A boss cantilevers off the side face of a top-pocketed cube (a non-convex solid). The
/// overhang Fuse now admits a non-convex solid when the contact face is convex; OCCT scores the
/// cantilever on the concave part.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn overhang_boss_on_non_convex_solid_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, pc) = nacre_ops::fixtures::pocketed_cube();
    // Boss on the +x side face, overhanging the bottom edge.
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.25, -0.25]),
        Point3::from_array([1.5, 0.75, 0.75]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Fuse, pc, boss).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, pc, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "overhang boss on non-convex solid: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
    // The boss touches the `x = 1` wall and overlaps nothing: `0.92 + 0.5·0.5·1`.
    let v = mass_props(&m, r).unwrap().volume;
    assert!(
        approx(v, 1.17),
        "overhang boss on non-convex solid: {v} vs hand 1.17"
    );
}

/// A boss raised on the top of a non-convex (L-prism) solid. OCCT scores the boss fused onto
/// the concave base.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn boss_onto_non_convex_solid_matches_occt() {
    use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
    let mut m = Model::new();
    let __w3 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: l, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w3,
            profile: Profile2d::polygon(
                [
                    [0.0, 0.0],
                    [2.0, 0.0],
                    [2.0, 1.0],
                    [1.0, 1.0],
                    [1.0, 2.0],
                    [0.0, 2.0],
                ]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
            )
            .unwrap(),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.3, 0.3, 1.0]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, boss).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, boss).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt.volume),
        "boss onto non-convex solid: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt.volume
    );
}

/// Two face-to-face cubes fuse into a clean box (cell fuse-coplanar-merge merges the
/// coplanar side faces and dissolves the interface corners), which then chains a Cut.
/// OCCT scores both the fuse and the chained cut on the exported b-rep.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn stacked_fuse_then_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
    let stack = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(
        approx(mass_props(&m, stack).unwrap().volume, occt_fuse.volume),
        "stacked fuse: {} vs {}",
        mass_props(&m, stack).unwrap().volume,
        occt_fuse.volume
    );
    // Chain a cut straddling the fused interface.
    let cutter = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    let occt_cut = occt_boolean_of(&m, OcctBool::Cut, stack, cutter).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, stack, cutter).unwrap();
    assert!(
        approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
        "stacked fuse then cut: {} vs {}",
        mass_props(&m, r).unwrap().volume,
        occt_cut.volume
    );
}

/// The two convex pokes, against OCCT. Both operands are convex, and both go down the seam
/// path. The notch is an edge crossed twice; the drill is a genus-1 solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn notch_cube_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([10.0; 3]),
    );
    let y = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([3.0, -1.0, -1.0]),
        Point3::from_array([7.0, 1.4, 1.2]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Cut, a, y).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, a, y).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "vol {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn notch_cube_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([10.0; 3]),
    );
    let y = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([3.0, -1.0, -1.0]),
        Point3::from_array([7.0, 1.4, 1.2]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Fuse, a, y).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, a, y).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "vol {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn drilled_cube_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([3.0; 3]),
    );
    let bar = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Cut, a, bar).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, a, bar).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "vol {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn drilled_cube_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([3.0; 3]),
    );
    let bar = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Fuse, a, bar).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, a, bar).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "vol {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The reflex-corner bite. The arc's single bend projects *outside*
/// its chord's endpoints, so a projection sort would order it only by luck.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn reflex_bite_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_reflex_box();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "reflex bite cut: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// The folded (staircase) arc, `Cut`.
///
/// **This is the only gate on the folded arc.** A folded arc mis-ordered would build a
/// self-intersecting face, and `validate` accepts one — still manifold, Euler holds.
/// Only an independent kernel's volume says the face is right.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn folded_arc_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_popup_box();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "folded arc cut: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// The `Fuse` counterpart of the folded arc.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn folded_arc_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_popup_box();
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "folded arc fuse: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// A stub standing in the L's top face, its footprint strictly inside that face.
/// The seam is a closed loop in the face interior, so the result has a face with
/// an inner loop.
fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
    l_prism_and_box([0.3, 0.3, 0.5], [0.7, 0.7, 1.5])
}

/// nacre's first boolean result carrying a *hole* vs OCCT: the
/// L with a blind pocket. Area is scored alongside volume here — the hole is an
/// area-visible feature, and nacre's own gates all read the same rings, so an
/// independent kernel is what makes the hole's size a real claim.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn blind_dimple_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_dimple();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "blind dimple cut volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "blind dimple cut area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The `Fuse` counterpart — a boss on the L. The same face gains the same hole,
/// so the two agree on area while their volumes straddle the L's own 3.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn blind_dimple_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_dimple();
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "blind dimple fuse volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "blind dimple fuse area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The L-prism and an L-shaped bar lying in its notch, `z ∈ [0.5, 1.5]`, its two arm
/// ends biting the cap's convex corners. Two chords on one face.
fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let raised = SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5]));
    let bar = extrude(
        &mut m,
        raised,
        &[
            [1.8, 0.8],
            [2.1, 0.8],
            [2.1, 2.1],
            [0.8, 2.1],
            [0.8, 1.8],
            [1.8, 1.8],
        ],
    );
    (m, l, bar)
}

/// The L-prism and a П-shaped staple drawn in the XZ plane and extruded along `−y`, so
/// the L's cap is parallel to the extrusion axis and the staple's section there falls in
/// two: a loop wholly inside the cap, and an arc wrapping the reflex corner.
fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
    use nacre_math::Vector3;
    let (mut m, l) = l_prism();
    let xz = SketchPlane::from_origin_normal(
        Point3::from_array([0.0, 1.3, 0.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
    )
    .expect("a unit normal");
    let st = extrude_dist(
        &mut m,
        xz,
        &[
            [0.1, 0.5],
            [0.6, 0.5],
            [0.6, 1.3],
            [0.8, 1.3],
            [0.8, 0.45],
            [1.4, 0.45],
            [1.4, 1.5],
            [0.1, 1.5],
        ],
        0.65,
    );
    (m, l, st)
}

/// The same loop is a **hole** here — inside the cap's kept region — and an **island**
/// under `Cut(staple, L)`, where the kept region is the corner bite alone. Containment
/// alone tells them apart; the three volumes close inclusion–exclusion on `V_∩ = 0.311`.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn l_staple_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, st) = l_and_staple();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, st).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, st).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "l staple cut volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "l staple cut area: {} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn l_staple_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, st) = l_and_staple();
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, st).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, st).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "l staple fuse volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "l staple fuse area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The island case, scored by a kernel that never heard of our containment test.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn staple_cut_by_l_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, st) = l_and_staple();
    let occt = occt_boolean_of(&m, OcctBool::Cut, st, l).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, st, l).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "staple cut by l volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "staple cut by l area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The U-prism (volume 5.3) and a slab shearing off both prong tops. The slab's
/// `y = 1.5` face carries **two** closed loops — one per prong — wholly inside it.
fn u_and_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let u = extrude(
        &mut m,
        SketchPlane::world_xy(),
        &[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 2.3],
            [2.0, 2.3],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ],
    );
    let slab = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5, 1.5, -0.5]),
        Point3::from_array([3.5, 2.5, 1.5]),
    );
    (m, u, slab)
}

/// Three diffs on one fixture, because the three reconstruct different things: two
/// **islands** from one face (`Cut(u, slab)`), the same face's two **holes** on B
/// (`Fuse`), and two holes on A with the U's caps split into cycles (`Cut(slab, u)`).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn u_slab_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, u, slab) = u_and_slab();
    let occt = occt_boolean_of(&m, OcctBool::Cut, u, slab).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, u, slab).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "u slab cut volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "u slab cut area: {} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn u_slab_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, u, slab) = u_and_slab();
    let occt = occt_boolean_of(&m, OcctBool::Fuse, u, slab).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, u, slab).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "u slab fuse volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "u slab fuse area: {} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_cut_by_u_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, u, slab) = u_and_slab();
    let occt = occt_boolean_of(&m, OcctBool::Cut, slab, u).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, slab, u).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "slab cut by u volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "slab cut by u area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The L-prism with an L-shaped stub standing wholly inside its cap. The blind pocket's
/// lid carries the suite's first **non-convex** inner loop.
fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let raised = SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5]));
    let stub = extrude(
        &mut m,
        raised,
        &[
            [0.2, 0.25],
            [0.85, 0.25],
            [0.85, 0.4],
            [0.35, 0.4],
            [0.35, 0.9],
            [0.2, 0.9],
        ],
    );
    (m, l, stub)
}

/// An independent kernel on a hole that is not a rectangle. Volume scores the pocket;
/// area scores its six walls, since the lid gives up exactly what the floor hands back.
/// Neither can see the loop's *winding* — both kernels sum unsigned areas — which is why
/// the winding has a second source and does not lean on this diff.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn ell_dimple_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, stub) = l_and_ell_stub();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, stub).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "ell dimple cut volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "ell dimple cut area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// nacre's first face carrying more than one chord vs OCCT, `Cut`.
/// The L's cap keeps a single ring that uses both arcs, while the bar's floor splits
/// into two faces — the two bites' floors. Volume scores the bites; area cannot, since
/// a corner cut hands back exactly the faces it removes (14.0 either way).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn notch_bar_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bar) = l_and_notch_bar();
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, bar).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bar).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "notch bar cut volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "notch bar cut area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The `Fuse` counterpart, and not a symmetry re-ask: it reconstructs different faces.
/// Keeping both outsides, the bar's floor stays one ring spanning both arcs instead of
/// splitting, and the area — `14 + 6.58 − 0.96` — finally has something to score.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn notch_bar_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bar) = l_and_notch_bar();
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bar).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bar).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "notch bar fuse volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "notch bar fuse area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// nacre's first face whose outer loop is *all* seam vs OCCT: the
/// stub cut by the L, leaving the `0.4 × 0.4 × 0.5` box above `z = 1`. Its floor is
/// that island face.
///
/// Both kernels take a face's normal from their own bookkeeping and sum unsigned
/// areas, so neither this diff nor any other can see the island wound backwards —
/// `validate` and the signed mesh volume do that. What an independent kernel scores
/// here is the face's *existence and extent*: get the ring's nodes wrong and the
/// polygon's area moves with it.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn island_cut_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, bx) = l_and_dimple();
    let occt = occt_boolean_of(&m, OcctBool::Cut, bx, l).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, bx, l).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "island cut volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "island cut area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// nacre's hole-aware classification vs OCCT: a pocketed cube
/// cut by a box sitting wholly inside the pocket void. The two solids are
/// disjoint, so the answer is the pocketed cube untouched — but only if the
/// classifier reads the lid's inner loop. Fanning the lid's outer ring alone
/// put the box's corners on both sides of the boundary.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocketed_cut_in_the_void_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, solid) = nacre_ops::fixtures::pocketed_cube();
    let bx = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.4, 0.4, 0.6]),
        Point3::from_array([0.6, 0.6, 0.9]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Cut, solid, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, solid, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "pocketed cut in the void: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.volume, 0.92),
        "pocketed cut in the void: {} vs the untouched 0.92",
        nacre.volume
    );
}

/// nacre's coincident-coplanar merge (M5-c5) vs OCCT: two cubes stacked on a
/// shared z=1 face fuse to a 1×1×2 box.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn stacked_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
}
