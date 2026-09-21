//! Booleans against OCCT: coplanar contact, shared planes, tunnels.

use super::*;

/// Two overlapping unit boxes A = [0,1]³, B = [0.5,1.5]³ (overlap [0.5,1]³ =
/// 0.125). Hand-computable boolean volumes: union 1.875, difference A−B
/// 0.875, intersection 0.125 — so they calibrate the OCCT boolean harness
/// end to end (helper commands + single-solid export + parse).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn boolean_volumes_match_occt() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
    let cut = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
    let common = occt_boolean_of(&m, OcctBool::Common, a, b).unwrap();
    assert!(approx(fuse.volume, 1.875), "fuse {}", fuse.volume);
    assert!(approx(cut.volume, 0.875), "cut {}", cut.volume);
    assert!(approx(common.volume, 0.125), "common {}", common.volume);
}

/// A top-flush tool punching clean through the base — the pocket becomes a bore. Cut leaves
/// 1 − 0.5·0.5·1 = 0.75, and Fuse keeps the stub that pokes out below for 1.125.
///
/// **Area is the point.** The claim being checked is topological — that the exit face came out
/// annular so the bore's walls have something to close against — and volume cannot see that.
/// Cut: 0.75 + 0.75 + 4 + 2.0 (bore walls) = 7.5, against 6.0 for the untouched cube. Fuse:
/// 1.0 + 0.75 + 4 + 1.0 + 0.25 = 7.0. Both used to be rejected outright.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn punch_through_bore_matches_occt() {
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let through = m.add_cuboid(
        Point3::from_array([0.25, 0.25, -0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let cut = occt_boolean_of(&m, OcctBool::Cut, base, through).unwrap();
    assert!(approx(cut.volume, 0.75), "bore cut volume {}", cut.volume);
    assert!(approx(cut.area, 7.5), "bore cut area {}", cut.area);
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, base, through).unwrap();
    assert!(
        approx(fuse.volume, 1.125),
        "stub fuse volume {}",
        fuse.volume
    );
    assert!(approx(fuse.area, 7.0), "stub fuse area {}", fuse.area);
}

/// The Cut twin of the corner-overhang boss: base [0,1]³ and a boss [0.5,1.5]²×[1,2] seated on
/// z=1 with an overhanging footprint. The boss lies entirely above the shared plane, so the cut
/// removes nothing and the base survives whole — volume 1.0. nacre used to reject this
/// (`coplanar_merge`) while building the Fuse of the very same pair, so the interesting claim is
/// "nothing was removed", which is exactly the kind of answer worth hearing from a second kernel.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn cut_by_a_corner_overhanging_boss_matches_occt() {
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let corner = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let cut = occt_boolean_of(&m, OcctBool::Cut, base, corner).unwrap();
    assert!(approx(cut.volume, 1.0), "seated-boss cut {}", cut.volume);
}

/// A coplanar contact on a **cavitied** operand: a hollow box ([0,3]³ minus a [1,2]³ void,
/// volume 26) with a top-flush boss on z=3 ([0.5,0.75]²×[3,4], volume 0.0625). The union keeps
/// the void ⇒ 26.0625. This is the independent check on the all-shell fix: the coplanar driver
/// used to emit outer-shell faces only, so the void vanished and the fuse read 27.0625 — the
/// un-hollowed cube plus the boss — with a clean `validate`, since what remained was still a
/// closed shell. Hand arithmetic and nacre agreeing would have proved nothing there; OCCT is a
/// second kernel. (`hollow_solid_step_volume_matches_occt` already pins that a void survives the
/// STEP round trip, so a failure here is about the boolean, not the transport.)
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn coplanar_boss_on_a_hollow_part_matches_occt() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 3.0]),
        Point3::from_array([0.75, 0.75, 4.0]),
    );
    m.rebuild_adjacency();
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, hollow, boss).unwrap();
    assert!(
        approx(fuse.volume, 26.0625),
        "hollow + coplanar boss fuse {}",
        fuse.volume
    );
}

/// A sever that leaves a surviving cavity, cross-checked against OCCT. A hollow box ([0,3]³
/// minus a [0.5,1.5]×[0.5,2.5]×[0.5,2.5] void, material 23) cut by a slab at x∈[2,2.2] severs
/// into two pieces (total material 21.2), the x<2 piece keeping the void. nacre assigns the
/// void by containment (`point_in_component`); OCCT is the second kernel confirming the total.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn severed_hollow_box_matches_occt() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 2.5, 2.5]),
    );
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let slab = m.add_cuboid(
        Point3::from_array([2.0, -1.0, -1.0]),
        Point3::from_array([2.2, 4.0, 4.0]),
    );
    m.rebuild_adjacency();
    let cut = occt_boolean_of(&m, OcctBool::Cut, hollow, slab).unwrap();
    assert!(
        approx(cut.volume, 21.2),
        "severed hollow box {}",
        cut.volume
    );
}

/// same_ground union: A = [0,1]³ and B = [0.5,1.5]²×[0,1] overlap in volume and share the
/// z=0 / z=1 planes with overlapping footprints. The union is an L-footprint prism: area
/// (1 + 1 − 0.25) × height 1 = 1.75. OCCT confirms it independently — the oracle for nacre's
/// general 2D coplanar merge (the driver's union-cell reconstruct + interpenetrating-wall clip).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn same_ground_union_matches_occt() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.0]),
        Point3::from_array([1.5, 1.5, 1.0]),
    );
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
    assert!(
        approx(fuse.volume, 1.75),
        "same_ground fuse {}",
        fuse.volume
    );
}

/// Cut/Common of the same_ground config (A=[0,1]³, B=[0.5,1.5]²×[0,1]). A−B removes the
/// overlap column [0.5,1]²×[0,1] ⇒ L-prism 0.75; A∩B is that overlap box ⇒ 0.25. OCCT confirms
/// both independently — the oracle for nacre's same_ground MinusQ/InterQ caps + wall clip.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn same_ground_cut_common_matches_occt() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.0]),
        Point3::from_array([1.5, 1.5, 1.0]),
    );
    let cut = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
    let common = occt_boolean_of(&m, OcctBool::Common, a, b).unwrap();
    assert!(approx(cut.volume, 0.75), "same_ground cut {}", cut.volume);
    assert!(
        approx(common.volume, 0.25),
        "same_ground common {}",
        common.volume
    );
}

/// Single-shared-plane Fuse: A=[0,1]³ and B=[0.5,1.5]²×[0,2] share only z=0 with overlapping
/// footprints; B is taller and pokes through A's top. The union is a stepped prism: z0–1 footprint
/// 1.75 + z1–2 tower 1.0 = 2.75. OCCT confirms it independently — the oracle for nacre's mixed
/// coplanar-cap + transversal-exit route.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn single_shared_plane_fuse_matches_occt() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
    assert!(
        approx(fuse.volume, 2.75),
        "single-shared fuse {}",
        fuse.volume
    );
}

/// Single-shared-plane Cut/Common (same config as the Fuse oracle above). b's tower above z=1 is
/// irrelevant to A−B and A∩B, so A−B = L-prism 0.75, A∩B = the overlap box 0.25 (same as
/// same_ground). OCCT confirms both independently — the oracle for nacre's tower-drop.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn single_shared_plane_cut_common_matches_occt() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let cut = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
    let common = occt_boolean_of(&m, OcctBool::Common, a, b).unwrap();
    assert!(approx(cut.volume, 0.75), "single-shared cut {}", cut.volume);
    assert!(
        approx(common.volume, 0.25),
        "single-shared common {}",
        common.volume
    );
}

/// Spanning-slab single-shared variant: base=[0,1]³, slab=[-0.5,1.5]×[0.4,0.6]×[0.5,1] shares the
/// base top (z=1, same-normal) and overhangs two opposite x edges (four crossings). OCCT confirms
/// all three: Fuse 1.1, Cut 0.9, Common 0.1 — the oracle for the broader single-shared route.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_overhang_matches_occt() {
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let slab = m.add_cuboid(
        Point3::from_array([-0.5, 0.4, 0.5]),
        Point3::from_array([1.5, 0.6, 1.0]),
    );
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, base, slab).unwrap();
    let cut = occt_boolean_of(&m, OcctBool::Cut, base, slab).unwrap();
    let common = occt_boolean_of(&m, OcctBool::Common, base, slab).unwrap();
    assert!(approx(fuse.volume, 1.1), "slab fuse {}", fuse.volume);
    assert!(approx(cut.volume, 0.9), "slab cut {}", cut.volume);
    assert!(approx(common.volume, 0.1), "slab common {}", common.volume);
}

/// A flush-edge pocket: the cutter sits flush on two adjacent base faces (top z=1 and
/// front y=0), so its walls are coplanar with the part's walls along the shared boundary edge.
/// OCCT confirms base − cutter = 0.92 independently.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn flush_edge_pocket_cut_matches_occt() {
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let cutter = m.add_cuboid(
        Point3::from_array([0.3, 0.0, 0.5]),
        Point3::from_array([0.7, 0.4, 1.0]),
    );
    let cut = occt_boolean_of(&m, OcctBool::Cut, base, cutter).unwrap();
    assert!(approx(cut.volume, 0.92), "flush cut {}", cut.volume);
}

/// A through-tunnel Cut: a cutter spanning the bar's full height (top and bottom both flush)
/// makes two parallel coplanar contacts, cut into a tunnel.
/// OCCT confirms bar − tunnel = 0.84 independently.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn through_tunnel_cut_matches_occt() {
    let mut m = Model::new();
    let bar = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let cutter = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.0]),
        Point3::from_array([0.7, 0.7, 1.0]),
    );
    let cut = occt_boolean_of(&m, OcctBool::Cut, bar, cutter).unwrap();
    assert!(approx(cut.volume, 0.84), "tunnel cut {}", cut.volume);
}

/// Disjoint boxes fuse to a compound whose total volume is the sum — a sanity
/// check that the harness handles a non-overlapping (compound) result.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn disjoint_fuse_sums_volumes() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([2.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );
    let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
    assert!(approx(fuse.volume, 2.0), "disjoint fuse {}", fuse.volume);
}

/// A severing `Cut` returns two nacre solids; OCCT returns a COMPOUND of two solids for the
/// same inputs. Their aggregate volume and area agree — nacre's multi-solid output
/// matches OCCT. The bar threads the cube and out both ends, leaving two 1×1×1 stubs.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn sever_cut_matches_occt_compound() {
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let bar = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    // OCCT ground truth from the inputs, before nacre supersedes them.
    let occt = occt_boolean_of(&m, OcctBool::Cut, bar, cube).unwrap();
    let solids = boolean(&mut m, BoolKind::Cut, bar, cube).unwrap();
    assert_eq!(solids.len(), 2, "nacre severs into two solids");
    let vol: f64 = solids
        .iter()
        .map(|&s| mass_props(&m, s).unwrap().volume)
        .sum();
    let area: f64 = solids
        .iter()
        .map(|&s| mass_props(&m, s).unwrap().area)
        .sum();
    assert!(
        approx(vol, occt.volume),
        "volume {vol} vs occt {}",
        occt.volume
    );
    assert!(approx(area, occt.area), "area {area} vs occt {}", occt.area);
}

/// A `Transform`-translated solid composes correctly into a boolean: translate a cube by a
/// rational offset, then `Cut` an overlapping cube from it. The moved geometry's boolean matches OCCT on the same inputs —
/// so the geometry-rewrite produced a boolean-valid solid, not just a
/// volume-invariant one.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn translated_solid_cut_matches_occt() {
    use nacre_exact::{Isometry, Rat};
    use nacre_ops::{BoolKind, Operation, apply, boolean};
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    // Translate A by [3/10, 2/5, 1/2] ⇒ A' = [0.3,1.3]×[0.4,1.4]×[0.5,1.5].
    let iso = Isometry::translation([
        Rat::new(3, 10).unwrap(),
        Rat::new(2, 5).unwrap(),
        Rat::new(1, 2).unwrap(),
    ]);
    let out = apply(
        &mut m,
        &Operation::Transform {
            solid: a,
            isometry: iso,
        },
    )
    .unwrap();
    let nacre_ops::OpOutput::Transform { solid: a2 } = out else {
        panic!("expected Transform output");
    };
    // B overlaps A' at a corner.
    let b = m.add_cuboid(Point3::from_array([0.8; 3]), Point3::from_array([1.8; 3]));
    // OCCT ground truth on the moved inputs, before nacre supersedes them.
    let occt = occt_boolean_of(&m, OcctBool::Cut, a2, b).unwrap();
    let solids = boolean(&mut m, BoolKind::Cut, a2, b).unwrap();
    let vol: f64 = solids
        .iter()
        .map(|&s| mass_props(&m, s).unwrap().volume)
        .sum();
    let area: f64 = solids
        .iter()
        .map(|&s| mass_props(&m, s).unwrap().area)
        .sum();
    assert!(
        approx(vol, occt.volume),
        "volume {vol} vs occt {}",
        occt.volume
    );
    assert!(approx(area, occt.area), "area {area} vs occt {}", occt.area);
}

/// A `Transform`-rotated solid is a well-formed b-rep:
/// rotate a cuboid 30° about Z through a rational axis point, then export just
/// that solid and ask OCCT for its volume/area. A rigid rotation leaves both
/// invariant, so OCCT must agree with nacre's `mass_props` — proving the
/// rotation rewrite produced a valid solid, not merely a volume-preserving
/// vertex shuffle. The same check on a rotated `Cut` result exercises the
/// `Discovered`-vertex rotation path. Single-solid export (`to_step_solid`)
/// avoids summing the superseded input still resident in the arena.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_solid_props_match_occt() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    use nacre_ops::{Operation, apply};
    let rot30 = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });

    // (a) rotate a plain cuboid.
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let out = apply(
        &mut m,
        &Operation::Transform {
            solid: c,
            isometry: rot30,
        },
    )
    .unwrap();
    let nacre_ops::OpOutput::Transform { solid: c2 } = out else {
        panic!("expected Transform output");
    };
    let occt = occt_props(&nacre_step::to_step_solid(&m, c2).unwrap()).unwrap();
    let nacre = mass_props(&m, c2).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "volume {} vs occt {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs occt {}",
        nacre.area,
        occt.area
    );

    // (b) rotate a Cut result (Discovered seam vertices → Rotated).
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    let cut = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    let out = apply(
        &mut m,
        &Operation::Transform {
            solid: cut,
            isometry: rot30,
        },
    )
    .unwrap();
    let nacre_ops::OpOutput::Transform { solid: cut2 } = out else {
        panic!("expected Transform output");
    };
    let occt = occt_props(&nacre_step::to_step_solid(&m, cut2).unwrap()).unwrap();
    let nacre = mass_props(&m, cut2).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "cut volume {} vs occt {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "cut area {} vs occt {}",
        nacre.area,
        occt.area
    );
}
