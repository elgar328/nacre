//! Non-convex operands — L-prisms with boxes biting convex corners, reflex
//! corners, containment, and closed seam loops. These are the cases a convex
//! half-space classifier gets wrong; only exact point-in-solid gets the volume
//! right.

use crate::common::{
    boolean_one, l_and_corner_box, l_and_inner_box, l_and_notch_bar, l_and_reflex_box, l_prism,
    u_and_slab,
};
use nacre_math::Point3;
use nacre_ops::{BoolKind, boolean};
use nacre_topo::Model;

#[test]
fn cut_non_convex_containment_makes_cavity() {
    // Cut(L − box) with box ⊂ L ⇒ a hollow L (outer L shell + box void).
    let (mut m, l, bx) = l_and_inner_box();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solids.get(r).cavities.len(), 1);
}

#[test]
fn fuse_non_convex_containment_is_container() {
    let (mut m, l, bx) = l_and_inner_box();
    let vol_l = nacre_props::mass_props(&m, l).unwrap().volume;
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_l).abs() < 1e-9);
    assert!(m.solids.get(r).cavities.is_empty());
}

#[test]
fn common_non_convex_containment_is_inner() {
    let (mut m, l, bx) = l_and_inner_box();
    let vol_bx = nacre_props::mass_props(&m, bx).unwrap().volume;
    let r = boolean_one(&mut m, BoolKind::Common, l, bx).unwrap();
    assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_bx).abs() < 1e-9);
}

#[test]
fn cut_box_inside_non_convex_is_empty() {
    // Cut(box − L): the box is wholly inside L ⇒ nothing remains. Nothing is an answer, so the
    // boolean succeeds with no solids — and consumes both operands like any other success.
    let (mut m, l, bx) = l_and_inner_box();
    assert!(boolean(&mut m, BoolKind::Cut, bx, l).unwrap().is_empty());
    assert!(m.live_solids.is_empty(), "both operands are consumed");
}

#[test]
fn disjoint_non_convex_operand() {
    // L and a far box (non-coplanar): Cut ⇒ L, Fuse ⇒ empty, Common ⇒ empty.
    let far = || Point3::from_array([10.0; 3]);
    let far_max = || Point3::from_array([11.0; 3]);
    let (mut m, l) = l_prism();
    let vol_l = nacre_props::mass_props(&m, l).unwrap().volume;
    let d = m.add_cuboid(far(), far_max());
    let r = boolean_one(&mut m, BoolKind::Cut, l, d).unwrap();
    assert!((nacre_props::mass_props(&m, r).unwrap().volume - vol_l).abs() < 1e-9);
    // Fusing things that never touch does not merge them — it keeps both, whole.
    let (mut m, l) = l_prism();
    let vol_d = 1.0; // far()..far_max() is the unit box
    let d = m.add_cuboid(far(), far_max());
    let both = boolean(&mut m, BoolKind::Fuse, l, d).unwrap();
    assert_eq!(both.len(), 2, "disjoint operands stay two solids");
    m.rebuild_adjacency();
    let mut vols: Vec<f64> = both
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .collect();
    vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
    assert!(
        (vols[0] - vol_d).abs() < 1e-9 && (vols[1] - vol_l).abs() < 1e-9,
        "{vols:?}"
    );
    // Their intersection, on the other hand, really is empty.
    let (mut m, l) = l_prism();
    let d = m.add_cuboid(far(), far_max());
    assert!(boolean(&mut m, BoolKind::Common, l, d).unwrap().is_empty());
}

#[test]
fn cut_non_convex_overlap_corner_bite() {
    // Cut(L − box): the corner bite carves 0.224 off the L.
    let (mut m, l, bx) = l_and_corner_box();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.224)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn fuse_non_convex_overlap_corner_bite() {
    // Fuse(L ∪ box): the protruding box adds (0.924 − 0.224) to the L.
    let (mut m, l, bx) = l_and_corner_box();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 + 0.924 - 0.224)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn common_non_convex_overlap_is_their_intersection() {
    let (mut m, l, bx) = l_and_corner_box();
    let r = boolean_one(&mut m, BoolKind::Common, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.224).abs() < 1e-9, "volume {vol}");
}

#[test]
fn a_common_can_leave_a_closed_seam_loop() {
    // A bar drilled through a cube: ∩ = [1,2]²×[0,3], the middle segment. On the
    // z=0/z=3 caps the kept square [1,2]² is bounded entirely by the cut (an island
    // face) — the loop is oriented with material inside it, a sign a hole never hits.
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let bar = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Common, cube, bar).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 3.0).abs() < 1e-9, "volume {vol}");
}

#[test]
fn cut_across_reflex_corner_bite() {
    // Only exact point_in_solid gets the volume right; a convex half-space test
    // would misclassify the box vertex sitting in the L's notch.
    let (mut m, l, bx) = l_and_reflex_box();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    // Overlap = xy(1.0 − notch 0.36 = 0.64) · z(0.8) = 0.512.
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn cut_u_by_slab() {
    // Two island loops on the slab's y=1.5 face: the two prong cross-sections become
    // the result's new end caps (one input face → two output faces, both flipped).
    let (mut m, u, slab) = u_and_slab();
    let r = boolean_one(&mut m, BoolKind::Cut, u, slab).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!((props.volume - 4.0).abs() < 1e-9, "volume {}", props.volume);
    assert!((props.area - 18.0).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn fuse_u_and_slab() {
    // The same face the other way: Fuse keeps both outsides, so both loops are holes.
    let (mut m, u, slab) = u_and_slab();
    let r = boolean_one(&mut m, BoolKind::Fuse, u, slab).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!((props.volume - 12.0).abs() < 1e-9, "volume {}", props.volume);
    assert!((props.area - 42.0).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn cut_slab_by_u() {
    // Operands swapped: the holed face is now on A, and the U's caps split into two
    // cycles apiece under flip — two blind pockets.
    let (mut m, u, slab) = u_and_slab();
    let r = boolean_one(&mut m, BoolKind::Cut, slab, u).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!((props.volume - 6.7).abs() < 1e-9, "volume {}", props.volume);
    assert!((props.area - 33.2).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn cut_notch_bar() {
    // Two chords on one face: the bar bites the cap's corners (2,1) and (1,2); the kept
    // region is the cap minus both, one ring using both arcs. 3 − 2·(0.2·0.2·0.5) = 2.96.
    let (mut m, l, bar) = l_and_notch_bar();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bar).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!((props.volume - 2.96).abs() < 1e-9, "volume {}", props.volume);
    assert!((props.area - 14.0).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn fuse_notch_bar() {
    // The Fuse counterpart, closing inclusion–exclusion: 3 + 0.69 − 0.04.
    let (mut m, l, bar) = l_and_notch_bar();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bar).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!((props.volume - (3.0 + 0.69 - 0.04)).abs() < 1e-9, "volume {}", props.volume);
    assert!((props.area - 19.62).abs() < 1e-9, "area {}", props.area);
}
