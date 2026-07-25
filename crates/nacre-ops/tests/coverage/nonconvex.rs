//! Non-convex operands — L/U prisms, reflex corners, notches, dimples, staples.

#![allow(unused_imports)]
use crate::common::*;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean,
    replay,
};
use nacre_scalar::Axis;
use nacre_store::Handle;
use nacre_topo::{Face, Loop, Model, Orientation, Shell, Solid, Vertex};

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
    assert!(
        (props.volume - 12.0).abs() < 1e-9,
        "volume {}",
        props.volume
    );
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
    assert!(
        (props.volume - 2.96).abs() < 1e-9,
        "volume {}",
        props.volume
    );
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
    assert!(
        (props.volume - (3.0 + 0.69 - 0.04)).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 19.62).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn cut_staircase_seam_arc() {
    // On the box's bottom face the seam runs (2,0.5) → (2,1) → (1,1) → (1,1.5): a
    // staircase whose two bends turn opposite ways. `strict` used to reject it,
    // unable to tell a reflex turn from an arc folded back on itself.
    //
    // The reconstructed face there is (0.5,0.5) → (0.5,1.5) → (1,1.5) → (1,1) →
    // (2,1) → (2,0.5): the overlap footprint, area 1.0, reflex at (1,1). A correct
    // simple polygon. `strict` was rejecting a right answer.
    //
    // Overlap = footprint 1.0 × z∈[0.2,1] = 0.8.
    let (mut m, l, bx) = l_and_popup_box();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.8)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn fuse_staircase_seam_arc() {
    // The `Fuse` counterpart, closing the inclusion–exclusion: V_L + V_box − 0.8.
    // `validate` cannot see a self-intersecting face (it stays manifold, Euler
    // holds), so the volume is what pins the folded arc — with OCCT alongside.
    let (mut m, l, bx) = l_and_popup_box();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 + 2.0 - 0.8)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn cut_ell_dimple() {
    // The blind pocket is the stub's L-shaped section, `0.65·0.15 + 0.15·0.5 = 0.1725`,
    // half a unit deep. Area `14 − 0.1725 + 2.6·0.5 + 0.1725`: the lid gives up exactly
    // what the floor hands back, so only the walls move it.
    let (mut m, l, stub) = l_and_ell_stub();
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!(
        (props.volume - (3.0 - 0.1725 * 0.5)).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 15.3).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn cut_l_staple() {
    // A loop beside an arc, resolved. The cap's kept region is the hexagon minus the
    // corner bite, and the near leg's rectangle sits inside it — so the loop is a
    // **hole** of that region. The winding says clockwise, which is only a cross-check
    // now: containment decided.
    //
    // `V_∩ = 0.5·0.5·0.65 + 0.27·0.55 = 0.311`, the two legs' parts inside the L.
    let (mut m, l, st) = l_and_staple();
    let r = boolean_one(&mut m, BoolKind::Cut, l, st).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!(
        (props.volume - 2.689).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 15.755).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn fuse_l_staple() {
    // `3 + 0.7605 − 0.311`. The staple measures `1.17` in section, `0.65` deep.
    let (mut m, l, st) = l_and_staple();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, st).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!(
        (props.volume - (3.0 + 0.7605 - 0.311)).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 16.72).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn cut_staple_by_l() {
    // ★ The same loop, on the same face, is now an **island**. Swap the operands and the
    // cap's kept region becomes the corner bite alone; the loop lies outside it, in a
    // dropped region, so its interior is what survives.
    //
    // Cell 3f-3's counterexample made flesh: a hole and an island wind oppositely without
    // nesting, so no winding could have told these two apart. Position did.
    //
    // `0.7605 − 0.311`, and the three volumes close inclusion–exclusion exactly.
    let (mut m, l, st) = l_and_staple();
    let r = boolean_one(&mut m, BoolKind::Cut, st, l).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!(
        (props.volume - (0.7605 - 0.311)).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 4.68).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn cut_the_stub_by_the_l_leaves_an_island_face() {
    // Swap the operands of `cut_blind_dimple` and the same loop lands on a face whose
    // boundary is *all* dropped: the kept region is the loop's interior alone. That
    // face has no `∂f` at all — its outer loop *is* the seam ring, four `Discovered`
    // vertices and nothing else. The answer is the `0.4 × 0.4 × 0.5` box above `z = 1`.
    //
    // This is where `flip` first meets a discovered hole ring. It is wound CCW about the
    // L's `+z`, keeping the material (inside the stub) on
    // its left; `flip` reverses it and the face becomes the box's downward-facing
    // floor. Nothing but `validate` and the signed mesh volume can see that go wrong.
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Cut, stub, l).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.08).abs() < 1e-9, "volume {vol}");
}

#[test]
fn fuse_the_stub_and_the_l() {
    // Not a new branch — the same hole, on the same face of the L, reached with the
    // operands the other way round. `Fuse` keeps both outsides and flips neither, so
    // what this pins is that the answer does not depend on which solid is `a`: the
    // hole now lands on B, and 3.08 is 3.08.
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Fuse, stub, l).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 + 0.16 - 0.08)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn cut_blind_dimple() {
    // The stub's footprint never reaches the L's top-face boundary, so the seam is
    // a closed ring in the face interior. It is the face's inner loop, and the
    // result is a blind pocket: `3 − 0.4² × 0.5`.
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.08)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn a_flipped_hole_loop_is_caught() {
    // `loop_orient_mismatch` cannot see a loop flipped as a whole, and `f2`'s golden
    // pins the derivation — but the two downstream detectors must actually fire on
    // *this* shape, not merely exist. They are different in kind: `validate` sees a
    // rim edge used twice the same way, `tessellate` sees a hole wound like its
    // outer ring. Volume and OCCT see neither: `props` sums `|area|`.
    //
    // `Store` is append-only, so the face cannot be edited. Push a replacement with
    // the hole reversed, swap it into a fresh shell and solid, and move the live
    // handle: the old face falls out of `reachable()`.
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();

    // The dimple is blind, so exactly one face carries a hole. Counted, not assumed:
    // its sibling in `pipeline.rs` said "the only face …" in prose, used `find`, and
    // silently flipped a different face once the engine stopped making the predicate
    // unique.
    let faces = m.shells.get(m.solids.get(r).outer).faces.clone();
    let holed_faces: Vec<_> = faces
        .iter()
        .copied()
        .filter(|&f| !m.faces.get(f).inner.is_empty())
        .collect();
    assert_eq!(holed_faces.len(), 1, "the L's top face carries the hole");
    let holed = holed_faces[0];
    let f = m.faces.get(holed).clone();
    let mut hole = f.inner[0].clone();
    hole.half_edges.reverse();
    for he in &mut hole.half_edges {
        he.forward = !he.forward;
    }
    let bad = m.faces.push(Face {
        inner: vec![hole],
        ..f
    });
    let swapped = faces
        .iter()
        .map(|&x| if x == holed { bad } else { x })
        .collect();
    let shell = m.shells.push(Shell { faces: swapped });
    let solid = m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    m.live_solids = vec![solid];
    m.rebuild_adjacency();

    let vs = nacre_validate::validate(&m);
    assert!(
        vs.iter()
            .any(|v| matches!(v, nacre_validate::Violation::NonOpposedEdge { .. })),
        "validate stayed quiet: {vs:?}"
    );
    assert!(matches!(
        nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()),
        Err(nacre_tess::TessError::HoleWinding)
    ));
}

#[test]
fn fuse_blind_dimple() {
    // The `Fuse` counterpart — a boss on the L — closing the inclusion–exclusion:
    // `V_L + V_stub − V_overlap`. Both put a hole in the same face of the L.
    let (mut m, l, stub) = l_and_dimple();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, stub).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 + 0.16 - 0.08)).abs() < 1e-9, "volume {vol}");
}

/// The first time a boolean result is fed back as an operand: the overlap box
/// (Discovered corners) stacked on a third box merges through the coincident-
/// interface path — which runs `is_convex` on that Discovered-cornered
/// operand. Volume is the sum and the shell stays closed.
#[test]
fn a_boolean_result_stacks_as_an_operand() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, 1.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let c = boolean_one(&mut m, BoolKind::Common, a, b).unwrap(); // [1,2]³
    m.rebuild_adjacency();
    let d = m.add_cuboid(
        Point3::from_array([1.0, 1.0, 2.0]),
        Point3::from_array([2.0, 2.0, 3.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, c, d).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 2.0).abs() < 1e-12, "volume {vol}");
    assert_eq!(m.solids.get(r).cavities.len(), 0);
}

#[test]
fn blind_hole_drills_into_a_void() {
    // A cut reaching *into* a void: a stub drilled from below the hollow box up into
    // its void. The void loses its enclosure and merges with the outer shell (an
    // open pocket, cavities 0). The seam machinery (all-shells) + the (5c) component
    // split reconstruct it exactly: remove the channel [1.4,1.6]²×[0,1]=0.04 through
    // the floor, void interior removes nothing ⇒ 26 − 0.04 = 25.96.
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    assert_eq!(m.solids.get(hollow).cavities.len(), 1);
    m.rebuild_adjacency();
    let stub = m.add_cuboid(
        Point3::from_array([1.4, 1.4, -0.5]),
        Point3::from_array([1.6, 1.6, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, hollow, stub).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 25.96).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solids.get(r).cavities.len(), 0); // the void opened to outside
}

#[test]
fn a_tunnel_drilled_through_a_void() {
    // A cut passing all the way through the void (bottom to top). Removes the floor
    // and ceiling channels [1.4,1.6]²×([0,1]∪[2,3]) = 0.08; the void interior removes
    // nothing ⇒ 26 − 0.08 = 25.92. Result is a genus-1 solid (a straight tunnel),
    // one shell, no cavity.
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let tunnel = m.add_cuboid(
        Point3::from_array([1.4, 1.4, -0.5]),
        Point3::from_array([1.6, 1.6, 3.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, hollow, tunnel).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 25.92).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solids.get(r).cavities.len(), 0);
}

#[test]
fn clockwise_input_is_auto_corrected() {
    // The square wound CW; auto-CCW makes it a valid cube anyway.
    let cw = Profile2d {
        points: vec![p2(0.0, 1.0), p2(1.0, 1.0), p2(1.0, 0.0), p2(0.0, 0.0)],
    };
    let m = replay(&[extrude_op(cw, 1.0)]).unwrap();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(m.faces.len(), 6);
}

#[test]
fn degenerate_inputs_are_rejected() {
    let plane = SketchPlane::world_xy();
    let two = Profile2d {
        points: vec![p2(0.0, 0.0), p2(1.0, 0.0)],
    };
    assert_eq!(
        apply(
            &mut Model::new(),
            &Operation::Extrude {
                plane,
                profile: two,
                dist: 1.0
            }
        ),
        Err(OpError::DegenerateProfile)
    );
    assert_eq!(
        apply(&mut Model::new(), &extrude_op(square(), 0.0)),
        Err(OpError::NonPositiveDistance)
    );
    let dup = Profile2d {
        points: vec![p2(0.0, 0.0), p2(0.0, 0.0), p2(1.0, 1.0)],
    };
    assert_eq!(
        apply(
            &mut Model::new(),
            &Operation::Extrude {
                plane,
                profile: dup,
                dist: 1.0
            }
        ),
        Err(OpError::DegenerateGeometry)
    );
}

/// The `Transform` op flows through `apply`, superseding via the op dispatch.
#[test]
fn transform_op_applies() {
    let (iso, _) = test_iso();
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let out = apply(
        &mut m,
        &Operation::Transform {
            solid: c,
            isometry: iso,
        },
    )
    .unwrap();
    match out {
        OpOutput::Transform { solid } => assert_eq!(m.live_solids, vec![solid]),
        other => panic!("expected Transform output, got {other:?}"),
    }
}

#[test]
fn cut_containment_makes_a_cavity() {
    // A = [0,3]³ (27) with B = [1,2]³ (1) strictly inside ⇒ A − B is a
    // hollow solid: volume 26, an outer + one void shell (V16/E24/F12/S2).
    let (mut m, a, b) = nested_boxes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 26.0).abs() < 1e-12, "volume {vol}");
    assert_eq!(m.solids.get(r).cavities.len(), 1);
    let reach = m.reachable();
    assert_eq!(reach.shells.len(), 2);
    assert_eq!(reach.faces.len(), 12);
    assert_eq!(reach.vertices.len(), 16);
    assert_eq!(reach.edges.len(), 24);
    assert_eq!(m.live_solids, vec![r]);
}

#[test]
fn cut_containment_off_center_cavity() {
    // The inner box need not be concentric — any strictly-interior B works.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 2.5, 3.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (64.0 - 6.0)).abs() < 1e-12, "volume {vol}"); // 4³ − 1·2·3
    assert_eq!(m.solids.get(r).cavities.len(), 1);
}

#[test]
fn fuse_containment_is_the_container() {
    // A ∪ B with B ⊂ A is just A (no cavity).
    let (mut m, a, b) = nested_boxes();
    let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
    let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - vol_a).abs() < 1e-12, "volume {vol}");
    assert!(m.solids.get(r).cavities.is_empty());
    assert_eq!(m.live_solids, vec![r]);
}

#[test]
fn containment_symmetric_when_a_inside_b() {
    // Arguments swapped: A = inner ⊂ B = outer.
    // Cut(inner − outer): inner is wholly removed ⇒ empty, which is an answer, not an error —
    // and a successful boolean consumes its operands, so the Fuse below needs a fresh model
    // (it used to reuse this one only because the empty Cut was an error that consumed nothing).
    let (mut m, outer, inner) = nested_boxes();
    assert!(
        boolean(&mut m, BoolKind::Cut, inner, outer)
            .unwrap()
            .is_empty()
    );
    assert!(m.live_solids.is_empty(), "both operands are consumed");

    // Fuse(inner ∪ outer) = outer.
    let (mut m, outer, inner) = nested_boxes();
    let vol_outer = nacre_props::mass_props(&m, outer).unwrap().volume;
    let r = boolean_one(&mut m, BoolKind::Fuse, inner, outer).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - vol_outer).abs() < 1e-12, "volume {vol}");
}

#[test]
fn common_containment_is_the_inner_solid() {
    // Either argument order ⇒ the intersection is the inner solid (handled
    // by the existing half-space enumeration path, no cavity code).
    for swap in [false, true] {
        let (mut m, outer, inner) = nested_boxes();
        let vol_inner = nacre_props::mass_props(&m, inner).unwrap().volume;
        let (x, y) = if swap { (inner, outer) } else { (outer, inner) };
        let r = boolean_one(&mut m, BoolKind::Common, x, y).unwrap();
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - vol_inner).abs() < 1e-9, "swap={swap} volume {vol}");
    }
}

#[test]
fn a_corner_cut_through_the_bottom() {
    // The corner prism pokes out the base's bottom, so the old blind gate declined it. The F2
    // collapse routes it to the unified driver: base 1.0 − corner column (x,y ∈ [0.5,1], full
    // height) = 1 − 0.25 = 0.75.
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let through = m.add_cuboid(
        Point3::from_array([0.5, 0.5, -0.5]),
        Point3::from_array([1.5, 1.5, 1.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, base, through).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
}

#[test]
fn a_fused_stack_chains_through_a_cut() {
    // The dissolved 1×1×2 box (cell fuse-coplanar-merge) feeds a second boolean. Before
    // the merge/dissolve this rejected — first as COPLANAR_PAIR (the flat edges), then as
    // LOOP_ORIENT_MISMATCH (the straight-angle interface corners). A clean box cuts.
    let (mut m, a, b) = stacked_cubes();
    let stack = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    // A cutter straddling z=1 (the fused interface) — the seam runs where the split
    // vertical edges used to be. Result: 2 − 0.5·0.5·1.0.
    let cutter = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, stack, cutter).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.75).abs() < 1e-12, "volume {vol}");
}

/// **The topology of a boolean does not depend on whether the coordinates are
/// f64-representable.** The same shape is built twice — once on tidy integers, once on
/// dimensions that are not exact binary fractions — and both must give the same b-rep counts,
/// with each volume matching its own formula.
///
/// This is the invariant family #3 restored. Plane identity used to be read from the faces'
/// *derived* coefficients, which are not exactly proportional for two differently-sized faces
/// on one plane, so one plane became two classes and the arrangement named one point twice —
/// but only when the arithmetic did not happen to cancel, which tidy coordinates hid
/// (measured before the fix: 200/200 random stacked pairs under-merged, 12/600 ops aborted).
#[test]
fn boolean_topology_is_the_same_on_untidy_coordinates() {
    let counts = |dx: f64, dy: f64, z0: f64, h1: f64, h2: f64, kind: BoolKind| {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, z0]),
            Point3::from_array([dx, dy, z0 + h1]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, z0 + h1]),
            Point3::from_array([dx, dy, z0 + h1 + h2]),
        );
        let r = boolean_one(&mut m, kind, a, b).expect("stacked boxes fuse/cut");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let s = m.solids.get(r);
        let sh = m.shells.get(s.outer);
        let faces = sh.faces.len();
        let mut edges = std::collections::HashSet::new();
        let mut verts = std::collections::HashSet::new();
        for &fh in &sh.faces {
            let f = m.faces.get(fh);
            for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &l.half_edges {
                    edges.insert(he.edge);
                    if let Some(bd) = m.edges.get(he.edge).bounds {
                        verts.extend(bd);
                    }
                }
            }
        }
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        ((faces, edges.len(), verts.len(), s.cavities.len()), vol)
    };
    // The untidy dimensions are the minimal case the proptest shrank to when the kernel aborted.
    let (tidy_dx, tidy_dy, tidy_z0, tidy_h1, tidy_h2) = (2.0, 0.5, 0.0, 0.5, 2.0);
    let (dx, dy, z0, h1, h2) = (
        1.628165457453874,
        0.5,
        0.11200046228159026,
        0.5,
        2.07926124157585,
    );
    for kind in [BoolKind::Fuse, BoolKind::Cut] {
        let (tidy_shape, tidy_vol) = counts(tidy_dx, tidy_dy, tidy_z0, tidy_h1, tidy_h2, kind);
        let (shape, vol) = counts(dx, dy, z0, h1, h2, kind);
        assert_eq!(tidy_shape, shape, "{kind:?}: same topology either way");
        let want = |a: f64, b: f64| match kind {
            BoolKind::Fuse => a + b,
            _ => a,
        };
        let (tw, w) = (
            want(tidy_dx * tidy_dy * tidy_h1, tidy_dx * tidy_dy * tidy_h2),
            want(dx * dy * h1, dx * dy * h2),
        );
        assert!((tidy_vol - tw).abs() < 1e-9, "{kind:?} tidy {tidy_vol}");
        assert!((vol - w).abs() < 1e-9 * w, "{kind:?} untidy {vol}");
    }
}

#[test]
fn boolean_rejects_non_live_input() {
    let (mut m, a, b) = two_boxes();
    m.live_solids.retain(|&s| s != b); // as if superseded
    assert_eq!(
        boolean_one(&mut m, BoolKind::Common, a, b),
        Err(BoolError::InputNotLive)
    );
}

#[test]
fn boolean_op_applies_and_wraps_error() {
    // A failing boolean's error is surfaced as `OpError::Boolean`. This used to be driven by a
    // disjoint `Common`, but that is no longer an error (it is an empty result, see
    // `boolean_op_passes_an_empty_result_through`), so the wrapping is exercised with a boolean
    // that genuinely fails: a handle that is not live.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
    m.live_solids.retain(|&s| s != b); // retire `b` behind the op's back
    assert_eq!(
        apply(
            &mut m,
            &Operation::Boolean {
                kind: BoolKind::Common,
                a,
                b
            }
        ),
        Err(OpError::Boolean(BoolError::InputNotLive))
    );
}

/// An empty boolean reaches the caller as an empty solid list, not an error — the op layer
/// passes the kernel's answer through rather than reinterpreting it.
#[test]
fn boolean_op_passes_an_empty_result_through() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    // Offset in all axes so no faces are coplanar with A.
    let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
    let out = apply(
        &mut m,
        &Operation::Boolean {
            kind: BoolKind::Common,
            a,
            b,
        },
    )
    .unwrap();
    assert_eq!(out, OpOutput::Boolean { solids: vec![] });
    assert!(m.live_solids.is_empty(), "both operands are consumed");
}
