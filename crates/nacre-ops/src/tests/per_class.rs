//! Unit tests of the per-class stage's own functions.

use super::*;
use nacre_exact::Rat;

/// **The net cannot measure an oblique section its pieces do not part, and says so.** The bored
/// plate with a blind bore along `x` — radius ½ about `y = 10, z = 5`, `x ∈ [13, 40]` — and the sunk
/// gusset, whose slanted plane `z = y/4 + 5/2` sections the `z` bore and runs through the `x` bore's
/// axis: the class carries the first's ellipse and the second's rulings. With the ruling's own ends
/// the two pieces part along `x` (`[8, 12]` against `[13, 40]`) and the net answers «no crossing».
/// With its ends moved out to the plate's wall `x = 0` they overlap along every direction the pair
/// rule asks — the two axes meet at `(10, 10, 5)` — and the strip's arc ends are a circle's, so the
/// answer is «not measured», whether the circle states no contribution or one whole one: refused as
/// the oblique coverage's end, not as the net's backstop.
#[test]
fn an_oblique_section_its_pieces_do_not_part_is_not_measured() {
    let mut m = Model::new();
    let plate = crate::tests::plate_with_bore(&mut m);
    let across = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([13.0, 10.0, 5.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        28.0,
    )
    .solid;
    let a = crate::boolean(&mut m, BoolKind::Cut, plate, across).expect("the blind bore")[0];
    let b = crate::tests::gusset(&mut m, crate::tests::SUNK, 30.0);
    let setup = crate::arrangement::plane_index_setup(&m, a, b).expect("the gate passes");
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let z = |v: i128| Rat::from_int(v);
    let class = |want: [Rat; 4]| {
        (0..setup.geom.len())
            .find(|&c| {
                jd.world_rat(c).is_some_and(|w| {
                    // The same plane up to a positive or negative scale.
                    let k = (0..4).find(|&i| want[i] != z(0)).unwrap();
                    let s = w[k]
                        .checked_mul(Rat::new(want[k].denom(), want[k].numer()).unwrap())
                        .unwrap();
                    (0..4).all(|i| w[i] == want[i].checked_mul(s).unwrap())
                })
            })
            .expect("the class is in the table")
    };
    let wc = class([
        z(0),
        Rat::new(-1, 4).unwrap(),
        z(1),
        Rat::new(-5, 2).unwrap(),
    ]);
    let (x0, x13, x40) = (
        class([z(1), z(0), z(0), z(0)]),
        class([z(1), z(0), z(0), z(-13)]),
        class([z(1), z(0), z(0), z(-40)]),
    );
    let cyl_along = |dir: [Rat; 3]| {
        setup
            .cyls
            .iter()
            .position(|cy| nacre_exact::parallel_rat(&cy.def.dir(), &dir))
            .expect("the cylinder is in the table")
    };
    let (upright, lying) = (cyl_along([z(0), z(0), z(1)]), cyl_along([z(1), z(0), z(0)]));
    assert!(
        setup.cyls[upright].sections_oblique(wc),
        "the gate records the section"
    );
    assert!(
        setup.cyls[lying].cuts(wc),
        "and the cut along the other axis"
    );
    let end = |c: usize| {
        crate::combinatorics::NodeId::pierce(wc, c, lying, crate::combinatorics::QuadRoot::Lo)
    };
    let ruling = |from: usize| MergedRuling {
        cyl: lying,
        side: 1,
        end: [end(from), end(x40)],
        merged: Vec::new(),
        orient: -1,
    };
    let whole = (
        crate::planes::SolidSide::A,
        SegKind::Transversal { mat: -1 },
        None,
    );
    for merged in [Vec::new(), vec![whole]] {
        let circle = MergedCircle {
            cyl: upright,
            merged,
        };
        assert_eq!(
            circle_crosses_ruling(&jd, &setup.cyls, wc, &circle, &ruling(x13)),
            Ok(false),
            "parted along x"
        );
        assert_eq!(
            circle_crosses_ruling(&jd, &setup.cyls, wc, &circle, &ruling(x0)),
            Err(RejectReason::ObliqueCylinderCut),
            "not parted, and an ellipse is not measured by a circle's strip — the oblique coverage's end"
        );
    }
}

/// A knife of `solid` whose two faces leave the plane toward `above`.
fn knife_of(
    solid: crate::planes::SolidSide,
    above: bool,
) -> Vec<(crate::planes::SolidSide, SegKind)> {
    vec![(solid, SegKind::Graze { body_above: above }); 2]
}

/// **A knife can part the material only where it can lie, and never under `Fuse`.** Its solid has
/// no face on the plane along it, so that solid's two bits in the label are equal; over every such
/// label, a one-solid knife can part the material under `Cut` (either solid's) and `Common`, and
/// under `Fuse` cannot — the one kept chamber beside a face keeps the wedge too. A knife two solids
/// graze at once is left whole. red: every label asked, the solid's two bits free → the label
/// `[1, 0, 0, 0]`, which [`knife_parts`]'s premise refuses (and which parts the material under
/// `Fuse`).
#[test]
fn a_knife_can_part_the_material_except_under_fuse() {
    use crate::planes::SolidSide::{A, B};
    for (solid, above) in [(A, true), (A, false), (B, true), (B, false)] {
        let m = knife_of(solid, above);
        assert!(!knife_can_part(BoolKind::Fuse, &m), "{solid:?} {above}");
        assert!(knife_can_part(BoolKind::Cut, &m), "{solid:?} {above}");
        assert!(knife_can_part(BoolKind::Common, &m), "{solid:?} {above}");
    }
    let mut two = knife_of(A, true);
    two.extend(knife_of(B, true));
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        assert!(!knife_can_part(kind, &two), "{kind:?}");
    }
}

/// A class split by hand: the skeleton is the squares named by `rings` (corners `4r..4r+4`), the
/// knives run between the named nodes. Only the ends and the contributions are read.
fn split_by_hand(rings: usize, knives: &[[usize; 2]]) -> ClassSplit {
    let node = |n: usize| NodeId::three_planes(Canon3::three([3 * n, 3 * n + 1, 3 * n + 2]));
    let piece = |[a, b]: [usize; 2], merged| Piece {
        on: OnCurve::Segment {
            wall: 0,
            end_h: [crate::combinatorics::EndPin::Class(0); 2],
            sense: None,
        },
        end: [node(a), node(b)],
        merged,
    };
    ClassSplit {
        pieces: (0..rings)
            .flat_map(|r| (0..4).map(move |i| [4 * r + i, 4 * r + (i + 1) % 4]))
            .map(|e| piece(e, Vec::new()))
            .collect(),
        circles: Vec::new(),
        cut_rims: Vec::new(),
        knives: knives
            .iter()
            .map(|&e| piece(e, knife_of(crate::planes::SolidSide::A, true)))
            .collect(),
    }
}

/// **A knife is an edge only where both its ends rest on edges, and it joins no two parts of the
/// skeleton.** A square `0-1-2-3` and a chain of two knives `0–9–2` across it: both are edges
/// together, neither alone (the other's end is free). A knife from the square to a second square
/// `4-5-6-7` joins the two and is no edge.
#[test]
fn a_knife_is_an_edge_where_it_ends_on_edges_and_bridges_nothing() {
    let chain = split_by_hand(1, &[[0, 9], [9, 2]]);
    assert_eq!(knives_as_edges(&chain, vec![0, 1]), vec![0, 1]);
    assert_eq!(knives_as_edges(&chain, vec![0]), Vec::<usize>::new());
    let bridge = split_by_hand(2, &[[0, 4]]);
    assert_eq!(knives_as_edges(&bridge, vec![0]), Vec::<usize>::new());
}

/// **The second walk carries the knives that part and can still be edges.** On the chain `0–9–2`
/// across a square, both knives part: the first walk stands. Only `0–9` parts: alone it ends at `9`
/// inside the face, so no knife is walked. red: the parting knives walked again without asking
/// whether they are edges → `Again([0])`, a dangling knife.
#[test]
fn a_parting_knife_left_with_a_free_end_is_not_walked_again() {
    let chain = split_by_hand(1, &[[0, 9], [9, 2]]);
    assert_eq!(rewalk(&chain, &[0, 1], &[true, true]), Rewalk::Kept);
    assert_eq!(rewalk(&chain, &[0, 1], &[true, false]), Rewalk::Plain);
    assert_eq!(rewalk(&chain, &[0, 1], &[false, false]), Rewalk::Plain);
    let two = split_by_hand(1, &[[0, 2], [1, 3]]);
    assert_eq!(
        rewalk(&two, &[0, 1], &[false, true]),
        Rewalk::Again(vec![1])
    );
}
