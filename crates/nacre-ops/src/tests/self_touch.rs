//! Unit tests of the tangency verdict's own functions. (`tests/suite/tangent_lumps.rs` drives the
//! verdict end to end through booleans.)

use super::*;
use crate::combinatorics::QuadRoot;
use crate::draft::{Bound, Ring};
use crate::planes::ClassIx;

/// **An edge on the tangent line counts only where the wall's faces and the lateral both end on
/// it** ([`edges_on_line`]). Cylinder `0` touches plane class `1` along its ruling between the
/// cylinder's meets with plane classes `2` and `3`; two wall faces end on that line from either
/// side. With the lateral cut there too, the two wall edges count. With the lateral running on
/// across the line — a whole band between two rims, no edge on it — neither does, whatever the
/// wall's faces hold. A wall edge on another cylinder's ruling is not this line's.
#[test]
fn a_line_is_an_edge_only_where_the_wall_and_the_lateral_both_end_on_it() {
    let (lo, hi) = (
        NodeId::pierce(1, 2, 0, QuadRoot::Double),
        NodeId::pierce(1, 3, 0, QuadRoot::Double),
    );
    let corner =
        |k: usize| NodeId::pierce(4, 5, 7, if k == 0 { QuadRoot::Lo } else { QuadRoot::Hi });
    let ruling = |cyl: usize, up: bool, plane: usize| Wall::Ruling {
        cyl,
        side: 0,
        up,
        plane,
    };
    let face = |surf, nodes: Vec<NodeId>, walls: Vec<Wall>| LocalFace {
        surf,
        outer: Bound::Ring(Ring::new(nodes, walls)),
        inner: Vec::new(),
        flip: false,
    };
    let walls = || {
        [
            face(
                ClassIx::Plane(1),
                vec![lo, hi, corner(0)],
                vec![ruling(0, true, 1), Wall::Plane(3), Wall::Plane(2)],
            ),
            face(
                ClassIx::Plane(1),
                vec![hi, lo, corner(1)],
                vec![ruling(0, false, 1), Wall::Plane(2), Wall::Plane(3)],
            ),
        ]
    };
    let cut = face(
        ClassIx::Cyl(0),
        vec![lo, hi, corner(0)],
        vec![ruling(0, true, 1), Wall::Plane(3), Wall::Plane(2)],
    );
    let whole = LocalFace {
        surf: ClassIx::Cyl(0),
        outer: Bound::Rim {
            plane: 2,
            ccw: true,
        },
        inner: vec![Bound::Rim {
            plane: 3,
            ccw: false,
        }],
        flip: false,
    };

    let with_cut: Vec<LocalFace> = walls().into_iter().chain([cut]).collect();
    let mut got = edges_on_line(&with_cut, 1, 0);
    got.sort();
    let mut want = vec![[lo, hi], [hi, lo]];
    want.sort();
    assert_eq!(got, want);

    let with_whole: Vec<LocalFace> = walls().into_iter().chain([whole]).collect();
    assert!(edges_on_line(&with_whole, 1, 0).is_empty());

    assert!(edges_on_line(&with_cut, 1, 1).is_empty());
}
