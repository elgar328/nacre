//! The folded map against the arithmetic it folds — one node at a time, then whole chains.

use super::*;
use proptest::prelude::*;

fn r(n: i128, d: i128) -> Rat {
    Rat::new(n, d).expect("a rational")
}

/// One chain node, as the existing per-motion arithmetic applies it — the independent oracle.
#[derive(Clone, Copy, Debug)]
enum Node {
    Translate([Rat; 3]),
    Rotate(Rotation),
    Mirror(Axis, Rat),
}

impl Node {
    fn fold(self) -> Option<AxisAffine> {
        match self {
            Node::Translate(t) => Some(AxisAffine::translation(t)),
            Node::Rotate(rot) => AxisAffine::rotation(rot),
            Node::Mirror(a, o) => Some(AxisAffine::mirror(a, o)),
        }
    }
    fn point(self, p: [Rat; 3]) -> Option<[Rat; 3]> {
        match self {
            Node::Translate(t) => Isometry::translation(t).point_rat(p),
            Node::Rotate(rot) => Isometry::rotation(rot).point_rat(p),
            Node::Mirror(a, o) => mirror_point_rat(p, a, o),
        }
    }
    fn plane(self, c: [Rat; 4]) -> Option<[Rat; 4]> {
        match self {
            Node::Translate(t) => Isometry::translation(t).plane_coeffs(c),
            Node::Rotate(rot) => Isometry::rotation(rot).plane_coeffs(c),
            Node::Mirror(a, o) => mirror_plane_coeffs(c, a, o),
        }
    }
    fn dir(self, d: [Rat; 3]) -> Option<[Rat; 3]> {
        match self {
            Node::Translate(_) => Some(d),
            Node::Rotate(rot) => Isometry::rotation(rot).dir_rat(d),
            Node::Mirror(a, _) => {
                let mut out = d;
                out[a.index()] = Rat::from_int(0).checked_sub(d[a.index()])?;
                Some(out)
            }
        }
    }
}

fn axis_of(k: u8) -> Axis {
    [Axis::X, Axis::Y, Axis::Z][k as usize % 3]
}

fn small() -> impl Strategy<Value = Rat> {
    (-12i128..=12, 1i128..=6).prop_map(|(n, d)| r(n, d))
}

fn point() -> impl Strategy<Value = [Rat; 3]> {
    proptest::array::uniform3(small())
}

/// A node that folds: a translation, a quarter turn about a rational pivot, or a reflection.
fn node() -> impl Strategy<Value = Node> {
    prop_oneof![
        point().prop_map(Node::Translate),
        (0u8..3, point(), 0i128..4).prop_map(|(a, pivot, q)| Node::Rotate(Rotation {
            axis: axis_of(a),
            pivot,
            angle: Angle::from_deg(Rat::from_int(90 * q)).expect("a quarter"),
        })),
        (0u8..3, small()).prop_map(|(a, o)| Node::Mirror(axis_of(a), o)),
    ]
}

/// An irrational turn has no rational map, and says so.
#[test]
fn a_turn_the_rationals_cannot_state_does_not_fold() {
    let rot = Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(37)).expect("an angle"),
    };
    assert_eq!(AxisAffine::rotation(rot), None);
}

/// A reflection is improper and a turn is not; composing multiplies.
#[test]
fn a_reflection_is_the_only_improper_fold() {
    let turn = AxisAffine::rotation(Rotation {
        axis: Axis::X,
        pivot: [r(1, 2); 3],
        angle: Angle::from_deg(Rat::from_int(90)).expect("an angle"),
    })
    .expect("a quarter turn folds");
    let flip = AxisAffine::mirror(Axis::Y, r(3, 4));
    assert_eq!(turn.det(), 1);
    assert_eq!(flip.det(), -1);
    assert_eq!(turn.then(&flip).det(), -1);
    assert_eq!(flip.then(&flip).det(), 1);
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// One node folded answers what that node's own arithmetic answers — points, planes and
    /// directions alike.
    #[test]
    fn one_node_folds_to_its_own_arithmetic(
        n in node(),
        p in point(),
        c in proptest::array::uniform4(small()),
    ) {
        let f = n.fold().expect("every generated node folds");
        prop_assert_eq!(f.point_rat(p), n.point(p));
        prop_assert_eq!(f.dir_rat(p), n.dir(p));
        prop_assume!(c[0] != Rat::from_int(0) || c[1] != Rat::from_int(0) || c[2] != Rat::from_int(0));
        prop_assert_eq!(f.plane_coeffs(c), n.plane(c));
    }

    /// A chain folded root first answers what applying its nodes in order answers.
    #[test]
    fn a_folded_chain_is_its_nodes_in_order(
        chain in proptest::collection::vec(node(), 1..6),
        p in point(),
        c in proptest::array::uniform4(small()),
    ) {
        let folded = chain
            .iter()
            .skip(1)
            .fold(chain[0].fold().expect("folds"), |acc, n| acc.then(&n.fold().expect("folds")));
        let stepped = chain.iter().try_fold(p, |q, n| n.point(q));
        prop_assert_eq!(folded.point_rat(p), stepped);
        let stepped_dir = chain.iter().try_fold(p, |q, n| n.dir(q));
        prop_assert_eq!(folded.dir_rat(p), stepped_dir);
        prop_assume!(c[0] != Rat::from_int(0) || c[1] != Rat::from_int(0) || c[2] != Rat::from_int(0));
        let stepped_plane = chain.iter().try_fold(c, |q, n| n.plane(q));
        prop_assert_eq!(folded.plane_coeffs(c), stepped_plane);
        let dir_f = folded.dir_f64(p.map(|x| x.to_f64()));
        if let Some(d) = folded.dir_rat(p) {
            prop_assert_eq!(dir_f, d.map(|x| x.to_f64()));
        }
    }

    /// A chain carried back by its inverse returns every point and direction where it started.
    #[test]
    fn a_folded_chain_undoes_itself(
        chain in proptest::collection::vec(node(), 1..6),
        p in point(),
    ) {
        let folded = chain
            .iter()
            .skip(1)
            .fold(chain[0].fold().expect("folds"), |acc, n| acc.then(&n.fold().expect("folds")));
        let back = folded.inverse();
        if let Some(q) = folded.point_rat(p) {
            prop_assert_eq!(back.point_rat(q), Some(p));
        }
        if let Some(d) = folded.dir_rat(p) {
            prop_assert_eq!(back.dir_rat(d), Some(p));
        }
        prop_assert_eq!(back.det(), folded.det());
    }

    /// A point on a plane lands on the plane's image — the two maps describe one motion.
    #[test]
    fn a_point_on_a_plane_stays_on_its_image(
        chain in proptest::collection::vec(node(), 1..5),
        p in point(),
        n in point(),
    ) {
        prop_assume!(n.iter().any(|x| *x != Rat::from_int(0)));
        let d = (0..3).try_fold(Rat::from_int(0), |acc, k| acc.checked_sub(n[k].checked_mul(p[k])?));
        let Some(d) = d else { return Ok(()) };
        let folded = chain
            .iter()
            .skip(1)
            .fold(chain[0].fold().expect("folds"), |acc, n| acc.then(&n.fold().expect("folds")));
        if let (Some(q), Some(c)) = (folded.point_rat(p), folded.plane_coeffs([n[0], n[1], n[2], d])) {
            let v = (0..3).try_fold(c[3], |acc, k| acc.checked_add(c[k].checked_mul(q[k])?));
            prop_assert_eq!(v, Some(Rat::from_int(0)));
        }
    }
}
