use super::*;

/// The canonical form is the sorted triple, and the constructor is the only way to get one —
/// so the same three planes named in any order are **one** name.
///
/// The expected value is written out rather than derived, because a test that re-derives it
/// through the constructor would agree with the constructor however the constructor behaved.
#[test]
fn three_planes_names_one_vertex_however_it_is_spelled() {
    let canonical = NodeId::ThreePlane([2, 5, 9]);
    for spelling in [
        [2, 5, 9],
        [2, 9, 5],
        [5, 2, 9],
        [5, 9, 2],
        [9, 2, 5],
        [9, 5, 2],
    ] {
        assert_eq!(
            NodeId::three_planes(Canon3::three(spelling)),
            canonical,
            "{spelling:?} names the same vertex as [2, 5, 9]"
        );
    }
    assert_eq!(
        three_plane_name(canonical),
        Some([2, 5, 9]),
        "and the door agrees"
    );
}

/// ★ **`Ord` is the bare triple's lexicographic order** — the proposition the whole migration
/// to this type stands on. Four rules read it and would answer differently if it moved:
/// `Aliases::union_point`'s "smallest name wins", `merge_component`'s sorted ring starts,
/// `reuse::canonical`'s ring rotation, and the reject witness's `Break::key`.
///
/// Each expected sign is written by hand from the pair, not computed from either operand, so
/// the test cannot agree with a wrong implementation by sharing its derivation.
#[test]
fn a_name_orders_like_the_triple_it_is() {
    use std::cmp::Ordering::{Equal, Greater, Less};
    for (a, b, want) in [
        ([0, 1, 2], [0, 1, 3], Less),    // last component decides
        ([0, 1, 3], [0, 1, 2], Greater), // …and antisymmetrically
        ([0, 2, 9], [0, 3, 4], Less),    // middle decides before last
        ([1, 0, 0], [0, 9, 9], Greater), // first decides before middle
        ([4, 4, 4], [4, 4, 4], Equal),
    ] {
        assert_eq!(
            NodeId::ThreePlane(a).cmp(&NodeId::ThreePlane(b)),
            want,
            "{a:?} vs {b:?}"
        );
    }
}

/// **A pierce point named either way round is one name** — the [`NodeId::three_planes`]
/// proposition for the second variant, where the canonicalization is bigger than a sort.
///
/// The expected values are written from the rule's sentence, not by calling the constructor a
/// second time. The rule itself is [`nacre_topo::QuadRoot::canonical`]'s, and its own locks
/// (in `nacre-topo`) are what say the *geometry* agrees; this only says the name reads it.
#[test]
fn a_pierce_is_one_name_however_the_pair_is_handed_in() {
    use nacre_topo::QuadRoot::{Double, Hi, Lo};
    for (root, other) in [(Lo, Hi), (Hi, Lo), (Double, Double)] {
        assert_eq!(
            NodeId::pierce(9, 2, 7, root),
            NodeId::pierce(2, 9, 7, other),
            "{root:?} handed in descending is {other:?} in stored order"
        );
    }
    // ★ And the two roots of one pair stay **two** names — a canonicalization that collapsed
    // them would make this pass by making everything equal.
    assert_ne!(NodeId::pierce(2, 9, 7, Lo), NodeId::pierce(2, 9, 7, Hi));
    // A different cylinder through the same two planes is a different point.
    assert_ne!(NodeId::pierce(2, 9, 7, Lo), NodeId::pierce(2, 9, 8, Lo));
}

/// ★ **The door tells the two variants apart, and `Ord` places them deterministically.**
///
/// The order itself is a *choice* — derived `Ord` puts every `ThreePlane` before every
/// `Pierce`, so "the smallest name wins" gains a systematic lean toward three-plane names at
/// the six rules that read it. Nothing can observe it yet (no ring holds a pierce node), which
/// is exactly why it is written down here rather than left to be discovered.
#[test]
fn the_two_variants_are_told_apart_and_ordered() {
    let three = NodeId::three_planes(Canon3::three([9, 2, 5]));
    let pierce = NodeId::pierce(2, 9, 7, nacre_topo::QuadRoot::Lo);
    assert_eq!(three_plane_name(three), Some([2, 5, 9]));
    assert_eq!(
        three_plane_name(pierce),
        None,
        "a pierce point has no triple"
    );
    assert!(
        three < pierce,
        "declaration order: ThreePlane before Pierce"
    );
    // The probe helper drops exactly the pierce node, and keeps the ring's order otherwise.
    assert_eq!(
        three_plane_probes([pierce, three, pierce]),
        vec![[2, 5, 9]],
        "a probe list may lose a member; that is its licence"
    );
}

/// ★★★★★ **The chart's shape, stated from the rule rather than from itself.**
///
/// One normal carries every clause: `n = 3e19 * (1, 2, 3)` is **wide** (a component of `9e19`,
/// whose square leaves `i128`, so an unrescaled `e2` is `None`) and **tilted** (`n_k != 0` for
/// the chosen `k`, without which two of the four clauses below cannot fail).
///
/// ⚠ **The expected axes are written out by hand from `e1 = ê_k × n`, `e2 = n × e1`, never by
/// calling `of_normal` a second time.** An earlier draft of this lock compared
/// `of_normal(K·n)` against `of_normal(n)` — which a chart that ignores magnitude entirely
/// (dropping a coordinate) passes on both sides while regressing the census by nine rows. An
/// oracle has to come from the inputs.
///
/// Cross-check on `e2`: `n × (ê₀ × n) = |n|²ê₀ − n₀n = 14(1,0,0) − (1,2,3) = (13, −2, −3)`,
/// the same vector the cross product gives — two derivations, one answer.
#[test]
fn the_chart_is_an_orthogonal_in_plane_frame_along_the_normals_own_directions() {
    use nacre_exact::{Orient, Rat, dot_sign_rat, parallel_rat};
    let q = Rat::from_int;
    let k = 30_000_000_000_000_000_000i128; // 3e19
    let n = [q(k), q(2 * k), q(3 * k)];
    let chart = Chart2dRat::of_normal(&n)
        .expect("a wide normal is a rescale away from a chart, not a refusal");
    let (e1, e2) = chart.axes();
    for (got, want, name) in [
        (e1, [q(0), q(-3), q(2)], "e1 = ê₀ × n"),
        (e2, [q(13), q(-2), q(-3)], "e2 = n × e1"),
    ] {
        assert!(
            parallel_rat(got, &want) && dot_sign_rat(got, &want) == Orient::Positive,
            "{name}: {got:?} is not along {want:?}"
        );
        assert_eq!(
            dot_sign_rat(got, &n),
            Orient::Zero,
            "{name} is a direction *in* the plane, and `ring_interior_candidates` walks it as one"
        );
    }
    assert_eq!(
        dot_sign_rat(e1, e2),
        Orient::Zero,
        "the frame is orthogonal (not orthonormal), which is what makes `axes`' sentence — \
             the parity walks its ray along e1 — true of the world and not just of the chart"
    );
}
