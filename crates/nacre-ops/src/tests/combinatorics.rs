use super::*;

/// The canonical form is the sorted triple, and the constructor is the only way to get one —
/// so the same three planes named in any order are **one** name.
///
/// The expected value is written out rather than derived, because a test that re-derives it
/// through the constructor would agree with the constructor however the constructor behaved.
#[test]
fn three_planes_names_one_vertex_however_it_is_spelled() {
    let canonical = NodeId(NodeKind::ThreePlane([2, 5, 9]));
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
/// to this type stands on. Rules read it that would answer differently if it moved:
/// `Aliases::union_point`'s "smallest name wins", `merge_component`'s sorted ring starts, and
/// `reuse::canonical`'s ring rotation.
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
            NodeId(NodeKind::ThreePlane(a)).cmp(&NodeId(NodeKind::ThreePlane(b))),
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

/// ★★ **A direction sign reads the truth, never the plane cache.** A prism over
/// `(0,0)·(1, 10⁻¹⁷)·(1,1)·(0,1)`: its wall from the origin runs along `(1, 10⁻¹⁷)`, so the
/// plane's world name is `(10⁻¹⁷·k, −k, 0, 0)` canonicalized — a first nonzero component of
/// `10⁻¹⁷` against the largest. The class's cache is then replaced by a rounded image of the same
/// plane with that component's sign flipped: an error of `2·10⁻¹⁷` against a unit normal, inside
/// one ulp of the large component, so nothing reading the cache as a rounding could object.
///
/// The class's facing ([`stored_coeffs_rat`]) and the canonical → outward sign
/// ([`loops::outward_fix`]) come out the same: they read the world name and its sense. Comparing
/// the name's first nonzero component with the cache's — the spelling `outward_fix` and two
/// siblings had — reads the flipped sign, which the contrast below states.
#[test]
fn a_direction_sign_does_not_read_the_plane_cache() {
    use crate::{OpOutput, Operation, SketchFrame, apply, from_rings};
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, nacre_exact::Axis::Z);
    let profile = from_rings(vec![vec![
        p2(0.0, 0.0),
        p2(1.0, 1e-17),
        p2(1.0, 1.0),
        p2(0.0, 1.0),
    ]])
    .expect("a quadrilateral")
    .remove(0);
    let Ok(OpOutput::Extrude { solid, .. }) = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    ) else {
        panic!("the prism")
    };
    m.rebuild_adjacency();
    let faces = crate::planes::collect_planes(&m, solid).expect("the face table");
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces));
    let (mut planes, _, _) = crate::planes::dense_planes(&faces, &canon);
    let zero = nacre_exact::Rat::from_int(0);
    let c = planes
        .iter()
        .position(|p| p.world_rat().is_some_and(|w| w[0] != zero && w[1] != zero))
        .expect("the slanted wall");
    let w = planes[c].world_rat().expect("a world name");
    assert!(
        w[0].to_f64().abs() < 1e-15 * w[1].to_f64().abs(),
        "the name's first nonzero component is tiny: {w:?}"
    );
    let answers = |planes: &[WorkingPlane]| {
        let jd = crate::planes::test_judge(planes);
        (stored_coeffs_rat(&jd, c), loops::outward_fix(&jd, c))
    };
    let truth = answers(&planes);
    assert!(truth.0.is_some() && truth.1.is_some(), "both answer");

    let n = planes[c].plane.normal().as_array();
    assert!(
        n[0] != 0.0 && n[0].abs() < 1e-15,
        "the cache carries the tiny component: {n:?}"
    );
    let lying = nacre_geom::Plane::from_point_normal(
        planes[c].plane.origin(),
        nacre_math::Vector3::from_array([-n[0], n[1], n[2]]),
    )
    .expect("a plane");
    planes[c].plane = lying;
    // The contrast: the first nonzero component of the name, compared with the cache's, now reads
    // the other way — the spelling this test retires would have flipped.
    assert_ne!(
        (w[0] > zero) == (lying.normal().as_array()[0] > 0.0),
        (w[0] > zero) == (n[0] > 0.0),
        "the planted cache flips the first-component reading"
    );
    assert_eq!(answers(&planes), truth, "the direction signs read no cache");
}

/// **An arc's span holds the seam point like any other** — `arc_span` against an oracle derived
/// apart from it: a point `P` of a circle is on the counter-clockwise arc `A → B` iff it lies to the
/// right of the chord, `orient2d(A, B, P) < 0`, viewed with the axis toward the viewer.
///
/// The circle is radius 5 about `+z` through the origin, and its twelve points with integer
/// coordinates (`(±5, 0)`, `(0, ±5)`, `(±3, ±4)`, `(±4, ±3)`) are every arc's ends and every
/// root. The seam is put on one of them three ways — and, the fourth, between them — so ends and
/// roots on the seam, arcs across it and arcs beside it all occur: every ordered pair of distinct
/// ends against every point, `AtLo`/`AtHi` at the ends.
#[test]
fn an_arc_span_reads_the_seam_point_as_any_other() {
    use nacre_exact::quad::{CylinderMeet, QuadVal, plane_plane_cylinder};
    use nacre_exact::{BigRat, Orient, Rat};
    let ri = Rat::from_int;
    let mut pts: Vec<(i128, i128)> = vec![(5, 0), (0, 5), (-5, 0), (0, -5)];
    for (a, b) in [(3, 4), (4, 3)] {
        for (sa, sb) in [(1, 1), (1, -1), (-1, 1), (-1, -1)] {
            pts.push((sa * a, sb * b));
        }
    }
    let r2 = BigRat::from(ri(25));
    for seam in [[1, 0], [0, -1], [3, 4], [1, 2]] {
        let def = nacre_topo::CylinderDef::new(
            [ri(0); 3],
            [ri(0), ri(0), ri(1)],
            [ri(seam[0]), ri(seam[1]), ri(0)],
            r2.clone(),
        )
        .expect("a cylinder");
        // The point `(x, y)` as the meet of `x = px` and `z = 0` and the root on its side of `y`.
        let at = |(px, py): (i128, i128)| {
            let meet = plane_plane_cylinder(
                &[ri(1), ri(0), ri(0), ri(-px)],
                &[ri(0), ri(0), ri(1), ri(0)],
                &def.origin(),
                &def.dir(),
                def.r2(),
            )
            .expect("fits");
            let (line, roots) = match meet {
                CylinderMeet::Pair { line, s } => (line, s.to_vec()),
                CylinderMeet::Tangent { line, s } => (line, vec![QuadVal::from_rat(s)]),
                other => panic!("({px}, {py}) is on the circle: {other:?}"),
            };
            let y_of = |s: &QuadVal| {
                QuadVal::from_rat(line.base()[1])
                    .checked_add(&s.checked_mul_rat(line.dir()[1]).unwrap())
                    .unwrap()
            };
            let s = *roots
                .iter()
                .find(|s| {
                    y_of(s)
                        .checked_sub(&QuadVal::from_rat(ri(py)))
                        .unwrap()
                        .sign()
                        == Orient::Zero
                })
                .expect("the root at the point");
            (line, s)
        };
        let mut decided = 0usize;
        for &a in &pts {
            for &b in &pts {
                if a == b {
                    continue;
                }
                for &p in &pts {
                    let want = if p == a {
                        ArcSpan::AtLo
                    } else if p == b {
                        ArcSpan::AtHi
                    } else {
                        let o = (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0);
                        if o < 0 {
                            ArcSpan::Inside
                        } else {
                            ArcSpan::Outside
                        }
                    };
                    let got = arc_span(&def, &at(a), &at(b), &at(p));
                    assert_eq!(
                        got,
                        Some(want),
                        "seam {seam:?}: arc {a:?} → {b:?}, point {p:?}"
                    );
                    decided += 1;
                }
            }
        }
        assert_eq!(decided, 12 * 11 * 12, "seam {seam:?}");
    }
}
