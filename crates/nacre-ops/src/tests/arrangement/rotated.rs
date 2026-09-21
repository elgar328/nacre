//! Rotation-generality: the trace engine's decisions are coordinate-free.

use super::*;

/// The boolean's error and the class audit's `failed_at` are the **same** reject.
///
/// They are two consumers of one `DeclineKind → RejectReason` mapping, and before
/// `decline_to_reject` they were two copies of it. A copy that drifts makes the audit — the
/// tool used to debug a reject — disagree with the reject being debugged, which is the worst
/// possible time to be lying.
///
/// The model is **chosen by measurement, not by taste**: a sweep over the fixture corpus found
/// no boolean that declines inside the arrangement at all (every fixture builds), so the
/// agreement had to be pinned on an input that still stops after the classes have run. If a
/// later capability makes it build, the fix is to take whatever still declines — not to weaken
/// the assertion.
#[test]
fn the_audit_does_not_invent_failures() {
    // The audit's scope is the per-class pipeline, and its duty is to run **the pipeline the
    // boolean runs** — with the alias fixpoint. Audited against an empty alias table it
    // reported `UnorderedEdges` for three classes of this input (names two classes discover
    // for each other were missing), failures the boolean never had: an instrument that
    // invents readings. The boolean's own reject here (`StraightAngle`) comes from the
    // assembly's vertex naming, after every class pipeline has run and outside the audit's
    // scope — so the audit's honest answer for this input is "no class failed".
    //
    // ★ **The fixture has moved twice, exactly as the note above prescribes**, and each move
    // is a capability the kernel gained. First it was the same pair at 30° with `Cut`,
    // rejecting `DegenerateWitness` — a reject from the component outwardness test, which the
    // nesting-parity label replaced. Then it was that pair at 60°, rejecting `StraightAngle`:
    // two unit cubes that only *touch*, which now come back as the two bodies they are
    // (measured: `Common` empty, fused volume 2.0, `validate` clean).
    //
    // So the fixture is now a pinch that **cannot** part: A and B meet only along the line
    // `x = 2, y = 2`, and a bridge overlapping both runs the material around the contact, so
    // cutting there leaves one piece. The closed-shell guard in `assembly::reconstruct` rejects
    // it — after every class pipeline has run, which is the property this test needs.
    //
    // (This test once asserted the audit reports the boolean's *class-level* reject, on a
    // fixture chosen as "some input that rejects" — a bar rotated through an L-shaped
    // target. The knife-edge fix turned that family, and every class-level-rejecting valid
    // input we could construct, into answers; the shared mapping the old test guarded,
    // `decline_to_reject`, is one function called by both consumers, so it cannot drift.)
    let build = || -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
            let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            s
        };
        let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
        let bridge = cub(&mut m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
        let b = cub(&mut m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
        let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
        m.rebuild_adjacency();
        (m, ab[0], b)
    };
    let (mut m, a, b) = build();
    let err = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::NonManifoldResultEdge,
                ..
            }
        ),
        "the fixture's premise: a reject from outside the class pipeline (got {err:?})"
    );
    let (m, a, b) = build();
    let audits = frame_audit(&m, BoolKind::Fuse, a, b).unwrap();
    let failed: Vec<RejectReason> = audits.iter().filter_map(|x| x.failed_at).collect();
    assert_eq!(
        failed,
        vec![],
        "every class runs clean under the boolean's own alias table — a failure invented \
             here is an artifact of running a different pipeline than the boolean runs"
    );
}

/// Rigidly rotating both operands (same single-Z tilt) leaves all three booleans' volumes
/// invariant — the axis values (hand-anchored by `end_to_end_overlapping_cubes_all_three`)
/// transfer to the rotated case. A coordinate-dependent decision that flipped under rotation
/// would add/drop a cell and move the volume by O(0.1), far past 1e-9.
#[test]
fn rotated_overlapping_cubes_all_three() {
    use nacre_exact::Axis;
    let vol_of = |kind: BoolKind| -> f64 {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let a = tilt(&mut m, a, &[Axis::Z]);
        let b = tilt(&mut m, b, &[Axis::Z]);
        let solids = boolean(&mut m, kind, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
        solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum()
    };
    assert!(
        (vol_of(BoolKind::Fuse) - 1.875).abs() < 1e-9,
        "rotated Fuse invariant 1.875"
    );
    assert!(
        (vol_of(BoolKind::Cut) - 0.875).abs() < 1e-9,
        "rotated Cut invariant 0.875"
    );
    assert!(
        (vol_of(BoolKind::Common) - 0.125).abs() < 1e-9,
        "rotated Common invariant 0.125"
    );
}

/// The through-tunnel Cut is invariant under every orientation: single Z, X, Y, and compound
/// Z∘X. The axis-DEPENDENCE was the tell of the bug (Y worked; Z/X declined before the
/// `Judge::orient3d` on-plane fix), so all four orientations returning 24 is the fix's direct
/// regression lock. The cube's own z=0/z=3 caps exercise the seated path under rotation.
#[test]
fn rotated_tunnel_cut_all_orientations() {
    use nacre_exact::Axis;
    let vol_of = |axes: &[Axis]| -> f64 {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let a = tilt(&mut m, a, axes);
        let b = tilt(&mut m, b, axes);
        let solids = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "manifold: {vs:?}");
        solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum()
    };
    for axes in [
        &[Axis::Z][..],
        &[Axis::X][..],
        &[Axis::Y][..],
        &[Axis::Z, Axis::X][..],
    ] {
        let v = vol_of(axes);
        assert!(
            (v - 24.0).abs() < 1e-9,
            "tunnel Cut vol 24 for {axes:?}, got {v}"
        );
    }
}

/// Compound-tilted (no face normal axis-aligned) tunnel Cut: AREA 64 is invariant (the bore
/// discriminator volume cannot see), and the pre-assembly face set is combinatorially identical
/// to the axis case (`tunnel_cut_emits_ten_faces_two_annular`) — 10 faces, 2 with an inner ring,
/// every undirected edge twice — proving the arrangement itself survived rotation.
#[test]
fn rotated_tunnel_area_and_faces() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let a = tilt(&mut m, a, &[Axis::Z, Axis::X]);
    let b = tilt(&mut m, b, &[Axis::Z, Axis::X]);

    // Combinatorial invariant (pre-assembly, A/B isolation).
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        plane_ix,
        class_owner,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let trace_in = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        Default::default(),
    );
    let (faces, _, _) = trace_result_faces(
        &m,
        BoolKind::Cut,
        a,
        b,
        &jd,
        &faces_tab,
        &plane_ix,
        &cyls,
        n_a,
        &class_owner,
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .unwrap();
    assert_eq!(faces.len(), 10, "rotated arrangement keeps 10 faces");
    assert_eq!(
        faces.iter().filter(|f| !f.inner.is_empty()).count(),
        2,
        "the two annular caps survive rotation"
    );
    let triples = |ns: &[combinatorics::NodeId]| -> Vec<[usize; 3]> {
        ns.iter()
            .map(|&n| three_plane_name(n).expect("a three-plane node"))
            .collect()
    };
    let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
    for f in &faces {
        for ring in f.poly_rings() {
            let ns = triples(ring);
            for w in ns
                .windows(2)
                .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
            {
                let key = if w[0] < w[1] {
                    (w[0], w[1])
                } else {
                    (w[1], w[0])
                };
                *count.entry(key).or_insert(0) += 1;
            }
        }
    }
    assert!(
        count.values().all(|&c| c == 2),
        "every edge used exactly twice: {:?}",
        count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
    );

    // Area (assembled).
    let solids = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let area: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().area)
        .sum();
    assert!(
        (area - 64.0).abs() < 1e-9,
        "rotated Cut area 64 (bore survives rotation), got {area}"
    );
}

/// The sharp tripwire: a compound-tilted tunnel must decline nothing, exactly like the axis
/// baseline `axis_aligned_cubes_decline_nothing`. This is the EXACT site the bug broke
/// (`coincident-features` from a `wall == wc` degenerate), so it fails immediately if the
/// `Judge::orient3d` on-plane fix is reverted.
#[test]
fn rotated_tunnel_declines_nothing() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let a = tilt(&mut m, a, &[Axis::Z, Axis::X]);
    let b = tilt(&mut m, b, &[Axis::Z, Axis::X]);
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    for wc in 0..planes.len() {
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        assert!(
            tr.declined.is_empty(),
            "rotated class {wc} declined: {tr:?}"
        );
    }
}
