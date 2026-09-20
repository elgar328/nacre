//! Contacts, arcs and chained booleans: slits, arc classes, subdivided twins, the band's loop,
//! contact cuts.

use super::*;

/// A minted vertex sits within `1e-12` of its definition — asserted from the cache where the
/// cache proves it, and **measured here** where it does not.
///
/// ★ The second half used to read a stored residual. The cache no longer stores one, and the
/// honest replacement is not "skip those" but "measure the same distance the residual was":
/// the point against the surfaces its definition names. Dropping to a bare `panic!` for the
/// unproven variants would have turned this from a measurement into a restatement of which
/// variant the funnel chose.
fn knowledge_is_tight(m: &Model, h: Handle<nacre_topo::Vertex>) {
    match m.vertex_cache(h) {
        nacre_topo::PointCache::Bounded { bound, .. } => assert!(
            bound
                .iter()
                .all(|b| nacre_exact::Mag::lt(*b, nacre_exact::Mag::of(1e-12))),
            "realized within rounding: {bound:?}"
        ),
        nacre_topo::PointCache::Ceiling { coord }
        | nacre_topo::PointCache::Unrealized { coord } => {
            for sh in m.vertex(h).carriers() {
                let d = m.surface_cache(sh).distance(*coord);
                assert!(
                    d < 1e-12,
                    "unrealized, but {d:e} from surface {}",
                    sh.index()
                );
            }
        }
    }
}

/// **The tangency is a *slit*, not a cut — and that is why no combinatorial check can see it.**
/// Two clauses, because either alone passes vacuously:
///
/// 1. **No vertex within `1e-9` of the touch.** The arrangement names the point exactly
///    (`NodeId::pierce(.., QuadRoot::Double)`) and then *discards* it — a tangency touches
///    without separating, so the circle keeps its closed cell. `nonmanifold_vertices` and the
///    Euler count read topology, and there is none here: **their silence is not evidence**
///    about this point, in either direction.
/// 2. **The circle that touches is still whole** — a closed `[v, v]` rim edge passing through
///    the point. Without this the first clause cannot tell "the circle was never cut" from
///    "the circle was cut and its vertex landed elsewhere".
///
/// What *is* the evidence that the solid is sound: the link of the boundary at the touch is a
/// single circle (the pinched face's two lobes are joined around through the neighbouring
/// curved face), the material is locally one piece, and a second kernel returns a body of the
/// same volume and area (`nacre-oracle`'s `a_segment_tangent_to_a_rim_is_a_body_to_occt` and
/// `a_rim_tangent_to_a_plate_top_is_a_body_to_occt`, measured against an ε-twin whose tangency
/// is broken). ⇒ **the surface is a 2-manifold there; only the face is pinched.**
fn the_touch_is_a_slit(m: &Model, s: Handle<Solid>, at: [f64; 3]) {
    let p = Point3::from_array(at);
    let sol = m.solid(s).clone();
    let (mut vertices_at, mut whole_circle_through) = (0usize, 0usize);
    let mut seen = std::collections::HashSet::new();
    for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
        for &fh in &m.shell(sh).faces {
            let face = m.face(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    if !seen.insert(he.edge) {
                        continue;
                    }
                    let e = m.edge(he.edge);
                    for &vh in e.vertices.iter() {
                        if m.vertex_point(vh).distance(p) <= 1e-9 {
                            vertices_at += 1;
                        }
                    }
                    // A whole circle is the closed `[v, v]` rim spelling; a cut one is arcs.
                    if e.vertices[0] == e.vertices[1] && m.edge_curve(he.edge).distance(p) <= 1e-9 {
                        whole_circle_through += 1;
                    }
                }
            }
        }
    }
    assert_eq!(vertices_at, 0, "a tangency mints no vertex: {at:?}");
    assert!(
        whole_circle_through > 0,
        "the touching circle is still whole (a closed rim edge through {at:?})"
    );
}

/// **A tangency builds.** The boss's `x = 11` edge is exactly tangent to the rim `(8,10)`,
/// `r = 3`, so the locator's quadratic has a double root — the third arm of `CylinderMeet`
/// that can reach here, and the only one whose point carries no radical (`(11, 10, 5)`,
/// rational).
///
/// ★★ **And a touch is not a break**: the circle keeps its closed cell, so nothing about the
/// arrangement had to change for this to work. It was refused only because the guard's
/// sentence was wider than its proposition — measured before the narrowing landed, with the
/// whole guard off, this already produced exactly the solid asserted below.
#[test]
fn a_segment_tangent_to_the_rim_builds() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = m.add_cylinder(
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        7.0,
    );
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    let boss = m.add_cuboid(
        Point3::from_array([11.0, 8.0, 5.0]),
        Point3::from_array([15.0, 12.0, 8.0]),
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("a tangent edge");
    assert_eq!(out.len(), 1, "the boss sits on material");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 40.0 * 20.0 * 5.0 - std::f64::consts::PI * 9.0 * 5.0 + 4.0 * 4.0 * 3.0;
    assert!(
        (v - want).abs() < 1e-9,
        "plate less bore plus boss: {v} vs {want}"
    );
    // The through bore survives the fuse, so the one shell still has one handle.
    let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
    assert_eq!(
        v_n - e_n + f_n - l_n,
        0,
        "genus 1: V{v_n} E{e_n} F{f_n} L{l_n}"
    );
    // ★ **This χ is the genus's lock, not the tangency's.** A pinch *that has a vertex* makes
    // χ odd — what `check_result_topology` reads — but a tangency mints no vertex at all
    // (below), so an even χ is compatible with a pinch and with none: it decides nothing
    // here. Deleting this assertion would lose the genus, not the manifold claim.
    m.rebuild_adjacency();
    the_touch_is_a_slit(&m, out[0], [11.0, 10.0, 5.0]);
    // ★★★★★ **And the solid meshes — across a bridge.** The boss's footprint
    // touches the bore's rim at exactly one point, so the plate's top face has two inner
    // loops meeting there: its interior is pinched. `SelfTouchingBoundary` here would be
    // a statement about the *tessellator* (`validate` is clean above and the volume is
    // exact). Instead the touching sample is put into the straight edge it lies on, the two
    // holes are spliced into one at that point, and the sweep orders the twins symbolically.
    // ★ Rebuilt first, on purpose: `boolean`'s own census meshes *before* the rebuild, and a
    // census that only agreed with itself would be measuring when it looks rather than what
    // came out. Both spellings say the same thing here.
    m.rebuild_adjacency();
    crate::tests::mesh_covers_faces("segment tangent to the rim", &m, &out);
}

/// ★★ **The fence post: a crossing that lands exactly on a segment's endpoint.** The boss's
/// corner `(8,13)` sits on the rim `(8,10)`, `r = 3`, so both edges leaving it meet the circle
/// *at* their own end — the one place the in-segment test's inequality can be open or closed,
/// and the two spellings give different answers.
///
/// **Measured, both ways:** counting `Zero` as inside names the corner `(8, 13, 5)`; excluding
/// it drops every root and falls back to the containment witness — a `Segment`, which would
/// say "this edge lies inside the circle" about an edge that runs *outward* from a single
/// touching point. Inclusive is the true sentence, and it is the one the arc split will need.
#[test]
fn a_crossing_on_a_segments_endpoint_is_inside_it() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = m.add_cylinder(
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        3.0,
        7.0,
    );
    m.rebuild_adjacency();
    let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    let boss = m.add_cuboid(
        Point3::from_array([8.0, 13.0, 5.0]),
        Point3::from_array([12.0, 17.0, 8.0]),
    );
    m.rebuild_adjacency();
    let err = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect_err("corner on rim");
    // ★★ **One point, two names** — and the split says so by name. `(8,13,5)` is the boss's
    // corner, so the arrangement already holds it as a three-plane vertex; the rim crossing
    // names the *same* point as a `NodeId::Pierce`. The DCEL keys vertices by name, so shipping
    // both would put two vertices where there is one — folding them is its own step, and until
    // then `CoincidentNodes` ("two names for one point") is the true sentence.
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::CoincidentNodes,
                ..
            }
        ),
        "{err:?}"
    );
}

/// ★★★ **The population where the naming rule actually fires — and the first non-`+Z`
/// cylinder in this repository.**
///
/// Every `add_cylinder*` call in the suite states the axis `[0, 0, 1]` (measured: 87 of them),
/// and that is not an accident of taste — `add_cuboid` pushes faces `[−Z, +Z, −Y, +Y, −X, +X]`
/// and `build_prism` pushes "base cap, top cap, then walls", so a `+Z` cylinder's circle always
/// lives on a class interned **before** the wall that crosses it. Measured over the whole ops
/// suite: 463 of the 2060 `(class, wall)` pairs the gate examines are descending, but **all
/// nine that reach the locator are ascending**. So `NodeId::pierce`'s canonicalization —
/// re-sorting the pair and restating the root with it — has never once been exercised by a
/// production path.
///
/// Turning the boss onto `+X` inverts it: the rim lives on the box's `x = 4` face (class 5,
/// interned last) and the segment it crosses is on the `z = 2` cap (class 1). The pair swaps.
///
/// ★ **A tangency, so there is exactly one root** — and that also means the
/// circle is not separated, so this **builds** rather than refusing. The single root is still
/// what makes the naming rule visible here (it is what the fixture was written for), and the
/// witness it once asserted is now the *split point that never happens*. That matters: the
/// straddling fixtures deliberately
/// accept either of their two crossings, and a fixture copied from that template would be green
/// whether or not the root rule is right. Here `disc = 0` — the rim (`x = 4`, centre
/// `(4, 2, 1.5)`, `r = 0.5`) touches `z = 2` at the single point `(4, 2, 2)`, solved from the
/// fixture's own numbers — and the name it must take is `QuadRoot::Double`, which a swap must
/// **not** toggle.
#[test]
fn a_turned_boss_tangent_to_the_plate_top_builds() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array([4.0, 2.0, 1.5]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        0.5,
        1.0,
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("a turned boss");
    assert_eq!(
        out.len(),
        1,
        "the boss's base disk sits on the plate's face"
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 32.0 + std::f64::consts::PI * 0.25;
    assert!(
        (v - want).abs() < 1e-9,
        "plate plus the whole boss: {v} vs {want}"
    );
    let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
    assert_eq!(
        v_n - e_n + f_n - l_n,
        2,
        "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
    );
    // ★ The genus's lock, not the tangency's — see the note in the fixture above.
    m.rebuild_adjacency();
    the_touch_is_a_slit(&m, out[0], [4.0, 2.0, 2.0]);
    // ★★★★★ **And the solid meshes — the same pinch, spelled inner-to-outer.** The boss's
    // base circle is tangent to the plate's `z = 2` edge at `(4, 2, 2)`, so the `x = 4` face's
    // hole touches its own outer ring at one point and the face's interior is pinched there.
    // `validate` is clean above and the volume is exact; the touch is bridged
    // and the face's triangles are held to its exact area.
    // ★ Rebuilt first, on purpose: `boolean`'s own census meshes *before* the rebuild, and a
    // census that only agreed with itself would be measuring when it looks rather than what
    // came out. Both spellings say the same thing here.
    m.rebuild_adjacency();
    crate::tests::mesh_covers_faces("turned boss tangent to the plate top", &m, &out);
}

/// **The turned boss builds — and its crescent is the winding witness's red switch.**
///
/// The near cap's circle is cut by two different plate edges (top and corner), so this
/// result carries the arc-dominated crescent face whose chord Newell reads **backwards**
/// (`cos = −1`, the measured wall) — `validate == []` here is what pins `loop_winding`'s
/// segment witnesses. ★ The old fence's proposition — a root that fails to follow its pair
/// through the sort names the wrong crossing — did not retire with the reject: the mint
/// fence asserts the pierce vertices sit **on the derived crossings**, and a wrong root
/// moves the minted point itself.
#[test]
fn a_turned_boss_over_the_plates_corner_builds() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array([4.0, 0.25, 2.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        0.5,
        1.0,
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
    assert_eq!(out.len(), 1, "one fused solid");
    m.rebuild_adjacency();
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 32.0 + std::f64::consts::PI * 0.25;
    assert!((v - want).abs() < 1e-12, "{v} vs {want}");
    let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
    assert_eq!(
        v_n - e_n + f_n - l_n,
        2,
        "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
    );
}

/// **The audit and the boolean say the same thing about an arc class — and it says what the
/// arrangement produced.**
///
/// ★★★★ Three things are pinned here and none had a fence before.
///
/// **One: the two copies of the per-class pipeline agree.** `frame_audit` re-runs pass A and B
/// inline rather than calling `arrange`, so a split hoisted into one and not the other compiles
/// perfectly and makes the audit report *no failure* for an input the boolean refuses —
/// `decline_to_reject`'s doc calls that "the worst possible time to be lying". The existing
/// guard (`the_audit_does_not_invent_failures`) cannot see it: its fixture is deliberately one
/// that fails **outside** the class pipeline.
///
/// **Two: what the arrangement produced on an arc class.** The stopper stands after every
/// arrangement stage and swallows their answer, so without this the whole arc population is
/// measured by a probe and the probe is deleted before the commit. The numbers are **derived,
/// not read back**: the
/// boss's rim crosses one edge of the plate's top face twice, cutting the circle into **2**
/// arcs and the plane into **4** cells — `rect−disk`, `rect∩disk`, `disk−rect` and the outside.
/// Only the outside winds `−1`, so there is **1** root group; and it shares nodes with all
/// three `+1` cells (the plate's corners, the two crossings), so `cell_in_cell` answers "not
/// comparable" every time and there are **0** holes. Exactly **one** class is cut: the boss
/// stands *on* the top face, so only its base circle lies in a plane of the plate.
///
/// **Three: the faces it emitted, as coordinates** — and this is the one that sees a wrong
/// `sense`. Flipping the sense the split carries onto its sub-segments attaches the arcs to
/// the wrong cells, and **every number in part Two stays put** (the cells still read four with
/// one `−1`, the roots one, the holes zero, the sorted labels identical) — that blindness was
/// derived before it was measured, and it held for four rungs. The emitted ring is where it
/// finally shows, as the **exact reverse**, which is why rotation is normalized here and
/// reversal is not. See `ClassAudit::outer_rings` for both arguments and for why the red probe
/// is read on the **turned** boss: the straddling one is blind to the flip *and* its 2-node
/// ring cannot express a reversal at all.
///
/// ★ **Two earlier claims here were wrong and are recorded rather than quietly dropped**: that
/// this fence is *green* under that flip (it is red — measured, three times), and that the
/// flip's first reader is `label_cells`' keep decision (refuted a rung earlier: every
/// order-independent summary of the labels is identical; it is `emit_faces`).
#[test]
fn the_audit_and_the_boolean_agree_about_an_arc_class() {
    // `y = 0` meets the turned boss's circle (centre `(y,z) = (0.25, 2)`, `r = 0.5`) at
    // `(z - 2)² = 0.1875`. The one irrational coordinate in either fixture, and the reason the
    // ring comparison is `near()` rather than `==`.
    let s = 2.0 - 3.0f64.sqrt() / 4.0;
    for (origin, axis, n_out, want_rings) in [
        (
            [4.0, 2.0, 2.0],
            [0.0, 0.0, 1.0],
            // The plate's top face carries this class, so the ring is CCW about `+z`.
            [0.0, 0.0, 1.0],
            vec![
                // `rect - disk`: the plate's top, three corners and the inner arc.
                vec![
                    [0.0, 0.0, 2.0],
                    [4.0, 0.0, 2.0],
                    [4.0, 1.5, 2.0],
                    [4.0, 2.5, 2.0],
                    [4.0, 4.0, 2.0],
                    [0.0, 4.0, 2.0],
                ],
                // `disk - rect`: the overhang's underside — chord and outer arc.
                vec![[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
            ],
        ),
        (
            [4.0, 0.25, 2.0],
            [1.0, 0.0, 0.0],
            // The plate's `x = 4` face carries this one, so CCW is about `+x`.
            [1.0, 0.0, 0.0],
            vec![
                // `rect - disk`: the plate's right face, the disk biting its top-left corner.
                vec![
                    [4.0, 0.0, 0.0],
                    [4.0, 4.0, 0.0],
                    [4.0, 4.0, 2.0],
                    [4.0, 0.75, 2.0],
                    [4.0, 0.0, s],
                ],
                // `disk - rect`: the boss's base cap outside the plate. The corner `(4,0,2)` is
                // `0.25` from the circle's centre, so it sits *inside* the disk and is one of
                // this ring's three nodes.
                vec![[4.0, 0.0, s], [4.0, 0.0, 2.0], [4.0, 0.75, 2.0]],
            ],
        ),
    ] {
        // The boolean takes `&mut Model` and the audit `&Model`; a rejected boolean leaves
        // arena residue, so each gets its own build of the same fixture.
        let build = || {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            (m, plate, boss)
        };
        let (mut m, plate, boss) = build();
        crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the arc class assembles");
        let (m, plate, boss) = build();
        let audits = crate::arrangement::frame_audit(&m, BoolKind::Fuse, plate, boss).unwrap();
        // The agreement, now that the population is green: the boolean built, and the audit's
        // replay of the same pipeline stops nowhere either.
        let stopped: Vec<_> = audits.iter().filter(|a| a.failed_at.is_some()).collect();
        assert!(stopped.is_empty(), "no class stops: {stopped:?}");
        // ★ The whole audit, not just `produced`: the ring coordinates are a sibling field, and
        // joining two filtered lists by index is the seam where they could come from different
        // classes. ★ By reference: `Produced` carries the labels now, so it is no longer `Copy`.
        let cut: Vec<_> = audits
            .iter()
            .filter(|a| a.produced.as_ref().is_some_and(|p| p.arcs > 0))
            .collect();
        assert_eq!(cut.len(), 1, "exactly one class has its circle cut");

        assert_eq!(
            *cut[0].produced.as_ref().unwrap(),
            crate::arrangement::Produced {
                cells: 4,
                arcs: 2,
                roots: 1,
                holes: 0,
                // ★★ **Derived from the geometry, not read back.** `[A_above, A_below,
                // B_above, B_below]` for the three `+1` regions, sorted: the boss stands *on*
                // the plate, so plate-minus-disk has the plate below it and nothing above;
                // disk-minus-plate — the overhang's underside — has the boss above it and
                // nothing below; and their intersection has both. The outside cell winds `-1`
                // and is not here.
                //
                // ★ **The same three come out for the turned boss, and the story is the same
                // one in its own frame** — this assertion is inside the two-fixture loop, so
                // reading it as a sentence about the straddling boss alone would be reading
                // half of what it checks. There the class is `x = 4`: the plate lies on the
                // `−x` side, the boss's base cap outside it on `+x`, and the overlap has both.
                pos_labels: vec![
                    [false, false, true, false],
                    [false, true, false, false],
                    [false, true, true, false],
                ],
            },
            "the arrangement walked the arcs, nested them, and labelled them"
        );
        // ★★★★★ **And the faces it emitted, as coordinates.** This is what a wrong `sense` on a
        // split segment moves and nothing before it does: the cell counts, the nesting and the
        // sorted labels above are all identical whichever way the rings run.
        //
        // ★★★ **The premise first.** `emit_faces` emits CCW about `n_out(wc)`, and which face
        // is the class root — hence which way `n_out` points — is a plane-table fact, not one
        // of the fixture's numbers. Pinning it here means a root flip fails *this* assertion,
        // with its own sentence, instead of silently reversing every ring below. Both come
        // from the plane table, which no arrangement stage can move.
        let got_n_out: Vec<f64> = cut[0]
            .root_normal
            .iter()
            .map(|c| c * cut[0].orient_sign as f64)
            .collect();
        assert!(
            got_n_out
                .iter()
                .zip(&n_out)
                .all(|(a, b)| (a - b).abs() < 1e-9),
            "the rings below are derived CCW about {n_out:?}, but the class faces {got_n_out:?}"
        );
        let got = cut[0]
            .outer_rings
            .as_ref()
            .expect("the per-class road emits the rings");
        assert_eq!(
            got.len(),
            want_rings.len(),
            "one outer ring per emitted face"
        );
        for (g, w) in got.iter().zip(&want_rings) {
            assert_eq!(g.len(), w.len(), "ring length: got {g:?}, want {w:?}");
            assert!(
                g.iter()
                    .zip(w)
                    .all(|(p, q)| (0..3).all(|i| (p[i] - q[i]).abs() < 1e-9)),
                "ring: got {g:?}, want {w:?}"
            );
        }
    }
}

/// **The seam realizes a pierce vertex and measures it — seen through the door production
/// uses, because nothing else can see it at all.**
///
/// ★★★ The deferred stopper intercepts the whole seam stretch, so with the pierce arm
/// deleted the seam fails `PierceVertexUnnamed`, the interception swallows it, and **every
/// boolean-level fence stays green** (measured — that probe is what forced this test). The
/// arm's only witness is a direct second consumer of `seam_table` on the very faces
/// production feeds it, so this walks production's stretch step for step: trace, clean,
/// append the bands, build the seam.
///
/// ★ The tolerance is asserted as a **bound**, never a copied value; the pierce coordinates
/// themselves are pinned by `ClassAudit::outer_rings` through the same `pierce_point` road,
/// so re-asserting them here would be a second copy of an existing lock — and a `Lo`/`Hi`
/// mix-up cannot hide behind the bound either, because both crossings lie on every defining
/// surface and `outer_rings` is what tells them apart.
#[test]
fn the_seam_realizes_a_pierce_vertex_and_measures_it() {
    for (origin, axis) in [
        ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
        ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let setup = plane_index_setup(&m, plate, boss).unwrap();
        let PlaneSetup {
            planes: faces_tab,
            geom,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
            cyls,
            ..
        } = &setup;
        let jd = Judge::new(geom, *standard, notes);
        let trace_in = crate::combinatorics::trace_input(
            &m,
            [(plate, inc_a), (boss, inc_b)],
            surf_ix,
            faces_tab.len(),
            &jd,
            plane_ix,
            cyls,
            Default::default(),
        );
        let (plane_faces, curved, deferred) = crate::arrangement::trace_result_faces_full_for_test(
            &m,
            BoolKind::Fuse,
            plate,
            boss,
            &jd,
            faces_tab,
            plane_ix,
            cyls,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the arc population traces");
        // The stopper socket is empty since the population went green — nothing defers.
        assert!(deferred.is_none(), "{deferred:?}");
        let faces =
            crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
        let mut faces = faces;
        faces.extend(
            crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                .expect("the lateral faces emit"),
        );
        let seam = crate::arrangement::seam_table(&faces, cyls, &jd)
            .expect("the seam realizes pierce nodes");
        let pierce: Vec<_> = seam
            .iter()
            .filter(|sv| crate::combinatorics::pierce_name(sv.triple).is_some())
            .collect();
        assert_eq!(pierce.len(), 2, "both crossings reach the seam, once each");
        for sv in pierce {
            assert!(
                sv.tol < 1e-12,
                "a pierce realization sits on everything that defines it: tol {}",
                sv.tol
            );
        }
    }
}

/// **Every ring node of the arc population's result faces earns a definition — measured
/// through the door production uses.**
///
/// ★★★ Same instrument shape as the seam fence, same reason: the deferred stopper stands
/// behind `name_result_vertices` (at the assembly's very end now) and intercepts
/// everything, so no reject name can testify the naming completed — a direct second consumer
/// is the only witness. This walks production's road (trace → clean → bands → seam → naming)
/// and asserts on its product.
///
/// ★★ Disabling the pierce arm reddens both fixtures (and the walls-fallback cannot fake a
/// pierce def past its arc-carrier guard). ★ The fallback itself went **zero-population** when
/// the split-twin subdivision landed — the bitten corner's twins match now, so its def comes
/// down the far-plane road and no probe reddens on the fallback alone; the subdivision has
/// its own fence (`a_subdivided_twin_matches_its_neighbour_edge_for_edge`).
#[test]
fn every_result_vertex_of_the_arc_population_is_named() {
    for (origin, axis, bites_corner) in [
        ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0], false),
        ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0], true),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let setup = plane_index_setup(&m, plate, boss).unwrap();
        let PlaneSetup {
            planes: faces_tab,
            geom,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
            cyls,
            ..
        } = &setup;
        let jd = Judge::new(geom, *standard, notes);
        let trace_in = crate::combinatorics::trace_input(
            &m,
            [(plate, inc_a), (boss, inc_b)],
            surf_ix,
            faces_tab.len(),
            &jd,
            plane_ix,
            cyls,
            Default::default(),
        );
        let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
            &m,
            BoolKind::Fuse,
            plate,
            boss,
            &jd,
            faces_tab,
            plane_ix,
            cyls,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the arc population traces");
        let faces =
            crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
        let mut faces = faces;
        faces.extend(
            crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                .expect("the lateral faces emit"),
        );
        let seam = crate::arrangement::seam_table(&faces, cyls, &jd).expect("seam");
        let named =
            crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                .expect("the naming stages run");
        let live: &[LocalFace] = named.per_solid.as_deref().unwrap_or(&faces);
        // Completeness: every ring node has a definition.
        let mut missing = Vec::new();
        for (fi, lf) in live.iter().enumerate() {
            for ring in lf.poly_rings() {
                for &n in ring.iter() {
                    if !named.defs.contains_key(&(named.group_of[fi], n)) {
                        missing.push(n);
                    }
                }
            }
        }
        assert!(missing.is_empty(), "def-less nodes: {missing:?}");
        let pierce_defs = named
            .defs
            .values()
            .filter(|d| matches!(d, crate::boolean::Def::Pierce { .. }))
            .count();
        assert_eq!(pierce_defs, 2, "both crossings are declared, once each");
        // The bitten corner's def names the right point: realize its three planes and land
        // on (4, 0, 2) — the fixture's own number, no class index copied.
        if bites_corner {
            let hit = named.defs.values().any(|d| {
                let crate::boolean::Def::Three(t) = d else {
                    return false;
                };
                let t = t.planes();
                nacre_geom::intersect::three_planes(
                    &geom[t[0]].plane,
                    &geom[t[1]].plane,
                    &geom[t[2]].plane,
                )
                .is_some_and(|p| {
                    let c = p.as_array();
                    (0..3).all(|i| (c[i] - [4.0, 0.0, 2.0][i]).abs() < 1e-9)
                })
            });
            assert!(hit, "the bitten corner's def realizes to (4, 0, 2)");
        }
    }
}

/// **After the split-twin subdivision, every segment edge has exactly one twin.**
///
/// ★★★ The subdivision cannot be locked through the naming fence — the walls-fallback
/// rescues the bitten corner whether or not the twins match, so that fence is green either
/// way. What the subdivision actually changes is the **edge-key census**: a neighbour's whole
/// edge and the arc class's subdivided pieces share `norm_edge` keys only once the whole edge
/// is cut at the same pierce nodes. So the proposition, with its one exception stated:
///
/// > every segment ring edge (a `Wall::Plane` carrier) has its key used by exactly two
/// > faces.
///
/// Arc edges are excluded because their far side is the band, which contributes no ring.
/// ★ Both-ends-pierce keys used to be excluded too — the chord and the two arcs between one
/// pierce pair folded into a single `norm_edge` key — but the carrier gave arcs their own
/// ordered key, so the chord's line key counts exactly its two coplanar faces now (measured:
/// the exclusion removed, both fixtures stay green).
///
/// red: with the subdivision disabled, the pre-split single-use keys come back.
#[test]
fn a_subdivided_twin_matches_its_neighbour_edge_for_edge() {
    for (origin, axis) in [
        ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
        ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let setup = plane_index_setup(&m, plate, boss).unwrap();
        let PlaneSetup {
            planes: faces_tab,
            geom,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
            cyls,
            ..
        } = &setup;
        let jd = Judge::new(geom, *standard, notes);
        let trace_in = crate::combinatorics::trace_input(
            &m,
            [(plate, inc_a), (boss, inc_b)],
            surf_ix,
            faces_tab.len(),
            &jd,
            plane_ix,
            cyls,
            Default::default(),
        );
        let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
            &m,
            BoolKind::Fuse,
            plate,
            boss,
            &jd,
            faces_tab,
            plane_ix,
            cyls,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the arc population traces");
        let faces =
            crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
        let mut faces = faces;
        faces.extend(
            crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                .expect("the lateral faces emit"),
        );
        let seam = crate::arrangement::seam_table(&faces, cyls, &jd).expect("seam");
        let named =
            crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                .expect("the naming stages run");
        let live: &[LocalFace] = named.per_solid.as_deref().unwrap_or(&faces);
        let mut uses: std::collections::HashMap<
            (crate::combinatorics::NodeId, crate::combinatorics::NodeId),
            usize,
        > = std::collections::HashMap::new();
        let mut plain: Vec<(crate::combinatorics::NodeId, crate::combinatorics::NodeId)> =
            Vec::new();
        for lf in live {
            for ring in lf.poly_rings() {
                let k = ring.nodes.len();
                for t in 0..k {
                    let (a, b) = (ring.nodes[t], ring.nodes[(t + 1) % k]);
                    if matches!(ring.walls[t], crate::combinatorics::Wall::Arc { .. }) {
                        continue;
                    }
                    let key = crate::boolean::norm_edge(a, b);
                    *uses.entry(key).or_insert(0) += 1;
                    if !plain.contains(&key) {
                        plain.push(key);
                    }
                }
            }
        }
        let odd: Vec<_> = plain
            .iter()
            .filter(|k| uses[*k] != 2)
            .map(|k| (*k, uses[k]))
            .collect();
        assert!(odd.is_empty(), "edges without exactly one twin: {odd:?}");
    }
}

/// **The pierce vertices are minted — canonical, measured, on the derived crossings.**
///
/// Before the boolean no `Vertex::Pierce` exists anywhere in the model, so a whole-store
/// filter is position-independent; a wrong `QuadRoot` canonicalization moves the minted
/// point itself, which is what keeps the old toggle-lock alive now that the reject (whose
/// witness once carried it) is gone.
///
/// ★ The coordinates are the fixtures' own crossing derivations (the same numbers the ring
/// and seam fences pin) — nothing here is copied from a run. The tolerance is a bound, and
/// ascending handle order is `Vertex::Pierce`'s own contract, minted through
/// `QuadRoot::canonical`'s second answer.
#[test]
fn a_pierce_vertex_is_minted_and_measured() {
    let s = 2.0 - 3.0f64.sqrt() / 4.0;
    for (origin, axis, crossings) in [
        (
            [4.0, 2.0, 2.0],
            [0.0, 0.0, 1.0],
            [[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
        ),
        (
            [4.0, 0.25, 2.0],
            [1.0, 0.0, 0.0],
            [[4.0, 0.75, 2.0], [4.0, 0.0, s]],
        ),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let before = m.live_solids().to_vec();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
            .expect("the cut-rim boolean builds");
        assert_eq!(out.len(), 1, "one fused solid");
        assert_ne!(m.live_solids().to_vec(), before, "the operands retired");
        assert_eq!(m.live_solids().to_vec(), out, "the result lives");
        let pierce: Vec<_> = (0..m.vertex_count() as u32)
            .filter_map(|i| m.vertex_handle_at(i))
            .map(|h| (h, m.vertex(h)))
            .filter(|(_, v)| matches!(**v, nacre_topo::Vertex::Pierce { .. }))
            .collect();
        assert_eq!(pierce.len(), 2, "both crossings minted, once each");
        for (h, v) in pierce {
            let nacre_topo::Vertex::Pierce { planes: [a, b], .. } = *v else {
                unreachable!("filtered above");
            };
            assert!(a < b, "planes in ascending handle order: {a:?} vs {b:?}");
            let p = m.vertex_point(h);
            assert!(
                crossings
                    .iter()
                    .any(|c| (0..3).all(|i| (p.as_array()[i] - c[i]).abs() < 1e-9)),
                "a minted pierce vertex sits on a derived crossing: {p:?}"
            );
            knowledge_is_tight(&m, h);
        }
    }
}

/// **The two complementary arcs are two edges, and the cut rim is none.**
///
/// ★★★ The observable form of the `[A, B]`-CCW convention (`derive_edge_curve`'s circle
/// arm): between one pair of pierce vertices a circle offers two pieces, the endpoints
/// alone cannot tell them apart, and the *vertex order* is the bit that does — so the store
/// must hold **two** circle-carrier edges whose vertex pairs are each other's reverse.
/// Erase the order from the welding key and they fold into one edge (the red probe this
/// fence was built against).
///
/// ★★ The rim skip's two sides, on one store: the **cut** circle mints no closed `[v, v]`
/// edge (before the skip, both fixtures minted one that only the reject discarded), while
/// the **uncut** far rim still mints exactly one — the skip's negative control, pinned to
/// the far cap's axis coordinate so a skip that turned into "skip every rim" reddens here.
///
/// ★ The two populations differ on the chord, deliberately: the straddling boss's pierce
/// pair is joined by the plate-top chord (welded with the subdivided middle piece into
/// **one** line edge used by both coplanar faces — the subdivision's promise realized
/// in the store), while the turned boss's pair sits across the plate corner, joined through
/// it by split boundary edges — no chord at all. Scoped to the edges the boolean minted
/// (a snapshot, not a whole-store filter: the input cylinder's own rims are `[v, v]` too).
#[test]
fn an_arc_and_its_complement_are_minted_as_two_ordered_edges() {
    let s = 2.0 - 3.0f64.sqrt() / 4.0;
    // `seam_split`: where θ = 0 sits. `None` = a pierce vertex lies on the seam generator
    // (the straddling boss — the split's own `SeamIncident` case), so no piece splits;
    // `Some(p)` = the seam vertex S is minted at `p` (centre + ref_dir·r, derived) and the
    // wrap arc is cut there into two pieces.
    for (origin, axis, height, crossings, chord_edges, seam_split) in [
        (
            [4.0, 2.0, 2.0],
            [0.0, 0.0, 1.0],
            1.0,
            [[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
            1usize,
            None,
        ),
        (
            [4.0, 0.25, 2.0],
            [1.0, 0.0, 0.0],
            1.0,
            [[4.0, 0.75, 2.0], [4.0, 0.0, s]],
            0usize,
            Some([4.0, 0.25, 1.5]),
        ),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            height,
        );
        m.rebuild_adjacency();
        let minted_from = m.edge_count();
        let vertices_from = m.vertex_count();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
            .expect("the cut-rim boolean builds");
        assert_eq!(out.len(), 1, "one fused solid");
        let minted: Vec<_> = (minted_from as u32..m.edge_count() as u32)
            .filter_map(|i| m.edge_handle_at(i))
            .map(|h| (h, m.edge(h)))
            .collect();
        let is_cyl = |sh| matches!(m.surface_cache(sh), nacre_geom::Surface::Cylinder(_));
        // A circle-carrier edge is a **mixed** pair (the cap plane and the lateral); the
        // band's seam edge is `[lat, lat]` — both carriers the cylinder — and is not a
        // piece of any circle.
        let on_circle = |e: &nacre_topo::Edge| is_cyl(e.surfaces[0]) != is_cyl(e.surfaces[1]);

        // The cut circle's pieces: their directed vertex pairs chain into **one** cycle —
        // the observable of the `[A, B]`-CCW convention (for two pieces the cycle *is* the
        // mutually-reversed pair the first version of this fence asserted), and of the seam
        // split (three pieces when S stands, still one circle).
        let arcs: Vec<_> = minted
            .iter()
            .filter(|(_, e)| on_circle(e) && e.vertices[0] != e.vertices[1])
            .collect();
        let expect_pieces = 2 + usize::from(seam_split.is_some());
        assert_eq!(
            arcs.len(),
            expect_pieces,
            "the cut circle's pieces: {arcs:?}"
        );
        let mut succ: std::collections::HashMap<_, _> = std::collections::HashMap::new();
        for (_, e) in &arcs {
            let prev = succ.insert(e.vertices[0], e.vertices[1]);
            assert!(prev.is_none(), "two pieces leave one vertex CCW: {arcs:?}");
        }
        let start = arcs[0].1.vertices[0];
        let (mut cur, mut steps) = (start, 0usize);
        loop {
            cur = succ[&cur];
            steps += 1;
            if cur == start || steps > expect_pieces {
                break;
            }
        }
        assert_eq!(steps, expect_pieces, "the pieces close one circle");
        // Endpoints: exactly two pierce vertices on the derived crossings, plus — when the
        // seam splits an arc — one OnSeam vertex at the derived seam point, with its
        // tolerance measured.
        let (mut pierce, mut on_seam) = (Vec::new(), Vec::new());
        for &v in succ.keys() {
            match *m.vertex(v) {
                nacre_topo::Vertex::Pierce { .. } => pierce.push(v),
                nacre_topo::Vertex::OnSeam(_) => on_seam.push(v),
                ref d => panic!("an arc endpoint is neither pierce nor seam: {d:?}"),
            }
        }
        assert_eq!(pierce.len(), 2, "one pierce pair");
        for &v in &pierce {
            let p = m.vertex_point(v);
            assert!(
                crossings
                    .iter()
                    .any(|c| (0..3).all(|i| (p.as_array()[i] - c[i]).abs() < 1e-9)),
                "an arc endpoint sits on a derived crossing: {p:?}"
            );
        }
        match seam_split {
            None => assert!(on_seam.is_empty(), "the seam is the pierce vertex itself"),
            Some(sp) => {
                assert_eq!(on_seam.len(), 1, "one seam vertex on the cut circle");
                let p = m.vertex_point(on_seam[0]);
                assert!(
                    (0..3).all(|i| (p.as_array()[i] - sp[i]).abs() < 1e-12),
                    "S sits on the derived seam point: {p:?}"
                );
                knowledge_is_tight(&m, on_seam[0]);
            }
        }
        // The minted OnSeam census: the uncut far rim's vertex, plus S when it stands —
        // and nothing else (a duplicate S at a seam-incident pierce vertex would show here).
        let minted_on_seam = (0..m.vertex_count() as u32)
            .filter_map(|i| m.vertex_handle_at(i))
            .map(|h| (h, m.vertex(h)))
            .skip(vertices_from)
            .filter(|(_, v)| matches!(**v, nacre_topo::Vertex::OnSeam(_)))
            .count();
        assert_eq!(
            minted_on_seam,
            1 + on_seam.len(),
            "far rim + S, nothing else"
        );

        // The rim skip's two sides: no closed edge on the cut circle, exactly one on the
        // uncut far rim (its axis coordinate is the far cap's, derived from the fixture).
        let closed: Vec<_> = minted
            .iter()
            .filter(|(_, e)| e.vertices[0] == e.vertices[1])
            .collect();
        assert_eq!(closed.len(), 1, "one uncut rim, no cut one: {closed:?}");
        let far = (0..3)
            .map(|i| (origin[i] + axis[i] * height) * axis[i])
            .sum::<f64>();
        let p = m.vertex_point(closed[0].1.vertices[0]).as_array();
        let along = (0..3).map(|i| p[i] * axis[i]).sum::<f64>();
        assert!(
            (along - far).abs() < 1e-12,
            "the closed rim is the far cap's: {along} vs {far}"
        );

        // The chord: welded into one line edge on the straddling boss, absent across the
        // turned boss's corner.
        let (ba, bb) = (pierce[0], pierce[1]);
        let chords = minted
            .iter()
            .filter(|(_, e)| !on_circle(e) && (e.vertices == [ba, bb] || e.vertices == [bb, ba]))
            .count();
        assert_eq!(chords, chord_edges, "the pierce pair's line edges");
    }
}

/// **The band's loop is one continuous cycle, and the shell closes over it.**
///
/// ★★★ The band assembles its cut rim from the arc chain; this fence counts the closure
/// directly on the store beside the production guard: every edge the boolean minted is used
/// exactly twice across its faces, and the band face's outer loop is one vertex-continuous
/// cycle of the derived length, with the seam edge traversed once in each sense.
///
/// ★ Three fixtures: the straddling boss (lo rim cut, seam ≡ pierce), the turned boss
/// (lo rim cut, seam splits the wrap arc — six half-edges), and the **hung** boss (the
/// straddling boss mirrored under the plate — the cut circle is the band's **hi** end, so
/// the chain is walked reversed; measured to pass the gate before this fence was written).
/// The hi-cut *and* seam-split combination has no fixture yet — the chain logic is shared,
/// and its population brings one when it arrives.
#[test]
fn the_bands_loop_is_one_continuous_cycle() {
    for (origin, axis, band_len) in [
        ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0], 5usize),
        ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0], 6usize),
        ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0], 5usize),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let faces_from = m.face_count();
        let before = m.live_solids().to_vec();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
            .expect("the cut-rim boolean builds");
        assert_eq!(out.len(), 1, "one fused solid");
        assert_ne!(m.live_solids().to_vec(), before, "the operands retired");
        assert_eq!(m.live_solids().to_vec(), out, "the result lives");
        let garbage: Vec<_> = (faces_from as u32..m.face_count() as u32)
            .filter_map(|i| m.face_handle_at(i))
            .map(|h| (h, m.face(h)))
            .collect();
        assert!(!garbage.is_empty(), "the face loop ran to completion");

        // Closure: every edge of the garbage faces is used exactly twice.
        let mut uses: std::collections::HashMap<_, usize> = std::collections::HashMap::new();
        for (_, f) in &garbage {
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    *uses.entry(he.edge).or_default() += 1;
                }
            }
        }
        let odd: Vec<_> = uses.iter().filter(|&(_, &n)| n != 2).collect();
        assert!(
            odd.is_empty(),
            "a closed shell uses every edge twice: {odd:?}"
        );

        // The band face: one, on the cylinder, its outer loop a continuous cycle of the
        // derived length, the seam edge once in each sense.
        let bands: Vec<_> = garbage
            .iter()
            .filter(|(_, f)| matches!(m.surface_cache(f.surface), nacre_geom::Surface::Cylinder(_)))
            .collect();
        assert_eq!(bands.len(), 1, "one band face");
        let lp = &bands[0].1.outer;
        assert_eq!(lp.half_edges.len(), band_len, "the derived loop length");
        let ends = |he: &nacre_topo::HalfEdge| {
            let [a, b] = m.edge(he.edge).vertices;
            if he.forward { (a, b) } else { (b, a) }
        };
        for w in 0..lp.half_edges.len() {
            let (_, e0) = ends(&lp.half_edges[w]);
            let (s1, _) = ends(&lp.half_edges[(w + 1) % lp.half_edges.len()]);
            assert_eq!(e0, s1, "the loop chains vertex to vertex at step {w}");
        }
        let mut seen: std::collections::HashMap<_, Vec<bool>> = std::collections::HashMap::new();
        for he in &lp.half_edges {
            seen.entry(he.edge).or_default().push(he.forward);
        }
        let twice: Vec<_> = seen.values().filter(|v| v.len() == 2).collect();
        assert_eq!(twice.len(), 1, "exactly one edge is walked twice: the seam");
        assert_ne!(twice[0][0], twice[0][1], "once in each sense");
    }
}

/// **The grouping joins across a cut rim — through the door production uses.**
///
/// ★★★ A cut circle bounds no whole disk, so the rim-key rule cannot see it, and the
/// unordered node rule cannot either: on a 2-node circle the chord and both complementary
/// arcs fold into one `norm_edge` pair (six users where each piece has two). The `JoinKey`'s
/// ordered arc pairs — "a line is unordered, a circle is ordered", third appearance — give
/// every piece exactly its cap and the band, so the result comes back as **one** component.
///
/// ★ Asserted on `name_result_vertices`' own product (the subdivision must run first for the
/// chord's line key to match), across all three fixtures. red: the band's registration
/// removed → n == 2 and the held grouping is the old `PierceVertexUnnamed`.
#[test]
fn the_grouping_joins_across_a_cut_rim() {
    for (origin, axis) in [
        ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
        ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
        ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let setup = plane_index_setup(&m, plate, boss).unwrap();
        let PlaneSetup {
            planes: faces_tab,
            geom,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
            cyls,
            ..
        } = &setup;
        let jd = Judge::new(geom, *standard, notes);
        let trace_in = crate::combinatorics::trace_input(
            &m,
            [(plate, inc_a), (boss, inc_b)],
            surf_ix,
            faces_tab.len(),
            &jd,
            plane_ix,
            cyls,
            Default::default(),
        );
        let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
            &m,
            BoolKind::Fuse,
            plate,
            boss,
            &jd,
            faces_tab,
            plane_ix,
            cyls,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the arc population traces");
        let faces =
            crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
        let mut faces = faces;
        faces.extend(
            crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                .expect("the lateral faces emit"),
        );
        let seam = crate::arrangement::seam_table(&faces, cyls, &jd)
            .expect("the seam realizes pierce nodes");
        let named =
            crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                .expect("the naming runs");
        let g = named
            .grouping
            .as_ref()
            .expect("the grouping joins across the cut rim");
        assert_eq!(g.n, 1, "one component");
        assert_eq!(g.positives, vec![0], "one material piece, no cavity");
    }
}

/// **A cut-rim boolean builds a complete solid, and the integrals pin it exactly.**
///
/// ★★★ The volume is derived (plate 4·4·2 = 32, boss π·0.25·1, zero overlap — every fixture
/// is a contact), so a wrong segment sign or scale moves an exact number; and because all
/// three fixtures share that one number, the straddling boss's two z = 2 faces additionally
/// split the segment **signs** (the plate top loses its half-disk bite, the digon *is* the
/// bite). The structure is pinned beside it: one solid, no cavities, its outer shell exactly
/// the boolean's minted faces.
#[test]
fn a_cut_rim_boolean_builds_a_complete_solid() {
    for (origin, axis) in [
        ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
        ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
        ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
    ] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(origin),
            Vector3::from_array(axis),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let (faces_from, solids_from) = (m.face_count(), m.solid_count());
        let before = m.live_solids().to_vec();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
            .expect("the cut-rim boolean builds");
        assert_eq!(out.len(), 1, "one fused solid");
        assert_ne!(m.live_solids().to_vec(), before, "the operands retired");
        assert_eq!(m.live_solids().to_vec(), out, "the result lives");
        let pushed: Vec<_> = (solids_from as u32..m.solid_count() as u32)
            .filter_map(|i| m.solid_handle_at(i))
            .map(|h| (h, m.solid(h)))
            .collect();
        let [(sh, solid)] = pushed[..] else {
            panic!("one result solid, got {}", pushed.len());
        };
        assert_eq!(sh, out[0], "the pushed solid is the returned one");
        assert!(solid.cavities.is_empty(), "one material piece, no cavity");
        let shell_faces: std::collections::HashSet<_> =
            m.shell(solid.outer).faces.iter().copied().collect();
        let minted: std::collections::HashSet<_> = (faces_from as u32..m.face_count() as u32)
            .filter_map(|i| m.face_handle_at(i))
            .collect();
        assert_eq!(
            shell_faces, minted,
            "the outer shell is exactly the boolean's minted faces"
        );
        // ★★ **The arc integrals answer exactly** —
        // the volume is derived (plate 4·4·2 = 32, boss π·0.25·1, zero overlap: every
        // fixture is a contact), so a wrong segment sign or scale moves an exact number.
        let props = nacre_props::mass_props(&m, sh).expect("the mixed loops integrate");
        let expect = 32.0 + std::f64::consts::PI / 4.0;
        assert!(
            (props.volume - expect).abs() < 1e-12,
            "volume {} vs derived {expect}",
            props.volume
        );
        // ★ All three fixtures share that one number, so a global scale error fools them
        // together — the straddling boss's two z = 2 faces split the segment **signs**: the
        // plate top loses its half-disk bite (16 − π/8), the overhang digon *is* the other
        // half-disk (π/8).
        if origin == [4.0, 2.0, 2.0] {
            let (mut top, mut digon) = (None, None);
            let mut i = faces_from as u32;
            while let Some(h) = m.face_handle_at(i) {
                i += 1;
                let f = m.face(h);
                let nacre_geom::Surface::Plane(p) = m.surface_cache(f.surface) else {
                    continue;
                };
                if (p.normal().as_array()[2] - 1.0).abs() > 1e-9 {
                    continue;
                }
                match f.outer.half_edges.len() {
                    2 => digon = Some(h),
                    6 => top = Some(h),
                    1 => {} // the boss's closed top cap
                    n => panic!("an unexpected z-normal face with {n} half-edges"),
                }
            }
            let area = |h| nacre_props::face_props(&m, h).expect("planar").area;
            let bite = std::f64::consts::PI / 8.0;
            assert!(
                (area(top.expect("plate top")) - (16.0 - bite)).abs() < 1e-12,
                "the plate top loses its bite"
            );
            assert!(
                (area(digon.expect("overhang digon")) - bite).abs() < 1e-12,
                "the digon is the bite"
            );
        }
    }
}

#[test]
#[ignore = "scratch OBJ dump for eyeballing — run on demand"]
fn measure_dump_straddling_boss_obj() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array([4.0, 2.0, 2.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        1.0,
    );
    m.rebuild_adjacency();
    crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("builds");
    m.rebuild_adjacency();
    let obj = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
        .expect("tessellates")
        .to_obj();
    // Written only on request — the `--ignored` sweep runs this test too, and a test that
    // writes outside the workspace on every sweep is a side effect nobody asked for.
    match std::env::var("OBJ_OUT") {
        Ok(path) => {
            std::fs::write(&path, obj).expect("write");
            println!("wrote {path}");
        }
        Err(_) => println!("set OBJ_OUT=<path> to write the OBJ ({} bytes)", obj.len()),
    }
}

/// The chained bore-then-boss body: `cut` a through-bore at `(12,12)`, then `fuse` a boss
/// whose rim the plate's edge cuts. The chain is what the milestone fence could not ask:
/// nesting must place the bore's rim circle inside a **bitten** top ring — two pierce
/// corners and an arc step — where the chart road's parity had no rational corners to read
/// (`WitnessNotRational`, chaining wall 3) and the mixed parity answers in ℚ(√c).
fn a_bored_plate_with_a_boss(boss_base: [f64; 3]) -> f64 {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 20.0]),
    );
    let bore = m.add_cylinder(
        Point3::from_array([12.0, 12.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        4.0,
        20.0,
    );
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
    let boss = m.add_cylinder(
        Point3::from_array(boss_base),
        Vector3::from_array([0.0, 0.0, 1.0]),
        5.0,
        10.0,
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, bored, boss).expect("the boss fuses");
    assert_eq!(out.len(), 1, "one solid");
    m.rebuild_adjacency();
    assert_eq!(
        m.live_solids().to_vec(),
        out,
        "the operands retired, the result lives"
    );
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
        .expect("the arcs tessellate");
    let mut uses: std::collections::HashMap<(u32, u32), usize> = std::collections::HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let open = uses.values().filter(|&&n| n != 2).count();
    assert_eq!(open, 0, "the mesh is watertight");
    nacre_props::mass_props(&m, out[0]).expect("props").volume
}

/// ★ **A bored plate takes a straddling boss — the first two-boolean chain through the arc
/// population is green.** Volume `40·40·20 − π·4²·20 + π·5²·10 = 32000 − 70π` (the boss
/// stands wholly above the plate top, so fuse adds its full cylinder).
#[test]
fn a_bored_plate_takes_a_straddling_boss() {
    let v = a_bored_plate_with_a_boss([40.0, 20.0, 20.0]);
    let want = 32000.0 - 70.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// The mirror: the boss hangs under the plate's bottom edge — same chain, same volume, the
/// cut circle at the band's hi end instead of lo.
#[test]
fn a_bored_plate_hangs_a_straddling_boss() {
    let v = a_bored_plate_with_a_boss([40.0, 20.0, -10.0]);
    let want = 32000.0 - 70.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★ **The chained contact-cut builds clean** (grouping-arm cell). Cut the same chain
/// instead of fusing: the boss only *touches* the bored plate's top, so the cut removes
/// nothing — and the result is the bored plate, exactly. It used to refuse
/// `VertexNamesAbsentSurface` here: the top ring's rim seam could not merge (the arc pair
/// collided in the merge's node-pair key), the unmerged ring kept pierce vertices whose
/// definitions name the boss's cylinder, and the result has no face on it. The two-pass
/// erase removes the seam — and the vertices with it — so the honest refusal became the
/// honest build. Volume and validate lock that the build is *right*, not merely green.
#[test]
fn a_chained_contact_cut_builds_clean() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 20.0]),
    );
    let bore = m.add_cylinder(
        Point3::from_array([12.0, 12.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        4.0,
        20.0,
    );
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
    let boss = m.add_cylinder(
        Point3::from_array([40.0, 20.0, 20.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        5.0,
        10.0,
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Cut, bored, boss).expect("the contact cut");
    assert_eq!(out.len(), 1, "one body");
    m.rebuild_adjacency();
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = 32000.0 - std::f64::consts::PI * 16.0 * 20.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **Two cylinders with coplanar caps fuse apart.** They stand `5` apart with their caps in
/// the same two planes, which is what once made them look like a seating problem; what was
/// actually hard was classifying **two curved bodies**, and both are now probed from a cap
/// disk's centre — neither has a vertex to probe from.
///
/// ★ Both volumes are the same number derived the same way (`π·0.5²·2`), so the assertion
/// cannot pass by matching one body against the other.
#[test]
fn two_cylinders_with_coplanar_caps_fuse_apart() {
    let mut m = Model::new();
    let a = m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    let b = m.add_cylinder(
        Point3::from_array([5.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, a, b).expect("two curved bodies");
    assert_eq!(out.len(), 2, "they do not touch");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let want = std::f64::consts::PI * 0.25 * 2.0;
    for &s in &out {
        let v = nacre_props::mass_props(&m, s).expect("props").volume;
        assert!((v - want).abs() < 1e-9, "an untouched cylinder: {v}");
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            2,
            "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
    }
}
