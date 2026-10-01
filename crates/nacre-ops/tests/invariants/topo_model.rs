//! The topology layer's own doors and live set, exercised on a box the product builds
//! (`nacre_ops::fixtures`): the write doors' refusals, the cross-model handle guard,
//! `reversed_shell`, seed interning and the live set's order. They live here because a solid
//! is needed and topo's own tests cannot build one through the operations.

use nacre_math::Point3;
use nacre_ops::fixtures::cuboid;
use nacre_topo::{Adjacency, Model, Solid};

/// The unit box, adjacency built.
fn unit_cube() -> Model {
    let mut m = Model::new();
    cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    m.rebuild_adjacency();
    m
}

/// ★★★ **A foreign handle does not read a cache** — the guard [`Store::get`] owns, kept on
/// the index-parallel cache reads too ([`Model::debug_guard`]).
///
/// ★ `Model::surface`, `vertex_point` and `edge_curve` index a `Vec`, so without the guard a
/// handle from another model reaches all three and answers with the wrong cell. The guard is
/// `cfg(debug_assertions)`, so the
/// test is too. Panic output is left unsuppressed on purpose: swapping the panic hook is
/// global state, and the suite runs tests in parallel.
#[test]
#[cfg(debug_assertions)]
fn a_foreign_handle_cannot_read_a_cache() {
    let mut a = Model::new();
    let mut b = Model::new();
    for m in [&mut a, &mut b] {
        cuboid(
            m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0; 3]),
        );
    }
    let surf = b.surface_handle_at(4).expect("b has surfaces");
    let vert = b.vertex_handle_at(3).expect("b has vertices");
    let edge = b.edge_handle_at(3).expect("b has edges");

    let refuses = |what: &str, call: &dyn Fn()| {
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(call));
        assert!(
            r.is_err(),
            "{what}: a handle minted by another model must not read this model's cell"
        );
    };
    refuses("surface_cache", &|| {
        let _ = a.surface_cache(surf);
    });
    refuses("surface", &|| {
        let _ = a.surface(surf);
    });
    refuses("vertex_point", &|| {
        let _ = a.vertex_point(vert);
    });
    refuses("vertex_cache", &|| {
        let _ = a.vertex_cache(vert);
    });
    refuses("edge_curve", &|| {
        let _ = a.edge_curve(edge);
    });
}

/// ★★★★ **The write doors bite** — the negative control for [`Model::push_face`] and
/// [`Model::push_shell`].
///
/// ★ It is measured that *every face the product builds closes its loop*: one temporary
/// assertion, the whole suite, 1324 tests, a single red — and that one was a fixture whose
/// own comment called it a franken-face. That is a **different proposition** from *the
/// assertion refuses a face that does not close*. The first says the population is clean;
/// the second says the door has teeth. Only this test says the second, and without it the
/// doors could assert nothing at all and every green would still be green.
///
/// Each case violates **exactly one** assertion and the panic **message is checked**, because
/// asking only "did it panic" lets a case go green for the wrong reason: a face
/// cloned out of a *throwaway* model makes walking its loop hit
/// `Store`'s cross-model handle guard before ever reaching the door's own assertion. The face
/// comes from the very model it is pushed into, and only the out-of-bounds handles are
/// strangers — those are read with `.index()`, which no guard sees. The assertions are
/// `debug_assert`, so this test is `cfg(debug_assertions)` — the shape the foreign-handle
/// lock above already uses. Panic output is left unsuppressed on purpose: swapping the
/// panic hook is global state, and the suite runs its tests in parallel.
#[test]
#[cfg(debug_assertions)]
fn the_write_doors_refuse_what_their_invariants_forbid() {
    use nacre_topo::{Loop, Shell};
    // The panic message is checked, not just the panic: the claim above is that each case
    // trips *its own* assertion, and only the message can say which one bit.
    fn refuses(what: &str, expect: &str, call: impl FnOnce()) {
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(call));
        let e = r
            .err()
            .unwrap_or_else(|| panic!("{what}: the door let an invalid cell into the arena"));
        let msg = e
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| e.downcast_ref::<&str>().copied())
            .unwrap_or("<non-string panic>");
        assert!(
            msg.contains(expect),
            "{what}: tripped a different assertion — wanted {expect:?}, got {msg:?}"
        );
    }
    let cube = unit_cube;
    let sample = |m: &Model| {
        m.face(m.face_handle_at(0).expect("a cuboid has faces"))
            .clone()
    };

    // A bigger model mints handles this one does not hold. It is the only way to name an
    // out-of-bounds cell: `handle_at` answers `None` past the end, by design.
    let mut big = cube();
    cuboid(
        &mut big,
        Point3::from_array([2.0, 2.0, 2.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let stranger_surface = big
        .surface_handle_at((big.surface_count() - 1) as u32)
        .expect("the bigger model has surfaces");
    let stranger_face = big
        .face_handle_at((big.face_count() - 1) as u32)
        .expect("the bigger model has faces");

    // (1) a face names a surface the arena does not hold
    let mut m = cube();
    let mut f = sample(&m);
    f.surface = stranger_surface;
    refuses(
        "push_face / surface in bounds",
        "a face names a surface the arena does not hold",
        || {
            m.push_face(f);
        },
    );

    // (2) a face's outer loop has no half-edges
    let mut m = cube();
    let mut f = sample(&m);
    f.outer.half_edges.clear();
    refuses(
        "push_face / outer loop non-empty",
        "a face's outer loop has no half-edges",
        || {
            m.push_face(f);
        },
    );

    // (3) a face's inner loop has no half-edges — the outer one is left intact so this
    //     case cannot be carried by (2).
    let mut m = cube();
    let mut f = sample(&m);
    f.inner.push(Loop {
        half_edges: Vec::new(),
    });
    refuses(
        "push_face / inner loop non-empty",
        "a face's inner loop has no half-edges",
        || {
            m.push_face(f);
        },
    );

    // (4) a face loop does not close: flipping one half-edge swaps its end for its start,
    //     which is precisely what walking the loop is there to catch.
    let mut m = cube();
    let mut f = sample(&m);
    f.outer.half_edges[0].forward = !f.outer.half_edges[0].forward;
    refuses(
        "push_face / loop closes",
        "a face loop does not close",
        || {
            m.push_face(f);
        },
    );

    // (5) a shell with no faces bounds nothing
    let mut m = cube();
    refuses(
        "push_shell / non-empty",
        "a shell with no faces bounds nothing",
        || {
            m.push_shell(Shell { faces: Vec::new() });
        },
    );

    // (6) a shell names a face the arena does not hold
    let mut m = cube();
    refuses(
        "push_shell / faces in bounds",
        "a shell names a face the arena does not hold",
        || {
            m.push_shell(Shell {
                faces: vec![stranger_face],
            });
        },
    );
}

#[test]
fn reversed_shell_toggles_orientation_and_reverses_loops() {
    let mut m = unit_cube();
    let outer = m.solid(m.live_solids()[0]).outer;
    let f0 = m.shell(outer).faces[0];
    let orig = m.face(f0).clone();

    let rev_shell = m.reversed_shell(outer);
    // Fresh cells (not reused faces), fresh shell.
    assert_ne!(rev_shell, outer);
    let rf0 = m.shell(rev_shell).faces[0];
    assert_ne!(rf0, f0);
    let rev = m.face(rf0);

    assert_eq!(rev.surface, orig.surface); // surface reused
    assert_eq!(rev.orientation, orig.orientation.flipped());
    let n = orig.outer.half_edges.len();
    assert_eq!(rev.outer.half_edges.len(), n);
    // Reversed winding: he[i] mirrors orig[n-1-i] with the edge reused and
    // the traversal direction flipped.
    for i in 0..n {
        let o = orig.outer.half_edges[n - 1 - i];
        let r = rev.outer.half_edges[i];
        assert_eq!(r.edge, o.edge);
        assert_ne!(r.forward, o.forward);
    }
}

#[test]
fn reversed_shell_is_a_valid_manifold() {
    // Reversing every face's winding preserves the b-rep manifold: each edge
    // is still used by exactly two faces with opposed half-edges. The
    // reversed shell reuses the cube's edges, but the source solid is
    // superseded, so `Adjacency` (reachable-scoped) counts only the reversed
    // faces — no 4-use false positive.
    let mut m = unit_cube();
    let outer = m.solid(m.live_solids()[0]).outer;
    let rev = m.reversed_shell(outer);
    let s = m.push_solid(Solid {
        outer: rev,
        cavities: vec![],
    });
    m.restore_live(vec![s]); // supersede the original cube
    let adj = Adjacency::rebuild(&m);
    assert_eq!(adj.edge_uses.len(), 12);
    for uses in adj.edge_uses.values() {
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0].1, uses[1].1); // opposite forward
    }
}

/// ★★ The seeds are the interning survivors — an origin box's floor, front and left faces
/// carry the seed handles, so "the world plane" and "that face's plane" are one surface,
/// stated once (the floor's datum interns onto the `z = 0` seed by name).
#[test]
fn an_origin_cuboids_axis_faces_intern_onto_the_seeds() {
    let mut m = Model::new();
    let s = cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let face_surfaces: Vec<_> = m
        .shell(m.solid(s).outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .collect();
    for axis in [
        nacre_exact::Axis::Z,
        nacre_exact::Axis::X,
        nacre_exact::Axis::Y,
    ] {
        assert!(
            face_surfaces.contains(&m.world_plane(axis)),
            "the {axis:?}-normal face at 0 must be the seed itself"
        );
    }
    // And nothing pointless was minted: 3 seeds + the 3 off-origin faces.
    assert_eq!(m.surface_count(), 6);
}

/// ★★★★ **`supersede_live` keeps the survivors in their original order — and nothing else in
/// this tree would notice if it did not.**
///
/// `nacre_step::to_step` exports the live set *in order*, the census reads the **arena** (it
/// never sees live order), and there is no golden STEP text anywhere. So a permutation here
/// would travel all the way out to the exported file unseen. The order is therefore locked at
/// the door itself, and again on the export side (`nacre-step`'s
/// `superseding_a_solid_leaves_the_export_order_alone`).
///
/// ⚠ The oracle is the **construction order**, not a second call to the same machinery —
/// comparing this against a hand-written `retain` would be one implementation checking itself.
#[test]
fn supersede_live_preserves_the_order_of_the_survivors() {
    let mut m = Model::new();
    let a = cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let b = cuboid(
        &mut m,
        Point3::from_array([10.0; 3]),
        Point3::from_array([11.0; 3]),
    );
    let c = cuboid(
        &mut m,
        Point3::from_array([20.0; 3]),
        Point3::from_array([21.0; 3]),
    );
    assert_eq!(
        m.live_solids(),
        [a, b, c].as_slice(),
        "the fixture's premise"
    );

    m.supersede_live(&[b]);
    assert_eq!(
        m.live_solids(),
        [a, c].as_slice(),
        "the middle solid goes and the order stays"
    );

    // It drops sets, not only singletons — that is what the 19 `retain` call sites ask for.
    m.supersede_live(&[a, c]);
    assert!(m.live_solids().is_empty());
}

/// `restore_live` is the rollback half: what a rejected operation puts back.
#[test]
fn restore_live_puts_the_snapshot_back() {
    let mut m = Model::new();
    let a = cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let b = cuboid(
        &mut m,
        Point3::from_array([10.0; 3]),
        Point3::from_array([11.0; 3]),
    );
    let snapshot = m.live_solids().to_vec();

    m.supersede_live(&[a]);
    assert_eq!(m.live_solids(), [b].as_slice());

    m.restore_live(snapshot);
    assert_eq!(m.live_solids(), [a, b].as_slice(), "order comes back too");
}
