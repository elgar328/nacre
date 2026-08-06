//! **The replay contract, measured.**
//!
//! `replay`'s doc (and `docs/design.md` §6) promises «the same log reproduces the same model
//! **down to handle indices**». Six of `Operation`'s seven variants carry a handle
//! (`PadOnFace`·`PocketOnFace`·`Boolean`·`Transform`·`Mirror`·`Copy`) — and for those the
//! promise is **false today**: `replay` starts from a fresh `Model::new()`, so the log's
//! handles carry the *original* model's `StoreId` and `Store::get`'s debug guard fires before
//! the index is ever used.
//!
//! Nothing noticed because nothing walks that path: all 46 `replay` call sites in the
//! workspace pass value-only `Extrude` logs. So today's replay determinism rests entirely on
//! `Extrude` stating its plane by value — and the migration is about to take that away.
//!
//! This file is where the contract is *measured* rather than asserted by documentation.

use nacre_math::Point2;
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply, replay};
use nacre_scalar::{Isometry, Rat};
use nacre_topo::{Model, VertexDef};

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
    .expect("a square is a fair profile")
}

fn extrude_op(a: f64, b: f64, dist: f64) -> Operation {
    Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: square(a, b),
        dist,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// The comparator
// ─────────────────────────────────────────────────────────────────────────────

/// One named, indexed observation of a model's arena.
///
/// A signature is a **vector of these**, not one string or hash: when two models differ, the
/// interesting question is *where first*, and a comparator that answers only "different" hands
/// that work back to a human. (`docs/dev-log.md`'s S9 ε-gate prints which line and which field
/// moved by how much, for the same reason.)
#[derive(Clone, PartialEq, Debug)]
struct SigItem {
    what: &'static str,
    at: usize,
    value: String,
}

/// **What "down to handle indices" actually means**, made observable.
///
/// Store lengths, then every cell in index order stated purely by *index* — a vertex's
/// defining surfaces, an edge's carriers and endpoints, a face's surface and half-edge
/// sequence, a shell's faces, a solid's shells, and the live set's order. Coordinates and
/// measured tolerances ride along as bits so a model that is topologically identical but
/// geometrically moved is still caught.
///
/// Surfaces are covered without iterating their (sealed) store: `surface_count` pins the
/// population, and every reference above pins identity — an orphan plane changes the count, a
/// mis-wired one changes a reference.
fn arena_sig(m: &Model) -> Vec<SigItem> {
    let mut out = Vec::new();
    let mut push =
        |what: &'static str, at: usize, value: String| out.push(SigItem { what, at, value });

    push("len.vertices", 0, m.vertices.len().to_string());
    push("len.edges", 0, m.edges.len().to_string());
    push("len.faces", 0, m.faces.len().to_string());
    push("len.shells", 0, m.shells.len().to_string());
    push("len.solids", 0, m.solids.len().to_string());
    push("len.surfaces", 0, m.surface_count().to_string());

    for (vh, v) in m.vertices.iter() {
        let i = vh.index() as usize;
        let p = m.vertex_point(vh).as_array();
        push(
            "vertex.coord",
            i,
            format!(
                "{:x},{:x},{:x}",
                p[0].to_bits(),
                p[1].to_bits(),
                p[2].to_bits()
            ),
        );
        push(
            "vertex.tol",
            i,
            match m.vertex_tol(vh) {
                Some(t) => format!("{:x}", t.to_bits()),
                None => "-".to_owned(),
            },
        );
        push(
            "vertex.def",
            i,
            match v.def {
                VertexDef::ThreePlane(s) => {
                    format!("3p[{},{},{}]", s[0].index(), s[1].index(), s[2].index())
                }
                VertexDef::OnSeam(s) => format!("seam[{},{}]", s[0].index(), s[1].index()),
            },
        );
    }
    for (eh, e) in m.edges.iter() {
        let i = eh.index() as usize;
        push(
            "edge.carriers",
            i,
            format!("{},{}", e.surfaces[0].index(), e.surfaces[1].index()),
        );
        push(
            "edge.vertices",
            i,
            format!("{},{}", e.vertices[0].index(), e.vertices[1].index()),
        );
    }
    for (fh, f) in m.faces.iter() {
        let i = fh.index() as usize;
        push("face.surface", i, f.surface.index().to_string());
        push("face.orientation", i, format!("{:?}", f.orientation));
        let loops: Vec<String> = std::iter::once(&f.outer)
            .chain(f.inner.iter())
            .map(|lp| {
                lp.half_edges
                    .iter()
                    .map(|he| format!("{}{}", he.edge.index(), if he.forward { "+" } else { "-" }))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect();
        push("face.loops", i, loops.join(" | "));
    }
    for (sh, shell) in m.shells.iter() {
        push(
            "shell.faces",
            sh.index() as usize,
            shell
                .faces
                .iter()
                .map(|h| h.index().to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
    }
    for (sh, solid) in m.solids.iter() {
        let i = sh.index() as usize;
        push("solid.outer", i, solid.outer.index().to_string());
        push(
            "solid.cavities",
            i,
            solid
                .cavities
                .iter()
                .map(|h| h.index().to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
    }
    push(
        "live_solids",
        0,
        m.live_solids
            .iter()
            .map(|h| h.index().to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    for &s in &m.live_solids {
        let v = nacre_props::mass_props(m, s)
            .map(|p| p.volume)
            .unwrap_or(f64::NAN);
        push("volume", s.index() as usize, format!("{:x}", v.to_bits()));
    }
    out
}

/// Assert two arenas agree, naming **the first place they do not**.
fn assert_same_arena(a: &Model, b: &Model, what: &str) {
    let (sa, sb) = (arena_sig(a), arena_sig(b));
    for (x, y) in sa.iter().zip(&sb) {
        assert_eq!(
            x, y,
            "{what}: arenas first differ at {}[{}] — {:?} vs {:?}",
            x.what, x.at, x.value, y.value
        );
    }
    assert_eq!(
        sa.len(),
        sb.len(),
        "{what}: arenas agree item-for-item but one has more ({} vs {})",
        sa.len(),
        sb.len()
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Fixtures — logs that carry a handle
// ─────────────────────────────────────────────────────────────────────────────

/// `[Extrude, PadOnFace]` — the log and a model built by applying it.
fn pad_log() -> (Vec<Operation>, Model) {
    let mut m = Model::new();
    let ex = extrude_op(0.0, 2.0, 1.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &ex).expect("extrude") else {
        unreachable!()
    };
    let pad = Operation::PadOnFace {
        face: faces[1], // top cap
        profile: square(0.5, 1.5),
        dist: 0.5,
    };
    apply(&mut m, &pad).expect("pad on the live model");
    m.rebuild_adjacency();
    (vec![ex, pad], m)
}

/// `[Extrude, Extrude, Boolean]`.
fn boolean_log() -> (Vec<Operation>, Model) {
    let mut m = Model::new();
    let (e1, e2) = (extrude_op(0.0, 2.0, 1.0), extrude_op(1.0, 3.0, 1.0));
    let OpOutput::Extrude { solid: a, .. } = apply(&mut m, &e1).expect("a") else {
        unreachable!()
    };
    let OpOutput::Extrude { solid: b, .. } = apply(&mut m, &e2).expect("b") else {
        unreachable!()
    };
    let bo = Operation::Boolean {
        kind: BoolKind::Fuse,
        a,
        b,
    };
    apply(&mut m, &bo).expect("fuse on the live model");
    m.rebuild_adjacency();
    (vec![e1, e2, bo], m)
}

/// `[Extrude, Transform]`.
fn transform_log() -> (Vec<Operation>, Model) {
    let mut m = Model::new();
    let ex = extrude_op(0.0, 2.0, 1.0);
    let OpOutput::Extrude { solid, .. } = apply(&mut m, &ex).expect("extrude") else {
        unreachable!()
    };
    let xf = Operation::Transform {
        solid,
        isometry: Isometry::translation([Rat::from_int(5), Rat::from_int(0), Rat::from_int(0)]),
    };
    apply(&mut m, &xf).expect("translate on the live model");
    m.rebuild_adjacency();
    (vec![ex, xf], m)
}

// ─────────────────────────────────────────────────────────────────────────────
// The witnesses — what today actually does
// ─────────────────────────────────────────────────────────────────────────────
//
// These three pin the defect as a *fact*, not as a description. Each dies in
// `nacre-store`'s `Store::get` guard (`"Handle was minted by a different Store"`), reached
// through the op's first dereference: pad at `ops.rs`'s `shells.get(sh)`, boolean at
// `planes.rs`'s `solids.get(solid)`, transform at `transform.rs`'s `solids.get(solid)`.
// Note what gets *past* first: `live_solids.contains(&h)` and `faces.contains(&face)` compare
// `Handle`s by index only, so a foreign handle answers "yes, I am live" — the guard is not the
// door, it is one step inside it.
//
// The repair commit deletes these three and asserts success in their place.

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "different Store")]
fn a_pad_log_cannot_replay_before_the_repair() {
    let (log, _) = pad_log();
    let _ = replay(&log);
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "different Store")]
fn a_boolean_log_cannot_replay_before_the_repair() {
    let (log, _) = boolean_log();
    let _ = replay(&log);
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "different Store")]
fn a_transform_log_cannot_replay_before_the_repair() {
    let (log, _) = transform_log();
    let _ = replay(&log);
}

/// ★ **The profiles disagree**, and that is the sharpest statement of the defect.
///
/// `StoreId` is `#[cfg(debug_assertions)]`, so in release the guard does not exist: the same
/// log replays fine, because the index was right all along. Debug panics; release is quietly
/// correct. The repair therefore does not *add* a capability — it makes debug agree with what
/// release already does, by construction rather than by luck.
///
/// **Run this by hand** (`cargo test -p nacre-ops --release --test replay`) and record the
/// numbers — it is a one-shot measurement, not a standing gate: the pre-commit hook runs debug
/// only, the same way census runs `--release` separately. Once the repair lands, both profiles
/// give the same answer and this twin is replaced by a profile-independent assertion.
#[cfg(not(debug_assertions))]
#[test]
fn in_release_a_handle_carrying_log_already_replays() {
    for (name, (log, scratch)) in [
        ("pad", pad_log()),
        ("boolean", boolean_log()),
        ("transform", transform_log()),
    ] {
        let replayed = replay(&log).unwrap_or_else(|e| panic!("{name}: replay failed: {e:?}"));
        assert_same_arena(&replayed, &scratch, name);
        println!("stat release_replay {name} ok");
    }
}

/// The comparator's own check: a value-only log (the population that *does* replay today)
/// reproduces its scratch-built twin item for item, and replaying twice is stable.
#[test]
fn the_comparator_agrees_on_a_value_only_log() {
    let log = vec![extrude_op(0.0, 2.0, 1.0), extrude_op(3.0, 4.0, 0.5)];
    let mut scratch = Model::new();
    for op in &log {
        apply(&mut scratch, op).expect("apply");
    }
    scratch.rebuild_adjacency();

    let a = replay(&log).expect("value-only logs replay today");
    let b = replay(&log).expect("value-only logs replay today");
    assert_same_arena(&a, &b, "replay twice");
    assert_same_arena(&a, &scratch, "replay vs scratch");
}

/// And it is not vacuous: two different models must differ, at a named place.
#[test]
fn the_comparator_notices_a_difference() {
    let a = replay(&[extrude_op(0.0, 2.0, 1.0)]).unwrap();
    let b = replay(&[extrude_op(0.0, 2.0, 2.0)]).unwrap();
    let (sa, sb) = (arena_sig(&a), arena_sig(&b));
    let first = sa
        .iter()
        .zip(&sb)
        .find(|(x, y)| x != y)
        .map(|(x, _)| (x.what, x.at));
    assert_eq!(
        first,
        Some(("vertex.coord", 4)),
        "a taller prism must first differ at a vertex coordinate"
    );
}
