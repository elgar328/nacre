//! **The replay contract, measured.**
//!
//! `replay`'s doc (and `docs/design.md` §6) promises «the same log reproduces the same model
//! **down to handle indices**». Six of `Operation`'s seven variants carry a handle
//! (`PadOnFace`·`PocketOnFace`·`Boolean`·`Transform`·`Mirror`·`Copy`) — and for those the
//! promise is **false today**: `replay` starts from a fresh `Model::new()`, so the log's
//! handles carry the *original* model's `StoreId` and `Store::get`'s debug guard fires before
//! the index is ever used.
//!
//! Nothing noticed because nothing walked that path: every `replay` call site in the workspace
//! passed value-only `Extrude` logs, so replay determinism rested entirely on `Extrude` stating
//! its plane by value.
//!
//! ★ **That sentence is now history.** S5(i)-b gave `Extrude` a `SketchFrame`, so **all seven
//! variants carry a handle** and there is no value-only operation left. The premise this file was
//! written to expose has been consumed; what it measures — that a log's indices are re-anchored
//! onto the arena being built — is now load-bearing for every log there is.
//!
//! This file is where the contract is *measured* rather than asserted by documentation.

use nacre_math::Point2;
use nacre_ops::SketchFrame;
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply, replay};
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

/// A world-XY extrude **for `m`** — an extrude names its plane now, and a handle only means
/// something in its own arena until `replay` re-anchors it.
fn extrude_op(m: &Model, a: f64, b: f64, dist: f64) -> Operation {
    Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
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
            "vertex.bound",
            i,
            match m.vertex_cache(vh).bound() {
                Some(b) => format!("{b:?}"),
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
                VertexDef::Pierce {
                    planes,
                    cylinder,
                    root,
                } => format!(
                    "pierce[{},{};{};{root:?}]",
                    planes[0].index(),
                    planes[1].index(),
                    cylinder.index()
                ),
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
    let ex = extrude_op(&m, 0.0, 2.0, 1.0);
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
    let (e1, e2) = (extrude_op(&m, 0.0, 2.0, 1.0), extrude_op(&m, 1.0, 3.0, 1.0));
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
    let ex = extrude_op(&m, 0.0, 2.0, 1.0);
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
// The contract, now kept
// ─────────────────────────────────────────────────────────────────────────────
//
// These three used to be `#[should_panic(expected = "different Store")]` witnesses: each log
// died in `Store::get`'s guard, reached through the op's first dereference (pad at `shells.get`,
// boolean and transform at `solids.get`). What got *past* first is worth remembering —
// `live_solids.contains(&h)` and `faces.contains(&f)` compare handles by index only, so a
// foreign handle answered "yes, I am live"; the guard was never the door, it was one step
// inside it.
//
// `replay` now re-anchors the log's indices onto the model it is building, so the same three
// logs reproduce their scratch-built twins — arena for arena, in both profiles.

#[test]
fn a_pad_log_replays() {
    let (log, scratch) = pad_log();
    let replayed = replay(&log).expect("a pad log replays");
    assert_same_arena(&replayed, &scratch, "pad");
}

#[test]
fn a_boolean_log_replays() {
    let (log, scratch) = boolean_log();
    let replayed = replay(&log).expect("a boolean log replays");
    assert_same_arena(&replayed, &scratch, "boolean");
}

#[test]
fn a_transform_log_replays() {
    let (log, scratch) = transform_log();
    let replayed = replay(&log).expect("a transform log replays");
    assert_same_arena(&replayed, &scratch, "transform");
}

/// ★ **The profiles agree now** — and that is the whole shape of the repair.
///
/// `StoreId` is `#[cfg(debug_assertions)]`, so before the repair release was quietly correct
/// (the index had been right all along) while debug panicked. Re-anchoring makes debug reach
/// the same answer *by construction* rather than by the guard happening to be compiled out, so
/// this assertion is profile-independent and replaces the release-only twin that measured the
/// divergence.
#[test]
fn both_profiles_replay_a_handle_carrying_log_the_same_way() {
    for (name, (log, scratch)) in [
        ("pad", pad_log()),
        ("boolean", boolean_log()),
        ("transform", transform_log()),
    ] {
        let replayed = replay(&log).unwrap_or_else(|e| panic!("{name}: replay failed: {e:?}"));
        assert_same_arena(&replayed, &scratch, name);
    }
}

/// A log naming a cell this model does not have is a **named reject**, not a panic — and the
/// same `Err` in both profiles. (Built with a handle taken from a *larger* model, which is how
/// a truncated or spliced log looks from the inside.)
#[test]
fn a_log_naming_a_cell_that_does_not_exist_is_rejected_by_name() {
    // A two-solid model, so its solid handle 1 exists…
    let mut big = Model::new();
    let e1 = extrude_op(&big, 0.0, 2.0, 1.0);
    let e2 = extrude_op(&big, 5.0, 6.0, 1.0);
    apply(&mut big, &e1).expect("first");
    let OpOutput::Extrude { solid: second, .. } = apply(&mut big, &e2).expect("second") else {
        unreachable!()
    };
    // …but a log with only the first extrude leaves the replayed model one solid short.
    let log = vec![
        e1,
        Operation::Transform {
            solid: second,
            isometry: Isometry::translation([Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)]),
        },
    ];
    match replay(&log) {
        Err(nacre_ops::OpError::LogHandleOutOfRange { cell, index }) => {
            assert_eq!(cell, nacre_ops::LogCell::Solid);
            assert_eq!(index, second.index());
        }
        other => panic!(
            "expected a named reject, got {:?}",
            other.map(|_| "a model")
        ),
    }
}

/// The comparator's own check: a value-only log (the population that *does* replay today)
/// reproduces its scratch-built twin item for item, and replaying twice is stable.
#[test]
fn the_comparator_agrees_on_a_value_only_log() {
    let mut scratch = Model::new();
    let log = vec![
        extrude_op(&scratch, 0.0, 2.0, 1.0),
        extrude_op(&scratch, 3.0, 4.0, 0.5),
    ];
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
    let a = replay(&[extrude_op(&Model::new(), 0.0, 2.0, 1.0)]).unwrap();
    let b = replay(&[extrude_op(&Model::new(), 0.0, 2.0, 2.0)]).unwrap();
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

// ─────────────────────────────────────────────────────────────────────────────
// The property that was promised — `docs/design.md` §7 (2)
// ─────────────────────────────────────────────────────────────────────────────

use nacre_scalar::{Angle, Axis, Rotation};
use nacre_store::Handle;
use nacre_topo::{Face, Solid};
use proptest::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([x0, y0]),
        Point2::from_array([x1, y0]),
        Point2::from_array([x1, y1]),
        Point2::from_array([x0, y1]),
    ])
    .expect("an axis-aligned rectangle is a fair profile")
}

/// A **choice**, not an operation.
///
/// A generator cannot name a handle: the cell an operation refers to does not exist until the
/// operations before it have run. So the generator emits ordinals — "the 2nd live solid, its
/// 3rd face" — and [`run_recipe`] resolves them against the model as it grows. This is the
/// same shape as a log: an index that only means something relative to the arena it is read in.
#[derive(Clone, Debug)]
enum Step {
    Extrude {
        x: i8,
        y: i8,
        w: u8,
        dist: u8,
    },
    Pad {
        solid: usize,
        face: usize,
        inset: u8,
        dist: u8,
    },
    Pocket {
        solid: usize,
        face: usize,
        inset: u8,
        dist: u8,
    },
    Boolean {
        kind: u8,
        a: usize,
        b: usize,
    },
    Translate {
        solid: usize,
        d: [i8; 3],
    },
    Rotate {
        solid: usize,
        axis: u8,
        deg: u8,
    },
    Mirror {
        solid: usize,
        axis: u8,
        offset: i8,
    },
    Copy {
        solid: usize,
    },
    /// ★★★★ **A datum that names vertices — the step this generator most needs.**
    ///
    /// Re-anchoring is a bounds check, and a bounds check cannot see an in-range index that names
    /// the *wrong* cell. Every earlier handle-carrying step names a face, a solid or a surface —
    /// few, and little disturbed by a reject. Vertices are neither: a late reject leaves 63–84
    /// cells behind (measured, R), and most of them are vertices, so a log recorded across one
    /// and replayed would re-anchor a datum onto three different corners with no error at all.
    /// Only a generated session that actually places this step after a reject can say it does
    /// not happen.
    DatumThroughVertices {
        solid: usize,
        a: usize,
        b: usize,
        c: usize,
    },
}

fn step_strategy() -> impl Strategy<Value = Step> {
    prop_oneof![
        // ★ Weighted like the other handle-carrying steps, and deliberately *not* filtered to
        // triples that succeed: its rejects are part of what the session must survive.
        3 => (0usize..4, 0usize..12, 0usize..12, 0usize..12)
            .prop_map(|(solid, a, b, c)| Step::DatumThroughVertices { solid, a, b, c }),
        4 => (-3i8..3, -3i8..3, 1u8..4, 1u8..4)
            .prop_map(|(x, y, w, dist)| Step::Extrude { x, y, w, dist }),
        3 => (0usize..4, 0usize..8, 0u8..3, 1u8..3)
            .prop_map(|(solid, face, inset, dist)| Step::Pad { solid, face, inset, dist }),
        3 => (0usize..4, 0usize..8, 0u8..3, 1u8..3)
            .prop_map(|(solid, face, inset, dist)| Step::Pocket { solid, face, inset, dist }),
        3 => (0u8..3, 0usize..4, 0usize..4)
            .prop_map(|(kind, a, b)| Step::Boolean { kind, a, b }),
        2 => (0usize..4, prop::array::uniform3(-3i8..3))
            .prop_map(|(solid, d)| Step::Translate { solid, d }),
        2 => (0usize..4, 0u8..3, 1u8..5)
            .prop_map(|(solid, axis, deg)| Step::Rotate { solid, axis, deg }),
        2 => (0usize..4, 0u8..3, -2i8..2)
            .prop_map(|(solid, axis, offset)| Step::Mirror { solid, axis, offset }),
        2 => (0usize..4).prop_map(|solid| Step::Copy { solid }),
    ]
}

fn axis_of(k: u8) -> Axis {
    match k % 3 {
        0 => Axis::X,
        1 => Axis::Y,
        _ => Axis::Z,
    }
}

/// Every face of a solid's outer shell, in shell order.
fn faces_of(m: &Model, s: Handle<Solid>) -> Vec<Handle<Face>> {
    m.shells.get(m.solids.get(s).outer).faces.clone()
}

/// Turn a [`Step`] into a concrete `Operation` against **this** model, or `None` if the model
/// cannot host it at all (no live solid to name). A step that resolves but will be *rejected*
/// still resolves — the reject is the operation's answer to give, not the generator's.
fn concretize(m: &Model, step: &Step) -> Option<Operation> {
    let live: Vec<Handle<Solid>> = m.live_solids.to_vec();
    let pick = |i: usize| live.get(i % live.len().max(1)).copied();
    Some(match *step {
        Step::DatumThroughVertices { solid, a, b, c } => {
            let s = pick(solid)?;
            let mut vs = Vec::new();
            let sol = m.solids.get(s);
            for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
                for &fh in &m.shells.get(sh).faces {
                    let f = m.faces.get(fh);
                    for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                        for &he in &lp.half_edges {
                            let vh = m.he_start(he);
                            if !vs.contains(&vh) {
                                vs.push(vh);
                            }
                        }
                    }
                }
            }
            if vs.len() < 3 {
                return None;
            }
            let g = |i: usize| vs[i % vs.len()];
            Operation::DatumPlane {
                def: nacre_ops::DatumDef::ThroughVertices([g(a), g(b), g(c)]),
            }
        }
        Step::Extrude { x, y, w, dist } => {
            let (x, y, w) = (f64::from(x), f64::from(y), f64::from(w));
            Operation::Extrude {
                frame: SketchFrame::world(m, Axis::Z),
                profile: rect(x, y, x + w, y + w),
                dist: f64::from(dist),
            }
        }
        Step::Pad {
            solid,
            face,
            inset,
            dist,
        }
        | Step::Pocket {
            solid,
            face,
            inset,
            dist,
        } => {
            let s = pick(solid)?;
            let faces = faces_of(m, s);
            let f = *faces.get(face % faces.len().max(1))?;
            let i = f64::from(inset) * 0.25;
            let profile = rect(i, i, i + 1.0, i + 1.0);
            let dist = f64::from(dist);
            if matches!(step, Step::Pad { .. }) {
                Operation::PadOnFace {
                    face: f,
                    profile,
                    dist,
                }
            } else {
                Operation::PocketOnFace {
                    face: f,
                    profile,
                    dist,
                }
            }
        }
        Step::Boolean { kind, a, b } => {
            let (ha, hb) = (pick(a)?, pick(b)?);
            if ha == hb {
                return None; // a boolean of a solid with itself is not a case, it is a typo
            }
            Operation::Boolean {
                kind: match kind % 3 {
                    0 => BoolKind::Fuse,
                    1 => BoolKind::Cut,
                    _ => BoolKind::Common,
                },
                a: ha,
                b: hb,
            }
        }
        Step::Translate { solid, d } => Operation::Transform {
            solid: pick(solid)?,
            isometry: Isometry::translation([
                Rat::from_int(i128::from(d[0])),
                Rat::from_int(i128::from(d[1])),
                Rat::from_int(i128::from(d[2])),
            ]),
        },
        Step::Rotate { solid, axis, deg } => Operation::Transform {
            solid: pick(solid)?,
            isometry: Isometry::rotation(Rotation {
                axis: axis_of(axis),
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(i128::from(deg) * 15))?,
            }),
        },
        Step::Mirror {
            solid,
            axis,
            offset,
        } => Operation::Mirror {
            solid: pick(solid)?,
            axis: axis_of(axis),
            offset: Rat::from_int(i128::from(offset)),
        },
        Step::Copy { solid } => Operation::Copy {
            solid: pick(solid)?,
        },
    })
}

/// What a generated session did, so the population can be reported instead of guessed.
#[derive(Default, Debug)]
struct Stats {
    attempted: usize,
    accepted: usize,
    rejected: usize,
    unresolvable: usize,
    handle_carrying: usize,
    /// ★ Vertex-naming datums accepted **after this session already rejected something** — the
    /// only shape in which "in range but the wrong vertex" could show itself, since a reject is
    /// what shifts the indices. Counted so the coverage cannot vanish silently.
    datum_after_reject: usize,
}

/// Interpret a recipe into `(log, model)` — **the discipline, executed.**
///
/// A rejected operation is not in the log, but it may already have pushed cells into the arena
/// (see [`a_late_reject_is_not_index_neutral`]). So after any reject the session throws its
/// model away and rebuilds it from the log. That is exactly the rule `docs/design.md` states
/// for a session that keeps recording after a reject, and it is what makes property (3)
/// (`replay(log)` reproduces the session) true rather than lucky.
fn run_recipe(steps: &[Step]) -> (Vec<Operation>, Model, Stats) {
    let mut model = Model::new();
    let mut log: Vec<Operation> = Vec::new();
    let mut st = Stats::default();

    // Seed: every session starts with something to name.
    let seed = extrude_op(&model, 0.0, 2.0, 1.0);
    apply(&mut model, &seed).expect("the seed extrude is unconditionally legal");
    log.push(seed);

    for step in steps {
        let Some(op) = concretize(&model, step) else {
            st.unresolvable += 1;
            continue;
        };
        st.attempted += 1;
        match apply(&mut model, &op) {
            Ok(_) => {
                // ★ Every variant names a cell now (S5(i)-b gave `Extrude` a frame), so this
                // counts accepted operations. It used to exclude `Extrude` because that was the
                // one variant carrying nothing.
                st.handle_carrying += 1;
                st.accepted += 1;
                if st.rejected > 0
                    && matches!(
                        op,
                        Operation::DatumPlane {
                            def: nacre_ops::DatumDef::ThroughVertices(_)
                        }
                    )
                {
                    st.datum_after_reject += 1;
                }
                log.push(op);
            }
            Err(_) => {
                st.rejected += 1;
                // Re-sync: the arena may carry residue from the declined operation.
                model = replay(&log).expect("a log of accepted operations replays");
            }
        }
    }

    // ★ Property (4) used to need a fallback here — if a recipe landed no handle-carrying
    // operation, the session closed with a `Copy` so the case proved something. Since S5(i)-b
    // there is no operation that carries no handle: `Extrude` names its plane, so the seed alone
    // satisfies it. The counter below stays because the *report* is still worth having; the
    // fallback is gone because it became unreachable.
    // The seed extrude names the world plane, so `handle_carrying >= 1` holds before any
    // generated step runs — property (4) is satisfied by construction rather than by a fallback.
    st.handle_carrying += 1;

    model.rebuild_adjacency();
    (log, model, st)
}

static CASES: AtomicUsize = AtomicUsize::new(0);
static ATTEMPTED: AtomicUsize = AtomicUsize::new(0);
static REJECTED: AtomicUsize = AtomicUsize::new(0);
static HANDLED: AtomicUsize = AtomicUsize::new(0);
static DATUM_AFTER_REJECT: AtomicUsize = AtomicUsize::new(0);
const CASE_COUNT: u32 = 32;

proptest! {
    #![proptest_config(ProptestConfig { cases: CASE_COUNT, ..ProptestConfig::default() })]

    /// ★ **The gate of this stage.** Four properties on one generated session:
    ///
    /// 1. the session model is a valid b-rep,
    /// 2. `replay` is **idempotent** — replaying twice lands on the same arena
    ///    (`docs/design.md` §7 (2), promised and until now unimplemented),
    /// 3. ★ `replay(log)` reproduces the session **down to handle indices** — the property the
    ///    contract advertised and no call site had ever exercised, because every log in the
    ///    workspace was value-only (which S5(i)-b then ended — every variant names a cell now),
    /// 4. every case carries at least one handle-bearing operation, guaranteed by construction
    ///    in [`run_recipe`], so (3) is never vacuously true.
    ///
    /// The population deliberately includes rotations and mirrors, i.e. widths that CIP has to
    /// judge: property (5) is that such a step is either a valid result or a **named reject** —
    /// never a panic, which the harness enforces for free.
    #[test]
    fn a_generated_session_replays_to_itself(steps in prop::collection::vec(step_strategy(), 1..5)) {
        let (log, session, st) = run_recipe(&steps);

        prop_assert!(st.handle_carrying >= 1, "a case with no handle-carrying op proves nothing");

        let violations = nacre_validate::validate(&session);
        prop_assert!(violations.is_empty(), "session model is invalid: {violations:?}");

        let once = replay(&log).expect("the accepted log replays");
        assert_same_arena(&once, &session, "replay(log) vs session");

        let twice = replay(&log).expect("replaying again is still a replay");
        assert_same_arena(&twice, &once, "replay idempotence");

        ATTEMPTED.fetch_add(st.attempted, Ordering::Relaxed);
        REJECTED.fetch_add(st.rejected, Ordering::Relaxed);
        HANDLED.fetch_add(st.handle_carrying, Ordering::Relaxed);
        DATUM_AFTER_REJECT.fetch_add(st.datum_after_reject, Ordering::Relaxed);
        if CASES.fetch_add(1, Ordering::Relaxed) + 1 == CASE_COUNT as usize {
            // A population that rejects everything measures nothing; print it, do not assume it.
            println!(
                "stat generated_population cases={CASE_COUNT} attempted={} rejected={} handle_carrying={} datum_after_reject={}",
                ATTEMPTED.load(Ordering::Relaxed),
                REJECTED.load(Ordering::Relaxed),
                HANDLED.load(Ordering::Relaxed),
                DATUM_AFTER_REJECT.load(Ordering::Relaxed),
            );
            // ★ Reported, not asserted: a random generator cannot be *made* to produce a shape,
            // and a suite that fails when it does not is flaky rather than strict. The shape this
            // number watches for is guaranteed deterministically instead, by
            // [`a_datum_naming_vertices_survives_a_session_that_rejected`].
        }
    }
}

/// Property (6), in plain text: **all six** handle-carrying variants in one log.
///
/// The generator reaches them by chance; this reaches them by name, so a variant that silently
/// stops being generated cannot take the coverage with it.
#[test]
fn a_log_using_every_handle_carrying_variant_replays() {
    let mut m = Model::new();
    let mut log = Vec::new();
    let mut run = |m: &mut Model, op: Operation| {
        let out = apply(m, &op).expect("each step of this fixture is legal");
        log.push(op);
        out
    };

    let op = extrude_op(&m, 0.0, 2.0, 2.0);
    let OpOutput::Extrude { faces, .. } = run(&mut m, op) else {
        unreachable!()
    };
    let top = faces[1];
    // Every edit but `Copy` supersedes its input, so the fixture threads the *new* handle on.
    let OpOutput::PadOnFace { solid: a, .. } = run(
        &mut m,
        Operation::PadOnFace {
            face: top,
            profile: rect(0.25, 0.25, 1.75, 1.75),
            dist: 1.0,
        },
    ) else {
        unreachable!()
    };
    let OpOutput::Copy { solid: b } = run(&mut m, Operation::Copy { solid: a }) else {
        unreachable!()
    };
    let OpOutput::Transform { solid: b } = run(
        &mut m,
        Operation::Transform {
            solid: b,
            isometry: Isometry::translation([Rat::from_int(5), Rat::from_int(0), Rat::from_int(0)]),
        },
    ) else {
        unreachable!()
    };
    let OpOutput::Mirror { solid: b } = run(
        &mut m,
        Operation::Mirror {
            solid: b,
            axis: Axis::X,
            offset: Rat::from_int(9),
        },
    ) else {
        unreachable!()
    };
    // The boss's top cap, found by geometry rather than by index — the solid has been copied,
    // translated and mirrored since it was built, so no remembered handle survives.
    let corners = |f: Handle<Face>| -> Vec<[f64; 3]> {
        m.faces
            .get(f)
            .outer
            .half_edges
            .iter()
            .map(|he| m.vertex_point(m.edges.get(he.edge).vertices[0]).as_array())
            .collect()
    };
    let cap = *faces_of(&m, b)
        .iter()
        .filter(|&&f| {
            let c = corners(f);
            let z = c[0][2];
            c.iter().all(|p| p[2] == z)
        })
        .max_by(|&&x, &&y| corners(x)[0][2].total_cmp(&corners(y)[0][2]))
        .expect("the boss has a top cap");
    let c = corners(cap);
    let (x0, x1) = (
        c.iter().map(|p| p[0]).fold(f64::MAX, f64::min),
        c.iter().map(|p| p[0]).fold(f64::MIN, f64::max),
    );
    let (y0, y1) = (
        c.iter().map(|p| p[1]).fold(f64::MAX, f64::min),
        c.iter().map(|p| p[1]).fold(f64::MIN, f64::max),
    );
    let (mx, my) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    run(
        &mut m,
        Operation::PocketOnFace {
            face: cap,
            profile: rect(mx - 0.25, my - 0.25, mx + 0.25, my + 0.25),
            dist: 0.5, // into a 1-thick boss: blind
        },
    );
    let live: Vec<_> = m.live_solids.to_vec();
    run(
        &mut m,
        Operation::Boolean {
            kind: BoolKind::Fuse,
            a: live[0],
            b: live[1],
        },
    );
    m.rebuild_adjacency();

    // All six carried a handle.
    let carried = log
        .iter()
        .filter(|o| !matches!(o, Operation::Extrude { .. }))
        .count();
    assert_eq!(carried, 6, "this fixture is supposed to exercise all six");

    let replayed = replay(&log).expect("a log of every handle-carrying variant replays");
    assert_same_arena(&replayed, &m, "all-six");
}

// ─────────────────────────────────────────────────────────────────────────────
// What a reject leaves behind
// ─────────────────────────────────────────────────────────────────────────────

fn arena_lengths(m: &Model) -> [usize; 5] {
    [
        m.vertices.len(),
        m.edges.len(),
        m.faces.len(),
        m.shells.len(),
        m.solids.len(),
    ]
}

/// ★ **A reject is atomic on the live model, but not on the arena.**
///
/// `ops.rs` says it in one line — "only `live_solids` is touched, the store stays append-only" —
/// and that line is the whole hazard: an operation that builds a tool prism and *then* declines
/// has already grown the stores. The declined operation is not in the log, so a session that
/// keeps recording afterwards hands `replay` indices it can never reproduce.
///
/// This measures the boundary rather than asserting it from memory: early rejects (decided from
/// the request alone) leave Δ = 0, late rejects (decided from geometry that had to be built)
/// leave Δ > 0. The numbers print so `docs/dev-log.md` quotes a measurement, not a belief.
#[test]
fn a_late_reject_is_not_index_neutral() {
    let mut early = 0usize;
    let mut late = 0usize;

    let probe = |name: &str, setup: &dyn Fn() -> (Model, Operation)| -> Option<usize> {
        let (mut m, op) = setup();
        let before = arena_lengths(&m);
        let live_before: Vec<_> = m.live_solids.to_vec();
        let err = apply(&mut m, &op).err()?;
        // Every reject, early or late, leaves the live model exactly as it was — this is the
        // half of the contract `ops.rs` names ("no reject-after-commit") and the half that is
        // repaired. What follows measures the half that is not.
        assert_eq!(
            m.live_solids.to_vec(),
            live_before,
            "{name}: a reject must not change the live model"
        );
        let after = arena_lengths(&m);
        let delta: usize = before.iter().zip(&after).map(|(b, a)| a - b).sum();
        let d: Vec<String> = before
            .iter()
            .zip(&after)
            .map(|(b, a)| (a - b).to_string())
            .collect();
        println!(
            "stat reject_delta {name} total={delta} v/e/f/sh/so=+{} err={err:?}",
            d.join("/+")
        );
        Some(delta)
    };

    let seeded = || {
        let mut m = Model::new();
        let __seed = extrude_op(&m, 0.0, 2.0, 1.0);
        let OpOutput::Extrude { solid, faces } = apply(&mut m, &__seed).expect("seed") else {
            unreachable!()
        };
        (m, solid, faces)
    };

    // ── Early: decided before any geometry is built ─────────────────────────
    for (name, setup) in [
        (
            "NonPositiveDistance",
            Box::new(|| {
                let m = Model::new();
                let op = extrude_op(&m, 0.0, 2.0, 0.0);
                (m, op)
            }) as Box<dyn Fn() -> _>,
        ),
        (
            "SolidNotLive",
            Box::new(|| {
                let (mut m, solid, _) = seeded();
                apply(&mut m, &Operation::Copy { solid }).expect("copy");
                apply(
                    &mut m,
                    &Operation::Transform {
                        solid,
                        isometry: Isometry::translation([Rat::from_int(1); 3]),
                    },
                )
                .expect("supersede it");
                (
                    m,
                    Operation::Transform {
                        solid, // now superseded
                        isometry: Isometry::translation([Rat::from_int(1); 3]),
                    },
                )
            }),
        ),
    ] {
        let d = probe(name, &setup).unwrap_or_else(|| panic!("{name} was supposed to reject"));
        assert_eq!(
            d, 0,
            "{name} is an early reject and must not grow the arena"
        );
        early += 1;
    }

    // ── Late: the geometry had to exist before the answer was known ─────────
    for (name, setup) in [
        (
            "PadMissesFace",
            Box::new(|| {
                let (m, _, faces) = seeded();
                let op = Operation::PadOnFace {
                    face: faces[1],
                    profile: rect(20.0, 20.0, 21.0, 21.0), // nowhere near the cap
                    dist: 1.0,
                };
                (m, op)
            }) as Box<dyn Fn() -> _>,
        ),
        (
            "PocketNotBlind",
            Box::new(|| {
                let (m, _, faces) = seeded();
                let op = Operation::PocketOnFace {
                    face: faces[1],
                    profile: rect(0.5, 0.5, 1.5, 1.5),
                    dist: 5.0, // straight through the 1-thick block
                };
                (m, op)
            }),
        ),
    ] {
        let d = probe(name, &setup).unwrap_or_else(|| panic!("{name} was supposed to reject"));
        assert!(
            d > 0,
            "{name} is a late reject and is expected to leave residue"
        );
        late += 1;
    }

    // A boolean reject too, so all three late families back the "no reject-after-commit"
    // claim rather than two of them. Corner coincidence is genuinely non-manifold and the
    // engine declines it (cf. `a_corner_coincident_cut_is_rejected_not_silently_wrong`).
    // `add_cuboid` puts cells in the arena outside any log, which is exactly why this probe is
    // not a replay test — it only asks what a reject does to the model in front of it.
    let boolean_setup = || {
        let mut m = Model::new();
        let a = m.add_cuboid(
            nacre_math::Point3::from_array([0.0; 3]),
            nacre_math::Point3::from_array([2.0; 3]),
        );
        let b = m.add_cuboid(
            nacre_math::Point3::from_array([1.0; 3]),
            nacre_math::Point3::from_array([3.0; 3]),
        );
        let OpOutput::Boolean { solids } = apply(
            &mut m,
            &Operation::Boolean {
                kind: BoolKind::Cut,
                a,
                b,
            },
        )
        .expect("the first cut is fine") else {
            unreachable!()
        };
        // C's + corner is R's concave corner: three shared planes meet there.
        let c = m.add_cuboid(
            nacre_math::Point3::from_array([-1.0; 3]),
            nacre_math::Point3::from_array([1.0; 3]),
        );
        (
            m,
            Operation::Boolean {
                kind: BoolKind::Cut,
                a: solids[0],
                b: c,
            },
        )
    };
    match probe("Boolean(corner-coincident Cut)", &boolean_setup) {
        Some(d) => {
            assert!(d > 0, "a boolean reject builds geometry before deciding");
            late += 1;
        }
        None => panic!("the corner-coincident cut was supposed to be declined"),
    }

    println!("stat reject_delta_summary early={early} late={late}");
}

/// ★★ **The hazard, executed** — and the reason the re-sync in [`run_recipe`] is a rule and not
/// a convenience.
///
/// A session that hits a late reject and keeps recording *without* rebuilding from its log is
/// naming cells whose indices include residue the log does not contain. Exactly two outcomes
/// are acceptable, and this test enumerates them: a **named reject**, or a **success whose arena
/// diverges**. A panic is not on the list — the kernel declines what it cannot do (see
/// `docs/overview.md`), and an off-by-a-few log is squarely "cannot do".
#[test]
fn a_session_that_keeps_recording_after_a_late_reject_diverges() {
    let mut m = Model::new();
    let mut log = Vec::new();

    let seed = extrude_op(&m, 0.0, 2.0, 1.0);
    let OpOutput::Extrude { solid, faces } = apply(&mut m, &seed).expect("seed") else {
        unreachable!()
    };
    log.push(seed);

    // A late reject: builds the tool prism, then declines. Not recorded — but its cells stay.
    let before = arena_lengths(&m);
    // `PadMissesFace` is used deliberately: it is the late reject that *does* restore
    // `live_solids` (`ops.rs`, "no reject-after-commit"), so what is left is arena residue and
    // nothing else. The reject that fails to restore is a separate defect, witnessed by
    // [`a_pocket_that_is_not_blind_rejects_after_committing`].
    let declined = Operation::PadOnFace {
        face: faces[1],
        profile: rect(20.0, 20.0, 21.0, 21.0),
        dist: 1.0,
    };
    let err = apply(&mut m, &declined).expect_err("this pad misses the face");
    let residue: usize = arena_lengths(&m)
        .iter()
        .zip(&before)
        .map(|(a, b)| a - b)
        .sum();
    assert!(
        residue > 0,
        "the premise of this test is residue; got none ({err:?})"
    );

    // The session keeps going *without* re-syncing, and records what follows.
    let cp = Operation::Copy { solid };
    let OpOutput::Copy { solid: copy } = apply(&mut m, &cp).expect("copy") else {
        unreachable!()
    };
    log.push(cp);
    let mv = Operation::Transform {
        solid: copy,
        isometry: Isometry::translation([Rat::from_int(5), Rat::from_int(0), Rat::from_int(0)]),
    };
    apply(&mut m, &mv).expect("move the copy");
    log.push(mv);
    m.rebuild_adjacency();

    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| replay(&log)));
    match outcome {
        Err(_) => panic!(
            "replaying a log recorded past a late reject PANICKED — the kernel must decline, \
             not abort (docs/overview.md). This is a defect, not an acceptable outcome."
        ),
        Ok(Err(e)) => println!("stat past_reject_outcome named_reject {e:?}"),
        Ok(Ok(r)) => {
            let (sa, sb) = (arena_sig(&r), arena_sig(&m));
            let first = sa.iter().zip(&sb).find(|(x, y)| x != y);
            match first {
                Some((x, y)) => println!(
                    "stat past_reject_outcome diverged at {}[{}] {:?} vs {:?}",
                    x.what, x.at, x.value, y.value
                ),
                None if sa.len() != sb.len() => {
                    println!(
                        "stat past_reject_outcome diverged in length {} vs {}",
                        sa.len(),
                        sb.len()
                    )
                }
                None => panic!(
                    "replay AGREED with a session that recorded past a late reject. That is a \
                     finding, not a pass: the residue ({residue} cells) did not move any index \
                     this log names. Narrow the claim in docs/design.md to what is measured."
                ),
            }
        }
    }
}

/// ★ **A reject does not commit — including the one that used to.**
///
/// `ops.rs` names this contract where `PadMissesFace` restores `live_solids`: "the model the
/// caller sees is the one it had before". `PocketNotBlind` used to break it — it was raised in
/// `pocket` *after* `extrude_and_boolean` returned `Ok`, by which point the boolean had already
/// retired the caller's solid and installed the through-cut result in its place, so the caller
/// got an `Err` **and** a different live model. Measured here first, then repaired by moving the
/// missing-cap verdict into the scope that still holds the result solids.
///
/// What a late reject still leaves is **arena residue** — the store is append-only and the
/// prism's cells stay. That is the honest remainder, and the reason a session must rebuild from
/// its log before recording again ([`a_session_that_keeps_recording_after_a_late_reject_diverges`]).
#[test]
fn a_pocket_that_is_not_blind_leaves_the_live_model_alone() {
    let mut m = Model::new();
    let __seed = extrude_op(&m, 0.0, 2.0, 1.0);
    let OpOutput::Extrude { solid, faces } = apply(&mut m, &__seed).expect("seed") else {
        unreachable!()
    };
    let live_before: Vec<_> = m.live_solids.to_vec();
    let arena_before = arena_lengths(&m);

    let err = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: faces[1],
            profile: rect(0.5, 0.5, 1.5, 1.5),
            dist: 5.0, // straight through the 1-thick block
        },
    )
    .expect_err("a through-cut is not a blind pocket");
    assert!(matches!(err, nacre_ops::OpError::PocketNotBlind));

    assert!(
        m.live_solids.contains(&solid),
        "the declined pocket must leave the caller's solid live"
    );
    assert_eq!(
        m.live_solids.to_vec(),
        live_before,
        "the live model after a reject is the one the caller had before it"
    );
    // The model still works: the solid can be pocketed for real.
    apply(
        &mut m,
        &Operation::PocketOnFace {
            face: faces[1],
            profile: rect(0.5, 0.5, 1.5, 1.5),
            dist: 0.5,
        },
    )
    .expect("a blind pocket on the restored solid");

    // The remainder that is *not* repaired, stated rather than implied.
    let grew: usize = arena_lengths(&m)
        .iter()
        .zip(&arena_before)
        .map(|(a, b)| a - b)
        .sum();
    assert!(
        grew > 0,
        "the append-only arena still keeps the declined prism's cells"
    );
}

/// ★ **The positive control for the repair: re-anchoring did not weaken the door.**
///
/// `replay` now launders a log's indices into its own model on purpose. The risk of a laundering
/// step is that it launders everything — so this checks the case it must *not* touch: a handle
/// from another model handed straight to `apply` still dies in `Store::get`'s cross-store guard.
/// `apply`'s model belongs to the caller, so its handles do too, and a re-anchor there would turn
/// a caller bug into a silently different answer.
///
/// Note what this test is really pinning: the guard is not the first thing the handle meets.
/// `live_solids.contains(&h)` compares by index and answers "yes, live" for the foreigner; the
/// panic comes one step later, at the first dereference. That asymmetry is deliberate and is
/// annotated at each `contains` site — `Handle` cannot carry the store in its `Eq` (it has no
/// `T: Eq` to lean on), and every hashed-handle map in `nacre-topo` depends on that.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "different Store")]
fn a_foreign_handle_still_dies_at_the_door() {
    let mut a = Model::new();
    let op_a = extrude_op(&a, 0.0, 2.0, 1.0);
    let OpOutput::Extrude { solid, .. } = apply(&mut a, &op_a).expect("a") else {
        unreachable!()
    };

    // `b` is built the same way, so `solid`'s *index* is live in `b` too — the handle is wrong
    // only in the one way that matters, and `contains` cannot see it.
    let mut b = Model::new();
    let op_b = extrude_op(&b, 0.0, 2.0, 1.0);
    apply(&mut b, &op_b).expect("b");
    assert!(
        b.live_solids.contains(&solid),
        "index-only equality: b agrees the foreign handle is live"
    );

    let _ = apply(
        &mut b,
        &Operation::Transform {
            solid, // minted by `a`
            isometry: Isometry::translation([Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)]),
        },
    );
}

/// ★ **An `Extrude` naming a surface that is not there is a named reject too.**
///
/// `Extrude` was the one variant `rebind` could borrow rather than re-anchor, because it carried
/// no handle. Since S5(i)-b it carries a `SketchFrame`, and forgetting to re-anchor it would have
/// put the defect R repaired back into the **most common** operation in any log. This is the test
/// that would have caught that.
#[test]
fn an_extrude_naming_a_surface_that_does_not_exist_is_rejected_by_name() {
    // A frame on a plane that only a datum creates — so a log that omits the datum names a
    // surface index the replayed model never reaches.
    let mut big = Model::new();
    let sp = nacre_ops::SketchPlane::from_origin_normal(
        nacre_math::Point3::from_array([0.0, 0.0, 1.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .expect("a plane");
    let datum = Operation::DatumPlane {
        def: nacre_ops::DatumDef::Stated(sp),
    };
    let OpOutput::DatumPlane { plane, frame } = apply(&mut big, &datum).expect("stated") else {
        unreachable!()
    };
    assert!(
        plane.index() >= 3,
        "the fixture needs a plane past the seeds, or the log would resolve by accident"
    );

    let orphan = vec![Operation::Extrude {
        frame,
        profile: square(0.0, 2.0),
        dist: 1.0,
    }];
    match replay(&orphan) {
        Err(nacre_ops::OpError::LogHandleOutOfRange { cell, index }) => {
            assert_eq!(cell, nacre_ops::LogCell::Surface);
            assert_eq!(index, plane.index());
        }
        other => panic!(
            "expected a named reject, got {:?}",
            other.map(|_| "a model")
        ),
    }

    // With the datum in front of it, the same frame replays.
    let full = vec![
        datum,
        Operation::Extrude {
            frame,
            profile: square(0.0, 2.0),
            dist: 1.0,
        },
    ];
    let replayed = replay(&full).expect("the whole log replays");
    assert_eq!(replayed.live_solids.len(), 1);
}

/// ★★★★★ **A vertex-naming datum, recorded after the session already rejected something.**
///
/// This is the one shape where re-anchoring could go wrong without any error: a late reject
/// leaves cells in the arena (measured here, and 63–84 in R's census), so a replay that rebuilt
/// from the log would find its vertices at *different indices* than the recording session did.
/// A bounds check cannot see that — the wrong index is still in range — so the only witness is a
/// session that contains both and an arena comparison afterwards.
///
/// The session follows the discipline (`docs/design.md`): after a reject it throws its model away
/// and rebuilds from the log, which is what makes the indices agree. The point of the test is
/// that a datum naming *vertices* is not an exception to it.
#[test]
fn a_datum_naming_vertices_survives_a_session_that_rejected() {
    let mut m = Model::new();
    let mut log: Vec<Operation> = Vec::new();

    let seed = extrude_op(&m, 0.0, 2.0, 1.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &seed).expect("seed") else {
        unreachable!()
    };
    log.push(seed);

    // A late reject — builds its tool prism, then declines. Its cells stay behind.
    let before = arena_lengths(&m);
    let err = apply(
        &mut m,
        &Operation::PadOnFace {
            face: faces[1],
            profile: rect(20.0, 20.0, 21.0, 21.0),
            dist: 1.0,
        },
    )
    .expect_err("this pad misses the face");
    let residue: usize = arena_lengths(&m)
        .iter()
        .zip(&before)
        .map(|(a, b)| a - b)
        .sum();
    assert!(
        residue > 0,
        "the premise of this test is residue; got none ({err:?})"
    );

    // The discipline: rebuild from the log, so the session's indices are the log's indices.
    // ★ The rebuilt model is a *different* arena, so `solid` from before it must not be reused —
    // the cross-store guard says so, loudly, and it is right to.
    m = replay(&log).expect("the log so far replays");
    let solid = m.live_solids[0];

    // Now name three corners of the seed solid — the step whose handles are vertices.
    let mut vs = Vec::new();
    let sol = m.solids.get(solid);
    for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
        for &fh in &m.shells.get(sh).faces {
            let f = m.faces.get(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for &he in &lp.half_edges {
                    let vh = m.he_start(he);
                    if !vs.contains(&vh) {
                        vs.push(vh);
                    }
                }
            }
        }
    }
    // Corners of one face are coplanar with it, so this datum interns onto a plane the box
    // already holds — which is fine here: what is under test is the *handles*, not the variant.
    let datum = Operation::DatumPlane {
        def: nacre_ops::DatumDef::ThroughVertices([vs[0], vs[1], vs[2]]),
    };
    apply(&mut m, &datum).expect("three corners of a box name a plane");
    log.push(datum);

    let once = replay(&log).expect("the log replays");
    assert_same_arena(&once, &m, "replay(log) vs the session that rejected");
}
