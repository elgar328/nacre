//! Fusing a part with rotated copies of itself — and the naming rule that made 30° fail.
//!
//! At 30° the arithmetic conspires: `sin 30° = ½` exactly, so a rotated copy's profile corner
//! lands *precisely* on one of the original's face planes, and four planes concur at that point.
//! The arrangement named it by its canonical (lexicographically first) triple, which included a
//! plane the union keeps no face on — the copy's wall, buried inside the material. The fuse then
//! succeeded while returning a solid that could not describe itself, and the failure surfaced two
//! operations later, when the next `rotateZ` refused to remap a definition pointing at a surface
//! this solid has no face on.
//!
//! So the propositions here are about **what a result says about itself**, not about volume.

use nacre_math::Point2;
use nacre_ops::{
    BoolError, BoolKind, Operation, Profile2d, RejectReason, SketchFrame, apply, boolean,
};
use nacre_scalar::{Angle, Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// The playground script's part: an L-ish plate with a notch, plus a bar across it.
fn unit(m: &mut Model) -> Handle<Solid> {
    let plate = [
        [0.0, 0.0],
        [50.0, 0.0],
        [50.0, 25.0],
        [38.0, 25.0],
        [38.0, 50.0],
        [50.0, 50.0],
        [50.0, 75.0],
        [0.0, 75.0],
    ];
    let p1 = prism(m, &plate, Axis::Z, 12.0);
    let bar = [[20.0, 12.0], [75.0, 12.0], [75.0, 37.0], [55.0, 37.0]];
    let p2 = prism(m, &bar, Axis::X, 25.0);
    let out = boolean(m, BoolKind::Fuse, p1, p2).expect("the part fuses")[0];
    m.rebuild_adjacency();
    out
}

fn prism(m: &mut Model, pts: &[[f64; 2]], axis: Axis, dist: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let op = Operation::Extrude {
        frame: SketchFrame::world(m, axis),
        profile,
        dist,
    };
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(m, &op).expect("extrude") else {
        panic!("extrude output")
    };
    m.rebuild_adjacency();
    solid
}

fn copy_of(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
    let nacre_ops::OpOutput::Copy { solid } =
        apply(m, &Operation::Copy { solid: s }).expect("copy")
    else {
        panic!("copy output")
    };
    m.rebuild_adjacency();
    solid
}

fn rot_z(m: &mut Model, s: Handle<Solid>, deg: i128) -> Result<Handle<Solid>, nacre_ops::OpError> {
    let out = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            }),
        },
    )?;
    m.rebuild_adjacency();
    match out {
        nacre_ops::OpOutput::Transform { solid } => Ok(solid),
        o => panic!("{o:?}"),
    }
}

/// Every vertex definition of `solid` that names a surface this solid has no face on — the
/// property `transform` needs and the one the old naming broke. Walked through the public model
/// so the test sees what any consumer would.
fn foreign_definitions(m: &Model, solid: Handle<Solid>) -> Vec<String> {
    use std::collections::HashSet;
    let s = m.solids.get(solid);
    let shells: Vec<_> = std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .collect();
    let mut mine: HashSet<_> = HashSet::new();
    for &sh in &shells {
        for &fh in &m.shells.get(sh).faces {
            mine.insert(m.faces.get(fh).surface);
        }
    }
    let mut bad = Vec::new();
    let mut seen: HashSet<u32> = HashSet::new();
    for &sh in &shells {
        for &fh in &m.shells.get(sh).faces {
            let face = m.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for vh in m.edges.get(he.edge).vertices.iter() {
                        if !seen.insert(vh.index()) {
                            continue;
                        }
                        let names: Vec<_> = m.vertices.get(*vh).def.carriers().collect();
                        let missing: Vec<u32> = names
                            .iter()
                            .filter(|s| !mine.contains(s))
                            .map(|s| s.index())
                            .collect();
                        if !missing.is_empty() {
                            let at = m.vertex_point(*vh).as_array();
                            bad.push(format!(
                                "vertex {} at {at:?} names surfaces {missing:?}",
                                vh.index()
                            ));
                        }
                    }
                }
            }
        }
    }
    bad
}

/// ① The fuse that started it: it always built, and what it built was unusable. Both halves are
/// asserted, because the second is the one that was silently false.
#[test]
fn a_fused_result_names_itself_in_its_own_surfaces() {
    let mut m = Model::new();
    let u = unit(&mut m);
    let c = copy_of(&mut m, u);
    let c = rot_z(&mut m, c, 30).expect("rotating a copy");
    let fused = boolean(&mut m, BoolKind::Fuse, u, c).expect("30° fuses")[0];
    m.rebuild_adjacency();
    assert_eq!(
        foreign_definitions(&m, fused),
        Vec::<String>::new(),
        "a result must be describable in the surfaces it actually has"
    );
    // …and the consequence, which is where the failure used to appear: the next rotation of that
    // result is exactly the operation that returned `OriginNotOnSolid`.
    rot_z(&mut m, fused, 30).expect("the result can be moved");
}

/// Fuse `unit` with a copy of itself every `step` degrees, all the way round, and say where it
/// stopped. The template stays live: a boolean supersedes its operands, so `u` itself must never
/// be one of them or the next copy has nothing to copy.
fn sweep(step: usize) -> Result<(), (usize, BoolError)> {
    let mut m = Model::new();
    let u = unit(&mut m);
    let mut part = copy_of(&mut m, u);
    for deg in (step..360).step_by(step) {
        let c = copy_of(&mut m, u);
        let c = rot_z(&mut m, c, deg as i128).expect("rotating a copy");
        match boolean(&mut m, BoolKind::Fuse, part, c) {
            Ok(out) => {
                // ★ Taking `out[0]` without asking how many there are would let a fold that starts
                // severing carry on against a *fragment* — green, and measuring something else.
                assert_eq!(out.len(), 1, "the {step}° fold split at the {deg}° copy");
                part = out[0];
                m.rebuild_adjacency();
            }
            Err(e) => return Err((deg, e)),
        }
    }
    Ok(())
}

/// ★ The band that one arbitrary probe vertex used to cost us.
///
/// "Is this hole inside that ring" is answered by casting a ray from a vertex, and a ring node on
/// the ray's line makes the parity ambiguous, so that candidate is dropped. The cleaning pass
/// probed the hole's *first* node and gave up if it was spoiled — while the arrangement, asking
/// the same question, tried every node. A small change of angle leaves the topology alone, so the
/// same unlucky first node persisted across a whole band of angles and took all of them down.
///
/// 41° is the cheapest angle that used to fail, so it is the one the default suite runs; the
/// rest of the band lives in the sweep below, which is slow enough to be opt-in.
#[test]
fn the_band_of_angles_a_single_probe_used_to_cost() {
    if let Err((deg, e)) = sweep(41) {
        panic!("the 41° sweep died at the {deg}° copy: {e:?}");
    }
}

/// ★ The angles the retry alone could not save — where **every** candidate ray was blocked.
///
/// A ring node sitting on the ray's line used to be abandoned rather than judged, and at 38–40°
/// the arrangement puts one on every line a probe can cast along: five isolated corners and one
/// whole edge, measured. The node's two off-line neighbours settle it — opposite sides is a
/// crossing, equal sides a touch — which is the rule the tracer had all along.
///
/// This carries the *second* fix the way the test above carries the first; they lock different
/// things, so both stay in the default suite.
#[test]
fn the_angles_where_no_ray_was_left_unblocked() {
    if let Err((deg, e)) = sweep(38) {
        panic!("the 38° sweep died at the {deg}° copy: {e:?}");
    }
}

/// The whole band, plus the boundaries that never failed and must not start (37°, 46°) and the
/// long 20° sweep — seventeen chained fuses, each on the result of the last.
///
/// `#[ignore]`: about eighty booleans on a solid that grows with every one, minutes in a debug
/// build. The population is what the two fixes were judged by, so it is written down; the default
/// suite carries one angle from each.
#[test]
#[ignore = "slow angle sweep (run with --ignored)"]
fn the_whole_band_of_angles() {
    let mut died = Vec::new();
    for step in [20, 37, 38, 39, 40, 41, 42, 43, 44, 46] {
        if let Err((deg, e)) = sweep(step) {
            died.push(format!("{step}° at the {deg}° copy: {e:?}"));
        }
    }
    assert!(died.is_empty(), "the failing population changed: {died:?}");
}

/// ② The whole 30° sweep — eleven chained fuses, each on the result of the last.
#[test]
fn the_thirty_degree_sweep_runs_to_completion() {
    let mut m = Model::new();
    let u = unit(&mut m);
    // The template stays live: a boolean supersedes its operands, so `u` itself must never be
    // one of them or the next copy has nothing to copy.
    let mut part = copy_of(&mut m, u);
    for deg in (30..360).step_by(30) {
        let c = copy_of(&mut m, u);
        let c = rot_z(&mut m, c, deg as i128).expect("rotating a copy");
        let out = boolean(&mut m, BoolKind::Fuse, part, c)
            .unwrap_or_else(|e| panic!("fusing the {deg}° copy: {e:?}"));
        assert_eq!(out.len(), 1, "the star stays one body at {deg}°");
        part = out[0];
        m.rebuild_adjacency();
        assert_eq!(
            foreign_definitions(&m, part),
            Vec::<String>::new(),
            "after the {deg}° fuse"
        );
    }
}

/// ③ ★ The negative control: a step angle that rejects for an *other* reason, which renaming
/// vertices must not touch. If it moves, the fix reached further than its argument says it does.
///
/// ★★ It used to carry `120°` as well, on the reading "two bodies meeting along one line — no
/// 2-manifold contains it". Two bodies meeting along one line are two bodies now, and that first
/// copy comes back as two of them (`the_hundred_and_twenty_degree_copy_is_two_bodies`). The
/// reading was right about the *single* body the reconstruction used to weld, and wrong about the
/// answer.
///
/// ★★★ **45° says its truth now** (2026-08-17): the arm landing coplanar on its own body is a
/// self-touch, and since the merge abstains on a pinching group (nothing to re-thread) and the
/// whole-result judgement runs before minting, the reject is `SelfTouchingResult` — `Impossible`,
/// with the touching edge itself as the witness — where for two days it was `CoplanarPinch`
/// (`NotSupported`, a capability-limit name raised by a site that could not see the whole result).
#[test]
fn the_other_rejections_are_untouched() {
    let cases: [(usize, RejectReason); 1] = [(45, RejectReason::SelfTouchingResult)];
    for (step, expected) in cases {
        let mut m = Model::new();
        let u = unit(&mut m);
        let mut part = copy_of(&mut m, u);
        let mut got = None;
        for deg in (step..360).step_by(step) {
            let c = copy_of(&mut m, u);
            let Ok(c) = rot_z(&mut m, c, deg as i128) else {
                panic!("{step}°: a copy stopped being movable — that is this fix's own subject")
            };
            match boolean(&mut m, BoolKind::Fuse, part, c) {
                Ok(out) => {
                    part = out[0];
                    m.rebuild_adjacency();
                }
                Err(e) => {
                    got = Some(e);
                    break;
                }
            }
        }
        let got = got.map(|e| match e {
            BoolError::Rejected { reason, at } => (reason, at),
            other => panic!("expected a named rejection, got {other:?}"),
        });
        let (reason, at) = got.expect("the sweep was expected to reject");
        assert_eq!(reason, expected, "the {step}° sweep's rejection changed");
        // The witness is the touching edge itself. Every fuse here rotates about Z, so the
        // self-touch is the vertical line the coplanar pinch extrudes along: its endpoints sit
        // on the two horizontal caps (z = 0 and z = 12, the fixture's own height) and share
        // their x/y — a wrongly-realized segment would break one of those. (Pinning the exact
        // x/y would re-derive the fold; the algebra pins the vertical-line shape only.)
        let Some(nacre_ops::RejectWhere::Segment([a, b])) = at else {
            panic!("the {step}° self-touch carries no witness segment: {at:?}");
        };
        let (mut lo, mut hi) = (a[2], b[2]);
        if lo > hi {
            std::mem::swap(&mut lo, &mut hi);
        }
        assert!(
            lo.abs() < 1e-9 && (hi - 12.0).abs() < 1e-9,
            "the {step}° touch segment does not span the caps: {a:?} {b:?}"
        );
        assert!(
            (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9,
            "the {step}° touch segment is not vertical: {a:?} {b:?}"
        );
    }
}

/// ④ ★ **The 120° copy, which used to be a reject.** It meets the part along one line and nowhere
/// else, so the fuse is the two bodies it was handed — on a real part, not a pair of cubes.
///
/// ★★ Written as its own statement rather than through `sweep`, because `sweep` is a *fold* and
/// asserts its result stays one solid. That assertion is what caught the first draft of this test,
/// which claimed the fold ran to completion: it does not, and the first step is why.
#[test]
fn the_hundred_and_twenty_degree_copy_is_two_bodies() {
    let mut m = Model::new();
    let u = unit(&mut m);
    let one = nacre_props::mass_props(&m, u).expect("props").volume;
    let part = copy_of(&mut m, u);
    let c = copy_of(&mut m, u);
    let c = rot_z(&mut m, c, 120).expect("rotating a copy");
    let out = boolean(&mut m, BoolKind::Fuse, part, c).expect("a line contact separates");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 2, "the part and its 120° copy only touch");
    let total: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    assert!(
        (total - 2.0 * one).abs() < 1e-9,
        "nothing was removed: {total} against {}",
        2.0 * one
    );
    assert!(nacre_validate::validate(&m).is_empty());
    for &s in &out {
        assert_eq!(foreign_definitions(&m, s), Vec::<String>::new());
    }
}

/// ★ One implementation, and it stays one.
///
/// The retry that answers "is this ring inside that one" used to live in the callers — spelled
/// two ways in the arrangement, missing entirely in the cleaning pass. It lives in
/// `ring_in_ring` now, and a caller that goes around it is a caller that will quietly lack the
/// retry again, so the source says so.
///
/// ★ **The exceptions are listed, not implied.** `ring_in_ring`'s value is picking *another node*
/// when one grazes, so the rule binds anything asking "is this ring inside that one". A caller
/// asking about **one named vertex** has no other node to offer, and the retry would have nothing
/// to retry with — `inside_trimmed_face` (the self-touch check) is that, and it is spelled out
/// below rather than left to a loose pattern.
#[test]
fn no_production_caller_reaches_past_the_shared_predicate() {
    let src = [
        "src/combinatorics.rs",
        "src/arrangement.rs",
        "src/boolean.rs",
    ];
    let mut offenders = Vec::new();
    for file in src {
        let text = std::fs::read_to_string(file).expect("source file");
        for (n, line) in text.lines().enumerate() {
            let calls = line.contains("point_in_ring(");
            let is_definition = line.contains("fn point_in_ring");
            let is_prose = line.trim_start().starts_with("//");
            if calls && !is_definition && !is_prose {
                offenders.push(format!("{file}:{}: {}", n + 1, line.trim()));
            }
        }
    }
    // `ring_in_ring` is the one ring-vs-ring caller, and `point_in_component` casts its own rays
    // in 3D. `inside_trimmed_face` asks about a single vertex — see this test's note.
    offenders.retain(|o| {
        !o.contains("point_in_ring(jd, p, v, outer)") && !o.contains("point_in_ring(jd, q, t, ")
    });
    assert_eq!(
        offenders,
        Vec::<String>::new(),
        "call `ring_in_ring` instead — the retry lives there"
    );
}

/// ★ …and so does the walk underneath it.
///
/// "Where does this ring meet that line, and does it cross or only touch" was answered in two
/// places — the tracer had the flank rule, the ray caster threw such a candidate away — which is
/// the same shape the retry above was in. It is `ring_against_plane` now. Reading `side_of` over a
/// **ring** anywhere else is a second walk being born, so the source says so.
#[test]
fn no_production_code_walks_a_ring_past_the_shared_walk() {
    let src = [
        "src/combinatorics.rs",
        "src/arrangement.rs",
        "src/boolean.rs",
    ];
    let mut offenders = Vec::new();
    for file in src {
        let text = std::fs::read_to_string(file).expect("source file");
        for (n, line) in text.lines().enumerate() {
            let calls = line.contains("side_of(");
            let is_definition = line.contains("fn side_of");
            let is_prose =
                line.trim_start().starts_with("//") || line.trim_start().starts_with("///");
            if calls && !is_definition && !is_prose {
                offenders.push(format!("{file}:{}: {}", n + 1, line.trim()));
            }
        }
    }
    // The walk itself, plus the two places that ask about a **single point** rather than a ring:
    // `point_on_ring` ("is `v` on this edge's line") and the alias seed ("does this vertex lie on
    // that class"). Neither reads a sign sequence, so neither is a walk.
    offenders.retain(|o| {
        !o.contains("side_of(jd, nodes[i], q)")
            && !o.contains("side_of(jd, v, r)")
            && !o.contains("side_of(&jd, sorted3(t), wc)")
    });
    assert_eq!(
        offenders,
        Vec::<String>::new(),
        "call `ring_against_plane` instead — the flank rule lives there"
    );
}
