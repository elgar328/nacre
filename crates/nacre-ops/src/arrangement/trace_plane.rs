use super::*;
/// Why a ring of names could not be carried.
///
/// **One variant.** A corner with no three-plane name and a carrier riding a cylinder are not
/// *decisions*: a ring spelled as plane ids cannot **describe** a curved ring, so refusing there
/// would dress a vocabulary limit up as a cause, and the ring vocabulary describes both.
///
/// What is left is a name that is degenerate **as a name**, which is a fact about the triple and
/// has nothing to do with cylinders — so it stands alone and says only itself.
#[derive(Debug)]
enum RingFail {
    /// Two of a vertex's three planes coincide, so the name denotes no point.
    Collapsed,
}

/// The decline kind for a [`RingFail`] at a caller that has no finer fact to add.
///
/// Two of the four callers *do* have one — they know which loop of the face failed and report
/// `OuterRing`/`HoleRing` — so this is not "the mapping", only the default.
fn decline_of(f: RingFail) -> DeclineKind {
    match f {
        RingFail::Collapsed => DeclineKind::CollapsedTriple,
    }
}

/// **The one place a ring of *names* becomes a ring of *plane data*** for the tracer, which asks
/// only "which side, which wall, which third plane" of each vertex. There is no sort here —
/// [`NodeId`]'s only constructor sorts — so what this does is the collapse check, and that is the
/// whole reason this function exists.
fn plane_ring(
    nr: &combinatorics::NamedRing,
) -> Result<(Vec<combinatorics::NodeId>, Vec<crate::combinatorics::Wall>), RingFail> {
    // **No refusal here, because the vessel can describe a curved ring.** Projected to
    // `[usize; 3]` corners and plane-class carriers, a ring a cylinder touches would not be
    // *described*, and a refusal here would be the type running out, dressed as a decision. What
    // a curved ring meets is the roads themselves, which answer or decline on their own terms.
    let ts: Vec<combinatorics::NodeId> = nr
        .triples
        .iter()
        .map(|&n| {
            // A collapsed name is still a collapsed name — but only a three-plane one can be.
            if let Some(c) = three_plane_name(n) {
                // Sorted by construction, so equal neighbours catch every duplicate.
                if c[0] == c[1] || c[1] == c[2] {
                    return Err(RingFail::Collapsed);
                }
            }
            Ok(n)
        })
        .collect::<Result<_, _>>()?;
    Ok((ts, nr.walls.clone()))
}

/// Trace one **transversal or mixed** face `f` (plane `fp ≠ W`) onto plane class `wc`, as the
/// segments where `f` meets `W`. This is 2D polygon-vs-line clipping on `L = W ∩ fp`, handling
/// on-line edges (a run of `side == 0` vertices) instead of rejecting them: the rings' features
/// become stations on `L` and [`sweep_line`] reads what the face leaves there. A **holed** face's
/// inner rings join the same station list, and each hole's two crossings of `L` toggle the face
/// off between them, carving the hole out of the emitted chord.
///
/// Every emitted segment rides `wall = fp` (its line is `W ∩ fp`); the side walls appear only as
/// the third plane naming each endpoint. `mat = orient_sign(fp)` is a per-face constant.
#[allow(clippy::too_many_arguments)]
fn trace_transversal_face(
    fp: usize,
    loops: &combinatorics::FaceLoops,
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    out: &mut Trace,
) {
    // `fp` names a *face* (`n_out`, `orient`, the declined log); `fc` names the *plane class* it
    // lies on (triples, comparisons, predicate arguments). Every ring below is in class form, so
    // the two must not be confused.
    let fc = plane_ix[fp].plane();
    // Each ring pairs its class-form triples with the **carried walls** the producer read off
    // the model's edges (`NamedRing`) — the Crossing arm below names an edge's wall from the
    // ride, never from the two endpoint names.
    // Every station this face leaves on `L`. A **circle** of the face leaves its chord's two ends
    // here first, for the outer and the holes alike; the rings' features join them below.
    let mut stations: Vec<ClassLineStation> = Vec::new();
    // The outer ring the scan walks — `None` for a circular outer, whose two chord ends are
    // already in `stations` and which has no polygon to walk.
    let outer: Option<(Vec<combinatorics::NodeId>, Vec<crate::combinatorics::Wall>)> =
        match loops.outer.as_ref().and_then(|r| r.poly()) {
            Some(nr) => match plane_ring(nr) {
                Ok(pair) => Some(pair),
                Err(f) => {
                    out.declined.push((fp, decline_of(f)));
                    return;
                }
            },
            // A **circular outer is a disk — and a disk can have holes.** A recorded wall
            // class cuts the circle in a chord (a diameter when it runs through the axis), and
            // the chord's two ends are stations of the same sweep every polygon face runs, so the
            // face's **holes carve the chord** exactly as they carve a polygon's
            // section: a tube's annular cap sectioned by a wall is two pieces, not the whole
            // chord across the bore. Emitting the whole chord and *returning before reading
            // the holes* («a circular outer is a disk face») is wrong for any solid with two
            // coaxial cylinders (annulus, tube): the tube's bore comes through as a chord piece
            // the cap does not cover and the class refuses `LabelConflict` (measured, the
            // annulus × box rows). One road for both.
            //
            // Any other class misses or is left silent: another ⊥ class is parallel to the disk's
            // own and meets it nowhere, a ∥-axis wall the gate did not record either clears the
            // disk (the gate's plane test) or cuts an **irrational** chord this road cannot state
            // yet (the face-cleared family — its wrongly-silent trace lands in cells no face
            // reaches, today's recorded state), and a class with no world statement is one the
            // gate proved clear of the whole cylinder. A disk the wall clears has holes the wall
            // clears too, so the silence covers the holes as well.
            //
            // `None` is also how an *unnamed* ring arrives, but not for a face the tracer
            // reaches: the only loops `loop_triples` leaves unnamed are the ones it rejects for,
            // and those are `Err` at the source. A circle is the one `None` that means "named,
            // and not a polygon".
            None if matches!(loops.outer, Some(combinatorics::LoopRing::Circle { .. })) => {
                let Some(combinatorics::LoopRing::Circle { cyl }) = loops.outer else {
                    unreachable!("the guard just matched a circle outer");
                };
                match chord_stations(jd, faces, plane_ix, wc, fc, cyl, cyls) {
                    Ok(Some(pair)) => stations.extend(pair),
                    Ok(None) => return,
                    Err(ChordUnstatable) => {
                        out.declined.push((fp, DeclineKind::Ruling));
                        return;
                    }
                }
                None
            }
            None => {
                out.declined.push((fp, DeclineKind::OuterRing));
                return;
            }
        };
    let mut holes: Vec<(Vec<combinatorics::NodeId>, Vec<crate::combinatorics::Wall>)> = Vec::new();
    // A hole whose ring cannot be named is not "no hole" — swallowing the error would trace the
    // face as solid where it is pierced, which is a silent wrong answer rather than a reject.
    let Some(raw_holes) = &loops.holes else {
        out.declined.push((fp, DeclineKind::HoleRing));
        return;
    };
    for r in raw_holes {
        match r {
            // A lateral's own rim is a cycle of a lateral row (`FaceLoops::cycles`), never a
            // plane face's hole: reaching here is a producer inconsistency, declined by the
            // loop it arrived as.
            combinatorics::LoopRing::Rim { .. } => {
                out.declined.push((fp, DeclineKind::HoleRing));
                return;
            }
            // A **circular hole meets the line in a chord, and then it toggles the face twice**.
            // The population gate admits a ∥-axis wall clear of the hole's cylinder
            // (skip: the circle cannot meet `L`) or **recorded** — within the radius, through
            // the axis or offset from it — and there the line cuts the hole in the chord's two
            // pierce points. Skipping them "by proof" is wrong: the proof covers the clear half
            // only, and the full-width chord then surfaces as `LabelConflict` on the
            // through-family × through-axis tool (measured). The two roots are pushed as
            // stations — the hole's own chord, named the same way a disk outer's is
            // ([`chord_stations`], the one spelling for both).
            //
            // **The gate's record, the same one the rulings road reads first.** A pair the
            // gate did not list was proven clear — no emitted face reaches the hole's footprint
            // on this class, and the rulings road is silent for it — so planting stations
            // would break the line against rulings that are (rightly) absent: measured, the
            // d = 0 boss-and-bore fixture walked its spur out and back and refused
            // `StraightAngle`. A cut pair gets the stations; any other pair is skipped, and the
            // two roads stay one rule.
            combinatorics::LoopRing::Circle { cyl } => {
                match chord_stations(jd, faces, plane_ix, wc, fc, *cyl, cyls) {
                    Ok(Some(pair)) => stations.extend(pair),
                    Ok(None) => {}
                    Err(_) => {
                        out.declined.push((fp, DeclineKind::HoleRing));
                        return;
                    }
                }
                continue;
            }
            combinatorics::LoopRing::Poly(r) => match plane_ring(r) {
                Ok(pair) => holes.push(pair),
                Err(f) => {
                    out.declined.push((fp, decline_of(f)));
                    return;
                }
            },
        }
    }

    // `L`'s third plane naming a point of the on-line edge: the vertex triple `{fc, W-class, r}`.
    //
    // Failing to name it has two *different* causes, and the caller must be able to tell them
    // apart. `t`'s three classes are distinct (`plane_ring` rejects a collapsed triple first) and
    // `fc != wc` is this function's precondition, so:
    //
    // - `wc ∉ t` — every caller passes a **run** vertex, whose point lies on `wc` (that is what
    //   `side == 0` says). A point on `wc` that `wc` does not name, plus `t`'s three distinct
    //   classes, is **four distinct planes through one point** — an identity the plane-triple
    //   substrate cannot express, whether or not `fc` is among them.
    // - `wc ∈ t` but `fc ∉ t` — the face's own plane does not name its own vertex. That is a
    //   naming anomaly, not a concurrency, and keeps the generic `RunName`.
    //
    // Counting to the end rather than returning early matters: a four-plane vertex has *two*
    // off-plane classes, so an early "two off-planes" bail would exit before `wc`'s absence is
    // ever noticed and report every such point as `RunName`.
    let third_on_l =
        |n: NodeId, out: &mut Trace| -> Result<(NodeId, combinatorics::EndPin), DeclineKind> {
            // **A corner a cylinder made keeps its own name**, and `wc` in that name means
            // there is no fourth-plane alias to record. What pins it on `L` is a function of the
            // name and the line ([`combinatorics::pin_for`]): the quadric where its pair *is*
            // `{wc, fc}`, and otherwise its other plane — a corner on three planes and the
            // cylinder is named by the pair its definition chose, which need not be this line's.
            let (n, t) = match three_plane_name(n) {
                Some(t) => (n, t),
                None => match combinatorics::pierce_name(n) {
                    Some((planes, _, _)) if planes.contains(&wc) => {
                        return combinatorics::pin_for(jd, wc, fc, n)
                            .map(|pin| (n, pin))
                            .ok_or(DeclineKind::NoPinOnLine);
                    }
                    // **On `wc` by coincidence**: a wall through a fillet's axis runs
                    // through the fillet's tangent corner, so the corner is on `L` while neither of
                    // its own planes is `wc`. Three planes pass through the point — its own two and
                    // `wc` — so it is a **three-plane point**, named by them like any other; the
                    // cylinder is a fourth carrier that names nothing here. (Restating it as a
                    // pierce of a new pair kept a root that pair does not have — measured, the
                    // fold's gusset beside the plate's fillet.)
                    Some((planes, _, _)) => {
                        // This corner has a third plane through it, so every name
                        // that set yields is one point — its own, this three-plane one, and the
                        // class's crossing of the ruling it sits on. The **table** holds that
                        // identity, and it learns it from the operands before any class is
                        // traced (`seed_from_operands` → `Aliases::record_on_cylinder`), not
                        // here: a round that declines returns, so a discovery made mid-round
                        // could arrive after the class that needed it had already refused.
                        let mut t = [planes[0], planes[1], wc];
                        t.sort_unstable();
                        (NodeId::three_planes(Canon3::three(t)), t)
                    }
                    None => return Err(DeclineKind::RunName),
                },
            };
            // The face's own plane need not appear in the name — a ring vertex is on
            // its face by topology, and a concurrency's canonical triple may name it by three
            // *other* planes. One question remains: does `wc` name it? If so the pin is whichever
            // plane of `t` cuts `L` (`pin_on_line`; with `t = {fc, wc, r}` that is `r`). If not,
            // `wc` is a further plane through the
            // point: record the concurrency, pin the same way, and call the point by the
            // canonical triple of what is now known through it — the rule the operand's ring
            // and the alias table use, so the three agree on the representative.
            if t.contains(&wc) {
                combinatorics::pin_on_line(jd, wc, fc, &t)
                    .map(|r| (n, combinatorics::EndPin::Class(r)))
                    .ok_or(DeclineKind::NoPinOnLine)
            } else {
                let mut set = t.to_vec();
                set.push(wc);
                set.sort_unstable();
                set.dedup();
                out.aliases.record(jd, &set);
                // A handle has a duty the identity does not: it must **cut** `L` — see
                // [`combinatorics::pin_on_line`], where that rule and its reason live.
                // `FourPlane` and not `NoPinOnLine`: this arm has already established that
                // `wc` is a fourth plane through the point, so the substrate limit is the cause
                // and a missing pin is its symptom.
                combinatorics::pin_on_line(jd, wc, fc, &t)
                    .ok_or(DeclineKind::FourPlane)
                    .map(|r| {
                        let name = combinatorics::canonical_triple(jd, &set)
                            .map(NodeId::three_planes)
                            .unwrap_or(n);
                        (name, combinatorics::EndPin::Class(r))
                    })
            }
        };

    // Phase A — scan the outer ring and every hole ring into ONE station list. A hole ring keeps
    // its stored CW winding, but the scan is winding-agnostic: the hole is carved because its two
    // crossings of `L` toggle the face back off between them. A ring that does not meet `W` (every
    // vertex one side) leaves no station.
    let mut declined: Option<DeclineKind> = None;
    'rings: for (ring, walls) in outer.iter().chain(holes.iter()) {
        let n = ring.len();
        // Where this ring meets `L`, and whether it crosses or only touches — the reading every
        // road that walks a plane face's ring along a line shares ([`combinatorics::ClassLine`]).
        // What is done with a feature is this function's own business: naming it, recording a
        // four-plane alias, deciding occupancy.
        // **What "on the meet" means here**: `W` cuts this *planar* face in a straight **line**.
        // (The lateral road's meet is a **circle**, where an arc *can* lie on it, so it tests the
        // carrier instead — see `cycle_on_class`.)
        let line = combinatorics::ClassLine::new(jd, cyls, fc, wc, ring, |i| walls[i].into());
        let features = match line.walk() {
            combinatorics::RingWalk::Met(f) => f,
            // Every vertex on `W`: a ring lying in the cut plane is degenerate here.
            combinatorics::RingWalk::AllOn => {
                declined = Some(DeclineKind::AllOnPlane);
                break;
            }
            // Not "a corner a cylinder made" — the walk reads those. This is the
            // walk's own `None`: a side it could not form exactly, an arc between two on-line
            // nodes tangent to the class at its end, or an arc whose crossings could not be
            // counted (`arc_crossings`' `None`, or a count its ends make impossible). Measured
            // 0 raises across the suite; the name is kept because the walk can still say it.
            combinatorics::RingWalk::Unnameable => {
                declined = Some(DeclineKind::PierceNode);
                break;
            }
        };
        for feature in features {
            match feature {
                combinatorics::Feature::Crossing { edge, nth, .. } => {
                    // The crossed edge's wall, **carried** from the producer — not re-derived
                    // from the two endpoint names, which is sound only while every vertex lies on
                    // exactly three planes and can hand back a plane the edge does not ride at a
                    // concurrency.
                    //
                    //
                    // A crossing on a **ruling** is a point on the cylinder — a pierce node
                    // `wc ∩ fc ∩ cyl`, pinned by the quadric ([`combinatorics::crossing_on_ruling`]); the second
                    // operation on a wall boss makes one wherever a ⊥ cap crosses the plate's
                    // wall face along the boss's rulings (the crossing census's mid slab). A
                    // crossing on an **arc** is the same pierce node, at the root the arc strictly
                    // contains — the `nth` in travel order where it contains both
                    // ([`combinatorics::arc_crossings`]).
                    // The gate records the pairs it could not prove clear, so a pair it proved
                    // clear has no root inside any arc — a crossing named here on an unrecorded pair
                    // would meet no ruling on the lateral's side.
                    if let crate::combinatorics::Wall::Arc { cyl, .. } = walls[edge]
                        && !cyls.get(cyl).is_some_and(|cy| cy.draws_rulings(wc))
                    {
                        declined = Some(DeclineKind::CurvedRingWall);
                        break 'rings;
                    }
                    let Some(at) = line.crossing(edge, nth) else {
                        // A wall this trace does not state.
                        declined = Some(DeclineKind::CurvedRingWall);
                        break 'rings;
                    };
                    let (id, pin) = match at {
                        combinatorics::PointOn::Class(w) => (
                            NodeId::three_planes(Canon3::three([wc, fc, w])),
                            combinatorics::EndPin::Class(w),
                        ),
                        combinatorics::PointOn::Pierce(id) => (id, combinatorics::EndPin::Cylinder),
                    };
                    #[cfg(test)]
                    if let crate::combinatorics::Wall::Ruling { cyl, side, .. } = walls[edge]
                        && let Some(wcy) = cyls.get(cyl)
                    {
                        crossing_probe::record(jd, &wcy.def, cyl, wc, fc, side, id);
                    }
                    // A crossing is one point strictly inside its edge where the boundary passes to
                    // the other side of `L` — twice for an arc crossed twice.
                    stations.push((LineEvent::Toggle, (id, pin)));
                }
                combinatorics::Feature::Run {
                    first,
                    len: m,
                    flanks_differ,
                    ..
                } => {
                    // A run's end may be a corner a cylinder made (a half-disk cap's chord);
                    // `third_on_l` names that one by itself, pinned by the quadric.
                    let name = |k: usize, out: &mut Trace| third_on_l(ring[(first + k) % n], out);
                    // Every on-line vertex is **named**, whether or not it becomes a station:
                    // naming is where a fourth plane through the point is recorded (`third_on_l`)
                    // and where a point that cannot be pinned is refused.
                    let names: Result<Vec<(NodeId, combinatorics::EndPin)>, DeclineKind> =
                        (0..m).map(|k| name(k, out)).collect();
                    let rs = match names {
                        Ok(rs) => rs,
                        Err(kind) => {
                            declined = Some(kind);
                            break 'rings;
                        }
                    };
                    if m == 1 {
                        // A single on-line vertex is a point: with differing flanks the boundary
                        // crosses there, with equal flanks it touches and turns back, which changes
                        // nothing on `L`.
                        if flanks_differ {
                            stations.push((LineEvent::Toggle, rs[0]));
                        }
                    } else {
                        // m >= 2: one on-line interval, an *edge* of `f` lying on `L`, so `f` is on
                        // one side of it whatever the ring does afterwards. Only its two ends
                        // become stations; the points between are names the arrangement may split
                        // at. `flanks_differ` says only whether the face goes on to the far side
                        // at the far end — not which side it occupies along the run, and gating
                        // the side on it left a crossing run classified as a straddling
                        // transversal.
                        let Some((first_is_lo, body_above)) =
                            run_body_above(jd, cyls, faces, wc, fc, fp, &rs)
                        else {
                            declined = Some(DeclineKind::PierceNode);
                            break 'rings;
                        };
                        // In ring order, as every other station is pushed: the sort below runs
                        // fewer comparisons on a list already in runs (lower end first costs the
                        // rotated fold 10% more `order_pinned` calls, measured).
                        let start = LineEvent::GrazeStart { body_above };
                        let end = LineEvent::GrazeEnd {
                            toggle: flanks_differ,
                        };
                        let (first, last) = if first_is_lo {
                            (start, end)
                        } else {
                            (end, start)
                        };
                        stations.push((first, rs[0]));
                        stations.push((last, rs[m - 1]));
                    }
                }
            }
            if declined.is_some() {
                break 'rings;
            }
        }
    }
    if let Some(reason) = declined {
        out.declined.push((fp, reason));
        return;
    }

    // Phase B — order the stations along `L`, and at one point by the order their events apply.
    // An order that cannot be formed declines the face — the sort stops at it
    // ([`combinatorics::try_sort_by`]) rather than reading it as a tie. What can fail has
    // narrowed: a cylinder-pinned node is *ordered* here rather than refused, and the only `None`
    // left is a description that could not be formed exactly. The neighbours compared after the
    // sort can fail too — the sort may have asked that pair the other way round.
    // Still unexercised in the corpus — measured through `reject_census`, which does not move.
    let order = |a: (NodeId, combinatorics::EndPin), b: (NodeId, combinatorics::EndPin)| {
        combinatorics::order_pinned(jd, cyls, wc, fc, a, b).map(|o| o.cmp(&0))
    };
    let same: Option<Vec<bool>> = combinatorics::try_sort_by(&mut stations, |a, b| {
        order(a.1, b.1).map(|o| o.then(a.0.cmp(&b.0))).ok_or(())
    })
    .ok()
    .and_then(|()| {
        stations
            .windows(2)
            .map(|w| order(w[0].1, w[1].1).map(|o| o == std::cmp::Ordering::Equal))
            .collect()
    });
    let Some(same) = same else {
        out.declined.push((fp, DeclineKind::PierceNode));
        return;
    };
    // Phase C — the sweep both line roads share.
    let mat = SegKind::Transversal {
        mat: faces[fp].plane().orient_sign,
    };
    match sweep_line(&stations, &same, |p| p.0, mat) {
        Ok(pieces) => out
            .segs
            .extend(pieces.into_iter().map(|([a, b], kind)| Seg {
                wall: plane_ix[fp].plane(),
                end: [a.0, b.0],
                end_h: [a.1, b.1],
                solid: which,
                kind,
            })),
        Err(kind) => out.declined.push((fp, kind)),
    }
}

/// Which side of `W` face `fp` occupies along an on-line run — **the left of the ring's travel**,
/// stated in the label frame (`true` = above) — and whether the run's first name is its lower end
/// along `L` (`t < 0` below: the one comparison answers both).
///
/// A run is an *edge* of `f`'s boundary lying on `L = W ∩ fp`, so along it `f` is on exactly one
/// side, and which side is the universal boundary convention: **material is on the left of the
/// ring's direction of travel**. That holds for an outer ring, a hole ring, a notch, and a reflex
/// corner alike — it does not care how the ring continues past the run.
///
/// Not off the run's *flank* (the neighbouring off-line vertex): that is only a proxy for "which
/// way the ring bulges", and it is **inverted wherever the
/// ring turns away from material**: on a hole ring the neighbours point into the hole, and on an
/// outer-ring notch they point across the notch — both void. Measured: over the whole corpus the
/// two agree on 251 of 252 runs, and the one disagreement (a fused boss's `z=1` annulus, where the
/// run is the hole's edge) is the one the flank gets wrong.
///
/// **Derivation.** `order_along` runs along `d = n_s(wc) × n_s(fc)` (stored normals — see
/// [`combinatorics::dir_sign`] and [`Judge::plane_pair_dir_sign`](nacre_judge::predicate::Judge)), and the label frame's "above" is
/// `n_s(wc)` (see [`combinatorics::side_of`]). With `t = order_along(wc, fc, first, last)` the travel
/// is `−t·d`, so material points along `n_out(fp) × (−t·d)`. Since `fp` is coplanar with its class
/// root, `n_out(fp) = σ·n_s(fc)`, and the triple product collapses to
/// `−t·σ·(1 − (n_s(fc)·n_s(wc))²)` — whose bracket is positive because the two planes are not
/// parallel. Hence `body_above = (t·σ < 0)`. Note `orient_sign(wc)` **cancels**: the ordering
/// direction and the label frame are defined by the same stored normal.
///
/// `σ` is an f64 dot of two **parallel** unit vectors (`fp` and `fc` are the same plane class), so
/// `|σ| ≈ 1` — a full unit from the sign boundary, the same robustness [`FaceInfo::orient_sign`]
/// relies on. Everything else here is exact.
#[allow(clippy::too_many_arguments)]
fn run_body_above(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    wc: usize,
    fc: usize,
    fp: usize,
    rs: &[(NodeId, combinatorics::EndPin)],
) -> Option<(bool, bool)> {
    let planes = jd.planes;
    // The run's own direction along `L`, through the one rule both roads share — a run whose
    // ends a cylinder pinned orders by the `a + b√c` tower, and by the integer predicates
    // otherwise. `None` is a missing description, which the caller turns into its own decline.
    let t = combinatorics::order_pinned(jd, cyls, wc, fc, rs[0], rs[rs.len() - 1])?;
    let sigma = planes[fc].frame_sign * crate::planes::face_facing_noted(jd, faces, fp, fc)?;
    Some((t < 0, (t < 0) == (sigma > 0)))
}

/// Trace one solid on plane class `wc` (a canon root, i.e. an index into `planes`). This brick:
/// seated faces → their boundary as segments; every other face → `declined`.
///
/// `faces_in` is that solid's `(slot, rings)` pairs — see [`combinatorics::TraceInput`]. The rings
/// are derived once, not per class here: nothing about them depends on `wc`.
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_one(
    faces_in: &[(usize, combinatorics::FaceLoops)],
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    aliases: &Aliases,
    out: &mut Trace,
) {
    let planes = jd.planes;
    // Every curved carrier a seated ring rode and this walk declined to re-emit, with the face
    // that rode it. Checked once at the end — a lateral face may be visited after the seated one,
    // so the question is only answerable when the solid's whole contribution is in.
    let mut curves_owed: Vec<(usize, crate::combinatorics::Wall)> = Vec::new();
    for (fp, fl) in faces_in {
        let fp = *fp;
        // A **lateral face** (a cylinder row): what it leaves on a ⊥ class is a **transversal**
        // circle where the class cuts it in two and a **graze** where the class carries one of its
        // rims — see [`circle_on_class`]. Beyond the span it leaves nothing. Non-⊥ classes cannot
        // carry a circle, and the population gate has already named every such interaction.
        if let ClassIx::Cyl(k) = plane_ix[fp] {
            let cf = match &faces[fp] {
                FaceRow::Cylinder(cf) => cf,
                FaceRow::Plane(_) => unreachable!("ClassIx::Cyl marks a cylinder row"),
            };
            match circle_on_class(jd, cyls, cf, fl, wc, k) {
                Ok(spans) if !spans.is_empty() => {
                    out.circles.extend(spans.into_iter().map(|s| CircleTrace {
                        cyl: k,
                        solid: which,
                        kind: match s.on {
                            CylOnClass::Crosses => SegKind::Transversal {
                                mat: cf.orient_sign,
                            },
                            CylOnClass::Grazes { body_above } => SegKind::Graze { body_above },
                        },
                        arc: s.arc,
                    }));
                }
                // No circle: a wall class **within the radius** leaves two rulings instead (the
                // rulings road, through the axis or offset from it); any other non-⊥ class still
                // leaves nothing, silently — the population gate names those interactions.
                Ok(_) => match rulings_on_class(jd, cf, fl, wc, k, which, cyls, aliases) {
                    Ok((v, grazes)) => {
                        out.rulings.extend(v);
                        // The lateral's side of a plane-pair line it touches.
                        out.segs.extend(grazes);
                    }
                    Err(kind) => out.declined.push((fp, kind)),
                },
                Err(kind) => out.declined.push((fp, kind)),
            }
            continue;
        }
        if plane_ix[fp].plane() != wc {
            trace_transversal_face(fp, fl, which, wc, jd, cyls, faces, plane_ix, out);
            continue;
        }
        // Seated: the face lies in W, so its whole boundary is trace. The body lies on one
        // side of W — the face's outward points away from the body, so the body is above W
        // exactly when that outward opposes the way W faces: `frame_sign · facing < 0`
        // ([`crate::planes::face_facing`]). An undecided facing reads as `false` here, as an
        // undecided sign reads as `0` everywhere — the judgement's evidence is recorded, and the
        // operation refuses by its name before this trace is used.
        let body_above = crate::planes::face_facing_noted(jd, faces, fp, wc)
            .is_some_and(|f| planes[wc].frame_sign * f < 0);
        let kind = SegKind::Seated { body_above };
        // Collect every ring in class form first: a collapsed name declines the whole face, and
        // deciding that before the emitting closure exists keeps the two borrows apart.
        // A **circular** boundary — a disk face's outer, or a circular drill hole — is a
        // closed trace element ([`CircleTrace`]) rather than a segment ring.
        // The **walls travel with the triples**. `plane_ring` hands both back, and this caller
        // keeps both.
        let mut rings: Vec<(Vec<NodeId>, Vec<crate::combinatorics::Wall>)> = Vec::new();
        match &fl.outer {
            Some(combinatorics::LoopRing::Poly(nr)) => match plane_ring(nr) {
                Ok(r) => rings.push(r),
                // `OuterRing` here, not `CollapsedTriple`: this caller reports *which loop*
                // failed, which is the finer fact when a face has several.
                Err(RingFail::Collapsed) => {
                    out.declined.push((fp, DeclineKind::OuterRing));
                    continue;
                }
            },
            Some(combinatorics::LoopRing::Circle { cyl }) => {
                out.circles.push(CircleTrace {
                    cyl: *cyl,
                    solid: which,
                    kind,
                    // A seated circle is a whole loop of the face lying in the class: the disk
                    // cap's own rim, or a circular hole through it. It has no partial extent.
                    arc: None,
                });
            }
            // A lateral's rim never names a plane face's outer loop (see `FaceLoops::cycles`).
            Some(combinatorics::LoopRing::Rim { .. }) | None => {
                out.declined.push((fp, DeclineKind::OuterRing));
                continue;
            }
        }
        let mut failed: Option<DeclineKind> = None;
        // As above: an unnameable hole is a reject, not "no hole".
        let Some(raw) = &fl.holes else {
            out.declined.push((fp, DeclineKind::HoleRing));
            continue;
        };
        for r in raw {
            match r {
                combinatorics::LoopRing::Circle { cyl } => out.circles.push(CircleTrace {
                    cyl: *cyl,
                    solid: which,
                    kind,
                    arc: None,
                }),
                combinatorics::LoopRing::Rim { .. } => failed = Some(DeclineKind::HoleRing),
                combinatorics::LoopRing::Poly(nr) => match plane_ring(nr) {
                    Ok(r) => rings.push(r),
                    Err(f) => failed = Some(decline_of(f)),
                },
            }
        }
        if let Some(kind) = failed {
            out.declined.push((fp, kind));
            continue;
        }
        let fc = plane_ix[fp].plane();
        let mut emit_ring = |ns: &[combinatorics::NodeId], ws: &[crate::combinatorics::Wall]| {
            // **The carrier is taken, not derived.** Re-reading it out of the two endpoint
            // *names* — "the class they share besides `fc`" — is sound only while every vertex
            // lies on exactly three planes, and at a four-plane concurrency can hand back a plane
            // the edge does not ride. `NamedRing.walls` is what the producer read off the **edge's
            // own two faces**, total even where the vertex names fall back. Measured across the
            // suite: the two agree 810,925 times and differ 0, and the carried one also answers 40
            // edges the derivation cannot name at all.
            //
            // The corner keeps its own name too: `NodeId::three_planes(t)` would be a round trip
            // through a projection, and `three_plane_name` is pure extraction of an already-sorted
            // triple, so the two are the same value.
            let n = ns.len();
            for i in 0..n {
                // **A curved edge emits nothing, and that is a refusal to say a thing
                // twice rather than a hole.** A [`Seg`] is a straight edge on `wc`; an arc or a
                // ruling is not that species. The element it needs is already on this class,
                // contributed by the **cylinder's own lateral face** of this same solid — the
                // circle where a ⊥ class cuts it, the two rulings where a class carries its axis.
                // Measured on every class where this fires: `circles = 1` on the cap classes and
                // `rulings = 2` on the wall class, all present before the seated face is reached.
                // The claim is checked at the end of this walk rather than trusted — see
                // `curves_owed` — because an element silently missing does not decline, it
                // corrupts every label on the class.
                // Except a **tangent** ruling (`side == 0`): the lateral touches this
                // class along it and crosses nowhere, so no ruling trace will ever back it — the
                // seated face is the only face on this class with that line as an edge, and it
                // states the piece itself, `Seated` like its straight edges. `up` says which end
                // is the lower one along the axis (a ruling [`Piece`]'s `end` ascends).
                if let crate::combinatorics::Wall::Ruling {
                    cyl, side: 0, up, ..
                } = ws[i]
                {
                    let (a, b) = (ns[i], ns[(i + 1) % n]);
                    // The axis side, in the label frame (`W`'s stored normal): the one bit a
                    // tangent edge flips ([`SegKind::Tangent`]).
                    let axis_above = (|| -> Option<bool> {
                        let w = jd.stored_world_rat(wc)?;
                        let o = cyls.get(cyl)?.def.origin();
                        let mut acc = w[3];
                        for k in 0..3 {
                            acc = acc.checked_add(w[k].checked_mul(o[k])?)?;
                        }
                        Some(acc > Rat::from_int(0))
                    })();
                    let Some(axis_above) = axis_above else {
                        out.declined.push((fp, DeclineKind::Ruling));
                        continue;
                    };
                    out.rulings.push(RulingTrace {
                        cyl,
                        side: 0,
                        end: if up { [a, b] } else { [b, a] },
                        solid: which,
                        kind: SegKind::Tangent { axis_above },
                        #[cfg(test)]
                        orient: faces[fp].plane().orient_sign,
                    });
                    curves_owed.push((fp, ws[i]));
                    continue;
                }
                let crate::combinatorics::Wall::Plane(wall) = ws[i] else {
                    curves_owed.push((fp, ws[i]));
                    continue;
                };
                // **This test comes first, and the order is what keeps two causes apart.**
                // A carrier that *is* this face's own class makes `fc ∩ wall` not a line, and
                // `plane_pair_dir_sign(fc, fc, ·)` is a determinant with two equal rows — zero
                // against every candidate. Ask the pin first and the pencil case arrives wearing
                // `NoPinOnLine`, which is the wrong sentence for it.
                // Measured 0 in the corpus. A derivation filtering `c != fc` cannot produce it, so
                // the guard is what carries that property. Measured: made `unreachable!()`, this
                // arm and the pin's below left the workspace suite and the ignored sweep green —
                // that is what "unexercised" means here, rather than an argument that they cannot
                // fire. (The curved arm above is *not* in that set: it fires, 20 times, and
                // `curves_owed` is what checks it.)
                if wall == fc {
                    out.declined.push((fp, DeclineKind::SeatedEdgeNaming));
                    continue;
                }
                // **What pins each end is a function of its name and this line**
                // ([`combinatorics::pin_for`] — the rule `draft::Ring::edges` and the scan road
                // read). A three-plane corner is pinned by a plane ([`combinatorics::pin_on_line`]).
                // A corner a cylinder made is pinned by the quadric where its pair is `{fc, wall}`
                // ([`combinatorics::EndPin::Cylinder`], the name beside it saying which root) —
                // and by its other plane where it is a corner on three planes whose name's pair is
                // another: read as the quadric there, the name's root would be taken for a root of
                // this line, the circle's far point.
                let pin = |k: usize| combinatorics::pin_for(jd, fc, wall, ns[k]);
                let (Some(h0), Some(h1)) = (pin(i), pin((i + 1) % n)) else {
                    out.declined.push((fp, DeclineKind::NoPinOnLine));
                    continue;
                };
                out.segs.push(Seg {
                    wall,
                    end: [ns[i], ns[(i + 1) % n]],
                    end_h: [h0, h1],
                    solid: which,
                    kind,
                });
            }
        };
        for (ring, walls) in &rings {
            emit_ring(ring, walls);
        }
    }
    // **The seated walk's silence over a curved edge, checked rather than argued.** It emits
    // no segment there because the cylinder's own lateral face already put the element on this
    // class. That follows from the population gate — a class meeting a cylinder is ⊥ to its axis
    // (a circle), carries its axis (two rulings), or was refused long before here — but an
    // argument is not a measurement, and an element that goes missing does not raise: it silently
    // relabels the class.
    //
    // Asked as narrowly as the carrier speaks. A `(cylinder, class)` pair has **two** rulings and
    // a ring rides one of them, so the side is part of the question; a plane meets a cylinder in
    // one circle, so there the cylinder alone names it.
    for (fp, w) in curves_owed {
        let backed = match w {
            crate::combinatorics::Wall::Ruling { cyl, side, .. } => out
                .rulings
                .iter()
                .any(|r| r.cyl == cyl && r.side == side && r.solid == which),
            crate::combinatorics::Wall::Arc { cyl, .. } => {
                out.circles.iter().any(|c| c.cyl == cyl && c.solid == which)
            }
            crate::combinatorics::Wall::Plane(_) => {
                unreachable!("only curved carriers are collected")
            }
        };
        if !backed {
            out.declined.push((fp, DeclineKind::SeatedCurveUnbacked));
        }
    }
}

#[cfg(test)]
#[path = "../tests/trace_plane.rs"]
mod tests;
