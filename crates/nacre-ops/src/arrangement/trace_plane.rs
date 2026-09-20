use super::*;
/// Why a ring of names could not be carried.
///
/// ★★★★★ **It had three variants and now has one, and that is the shape of what changed.** The
/// other two — a corner with no three-plane name, a carrier riding a cylinder — were never
/// *decisions*: the ring came out as plane ids and a curved ring could not be **described** in it,
/// so the type ran out and the refusal was dressed up as a cause. Two rungs widened the vessel and
/// both refusals went with the projections they were standing in for.
///
/// What is left is a name that is degenerate **as a name**, which is a fact about the triple and
/// has nothing to do with cylinders — so it survives, alone, and says only itself.
#[derive(Debug)]
pub(super) enum RingFail {
    /// Two of a vertex's three planes coincide, so the name denotes no point.
    Collapsed,
}

/// The decline kind for a [`RingFail`] at a caller that has no finer fact to add.
///
/// ★ Two of the four callers *do* have one — they know which loop of the face failed and report
/// `OuterRing`/`HoleRing` — so this is not "the mapping", only the default.
fn decline_of(f: RingFail) -> DeclineKind {
    match f {
        RingFail::Collapsed => DeclineKind::CollapsedTriple,
    }
}

/// ★ **The one place a ring of *names* becomes a ring of *plane data*** for the tracer, which asks
/// only "which side, which wall, which third plane" of each vertex. The sort that used to stand
/// here is gone — [`NodeId`]'s only constructor sorts — so what remains is the collapse check, and
/// that is the whole reason this function exists.
pub(super) fn plane_ring(
    nr: &combinatorics::NamedRing,
) -> Result<(Vec<combinatorics::NodeId>, Vec<crate::boolean::Wall>), RingFail> {
    // ★★★★★ **Both refusals are gone, and the vessel is why they could go.** The corner used to
    // be projected to `[usize; 3]` and the carrier to a plane class, so a ring a cylinder touched
    // could not be *described* — the refusals here were the type running out, dressed as
    // decisions. Two rungs widened the vessel; this one stops asking. What a curved ring now meets
    // is the roads themselves, which answer or decline on their own terms.
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
/// on-line edges (a run of `side == 0` vertices) instead of rejecting them. A **holed** face is
/// handled by scanning its inner rings into the same node list: each hole's two crossings of `L` toggle the parity sweep back to void between them,
/// carving the hole out of the emitted chord.
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
    crossings: &std::collections::HashSet<(usize, usize)>,
    out: &mut Trace,
) {
    // `fp` names a *face* (`n_out`, `orient`, the declined log); `fc` names the *plane class* it
    // lies on (triples, comparisons, predicate arguments). Every ring below is in class form, so
    // the two must not be confused — see `canon_ring`.
    let fc = plane_ix[fp].plane();
    // Each ring pairs its class-form triples with the **carried walls** the producer read off
    // the model's edges (`NamedRing`) — the Crossing arm below names an edge's wall from the
    // ride, never from the two endpoint names.
    // The flip nodes a **circle** of this face leaves on `L` — its chord's two ends — for the
    // outer and the holes alike, joined to the
    // scan's node list below.
    let mut circle_nodes: Vec<Node> = Vec::new();
    // The outer ring the scan walks — `None` for a circular outer, whose two chord ends are
    // already in `circle_nodes` and which has no polygon to walk.
    let outer: Option<(Vec<combinatorics::NodeId>, Vec<crate::boolean::Wall>)> =
        match loops.outer.as_ref().and_then(|r| r.poly()) {
            Some(nr) => match plane_ring(nr) {
                Ok(pair) => Some(pair),
                Err(f) => {
                    out.declined.push((fp, decline_of(f)));
                    return;
                }
            },
            // ★ A **circular outer is a disk — and a disk can have holes.** A recorded wall
            // class cuts the circle in a chord (a diameter when it runs through the axis), and
            // the chord's two ends are flip nodes of the same parity sweep every polygon face
            // runs, so the face's **holes carve the chord** exactly as they carve a polygon's
            // section: a tube's annular cap sectioned by a wall is two pieces, not the whole
            // chord across the bore. ★★ This arm used to emit the whole chord on a vessel of its
            // own (`ChordTrace`) and *return before reading the holes* — "a circular outer is a
            // disk face", true only while the gate refused every solid with two coaxial cylinders
            // (annulus, tube); the day the pair rule read faces, the tube's bore came through as
            // a chord piece the cap does not cover and the class refused `LabelConflict`
            // (measured, the annulus × box rows). One road now; the vessel is gone.
            //
            // Any other class misses or is left silent: another ⊥ class is parallel to the disk's
            // own and meets it nowhere, and a ∥-axis wall the gate did not record either clears
            // the disk (the gate's plane test) or cuts an **irrational** chord this road cannot
            // state yet (the face-cleared family — its wrongly-silent trace lands in cells no
            // face reaches, today's recorded state). A disk the wall clears has holes the wall
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
                match chord_nodes(jd, faces, plane_ix, wc, fc, cyl, crossings) {
                    Ok(Some(pair)) => circle_nodes.extend(pair),
                    Ok(None) | Err(ChordFail::NoCoefficients) => return,
                    Err(ChordFail::Unstatable) => {
                        out.declined.push((fp, DeclineKind::Ruling));
                        return;
                    }
                }
                None
            }
            None => {
                #[cfg(test)]
                decline_probe::mark(decline_probe::Site::NoOuter);
                out.declined.push((fp, DeclineKind::OuterRing));
                return;
            }
        };
    let mut holes: Vec<(Vec<combinatorics::NodeId>, Vec<crate::boolean::Wall>)> = Vec::new();
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
            // ★ A **circular hole meets the line in a chord, and then it flips parity twice**.
            // The population gate admits a ∥-axis wall clear of the hole's cylinder
            // (skip: the circle cannot meet `L`) or **recorded** — within the radius, through
            // the axis or offset from it — and there the line cuts the hole in the
            // chord's two pierce points, which used to be skipped "by proof": the proof covered
            // the clear half only, and the planted full-width chord surfaced as `LabelConflict`
            // on the through-family × through-axis tool (measured). The two roots are pushed as
            // flip nodes — the hole's own chord, named the same way a disk outer's is
            // ([`chord_nodes`], the one spelling for both).
            //
            // ★ **The gate's record, the same one the rulings road reads first.** A pair the
            // gate did not list was proven clear — no emitted face reaches the hole's footprint
            // on this class, and the rulings road is silent for it — so planting flip nodes
            // would break the line against rulings that are (rightly) absent: measured, the
            // d = 0 boss-and-bore fixture walked its spur out and back and refused
            // `StraightAngle`. Listed pairs get the nodes; unlisted pairs keep today's skip, and
            // the two roads stay one rule.
            combinatorics::LoopRing::Circle { cyl } => {
                match chord_nodes(jd, faces, plane_ix, wc, fc, *cyl, crossings) {
                    Ok(Some(pair)) => circle_nodes.extend(pair),
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
            // ★★ **A corner a cylinder made needs no third plane, and has none.** It is already the
            // point's own name, and what pins it on `L` is the quadric — measured, every such corner
            // that lands on a cut line is pinned by that very cut plane (`far == wc`), so `wc` is in
            // its name and there is no fourth-plane alias to record either.
            let (n, t) = match three_plane_name(n) {
                Some(t) => (n, t),
                None => match combinatorics::pierce_name(n) {
                    // Pinned by the quadric: its own pair is this line's.
                    Some((planes, _, _)) if planes.contains(&wc) => {
                        return Ok((n, combinatorics::EndPin::Cylinder));
                    }
                    // ★ **On `wc` by coincidence**: a wall through a fillet's axis runs
                    // through the fillet's tangent corner, so the corner is on `L` while neither of
                    // its own planes is `wc`. Three planes pass through the point — its own two and
                    // `wc` — so it is a **three-plane point**, named by them like any other; the
                    // cylinder is a fourth carrier that names nothing here. (Restating it as a
                    // pierce of a new pair kept a root that pair does not have — measured, the
                    // fold's gusset beside the plate's fillet.)
                    Some((planes, _, _)) => {
                        // ★ This corner has a third plane through it, so every name
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
            // ★ The face's own plane need not appear in the name — a ring vertex is on
            // its face by topology, and a concurrency's canonical triple may name it by three
            // *other* planes. One question remains: does `wc` name it? If so the pin is whichever
            // plane of `t` cuts `L` (`pin_on_line`; with `t = {fc, wc, r}` that is `r`, the old
            // ordinary arm, measured identical). If not, `wc` is a further plane through the
            // point: record the concurrency, pin the same way, and call the point by the
            // canonical triple of what is now known through it — the rule the operand's ring
            // and the alias table use, so the three agree on the representative.
            if t.contains(&wc) {
                combinatorics::pin_on_line(jd, wc, fc, t)
                    .map(|r| (n, combinatorics::EndPin::Class(r)))
                    .ok_or(DeclineKind::NoPinOnLine)
            } else {
                let mut set = t.to_vec();
                set.push(wc);
                set.sort_unstable();
                set.dedup();
                out.aliases.record(jd, &set);
                // ★ A handle has a duty the identity does not: it must **cut** `L` — see
                // [`combinatorics::pin_on_line`], where that rule and its reason live.
                // ★ `FourPlane` and not `NoPinOnLine`: this arm has already established that
                // `wc` is a fourth plane through the point, so the substrate limit is the cause
                // and a missing pin is its symptom.
                combinatorics::pin_on_line(jd, wc, fc, t)
                    .ok_or(DeclineKind::FourPlane)
                    .map(|r| {
                        let name = combinatorics::canonical_triple(jd, &set)
                            .map(NodeId::three_planes)
                            .unwrap_or(n);
                        (name, combinatorics::EndPin::Class(r))
                    })
            }
        };

    // Phase A — scan the outer ring and every hole ring, collecting feature nodes into ONE list.
    // A hole ring keeps its stored CW winding, but the Phase-A predicates are winding-agnostic; the
    // parity sweep (Phase C) then carves the hole because its two crossings of `L` toggle parity
    // back to void between them. A ring that does not meet `W` (every vertex one side) yields no
    // nodes; `run_counter` is shared so run ids stay unique across rings.
    let mut nodes: Vec<Node> = Vec::new();
    nodes.append(&mut circle_nodes);
    let mut run_counter = 0usize;
    let mut declined: Option<DeclineKind> = None;
    'rings: for (ring, walls) in outer.iter().chain(holes.iter()) {
        let n = ring.len();
        // Where this ring meets `L`, and whether it crosses or only touches — the walk the ray
        // caster shares (`combinatorics::ring_against_plane`). What is done with a feature is this
        // function's own business: naming it, recording a four-plane alias, deciding occupancy.
        // ★ **What "on the meet" means here**: `W` cuts this *planar* face in a straight **line**,
        // and two points fix a line — so a straight edge between two on-line nodes is on it and a
        // curved one is not. (The lateral road's meet is a **circle**, where an arc *can* lie on it,
        // so it tests the carrier instead — see `cycle_on_class`.)
        // An arc between two on-line nodes leaves the line to the side its tangent points
        // ([`combinatorics::arc_departure_side`]); a straight edge stays on it.
        let on_meet = |i: usize| match walls[i] {
            crate::boolean::Wall::Arc { cyl, ccw } => {
                combinatorics::arc_departure_side(jd, cyls, ring[i], wc, cyl, ccw)
                    .map(combinatorics::EdgeMeet::Departs)
            }
            crate::boolean::Wall::Plane(_) | crate::boolean::Wall::Ruling { .. } => {
                Some(combinatorics::EdgeMeet::On)
            }
        };
        let features = match combinatorics::ring_against_plane(jd, cyls, ring, wc, on_meet) {
            combinatorics::RingWalk::Met(f) => f,
            // Every vertex on `W`: a ring lying in the cut plane is degenerate here.
            combinatorics::RingWalk::AllOn => {
                declined = Some(DeclineKind::AllOnPlane);
                break;
            }
            // ★ Not "a corner a cylinder made" any more — the walk reads those. This is the
            // walk's own `None`: a side it could not form exactly, or an arc tangent to the
            // class at its end. ☑ Measured 0 raises across the workspace suite; the name
            // is kept because the walk can still say it.
            combinatorics::RingWalk::Unnameable => {
                declined = Some(DeclineKind::PierceNode);
                break;
            }
        };
        for feature in features {
            match feature {
                combinatorics::Feature::Crossing { edge, .. } => {
                    // The crossed edge's wall, **carried** from the producer — it used to be
                    // re-derived from the two endpoint names (`ring_from_names`), which is
                    // sound only while every vertex lies on exactly three planes and could
                    // hand back a plane the edge does not ride at a concurrency.
                    //
                    //
                    // ★ A crossing on a **ruling** is a point on the cylinder — a pierce node
                    // `wc ∩ fc ∩ cyl`, pinned by the quadric ([`crossing_on_ruling`]); the second
                    // operation on a wall boss makes one wherever a ⊥ cap crosses the plate's
                    // wall face along the boss's rulings (the crossing census's mid slab). A
                    // crossing on an **arc** still declines: the lateral's ruling road cannot
                    // yet cut a ruling at a hole it does not carry on the class,
                    // so naming the arc's crossing here would meet a phantom ruling there.
                    let (id, pin) = match walls[edge] {
                        crate::boolean::Wall::Plane(w) => (
                            NodeId::three_planes(Canon3::three([wc, fc, w])),
                            combinatorics::EndPin::Class(w),
                        ),
                        crate::boolean::Wall::Ruling { cyl, side, .. } => {
                            let Some(wcy) = cyls.get(cyl) else {
                                declined = Some(DeclineKind::CurvedRingWall);
                                break 'rings;
                            };
                            match crossing_on_ruling(jd, &wcy.def, wc, fc, cyl, side) {
                                Ok(id) => {
                                    crossing_probe::record(
                                        jd, &wcy.def, cyl, wc, fc, side, false, id,
                                    );
                                    (id, combinatorics::EndPin::Cylinder)
                                }
                                Err(d) => {
                                    declined = Some(d);
                                    break 'rings;
                                }
                            }
                        }
                        crate::boolean::Wall::Arc { cyl, ccw } => {
                            let Some(wcy) = cyls.get(cyl) else {
                                declined = Some(DeclineKind::CurvedRingWall);
                                break 'rings;
                            };
                            let b = ring[(edge + 1) % ring.len()];
                            match crossing_on_arc(jd, &wcy.def, wc, fc, cyl, ccw, ring[edge], b) {
                                Ok(id) => {
                                    crossing_probe::record(jd, &wcy.def, cyl, wc, fc, 0, true, id);
                                    (id, combinatorics::EndPin::Cylinder)
                                }
                                Err(d) => {
                                    declined = Some(d);
                                    break 'rings;
                                }
                            }
                        }
                    };
                    // A crossing is strict on both ends and its carrier is a line or a ruling —
                    // straight either way — so it is met exactly once: parity flips.
                    nodes.push(Node {
                        id,
                        pin,
                        flip: true,
                        run: None,
                        flanks_differ: false,
                        single_touch: false,
                        graze_above: None,
                    });
                }
                combinatorics::Feature::Run {
                    first,
                    len: m,
                    flanks_differ,
                    ..
                } => {
                    // ★ A run's end may be a corner a cylinder made (a half-disk cap's chord);
                    // `third_on_l` names that one by itself, pinned by the quadric.
                    let name = |k: usize, out: &mut Trace| third_on_l(ring[(first + k) % n], out);
                    if m == 1 {
                        match name(0, out) {
                            Ok((id, pin)) => nodes.push(Node {
                                id,
                                pin,
                                flip: flanks_differ,
                                run: None,
                                flanks_differ,
                                single_touch: !flanks_differ,
                                // A single-vertex tangential touch is a point (`touches`), never a
                                // graze segment; a flank-differing single vertex is a strict
                                // crossing.
                                graze_above: None,
                            }),
                            Err(kind) => declined = Some(kind),
                        }
                    } else {
                        // m >= 2: one on-line interval. Only its two ends and the flanks decide
                        // anything; the interior points are names the arrangement may split at.
                        let id = run_counter;
                        run_counter += 1;
                        let names: Result<Vec<(NodeId, combinatorics::EndPin)>, DeclineKind> =
                            (0..m).map(|k| name(k, out)).collect();
                        match names {
                            Ok(rs) => {
                                // Every run is one-sided: it is an *edge* of `f` lying on `L`, so
                                // `f` is on one side of it whatever the ring does afterwards.
                                // `flanks_differ` says only whether the sweep's parity toggles here
                                // (Phase B gives it to `flip`) — it is not an occupancy fact, and
                                // gating the side on it left a crossing run classified as a
                                // straddling transversal.
                                let Some(ba) = run_body_above(jd, cyls, faces, wc, fc, fp, &rs)
                                else {
                                    declined = Some(DeclineKind::PierceNode);
                                    break 'rings;
                                };
                                let graze_above = Some(ba);
                                for (nid, pin) in rs {
                                    nodes.push(Node {
                                        id: nid,
                                        pin,
                                        flip: false,
                                        run: Some(id),
                                        flanks_differ,
                                        single_touch: false,
                                        graze_above,
                                    });
                                }
                            }
                            Err(kind) => declined = Some(kind),
                        }
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

    // Phase B — order the nodes along L and fix run structure.
    // ★ A comparator cannot decline, so it raises a flag the caller reads. What can make it fire
    // has narrowed: the second road landed, so a cylinder-pinned node is *ordered* here rather
    // than refused, and the only `None` left is a description that could not be formed exactly.
    // ☑ Still unexercised in the corpus — measured through `reject_census`, which does not move.
    let mut unordered = false;
    let order = |a: &Node, b: &Node, unordered: &mut bool| -> std::cmp::Ordering {
        let Some(o) = combinatorics::order_pinned(jd, cyls, wc, fc, (a.id, a.pin), (b.id, b.pin))
        else {
            *unordered = true;
            return std::cmp::Ordering::Equal;
        };
        match o {
            -1 => std::cmp::Ordering::Less,
            1 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        }
    };
    nodes.sort_by(|a, b| order(a, b, &mut unordered));
    if unordered {
        out.declined.push((fp, DeclineKind::PierceNode));
        return;
    }
    for w in nodes.windows(2) {
        if order(&w[0], &w[1], &mut unordered) == std::cmp::Ordering::Equal {
            out.declined.push((fp, DeclineKind::CoincidentFeatures));
            return;
        }
    }
    // Each run's two nodes must be adjacent after the sort; the hi node carries the flip.
    for k in 0..nodes.len() {
        if let Some(id) = nodes[k].run {
            let prev = k > 0 && nodes[k - 1].run == Some(id);
            let next = k + 1 < nodes.len() && nodes[k + 1].run == Some(id);
            if !prev && !next {
                out.declined.push((fp, DeclineKind::RunSplit));
                return;
            }
            // The whole run is one interval, so its flip happens once — at the far end.
            if prev && !next {
                nodes[k].flip = nodes[k].flanks_differ;
            }
        }
    }

    // Phase C — sweep left to right; covered gaps merge into segments.
    let mut parity = 0i8;
    // `seg_start` also carries the graze side of the segment being built: `Some(ba)` when the
    // segment is a pure graze gap (opened by a forced run while outside material), else `None`.
    let mut seg_start: Option<((NodeId, combinatorics::EndPin), Option<bool>)> = None;
    let emit = |a: (NodeId, combinatorics::EndPin),
                b: (NodeId, combinatorics::EndPin),
                graze: Option<bool>,
                out: &mut Trace| {
        out.segs.push(Seg {
            wall: plane_ix[fp].plane(),
            end: [a.0, b.0],
            end_h: [a.1, b.1],
            solid: which,
            kind: match graze {
                Some(body_above) => SegKind::Graze { body_above },
                None => SegKind::Transversal {
                    mat: faces[fp].plane().orient_sign,
                },
            },
        });
    };
    for k in 0..nodes.len() {
        if nodes[k].single_touch && parity == 0 {
            out.touches.push(nodes[k].id);
        }
        if nodes[k].flip {
            parity ^= 1;
        }
        let forced =
            k + 1 < nodes.len() && nodes[k].run.is_some() && nodes[k].run == nodes[k + 1].run;
        let covered = parity == 1 || forced;
        // The kind of the gap *after* node `k`, as a fact about that gap alone: riding a run means
        // `L` is on `f`'s boundary there, so `f` occupies one side; any other covered gap has `L`
        // in `f`'s interior, so `f` straddles. `parity` says whether the gap is covered, not which
        // arcs it fills — reading occupancy off it classified a hole's edge as a straddle.
        let gap_graze = forced.then(|| {
            nodes[k]
                .graze_above
                .expect("every run node carries its occupied side")
        });
        match (seg_start, covered) {
            (None, true) => seg_start = Some(((nodes[k].id, nodes[k].pin), gap_graze)),
            // The kind changes here, so close and reopen: a segment must be homogeneous, because
            // `split_at_crossings` subdivides it later and every piece inherits its kind.
            (Some((a, graze)), true) if graze != gap_graze => {
                emit(a, (nodes[k].id, nodes[k].pin), graze, out);
                seg_start = Some(((nodes[k].id, nodes[k].pin), gap_graze));
            }
            (Some((a, graze)), false) => {
                emit(a, (nodes[k].id, nodes[k].pin), graze, out);
                seg_start = None;
            }
            _ => {}
        }
    }
    if seg_start.is_some() || parity != 0 {
        // A line enters and leaves a bounded region equally; an unbalanced sweep is degenerate.
        out.declined.push((fp, DeclineKind::OddParity));
    }
}

/// Which side of `W` face `fp` occupies along an on-line run — **the left of the ring's travel**,
/// stated in the label frame (`true` = above).
///
/// A run is an *edge* of `f`'s boundary lying on `L = W ∩ fp`, so along it `f` is on exactly one
/// side, and which side is the universal boundary convention: **material is on the left of the
/// ring's direction of travel**. That holds for an outer ring, a hole ring, a notch, and a reflex
/// corner alike — it does not care how the ring continues past the run.
///
/// The earlier version read the side off the run's *flank* (the neighbouring off-line vertex)
/// instead. That is only a proxy for "which way the ring bulges", and it is **inverted wherever the
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
/// `|σ| ≈ 1` — a full unit from the sign boundary, the same robustness [`FaceInfo::orient_sign`] and
/// `trace_seated_face` already rely on. Everything else here is exact.
#[allow(clippy::too_many_arguments)]
fn run_body_above(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    wc: usize,
    fc: usize,
    fp: usize,
    rs: &[(NodeId, combinatorics::EndPin)],
) -> Option<bool> {
    let planes = jd.planes;
    // ★ The run's own direction along `L`, through the one rule both roads share — a run whose
    // ends a cylinder pinned orders by the `a + b√c` tower, and by the integer predicates
    // otherwise. `None` is a missing description, which the caller turns into its own decline.
    let t = combinatorics::order_pinned(jd, cyls, wc, fc, rs[0], rs[rs.len() - 1])?;
    let sigma = faces[fp].plane().n_out.dot(planes[fc].plane.normal());
    Some((t < 0) == (sigma > 0.0))
}

/// Trace one solid on plane class `wc` (a canon root, i.e. an index into `planes`). This brick:
/// seated faces → their boundary as segments; every other face → `declined`.
///
/// `faces_in` is that solid's `(slot, rings)` pairs — see [`combinatorics::TraceInput`]. This walk
/// used to read the shells out of the `Model` and derive the rings here, once per class; nothing
/// about either depends on `wc`.
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_one(
    faces_in: &[(usize, combinatorics::FaceLoops)],
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    crossings: &std::collections::HashSet<(usize, usize)>,
    aliases: &Aliases,
    out: &mut Trace,
) {
    let planes = jd.planes;
    let w_normal = planes[wc].plane.normal();
    // ★ Every curved carrier a seated ring rode and this walk declined to re-emit, with the face
    // that rode it. Checked once at the end — a lateral face may be visited after the seated one,
    // so the question is only answerable when the solid's whole contribution is in.
    let mut curves_owed: Vec<(usize, crate::boolean::Wall)> = Vec::new();
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
                Ok(_) => match rulings_on_class(jd, cf, fl, wc, k, which, crossings, aliases) {
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
            trace_transversal_face(fp, fl, which, wc, jd, cyls, faces, plane_ix, crossings, out);
            continue;
        }
        // Seated: the face lies in W, so its whole boundary is trace. The body lies on one
        // side of W — `n_out` points away from the body, so the body is above W exactly when
        // `n_out · n_W < 0`. The f64 sign is robust even rotated: seated means `canon[fp]==wc`,
        // so `n_out ∥ n_W` (both unit) and the dot is ≈ ±1, a full unit from the sign boundary
        // (the rotated-tunnel tests exercise this seated path through the cube's own caps).
        let body_above = faces[fp].plane().n_out.dot(w_normal) < 0.0;
        let kind = SegKind::Seated { body_above };
        // Collect every ring in class form first: a collapsed name declines the whole face, and
        // deciding that before the emitting closure exists keeps the two borrows apart.
        // A **circular** boundary — a disk face's outer, or a circular drill hole — is a
        // closed trace element ([`CircleTrace`]) rather than a segment ring.
        // ★ The **walls travel with the triples**. `plane_ring` has handed both back since the
        // vessel widened, and this caller dropped one of them on the floor.
        let mut rings: Vec<(Vec<NodeId>, Vec<crate::boolean::Wall>)> = Vec::new();
        match &fl.outer {
            Some(combinatorics::LoopRing::Poly(nr)) => match plane_ring(nr) {
                Ok(r) => rings.push(r),
                // ★ `OuterRing` here, not `CollapsedTriple`: this caller reports *which loop*
                // failed, which is the finer fact when a face has several. The two causes that
                // used to keep their own names here were the vessel's, and went with it.
                Err(RingFail::Collapsed) => {
                    #[cfg(test)]
                    decline_probe::mark(decline_probe::Site::PolyCollapsed);
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
                #[cfg(test)]
                decline_probe::mark(decline_probe::Site::RimOrNone);
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
        let mut emit_ring = |ns: &[combinatorics::NodeId], ws: &[crate::boolean::Wall]| {
            // ★★★★★ **The carrier is taken, not derived.** It used to be re-read out of the two
            // endpoint *names* — "the class they share besides `fc`" — which is sound only while
            // every vertex lies on exactly three planes, and at a four-plane concurrency can hand
            // back a plane the edge does not ride.
            // `NamedRing.walls` is what the producer read off the **edge's own two faces**, total
            // even where the vertex names fall back. Measured across the suite before the swap:
            // the two agree 810,925 times and differ 0, and the carried one additionally answers
            // 40 edges the derivation cannot name at all.
            //
            // ★ The corner keeps its own name too: `NodeId::three_planes(t)` was a round trip
            // through a projection, and `three_plane_name` is pure extraction of an already-sorted
            // triple, so the two are the same value.
            let n = ns.len();
            for i in 0..n {
                // ★★★★★ **A curved edge emits nothing, and that is a refusal to say a thing
                // twice rather than a hole.** A [`Seg`] is a straight edge on `wc`; an arc or a
                // ruling is not that species. The element it needs is already on this class,
                // contributed by the **cylinder's own lateral face** of this same solid — the
                // circle where a ⊥ class cuts it, the two rulings where a class carries its axis.
                // Measured on every class where this fires: `circles = 1` on the cap classes and
                // `rulings = 2` on the wall class, all present before the seated face is reached.
                // ★ The claim is checked at the end of this walk rather than trusted — see
                // `curves_owed` — because an element silently missing does not decline, it
                // corrupts every label on the class.
                // ★ Except a **tangent** ruling (`side == 0`): the lateral touches this
                // class along it and crosses nowhere, so no ruling trace will ever back it — the
                // seated face is the only face on this class with that line as an edge, and it
                // states the piece itself, `Seated` like its straight edges. `up` says which end
                // is the lower one along the axis (`MergedRuling::end` ascends).
                if let crate::boolean::Wall::Ruling { cyl, side: 0, up } = ws[i] {
                    let (a, b) = (ns[i], ns[(i + 1) % n]);
                    // The axis side, in the label frame (`W`'s stored normal): the one bit a
                    // tangent edge flips ([`SegKind::Tangent`]).
                    let axis_above = (|| -> Option<bool> {
                        let w = combinatorics::stored_coeffs_rat(jd, wc)?;
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
                let crate::boolean::Wall::Plane(wall) = ws[i] else {
                    curves_owed.push((fp, ws[i]));
                    continue;
                };
                // ★★★★★ **This test comes first, and the order is what keeps two causes apart.**
                // A carrier that *is* this face's own class makes `fc ∩ wall` not a line, and
                // `plane_pair_dir_sign(fc, fc, ·)` is a determinant with two equal rows — zero
                // against every candidate. Ask the pin first and the pencil case arrives wearing
                // `NoPinOnLine`, which is the wrong sentence for it.
                // ☑ Measured 0 in the corpus. The derivation this replaced could not produce it
                // (it filtered `c != fc`), so the guard is what carries that property across.
                // ☑ This arm and the pin's below were both made `unreachable!()` with the workspace
                // suite and the ignored sweep green — that is what "unexercised" means here, rather
                // than an argument that they cannot fire. (The curved arm above is *not* in that
                // set: it fires, 20 times, and `curves_owed` is what checks it.)
                if wall == fc {
                    out.declined.push((fp, DeclineKind::SeatedEdgeNaming));
                    continue;
                }
                // ★★ **What pins each end, in the vocabulary that names it.** A three-plane
                // corner is pinned by a plane — [`combinatorics::pin_on_line`], the same rule the
                // scan road and the ray caster ask. A corner a cylinder made has no third plane
                // and needs none: the quadric pins it, and [`combinatorics::EndPin::Cylinder`]
                // says so while the name beside it says which point.
                let pin = |k: usize| -> Option<combinatorics::EndPin> {
                    match three_plane_name(ns[k]) {
                        Some(t) => combinatorics::pin_on_line(jd, fc, wall, t)
                            .map(combinatorics::EndPin::Class),
                        None => Some(combinatorics::EndPin::Cylinder),
                    }
                };
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
    // ★★★★ **The seated walk's silence over a curved edge, checked rather than argued.** It emits
    // no segment there because the cylinder's own lateral face already put the element on this
    // class. That follows from the population gate — a class meeting a cylinder is ⊥ to its axis
    // (a circle), carries its axis (two rulings), or was refused long before here — but an
    // argument is not a measurement, and an element that goes missing does not raise: it silently
    // relabels the class.
    //
    // ★ Asked as narrowly as the carrier speaks. A `(cylinder, class)` pair has **two** rulings and
    // a ring rides one of them, so the side is part of the question; a plane meets a cylinder in
    // one circle, so there the cylinder alone names it.
    for (fp, w) in curves_owed {
        let backed = match w {
            crate::boolean::Wall::Ruling { cyl, side, .. } => out
                .rulings
                .iter()
                .any(|r| r.cyl == cyl && r.side == side && r.solid == which),
            crate::boolean::Wall::Arc { cyl, .. } => {
                out.circles.iter().any(|c| c.cyl == cyl && c.solid == which)
            }
            crate::boolean::Wall::Plane(_) => unreachable!("only curved carriers are collected"),
        };
        if !backed {
            out.declined.push((fp, DeclineKind::SeatedCurveUnbacked));
        }
    }
}
