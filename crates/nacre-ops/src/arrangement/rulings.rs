use super::*;
/// **What a lateral face leaves on a recorded wall class** — the ∥ sibling of
/// [`circle_on_class`], and like it, an answer **per extent** rather than one for the whole
/// ruling.
///
/// Asked only after the circle road answered "no circle"; **empty** (silently, today's state)
/// unless the gate **recorded** the pair (the wall within the radius — through the axis or
/// offset from it), and declining rather than contributing
/// *partially* when it does but a piece cannot be stated — a half-contributed rectangle would
/// leave the class's 1-skeleton dangling, which is a worse lie than an honest incomplete-trace
/// mark.
///
/// The class cuts **two** rulings, and each is [`ruling_sweep`]'s: the face's boundary cycles
/// against the line `θ = θ_side` of the `(θ, z)` chart, swept along the axis — `Transversal`
/// where the solid straddles the wall, `Graze` where a cycle's own edge lies on the ruling, and
/// nothing where the face is not there (a panel, a chain rim, a band with holes, all by the
/// one rule).
#[allow(clippy::too_many_arguments)]
pub(super) fn rulings_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    k: usize,
    which: SolidSide,
    crossings: &std::collections::HashSet<(usize, usize)>,
    aliases: &Aliases,
) -> Result<(Vec<RulingTrace>, Vec<Seg>), DeclineKind> {
    // ★ **The gate's answer, first** — see [`combinatorics::TraceInput::crossings`]. A pair not
    // listed there was proven clear (or never in question), and contributing its rectangle
    // anyway breaks the straddling family: the boss's own wall class is a d=0 pair, and its
    // coplanar bottom cap's chord lands on a line the plate already traces. Unlisted pairs are
    // silent; the gate lists exactly the pairs that need the rectangle.
    if !crossings.contains(&(wc, k)) {
        return Ok((Vec::new(), Vec::new()));
    }
    let Some(w) = combinatorics::class_coeffs_rat(jd, wc) else {
        return Ok((Vec::new(), Vec::new())); // no exact description: the circle road's decline covers ⊥ classes
    };
    let Some(def) = cf.def.as_ref() else {
        return Err(DeclineKind::Ruling);
    };
    // ★ The record is the crossing statement: a listed pair's plane runs within the
    // radius, so the wall meets the lateral in two rulings; the sweep reads the face's cycles
    // and says where. The premise is the gate's, asserted rather than re-derived — and it survived
    // the tangent wall opening precisely because a tangency is written to `tangencies`
    // and **not** to `crossings`: "within" still means within.
    debug_assert_eq!(
        nacre_exact::point_plane_clearance_rat(&w, &def.origin(), def.r2()),
        nacre_exact::Orient::Negative,
        "a recorded pair's plane runs within the radius"
    );
    // The face's shape from its cycles — the same reading the ⊥ road makes.
    let LateralShape { rims, cycles, .. } = lateral_shape(jd, cf, fl, def)?;
    let mat = SegKind::Transversal {
        mat: cf.orient_sign,
    };
    let mut out = Vec::with_capacity(2);
    let mut grazes = Vec::new();
    for side in [1i8, -1i8] {
        let (pieces, g) = ruling_sweep(
            jd, cf, def, &w, wc, k, side, &rims, &cycles, mat, which, aliases,
        )?;
        grazes.extend(g);
        #[cfg(test)]
        if pieces
            .iter()
            .any(|(_, kind)| matches!(kind, SegKind::Graze { .. }))
        {
            let o = def.origin().map(|x| x.to_f64());
            let at = |n: NodeId| node_axis_param(jd, def, n).map_or(f64::NAN, |t| t.to_f64());
            ruling_probe::CARVED.push(ruling_probe::Carved {
                origin: o,
                span: [at(pieces[0].0[0]), at(pieces[pieces.len() - 1].0[1])],
                kinds: pieces.iter().map(|&(_, k)| k).collect(),
            });
        }
        for (end, kind) in pieces {
            out.push(RulingTrace {
                cyl: k,
                side,
                end,
                solid: which,
                kind,
                #[cfg(test)]
                orient: cf.orient_sign,
            });
        }
    }
    Ok((out, grazes))
}

/// The θ side a cycle's arc lies on, seen from the node where it meets a ruling: an arc arriving
/// counter-clockwise came from smaller θ (`−1`), one departing counter-clockwise goes to larger θ
/// (`+1`) — the flank of an on-line run, read from the arc's own sense (`Wall::Arc`'s `ccw`)
/// rather than from a coordinate, because a global "side" does not exist on a cylinder.
fn arc_flank(ccw: bool, arriving: bool) -> i8 {
    if ccw == arriving { -1 } else { 1 }
}

/// **Is `node` strictly inside the counter-clockwise arc from `lo` to `hi`?** — cyclic order about
/// the axis, by name ([`circular_order`]). `Err` when the order cannot be formed or two of the
/// three are one point wearing two names — the caller's population, not a tie to shrug at.
pub(super) fn theta_between(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    lo: NodeId,
    hi: NodeId,
    node: NodeId,
) -> Result<bool, DeclineKind> {
    let nodes = [lo, hi, node];
    let order = circular_order(jd, k, def, &nodes).map_err(|_| DeclineKind::Ruling)?;
    let pos = |x: usize| {
        order
            .iter()
            .position(|&i| i == x)
            .expect("a permutation of the three")
    };
    let (p_lo, p_hi, p_n) = (pos(0), pos(1), pos(2));
    Ok((p_n + 3 - p_lo) % 3 < (p_hi + 3 - p_lo) % 3)
}

/// What one side's sweep states: the ruling's pieces (its own vocabulary) and, for a plane-pair
/// line the face touches on this ruling, the face's side of that line (the plane vocabulary).
type SweepOut = (Vec<([NodeId; 2], SegKind)>, Vec<Seg>);

/// **One ruling of a through-axis class against the face's boundary cycles, swept along the
/// axis** — the ∥ twin of the circle road's carving, and the polygon-against-a-line rule
/// of `ring_against_plane` stated on the `(θ, z)` chart for the line `θ = θ_side`.
///
/// Every cycle contributes **stations** on the ruling, each with an event:
/// - a whole **rim** toggles the face (it starts or ends there);
/// - an **arc** that strictly contains the ruling's θ ([`theta_between`]) toggles it — the boundary
///   crosses the line there, at the pierce node `{arc's plane, wc}` ([`crossing_on_ruling`], the
///   one door every road names such a point through);
/// - a **run** of the cycle's own edges lying on the ruling is a `Graze` over its extent, the side
///   the face occupies along it `σ·τ·κ·side` (the rulings ladder's derivation, unchanged), and it
///   toggles the face exactly when the arcs at its two ends lie on different θ flanks
///   ([`arc_flank`]) — a collinear piece the boundary passes *through* rather than touches, the
///   `flanks_differ` of an on-line run. A panel's arcs both lie on its own side (no toggle: absent
///   above and below); a chain's arrive from the outer half and leave into the inner (toggle:
///   transversal below, absent above); a band's notch hole has both on the hole's side (no toggle).
///
/// Then the sweep: absent below the first station, a piece per gap between stations — `Graze`
/// inside a run, `Transversal` while the face is present, nothing while it is not — and adjacent
/// equal pieces merged. Declined by name: a run whose neighbours are not arcs, an arc ending on the
/// ruling where no run of its cycle sits (the boundary would cross at a vertex — a shape
/// `loop_triples` does not produce), two names at one station, a run inside a run, a toggle inside
/// a run (the successor of "a graze past a rim"), and a face still present past the last station.
#[allow(clippy::too_many_arguments)]
fn ruling_sweep(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    w: &[Rat; 4],
    wc: usize,
    k: usize,
    side: i8,
    rims: &[(Rat, usize)],
    cycles: &[(bool, combinatorics::LoopRing)],
    mat: SegKind,
    which: SolidSide,
    aliases: &Aliases,
) -> Result<SweepOut, DeclineKind> {
    use crate::combinatorics::Wall;
    // ★ **The lateral's half of a plane-pair line.** An edge of this cycle can lie on
    // this ruling with the face *across* it being another plane `t` — the fillet's tangent edge,
    // when `wc` passes through the axis and so contains the tangent ruling. That line is a
    // plane-pair line (`wc ∩ t`), and the lateral **ends** there, so the plane vocabulary states
    // it — the half of the one rule for a line that is a ruling and a plane-pair line at once
    // (the other half: where the lateral runs on across the line, the line is the other solid's
    // edge and the ruling vocabulary states it — `split_rulings`' fold). The
    // tangent wall's own run already does (`Feature::Run` → `Graze`), from *its* side. But the
    // solid's material near that line lies on **both** sides of `wc` — the tangent wall's below,
    // the fillet face's above — and an edge's mask is the sum of its faces' statements
    // (`edge_mask`: two grazes of one solid, one per side, flip both bits). So the lateral states
    // its side too — as a `Seg` on that line, in the plane vocabulary, so `merge_coincident`
    // joins it with the wall's — and states nothing in the ruling vocabulary (`quiet_nodes`).
    // The side is read exactly as a run's on this ruling is (`body_side`): the same cycle, the
    // same ruling, the same traversal.
    let mut grazes: Vec<Seg> = Vec::new();
    // ★ **Identity is the table's** — a station made here (`wc`'s crossing with a cap)
    // and a corner the cycle carries can be one point under two names, and the plane roads have
    // already told the table so (`Aliases::record_on_cylinder`). Every "same node?" below asks
    // through it; the sweep never decides coincidence on its own.
    let same = |x: NodeId, y: NodeId| aliases.canon_point(x) == aliases.canon_point(y);
    #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    enum Event {
        // Ordered as they apply at one station: a run ends, the face toggles, a run starts.
        GrazeEnd { toggle: bool },
        Toggle,
        GrazeStart { body_above: bool },
    }
    let on_ruling =
        |c: usize| crossing_on_ruling(jd, def, c, wc, k, side).map_err(|_| DeclineKind::Ruling);
    let mut stations: Vec<(Rat, Event, NodeId)> = Vec::new();
    for &(t, c) in rims {
        stations.push((t, Event::Toggle, on_ruling(c)?));
    }
    // The side of `wc` the face's body lies on along an edge of this cycle that runs on the
    // ruling from `start` to `end` — see [`plus_theta_is_above`], the one spelling of the product.
    let body_side = |start: NodeId, end: NodeId| -> Result<(Rat, Rat, i8, bool), DeclineKind> {
        let (Some(ta), Some(tb)) = (
            node_axis_param(jd, def, start),
            node_axis_param(jd, def, end),
        ) else {
            return Err(DeclineKind::Ruling);
        };
        if ta == tb {
            return Err(DeclineKind::CylHoleFeature); // an extent of no length
        }
        let tau: i8 = if tb > ta { 1 } else { -1 };
        let theta_dot_stored: i8 =
            if plus_theta_is_above(jd, wc, side).ok_or(DeclineKind::Ruling)? {
                1
            } else {
                -1
            };
        let body_above =
            i32::from(material_theta_sign(cf.orient_sign, tau)) * i32::from(theta_dot_stored) > 0;
        Ok((ta, tb, tau, body_above))
    };
    for (inner, ring) in cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair.
        let Some(nr) = ring.poly() else {
            return Err(DeclineKind::CylFaceHole);
        };
        let n = nr.triples.len();
        #[cfg(not(test))]
        let _ = inner;
        // Which edges lie on this ruling. Two shapes, one rule — **the line's carrier says who
        // states it**:
        // - a ruling held by `wc` itself (its far face is *seated* on the class): the
        //   line is the cylinder's alone, and this sweep states it as a run (`on`);
        // - a ruling held by another plane `t` whose two ends are this ruling's cap
        //   crossings (by name, through the table): the line is the plane pair `wc ∩ t`, and
        //   `t`'s own trace states it as a graze run on this class — this sweep stays **silent**
        //   there (`quiet`), and only remembers the ends so an arc ending on them is "an arc that
        //   ends on this ruling", not a crossing.
        let mut on = vec![false; n];
        let mut quiet_nodes: Vec<NodeId> = Vec::new();
        for (i, slot) in on.iter_mut().enumerate() {
            let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
            match nr.walls[i] {
                // The two ends are asked again rather than the carried `side` read: that is one
                // end's, and an edge whose ends answer differently is no ruling.
                Wall::Ruling { plane: t, .. } if t == wc => {
                    let (Some(sa), Some(sb)) = (
                        node_ruling_side(jd, def, w, a),
                        node_ruling_side(jd, def, w, b),
                    ) else {
                        return Err(DeclineKind::Ruling);
                    };
                    // An edge whose two ends answer differently is not a ruling at all — it
                    // would have to cross the plane through the axis. Named rather than silently
                    // mis-placed.
                    if sa != sb {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    *slot = sa == side;
                }
                Wall::Ruling { plane: t, .. } => {
                    // Both ends this ruling's crossings with the caps the ends name? Then the edge
                    // lies on the ruling and `wc ∩ t` carries it.
                    let on_ruling_end = |x: NodeId| -> bool {
                        combinatorics::pierce_name(x).is_some_and(|(pl, _, _)| {
                            pl.iter().any(|&c| {
                                c != t
                                    && crossing_on_ruling(jd, def, c, wc, k, side)
                                        .is_ok_and(|id| same(id, x))
                            })
                        })
                    };
                    if on_ruling_end(a) && on_ruling_end(b) {
                        quiet_nodes.push(a);
                        quiet_nodes.push(b);
                        let (_, _, _, body_above) = body_side(a, b)?;
                        let pin = |x: NodeId| {
                            combinatorics::pin_for(jd, wc, t, x).ok_or(DeclineKind::Ruling)
                        };
                        grazes.push(Seg {
                            wall: t,
                            end: [a, b],
                            end_h: [pin(a)?, pin(b)?],
                            solid: which,
                            kind: SegKind::Graze { body_above },
                        });
                    }
                }
                _ => {}
            }
        }
        if n > 0 && on.iter().all(|&x| x) {
            return Err(DeclineKind::Ruling); // a cycle lying wholly on one ruling is no face
        }
        let mut run_nodes: Vec<NodeId> = Vec::new();
        for i0 in 0..n {
            if !on[i0] || on[(i0 + n - 1) % n] {
                continue; // not the first edge of a run
            }
            let mut i1 = i0;
            while on[(i1 + 1) % n] {
                i1 = (i1 + 1) % n;
            }
            let (start, end) = (nr.triples[i0], nr.triples[(i1 + 1) % n]);
            let (prev, next) = ((i0 + n - 1) % n, (i1 + 1) % n);
            let (Wall::Arc { ccw: nu_in, .. }, Wall::Arc { ccw: nu_out, .. }) =
                (nr.walls[prev], nr.walls[next])
            else {
                return Err(DeclineKind::Ruling); // a run's neighbours are arcs
            };
            let (ta, tb, tau, body_above) = body_side(start, end)?;
            let toggle = arc_flank(nu_in, true) != arc_flank(nu_out, false);
            let ((t_lo, n_lo), (t_hi, n_hi)) = if tau > 0 {
                ((ta, start), (tb, end))
            } else {
                ((tb, end), (ta, start))
            };
            #[cfg(test)]
            ruling_probe::GRAZE_SIDE.push(ruling_probe::GrazeSide {
                inner: *inner,
                body_above,
                ny: jd.planes[wc].plane.normal().as_array()[1],
                origin: def.origin().map(|x| x.to_f64()),
                span: [t_lo.to_f64(), t_hi.to_f64()],
            });
            stations.push((t_lo, Event::GrazeStart { body_above }, n_lo));
            stations.push((t_hi, Event::GrazeEnd { toggle }, n_hi));
            run_nodes.push(start);
            run_nodes.push(end);
        }
        // Arcs strictly containing the ruling's θ: crossings of the line.
        for i in 0..n {
            let Wall::Arc {
                ccw: nu, plane: c, ..
            } = nr.walls[i]
            else {
                continue;
            };
            let node = on_ruling(c)?;
            let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
            if same(node, a) || same(node, b) {
                // The arc ends on this ruling: an edge of this cycle must lie there — a run this
                // sweep states, or a plane-pair line another face states (`quiet_nodes`) — or the
                // boundary crosses the ruling at a vertex: not a shape the ring producer makes,
                // and not one to answer "no" to silently.
                if !run_nodes
                    .iter()
                    .chain(quiet_nodes.iter())
                    .any(|&r| same(r, node))
                {
                    return Err(DeclineKind::Ruling);
                }
                continue;
            }
            let (lo, hi) = if nu { (a, b) } else { (b, a) };
            if theta_between(jd, k, def, lo, hi, node)? {
                let t = crate::bands::param_opt(jd, c, def).ok_or(DeclineKind::Ruling)?;
                stations.push((t, Event::Toggle, node));
            }
        }
    }
    stations.sort_by(|x, y| x.0.cmp(&y.0).then(x.1.cmp(&y.1)));
    for pair in stations.windows(2) {
        if pair[0].0 == pair[1].0 && pair[0].2 != pair[1].2 {
            return Err(DeclineKind::CylHoleFeature); // two names for one station
        }
    }
    let mut pieces: Vec<([NodeId; 2], SegKind)> = Vec::new();
    let (mut present, mut graze): (bool, Option<bool>) = (false, None);
    let mut i = 0;
    while i < stations.len() {
        let (t, _, node) = stations[i];
        let mut j = i;
        while j < stations.len() && stations[j].0 == t {
            match stations[j].1 {
                Event::GrazeEnd { toggle } => {
                    if graze.take().is_none() {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    if toggle {
                        present = !present;
                    }
                }
                Event::Toggle => {
                    if graze.is_some() {
                        return Err(DeclineKind::CylHoleFeature); // a toggle inside a run
                    }
                    present = !present;
                }
                Event::GrazeStart { body_above } => {
                    if graze.replace(body_above).is_some() {
                        return Err(DeclineKind::CylHoleFeature); // a run inside a run
                    }
                }
            }
            j += 1;
        }
        if j < stations.len() {
            let next = stations[j].2;
            let kind = match graze {
                Some(body_above) => Some(SegKind::Graze { body_above }),
                None if present => Some(mat),
                None => None,
            };
            if let Some(kind) = kind {
                match pieces.last_mut() {
                    Some((end, k)) if *k == kind && end[1] == node => end[1] = next,
                    _ => pieces.push(([node, next], kind)),
                }
            }
        }
        i = j;
    }
    if present || graze.is_some() {
        return Err(DeclineKind::CylHoleFeature); // an open boundary
    }
    Ok((pieces, grazes))
}

/// Which of the two rulings of `w` this pierce node sits on, `0` being a **tangent** wall's
/// single ruling — [`crate::combinatorics::ruling_side_signed`] asked of the point the name
/// denotes. The one spelling:
/// the chart's `ruling_name` reads it too.
///
/// ★ The `0` is asked of the point, not read off the name's root (`Double`): see
/// [`crate::combinatorics::ruling_side_signed`] for why the two do not agree.
pub(crate) fn node_ruling_side(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    w: &[nacre_exact::Rat; 4],
    n: NodeId,
) -> Option<i8> {
    let (_, cyl, _) = combinatorics::pierce_name(n)?;
    let (line, s) = combinatorics::pierce_meet(jd, cyl, def, n)?;
    combinatorics::ruling_side_signed(w, def, (&line, &s))
}

/// A pierce point's axis parameter: one of the two planes in its name is ⊥ the axis (a cap, a
/// rim), and **that class's** parameter is the point's — a rational, so the ordering along a
/// ruling needs no quadratic comparison at all.
///
/// ★ The name need not carry `wc` (a station the sweep minted does). A corner the
/// operand named (`Pierce{[cap, t], Double}`, a fillet's tangent corner) lies on `wc`'s ruling
/// too when `wc` contains it, and its parameter is read the same way: from its cap.
pub(super) fn node_axis_param(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<nacre_exact::Rat> {
    let (planes, _, _) = combinatorics::pierce_name(n)?;
    planes.iter().find_map(|&c| {
        crate::planes::axis_param_of_plane(&combinatorics::class_coeffs_rat(jd, c)?, def)
    })
}
