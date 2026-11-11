use super::*;
/// **What a lateral face leaves on a recorded wall class** — the ∥ sibling of
/// [`circle_on_class`], and like it, an answer **per extent** rather than one for the whole
/// ruling.
///
/// Asked only after the circle road answered "no circle"; **empty** (silently) unless the gate
/// **recorded** the pair — the wall within the radius (through the axis or offset from it), or
/// touching the cylinder along one ruling, where the face states only its half of a plane-pair
/// line on it ([`lateral_half_of_plane_pair_lines`]) — and declining rather than contributing
/// *partially* when it does but a piece cannot be stated — a half-contributed rectangle would
/// leave the class's 1-skeleton dangling, which is a worse lie than an honest incomplete-trace
/// mark.
///
/// The class cuts **two** rulings. Each boundary cycle is read against the class once
/// ([`cycle_on_wall`]) and leaves stations on both; each ruling is then swept along the axis
/// ([`sweep_line`]) over its stations and the face's whole rims — `Transversal` where the solid
/// straddles the wall, `Graze` where a cycle's own edge lies on the ruling, and nothing where the
/// face is not there (a panel, a chain rim, a band with holes, all by the one rule).
#[allow(clippy::too_many_arguments)]
pub(super) fn rulings_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    k: usize,
    which: SolidSide,
    cyls: &[crate::planes::WorkingCyl],
    aliases: &Aliases,
    found: &mut Aliases,
) -> Result<(Vec<RulingTrace>, Vec<Seg>), DeclineKind> {
    // **The gate's answer, first** — a class touching the cylinder
    // ([`WorkingCyl::touches`](crate::planes::WorkingCyl::touches)) gets the lateral's half of a
    // plane-pair line and no rectangle (below); otherwise
    // [`WorkingCyl::draws_rulings`](crate::planes::WorkingCyl::draws_rulings) — a class that cuts
    // the lateral, or one with no world statement that sections it on this its own arrangement
    // only, [`Silence::Sectioned`]. A pair without either record was proven clear (or never in
    // question), and contributing its rectangle anyway breaks the straddling family: the boss's own
    // wall class is a d=0 pair, and its coplanar bottom cap's chord lands on a line the plate
    // already traces. The rectangle goes exactly to the pairs whose record draws rulings.
    let sections = cyls.get(k).is_some_and(|cy| cy.sections_along(wc));
    // **A class touching the cylinder along one ruling** draws no ruling piece — the tangent
    // divides nothing — but where the face's own edge lies on that ruling with another plane
    // across it, the lateral states its half of that plane-pair line, by the rule the secant
    // rulings use, on the class's one ruling (`side = 0`).
    if cyls.get(k).is_some_and(|cy| cy.touches(wc)) {
        let Some(def) = cf.def.as_ref() else {
            return Err(DeclineKind::Ruling);
        };
        let LateralShape { cycles, .. } = lateral_shape(jd, fl, def, k)?;
        let above = axis_above(jd, wc, def).ok_or(DeclineKind::Ruling)?;
        let mut grazes = Vec::new();
        for (inner, ring) in &cycles {
            // A loop on another lateral states no plane-pair line (two laterals meeting are the
            // gate's to refuse, by name).
            let Some(nr) = ring.poly() else {
                continue;
            };
            let read = cycle_on_wall(jd, cf, def, wc, k, true, *inner, nr, aliases, cyls, found)?;
            let edges: Vec<usize> = read.quiet.iter().map(|&(_, i)| i).collect();
            lateral_half_of_plane_pair_lines(
                jd,
                wc,
                nr,
                &edges,
                which,
                &|_, _| Ok(above),
                &mut grazes,
            )?;
        }
        return Ok((Vec::new(), grazes));
    }
    if !cyls.get(k).is_some_and(|cy| cy.draws_rulings(wc)) {
        return Ok((Vec::new(), Vec::new()));
    }
    let Some(def) = cf.def.as_ref() else {
        return Err(DeclineKind::Ruling);
    };
    // The record is the crossing statement: a cut pair's plane runs within the
    // radius, so the wall meets the lateral in two rulings; the cycles' reading says where. The
    // premise is the gate's, asserted rather than re-derived — and it survived
    // the tangent wall opening precisely because a tangency is written to `tangencies`
    // and **not** as a cut: "within" still means within. A recorded pair's class has world
    // coefficients (the gate records a crossing only from them); a sectioned class's premise is
    // the gate's interval reading of the same inequality.
    debug_assert!(
        sections
            || jd.world_rat(wc).is_some_and(|w| {
                nacre_exact::point_plane_clearance_rat(&w, &def.origin(), def.r2())
                    == nacre_exact::Orient::Negative
            }),
        "a recorded pair's plane runs within the radius"
    );
    // The face's shape from its cycles — the same reading the ⊥ road makes.
    let LateralShape { rims, cycles, .. } = lateral_shape(jd, fl, def, k)?;
    let mut read: Vec<(&combinatorics::NamedRing, WallCycle)> = Vec::with_capacity(cycles.len());
    for (inner, ring) in &cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair.
        let Some(nr) = ring.poly() else {
            return Err(DeclineKind::CylFaceHole);
        };
        let r = cycle_on_wall(jd, cf, def, wc, k, false, *inner, nr, aliases, cyls, found)?;
        read.push((nr, r));
    }
    let mat = SegKind::Transversal {
        mat: cf.orient_sign,
    };
    // **Where two points stand on a ruling** — their sections' order on its generator, the one
    // spelling ([`combinatorics::height_order_on_generator`]); a point it cannot place declines.
    let order = |a: NodeId, b: NodeId| {
        combinatorics::height_order_on_generator(jd, k, def, a, b).ok_or(DeclineKind::Ruling)
    };
    let mut out = Vec::with_capacity(2);
    let mut grazes = Vec::new();
    for side in [1i8, -1i8] {
        let mut stations: Vec<(LineEvent, NodeId)> = Vec::new();
        for &(c, _) in &rims {
            let node = combinatorics::crossing_on_ruling(jd, def, c, wc, k, side)
                .map_err(|_| DeclineKind::Ruling)?;
            stations.push((LineEvent::Toggle, node));
        }
        for (nr, r) in &read {
            // **The lateral's half of a plane-pair line** on this ruling is stated in the plane
            // vocabulary; its side is read exactly as a run's on this ruling is (`body_side`): the
            // same cycle, the same ruling, the same traversal.
            let edges: Vec<usize> = r
                .quiet
                .iter()
                .filter(|&&(s, _)| s == side)
                .map(|&(_, i)| i)
                .collect();
            lateral_half_of_plane_pair_lines(
                jd,
                wc,
                nr,
                &edges,
                which,
                &|a, b| body_side(jd, cf, def, wc, k, side, a, b).map(|(_, above)| above),
                &mut grazes,
            )?;
            stations.extend(
                r.stations
                    .iter()
                    .filter(|&&(s, ..)| s == side)
                    .map(|&(_, e, x)| (e, x)),
            );
        }
        combinatorics::try_sort_by(&mut stations, |x, y| {
            Ok(order(x.1, y.1)?.then(x.0.cmp(&y.0)))
        })?;
        let same: Vec<bool> = stations
            .windows(2)
            .map(|w| Ok(order(w[0].1, w[1].1)? == core::cmp::Ordering::Equal))
            .collect::<Result<_, DeclineKind>>()?;
        // Two names at one height on this ruling are one point — the order is exact. Recorded for
        // the next round and read under the first name, as the plane road's sweep does.
        for i in 0..same.len() {
            if same[i] && stations[i].1 != stations[i + 1].1 {
                found.record_same_point_traced(jd, cyls, stations[i].1, stations[i + 1].1)?;
                stations[i + 1].1 = stations[i].1;
            }
        }
        let pieces = sweep_line(&stations, &same, |n| n, mat)?;
        #[cfg(test)]
        if pieces
            .iter()
            .any(|(_, kind)| matches!(kind, SegKind::Graze { .. }))
        {
            let o = def.origin().map(|x| x.to_f64());
            let at = |n: NodeId| ruling_probe::height(jd, def, n);
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

/// What one boundary cycle of a lateral leaves on a wall class ([`cycle_on_wall`]).
struct WallCycle {
    /// Stations on the class's rulings, `(side, event, name)` — `side` the ruling's identity
    /// ([`node_ruling_side`]).
    stations: Vec<(i8, LineEvent, NodeId)>,
    /// The cycle's edges lying on a ruling of the class with another plane across them, `(side,
    /// edge)` — plane-pair lines the lateral ends on ([`lateral_half_of_plane_pair_lines`]).
    quiet: Vec<(i8, usize)>,
}

/// Where an arc of a cycle meets a wall class's rulings: the points strictly inside it, and the
/// point at each of its two ends (`at[0]` the edge's start, `at[1]` its end) — each with the
/// ruling it lies on.
#[derive(Clone, Default)]
struct ArcMeet {
    inside: Vec<(i8, NodeId)>,
    at: [Option<(i8, NodeId)>; 2],
}

/// **One boundary cycle of a lateral read against a wall class** — a class along the axis, which
/// meets the lateral in its rulings: two for a secant class, one for a tangent one (`tangent`).
///
/// **Whether the boundary crosses the class where it meets it is the walk's**
/// ([`combinatorics::ring_against_plane`]), the one rule every line road reads a ring by: a
/// contact — one corner on a ruling, or a run of the cycle's edges along it — is a crossing
/// exactly when its two off-line neighbours lie on opposite sides of the class. A ruling edge
/// with a point on the class lies in it (both run along the axis), so it is `On`; an arc leaving
/// a node on the class leaves to the side its tangent points ([`combinatorics::arc_departure`]).
///
/// **Where an arc meets a ruling is read once per arc and side**, before the walk: the class
/// crosses the arc's section at one point on each of its two rulings
/// ([`combinatorics::crossing_on_ruling`], the name every road gives that point), which the arc
/// holds strictly inside, at an end, or not at all ([`combinatorics::theta_place`]). An end there
/// is one point under two names, an exact proof recorded for the next round
/// ([`Aliases::record_same_point`]). The walk counts an arc's crossings off this reading — an end
/// is never inside, so the count needs no end the walk found — and a `Crossing` is named by it. A
/// tangent class meets each section's plane in a line that touches the cylinder or misses it, so
/// no arc crosses it and nothing is read.
///
/// The walk's features become stations on the ruling they lie on ([`node_ruling_side`]):
/// - a **crossing** inside an arc toggles the face there;
/// - a **corner** on a ruling — one node, two arcs meeting there — toggles it exactly when the
///   walk's flanks differ, named by the arriving arc's end;
/// - a **run** of the cycle's own edges held by the class is a `Graze` over its extent, the side
///   the face occupies along it `σ·τ·κ·side` (`body_side`), toggling the face exactly when its
///   flanks differ — a panel's arcs both lie on its own side (no toggle: absent above and below);
///   a chain's arrive from the outer half and leave into the inner (toggle: transversal below,
///   absent above); a band's notch hole has both on the hole's side (no toggle);
/// - a **run** held by another plane `t` is the plane pair `wc ∩ t`, which `t`'s own trace states
///   as a graze run on this class: no station, and its edges go out as `quiet` for the lateral's
///   half in the plane vocabulary.
///
/// Declined by name: wherever a point, an order or a side cannot be formed exactly — an arc's
/// meeting with a ruling, a node's ruling, two points on a generator or about the axis — and a
/// cycle the walk cannot read or that lies wholly in the class (`Ruling`); a run whose nodes stand
/// on two rulings, whose edges mix the two holders, of no length, or held by another plane with
/// flanks that differ — what the sweep's events cannot say over that extent (`CylHoleFeature`).
#[allow(clippy::too_many_arguments)]
fn cycle_on_wall(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    k: usize,
    tangent: bool,
    inner: bool,
    nr: &combinatorics::NamedRing,
    aliases: &Aliases,
    cyls: &[crate::planes::WorkingCyl],
    found: &mut Aliases,
) -> Result<WallCycle, DeclineKind> {
    use crate::combinatorics::{EdgeMeet, Feature, RingWalk, ThetaPlace, Wall};
    #[cfg(not(test))]
    let _ = inner;
    let ring = &nr.triples;
    let n = ring.len();
    // **Identity is the table's** — a crossing named here and a corner the cycle carries can be
    // one point under two names, and the seeds have usually told the table so
    // (`Aliases::record_on_cylinder`). Where the geometry proves two names one point the tables do
    // not hold yet (one angle on an arc's circle, one height on a ruling), the proof is recorded in
    // `found` for the next round, not refused.
    let same =
        |found: &Aliases, x: NodeId, y: NodeId| aliases.one_point(x, y) || found.one_point(x, y);
    let mut meets: Vec<ArcMeet> = vec![ArcMeet::default(); n];
    if !tangent {
        for (i, meet) in meets.iter_mut().enumerate() {
            let Wall::Arc {
                ccw: nu, plane: c, ..
            } = nr.walls[i]
            else {
                continue;
            };
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            let (lo, hi) = if nu { (a, b) } else { (b, a) };
            for side in [1i8, -1i8] {
                let node = combinatorics::crossing_on_ruling(jd, def, c, wc, k, side)
                    .map_err(|_| DeclineKind::Ruling)?;
                let end = if same(found, node, a) {
                    Some(a)
                } else if same(found, node, b) {
                    Some(b)
                } else {
                    match combinatorics::theta_place(jd, k, def, lo, hi, node)
                        .ok_or(DeclineKind::Ruling)?
                    {
                        ThetaPlace::At(e) => {
                            found.record_same_point_traced(jd, cyls, [lo, hi][e], node)?;
                            Some([lo, hi][e])
                        }
                        ThetaPlace::Inside => {
                            meet.inside.push((side, node));
                            None
                        }
                        ThetaPlace::Outside => None,
                    }
                };
                // A closed arc's two ends are one node.
                if end == Some(a) {
                    meet.at[0] = Some((side, node));
                }
                if end == Some(b) {
                    meet.at[1] = Some((side, node));
                }
            }
        }
    }
    let features = match combinatorics::ring_against_plane(
        jd,
        cyls,
        ring,
        wc,
        |i| match nr.walls[i] {
            Wall::Ruling { .. } => Some(EdgeMeet::On),
            Wall::Arc { cyl: c, ccw, plane } => {
                combinatorics::arc_departure(jd, cyls, ring[i], plane, wc, c, ccw)
            }
            Wall::Plane(_) => None,
        },
        |i, _| u8::try_from(meets[i].inside.len()).ok(),
    ) {
        RingWalk::Met(f) => f,
        // A cycle lying wholly in the class is no face; a node the walk cannot name is the
        // refusal every other consumer makes of it.
        RingWalk::AllOn | RingWalk::Unnameable => return Err(DeclineKind::Ruling),
    };
    let held = |i: usize| matches!(nr.walls[i], Wall::Ruling { plane: t, .. } if t == wc);
    // The ruling a run stands on: every node of it on one.
    let run_side = |nodes: &[NodeId]| -> Result<i8, DeclineKind> {
        let mut out = None;
        for &x in nodes {
            let s = node_ruling_side(jd, def, wc, x).ok_or(DeclineKind::Ruling)?;
            match out {
                None => out = Some(s),
                // An edge whose two ends answer differently is not a ruling at all — it
                // would have to cross the plane through the axis. Named rather than silently
                // mis-placed.
                Some(p) if p != s => return Err(DeclineKind::CylHoleFeature),
                Some(_) => {}
            }
        }
        out.ok_or(DeclineKind::Ruling)
    };
    let mut read = WallCycle {
        stations: Vec::new(),
        quiet: Vec::new(),
    };
    let mut run_ends: Vec<(i8, NodeId)> = Vec::new();
    for f in features {
        match f {
            Feature::Crossing { edge, nth, .. } => {
                let Wall::Arc { ccw, .. } = nr.walls[edge] else {
                    unreachable!(
                        "a ruling runs along the axis as the class does: it lies in the class or \
                         wholly on one side, so its two ends never straddle it"
                    );
                };
                let (side, node) = match meets[edge].inside[..] {
                    [x] => x,
                    [x, y] => {
                        let x_first =
                            combinatorics::met_first(jd, k, def, ccw, ring[edge], x.1, y.1)
                                .ok_or(DeclineKind::Ruling)?;
                        if (nth == 0) == x_first { x } else { y }
                    }
                    // Where the walk asked, it was given this count. Where it did not — ends on
                    // opposite sides — the arc's section meets the class once on each ruling, and
                    // an arc between the two sides holds exactly one of those points strictly
                    // inside, and the ends' sides and the points' places are both proven readings.
                    _ => unreachable!("the walk counts an arc's crossings off its reading"),
                };
                read.stations.push((side, LineEvent::Toggle, node));
            }
            Feature::Run { first, len: 1, .. } => {
                if tangent {
                    continue;
                }
                let _ = first;
                return Err(DeclineKind::Ruling);
            }
            Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } => {
                let edges: Vec<usize> = (0..len - 1).map(|j| (first + j) % n).collect();
                let nodes: Vec<NodeId> = (0..len).map(|j| ring[(first + j) % n]).collect();
                let (start, end) = (nodes[0], nodes[len - 1]);
                if edges.iter().all(|&i| held(i)) {
                    // A tangent class seated on the lateral states its edge itself
                    // (`trace_plane`'s tangent ruling); nothing to add here.
                    if tangent {
                        continue;
                    }
                    let side = run_side(&nodes)?;
                    let (tau, body_above) = body_side(jd, cf, def, wc, k, side, start, end)?;
                    let (n_lo, n_hi) = if tau > 0 { (start, end) } else { (end, start) };
                    #[cfg(test)]
                    ruling_probe::GRAZE_SIDE.push(ruling_probe::GrazeSide {
                        inner,
                        body_above,
                        ny: jd.planes[wc].plane.normal().as_array()[1],
                        origin: def.origin().map(|x| x.to_f64()),
                        span: [
                            ruling_probe::height(jd, def, n_lo),
                            ruling_probe::height(jd, def, n_hi),
                        ],
                    });
                    read.stations
                        .push((side, LineEvent::GrazeStart { body_above }, n_lo));
                    read.stations.push((
                        side,
                        LineEvent::GrazeEnd {
                            toggle: flanks_differ,
                        },
                        n_hi,
                    ));
                    run_ends.extend([(side, start), (side, end)]);
                } else if edges.iter().any(|&i| held(i)) || flanks_differ {
                    return Err(DeclineKind::CylHoleFeature);
                } else {
                    let side = if tangent { 0 } else { run_side(&nodes)? };
                    read.quiet.extend(edges.iter().map(|&i| (side, i)));
                    run_ends.extend([(side, start), (side, end)]);
                }
            }
        }
    }
    // A run's end at an arc's meeting height on its ruling is that meeting too: one height on one
    // generator, an exact proof recorded for the next round.
    for meet in &meets {
        for &(side, node) in meet.at.iter().flatten() {
            for &(s, r) in &run_ends {
                if s == side
                    && !same(found, r, node)
                    && combinatorics::height_order_on_generator(jd, k, def, r, node)
                        == Some(core::cmp::Ordering::Equal)
                {
                    found.record_same_point_traced(jd, cyls, r, node)?;
                }
            }
        }
    }
    Ok(read)
}

/// The side of `wc` the face's body lies on along an edge of its cycle that runs on ruling `side`
/// from `start` to `end`, with the ring's travel along the axis `τ` — see [`plus_theta_is_above`],
/// the one spelling of the product.
#[allow(clippy::too_many_arguments)]
fn body_side(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    k: usize,
    side: i8,
    start: NodeId,
    end: NodeId,
) -> Result<(i8, bool), DeclineKind> {
    let tau: i8 = match combinatorics::height_order_on_generator(jd, k, def, end, start)
        .ok_or(DeclineKind::Ruling)?
    {
        core::cmp::Ordering::Greater => 1,
        core::cmp::Ordering::Less => -1,
        core::cmp::Ordering::Equal => return Err(DeclineKind::CylHoleFeature), // no length
    };
    let theta_dot_stored: i8 = if plus_theta_is_above(jd, wc, side).ok_or(DeclineKind::Ruling)? {
        1
    } else {
        -1
    };
    let body_above =
        i32::from(material_theta_sign(cf.orient_sign, tau)) * i32::from(theta_dot_stored) > 0;
    Ok((tau, body_above))
}

/// Which of the two rulings of `wall` this pierce node sits on, `0` being a **tangent** wall's
/// single ruling — [`crate::combinatorics::ruling_side_signed`] asked of the point the name
/// denotes. The one spelling:
/// the chart's `ruling_name` reads it too.
///
/// The `0` is asked of the point, not read off the name's root (`Double`): see
/// [`crate::combinatorics::ruling_side_signed`] for why the two do not agree.
pub(crate) fn node_ruling_side(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wall: usize,
    n: NodeId,
) -> Option<i8> {
    let (_, cyl, _) = combinatorics::pierce_name(n)?;
    combinatorics::pierce_meet(jd, cyl, def, n)?.ruling_side(jd, wall, def)
}

/// **The lateral's half of a plane-pair line on a ruling of class `wc`** — the cycle's `edges`
/// lying there with another plane `t` across them ([`cycle_on_wall`]'s `quiet`).
///
/// Such an edge lies on the ruling with the face *across* it being another plane `t` —
/// the fillet's tangent edge, when `wc` passes through the axis and so contains the tangent
/// ruling; a half cylinder's flat edge, when `wc` is a wall tangent to it there. That line is a
/// plane-pair line (`wc ∩ t`), and the lateral **ends** there, so the plane vocabulary states it —
/// the half of the one rule for a line that is a ruling and a plane-pair line at once (the other
/// half: where the lateral runs on across the line, the line is the other solid's edge and the
/// ruling vocabulary states it — `split_rulings`' fold). The face `t` states the line from its
/// side (a graze run). But the solid's material near that line lies on **both** sides of `wc`, or
/// on one side twice — and an edge's mask is the sum of its faces' statements (`edge_mask`: two
/// grazes of one solid flip both bits when they lie on opposite sides, none when on one). So the
/// lateral states its side too — as a `Seg` on that line, so `merge_coincident` joins it with
/// `t`'s — and nothing in the ruling vocabulary.
///
/// The side, `body_above`, is the caller's: on a secant ruling the face leaves the line along
/// `±θ` and crosses `wc` there (first order — `body_side`); on a tangent one it leaves along `wc`
/// itself and bends to the axis on both sides of the line, so its side is the axis side
/// ([`axis_above`] — second order, a boss's and a bore's alike).
///
/// The two orders are the ladder CGAL's conic arrangement climbs to order two curves that touch:
/// by the first derivative, then the second.
fn lateral_half_of_plane_pair_lines(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    nr: &combinatorics::NamedRing,
    edges: &[usize],
    which: SolidSide,
    body_above: &dyn Fn(NodeId, NodeId) -> Result<bool, DeclineKind>,
    grazes: &mut Vec<Seg>,
) -> Result<(), DeclineKind> {
    let n = nr.triples.len();
    for &i in edges {
        let combinatorics::Wall::Ruling { plane: t, .. } = nr.walls[i] else {
            unreachable!("a run on a ruling is made of ruling edges");
        };
        let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
        let above = body_above(a, b)?;
        let pin = |x: NodeId| combinatorics::pin_for(jd, wc, t, x).ok_or(DeclineKind::Ruling);
        grazes.push(Seg {
            wall: t,
            end: [a, b],
            end_h: [pin(a)?, pin(b)?],
            solid: which,
            kind: SegKind::Graze { body_above: above },
        });
    }
    Ok(())
}

/// Whether the axis of the cylinder `def` lies above class `wc` in the label frame (`wc`'s stored
/// normal) — the one side a lateral tangent to `wc` lies on, the cylinder bending to its axis on
/// both sides of the tangent line. `None` where the class has no world statement.
pub(super) fn axis_above(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    def: &nacre_topo::CylinderDef,
) -> Option<bool> {
    let w = jd.planes[wc].world.as_ref()?;
    Some(w.sense.sign() * nacre_exact::plane_residual_sign(&w.name, def.origin()) > 0)
}
