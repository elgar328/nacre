use super::*;
/// What a lateral face leaves on a ⊥ plane class — see [`circle_on_class`], which answers it
/// **per angular extent** ([`CircleSpan`]) rather than once for the whole circle.
///
/// The empty answer is "this class leaves nothing here": a non-⊥ class carries no circle at all,
/// and a ⊥ class beyond both rims never meets the face. `Err` is "the face cannot answer" (no rim
/// span, or a class with no exact description the gate would already have refused) — declining
/// beats a silently missing circle, which would corrupt every label on the class.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CylOnClass {
    /// The class cuts the face in two: the solid straddles it here.
    Crosses,
    /// The class carries one of the face's **rims**: the wall touches it along that circle with
    /// its body to one side. `body_above` is that side, about the class's **stored** normal —
    /// the frame every cell label is written in.
    Grazes { body_above: bool },
}

pub(super) fn circle_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    cyl: usize,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    // The world description — the cylinder's statement is world, and a comparison across two
    // frames is a silently wrong answer, not a slow one.
    let Some(world) = jd.planes[wc].world.as_ref() else {
        return Err(DeclineKind::CylSpan);
    };
    let Some(coeffs) = world.name.narrow().copied() else {
        return Err(DeclineKind::CylSpan);
    };
    // No world statement for this lateral (a rotated or frame-borne truth): the same decline
    // the class side makes — both sides of this comparison have to speak about the world.
    let Some(def) = cf.def.as_ref() else {
        return Err(DeclineKind::CylSpan);
    };
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let m = def.dir();
    // ⊥ to the axis, decided **totally** (`parallel_rat` clears denominators and answers in
    // integers, so this cannot decline for want of bits). A non-⊥ class carries no circle at
    // all — a miss, like a parallel plane, not a decline.
    if !nacre_exact::parallel_rat(&n, &m) {
        return Ok(Vec::new());
    }
    // ★★★★★ **The face's boundary as the tracer names it — its cycles — and nothing else.** A
    // band's two rims and its holes come from the outer loop cut at its slits
    // (`combinatorics::lateral_cycles`), so a hole the assembly spliced into the outer walk is
    // a hole here like any other; a panel or a chain rim is a cycle like a hole, carved
    // out of an outer answer that may be **absent** rather than out of a whole circle.
    let LateralShape {
        range,
        rims,
        cycles,
    } = lateral_shape(jd, cf, fl, def)?;
    let Some(t) = crate::planes::axis_param_of_plane(&coeffs, def) else {
        return Err(DeclineKind::CylSpan);
    };
    // The **outer** answer, which the holes below carve out of. Unchanged from when it was the
    // whole answer:
    //
    // ★★ **A rim is a graze, not a miss.** This used to answer "the lateral does not reach this
    // plane", on the premise that coplanarity would have folded a rim into a seated class — but a
    // *cylinder* face never becomes seated on a plane class, so the touch was simply dropped.
    // Every planar wall meeting a class along an edge contributes a `Graze`; the lateral
    // contributed nothing, and `edge_mask`'s precedence (`Graze > Seated`) exists exactly for the
    // corner where the two disagree. On a convex cap they agree and the omission is invisible; at
    // the **reflex dihedral of a blind bore's ceiling** they disagree, the seated rule then flipped
    // the wrong label bit, and the *next* boolean on that solid came back
    // `CylinderGateUndecided` — while a through bore, which has no ceiling, was fine.
    //
    // Which side the body is on: the face runs from `span[0]` toward `span[1]`, so at the low rim
    // it lies toward `+t` and at the high rim toward `−t`.
    let up = crate::planes::plus_t_is_above(world, def);
    if t < range[0] || range[1] < t {
        // Beyond the face: this class does not meet it at all, and every cycle lives inside the
        // range, so there is nothing for the walk below to find either.
        return Ok(Vec::new());
    }
    // ★★★★★ **The outer answer is optional.** A whole rim at this station grazes as a
    // whole circle; strictly inside the range the face crosses the class *somewhere* — a
    // connected face covers every interior station at some θ — and the cycles below carve out
    // where it does not; at an end of the range that is **no** rim (a panel's arc, a chain's top)
    // nothing stands: only the runs the cycles lay on the class speak, and the rest is absent.
    // ★ `Crosses` is exact only because a run whose flanks differ is declined below
    // (`cycle_on_class`): such a run is a parity boundary in θ, and the face beyond it at this
    // station would be absent with no crossing to carve it — a staircase chain, not today's
    // population. Pairing run ends with crossings is that road's generalization, not this one.
    let outer: Option<CylOnClass> = if rims.iter().any(|(s, _)| *s == t) {
        Some(CylOnClass::Grazes {
            body_above: if t == range[0] { up } else { !up },
        })
    } else if range[0] < t && t < range[1] {
        Some(CylOnClass::Crosses)
    } else {
        None
    };
    // ★★★★★ **The hole is read by the walk every other face's boundary is read by** — the face's
    // own loop, against this class, through [`combinatorics::ring_against_plane`]. An interval
    // derived from the loop's ⊥ carriers would be exact for a chart rectangle and a *premise* for
    // anything else; the walk asks the ring and needs no premise about the hole's shape.
    let mut carved: Vec<Carved> = Vec::new();
    for (kind, ring) in &cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair, not
        // this road's. It has no ring to walk, so it declines rather than passing unread.
        let Some(nr) = ring.poly() else {
            return Err(DeclineKind::CylFaceHole);
        };
        cycle_on_class(jd, cyls, cf, def, wc, cyl, up, *kind, nr, &mut carved)?;
    }
    let spans = if carved.is_empty() {
        match outer {
            Some(on) => vec![CircleSpan { arc: None, on }],
            None => Vec::new(),
        }
    } else {
        assemble_spans(jd, cyl, &cyls[cyl].def, &carved, outer)?
    };
    #[cfg(test)]
    cycle_probe::record(t, outer, rims.len(), &cycles, carved.len(), spans.len());
    Ok(spans)
}

/// **A lateral face's shape as its cycles state it** — what both lateral roads read of a
/// lateral's loops: its axial `range`, its whole rims (station and plane class, in station order),
/// and every other cycle — a panel, a chain rim, a hole — with its kind. The stations are read
/// from the planes' world coefficients, so `range` is a second derivation of
/// `CylFaceInfo::footprint.span` (the model side's: the outer loop's ⊥ carriers), and the two are held
/// equal where they meet.
///
/// Declines by the one name: `None` cycles (the outer loop could not be cut or named), a rim kind
/// that is not a whole circle, a class with no exact station, more than two rims, two rims at one
/// station, or a rim that is not at an end of the range — a whole circle bounds the face, so it is
/// an extreme of it.
pub(super) fn lateral_shape(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    def: &nacre_topo::CylinderDef,
) -> Result<LateralShape, DeclineKind> {
    use combinatorics::{CycleKind, LoopRing};
    let Some(cycles) = &fl.cycles else {
        return Err(DeclineKind::CylSpan);
    };
    // One spelling of "a class's axis station" (`bands::param_opt`), not a third.
    let station = |c: usize| crate::bands::param_opt(jd, c, def).ok_or(DeclineKind::CylSpan);
    let mut rims: Vec<(Rat, usize)> = Vec::new();
    let mut rest: Vec<(CycleKind, LoopRing)> = Vec::new();
    let mut range: Option<[Rat; 2]> = None;
    let mut widen = |t: Rat| {
        range = Some(match range {
            None => [t, t],
            Some([lo, hi]) => [if t < lo { t } else { lo }, if hi < t { t } else { hi }],
        });
    };
    for (kind, ring) in cycles {
        match (kind, ring) {
            (CycleKind::Rim, LoopRing::Rim { plane }) => {
                let t = station(*plane)?;
                widen(t);
                rims.push((t, *plane));
            }
            (CycleKind::Rim, _) => return Err(DeclineKind::CylSpan),
            (_, ring) => {
                // The cycle's arcs (an edge with a stated sense) each lie on a ⊥ class; its
                // rulings have no station of their own.
                if let Some(nr) = ring.poly() {
                    for (i, w) in nr.walls.iter().enumerate() {
                        if nr.arc_ccw[i].is_none() {
                            continue;
                        }
                        let crate::combinatorics::Wall::Plane(c) = *w else {
                            return Err(DeclineKind::CylSpan);
                        };
                        widen(station(c)?);
                    }
                }
                rest.push((*kind, ring.clone()));
            }
        }
    }
    let Some(range) = range else {
        return Err(DeclineKind::CylSpan);
    };
    if range[0] == range[1] || rims.len() > 2 {
        return Err(DeclineKind::CylSpan);
    }
    rims.sort_by_key(|r| r.0);
    if rims.windows(2).any(|w| w[0].0 == w[1].0)
        || rims.iter().any(|(t, _)| *t != range[0] && *t != range[1])
    {
        return Err(DeclineKind::CylSpan);
    }
    debug_assert_eq!(
        Some(range),
        cf.footprint.span,
        "the cycles' stations by name and the face's range by model disagree"
    );
    Ok(LateralShape {
        range,
        rims,
        cycles: rest,
    })
}

/// **Lay the carved extents over the outer answer and hand back *exclusive* spans.**
///
/// ★★★★★ **Exclusive, not "the whole circle plus overrides".** The cheaper spelling would emit the
/// outer answer over the whole circle and let the hole's arcs sit on top, leaning on `edge_mask`'s
/// `Graze > Transversal` precedence to pick the winner. That precedence is for **different faces
/// meeting on one edge**, not for one face's own extent — borrowing it would be a category error,
/// and it cannot express "absent" at all. So the circle is cut at every boundary and each piece
/// carries exactly one answer.
///
/// Adjacent pieces with the same answer are **merged**, so a boundary where nothing changes leaves
/// no cut behind: an extra split point would be a degree-2 vertex the winding walk has to have a
/// turn for.
fn assemble_spans(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    carved: &[Carved],
    outer: Option<CylOnClass>,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    let mut nodes: Vec<NodeId> = carved.iter().flat_map(|c| c.arc).collect();
    nodes.sort_unstable();
    nodes.dedup();
    let (order, _) = circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
        CircleOrderFail::Undecided => DeclineKind::CylSpan,
        CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
    })?;
    let n = order.len();
    // θ position of each node, keyed by its index in the sorted `nodes`.
    let mut at = vec![0usize; n];
    for (k, &i) in order.iter().enumerate() {
        at[i] = k;
    }
    let place = |x: NodeId| {
        at[nodes
            .binary_search(&x)
            .expect("every boundary node is in the set")]
    };
    // The answer over elementary arc `k` (from `order[k]` to `order[k+1]`), or `None` where no
    // cycle reaches: there the outer answer stands — which may itself be "absent".
    let mut ans: Vec<Option<Option<CylOnClass>>> = vec![None; n];
    for c in carved {
        let (a, b) = (place(c.arc[0]), place(c.arc[1]));
        if a == b {
            return Err(DeclineKind::CylHoleFeature); // an extent of no length, or of the whole circle
        }
        let mut k = a;
        while k != b {
            if ans[k].is_some() {
                return Err(DeclineKind::CylHoleFeature); // two cycles claiming one extent
            }
            ans[k] = Some(c.on);
            k = (k + 1) % n;
        }
    }
    let answer = |k: usize| ans[k].unwrap_or(outer);
    // One answer everywhere: the circle was never divided, whatever the holes touched.
    let Some(start) = (0..n).find(|&k| answer(k) != answer((k + n - 1) % n)) else {
        return Ok(match answer(0) {
            Some(on) => vec![CircleSpan { arc: None, on }],
            None => Vec::new(),
        });
    };
    let mut out = Vec::new();
    let mut k = 0;
    while k < n {
        let i = (start + k) % n;
        let a = nodes[order[i]];
        let this = answer(i);
        let mut len = 1;
        while k + len < n && answer((start + k + len) % n) == this {
            len += 1;
        }
        if let Some(on) = this {
            out.push(CircleSpan {
                arc: Some([a, nodes[order[(start + k + len) % n]]]),
                on,
            });
        }
        k += len;
    }
    Ok(out)
}

/// **Which way round a lateral face's material lies, for a boundary edge travelling along the
/// axis** — the one place that sign is made, and the third application of one convention.
///
/// *Material is on the left of the ring's direction of travel* — the sentence [`run_body_above`]
/// states for a planar face's on-line run. On a lateral face, with `r̂` the outward radial
/// direction, `m̂` the axis and `θ̂` increasing θ (counter-clockwise about `m̂`), cylindrical
/// coordinates are right-handed so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The face's outward is `σ·r̂`
/// with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], travel is `τ·m̂`, and left of travel is
///
/// ```text
///   n_out × travel = (σ·r̂) × (τ·m̂) = σ·τ·(r̂ × m̂) = −σ·τ·θ̂
/// ```
///
/// so the answer is `−σ·τ`: the θ direction the face occupies. Two readers need it — the ⊥ road,
/// where a ruling edge crossing the circle puts the **hole** in the opposite θ direction, and the
/// ∥ road, where a ruling edge lying on the class puts the **face** to one side of the wall — and
/// they must not spell it twice.
pub(super) fn material_theta_sign(orient_sign: i8, travel_up: i8) -> i8 {
    -(orient_sign * travel_up)
}

/// **Which way `+θ̂` points across a ∥ wall, in the frame a label is written in** — `+1` when
/// leaving a ruling counter-clockwise about the axis enters the side the wall's plane faces.
///
/// ★★★ **One atom, three readers.** [`combinatorics::RulingCarrier::side`] is
/// `sign((x − o) · (m̂ × n̂_r))` against the class's *rational* name, and the scalar triple product
/// gives `(x − o) · (m̂ × n̂) = n̂ · ((x − o) × m̂) = −r·(n̂ · θ̂)`, so `sign(n̂_r · θ̂) = −side`;
/// the world name's sense (`κ`, [`crate::planes::WorldName`] — the truth's, since the name alone
/// carries no direction and the plane cache is a rounded image) carries that to the plane's own
/// facing. The three consumers are the
/// ruling sweep's graze side, [`ruling_interior_is_even`], and the chart's vertical read — and the
/// whole point of naming it is that none of them spells the product a second time. ★ Two of them
/// read a **label** and stay in the stored frame; the one that picks a **cell** —
/// [`ruling_interior_is_even`] — crosses into the chart's frame with `frame_sign`, and it is the
/// only reader that does.
pub(crate) fn plus_theta_is_above(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    side: i8,
) -> Option<bool> {
    // A tangent ruling (`side == 0`) is never a recorded crossing, so nothing labelled by this
    // sign reaches it; the day one does, the answer is a measurement, not a sign.
    debug_assert_ne!(side, 0, "a tangent ruling has no +θ side to be above");
    let kappa = jd.planes[wc].world.as_ref()?.sense.sign();
    Some(-(kappa * side) == 1)
}

/// **One boundary cycle — a hole, a panel, a chain rim — read against one ⊥ class**: the walk's
/// features turned into extents. A hole is carved out of a whole-circle outer answer; a panel or a
/// chain is carved out of one that may be absent — the same arms, because neither reads
/// "inside" or "outside": only the ring's winding, below.
///
/// ★★★★★ **Which way round is decided by the ring's own winding, and by nothing else.** A plane
/// cuts a circle in two points, and "which of the two arcs is the hole" is the whole difficulty:
/// a hole whose two ruling edges lie on **one** plane cannot be told apart by that plane's sides.
/// The universal boundary convention answers it locally instead — *material is on the left of the
/// ring's direction of travel* — the same sentence [`run_body_above`] states for a planar face's
/// on-line run, and it holds for an outer ring, a hole ring, a notch and a reflex corner alike.
///
/// **Derivation.** Write `r̂` for the cylinder's outward radial direction, `m̂` for its axis and `θ̂`
/// for increasing θ (counter-clockwise about `m̂`, which is what
/// `nacre_exact::quad::circular_order_about_seam` ranks and what [`MergedArc::end`] runs along).
/// Cylindrical coordinates are right-handed, so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The face's outward
/// is `σ·r̂` with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], and "left of travel" is
/// `n_out × travel`.
///
/// - **A ruling edge crossing the class** travels `τ·m̂`. Material lies along
///   `(σ·r̂) × (τ·m̂) = −σ·τ·θ̂`, so the **hole** lies along `+σ·τ·θ̂`: the hole runs counter-clockwise
///   from the crossing exactly when `σ·τ = +1`.
/// - **An arc edge lying on the class** (one of the hole's own rims) travels `ν·θ̂`. Material lies
///   along `(σ·r̂) × (ν·θ̂) = σ·ν·m̂`, so `ν = σ·μ` where `μ` is the axis direction the face occupies
///   — which is the opposite of the side the run's flanks sit on. That fixes the rim arc's CCW
///   orientation without reading a coordinate.
///
/// `up` is `plus_t_is_above`: whether the class's **stored** normal points along `+m̂`. It is the
/// bridge between the two frames here, because [`combinatorics::side_of`] answers in the class
/// root's **outward** frame and a label's "above" is the stored one — the correction is
/// [`crate::planes::WorkingPlane::frame_sign`], and forgetting it is a silently mirrored answer.
#[allow(clippy::too_many_arguments)]
fn cycle_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    cyl: usize,
    up: bool,
    kind: combinatorics::CycleKind,
    nr: &combinatorics::NamedRing,
    out: &mut Vec<Carved>,
) -> Result<(), DeclineKind> {
    let _ = kind; // the ledger's, under test
    // ★★ **`true` here is not a shrug — the meet on this road is a *circle*, and an arc can lie on
    // one.** A rim arc of the hole at this very axis parameter *is* the class's meet, so "not an
    // arc" would cut runs that are genuinely continuous. What decides it is the carrier, and the
    // `Run` arm below asks that directly (declining, rather than splitting, is this road's scope).
    // No edge of a lateral's ring crosses a ⊥ class strictly inside: its arcs lie in ⊥ planes
    // parallel to `wc` and its rulings are straight — so every edge answers `0` (the ring's arcs
    // carry `Wall::Plane`, so the wall could not tell them apart anyway).
    let features = match combinatorics::ring_against_plane(
        jd,
        cyls,
        &nr.triples,
        wc,
        |_| Some(combinatorics::EdgeMeet::On),
        |_, _| Some(0),
    ) {
        combinatorics::RingWalk::Met(f) => f,
        // A hole ring lying wholly in the class has no thickness to bound anything with, and a
        // node the walk cannot name is the same refusal every other consumer makes of it.
        combinatorics::RingWalk::AllOn | combinatorics::RingWalk::Unnameable => {
            return Err(DeclineKind::CylFaceHole);
        }
    };
    let ring = &nr.triples;
    let n = ring.len();
    let sigma = i32::from(cf.orient_sign);
    let fs = i32::from(jd.planes[wc].frame_sign);
    let axis_of = |stored_side: i32| if up { stored_side } else { -stored_side };
    // Crossings, each with the direction the hole runs in from it. Paired below, after the θ sort.
    let mut cuts: Vec<(NodeId, bool)> = Vec::new();
    for f in features {
        match f {
            combinatorics::Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } => {
                // A run whose flanks differ is a rim the hole continues *through* — it is both an
                // extent and a parity toggle, and no fixture makes one. Declined and counted
                // rather than guessed at.
                if flanks_differ {
                    return Err(DeclineKind::CylHoleFeature);
                }
                // A single on-line vertex with equal flanks is the hole's corner touching the
                // circle at a point. A point takes no extent away, so there is nothing to carve.
                if len < 2 {
                    continue;
                }
                // ★★★★★ **Which way the arc runs, and so which side the face is on, is the
                // producer's to say — not the flank's.** The flank (the side the ring's
                // off-class neighbours sit on) gives the ring's *interior* side — right for a
                // convex hole and wrong for a wrapping rim, whose highest arc has both
                // neighbours below it and the face below too. The arc's own direction about the
                // axis (`NamedRing::arc_ccw`, the
                // producer's convention) settles both: walked with sense `ν`, material lies
                // along `σ·ν·m̂`, so the face is above the class exactly when that agrees with
                // the class's stored normal (`up`), and the carved extent is the arc as walked.
                // ☑ Measured: this reading and the flank's are equal on every on-class arc of
                // today's holes (50 of 50).
                for k in 0..len - 1 {
                    // ★★★★★ **"Both ends on the class" is not "the edge is on the circle."** The
                    // walk reads a run off its *nodes* and this road's blanket `On`; on a plane the
                    // edge between two on-line nodes
                    // follows, and on a **cylinder** it does not — a tilted carrier meets the
                    // lateral in an ellipse, which can cross this class at both ends without
                    // lying on it. Taking that for a rim arc would state an extent along a curve
                    // that is not there. The carrier says it directly: an edge on the circle
                    // `wc ∩ cylinder` lies in `wc`, so its far face is of that class.
                    // ★ The check belongs here, per edge, where the answer is used — not on the
                    // row ("neither ⊥ nor ∥: an ellipse, outside this vocabulary").
                    let edge = (first + k) % n;
                    if nr.walls[edge] != crate::combinatorics::Wall::Plane(wc) {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    // An arc on the class with no stated sense is not a lateral's arc at all.
                    let Some(nu) = nr.arc_ccw[edge] else {
                        return Err(DeclineKind::CylHoleFeature);
                    };
                    let body_above = (sigma * if nu { 1 } else { -1 } > 0) == up;
                    let a = ring[edge];
                    let b = ring[(edge + 1) % n];
                    let arc = if nu { [a, b] } else { [b, a] };
                    #[cfg(test)]
                    arc_probe::record(jd, cyl, def, kind, arc);
                    out.push(Carved {
                        arc,
                        on: Some(CylOnClass::Grazes { body_above }),
                    });
                }
            }
            combinatorics::Feature::Crossing { edge, from, .. } => {
                // The crossed edge's carrier, **carried** from the producer rather than re-derived
                // from the two endpoint names.
                let crate::combinatorics::Wall::Plane(j) = nr.walls[edge] else {
                    return Err(DeclineKind::CylHoleFeature);
                };
                let start = ring[edge];
                let Some((planes, ncyl, _)) = combinatorics::pierce_name(start) else {
                    return Err(DeclineKind::CylHoleFeature);
                };
                if !planes.contains(&j) {
                    return Err(DeclineKind::CylHoleFeature);
                }
                // ★★★★★ **The crossing is the *same ruling*, restated for this class — through
                // the one door the planar scan uses** ([`crossing_on_ruling`]): the ruling's side,
                // measured against the wall `j` at the hole's own corner ([`node_ruling_side`],
                // the predicate `curved_wall` carries on a ruling edge), picks the root of
                // `{wc, j}` on the cylinder. Its previous spelling restated the corner's *root* by
                // the axis senses of the two ⊥ classes (`ε·sign(k)`, the order of the roots along
                // `ℓ`) — the same point derived the other way round, and two spellings of one
                // rule were one too many. `j` must hold a ruling (a wall plane within the radius)
                // and `wc` be ⊥, which the door checks; anything else is not this feature.
                let cj = combinatorics::class_coeffs_rat(jd, j).ok_or(DeclineKind::CylSpan)?;
                let wdef = &cyls.get(ncyl).ok_or(DeclineKind::CylHoleFeature)?.def;
                let side =
                    node_ruling_side(jd, wdef, &cj, start).ok_or(DeclineKind::CylHoleFeature)?;
                let cut = crossing_on_ruling(jd, wdef, wc, j, ncyl, side)
                    .map_err(|_| DeclineKind::CylHoleFeature)?;
                // Which way the hole runs from here: the travel's axis sense, then the winding.
                // The side the edge leaves comes from the walk, for the same reason the run's
                // flank does.
                let s0 = i32::from(from);
                let tau = axis_of(if s0 * fs < 0 { 1 } else { -1 });
                // The face occupies `material_theta_sign`; the hole is the other way.
                let hole_ccw = material_theta_sign(cf.orient_sign, tau as i8) < 0;
                cuts.push((cut, hole_ccw));
            }
        }
    }
    // A closed ring enters and leaves the circle equally often.
    if cuts.len() % 2 != 0 {
        return Err(DeclineKind::CylHoleFeature);
    }
    // The crossings alternate in θ, so sorting them and pairing each "hole ahead" with the next
    // boundary is the whole assignment — the circle twin of the parity sweep
    // `trace_transversal_face` runs along a line.
    if !cuts.is_empty() {
        let nodes: Vec<NodeId> = cuts.iter().map(|c| c.0).collect();
        let (order, _) = circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
            CircleOrderFail::Undecided => DeclineKind::CylSpan,
            CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
        })?;
        for k in 0..order.len() {
            let (a, ahead) = cuts[order[k]];
            let (b, next_ahead) = cuts[order[(k + 1) % order.len()]];
            if ahead == next_ahead {
                return Err(DeclineKind::CylHoleFeature); // the boundaries did not alternate
            }
            if ahead {
                let arc = [a, b];
                #[cfg(test)]
                arc_probe::record(jd, cyl, def, kind, arc);
                out.push(Carved { arc, on: None });
            }
        }
    }
    Ok(())
}
