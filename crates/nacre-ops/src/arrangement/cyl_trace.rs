use super::*;
/// What a lateral face leaves on a ⊥ plane class — see [`circle_on_class`], which answers it
/// **per angular extent** ([`CircleSpan`]) rather than once for the whole circle.
///
/// The empty answer is "this class leaves nothing here": a non-⊥ class carries no circle at all,
/// a ⊥ class beyond both rims never meets the face, and a class the gate proved silent about the
/// cylinder stands clear of it. `Err` is "the face cannot answer" (no rim span) — declining beats a
/// silently missing circle, which would corrupt every label on the class.
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
    // A class the gate proved clear of this cylinder carries none of it — the gate's record, not a
    // question this road could ask of a class with no world statement.
    if cyls.get(cyl).is_some_and(|cy| cy.is_silent(wc)) {
        return Ok(Vec::new());
    }
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
    // ⊥ to the axis, decided **totally** (`class_carries_circle` clears denominators and answers
    // in integers, so this cannot decline for want of bits). A non-⊥ class carries no circle at
    // all — a miss, like a parallel plane, not a decline.
    if !class_carries_circle(&n, &m) {
        return Ok(Vec::new());
    }
    // **The face's boundary as the tracer names it — its cycles — and nothing else.** A
    // band's two rims and its holes are its loops (`combinatorics::lateral_cycles`); a panel or a
    // chain rim is a cycle like a hole, carved out of an outer answer that may be **absent**
    // rather than out of a whole circle.
    let LateralShape {
        range,
        rims,
        cycles,
    } = lateral_shape(jd, cf, fl, def)?;
    let Some(t) = crate::planes::axis_param_of_plane(&coeffs, def) else {
        return Err(DeclineKind::CylSpan);
    };
    // The **outer** answer, which the holes below carve out of:
    //
    // **A rim is a graze, not a miss.** Answering "the lateral does not reach this plane", on the
    // premise that coplanarity would fold a rim into a seated class, would drop the touch: a
    // *cylinder* face never becomes seated on a plane class. Every planar wall meeting a class
    // along an edge contributes a `Graze`, and `edge_mask`'s precedence (`Graze > Seated`) exists
    // exactly for the corner where the two disagree. On a convex cap they agree and the omission
    // would be invisible; at the **reflex dihedral of a blind bore's ceiling** they disagree, the
    // seated rule would flip the wrong label bit, and the *next* boolean on that solid would come
    // back `CylinderGateUndecided` — while a through bore, which has no ceiling, would be fine.
    //
    // Which side the body is on: the face runs from `span[0]` toward `span[1]`, so at the low rim
    // it lies toward `+t` and at the high rim toward `−t`.
    let up = crate::planes::plus_t_is_above(world, def);
    if t < range[0] || range[1] < t {
        // Beyond the face: this class does not meet it at all, and every cycle lives inside the
        // range, so there is nothing for the walk below to find either.
        return Ok(Vec::new());
    }
    // **The outer answer is optional.** A whole rim at this station grazes as a
    // whole circle; strictly inside the range the face crosses the class *somewhere* — a
    // connected face covers every interior station at some θ — and the cycles below carve out
    // where it does not; at an end of the range that is **no** rim (a panel's arc, a chain's top)
    // nothing stands: only the runs the cycles lay on the class speak, and the rest is absent.
    // `Crosses` is exact only because a run whose flanks differ is declined below
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
    // **The hole is read by the walk every other face's boundary is read by** — the face's
    // own loop, against this class, through [`combinatorics::ring_against_plane`]. An interval
    // derived from the loop's ⊥ carriers would be exact for a chart rectangle and a *premise* for
    // anything else; the walk asks the ring and needs no premise about the hole's shape.
    let mut carved: Vec<Carved> = Vec::new();
    for (inner, ring) in &cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair, not
        // this road's. It has no ring to walk, so it declines rather than passing unread.
        let Some(nr) = ring.poly() else {
            return Err(DeclineKind::CylFaceHole);
        };
        cycle_on_class(jd, cyls, cf, def, wc, cyl, up, *inner, nr, &mut carved)?;
    }
    let spans = if carved.is_empty() {
        match outer {
            Some(on) => vec![CircleSpan { arc: None, on }],
            None => Vec::new(),
        }
    } else {
        assemble_spans(jd, cyl, &cyls[cyl].def, &carved, outer)?
    };
    Ok(spans)
}

/// **A lateral face's shape as its cycles state it** — what both lateral roads read of a
/// lateral's loops: its axial `range`, its rims — each a closed rim edge or a rim split into arcs
/// ([`arcs_on_one_class`]) — (station and plane class, in station order), and every other cycle —
/// a panel, a chain rim, a hole — with whether it is an inner loop. The
/// stations are read
/// from the planes' world coefficients, so `range` is a second derivation of
/// `CylFaceInfo::footprint.span` (the model side's: the outer loop's ⊥ carriers), and the two are held
/// equal where they meet.
///
/// Declines by the one name: `None` cycles (a loop could not be named), a class with no exact
/// station, more than two rims, two rims at one
/// station, or a rim that is not at an end of the range — a whole circle bounds the face, so it is
/// an extreme of it.
pub(super) fn lateral_shape(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    def: &nacre_topo::CylinderDef,
) -> Result<LateralShape, DeclineKind> {
    use combinatorics::LoopRing;
    let Some(cycles) = &fl.cycles else {
        return Err(DeclineKind::CylSpan);
    };
    // One spelling of "a class's axis station" (`bands::param_opt`), not a third.
    let station = |c: usize| crate::bands::param_opt(jd, c, def).ok_or(DeclineKind::CylSpan);
    let mut rims: Vec<(Rat, usize)> = Vec::new();
    let mut rest: Vec<(bool, LoopRing)> = Vec::new();
    let mut range: Option<[Rat; 2]> = None;
    let mut widen = |t: Rat| {
        range = Some(match range {
            None => [t, t],
            Some([lo, hi]) => [if t < lo { t } else { lo }, if hi < t { t } else { hi }],
        });
    };
    for (i, ring) in cycles.iter().enumerate() {
        match ring {
            LoopRing::Rim { plane } => {
                let t = station(*plane)?;
                widen(t);
                rims.push((t, *plane));
            }
            ring => {
                // **A rim split into arcs is a rim.** Corners on the rim's own plane (a pin over a
                // plate's edge: the edge's wall meets the rim there) cut the closed rim edge into
                // arcs, and the loop arrives as a ring — but a simple loop of arcs of one circle
                // goes once round it, so it is the whole circle at that station, as
                // `LoopRing::Rim` is. Read as a hole or a panel it lies wholly on its own class
                // and has nothing to bound with. The ring stays as the model has it
                // (`combinatorics::FaceLoops::cycles`): «is this a rim» is a question about the
                // face's shape, asked here, where both lateral roads read it — a different question
                // from the gate's whole rim (`planes::cyl_gate::only_whole_rims`).
                if let Some(nr) = ring.poly()
                    && let Some(plane) = arcs_on_one_class(&nr.walls)
                {
                    let t = station(plane)?;
                    widen(t);
                    rims.push((t, plane));
                    continue;
                }
                // The cycle's arcs each lie on their rim's ⊥ class; its rulings have no station of
                // their own.
                if let Some(nr) = ring.poly() {
                    for w in &nr.walls {
                        if let crate::combinatorics::Wall::Arc { plane, .. } = *w {
                            widen(station(plane)?);
                        }
                    }
                }
                rest.push((i > 0, ring.clone()));
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

/// The one plane class every wall of a ring is an arc on, `None` when a wall is not an arc or two
/// arcs ride different classes. A lateral's arcs are its own cylinder's, so they lie on one
/// circle, and a simple loop of arcs of one circle runs one way round it — asserted, as the
/// producer's premise.
fn arcs_on_one_class(walls: &[crate::combinatorics::Wall]) -> Option<usize> {
    use crate::combinatorics::Wall;
    let Some(&Wall::Arc { plane, .. }) = walls.first() else {
        return None;
    };
    if !walls
        .iter()
        .all(|w| matches!(*w, Wall::Arc { plane: c, .. } if c == plane))
    {
        return None;
    }
    debug_assert!(
        walls.windows(2).all(|w| match (w[0], w[1]) {
            (Wall::Arc { ccw: d0, .. }, Wall::Arc { ccw: d1, .. }) => d0 == d1,
            _ => false,
        }),
        "a rim split into arcs runs one way round: {walls:?}"
    );
    Some(plane)
}

/// **Lay the carved extents over the outer answer and hand back *exclusive* spans.**
///
/// **Exclusive, not "the whole circle plus overrides".** The cheaper spelling would emit the
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
    let order = combinatorics::circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
        combinatorics::CircleOrderFail::Undecided => DeclineKind::CylSpan,
        combinatorics::CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
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
    // One answer everywhere: the circle is undivided, whatever the holes touch.
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
/// *Material is on the left of the ring's direction of travel* — the sentence
/// `trace_plane::run_body_above` states for a planar face's on-line run. On a lateral face, with
/// `r̂` the outward radial direction, `m̂` the axis and `θ̂` increasing θ (counter-clockwise about
/// `m̂`), cylindrical coordinates are right-handed so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The
/// face's outward is `σ·r̂` with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], travel is
/// `τ·m̂`, and left of travel is
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
/// **One atom, three readers.** [`combinatorics::RulingCarrier::side`] is
/// `sign((x − o) · (m̂ × n̂_r))` against the class's *rational* name, and the scalar triple product
/// gives `(x − o) · (m̂ × n̂) = n̂ · ((x − o) × m̂) = −r·(n̂ · θ̂)`, so `sign(n̂_r · θ̂) = −side`;
/// the world name's sense (`κ`, [`crate::planes::WorldName`] — the truth's, since the name alone
/// carries no direction and the plane cache is a rounded image) carries that to the plane's own
/// facing. The three consumers are the
/// ruling sweep's graze side, [`ruling_interior_is_even`], and the chart's vertical read — and the
/// whole point of naming it is that none of them spells the product a second time. Two of them
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
    let kappa = jd.canon_sense(wc);
    Some(-(kappa * side) == 1)
}

/// **One boundary cycle — a hole, a panel, a chain rim — read against one ⊥ class**: the walk's
/// features turned into extents. A hole is carved out of a whole-circle outer answer; a panel or a
/// chain is carved out of one that may be absent — the same arms, because neither reads
/// "inside" or "outside": only the ring's winding, below.
///
/// **Which way round is decided by the ring's own winding, and by nothing else.** A plane cuts a
/// circle in two points, and "which of the two arcs is the hole" is the whole difficulty: a hole
/// whose two ruling edges lie on **one** plane cannot be told apart by that plane's sides. The
/// universal boundary convention answers it locally instead — *material is on the left of the
/// ring's direction of travel* — the same sentence `trace_plane::run_body_above` states for a
/// planar face's on-line run, and it holds for an outer ring, a hole ring, a notch and a reflex
/// corner alike.
///
/// **Derivation.** Write `r̂` for the cylinder's outward radial direction, `m̂` for its axis and `θ̂`
/// for increasing θ (counter-clockwise about `m̂`, which is what
/// `nacre_exact::quad::circular_order_about_seam` ranks and what an arc [`Piece`]'s `end` runs along).
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
    inner: bool,
    nr: &combinatorics::NamedRing,
    out: &mut Vec<Carved>,
) -> Result<(), DeclineKind> {
    let _ = inner; // the ledger's, under test
    #[cfg(test)]
    let start = out.len();
    // **`true` here is not a shrug — the meet on this road is a *circle*, and an arc can lie on
    // one.** A rim arc of the hole at this very axis parameter *is* the class's meet, so "not an
    // arc" would cut runs that are genuinely continuous. What decides it is the carrier, and the
    // `Run` arm below asks that directly (declining, rather than splitting, is this road's scope).
    // No edge of a lateral's ring crosses a ⊥ class strictly inside: its arcs lie in ⊥ planes
    // parallel to `wc` and its rulings are straight — so every edge answers `0`, whichever the
    // wall says it is.
    let features = match combinatorics::ring_against_plane(
        jd,
        cyls,
        &nr.triples,
        wc,
        |_| Some(combinatorics::EdgeMeet::On),
        |_, _| Some(0),
    ) {
        combinatorics::RingWalk::Met(f) => f,
        // A hole ring lying wholly in the class has no thickness to bound anything with (a rim
        // split into arcs lies wholly in its class too, and `lateral_shape` reads it as a rim
        // before this road), and a node the walk cannot name is the same refusal every other
        // consumer makes of it.
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
                // **Which way the arc runs, and so which side the face is on, is the
                // producer's to say — not the flank's.** The flank (the side the ring's
                // off-class neighbours sit on) gives the ring's *interior* side — right for a
                // convex hole and wrong for a wrapping rim, whose highest arc has both
                // neighbours below it and the face below too. The arc's own direction about the
                // axis (`Wall::Arc`'s `ccw`, the
                // producer's convention) settles both: walked with sense `ν`, material lies
                // along `σ·ν·m̂`, so the face is above the class exactly when that agrees with
                // the class's stored normal (`up`), and the carved extent is the arc as walked.
                // Measured: this reading and the flank's are equal on every on-class arc of
                // today's holes (50 of 50).
                for k in 0..len - 1 {
                    // **"Both ends on the class" is not "the edge is on the circle."** The
                    // walk reads a run off its *nodes* and this road's blanket `On`; on a plane the
                    // edge between two on-line nodes
                    // follows, and on a **cylinder** it does not — a tilted carrier meets the
                    // lateral in an ellipse, which can cross this class at both ends without
                    // lying on it. Taking that for a rim arc would state an extent along a curve
                    // that is not there. The carrier says it directly: an edge on the circle
                    // `wc ∩ cylinder` is an arc whose rim's plane is that class.
                    // The check belongs here, per edge, where the answer is used — not on the
                    // row ("neither ⊥ nor ∥: an ellipse, outside this vocabulary").
                    let edge = (first + k) % n;
                    let crate::combinatorics::Wall::Arc { ccw: nu, plane, .. } = nr.walls[edge]
                    else {
                        return Err(DeclineKind::CylHoleFeature);
                    };
                    if plane != wc {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    let body_above = (sigma * if nu { 1 } else { -1 } > 0) == up;
                    let a = ring[edge];
                    let b = ring[(edge + 1) % n];
                    let arc = if nu { [a, b] } else { [b, a] };
                    out.push(Carved {
                        arc,
                        on: Some(CylOnClass::Grazes { body_above }),
                    });
                    #[cfg(test)]
                    arc_probe::record(jd, cyl, def, inner, &out[start..]);
                }
            }
            combinatorics::Feature::Crossing { edge, from, .. } => {
                // The crossed edge's carrier, **carried** from the producer rather than re-derived
                // from the two endpoint names.
                // The ruling is read off the wall alone — its cylinder, its side and the plane the
                // side was measured against — not off the corner's name: a corner on three planes
                // is named by the pair its definition chose, which need not hold `j`.
                let crate::combinatorics::Wall::Ruling {
                    cyl: ncyl,
                    side,
                    plane: j,
                    ..
                } = nr.walls[edge]
                else {
                    return Err(DeclineKind::CylHoleFeature);
                };
                // **The crossing is the *same ruling*, restated for this class — through
                // the one door the planar scan uses** ([`combinatorics::crossing_on_ruling`]): the ruling's side,
                // measured against its wall `j` at the hole's own corner — the side the wall
                // carries (`curved_wall`) — picks the root of `{wc, j}` on the cylinder; restating
                // the corner's root by the axis senses of the two ⊥ classes would be a second
                // spelling of the same point. `j` must hold a ruling (a wall plane within the
                // radius) and `wc` be ⊥, which the door checks; anything else is not this feature.
                let wdef = &cyls.get(ncyl).ok_or(DeclineKind::CylHoleFeature)?.def;
                let cut = combinatorics::crossing_on_ruling(jd, wdef, wc, j, ncyl, side)
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
    // boundary is the whole assignment — the circle twin of the sweep along a line
    // ([`sweep_line`]).
    if !cuts.is_empty() {
        let nodes: Vec<NodeId> = cuts.iter().map(|c| c.0).collect();
        let order = combinatorics::circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
            combinatorics::CircleOrderFail::Undecided => DeclineKind::CylSpan,
            combinatorics::CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
        })?;
        for k in 0..order.len() {
            let (a, ahead) = cuts[order[k]];
            let (b, next_ahead) = cuts[order[(k + 1) % order.len()]];
            if ahead == next_ahead {
                return Err(DeclineKind::CylHoleFeature); // the boundaries did not alternate
            }
            if ahead {
                out.push(Carved {
                    arc: [a, b],
                    on: None,
                });
                #[cfg(test)]
                arc_probe::record(jd, cyl, def, inner, &out[start..]);
            }
        }
    }
    Ok(())
}
