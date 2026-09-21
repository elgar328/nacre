use super::*;
/// A ring the rational chart road cannot name: a pierce corner (no rational coordinates) or an
/// arc step (no straight chart image). Such a ring takes [`point_in_mixed_ring`], which walks it
/// step by step in ℚ(√c) — the one predicate, asked by the circle arm of the arrangement's `cell_in_cell` of a
/// centre and by its polygon arm of a node.
pub(crate) fn ring_is_mixed(ring: &[RingEdge]) -> bool {
    ring.iter()
        .any(|e| matches!(e.carrier, Carrier::Arc(_)) || pierce_name(e.node).is_some())
}

/// **Parity of a rational point against a ring with pierce corners and arc steps** — the mixed
/// sibling of the chart road, asked only when `node_coords_rat` cannot name every corner.
///
/// The ray runs along the chart's own first axis (`Chart2dRat::axes` — one decision rule, not a
/// second basis spelling): in chart coordinates it is `{ y = q.y, x > q.x }`. Each ring step
/// answers by its carrier:
///
/// * a **line step** compares in ℚ(√c): the corners' chart coordinates are `QuadVal`s (a
///   rational corner lifted by `from_rat`, a pierce corner evaluated along its canonical meet
///   line — [`pierce_meet`]), the y-straddle is two signs, and "right of the
///   probe" is the 2-D orientation `(b−a) × (q−a)` — products stay in one radical because a
///   step carries at most one circle's corners; two *different* circles' corners on one step
///   make `checked_mul` refuse the radical mismatch and the whole answer abstains honestly;
/// * an **arc step** solves ray × circle exactly — the ray's plane `{ e2·p = q.y }` against the
///   class plane and the cylinder is [`nacre_exact::quad::plane_plane_cylinder`], the pierce
///   shape — and asks each root: right of the probe (chart-x as a `QuadVal`), and inside the
///   arc's CCW span (`circular_order_about_seam` on the carrier's own `end` pair, cyclic with
///   the wrap arm).
///
/// **A corner on the ray is a decision, not a tie**: the rule is the planar roads'
/// half-open one, spelled once in [`nacre_geom::intersect::ray_step_crossing`] — the corner is
/// counted by the step that leaves it upward. A line step reads that off its other end's sign;
/// an arc whose root is its own end reads it off its tangent there (the CCW tangent's side of
/// the ray's plane is minus [`ruling_side`], `arc_departure_side`'s
/// convention). Abstaining at the corner in both arms is not harmless: eight `NoClearRay`
/// cells of the crossing census were exactly that abstention exhausting every probe.
///
/// `None` is an honest abstention — the probe *on* the ring (at a corner, on a step along or
/// across the ray, at an arc root), a tangent ray, a horizontal tangent at an arc end the root
/// lands on, a seam-incident root, checked-`Rat` overflow — and the caller keeps its
/// `WitnessNotRational`. The kinds are counted under `tie_probe` in tests.
///
/// ★ **A radical mismatch is not among them, whatever the line above used to say.**
/// `QuadVal::common_radical` returns `None` there, but only after a `debug_assert!(false)` — so in
/// a test or debug build it **panics** rather than abstaining. The contract it states is
/// same-radical arithmetic, and a caller that could mix two must not reach it.
/// Where a point of an arc's circle sits against the arc's CCW span `lo → hi`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ArcSpan {
    Inside,
    Outside,
    /// Exactly at the `lo` end — the end the arc leaves counter-clockwise.
    AtLo,
    /// Exactly at the `hi` end — the end the arc leaves clockwise.
    AtHi,
}

/// **Is the circle point `root` in the arc's CCW span `lo → hi`, or at one of its ends?** — the
/// one spelling of the span question, read by the planar ring parity (an arc step of a mixed
/// ring, [`point_in_mixed_ring`]) and by the lateral face parity ([`loop_parity`]).
///
/// Cyclic in the seam chart. A seam-incident **end** is information, not a tie (the straddling
/// boss's alias corner sits exactly there): its θ is the chart boundary, so the span test
/// collapses to one comparison against the other end. Only a seam-incident **root** — the
/// crossing at the joint itself — abstains, as do two seam ends (one point twice; upstream
/// refuses it) and a zero-span arc. `None` for those and for checked arithmetic running out;
/// `tie_probe` says which.
pub(crate) fn arc_span(
    def: &nacre_topo::CylinderDef,
    e_lo: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    e_hi: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    root: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
) -> Option<ArcSpan> {
    use core::cmp::Ordering;
    use nacre_exact::quad::{MeetLine, QuadVal, SeamOrder, circular_order_about_seam};
    let (o, m, rd) = (def.origin(), def.dir(), def.ref_dir());
    let on_seam = |p: &(MeetLine, QuadVal)| -> Option<bool> {
        match circular_order_about_seam(&o, &m, &rd, (&p.0, &p.1), (&p.0, &p.1))? {
            SeamOrder::SeamIncident { first, .. } => Some(first),
            SeamOrder::Ordered(_) => Some(false),
        }
    };
    let ord = |p: &(MeetLine, QuadVal), qq: &(MeetLine, QuadVal)| -> Option<Ordering> {
        match circular_order_about_seam(&o, &m, &rd, (&p.0, &p.1), (&qq.0, &qq.1))? {
            SeamOrder::Ordered(o) => Some(o),
            SeamOrder::SeamIncident { .. } => None,
        }
    };
    if on_seam(root)? {
        // the root is the joint itself — a tie
        #[cfg(test)]
        tie_probe::mark(tie_probe::Tie::SeamRoot);
        return None;
    }
    Some(match (on_seam(e_lo)?, on_seam(e_hi)?) {
        // Two seam ends would be one point twice — upstream refuses it.
        (true, true) => {
            #[cfg(test)]
            tie_probe::mark(tie_probe::Tie::TwoSeamEnds);
            return None;
        }
        // From the seam CCW to `hi`: chart order θ ∈ (0, θ_hi). Only the seam end is unreachable
        // (a seam-incident root already returned above), so the **other** end is exactly what
        // can coincide, and `ord` answers it totally. **Three arms, one convention.**
        (true, false) => match ord(root, e_hi)? {
            Ordering::Equal => ArcSpan::AtHi,
            Ordering::Less => ArcSpan::Inside,
            Ordering::Greater => ArcSpan::Outside,
        },
        // From `lo` CCW back to the seam: θ ∈ (θ_lo, 2π).
        (false, true) => match ord(root, e_lo)? {
            Ordering::Equal => ArcSpan::AtLo,
            Ordering::Greater => ArcSpan::Inside,
            Ordering::Less => ArcSpan::Outside,
        },
        (false, false) => {
            let x0 = ord(root, e_lo)?;
            let x1 = ord(root, e_hi)?;
            match (x0 == Ordering::Equal, x1 == Ordering::Equal) {
                (true, true) => {
                    // zero-span arc cannot stand
                    #[cfg(test)]
                    tie_probe::mark(tie_probe::Tie::ZeroSpanArc);
                    return None;
                }
                (true, false) => ArcSpan::AtLo,
                (false, true) => ArcSpan::AtHi,
                (false, false) => {
                    let inside = match ord(e_lo, e_hi)? {
                        Ordering::Less => x0 == Ordering::Greater && x1 == Ordering::Less,
                        Ordering::Greater => x0 == Ordering::Greater || x1 == Ordering::Less,
                        Ordering::Equal => {
                            // zero-span arc cannot stand
                            #[cfg(test)]
                            tie_probe::mark(tie_probe::Tie::ZeroSpanArc);
                            return None;
                        }
                    };
                    if inside {
                        ArcSpan::Inside
                    } else {
                        ArcSpan::Outside
                    }
                }
            }
        }
    })
}

pub(crate) fn point_in_mixed_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc_coeffs: &[nacre_exact::Rat; 4],
    probe: &[nacre_exact::Rat; 3],
    ring: &[RingEdge],
) -> Option<bool> {
    #[cfg(test)]
    tie_probe::begin();
    let out = point_in_mixed_ring_inner(jd, cyls, wc_coeffs, probe, ring);
    #[cfg(test)]
    if out.is_none() {
        tie_probe::abstained();
    }
    out
}

fn point_in_mixed_ring_inner(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc_coeffs: &[nacre_exact::Rat; 4],
    probe: &[nacre_exact::Rat; 3],
    ring: &[RingEdge],
) -> Option<bool> {
    use nacre_exact::Orient;
    use nacre_exact::Rat;
    use nacre_exact::quad::{CylinderMeet, MeetLine, QuadVal};
    use nacre_geom::intersect::{ray_step_crossing, ray_straddle};
    let n = [wc_coeffs[0], wc_coeffs[1], wc_coeffs[2]];
    let chart = Chart2dRat::of_normal(&n)?;
    let (e1, e2) = chart.axes();
    let q = chart.project(probe)?;
    let (qx, qy) = (q[0], q[1]);
    // The ring's cylinders, for evaluating pierce corners — an arc or a ruling both carry theirs.
    // The class table is the identity's source; a chord-bounded ring carries no cylinder of its
    // own, so the ring's carriers cannot be the door.
    let def_of = |cyl: usize| cyls.get(cyl).map(|c| &c.def);
    let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    // A corner's chart coordinates, as `QuadVal`s.
    let corner = |nd: NodeId| -> Option<[QuadVal; 2]> {
        if let Some((_, cyl, _)) = pierce_name(nd) {
            let def = def_of(cyl)?;
            let (line, s) = pierce_meet(jd, cyl, def, nd)?;
            let (b, d) = (line.base(), line.dir());
            let coord = |e: &[Rat; 3]| -> Option<QuadVal> {
                QuadVal::from_rat(dot(&b, e)?).checked_add(&s.checked_mul_rat(dot(&d, e)?)?)
            };
            Some([coord(e1)?, coord(e2)?])
        } else {
            let p = node_coords_rat(jd, nd)?;
            let pr = chart.project(&p)?;
            Some([QuadVal::from_rat(pr[0]), QuadVal::from_rat(pr[1])])
        }
    };
    let k = ring.len();
    let mut inside = false;
    for i in 0..k {
        let (na, nb) = (ring[i].node, ring[(i + 1) % k].node);
        match &ring[i].carrier {
            // A ruling is a straight step like a plane step — same y-straddle, same orient2d;
            // its corners are pierce points, which `corner` already evaluates exactly.
            Carrier::Plane { .. } | Carrier::Ruling(_) => {
                let a = corner(na)?;
                let b = corner(nb)?;
                let ya = a[1].checked_sub(&QuadVal::from_rat(qy))?.sign();
                let yb = b[1].checked_sub(&QuadVal::from_rat(qy))?.sign();
                // A corner on the ray is the boundary only when the probe *is* that corner —
                // corner against the rational probe, never corner against corner. Every corner
                // is `a` of exactly one step, so this asks each corner once.
                if ya == Orient::Zero
                    && a[0].checked_sub(&QuadVal::from_rat(qx))?.sign() == Orient::Zero
                {
                    #[cfg(test)]
                    tie_probe::mark(tie_probe::Tie::ProbeAtCorner);
                    return None;
                }
                // A step lying along the ray straddles nothing (half-open: neither end is
                // above); it is the boundary iff the probe lies between its ends.
                if ya == Orient::Zero && yb == Orient::Zero {
                    let sa = a[0].checked_sub(&QuadVal::from_rat(qx))?.sign();
                    let sb = b[0].checked_sub(&QuadVal::from_rat(qx))?.sign();
                    if sa != sb {
                        #[cfg(test)]
                        tie_probe::mark(if sb == Orient::Zero {
                            tie_probe::Tie::ProbeAtCorner
                        } else {
                            tie_probe::Tie::ProbeOnStep
                        });
                        return None;
                    }
                    continue;
                }
                // The half-open rule (`ray_step_crossing`), read lazily: only a straddling step
                // pays for `orient2d`.
                if ray_straddle(ya, yb).is_none() {
                    continue;
                }
                // orient2d(a, b, q) = (b−a) × (q−a), all in one radical (or an honest None).
                let (qxv, qyv) = (QuadVal::from_rat(qx), QuadVal::from_rat(qy));
                let o = b[0]
                    .checked_sub(&a[0])?
                    .checked_mul(&qyv.checked_sub(&a[1])?)?
                    .checked_sub(
                        &b[1]
                            .checked_sub(&a[1])?
                            .checked_mul(&qxv.checked_sub(&a[0])?)?,
                    )?;
                match ray_step_crossing(ya, yb, o.sign()) {
                    Some(true) => inside = !inside,
                    Some(false) => {}
                    None => {
                        // probe on the step
                        #[cfg(test)]
                        tie_probe::mark(tie_probe::Tie::ProbeOnStep);
                        return None;
                    }
                }
            }
            Carrier::Arc(arc) => {
                // The ray's own plane: e2·p − qy = 0 (rational).
                let ray_plane = [e2[0], e2[1], e2[2], Rat::from_int(0).checked_sub(qy)?];
                let (o, m, r2) = (arc.def.origin(), arc.def.dir(), arc.def.r2());
                let roots = match nacre_exact::quad::plane_plane_cylinder(
                    wc_coeffs, &ray_plane, &o, &m, r2,
                )? {
                    CylinderMeet::Pair { line, s } => Some((line, s)),
                    CylinderMeet::Miss(_) | CylinderMeet::AxisParallelMiss(_) => None,
                    // A tangent ray, a ruling, or degenerate planes: ties and shapes the parity
                    // cannot count — abstain.
                    _ => {
                        #[cfg(test)]
                        tie_probe::mark(tie_probe::Tie::TangentRay);
                        return None;
                    }
                };
                let Some((line, s)) = roots else { continue };
                // The arc's CCW span: the step's ends oriented by the carried `ccw` bit — the
                // same convention every arc consumer reads (membership is direction-agnostic,
                // so the *set* is what the CCW pair names).
                let (lo_nd, hi_nd) = if arc.ccw { (na, nb) } else { (nb, na) };
                let e_lo = pierce_meet(jd, arc.cyl, &arc.def, lo_nd)?;
                let e_hi = pierce_meet(jd, arc.cyl, &arc.def, hi_nd)?;
                for root in s {
                    // Right of the probe along the ray: chart-x of the root.
                    let (bse, dir) = (line.base(), line.dir());
                    let x = QuadVal::from_rat(dot(&bse, e1)?)
                        .checked_add(&root.checked_mul_rat(dot(&dir, e1)?)?)?;
                    let xsign = x.checked_sub(&QuadVal::from_rat(qx))?.sign();
                    match xsign {
                        Orient::Zero => {
                            // root exactly at the probe
                            #[cfg(test)]
                            tie_probe::mark(tie_probe::Tie::ArcRootAtProbe);
                            return None;
                        }
                        Orient::Negative => continue,
                        Orient::Positive => {}
                    }
                    // Inside the CCW span end[0] → end[1]? Cyclic in the seam chart. A
                    // seam-incident **end** is information, not a tie (the straddling boss's
                    // alias corner sits exactly there): its θ is the chart boundary, so the
                    // span test collapses to one comparison against the other end. Only a
                    // seam-incident **root** — the crossing at the joint itself — abstains.
                    let rootp = (line.clone(), root);
                    // ★ **A root at the arc's own end is the corner on the ray in the arc's
                    // clothing**, and it takes the rule the line steps and the planar roads take
                    // (`ray_step_crossing`): the corner is counted by the step that leaves it
                    // upward. The arc leaves an end along its tangent, whose side of the ray's
                    // plane is minus `ruling_side` for the CCW tangent (`arc_departure_side`'s
                    // convention): `lo` departs CCW at `−rs`, `hi` — walked backwards — at `+rs`.
                    // With the end on the ray and right of the probe, that sign is the virtual
                    // tangent step's `side` too (`t × (q − E) = ty·(Ex − qx)`), so the call is
                    // `ray_step_crossing(Zero, ty, ty)` and no second cross product is spelled.
                    // A horizontal tangent (`rs` zero) is the genuine second-order tie.
                    let departs_across = |end: &(MeetLine, QuadVal), ccw: bool| -> Option<bool> {
                        let Some(rs) = ruling_side(&ray_plane, &arc.def, (&end.0, &end.1)) else {
                            #[cfg(test)]
                            tie_probe::mark(tie_probe::Tie::TangentAtEnd);
                            return None;
                        };
                        let ty = if (if ccw { -rs } else { rs }) > 0 {
                            Orient::Positive
                        } else {
                            Orient::Negative
                        };
                        #[cfg(test)]
                        tie_probe::arc_end_decided();
                        Some(ray_step_crossing(Orient::Zero, ty, ty) == Some(true))
                    };
                    // ★★★★★ **`Equal` used to be read as "outside the span"** — a root on the
                    // arc's own end counted nothing, silently: a confident wrong answer, not an
                    // abstention. Then the three arms abstained on it alike; now they decide it
                    // alike (`departs_across`), and the span itself is one spelling
                    // (`arc_span`) the lateral road reads too.
                    let contained = match arc_span(&arc.def, &e_lo, &e_hi, &rootp)? {
                        ArcSpan::Inside => true,
                        ArcSpan::Outside => false,
                        ArcSpan::AtLo => departs_across(&e_lo, true)?,
                        ArcSpan::AtHi => departs_across(&e_hi, false)?,
                    };
                    if contained {
                        inside = !inside;
                    }
                }
            }
        }
    }
    Some(inside)
}
