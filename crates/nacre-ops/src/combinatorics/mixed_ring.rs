use super::*;
/// A ring the rational chart road cannot name: a pierce corner (no rational coordinates) or an
/// arc step (no straight chart image). Such a ring takes [`point_in_mixed_ring`], which walks it
/// step by step in quadratic radicals — the one predicate, asked by the circle arm of the
/// arrangement's `cell_in_cell` of a centre and by its polygon arm of a node.
pub(crate) fn ring_is_mixed(ring: &[RingEdge]) -> bool {
    ring.iter()
        .any(|e| matches!(e.carrier, Carrier::Arc(_)) || pierce_name(e.node).is_some())
}

/// Where a point of an arc's circle sits against the arc's CCW span `lo → hi`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ArcSpan {
    Inside,
    Outside,
    /// Exactly at the `lo` end — the end the arc leaves counter-clockwise.
    AtLo,
    /// Exactly at the `hi` end — the end the arc leaves clockwise.
    AtHi,
    /// At the one node of a **closed** arc ([`ArcEnds::Closed`]), which is both its ends: the arc
    /// leaves it counter-clockwise and arrives back at it. ★ Its own answer, not `AtLo`: a reader
    /// that counts an end by the way the arc leaves it (the planar ring parity) has to count both
    /// here, and one folded into `AtLo` would count half and answer wrong without a refusal.
    AtBothEnds,
}

/// The ends of an arc as [`arc_span`] reads them — by **topology**, which the caller knows: an
/// arc whose two ends are one node is the circle that node cuts, its whole turn, while two ends
/// that merely stand at one place under two names are a zero-span arc. Coordinates alone cannot
/// tell the two apart.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ArcEnds<'a> {
    /// Two nodes, `lo` then `hi` counter-clockwise.
    Open(
        &'a (nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
        &'a (nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    ),
    /// One node, the arc running its whole turn from it back to it.
    Closed(&'a (nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)),
}

/// **Is the circle point `root` in the arc's CCW span `lo → hi`, or at one of its ends?** — the
/// one spelling of the span question, read by the planar ring parity (an arc step of a mixed
/// ring, [`point_in_mixed_ring`]) and by the lateral face parity ([`loop_parity`]).
///
/// Cyclic, with the seam point first ([`nacre_exact::quad::SeamOrder::seam_first`]) — an end or
/// a root on the seam is ordered like any other point, so a span that runs across the seam holds
/// a root there and one that does not leaves it out. A seam-incident end is how the straddling
/// boss's alias corner arrives. A closed arc holds every point of its circle but its node, which
/// is [`ArcSpan::AtBothEnds`]. `None` for a zero-span arc (two distinct ends at one place, the
/// seam's two names for one point included — upstream refuses it) and for checked arithmetic
/// running out.
pub(crate) fn arc_span(
    def: &nacre_topo::CylinderDef,
    ends: ArcEnds<'_>,
    root: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
) -> Option<ArcSpan> {
    use core::cmp::Ordering;
    use nacre_exact::quad::{MeetLine, QuadVal, circular_order_about_seam};
    let (o, m, rd) = (def.origin(), def.dir(), def.ref_dir());
    let ord = |p: &(MeetLine, QuadVal), qq: &(MeetLine, QuadVal)| -> Option<Ordering> {
        Some(circular_order_about_seam(&o, &m, &rd, (&p.0, &p.1), (&qq.0, &qq.1))?.seam_first())
    };
    let (e_lo, e_hi) = match ends {
        ArcEnds::Closed(node) => {
            return Some(if ord(root, node)? == Ordering::Equal {
                ArcSpan::AtBothEnds
            } else {
                ArcSpan::Inside
            });
        }
        ArcEnds::Open(e_lo, e_hi) => (e_lo, e_hi),
    };
    let x0 = ord(root, e_lo)?;
    let x1 = ord(root, e_hi)?;
    Some(match (x0 == Ordering::Equal, x1 == Ordering::Equal) {
        (true, true) => {
            // zero-span arc cannot stand
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
                    return None;
                }
            };
            if inside {
                ArcSpan::Inside
            } else {
                ArcSpan::Outside
            }
        }
    })
}

/// A line step's corner in chart coordinates, `r₀ + r₁·√c` — one radical for both coordinates.
struct ChartCorner {
    r0: [nacre_exact::Rat; 2],
    r1: [nacre_exact::Rat; 2],
    c: nacre_exact::Rat,
}

impl ChartCorner {
    /// The radical the corner carries — `0` for a rational corner, whatever `c` its root spelled.
    fn radical(&self) -> nacre_exact::Rat {
        let zero = nacre_exact::Rat::from_int(0);
        if self.r1 == [zero, zero] {
            zero
        } else {
            self.c
        }
    }

    /// The two coordinates as `QuadVal`s, for the one-radical questions a corner is asked alone.
    fn values(&self) -> Option<[nacre_exact::QuadVal; 2]> {
        let v = |k: usize| nacre_exact::QuadVal::new(self.r0[k], self.r1[k], self.c);
        Some([v(0)?, v(1)?])
    }
}

/// **Parity of a rational point against a ring with pierce corners and arc steps** — the mixed
/// sibling of the chart road, asked only when `node_coords_rat` cannot name every corner.
///
/// The ray runs along the chart's own first axis (`Chart2dRat::axes` — one decision rule, not a
/// second basis spelling): in chart coordinates it is `{ y = q.y, x > q.x }`. Each ring step
/// answers by its carrier:
///
/// * a **line step** reads its corners as chart coordinates `r₀ + r₁·√c` (a rational corner with
///   `r₁ = 0`, a pierce corner evaluated along its canonical meet line — [`pierce_meet`]), the
///   y-straddle is two signs, and "right of the probe" is the 2-D orientation `(b−a) × (q−a)`.
///   Its two corners need not spell one radical — the ends of one ruling named by two plane
///   pairs carry `c` and `k²c`, and two circles' corners carry two radicals — so the orientation
///   is one `QuadVal` expression where they do and the two-radical sign (`biquad_sign`) where
///   they do not;
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
/// across the ray, at an arc root inside the arc's closed span), a tangent ray, a horizontal
/// tangent at an arc end the root lands on, checked-`Rat` overflow — and the caller keeps its
/// `WitnessNotRational`. A root on the seam is
/// none of them: the seam point is the cyclic order's first ([`arc_span`]).
///
/// ★ **A radical mismatch is not among them**, and nothing here can reach one: every value two
/// corners are combined into goes through the line step's orientation above, which asks
/// `QuadVal` only when the two share a spelling. `QuadVal`'s same-radical arithmetic
/// `debug_assert`s on a mix, so a mix would be a panic in a test build, not an abstention.
pub(crate) fn point_in_mixed_ring(
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
    // A corner's chart coordinates `r₀ + r₁·√c`, with the radical read **once, off the root** — a
    // coordinate's own `QuadVal` folds `c` to 0 where its radical part vanishes, so the two
    // coordinates of one corner can spell different `c()`s, and the radical read off either one
    // would drop the other's `√c` term.
    let corner = |nd: NodeId| -> Option<ChartCorner> {
        if let Some((_, cyl, _)) = pierce_name(nd) {
            let def = def_of(cyl)?;
            let (line, s) = pierce_meet(jd, cyl, def, nd)?.into_exact()?;
            let (b, d) = (line.base(), line.dir());
            let coord = |e: &[Rat; 3]| -> Option<(Rat, Rat)> {
                let along = dot(&d, e)?;
                Some((
                    dot(&b, e)?.checked_add(s.a().checked_mul(along)?)?,
                    s.b().checked_mul(along)?,
                ))
            };
            let ((x0, x1), (y0, y1)) = (coord(e1)?, coord(e2)?);
            Some(ChartCorner {
                r0: [x0, y0],
                r1: [x1, y1],
                c: s.c(),
            })
        } else {
            let p = node_coords_rat(jd, nd)?;
            let pr = chart.project(&p)?;
            let zero = Rat::from_int(0);
            Some(ChartCorner {
                r0: pr,
                r1: [zero, zero],
                c: zero,
            })
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
                let (ca, cb) = (corner(na)?, corner(nb)?);
                let (a, b) = (ca.values()?, cb.values()?);
                let ya = a[1].checked_sub(&QuadVal::from_rat(qy))?.sign();
                let yb = b[1].checked_sub(&QuadVal::from_rat(qy))?.sign();
                // A corner on the ray is the boundary only when the probe *is* that corner —
                // corner against the rational probe, never corner against corner. Every corner
                // is `a` of exactly one step, so this asks each corner once.
                if ya == Orient::Zero
                    && a[0].checked_sub(&QuadVal::from_rat(qx))?.sign() == Orient::Zero
                {
                    return None;
                }
                // A step lying along the ray straddles nothing (half-open: neither end is
                // above); it is the boundary iff the probe lies between its ends.
                if ya == Orient::Zero && yb == Orient::Zero {
                    let sa = a[0].checked_sub(&QuadVal::from_rat(qx))?.sign();
                    let sb = b[0].checked_sub(&QuadVal::from_rat(qx))?.sign();
                    if sa != sb {
                        return None;
                    }
                    continue;
                }
                // The half-open rule (`ray_step_crossing`), read lazily: only a straddling step
                // pays for `orient2d`.
                if ray_straddle(ya, yb).is_none() {
                    continue;
                }
                // orient2d(a, b, q) = (b−a) × (q−a). In one radical — the two corners spell the
                // same `c`, or one of them is rational — it is `QuadVal` arithmetic. Two spellings
                // (the corners of one ruling named by two plane pairs carry `c` and `k²c`) or two
                // circles' radicals take the form every two-radical sign in the kernel takes: the
                // same value is `(a−q) × (b−q)`, bilinear in `a − q ∈ ℚ(√u)²` and `b − q ∈ ℚ(√v)²`,
                // so its four terms go to `biquad_sign`.
                let o = if ca.radical() == cb.radical()
                    || ca.radical() == Rat::from_int(0)
                    || cb.radical() == Rat::from_int(0)
                {
                    let (qxv, qyv) = (QuadVal::from_rat(qx), QuadVal::from_rat(qy));
                    b[0].checked_sub(&a[0])?
                        .checked_mul(&qyv.checked_sub(&a[1])?)?
                        .checked_sub(
                            &b[1]
                                .checked_sub(&a[1])?
                                .checked_mul(&qxv.checked_sub(&a[0])?)?,
                        )?
                        .sign()
                } else {
                    let ra = [ca.r0[0].checked_sub(qx)?, ca.r0[1].checked_sub(qy)?];
                    let rb = [cb.r0[0].checked_sub(qx)?, cb.r0[1].checked_sub(qy)?];
                    let cross = |x: &[Rat; 2], y: &[Rat; 2]| -> Option<Rat> {
                        x[0].checked_mul(y[1])?.checked_sub(x[1].checked_mul(y[0])?)
                    };
                    nacre_exact::biquad_sign(
                        cross(&ra, &rb)?,
                        cross(&ca.r1, &rb)?,
                        cross(&ra, &cb.r1)?,
                        cross(&ca.r1, &cb.r1)?,
                        ca.c,
                        cb.c,
                    )?
                };
                match ray_step_crossing(ya, yb, o) {
                    Some(true) => inside = !inside,
                    Some(false) => {}
                    None => {
                        // probe on the step
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
                        return None;
                    }
                };
                let Some((line, s)) = roots else { continue };
                // The arc's CCW span: the step's ends oriented by the carried `ccw` bit — the
                // same convention every arc consumer reads (membership is direction-agnostic,
                // so the *set* is what the CCW pair names).
                let (lo_nd, hi_nd) = if arc.ccw { (na, nb) } else { (nb, na) };
                let e_lo = pierce_meet(jd, arc.cyl, &arc.def, lo_nd)?.into_exact()?;
                let e_hi = pierce_meet(jd, arc.cyl, &arc.def, hi_nd)?.into_exact()?;
                let ends = if lo_nd == hi_nd {
                    ArcEnds::Closed(&e_lo)
                } else {
                    ArcEnds::Open(&e_lo, &e_hi)
                };
                for root in s {
                    // Right of the probe along the ray: chart-x of the root.
                    let (bse, dir) = (line.base(), line.dir());
                    let x = QuadVal::from_rat(dot(&bse, e1)?)
                        .checked_add(&root.checked_mul_rat(dot(&dir, e1)?)?)?;
                    let xsign = x.checked_sub(&QuadVal::from_rat(qx))?.sign();
                    let rootp = (line.clone(), root);
                    match xsign {
                        // ★ **A root at the probe is the probe itself** (it is on the ray's
                        // plane and the class's), so the probe is on this circle — and on this
                        // **ring's boundary** only where the arc is: within its closed span. Off
                        // the span it is another part of the circle, and the root, not on the
                        // open ray, counts nothing. The lateral road reads the same rule
                        // (`loop_parity`'s arc at the point's own height), and so do the line
                        // steps (a corner is the boundary only when the probe *is* that corner).
                        // The span is asked here alone: a root left of the probe never needs it.
                        Orient::Zero => {
                            if arc_span(&arc.def, ends, &rootp) == Some(ArcSpan::Outside) {
                                continue;
                            }
                            return None;
                        }
                        Orient::Negative => continue,
                        Orient::Positive => {}
                    }
                    // Inside the CCW span end[0] → end[1]? Cyclic, with the seam point first
                    // (`arc_span`): an end or a root on the seam is ordered like any other.
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
                        let rs = ruling_side(&ray_plane, &arc.def, (&end.0, &end.1))?;
                        let ty = if (if ccw { -rs } else { rs }) > 0 {
                            Orient::Positive
                        } else {
                            Orient::Negative
                        };
                        #[cfg(test)]
                        tie_probe::arc_end_decided();
                        Some(ray_step_crossing(Orient::Zero, ty, ty) == Some(true))
                    };
                    // ★★★★★ **`Equal` is not "outside the span"** — read so, a root on the arc's
                    // own end counts nothing, silently: a confident wrong answer, not an
                    // abstention. Both ends decide it alike (`departs_across`), and the span
                    // itself is one spelling (`arc_span`) the lateral road reads too.
                    // A closed arc's node is both its ends, and each end counts as the step it
                    // is: the arc leaving counter-clockwise and the arc arriving, walked back.
                    let contained = match arc_span(&arc.def, ends, &rootp)? {
                        ArcSpan::Inside => true,
                        ArcSpan::Outside => false,
                        ArcSpan::AtLo => departs_across(&e_lo, true)?,
                        ArcSpan::AtHi => departs_across(&e_hi, false)?,
                        ArcSpan::AtBothEnds => {
                            departs_across(&e_lo, true)? != departs_across(&e_hi, false)?
                        }
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
