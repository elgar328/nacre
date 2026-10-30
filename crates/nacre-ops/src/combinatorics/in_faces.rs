use super::*;
use nacre_geom::intersect::RingSide;
/// **Is the point `p` inside this component**, asked along the ray `p + t·dir` — the coordinate
/// twin of [`point_in_component`].
///
/// **Two roads, one crossing rule.** The two differ only in how the *point* arrives: a
/// three-plane name (exact without coordinates, so it survives a rotated class) or rational
/// coordinates (the only thing a curved component can offer, since none of its faces carries a
/// vertex). Everything they ask about a **face** — a circle's radial side, a band's roots and
/// axial span, the abandon-on-boundary policy — is the same rule and is called, not restated.
/// Restating it is how the two would come to disagree about a graze.
///
/// **What is measured, and what is owed.** `an_enclosed_cylindrical_void_is_a_cavity` is the
/// fixture that makes this road say `true`, and it goes red if the road is stubbed to `false` —
/// two *disjoint* bodies would not, since their answer is "outside" whatever the road does.
///
/// The **curved** arm is barely loaded from here: this road is entered a handful of times in
/// the suite and mostly looks at no cylinder face at all (the other component is all planes —
/// the void fixture's is a box); where it does (`two_cylinders_with_coplanar_caps_fuse_apart`)
/// the ray misses. Nor can a fixture with parallel axes do better: the gate refuses any boolean
/// whose two cylinders' faces it cannot prove apart (spans and arcs), so a ray from one
/// cap's centre crosses the other lateral **0 or 2 times** and the parity is the same. What loads the lateral
/// road ([`lateral_face_crossings`]) is the **named** road's probes — the crossing
/// census's corner Commons, whose vertex rays cross the other half's panel — and the lattice
/// oracle on the through-boss; a real `k = 1` from a coordinate probe wants a ∥ wall inside the
/// strip — a population the gate **serves** rather than refuses (a crossing or a
/// tangency), so what owes this arm a fixture is the rulings road, not a refusal.
///
/// [`Said::In`] or [`Said::Out`] decide; [`Said::Abstain`] is a ray that grazed, and
/// [`Said::Silent`] one stopped by a class the gate passed silent — a face on such a class it could
/// not show it misses, or a crossing a ring that class names could not show outside by its box.
/// In both the caller has other directions to try.
pub(crate) fn point_in_faces_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: &[nacre_exact::Rat; 3],
    dir: &[nacre_exact::Rat; 3],
    faces: &[CompFace],
) -> Result<Said, BoolError> {
    use nacre_exact::Rat;
    let not_rational = || reject(RejectReason::WitnessNotRational);
    let zero = Rat::from_int(0);
    if dir.iter().all(|c| *c == zero) {
        return Err(reject(RejectReason::DegenerateWitness));
    }
    // The half this road counts is **ahead** of the origin, so its cut plane's negative side is
    // `dir·(x − p) > 0` — the mirror of the named road, which counts behind. Either half has the
    // same parity; what matters is that one road picks one.
    let ahead = (|| {
        let n = [
            zero.checked_sub(dir[0])?,
            zero.checked_sub(dir[1])?,
            zero.checked_sub(dir[2])?,
        ];
        Some([
            n[0],
            n[1],
            n[2],
            zero.checked_sub(nacre_exact::dot3_rat(&n, p)?)?,
        ])
    })()
    .ok_or_else(not_rational)?;
    let line = planes_through_line(p, dir).ok_or_else(not_rational)?;

    let mut count = 0usize;
    for f in faces {
        let q = match &f.surf {
            CompSurf::Cylinder(def) => {
                // A circle or polygon outer on a cylinder face has no producer; refusing to
                // guess costs the caller another direction, never a wrong answer.
                let BoundEdges::Lateral(loops) = &f.outer else {
                    return Ok(Said::Abstain);
                };
                match lateral_face_crossings(jd, [&line[0], &line[1]], def, loops, &ahead) {
                    Some(CurvedHit::Counted(k)) => {
                        count += k;
                        continue;
                    }
                    Some(CurvedHit::Graze) | None => return Ok(Said::Abstain),
                }
            }
            CompSurf::Plane(q) => *q,
        };
        // **The face is read by its cause.** A class with no world statement — one the gate proved
        // silent about the cylinders — has no rational crossing with the ray, so the one question
        // its interval answers is whether the ray misses the face: every corner strictly to one
        // side of a plane through the ray, or strictly behind its origin ([`ray_misses_face`]).
        // What that cannot show — the ray may well cross the face — this ray cannot count, and
        // another direction may miss it. A face on a stated plane whose outline a silent class
        // names has an exact crossing, and is read like any stated face below; only the ring that
        // class names answers by its box, inside the one coordinate reading of a ring.
        let Some(coeffs) = jd.world_rat(q) else {
            let unstated = crate::planes::silent_among(cyls, &[q])
                || std::iter::once(&f.outer)
                    .chain(&f.inner)
                    .any(|b| matches!(b, BoundEdges::Ring(r) if names_a_silent_class(cyls, r)));
            if !unstated {
                return Err(not_rational());
            }
            if ray_misses_face(jd, cyls, q, &line, &ahead, f) {
                continue;
            }
            return Ok(Said::Silent);
        };
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        let nd = nacre_exact::dot3_rat(&n, dir).ok_or_else(not_rational)?;
        let residual = nacre_exact::dot3_rat(&n, p)
            .and_then(|v| v.checked_add(coeffs[3]))
            .ok_or_else(not_rational)?;
        // **Parallel is the origin's question, not the ray's.** `nd == 0` and a nonzero
        // residual is a plain miss. `nd == 0` with a zero residual is the ray lying **in** this
        // plane, where it never passes from one side of the face to the other — zero crossings —
        // *unless* the origin is on the face itself, and then "is the origin inside the component"
        // has no answer. That is exactly the `t == 0` policy below, so the two are one arm: a
        // coplanar ray is the origin's own crossing.
        let t = if nd == zero {
            if residual != zero {
                continue;
            }
            zero
        } else {
            (|| -> Option<Rat> {
                zero.checked_sub(residual)?
                    .checked_mul(Rat::new(nd.denom(), nd.numer())?)
            })()
            .ok_or_else(not_rational)?
        };
        if t < zero {
            continue; // behind the origin
        }
        let x = (|| -> Option<[Rat; 3]> {
            let mut x = *p;
            for k in 0..3 {
                x[k] = x[k].checked_add(t.checked_mul(dir[k])?)?;
            }
            Some(x)
        })()
        .ok_or_else(not_rational)?;
        // The face is read by the one coordinate reading of a face ([`rational_point_in_face`]); a
        // ring a silent class names answers by its box there, and what the box cannot show leaves
        // this ray to the gate's undecidedness. A value that reading could not form is raised here,
        // not abstained on — the coordinate ray's rule.
        let material = match rational_point_in_face(jd, cyls, q, &x, f)? {
            Said::In => true,
            Said::Out => false,
            Said::Unformed => return Err(not_rational()),
            Said::Silent => return Ok(Said::Silent),
            Said::Abstain | Said::Declines => return Ok(Said::Abstain),
        };
        if t == zero {
            // The crossing is the ray's own origin — either the ray meets this plane there, or it
            // lies in it. On the face's material the query sits on the component's boundary, where
            // "inside" has no answer; off it, the plane is touched (or run along) at points
            // outside the face and counts for nothing.
            if material {
                return Ok(Said::Abstain);
            }
            continue;
        }
        if material {
            count += 1;
        }
    }
    Ok(if count % 2 == 1 { Said::In } else { Said::Out })
}

/// **Does the ray miss this face, provably?** — the points whose box holds the face's outline, in
/// the judge's interval arithmetic ([`ring_box_points_hp`]), all strictly on one side of either
/// plane through the ray (`line`: the ray's line misses the outline's hull), or all strictly on the
/// positive side of `ahead` (behind the origin, the half this road does not count). `false` — not
/// shown, not "crossed" — for an outline that is not a ring or a point the precision could not
/// separate. `q` is the face's class.
fn ray_misses_face(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    line: &[[nacre_exact::Rat; 4]; 2],
    ahead: &[nacre_exact::Rat; 4],
    f: &CompFace,
) -> bool {
    let BoundEdges::Ring(r) = &f.outer else {
        return false;
    };
    let Some(corners) = ring_box_points_hp(jd, cyls, q, r) else {
        return false;
    };
    let prec = jd.standard.prec;
    let one_side = |plane: &[nacre_exact::Rat; 4]| {
        let signs: Vec<Option<bool>> = corners
            .iter()
            .map(|x| plane_sign_hp(plane, x, prec))
            .collect();
        signs.iter().all(|s| *s == Some(true)) || signs.iter().all(|s| *s == Some(false))
    };
    one_side(&line[0])
        || one_side(&line[1])
        || corners
            .iter()
            .all(|x| plane_sign_hp(ahead, x, prec) == Some(true))
}

/// **Where `p` is against this circle** — inside its disk, outside, or on the circle.
///
/// It takes **coordinates**, so both roads ask it: the named probe realizes its three-plane
/// crossing first, the coordinate road already has one. The rule must not be written twice.
///
/// The rule is `cylinder_radial_side`'s, the one the nesting engine reads for a disk target —
/// a circle bound is `cylinder ∩ plane`, so "inside the disk" is "inside the cylinder's radius".
pub(crate) fn point_in_disk(p: &[nacre_exact::Rat; 3], def: &nacre_topo::CylinderDef) -> RingSide {
    match nacre_exact::cylinder_radial_side(p, &def.origin(), &def.dir(), def.r2()) {
        nacre_exact::Orient::Negative => RingSide::Inside,
        nacre_exact::Orient::Positive => RingSide::Outside,
        nacre_exact::Orient::Zero => RingSide::OnBoundary,
    }
}

/// **Is the rational point `p` on class `wc` inside this ring?** — the one exact reading of a
/// coordinate against a ring. The 2-D coordinate witness (`nesting`) asks it of a cell, the one
/// coordinate reading of a face ([`rational_point_in_face`]) of each of the face's rings, and the
/// named road ([`named_point_in_face`]) of a mixed ring at the name's rational coordinates.
///
/// The answer is what one witness says ([`Said`]): `In` or `Out`; [`Said::Abstain`] where this
/// point cannot rank itself — on the ring, or where the mixed walk abstains; [`Said::Silent`] for a
/// ring whose corners name a class the gate passed silent, when the ring's box does not show the
/// point outside; [`Said::Unformed`] for a value that could not be formed — a class with no rational
/// description, or a corner or the chart past `Rat`. It never declines: a coordinate reads every
/// ring kind. `Err` is a producer defect, a straight ring of fewer than three edges.
///
/// **The order is the rule; the cause that knows itself speaks first.**
/// 1. A ring a silent class names has no rational corners, for the gate's own reason, so its box
///    answers ([`ring_box_says`]) before anything reads a corner.
/// 2. A straight ring of fewer than three edges is no boundary: the polygon walk would answer
///    `Outside` by contract and make it invisible. A mixed ring that short is a lens — a chord and
///    its arc, a half-disc cap — and is read below: asked of every ring, this refuses six census
///    cells for being what they are, a boss straddling a bore and an annulus profile in each of the
///    three operations (measured).
/// 3. A three-plane corner with no rational meet is no circumstance of this point: every other
///    witness of the ring meets it too, so it is `Unformed` and not an abstention of the mixed walk
///    — which is why it is asked before the fork. The chart is formed here, once per ring, so its
///    overflow has one name on both roads.
/// 4. A ring with a pierce corner or an arc step takes [`point_in_mixed_ring`].
/// 5. Every other ring is a polygon of rational corners in the class's chart
///    ([`nacre_geom::intersect::point_in_ring_2d_rat`] — half-open in y, so a corner on the ray is
///    decided rather than abstained on).
pub(crate) fn rational_point_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    p: &[nacre_exact::Rat; 3],
    ring: &[RingEdge],
) -> Result<Said, BoolError> {
    if names_a_silent_class(cyls, ring) {
        let prec = jd.standard.prec;
        let ph = p.map(|v| nacre_exact::HpBounded::of_rat(v, prec));
        return Ok(ring_box_says(jd, cyls, wc, &ph, ring));
    }
    let mixed = ring_is_mixed(ring);
    if !mixed && ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    let Some(coeffs) = jd.world_rat(wc) else {
        return Ok(Said::Unformed);
    };
    let Some((chart, p2)) = Chart2dRat::of_normal(&[coeffs[0], coeffs[1], coeffs[2]])
        .and_then(|chart| chart.project(p).map(|p2| (chart, p2)))
    else {
        return Ok(Said::Unformed);
    };
    if ring
        .iter()
        .any(|e| three_plane_name(e.node).is_some() && node_coords_rat(jd, e.node).is_none())
    {
        return Ok(Said::Unformed);
    }
    if mixed {
        #[expect(
            clippy::disallowed_methods,
            reason = "the one coordinate reading of a ring: its mixed arm"
        )]
        let inside = point_in_mixed_ring(jd, cyls, &coeffs, p, ring);
        return Ok(match inside {
            Some(true) => Said::In,
            Some(false) => Said::Out,
            None => Said::Abstain,
        });
    }
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let Some(ring2) = chart.ring(jd, &nodes) else {
        return Ok(Said::Unformed);
    };
    Ok(
        match nacre_geom::intersect::point_in_ring_2d_rat(p2, &ring2) {
            RingSide::Inside => Said::In,
            RingSide::Outside => Said::Out,
            RingSide::OnBoundary => Said::Abstain,
        },
    )
}

/// **Where the rational point `x` on class `q` is against the plane face `f`** — the one coordinate
/// reading of a face. Each bound answers by its kind — a ring by the one coordinate reading of a
/// ring ([`rational_point_in_ring`]), a circle by its radial side ([`point_in_disk`]) — and
/// [`material_of`] combines them. The 3-D coordinate ray ([`point_in_faces_rat`]) asks it of the
/// face its crossing lands on, and the coordinate probe supply ([`coord_probes`]) of every
/// candidate point, taking only [`Said::In`].
///
/// [`Said::In`] or [`Said::Out`] is the face's material. Otherwise the answer is what the first
/// bound that could not tell said — [`Said::Abstain`] for a point on a bound or a lateral loop
/// (loops bound a cylinder, never a plane), and a ring's own `Abstain`, `Silent` or `Unformed` as
/// it said it; what each costs is the reader's to decide. `Err` is the ring reading's, a producer
/// defect.
pub(crate) fn rational_point_in_face(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    x: &[nacre_exact::Rat; 3],
    f: &CompFace,
) -> Result<Said, BoolError> {
    let mut why = Said::Abstain;
    let side = material_of(f, |b| {
        Ok(match b {
            BoundEdges::Ring(r) => match rational_point_in_ring(jd, cyls, q, x, r)? {
                Said::In => Some(RingSide::Inside),
                Said::Out => Some(RingSide::Outside),
                said => {
                    why = said;
                    None
                }
            },
            BoundEdges::Circle(def) => Some(point_in_disk(x, def)),
            BoundEdges::Lateral(_) => None,
        })
    })?;
    Ok(match decided(side) {
        Some(true) => Said::In,
        Some(false) => Said::Out,
        None => why,
    })
}

/// **Where the point named `t` is against the disk `def` bounds**, on a plane through both.
///
/// **On a class the gate proved this cylinder beside, outside — exactly.** Every lateral face
/// of the cylinder lies strictly on one side of that class's plane, so its rims and the disks
/// they bound do too, and the point, on the plane, is off the disk whatever else names it. This
/// needs nothing of where the point came from, which is why every reader of a named point takes
/// it. Otherwise the point's rational coordinates answer by their radial side, and with none the
/// answer is unknown.
pub(crate) fn named_point_in_disk(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    t: [usize; 3],
    def: &nacre_topo::CylinderDef,
) -> Option<RingSide> {
    if cyls
        .iter()
        .find(|cy| cy.def == *def)
        .is_some_and(|cy| t.iter().any(|&k| cy.beside(k).is_some()))
    {
        return Some(RingSide::Outside);
    }
    node_coords_rat(jd, NodeId::three_planes(Canon3::three(t))).map(|p| point_in_disk(&p, def))
}

/// [`lateral_face_crossings`] for a probe ray named by three plane classes: it states the ray's
/// two planes and the "behind" cut plane from the class table and hands them over.
///
/// `Err` is the abandon-this-ray answer, and says why: [`Said::Silent`] for a ray on a class the
/// gate passed silent, which has no statement to count with, [`Said::Abstain`] for the rest —
/// a graze, a `Wide` name, arithmetic that cannot answer.
fn curved_count(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    a: usize,
    b: usize,
    c: usize,
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
) -> Result<usize, Said> {
    // **A ray on a plane the gate proved the whole cylinder beside crosses none of it.** The ray
    // lies in both `a` and `b`; if the gate recorded this cylinder beside either
    // ([`crate::planes::WorkingCyl::beside`] — by its boxes or by the axis's distance), the count is
    // zero. A class with no world statement reaches here only that way: nothing below can state
    // it. The record and not a second asking of the boxes: those exist only when the boolean has a
    // nameless class, so a shortcut on them would answer for a named ray depending on an unrelated
    // class, and would miss the axis-distance proof.
    if let Some(cy) = cyls.iter().find(|cy| cy.def == *def)
        && (cy.beside(a).is_some() || cy.beside(b).is_some())
    {
        return Ok(0);
    }
    // A class with no rational description and no such proof: this ray cannot be counted, and the
    // caller casts the next one.
    let Some(((ca, cb), cc)) = jd.world_rat(a).zip(jd.world_rat(b)).zip(jd.world_rat(c)) else {
        return Err(if crate::planes::silent_among(cyls, &[a, b, c]) {
            Said::Silent
        } else {
            Said::Abstain
        });
    };
    // The cut plane through the ray's origin, normal `n_a × n_b` — the very direction
    // `plane_plane_cylinder` gives its meet line, so "negative side" is "behind the query".
    let Some(half) = (|| {
        let d = nacre_exact::cross3_rat(&[ca[0], ca[1], ca[2]], &[cb[0], cb[1], cb[2]])?;
        let at = nacre_exact::three_planes_rat([ca, cb, cc])?;
        let d0 = nacre_exact::Rat::from_int(0).checked_sub(nacre_exact::dot3_rat(&d, &at)?)?;
        Some([d[0], d[1], d[2], d0])
    })() else {
        return Err(Said::Abstain);
    };
    match lateral_face_crossings(jd, [&ca, &cb], def, loops, &half) {
        Some(CurvedHit::Counted(k)) => Ok(k),
        Some(CurvedHit::Graze) | None => Err(Said::Abstain),
    }
}

/// A rim's plane restated with the **axis** as its normal — `[m, −m·(o + t·m)]`, so
/// [`nacre_exact::quad::plane_side`] at a point of the cylinder is the sign of `z − t` along
/// the axis. The loops road builds it for every ⊥ class it meets.
fn rim_plane(def: &nacre_topo::CylinderDef, t: nacre_exact::Rat) -> Option<[nacre_exact::Rat; 4]> {
    let (o, m) = (def.origin(), def.dir());
    let mut d0 = nacre_exact::Rat::from_int(0);
    for k in 0..3 {
        let at = o[k].checked_add(t.checked_mul(m[k])?)?;
        d0 = d0.checked_sub(m[k].checked_mul(at)?)?;
    }
    Some([m[0], m[1], m[2], d0])
}

/// **Is the cylinder point `x` on the lateral face these loops bound?** — the parity of the
/// ray up the axis from `x` against the boundary loops, on the chart `(θ, z)`.
///
/// The rule is the planar rings' half-open rule read on this chart. A whole-circle rim is
/// crossed iff it is above `x`. A ring's **arc** (z = const, a CCW span in θ) is crossed iff
/// it is above `x` and `θ_x` is in its span — with a span end exactly at `θ_x` counted by the
/// end the arc *leaves upward* (`+θ`): `lo` counts, `hi` does not
/// ([`nacre_geom::intersect::ray_step_crossing`]`(Zero, ±, ±)`), and the ruling at that corner
/// runs along the ray and never counts. A **ruling** piece is crossed by nothing; `x` on it is
/// the boundary. Every comparison is one the arrangement already owns: `z` by
/// [`nacre_exact::quad::plane_side`] against [`rim_plane`], `θ` by [`arc_span`], a ruling by
/// its own carrier (`plane_side(plane) == 0 ∧ ruling_side_signed(plane) == side`, against the
/// plane its side was measured against — the predicate `crossing_on_ruling` names stations with).
///
/// `None` = `x` on the boundary (a rim, an arc, a ruling), a zero-span
/// arc inside `arc_span`, a loop the road cannot read (a tilted arc, a plane carrier), or checked
/// arithmetic running out. `z` is asked before `θ`: an arc below `x` needs no order.
pub(crate) fn loop_parity(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
    meet: &nacre_exact::quad::MeetLine,
    s: &nacre_exact::quad::QuadVal,
) -> Option<bool> {
    use nacre_exact::Orient;
    use nacre_exact::quad::plane_side;
    let above =
        |t: nacre_exact::Rat| -> Option<Orient> { Some(plane_side(&rim_plane(def, t)?, meet, s)) };
    let mut crossings = 0usize;
    for lp in loops {
        match lp {
            LateralLoop::Rim(c) => {
                let t = crate::planes::axis_param_of_plane(&jd.world_rat(*c)?, def)?;
                match above(t)? {
                    // `z_x < t`: the rim is above the point, and the ray up the axis crosses it.
                    Orient::Negative => crossings += 1,
                    Orient::Positive => {}
                    Orient::Zero => {
                        return None;
                    }
                }
            }
            LateralLoop::Ring(edges) => {
                for e in edges {
                    match &e.carrier {
                        Carrier::Arc(arc) => {
                            let t = pierce_axis_param(jd, def, e.node)?;
                            let side = above(t)?;
                            if side == Orient::Positive {
                                continue; // the arc is below the point
                            }
                            let (lo_nd, hi_nd) = if arc.ccw {
                                (e.node, e.to)
                            } else {
                                (e.to, e.node)
                            };
                            let e_lo = pierce_meet(jd, arc.cyl, &arc.def, lo_nd)?.into_exact()?;
                            let e_hi = pierce_meet(jd, arc.cyl, &arc.def, hi_nd)?.into_exact()?;
                            let ends = if lo_nd == hi_nd {
                                ArcEnds::Closed(&e_lo)
                            } else {
                                ArcEnds::Open(&e_lo, &e_hi)
                            };
                            let Some(span) = arc_span(&arc.def, ends, &(meet.clone(), *s)) else {
                                // A zero-span arc (two ends at one place) or arithmetic out.
                                return None;
                            };
                            if side == Orient::Zero {
                                // At the arc's own z: on the arc iff within its closed span.
                                if span == ArcSpan::Outside {
                                    continue;
                                }
                                return None;
                            }
                            // Half-open `[lo, hi)`: a closed arc's node is its `lo`, counted once.
                            if matches!(span, ArcSpan::Inside | ArcSpan::AtLo | ArcSpan::AtBothEnds)
                            {
                                crossings += 1;
                            }
                        }
                        Carrier::Ruling(rl) => {
                            let (Some(t0), Some(t1)) = (
                                pierce_axis_param(jd, def, e.node),
                                pierce_axis_param(jd, def, e.to),
                            ) else {
                                return None;
                            };
                            let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
                            // Outside the piece's closed axial range: not on it, whatever θ.
                            if above(lo)? == Orient::Negative || above(hi)? == Orient::Positive {
                                continue;
                            }
                            // The plane the ruling's `side` was measured against, carried — not
                            // the corner name's ∥ class, which can be another plane through the
                            // same line — and the signed side, which answers `0` for a tangent
                            // one; arithmetic out abstains like everywhere here.
                            let w = jd.world_rat(rl.plane)?;
                            if plane_side(&w, meet, s) == Orient::Zero
                                && ruling_side_signed(&w, def, (meet, s))? == rl.side
                            {
                                return None;
                            }
                        }
                        Carrier::Plane { .. } => {
                            return None;
                        }
                    }
                }
            }
        }
    }
    Some(crossings % 2 == 1)
}

/// **Where the line `line[0] ∩ line[1]` crosses one lateral face, and how many of those lie on
/// the half the caller is counting.**
///
/// **One copy, two roads.** The named-point probe and the coordinate-point road both state
/// their ray as *two rational planes*, so both ask this. Writing the arms twice is how the two
/// would come to disagree about a graze. A `half` is a rational plane through the ray's origin,
/// and a crossing counts exactly when it is on that plane's **negative** side — so each road
/// states its own half and nothing here has a front or a back: the named road's normal is the
/// line's own direction `n_a × n_b` (which counts **behind** the query — `plane_plane_cylinder`
/// builds the meet line's `dir` as `cross3(n_a, n_b)`, the very `d` that `order_along`'s
/// `fwd == 1` measures against), the coordinate road's is `−dir` (**ahead** of the origin).
///
/// Each root asks the face by [`loop_parity`] **before** the half — a root off the face is not
/// a crossing whichever side it is on, and a root on the face at the ray's own origin is the
/// query on the other component's surface (`Graze`). Read as loops, not as a band's two rims
/// taken for an axial span (the crossing between them iff on opposite sides of the two rim
/// planes): the two say the same thing on a band — measured identical on 2,180 lattice rays,
/// every root and every graze — and the loops also read a panel, a chain rim or a hole, where the
/// span has nothing to say.
///
/// `None` is checked-`Rat` arithmetic that could not answer — an honest decline, never a guess.
pub(crate) fn lateral_face_crossings(
    jd: &Judge<'_, WorkingPlane>,
    line: [&[nacre_exact::Rat; 4]; 2],
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
    half: &[nacre_exact::Rat; 4],
) -> Option<CurvedHit> {
    use nacre_exact::Orient;
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let (meet, roots): (_, [QuadVal; 2]) =
        match nacre_exact::quad::plane_plane_cylinder(line[0], line[1], &o, &m, r2)? {
            CylinderMeet::Pair { line, s } => (line, s),
            CylinderMeet::Tangent { .. } | CylinderMeet::OnRuling(_) => {
                return Some(CurvedHit::Graze);
            }
            CylinderMeet::AxisParallelMiss(_) | CylinderMeet::Miss(_) => {
                return Some(CurvedHit::Counted(0));
            }
            CylinderMeet::CoincidentPlanes | CylinderMeet::ParallelPlanes => {
                return Some(CurvedHit::Graze);
            }
        };
    let mut count = 0usize;
    for s in &roots {
        match loop_parity(jd, def, loops, &meet, s) {
            // On the boundary, a zero-span arc, an unreadable loop, arithmetic out: this ray
            // cannot count this face — the banded arm's `Graze` on a rim, generalized.
            None => return Some(CurvedHit::Graze),
            Some(false) => continue,
            Some(true) => match nacre_exact::quad::plane_side(half, &meet, s) {
                Orient::Zero => return Some(CurvedHit::Graze),
                Orient::Negative => count += 1,
                Orient::Positive => {}
            },
        }
    }
    Some(CurvedHit::Counted(count))
}

/// **Does the open segment of the line `w[0] ∩ w[1]`, a ruling of `def`, between the points
/// `ends` meet the interior of the lateral face `loops` bounds?** — the self-touch question for
/// an edge two plane faces of one solid share along that solid's own lateral.
///
/// Along a ruling only the face's arcs and rims can change whether a point is on the face (a
/// ruling of the face at that angle makes the whole run its boundary), so the segment is cut at
/// their axis parameters and one point per piece — at the piece's middle height, the root of the
/// secant wall and that height's rim plane that lies on the other wall too — is asked by
/// [`loop_parity`]. A piece inside answers for the edge; a piece on the boundary abstains, and
/// so does the edge unless another piece is inside.
///
/// **What a `true` stands on** is that one point: it lies on the edge and strictly inside the
/// face, which is the proposition whole. So the ends' heights must be exact — they are read
/// where a plane of each end's name crosses **this ruling**, not the axis (a plane naming the end
/// may be slanted, and where it meets the axis says nothing of the end; ☐ no input names one so —
/// over the suite's 1,108 ends the two readings agree, and the census asks none) — while the cut heights,
/// read off the arcs' corners as [`loop_parity`] reads them, only decide which points are asked:
/// wrong there, a touch can be missed, never invented. ☐ No input yet where the cut changes the
/// answer — an arc of the face crossing the edge's line meets the two walls there too, and the
/// arrangement makes that point a vertex, so the edge already ends at it; asking a piece at a time
/// keeps the question «does the open segment meet the interior» rather than leaning on that.
///
/// `None` is a point on the face's boundary, an end with no plane across the ruling, a loop the
/// road cannot read, or arithmetic out — the guard calling this abstains then, it does not refuse.
pub(crate) fn ruling_meets_lateral(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
    w: [usize; 2],
    ends: [NodeId; 2],
) -> Option<bool> {
    use nacre_exact::quad::{CylinderMeet, plane_side};
    use nacre_exact::{Orient, Rat, dot3_rat};
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let walls = [jd.world_rat(w[0])?, jd.world_rat(w[1])?];
    let CylinderMeet::OnRuling(ruling) =
        nacre_exact::quad::plane_plane_cylinder(&walls[0], &walls[1], &o, &m, r2)?
    else {
        return None;
    };
    // The ruling's base, as an axis parameter: `m·(base − o) / m·m`.
    let base = ruling.base();
    let mm = dot3_rat(&m, &m)?;
    let mut d = [Rat::from_int(0); 3];
    for k in 0..3 {
        d[k] = base[k].checked_sub(o[k])?;
    }
    let at_base = dot3_rat(&m, &d)?.checked_mul(Rat::new(mm.denom(), mm.numer())?)?;
    // An end's height: a plane of its name across the ruling meets it at the end.
    let z_of = |n: NodeId| -> Option<Rat> {
        let planes: Vec<usize> = match (three_plane_name(n), pierce_name(n)) {
            (Some(t), _) => t.to_vec(),
            (None, Some((p, _, _))) => p.to_vec(),
            (None, None) => return None,
        };
        planes.into_iter().find_map(|c| {
            let t = nacre_exact::axis_param_of_plane(&jd.world_rat(c)?, &base, &m)?;
            at_base.checked_add(t)
        })
    };
    let (a, b) = (z_of(ends[0])?, z_of(ends[1])?);
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let mut ts = vec![lo, hi];
    for lp in loops {
        match lp {
            LateralLoop::Rim(c) => {
                ts.push(crate::planes::axis_param_of_plane(&jd.world_rat(*c)?, def)?)
            }
            LateralLoop::Ring(edges) => {
                for e in edges {
                    if matches!(e.carrier, Carrier::Arc(_)) {
                        ts.push(pierce_axis_param(jd, def, e.node)?);
                    }
                }
            }
        }
    }
    ts.retain(|t| *t >= lo && *t <= hi);
    ts.sort();
    ts.dedup();
    let half = Rat::new(1, 2)?;
    let mut abstained = false;
    for piece in ts.windows(2) {
        let mid = piece[0].checked_add(piece[1])?.checked_mul(half)?;
        let rim = rim_plane(def, mid)?;
        // A secant wall meets the cylinder in two rulings at this height; the point is the one on
        // the other wall too. (A tangent wall meets it in its own ruling alone — the other wall is
        // then the secant one.)
        let point = (0..2).find_map(|i| {
            match nacre_exact::quad::plane_plane_cylinder(&walls[i], &rim, &o, &m, r2)? {
                CylinderMeet::Pair { line, s } => s
                    .into_iter()
                    .find(|x| plane_side(&walls[1 - i], &line, x) == Orient::Zero)
                    .map(|x| (line, x)),
                _ => None,
            }
        });
        match point.and_then(|(line, x)| loop_parity(jd, def, loops, &line, &x)) {
            Some(true) => return Some(true),
            Some(false) => {}
            None => abstained = true,
        }
    }
    (!abstained).then_some(false)
}

/// A component as its faces, each boundary already carrying its edges' walls.
pub(crate) type ComponentFaces = Vec<CompFace>;

/// **Where the point named `x` is against the plane face `f` on class `q`** — inside it, outside
/// it, or on its boundary; `None` where a bound could not tell. `x` lies on `q` and its name need
/// not carry `q`.
///
/// Each bound answers by its own kind — a polygon by a ring walk, a circle by its radial side, a
/// mixed ring (an arc step, a pierce corner) by the one coordinate reading of a ring
/// ([`rational_point_in_ring`]) at the name's rational coordinates. The ray road's walk cannot read
/// a pierce corner against a class with a world statement (it sides nodes with no cylinder table
/// there) — letting it try would refuse with `PierceVertexUnnamed`. That reading's abstention does
/// not tell a point on its ring from a graze, so a point on a mixed ring is unknown.
///
/// **A name holding a class the gate passed silent** has no rational coordinates, so where a ring
/// bound says nothing it is asked the one question its interval answers — whether it stands
/// outside the ring's box ([`silent_name_by_box`], the rule the 2-D witnesses read).
pub(crate) fn named_point_in_face(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    x: [usize; 3],
    f: &CompFace,
) -> Result<Option<RingSide>, BoolError> {
    let at = || node_coords_rat(jd, NodeId::three_planes(Canon3::three(x)));
    let inside = |b: &BoundEdges| -> Result<Option<RingSide>, BoolError> {
        match b {
            BoundEdges::Ring(r) => {
                if ring_is_mixed(r) {
                    let Some(px) = at() else {
                        return Ok(None); // no rational coordinate: abstain
                    };
                    return Ok(match rational_point_in_ring(jd, cyls, q, &px, r)? {
                        Said::In => Some(RingSide::Inside),
                        Said::Out => Some(RingSide::Outside),
                        _ => None,
                    });
                }
                if point_on_ring(jd, q, x, r)? {
                    return Ok(Some(RingSide::OnBoundary));
                }
                Ok(
                    every_ray(jd, cyls, q, NodeId::three_planes(Canon3::three(x)), r)?
                        .first()
                        .map(|&b| {
                            if b {
                                RingSide::Inside
                            } else {
                                RingSide::Outside
                            }
                        }),
                )
            }
            BoundEdges::Circle(def) => Ok(named_point_in_disk(jd, cyls, x, def)),
            // Loops bound a cylinder, never a plane — a producer error, not an input.
            BoundEdges::Lateral(_) => Ok(None),
        }
    };
    let silent = crate::planes::silent_among(cyls, &x);
    material_of(f, |b| match (inside(b)?, b) {
        (None, BoundEdges::Ring(r)) if silent => Ok(
            match silent_name_by_box(jd, cyls, q, NodeId::three_planes(Canon3::three(x)), r) {
                Said::Out => Some(RingSide::Outside),
                _ => None,
            },
        ),
        (side, _) => Ok(side),
    })
}

pub(crate) fn point_in_component(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    query: [usize; 3],
    faces: &[CompFace],
) -> Result<Said, BoolError> {
    let mut vplanes = query.to_vec();
    vplanes.sort_unstable();
    vplanes.dedup();
    // A face reading that said nothing at a point whose name holds a class the gate passed silent
    // is the gate's undecidedness; anything else that stops a ray is a graze.
    let blocked = |x: &[usize]| {
        if crate::planes::silent_among(cyls, x) {
            Said::Silent
        } else {
            Said::Abstain
        }
    };

    // One ray attempt along `L = a ∩ b`, located by the query point `{a,b,c}`. `In`/`Out` decide;
    // anything else is why the line was given up, and the caller tries the next plane pair.
    let attempt = |a: usize, b: usize, c: usize| -> Result<Said, BoolError> {
        let mut count = 0usize;
        for f in faces {
            // A cylindrical face is counted by its own arm — the crossings are roots of a
            // quadratic, not three-plane points, and "inside the face" is an axial span when
            // the face states one, a bare miss-oracle when it does not.
            if let CompSurf::Cylinder(def) = &f.surf {
                // A circle or polygon outer on a cylinder face has no producer; refusing
                // to guess costs the caller another node, never a wrong answer.
                let BoundEdges::Lateral(loops) = &f.outer else {
                    return Ok(Said::Abstain);
                };
                match curved_count(jd, cyls, a, b, c, def, loops) {
                    Ok(k) => {
                        count += k;
                        continue;
                    }
                    Err(said) => return Ok(said),
                }
            }
            let CompSurf::Plane(q) = &f.surf else {
                return Ok(Said::Abstain);
            };
            let q = *q;
            // Inside `q`'s material at the point `x`: inside the outer bound, outside every hole.
            // On the face's boundary is no answer here — the same non-generic abandon
            // `point_on_ring` has always raised.
            let material = |x: [usize; 3]| named_point_in_face(jd, cyls, q, x, f).map(decided);
            if jd.plane_pair_dir_sign(a, b, q) == 0 {
                // `L` is parallel to `q` — and possibly **in** it, which is not the same thing.
                // A coplanar ray never passes from one side of this face to the other, so zero
                // crossings is the right count. What a bare `continue` would also swallow is the
                // *query*: if it lies on this face, "is the query inside the component" has no
                // answer at all, and skipping the face would answer it anyway. That is the same
                // proposition the `fwd == 0` arm below abandons for a transversal plane — one rule,
                // both spellings.
                //
                // A query on the face's material is a body standing on another's face: a
                // corner on it, or an edge resting along it
                // (`two_bodies_touching_at_a_corner_or_an_edge_build_in_either_order`). Such a
                // query's name need not carry `q` — it is the other body's corner — so `material`
                // reads it with a name that misses `q`, and the ray orders it by a plane of that
                // name that cuts the line (`pin_on_line`).
                let mut vq = [query[0], query[1], query[2]];
                vq.sort_unstable();
                #[expect(
                    clippy::disallowed_methods,
                    reason = "a point question: does the query lie on the face's plane `q`"
                )]
                let on_plane =
                    side_of(jd, &[], NodeId::three_planes(Canon3::three(vq)), q) == Some(0);
                if on_plane {
                    match material(vq)? {
                        Some(false) => {}
                        Some(true) => return Ok(Said::Abstain),
                        None => return Ok(blocked(&vq)),
                    }
                }
                continue;
            }
            let mut x = [a, b, q];
            x.sort_unstable();
            let Some(in_g) = material(x)? else {
                return Ok(blocked(&x));
            };
            let fwd = order_along(jd, a, b, c, q);
            if fwd == 0 {
                // The crossing is the ray origin itself (query on plane `q`). Inside `q`'s
                // material ⇒ query on the component's surface ⇒ undecidable → abandon.
                if in_g {
                    return Ok(Said::Abstain);
                }
                continue;
            }
            if fwd == 1 && in_g {
                count += 1; // one side of the line — parity is the same on either half
            }
        }
        Ok(if count % 2 == 1 { Said::In } else { Said::Out })
    };

    // Every plane pair of this query blocked: the node abstains, and says [`Said::Silent`] if any
    // pair was stopped by the gate's undecidedness — the name the caller's list runs out by.
    let mut said = Said::Abstain;
    for i in 0..vplanes.len() {
        for j in (i + 1)..vplanes.len() {
            let (a, b) = (vplanes[i], vplanes[j]);
            // Locator `c`: a plane of the query off the line `a ∩ b`, so `{a,b,c}` is the query.
            let Some(c) = pin_on_line(jd, a, b, &vplanes) else {
                continue;
            };
            match attempt(a, b, c)? {
                decided @ (Said::In | Said::Out) => return Ok(decided),
                Said::Silent => said = Said::Silent,
                _ => {}
            }
        }
    }
    Ok(said)
}
