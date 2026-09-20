use super::*;
/// **Where a segment crosses a cylinder's rulings on a ∥ class** — the straight sibling of
/// [`circle_crossings`], same machinery ([`nacre_exact::quad::plane_plane_cylinder`], the same
/// fences, the same names), different reachable arms: on a class **parallel** to the axis the
/// meet can be tangent to the lateral — **skipped**, because a touch divides nothing and cutting
/// there would make a zero-length piece — or lie on it, which the ladder refuses by name rather
/// than reasoning it unreachable.
///
/// Returns each crossing with the exact point it was solved at (`(line, s)` — the side and
/// extent tests below read it without re-solving). `None` is overflow.
fn lateral_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    w: &[Rat; 4],
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    sg: &MergedSeg,
) -> Result<Option<Vec<LateralCrossing>>, BoolError> {
    use nacre_exact::quad::CylinderMeet;
    use nacre_topo::QuadRoot;
    let Some(v) = combinatorics::class_coeffs_rat(jd, sg.wall) else {
        return Ok(None);
    };
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let (line, roots) = match nacre_exact::quad::plane_plane_cylinder(w, &v, &o, &m, r2) {
        Some(CylinderMeet::Pair { line, s }) => {
            (line, vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])])
        }
        // A tangent line touches without separating; cutting there would make a zero-length
        // piece — the same skip the circle side's `Double` arm makes.
        // ★ A **tangent** class touches the cylinder along one ruling: the
        // segment's line meets it in one point, root `Double` — the crossing the seated face's own
        // tangent piece (`side == 0`) is cut at, named as the scan names it.
        Some(CylinderMeet::Tangent { line, s }) => (
            line,
            vec![(QuadRoot::Double, nacre_exact::quad::QuadVal::from_rat(s))],
        ),
        Some(CylinderMeet::Miss(_)) | Some(CylinderMeet::AxisParallelMiss(_)) => {
            return Ok(Some(Vec::new()));
        }
        // The segment's line lies on the lateral (a wall meeting the cylinder exactly along a
        // ruling), or the planes degenerate: shapes this ladder does not arrange yet.
        Some(CylinderMeet::OnRuling(_))
        | Some(CylinderMeet::CoincidentPlanes)
        | Some(CylinderMeet::ParallelPlanes) => {
            return Err(reject(RejectReason::RulingBoundNotYet));
        }
        None => return Ok(None),
    };
    // ★ Whether a crossing is **on this segment** is the caller's question now, asked in the one
    // vocabulary that answers it for both kinds of end — the same move the circle side makes.
    Ok(Some(
        roots
            .into_iter()
            .map(|(root, s)| {
                (
                    combinatorics::NodeId::pierce(wc, sg.wall, cyl, root),
                    line.clone(),
                    s,
                )
            })
            .collect(),
    ))
}

/// One crossing [`lateral_crossings`] found: its Pierce name, and the exact `(line, s)` it was
/// solved at (the side and extent tests read it without re-solving).
type LateralCrossing = (
    combinatorics::NodeId,
    nacre_exact::quad::MeetLine,
    nacre_exact::quad::QuadVal,
);

/// `split_rulings`' product — `None` when nothing crossed.
type SplitRulings = Option<(Vec<MergedSeg>, Vec<MergedRuling>)>;

/// **Cut every ruling a segment crosses into pieces, and the segments with it** — the straight
/// twin of [`split_circles`], run after it.
///
/// ★★ **It cuts rulings against *segments*, and deliberately not against arcs** — as its twin cuts
/// circles against segments and not against rulings. One class *can* carry both populations
/// (there is no blanket refusal that says otherwise), so what keeps the pair of splits
/// complete is not «they never share a class» but the sharper fact that **a class's edge is
/// some face's boundary**: two edges of different cylinders meeting would put two lateral faces on
/// one point, which `crate::planes::lateral_faces_clear` denies. That one sentence covers
/// every pair no split cuts — circle×circle and ruling×ruling as well as this one — and
/// `ClassEdges::of` carries a shipped backstop over the pair that has an exact predicate
/// ([`crate::RejectReason::CircleCrossesRuling`]).
///
/// The ordering vocabulary is the established one: along a segment, [`cmp_along`] on the
/// canonical meet line; along a ruling, the same comparator on the **axis coordinate**
/// (`m · x`, evaluated exactly from each node's re-solved name).
pub(super) fn split_rulings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    rulings: &[MergedRuling],
    aliases: &Aliases,
) -> Result<SplitRulings, BoolError> {
    use nacre_exact::quad::QuadVal;
    let undecided = || reject(RejectReason::WitnessNotRational);
    if rulings.is_empty() {
        return Ok(None);
    }
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
    let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    // A point's axis coordinate, from the exact `(line, s)` it was solved at.
    let coord_at = |def: &nacre_topo::CylinderDef,
                    line: &nacre_exact::quad::MeetLine,
                    s: &QuadVal|
     -> Option<QuadVal> {
        let m = def.dir();
        let (b, d) = (line.base(), line.dir());
        QuadVal::from_rat(dot(&b, &m)?).checked_add(&s.checked_mul_rat(dot(&d, &m)?)?)
    };
    // Each ruling's end coordinates, re-solved from the names once.
    let ruling_ends: Vec<[QuadVal; 2]> = rulings
        .iter()
        .map(|r| {
            let coord = |n| {
                let (line, s) = combinatorics::pierce_meet(jd, r.cyl, &r.def, n)?;
                coord_at(&r.def, &line, &s)
            };
            match (coord(r.end[0]), coord(r.end[1])) {
                (Some(a), Some(b)) => Ok([a, b]),
                _ => Err(undecided()),
            }
        })
        .collect::<Result<_, BoolError>>()?;
    // Distinct cylinders carrying rulings here, each solved once per segment.
    let mut cyl_list: Vec<(usize, &nacre_topo::CylinderDef)> = Vec::new();
    for r in rulings {
        if !cyl_list.iter().any(|&(c, _)| c == r.cyl) {
            cyl_list.push((r.cyl, &r.def));
        }
    }
    let mut on_ruling: Vec<Vec<combinatorics::NodeId>> = vec![Vec::new(); rulings.len()];
    let mut on_seg: Vec<Vec<combinatorics::NodeId>> = vec![Vec::new(); segs.len()];
    for (si, sg) in segs.iter().enumerate() {
        // ★★★★★ **The skip is gone with the reason it gave.** It said a segment whose end a
        // cylinder pinned "carries no fence planes to test extent with" — true of the fences, which
        // are retired; the extent question is `closed_contains` now and it answers for both kinds of
        // end. The rest of that sentence ("such a segment shares only endpoints with a ruling") was
        // an argument about one narrow population, never a measurement, and it is exactly the sort
        // of claim that goes quietly false when the population widens. So it is asked instead of
        // assumed.
        // ☑ Measured: **200** segments now reach here that the skip dropped, and **392** crossings
        // on them are examined. The census is bit-identical, so the skip's *conclusion* held — none
        // of those crossings splits a ruling. It is a measurement now rather than an argument.
        //
        // A segment whose two ends are one point cuts nothing — asked by **name**, the identity,
        // rather than by a coordinate a pierce end does not have. ☑ Measured unexercised.
        if sg.end[0] == sg.end[1] {
            continue;
        }
        // Both ends on this segment's own line, once — the pair `closed_contains` is asked with.
        let mut it = (0..2).map(|k| {
            combinatorics::on_line(
                jd,
                cyls,
                wc,
                sg.wall,
                Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
            )
            .ok_or_else(undecided)
        });
        let seg_ends = [it.next().unwrap()?, it.next().unwrap()?];
        for &(cyl, def) in &cyl_list {
            let Some(xs) = lateral_crossings(jd, wc, &w, cyl, def, sg)? else {
                return Err(undecided());
            };
            for (n, line, s) in xs {
                // ★ On this segment? The caller's question since the fences went — same rule, same
                // vocabulary as the circle side.
                let at =
                    combinatorics::locate(jd, cyls, wc, sg.wall, combinatorics::PointOn::Pierce(n))
                        .ok_or_else(undecided)?;
                if !combinatorics::closed_contains(jd, wc, sg.wall, &at, &seg_ends)
                    .ok_or_else(undecided)?
                {
                    continue;
                }
                // Which of `wc`'s rulings the crossing is on, `0` being `wc`'s tangent one —
                // asked of the point ([`ruling_side_signed`]), not read off the
                // crossing's own root (`Double`, i.e. the *line* `wc ∩ wall` touches the
                // cylinder), which is the same statement only while `wc` is the wall that
                // touches: a line tangent to the cylinder can cross a `wc` that cuts it, and
                // then the root's `0` matches no ruling of `wc` at all.
                let side = ruling_side_signed(&w, def, (&line, &s)).ok_or_else(undecided)?;
                let Some(c) = coord_at(def, &line, &s) else {
                    return Err(undecided());
                };
                // Which piece of that ruling (several faces of one solid can put several
                // pieces on one side): closed-extent test on the axis coordinate.
                for (ri, r) in rulings.iter().enumerate() {
                    if r.cyl != cyl || r.side != side {
                        continue;
                    }
                    let [lo, hi] = &ruling_ends[ri];
                    let (a, b) = (
                        cmp_along(&c, lo).ok_or_else(undecided)?,
                        cmp_along(&c, hi).ok_or_else(undecided)?,
                    );
                    use core::cmp::Ordering::*;
                    match (a, b) {
                        (Greater, Less) => {
                            // Strictly inside: the ruling splits here, and the segment with it.
                            on_ruling[ri].push(aliases.canon_point(n));
                            on_seg[si].push(aliases.canon_point(n));
                        }
                        (Equal, _) | (_, Equal) => {
                            // On an end: the T-junction. One point, and it must wear **one**
                            // name — the crossing's `{wc, wall}` pair must be the very pair
                            // that named the ruling's end, or the table must know the
                            // two names as one point — else two names share the point.
                            let end_n = if a == Equal { r.end[0] } else { r.end[1] };
                            if aliases.canon_point(end_n) != aliases.canon_point(n) {
                                return Err(reject(RejectReason::CoincidentNodes));
                            }
                            on_seg[si].push(aliases.canon_point(n));
                        }
                        _ => {} // beyond this piece: the surface continues, the face does not
                    }
                }
            }
        }
    }
    if on_ruling.iter().all(|v| v.is_empty()) && on_seg.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    // ---- segments → sub-segments, the arc split's own idiom ----
    let out_segs = split_segments_at(jd, cyls, wc, segs, &mut on_seg, aliases)?;
    // ---- rulings → pieces, in axis order ----
    let mut out_rulings = Vec::with_capacity(rulings.len());
    for (ri, r) in rulings.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_ruling[ri]);
        if nodes.is_empty() {
            out_rulings.push(r);
            continue;
        }
        nodes.sort_unstable();
        nodes.dedup();
        let mut keyed: Vec<(QuadVal, combinatorics::NodeId)> = Vec::new();
        for &n in &nodes {
            let (line, s) =
                combinatorics::pierce_meet(jd, r.cyl, &r.def, n).ok_or_else(undecided)?;
            keyed.push((coord_at(&r.def, &line, &s).ok_or_else(undecided)?, n));
        }
        let [lo, hi] = ruling_ends[ri];
        keyed.push((lo, r.end[0]));
        keyed.push((hi, r.end[1]));
        let mut order: Vec<usize> = (0..keyed.len()).collect();
        let mut bad = false;
        order.sort_by(|&i, &j| match cmp_along(&keyed[i].0, &keyed[j].0) {
            Some(o) => o,
            None => {
                bad = true;
                core::cmp::Ordering::Equal
            }
        });
        if bad {
            return Err(undecided());
        }
        for pair in order.windows(2) {
            if cmp_along(&keyed[pair[0]].0, &keyed[pair[1]].0) == Some(core::cmp::Ordering::Equal) {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        }
        // Ascending axis coordinate **is** the stored convention (`end[0] → end[1]` ascends
        // `+m`), so each window is already a well-formed piece.
        for pair in order.windows(2) {
            out_rulings.push(MergedRuling {
                cyl: r.cyl,
                def: r.def.clone(),
                side: r.side,
                end: [keyed[pair[0]].1, keyed[pair[1]].1],
                merged: r.merged.clone(),
                #[cfg(test)]
                orient: r.orient,
            });
        }
    }
    Ok(Some((out_segs, out_rulings)))
}
