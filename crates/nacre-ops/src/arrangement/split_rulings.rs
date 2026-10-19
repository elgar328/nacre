use super::*;
/// **Where a segment crosses a cylinder's rulings on a ∥ class** — the straight sibling of
/// [`circle_crossings`], same machinery ([`nacre_exact::quad::plane_plane_cylinder`], the same
/// fences, the same names), different reachable arms: on a class **parallel** to the axis the
/// meet can be tangent to the lateral — **skipped**, because a touch divides nothing and cutting
/// there would make a zero-length piece — or lie on it, which the ladder refuses by name rather
/// than reasoning it unreachable.
///
/// Returns each crossing with the meet it was solved as ([`combinatorics::PierceMeet`] — exact, or
/// an interval box where a class has no world statement; the side and extent tests below read it
/// without re-solving). `None` is overflow.
fn lateral_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    sg: &MergedSeg,
    on_a_shared_line: bool,
) -> Result<Option<Vec<LateralCrossing>>, BoolError> {
    use nacre_exact::quad::CylinderMeet;
    use nacre_topo::QuadRoot;
    let Some(v) = combinatorics::class_coeffs_rat(jd, sg.wall) else {
        return Ok(None);
    };
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    // **A class with no world statement here sections the cylinder** — it carries rulings, so the
    // gate found it along the axis (`Silence::Sectioned`). Against a wall along the axis too, the
    // line between them is exactly parallel to it (both normals are perpendicular to it) and
    // crosses no ruling: it misses the surface, which `C ≠ 0` shows, or lies on it, which no record
    // states here. Against any other wall the meet's shape is read off the boxes.
    if jd.planes[wc].world.is_none() {
        let prec = jd.standard.prec;
        let undecided = || reject(RejectReason::CylinderGateUndecided);
        let line = nacre_exact::quad::plane_plane_cylinder_hp(
            &jd.planes[wc].canon_hp(prec),
            &jd.planes[sg.wall].canon_hp(prec),
            &o,
            &m,
            r2,
            prec,
        )
        .ok_or_else(undecided)?;
        if nacre_exact::dot_sign_rat(&[v[0], v[1], v[2]], &m) == nacre_exact::Orient::Zero {
            return match line.c.sign() {
                Some(_) => Ok(Some(Vec::new())),
                None => Err(undecided()),
            };
        }
        if line.a.sign() != Some(true) || !line.disc(prec).sign().ok_or_else(undecided)? {
            return match line.a.sign() {
                Some(true) => Ok(Some(Vec::new())), // the discriminant is negative: a miss
                _ => Err(undecided()),
            };
        }
        let mut out = Vec::with_capacity(2);
        for root in [QuadRoot::Lo, QuadRoot::Hi] {
            let n = combinatorics::NodeId::pierce(wc, sg.wall, cyl, root);
            out.push((
                n,
                combinatorics::pierce_meet(jd, cyl, def, n).ok_or_else(undecided)?,
            ));
        }
        return Ok(Some(out));
    }
    let w = combinatorics::class_coeffs_rat(jd, wc)
        .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
    let (line, roots) = match nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2) {
        Some(CylinderMeet::Pair { line, s }) => {
            (line, vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])])
        }
        // A tangent line touches without separating; cutting there would make a zero-length
        // piece — the same skip the circle side's `Double` arm makes.
        // A **tangent** class touches the cylinder along one ruling: the
        // segment's line meets it in one point, root `Double` — the crossing the seated face's own
        // tangent piece (`side == 0`) is cut at, named as the scan names it.
        Some(CylinderMeet::Tangent { line, s }) => (
            line,
            vec![(QuadRoot::Double, nacre_exact::quad::QuadVal::from_rat(s))],
        ),
        Some(CylinderMeet::Miss(_)) | Some(CylinderMeet::AxisParallelMiss(_)) => {
            return Ok(Some(Vec::new()));
        }
        // The segment's line lies on the lateral: a line the gate recorded
        // ([`crate::planes::SharedRuling`]) crosses no ruling — it *is* one, and the fold in
        // [`split_rulings`] states it there. Anywhere else, and where the planes degenerate, the
        // ladder does not arrange the shape yet.
        Some(CylinderMeet::OnRuling(_)) if on_a_shared_line => return Ok(Some(Vec::new())),
        Some(CylinderMeet::OnRuling(_))
        | Some(CylinderMeet::CoincidentPlanes)
        | Some(CylinderMeet::ParallelPlanes) => {
            return Err(reject(RejectReason::RulingBoundNotYet));
        }
        None => return Ok(None),
    };
    // Whether a crossing is **on this segment** is the caller's question now, asked in the one
    // vocabulary that answers it for both kinds of end — the same move the circle side makes.
    Ok(Some(
        roots
            .into_iter()
            .map(|(root, s)| {
                (
                    combinatorics::NodeId::pierce(wc, sg.wall, cyl, root),
                    combinatorics::PierceMeet::Exact(line.clone(), s),
                )
            })
            .collect(),
    ))
}

/// One crossing [`lateral_crossings`] found: its Pierce name, and the meet it was solved as (the
/// side and extent tests read it without re-solving).
type LateralCrossing = (combinatorics::NodeId, combinatorics::PierceMeet);

/// **A point's coordinate along a cylinder's axis** (`m · x`) — exact (`a + b√c`) where the point has
/// an exact meet or lies on a plane across the axis, else an interval. Two exact ones compare
/// exactly ([`nacre_exact::QuadVal::cmp_value`]); a pair with an interval compares by the boxes
/// and is never a tie: `None` where they cannot be separated.
#[derive(Clone)]
enum AxisCoord {
    Exact(nacre_exact::quad::QuadVal),
    Bounded(nacre_exact::HpBounded),
}

impl AxisCoord {
    /// The coordinate of pierce point `n` of cylinder `cyl`.
    fn of(
        jd: &Judge<'_, WorkingPlane>,
        cyl: usize,
        def: &nacre_topo::CylinderDef,
        n: combinatorics::NodeId,
    ) -> Option<Self> {
        Self::of_meet(jd, def, n, &combinatorics::pierce_meet(jd, cyl, def, n)?)
    }

    fn of_meet(
        jd: &Judge<'_, WorkingPlane>,
        def: &nacre_topo::CylinderDef,
        n: combinatorics::NodeId,
        meet: &combinatorics::PierceMeet,
    ) -> Option<Self> {
        use nacre_exact::quad::QuadVal;
        let m = def.dir();
        let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
            x[0].checked_mul(y[0])?
                .checked_add(x[1].checked_mul(y[1])?)?
                .checked_add(x[2].checked_mul(y[2])?)
        };
        if let Some((line, s)) = meet.exact() {
            let (b, d) = (line.base(), line.dir());
            return Some(Self::Exact(
                QuadVal::from_rat(dot(&b, &m)?).checked_add(&s.checked_mul_rat(dot(&d, &m)?)?)?,
            ));
        }
        // A plane of the name across the axis holds every point of it at one axis coordinate,
        // `m·o + t·|m|²` with `t` where the axis meets it — stated exactly, so two points on one
        // cap compare as equal.
        let (planes, _, _) = combinatorics::pierce_name(n)?;
        for c in planes {
            if let Some(w) = combinatorics::class_coeffs_rat(jd, c)
                && nacre_exact::parallel_rat(&[w[0], w[1], w[2]], &m)
            {
                let t = crate::planes::axis_param_of_plane(&w, def)?;
                return Some(Self::Exact(QuadVal::from_rat(
                    dot(&def.origin(), &m)?.checked_add(t.checked_mul(dot(&m, &m)?)?)?,
                )));
            }
        }
        let prec = jd.standard.prec;
        let at = meet.point_hp(jd)?;
        let mh = m.map(|x| nacre_exact::HpBounded::of_rat(x, prec));
        Some(Self::Bounded(
            at[0]
                .mul(&mh[0], prec)
                .add(&at[1].mul(&mh[1], prec), prec)
                .add(&at[2].mul(&mh[2], prec), prec),
        ))
    }

    fn cmp(&self, other: &Self, prec: usize) -> Option<core::cmp::Ordering> {
        use core::cmp::Ordering;
        let hp = |c: &Self| match c {
            Self::Exact(q) => nacre_exact::quad_hp(q, prec),
            Self::Bounded(b) => Some(b.clone()),
        };
        match (self, other) {
            (Self::Exact(a), Self::Exact(b)) => a.cmp_value(b),
            _ => match hp(self)?.sub(&hp(other)?, prec).sign()? {
                true => Some(Ordering::Greater),
                false => Some(Ordering::Less),
            },
        }
    }
}

/// `split_rulings`' product — `None` when nothing crossed.
type SplitRulings = Option<(Vec<MergedSeg>, Vec<MergedRuling>)>;

/// **Cut every ruling a segment crosses into pieces, and the segments with it** — the straight
/// twin of [`split_circles`], run after it.
///
/// **It cuts rulings against *segments*, and deliberately not against arcs** — as its twin cuts
/// circles against segments and not against rulings. One class *can* carry both populations
/// (there is no blanket refusal that says otherwise), so what keeps the pair of splits
/// complete is not «they never share a class» but the sharper fact that **a class's edge is
/// some face's boundary**: two edges of different cylinders meeting would put two lateral faces on
/// one point, which `crate::planes::lateral_faces_clear` denies. That one sentence covers
/// every pair no split cuts — circle×circle and ruling×ruling as well as this one — and
/// `ClassEdges::of` carries a shipped backstop over the pair that has an exact predicate
/// ([`crate::RejectReason::CircleCrossesRuling`]).
///
/// The ordering vocabulary is the established one: along a segment,
/// [`nacre_exact::QuadVal::cmp_value`] on the canonical meet line; along a ruling, the same
/// comparator on the **axis coordinate** (`m · x`, evaluated exactly from each node's re-solved
/// name).
pub(super) fn split_rulings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    rulings: &[MergedRuling],
    aliases: &Aliases,
) -> Result<SplitRulings, BoolError> {
    // A class with no world statement (it sections the cylinder) refuses by the gate's name.
    let undecided = || reject(combinatorics::undecided_reason(cyls, &[wc]));
    if rulings.is_empty() {
        return Ok(None);
    }
    let prec = jd.standard.prec;
    // Each ruling's end coordinates, re-solved from the names once.
    let ruling_ends: Vec<[AxisCoord; 2]> = rulings
        .iter()
        .map(|r| {
            let coord = |n| AxisCoord::of(jd, r.cyl, &r.def, n);
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
    // The lines on this class that lie on a cylinder and that this class states as a ruling —
    // `(the other class, cylinder, which of this class's two rulings)`, the pair read off the
    // gate's record, never decided here, and the side asked of the line once
    // ([`combinatorics::shared_line`], which the straight-angle cleaning reads too): it is the
    // identity a ruling piece carries (`MergedRuling::side`).
    let mut shared_here: Vec<(usize, usize, i8)> = Vec::new();
    for (ci, cy) in cyls.iter().enumerate() {
        for sr in cy.shared.iter().filter(|sr| sr.stated_by.contains(&wc)) {
            let l = combinatorics::shared_line(jd, cyls, ci, sr, wc).ok_or_else(undecided)?;
            shared_here.push((l.other, ci, l.side));
        }
    }
    let mut on_ruling: Vec<Vec<combinatorics::NodeId>> = vec![Vec::new(); rulings.len()];
    let mut on_seg: Vec<Vec<combinatorics::NodeId>> = vec![Vec::new(); segs.len()];
    for (si, sg) in segs.iter().enumerate() {
        // **No skip for a segment whose end a cylinder pins.** The extent question is
        // `closed_contains`, which answers for both kinds of end, and «such a segment shares only
        // endpoints with a ruling» is an argument about one narrow population — the sort of claim
        // that goes quietly false when the population widens — so it is asked instead of assumed.
        // Measured: **200** such segments reach here and **392** crossings on them are examined;
        // none splits a ruling.
        //
        // A segment whose two ends are one point cuts nothing — asked by **name**, the identity,
        // rather than by a coordinate a pierce end does not have. Measured unexercised.
        if sg.end[0] == sg.end[1] {
            continue;
        }
        // A segment on a class the gate proved clear of a cylinder cannot reach its rulings — the
        // gate's record, asked before anything this segment's class would have to describe.
        let here: Vec<(usize, &nacre_topo::CylinderDef)> = cyl_list
            .iter()
            .copied()
            .filter(|&(cyl, _)| !cyls.get(cyl).is_some_and(|cy| cy.is_silent(sg.wall)))
            .collect();
        if here.is_empty() {
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
        for &(cyl, def) in &here {
            let shared = shared_here
                .iter()
                .any(|&(t, c, _)| t == sg.wall && c == cyl);
            let Some(xs) = lateral_crossings(jd, wc, cyl, def, sg, shared)? else {
                return Err(undecided());
            };
            for (n, meet) in xs {
                // On this segment? The caller's question since the fences went — same rule, same
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
                let side = meet.ruling_side(jd, wc, def).ok_or_else(undecided)?;
                let Some(c) = AxisCoord::of_meet(jd, def, n, &meet) else {
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
                        c.cmp(lo, prec).ok_or_else(undecided)?,
                        c.cmp(hi, prec).ok_or_else(undecided)?,
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
    let folds = segs
        .iter()
        .any(|sg| shared_here.iter().any(|&(t, _, _)| t == sg.wall));
    if on_ruling.iter().all(|v| v.is_empty()) && on_seg.iter().all(|v| v.is_empty()) && !folds {
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
        let mut keyed: Vec<(AxisCoord, combinatorics::NodeId)> = Vec::new();
        for &n in &nodes {
            keyed.push((
                AxisCoord::of(jd, r.cyl, &r.def, n).ok_or_else(undecided)?,
                n,
            ));
        }
        let [lo, hi] = ruling_ends[ri].clone();
        keyed.push((lo, r.end[0]));
        keyed.push((hi, r.end[1]));
        let mut order: Vec<usize> = (0..keyed.len()).collect();
        combinatorics::try_sort_by(&mut order, |&i, &j| {
            keyed[i].0.cmp(&keyed[j].0, prec).ok_or(())
        })
        .map_err(|()| undecided())?;
        for pair in order.windows(2) {
            match keyed[pair[0]].0.cmp(&keyed[pair[1]].0, prec) {
                Some(core::cmp::Ordering::Equal) => {
                    return Err(reject(RejectReason::CoincidentNodes));
                }
                Some(_) => {}
                // The sort may have asked this pair the other way round; a failure is no answer.
                None => return Err(undecided()),
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
    let out_segs = if folds {
        fold_onto_rulings(jd, cyls, out_segs, &mut out_rulings, &shared_here, aliases)?
    } else {
        out_segs
    };
    Ok(Some((out_segs, out_rulings)))
}

/// **A segment on a line the gate recorded is that ruling's edge** — so it is stated once, in
/// the ruling vocabulary: its `(solid, kind)` contributions join the ruling piece with the same
/// two ends, and the segment goes. The chart's vertical lines come from rulings alone, and the
/// lateral continues across such a line (the prism's corner is the *other* solid's edge), so the
/// station has to be there; stated as a plane segment as well, one line would reach the walk as
/// two edges.
///
/// The splits above cut both at every point the other crosses, so a piece either matches a
/// ruling piece end for end or shares no stretch with one. No stretch shared is a segment off the
/// lateral — past its cap, or where the ruling is silent because the line is the lateral face's
/// own edge (`ruling_sweep`'s `quiet_nodes`: there the plane vocabulary states it, the other half
/// of the one rule) — and it stays a segment. A stretch shared without a match is a shape the
/// splits did not cut, and is refused by name rather than stated twice.
///
/// `Seated`, `Graze` and `Tangent` say which side of the class plane a face fills, so they mean
/// the same on a ruling. `Transversal` says which side of the *directed line*, and a ruling runs
/// up the axis rather than along `wc × wall`: it is refused by name, not carried across frames.
fn fold_onto_rulings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    segs: Vec<MergedSeg>,
    rulings: &mut [MergedRuling],
    shared_here: &[(usize, usize, i8)],
    aliases: &Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    use nacre_exact::quad::QuadVal;
    let undecided = || reject(RejectReason::WitnessNotRational);
    let canon = |e: [combinatorics::NodeId; 2]| {
        let mut x = [aliases.canon_point(e[0]), aliases.canon_point(e[1])];
        x.sort_unstable();
        x
    };
    // A point's coordinate along the cylinder's axis, exact for either kind of name.
    let axis_at = |cyl: usize, n: combinatorics::NodeId| -> Option<QuadVal> {
        let m = cyls.get(cyl)?.def.dir();
        let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
            x[0].checked_mul(y[0])?
                .checked_add(x[1].checked_mul(y[1])?)?
                .checked_add(x[2].checked_mul(y[2])?)
        };
        if let Some(p) = combinatorics::node_coords_rat(jd, n) {
            return Some(QuadVal::from_rat(dot(&p, &m)?));
        }
        let (_, own, _) = combinatorics::pierce_name(n)?;
        let (line, s) =
            combinatorics::pierce_meet(jd, own, &cyls.get(own)?.def, n)?.into_exact()?;
        QuadVal::from_rat(dot(&line.base(), &m)?)
            .checked_add(&s.checked_mul_rat(dot(&line.dir(), &m)?)?)
    };
    let mut kept = Vec::with_capacity(segs.len());
    for sg in segs {
        let Some(&(_, cyl, side)) = shared_here.iter().find(|&&(t, _, _)| t == sg.wall) else {
            kept.push(sg);
            continue;
        };
        let key = canon(sg.end);
        if let Some(r) = rulings
            .iter_mut()
            .find(|r| r.cyl == cyl && canon(r.end) == key)
        {
            if sg
                .merged
                .iter()
                .any(|(_, k)| matches!(k, SegKind::Transversal { .. }))
            {
                return Err(reject(RejectReason::RulingBoundNotYet));
            }
            r.merged.extend(sg.merged.iter().copied());
            continue;
        }
        // No piece end for end: off the lateral, or a stretch the splits left uncut.
        let (Some(a), Some(b)) = (axis_at(cyl, sg.end[0]), axis_at(cyl, sg.end[1])) else {
            return Err(undecided());
        };
        let (a, b) = match a.cmp_value(&b).ok_or_else(undecided)? {
            core::cmp::Ordering::Greater => (b, a),
            _ => (a, b),
        };
        for r in rulings.iter().filter(|r| r.cyl == cyl && r.side == side) {
            let (Some(lo), Some(hi)) = (axis_at(cyl, r.end[0]), axis_at(cyl, r.end[1])) else {
                return Err(undecided());
            };
            use core::cmp::Ordering::Less;
            if a.cmp_value(&hi).ok_or_else(undecided)? == Less
                && lo.cmp_value(&b).ok_or_else(undecided)? == Less
            {
                return Err(reject(RejectReason::RulingBoundNotYet));
            }
        }
        kept.push(sg);
    }
    Ok(kept)
}
