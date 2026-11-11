use super::*;
/// **A split point on one wall's line, in sixteen bytes.**
///
/// The line is `wc ∩ w`, so a plane split point is named by the third plane alone and its full
/// [`NodeId`] is derived by [`Split::name`] where it is needed — the same trade
/// [`combinatorics::OnLine::Class`] makes, and for the same measured reason: these vectors are
/// rebuilt once per wall and a rotated fold spends most of itself allocating. A pierce split point
/// carries its own name, pair and all, narrowed to fit the same size.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(super) enum Split {
    /// First, so `Ord` still runs ascending in the class — the tie groups below hand
    /// `group.first()` to the output as the point's name, and that is the "smallest name wins"
    /// rule the arrangement replays on.
    Class(usize),
    /// The pierce name itself, its pair included: the pair need not be the line's. A pencil — three
    /// or more classes holding the line — names one crossing by any pair of its classes, and the
    /// root is along *that* pair's direction, which can be the line's reversed; restating the name
    /// on the line's pair would keep the root and so name the other crossing. The indices are
    /// narrowed so the whole name fits where a class does (the size lock below).
    Pierce {
        planes: [u32; 2],
        cyl: u32,
        root: crate::combinatorics::QuadRoot,
    },
}

// **Sixteen bytes, as a plane split point is.** These vectors are rebuilt per wall, and widening
// the element to a full `NodeId` (32 bytes) cost the collecting loop 36% on a rotated fold — see
// `split_at_crossings`.
const _: () = assert!(std::mem::size_of::<Split>() == 16);

/// A class or cylinder index, narrowed for [`Split::Pierce`]; a table of `u32::MAX` classes is
/// out of every reach.
fn narrow(i: usize) -> u32 {
    u32::try_from(i).expect("a class index fits in u32")
}

impl Split {
    /// **The pierce arm keeps the name the producer wrote.** A segment on `w` in class `wc` has a
    /// pierce end pinned by the cylinder where both planes of the name's pair hold the line `wc ∩ w`
    /// ([`combinatorics::pin_for`] — a corner whose pair is another line is pinned by its other
    /// plane and arrives here as [`Split::Class`]), and the pair is `{wc, w}` itself or another pair
    /// of the classes holding it; [`Split::name`] puts the name back verbatim.
    /// **A cylinder pin arrives with a pierce name, and that is a producer's invariant, so
    /// this panics rather than naming a refusal.** Every site that writes an `EndPin::Cylinder`
    /// writes the `NodeId::Pierce` beside it — the seated walk, the arc split, the chord pass — so
    /// the two disagreeing is a defect in *this* kernel, not a property of the model. The
    /// alternative would be `PierceVertexUnnamed`, whose sentence is the mirror image of this case:
    /// *"the point is exactly named; what is missing is that these paths have no other name to
    /// carry it by."* Here the point is **not** exactly named — its two halves contradict. A false
    /// sentence in a reject is worse than a loud stop, which is what `ClassIx::plane` says one
    /// door over: *"a loud panic beats a silently wrong plane."*
    /// Measured unexercised over the whole suite and the ignored sweep.
    pub(super) fn of(name: NodeId, pin: combinatorics::EndPin) -> Split {
        match (name.kind(), pin) {
            (_, combinatorics::EndPin::Class(r)) => Split::Class(r),
            (NodeKind::Pierce { planes, cyl, root }, combinatorics::EndPin::Cylinder) => {
                Split::Pierce {
                    planes: planes.map(narrow),
                    cyl: narrow(cyl),
                    root,
                }
            }
            (n, combinatorics::EndPin::Cylinder) => unreachable!(
                "a cylinder pin was written beside a three-plane name: {n:?} — the producer that \
                 made this segment set its two halves apart"
            ),
        }
    }

    /// The point in the form the ordering rule takes, on the line `p ∩ q` it was collected for.
    /// A plane split point costs nothing here — its name *is* its pin, and the rule reads the pin.
    pub(super) fn on(self, p: usize, q: usize) -> combinatorics::PointOn {
        match self {
            Split::Class(r) => combinatorics::PointOn::Class(r),
            Split::Pierce { .. } => combinatorics::PointOn::Pierce(self.name(p, q)),
        }
    }

    /// The name this point ships under.
    fn name(self, p: usize, q: usize) -> NodeId {
        match self {
            Split::Class(r) => NodeId::three_planes(Canon3::three([p, q, r])),
            // The pair came out of a name, so it is sorted and `QuadRoot::canonical` finds nothing
            // to do: the name goes back as it came.
            Split::Pierce { planes, cyl, root } => {
                let [p, q] = planes.map(|i| i as usize);
                NodeId::pierce(p, q, cyl as usize, root)
            }
        }
    }
}

/// **The walls of one class, in first-appearance order, each with its direction family** and the
/// segments riding it — what [`split_at_crossings`] partitions before it asks anything. See
/// [`Wall`].
///
/// **The family replaces a predicate per wall pair, and that is `|W|²` of them.** Two walls'
/// lines on `wc` meet in a point iff `plane_pair_dir_sign(wc, w, r) != 0`, which is the sign of
/// `det(n_wc, n_w, n_r)` — zero exactly when `n_wc × n_w` lies in `r`, i.e. when `wc ∩ w` is
/// **parallel** to `wc ∩ r`. Parallelism of two lines in one plane is an **equivalence relation**,
/// so the walls partition, and each wall needs only to find its family: ask the representatives
/// until one matches. `|W| × families` questions instead of `|W|²`, and the collector then asks
/// none at all.
///
/// The relation needs `wc ∩ w` to *be* a line. It is: a wall comes from a segment that rides it,
/// so the meet carries a segment by construction. (Were `n_w ∥ n_wc` the determinant would vanish
/// against every `r` and the partition would collapse to one family.)
///
/// **Not union-find**, though `Aliases` next door is one. That structure merges *given* pairs;
/// the whole point here is to **not ask** the pairs — transitivity is what lets one question per
/// family stand in for all of them.
///
/// It takes the walls the segments ride rather than the segments: the family is a fact about
/// `s.wall` alone, and the two callers hold different segment types (`MergedSeg` here, `Seg` in
/// the audit that checks this partition against the predicate pair by pair).
pub(super) fn direction_partition(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    seg_walls: impl IntoIterator<Item = usize>,
) -> Vec<Wall> {
    watch!(S_PART);
    let mut walls: Vec<Wall> = Vec::new();
    // Class → its slot in `walls`, for construction only. Nothing downstream reads it: a wall's
    // segments and family travel with the wall, so the loops have one index space.
    let mut slot: HashMap<usize, usize> = HashMap::new();
    // One representative class per direction family, in first-appearance order.
    let mut reps: Vec<usize> = Vec::new();
    for (i, wall) in seg_walls.into_iter().enumerate() {
        let at = *slot.entry(wall).or_insert_with(|| {
            let dir = reps
                .iter()
                .position(|&rep| combinatorics::parallel_carriers(jd, wc, rep, wall))
                .unwrap_or_else(|| {
                    reps.push(wall);
                    reps.len() - 1
                });
            walls.push(Wall {
                class: wall,
                dir,
                segs: Vec::new(),
            });
            walls.len() - 1
        });
        walls[at].segs.push(i);
    }
    walls
}

/// Split the segments on plane class `wc` into the arrangement's 1-skeleton by a **per-wall 1D
/// interval overlay**: on each wall `W`, the line `wc ∩ W` is cut at the sorted union of every
/// W-segment's endpoints and every different-wall crossing, and one `MergedSeg` is emitted per
/// non-empty sub-interval carrying the **union** of the contributions of the W-segments whose
/// closed extent covers it. A crossing of two segments riding `W`/`o.wall` is the plane triple
/// `{wc, W, o.wall}` — no new point species.
///
/// This resolves same-wall partial overlap: where two segments overlap, the shared sub-interval
/// carries both their contributions, which the label brick reads per solid. Sub-intervals run
/// between **distinct** points — names the order ties are kept once and recorded as one point — so
/// no piece of no length is emitted, and every emitted segment's two ends are two points of its
/// line, the contract the splits after this one read. A segment that covers no sub-interval has its
/// two ends at one place and is refused (`CoincidentNodes`). A clean arrangement — no same-wall
/// overlap — has exactly one segment covering each sub-interval, so the output is identical to the
/// naive per-segment split. No re-merge is needed: each `(wall, sub-interval)` is emitted once with
/// its full union, so no two outputs can coincide.
///
/// `Err` (honest reject) when an end or a crossing cannot be located on its line or ordered there
/// (named by [`combinatorics::undecided_reason`] — the gate's or the arithmetic's; locating a
/// piece's ends or asking whether a segment covers them, `WitnessNotRational`), when the alias
/// table refuses to join the names a tie or a line's plane family puts at one point (two distinct
/// plane classes coincident on a wall's line is a four-plane concurrency `{wc, W, a, b}` the
/// 3-plane DCEL cannot name, `FourPlane`), when a piece's end has no pin on the line
/// (`RingNaming`), and on a segment of no length (`CoincidentNodes`).
pub(super) fn split_at_crossings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    aliases: &mut Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    // **The direction sign of each endpoint, once per segment.** `order_along(wc, wall, i, j)`
    // factors into `orient3d(wc, wall, i, j) × dir_sign(wc, wall, j)`, and the second factor does
    // **not** mention `i` — for a segment's endpoints it is a property of the segment alone. The
    // containment test below sweeps `i` over every wall, so asking it there would ask the same
    // question `|walls|` times over (measured: 1.5M `dir_sign` calls where 113k are distinct).
    // **This pass is not plane-only.** A split point is
    // a [`Split`] — a plane class **or** which root of which cylinder — and the order is one rule
    // that takes both ([`combinatorics::order_on`]; here [`combinatorics::closed_contains`] reads
    // [`combinatorics::order_located`] with the hoists handed in). The endpoints are not read as
    // class ids up front either: they are located
    // as [`combinatorics::OnLine`]s, so an end a cylinder pinned
    // has a form to be compared by, and nothing here refuses: [`Split::of`] and
    // `combinatorics::PointOn::of` **panic** on a cylinder pin beside a three-plane name, because
    // those two halves disagreeing is this kernel's defect and not a shape a model can have.
    // **Both ends of every segment, located on that segment's own line, once.** A
    // [`combinatorics::Located`] carries each endpoint's pin and its `dir_sign`, built for the
    // pair `(wc, s.wall)` it will be asked
    // about. Building it here rather than per question is the hoist the containment loop below
    // rests on — measured, that loop runs 1.4–3.2M times per 60-fin fold.
    // **What a failure here is called depends on whose description failed.** A class the gate
    // passed as silent about the cylinders has no world statement, and a question this pass cannot
    // answer without one refuses by the gate's name — the cause is the gate's; anything else is
    // the arithmetic's.
    let cause = |classes: &[usize]| reject(combinatorics::undecided_reason(cyls, classes));
    let end_l: Vec<[combinatorics::OnLine; 2]> = segs
        .iter()
        .map(|s| {
            let mut it = (0..2).map(|k| {
                let on = Split::of(s.end[k], s.end_h[k]).on(wc, s.wall);
                combinatorics::on_line(jd, cyls, wc, s.wall, on).ok_or_else(|| cause(&[s.wall]))
            });
            Ok([it.next().unwrap()?, it.next().unwrap()?])
        })
        .collect::<Result<_, BoolError>>()?;

    // The extent question is [`combinatorics::closed_contains`]; what stays here is *which* pair to
    // ask it for.
    // **The line is `wc ∩ segs[si].wall`, not `wc ∩` the wall being processed.** The
    // collector below asks *other* walls' segments about this wall's crossing, and each of those
    // segments rides its own line — which is also the line `end_l[si]` was located on, and which
    // the rule requires the pair to be.
    // `q` is a **call-site constant**, not a per-segment lookup: the collector asks one
    // wall's crossing of every segment on **one other wall**, and the cover loop asks about
    // segments on **this** wall. Reading `segs[si].wall` at the question instead cost a touch of a
    // wide `MergedSeg` on a loop that runs millions of times — measured, and this is the fix.
    // The alias table's representative comes with each question, not with the closure: the loops
    // below add to the table between questions, and a containment reads what is known by then.
    let closed_contains = |si: usize,
                           at: &combinatorics::Located<'_>,
                           q: usize,
                           canon: &dyn Fn(NodeId) -> NodeId|
     -> Option<bool> {
        combinatorics::closed_contains(jd, wc, q, at, &end_l[si], &canon)
    };

    let walls: Vec<Wall> = direction_partition(jd, wc, segs.iter().map(|s| s.wall));

    #[cfg(test)]
    {
        phase::scale::add(&phase::scale::SEGS, segs.len());
        phase::scale::add(&phase::scale::WALLS, walls.len());
    }

    let mut out = Vec::new();
    for wall in &walls {
        let w = wall.class;
        // Split points on W's line: every W-segment endpoint, plus every real different-wall
        // crossing (a segment on `o.wall` whose closed extent reaches W's line).
        // **A point, not a plane class.** A `usize` — the third plane naming the point —
        // cannot say "a cylinder pins this one", and a `Vec<usize>` silently drops such an
        // endpoint: the two pieces either side of a boss then cover no interval and vanish,
        // leaving the class's boundary open (measured).
        //
        // **`Split`, not a `NodeId` — sixteen bytes.** On *this* line a plane point is a third
        // plane, `{wc, w}` common to every one of them, and storing the pair 114k–134k times is
        // waste. Carrying the full name instead cost the collecting loop 36% on a rotated fold
        // (75.8→103.5ms, same trip counts): these vectors are rebuilt per wall, and allocation is
        // where that fixture spends most of its time. The same lesson as
        // [`combinatorics::OnLine::Class`] dropping its `name`, one file over. A pierce point keeps
        // its pair — a pencil's need not be `{wc, w}` — in narrowed indices, so the element does
        // not grow.
        let mut pts: Vec<Split> = Vec::new();
        {
            watch!(COLLECT);
            #[cfg(test)]
            phase::scale::add(&phase::scale::COLLECT_TRIPS, segs.len());
            // A segment's own endpoints are split points on its line — whatever pins them.
            for &i in &wall.segs {
                for k in 0..2 {
                    pts.push(Split::of(segs[i].end[k], segs[i].end_h[k]));
                }
            }
            // **Wall-major, because the question is about walls.** What lands in `pts` is a *wall*
            // — the class naming the crossing — so sweeping every segment would ask "is `w`
            // parallel to this segment's wall?" once per segment when the answer depends only on
            // the pair. And once one segment on `r` reaches `w`'s line, the rest cannot add
            // anything: `r` is already a split point.
            for other in &walls {
                // Same family ⇒ the two lines are parallel and meet in no point. This also
                // subsumes `other.class == w`: a wall is in its own family, so one test does the
                // work of an identity check and a predicate call.
                if other.dir == wall.dir {
                    continue;
                }
                let r = other.class;
                // One handle per wall pair — `r`'s segments all ask about this same crossing,
                // and they all ride `r`'s line, which is the pair this point is located for.
                let at = combinatorics::locate(jd, cyls, wc, r, combinatorics::PointOn::Class(w))
                    .ok_or_else(|| cause(&[w, r]))?;
                let mut reaches = false;
                for &i in &other.segs {
                    match closed_contains(i, &at, r, &|n| aliases.canon_point(n)) {
                        Some(true) => {
                            reaches = true;
                            break;
                        }
                        Some(false) => {}
                        // A pierce end against a class the gate passed as silent is ordered by the
                        // side it recorded ([`combinatorics::closed_contains`]), so what is left is
                        // a description that could not be formed.
                        None => return Err(cause(&[w, r])),
                    }
                }
                if reaches {
                    // A crossing of two plane lines is a three-plane point; only a segment's own
                    // endpoint can be one a cylinder pins.
                    pts.push(Split::Class(r));
                }
            }
        }
        let pts = {
            watch!(SORT);
            // Distinct names (same name = same point), then ordered along the line.
            // Sorting by `NodeId` orders by `r` — `{wc, w, r}`'s sorted triple is monotone in
            // `r` — so a group of tied points hands its **smallest class** to
            // `group.first()` below, which is the representative and so the emitted name.
            pts.sort_unstable();
            pts.dedup();
            // An order this line cannot form is named by the classes on it — which the gate's,
            // which the arithmetic's ([`combinatorics::undecided_reason`]) — read off the points.
            let undecided = |pts: &[Split]| {
                let classes: Vec<usize> = std::iter::once(w)
                    .chain(pts.iter().filter_map(|x| match x {
                        Split::Class(r) => Some(*r),
                        _ => None,
                    }))
                    .collect();
                cause(&classes)
            };
            // A comparison that cannot be formed stops the sort rather than reading as a tie
            // ([`combinatorics::try_sort_by`]); a stopped sort still holds every point.
            let sorted = combinatorics::try_sort_by(&mut pts, |&x, &y| {
                combinatorics::order_on(jd, cyls, wc, w, x.on(wc, w), y.on(wc, w))
                    .map(|o| o.cmp(&0))
                    .ok_or(())
            });
            if sorted.is_err() {
                return Err(undecided(&pts));
            }
            // Two DISTINCT points ordering equal are one point wearing two names — the order is
            // exact where it says `0`. Record it, and keep one representative as a split point:
            // splitting at both would emit a zero-length piece between them. A tie of plane names is
            // a concurrency `{wc, w, ·, ·}` and records as one; a tie with a cylinder-pinned point
            // records every name the point has ([`Aliases::record_same_point`]) — this is a round of
            // the fixed-point loop, and the next round reads the point under one name, which the
            // pieces below ship (`sorted`). Two pierce names on two cylinders are where two
            // cylinders meet, a point no road names, and are refused.
            let mut reps: Vec<Split> = Vec::with_capacity(pts.len());
            let mut group: Vec<Split> = Vec::new();
            let mut refused: Option<RejectReason> = None;
            let flush = |group: &mut Vec<Split>,
                         reps: &mut Vec<Split>,
                         al: &mut Aliases,
                         refused: &mut Option<RejectReason>| {
                if let Some(&rep) = group.first() {
                    reps.push(rep);
                }
                if group.len() > 1 {
                    let mut set: Vec<usize> = vec![wc, w];
                    set.extend(group.iter().filter_map(|s| match s {
                        Split::Class(r) => Some(*r),
                        Split::Pierce { .. } => None,
                    }));
                    set.sort_unstable();
                    set.dedup();
                    if let Err(e) = al.record(jd, &set) {
                        refused.get_or_insert(e.into());
                    }
                    let first = group[0].name(wc, w);
                    for g in &group[1..] {
                        if let Err(r) = al.record_same_point(jd, cyls, first, g.name(wc, w)) {
                            refused.get_or_insert(r);
                        }
                    }
                }
                group.clear();
            };
            for &r in &pts {
                // The group's first point against this one is not always a pair the sort asked,
                // so it can fail here too — and a failure is no answer, not «another point»: read
                // as one, it would split a group of tied names in two.
                let same = match group.first() {
                    None => false,
                    Some(&g) => {
                        match combinatorics::order_on(jd, cyls, wc, w, g.on(wc, w), r.on(wc, w)) {
                            Some(o) => o == 0,
                            None => return Err(undecided(&pts)),
                        }
                    }
                };
                if !same {
                    flush(&mut group, &mut reps, aliases, &mut refused);
                }
                group.push(r);
            }
            flush(&mut group, &mut reps, aliases, &mut refused);
            if let Some(r) = refused {
                return Err(reject(r));
            }
            reps
        };
        // **The fourth plane that *carries* this line, rather than crossing it.**
        //
        // The report above sees a concurrency when two classes **cross** `W`'s line at one point:
        // they order equal, so one point wears two handles. A plane that *contains* the line cannot
        // be seen that way — it is in `W`'s own direction family, which the collector skips
        // (`other.dir == wall.dir`), so it never becomes a handle and never orders against anything.
        //
        // But `Aliases::canon_wall` already knows it: folding two walls onto one line is exactly the
        // statement that a second plane carries it. Every split point on this line therefore lies on
        // every plane of the family, and that is a four-plane concurrency the point fold must hear
        // about — otherwise one class names the point with the canonical wall and another names it
        // with a plane that cuts the line, and the two names reach the seam table unlinked
        // (`SeamAlias`, measured: `[3,11,12]` vs `[3,6,11]` at one point).
        //
        // The report happens on the class whose wall was **folded**, and the class that produces
        // the other name reads the alias afterwards — the fixpoint loop above runs another round
        // whenever this grows the table, which is what carries it across.
        {
            let family = aliases.wall_family(wc, w);
            if family.len() > 1 {
                // Only a point three planes name joins a plane-class fold — a pierce name has
                // no third class to contribute to the set.
                for r in pts.iter().filter_map(|s| match s {
                    Split::Class(r) => Some(*r),
                    Split::Pierce { .. } => None,
                }) {
                    let mut set: Vec<usize> = vec![wc, r];
                    set.extend(family.iter().copied());
                    set.sort_unstable();
                    set.dedup();
                    aliases.record(jd, &set).map_err(|e| reject(e.into()))?;
                }
            }
        }
        // A split point carries its own name — folded, because a point rebuilt from a handle would
        // otherwise re-introduce the very alias `merge_coincident` just removed. Canonicalizing
        // here and in the merge means every name **downstream** is already the canonical one, and
        // no later stage has to know the table exists.
        // The **raw** name is what `dedup` and the tie groups above key on, and the canonical one
        // is what ships: fold first and a four-plane concurrency would collapse before it is
        // recorded, and the symptom surfaces much later as `SeamAlias`.
        let sorted = |s: Split, al: &Aliases| al.canon_point(s.name(wc, w));
        // Each sub-interval [p, q] carries the union of the W-segments that cover it. Every segment
        // endpoint is itself a split point, so "covers both ends" means "spans the whole interval".
        watch!(COVER);
        #[cfg(test)]
        {
            phase::scale::add(&phase::scale::PTS, pts.len());
            phase::scale::add(
                &phase::scale::COVER_TRIPS,
                pts.len().saturating_sub(1) * wall.segs.len(),
            );
        }
        let mut covers = vec![false; wall.segs.len()];
        for pair in pts.windows(2) {
            let (p, q) = (pair[0], pair[1]);
            let mut merged: Vec<(SolidSide, SegKind)> = Vec::new();
            // The same hoist as the collector: every `w`-segment asks about these two points, and
            // they all ride `w`'s line.
            let loc = |s: Split| {
                combinatorics::locate(jd, cyls, wc, w, s.on(wc, w))
                    .ok_or_else(|| reject(RejectReason::WitnessNotRational))
            };
            let (at_p, at_q) = (loc(p)?, loc(q)?);
            for (x, &i) in wall.segs.iter().enumerate() {
                // **Short-circuiting, deliberately.** `A && B` never asks the second
                // question when the first says no, and this loop runs 1.6M times per fold —
                // evaluating both eagerly doubled it (48.6→89.9ms, same trip counts). It also
                // keeps the *declines* identical: a second question that cannot be formed is never
                // reached when the first already ruled the segment out.
                let undecided = || reject(RejectReason::WitnessNotRational);
                if !closed_contains(i, &at_p, w, &|n| aliases.canon_point(n))
                    .ok_or_else(undecided)?
                {
                    continue;
                }
                if closed_contains(i, &at_q, w, &|n| aliases.canon_point(n))
                    .ok_or_else(undecided)?
                {
                    covers[x] = true;
                    merged.extend(segs[i].merged.iter().copied());
                }
            }
            if !merged.is_empty() {
                // The pin is a fact about the name beside it (`pin_for`): a handle whose
                // point the table represents under another name — a class's root that is an
                // operand's corner — is pinned as that name is on this line.
                let (np, nq) = (sorted(p, aliases), sorted(q, aliases));
                let pin = |n: NodeId| {
                    combinatorics::pin_for(jd, wc, w, n)
                        .ok_or_else(|| reject(RejectReason::RingNaming))
                };
                out.push(MergedSeg {
                    wall: w,
                    end: [np, nq],
                    end_h: [pin(np)?, pin(nq)?],
                    merged,
                    sense: None,
                });
            }
        }
        // **A segment's two ends at one place are an edge of no length, refused here, once.**
        // Both ends are split points and the pieces run between neighbouring points the order
        // tells apart, so a segment covers no piece exactly when its ends are one place — one
        // name, two the table joined (the merge folds them), or two the order above has just tied.
        // Dropped instead, it vanishes without a word: measured, all three spellings left the
        // class's pieces as they are without the segment. The splits after this one read the
        // contract (`split_segments_at`); arcs and rulings ask the same question where their ends
        // are first ordered.
        if covers.contains(&false) {
            return Err(reject(RejectReason::CoincidentNodes));
        }
    }
    Ok(out)
}

/// One face of the arrangement: the cyclic list of half-edges bounding it, and its winding
/// (`+1` a bounded island, `-1` the unbounded outer contour).
#[derive(Clone, Debug)]
pub(crate) struct Cell {
    pub half_edges: Vec<usize>,
    pub winding: i8,
}

/// Number of connected components of the 1-skeleton (union-find over vertex names joined by each
/// piece — segment, arc, or ruling). The `-1`-winding face count must equal this.
pub(super) fn component_count(pieces: &[Piece]) -> usize {
    let mut idx: HashMap<NodeId, usize> = HashMap::new();
    let mut id = |t: NodeId, parent: &mut Vec<usize>| -> usize {
        let n = idx.len();
        *idx.entry(t).or_insert_with(|| {
            parent.push(n);
            n
        })
    };
    let mut parent: Vec<usize> = Vec::new();
    // **The arcs join too.** The `-1` count the walk compares against is a count of the
    // 1-skeleton's components, and an arc is an edge of that skeleton like any other — leaving them
    // out would make a crossed circle look like a component of its own and the walk's check would
    // then be off by one wherever a circle was cut.
    for e in pieces.iter().map(|p| p.end) {
        let (a, b) = (id(e[0], &mut parent), id(e[1], &mut parent));
        let (ra, rb) = (uf_find(&mut parent, a), uf_find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..parent.len())
        .filter(|&i| uf_find(&mut parent, i) == i)
        .count()
}
