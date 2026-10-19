use super::*;
/// **A split point on one wall's line, named by what is common to every point there.**
///
/// The line is `wc ∩ w`, so those two planes are the same for all of them and only the third thing
/// differs: a plane class that cuts the line, or which root of which cylinder. The full
/// [`NodeId`] is derived by [`Split::name`] where it is needed — the same trade
/// [`combinatorics::OnLine::Class`] makes, and for the same measured reason: these vectors are
/// rebuilt once per wall and a rotated fold spends most of itself allocating.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(super) enum Split {
    /// First, so `Ord` still runs ascending in the class — the tie groups below hand
    /// `group.first()` to the output as the point's name, and that is the "smallest name wins"
    /// rule the arrangement replays on.
    Class(usize),
    Pierce {
        cyl: usize,
        /// **Already canonical**, because it was read out of a [`NodeId`] rather than solved
        /// here — which is why [`Split::name`] hands the door a *sorted* pair, where
        /// `QuadRoot::canonical` has nothing to do. Passing `(p, q)` in caller order instead
        /// would flip a root that is already right.
        root: nacre_topo::QuadRoot,
    },
}

impl Split {
    /// **The pierce arm keeps the root the producer wrote, not a re-canonicalized one.** A
    /// segment on `w` in class `wc` has a pierce end pinned by the cylinder only where the name's
    /// plane pair *is* `{wc, w}` ([`combinatorics::pin_for`] — a corner whose pair is another line
    /// is pinned by its other plane and arrives here as [`Split::Class`]), so the stored name is
    /// already canonical for it and [`Split::name`] can put it back verbatim. Restating it
    /// through `NodeId::pierce` would re-run the pair ordering and could flip the root.
    /// **A cylinder pin arrives with a pierce name, and that is a producer's invariant, so
    /// this panics rather than naming a refusal.** Every site that writes an `EndPin::Cylinder`
    /// writes the `NodeId::Pierce` beside it — the seated walk, the arc split, the chord pass — so
    /// the two disagreeing is a defect in *this* kernel, not a property of the model. The
    /// alternative was `PierceVertexUnnamed`, and its sentence is the mirror image of this case:
    /// *"the point is exactly named; what is missing is that these paths have no other name to
    /// carry it by."* Here the point is **not** exactly named — its two halves contradict. A false
    /// sentence in a reject is worse than a loud stop, which is what `ClassIx::plane` says one
    /// door over: *"a loud panic beats a silently wrong plane."*
    /// Measured unexercised over the whole suite and the ignored sweep before it was made loud.
    pub(super) fn of(name: NodeId, pin: combinatorics::EndPin) -> Split {
        match (name.kind(), pin) {
            (_, combinatorics::EndPin::Class(r)) => Split::Class(r),
            (NodeKind::Pierce { cyl, root, .. }, combinatorics::EndPin::Cylinder) => {
                Split::Pierce { cyl, root }
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
            // The pair is sorted before the door, not after: `root` came out of a name and is
            // already canonical against the sorted pair, so `QuadRoot::canonical` must find
            // nothing to do. Handing it `(p, q)` in caller order would flip a root that is
            // already right — and the classes do arrive out of order.
            Split::Pierce { cyl, root } => NodeId::pierce(p.min(q), p.max(q), cyl, root),
        }
    }

    /// What pins it — the arrangement ships this beside the name.
    fn pin(self) -> combinatorics::EndPin {
        match self {
            Split::Class(r) => combinatorics::EndPin::Class(r),
            Split::Pierce { .. } => combinatorics::EndPin::Cylinder,
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
    let closed_contains = |si: usize, at: &combinatorics::Located<'_>, q: usize| -> Option<bool> {
        combinatorics::closed_contains(jd, wc, q, at, &end_l[si])
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
        // **`Split`, not a `NodeId` — the name is derived.** On *this* line a point is either
        // a third plane or a root of a cylinder, so `{wc, w}` is common to every one of them and
        // storing it 114k–134k times is waste. Carrying the full name instead cost the collecting
        // loop 36% on a rotated fold (75.8→103.5ms, same trip counts): these vectors are rebuilt
        // per wall, and allocation is where that fixture spends most of its time. The same lesson
        // as [`combinatorics::OnLine::Class`] dropping its `name`, one file over.
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
                    match closed_contains(i, &at, r) {
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
            // Two DISTINCT points ordering equal are one point wearing two names — a four-plane
            // concurrency `{wc, w, ·, ·}`. Record it, and keep one representative as a split point:
            // splitting at both would emit a zero-length piece between them.
            // **A tie that involves a cylinder-pinned point is one the alias table must
            // already know.** The plane fold below is a statement about *plane classes* — it
            // hands `Aliases` a set of them — and a pierce name has none to contribute; nor could
            // a fold made here be trusted, since the DCEL keys vertices by name and the two
            // handles would still ship two names. The table *does* know such a point
            // when it is an operand's corner a class passes through (the seed,
            // `Aliases::record_on_cylinder`): every handle's name canonicalizes to one
            // representative, the pieces below are emitted under it (`sorted`), and one handle
            // is kept as the split point. A tie the table does not know is refused — the honest
            // floor, as `split_segments_at` and `split_circles` refuse the same shape.
            // Measured: the refusal is exercised
            // only with the slab as the first operand, where the corner's own name is not the
            // representative and the tie is between the three-plane name and the class's root.
            let mut reps: Vec<Split> = Vec::with_capacity(pts.len());
            let mut group: Vec<Split> = Vec::new();
            let mut tied_pierce = false;
            let flush = |group: &mut Vec<Split>,
                         reps: &mut Vec<Split>,
                         al: &mut Aliases,
                         tied_pierce: &mut bool| {
                if group.len() > 1 && group.iter().any(|s| matches!(s, Split::Pierce { .. })) {
                    let rep = al.canon_point(group[0].name(wc, w));
                    if group.iter().all(|g| al.canon_point(g.name(wc, w)) == rep) {
                        reps.push(group[0]);
                    } else {
                        *tied_pierce = true;
                    }
                    group.clear();
                    return;
                }
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
                    al.record(jd, &set);
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
                    flush(&mut group, &mut reps, aliases, &mut tied_pierce);
                }
                group.push(r);
            }
            flush(&mut group, &mut reps, aliases, &mut tied_pierce);
            if tied_pierce {
                return Err(reject(RejectReason::CoincidentNodes));
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
        // But `Aliases::wall` already knows it: folding two walls onto one line is exactly the
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
                    aliases.record(jd, &set);
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
            for &i in &wall.segs {
                // **Short-circuiting, deliberately.** `A && B` never asks the second
                // question when the first says no, and this loop runs 1.6M times per fold —
                // evaluating both eagerly doubled it (48.6→89.9ms, same trip counts). It also
                // keeps the *declines* identical: a second question that cannot be formed is never
                // reached when the first already ruled the segment out.
                let undecided = || reject(RejectReason::WitnessNotRational);
                if !closed_contains(i, &at_p, w).ok_or_else(undecided)? {
                    continue;
                }
                if closed_contains(i, &at_q, w).ok_or_else(undecided)? {
                    merged.extend(segs[i].merged.iter().copied());
                }
            }
            if !merged.is_empty() {
                // The pin is a fact about the name beside it (`pin_for`): a handle whose
                // point the table represents under another name — a class's root that is an
                // operand's corner — is pinned as that name is on this line.
                let (np, nq) = (sorted(p, aliases), sorted(q, aliases));
                let pin = |s: Split, n: NodeId| {
                    combinatorics::pin_for(jd, wc, w, n).unwrap_or_else(|| s.pin())
                };
                out.push(MergedSeg {
                    wall: w,
                    end: [np, nq],
                    end_h: [pin(p, np), pin(q, nq)],
                    merged,
                    sense: None,
                });
            }
        }
    }
    Ok(out)
}

/// Keep only the edges that carry news: an edge whose [`edge_mask`] is all-false changes no
/// label when crossed, so the cells on its two sides are the same region — of the arrangement
/// *and* of the result, since equal labels mean an equal keep decision. Left in the skeleton it
/// only damages structure: a knife-edge tangency far from everything else is a dangling edge the
/// face walk cannot close (a two-half-edge orbit), which surfaced as `RingOrientation` two layers
/// from its cause.
///
/// A single-point tangency on a class line becomes no segment at all, and a tangency along a
/// line is likewise an edge that never enters the DCEL.
/// Where a knife edge coincides with the *other* solid's real boundary, that solid's
/// contributions make the mask non-false and the edge stays — a true boundary cannot be dropped
/// by this filter.
///
/// Mask evaluation runs ahead of the walk here, so a conflict (`EdgeOccupancyConflict`) surfaces
/// before any walk-stage reject — nearer its cause.
pub(super) fn drop_newsless(segs: Vec<MergedSeg>) -> Result<Vec<MergedSeg>, BoolError> {
    let mut kept = Vec::with_capacity(segs.len());
    for s in segs {
        if edge_mask(&s.merged)? != [false; 4] {
            kept.push(s);
            continue;
        }
        // The only way to all-false is an even count of same-side grazes per solid (seated flips
        // one bit, a transversal flips two, and an edge has at least one contribution) — anything
        // else dropped here would be a hole in that argument, so it is checked, not assumed.
        debug_assert!(
            s.merged
                .iter()
                .all(|(_, k)| matches!(k, SegKind::Graze { .. })),
            "a newsless edge that is not a graze pair: {:?}",
            s.merged
        );
    }
    Ok(kept)
}

/// One face of the arrangement: the cyclic list of half-edges bounding it, and its winding
/// (`+1` a bounded island, `-1` the unbounded outer contour).
#[derive(Clone, Debug)]
pub(crate) struct Cell {
    pub half_edges: Vec<usize>,
    pub winding: i8,
}

/// Number of connected components of the 1-skeleton (union-find over vertex names joined by each
/// edge — segment, arc, or ruling). The `-1`-winding face count must equal this.
pub(super) fn component_count(
    segs: &[MergedSeg],
    arcs: &[MergedArc],
    rulings: &[MergedRuling],
) -> usize {
    let mut idx: HashMap<NodeId, usize> = HashMap::new();
    let mut id = |t: NodeId, parent: &mut Vec<usize>| -> usize {
        let n = idx.len();
        *idx.entry(t).or_insert_with(|| {
            parent.push(n);
            n
        })
    };
    let mut parent: Vec<usize> = Vec::new();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    // **The arcs join too.** The `-1` count the walk compares against is a count of the
    // 1-skeleton's components, and an arc is an edge of that skeleton like any other — leaving them
    // out would make a crossed circle look like a component of its own and the walk's check would
    // then be off by one wherever a circle was cut.
    let ends = segs
        .iter()
        .map(|s| s.end)
        .chain(arcs.iter().map(|a| a.end))
        .chain(rulings.iter().map(|r| r.end));
    for e in ends {
        let (a, b) = (id(e[0], &mut parent), id(e[1], &mut parent));
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..parent.len())
        .filter(|&i| find(&mut parent, i) == i)
        .count()
}
