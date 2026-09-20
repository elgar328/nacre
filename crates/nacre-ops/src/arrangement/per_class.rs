use super::*;
/// **The per-class stages, in one place** — the walk and the nesting.
///
/// ★★★★ **It exists because the sequence was written twice.** `trace_result_faces`' `arrange` and
/// `frame_audit`'s replay ran the same stages against the same edges, and the second is what makes
/// a reject debuggable — so a drift between them is, in `decline_to_reject`'s words, *"the worst
/// possible time to be lying"*. One copy, and the audit is auditing what the boolean ran.
///
/// ★★★ **And it puts every stage behind one `Result`, which is the failure this shape is really
/// for.** While a stopper occupies the caller's socket, it must run before this `Result` is
/// unwrapped, so an out-of-coverage class carries the same name out however far the pipeline
/// got. With the stages held separately that order had to be repeated per stage, and getting it
/// wrong is silent — the suite stays green and only the reject's *name* changes, which is
/// exactly what happened once (`f8eb935`). One `Result` leaves one place to put the `?`, and it
/// is after the socket.
///
/// ★ **`emit_faces` is the last stage, and it cannot fail.** Every other stage returns a
/// `Result`; this one returns its product outright — so a probe of the socket's interception
/// cannot be planted from outside by failing this stage.
pub(super) fn per_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    kind: BoolKind,
    wc: usize,
    edges: &ClassEdges<'_>,
) -> Result<Staged, BoolError> {
    let (cells, face_of) = timed!(C_EXTRACT, walk_cells(jd, cyls, wc, edges))?;
    let nesting = timed!(C_NEST, nest_cells(jd, cyls, wc, &cells, edges))?;
    // ★ **The seed is `[false; 4]`, and the argument is why it stays an argument.** The
    // arrangement covers all of space, so its unbounded cells reach infinity, where neither solid
    // is. That is a fact about arranging the *whole* model — restrict the input to a region of
    // space and the unbounded cells become an artifact of the restriction, which is what the
    // parameter records.
    let labels = timed!(
        C_LABEL,
        label_cells(&cells, &face_of, edges, &nesting, [false; 4])
    )?;
    let (faces, disk_labels, arc_labels) = timed!(
        C_EMIT,
        emit_faces(kind, &labels, &cells, edges, jd, wc, &nesting.holes)
    );
    // The cylinder chart's **vertical** lines, read off the pieces this class already made. Both
    // ends carry an axis parameter or the piece is not stated: a ruling whose end has no ⊥ partner
    // has no place on a chart's `z` axis, and inventing one would be worse than leaving it out —
    // so it is dropped here and **counted** by the chart's own census rather than guessed at.
    // ★ Walked by index rather than `enumerate`: the half-edge number below is `2·(ns+na) + 2i`,
    // and indexing keeps `i` used in every build (a discarded `enumerate` index is a lint, and an
    // `allow` for it would be a second thing to keep true).
    let ruling_extents: Vec<(usize, RulingExtent)> = (0..edges.rulings.len())
        .filter_map(|i| {
            let r = &edges.rulings[i];
            let z = [
                node_axis_param(jd, &r.def, r.end[0])?,
                node_axis_param(jd, &r.def, r.end[1])?,
            ];
            // ★★ The chart's vertical answer — see [`RulingExtent::label`]. It is an `Option`
            // rather than a `?` deliberately: dropping the piece here would make the recorded
            // population differ between a test build and a release one, and then the census would
            // be measuring a corpus production never sees.
            let label = {
                let base = 2 * (edges.segs.len() + edges.arcs.len());
                let cell_at = |even: bool| -> Option<Label> {
                    Some(labels[*face_of.get(&(base + 2 * i + usize::from(!even)))?])
                };
                // ★ A **tangent** ruling (`side == 0`) has no interior side on this class:
                // the plane touches the cylinder, and both cells beside the line lie outside its
                // disk. Its extent carries no label — the chart reads it as a station where the
                // face ends, not as a wall with a chamber behind it.
                let inside = if r.side == 0 {
                    None
                } else {
                    ruling_interior_is_even(jd, wc, r.side)
                };
                let out = inside.and_then(cell_at);
                #[cfg(test)]
                {
                    // ★ The postcondition, checked rather than assumed — see `ruling_probe::SIDE_CHECK`.
                    // ★★ It asks whether a **cell** carries the lateral's solid, not whether the two
                    // cells *differ*: crossing a ruling on this wall crosses the **lateral**, so they
                    // differ always (☑ measured while designing this rung: every piece, without
                    // exception) and that says nothing about the side.
                    // ★ Recorded only where a label was formed, so `None` means one thing — the content
                    // did not distinguish — and never "there was nothing to check".
                    if let (Some(inside), Some(s)) = (inside, r.merged.first().map(|m| m.0)) {
                        let has = |l: Label| {
                            let b = match s {
                                crate::planes::SolidSide::A => 0,
                                crate::planes::SolidSide::B => 2,
                            };
                            l[b] || l[b + 1]
                        };
                        if let (Some(i_l), Some(e_l)) = (cell_at(inside), cell_at(!inside)) {
                            // ★ The lateral's material is inside the cylinder for a boss and
                            // **outside** for a bore or a notch (`orient` — re-operated Cut
                            // results are such a population on a through-axis class), so
                            // the content agrees with the derivation when the cell on the material's
                            // side has the solid and the other does not.
                            let material_inside = r.orient > 0;
                            let verdict = match (has(i_l), has(e_l)) {
                                (a, b) if a != b => Some(a == material_inside),
                                _ => None,
                            };
                            // ★★★★★ **Asserted where the fact is made**, so the coverage is total and
                            // the panic names the offending test — a test reading the ledger afterwards
                            // sees only what ran before it (this ladder measured that hole twice).
                            assert_ne!(
                                verdict,
                                Some(false),
                                "{}: class {wc}, side {}",
                                ruling_probe::WRONG_SIDE,
                                r.side
                            );
                            ruling_probe::SIDE_CHECK
                                .lock()
                                .expect("the probe's lock is never held across a panic")
                                .push(verdict);
                        }
                    }
                    // A tangent piece (`side == 0`) carries no label by design; the
                    // ledger's proposition is about the rulings a chamber lies behind.
                    if r.side != 0 {
                        ruling_probe::LABELLED
                            .lock()
                            .expect("the probe's lock is never held across a panic")
                            .push(out.is_some());
                    }
                }
                out
            };
            Some((
                r.cyl,
                RulingExtent {
                    wall: wc,
                    side: r.side,
                    end: r.end,
                    z,
                    label,
                    #[cfg(test)]
                    marks: r.merged.clone(),
                },
            ))
        })
        .collect();
    Ok(Staged {
        #[cfg(test)]
        cells,
        #[cfg(test)]
        nesting,
        #[cfg(test)]
        labels,
        faces,
        disk_labels,
        arc_labels,
        ruling_extents,
    })
}

/// **Does this circle actually meet that ruling?** — asked of the class's own tables.
///
/// ★ This is the check `RejectReason::CircleCrossesRuling` is raised from, so it ships.
///
/// ★★ The comparison is exactly
/// [`nacre_scalar::cylinder_ruling_reached`] — a disk of the circle's radius about its centre,
/// against the one named ruling — so the net asks the predicate the gate asks, and there is
/// one derivation rather than two. `None` stays "could not be measured", never "no".
fn circle_crosses_ruling(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    circle: &MergedCircle,
    ruling: &MergedRuling,
) -> Option<bool> {
    let coeffs = combinatorics::class_coeffs_rat(jd, wc)?;
    let centre = combinatorics::circle_centre_rat(jd, wc, &circle.def)?;
    let (o, m, r2) = (ruling.def.origin(), ruling.def.dir(), ruling.def.r2());
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    // The predicate's own precondition, asked rather than assumed: a class that is not parallel
    // to the axis carries no ruling of it, and the decomposition would be about other geometry.
    if nacre_scalar::dot_sign_rat(&n, &m) != nacre_scalar::Orient::Zero {
        return None;
    }
    // The strip runs along `e = n × m`, so that is the direction a piece states its reach in.
    let e = combinatorics::cross3_rat(&n, &m)?;
    let ask = |arc: Option<&crate::planes::RimArc>| -> Option<bool> {
        let (lo, hi) =
            crate::planes::arc_ends_along(&centre, circle.def.r2(), &circle.def.dir(), arc, &e)?;
        Some(nacre_scalar::cylinder_ruling_reached_extent(
            &coeffs,
            &nacre_scalar::StripReach {
                lo: (&nacre_scalar::MeetPoint::Narrow(lo.0), &lo.1),
                hi: Some((&nacre_scalar::MeetPoint::Narrow(hi.0), &hi.1)),
            },
            &o,
            &m,
            r2,
            // ★ **The two vocabularies are opposite.** `ruling.side` is measured against `m × n`
            // (`ruling_side_signed`), and the scalar family writes its sides about `n × m`. The
            // restatement is one negation, and it belongs here — at the boundary between the two
            // index spaces — rather than inside a predicate that would then have to guess.
            -ruling.side,
            // ★★★ **A touch is not what this asks about**. The node no road mints is a
            // **crossing**; an arc tangent to a ruling divides nothing, which is the same sentence
            // as for a tangency at a vertex. Measured: a slab wall at `x = ±12` puts a
            // `d 30` cylinder's rulings exactly on a `d 18` tool's circle, and the boolean that
            // follows validates clean.
            false,
        ))
    };
    // ★★★★★ **The question is about the arc, not the circle**. A contribution states the
    // angular extent its face covers, and *that* is the edge — the rest of the circle is a
    // continuation no face uses. Asking the whole circle refuses a plate whose corner
    // fillets never come near the drill's rulings.
    //
    // An extent that cannot be realized, and a circle with nothing partial to read, fall back to
    // the whole circle: a superset of every arc on it, so the answer stays sound and is the one
    // this predicate always gave.
    if circle.merged.is_empty() {
        return ask(None);
    }
    for (_, _, extent) in &circle.merged {
        let hit = match extent {
            None => ask(None)?,
            Some(ends) => match rim_arc_of(jd, cyls, &centre, *ends) {
                Some(arc) => ask(Some(&arc))?,
                None => ask(None)?,
            },
        };
        if hit {
            return Some(true);
        }
    }
    Some(false)
}

/// **A contribution's angular extent as the two radial vectors [`crate::planes::RimArc`] speaks**
/// — its ends realized and taken from the circle's centre.
///
/// The order is the contribution's own, which [`CircleTrace::arc`] states is counter-clockwise
/// about the cylinder's axis — the same convention `RimArc` carries, so no reordering. A node with
/// no rational coordinate (a `Wide` meet, a pierce whose root is irrational) gives `None`, and the
/// caller then reads the whole circle rather than guess.
fn rim_arc_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    centre: &[nacre_scalar::Rat; 3],
    ends: [NodeId; 2],
) -> Option<crate::planes::RimArc> {
    let radial = |node: NodeId| -> Option<[nacre_scalar::Rat; 3]> {
        let p = combinatorics::node_coords_rat(jd, node)
            .or_else(|| combinatorics::pierce_coords_rat(jd, cyls, node))?;
        let mut v = p;
        for k in 0..3 {
            v[k] = v[k].checked_sub(centre[k])?;
        }
        Some(v)
    };
    Some(crate::planes::RimArc {
        from: radial(ends[0])?,
        to: radial(ends[1])?,
    })
}

/// **Does the circle lie wholly between the two rulings?** — the sibling population, counted
/// because it is the one shape the walk sees nowhere else (a disk afloat in a
/// ruling-bounded region). Same door, read for its third answer.
#[cfg(debug_assertions)]
fn circle_inside_strip(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    circle: &MergedCircle,
    ruling: &MergedRuling,
) -> Option<bool> {
    let coeffs = combinatorics::class_coeffs_rat(jd, wc)?;
    let centre = combinatorics::circle_centre_rat(jd, wc, &circle.def)?;
    let (o, m, r2) = (ruling.def.origin(), ruling.def.dir(), ruling.def.r2());
    if nacre_scalar::dot_sign_rat(&[coeffs[0], coeffs[1], coeffs[2]], &m)
        != nacre_scalar::Orient::Zero
    {
        return None;
    }
    Some(matches!(
        nacre_scalar::cylinder_strip_side_margin(
            &coeffs,
            &nacre_scalar::MeetPoint::Narrow(centre),
            circle.def.r2(),
            &o,
            &m,
            r2,
        ),
        nacre_scalar::StripSide::Inside
    ))
}

/// How many `(circle, ruling)` pairs the net above looked at, and how many it **could not
/// measure** — the denominator of "zero crossings", without which the zero could mean
/// "the instrument saw nothing".
#[cfg(debug_assertions)]
pub(crate) struct MixedClassAudit {
    pub(crate) pairs: std::sync::atomic::AtomicUsize,
    pub(crate) unmeasured: std::sync::atomic::AtomicUsize,
    /// Circles that lie wholly inside a ruling cylinder's strip — the walk's other first-time
    /// shape, counted so "never exercised" is a number rather than a guess.
    pub(crate) inside: std::sync::atomic::AtomicUsize,
}

#[cfg(debug_assertions)]
pub(crate) static MIXED_CLASS_AUDIT: MixedClassAudit = MixedClassAudit {
    pairs: std::sync::atomic::AtomicUsize::new(0),
    unmeasured: std::sync::atomic::AtomicUsize::new(0),
    inside: std::sync::atomic::AtomicUsize::new(0),
};

impl<'a> ClassEdges<'a> {
    /// Run the arc split — and, when the class carries them, the chord injection and the ruling
    /// split — and hold the result: **the one place a class's edges are assembled.**
    pub(super) fn of(
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        wc: usize,
        segs: &'a [MergedSeg],
        circles: &'a [MergedCircle],
        rulings: &'a [MergedRuling],
        aliases: &Aliases,
    ) -> Result<Self, BoolError> {
        use std::borrow::Cow;
        // ★★★★★ **A class may carry a circle (⊥ one cylinder's axis) and rulings (∥ another's)
        // at once, and that is not a refusal**. What would be unarrangeable is the two
        // **meeting**: neither split cuts the other's kind — [`split_circles`] cuts circles
        // against segments, [`split_rulings`] cuts rulings against segments — so a crossing they
        // both walked past would be a node no road ever mints, and the point itself is
        // `plane ∩ cylinder ∩ cylinder`, a nested radical with no name.
        //
        // ★★★ **That population cannot reach here, by two layers.** A circle only becomes an edge
        // *whole* (a partial trace forces its own cut, and `PartialCircleUncut` names the
        // alternative) and a ruling is a segment of a **face**, so two such edges meeting means
        // two lateral **faces** share a point — exactly what [`crate::planes::lateral_faces_clear`]
        // is written to deny ("the two classes share no face — the proposition the arrangement
        // needs"). A pair not shown to clear is `CylinderPairContact`, raised before any
        // arrangement runs. That covers pairs from **different** solids; the pair loop skips
        // same-solid pairs on purpose. But a same-solid pair cannot be built either: two cylinder
        // faces whose axes are *not* parallel need two extrude directions, hence a boolean, and
        // at that boolean they are on opposite sides and meet the very gate above. (Two overlapping
        // circles in one sketch do make one body with intersecting cylinders — with **parallel**
        // axes, which never put a circle and a ruling on one class.)
        //
        // ★ So the check below is a **net over an argument, not a filter over a population**, and
        // it lived under `debug_assertions` for that reason: shipping it in release would have been
        // a device behind a wall nothing could reach. The day the cylinder-pair refusal opens, this
        // population becomes real — and the obligation to carry a *shipped* check then was written
        // at that refusal, where it would be read.
        //
        // ★★★ **The check is here, and the argument above is why it never fires**.
        // Read about the whole **circle** — which is not an edge of any face
        // when the face uses a quarter of it — it refuses a plate whose corner fillets never come
        // near a crosswise drill. The question is asked of the arc a contribution states,
        // and with that reading the argument holds: an arc reaching a ruling is inside the other
        // cylinder's own reach (a ruling sits at `√(r² − h²) < r` from its axis), so the pair rule
        // speaks first.
        //
        // ★★ **One sentence covers all three uncut pairs.** `split_at_crossings` cuts segments
        // against segments, [`split_circles`] circles against segments, [`split_rulings`] rulings
        // against segments — so circle×circle, circle×ruling and ruling×ruling are never cut by
        // anything. All three are safe for the same reason: **a class's edge is some face's
        // boundary**, so two of them from different cylinders meeting puts two lateral faces on
        // one point. Only this pair carries a check, because only this pair has a cheap exact
        // predicate; the other two rest on the sentence alone, and building detectors for
        // populations that cannot arrive measures nothing.
        //
        // `None` — the decomposition could not be stated — refuses with the same name:
        // honest-reject over silent-wrong, and measured to fire on nothing today (`unmeasured` is
        // 0 across the corpus, the census and the kit).
        if !rulings.is_empty() && !circles.is_empty() {
            for c in circles {
                for ru in rulings {
                    #[cfg(debug_assertions)]
                    {
                        use std::sync::atomic::Ordering;
                        MIXED_CLASS_AUDIT.pairs.fetch_add(1, Ordering::Relaxed);
                        if circle_crosses_ruling(jd, cyls, wc, c, ru).is_none() {
                            MIXED_CLASS_AUDIT.unmeasured.fetch_add(1, Ordering::Relaxed);
                        }
                        if circle_inside_strip(jd, wc, c, ru) == Some(true) {
                            MIXED_CLASS_AUDIT.inside.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    if circle_crosses_ruling(jd, cyls, wc, c, ru) != Some(false) {
                        return Err(reject(RejectReason::CircleCrossesRuling));
                    }
                }
            }
        }
        type AfterCircles<'a> = (
            Cow<'a, [MergedSeg]>,
            Cow<'a, [MergedArc]>,
            Cow<'a, [MergedCircle]>,
            Vec<(usize, CutRim)>,
        );
        let (segs2, arcs, circles2, cut_rims): AfterCircles<'a> =
            match split_circles(jd, cyls, wc, segs, circles, aliases)? {
                Some((s, c, a, r)) => {
                    // ★ The names below rest on this: `split_circles` answers `Some` only when some
                    // circle collected a crossing, and a crossing yields at least one arc (a tangency
                    // is skipped before it is collected). So "was split" and "has arcs" are one fact.
                    debug_assert!(!a.is_empty(), "a split that cut no circle into arcs");
                    (Cow::Owned(s), Cow::Owned(a), Cow::Owned(c), r)
                }
                None => (
                    Cow::Borrowed(segs),
                    Cow::Borrowed(&[][..]),
                    Cow::Borrowed(circles),
                    Vec::new(),
                ),
            };

        let (segs3, rulings2) = match split_rulings(jd, cyls, wc, &segs2, rulings, aliases)? {
            Some((s, r)) => (Cow::Owned(s), Cow::Owned(r)),
            None => (segs2, Cow::Borrowed(rulings)),
        };
        Ok(ClassEdges {
            segs: segs3,
            arcs,
            rulings: rulings2,
            circles: circles2,
            cut_rims,
        })
    }

    /// Where the walk's half-edges end and the circles' pseudo-half-edges begin.
    pub(super) fn he_count(&self) -> usize {
        2 * (self.segs.len() + self.arcs.len() + self.rulings.len())
    }

    /// **The one classifier.** Everything that used to compare against `2 * segs.len()` asks this.
    pub(super) fn kind(&self, he: usize) -> HalfEdgeKind {
        let ns = self.segs.len();
        let na = self.arcs.len();
        if he < 2 * ns {
            HalfEdgeKind::Seg(he / 2)
        } else if he < 2 * (ns + na) {
            HalfEdgeKind::Arc((he - 2 * ns) / 2)
        } else if he < self.he_count() {
            HalfEdgeKind::Ruling((he - 2 * (ns + na)) / 2)
        } else {
            HalfEdgeKind::Circle((he - self.he_count()) / 2)
        }
    }

    /// The vertex a half-edge leaves. `he % 2 == 0` takes `end[0]` — one rule, all ranges.
    pub(super) fn origin(&self, he: usize) -> NodeId {
        match self.kind(he) {
            HalfEdgeKind::Seg(i) => self.segs[i].end[he % 2],
            HalfEdgeKind::Arc(i) => self.arcs[i].end[he % 2],
            HalfEdgeKind::Ruling(i) => self.rulings[i].end[he % 2],
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge has no vertex to leave")
            }
        }
    }

    /// **The one place a `RingEdge` is made from a half-edge.**
    ///
    /// ★★★★ It used to be written twice, and the two spellings **differed**: the walk passed the
    /// sense a split carried onto its sub-segments, and `nest_cells`' ring builder passed
    /// `Carrier::plane(wall)` — that is, `sense: None`. The answers agreed only because the second
    /// never saw a split segment; the day it does, it would drop the one fact the endpoints can no
    /// longer supply. One spelling, so there is nothing to drift.
    pub(super) fn edge_at(&self, he: usize) -> combinatorics::RingEdge {
        let carrier = match self.kind(he) {
            // ★ `MergedArc::end` runs counter-clockwise about the axis, so the even half-edge
            // travels that way and its twin the other. The convention is stated once, at the
            // split; this is the only place it is read.
            HalfEdgeKind::Arc(i) => {
                let a = &self.arcs[i];
                combinatorics::Carrier::Arc(Box::new(combinatorics::ArcCarrier {
                    cyl: a.cyl,
                    def: a.def.clone(),
                    ccw: he % 2 == 0,
                }))
            }
            HalfEdgeKind::Seg(i) => combinatorics::Carrier::Plane {
                wall: self.segs[i].wall,
                sense: self.segs[i].sense.map(|s| if he % 2 == 0 { s } else { -s }),
            },
            // ★ `MergedRuling::end` runs along `+m`, so the even half-edge travels up and its
            // twin down — the straight reading of the arc convention above.
            HalfEdgeKind::Ruling(i) => {
                let r = &self.rulings[i];
                combinatorics::Carrier::Ruling(Box::new(combinatorics::RulingCarrier {
                    cyl: r.cyl,
                    def: r.def.clone(),
                    side: r.side,
                    up: he % 2 == 0,
                }))
            }
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge is not a ring edge")
            }
        };
        // The endpoints as handles on this edge's line, carried by the producer. Not recovered from
        // the names: a canonical name need not mention `wc` or the wall (see
        // `combinatorics::RingEdge`). ★ An arc's ends are pierce points by construction — and a
        // ruling's too — which is what `EndPin::Cylinder` says.
        let (from_h, to_h) = match self.kind(he) {
            HalfEdgeKind::Seg(i) => (self.segs[i].end_h[he % 2], self.segs[i].end_h[1 - he % 2]),
            HalfEdgeKind::Arc(_) | HalfEdgeKind::Ruling(_) | HalfEdgeKind::Circle(_) => (
                combinatorics::EndPin::Cylinder,
                combinatorics::EndPin::Cylinder,
            ),
        };
        combinatorics::RingEdge {
            node: self.origin(he),
            to: self.origin(he ^ 1),
            carrier,
            from_h,
            to_h,
        }
    }
}
