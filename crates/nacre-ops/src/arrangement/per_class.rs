use super::*;
/// **The per-class stages, in one place** — the walk and the nesting.
///
/// **It exists because the sequence would otherwise be written twice.** `trace_result_faces`'
/// `arrange` and `frame_audit`'s replay run the same stages against the same edges, and the second
/// is what makes a reject debuggable — so a drift between them would be, in `decline_to_reject`'s
/// words, *"the worst possible time to be lying"*. One copy, and the audit is auditing what the
/// boolean ran.
///
/// **And it puts every stage behind one `Result`, which is the failure this shape is really
/// for.** While a stopper occupies the caller's socket, it must run before this `Result` is
/// unwrapped, so an out-of-coverage class carries the same name out however far the pipeline
/// got. With the stages held separately that order would be repeated per stage, and getting it
/// wrong is silent — the suite stays green and only the reject's *name* changes. One `Result`
/// leaves one place to put the `?`, and it is after the socket.
///
/// **`emit_faces` is the last stage, and it cannot fail.** Every other stage returns a
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
    // **The seed is `[false; 4]`, and the argument is why it stays an argument.** The
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
    // **A sectioned class's rulings close its own section and nothing else.** The gate proved
    // every face on the class apart from the lateral boxes (`Silence::Sectioned`), so no face this
    // class emits borders one of those rulings — and the lateral chart is not told of them below,
    // so a face that did would have no lateral partner. Checked here, where the stage can refuse.
    let sectioned = |cyl: usize| cyls.get(cyl).is_some_and(|cy| cy.sections(wc));
    if faces.iter().any(|f| {
        f.poly_rings().any(|r| {
            r.walls.iter().any(
                |w| matches!(*w, crate::combinatorics::Wall::Ruling { cyl, .. } if sectioned(cyl)),
            )
        })
    }) {
        return Err(reject(RejectReason::CylinderStagesDisagree));
    }
    // The cylinder chart's **vertical** lines, read off the pieces this class already made. Both
    // ends carry an axis parameter or the piece is not stated: a ruling whose end has no ⊥ partner
    // has no place on a chart's `z` axis, and inventing one would be worse than leaving it out —
    // so it is dropped here and **counted** by the chart's own census rather than guessed at.
    let ruling_extents: Vec<(usize, RulingExtent)> = edges
        .pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            let OnCurve::Ruling { cyl, side, .. } = p.on else {
                return None;
            };
            let def = &cyls[cyl].def;
            // A sectioned class's ruling is not a line of the lateral's chart: no face of the class
            // reaches the lateral, so it splits no cell there (see the check above).
            if sectioned(cyl) {
                return None;
            }
            let z = [
                node_axis_param(jd, def, p.end[0])?,
                node_axis_param(jd, def, p.end[1])?,
            ];
            // The chart's vertical answer — see [`RulingExtent::label`] and
            // [`RulingExtent::outer`], the two cells beside the ruling. Each is an `Option` rather
            // than a `?` deliberately: dropping the piece here would make the recorded population
            // differ between a test build and a release one, and then the census would be measuring
            // a corpus production never sees.
            let (label, outer) = {
                let cell_at = |even: bool| -> Option<Label> {
                    Some(labels[*face_of.get(&(2 * i + usize::from(!even)))?])
                };
                // A **tangent** ruling (`side == 0`) has no interior side on this class:
                // the plane touches the cylinder, and both cells beside the line lie outside its
                // disk. Its extent carries no label — the chart reads it as a station where the
                // face ends, not as a wall with a chamber behind it.
                let inside = if side == 0 {
                    None
                } else {
                    ruling_interior_is_even(jd, wc, side)
                };
                let out = inside.and_then(cell_at);
                let outer = inside.and_then(|inside| cell_at(!inside));
                #[cfg(test)]
                {
                    // The postcondition, checked rather than assumed — see `ruling_probe::SIDE_CHECK`.
                    // It asks whether a **cell** carries the lateral's solid, not whether the two
                    // cells *differ*: crossing a ruling on this wall crosses the **lateral**, so they
                    // differ always (measured: every piece, without
                    // exception) and that says nothing about the side.
                    // Recorded only where a label was formed, so `None` means one thing — the content
                    // did not distinguish — and never "there was nothing to check".
                    if let Some(s) = p.merged.first().map(|m| m.0) {
                        let has = |l: Label| {
                            let b = match s {
                                crate::planes::SolidSide::A => 0,
                                crate::planes::SolidSide::B => 2,
                            };
                            l[b] || l[b + 1]
                        };
                        if let (Some(i_l), Some(e_l)) = (out, outer) {
                            // The lateral's material is inside the cylinder for a boss and
                            // **outside** for a bore or a notch (`orient` — re-operated Cut
                            // results are such a population on a through-axis class), so
                            // the content agrees with the derivation when the cell on the material's
                            // side has the solid and the other does not.
                            let material_inside =
                                matches!(p.on, OnCurve::Ruling { orient, .. } if orient > 0);
                            let verdict = match (has(i_l), has(e_l)) {
                                (a, b) if a != b => Some(a == material_inside),
                                _ => None,
                            };
                            // **Asserted where the fact is made**, so the coverage is total and
                            // the panic names the offending test — a test reading the ledger afterwards
                            // sees only what ran before it.
                            assert_ne!(
                                verdict,
                                Some(false),
                                "{}: class {wc}, side {}",
                                ruling_probe::WRONG_SIDE,
                                side
                            );
                            ruling_probe::SIDE_CHECK.push(verdict);
                        }
                    }
                    // A tangent piece (`side == 0`) carries no label by design; the
                    // ledger's proposition is about the rulings a chamber lies behind.
                    if side != 0 {
                        ruling_probe::LABELLED.push(out.is_some());
                    }
                }
                (out, outer)
            };
            Some((
                cyl,
                RulingExtent {
                    wall: wc,
                    side,
                    end: p.end,
                    z,
                    label,
                    outer,
                    #[cfg(test)]
                    marks: p.merged.clone(),
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
/// This is the check `RejectReason::CircleCrossesRuling` is raised from, so it ships.
///
/// **First, the pair rule's own question, of the two edges as pieces**
/// ([`crate::planes::pieces_clear`]): an arc is the piece `[t, t]` × its arc on its cylinder, a
/// ruling its axis span × the whole circle on its own (one ruling's angle is seldom rational, and
/// the whole circle is wider, never narrower). So the net and the gate ask one predicate, and
/// what it is asked about is the bounded piece — the ruling's segment, not the line it lies on.
/// Of its three witnesses only two can speak here: a circle and a ruling stand on one class only
/// when the two axes are perpendicular, so the parallel cross-section never runs, and the common
/// perpendicular is the strip direction below, where the ruling's whole-circle reach says less
/// than the strip does. What this guards, then, is that the edges the arrangement drew lie in the
/// faces the gate proved apart (a whole circle drawn where a face used a quarter of it is wider
/// than the face).
///
/// Then the strip: [`nacre_exact::cylinder_ruling_reached_extent`] — the arc's reach across the
/// strip, against the one named ruling's line — which reads a tangency as no crossing, as this net
/// must (a disjointness proof is about closed sets and cannot). `None` stays "could not be
/// measured", never "no"; a piece that cannot be stated only skips the first question.
fn circle_crosses_ruling(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    circle: &MergedCircle,
    ruling: &MergedRuling,
) -> Option<bool> {
    let coeffs = jd.world_rat(wc)?;
    let (circle_def, ruling_def) = (&cyls[circle.cyl].def, &cyls[ruling.cyl].def);
    let centre = combinatorics::circle_centre_rat(jd, wc, circle_def)?;
    let (o, m, r2) = (ruling_def.origin(), ruling_def.dir(), ruling_def.r2());
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    // The predicate's own precondition, asked rather than assumed: a class that is not parallel
    // to the axis carries no ruling of it, and the decomposition would be about other geometry.
    if nacre_exact::dot_sign_rat(&n, &m) != nacre_exact::Orient::Zero {
        return None;
    }
    // The ruling as a piece: its two ends' axis parameters (rational — each end's other plane is
    // ⊥ the axis), the whole circle round. The arc's axis parameter is the class's own.
    let ruling_piece = (|| {
        let t0 = node_axis_param(jd, ruling_def, ruling.end[0])?;
        let t1 = node_axis_param(jd, ruling_def, ruling.end[1])?;
        Some(crate::planes::Footprint {
            span: Some([t0.min(t1), t0.max(t1)]),
            theta: crate::planes::Theta::Whole,
        })
    })();
    let circle_t = crate::planes::axis_param_of_plane(&coeffs, circle_def);
    let apart = |theta: crate::planes::Theta| -> bool {
        match (&ruling_piece, circle_t) {
            (Some(rp), Some(t)) => crate::planes::pieces_clear(
                circle_def,
                &crate::planes::Footprint {
                    span: Some([t, t]),
                    theta,
                },
                ruling_def,
                rp,
            ),
            _ => false,
        }
    };
    // The strip runs along `e = n × m`, so that is the direction a piece states its reach in.
    let e = nacre_exact::cross3_rat(&n, &m)?;
    let ask = |theta: &crate::planes::Theta| -> Option<bool> {
        let (lo, hi) =
            crate::planes::arc_ends_along(&centre, circle_def.r2(), &circle_def.dir(), theta, &e)?;
        Some(nacre_exact::cylinder_ruling_reached_extent(
            &coeffs,
            &nacre_exact::StripReach {
                lo: (&nacre_exact::MeetPoint::Narrow(lo.0), &lo.1),
                hi: Some((&nacre_exact::MeetPoint::Narrow(hi.0), &hi.1)),
            },
            &o,
            &m,
            r2,
            // **The two vocabularies are opposite.** `ruling.side` is measured against `m × n`
            // (`ruling_side_signed`), and the scalar family writes its sides about `n × m`. The
            // restatement is one negation, and it belongs here — at the boundary between the two
            // index spaces — rather than inside a predicate that would then have to guess.
            -ruling.side,
            // **A touch is not what this asks about**. The node no road mints is a
            // **crossing**; an arc tangent to a ruling divides nothing, which is the same sentence
            // as for a tangency at a vertex. Measured: a slab wall at `x = ±12` puts a
            // `d 30` cylinder's rulings exactly on a `d 18` tool's circle, and the boolean that
            // follows validates clean.
            false,
        ))
    };
    // **The question is about the arc, not the circle**. A contribution states the
    // angular extent its face covers, and *that* is the edge — the rest of the circle is a
    // continuation no face uses. Asking the whole circle refuses a plate whose corner
    // fillets never come near the drill's rulings.
    //
    // An extent that cannot be realized, and a circle with nothing partial to read, fall back to
    // the whole circle: a superset of every arc on it, so the answer stays sound and is the one
    // this predicate always gave.
    if circle.merged.is_empty() {
        let whole = crate::planes::Theta::Whole;
        return if apart(whole) {
            Some(false)
        } else {
            ask(&whole)
        };
    }
    for (_, _, extent) in &circle.merged {
        let theta = extent.map_or(crate::planes::Theta::Whole, |ends| {
            combinatorics::rim_arc_of(jd, cyls, &centre, ends)
        });
        if apart(theta) {
            continue;
        }
        if ask(&theta)? {
            return Some(true);
        }
    }
    Some(false)
}

/// How many `(circle, ruling)` pairs the net above looked at, and how many it **could not
/// measure** — the denominator of "zero crossings", without which the zero could mean
/// "the instrument saw nothing".
#[cfg(debug_assertions)]
pub(crate) struct MixedClassAudit {
    pub(crate) pairs: std::sync::atomic::AtomicUsize,
    pub(crate) unmeasured: std::sync::atomic::AtomicUsize,
}

#[cfg(debug_assertions)]
pub(crate) static MIXED_CLASS_AUDIT: MixedClassAudit = MixedClassAudit {
    pairs: std::sync::atomic::AtomicUsize::new(0),
    unmeasured: std::sync::atomic::AtomicUsize::new(0),
};

impl<'a> ClassEdges<'a> {
    /// Run the arc split and, when the class carries rulings, the ruling split, and hold the result
    /// as one list of pieces: **the one place a class's edges are assembled.**
    #[allow(clippy::too_many_arguments)]
    pub(super) fn of(
        jd: &Judge<'_, WorkingPlane>,
        cyls: &'a [crate::planes::WorkingCyl],
        faces: &[FaceRow],
        wc: usize,
        segs: &'a [MergedSeg],
        circles: &'a [MergedCircle],
        rulings: &'a [MergedRuling],
        aliases: &Aliases,
    ) -> Result<Self, BoolError> {
        use std::borrow::Cow;
        // **A class may carry a circle (⊥ one cylinder's axis) and rulings (∥ another's)
        // at once, and that is not a refusal**. What would be unarrangeable is the two
        // **meeting**: neither split cuts the other's kind — [`split_circles`] cuts circles
        // against segments, [`split_rulings`] cuts rulings against segments — so a crossing they
        // both walked past would be a node no road ever mints, and the point itself is
        // `plane ∩ cylinder ∩ cylinder`, a nested radical with no name.
        //
        // **That population cannot reach here, by two layers.** A circle only becomes an edge
        // *whole* (a partial trace forces its own cut, and `PartialCircleUncut` names the
        // alternative) and a ruling is a segment of a **face**, so two such edges meeting means
        // two lateral **faces** share a point — exactly what [`crate::planes::lateral_faces_clear`]
        // is written to deny ("the two classes share no face — the proposition the arrangement
        // needs"). A pair not shown to clear is `CylinderPairContact`, raised before any
        // arrangement runs. That covers pairs from **different** solids; the pair loop skips
        // same-solid pairs on purpose. But a same-solid pair cannot be built either: two cylinder
        // faces whose axes are *not* parallel need two extrude directions, hence a boolean, and
        // at that boolean they are on opposite sides and meet the very gate above. (Two arcs of one
        // sketch do make one body whose cylinders meet — with **parallel** axes, which never put a
        // circle and a ruling on one class, and such an operand is refused at the boolean's entry
        // anyway, `CylinderPairEdge`.) A circle and rulings that do not meet make an ordinary
        // class: the cap of a pin nested crosswise in a larger cylinder holds the pin's circle
        // between the larger one's rulings, and the cap a crosswise stud stands proud of holds the
        // stud's rulings inside its circle.
        //
        // **So the check below is a net over an argument, not a filter over a population**,
        // and it ships rather than sitting under `debug_assertions`: the day the cylinder-pair
        // refusal opens, this population becomes real, and the obligation to carry a shipped check
        // then is written at that refusal, where it would be read.
        //
        // **The question is asked of the two edges as pieces** — the arc a contribution states,
        // not the whole circle (which is not an edge of any face when the face uses a quarter of
        // it, and reading it that way refuses a plate whose corner fillets never come near a
        // crosswise drill), and the ruling's segment, not its line (whose line runs past the end
        // of a short crosswise cylinder and through a cap circle its segment never reaches). Asked
        // so, it is the pair rule's own question ([`crate::planes::pieces_clear`]), and the
        // argument holds edge for edge: a pair the gate proved apart has its edges apart.
        //
        // **One sentence covers all three uncut pairs.** Circle×circle, circle×ruling and
        // ruling×ruling are cut by nothing, and all three are safe for the same reason: **a
        // class's edge is some face's boundary**, so two of them from different cylinders meeting
        // puts two lateral faces on one point. Only this pair carries a check, because only this
        // pair has a cheap exact predicate; building detectors for populations that cannot arrive
        // measures nothing.
        //
        // `None` — the decomposition could not be stated — refuses with the same name:
        // honest-reject over silent-wrong, and measured to fire on nothing today (`unmeasured` is
        // 0 across the corpus, the census and the kit).
        if !rulings.is_empty() && !circles.is_empty() {
            for c in circles {
                for ru in rulings {
                    let crosses = circle_crosses_ruling(jd, cyls, wc, c, ru);
                    #[cfg(debug_assertions)]
                    {
                        use std::sync::atomic::Ordering;
                        MIXED_CLASS_AUDIT.pairs.fetch_add(1, Ordering::Relaxed);
                        if crosses.is_none() {
                            MIXED_CLASS_AUDIT.unmeasured.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    if crosses != Some(false) {
                        return Err(reject(RejectReason::CircleCrossesRuling));
                    }
                }
            }
        }
        type AfterCircles<'a> = (
            Cow<'a, [MergedSeg]>,
            Vec<Piece>,
            Cow<'a, [MergedCircle]>,
            Vec<(usize, CutRim)>,
        );
        let (segs2, arcs, circles2, cut_rims): AfterCircles<'a> =
            match split_circles(jd, cyls, wc, segs, circles, aliases)? {
                Some((s, c, a, r)) => {
                    // The names below rest on this: `split_circles` answers `Some` only when some
                    // circle collected a node — an arc's end or a crossing — and either yields at
                    // least one arc (a tangency is skipped before it is collected). So "was split"
                    // and "has arcs" are one fact.
                    debug_assert!(!a.is_empty(), "a split that cut no circle into arcs");
                    (Cow::Owned(s), a, Cow::Owned(c), r)
                }
                None => (
                    Cow::Borrowed(segs),
                    Vec::new(),
                    Cow::Borrowed(circles),
                    Vec::new(),
                ),
            };

        let (segs3, rulings2) = match split_rulings(jd, cyls, faces, wc, &segs2, rulings, aliases)?
        {
            Some((s, r)) => (Cow::Owned(s), Cow::Owned(r)),
            None => (segs2, Cow::Borrowed(rulings)),
        };
        // Segments, then arcs, then ruling pieces — the order [`ClassEdges`] numbers them in.
        let pieces = segs3
            .into_owned()
            .into_iter()
            .map(Piece::from)
            .chain(arcs)
            .chain(rulings2.into_owned().into_iter().map(Piece::from))
            .collect();
        Ok(ClassEdges {
            pieces,
            cyls,
            circles: circles2,
            cut_rims,
            wc,
        })
    }

    /// Where the walk's half-edges end and the circles' pseudo-half-edges begin.
    pub(super) fn he_count(&self) -> usize {
        2 * self.pieces.len()
    }

    /// **The one classifier.** Every reader of a half-edge's kind asks this.
    pub(super) fn kind(&self, he: usize) -> HalfEdgeKind {
        if he < self.he_count() {
            HalfEdgeKind::Piece(he / 2)
        } else {
            HalfEdgeKind::Circle((he - self.he_count()) / 2)
        }
    }

    /// The vertex a half-edge leaves. `he % 2 == 0` takes `end[0]` — one rule, all pieces.
    pub(super) fn origin(&self, he: usize) -> NodeId {
        match self.kind(he) {
            HalfEdgeKind::Piece(i) => self.pieces[i].end[he % 2],
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge has no vertex to leave")
            }
        }
    }

    /// **The one place a `RingEdge` is made from a half-edge.** The carrier, sense included, is the
    /// piece's ([`Piece::carrier`]) — a second spelling that passed `Carrier::plane(wall)`, that is
    /// `sense: None`, would drop the sense a split carries onto its sub-segments, the one fact the
    /// endpoints cannot supply.
    pub(super) fn edge_at(&self, he: usize) -> combinatorics::RingEdge {
        let HalfEdgeKind::Piece(i) = self.kind(he) else {
            unreachable!("a circle's pseudo-half-edge is not a ring edge")
        };
        let p = &self.pieces[i];
        // The endpoints as handles on this edge's line, carried by the producer. Not recovered from
        // the names: a canonical name need not mention `wc` or the wall (see
        // `combinatorics::RingEdge`).
        let pins = p.pins();
        combinatorics::RingEdge {
            node: self.origin(he),
            to: self.origin(he ^ 1),
            carrier: p.carrier(self.cyls, he % 2 == 0, self.wc),
            from_h: pins[he % 2],
            to_h: pins[1 - he % 2],
        }
    }
}
